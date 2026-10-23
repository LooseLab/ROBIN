"""Replay archived completed BAMs through ROBIN's native watcher and Ray queues."""
import argparse
import asyncio
import contextlib
import csv
import hashlib
import json
import math
import os
import time
from pathlib import Path


def read_manifest(path, limit=None):
    with open(path) as handle:
        rows=[dict(bam=str(Path(r['bam']).resolve()),available_at=float(r['available_at']))
              for r in csv.DictReader(handle)]
    rows.sort(key=lambda r:(r['available_at'],r['bam']))
    if limit:rows=rows[:limit]
    if not rows or any(not math.isfinite(r['available_at']) or r['available_at']<0 for r in rows):
        raise ValueError('Invalid or empty replay manifest')
    if len({r['bam'] for r in rows})!=len(rows):raise ValueError('Duplicate BAM input')
    return rows


def publish(path, incoming):
    source=Path(path)
    folder=incoming/source.parent.name
    folder.mkdir(exist_ok=True)
    dest=folder/source.name
    temporary=folder/(source.name+'.pending')
    temporary.symlink_to(source)
    temporary.replace(dest)
    return dest


def processes():
    import psutil
    proc=psutil.Process()
    output=[]
    for p in [proc,*proc.children(recursive=True)]:
        try:
            cpu=p.cpu_times();io=p.io_counters()
            output.append(dict(pid=p.pid,name=p.name(),rss=p.memory_info().rss,
                cpu_seconds=cpu.user+cpu.system,read_bytes=io.read_bytes,write_bytes=io.write_bytes))
        except (psutil.NoSuchProcess,psutil.AccessDenied):pass
    return output


async def simulate(args):
    import ray
    from watchdog.observers import Observer
    from robin import workflow_ray as wr
    from robin.sv_replay.replay import export
    import sqlite3
    import zipfile
    from robin import models
    model=models.DIR/"general.zip"
    if not zipfile.is_zipfile(model):
        raise ValueError("Sturgeon model archive is missing or incomplete")

    root=args.output.resolve()
    if root.exists() and any(root.iterdir()):raise ValueError('Output directory must be new and empty')
    incoming=root/'incoming';incoming.mkdir(parents=True)
    work=root/'robin';work.mkdir()
    rows=read_manifest(args.manifest,args.limit)
    stages=['preprocessing:preprocessing','cnv:cnv','bed_conversion:bed_conversion',
            'classification:sturgeon','classification:nanodx','mgmt:mgmt','fusion:fusion','sv:sv']
    config=dict(manifest=str(Path(args.manifest).resolve()),
        manifest_sha256=hashlib.sha256(Path(args.manifest).read_bytes()).hexdigest(),
        input_batches=len(rows),speed=args.speed,reference=str(args.reference.resolve()),
        cpus=args.cpus,preset='standard',target_panel=args.panel,stages=stages,
        model_assets=[dict(path=str(p),size=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
          for p in [model,models.DIR/'Capper_et_al_NN_v2.pkl']],
        timing_basis='original BAM latest read completion; excludes historical processing/publication delays',
        source_interval_seconds=rows[-1]['available_at'],
        pipeline='native ROBIN Ray coordinator, preprocessing, watcher, batching and handlers; headless',
        published_as='atomic symlinks to immutable archived aligned BAMs',
        source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
         [Path(wr.__file__),Path(__file__),Path(__file__).parent.parent/'analysis/sv_analysis.py']})
    (root/'run.json').write_text(json.dumps(config,indent=2)+'\n')
    temp=f'/tmp/rr{os.getpid()}'
    ray.init(address='local',namespace=f'robin_replay_{os.getpid()}',num_cpus=args.cpus,
             include_dashboard=False,_temp_dir=temp,object_store_memory=1024**3,
             log_to_driver=True)
    coord=wr.Coordinator.options(name='robin_coordinator',num_cpus=0,max_concurrency=32).remote(
        target_panel=args.panel,analysis_workers=1,preprocessing_workers=1,bed_workers=1,
        preset='standard',reference=str(args.reference.resolve()),enable_batching=True)
    observer=None;telemetry=None;last_stats={};drained=False
    try:
        await coord.setup.remote();await coord.set_work_dir.remote(str(work))
        required=[wr._cnv_handler,wr._sturgeon_handler,wr._nanodx_handler,wr._fusion_handler,wr._sv_handler]
        if any(h is None for h in required):raise RuntimeError('Required ROBIN handler unavailable')
        watcher=wr.RayFileWatcher(coord,stages,args.panel,patterns=['*.bam'],
                 work_dir=str(work),reference=str(args.reference.resolve()))
        observer=Observer();observer.schedule(watcher,str(incoming),recursive=True);observer.start()
        start=time.monotonic();epoch=time.time()
        (root/'epoch.json').write_text(json.dumps(dict(wall_start_epoch=epoch))+'\n')
        async def sample():
            nonlocal last_stats
            with (root/'telemetry.jsonl').open('w') as handle:
                tick=0
                while True:
                    last_stats=await coord.stats.remote()
                    proc=processes()
                    record=dict(elapsed=time.monotonic()-start,wall_time=time.time(),
                        process_tree_rss_sum=sum(p['rss'] for p in proc),processes=proc,coordinator=last_stats)
                    if tick % 12 == 0:
                        record['output_disk_bytes']=sum(p.stat().st_size for p in work.rglob('*') if p.is_file() and not p.is_symlink())
                    tick+=1
                    handle.write(json.dumps(record,default=str)+'\n');handle.flush()
                    await asyncio.sleep(5)
        telemetry=asyncio.create_task(sample())
        with (root/'arrivals.jsonl').open('w') as arrivals:
            for i,row in enumerate(rows):
                while time.monotonic()-start<row['available_at']/args.speed:
                    if telemetry.done():await telemetry
                    await asyncio.sleep(min(1,max(0,row['available_at']/args.speed-(time.monotonic()-start))))
                dest=publish(row['bam'],incoming)
                arrivals.write(json.dumps(dict(batch=i,source_path=row['bam'],published_path=str(dest),
                    source_available_at=row['available_at'],published_epoch=time.time(),
                    elapsed=time.monotonic()-start))+'\n');arrivals.flush()
                if i%50==0:print(f'Published {i+1}/{len(rows)} batches; elapsed {time.monotonic()-start:.1f}s',flush=True)
        # Watchdog dispatch is asynchronous. Wait for its last events before flushing.
        for _ in range(60):
            if len(watcher.processed)>=len(rows):break
            await asyncio.sleep(1)
        await asyncio.to_thread(watcher.flush_remaining)
        if len(watcher.processed)!=len(rows):raise RuntimeError('Not every published BAM reached the file watcher')
        idle=0;deadline=time.monotonic()+args.drain_timeout
        while time.monotonic()<deadline:
            stats=await coord.stats.remote();last_stats=stats
            pending=sum(int(stats.get(k,0) or 0) for k in ('active_count','queued_count','waiting_global'))
            idle=idle+1 if pending==0 else 0
            if idle>=6:drained=True;break
            if telemetry.done():await telemetry
            await asyncio.sleep(5)
        if not drained:raise RuntimeError('Queue drain exceeded limit')
        # Flush native duration CSVs / failed-job log before checking outcomes.
        await coord.shutdown.remote()
        for dbpath in work.glob('*/sv_live/state.sqlite'):
            db=sqlite3.connect(dbpath)
            try:export(db,dbpath.parent,3)
            finally:db.close()
        expected=['preprocessing','cnv','bed_conversion','sturgeon','nanodx','mgmt','fusion','sv']
        completed=last_stats.get('completed_by_type',{})
        missing=[stage for stage in expected if not completed.get(stage)]
        artifact_checks={name:any(p.stat().st_size>100 for p in work.glob('*/'+name))
                         for name in ['CNV.npy','NanoDX_scores.csv','sturgeon_scores.csv']}
        summary=dict(artifact_checks=artifact_checks,drained=drained,input_batches=len(rows),elapsed_seconds=time.monotonic()-start,
            failed=last_stats.get('failed',0),completed_by_type=completed,
            failed_by_type=last_stats.get('failed_by_type',{}),missing_completed_stages=missing,
            watcher_batches=len(watcher.processed),stats=last_stats)
        (root/'summary.json').write_text(json.dumps(summary,indent=2,default=str)+'\n')
        print(json.dumps({k:v for k,v in summary.items() if k!='stats'},indent=2),flush=True)
        if summary['failed'] or missing or not all(artifact_checks.values()):raise RuntimeError('ROBIN reported failed or missing stages; see summary.json')
    finally:
        if observer:observer.stop();observer.join(timeout=10)
        if telemetry:
            telemetry.cancel()
            with contextlib.suppress(asyncio.CancelledError):await telemetry
        with contextlib.suppress(Exception):await coord.shutdown.remote()
        ray.shutdown()


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest',required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--reference',type=Path,required=True,help='Reference FASTA matching the archived alignments')
    p.add_argument('--panel',default='rCNS2');p.add_argument('--speed',type=float,default=1)
    p.add_argument('--limit',type=int);p.add_argument('--cpus',type=int,default=6)
    p.add_argument('--drain-timeout',type=float,default=14400)
    args=p.parse_args(argv)
    import math
    if not math.isfinite(args.speed) or args.speed<=0 or args.cpus<1 or (args.limit is not None and args.limit<1):p.error('Invalid speed, CPU count or batch limit')
    asyncio.run(simulate(args))


if __name__=='__main__':main()
