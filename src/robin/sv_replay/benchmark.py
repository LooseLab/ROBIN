"""Paced completed-BAM benchmark for an isolated incremental SV worker."""
import argparse
import collections
import csv
import hashlib
import json
import math
import os
import resource
import sqlite3
import time
from pathlib import Path

import pysam

from .interpretation import annotate
from .scoring import confidence
from .replay import connect, segment, signatures, add_signature, export
from .regions import RESOURCE_DIR, RegionContext, call_status


def prepare(summary, directories, duration, output):
    paths = {}
    for directory in directories:
        for path in Path(directory).glob('*.bam'):
            if path.name in paths:
                raise ValueError(f'Duplicate BAM basename: {path.name}')
            paths[path.name] = path.resolve()
    arrivals = {}
    with open(summary) as handle:
        header = handle.readline().rstrip().split('\t')
        fi, si, di = [header.index(k) for k in ('filename_bam', 'start_time', 'duration')]
        limit = max(fi, si, di) + 1
        for line in handle:
            row = line.rstrip('\n').split('\t', limit)
            name = row[fi]
            if name not in paths:
                continue
            at = float(row[si]) + float(row[di])
            if not math.isfinite(at) or at < 0:
                raise ValueError(f'Invalid completion time for {name}')
            arrivals[name] = max(arrivals.get(name, 0), at)
    selected = sorted((at, str(paths[name])) for name, at in arrivals.items() if at <= duration)
    if not selected:
        raise ValueError('No completed BAM batches selected')
    with open(output, 'w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['bam', 'available_at'])
        writer.writerows((path, at) for at, path in selected)
    Path(str(output) + '.json').write_text(json.dumps(dict(
        timing_basis='latest read completion per original BAM; excludes processing/publication delay',
        summary=str(Path(summary).resolve()), directories=[str(Path(p).resolve()) for p in directories],
        duration_seconds=duration, batches=len(selected),
        first_completion_seconds=selected[0][0], last_completion_seconds=selected[-1][0]), indent=2)+'\n')
    print(f'Prepared {len(selected)} original BAM batches', flush=True)


def schedule(arrival, available, service):
    """Single-worker FIFO model in source-run seconds (no artificial pacing)."""
    start = max(arrival, available)
    return start, start + service


class IncrementalWorker:
    """Single writer; re-extract changed reads, including late supplementary records."""
    def __init__(self, path, lengths, collect_metrics=True):
        self.collect_metrics = collect_metrics
        self.db = connect(path)
        self.db.execute('CREATE TABLE IF NOT EXISTS mapping_summary(call_id INTEGER PRIMARY KEY, data TEXT)')
        self.context = RegionContext(lengths, RESOURCE_DIR/'sv_grch38_tandem_repeats.bed',
                                     RESOURCE_DIR/'cenSatRegions.bed')
        self.db.execute('CREATE TABLE IF NOT EXISTS delivered_quality(id INTEGER PRIMARY KEY, support INTEGER)')
        self.db.execute('CREATE TABLE IF NOT EXISTS delivered(id TEXT PRIMARY KEY, support INTEGER, status TEXT, mask_flags TEXT)')

    def process(self, path, at):
        db = self.db
        changed = set()
        with db:
            with pysam.AlignmentFile(path, 'rb') as bam:
                if dict(zip(bam.references, bam.lengths)) != self.context.lengths:
                    raise ValueError('BAM reference dictionaries differ')
                records = 0
                for r in bam.fetch(until_eof=True):
                    records += 1
                    if r.is_unmapped or r.is_secondary or r.is_duplicate or r.is_qcfail:
                        continue
                    rid = r.query_name
                    db.execute('INSERT OR IGNORE INTO reads VALUES (?,?)', (rid, at))
                    segs = [segment(r.reference_name, r.reference_start,
                                    '-' if r.is_reverse else '+', r.cigarstring, r.mapping_quality,
                                    r.get_tag('NM') if self.collect_metrics and r.has_tag('NM') else None,actual=self.collect_metrics)]
                    if r.has_tag('SA'):
                        for entry in r.get_tag('SA').rstrip(';').split(';'):
                            chrom, pos, strand, cigar, mq, nm = entry.split(',')
                            segs.append(segment(chrom, int(pos)-1, strand, cigar, int(mq), int(nm) if self.collect_metrics else None))
                    for s in segs:
                        if db.execute('INSERT OR IGNORE INTO segments VALUES (?,?)',
                                      (rid, json.dumps(s, sort_keys=True))).rowcount:
                            changed.add(rid)
            touched = set()
            for rid in sorted(changed):
                segs = [json.loads(r[0]) for r in db.execute('SELECT data FROM segments WHERE read_id=?', (rid,))]
                for sig in signatures(segs, include_metrics=self.collect_metrics):
                    cid = add_signature(db, rid, at, sig, 100, self.context)
                    touched.add(cid)
                    db.executemany('INSERT OR IGNORE INTO call_flags VALUES (?,?)',
                                   [(cid, flag) for flag in sorted(self.context.flags(sig))])
            updates = []
            for cid in sorted(touched):
                support = db.execute('SELECT count(*) FROM evidence WHERE call_id=?', (cid,)).fetchone()[0]
                flags = [r[0] for r in db.execute('SELECT flag FROM call_flags WHERE call_id=? ORDER BY flag', (cid,))]
                status = call_status(support, flags, 3)
                flags = ','.join(flags)
                db.execute('INSERT OR REPLACE INTO history VALUES (?,?,?,?,?)', (cid, at, support, status, flags))
                old = db.execute('SELECT support,status,mask_flags FROM delivered WHERE id=?', (cid,)).fetchone()
                hq = db.execute('SELECT coalesce(sum(high_quality),0) FROM evidence_quality WHERE call_id=?',(cid,)).fetchone()[0]
                old_hq = db.execute('SELECT support FROM delivered_quality WHERE id=?',(cid,)).fetchone()
                if old != (support, status, flags) or old_hq != (hq,) or self.collect_metrics:
                    db.execute('INSERT OR REPLACE INTO delivered_quality VALUES (?,?)',(cid,hq))
                    db.execute('INSERT OR REPLACE INTO delivered VALUES (?,?,?,?)', (cid, support, status, flags))
                    position = db.execute('SELECT kind,size,c1,p1,o1,c2,p2,o2 FROM calls WHERE id=?', (cid,)).fetchone()
                    updates.append(dict(call_id=cid, support=support, status=status, mask_flags=flags,
                                        **dict(zip(('kind','size','c1','p1','o1','c2','p2','o2'),position))))
        for update in updates:
            update.update(annotate(update))
            if self.collect_metrics:
                summary=summarize_mapping(db,update['call_id'])
                update.update(summary)
                with db:db.execute('INSERT OR REPLACE INTO mapping_summary VALUES (?,?)',(update['call_id'],json.dumps(summary)))
            update['high_quality_support'] = db.execute('SELECT coalesce(sum(high_quality),0) FROM evidence_quality WHERE call_id=?',(update['call_id'],)).fetchone()[0]
            update.update(confidence(update))
        return dict(alignment_records=records, changed_read_ids=len(changed), updates=updates)


def disk_bytes(output):
    return sum(p.stat().st_size for p in output.iterdir() if p.is_file())


def run(manifest, output, mode, speed):
    if output.exists() and any(output.iterdir()):
        raise ValueError('Use a new empty benchmark directory; runs are not resumable')
    output.mkdir(parents=True, exist_ok=True)
    batches = []
    with open(manifest) as handle:
        for row in csv.DictReader(handle):
            at = float(row['available_at'])
            if not math.isfinite(at) or at < 0:
                raise ValueError('Invalid arrival time')
            batches.append((at, Path(row['bam']).resolve()))
    batches.sort(key=lambda r:(r[0],str(r[1])))
    if not batches:
        raise ValueError('Empty manifest')
    with pysam.AlignmentFile(batches[0][1], 'rb') as bam:
        lengths = dict(zip(bam.references, bam.lengths))
    worker = IncrementalWorker(output/'benchmark.sqlite', lengths)
    config = dict(mode=mode, speed=speed, batches=len(batches),
        timing_basis='manifest availability; see manifest provenance',
        manifest=str(Path(manifest).resolve()), manifest_sha256=hashlib.sha256(Path(manifest).read_bytes()).hexdigest(),
        inputs=[dict(path=str(p), bytes=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns,arrival=at) for at,p in batches],
        comparison_scope='isolated SV worker; no concurrent ROBIN analyses',
        masks=[dict(path=str(RESOURCE_DIR/name),sha256=hashlib.sha256((RESOURCE_DIR/name).read_bytes()).hexdigest())
               for name in ('sv_grch38_tandem_repeats.bed','cenSatRegions.bed')])
    (output/'run.json').write_text(json.dumps(config,indent=2)+'\n')
    wall_start=time.monotonic(); available=0; metrics=[]
    arrivals=[at/speed for at,_ in batches]
    try:
        with (output/'updates.jsonl').open('w') as updates, (output/'metrics.jsonl').open('w') as stream:
            for i,(source_at,path) in enumerate(batches):
                due=source_at/speed
                if mode=='paced':
                    while time.monotonic()-wall_start<due:
                        time.sleep(max(0,min(1,due-(time.monotonic()-wall_start))))
                actual_start=time.monotonic()-wall_start
                cpu_start=time.process_time(); started=time.monotonic()
                result=worker.process(path,source_at)
                # Publication cost is included in service time, including flush.
                updates.write(json.dumps(dict(batch=i,available_at=source_at,updates=result.pop('updates')))+'\n')
                updates.flush()
                service=time.monotonic()-started
                model_start, model_end=schedule(due,available,service);available=model_end
                end=time.monotonic()-wall_start
                row=dict(batch=i,bam=str(path),source_arrival=source_at,scheduled_arrival=due,
                    service_seconds=service,cpu_seconds=time.process_time()-cpu_start,
                    actual_start=actual_start,actual_end=end,
                    model_start=model_start,model_end=model_end,
                    model_latency=model_end-due,actual_latency=end-due if mode=='paced' else None,
                    queued_on_arrival=max(0,i-sum(m['model_end']<=due for m in metrics)),
                    model_backlog_at_finish=sum(t<=model_end for t in arrivals[i+1:]),
                    actual_backlog_at_finish=sum(t<=end for t in arrivals[i+1:]) if mode=='paced' else None,
                    peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
                    output_bytes=disk_bytes(output),input_bytes=path.stat().st_size,**result)
                metrics.append(row);stream.write(json.dumps(row)+'\n');stream.flush()
                print(f'Batch {i+1}/{len(batches)}: service={service:.2f}s model_latency={row["model_latency"]:.2f}s backlog={row["model_backlog_at_finish"]}',flush=True)
        export(worker.db,output,3)
        latencies=sorted(r['actual_latency'] if mode=='paced' else r['model_latency'] for r in metrics)
        source_span=arrivals[-1]-arrivals[0]
        summary=dict(mode=mode,batches=len(metrics),alignment_records=sum(r['alignment_records'] for r in metrics),
            total_service_seconds=sum(r['service_seconds'] for r in metrics),source_arrival_span_seconds=source_span,
            service_to_arrival_span_ratio=sum(r['service_seconds'] for r in metrics)/source_span if source_span else None,
            median_batch_latency_seconds=latencies[len(latencies)//2],
            p95_batch_latency_seconds=latencies[min(len(latencies)-1,math.ceil(len(latencies)*.95)-1)],
            maximum_batch_latency_seconds=max(latencies),max_modeled_backlog=max(r['model_backlog_at_finish'] for r in metrics),
            max_observed_backlog=max(r['actual_backlog_at_finish'] for r in metrics) if mode=='paced' else None,
            final_drain_seconds=max(0,(metrics[-1]['actual_end'] if mode=='paced' else available)-arrivals[-1]),
            peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,output_bytes=disk_bytes(output),
            wall_seconds_including_final_export=time.monotonic()-wall_start,
            qualification='Isolated worker benchmark. Does not measure contention or latency of other ROBIN analyses.')
        (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2),flush=True)
    finally:
        worker.db.close()


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    prep=sub.add_parser('prepare');prep.add_argument('--summary',required=True);prep.add_argument('--bam-dir',action='append',required=True)
    prep.add_argument('--duration',type=float,default=3600);prep.add_argument('--output',type=Path,required=True)
    test=sub.add_parser('run');test.add_argument('--manifest',required=True);test.add_argument('--output',type=Path,required=True)
    test.add_argument('--mode',choices=['modeled','paced'],default='modeled');test.add_argument('--speed',type=float,default=1)
    args=parser.parse_args(argv)
    if args.command=='prepare':
        if not math.isfinite(args.duration) or args.duration<=0:parser.error('duration must be positive and finite')
        prepare(args.summary,args.bam_dir,args.duration,args.output)
    else:
        if not math.isfinite(args.speed) or args.speed<=0:parser.error('speed must be positive and finite')
        run(args.manifest,args.output,args.mode,args.speed)




def summarize_mapping(db, call_id):
    import statistics
    rows=db.execute('SELECT m.mapq1,m.mapq2,m.error1,m.error2,m.event_size,e.p1,e.p2 '
                    'FROM evidence_mapping m JOIN evidence e USING(call_id,read_id) WHERE m.call_id=?',(call_id,)).fetchall()
    if not rows:return dict(mapping_measured_ids=0)
    def sd(index):
        values=[r[index] for r in rows if r[index] is not None]
        return statistics.pstdev(values) if values else None
    errors=[e for row in rows for e in row[2:4] if e is not None]
    return dict(mapping_measured_ids=len(rows),mapq30_support=sum(min(r[:2])>=30 for r in rows),
                mean_min_mapq=statistics.mean(min(r[:2]) for r in rows),
                substitution_rate_mean=statistics.mean(errors) if errors else None,
                substitution_measured_segments=len(errors),
                breakpoint1_sd=sd(5),breakpoint2_sd=sd(6),event_size_sd=sd(4),
                event_size_mean=statistics.mean(r[4] for r in rows if r[4] is not None) if any(r[4] is not None for r in rows) else None)


if __name__=='__main__':main()
