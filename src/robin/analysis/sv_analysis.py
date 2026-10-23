"""Experimental incremental SV stage for ROBIN's native workflow."""
import fcntl
import json
import time
from pathlib import Path

import pysam

from robin.sv_replay.benchmark import IncrementalWorker


def sv_handler(job, work_dir=None):
    context = job.context
    work_dir = work_dir or context.metadata.get('work_dir')
    if not work_dir:
        raise ValueError('SV stage requires work_dir')
    sample = context.get_sample_id()
    if not sample or Path(sample).name != sample or sample in ('.', '..'):
        raise ValueError('Invalid SV sample ID')
    root = Path(work_dir) / sample / 'sv_live'
    root.mkdir(parents=True, exist_ok=True)
    batch = context.metadata.get('_batched_job')
    paths = batch.get_filepaths() if batch else [context.filepath]
    stats = []
    # Durable shared state, serial across Ray worker processes and task retries.
    with (root / 'writer.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        with pysam.AlignmentFile(paths[0], 'rb') as bam:
            lengths = dict(zip(bam.references, bam.lengths))
        worker = IncrementalWorker(root / 'state.sqlite', lengths)
        try:
            with (root / 'updates.jsonl').open('a') as updates, (root / 'batches.jsonl').open('a') as metrics:
                for path in paths:
                    started = time.time()
                    result = worker.process(path, started)
                    changed = result.pop('updates')
                    updates.write(json.dumps(dict(filepath=str(path),published_at=time.time(),updates=changed))+'\n')
                    updates.flush()
                    record = dict(filepath=str(path),started_at=started,finished_at=time.time(),
                                  changed_calls=len(changed),**result)
                    metrics.write(json.dumps(record)+'\n');metrics.flush();stats.append(record)
        finally:
            worker.db.close()
    context.add_result('sv', dict(status='success',output_dir=str(root),batches=stats))
