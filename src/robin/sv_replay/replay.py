"""Persistent offline SV replay. Coordinates are zero-based interbase positions."""
import argparse
import csv
import hashlib
import json
import math
import re
import sqlite3
from pathlib import Path

import pysam

from .interpretation import annotate
from .regions import RESOURCE_DIR, RegionContext, call_status

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def segment(chrom, start, strand, cigar, mapq, nm=None, actual=False):
    ops = [(int(n), op) for n, op in CIGAR.findall(cigar)]
    if not ops or ''.join(f'{n}{op}' for n, op in ops) != cigar:
        raise ValueError(f'Invalid CIGAR: {cigar}')
    total = sum(n for n, op in ops if op in 'MISH=X')
    leading = 0
    for n, op in ops:
        if op not in 'SH':
            break
        leading += n
    aligned = sum(n for n, op in ops if op in 'MI=X')
    qstart, qend = leading, leading + aligned
    if strand == '-':
        qstart, qend = total - qend, total - qstart
    result = dict(chrom=chrom, start=start,
                end=start + sum(n for n, op in ops if op in 'MDN=X'),
                strand=strand, cigar=cigar, mapq=mapq,
                qs=qstart, qe=qend)
    if nm is not None:result["nm"] = int(nm)
    if actual:result["actual_record"] = True
    return result


def canonical(a, b, kind, size=0, anchors=None, metrics=None):
    swapped = tuple(a) > tuple(b)
    a, b = sorted((tuple(a), tuple(b)))
    result = dict(a=a, b=b, kind=kind, size=size)
    if anchors is not None:
        result['anchors'] = list(reversed(anchors)) if swapped else list(anchors)
    if metrics is not None:
        result['mapping_metrics'] = list(reversed(metrics)) if swapped else list(metrics)
    return result


def local_anchor(ops):
    """Consecutive aligned CIGAR bases next to a junction; M is not proven identity."""
    count = 0
    for n, op in ops:
        if op not in 'M=':
            break
        count += n
    return count


def segment_anchor(s, exit):
    ops = [(int(n), op) for n, op in CIGAR.findall(s['cigar'])]
    # Reference-right is query exit on +, query entry on -.
    if (s['strand'] == '+') == exit:
        ops.reverse()
    while ops and ops[0][1] in 'SH':
        ops.pop(0)
    return local_anchor(ops)


def mapping_metrics(s):
    ops = [(int(n), op) for n, op in CIGAR.findall(s['cigar'])]
    aligned = sum(n for n, op in ops if op in 'M=X')
    indels = sum(n for n, op in ops if op in 'ID')
    nm = s.get('nm')
    # Whole-segment estimate, not breakpoint-local identity. Remove indels,
    # including the SV itself, so a genuine deletion is not a mismatch penalty.
    rate = (nm-indels)/aligned if nm is not None and indels <= nm <= indels+aligned and aligned else None
    return dict(mapq=s['mapq'], substitution_rate=rate)


def signatures(segments, min_mapq=20, min_anchor=100, min_sv=50,
               max_overlap=100, max_gap=100, include_metrics=False):
    """Combine SA/record segments once, then emit CIGAR and query-adjacent junctions."""
    quality_records = {}
    if include_metrics:
        for item in segments:
            if item.get('actual_record'):
                key=(item['chrom'],item['start'],item['strand'])
                if key not in quality_records or item['mapq']>quality_records[key]['mapq']:
                    quality_records[key]=item
    def metric_for(item):
        return mapping_metrics(quality_records.get((item['chrom'],item['start'],item['strand']),item))
    unique = {}
    for s in segments:
        key = (s['chrom'], s['start'], s['strand'], s['cigar'].replace('H', 'S'))
        if key not in unique or s['mapq'] > unique[key]['mapq']:
            unique[key] = s
    good = [s for s in unique.values()
            if s['mapq'] >= min_mapq and s['qe'] - s['qs'] >= min_anchor]
    out = []
    for s in good:
        ops = [(int(n), op) for n, op in CIGAR.findall(s['cigar'])]
        rpos = s['start']
        aligned_total = sum(n for n, op in ops if op in 'M=X')
        left = 0
        for index, (n, op) in enumerate(ops):
            # Clips and insertions do not contribute to alignment anchors.
            right = aligned_total - left
            if op in 'ID' and n >= min_sv and min(left, right) >= min_anchor:
                end = rpos + n if op == 'D' else rpos
                out.append(canonical((s['chrom'], rpos, '+'),
                                     (s['chrom'], end, '-'),
                                     'DEL' if op == 'D' else 'INS', n,
                                     anchors=(local_anchor(ops[:index][::-1]), local_anchor(ops[index+1:])),
                                     metrics=[metric_for(s)]*2 if include_metrics else None))
            if op in 'MDN=X':
                rpos += n
            if op in 'M=X':
                left += n
    chain = sorted(good, key=lambda s: (s['qs'], s['qe'], s['chrom'], s['start']))
    for x, y in zip(chain, chain[1:]):
        gap = y['qs'] - x['qe']
        if not -max_overlap <= gap <= max_gap:
            continue
        a = (x['chrom'], x['end'] if x['strand'] == '+' else x['start'],
             '+' if x['strand'] == '+' else '-')
        b = (y['chrom'], y['start'] if y['strand'] == '+' else y['end'],
             '-' if y['strand'] == '+' else '+')
        if a[0] == b[0] and a[2] != b[2] and abs(a[1] - b[1]) < min_sv:
            continue
        out.append(canonical(a, b, 'BND', anchors=(segment_anchor(x, exit=True), segment_anchor(y, exit=False)),
                             metrics=[metric_for(x),metric_for(y)] if include_metrics else None))
    # One molecule can carry multiple distinct junctions, never duplicate signatures.
    return list({json.dumps(s, sort_keys=True): s for s in out}.values())


def connect(path):
    db = sqlite3.connect(path)
    db.execute('PRAGMA journal_mode=WAL')
    db.executescript('''
    CREATE TABLE IF NOT EXISTS meta(key TEXT PRIMARY KEY, value TEXT NOT NULL);
    CREATE TABLE IF NOT EXISTS reads(id TEXT PRIMARY KEY, at REAL NOT NULL);
    CREATE TABLE IF NOT EXISTS segments(read_id TEXT, data TEXT,
      PRIMARY KEY(read_id,data));
    CREATE TABLE IF NOT EXISTS processed(read_id TEXT PRIMARY KEY);
    CREATE TABLE IF NOT EXISTS calls(id INTEGER PRIMARY KEY, kind TEXT, size INTEGER,
      c1 TEXT,p1 INTEGER,o1 TEXT,c2 TEXT,p2 INTEGER,o2 TEXT,first REAL,last REAL);
    CREATE INDEX IF NOT EXISTS junction_lookup ON calls(kind,c1,o1,c2,o2,p1);
    CREATE TABLE IF NOT EXISTS evidence(call_id INTEGER,read_id TEXT,at REAL,
      p1 INTEGER,p2 INTEGER, PRIMARY KEY(call_id,read_id));
    CREATE TABLE IF NOT EXISTS history(call_id INTEGER,at REAL,support INTEGER,
      status TEXT, mask_flags TEXT DEFAULT '', PRIMARY KEY(call_id,at));
    CREATE TABLE IF NOT EXISTS call_flags(call_id INTEGER, flag TEXT,
      PRIMARY KEY(call_id,flag));
    ''')
    db.execute('CREATE TABLE IF NOT EXISTS evidence_quality(call_id INTEGER, read_id TEXT, anchor1 INTEGER, anchor2 INTEGER, high_quality INTEGER, repeat_context TEXT, PRIMARY KEY(call_id,read_id))')
    db.execute('CREATE TABLE IF NOT EXISTS evidence_mapping(call_id INTEGER, read_id TEXT, mapq1 INTEGER, mapq2 INTEGER, error1 REAL, error2 REAL, event_size INTEGER, PRIMARY KEY(call_id,read_id))')
    return db


def load_times(path):
    if path is None:
        return None
    times = {}
    with open(path, newline='') as handle:
        for row in csv.DictReader(handle):
            rid, at = row['read_id'], float(row['available_at'])
            if at < 0 or not math.isfinite(at):
                raise ValueError(f'Invalid available_at for {rid}')
            if rid in times and times[rid] != at:
                raise ValueError(f'Conflicting availability for {rid}')
            times[rid] = at
    return times


def stage(db, bams, times, batch_seconds):
    """Disk staging reunites primary/supplementary records even across BAM batches."""
    with db:
        for batch, path in enumerate(bams):
            with pysam.AlignmentFile(path, 'rb') as bam:
                for r in bam.fetch(until_eof=True):
                    if r.is_unmapped or r.is_secondary or r.is_duplicate or r.is_qcfail:
                        continue
                    rid = r.query_name
                    if times is not None and rid not in times:
                        raise ValueError(f'Missing read availability: {rid}')
                    at = times[rid] if times is not None else batch * batch_seconds
                    db.execute('INSERT INTO reads VALUES (?,?) ON CONFLICT(id) '
                               'DO UPDATE SET at=max(at,excluded.at)', (rid, at))
                    segs = [segment(r.reference_name, r.reference_start,
                                    '-' if r.is_reverse else '+',
                                    r.cigarstring, r.mapping_quality)]
                    if r.has_tag('SA'):
                        for entry in r.get_tag('SA').rstrip(';').split(';'):
                            chrom, pos, strand, cigar, mq, _nm = entry.split(',')
                            segs.append(segment(chrom, int(pos) - 1, strand, cigar, int(mq)))
                    db.executemany('INSERT OR IGNORE INTO segments VALUES (?,?)',
                                   [(rid, json.dumps(s, sort_keys=True)) for s in segs])


def add_signature(db, rid, at, sig, tolerance, context=None):
    a, b = sig['a'], sig['b']
    # Fixed seed coordinates bound cluster drift; both oriented ends must agree.
    matches = db.execute('''SELECT id,p1,p2,size FROM calls WHERE kind=?
      AND c1=? AND o1=? AND c2=? AND o2=? AND p1 BETWEEN ? AND ?
      AND p2 BETWEEN ? AND ? ORDER BY id''',
      (sig['kind'], a[0], a[2], b[0], b[2], a[1]-tolerance,
       a[1]+tolerance, b[1]-tolerance, b[1]+tolerance)).fetchall()
    matches = [m for m in matches if sig['kind'] != 'INS' or
               min(m[3], sig['size']) / max(m[3], sig['size']) >= .7]
    if matches:
        cid = min(matches, key=lambda m: (abs(m[1]-a[1])+abs(m[2]-b[1]), m[0]))[0]
    else:
        cid = db.execute('INSERT INTO calls(kind,size,c1,p1,o1,c2,p2,o2,first,last) '
                         'VALUES (?,?,?,?,?,?,?,?,?,?)',
                         (sig['kind'], sig['size'], *a, *b, at, at)).lastrowid
    db.execute('INSERT OR IGNORE INTO evidence VALUES (?,?,?,?,?)',
               (cid, rid, at, a[1], b[1]))
    if 'anchors' in sig:
        anchors = sig['anchors']
        repeat_context = []
        for chrom, pos, orientation in (a, b):
            base = max(0, pos-1 if orientation == '+' else pos)
            repeat_context.append(context.indexes[0][1].overlap_context(chrom, base) if context else None)
        # Experimental local-anchor support, not independent unique-mapping validation.
        high = int(min(anchors) >= 100)
        db.execute('INSERT INTO evidence_quality VALUES (?,?,?,?,?,?) ON CONFLICT(call_id,read_id) DO UPDATE SET anchor1=excluded.anchor1,anchor2=excluded.anchor2,high_quality=excluded.high_quality,repeat_context=excluded.repeat_context',
                   (cid, rid, *anchors, high, json.dumps(repeat_context) if context else None))
    if 'mapping_metrics' in sig:
        m1,m2=sig['mapping_metrics']
        event_size=sig['size'] if sig['kind']=='INS' else abs(b[1]-a[1]) if a[0]==b[0] else None
        db.execute('INSERT OR REPLACE INTO evidence_mapping VALUES (?,?,?,?,?,?,?)',
                   (cid,rid,m1['mapq'],m2['mapq'],m1['substitution_rate'],m2['substitution_rate'],event_size))
    db.execute('UPDATE calls SET last=max(last,?) WHERE id=?', (at, cid))
    return cid


def replay(db, args, context):
    cursor = db.execute('SELECT r.id,r.at FROM reads r LEFT JOIN processed p '
                        'ON r.id=p.read_id WHERE p.read_id IS NULL ORDER BY r.at,r.id')
    with db:
        for count, (rid, at) in enumerate(cursor, 1):
            segs = [json.loads(row[0]) for row in db.execute(
                'SELECT data FROM segments WHERE read_id=?', (rid,))]
            touched = set()
            for sig in signatures(segs, args.min_mapq, args.min_anchor,
                                  args.min_sv, args.max_overlap, args.max_gap):
                cid = add_signature(db, rid, at, sig, args.cluster_distance, context)
                touched.add(cid)
                db.executemany('INSERT OR IGNORE INTO call_flags VALUES (?,?)',
                               [(cid, flag) for flag in sorted(context.flags(sig))])
            for cid in touched:
                support = db.execute('SELECT count(*) FROM evidence WHERE call_id=?',
                                     (cid,)).fetchone()[0]
                flags = [r[0] for r in db.execute(
                    'SELECT flag FROM call_flags WHERE call_id=? ORDER BY flag', (cid,))]
                status = call_status(support, flags, args.min_support)
                db.execute('INSERT OR REPLACE INTO history VALUES (?,?,?,?,?)',
                           (cid, at, support, status, ','.join(flags)))
            db.execute('INSERT INTO processed VALUES (?)', (rid,))
            if count % 1000 == 0:
                db.commit()


def export(db, output, min_support):
    queries = {
        'calls.tsv': ('''SELECT c.id,c.kind,c.c1,c.p1,c.o1,c.c2,c.p2,c.o2,c.size,
          count(e.read_id) AS support,c.first,c.last,
          min(e.p1) AS p1_min,max(e.p1) AS p1_max,
          min(e.p2) AS p2_min,max(e.p2) AS p2_max,
          CASE WHEN EXISTS(SELECT 1 FROM call_flags f WHERE f.call_id=c.id
             AND f.flag='mitochondrial_reference_boundary') THEN 'reference_boundary'
           WHEN count(e.read_id)<? THEN 'candidate'
           WHEN EXISTS(SELECT 1 FROM call_flags f WHERE f.call_id=c.id)
             THEN 'masked_candidate' ELSE 'provisional' END AS status,
          COALESCE((SELECT group_concat(flag, ',') FROM
            (SELECT flag FROM call_flags WHERE call_id=c.id ORDER BY flag)), '') AS mask_flags
          FROM calls c JOIN evidence e ON c.id=e.call_id GROUP BY c.id ORDER BY c.id''',
          (min_support,)),
        'evidence.tsv': ('SELECT * FROM evidence ORDER BY at,read_id,call_id', ()),
        'history.tsv': ('SELECT * FROM history ORDER BY at,call_id', ()),
    }
    for filename, (query, params) in queries.items():
        cur = db.execute(query, params)
        temp = output / (filename + '.tmp')
        with temp.open('w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t')
            fields = [c[0] for c in cur.description]
            if filename == 'calls.tsv':
                writer.writerow(fields + ['topology', 'interval_bp', 'high_quality_support'])
                for row in cur:
                    annotation = annotate(dict(zip(fields, row)))
                    writer.writerow([*row, annotation['topology'], annotation['interval_bp'], db.execute('SELECT coalesce(sum(high_quality),0) FROM evidence_quality WHERE call_id=?',(row[0],)).fetchone()[0]])
            else:
                writer.writerow(fields)
                writer.writerows(cur)
        temp.replace(output / filename)


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bam', action='append', required=True,
                   help='Completed aligned BAM; repeat in arrival order.')
    p.add_argument('--read-times', help='CSV: read_id,available_at (seconds from run start).')
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--tandem-repeats', type=Path,
                   default=RESOURCE_DIR / 'sv_grch38_tandem_repeats.bed')
    p.add_argument('--centromeres', type=Path, default=RESOURCE_DIR / 'cenSatRegions.bed')
    p.add_argument('--mask-bed', type=Path, action='append', default=[],
                   help='Additional BED, e.g. low-complexity intervals; repeatable.')
    p.add_argument('--chromosome-end-margin', type=int, default=1_000_000)
    p.add_argument('--mitochondrial-origin-margin', type=int, default=100)
    for name, default in [('min-mapq',20), ('min-anchor',100), ('min-sv',50),
                          ('max-overlap',100), ('max-gap',100),
                          ('cluster-distance',100), ('min-support',3),
                          ('batch-seconds',60)]:
        p.add_argument('--' + name, type=int, default=default)
    return p


def main(argv=None):
    p = parser()
    args = p.parse_args(argv)
    for name in ('min_mapq','min_anchor','min_sv','max_overlap','max_gap',
                 'cluster_distance','min_support','batch_seconds',
                 'chromosome_end_margin','mitochondrial_origin_margin'):
        if getattr(args, name) < (0 if name in ('min_mapq','max_overlap','max_gap','cluster_distance',
                                                  'chromosome_end_margin','mitochondrial_origin_margin') else 1):
            p.error(f'Invalid {name}')
    lengths = {}
    for path in args.bam:
        with pysam.AlignmentFile(path, 'rb') as bam:
            current = dict(zip(bam.references, bam.lengths))
            if lengths and current != lengths:
                raise ValueError('BAM reference dictionaries differ.')
            lengths = current
    context = RegionContext(lengths, args.tandem_repeats, args.centromeres,
                            args.mask_bed, args.chromosome_end_margin,
                            args.mitochondrial_origin_margin)
    args.output.mkdir(parents=True, exist_ok=True)
    # Bind a checkpoint to immutable inputs/configuration; no silent mixed analyses.
    config = {k: v for k,v in vars(args).items() if k not in ('output','bam','read_times','tandem_repeats','centromeres','mask_bed')}
    config['inputs'] = [(str(Path(f).resolve()), Path(f).stat().st_size,
                         Path(f).stat().st_mtime_ns) for f in args.bam]
    config['timing_mode'] = 'read_availability' if args.read_times else 'bam_batch_order'
    config['read_times_sha256'] = (hashlib.sha256(Path(args.read_times).read_bytes()).hexdigest()
                                   if args.read_times else None)
    config['schema_version'] = 2
    config['region_masks'] = [dict(role=role, path=str(path.resolve()),
        sha256=hashlib.sha256(path.read_bytes()).hexdigest())
        for role, path in [('tandem_repeat', args.tandem_repeats),
                           ('centromere_satellite', args.centromeres)] +
                          [('custom_mask', p) for p in args.mask_bed]]
    encoded = json.dumps(config, sort_keys=True)
    db = connect(args.output / 'replay.sqlite')
    try:
        previous = db.execute("SELECT value FROM meta WHERE key='config'").fetchone()
        if previous and previous[0] != encoded:
            raise ValueError('Checkpoint inputs/settings differ; use a new output directory.')
        with db:
            db.execute("INSERT OR IGNORE INTO meta VALUES ('config',?)", (encoded,))
        if not db.execute("SELECT 1 FROM meta WHERE key='staged'").fetchone():
            stage(db, args.bam, load_times(args.read_times), args.batch_seconds)
            with db:
                db.execute("INSERT INTO meta VALUES ('staged','true')")
        replay(db, args, context)
        export(db, args.output, args.min_support)
        (args.output / 'run.json').write_text(json.dumps(config, indent=2) + '\n')
        print(f"Replay complete: {db.execute('SELECT count(*) FROM processed').fetchone()[0]} reads; "
              f"{db.execute('SELECT count(*) FROM calls').fetchone()[0]} candidate clusters. "
              f"Outputs: {args.output}")
    finally:
        db.close()


if __name__ == '__main__':
    main()
