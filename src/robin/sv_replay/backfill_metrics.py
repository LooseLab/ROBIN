"""Backfill mapping metrics from archived batches without changing calls/support."""
import argparse,fcntl,json,sqlite3,time
from pathlib import Path
import pysam
from .replay import segment,signatures
from .benchmark import summarize_mapping

def run(root,follow=False):
    path=root/'state.sqlite'
    db=sqlite3.connect(path,timeout=30)
    db.executescript('CREATE TABLE IF NOT EXISTS evidence_mapping(call_id INTEGER,read_id TEXT,mapq1 INTEGER,mapq2 INTEGER,error1 REAL,error2 REAL,event_size INTEGER,PRIMARY KEY(call_id,read_id)); CREATE TABLE IF NOT EXISTS mapping_summary(call_id INTEGER PRIMARY KEY,data TEXT); CREATE TABLE IF NOT EXISTS metrics_backfill_batches(at REAL PRIMARY KEY, evidence_count INTEGER);')
    while True:
        batches={r['started_at']:r['filepath'] for r in (json.loads(s) for s in (root/'batches.jsonl').read_text().splitlines(keepends=True) if s.endswith('\n'))}
        pending=db.execute('SELECT DISTINCT e.at FROM evidence e LEFT JOIN evidence_mapping m USING(call_id,read_id) WHERE m.call_id IS NULL').fetchall()
        priority={r[0] for r in db.execute("SELECT DISTINCT e.at FROM evidence e JOIN calls c ON e.call_id=c.id WHERE c.id=376 OR (c.c1='chr9' AND abs(c.p1-20830605)<100 AND abs(c.p2-22875607)<100)")}
        pending.sort(key=lambda r:(r[0] not in priority,r[0]))
        processed=0
        for at, in pending:
            if at not in batches:continue
            targets=db.execute('SELECT e.call_id,e.read_id,c.kind,c.size,c.c1,c.p1,c.o1,c.c2,c.p2,c.o2 FROM evidence e JOIN calls c ON c.id=e.call_id WHERE e.at=?',(at,)).fetchall()
            previous=db.execute('SELECT evidence_count FROM metrics_backfill_batches WHERE at=?',(at,)).fetchone()
            if previous and previous[0]==len(targets):continue
            wanted={r[1] for r in targets};segs={rid:[] for rid in wanted}
            with pysam.AlignmentFile(batches[at],'rb') as bam:
                for r in bam.fetch(until_eof=True):
                    if r.query_name not in wanted or r.is_unmapped or r.is_secondary or r.is_duplicate or r.is_qcfail:continue
                    segs[r.query_name].append(segment(r.reference_name,r.reference_start,'-' if r.is_reverse else '+',r.cigarstring,r.mapping_quality,r.get_tag('NM') if r.has_tag('NM') else None,actual=True))
                    if r.has_tag('SA'):
                        for entry in r.get_tag('SA').rstrip(';').split(';'):
                            chrom,pos,strand,cigar,mq,nm=entry.split(',')
                            segs[r.query_name].append(segment(chrom,int(pos)-1,strand,cigar,int(mq),int(nm)))
            values=[]
            for cid,rid,kind,size,c1,p1,o1,c2,p2,o2 in targets:
                if not segs[rid]:continue
                sigs=signatures(segs[rid],include_metrics=True)
                matches=[s for s in sigs if s['kind']==kind and s['a'][0]==c1 and s['a'][2]==o1 and s['b'][0]==c2 and s['b'][2]==o2 and abs(s['a'][1]-p1)<=100 and abs(s['b'][1]-p2)<=100 and (kind!='INS' or min(size,s['size'])/max(size,s['size'])>=.7)]
                if not matches:continue
                s=min(matches,key=lambda s:abs(s['a'][1]-p1)+abs(s['b'][1]-p2));a,b=s['mapping_metrics']
                values.append((cid,rid,a['mapq'],b['mapq'],a['substitution_rate'],b['substitution_rate'],s['size'] if kind=='INS' else abs(s['b'][1]-s['a'][1]) if c1==c2 else None))
            # Match native worker locking; hold only for short write/summary transactions.
            with (root/'writer.lock').open('a') as lock:
                fcntl.flock(lock,fcntl.LOCK_EX)
                with db:
                    db.execute('INSERT OR REPLACE INTO metrics_backfill_batches VALUES (?,?)',(at,len(targets)))
                    db.executemany('INSERT OR REPLACE INTO evidence_mapping VALUES (?,?,?,?,?,?,?)',values)
                    for cid in {v[0] for v in values}:
                        db.execute('INSERT OR REPLACE INTO mapping_summary VALUES (?,?)',(cid,json.dumps(summarize_mapping(db,cid))))
            processed+=len(values)
            print(json.dumps(dict(batch=batches[at],measured=len(values),total_mapping=db.execute('SELECT count(*) FROM evidence_mapping').fetchone()[0])),flush=True)
        if not follow:break
        time.sleep(30 if processed else 60)
    db.close()

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('root',type=Path);p.add_argument('--follow',action='store_true');a=p.parse_args();run(a.root,a.follow)
