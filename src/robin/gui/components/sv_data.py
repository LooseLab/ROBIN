"""Read-only paged access to the experimental SV worker checkpoint."""
import sqlite3
import json
from robin.sv_replay.scoring import confidence
from robin.sv_replay.interpretation import annotate
from pathlib import Path

STATUSES = {'all', 'candidate', 'provisional', 'masked_candidate', 'reference_boundary'}
MASKS = {'all', 'unmasked', 'tandem_repeat', 'centromere_satellite',
         'chromosome_end_margin', 'custom_mask', 'mitochondrial_reference_boundary'}


def database_path(sample_dir):
    return Path(sample_dir) / 'sv_live' / 'state.sqlite'


def _connect(path):
    db = sqlite3.connect(Path(path).resolve().as_uri() + '?mode=ro', uri=True, timeout=1)
    db.row_factory = sqlite3.Row
    db.execute('PRAGMA query_only=ON')
    db.execute('BEGIN')
    return db


def load_page(sample_dir, *, status='all', mask='all', search='', min_support=3,
              page=1, page_size=50, sort_by='support', descending=True, show_low_mapq=False):
    if status not in STATUSES or mask not in MASKS:
        raise ValueError('Invalid SV filter')
    if page < 1 or not 1 <= page_size <= 200 or min_support < 0:
        raise ValueError('Invalid SV pagination/support')
    sorts = {'id':'c.id', 'size':'c.size', 'support':'d.support',
             'interval_bp':"CASE WHEN c.c1=c.c2 THEN abs(c.p2-c.p1) END",
             'high_quality_support':'hq_support', 'confidence_score':'score'}
    for key in ['mapq30_support','mean_min_mapq','substitution_rate_mean','mapping_measured_ids','substitution_measured_segments','breakpoint1_sd','breakpoint2_sd','event_size_sd']:
        sorts[key]="json_extract(mapping_json,'$."+key+"')"
    if sort_by not in sorts:raise ValueError('Invalid SV sort column')
    path = database_path(sample_dir)
    if not path.is_file():
        return dict(available=False, rows=[], total=0, counts={}, page=1, pages=1)
    db = _connect(path)
    try:
        has_summary=db.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='mapping_summary'").fetchone()
        metrics_sql='(SELECT data FROM mapping_summary m WHERE m.call_id=c.id)' if has_summary else 'NULL'
        clauses=['d.support>=?'];params=[int(min_support)]
        if not show_low_mapq:
            mean_mapq_sql="json_extract("+metrics_sql+",'$.mean_min_mapq')"
            clauses.append('('+mean_mapq_sql+' IS NULL OR '+mean_mapq_sql+'>=55)')
        if status != 'all':clauses.append('d.status=?');params.append(status)
        if mask == 'unmasked':clauses.append("d.mask_flags=''")
        elif mask != 'all':clauses.append("instr(','||d.mask_flags||',',?)>0");params.append(','+mask+',')
        if search.strip():
            clauses.append('(c.c1 LIKE ? OR c.c2 LIKE ? OR CAST(c.id AS TEXT)=?)')
            params.extend(['%'+search.strip()+'%','%'+search.strip()+'%',search.strip()])
        source=' FROM calls c JOIN delivered d ON d.id=CAST(c.id AS TEXT) WHERE '+' AND '.join(clauses)
        total=db.execute('SELECT count(*)'+source,params).fetchone()[0]
        pages=max(1,(total+page_size-1)//page_size);page=min(page,pages)
        has_quality = db.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='evidence_quality'").fetchone()
        hq_sql = '(SELECT sum(q.high_quality) FROM evidence_quality q WHERE q.call_id=c.id)' if has_quality else 'NULL'
        has_summary=db.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='mapping_summary'").fetchone()
        metrics_sql='(SELECT data FROM mapping_summary m WHERE m.call_id=c.id)' if has_summary else 'NULL'
        db.create_function('sv_score',5,lambda n,hq,status,flags,data:confidence(dict(support=n,high_quality_support=hq,status=status,mask_flags=flags,**(json.loads(data) if data else {})))['confidence_score'])
        score_sql='sv_score(d.support,'+hq_sql+',d.status,d.mask_flags,'+metrics_sql+')'
        order = sorts[sort_by] + (' DESC' if descending else ' ASC') + ',c.id ASC'
        rows=[dict(r) for r in db.execute('SELECT c.id,c.kind,c.size,c.c1,c.p1,c.o1,c.c2,c.p2,c.o2,'
            'd.support,d.status,d.mask_flags,'+hq_sql+' AS hq_support,'+score_sql+' AS score,'+metrics_sql+' AS mapping_json'+source+' ORDER BY '+order+' LIMIT ? OFFSET ?',
            [*params,page_size,(page-1)*page_size])]
        counts=dict(db.execute('SELECT status,count(*) FROM delivered GROUP BY status').fetchall())
        has_quality = db.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='evidence_quality'").fetchone()
        for row in rows:
            row.pop('score')
            row.update(json.loads(row.pop('mapping_json') or '{}'))
            row['high_quality_support'] = row.pop('hq_support')
            row.update(confidence(row))
            row.update(annotate(row))
            row['breakend1']=f"{row['c1']}:{row['p1']} ({row['o1']})"
            row['breakend2']=f"{row['c2']}:{row['p2']} ({row['o2']})"
            row['mask_flags']=row['mask_flags'] or 'None'
        return dict(available=True,rows=rows,total=total,counts=counts,page=page,pages=pages)
    finally:db.close()


def load_evidence(sample_dir, call_id, limit=100):
    if not 1 <= limit <= 1000:raise ValueError('Invalid evidence limit')
    db=_connect(database_path(sample_dir))
    try:
        count=db.execute('SELECT count(*) FROM evidence WHERE call_id=?',(int(call_id),)).fetchone()[0]
        has_quality = db.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='evidence_quality'").fetchone()
        hq_sql = '(SELECT sum(q.high_quality) FROM evidence_quality q WHERE q.call_id=c.id)' if has_quality else 'NULL'
        rows=[dict(r) for r in db.execute('SELECT read_id,p1,p2 FROM evidence WHERE call_id=? '
            'ORDER BY read_id LIMIT ?',(int(call_id),limit))]
        if db.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='evidence_quality'").fetchone():
            for row in rows:
                quality=db.execute('SELECT anchor1,anchor2,high_quality,repeat_context FROM evidence_quality WHERE call_id=? AND read_id=?',(int(call_id),row['read_id'])).fetchone()
                if quality:row.update(dict(quality))
        if db.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='evidence_mapping'").fetchone():
            for row in rows:
                mapping=db.execute('SELECT mapq1,mapq2,error1,error2 FROM evidence_mapping WHERE call_id=? AND read_id=?',(int(call_id),row['read_id'])).fetchone()
                if mapping:row.update(dict(mapping))
        return dict(total=count,rows=rows)
    finally:db.close()
