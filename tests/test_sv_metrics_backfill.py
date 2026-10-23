import json,sqlite3,tempfile,unittest
from pathlib import Path
import pysam
from robin.sv_replay.benchmark import IncrementalWorker
from robin.sv_replay.backfill_metrics import run
from robin.gui.components.sv_data import load_page
class BackfillTests(unittest.TestCase):
    def test_backfill_preserves_calls_support_and_score_sort(self):
        with tempfile.TemporaryDirectory() as tmp:
            sample=Path(tmp);root=sample/'sv_live';root.mkdir();bam=sample/'input.bam'
            with pysam.AlignmentFile(bam,'wb',header={'HD':{'VN':'1.6'},'SQ':[{'SN':'chr1','LN':10000000}]}) as out:
                for i in range(5):
                    r=pysam.AlignedSegment();r.query_name=str(i);r.query_sequence='A'*400;r.flag=0;r.reference_id=0;r.reference_start=2000000+i;r.mapping_quality=60;r.cigarstring='200M50D200M';r.set_tag('NM',54);out.write(r)
            worker=IncrementalWorker(root/'state.sqlite',{'chr1':10000000},collect_metrics=False)
            worker.process(bam,0);before=worker.db.execute('SELECT * FROM evidence ORDER BY read_id').fetchall();worker.db.close()
            (root/'batches.jsonl').write_text(json.dumps(dict(started_at=0,filepath=str(bam)))+'\n')
            run(root);run(root)
            db=sqlite3.connect(root/'state.sqlite')
            self.assertEqual(before,db.execute('SELECT * FROM evidence ORDER BY read_id').fetchall())
            self.assertEqual(db.execute('SELECT count(*) FROM evidence_mapping').fetchone()[0],5)
            db.close()
            page=load_page(sample,sort_by='confidence_score');row=page['rows'][0]
            self.assertEqual(row['mapq30_support'],5)
            self.assertAlmostEqual(row['substitution_rate_mean'],.01)
            self.assertGreaterEqual(row['confidence_score'],70)
            self.assertEqual(load_page(sample,sort_by='mapq30_support')['rows'][0]['id'],row['id'])
