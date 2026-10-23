import tempfile
import unittest
from pathlib import Path
from robin.sv_replay.replay import segment, signatures, connect, add_signature
from robin.sv_replay.regions import IntervalIndex

class QualityTests(unittest.TestCase):
    def test_local_anchor_stops_at_indel(self):
        sig=signatures([segment('chr9',100,'+','200M1I20M100D200M',60)])[0]
        self.assertEqual(sig['anchors'],[20,200])
        with tempfile.TemporaryDirectory() as d:
            db=connect(Path(d)/'test.sqlite')
            cid=add_signature(db,'r',0,sig,100)
            add_signature(db,'r',1,sig,100)
            self.assertEqual(db.execute('SELECT anchor1,anchor2,high_quality FROM evidence_quality').fetchall(),[(20,200,0)])
            self.assertEqual(db.execute('SELECT count(*) FROM evidence').fetchone()[0],1)
            db.close()

    def test_reverse_split_anchors(self):
        sig=signatures([segment('chr9',100,'-','200S150M',60),segment('chr10',100,'+','150S200M',60)])[0]
        self.assertEqual(sorted(sig['anchors']),[150,200])

    def test_repeat_union_edges(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'mask.bed';p.write_text('chr9\t10\t20\nchr9\t15\t30\n')
            index=IntervalIndex(p)
            self.assertEqual(index.overlap_context('chr9',18),dict(start=10,end=30,edge_distance=8))
            self.assertIsNone(index.overlap_context('chr9',30))
