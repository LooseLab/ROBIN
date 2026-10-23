import tempfile
import unittest
from pathlib import Path
from robin.sv_replay.replay import segment,signatures,connect,add_signature,mapping_metrics
from robin.sv_replay.benchmark import summarize_mapping

class MappingTests(unittest.TestCase):
    def test_nm_subtracts_sv_and_missing_stays_unknown(self):
        self.assertEqual(mapping_metrics(segment('chr1',0,'+','100M50D100M',60,54))['substitution_rate'],.02)
        self.assertIsNone(mapping_metrics(segment('chr1',0,'+','200M',60))['substitution_rate'])
        self.assertIsNone(mapping_metrics(segment('chr1',0,'+','100M50D100M',60,2))['substitution_rate'])

    def test_malformed_sa_nm_and_actual_record_preference(self):
        bad=segment('chr1',0,'+','200M200S',60,758331583)
        self.assertIsNone(mapping_metrics(bad)['substitution_rate'])
        actual=segment('chr1',0,'+','100M1I99M200S',60,5,actual=True)
        other=segment('chr2',500,'+','200S200M',60,2)
        sigs=signatures([bad,actual,other],include_metrics=True)
        bnd=next(s for s in sigs if s['kind']=='BND')
        self.assertAlmostEqual(bnd['mapping_metrics'][0]['substitution_rate'],4/199)

    def test_canonical_metrics_follow_breakends(self):
        sig=signatures([segment('chr9',100,'+','200M200S',30,4),segment('chr1',500,'+','200S200M',60,2)],include_metrics=True)[0]
        self.assertEqual([m['mapq'] for m in sig['mapping_metrics']],[60,30])

    def test_unique_read_consistency_stats(self):
        with tempfile.TemporaryDirectory() as d:
            db=connect(Path(d)/'state.sqlite')
            for read,start,length,mq in [('a',100,50,25),('b',110,70,60)]:
                sig=signatures([segment('chr1',start,'+','200M'+str(length)+'D200M',mq,length+4)],include_metrics=True)[0]
                cid=add_signature(db,read,0,sig,100)
                add_signature(db,read,1,sig,100)
            stats=summarize_mapping(db,cid)
            self.assertEqual(stats['mapping_measured_ids'],2)
            self.assertEqual(stats['mapq30_support'],1)
            self.assertEqual(stats['breakpoint1_sd'],5)
            self.assertEqual(stats['breakpoint2_sd'],15)
            self.assertEqual(stats['event_size_sd'],10)
            self.assertAlmostEqual(stats['substitution_rate_mean'],.01)
            db.close()
