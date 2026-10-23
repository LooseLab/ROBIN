import unittest
from robin.gui.components.sv_data import confidence,load_page
class ConfidenceTests(unittest.TestCase):
    def row(self,n=10,**kwargs):
        row=dict(support=n,high_quality_support=n,status='provisional',mask_flags='',mapping_measured_ids=n,mapq30_support=n,substitution_measured_segments=2*n,substitution_rate_mean=.01,breakpoint1_sd=2,breakpoint2_sd=2,event_size_sd=2,event_size_mean=1000)
        row.update(kwargs);return row
    def test_metrics_influence_score(self):
        base=confidence(self.row())['confidence_score']
        self.assertGreaterEqual(base,90)
        for changes in [dict(mapq30_support=0),dict(substitution_rate_mean=.08),dict(breakpoint1_sd=100),dict(high_quality_support=0)]:
            self.assertLess(confidence(self.row(**changes))['confidence_score'],base)
    def test_annotations_do_not_penalize(self):
        base=confidence(self.row())['confidence_score']
        self.assertEqual(confidence(self.row(mask_flags='centromere_satellite,tandem_repeat'))['confidence_score'],base)
    def test_unknown_and_low_support_caps(self):
        self.assertLessEqual(confidence(self.row(1))['confidence_score'],25)
        self.assertLessEqual(confidence(self.row(2))['confidence_score'],39)
        self.assertLessEqual(confidence(self.row(substitution_measured_segments=0,substitution_rate_mean=None))['confidence_score'],69)
        self.assertLessEqual(confidence(self.row(status='reference_boundary'))['confidence_score'],5)
    def test_sort_allowlist(self):
        with self.assertRaises(ValueError):load_page('/tmp/absent',sort_by='id;DROP TABLE calls')
