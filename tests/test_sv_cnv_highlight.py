import unittest
from robin.gui.components.sv_cnv import boundary_series
class BoundaryTests(unittest.TestCase):
    def test_chromosome_coordinates(self):
        row=dict(id=376,c1='chr7',p1=51339887,c2='chr7',p2=55967606,bin_width=100)
        lines=boundary_series(row,{'chr7':[1,2]},'chr7')['markArea']['data']
        self.assertEqual([r['xAxis'] for r in lines[0]],[51339887,55967606])
    def test_genome_offsets_match_bins_and_skip_alt_contigs(self):
        row=dict(id=1,c1='chr2',p1=25,c2='chr10',p2=50,bin_width=100)
        lines=boundary_series(row,{'chr10':[1],'chr2':[1,2],'chr1':[1,2,3],'chr1_alt':[0]*99},'All')['markArea']['data']
        self.assertEqual([(a['xAxis']+b['xAxis'])/2 for a,b in lines],[325,550])
        self.assertEqual(len(boundary_series(row,{'chr2':[1]},'chr2')['markArea']['data']),1)

    def test_translucency_and_no_opaque_edges(self):
        overlay=boundary_series(dict(id=1,c1='chr7',p1=100,c2='chr7',p2=200,bin_width=100),{'chr7':[2]},'chr7')
        self.assertNotIn('markLine',overlay)
        self.assertEqual(overlay['markArea']['itemStyle']['opacity'],.22)
        self.assertEqual(overlay['markArea']['itemStyle']['borderWidth'],0)

    def test_small_sv_midpoint_and_exact_10kb_boundary(self):
        row=dict(id=2,c1='chr7',c2='chr7',p1=100000,p2=109999,bin_width=1000)
        small=boundary_series(row,{'chr7':[2]*200},'chr7')
        self.assertTrue(small['midpoint_line'])
        endpoints=small['markArea']['data'][0]
        self.assertEqual([p['xAxis'] for p in endpoints],[104999.5,104999.5])
        large=boundary_series(dict(row,p2=110000),{'chr7':[2]*200},'chr7')
        self.assertFalse(large['midpoint_line'])
        self.assertEqual([p['xAxis'] for p in large['markArea']['data'][0]],[100000,110000])
