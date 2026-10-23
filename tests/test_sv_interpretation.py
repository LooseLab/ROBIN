import unittest
from robin.sv_replay.interpretation import annotate

class InterpretationTests(unittest.TestCase):
    def call(self, a, b, c2='chr9', p1=100, p2=200, kind='BND'):
        return annotate(dict(kind=kind,c1='chr9',p1=p1,o1=a,c2=c2,p2=p2,o2=b))

    def test_topology_and_reversed_order(self):
        for a,b,expected in [('+','+','inversion_like'),('-','-','inversion_like'),
                             ('-','+','duplication_like'),('+','-','deletion_like')]:
            self.assertEqual(self.call(a,b),dict(topology=expected,interval_bp=100))
            self.assertEqual(self.call(b,a,p1=200,p2=100),dict(topology=expected,interval_bp=100))

    def test_no_overinterpretation(self):
        self.assertEqual(self.call('+','+',c2='chr10'),dict(topology='interchromosomal',interval_bp=None))
        self.assertEqual(self.call('-','+',p2=100)['topology'],'unresolved')
        self.assertEqual(self.call('+','-',kind='DEL')['topology'],'del')
        self.assertEqual(self.call('+','+',kind='INS')['topology'],'ins')
