import sqlite3
import tempfile
import unittest
from pathlib import Path

from robin.gui.components.sv_data import load_page, load_evidence


class SVViewerTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name)
        folder=self.root/'sv_live';folder.mkdir();self.path=folder/'state.sqlite'
        db=sqlite3.connect(self.path)
        db.executescript('''
          CREATE TABLE calls(id INTEGER PRIMARY KEY,kind TEXT,size INTEGER,c1 TEXT,p1 INT,o1 TEXT,c2 TEXT,p2 INT,o2 TEXT);
          CREATE TABLE delivered(id TEXT PRIMARY KEY,support INT,status TEXT,mask_flags TEXT);
          CREATE TABLE evidence(call_id INT,read_id TEXT,p1 INT,p2 INT);
        ''')
        db.executemany('INSERT INTO calls VALUES (?,?,?,?,?,?,?,?,?)',[
            (1,'DEL',80,'chr1',100,'+','chr1',180,'-'),
            (2,'BND',0,'chr9',200,'-','chr22',300,'+'),
            (3,'INS',60,'chr1',500,'+','chr1',500,'-')])
        db.executemany('INSERT INTO delivered VALUES (?,?,?,?)',[
            ('1',3,'provisional',''),('2',5,'masked_candidate','centromere_satellite,tandem_repeat'),
            ('3',1,'candidate','')])
        db.executemany('INSERT INTO evidence VALUES (?,?,?,?)',[(1,f'r{i}',100+i,180+i) for i in range(3)])
        db.commit();db.close()

    def tearDown(self):self.temp.cleanup()

    def test_default_filters_pagination_and_read_only(self):
        before=self.path.read_bytes()
        data=load_page(self.root,page_size=1)
        self.assertEqual(data['total'],2);self.assertEqual(data['pages'],2)
        self.assertEqual(data['rows'][0]['id'],2)
        self.assertEqual(load_page(self.root,page=99,page_size=1)['rows'][0]['id'],1)
        self.assertEqual(load_page(self.root,status='candidate',min_support=0)['rows'][0]['id'],3)
        self.assertEqual(load_page(self.root,mask='unmasked')['rows'][0]['id'],1)
        self.assertEqual(load_page(self.root,mask='tandem_repeat')['rows'][0]['id'],2)
        self.assertEqual(before,self.path.read_bytes())

    def test_numeric_sort_across_pages(self):
        self.assertEqual(load_page(self.root,sort_by='support',descending=False,page_size=1)['rows'][0]['id'],1)
        self.assertEqual(load_page(self.root,sort_by='size',descending=True,page_size=1)['rows'][0]['id'],1)
        self.assertEqual(load_page(self.root,sort_by='size',descending=True,page_size=1,page=2)['rows'][0]['id'],2)
        self.assertIsNone(load_page(self.root)['rows'][0]['high_quality_support'])

    def test_mapq_filter_boundary_unknown_and_pagination(self):
        import json
        db=sqlite3.connect(self.path)
        db.execute('CREATE TABLE mapping_summary(call_id INTEGER PRIMARY KEY,data TEXT)')
        db.executemany('INSERT INTO mapping_summary VALUES (?,?)',[(1,json.dumps({'mean_min_mapq':55})),(2,json.dumps({'mean_min_mapq':54.99}))])
        db.commit();db.close()
        page=load_page(self.root,page_size=1)
        self.assertEqual(page['total'],1);self.assertEqual(page['rows'][0]['id'],1)
        self.assertEqual(load_page(self.root,show_low_mapq=True)['total'],2)
        self.assertEqual(load_page(self.root,search='chr22')['total'],0)
        self.assertEqual(load_page(self.root,min_support=0)['total'],2)  # unknown MAPQ retained
        self.assertEqual(load_page(self.root,show_low_mapq=True,sort_by='support')['rows'][0]['id'],2)

    def test_both_breakends_search_and_query_safety(self):
        self.assertEqual(load_page(self.root,search='chr22')['rows'][0]['id'],2)
        self.assertEqual(load_page(self.root,search='1')['rows'][0]['id'],1)

    def test_evidence_missing_data_and_validation(self):
        data=load_evidence(self.root,1,limit=2)
        self.assertEqual(data['total'],3);self.assertEqual(len(data['rows']),2)
        self.assertEqual(load_page(self.root/'missing')['available'],False)
        self.assertFalse((self.root/'missing').exists())
        self.assertEqual(load_page(self.root,search="' OR 1=1 --")['total'],0)
        with self.assertRaises(ValueError):load_page(self.root,status='bad')
        with self.assertRaises(ValueError):load_page(self.root,mask='bad')
        with self.assertRaises(ValueError):load_page(self.root,page_size=100000)


if __name__=='__main__':unittest.main()
