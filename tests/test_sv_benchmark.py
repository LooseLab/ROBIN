import csv
import json
import tempfile
import unittest
from pathlib import Path

import pysam
from robin.sv_replay.benchmark import IncrementalWorker, prepare, run, schedule


class BenchmarkTests(unittest.TestCase):
    def test_fifo_queue(self):
        self.assertEqual(schedule(5, 0, 2), (5, 7))
        self.assertEqual(schedule(6, 7, 2), (7, 9))
        self.assertEqual(schedule(20, 9, 1), (20, 21))

    def test_late_supplementary_and_repeat_batch(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            paths=[]
            for number,(chrom,cigar,flag) in enumerate([(0,'500M500H',0),(1,'500H500M',2048)]):
                path=root/f'{number}.bam';paths.append(path)
                with pysam.AlignmentFile(path,'wb',header={'SQ':[{'SN':'chr1','LN':10000000},{'SN':'chr2','LN':10000000}]}) as bam:
                    r=pysam.AlignedSegment();r.query_name='same';r.reference_id=chrom
                    r.reference_start=5000000;r.mapping_quality=60;r.flag=flag;r.cigarstring=cigar
                    r.query_sequence='A'*500;bam.write(r)
            worker=IncrementalWorker(root/'state.sqlite',{'chr1':10000000,'chr2':10000000})
            try:
                self.assertEqual(worker.process(paths[0],1)['updates'],[])
                second=worker.process(paths[1],2)
                self.assertEqual(len(second['updates']),1)
                self.assertEqual(second['updates'][0]['support'],1)
                self.assertEqual(second['updates'][0]['c2'],'chr2')
                self.assertEqual(worker.process(paths[1],3)['updates'],[])
                self.assertEqual(worker.db.execute('SELECT count(*) FROM evidence').fetchone()[0],1)
            finally:worker.db.close()

    def test_completion_manifest_and_run(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            bam=root/'batch.bam'
            with pysam.AlignmentFile(bam,'wb',header={'SQ':[{'SN':'chr1','LN':10000000}]}) as handle:
                for i in range(3):
                    r=pysam.AlignedSegment();r.query_name=f'r{i}';r.reference_id=0;r.reference_start=5000000
                    r.mapping_quality=60;r.cigarstring='200M80D200M';r.query_sequence='A'*400;handle.write(r)
            summary=root/'summary.tsv'
            summary.write_text('filename_bam\tstart_time\tduration\nbatch.bam\t0\t2\nbatch.bam\t3\t1\n')
            manifest=root/'manifest.csv';prepare(summary,[root],10,manifest)
            with manifest.open() as f:self.assertEqual(list(csv.DictReader(f))[0]['available_at'],'4.0')
            out=root/'modeled';run(manifest,out,'modeled',1)
            with (out/'calls.tsv').open() as f:
                call=list(csv.DictReader(f,delimiter='\t'))[0]
                self.assertEqual(call['support'],'3');self.assertEqual(call['status'],'provisional')
            self.assertEqual(json.loads((out/'summary.json').read_text())['alignment_records'],3)
            self.assertIn('p1',json.loads((out/'updates.jsonl').read_text())['updates'][0])
            with self.assertRaises(ValueError):run(manifest,out,'modeled',1)
            manifest.write_text('bam,available_at\n'+str(bam)+',0\n')
            paced=root/'paced';run(manifest,paced,'paced',1)
            report=json.loads((paced/'summary.json').read_text())
            self.assertEqual(report['max_observed_backlog'],0)
            self.assertGreaterEqual(report['maximum_batch_latency_seconds'],0)


if __name__=='__main__':unittest.main()
