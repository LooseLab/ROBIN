import json
from pathlib import Path
import sqlite3
import tempfile
from types import SimpleNamespace
import unittest

import pysam
from robin.sv_replay.full_run import read_manifest, publish
from robin.analysis.sv_analysis import sv_handler


class FullReplayTests(unittest.TestCase):
    def test_manifest_and_atomic_publication(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);source=root/'bam_pass';source.mkdir()
            bam=source/'a.bam';bam.write_bytes(b'original')
            manifest=root/'input.csv';manifest.write_text(f'bam,available_at\n{bam},2\n')
            self.assertEqual(read_manifest(manifest)[0]['available_at'],2)
            incoming=root/'incoming';incoming.mkdir()
            published=publish(bam,incoming)
            self.assertTrue(published.is_symlink());self.assertEqual(published.read_bytes(),b'original')
            self.assertFalse(list(incoming.rglob('*.pending')))
            manifest.write_text(f'bam,available_at\n{bam},nan\n')
            with self.assertRaises(ValueError):read_manifest(manifest)

    def test_native_handler_retries_do_not_duplicate_support(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);bam=root/'source.bam'
            with pysam.AlignmentFile(bam,'wb',header={'SQ':[{'SN':'chr1','LN':10000000}]}) as handle:
                for i in range(3):
                    r=pysam.AlignedSegment();r.query_name=f'r{i}';r.reference_id=0;r.reference_start=5000000
                    r.mapping_quality=60;r.cigarstring='200M80D200M';r.query_sequence='A'*400;handle.write(r)
            results={}
            ctx=SimpleNamespace(metadata={'work_dir':str(root/'output')},filepath=str(bam),
                  get_sample_id=lambda:'sample',add_result=lambda key,value:results.update({key:value}))
            job=SimpleNamespace(context=ctx)
            sv_handler(job);sv_handler(job)
            db=sqlite3.connect(root/'output/sample/sv_live/state.sqlite')
            self.assertEqual(db.execute('SELECT count(*) FROM evidence').fetchone()[0],3);db.close()
            with (root/'output/sample/sv_live/updates.jsonl').open() as handle:
                updates=[json.loads(l) for l in handle]
            self.assertEqual(updates[0]['updates'][0]['support'],3)
            self.assertEqual(updates[1]['updates'],[])
            self.assertEqual(results['sv']['status'],'success')


if __name__=='__main__':unittest.main()
