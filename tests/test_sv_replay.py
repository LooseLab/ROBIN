"""Exercise extraction, availability, molecule counting and checkpoint recovery."""
import csv
import sqlite3
from pathlib import Path
import tempfile
import unittest

import pysam

from robin.sv_replay.replay import main as replay_main, segment, signatures, connect


def main(argv):
    # Tiny synthetic chromosomes would all be inside the production 1 Mb margin.
    return replay_main(argv + ['--chromosome-end-margin', '0'])


def test_reverse_hard_clips_and_query_adjacency():
    a = segment('chr1', 1000, '+', '500M500H', 60)
    b = segment('chr2', 2000, '-', '500M500H', 60)
    assert (b['qs'], b['qe']) == (500, 1000)
    sigs = signatures([a, b, a])
    assert sigs == [dict(a=('chr1', 1500, '+'), b=('chr2', 2500, '+'), kind='BND', size=0, anchors=[500,500])]
    assert signatures([a, dict(b, mapq=0)]) == []
    assert signatures([a, dict(b, qs=300)]) == []


def test_cigar_indels_and_anchors():
    sigs = signatures([segment('chr1', 100, '+', '200M80D200M60I200M', 60)])
    assert [(s['kind'], s['a'][1], s['b'][1], s['size']) for s in sigs] == [
        ('DEL', 300, 380, 80), ('INS', 580, 580, 60)]
    assert signatures([segment('chr1', 0, '+', '20M100D300M', 60)]) == []
    assert signatures([segment('chr1', 0, '+', '200M100N300M', 60)]) == []


def write_bam(path, ids, start=100):
    with pysam.AlignmentFile(path, 'wb', header={
        'HD': {'VN': '1.6'}, 'SQ': [{'SN': 'chr1', 'LN': 100000}]}) as bam:
        for rid in ids:
            r = pysam.AlignedSegment()
            r.query_name = rid
            r.query_sequence = 'A' * 400
            r.reference_id = 0
            r.reference_start = start
            r.mapping_quality = 60
            r.cigarstring = '200M80D200M'
            bam.write(r)


def rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def test_end_to_end_order_duplicates_resume_and_config(tmp_path):
    bam = tmp_path / 'reads.bam'
    write_bam(bam, ['r3', 'r1', 'r2', 'r1'])
    times = tmp_path / 'times.csv'
    times.write_text('read_id,available_at\nr1,10\nr2,20\nr3,30\n')
    out = tmp_path / 'output'
    argv = ['--bam', str(bam), '--read-times', str(times), '--output', str(out)]
    main(argv)
    calls = rows(out / 'calls.tsv')
    assert len(calls) == 1
    assert calls[0]['support'] == '3'
    assert calls[0]['first'] == '10.0'
    assert calls[0]['status'] == 'provisional'
    assert [r['support'] for r in rows(out / 'history.tsv')] == ['1', '2', '3']
    original = {p.name: p.read_bytes() for p in out.glob('*.tsv')}
    main(argv)
    assert original == {p.name: p.read_bytes() for p in out.glob('*.tsv')}
    # Simulate stopping after the first read's committed transaction.
    db = sqlite3.connect(out / 'replay.sqlite')
    with db:
        db.execute("DELETE FROM processed WHERE read_id!='r1'")
        db.execute("DELETE FROM evidence WHERE read_id!='r1'")
        db.execute('DELETE FROM history WHERE at>10')
        db.execute('UPDATE calls SET last=10')
    db.close()
    main(argv)
    assert original == {p.name: p.read_bytes() for p in out.glob('*.tsv')}
    with unittest.TestCase().assertRaisesRegex(ValueError, 'settings differ'):
        main(argv + ['--min-support', '4'])


def test_both_breakends_and_orientation_bound_clustering(tmp_path):
    from robin.sv_replay.replay import add_signature, canonical
    db = connect(tmp_path / 'test.sqlite')
    a = canonical(('chr1', 1000, '+'), ('chr2', 2000, '-'), 'BND')
    b = canonical(('chr1', 1010, '+'), ('chr2', 3000, '-'), 'BND')
    c = canonical(('chr1', 1010, '-'), ('chr2', 2000, '-'), 'BND')
    with db:
        ids = [add_signature(db, str(i), i, s, 100) for i, s in enumerate([a,b,c])]
    assert len(set(ids)) == 3
    db.close()


def test_missing_times_roll_back_and_late_batch_records(tmp_path):
    bam = tmp_path / 'one.bam'
    later = tmp_path / 'two.bam'
    write_bam(bam, ['r1'])
    write_bam(later, ['r1', 'r2'])
    times = tmp_path / 'times.csv'
    times.write_text('read_id,available_at\nr1,0\n')
    out = tmp_path / 'bad'
    with unittest.TestCase().assertRaisesRegex(ValueError, 'Missing read availability'):
        main(['--bam', str(later), '--read-times', str(times), '--output', str(out)])
    db = sqlite3.connect(out / 'replay.sqlite')
    assert db.execute('SELECT count(*) FROM reads').fetchone()[0] == 0
    db.close()
    out = tmp_path / 'batches'
    main(['--bam', str(bam), '--bam', str(later), '--output', str(out)])
    assert all(r['at'] == '60.0' for r in rows(out / 'evidence.tsv'))


def test_sa_and_explicit_supplementary_count_once(tmp_path):
    bam = tmp_path / 'split.bam'
    header = {'HD': {'VN': '1.6'}, 'SQ': [
        {'SN': 'chr1', 'LN': 100000}, {'SN': 'chr2', 'LN': 100000}]}
    with pysam.AlignmentFile(bam, 'wb', header=header) as handle:
        for i in range(3):
            r = pysam.AlignedSegment()
            r.query_name = f'r{i}'
            r.query_sequence = 'A' * 1000
            r.reference_id = 0
            r.reference_start = 1000
            r.mapping_quality = 60
            r.cigarstring = '500M500S'
            r.set_tag('SA', 'chr2,2001,+,500S500M,60,0;')
            handle.write(r)
            r = pysam.AlignedSegment()
            r.query_name = f'r{i}'
            r.query_sequence = 'A' * 500
            r.reference_id = 1
            r.reference_start = 2000
            r.mapping_quality = 60
            r.flag = 2048
            r.cigarstring = '500H500M'
            handle.write(r)
    out = tmp_path / 'output'
    main(['--bam', str(bam), '--output', str(out)])
    calls = rows(out / 'calls.tsv')
    assert len(calls) == 1
    assert calls[0]['kind'] == 'BND'
    assert calls[0]['support'] == '3'


def test_mask_boundaries_nested_intervals_and_mitochondria(tmp_path):
    from robin.sv_replay.regions import IntervalIndex, RegionContext, call_status
    from robin.sv_replay.replay import canonical
    bed = tmp_path / 'repeat.bed'
    bed.write_text('chr1\t100\t200\nchr1\t120\t130\n')
    empty = tmp_path / 'empty.bed'
    empty.write_text('')
    index = IntervalIndex(bed)
    assert index.contains('chr1', 199)
    assert not index.contains('chr1', 200)
    context = RegionContext({'chr1': 1000, 'chrM': 16569}, bed, empty,
                            end_margin=10, origin_margin=5)
    at_end = canonical(('chr1', 200, '+'), ('chr1', 500, '-'), 'DEL')
    assert context.flags(at_end) == {'tandem_repeat'}
    outside = canonical(('chr1', 200, '-'), ('chr1', 500, '+'), 'BND')
    assert not context.flags(outside)
    spanning = canonical(('chr1', 50, '+'), ('chr1', 500, '-'), 'DEL')
    assert not context.flags(spanning)  # Crossing a repeat does not mask its unique ends.
    mt = canonical(('chrM', 0, '-'), ('chrM', 16569, '+'), 'BND')
    assert context.flags(mt) == {'mitochondrial_reference_boundary'}
    deletion = canonical(('chrM', 0, '+'), ('chrM', 16569, '-'), 'DEL')
    assert not context.flags(deletion)
    assert call_status(1, context.flags(mt), 3) == 'reference_boundary'
    assert call_status(3, {'tandem_repeat'}, 3) == 'masked_candidate'
    assert call_status(2, {'tandem_repeat'}, 3) == 'candidate'
    assert call_status(3, set(), 3) == 'provisional'


def test_masked_evidence_retained_history_and_annotation_checkpoint(tmp_path):
    bam = tmp_path / 'test.bam'
    write_bam(bam, ['r1', 'r2', 'r3'])
    bed = tmp_path / 'mask.bed'
    bed.write_text('chr1\t299\t300\n')
    out = tmp_path / 'output'
    argv = ['--bam', str(bam), '--output', str(out), '--mask-bed', str(bed)]
    main(argv)
    call = rows(out / 'calls.tsv')[0]
    assert call['status'] == 'masked_candidate'
    assert call['mask_flags'] == 'custom_mask'
    assert call['support'] == '3'
    assert len(rows(out / 'evidence.tsv')) == 3
    assert rows(out / 'history.tsv')[-1]['status'] == 'masked_candidate'
    assert rows(out / 'history.tsv')[-1]['mask_flags'] == 'custom_mask'
    before = (out / 'calls.tsv').read_bytes()
    main(argv)
    assert (out / 'calls.tsv').read_bytes() == before
    bed.write_text('chr1\t500\t501\n')
    with unittest.TestCase().assertRaisesRegex(ValueError, 'settings differ'):
        main(argv)


def test_later_masked_signature_retracts_promotion(tmp_path):
    first, later = tmp_path / 'first.bam', tmp_path / 'later.bam'
    write_bam(first, ['r1', 'r2', 'r3'])
    write_bam(later, ['r4'], start=150)
    bed = tmp_path / 'mask.bed'
    bed.write_text('chr1\t349\t350\n')
    times = tmp_path / 'times.csv'
    times.write_text('read_id,available_at\nr1,1\nr2,2\nr3,3\nr4,4\n')
    out = tmp_path / 'output'
    main(['--bam', str(first), '--bam', str(later), '--mask-bed', str(bed),
          '--read-times', str(times), '--output', str(out)])
    call = rows(out / 'calls.tsv')[0]
    assert call['support'] == '4'
    assert call['status'] == 'masked_candidate'
    history = rows(out / 'history.tsv')
    assert [r['status'] for r in history] == [
        'candidate', 'candidate', 'provisional', 'masked_candidate']
    assert history[2]['mask_flags'] == ''
    assert history[3]['mask_flags'] == 'custom_mask'


def load_tests(loader, tests, pattern):
    suite = unittest.TestSuite()
    for name, fn in sorted(globals().items()):
        if name.startswith('test_') and callable(fn):
            def run(fn=fn):
                with tempfile.TemporaryDirectory() as directory:
                    if fn.__code__.co_argcount:
                        fn(Path(directory))
                    else:
                        fn()
            suite.addTest(unittest.FunctionTestCase(run, description=name))
    return suite


if __name__ == '__main__':
    unittest.main()
