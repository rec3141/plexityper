"""Fixtures for scratch/p6_field_count.py (Astra 251098919257686016, 251099116675186688)."""
import os, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p6_field_count as fc

CDS = "cds_id\tcontig\tstart\tend\nA\tc1\t101\t200\nB\tc1\t301\t400\nC\tc2\t1\t100\nD\tc2\t101\t200\n"

def run(sam, mapq=0):
    with tempfile.NamedTemporaryFile('w', suffix='.tsv', delete=False) as fh:
        fh.write(CDS); p = fh.name
    try:
        return fc.count(sam.strip().split('\n'), fc.load_cds(p), mapq)
    finally:
        os.unlink(p)

def rec(q, flag, ctg, pos, cig, mapq=60):
    return f"{q}\t{flag}\t{ctg}\t{pos}\t{mapq}\t{cig}\t*\t0\t0\t*\t*"

def test_30_vs_90():
    # SE read r1: 30 bp in A; SE read r2: 90 bp in B -> RNA 1:1, DNA bp 30:90
    sam = "\n".join([rec('r1', 0, 'c1', 171, '30M'), rec('r2', 0, 'c1', 311, '90M')])
    frag, bp, fate = run(sam)
    assert frag['A'] == 1 and frag['B'] == 1, frag
    assert bp['A'] == 30 and bp['B'] == 90, bp

def test_tie_dropped():
    # 50 bp in C (51..100) and 50 bp in D (101..150): exact tie -> unassigned
    frag, bp, fate = run(rec('t', 0, 'c2', 51, '100M'))
    assert not frag and fate['tie'] == 1, (frag, fate)

def test_strictly_greater():
    frag, bp, fate = run(rec('g', 0, 'c2', 52, '100M'))   # 49 in C, 51 in D
    assert frag == {'D': 1} and bp['D'] == 51, (frag, bp)

def test_flags_dropped():
    sam = "\n".join([rec('u', 4, '*', 0, '*', 0), rec('s', 256, 'c1', 101, '50M'),
                     rec('x', 2048, 'c1', 101, '50M')])
    frag, bp, fate = run(sam)
    assert not frag and fate['dropped'] == 3, (frag, fate)

def test_primary_plus_supplementary():
    # primary 20 bp in A, supplementary 80 bp in B: supplementary ignored -> A
    sam = "\n".join([rec('p', 0, 'c1', 181, '20M'), rec('p', 2048, 'c1', 301, '80M')])
    frag, bp, fate = run(sam)
    assert frag == {'A': 1} and bp['A'] == 20, frag

def test_pe_one_fragment_union():
    # mates overlap each other on A (101..160 and 141..200): union 101..200 = 100 bp, count once
    sam = "\n".join([rec('pe', 1 | 64, 'c1', 101, '60M'), rec('pe', 1 | 128, 'c1', 141, '60M')])
    frag, bp, fate = run(sam)
    assert frag == {'A': 1} and bp['A'] == 100, (frag, bp)

def test_pe_mates_on_two_cds():
    # mate1 60 bp in A, mate2 40 bp in B -> A wins (fragment once)
    sam = "\n".join([rec('m', 65, 'c1', 141, '60M'), rec('m', 129, 'c1', 301, '40M')])
    frag, bp, fate = run(sam)
    assert frag == {'A': 1} and bp['A'] == 60, frag

def test_se_vs_pe_same_count():
    se = rec('a', 0, 'c1', 101, '100M')
    pe = "\n".join([rec('b', 65, 'c1', 101, '50M'), rec('b', 129, 'c1', 151, '50M')])
    f1, b1, _ = run(se); f2, b2, _ = run(pe)
    assert f1 == f2 == {'A': 1} and b1['A'] == b2['A'] == 100

def test_deletion_not_counted():
    frag, bp, fate = run(rec('d', 0, 'c1', 101, '40M20D40M'))   # 80 aligned bases in A
    assert bp['A'] == 80, bp

def test_mapq20():
    sam = "\n".join([rec('lo', 0, 'c1', 101, '50M', 10), rec('hi', 0, 'c1', 301, '50M', 30)])
    f0, _, _ = run(sam, 0); f20, _, fate = run(sam, 20)
    assert f0 == {'A': 1, 'B': 1} and f20 == {'B': 1} and fate['dropped'] == 1

def test_noncontiguous_rejected():
    sam = "\n".join([rec('q', 65, 'c1', 101, '50M'), rec('z', 0, 'c1', 301, '50M'), rec('q', 129, 'c1', 151, '50M')])
    try:
        run(sam); raise AssertionError('should reject')
    except ValueError:
        pass

if __name__ == '__main__':
    n = 0
    for k, v in sorted(globals().items()):
        if k.startswith('test_'):
            v(); n += 1
    print(f'OK {n} fixtures')
