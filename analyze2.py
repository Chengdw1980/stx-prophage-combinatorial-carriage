import csv, re, sys, statistics as st
from collections import Counter, defaultdict

PATH = sys.argv[1]
csv.field_size_limit(10**9)
ASSIGNED = re.compile(r'^stx([12][a-z])_operon$')
UNASSIGNED = re.compile(r'^stx[12]_operon(=.+)?$')
MARKERS = ['eae', 'ehxA', 'espK', 'nleB', 'toxB', 'katP', 'espP', 'iha', 'saa', 'subA']
COMBINING = ['stx2a', 'stx2b', 'stx2c', 'stx2d']
SOLITARY = ['stx2e', 'stx2f', 'stx2i', 'stx2j', 'stx2k', 'stx2l', 'stx2n', 'stx2o']

rows = []
for row in csv.DictReader(open(PATH, newline='', encoding='utf-8', errors='replace'), delimiter='\t'):
    vg = row.get('virulence_genotypes') or ''
    if 'stx' not in vg: continue
    toks = [t.strip() for t in vg.split(',') if t.strip()]
    assigned, unassigned = set(), 0
    for t in toks:
        m = ASSIGNED.match(t)
        if m: assigned.add('stx' + m.group(1))
        elif UNASSIGNED.match(t): unassigned += 1
    if not assigned or unassigned: continue
    base = set(t.split('=')[0] for t in toks)
    def i(k):
        try: return int(row.get(k) or 0)
        except ValueError: return 0
    name = (row.get('scientific_name') or '') + ' ' + (row.get('serovar') or '')
    rows.append({'s2': set(s for s in assigned if s.startswith('stx2')),
                 's1': set(s for s in assigned if s.startswith('stx1')),
                 'mk': set(m for m in MARKERS if m in base),
                 'nc': i('asm_stats_n_contig'), 'nvf': i('number_virulence_genes'),
                 'len': i('asm_stats_length_bp'), 'o157': 'O157' in name})

def med(v): return st.median(v) if v else float('nan')

for label, excl in [('A: single stx2 AND no stx1  (what Methods says)', True),
                    ('B: single stx2, stx1 allowed (what Table S1 counts)', False)]:
    single = defaultdict(list)
    for r in rows:
        if len(r['s2']) == 1 and (not r['s1'] if excl else True):
            single[next(iter(r['s2']))].append(r)
    print('=' * 78)
    print(label)
    print('  %-7s %8s %7s %7s %7s %7s %7s' % ('subtype', 'n', 'medVF', 'eae%', 'ehxA%', 'katP%', 'iha%'))
    for s in COMBINING + ['stx2g'] + SOLITARY:
        g = single.get(s, [])
        if not g: continue
        p = lambda m: 100.0 * sum(m in r['mk'] for r in g) / len(g)
        print('  %-7s %8d %7.0f %7.1f %7.1f %7.1f %7.1f' % (
            s, len(g), med([r['nvf'] for r in g]), p('eae'), p('ehxA'), p('katP'), p('iha')))
    comb = [r for s in COMBINING for r in single.get(s, [])]
    soli = [r for s in SOLITARY for r in single.get(s, [])]
    print('  POOLED median VF: combining=%.0f (n=%d)  solitary=%.0f (n=%d)' % (
        med([r['nvf'] for r in comb]), len(comb), med([r['nvf'] for r in soli]), len(soli)))
    # stratified
    print('  --- within contig strata ---')
    for nm, lo, hi in [('<=50', 1, 50), ('51-200', 51, 200), ('>200', 201, 10**9)]:
        c = [r for r in comb if lo <= r['nc'] <= hi]
        s = [r for r in soli if lo <= r['nc'] <= hi]
        if len(c) < 10 or len(s) < 10: continue
        pc = lambda g, m: 100.0 * sum(m in r['mk'] for r in g) / len(g)
        print('   %-8s comb n=%6d med=%2.0f ehxA=%5.1f | soli n=%5d med=%2.0f ehxA=%4.1f | gap=%2.0f'
              % (nm, len(c), med([r['nvf'] for r in c]), pc(c, 'ehxA'),
                 len(s), med([r['nvf'] for r in s]), pc(s, 'ehxA'),
                 med([r['nvf'] for r in c]) - med([r['nvf'] for r in s])))
    print()

# ---- what does "closed" mean? O157, stx2a carriers, by contig cut ----
print('=' * 78)
print('O157 genomes carrying stx2a: co-carriage of stx2c and assembly length, by contig cut')
print('  %-10s %7s %7s %9s %14s %14s' % ('cut', 'n', 'both', 'a-only', 'med len both', 'med len a-only'))
for nm, lo, hi in [('==1', 1, 1), ('<=2', 1, 2), ('<=3', 1, 3), ('<=5', 1, 5)]:
    g = [r for r in rows if r['o157'] and 'stx2a' in r['s2'] and lo <= r['nc'] <= hi]
    both = [r for r in g if 'stx2c' in r['s2']]
    only = [r for r in g if 'stx2c' not in r['s2']]
    if not g: continue
    print('  %-10s %7d %7d %9d %14s %14s' % (
        nm, len(g), len(both), len(only),
        ('%d' % med([r['len'] for r in both])) if both else '-',
        ('%d' % med([r['len'] for r in only])) if only else '-'))
