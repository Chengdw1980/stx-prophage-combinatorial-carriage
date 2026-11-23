import csv, re, sys, statistics as st
from collections import Counter, defaultdict

PATH = sys.argv[1]
csv.field_size_limit(10**9)

ASSIGNED = re.compile(r'^stx([12][a-z])_operon$')
UNASSIGNED = re.compile(r'^stx[12]_operon(=.+)?$')
MARKERS = ['eae', 'ehxA', 'espK', 'nleB', 'toxB', 'katP', 'espP', 'iha', 'saa', 'subA']

COMBINING = ['stx2a', 'stx2b', 'stx2c', 'stx2d']
SOLITARY = ['stx2e', 'stx2f', 'stx2i', 'stx2j', 'stx2k', 'stx2l', 'stx2n', 'stx2o']
INTERMEDIATE = ['stx2g']

def bin_contigs(n):
    if n is None: return None
    if n <= 5: return '1-5'
    if n <= 50: return '6-50'
    if n <= 100: return '51-100'
    if n <= 200: return '101-200'
    if n <= 500: return '201-500'
    return '>500'

COARSE = [('<=50', lambda n: n <= 50),
          ('51-200', lambda n: 50 < n <= 200),
          ('>200', lambda n: n > 200)]

def coarse(n):
    if n is None: return None
    for name, f in COARSE:
        if f(n): return name
    return None

rows = []
n_stx = 0
f = open(PATH, newline='', encoding='utf-8', errors='replace')
for row in csv.DictReader(f, delimiter='\t'):
    vg = row.get('virulence_genotypes') or ''
    if 'stx' not in vg:
        continue
    toks = [t.strip() for t in vg.split(',') if t.strip()]
    assigned, unassigned = set(), 0
    for t in toks:
        m = ASSIGNED.match(t)
        if m:
            assigned.add('stx' + m.group(1))
        elif UNASSIGNED.match(t):
            unassigned += 1
    if not assigned and unassigned == 0:
        continue          # stxA/stxB subunit hits only, no operon call
    n_stx += 1
    if not assigned or unassigned:
        continue          # not fully resolved
    # markers: base name present, ignoring =STATUS suffix
    base = set(t.split('=')[0] for t in toks)
    try:
        nc = int(row.get('asm_stats_n_contig') or 0) or None
    except ValueError:
        nc = None
    try:
        nvf = int(row.get('number_virulence_genes') or -1)
    except ValueError:
        nvf = -1
    try:
        alen = int(row.get('asm_stats_length_bp') or 0) or None
    except ValueError:
        alen = None
    name = (row.get('scientific_name') or '') + ' ' + (row.get('serovar') or '')
    rows.append({
        'sub': assigned,
        'stx1': set(s for s in assigned if s.startswith('stx1')),
        'stx2': set(s for s in assigned if s.startswith('stx2')),
        'mk': set(m for m in MARKERS if m in base),
        'nc': nc, 'nvf': nvf, 'len': alen,
        'o157': 'O157' in name,
    })

print('stx-positive records (operon called): %d' % n_stx)
print('fully resolved genomes:              %d' % len(rows))
print()

# ---- Table 1 reproduction ----
cnt = Counter()
for r in rows:
    for s in r['sub']:
        cnt[s] += 1
print('--- Subtype prevalence among fully resolved (cf. Table 1) ---')
for s, c in sorted(cnt.items(), key=lambda kv: -kv[1]):
    print('  %-7s %7d  %5.2f%%' % (s, c, 100.0 * c / len(rows)))
print()

# ---- single-stx2, no stx1 set (manuscript's accessory-content set) ----
single = defaultdict(list)
for r in rows:
    if len(r['stx2']) == 1 and not r['stx1']:
        single[next(iter(r['stx2']))].append(r)

def med(v):
    return st.median(v) if v else float('nan')

print('--- Accessory content, single stx2 subtype and no stx1 (cf. Table S2 / Figure 2) ---')
print('  %-7s %7s %8s %8s %8s' % ('subtype', 'n', 'medVF', 'ehxA%', 'katP%'))
for s in COMBINING + INTERMEDIATE + SOLITARY:
    g = single.get(s, [])
    if not g: continue
    vf = [r['nvf'] for r in g if r['nvf'] >= 0]
    print('  %-7s %7d %8.0f %8.1f %8.1f' % (
        s, len(g), med(vf),
        100.0 * sum('ehxA' in r['mk'] for r in g) / len(g),
        100.0 * sum('katP' in r['mk'] for r in g) / len(g)))
print()

comb = [r for s in COMBINING for r in single.get(s, [])]
soli = [r for s in SOLITARY for r in single.get(s, [])]

def vfs(g): return [r['nvf'] for r in g if r['nvf'] >= 0]
print('POOLED medians  combining=%.1f (n=%d)   solitary=%.1f (n=%d)' % (
    med(vfs(comb)), len(comb), med(vfs(soli)), len(soli)))
print()

# ---- assembly quality distribution of the two classes ----
print('--- Assembly contiguity of the two classes (the confounder) ---')
print('  %-12s %8s %8s %8s %8s' % ('stratum', 'comb n', 'comb %', 'soli n', 'soli %'))
for name, _ in COARSE:
    cn = sum(1 for r in comb if coarse(r['nc']) == name)
    sn = sum(1 for r in soli if coarse(r['nc']) == name)
    print('  %-12s %8d %7.1f%% %8d %7.1f%%' % (
        name, cn, 100.0 * cn / len(comb), sn, 100.0 * sn / len(soli)))
print()

# ---- THE TEST: accessory content within contig strata ----
print('=== Accessory virulence gene content WITHIN assembly-quality strata ===')
print('  %-10s %9s %7s %9s %7s %8s' % ('stratum', 'comb n', 'medVF', 'soli n', 'medVF', 'gap'))
for name, _ in COARSE:
    c = [r for r in comb if coarse(r['nc']) == name]
    s = [r for r in soli if coarse(r['nc']) == name]
    if len(c) < 10 or len(s) < 10:
        print('  %-10s   (insufficient: comb=%d soli=%d)' % (name, len(c), len(s)))
        continue
    mc, ms = med(vfs(c)), med(vfs(s))
    print('  %-10s %9d %7.0f %9d %7.0f %8.0f' % (name, len(c), mc, len(s), ms, mc - ms))
print()

print('=== Marker prevalence WITHIN strata (ehxA / katP) ===')
print('  %-10s %18s %18s' % ('stratum', 'combining', 'solitary'))
for name, _ in COARSE:
    c = [r for r in comb if coarse(r['nc']) == name]
    s = [r for r in soli if coarse(r['nc']) == name]
    if len(c) < 10 or len(s) < 10: continue
    f1 = lambda g, m: 100.0 * sum(m in r['mk'] for r in g) / len(g)
    print('  %-10s   ehxA %5.1f katP %5.1f   ehxA %5.1f katP %5.1f' % (
        name, f1(c, 'ehxA'), f1(c, 'katP'), f1(s, 'ehxA'), f1(s, 'katP')))
