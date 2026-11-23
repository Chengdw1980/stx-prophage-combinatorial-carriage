import csv, re, sys, json, math, statistics as st
from collections import Counter, defaultdict

PATH = sys.argv[1]
csv.field_size_limit(10**9)
A = re.compile(r'^stx([12][a-z])_operon$')
U = re.compile(r'^stx[12]_operon(=.+)?$')
MARKERS = ['eae', 'ehxA', 'espK', 'nleB', 'toxB', 'katP', 'espP', 'iha', 'saa', 'subA']
COMB = ['stx2a', 'stx2b', 'stx2c', 'stx2d']
INTER = ['stx2g']
SOLI = ['stx2e', 'stx2f', 'stx2i', 'stx2j', 'stx2k', 'stx2l', 'stx2n', 'stx2o']

rows = []
n_pos = 0
for r in csv.DictReader(open(PATH, newline='', encoding='utf-8', errors='replace'), delimiter='\t'):
    vg = r.get('virulence_genotypes') or ''
    if 'stx' not in vg: continue
    toks = [t.strip() for t in vg.split(',') if t.strip()]
    asg, un = set(), 0
    for t in toks:
        m = A.match(t)
        if m: asg.add('stx' + m.group(1))
        elif U.match(t): un += 1
    if not asg and un == 0: continue
    n_pos += 1
    if not asg or un: continue
    base = set(t.split('=')[0] for t in toks)
    def I(k):
        try: return int(r.get(k) or 0)
        except ValueError: return 0
    sname = (r.get('scientific_name') or '').strip() or 'NULL'
    serov = (r.get('serovar') or '').strip()
    rows.append({
        's1': frozenset(s for s in asg if s.startswith('stx1')),
        's2': frozenset(s for s in asg if s.startswith('stx2')),
        'mk': frozenset(m for m in MARKERS if m in set(toks)),
        'nc': I('asm_stats_n_contig'), 'vf': I('number_virulence_genes'),
        'len': I('asm_stats_length_bp'), 'lvl': (r.get('asm_level') or '').strip(),
        'sn': sname, 'o157': ('O157' in sname or 'O157' in serov),
    })

N = len(rows)
out = {'n_stx_positive': n_pos, 'n_resolved': N}

def wald(k, n):
    if n == 0: return (float('nan'),) * 3
    p = k / n
    h = 1.96 * math.sqrt(p * (1 - p) / n)
    return 100 * p, max(0.0, 100 * (p - h)), min(100.0, 100 * (p + h))

def med(v): return st.median(v) if v else float('nan')

# ---------- Table 1 ----------
cnt = Counter()
for r in rows:
    for s in (r['s1'] | r['s2']): cnt[s] += 1
t1 = [(s, c, 100 * c / N) for s, c in sorted(cnt.items(), key=lambda kv: -kv[1])]
out['table1'] = t1

# ---------- Table S1 / Figure 1 : single stx2 subtype, does it carry stx1 ----------
only2 = defaultdict(list)
for r in rows:
    if len(r['s2']) == 1: only2[next(iter(r['s2']))].append(r)
s1tab = []
for s, g in only2.items():
    k = sum(1 for r in g if r['s1'])
    p, lo, hi = wald(k, len(g))
    s1tab.append((s, len(g), k, p, lo, hi))
s1tab.sort(key=lambda x: -x[1])
out['tableS1'] = s1tab

# ---------- Table S2 / Figure 2 ----------
s2tab = []
for s in COMB + INTER + SOLI:
    g = only2.get(s, [])
    if not g: continue
    row = {'sub': s, 'n': len(g), 'medVF': med([r['vf'] for r in g])}
    for m in MARKERS:
        row[m] = 100 * sum(m in r['mk'] for r in g) / len(g)
    s2tab.append(row)
out['tableS2'] = s2tab
combg = [r for s in COMB for r in only2.get(s, [])]
solig = [r for s in SOLI for r in only2.get(s, [])]
out['pooled'] = {
    'comb_med': med([r['vf'] for r in combg]), 'comb_n': len(combg),
    'soli_med': med([r['vf'] for r in solig]), 'soli_n': len(solig),
    'soli_ehxA': 100 * sum('ehxA' in r['mk'] for r in solig) / len(solig),
    'soli_katP': 100 * sum('katP' in r['mk'] for r in solig) / len(solig),
}

# ---------- stratified robustness (Figure 1 separation + Figure 2 gap) ----------
def stratum(n):
    if n <= 50: return '<=50'
    if n <= 200: return '51-200'
    return '>200'
strat = {}
for s in ['stx2e', 'stx2j', 'stx2k', 'stx2n', 'stx2a', 'stx2c']:
    d = {}
    for st_ in ['<=50', '51-200', '>200']:
        g = [r for r in only2.get(s, []) if stratum(r['nc']) == st_]
        if g: d[st_] = (len(g), 100 * sum(1 for r in g if r['s1']) / len(g))
    strat[s] = d
out['strat_partnership'] = strat
gap = {}
for st_ in ['<=50', '51-200', '>200']:
    c = [r for r in combg if stratum(r['nc']) == st_]
    s_ = [r for r in solig if stratum(r['nc']) == st_]
    gap[st_] = {'cn': len(c), 'cmed': med([r['vf'] for r in c]),
                'sn': len(s_), 'smed': med([r['vf'] for r in s_]),
                'cehxA': 100 * sum('ehxA' in r['mk'] for r in c) / len(c),
                'sehxA': 100 * sum('ehxA' in r['mk'] for r in s_) / len(s_),
                'ckatP': 100 * sum('katP' in r['mk'] for r in c) / len(c),
                'skatP': 100 * sum('katP' in r['mk'] for r in s_) / len(s_)}
out['strat_accessory'] = gap

# ---------- Table S4 species ----------
s4 = []
for s in COMB + INTER + SOLI:
    g = only2.get(s, [])
    if not g: continue
    c = Counter(r['sn'] for r in g)
    top = [(nm, 100 * k / len(g)) for nm, k in c.most_common(3)]
    ecoli = 100 * sum(1 for r in g if r['sn'].startswith('Escherichia coli')) / len(g)
    other = [(nm, 100 * k / len(g)) for nm, k in c.most_common()
             if not nm.startswith('Escherichia coli')]
    s4.append({'sub': s, 'n': len(g), 'top': top, 'ecoli_all': ecoli,
               'top_other': other[0] if other else None})
out['tableS4'] = s4

# ---------- Table 2 / S5 pairing ----------
pairset = [r for r in rows if len(r['s1']) == 1 and len(r['s2']) == 1]
M = len(pairset)
m1 = Counter(next(iter(r['s1'])) for r in pairset)
m2 = Counter(next(iter(r['s2'])) for r in pairset)
obs = Counter((next(iter(r['s1'])), next(iter(r['s2']))) for r in pairset)
s5 = []
for (a, b), o in sorted(obs.items(), key=lambda kv: (kv[0][0], -kv[1])):
    e = m1[a] * m2[b] / M
    s5.append({'s1': a, 's2': b, 'obs': o, 'exp': e, 'oe': o / e,
               'pct': 100 * o / m1[a]})
out['tableS5'] = s5
out['pair_total'] = M
out['pair_marg1'] = dict(m1)

# ---------- Table S3 / Figure 3 ----------
def cbin(n):
    if n <= 5: return '1-5'
    if n <= 50: return '6-50'
    if n <= 100: return '51-100'
    if n <= 200: return '101-200'
    if n <= 500: return '201-500'
    return '>500'
BINS = ['1-5', '6-50', '51-100', '101-200', '201-500', '>500']
s3 = []
for lin, f in [('O157', lambda r: r['o157']), ('non-O157', lambda r: not r['o157'])]:
    for b in BINS:
        g = [r for r in rows if f(r) and 'stx2a' in r['s2'] and cbin(r['nc']) == b]
        if not g: continue
        k = sum(1 for r in g if 'stx2c' in r['s2'])
        p, lo, hi = wald(k, len(g))
        s3.append({'lin': lin, 'bin': b, 'n': len(g), 'k': k, 'p': p, 'lo': lo, 'hi': hi})
out['tableS3'] = s3
agg = {}
for lin, f in [('O157', lambda r: r['o157']), ('non-O157', lambda r: not r['o157'])]:
    g = [r for r in rows if f(r) and 'stx2a' in r['s2'] and r['nc'] > 100]
    k = sum(1 for r in g if 'stx2c' in r['s2'])
    agg[lin] = {'n': len(g), 'k': k, 'p': 100 * k / len(g)}
out['gt100'] = agg

# ---------- closed-genome assembly length ----------
closed = [r for r in rows if r['o157'] and r['lvl'] == 'Complete Genome']
both = [r['len'] for r in closed if {'stx2a', 'stx2c'} <= set(r['s2'])]
alone = [r['len'] for r in closed if set(r['s2']) == {'stx2a'} and not r['s1']]
n1, n2 = len(both), len(alone)
Uv = sum((1.0 if x > y else 0.5 if x == y else 0.0) for x in both for y in alone)
mu = n1 * n2 / 2.0
allv = both + alone
tc = sum(t ** 3 - t for t in Counter(allv).values())
n = n1 + n2
sig = math.sqrt(n1 * n2 / 12.0 * ((n + 1) - tc / (n * (n - 1.0))))
z = (Uv - mu) / sig
p = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))
out['closed'] = {'n_both': n1, 'med_both': med(both), 'n_alone': n2,
                 'med_alone': med(alone), 'diff': med(both) - med(alone),
                 'U': Uv, 'z': z, 'p': p,
                 'distinct_both': len(set(both)) == len(both)}

json.dump(out, open('results.json', 'w'), indent=1, default=str)
print(json.dumps(out, indent=1, default=str)[:400])
print('...\nwrote results.json')
