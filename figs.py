import json, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

d = json.load(open('/home/claude/results.json'))
OUT = '/mnt/user-data/outputs/'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9,
                     'axes.linewidth': 0.8, 'savefig.bbox': 'tight'})
BLUE, GREEN, ORANGE = '#1f77b4', '#2a9d6a', '#e8761a'
COMB = ['stx2a', 'stx2b', 'stx2c', 'stx2d']
SOLI = ['stx2e', 'stx2f', 'stx2i', 'stx2j', 'stx2k', 'stx2l', 'stx2n', 'stx2o']

def ital(s):   # stx2b -> italic stx + roman 2b
    return r'$\it{stx}$' + s[3:]

def save(fig, name):
    fig.savefig(OUT + name + '.pdf')
    fig.savefig(OUT + name + '.tiff', dpi=300, pil_kwargs={'compression': 'tiff_lzw'})
    plt.close(fig)
    print('  wrote', name)

# ---------------- Figure 1 ----------------
rowsS1 = [r for r in d['tableS1'] if r[1] >= 20]          # drop stx2m (8), stx2h (7)
rowsS1.sort(key=lambda r: -r[3])
fig, ax = plt.subplots(figsize=(6.6, 4.3))
ax.axvspan(0, 5, color='0.93', zorder=0)
y = np.arange(len(rowsS1))[::-1]
for yi, (s, n, k, p, lo, hi) in zip(y, rowsS1):
    if s in COMB:   c, m, fc = BLUE, 'o', BLUE
    elif s in SOLI: c, m, fc = ORANGE, '^', 'none'
    else:           c, m, fc = GREEN, 's', GREEN
    ax.plot([lo, hi], [yi, yi], color=c, lw=1.1, zorder=2)
    ax.plot([p], [yi], marker=m, color=c, mfc=fc, ms=6.5, zorder=3)
    ax.text(1.013, yi, format(n, ','), transform=ax.get_yaxis_transform(),
            ha='right', va='center', fontsize=8)
ax.set_yticks(y); ax.set_yticklabels([ital(r[0]) for r in rowsS1])
ax.set_xlabel('Genomes also carrying an ' + r'$\it{stx}$' + '1 subtype (%)')
ax.text(1.013, len(rowsS1) - 0.35, '$n$', transform=ax.get_yaxis_transform(),
        ha='right', va='center', fontsize=8, style='italic')
ax.set_xlim(-2, 72); ax.set_ylim(-0.7, len(rowsS1) - 0.3)
ax.grid(axis='x', color='0.9', lw=0.6); ax.set_axisbelow(True)
for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
ax.legend(handles=[Line2D([], [], marker='o', color=BLUE, ls='', ms=6.5, label='Combining class'),
                   Line2D([], [], marker='s', color=GREEN, ls='', ms=6.5, label='Intermediate'),
                   Line2D([], [], marker='^', color=ORANGE, mfc='none', ls='', ms=6.5, label='Solitary class')],
          loc='upper center', bbox_to_anchor=(0.5, -0.14), ncol=3, frameon=False, fontsize=8.5)
save(fig, 'Figure1_partner_acceptance')

# ---------------- Figure 2 ----------------
MK = ['eae', 'ehxA', 'espK', 'nleB', 'toxB', 'katP', 'espP', 'iha', 'subA']
order = ['stx2c', 'stx2a', 'stx2b', 'stx2d', 'stx2g', 'stx2f', 'stx2e', 'stx2j', 'stx2k', 'stx2n']
S2 = {r['sub']: r for r in d['tableS2']}
M = np.array([[S2[s][m] for m in MK] for s in order])
med = [float(S2[s]['medVF']) for s in order]
fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.9, 4.5), sharey=True,
                             gridspec_kw={'width_ratios': [3.3, 1.0], 'wspace': 0.06})
im = ax.imshow(M, cmap='Blues', vmin=0, vmax=100, aspect='auto')
ax.set_xticks(range(len(MK)))
ax.set_xticklabels([r'$\it{%s}$' % m for m in MK], rotation=40, ha='right')
ax.set_yticks(range(len(order))); ax.set_yticklabels([ital(s) for s in order])
for i in range(len(order)):
    for j in range(len(MK)):
        ax.text(j, i, '%.0f' % M[i, j], ha='center', va='center', fontsize=7.5,
                color='white' if M[i, j] > 55 else '0.15')
ax.axhline(4.5, color='black', lw=1.2)
ax.set_xticks(np.arange(-.5, len(MK), 1), minor=True)
ax.set_yticks(np.arange(-.5, len(order), 1), minor=True)
ax.grid(which='minor', color='white', lw=0.8); ax.tick_params(which='minor', length=0)
bx.barh(range(len(order)), med, color=[BLUE] * 5 + [ORANGE] * 5, height=0.62)
for i, v in enumerate(med):
    bx.text(v + 0.7, i, '%.0f' % v, va='center', fontsize=8)
bx.invert_yaxis(); ax.invert_yaxis()
bx.set_xlim(0, max(med) * 1.35); bx.set_xlabel('Median\nvirulence genes', fontsize=8.5)
for sp in ('top', 'right', 'left'): bx.spines[sp].set_visible(False)
bx.tick_params(axis='y', length=0)
cax = fig.add_axes([0.13, -0.04, 0.45, 0.028])
cb = fig.colorbar(im, cax=cax, orientation='horizontal')
cb.set_label('Genomes carrying marker (%)', fontsize=8.5); cb.ax.tick_params(labelsize=8)
save(fig, 'Figure2_accessory_repertoire')

# ---------------- Figure 3 ----------------
BINS = ['1-5', '6-50', '51-100', '101-200', '201-500', '>500']
fig, ax = plt.subplots(figsize=(6.4, 4.0))
for lin, col, mk, ls, fc in [('O157', BLUE, 'o', '-', BLUE),
                             ('non-O157', ORANGE, 's', '--', 'none')]:
    pts = [(BINS.index(r['bin']), r['p'], r['lo'], r['hi'])
           for r in d['tableS3'] if r['lin'] == lin and r['n'] >= 15]
    pts.sort()
    x = [p[0] for p in pts]; yv = [p[1] for p in pts]
    err = [[p[1] - p[2] for p in pts], [p[3] - p[1] for p in pts]]
    ax.errorbar(x, yv, yerr=err, color=col, marker=mk, mfc=fc, ls=ls, lw=1.4,
                ms=6, capsize=3, label='%s lineages' % lin)
ax.set_xticks(range(len(BINS))); ax.set_xticklabels(BINS)
ax.set_xlabel('Assembly contig count')
ax.set_ylabel(r'$\it{stx}$' + '2a genomes also carrying ' + r'$\it{stx}$' + '2c (%)')
ax.annotate('Contiguous assemblies', xy=(0, 50), xytext=(0.6, 57), fontsize=8.5,
            arrowprops=dict(arrowstyle='-', lw=0.7, color='0.5'))
ax.annotate('Fragmented assemblies', xy=(3, 1.5), xytext=(2.0, 17), fontsize=8.5,
            arrowprops=dict(arrowstyle='-', lw=0.7, color='0.5'))
ax.set_ylim(-3, 62)
ax.grid(axis='y', color='0.9', lw=0.6); ax.set_axisbelow(True)
for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
ax.legend(frameon=False, fontsize=8.5, loc='upper right')
save(fig, 'Figure3_stx2a_stx2c_contiguity')
