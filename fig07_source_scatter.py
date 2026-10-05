#!/usr/bin/env python3
"""
Figure 7 v2 — Source-plane scatter.
- Alternating label offsets in panel (c) to prevent overlap at 49 systems.
- All three panels at the same axes width, aligned flush at left and right.
"""
import json
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

with open('/home/claude/bullet/final_results.json') as fp:
    r = json.load(fp)
bp = pd.read_csv('/home/claude/bullet/final_per_image_backprojection.csv')

def per_sys(cols):
    vals, sysnames, Ns, zs_sys = [], [], [], []
    for s, g in bp.groupby('sys'):
        if len(g) < 2: continue
        bx = g[cols[0]].values; by = g[cols[1]].values
        d2 = (bx-bx.mean())**2 + (by-by.mean())**2
        vals.append(np.sqrt(d2.mean()))
        sysnames.append(s); Ns.append(len(g)); zs_sys.append(g.zs.iloc[0])
    return sysnames, np.array(vals), np.array(Ns), np.array(zs_sys)

sysnames, lc_vals, Ns, zs_sys = per_sys(('bxL','byL'))
_, ct_vals, _, _               = per_sys(('bxC','byC'))

FS=8; FST=8
plt.rcParams.update({
    'font.size':FS,'axes.labelsize':FS,'axes.titlesize':FST,
    'xtick.labelsize':7.5,'ytick.labelsize':7.5,'legend.fontsize':7,
    'figure.dpi':300,'savefig.dpi':300,'lines.linewidth':1.0,
})

# ─ Figure layout: three stacked panels, same axes width ───────────────────
# Use GridSpec with uniform width_ratios; attach colorbars via divider so
# every main axes has the same horizontal extent regardless of whether a
# colorbar is drawn.
# Full-page-width figure (declared \begin{figure*} in LaTeX).  MDPI
# text-column width is 17 cm ~ 6.7 in; use that so all 49 system labels
# can be read at a comfortable 7pt.
fig = plt.figure(figsize=(6.7, 7.5))
gs  = fig.add_gridspec(3, 1, hspace=0.55, left=0.09, right=0.90,
                       top=0.985, bottom=0.07)
axA = fig.add_subplot(gs[0])
axB = fig.add_subplot(gs[1])
axC = fig.add_subplot(gs[2])

def _reserve_colorbar_space(ax):
    """Steal a thin invisible strip on the right so axis width matches
    panels that have a visible colorbar."""
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.08)
    cax.set_visible(False)
    return cax

# ─ (a) aggregate bars with mass-band error bars ───────────────────────────
central, lo, hi = r['central'], r['M-30'], r['M+30']
metric_keys = ['RMS','Mean','Med']
metric_labels = ['RMS','Mean/system','Median/system']
x = np.arange(3); w = 0.38
LC_K='LC'; CT_K='CTL'

lc_c  = [central[LC_K][k]  for k in metric_keys]
lc_lo = [min(lo[LC_K][k], hi[LC_K][k]) for k in metric_keys]
lc_hi = [max(lo[LC_K][k], hi[LC_K][k]) for k in metric_keys]
lc_err = np.array([[c-l for c,l in zip(lc_c,lc_lo)],
                   [h-c for c,h in zip(lc_c,lc_hi)]])

ct_c  = [central[CT_K][k] for k in metric_keys]
ct_lo = [min(lo[CT_K][k], hi[CT_K][k]) for k in metric_keys]
ct_hi = [max(lo[CT_K][k], hi[CT_K][k]) for k in metric_keys]
ct_err = np.array([[c-l for c,l in zip(ct_c,ct_lo)],
                   [h-c for c,h in zip(ct_c,ct_hi)]])

axA.bar(x - w/2, lc_c, w, color='steelblue', alpha=0.85, edgecolor='white',
        yerr=lc_err, capsize=4, ecolor='navy',
        label=r'Unoptimised $\Lambda$CDM (k=2)')
axA.bar(x + w/2, ct_c, w, color='crimson', alpha=0.85, edgecolor='white',
        yerr=ct_err, capsize=4, ecolor='darkred', label='CCC+TL (k=2)')
axA.axhline(8.2801, ls=':', color='royalblue', lw=1.3,
            label=r'Lenstool $\Lambda$CDM optimised (k=32)')

axA.set_xticks(x); axA.set_xticklabels(metric_labels, fontsize=FS)
axA.set_ylabel('Source-plane scatter [arcsec]', fontsize=FS)
axA.legend(fontsize=6.5, loc='upper left')
axA.grid(True, alpha=0.3, axis='y')
axA.annotate(f'CTL/$\\Lambda$CDM\nRMS = {central[CT_K]["RMS"]/central[LC_K]["RMS"]:.3f}\n'
             f'Mean = {central[CT_K]["Mean"]/central[LC_K]["Mean"]:.3f}\n'
             f'Error bars: $\\pm 30\\%$ $M_{{200}}$',
             xy=(0.62,0.03), xycoords='axes fraction', fontsize=6,
             ha='left', va='bottom',
             bbox=dict(boxstyle='round', fc='lightyellow', ec='orange', alpha=0.9))
_reserve_colorbar_space(axA)      # match panels (b),(c) axis width

# ─ (b) CTL vs LCDM scatter plot ──────────────────────────────────────────
sc = axB.scatter(lc_vals, ct_vals, c=zs_sys, cmap='plasma',
                 s=32, alpha=0.85, zorder=5, vmin=0.9, vmax=7.0,
                 edgecolors='k', linewidths=0.4)
divB = make_axes_locatable(axB)
caxB = divB.append_axes('right', size='5%', pad=0.08)
cbB = plt.colorbar(sc, cax=caxB); cbB.set_label('Source redshift', fontsize=FS)
cbB.ax.tick_params(labelsize=7)

vmax = max(lc_vals.max(), ct_vals.max()) * 1.1
axB.plot([0,vmax],[0,vmax],'k--',lw=1.5,label='CTL = $\\Lambda$CDM')
axB.fill_between([0,vmax],[0,0.7*vmax],[0,1.3*vmax],alpha=0.10,color='green',
                 label='$\\pm 30\\%$ $M_{200}$ band')
axB.set_xlabel(r'Unoptimised $\Lambda$CDM scatter [arcsec]', fontsize=FS)
axB.set_ylabel('CCC+TL scatter [arcsec]', fontsize=FS)
axB.legend(fontsize=6); axB.grid(True, alpha=0.3)
axB.set_xlim(0,vmax); axB.set_ylim(0,vmax)
# aspect='equal' not used: panel (b) spans the full page width to match (a)
# and (c); the diagonal line y=x is still shown.
axB.annotate(f'All-system RMS\n$\\Lambda$CDM = {central[LC_K]["RMS"]:.2f}\"\n'
             f'CCC+TL = {central[CT_K]["RMS"]:.2f}\"',
             xy=(0.03,0.95), xycoords='axes fraction', fontsize=6.5, va='top',
             bbox=dict(boxstyle='round',fc='w',alpha=0.85))

# ─ (c) per-system ratio bars, staggered labels ───────────────────────────
ratios = sorted(zip(sysnames, ct_vals/np.maximum(lc_vals,1e-6), zs_sys),
                key=lambda x: x[1])
snames = [t[0] for t in ratios]
rvals  = [t[1] for t in ratios]
zvals  = [t[2] for t in ratios]
norm = plt.Normalize(0.9, 7.0)
cols = [plt.cm.plasma(norm(z)) for z in zvals]

axC.bar(range(len(rvals)), rvals, color=cols, alpha=0.85, edgecolor='white')
axC.axhline(1.0, color='k', lw=1.5, ls='--', label='CTL = $\\Lambda$CDM')
rms_ratio = central[CT_K]['RMS']/central[LC_K]['RMS']
axC.axhline(rms_ratio, color='crimson', lw=1.3, ls=':',
            label=f'Overall RMS ratio = {rms_ratio:.3f}')
axC.fill_between([-0.5,len(rvals)-0.5],[0.7],[1.3],alpha=0.10,color='green',
                 label='$\\pm 30\\%$ $M_{200}$ band')

# Staggered x-tick labels: alternating two-row placement with vertical offsets.
# With the figure now at full-page width (~6.7 in) and 49 systems, each label
# gets ~0.13 in of horizontal space — comfortable at 6.5pt with a two-row
# stagger.  Leader lines connect the lower-row labels to their bar centres.
axC.set_xticks(range(len(snames)))
axC.set_xticklabels([''] * len(snames))     # hide default labels
y_lo = -0.03
y_hi = -0.10
trans = axC.get_xaxis_transform()
for i, name in enumerate(snames):
    y_text = y_lo if i % 2 == 0 else y_hi
    if i % 2 == 1:                           # thin leader for lower-row labels
        axC.plot([i, i], [0, y_hi + 0.015], color='0.5', lw=0.4,
                 transform=trans, clip_on=False)
    axC.text(i, y_text, name, rotation=0, ha='center', va='top',
             fontsize=6.5, transform=trans, clip_on=False)

axC.set_ylabel('CCC+TL / $\\Lambda$CDM scatter ratio', fontsize=FS)
axC.legend(fontsize=6, loc='upper left'); axC.grid(True, alpha=0.3, axis='y')
axC.set_xlim(-0.6, len(rvals)-0.4)
axC.margins(x=0.005)

divC = make_axes_locatable(axC)
caxC = divC.append_axes('right', size='5%', pad=0.08)
sm = plt.cm.ScalarMappable(cmap='plasma', norm=norm); sm.set_array([])
cbC = plt.colorbar(sm, cax=caxC); cbC.set_label('Source redshift', fontsize=FS)
cbC.ax.tick_params(labelsize=7)

fig.savefig('/home/claude/bullet/fig07_v2.pdf', bbox_inches='tight', dpi=300)
plt.close(fig)
print(f"fig07_v2.pdf written. N_sys={len(snames)}, label-stagger alternating.")
