"""How the OLD flux-flux host recovery changes with the probed variation %.

The old method anchors the recovered host at the driver value where the bluest
band's linear fit reaches zero total flux (nonvar_anchor). CAMP_FRAC sets how
much of the lamp's full amplitude the "observed campaign" spans (currently
0.20). Changing it changes which states are observed, the driver normalization,
the per-band linear fits, and the anchor -- so the recovered host moves.

Here we recompute the band SEDs once on a DENSE uniform lamp grid (so every
campaign window is sampled evenly), then for a range of CAMP_FRAC values redo
the driver / fit / anchor and read off the recovered host. We plot the ratio
recovered/true host versus the probed percentage, for the pillar and bowl discs.

Imports the built, rescaled discs and helpers from test_fluxflux_mock.
"""
import numpy as np
import matplotlib.pyplot as plt
import test_fluxflux_mock as m

S_FULL = m.S_FULL

# --- dense uniform lamp grid and band SEDs (computed once per disc) ---
s_dense = np.linspace(0.0, S_FULL, 40)
print("computing dense-grid SEDs (pillar)...")
Fp, _ = m.fluxflux(m.disk, m.bands, s_dense)
print("computing dense-grid SEDs (bowl)...")
Fb, _ = m.fluxflux(m.bowl, m.bands, s_dense)
Fp *= m.FLUX_SCALE
Fb *= m.FLUX_SCALE
Fp_h = Fp + m.HOST                       # + host
Fb_h = Fb + m.HOST

true_p = m.HOST + m.F_const_b            # true host = host + viscous floor
true_b = m.HOST + m.F_const_bowl


def driving_lc(F, camp):
    B = F[camp].mean(axis=0)
    A = F[camp].std(axis=0)
    X = ((F - B) / A).mean(axis=1)
    return (X - X[camp].mean()) / X[camp].std()


def recovered_host(F_h, camp):
    X = driving_lc(F_h, camp)
    fit = np.array([np.polyfit(X[camp], F_h[camp, k], 1)
                    for k in range(len(m.bands))])
    anchor = np.max(-fit[:, 1] / fit[:, 0])      # nonvar_anchor
    return fit[:, 0] * anchor + fit[:, 1]


fracs = np.linspace(0.1, 0.95, 25)
rec_p = np.zeros((len(fracs), len(m.bands)))
rec_b = np.zeros((len(fracs), len(m.bands)))
for i, fr in enumerate(fracs):
    S_LO = S_FULL * (1.0 - fr)
    camp = s_dense >= S_LO
    rec_p[i] = recovered_host(Fp_h, camp)
    rec_b[i] = recovered_host(Fb_h, camp)

ratio_p = rec_p / true_p                  # recovered / true, per band
ratio_b = rec_b / true_b

# print a table at a few fractions for the bands shown
show = [int(np.argmin(np.abs(m.bands - w))) for w in (3465, 5468, 8700)]
print("\nrecovered/true host (pillar | bowl) at selected bands:")
print(f"{'frac%':>6} " + " ".join(f"{m.bands[k]:>6.0f}A" for k in show) + "  |  "
      + " ".join(f"{m.bands[k]:>6.0f}A" for k in show))
for i, fr in enumerate(fracs):
    if i % 4 == 0 or i == len(fracs) - 1:
        pp = " ".join(f"{ratio_p[i, k]:7.3f}" for k in show)
        bb = " ".join(f"{ratio_b[i, k]:7.3f}" for k in show)
        print(f"{fr*100:6.0f} {pp}  |  {bb}")

# --- plot ---
plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

band_lab = {3465: r'$3465~\AA$', 5468: r'$5468~\AA$', 8700: r'$8700~\AA$'}
band_clr = {3465: 'royalblue', 5468: 'seagreen', 8700: 'orangered'}

fig, ax = plt.subplots(figsize=(9, 7))
for k, w in zip(show, (3465, 5468, 8700)):
    ax.plot(fracs * 100, ratio_p[:, k], '-', color=band_clr[w], lw=3, alpha=0.9,
            label=r'$\rm pillar,~' + band_lab[w].strip('$') + r'$')
    ax.plot(fracs * 100, ratio_b[:, k], '--', color=band_clr[w], lw=2.5, alpha=0.8,
            label=r'$\rm bowl,~' + band_lab[w].strip('$') + r'$')
ax.axhline(1.0, color='gray', ls=':', lw=2, alpha=0.8)        # perfect recovery
ax.axvline(20.0, color='k', ls=':', lw=1.5, alpha=0.6)        # current value
ax.text(20.5, ax.get_ylim()[0], r'$\rm current$', rotation=90,
        va='bottom', ha='left', fontsize=12)
ax.set_xlabel(r'$\rm probed~variation~amplitude~[\%]$', fontsize=18)
ax.set_ylabel(r'$f_\nu^{\rm recovered}/f_\nu^{\rm true}~\rm (host)$', fontsize=18)
ax.tick_params(which='major', direction='in', length=8, width=1.5,
               top=True, right=True, labelsize=14)
ax.tick_params(which='minor', direction='in', length=4, width=1.0,
               top=True, right=True)
ax.minorticks_on()
ax.legend(fontsize=12, frameon=False, ncol=2, loc='best')
plt.tight_layout()
plt.savefig('plots/test_campfrac_scan.png', dpi=200, bbox_inches='tight')
plt.close()
print("\nSaved plots/test_campfrac_scan.png")
