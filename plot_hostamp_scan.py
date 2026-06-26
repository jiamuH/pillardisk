"""How does the assumed galaxy-light amplitude affect the host overestimate?

Cai+24 (their Fig 3): the flux-variation-gradient overestimate grows as the
host fraction shrinks. Test it here in the current pillar regime (ETA_LP=50,
50 r_g, thin pillars). The disc fluxes are fixed; only the added host changes,
and the host cancels in the driver, so we just rescale HOST and redo the
recovery. x-axis = host fraction at V; y = recovered/true host. Pillar solid,
bowl dashed.
"""
import numpy as np
import matplotlib.pyplot as plt
import test_fluxflux_mock as m

bands = m.bands
camp = m.camp
iV = m.iV
HOST_BETA = m.HOST_BETA
Fp, Fb = m.F_states, m.F_bowl              # AGN band fluxes (scaled), pillar/bowl
floor_p, floor_b = m.F_const_b, m.F_const_bowl   # viscous floors (scaled)
AGN_V = Fp[:, iV].max()                     # bright-state AGN V flux (=9 by scaling)


def driving_lc(F):
    B = F[camp].mean(axis=0); A = F[camp].std(axis=0)
    X = ((F - B) / A).mean(axis=1)
    return (X - X[camp].mean()) / X[camp].std()


def recovered(F_h):
    X = driving_lc(F_h)
    fit = np.array([np.polyfit(X[camp], F_h[camp, k], 1) for k in range(len(bands))])
    anchor = np.max(-fit[:, 1] / fit[:, 0])
    return fit[:, 0] * anchor + fit[:, 1]


h5_vals = np.array([0.5, 1.0, 2.0, 3.5, 6.0, 10.0, 15.0, 25.0, 40.0])
frac_V = np.zeros(len(h5_vals))
rp = np.zeros((len(h5_vals), len(bands)))
rb = np.zeros((len(h5_vals), len(bands)))
for i, h5 in enumerate(h5_vals):
    HOST = h5 * (bands / 5100.0) ** HOST_BETA
    frac_V[i] = HOST[iV] / (HOST[iV] + AGN_V)
    rp[i] = recovered(Fp + HOST) / (HOST + floor_p)
    rb[i] = recovered(Fb + HOST) / (HOST + floor_b)

show = [int(np.argmin(np.abs(bands - w))) for w in (5468, 7545, 8700)]
print(f"{'host5100':>8} {'fracV':>6}  pillar(5468/7545/8700) | bowl")
for i, h5 in enumerate(h5_vals):
    pp = "/".join(f"{rp[i, k]:.2f}" for k in show)
    bb = "/".join(f"{rb[i, k]:.2f}" for k in show)
    print(f"{h5:8.1f} {frac_V[i]:6.2f}  {pp} | {bb}")

# --- plot ---
plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

clr = {5468: 'seagreen', 7545: 'darkorange', 8700: 'orangered'}
fig, ax = plt.subplots(figsize=(9, 7))
for k, w in zip(show, (5468, 7545, 8700)):
    ax.plot(frac_V, rp[:, k], '-o', color=clr[w], lw=3, ms=7,
            label=r'$\rm pillar,~' + f'{w}' + r'~\AA$')
    ax.plot(frac_V, rb[:, k], '--s', color=clr[w], lw=2.2, ms=6, alpha=0.7,
            label=r'$\rm bowl,~' + f'{w}' + r'~\AA$')
ax.axhline(1.0, color='gray', ls=':', lw=2, alpha=0.8)
ax.set_xlabel(r'$\rm host~fraction~at~V$', fontsize=18)
ax.set_ylabel(r'$f_\nu^{\rm recovered}/f_\nu^{\rm true}~\rm (host)$', fontsize=18)
ax.tick_params(which='major', direction='in', length=8, width=1.5,
               top=True, right=True, labelsize=14)
ax.tick_params(which='minor', direction='in', length=4, width=1.0,
               top=True, right=True)
ax.minorticks_on()
ax.legend(fontsize=12, frameon=False, ncol=2, loc='best')
plt.tight_layout()
plt.savefig('plots/test_hostamp_scan.png', dpi=200, bbox_inches='tight')
plt.close()
print("\nSaved plots/test_hostamp_scan.png")
