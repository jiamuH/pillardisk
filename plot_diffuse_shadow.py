"""Does a wave-independent diffuse glow filling the shadows feed the z-band host?

Changes from the first version, per discussion:
- Pillars moved in to r_min = 5 light days (hotter shadows).
- Diffuse glow is FLAT in F_nu (wavelength-independent, scattered-light-like),
  carrying a fraction f_diff of the local direct-irradiation power spread flat
  over the observed band range, with a Balmer bound-free rise (x1.5) blueward of
  3646 A. It is constant (lamp-independent) and fills the shadowed fraction, so
  it survives at lamp-off and is pillar-specific (the bowl has no shadows).
- Viscous floor grown: lamp set to eta_LP = 50 with the original viscous, then
  the viscous temperature raised x1.5 (lamp/viscous power ratio drops to ~10).
- Host overestimate measured in z band (8700 A), scored against the stellar host.
"""
import numpy as np
import matplotlib.pyplot as plt
from pillardisk.pillar_disk import C, H, K, ANGSTROM
from test_lag_spectrum import build_disk, add_pillars

H_P, SIGMA_R, SIGMA_PHI = 0.5, 0.25, 0.3
HLAMP = 50 * 0.00399
ETA_LP = 50.0
TVISC_BOOST = 1.5
R_MIN = 5.0
HOST_BETA, HOST_5100 = 3.0, 1.825
bands = np.array([1928., 2246., 2600., 3465., 4392., 5468., 6215., 7545., 8700.])
iV = int(np.argmin(np.abs(bands - 5100.0)))        # for the flux normalization
iZ = int(np.argmin(np.abs(bands - 8700.0)))        # z band, where we report
HOST = HOST_5100 * (bands / 5100.0) ** HOST_BETA
s_vals = np.linspace(1.08, 1.35, 12)               # 20% campaign window
fdiff_vals = np.array([0.0, 0.05, 0.1, 0.2, 0.3, 0.5])

# flat-F_nu diffuse normalization: bolometric blackbody intensity = SIGMA_I * T^4,
# spread flat over the observed frequency range DNU.
SIGMA_I = 2.0 * np.pi**4 * K**4 / (15.0 * H**3 * C**2)
nu_bands = C / (bands * ANGSTROM)
DNU = nu_bands.max() - nu_bands.min()
balmer = np.where(bands < 3646.0, 1.5, 1.0)        # Balmer bound-free rise blueward of 3646 A


def b_nu(lam_A, T):
    nu = C / (lam_A * ANGSTROM)
    x = np.clip(H * nu / (K * np.maximum(T, 1.0)), 1e-10, 700.0)
    return 2.0 * H * nu**3 / C**2 / (np.exp(x) - 1.0)


def eta_lp(dk):
    r2, p2 = np.meshgrid(dk.r, dk.phi, indexing='ij')
    h2 = dk.get_height(r2, p2)
    ds = np.sqrt(np.diff(dk.r)[:, None]**2 + (h2[1:, :] - h2[:-1, :])**2)
    da = ds * dk.r[:-1, None] * dk.dphi
    tv2 = np.interp(r2.ravel(), dk.r, dk.tv_base).reshape(r2.shape)
    T = dk.get_temperature(r2, p2, compute_shadows=True)
    tx4 = np.maximum(T**4 - tv2**4, 0.0)
    return np.sum(tx4[:-1, :] * da) / np.sum(tv2[:-1, :]**4 * da)


def setup(dk):
    dk.tx_base *= (ETA_LP / eta_lp(dk))**0.25      # lamp to eta=50 with original viscous
    dk.tv_base *= TVISC_BOOST                       # grow the viscous floor (ratio drops to ~10)


def precompute(dk):
    r2, p2 = np.meshgrid(dk.r, dk.phi, indexing='ij')
    h2 = dk.get_height(r2, p2)
    Tvisc = np.interp(r2.ravel(), dk.r, dk.tv_base).reshape(r2.shape)
    Tns = dk.get_temperature(r2, p2, compute_shadows=False)
    Tdir_nom = np.maximum(Tns**4 - Tvisc**4, 0.0)**0.25
    if len(dk.pillars) > 0:
        smask = dk._compute_shadow_mask(r2, p2, h2)
    else:
        smask = np.ones_like(r2)
    dh = h2[1:, :] - h2[:-1, :]
    ds = np.sqrt(np.diff(dk.r)[:, None]**2 + dh**2)
    px = -(dh / (ds + 1e-10)) * np.cos(dk.phi)[None, :]
    pz = np.diff(dk.r)[:, None] / (ds + 1e-10)
    dot = np.where(dk.ex * px + dk.ez * pz > 0, dk.ex * px + dk.ez * pz, 0.0)
    wgeo = ds * dk.r[:-1, None] * dk.dphi * dot
    return Tvisc, Tdir_nom, smask, wgeo, dk.fcol


def band_flux(prec, s, fdiff):
    Tvisc, Tdir_nom, smask, wgeo, fcol = prec
    Tdir = Tdir_nom * s
    T = (Tvisc**4 + (Tdir * smask)**4)**0.25                 # disc thermal (direct blocked by shadow)
    diff_amp = fdiff * SIGMA_I * Tdir_nom**4 / DNU           # flat-F_nu diffuse level per cell (constant)
    out = np.empty(len(bands))
    for k, lam in enumerate(bands):
        therm = b_nu(lam, T[:-1, :] * fcol) / fcol**4
        diff = diff_amp[:-1, :] * balmer[k] * (1.0 - smask[:-1, :])   # only in shadowed fraction
        out[k] = np.sum((therm + diff) * wgeo)
    return out


def driving_lc(F):
    B = F.mean(axis=0); A = F.std(axis=0)
    X = ((F - B) / A).mean(axis=1)
    return (X - X.mean()) / X.std()


def recovered(F_h):
    X = driving_lc(F_h)
    fit = np.array([np.polyfit(X, F_h[:, k], 1) for k in range(len(bands))])
    cross = np.max(-fit[:, 1] / fit[:, 0])
    return fit[:, 0] * cross + fit[:, 1]


dp = build_disk(hlamp=HLAMP); add_pillars(dp, 100, H_P, SIGMA_R, SIGMA_PHI, r_min=R_MIN)
setup(dp)
db = build_disk(hlamp=HLAMP); setup(db)
print(f"after raising viscous: lamp/viscous power ratio (pillar) = {eta_lp(dp):.1f}")
prec_p, prec_b = precompute(dp), precompute(db)
shadowed = prec_p[2] < 0.5 * (1.0 + dp.f_trans)
belt = dp.r >= R_MIN
print(f"mean shadow covering fraction over belt (r >= {R_MIN} ld) = "
      f"{shadowed.mean(axis=1)[belt].mean():.3f}")

ideal_p = np.zeros(len(fdiff_vals)); ideal_b = np.zeros(len(fdiff_vals))
act_p = np.zeros(len(fdiff_vals)); act_b = np.zeros(len(fdiff_vals))
print("All ratios are at the z band (8700 A).")
print("'true constant ratio' = (stellar host + constant disc flux) / stellar host")
print("'recovered ratio'      = recovered host / stellar host")
for i, fd in enumerate(fdiff_vals):
    Fp = np.array([band_flux(prec_p, s, fd) for s in s_vals])
    Fb = np.array([band_flux(prec_b, s, fd) for s in s_vals])
    cp = band_flux(prec_p, 0.0, fd)            # lamp-off constant (viscous + diffuse-in-shadow)
    cb = band_flux(prec_b, 0.0, fd)
    FS = 9.0 / Fp[:, iV].max()
    Fp *= FS; Fb *= FS; cp *= FS; cb *= FS
    ideal_p[i] = (HOST[iZ] + cp[iZ]) / HOST[iZ]
    ideal_b[i] = (HOST[iZ] + cb[iZ]) / HOST[iZ]
    act_p[i] = recovered(Fp + HOST)[iZ] / HOST[iZ]
    act_b[i] = recovered(Fb + HOST)[iZ] / HOST[iZ]
    print(f"  diffuse fraction = {fd:.2f}:  "
          f"true constant ratio  pillar = {ideal_p[i]:.2f}  bowl = {ideal_b[i]:.2f}   |   "
          f"recovered ratio  pillar = {act_p[i]:.2f}  bowl = {act_b[i]:.2f}")

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

fig, ax = plt.subplots(figsize=(9, 7))
ax.plot(fdiff_vals, ideal_p, '-o', color='purple', lw=3, ms=7,
        label=r'$\rm pillar~disc:~(stellar~host + constant~disc~flux)/stellar~host$')
ax.plot(fdiff_vals, ideal_b, '--s', color='purple', lw=2.2, ms=6, alpha=0.6,
        label=r'$\rm bowl~disc:~(stellar~host + constant~disc~flux)/stellar~host$')
ax.plot(fdiff_vals, act_p, '-o', color='orangered', lw=3, ms=7,
        label=r'$\rm pillar~disc:~recovered~host/stellar~host$')
ax.plot(fdiff_vals, act_b, '--s', color='orangered', lw=2.2, ms=6, alpha=0.6,
        label=r'$\rm bowl~disc:~recovered~host/stellar~host$')
ax.axhline(1.0, color='gray', ls=':', lw=2, alpha=0.8)
ax.set_xlabel(r'$\rm diffuse~floor~/~local~direct~irradiation~power$', fontsize=18)
ax.set_ylabel(r'$\rm recovered~host~/~stellar~host~~(z~band)$', fontsize=18)
ax.tick_params(which='major', direction='in', length=8, width=1.5,
               top=True, right=True, labelsize=14)
ax.tick_params(which='minor', direction='in', length=4, width=1.0,
               top=True, right=True)
ax.minorticks_on()
ax.legend(fontsize=10, frameon=False, loc='best')
plt.tight_layout()
plt.savefig('plots/test_diffuse_shadow.png', dpi=200, bbox_inches='tight')
plt.close()
print("\nSaved plots/test_diffuse_shadow.png")
