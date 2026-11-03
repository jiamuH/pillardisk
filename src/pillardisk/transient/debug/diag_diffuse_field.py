"""diag_diffuse_field.py - how strong is the diffuse ionising field that the
1D per-ray Cloudy models leave out? (docs/NOTE_diffuse_field.md)

Each Cloudy ray's reflected ionising photons (.conA reflc + reflin, all
wavelengths < 911.76 A) are emitted isotropically from its skin point
(spread over 6 points across its phi wedge) and received (a) at every other
ray's skin point and (b) on the disc's upper surface (the last cell above
the first n_H > N_SURF cell from the top, every 8th phi and 4th r column),
wherever the straight path does not cross a cell with n_H > N_SURF. Printed
and saved: diffuse / direct lamp photon flux ((b): the lamp flux the point
would get if unblocked; lit = the lamp sees it).

Run: python3 -m pillardisk.transient.debug.diag_diffuse_field <run> <out.npz>
(~30 s for a 300-ray run; scales as the number of rays); plot with
diag_diffuse_field_plot <dir holding diffuse_<h>_surface.npz>."""
import numpy as np, sys, time, pandas as pd
from pillardisk.transient.cloudy_rays import (load_sim, to_physical, DUMP,
    N_SURF, ray_geometry, density_interp)
from pillardisk.transient.ray_emulator import read_deck, ray_features
from pillardisk.transient.sim_spectrum import CONFIG, LD_CM
H, C = 6.626e-27, 2.998e10
t0 = time.time()
src = sys.argv[1]; nsub = 6
r0 = CONFIG['R0_LD'] * LD_CM
sim = load_sim(DUMP); phys = to_physical(sim, CONFIG)
dense3 = phys['nH'] > N_SURF            # (nph, nth, nr), nearest cell
def isdense(P):
    r = np.linalg.norm(P, axis=-1)
    t = np.arccos(np.clip(P[..., 2] / r, -1, 1))
    p = np.arctan2(P[..., 1], P[..., 0]) % (2 * np.pi)
    ir = np.clip(np.searchsorted(sim['rf'], r) - 1, 0, len(sim['r']) - 1)
    it = np.clip(np.searchsorted(sim['thf'], t) - 1, 0, len(sim['th']) - 1)
    ip = np.clip(np.searchsorted(sim['phf'], p) - 1, 0, len(sim['ph']) - 1)
    inside = (r >= sim['rf'][0]) & (r <= sim['rf'][-1])
    inside = inside & (t >= sim['thf'][0]) & (t <= sim['thf'][-1])
    return inside & dense3[ip, it, ir]
geo, _, h = ray_geometry(src, sim)
m = np.load(f'{src}/rays_meta.npz')
rows = []
for nm, plo, phi_hi in zip(m['name'].astype(str), m['phi_lo'], m['phi_hi']):
    try:
        a = pd.read_csv(f'{src}/rays/{nm}.conA', sep='	', usecols=[0, 1, 5, 7],
                        skiprows=1, header=None).values[::-1]
    except (OSError, pd.errors.EmptyDataError):
        continue
    if a.ndim < 2 or len(a) < 10:
        continue
    ion = a[:, 0] < 911.76
    w = (a[:, 0] * 1e-8 / (H * C) * np.gradient(np.log(a[:, 0])))[ion]
    qi = np.sum(a[ion, 1] * w); qr = np.sum((a[ion, 2] + a[ion, 3]) * w)
    lphi, s, n = read_deck(f'{src}/rays/{nm}.in')
    f, sk = ray_features(lphi, s, n, r0)
    s0, sa, ca, dom = geo(nm)
    rows.append((nm, plo, phi_hi, sa, ca, sk / r0, s0, dom, qr / qi,
                 10 ** lphi, 10 ** f['lphi_skin']))
nm = [r[0] for r in rows]
plo, phh, sa, ca, sk, s0, dom, alb, phi0, phisk = (np.array([r[k] for r in rows])
                                                    for k in range(1, 11))
ph = 0.5 * (plo + phh)
pos = lambda p: np.stack([sk * sa * np.cos(p), sk * sa * np.sin(p), h + sk * ca], -1)
X = pos(ph)                                         # receivers (n, 3) [r0]
Nph = alb * phi0 * s0 ** 2 * dom                    # reflected photons/s
# sources: each ray spread over nsub points across its phi wedge
fr = (np.arange(nsub) + 0.5) / nsub
PS = np.stack([sk[:, None] * sa[:, None] * np.cos(plo[:, None] + fr[None, :] * (phh - plo)[:, None]),
               sk[:, None] * sa[:, None] * np.sin(plo[:, None] + fr[None, :] * (phh - plo)[:, None]),
               np.broadcast_to((h + sk * ca)[:, None], (len(sk), nsub))], -1)
soft2 = (sk ** 2 * dom)                             # patch area [r0^2]
recv = np.zeros(len(nm)); recv_vis = np.zeros(len(nm))
lt = np.linspace(0.03, 0.97, 60)
for j in range(len(nm)):
    d = PS - X[j][None, None]                       # (n, nsub, 3)
    d2 = np.sum(d ** 2, -1) + soft2[:, None]
    ok = np.ones(d2.shape, bool); ok[j] = False
    # occlusion: points along each segment, dense gas blocks
    P = X[j][None, None, None] + lt[None, None, :, None] * d[:, :, None, :]
    dense = isdense(P)
    vis = ok & ~dense.any(-1)
    flux = (Nph[:, None] / nsub) / (4 * np.pi * d2 * r0 ** 2)
    recv[j] = np.sum(np.where(ok, flux, 0))
    recv_vis[j] = np.sum(np.where(vis, flux, 0))
R = recv_vis / phisk
R0 = recv / phisk
phd = np.degrees(ph) % 360
np.savez(sys.argv[2], name=nm, ph=phd, R_skin=sk, ratio=R, ratio_noocc=R0,
         albedo=alb, h=h)
print(f'{src}: {len(nm)} rays, h = {h} r0, {time.time()-t0:.0f} s')
print(f'reflected/incident ionising photons: median {np.median(alb):.2f}')
for lab, m_ in (('all', np.ones(len(nm), bool)), ('fan 240-300', (phd >= 240) & (phd < 300))):
    print(f'{lab:12s} diffuse/direct at skin (occluded): median {np.median(R[m_]):.3f}, '
          f'84% {np.percentile(R[m_], 84):.3f}, >0.3: {np.mean(R[m_] > 0.3):.0%} | no occlusion median {np.median(R0[m_]):.3f}')

# ---- receivers on the whole dense surface (lit and shadowed) ----
t1 = time.time()
nH = phys['nH']                                   # (nph, nth, nr)
ths, rs, phs = sim['th'], sim['r'], sim['ph']
step_p, step_r = 8, 4
surf = []
for ip in range(0, len(phs), step_p):
    for ir in range(0, len(rs), step_r):
        col = nH[ip, :, ir]
        k = np.nonzero(col > N_SURF)[0]
        if len(k) == 0:
            continue
        it = max(k[0] - 1, 0)                      # last thin cell above it
        surf.append((ip, ir, rs[ir], ths[it], phs[ip]))
surf = np.array(surf)
rr, tt, pp = surf[:, 2], surf[:, 3], surf[:, 4]
Y = np.stack([rr * np.sin(tt) * np.cos(pp), rr * np.sin(tt) * np.sin(pp),
              rr * np.cos(tt)], 1)
# step just above the surface (toward the lamp side): half a theta cell up
Y = Y * 1.0
lamp = np.array([0, 0, h])
Qsr = np.median(phi0 * s0 ** 2)                   # lamp photons/s/sr
lt2 = np.linspace(0.03, 0.97, 30)
dl = Y - lamp
dlamp2 = np.sum(dl ** 2, 1)
Pl = lamp[None, None] + lt2[None, :, None] * dl[:, None, :]
lit = ~isdense(Pl).any(-1)
direct = Qsr / (dlamp2 * r0 ** 2)
PSf = PS.reshape(-1, 3); Nf = np.repeat(Nph / nsub, nsub)
soft_f = np.repeat(soft2, nsub)
diff = np.zeros(len(Y))
for j in range(len(Y)):
    d = PSf - Y[j]
    P = Y[j][None, None] + lt2[None, :, None] * d[:, None, :]
    vis = ~isdense(P).any(-1)
    diff[j] = np.sum(np.where(vis, Nf / (4 * np.pi * (np.sum(d ** 2, 1) + soft_f) * r0 ** 2), 0))
ratio = diff / direct          # vs the lamp flux the point would get unblocked
np.savez(sys.argv[2].replace('.npz', '_surface.npz'), ip=surf[:, 0], ir=surf[:, 1],
         r=rr, th=tt, ph=pp, lit=lit, diffuse=diff, direct=direct, h=h)
print(f'surface: {len(Y)} points ({lit.mean():.0%} lit), {time.time()-t1:.0f} s')
for lab, m_ in (('lit', lit), ('shadowed', ~lit)):
    if m_.any():
        print(f'{lab:9s} diffuse / unblocked-lamp flux: median {np.median(ratio[m_]):.3f}, '
              f'16-84% {np.percentile(ratio[m_], 16):.3f}-{np.percentile(ratio[m_], 84):.3f}, >0.1: {np.mean(ratio[m_] > 0.1):.0%}')
