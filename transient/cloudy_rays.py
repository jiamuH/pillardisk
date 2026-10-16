#!/usr/bin/env python3
"""cloudy_rays.py - one Cloudy model per lamp ray, following the sim's own
density profile along the ray (instead of one uniform slab per ray read
off the arm_column2 grid).

The rays are the sim's radial columns: the lamp sits at the origin, so
every (phi, theta) column r = rin..rout is a straight ray from the lamp.
Each deck (settings as transient/cloudy/arm_column2.in, the grid the
pipeline uses):
  - the same incident SED (interpolate table), init "c84.ini", solar,
    turbulence 100 km/s, stop temperature 3e3 K, continuum resolution 0.2
  - radius = the ray's inner and outer sim boundary, phi(H) = Q_ION /
    4 pi r_in^2 at the inner face; Cloudy dilutes it as 1/r^2 along the
    ray (open geometry), so no foreshortening factor: a tilted surface is
    a gradual density rise along the radial ray, which Cloudy resolves
  - dlaw table radius = the sim's log n_H at the cell centres (log-log
    interpolation), padded past both ends
  - stop column density STOP_COL (as the grid's maximum)
Shadowing is built in: each ray carries all the gas between the lamp and
any depth. Still 1D: no diffuse light exchanged between rays.

Outputs per ray in <out>/rays/<name>.{in,out,conA,ems}; `extract`
gathers the observer-relevant spectra into <out>/rays.npz: per ray
refl/out line and continuum nu F_nu [erg/s/cm^2 at the inner face] on a
common wavelength grid, plus (iph, ith, r_in) - a ray's luminosity per
sr is F r_in^2.

Steps (all with --out DIR):
  write    decks for the upper-hemisphere rays (every --dph-th phi)
  run      run the decks in parallel (--ncpu), skipping finished rays
  extract  collect the spectra into rays.npz
  plot     finished rays of one phi column (--iph): their paths over the
           sim's n_H in the meridional plane, and per ray the depth
           profiles of the H alpha, Mg II, C IV emissivities (red, green,
           blue, as the RGB line maps) with n_H and T_e

Run:  python3 -m transient.cloudy_rays write --out /data2/jhuang/runs/cloudy_rays/edd10
      python3 -m transient.cloudy_rays run --out ... --ncpu 90
      python3 -m transient.cloudy_rays extract --out ...
"""

import argparse
import glob
import os
import subprocess
import time
from multiprocessing import Pool

import numpy as np

from transient.sim_spectrum import load_sim, to_physical, CONFIG, LD_CM
from transient.sim_cloudy_spectrum import Q_ION

HERE = os.path.dirname(os.path.abspath(__file__))
DUMP = os.path.join(HERE, 'data', 'sim', 'disk.out1.00012.athdf')
CLOUDY_DIR = '/data2/jhuang/repos/cloudy'
EXE = os.path.join(CLOUDY_DIR, 'source', 'cloudy.exe')
STOP_COL = 24.5           # log N_H [cm^-2], the arm_column2 grid maximum
WMIN, WMAX = 300.0, 12000.0

DECK = '''title pillardisk lamp ray {name} (phi {ph:.4f}, theta {th:.4f})
init "c84.ini"
# incident AGN SED: as transient/cloudy/arm_column2.in
interpolate (-8.0 -13.6) (-1.7 -1.0) (-1.3 1.48) (-1.0 2.0)
continue    (-0.19 1.65) (0.4 1.32) (1.87  -1.87)  (2.43 -2.37)
continue    (2.98 -2.79) (3.3 -2.96) (3.7 -3.38) (4.18 -4.05)
continue    (4.99 -5.59) (5.1 -6.67) (7.0 -12.37)
phi(H) {logphi:.4f}
radius {logr_in:.5f} {logr_out:.5f}
dlaw table radius
{dlaw}end of dlaw
turbulence = 100 km/s
stop column density {stop_col:.2f}
stop temperature 3e3 K
set continuum resolution 0.2
print last
save continuum "{name}.conA" units Angstroms last
save overview "{name}.ovr" last
save lines emissivity "{name}.ems" last
H  1  6562.80A
H  1  4861.32A
C  4  1548.19A
C  4  1550.77A
Mg 2  2795.53A
Mg 2  2802.71A
end of lines
'''


def ray_decks(sim, phys, dph):
    """(name, iph, ith, r_in [cm], deck text) for the upper-hemisphere
    rays, every dph-th phi."""
    r0 = CONFIG['R0_LD'] * LD_CM
    rf = sim['rf'] * r0
    logr = np.log10(sim['r'] * r0)
    up = np.nonzero(sim['th'] <= np.pi / 2)[0]
    out = []
    for iph in range(0, len(sim['ph']), dph):
        for ith in up:
            logn = np.log10(phys['nH'][iph, ith, :])
            # pad past both faces so Cloudy never extrapolates
            pts = ([(np.log10(rf[0]) - 0.01, logn[0])]
                   + list(zip(logr, logn))
                   + [(np.log10(rf[-1]) + 0.01, logn[-1])])
            dlaw = ''.join(f'continue {a:.6f} {b:.4f}\n' for a, b in pts)
            name = f'ray_p{iph:03d}_t{ith:02d}'
            deck = DECK.format(
                name=name, ph=sim['ph'][iph], th=sim['th'][ith],
                logphi=np.log10(Q_ION / (4 * np.pi * rf[0] ** 2)),
                logr_in=np.log10(rf[0]), logr_out=np.log10(rf[-1]),
                dlaw=dlaw, stop_col=STOP_COL)
            out.append((name, iph, ith, rf[0], deck))
    return out


def run_one(args):
    rdir, name = args
    if os.path.exists(os.path.join(rdir, name + '.conA')) and \
            'Cloudy exited OK' in open(os.path.join(rdir, name + '.out'),
                                       errors='ignore').read()[-3000:]:
        return name, 0.0, 'skip'
    t = time.time()
    with open(os.path.join(rdir, name + '.in')) as fi, \
            open(os.path.join(rdir, name + '.out'), 'w') as fo:
        rc = subprocess.run([EXE], stdin=fi, stdout=fo, cwd=rdir,
                            env=dict(os.environ,
                                     CLOUDY_DATA_PATH=os.path.join(
                                         CLOUDY_DIR, 'data'))).returncode
    return name, time.time() - t, 'ok' if rc == 0 else f'rc={rc}'


def read_conA(path):
    """[wave, DiffOut, reflc, reflin, outlin] rows in WMIN..WMAX, ascending
    wavelength (columns as transient.extract_cloudy_arm)."""
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.split('\t')
            w = float(f[0])
            if WMIN <= w <= WMAX:
                rows.append((w, float(f[3]), float(f[5]), float(f[7]),
                             float(f[8])))
    return np.array(rows)[::-1]


def plot_rays(out, iph):
    """Meridional n_H slice at phi index iph with the finished rays, and
    per ray the line emissivity, n_H and T_e against depth."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                         'font.family': 'serif', 'font.weight': 'heavy',
                         'font.size': 18})
    plt.rcParams['text.latex.preamble'] = \
        r'\usepackage{amsmath} \usepackage{bm} \boldmath'
    rdir = os.path.join(out, 'rays')
    names = sorted(os.path.basename(p)[:-4] for p in glob.glob(
        os.path.join(rdir, f'ray_p{iph:03d}_t*.ems')))
    names = [n for n in names if 'Cloudy exited OK' in open(
        os.path.join(rdir, n + '.out'), errors='ignore').read()[-3000:]]
    if not names:
        raise SystemExit(f'no finished rays at iph = {iph}')
    sim = load_sim(DUMP)
    phys = to_physical(sim, CONFIG)
    r0 = CONFIG['R0_LD'] * LD_CM
    fig = plt.figure(figsize=(4.6 * len(names) + 1, 15))
    gs = fig.add_gridspec(3, len(names), height_ratios=[1.25, 1, 1],
                          wspace=0.42, hspace=0.25)
    ax = fig.add_subplot(gs[0, :])
    R, T = np.meshgrid(sim['rf'], sim['thf'])
    nup = int((sim['th'] <= np.pi / 2).sum())        # upper cells
    pc = ax.pcolormesh((R * np.sin(T))[:nup + 1], (R * np.cos(T))[:nup + 1],
                       phys['nH'][iph][:nup],
                       norm=LogNorm(), cmap='Greys', shading='flat',
                       rasterized=True)
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    cax = make_axes_locatable(ax).append_axes('right', size='4%', pad=0.08)
    fig.colorbar(pc, cax=cax).set_label(r'$n_{\rm H}\,[\rm cm^{-3}]$')
    # .ems columns: depth, Ha, Hb, C IV 1548, 1551, Mg II 2796, 2803
    lines = [((1,), 'crimson', r'$\rm H\alpha$'),
             ((5, 6), 'seagreen', r'$\rm Mg\,II$'),
             ((3, 4), 'royalblue', r'$\rm C\,IV$')]
    # RGB code of each ray point: R = H alpha, G = Mg II, B = C IV, each
    # line's emissivity as log over RGB_DEX below its maximum over all the
    # plotted rays; opacity = the brightest channel, so non-emitting gas
    # (the hot, thin layers) is not drawn
    RGB_DEX = 4.0
    # zoom on the inner rim, where the dense rays emit
    axz = ax.inset_axes([0.6, 0.38, 0.38, 0.58])
    axz.pcolormesh((R * np.sin(T))[:nup + 1], (R * np.cos(T))[:nup + 1],
                   phys['nH'][iph][:nup], norm=pc.norm, cmap='Greys',
                   shading='flat', rasterized=True)
    axz.set_xlim(0.47, 0.62)
    axz.set_ylim(-0.005, 0.09)
    axz.set_aspect('equal')
    axz.tick_params(which='both', direction='in', labelsize=10)
    ax.indicate_inset_zoom(axz, edgecolor='purple')
    allems = {n: np.loadtxt(os.path.join(rdir, n + '.ems')) for n in names}
    jmax = [max(e[:, list(col)].sum(axis=1).max() for e in allems.values())
            for col, _, _ in lines]
    for k, name in enumerate(names):
        ith = int(name[-2:])
        th = sim['th'][ith]
        ovr = np.loadtxt(os.path.join(rdir, name + '.ovr'),
                         usecols=(0, 1, 3))           # depth, Te, hden
        depth, Te, nh = ovr.T
        r_end = (sim['rf'][0] * r0 + depth[-1]) / r0
        ax.plot([sim['rf'][0] * np.sin(th), r_end * np.sin(th)],
                [sim['rf'][0] * np.cos(th), r_end * np.cos(th)], '-',
                color='purple', lw=0.8)
        e = allems[name]
        rgb = np.stack([np.clip(1 + np.log10(
            e[:, list(col)].sum(axis=1) / jm + 1e-300) / RGB_DEX, 0, 1)
            for (col, _, _), jm in zip(lines, jmax)], axis=1)
        alpha = rgb.max(axis=1)
        m = alpha > 0
        rr = (sim['rf'][0] * r0 + e[m, 0]) / r0
        o = np.argsort(alpha[m])                  # brightest on top
        cols = np.column_stack([rgb[m][o] / alpha[m][o, None], alpha[m][o]])
        for x in (ax, axz):
            x.scatter((rr * np.sin(th))[o], (rr * np.cos(th))[o], s=90,
                      c=cols, edgecolors='none', zorder=5)
        axz.plot([sim['rf'][0] * np.sin(th), r_end * np.sin(th)],
                 [sim['rf'][0] * np.cos(th), r_end * np.cos(th)], '-',
                 color='purple', lw=0.8)
        ax.text(r_end * np.sin(th) + 0.03, r_end * np.cos(th),
                rf'$\rm {ith}$', color='purple', fontsize=14,
                va='center')
        ems = np.loadtxt(os.path.join(rdir, name + '.ems'))
        a = fig.add_subplot(gs[1, k])
        b = fig.add_subplot(gs[2, k], sharex=a)
        d = ems[:, 0]
        for col, c, lab in lines:
            a.plot(d, ems[:, list(col)].sum(axis=1), '-', color=c, lw=2,
                   label=lab)
        a.set_xscale('log')
        a.set_yscale('log')
        top = max(ems[:, list(col)].sum(axis=1).max()
                  for col, _, _ in lines)
        a.set_ylim(top * 1e-4, top * 3)
        plt.setp(a.get_xticklabels(), visible=False)
        if k == 0:
            a.set_ylabel(r'$4\pi j\,[\rm erg\,s^{-1}\,cm^{-3}]$')
            a.legend(fontsize=13, frameon=False, loc='lower left')
        a.set_title(rf'$\rm ray~{ith}:~\theta={np.degrees(th):.1f}^\circ$',
                    fontsize=16)
        b.plot(depth, nh, '-', color='0.3', lw=2)
        b.set_yscale('log')
        b.set_xlabel(r'$\rm depth\,[cm]$')
        c = b.twinx()
        c.plot(depth, Te, '-', color='darkorange', lw=2)
        c.set_yscale('log')
        if k == 0:
            b.set_ylabel(r'$n_{\rm H}\,[\rm cm^{-3}]$')
        if k == len(names) - 1:
            c.set_ylabel(r'$T_e\,[\rm K]$', color='darkorange')
        for x in (a, b, c):
            x.tick_params(which='both', direction='in', labelsize=12)
        c.tick_params(axis='y', colors='darkorange')
    for i, (_, c, lab) in enumerate(lines):
        ax.text(0.03, 0.95 - 0.08 * i, lab, color=c, transform=ax.transAxes,
                va='top', fontsize=15)
    ax.set_aspect('equal')
    ax.set_xlabel(r'$R\,[r_0]$')
    ax.set_ylabel(r'$z\,[r_0]$')
    ax.tick_params(which='both', direction='in', labelsize=12)
    pdir = os.path.join(out, 'plots')
    os.makedirs(pdir, exist_ok=True)
    path = os.path.join(pdir, f'rays_phi{iph:03d}.png')
    fig.savefig(path, dpi=130, bbox_inches='tight')
    plt.close(fig)
    print(f"wrote {path} ({len(names)} rays)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('step', choices=['write', 'run', 'extract', 'plot'])
    ap.add_argument('--iph', type=int, default=0,
                    help='plot: phi index of the rays')
    ap.add_argument('--out', required=True)
    ap.add_argument('--dph', type=int, default=2,
                    help='take every dph-th phi (sim has 256)')
    ap.add_argument('--ncpu', type=int, default=90)
    ap.add_argument('--only', type=int, default=0,
                    help='run only the first N decks (timing tests)')
    a = ap.parse_args()
    rdir = os.path.join(a.out, 'rays')
    os.makedirs(rdir, exist_ok=True)

    if a.step == 'write':
        sim = load_sim(DUMP)
        phys = to_physical(sim, CONFIG)
        decks = ray_decks(sim, phys, a.dph)
        meta = []
        for name, iph, ith, r_in, deck in decks:
            with open(os.path.join(rdir, name + '.in'), 'w') as fh:
                fh.write(deck)
            meta.append((iph, ith, r_in))
        np.save(os.path.join(a.out, 'rays_meta.npy'), np.array(meta))
        print(f"wrote {len(decks)} decks to {rdir} (every {a.dph}. phi, "
              f"upper hemisphere; Q_ION = {Q_ION:.3g}/s)")

    elif a.step == 'run':
        names = sorted(os.path.basename(p)[:-3]
                       for p in glob.glob(os.path.join(rdir, 'ray_*.in')))
        if a.only:
            names = names[:a.only]
        print(f"running {len(names)} rays on {a.ncpu} cores with {EXE}",
              flush=True)
        t0, done = time.time(), 0
        with Pool(a.ncpu) as pool:
            for name, dt, st in pool.imap_unordered(
                    run_one, [(rdir, n) for n in names]):
                done += 1
                if st != 'skip' and (done % 50 == 0 or st != 'ok'
                                     or done == len(names)):
                    el = time.time() - t0
                    print(f"  {done}/{len(names)} {name} {st} {dt:.0f} s; "
                          f"elapsed {el / 60:.1f} min, eta "
                          f"{el / done * (len(names) - done) / 60:.0f} min",
                          flush=True)

    elif a.step == 'plot':
        plot_rays(a.out, a.iph)

    else:
        meta = np.load(os.path.join(a.out, 'rays_meta.npy'))
        spec, keep = [], []
        for k, (iph, ith, r_in) in enumerate(meta):
            name = f'ray_p{int(iph):03d}_t{int(ith):02d}'
            p = os.path.join(rdir, name + '.conA')
            if not os.path.exists(p):
                continue
            spec.append(read_conA(p))
            keep.append(k)
        wave = spec[0][:, 0]
        for s in spec:
            assert len(s) == len(wave) and np.allclose(s[:, 0], wave), \
                'wavelength mesh differs between rays'
        S = np.array(spec)
        out = os.path.join(a.out, 'rays.npz')
        np.savez_compressed(out, wave=wave, out_cont=S[:, :, 1],
                            refl_cont=S[:, :, 2], refl_line=S[:, :, 3],
                            out_line=S[:, :, 4], iph=meta[keep, 0].astype(int),
                            ith=meta[keep, 1].astype(int),
                            r_in=meta[keep, 2])
        print(f"saved {out}: {len(keep)}/{len(meta)} rays, "
              f"{len(wave)} wavelengths")


if __name__ == '__main__':
    main()
