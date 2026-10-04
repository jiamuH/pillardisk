"""Compare the single-model Cloudy runs with the Nagao (2006) SED and the
Cloudy-shipped NGC5548.sed at one grid point (log Phi_H = 19, log n_H = 10,
solar metallicity, same stopping criteria as the paper grid).

Prints emergent line fluxes, equivalent widths referenced to the incident
continuum at 1215 A (Korista & Goad 2000 convention used in the paper), and
line-to-Hbeta ratios for both SEDs.

Run:  python3 compare_cloudy_sed_test.py
"""
import os
import numpy as np

D = '/Users/jiamuh/c23.01/my_models/sed_test_ngc5548'
RUNS = {'Nagao 2006': 'nagao_n10_phi19', 'NGC5548.sed': 'ngc5548_n10_phi19'}
# shadow-regime pair (log Phi_H = 17.5) for the illuminated-to-shadow contrast
RUNS_SHADOW = {'Nagao 2006': 'nagao_n10_phi17.5', 'NGC5548.sed': 'ngc5548_n10_phi17.5'}
LINES = [('Ly alpha', 'h  1 1215.67A'), ('C IV 1549', 'blnd 1549.00A'),
         ('He II 1640', 'he 2 1640.41A'), ('C III] 1909', 'blnd 1909.00A'),
         ('Mg II 2798', 'blnd 2798.00A'), ('H beta', 'h  1 4861.32A'),
         ('H alpha', 'h  1 6562.80A')]


def read_linelist(prefix):
    path = os.path.join(D, prefix + '_LineList_BLR_Fe2_flux.txt')
    with open(path) as f:
        head = f.readline().rstrip('\n').split('\t')[1:]
        vals = f.readline().rstrip('\n').split('\t')[1:]
    return dict(zip([h.strip() for h in head], [float(v) for v in vals]))


def incident_nufnu_1215(prefix):
    """incident nu F_nu at 1215 A from the save continuum (units Angstroms)."""
    path = os.path.join(D, prefix + '_SED.conA')
    a = np.loadtxt(path, comments='#', usecols=(0, 1))
    lam, inc = a[:, 0], a[:, 1]
    k = np.argmin(np.abs(lam - 1215.67))
    return inc[k]


def find_key(table, key):
    if key in table:
        return key
    # tolerate small differences in spacing / wavelength formatting
    tgt = key.split()
    for k in table:
        p = k.split()
        if p[0].lower() == tgt[0].lower() and p[-1] == tgt[-1]:
            return k
    raise KeyError(key)


def main():
    tabs = {n: read_linelist(p) for n, p in RUNS.items()}
    cont = {n: incident_nufnu_1215(p) for n, p in RUNS.items()}
    names = list(RUNS)
    print(f"incident nu F_nu(1215 A) [erg s^-1 cm^-2]: "
          + ", ".join(f"{n}: {cont[n]:.3e}" for n in names))
    hb = {n: tabs[n][find_key(tabs[n], 'h  1 4861.32A')] for n in names}
    print()
    print(f"{'line':<12}{'log F (Nagao)':>15}{'log F (5548)':>14}{'dlogF':>8}"
          f"{'log EW (Nagao)':>16}{'log EW (5548)':>15}{'dlogEW':>8}"
          f"{'log X/Hb (Nagao)':>18}{'log X/Hb (5548)':>17}{'d':>7}")
    for label, key in LINES:
        F = {n: tabs[n][find_key(tabs[n], key)] for n in names}
        ew = {n: F[n] * 1215.67 / cont[n] for n in names}
        rat = {n: F[n] / hb[n] for n in names}
        lf = [np.log10(F[n]) for n in names]
        le = [np.log10(ew[n]) for n in names]
        lr = [np.log10(rat[n]) for n in names]
        print(f"{label:<12}{lf[0]:>15.3f}{lf[1]:>14.3f}{lf[1]-lf[0]:>8.3f}"
              f"{le[0]:>16.3f}{le[1]:>15.3f}{le[1]-le[0]:>8.3f}"
              f"{lr[0]:>18.3f}{lr[1]:>17.3f}{lr[1]-lr[0]:>7.3f}")


def contrast():
    """change of log flux, log EW and log(X/Hbeta) from the illuminated point
    (log Phi_H = 19) to the shadow point (log Phi_H = 17.5), for each SED."""
    names = list(RUNS)
    lit = {n: read_linelist(RUNS[n]) for n in names}
    sha = {n: read_linelist(RUNS_SHADOW[n]) for n in names}
    c_lit = {n: incident_nufnu_1215(RUNS[n]) for n in names}
    c_sha = {n: incident_nufnu_1215(RUNS_SHADOW[n]) for n in names}
    hb_l = {n: lit[n][find_key(lit[n], 'h  1 4861.32A')] for n in names}
    hb_s = {n: sha[n][find_key(sha[n], 'h  1 4861.32A')] for n in names}
    print()
    print("Illuminated (log Phi_H = 19) to shadow (log Phi_H = 17.5): change in dex")
    print(f"{'line':<12}"
          f"{'flux change, Nagao SED':>24}{'flux change, NGC 5548 SED':>27}"
          f"{'EW change, Nagao SED':>22}{'EW change, NGC 5548 SED':>25}"
          f"{'X/Hbeta change, Nagao SED':>27}{'X/Hbeta change, NGC 5548 SED':>30}")
    for label, key in LINES:
        row = f"{label:<12}"
        cols = []
        for q in ('flux', 'ew', 'ratio'):
            for n in names:
                Fl = lit[n][find_key(lit[n], key)]
                Fs = sha[n][find_key(sha[n], key)]
                if q == 'flux':
                    d = np.log10(Fs / Fl)
                elif q == 'ew':
                    d = np.log10((Fs / c_sha[n]) / (Fl / c_lit[n]))
                else:
                    d = np.log10((Fs / hb_s[n]) / (Fl / hb_l[n]))
                cols.append(d)
        widths = (24, 27, 22, 25, 27, 30)
        for d, w in zip(cols, widths):
            row += f"{d:>{w}.3f}"
        print(row)


if __name__ == '__main__':
    main()
    contrast()
