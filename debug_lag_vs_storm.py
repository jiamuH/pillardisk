"""Referee reply, Section 8 bullet 4: how do the cached model lag spectra compare
with the NGC 5548 STORM lags (Fausnaugh et al. 2016)?

For the bare bowl, the ripple model and every pillar model in the four scans:
  - chi^2 against the 17 UV/optical bands (errors symmetrized), and against the
    15 bands left when U and u are excluded;
  - the best-fitting alpha when the model curve is compared with
    tau = alpha [(lambda/1367)^(4/3) - 1]  (thin disc: alpha = 0.14 day).

Run:  python3 -m pillardisk.debug_lag_vs_storm
"""
import numpy as np

from pillardisk.test_lag_spectrum import (ALPHA_THIN_DISC, LAMBDA_REF, load_data,
                                          load_or_compute_ripple, load_storm_lags,
                                          scans)

data = load_data()
wavelengths = data['wavelengths']
band, lam, lag, err_plus, err_minus = load_storm_lags()
err = 0.5 * (err_plus + err_minus)
no_balmer = ~np.isin(band, ['U', 'u'])
shape = (wavelengths / LAMBDA_REF) ** (4. / 3.) - 1.
fit_range = (wavelengths > 1100.) & (wavelengths < 9500.)


def summarize(label, tau):
    rel = tau - np.interp(LAMBDA_REF, wavelengths, tau)
    model = np.interp(lam, wavelengths, rel)
    chi2_all = np.sum(((lag - model) / err) ** 2)
    chi2_nb = np.sum((((lag - model) / err) ** 2)[no_balmer])
    alpha = np.sum(rel[fit_range] * shape[fit_range]) / np.sum(shape[fit_range] ** 2)
    print('%-26s chi2 = %6.1f (17 bands)  %6.1f (15 bands, no U/u)   alpha = %.2f d'
          % (label, chi2_all, chi2_nb, alpha))


print('thin disc reference: alpha = %.2f d' % ALPHA_THIN_DISC)
summarize('thin disc (alpha = 0.14)', ALPHA_THIN_DISC * shape)
summarize('bowl (h_LP = 0.5)', data['tau_bowl_0.5'])
summarize('ripple', load_or_compute_ripple(hlamp=0.5)['tau_ripple'])
for scan_name, scan in scans.items():
    print('\n' + scan_name)
    for i, value in enumerate(scan['values']):
        summarize('  %s = %s' % (scan_name, value), data['tau_%s_%d' % (scan_name, i)])
