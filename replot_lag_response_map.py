"""Replot the continuum response-map figure for the sigma_phi scan (paper
Figure 11, plots/test_lag_spectrum_psi_vary_sp_inv.png) from the cached scan
data. The quick line plots (Figures 9 and 10) are remade as well.

Run from the project directory:
    python3 -m pillardisk.replot_lag_response_map
"""
from pillardisk.test_lag_spectrum import load_data, plot_all

if __name__ == '__main__':
    plot_all(load_data(), response_map_scans=['vary_sp'])
