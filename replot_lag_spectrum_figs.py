"""Replot the lag-spectrum and T(r) scan figures (paper Figures 9 and 10) from
the cached scan data, skipping the slow per-scan response-map figures.

Run from the project directory:
    python3 -m pillardisk.replot_lag_spectrum_figs
"""
from pillardisk.test_lag_spectrum import load_data, plot_all

if __name__ == '__main__':
    plot_all(load_data(), plot_response_maps=False)
