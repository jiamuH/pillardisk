"""Replot the manual-pillar (two-pillar) H-alpha velocity-delay figure from the
cached maps written by regenerate_figures._run_cloudy_vdm, without recomputing.

Run:  python3 -m pillardisk.replot_vdm_manual_pillars [cache.npz] [output.png]
"""
import sys

import numpy as np

from pillardisk.pillar_line import plot_velocity_delay_map


def main(cache_file='plots/vdm_manual_pillars_maps.npz',
         filename='plots/velocity_delay_map_replot.png'):
    d = np.load(cache_file)
    plot_velocity_delay_map(d['lam'], d['tau'], d['psi'], float(d['lambda0']),
                            psi_map_no_pillars=d['psi_no'],
                            psi_per_pillar=list(d['psi_per_pillar']),
                            filename=filename)


if __name__ == '__main__':
    main(*sys.argv[1:3])
