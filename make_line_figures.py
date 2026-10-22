#!/usr/bin/env python3
"""Make every paper figure that uses the Cloudy line emissivity. The disc
figures without Cloudy (Figures 1 and 2) are made by make_disc_figures.py.
Versions made with the old ionizing-flux floor (log Phi_H = 17) are in
plots/archive/archive_logphi_floor17_2026-10-03/.

Runs, in order (each step is slow):
  1. Figure 5   make_vdm_1pillar            -> tagged figure folder
  2. Figure 6   make_vdm_frame              -> tagged figure folder
  3. Figures 7, 8 make_time_evolving_figures -> tagged figure folder
  4. Figure C1  make_vdm_weighting_compare  -> plots/
  5. Appendix   plot_inclination_comparison -> plots/
  6. Appendix   plot_height_comparison      -> plots/
  7. Appendix   test_ftrans_vdm             -> plots/
  8. RGB emission map, 100 pillars (Section 6.1) -> tagged figure folder

Nothing is copied into pillar_disc_draft/figures/; that is done after checking.

Usage:
    python3 -m pillardisk.make_line_figures [--start N] [--steps 4,8]
"""

import argparse
import os
import runpy


def make_rgb_map(config_file='config_line.yaml'):
    """100-pillar RGB emission-line map (moved here from make_disc_figures)."""
    from pillardisk.make_disc_figures import _make_figdir, make_config, move
    from pillardisk.pillar_line_cloudy import main as cloudy_main
    figdir = _make_figdir(config_file)
    tmp = make_config(config_file, make_many=True)
    try:
        cloudy_main(config_file=tmp, use_absolute_flux=False)
    finally:
        os.unlink(tmp)
    for suffix in ['_xy', '_ionizing_flux', '_logF_logphi', '_EW_xy']:
        move(f'plots/line_intensity_rgb{suffix}.png', f'{figdir}/fig_rgb_N100{suffix}.png')


def main(start=1, steps_wanted=None):
    from pillardisk import (make_vdm_1pillar, make_vdm_frame,
                            make_time_evolving_figures, make_vdm_weighting_compare,
                            plot_inclination_comparison, plot_height_comparison)
    steps = [
        ('Figure 5', lambda: make_vdm_1pillar.main()),
        ('Figure 6', lambda: make_vdm_frame.main(recompute=True)),
        ('Figures 7 and 8', lambda: make_time_evolving_figures.main(recompute=True)),
        ('Figure C1', lambda: make_vdm_weighting_compare.main()),
        ('inclination appendix figure', lambda: plot_inclination_comparison.main()),
        ('height appendix figure', lambda: plot_height_comparison.main()),
        ('transparency appendix figure',
         lambda: runpy.run_module('pillardisk.test_ftrans_vdm', run_name='__main__')),
        ('RGB emission map', make_rgb_map),
    ]
    for i, (name, run) in enumerate(steps, start=1):
        if i < start or (steps_wanted and i not in steps_wanted):
            continue
        print('\n' + '#' * 60 + f'\n# [{i}/{len(steps)}] {name}\n' + '#' * 60)
        run()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--start', type=int, default=1,
                        help='step to start from (to resume after a failure)')
    parser.add_argument('--steps', type=str, default=None,
                        help='comma-separated step numbers to run, e.g. 4,8')
    args = parser.parse_args()
    wanted = {int(x) for x in args.steps.split(',')} if args.steps else None
    main(start=args.start, steps_wanted=wanted)
