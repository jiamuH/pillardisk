#!/usr/bin/env python3
"""Regenerate every paper figure that uses the Cloudy line emissivity, after
the ionizing-flux floor was lowered from log Phi_H = 17 to 15 (2026-10-03).
The previous versions are in plots/archive_logphi_floor17_2026-10-03/.

Runs, in order (each step is slow):
  1. Figure 5   make_vdm_1pillar            -> tagged figure folder
  2. Figure 6   make_vdm_frame              -> tagged figure folder
  3. Figures 7, 8 make_time_evolving_figures -> tagged figure folder
  4. Figure C1  make_vdm_weighting_compare  -> plots/
  5. Appendix   plot_inclination_comparison -> plots/
  6. Appendix   plot_height_comparison      -> plots/
  7. Appendix   test_ftrans_vdm             -> plots/

Nothing is copied into pillar_disc_draft/figures/; that is done after checking.

Usage:
    python3 -m pillardisk.make_line_figures [--start N]
"""

import argparse
import runpy


def main(start=1):
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
    ]
    for i, (name, run) in enumerate(steps, start=1):
        if i < start:
            continue
        print('\n' + '#' * 60 + f'\n# [{i}/{len(steps)}] {name}\n' + '#' * 60)
        run()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--start', type=int, default=1,
                        help='step to start from (to resume after a failure)')
    main(start=parser.parse_args().start)
