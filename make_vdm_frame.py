#!/usr/bin/env python3
"""Regenerate just fig_vdm_frame.png (responsivity-weighted multi-line VDM).

Equivalent to step 8 of regenerate_figures.py, but skips steps 1-7
when only the VDM colormap needs refreshing.

The 100 pillars are drawn with a fixed random seed, and the maps are cached
in plots/vdm_frame_maps.npz, so a rerun after a plotting change only replots.

Usage:
    python3 -m pillardisk.make_vdm_frame [config_line.yaml] [--seed N] [--recompute]
"""

import argparse
import os
import shutil

from pillardisk.regenerate_figures import _make_figdir, make_config, move


SEED = 1
CACHE_FILE = 'plots/vdm_frame_maps.npz'


def main(config_file='config_line.yaml', seed=SEED, recompute=False):
    from pillardisk.pillar_line_time_cloudy import generate_movie_frames

    figdir = _make_figdir(config_file)
    os.makedirs(figdir, exist_ok=True)
    print(f'Output directory: {figdir}')

    out_dir = 'plots/_tmp_movie'
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)

    if recompute and os.path.exists(CACHE_FILE):
        os.unlink(CACHE_FILE)

    tmp = make_config(config_file, make_many=True)
    try:
        generate_movie_frames(config_file=tmp, n_frames=1,
                              output_dir=out_dir, skip_geometry=True,
                              weighting='responsivity', seed=seed,
                              cache_file=CACHE_FILE)
    finally:
        os.unlink(tmp)

    move(f'{out_dir}/frame_0000.png', f'{figdir}/fig_vdm_frame.png')
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)

    print(f' Done. Figure saved to {figdir}/fig_vdm_frame.png')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config', nargs='?', default='config_line.yaml')
    parser.add_argument('--seed', type=int, default=SEED,
                        help='random seed for the pillar positions')
    parser.add_argument('--recompute', action='store_true',
                        help='ignore the map cache and recompute the maps')
    args = parser.parse_args()
    main(config_file=args.config, seed=args.seed, recompute=args.recompute)
