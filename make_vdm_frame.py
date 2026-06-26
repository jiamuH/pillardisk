#!/usr/bin/env python3
"""Regenerate just fig_vdm_frame.png (emissivity-weighted multi-line VDM).

Equivalent to step 8 of regenerate_figures.py, but skips steps 1-7
when only the VDM colormap needs refreshing.

Usage:
    python3 -m pillardisk.make_vdm_frame [config_line.yaml]
"""

import argparse
import os
import shutil

from pillardisk.regenerate_figures import _make_figdir, make_config, move


def main(config_file='config_line.yaml'):
    from pillardisk.pillar_line_time_cloudy import generate_movie_frames

    figdir = _make_figdir(config_file)
    os.makedirs(figdir, exist_ok=True)
    print(f'Output directory: {figdir}')

    out_dir = 'plots/_tmp_movie'
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)

    tmp = make_config(config_file, make_many=True)
    try:
        generate_movie_frames(config_file=tmp, n_frames=1,
                              output_dir=out_dir, skip_geometry=True,
                              weighting='responsivity')
    finally:
        os.unlink(tmp)

    move(f'{out_dir}/frame_0000.png', f'{figdir}/fig_vdm_frame.png')
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)

    print(f' Done. Figure saved to {figdir}/fig_vdm_frame.png')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config', nargs='?', default='config_line.yaml')
    args = parser.parse_args()
    main(config_file=args.config)
