#!/usr/bin/env python3
r"""
Generate the responsivity-weighted velocity-delay-map figure.

This is the companion to the emissivity-weighted ``fig_vdm_frame.png``
produced by ``regenerate_figures.py`` step 8. It reuses the same
``generate_movie_frames`` pipeline but with ``weighting='responsivity'``,
so each cell's per-line weight is

    F_line(log Phi_H) * eta,    where eta = d log F_line / d log Phi_H,

i.e. the local logarithmic slope of the Cloudy F-L curve. The resulting
$\Psi(\lambda, \tau)$ is the standard Blandford-McKee linear response
function (up to the global $1/L_{\rm ion}$ factor that is constant
across cells).

Usage:
    python3 -m pillardisk.make_vdm_responsivity [config_line.yaml]
"""

import argparse
import os
import shutil
import tempfile

import yaml

from pillardisk.regenerate_figures import _make_figdir, move


def make_n100_config(config_file):
    """Force ``make_many=True`` and disable 3D geometry; return temp path."""
    with open(config_file) as f:
        config = yaml.safe_load(f)
    config.setdefault('pillars', {})
    config['pillars']['make_many'] = True
    config.setdefault('plotting', {})['plot_3d_geometry'] = False
    tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False)
    yaml.dump(config, tmp, default_flow_style=False)
    tmp.close()
    return tmp.name


def main(config_file='config_line.yaml'):
    from pillardisk.pillar_line_time_cloudy import generate_movie_frames

    figdir = _make_figdir(config_file)
    os.makedirs(figdir, exist_ok=True)
    print(f'Output directory: {figdir}')

    print('\n' + '=' * 50)
    print(' Responsivity-weighted multi-line VDM frame')
    print('=' * 50)

    tmp_cfg = make_n100_config(config_file)
    out_dir = 'plots/_tmp_movie_resp'
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
    try:
        generate_movie_frames(
            config_file=tmp_cfg,
            n_frames=1,
            output_dir=out_dir,
            skip_geometry=True,
            weighting='responsivity',
        )
    finally:
        os.unlink(tmp_cfg)

    src = os.path.join(out_dir, 'frame_0000.png')
    dst = f'{figdir}/fig_vdm_frame_resp.png'
    move(src, dst)
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)

    print('\n' + '=' * 50)
    print(f' Done. Figure saved to {dst}')
    print('=' * 50)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config', nargs='?', default='config_line.yaml',
                        help='Config YAML file (default: config_line.yaml)')
    args = parser.parse_args()
    main(config_file=args.config)
