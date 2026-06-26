#!/usr/bin/env python3
r"""
Single-pillar Cloudy-weighted H-alpha velocity-delay map figure.

Runs the same computation that step 6 of ``regenerate_figures.py`` did:
build the disk with one manual pillar, compute the H-alpha
\Psi(\lambda, \tau) with and without the pillar via Cloudy LOC weights,
and save the resulting figure into the parameter-tagged
pillar_disc_draft/figures/ directory.

Usage:
    python3 -m pillardisk.make_vdm_1pillar [config_line.yaml]
"""

import argparse
import os

from pillardisk.regenerate_figures import (
    _make_figdir,
    _run_cloudy_vdm,
    make_config,
    move,
)


def main(config_file='config_line.yaml'):
    figdir = _make_figdir(config_file)
    os.makedirs(figdir, exist_ok=True)
    print(f'Output directory: {figdir}')

    print('\n' + '=' * 50)
    print(' Manual-pillar Cloudy H-alpha VDM (2 pillars from config)')
    print('=' * 50)

    tmp = make_config(config_file, make_many=False)
    try:
        _run_cloudy_vdm(tmp, filename='plots/velocity_delay_map.png')
    finally:
        os.unlink(tmp)
    move('plots/velocity_delay_map.png', f'{figdir}/fig_vdm_1pillar.png')

    print('\n' + '=' * 50)
    print(f' Done. Figure saved to {figdir}/fig_vdm_1pillar.png')
    print('=' * 50)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config', nargs='?', default='config_line.yaml',
                        help='Config YAML file (default: config_line.yaml)')
    args = parser.parse_args()
    main(config_file=args.config)
