#!/usr/bin/env python3
"""
Time-evolving spectra figures (1-pillar and 100-pillar setups).

Runs the slow time-evolving Cloudy-weighted velocity-map computation for
both a single co-rotating pillar AND the 100-pillar random ensemble, saving
each set of figures separately into the parameter-tagged output directory
under ``pillar_disc_draft/figures``.

This is split out from ``regenerate_figures.py`` so the fast plots there can
be regenerated without paying the cost of the time-evolving computation.

Usage:
    python3 -m pillardisk.make_time_evolving_figures [config_line.yaml]
"""

import argparse
import os
import tempfile

import yaml

from pillardisk.regenerate_figures import (
    _make_figdir,
    cleanup_geometry_snapshots,
    move,
)


FIXED_TMAX_DAYS = 1200.0


def _force_common_tmax(config):
    """Set ``computation.tmax`` and disable the auto-rounding so both 1-pillar
    and 100-pillar runs share the same time axis."""
    config.setdefault('computation', {})
    config['computation']['tmax'] = FIXED_TMAX_DAYS
    config['computation']['round_tmax'] = False


def make_one_pillar_config(config_file):
    """Load ``config_file``, force ``make_many=False``, and trim every
    list-style pillar parameter to length 1 so only a single pillar is used.
    Returns the path to a temporary YAML file."""
    with open(config_file) as f:
        config = yaml.safe_load(f)

    config.setdefault('pillars', {})
    config['pillars']['make_many'] = False
    config.setdefault('plotting', {})['plot_3d_geometry'] = False
    _force_common_tmax(config)

    pillar_list_keys = ['r_pillar', 'phi_pillar', 'height',
                        'sigma_r', 'sigma_phi',
                        'modify_height', 'modify_temp', 'temp_factor']
    for key in pillar_list_keys:
        val = config['pillars'].get(key)
        if isinstance(val, list) and len(val) > 1:
            config['pillars'][key] = val[:1]

    tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False)
    yaml.dump(config, tmp, default_flow_style=False)
    tmp.close()
    return tmp.name


def make_many_pillars_config(config_file):
    """Load ``config_file`` and force ``make_many=True`` so the random
    100-pillar ensemble (controlled by N_pillar / r_mean / sig_r in the
    config) is used. Returns the path to a temporary YAML file."""
    with open(config_file) as f:
        config = yaml.safe_load(f)

    config.setdefault('pillars', {})
    config['pillars']['make_many'] = True
    config.setdefault('plotting', {})['plot_3d_geometry'] = False
    _force_common_tmax(config)

    tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False)
    yaml.dump(config, tmp, default_flow_style=False)
    tmp.close()
    return tmp.name


def _move_outputs(figdir, suffix):
    """Move standard time-evolving plot outputs to ``figdir`` with the
    given filename suffix (e.g. ``'1pillar'`` or ``'N100'``)."""
    pairs = [
        ('plots/velocity_time_cloudy.png',
         f'{figdir}/fig_velocity_time_{suffix}.png'),
        ('plots/velocity_time_cloudy_residual.png',
         f'{figdir}/fig_barberpole_{suffix}.png'),
        ('plots/velocity_time_cloudy_diff.png',
         f'{figdir}/fig_diff_{suffix}.png'),
        ('plots/velocity_time_cloudy_combined.png',
         f'{figdir}/fig_velocity_time_combined_{suffix}.png'),
        ('plots/velocity_time_cloudy_lightcurves.png',
         f'{figdir}/fig_lightcurves_{suffix}.png'),
        ('plots/geometry_t0_faceon.png',
         f'{figdir}/fig_ionflux_{suffix}.png'),
    ]
    for src, dst in pairs:
        move(src, dst)
    cleanup_geometry_snapshots()


def main(config_file='config_line.yaml'):
    from pillardisk.pillar_line_time_cloudy import main as time_cloudy_main

    figdir = _make_figdir(config_file)
    os.makedirs(figdir, exist_ok=True)
    print(f'Output directory: {figdir}')

    # ── Single-pillar variant ──────────────────────────────
    print('\n' + '=' * 50)
    print(' [1/2] Single-pillar time-evolving maps')
    print('=' * 50)

    tmp = make_one_pillar_config(config_file)
    try:
        time_cloudy_main(config_file=tmp)
    finally:
        os.unlink(tmp)
    _move_outputs(figdir, '1pillar')

    # ── 100-pillar variant ─────────────────────────────────
    print('\n' + '=' * 50)
    print(' [2/2] 100-pillar time-evolving maps')
    print('=' * 50)

    tmp = make_many_pillars_config(config_file)
    try:
        time_cloudy_main(config_file=tmp)
    finally:
        os.unlink(tmp)
    _move_outputs(figdir, 'N100')

    print('\n' + '=' * 50)
    print(f' Done. Figures saved to {figdir}/')
    print('=' * 50)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('config', nargs='?', default='config_line.yaml',
                        help='Config YAML file (default: config_line.yaml)')
    args = parser.parse_args()
    main(config_file=args.config)
