"""Low-resolution version of the two-pillar velocity-delay maps, for testing the
figure layout quickly. Same config as the paper figure but with a coarse disc
grid (nr = 150, nphi = 90, 16 times fewer cells), written to a separate cache so
it never mixes with the full-resolution run.

Run:  python3 -m pillardisk.debug_vdm_lowres
"""
import os
import tempfile

import yaml

from pillardisk.regenerate_figures import _run_cloudy_vdm

CONFIG = 'config_line.yaml'
CACHE = 'plots/vdm_manual_pillars_maps_lowres.npz'
FIGURE = 'plots/velocity_delay_map_lowres.png'

with open(CONFIG) as f:
    config = yaml.safe_load(f)
config['pillars']['make_many'] = False
config['disk']['nr'] = 150
config['disk']['nphi'] = 90
config.setdefault('plotting', {})['plot_3d_geometry'] = False

tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False)
yaml.dump(config, tmp, default_flow_style=False)
tmp.close()
try:
    _run_cloudy_vdm(tmp.name, filename=FIGURE, cache_file=CACHE)
finally:
    os.unlink(tmp.name)
