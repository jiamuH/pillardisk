"""Quick check of the seed and map cache in generate_movie_frames (Figure 6),
on a very coarse grid. Runs it twice with the same seed: the first run computes
and writes the cache, the second must read it and give the same figure.

Run:  python3 -m pillardisk.debug_vdm_frame_cache
"""
import os
import shutil
import tempfile

import yaml

from pillardisk.pillar_line_time_cloudy import generate_movie_frames

CACHE = 'plots/vdm_frame_maps_test.npz'
OUT = 'plots/_tmp_movie_test'

with open('config_line.yaml') as f:
    config = yaml.safe_load(f)
config['pillars']['make_many'] = True
config['disk']['nr'] = 40
config['disk']['nphi'] = 36
config.setdefault('computation', {})['nlambda'] = 40
config['computation']['ntau_line'] = 30
tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False)
yaml.dump(config, tmp, default_flow_style=False)
tmp.close()

if os.path.exists(CACHE):
    os.unlink(CACHE)
try:
    for run in (1, 2):
        print(f'===== run {run} =====')
        generate_movie_frames(config_file=tmp.name, n_frames=1, output_dir=OUT,
                              skip_geometry=True, weighting='responsivity',
                              seed=1, cache_file=CACHE)
        shutil.move(f'{OUT}/frame_0000.png', f'plots/vdm_frame_cache_test_run{run}.png')
finally:
    os.unlink(tmp.name)
    shutil.rmtree(OUT, ignore_errors=True)
