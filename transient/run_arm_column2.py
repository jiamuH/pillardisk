#!/usr/bin/env python3
"""run_arm_column2.py - one-command driver for the extended Cloudy grid.

Copies the deck (transient/cloudy/arm_column2.in) into the Cloudy model
directory, runs the 196-model grid (Cloudy forks across all cores on its
own, ~20-45 min wall-clock), and extracts the spectra cache. The
pipeline then adopts the extended grid automatically.

Run:  python3 -m transient.run_arm_column2
"""

import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DECK = os.path.join(HERE, 'cloudy', 'arm_column2.in')
BASE = '/Users/jiamuh/c23.01/my_models/arm_column2'
EXE = '/Users/jiamuh/c23.01/source/cloudy.exe'


def main():
    os.makedirs(BASE, exist_ok=True)
    shutil.copy2(DECK, BASE)
    print(f"deck copied to {BASE}")
    print("running the Cloudy grid (196 models; Cloudy forks across "
          "cores)...")
    rc = subprocess.run([EXE, '-r', 'arm_column2'], cwd=BASE).returncode
    if rc != 0:
        sys.exit(f"Cloudy exited with code {rc}; "
                 f"check {BASE}/arm_column2.out")
    print("extracting the spectra cache...")
    rc = subprocess.run(
        [sys.executable, os.path.join(HERE, 'extract_cloudy_arm.py'),
         'arm_column2']).returncode
    if rc == 0:
        print("done - the pipeline now uses the extended grid "
              "(data/cloudy_arm_spectra_arm_column2.npz)")
    sys.exit(rc)


if __name__ == '__main__':
    main()
