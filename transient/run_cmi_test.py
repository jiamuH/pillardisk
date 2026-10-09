#!/usr/bin/env python3
"""run_cmi_test.py - run the minimal CMacIonize photoionization test on
the converted disk snapshot (inputs written by transient.cmi_convert).

Runs the task-based Monte Carlo ionization solver in
transient/data/cmi_test, logging to cmi_run.log there. The final
snapshot (disk_020.hdf5 for 20 iterations) holds the neutral hydrogen
fraction per Cartesian cell.

Run:  python3 -m transient.run_cmi_test [--threads 16]
"""

import argparse
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
EXE = '/Users/jiamuh/codes/CMacIonize/build/rundir/CMacIonize'
RUNDIR = os.path.join(HERE, 'data', 'cmi_test')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--threads', type=int, default=16)
    a = ap.parse_args()
    if not os.path.exists(os.path.join(RUNDIR, 'density.txt')):
        sys.exit("no inputs found: run  python3 -m transient.cmi_convert  "
                 "first")
    log = os.path.join(RUNDIR, 'cmi_run.log')
    print(f"running CMacIonize in {RUNDIR} with {a.threads} threads "
          f"(log: {log})...")
    with open(log, 'w') as fh:
        rc = subprocess.run([EXE, '--params', 'test.param', '--task-based',
                             '--dirty', '--threads', str(a.threads)],
                            cwd=RUNDIR, stdout=fh,
                            stderr=subprocess.STDOUT).returncode
    print("done" if rc == 0 else f"CMacIonize exited with code {rc}; "
          f"see {log}")
    sys.exit(rc)


if __name__ == '__main__':
    main()
