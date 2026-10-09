#!/usr/bin/env python3
"""moc_server.py - run the MOCASSIN disk photoionization test on the igm
server (128 cores, ~2 TB RAM), one short command per step.

Steps (run in order):
  push    copy the patched MOCASSIN source and the run directory
          (transient/data/moc_test, written by transient.moc_convert)
          to the server; link the atomic data into the run directory
  build   compile MOCASSIN on the server with the system OpenMPI
          (gfortran 8.5 needs none of the compatibility flags used on
          the Mac)
  run     start the run in the background on the server (survives
          logging out); --np sets the number of MPI processes
  status  show the last progress lines of the server log
  pull    copy the results back to transient/data/moc_test/output_server

Remote layout:  ~/codes/MOCASSIN-2.0   and   ~/mocassin_runs/moc_test

Run:  python3 -m transient.moc_server push
      python3 -m transient.moc_server build
      python3 -m transient.moc_server run --np 32
      python3 -m transient.moc_server status
      python3 -m transient.moc_server pull
"""

import argparse
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
HOST = 'jhuang@igm.strw.leidenuniv.nl'
JUMP = 'jhuang@ssh.strw.leidenuniv.nl'
SSH = f'ssh -o ProxyJump={JUMP}'
LOCAL_SRC = '/Users/jiamuh/codes/MOCASSIN-2.0/'
LOCAL_RUN = os.path.join(HERE, 'data', 'moc_test') + '/'
REMOTE_SRC = 'codes/MOCASSIN-2.0'
REMOTE_RUN = 'mocassin_runs/moc_test'
MPI = '/usr/lib64/openmpi/bin'


def sh(cmd):
    print(f"$ {cmd}")
    rc = subprocess.run(cmd, shell=True).returncode
    if rc != 0:
        sys.exit(f"step failed (exit code {rc})")


def remote(cmd):
    # the login shell on igm is csh: run everything through bash
    sh(f"{SSH} {HOST} \"bash -lc '{cmd}'\"")


def push():
    remote(f'mkdir -p ~/{REMOTE_SRC} ~/{REMOTE_RUN}/output')
    excl = ' '.join(f"--exclude '{e}'" for e in (
        '.git', '*.o', '*.mod', 'mocassin', 'mocassin_debug', '*.dSYM'))
    sh(f'rsync -avz --progress -e "{SSH}" {excl} {LOCAL_SRC} '
       f'{HOST}:{REMOTE_SRC}/')
    # the local data/dustData entries are links to the Mac source tree:
    # skip them and link the server's copy instead
    sh(f'rsync -avz --progress -e "{SSH}" --exclude data --exclude '
       f'dustData --exclude output {LOCAL_RUN} {HOST}:{REMOTE_RUN}/')
    remote(f'cd ~/{REMOTE_RUN} && ln -sfn ~/{REMOTE_SRC}/data data && '
           f'ln -sfn ~/{REMOTE_SRC}/dustData dustData && ls -l')


def build():
    remote(f'cd ~/{REMOTE_SRC} && make -B mocassin F90={MPI}/mpif90 '
           f'OPT1=\\\"-fno-range-check -O2\\\" > build.log 2>&1; '
           f'tail -3 build.log; ls -l mocassin')


def run(np_):
    remote(f'cd ~/{REMOTE_RUN} && setsid nohup {MPI}/mpirun -np {np_} '
           f'~/{REMOTE_SRC}/mocassin > run.log 2>&1 < /dev/null & '
           f'sleep 2; echo started; tail -2 ~/{REMOTE_RUN}/run.log')


def status():
    remote(f'cd ~/{REMOTE_RUN} && ls -l output | tail -5; '
           f'grep -E \\\"iterateMC: (Starting|updateCell out)|convergence\\\" '
           f'run.log | tail -6; tail -2 run.log; '
           f'echo mocassin processes running: \\$(pgrep -c -x mocassin)')


def pull():
    dst = os.path.join(LOCAL_RUN, 'output_server') + '/'
    os.makedirs(dst, exist_ok=True)
    sh(f'rsync -avz --progress -e "{SSH}" {HOST}:{REMOTE_RUN}/output/ '
       f'{dst}')
    sh(f'rsync -avz -e "{SSH}" {HOST}:{REMOTE_RUN}/run.log {dst}')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('step', choices=['push', 'build', 'run', 'status',
                                     'pull'])
    ap.add_argument('--np', type=int, default=32,
                    help='MPI processes for the run step (default 32)')
    a = ap.parse_args()
    {'push': push, 'build': build, 'status': status, 'pull': pull,
     'run': lambda: run(a.np)}[a.step]()


if __name__ == '__main__':
    main()
