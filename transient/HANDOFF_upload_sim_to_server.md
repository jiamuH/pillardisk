# Handoff (Mac -> igm server): upload the Athena++ snapshot and caches

## Goal

The server session will build MOCASSIN line profiles and a predicted
observed spectrum (observer at i = 45 deg, as `INCL_DEG` in
`transient/sim_spectrum.py`). It needs the gas velocities, which live only
in the Athena++ snapshot on the Mac. Your only job: copy the files below
to the server, verify them, and report back. Do not run anything else on
the server.

## Files to upload (all under `transient/data/` on the Mac, untracked)

1. `transient/data/sim/disk.out1.00012.athdf` (the snapshot; also its
   `.xdmf` file if one exists next to it)
2. `transient/data/sim/athena_read.py` (`sim_spectrum.py` imports it from
   this folder; the server has no other copy)
3. The Cloudy caches the pipeline's line profiles need, so the server can
   overplot the pipeline spectrum for comparison:
   - `transient/data/cloudy_arm_spectra_arm_column.npz`
   - `transient/data/cloudy_arm_spectra_arm_column2.npz`
   - `transient/data/cloudy_arm_ems_arm_column2.npz`

Destination: the same relative paths in the server's clone,
`/data2/jhuang/repos/pillardisk/transient/data/` (the folder does not
exist on the server yet; create it).

## Steps

1. `git pull` (this handoff arrives that way).
2. List the files with sizes (`ls -la`) and show them to the user. If the
   total is more than ~20 GB, ask before uploading. The server has 1.9 TB
   free on `/data2`.
3. Upload with rsync over the same SSH route `transient/moc_server.py`
   uses (`HOST`, `JUMP` there):

       ssh -o ProxyJump=jhuang@ssh.strw.leidenuniv.nl jhuang@igm.strw.leidenuniv.nl "mkdir -p /data2/jhuang/repos/pillardisk/transient/data/sim"
       rsync -avz --progress -e "ssh -o ProxyJump=jhuang@ssh.strw.leidenuniv.nl" transient/data/sim/disk.out1.00012.athdf* transient/data/sim/athena_read.py jhuang@igm.strw.leidenuniv.nl:/data2/jhuang/repos/pillardisk/transient/data/sim/
       rsync -avz --progress -e "ssh -o ProxyJump=jhuang@ssh.strw.leidenuniv.nl" transient/data/cloudy_arm_spectra_arm_column.npz transient/data/cloudy_arm_spectra_arm_column2.npz transient/data/cloudy_arm_ems_arm_column2.npz jhuang@igm.strw.leidenuniv.nl:/data2/jhuang/repos/pillardisk/transient/data/

   The upload may take minutes: per the global rules, give the user the
   commands (one line each) rather than running a multi-minute transfer
   yourself, or run it in the background with progress updates if the
   user asks you to.
4. Verify: compare `md5sum` (server) with `md5 -r` (Mac) for each file.
5. Report the md5 check and sizes to the user.

## Rules

- NEVER write to `/home/jhuang` on the server (tiny quota). Only
  `/data2/jhuang/...`.
- Do NOT run `python3 -m transient.moc_server push`: it rsyncs the Mac's
  MOCASSIN source over the server's and would wipe server-only fixes
  (below). The `push` step is also obsolete for the run inputs.
- Ask the user when unsure.

## Context: what changed on the server (2026-10-09), for the Mac side

- MOCASSIN source fixes 4 and 5 (diagnostic prints only) were made on the
  server; see `transient/HANDOFF_mocassin_server.md`. They are not yet in
  the Mac's MOCASSIN tree.
- Serious bug found in the inputs: with `LPhot` and a file spectrum,
  MOCASSIN sets the luminosity from a blackbody formula
  (`LStar = 4 pi R^2 sigma TStellar^4`, `continuum_mod.f90:303-312`),
  giving 9e58 erg/s instead of 7.2e42 erg/s (lamp ~1e16 times too bright;
  T(H+) came out ~1e5 K). Fixed on the server by using `LStar 7.2400e+06`
  (1e36 erg/s units; L/Q of our SED = 8.5e-11 erg per ionizing photon).
  `transient/moc_convert.py` still writes `LPhot`: it must write `LStar`
  before any inputs are regenerated (the server session will fix it).
- The `output` keyword (per-iteration outputGas, serial, ~80% of the run
  time) was dropped; outputs are written once at the end.
- Coarse 33x33x17 run converged to 84% in 20 iterations (~12 min on 64
  cores); T(H+) = 14,600 K. Maps: `transient/moc_maps.py`.
