# Handoff: MOCASSIN photoionization run on the igm server

## Goal

Validate our photoionization post-processing of an Athena++ disk snapshot
with an independent three-dimensional Monte Carlo photoionization code.
Our pipeline marches straight rays from a central ionizing lamp and counts
photons against recombinations (case B, on-the-spot), then looks up
precomputed Cloudy slab models. A first check with CMacIonize (hydrogen
only, diffuse field on) already agrees: 97% of cells classified the same,
ionization front at most one cell (0.02 r0) shallower in CMacIonize. This
run adds what CMacIonize cannot do: the real AGN spectrum including
X-rays, helium and metals, and a self-consistent temperature. It tests the
partially ionized zone behind the hydrogen front and the placement of
C IV and Mg II.

Done means a converged MOCASSIN run whose `output/` folder is copied
back to the Mac for analysis (`python3 -m transient.moc_server pull`
on the Mac).

## Where things are on the server

- `/data2/jhuang/repos/MOCASSIN-2.0` : MOCASSIN 2.0 source (cloned from
  github.com/mocassin/MOCASSIN-2.0, with five local fixes below).
- `/data2/jhuang/runs/mocassin/moc_test` : the run directory.
  - `input/input.in` : MOCASSIN input file.
  - `input/density.dat` : hydrogen density on a Cartesian grid of
    65 x 65 x 33 cells (x, y, z in cm, n_H in cm^-3), resampled from
    Athena++ dump `disk.out1.00012.athdf`. 76,852 cells contain gas; cells
    outside the simulated volume have zero density (MOCASSIN treats them
    as empty). Odd cell counts so one cell is centred on the lamp.
  - `input/agn_sed.dat` : the AGN spectrum from our Cloudy decks
    (wavelength in Angstrom, f_lambda; MOCASSIN rescales it to `LPhot`).
  - `input/abun.in` : solar abundances (Asplund et al. 2009) for H, He,
    C, N, O, Ne, Mg, Si, S. Iron is switched OFF (see below).
  - `data`, `dustData` : links to the atomic data in the source folder.
  - `output/` : results.

## Input settings (input/input.in)

- Lamp at the origin, `LPhot 8.5e16` (units of 1e36 photons per second, so
  Q = 8.5e52 ionizing photons per second).
- Frequency grid `nuMin 1.001e-5` to `nuMax 300` Rydberg (about 4 keV),
  1500 bins.
- `nPhotons 10000000`, `autoPackets 0.20 2. 80000000` (photon count
  doubles automatically when convergence stalls, up to 8e7).
- `maxIterateMC 20 95.` (stop after 20 iterations or when 95% of cells
  have converged), `convLimit 0.05`, `nstages 7`, `Rin` = 0.48 r0.
- `output` (write output files after every iteration).
- Keyword spelling matters: the parser is case-sensitive (`nstages`, not
  `nStages`), and `Rin` is required.

## Why iron is off

With iron, MOCASSIN stops in the line-emission step: the Fe VI atomic
data has 19 energy levels, but the compile-time limit `nForLevels`
(`source/constants_mod.f90`) is 17. Raising it is what the error message
suggests, but that limit also sizes two arrays stored per grid cell (line
packets and line probabilities, one entry per possible line, growing as
the limit squared, with Fe II alone about 20,000 lines). Each MPI process
holds a full grid copy, so on the Mac (137 GB) this was not feasible. The
server's ~2 TB of memory makes an iron run possible later: raise
`nForLevels` (for example to 25), rebuild, and regenerate the inputs on
the Mac with `python3 -m transient.moc_convert --with-iron`. Do the
iron-free run first.

## Local source fixes

Fixes 1 to 3 are in the pushed source; fixes 4 and 5 were made on the
server (2026-10-09) and must be copied to the Mac's source tree.

1. `source/photon_mod.f90`, line 32: argument `gpLoc` of
   `energyPacketDriver` changed from `intent(inout)` to `intent(in)`. It is
   only read inside the subroutine, but a loop index is passed to it,
   which modern gfortran rejects.
2. `source/grid_mod.f90`, start of `setMotherGrid`: added
   `nullify(MdMg, HdenTemp, NdustTemp, dustAbunIndexTemp, twoDscaleJTemp)`.
3. `source/grid_mod.f90`, start of the sub-grid routine:
   added `nullify(HdenTemp, NdustTemp, dustAbunIndexTemp)`.

Fixes 2 and 3: these pointers were never initialized, so `associated()`
returned garbage, and the code freed memory it never allocated ("pointer
being freed was not allocated", found with lldb on a debug build).

4. `source/photon_mod.f90`, line 22: `totalEscaped` changed from `integer`
   to `real`. It sums the real-valued `escapedPackets`, so the integer
   overflowed and every iteration printed
   `total Escaped Packets : -2147483648`.
5. `source/set_input_mod.f90`, keywords `LStar` and `LPhot`: added
   `Lstar = 0.` right after `allocate(Lstar(0:1))`. Only `Lstar(1)` is
   read, so `Lstar(0)` was uninitialized and `writeSED` printed garbage
   (e.g. `-1.04E+38`) in "Total energy radiated out of the nebula".

Fixes 4 and 5 only affect printed diagnostics, not the physics. The first
server run (started 2026-10-09 04:23) used the binary built before them.

## Building

On the Mac (gfortran 16) the build needed
`-fallow-argument-mismatch -fallow-invalid-boz` (new strictness about
old Fortran) and `make -B` (the makefile rule has no dependencies, so
plain `make` skips rebuilding when a binary exists). The server's
gfortran 8.5 should not need those two flags, and they do not exist in
gfortran 8, so do not pass them. The build command:

    cd /data2/jhuang/repos/MOCASSIN-2.0 && make -B mocassin F90=/usr/lib64/openmpi/bin/mpif90 OPT1="-fno-range-check -O2" > build.log 2>&1

The system OpenMPI is version 4.1.1 in `/usr/lib64/openmpi/bin` (not on
the default PATH). The login shell is csh: for redirections like `2>&1`,
run commands through `bash -lc '...'` or start bash first.

## Running

    cd /data2/jhuang/runs/mocassin/moc_test && setsid nohup /usr/lib64/openmpi/bin/mpirun -np 32 /data2/jhuang/repos/MOCASSIN-2.0/mocassin > run.log 2>&1 < /dev/null &

MOCASSIN must be started from the run directory (it reads `input/` and
`data/` relative to the current directory). Choose `-np` with the shared
machine in mind: igm has 128 cores and about 2 TB RAM (about 184 GB in
use when last checked). 32 is the default suggestion.

## Timing so far (Mac, 4 processes, 10,000 photons, 1 iteration)

- Setup and photon transfer: about 4 minutes.
- Cell update (ionization and thermal balance for 77,000 cells): about
  86 minutes. This dominates and does NOT shrink with fewer photons.
- Writing outputs: about 10 minutes (`grid2.out` alone is 129 MB).

Scaled to 32 processes, roughly 10 to 15 minutes of cell update per
iteration, but the server's processors differ: time the first iteration.
If writing outputs every iteration turns out slow, consider dropping the
`output` keyword and recovering the output at the end with the
`mocassinOutput` driver (`make -B mocassinOutput`, same flags).

## Rules for the session

- NEVER write anything to the home directory (/home/jhuang, small quota). All code goes in /data2/jhuang/repos/, runs and outputs under /data2/jhuang/ (see /data2/jhuang/CLAUDE.md).

- Anything expected to run longer than about one minute is the user's to
  start: hand over a one-line command. A small photon count does not make
  a MOCASSIN run short (the per-cell update cost is fixed).
- Never chain `rm` (especially with wildcards) after a `cd`; MOCASSIN
  overwrites its own outputs, so cleanup is not needed.
- Ask when unsure rather than picking an interpretation.

## Checking progress

    grep -E "iterateMC: (Starting|updateCell out)|convergence" /data2/jhuang/runs/mocassin/moc_test/run.log | tail
    pgrep -c -x mocassin

From the Mac the same is `python3 -m transient.moc_server status`.
