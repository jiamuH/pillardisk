# Note: line escape, turbulence and velocities in the per-ray Cloudy runs

Context: `transient/cloudy_rays.py`, one Cloudy c23.01 model per lamp ray
through the Athena++ snapshot (`disk.out1.00012.athdf`), compared with
MOCASSIN (`/data2/jhuang/runs/mocassin/moc_edd10`). Written 2026-10-17.

## Sim velocity structure (`cloudy_rays vmaps`)

- Command: `python3 -m transient.cloudy_rays vmaps --out
  /data2/jhuang/runs/cloudy_rays/edd10`; writes `plots/vel_{grad,shear,
  disp,normal}.png`, meridional cuts through the star (phi 335, index
  238), the arm (phi 70, index 49, strongest Sigma excess away from the
  star, 3.6x the phi-median at R = 1) and the quietest column (phi 199,
  index 141). Black contours: n_H = 1e9, 1e10, 1e11 cm^-3 (thin to
  thick). The lit skin is between the 1e10 and 1e11 contours.
- Fields (`velocity_fields`):
  - `grad`: dv_r/dr, along an in-plane lamp ray [km/s per 0.01 r0].
  - `shear`: |sigma| = sqrt(2 sigma_ij sigma_ij), traceless strain rate;
    dominated by Keplerian r dOmega/dr; excludes the divergence.
  - `disp`: 1D dispersion of (v_r, v_th, v_ph) over 3x3x3 cells [km/s].
    Not sub-grid turbulence: the resolved gradients over ~1e14 cm.
  - `normal`: n.e.n, gradient of the normal velocity along the surface
    normal n = -grad n_H / |grad n_H| (full strain rate, trace
    included). Red = layer stretched (outflow base), blue = squeezed.
- Grid: 128 r x 64 th x 256 ph; dr = 0.0117 r0, r dth = 0.005-0.02 r0,
  r dph = 0.012-0.05 r0. "Per 0.01 r0" is ~1 pixel at R ~ 1.
- Lit skin numbers (16/50/84%): |dv_r/dr| 0.3/2.8/26, |n.e.n| 0.15/1.4/21
  km/s per 0.01 r0; |n.e.n| x 1e12 cm = 0.003/0.03/0.4 km/s; skin normal
  is vertical (|n_r| ~ 0.04). Away from the star the skin is nearly
  hydrostatic (v_z ~ -1 +- 30 km/s). Above the star: squeezed (n.e.n ~
  -9). Envelope around the star: X pattern, rising (v_z ~ +14), T ~
  6500 K. Hot layer (n < 1e9): +-300 km/s streams; n is ill defined
  there (checkerboard), so only coherent patches mean anything.
- Conclusion: resolved velocities shift a photon by << the thermal width
  (Mg II 1.7, H 13 km/s) across a ~1e12 cm emitting layer; only the
  star envelope / arm reach ~0.5-2 km/s. The sim says nothing about
  sub-cell (< 1e14 cm) motions.

## Turbulence in Cloudy

- Stock c23.01: `turbulence` (constant), `vlaw alpha` (power law in
  radius, alpha must be < 0), `turbulence ... dissipate` (exp decay, adds
  wave heating). No table form. `wind` solves its own dynamics (LVG), so
  an imposed v(s) is impossible.
- A `turbulence table depth` patch would be easy: `vlaw` resets
  `DoppVel.TurbVel` in `ZoneStart` (`source/zone_startend.cpp:55`) and
  every line width comes from it via `GetDopplerWidth`. Dropped: the sim
  velocities give vturb ~ 0 anyway.
- The current 100 km/s is supersonic and not supported by the sim;
  sub-cell turbulence, if any, is plausibly subsonic (c_s ~ 10 km/s at
  1e4 K). MOCASSIN has no turbulence (thermal widths only).

## Line escape geometry in Cloudy (important)

- 1D two-stream escape: P = [P(tau_in) + P(T - tau_in)] / 2
  (`docs/latex/hazy2/linedetl.tex:238`): toward the lamp, and deeper
  along the ray. The deeper direction is the second term; its velocity
  gradient is dv_r/dr, not n.e.n. The real shortest escape route
  (vertical, through the skin) does not exist in the 1D model.
- Deep route is closed physically: Sobolev length along the ray ~3e13 cm
  in the skin, Mg+ column ~1e20 over it, Mg II tau0 ~ 1e8 even then.
- OUR DECKS RUN ONE ITERATION (no `iterate`; Cloudy warns "This is the
  last iteration ... Use the ITERATE command"). On iteration 1, T is
  unknown and starts ~0, so tau_out = |T - tau_in| ~ tau_in
  (`tau_from_here`, `source/rt_escprob.cpp:47`): the deep side is as
  open as the lamp side. Only H Lya gets a guessed T (from `stop neutral
  column density`, `source/rt_tau_init.cpp:325`). Effect: thick lines
  (Mg II, Ha, C IV) escape up to 2x too easily.

## Horne & Marsh (1986, MNRAS 218, 761): shear and U vs V profiles

- Shear broadening across an emitting layer of thickness H:
  dv = dv_th [1 + Q^2 sin^2 i tan^2 i sin^2 2theta]^1/2,
  Q = (3/4)(V_K / V_th)(H/R) (eqs. 20-21); theta = disk azimuth from the
  line of sight; largest at theta = 45, 135 deg.
- Optically thin lines: isotropic, sharp peaks, shallow U-shaped valley.
  Optically thick lines: intensity ~ S_L dv (8 ln tau0)^1/2 (eq. 25),
  so they leak out along the high-shear 45 deg directions -> broader
  peaks, deeper V-shaped valley (once shear > thermal, sin i tan i >~ 1).
- Here (i = 45, sin i tan i = 0.71; R = 1, V_K ~ 5060 km/s; layer 1e12-
  1e13 cm): Q(Mg II) ~ 0.5-4.5 (dv up to ~3x thermal), Q(Ha) ~ 0.06-0.6
  (<~10%), Q with vturb = 100 km/s ~ 0.01-0.08 (no effect). Neither
  MOCASSIN (isotropic) nor our profile code includes this.
- The sim's real velocities belong here, in the profile synthesis, not
  in Cloudy: weight each emitting element by its line-of-sight width
  from |o.e.o| (strain tensor) x path through the layer; the arm and
  star then get their actual non-Keplerian shear.

## The Mg II gap (per-ray Cloudy ~100x below MOCASSIN)

- More turbulence makes thick-line flux larger (eq. 25), so vturb = 100
  already flatters Cloudy's Mg II; vturb ~ 0 makes it weaker. Turbulence
  cannot close the gap.
- Per-ray Cloudy Mg II/Ha ~ 0.4 is normal for AGN (Mg II/Hb ~ 1-2,
  Ha/Hb ~ 3-4). MOCASSIN gives ~4-5 (peak heights of its profiles, per-A
  units corrected). Its excess comes from the deep neutral midplane
  layer (|z| < 0.04, R 0.7-1.2, x(H+) ~ 0-0.01, 4000-8000 K), where
  MOCASSIN treats Mg II as optically thin (no trapping/thermalization).
- Working hypothesis: MOCASSIN, not Cloudy, is off for Mg II.

## Plan (in order)

1. MOCASSIN line totals excluding cells with x(H+) < 0.1: if Mg II drops
   toward Cloudy, the gap is a MOCASSIN artifact.
2. ~10 test rays (rim, ray 16, skin, star, arm, quiet) with `iterate`, at
   vturb = 0, 10, 100 km/s (~30 runs, user runs); compare with the
   1-iteration 100 km/s outputs; record tau0 per line.
3. Horne & Marsh anisotropic weighting in `ray_spectra` for thick lines,
   using the sim strain tensor; needs the vturb ~ 0 outputs from step 2.
