# Note: the diffuse ionising field between rays (missing from 1D Cloudy)

Context: each lamp ray is a separate 1D Cloudy model (`cloudy_rays.py`):
- open geometry, thick shell, `iterate 2`;
- Cloudy follows the diffuse field only along its own ray;
- ionising photons re-emitted back out of a ray's lit face leave the model.

In 3D those photons cross the cavity and land on other gas, both other
rays' gas and gas the lamp never sees. The emulator learns the 1D answer,
so it carries the same omission. A fix belongs on the Cloudy side, not in
the emulator or the maps. Written 2026-11-03.

## How much Cloudy reflects (`.conA`, reflected face)

Per ray, reflected / incident (medians, train/h0.0, ~30 rays; h0.5 alike):

| band | share of incident photons | reflected continuum / incident | reflected lines / incident |
|---|---|---|---|
| H, 13.6-24.6 eV | 0.51 | 0.78 | 0.08 |
| He I, 24.6-54.4 eV | 0.39 | 1.02 | 0.21 |
| He II, 54-250 eV | 0.10 | 0.63 | 0.11 |
| > 250 eV | 0.004 | 0.46 | 0.15 |

- All ionising photons together: reflected / incident ~ 0.9-1.0.
- In energy, everything reflected / incident is 1.0-1.1, and the ionising
  part alone 0.7-0.98.
  - The thick shell lets almost nothing through, so what is absorbed comes
    back out the lit face, and most of it stays ionising.
  - The ~10% excess over 1 is unexplained (shell normalisation?): treat
    these numbers as +-10%.
- Much higher than the textbook ~20% (40% of H recombinations go to the
  ground state, half of those outward). Likely causes:
  - the hot gas in front (T_e ~ 1e5-3e6 K at the ray entry) re-emits EUV;
  - higher energy photons split into several lower energy ones (He I
    band > 1).

## Where it lands (`debug/diag_diffuse_field.py`)

Method:
- Each ray's reflected photons, N = (reflected / incident) x phi s0^2
  dOmega, are emitted isotropically from its skin point.
- A path is blocked if it crosses any sim cell with n_H > N_SURF.
- The ratio is diffuse / direct lamp photon flux (photon density, as
  Cloudy sees it).

Outputs: `/data2/jhuang/runs/cloudy_rays/diffuse/` (npz per run;
`diffuse_surface.png`).

| run | at lit rays' skins: median [84%] | lit upper surface: median (16-84%) | shadowed upper surface: median (16-84%) | lit share of upper surface |
|---|---|---|---|---|
| h0.0 | 0.09 [0.11] | 0.11 (0.08-0.14) | 0.04 (0.01-0.08) | 19% |
| h0.2 | 0.08 [0.14] | 0.13 (0.08-0.18) | 0.03 (0.01-0.06) | 63% |
| h1.0 | 0.06 [0.16] | 0.14 (0.07-0.26) | 0.03 (0.01-0.07) | 95% |

- **Occlusion matters.** Without it, the lit-skin ratio would be about
  0.45-0.6.
- **Shadowed gas** is the "unblocked lamp flux" column: the diffuse light
  relative to what the lamp would give there if nothing blocked it.

Reading:
- **Lit gas:** the diffuse field adds ~10% (up to ~30% on h1.0's outer
  disc) to the lamp's ionising flux. That is a modest correction in U,
  about 0.05 dex.
- **Shadowed gas:** it gets ~1-8% of an unblocked lamp's flux instead of
  zero. Our maps have it black; with this field it would glow roughly
  1-2 dex below lit gas. This matters most at low lamps, where 81% of the
  upper surface is in shadow at h = 0.
- **Spectrum:** the diffuse field is softer (recombination photons near
  13.6 and 24.6 eV), so it ionises a thinner layer. The effect on H and
  Mg II should be larger than on C IV and He II.

Caveats:
- Emission is taken as isotropic. If it is beamed back toward the lamp
  (Lambertian about the lit face), shadowed gas gets less, by a factor of
  about 2.
- Each source sits at its ray's skin point, but part of the reflected
  light comes from the hot gas in front.
- Only one bounce is counted (the albedo ~1 makes more bounces add ~10%).
- Receivers lie on the upper surface only. At h = 0 the lit gas is
  largely the walls that face the lamp, so the lit share there is not the
  face-on lit area.
- Grids are coarse: 6 sources per ray wedge; receivers on every 8th phi
  and 4th r column.

## Next, if pursued

- Check against MOCASSIN (`/data2/jhuang/runs/mocassin/moc_edd10`), which
  includes the diffuse field in 3D: does gas in the shadowed fan or behind
  the rim come out ionised there?
- Or add the diffuse field to Cloudy:
  - lit rays: raise the incident flux by the local diffuse/direct ratio;
  - shadowed gas: give it its own Cloudy rays with a soft, recombination-
    like SED of the estimated flux.

## MOCASSIN check (2026-11-03, `debug/diag_moc_diffuse.py`)

Run: `moc_edd10`, lamp at the origin (h = 0), on a 49 x 49 x 65 grid of
~1e14 cm cells.
- "Direct front": a straight-ray case-B march on MOCASSIN's own grid
  (`march_mask`, the photon budget Q/4pi), i.e. what the lamp alone
  ionises.
- Cells beyond that front, binned by their distance from the nearest
  directly reached cell (upper half):

  | cells beyond | n | x(H+) median | x > 0.1 | x > 0.5 | Ha share | Mg II share | C IV share |
  |---|---|---|---|---|---|---|---|
  | 1 | 3457 | 1e-14 | 17% | 12% | 30% | 60% | 0% |
  | 2 | 2519 | 4e-15 | 2.5% | 1.3% | 0% | 0% | 0% |
  | 3-5 | 6645 | 2e-15 | 0% | 0% | 0.6% | 0% | 0% |
  | > 5 | 9682 | 1e-15 | 0% | 0% | 0% | 0% | 0% |

- **Top dense cell of each column (n_H > N_SURF):**
  - reached directly: 139 columns, x(H+) median 0.32 (16-84%:
    0.07-0.95);
  - beyond the front: 1693 columns, x(H+) ~ 2e-15, at the floor.
- **MOCASSIN shows no diffuse ionisation of shadowed gas.** Light beyond
  the direct front comes only from the first cell behind it, i.e. cells
  that straddle the front (or a slightly deeper front). Shadowed surfaces
  get no ionising packets at all.
- The estimate above predicted 1-8% of an unblocked lamp there. At
  n ~ 1e10 that gives x(H+) ~ 1e-2 for a whole 1e14 cm cell, which would
  be visible.

**MOCASSIN cannot settle it, because of resolution.** In Cloudy the
reflected ionising light leaves a ~1e12 cm ionised skin. In MOCASSIN that
skin is smeared over a 1e14 cm cell with x(H+) ~ 0.3. A recombination
photon emitted there faces tau_912 ~ 6e-18 x 7e9 x 1e14 ~ 4e6 inside its
own cell, so it is absorbed on the spot and never reaches other gas.
MOCASSIN's diffuse field therefore comes only from fully ionised, low
density cells, and is weak. So:
- **MOCASSIN does not support the estimate.**
- **MOCASSIN cannot rule it out either.** Its cells suppress the very
  escape that the estimate rests on.

Ways to decide:
- **Test Cloudy's reflection.** Rerun one ray with the hot entry gas cut
  (deck starting at the skin) and see whether the reflected ionising
  fraction stays ~1. That separates skin recombination photons from the
  hot gas's EUV.
- **Check how Cloudy's own escape depends on the skin.** If the reflected
  Lyman continuum comes from a skin much thinner than MOCASSIN's cells,
  that confirms that the resolution blocks it.
- **Run MOCASSIN at finer resolution** near the lit surface (expensive).

## Skin-start Cloudy test (set up 2026-11-03, `debug/diag_skin_ray.py`)

- **Test rays:** ray_p178_a09 and ray_p009_a00 from train/h0.0. Both
  enter through thin hot gas:
  - entry n_H 2e9 and 6e8, T_e 1e6 and 3e6 K;
  - ray_p009_a00's skin lies 1.2e15 cm in.
- **What changes:** each is rerun starting at its skin (the first
  n_H >= N_SURF point), with phi(H) diluted to there.
- **Question:** does the reflected ionising fraction (~0.95) survive
  without the hot gas in front?
- **Already known:** rays that enter dense gas, with nothing to cut,
  reflect almost as much. ray_p215_a06 (n_H 1.3e10, T_e 2e5 K at its
  face) reflects 0.94 of the ionising energy; ray_p044_a07 (3e11, 4e4 K)
  reflects 0.70. So the lit face itself seems to reflect most of it.
- **Run:** `python3 -m pillardisk.transient.debug.diag_skin_ray all`
  (~40 min). Output: `/data2/jhuang/runs/cloudy_rays/diffuse/skin_test/`,
  holding `skin_test.npz` and `skin_test.png`.
