"""Shadow covering fraction versus radius for the current pillar disc.

Covering fraction at radius r = the fraction of the azimuth (the full 2*pi ring
at that radius) that is blocked from the direct lamp by a pillar. Computed as
the azimuthal mean of the shadow strength, where shadow strength is 1 where a
pillar fully blocks the lamp and 0 in full light. Current geometry: lamp at
50 r_g, thin pillars (height 0.5, radial width 0.25, azimuthal width 0.3 rad).
"""
import numpy as np
import matplotlib.pyplot as plt
from test_lag_spectrum import build_disk, add_pillars

H_P, SIGMA_R, SIGMA_PHI = 0.5, 0.25, 0.3
HLAMP = 50 * 0.00399

disk = build_disk(hlamp=HLAMP)
add_pillars(disk, 100, H_P, SIGMA_R, SIGMA_PHI, r_min=10.0)

r2, p2 = np.meshgrid(disk.r, disk.phi, indexing='ij')
h2 = disk.get_height(r2, p2)
mask = disk._compute_shadow_mask(r2, p2, h2)             # 1 = lit, f_trans in full shadow
# a cell counts as shadowed if it is more than half-blocked (handles overlapping
# pillar shadows without double counting, so the fraction stays in [0, 1])
shadowed = mask < 0.5 * (1.0 + disk.f_trans)
covering = shadowed.mean(axis=1)                         # fraction of azimuth in shadow

mean_cov = covering[disk.r >= 10.0].mean()
print(f"peak covering fraction = {covering.max():.3f} at r = {disk.r[np.argmax(covering)]:.1f} ld")
print(f"mean covering fraction over the pillar belt (r >= 10 ld) = {mean_cov:.3f}")
print(f"mean covering fraction over the whole disc = {covering.mean():.3f}")

plt.rcParams.update({'text.usetex': True, 'axes.linewidth': 2,
                     'font.family': 'serif', 'font.weight': 'heavy',
                     'font.size': 20})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{bm} \boldmath'

fig, ax = plt.subplots(figsize=(9, 7))
ax.plot(disk.r, covering, '-', color='navy', lw=3)
ax.set_xlim(0, disk.rout)
ax.set_ylim(bottom=0)
ax.set_xlabel(r'$r~\rm [light~days]$', fontsize=20)
ax.set_ylabel(r'$\rm shadow~covering~fraction$', fontsize=20)
ax.tick_params(which='major', direction='in', length=8, width=1.5,
               top=True, right=True, labelsize=14)
ax.tick_params(which='minor', direction='in', length=4, width=1.0,
               top=True, right=True)
ax.minorticks_on()
plt.tight_layout()
plt.savefig('plots/test_shadow_covering.png', dpi=200, bbox_inches='tight')
plt.close()
print("Saved plots/test_shadow_covering.png")
