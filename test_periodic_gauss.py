"""Compare single-wrap vs periodic azimuthal Gaussian for pillar height."""
import numpy as np
import matplotlib.pyplot as plt
from pillardisk.pillar_disk import PillarDisk

phi = np.linspace(0, 2*np.pi, 1000)
phi_pillar = 0.0
dphi = np.mod(phi - phi_pillar + np.pi, 2*np.pi) - np.pi

height = 0.15  # pillar height in ld

fig, axes = plt.subplots(2, 2, figsize=(13, 10))

# Top row: 1D height profiles at r = r_pillar
for idx, sig in enumerate([1.0, 2.0]):
    ax = axes[0, idx]
    h_single = height * np.exp(-0.5 * (dphi / sig)**2)
    h_periodic = height * (np.exp(-0.5 * (dphi / sig)**2)
                         + np.exp(-0.5 * ((dphi - 2*np.pi) / sig)**2)
                         + np.exp(-0.5 * ((dphi + 2*np.pi) / sig)**2))
    ax.plot(np.degrees(phi), h_single, '--', color='orangered', lw=3, alpha=0.8,
            label=r'$\rm non$-$\rm periodic$')
    ax.plot(np.degrees(phi), h_periodic, '-', color='royalblue', lw=3, alpha=0.8,
            label=r'$\rm periodic$')
    ax.set_xlabel(r'$\phi~\rm [deg]$', fontsize=16)
    ax.set_ylabel(r'$h_{\rm bump}~\rm [ld]$', fontsize=16)
    ax.set_title(rf'$\sigma_\phi = {sig}~\rm rad$', fontsize=16)
    ax.legend(fontsize=13)
    ax.set_xlim(0, 360)
    ax.tick_params(direction='in', which='major', length=8, width=1.5, labelsize=13)
    ax.tick_params(direction='in', which='minor', length=5, width=1.0)
    ax.minorticks_on()

# Bottom row: 2D polar views
r = np.linspace(1, 20, 200)
phi_grid = np.linspace(0, 2*np.pi, 360)

for idx, sig_phi in enumerate([1.0, 2.0]):
    ax = fig.add_subplot(2, 2, idx + 3, projection='polar')
    disk = PillarDisk(rin=0.5, rout=20, nr=200, nphi=360,
                      h1=0.05, r0=20, beta=10, tv1=446, hlamp=0.05)
    disk.add_pillar(r_pillar=4.0, phi_pillar=0.0, height=height,
                    sigma_r=0.5, sigma_phi=sig_phi, modify_height=True)

    R, PHI = np.meshgrid(r, phi_grid, indexing='ij')
    H = disk.get_height(R, PHI)
    H_base = np.interp(R.flatten(), disk.r, disk.h_base).reshape(R.shape)
    bump = H - H_base

    c = ax.pcolormesh(PHI, R, bump, shading='auto', cmap='inferno')
    cbar = fig.colorbar(c, ax=ax, pad=0.1, shrink=0.8, aspect=20)
    cbar.set_label(r'$h_{\rm bump}~\rm [ld]$', fontsize=14)
    ax.grid(True, alpha=0.3)
    # Bold angle labels
    ticks = ax.get_xticks()
    ax.set_xticks(ticks)
    ax.set_xticklabels([f'{int(np.degrees(t))}°' for t in ticks],
                       fontsize=14, fontweight='bold')

axes[1, 0].remove()
axes[1, 1].remove()

plt.tight_layout()
plt.savefig('plots/test_periodic_gauss.png', dpi=150, bbox_inches='tight')
plt.close()
print("Saved plots/test_periodic_gauss.png")
