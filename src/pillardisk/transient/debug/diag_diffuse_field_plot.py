import numpy as np, sys
from pillardisk.transient.cloudy_rays import plot_style
plt = plot_style()
S = sys.argv[1]
runs = (('h0.0', 0.0), ('h0.2', 0.2), ('h1.0', 1.0))
fig, axs = plt.subplots(1, 3, figsize=(19, 6.6), sharex=True, sharey=True)
for ax, (r, hld) in zip(axs, runs):
    d = np.load(f'{S}/diffuse_{r}_surface.npz')
    x, y = d['r'] * np.sin(d['th']) * np.cos(d['ph']), d['r'] * np.sin(d['th']) * np.sin(d['ph'])
    lr = np.log10(np.maximum(d['diffuse'] / d['direct'], 1e-4))
    for m, mk in ((d['lit'], 'o'), (~d['lit'], 's')):
        sc = ax.scatter(x[m], y[m], c=lr[m], vmin=-2.5, vmax=-0.5, cmap='viridis',
                        marker=mk, s=26, edgecolors='k' if mk == 's' else 'none',
                        linewidths=0.4)
    ax.set_aspect('equal'); ax.set_xlim(-2.1, 2.1); ax.set_ylim(-2.1, 2.1)
    ax.set_title(rf'$\bm{{h = {hld:g}}}$\,\textbf{{ld}}', fontsize=24)
    ax.set_xlabel(r'$\bm{x\,[r_0]}$')
    ax.tick_params(which='both', direction='in', labelsize=19)
axs[0].set_ylabel(r'$\bm{y\,[r_0]}$')
axs[0].text(0.03, 0.97, r'\textbf{circles: lit}' '\n' r'\textbf{squares: shadowed}',
            transform=axs[0].transAxes, va='top', fontsize=16)
cb = fig.colorbar(sc, ax=axs, pad=0.01, fraction=0.02)
cb.set_label(r'\textbf{log diffuse / unblocked lamp flux}', fontsize=20)
fig.savefig(f'{S}/diffuse_surface.png', dpi=110, bbox_inches='tight')
print('ok')
