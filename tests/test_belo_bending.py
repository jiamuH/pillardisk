"""Sanity checks for the Beloborodov direct-ray bending implementation."""
import numpy as np

# Check 1: r_s -> 0 recovers flat space (alpha = psi, cos_e = s_t)
psi = np.linspace(0.01, np.pi / 2, 50)
cos_psi = np.cos(psi)
for rs_over_R in [0.0, 0.3, 0.5]:
    cos_alpha = 1.0 - (1.0 - cos_psi) * (1.0 - rs_over_R)
    sin_alpha = np.sqrt(1.0 - cos_alpha**2)
    # in-plane target at elevation eta_t: psi = eta_t, s_t = sin(psi)
    s_t = np.sin(psi)
    cos_e = sin_alpha * s_t / np.sin(psi)   # = sin_alpha
    if rs_over_R == 0.0:
        assert np.allclose(cos_e, np.sin(psi), atol=1e-12), 'flat limit fails'
        print('flat-space limit: cos_e == sin(psi)  OK')

# Check 2: grazing suppression ~ (1 - rs/R)^{3/2}
print('\n grazing-ray kernel ratio vs (1-rs/R)^1.5:')
psi0 = 0.01  # near-plane target
for rs_over_R in [0.1, 0.3, 0.5]:
    cos_psi0 = np.cos(psi0)
    cos_alpha0 = 1.0 - (1.0 - cos_psi0) * (1.0 - rs_over_R)
    sin_alpha0 = np.sqrt(1.0 - cos_alpha0**2)
    kernel_bend = sin_alpha0 * (1.0 - rs_over_R)   # cos_e_local * Jacobian
    kernel_flat = np.sin(psi0)
    print(f'  rs/R={rs_over_R}: ratio={kernel_bend/kernel_flat:.4f}, '
          f'(1-rs/R)^1.5={(1-rs_over_R)**1.5:.4f}')
