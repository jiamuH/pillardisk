"""What the embedded-star budget demands, and whether subdividing the mass is free.

Only the total embedded mass per pillar enters the flux-flux model, because
L_Edd is proportional to M. This checks the global budget that total implies,
and the point at which splitting it into smaller objects stops being free:
an Eddington-limited accretor has L proportional to m, so subdividing preserves
the total, but a Bondi-limited accretor has Mdot proportional to m^2, so the
total luminosity then falls linearly with the individual mass.
"""
import numpy as np

L_EDD_PER_MSUN = 1.26e38          # erg/s per solar mass
C_LIGHT = 2.998e10
M_SUN_G = 1.989e33
YR = 3.156e7
L_SUN = 3.828e33

M_BH = 7.0e7                      # solar masses
M_PER_PILLAR = 2.9e4              # solar masses, total, from the RHD run
N_PILLARS = 100
ETA = 0.1                         # radiative efficiency
L_AGN_BOL = 1.0e45                # erg/s, rough bolometric for PG 0844

M_TOT = M_PER_PILLAR * N_PILLARS
L_TOT = L_EDD_PER_MSUN * M_TOT

print(f"Total embedded mass  = {M_TOT:.3g} Msun "
      f"({M_TOT/M_BH*100:.2f}% of M_BH = {M_BH:.1g} Msun)")
print(f"Total Eddington L    = {L_TOT:.3g} erg/s")
print(f"  as a fraction of the black hole's own L_Edd "
      f"({L_EDD_PER_MSUN*M_BH:.3g}): {L_TOT/(L_EDD_PER_MSUN*M_BH)*100:.2f}%")
print(f"  as a fraction of the AGN bolometric ({L_AGN_BOL:.1g}): "
      f"{L_TOT/L_AGN_BOL*100:.1f}%")

mdot = L_TOT / (ETA * C_LIGHT**2)
mdot_agn = L_AGN_BOL / (ETA * C_LIGHT**2)
print(f"\nAccretion rate to sustain it at eta={ETA}: "
      f"{mdot:.3g} g/s = {mdot*YR/M_SUN_G:.4f} Msun/yr")
print(f"  the SMBH itself needs {mdot_agn*YR/M_SUN_G:.4f} Msun/yr, so the "
      f"embedded population would consume {mdot/mdot_agn*100:.0f}% as much")

print("\nSame total mass, different individual objects (all at Eddington):")
print(f"{'m_star':>10}{'N per pillar':>15}{'N total':>12}")
for m in [200.0, 50.0, 10.0, 1.0]:
    print(f"{m:>10.0f}{M_PER_PILLAR/m:>15,.0f}{M_TOT/m:>12,.0f}")

print("\nWhy 'Eddington' cannot mean ordinary stellar luminosity:")
print(f"{'m_star':>10}{'L_MS ~ m^3.5':>16}{'L_Edd':>12}{'L_MS/L_Edd':>13}"
      f"{'mass needed':>14}")
for m in [200.0, 50.0, 10.0, 1.0]:
    l_ms = L_SUN * m**3.5
    l_edd = L_EDD_PER_MSUN * m
    m_needed = L_TOT / l_ms * m            # total mass if each shines at L_MS
    print(f"{m:>10.0f}{l_ms:>16.3g}{l_edd:>12.3g}{l_ms/l_edd:>13.2e}"
          f"{m_needed:>14.3g}")
print("  (last column = total stellar mass required if the objects shine at "
      "their\n   main-sequence luminosity instead of Eddington)")

print("\nSubdivision scaling, fixed total mass M_tot = N m:")
print("  Eddington-limited:  L_tot = N L_Edd(m) ~ N m     = M_tot      (free)")
print("  Bondi-limited:      L_tot ~ N m^2               = M_tot * m  (falls)")
print("  Mdot_Bondi/Mdot_Edd ~ m, so the small objects are the ones that drop")
print("  below Eddington first. The critical mass needs a disc surface density,")
print("  which this model does not carry, so it cannot be evaluated here.")
