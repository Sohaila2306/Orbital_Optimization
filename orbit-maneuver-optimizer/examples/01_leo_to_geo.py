"""LEO to GEO, the bread-and-butter mission.

A 5-tonne satellite starts in a 300 km parking orbit at 28.5 degrees (a Cape
Canaveral launch) and needs to reach the equatorial geostationary ring.

Run:  python examples/01_leo_to_geo.py
"""

import math

from orbit_maneuvers import compare_methods, get_body, simulate_transfer

earth = get_body("earth")
r_leo = earth.altitude_to_radius(300e3)
r_geo = 42_164e3

result = compare_methods(earth, r_leo, r_geo, plane_change=math.radians(28.5),
                         initial_mass=5000, isp=320)
print(result.summary())

print("\nBurn by burn:\n")
print(result.best_by_delta_v().burn_report())

# Don't take the formulas' word for it. Fly it.
print()
print(simulate_transfer(result.best_by_delta_v()).verification_report())
