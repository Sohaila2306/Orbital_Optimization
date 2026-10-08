"""Playing with the plane change: where should it happen, and when is it worth
swinging way out to do it?

Run:  python examples/04_plane_change_tradeoffs.py
"""

import math

from orbit_maneuvers import hohmann, plane_change_sweep, simple_plane_change

r_leo, r_geo = 6678e3, 42164e3
tilt = math.radians(28.5)

print("LEO -> GEO with a 28.5 deg plane change, Hohmann transfer")
for mode in ("departure", "arrival", "optimal"):
    t = hohmann("earth", r_leo, r_geo, tilt, mode)
    split = ", ".join(f"{math.degrees(b.plane_change):5.2f} deg" for b in t.burns)
    print(f"  {mode:<10} {t.total_delta_v / 1e3:.3f} km/s   (plane change per burn: {split})")

print("\nPure plane change in a 7000 km circular orbit (best bi-elliptic vs doing it directly):")
v = math.sqrt(3.986004418e14 / 7000e3)
sweep = plane_change_sweep("earth", 7000e3, 7000e3 * 1.0001, [10, 30, 50, 70, 90])
for ang, bi in zip(sweep["angle_deg"], sweep["bielliptic"]):
    direct = simple_plane_change(v, math.radians(ang))
    print(f"  {ang:>3.0f} deg   direct {direct / 1e3:6.3f} km/s   bi-elliptic {bi / 1e3:6.3f} km/s")
