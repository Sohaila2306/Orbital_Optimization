"""Hunting for the point where bi-elliptic starts to pay off.

For planar transfers the answer is purely a function of r2/r1: 11.94 and 15.58.
This script prints those, then tries a few real-ish missions either side.

Run:  python examples/03_when_does_bielliptic_win.py
"""

from orbit_maneuvers import compare_methods, find_crossover_ratios

lo, hi = find_crossover_ratios()
print(f"Crossover ratios: {lo:.2f} and {hi:.2f}\n")

r1 = 7000e3   # a 7000 km radius Earth orbit
for ratio in (6, 12, 14, 18, 30):
    result = compare_methods("earth", r1, ratio * r1, rb_max=400_000e3)
    h, b = result.transfers
    verdict = "bi-elliptic" if b.total_delta_v < h.total_delta_v else "Hohmann"
    print(f"r2/r1 = {ratio:>2}:  Hohmann {h.total_delta_v:7.1f} m/s   "
          f"bi-elliptic {b.total_delta_v:7.1f} m/s   -> {verdict}"
          f"   (trip time {b.total_time / h.total_time:.0f}x longer)")
