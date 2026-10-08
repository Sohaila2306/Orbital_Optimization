"""Regenerate every figure used in the README (saved to docs/images/).

Run:  python examples/make_readme_figures.py
"""

import math
import pathlib

import matplotlib

matplotlib.use("Agg")

from orbit_maneuvers import compare_methods, plane_change_sweep, simulate_transfer
from orbit_maneuvers import visualization as viz

out = pathlib.Path(__file__).resolve().parent.parent / "docs" / "images"
out.mkdir(parents=True, exist_ok=True)

# 1. planar transfer well past the crossover ratio (r2/r1 ~ 18): the bi-elliptic wins
far = compare_methods("earth", 6678e3, 120_000e3, rb=400_000e3, initial_mass=2000, isp=450)
sims = [simulate_transfer(t, final_orbit_fraction=0.0) for t in far.transfers]
viz.plot_trajectories(sims, extent=430_000, title="LEO to 120,000 km: Hohmann vs bi-elliptic",
                      save_path=out / "trajectories_planar.png")
viz.plot_delta_v_breakdown(far, save_path=out / "delta_v_breakdown.png")
viz.plot_tradeoff(far, save_path=out / "fuel_vs_time.png")

# 2. LEO to GEO with the 28.5 degree plane change, in 3D
geo = compare_methods("earth", 6678e3, 42164e3, math.radians(28.5))
viz.plot_comparison_trajectories(geo, save_path=out / "trajectories_3d.png")
viz.plot_plane_change_split("earth", 6678e3, 42164e3, math.radians(28.5),
                            save_path=out / "plane_change_split.png")

# 3. the classic ratio chart
viz.plot_delta_v_vs_ratio(save_path=out / "delta_v_vs_ratio.png")

# 4. pure plane change sweep
sweep = plane_change_sweep("earth", 7000e3, 7000e3 * 1.0001, range(0, 91, 5))
viz.plot_plane_change_sweep(sweep, save_path=out / "plane_change_sweep.png")

print(f"Saved figures to {out}")
