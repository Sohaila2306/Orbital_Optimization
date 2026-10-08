# Changelog

## 0.1.0

First public version.

- Hohmann and bi-elliptic transfers between circular orbits, planar or with a plane change
- Plane change placement: departure, arrival, apoapsis, or an optimised split across all burns
- Bi-elliptic apoapsis optimiser
- Rocket-equation propellant budgets
- Phase angle, synodic period and launch-window wait time
- RK4 two-body simulation that flies each transfer and checks the numbers
- Matplotlib plots (2D/3D trajectories, delta-v vs ratio, plane-change split, fuel vs time)
- `orbit-maneuvers` command line tool
