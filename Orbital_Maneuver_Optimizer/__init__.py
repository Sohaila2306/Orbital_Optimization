"""orbit_maneuvers: compare Hohmann and bi-elliptic transfers (planar and non-coplanar).

Quick taste::

    import math
    from orbit_maneuvers import compare_methods

    result = compare_methods("earth", 6678e3, 42164e3, plane_change=math.radians(28.5),
                             initial_mass=5000, isp=320)
    print(result.summary())

Plotting lives in ``orbit_maneuvers.visualization`` and is imported separately so that
the core stays usable without matplotlib doing anything at import time.
"""

from .analysis import (bi_parabolic_dv_norm, bielliptic_dv_norm, find_crossover_ratios,
                       hohmann_dv_norm, plane_change_split_curve, plane_change_sweep)
from .bodies import BODIES, CentralBody, get_body
from .compare import Comparison, compare_methods
from .fuel import FuelBudget, delta_v_from_masses, fuel_budget, mass_ratio, propellant_mass
from .plane_change import burn_delta_v, optimal_split, plane_angle, simple_plane_change
from .simulation import SimulationResult, simulate_transfer
from .timing import hohmann_phase_angle, synodic_period, wait_time
from .transfers import Burn, Transfer, bielliptic, hohmann, optimal_bielliptic_rb

__version__ = "0.1.0"

__all__ = [
    "BODIES", "Burn", "CentralBody", "Comparison", "FuelBudget", "SimulationResult", "Transfer",
    "bi_parabolic_dv_norm", "bielliptic", "bielliptic_dv_norm", "burn_delta_v", "compare_methods",
    "delta_v_from_masses", "find_crossover_ratios", "fuel_budget", "get_body", "hohmann",
    "hohmann_dv_norm", "hohmann_phase_angle", "mass_ratio", "optimal_bielliptic_rb",
    "optimal_split", "plane_angle", "plane_change_split_curve", "plane_change_sweep",
    "propellant_mass", "simple_plane_change", "simulate_transfer", "synodic_period", "wait_time",
]
