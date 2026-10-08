"""Basic two-body relations. All inputs/outputs in SI units."""

from __future__ import annotations

import math


def circular_speed(mu: float, r: float) -> float:
    return math.sqrt(mu / r)


def vis_viva(mu: float, r: float, a: float) -> float:
    """Speed at radius r on an orbit with semi-major axis a."""
    return math.sqrt(mu * (2.0 / r - 1.0 / a))


def orbital_period(mu: float, a: float) -> float:
    return 2.0 * math.pi * math.sqrt(a**3 / mu)


def half_period(mu: float, a: float) -> float:
    """Time from periapsis to apoapsis."""
    return math.pi * math.sqrt(a**3 / mu)


def mean_motion(mu: float, a: float) -> float:
    """Angular rate of a circular orbit of radius a, rad/s."""
    return math.sqrt(mu / a**3)
