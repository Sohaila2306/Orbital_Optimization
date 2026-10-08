"""Phasing and launch-window helpers for circular, coplanar orbits.

These matter when you're transferring to a *moving* target (a station, a
planet): the transfer is only useful if the target arrives at the right spot
at the right time.
"""

from __future__ import annotations

import math

from .orbits import half_period, mean_motion


def hohmann_phase_angle(mu: float, r1: float, r2: float) -> float:
    """Where the target must be at departure, in radians, measured ahead of the
    spacecraft along its orbit. Negative means the target must be *behind*
    (typical for inward transfers).

    The target keeps moving during the coast, so it has to start
    ``pi - n2 * t_transfer`` ahead of us.
    """
    a = 0.5 * (r1 + r2)
    t = half_period(mu, a)
    n2 = mean_motion(mu, r2)
    phase = math.pi - n2 * t
    # wrap to (-pi, pi]
    return (phase + math.pi) % (2 * math.pi) - math.pi


def synodic_period(mu: float, r1: float, r2: float) -> float:
    """Time between repeats of the same relative geometry (i.e. between windows)."""
    n1, n2 = mean_motion(mu, r1), mean_motion(mu, r2)
    if math.isclose(n1, n2):
        return math.inf
    return 2 * math.pi / abs(n1 - n2)


def wait_time(mu: float, r1: float, r2: float, current_phase: float,
              required_phase: float | None = None) -> float:
    """How long until the departure window opens (seconds).

    ``current_phase`` is the target's present angle ahead of the spacecraft.
    ``required_phase`` defaults to the Hohmann phase angle.
    """
    if required_phase is None:
        required_phase = hohmann_phase_angle(mu, r1, r2)
    n1, n2 = mean_motion(mu, r1), mean_motion(mu, r2)
    rate = n2 - n1                    # how fast the phase angle changes
    if math.isclose(rate, 0.0, abs_tol=1e-15):
        raise ValueError("Same orbital rate: phase never changes, no window to wait for")
    two_pi = 2 * math.pi
    if rate < 0:
        return ((current_phase - required_phase) % two_pi) / -rate
    return ((required_phase - current_phase) % two_pi) / rate
