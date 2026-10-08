"""Transfer definitions: Hohmann, bi-elliptic, and plain plane changes.

Every transfer here goes between two *circular* orbits around the same body,
using impulsive burns at apsides. The orbits can differ in inclination; the
plane change is folded into the burns, and you choose where it happens
(``plane_change_mode``):

    "optimal"    split it across all burns in whatever way is cheapest
    "departure"  all of it in the first burn
    "arrival"    all of it in the last burn
    "apoapsis"   bi-elliptic only: all of it in the middle burn, where the
                 spacecraft is slowest (this is the textbook version)

Angles are radians, distances metres, speeds m/s, times seconds.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

import numpy as np

from .bodies import CentralBody, get_body
from .orbits import circular_speed, half_period, vis_viva
from .plane_change import burn_delta_v, optimal_split

_HOHMANN_MODES = ("optimal", "departure", "arrival")
_BIELLIPTIC_MODES = ("optimal", "departure", "apoapsis", "arrival")


@dataclass(frozen=True)
class Burn:
    label: str
    radius: float            # where it happens (distance from body centre)
    speed_before: float
    speed_after: float
    plane_change: float = 0.0   # how much of the plane change this burn does, rad

    @property
    def delta_v(self) -> float:
        return float(burn_delta_v(self.speed_before, self.speed_after, self.plane_change))

    @property
    def planar_delta_v(self) -> float:
        """What the burn would cost with no plane change at all."""
        return abs(self.speed_after - self.speed_before)


@dataclass(frozen=True)
class Transfer:
    method: str
    body: CentralBody
    r1: float
    r2: float
    plane_change: float                 # total, rad
    burns: tuple[Burn, ...]
    coast_times: tuple[float, ...]      # time between consecutive burns
    rb: float | None = None             # intermediate apoapsis (bi-elliptic only)
    plane_change_mode: str = "none"

    @property
    def total_delta_v(self) -> float:
        return sum(b.delta_v for b in self.burns)

    @property
    def total_time(self) -> float:
        return sum(self.coast_times)

    @property
    def plane_change_penalty(self) -> float:
        """Extra delta-v you pay because of the plane change."""
        return self.total_delta_v - sum(b.planar_delta_v for b in self.burns)

    @property
    def label(self) -> str:
        if self.rb is None:
            return self.method
        return f"{self.method} (rb = {self.body.format_distance(self.rb)})"

    def burn_report(self) -> str:
        """Burn-by-burn breakdown as text."""
        lines = [f"{self.label}: total {self.total_delta_v / 1e3:.3f} km/s over "
                 f"{self.total_time / 3600:.2f} h"]
        for b in self.burns:
            extra = f", plane change {math.degrees(b.plane_change):.2f} deg" if b.plane_change else ""
            lines.append(f"  - {b.label:<14} {b.delta_v / 1e3:7.3f} km/s  "
                         f"at r = {self.body.format_distance(b.radius)}{extra}")
        return "\n".join(lines)


# --------------------------------------------------------------------------- helpers

def _check_radii(body: CentralBody, *radii: float) -> None:
    for r in radii:
        if r <= 0:
            raise ValueError("Orbit radii must be positive")
        if r < body.radius:
            raise ValueError(
                f"Radius {r / 1e3:,.0f} km is below the surface of {body.name} "
                f"({body.radius / 1e3:,.0f} km). Did you give an altitude instead of a radius?")


def _angles(mode: str, valid: Sequence[str], pairs, total: float, apoapsis_index=None) -> list[float]:
    if mode not in valid:
        raise ValueError(f"plane_change_mode {mode!r} not valid here; choose from {valid}")
    n = len(pairs)
    if total == 0:
        return [0.0] * n
    if mode == "optimal":
        return optimal_split(pairs, total)
    idx = {"departure": 0, "arrival": n - 1, "apoapsis": apoapsis_index}[mode]
    out = [0.0] * n
    out[idx] = total
    return out


# --------------------------------------------------------------------------- transfers

def hohmann(body, r1: float, r2: float, plane_change: float = 0.0,
            plane_change_mode: str = "optimal") -> Transfer:
    """Two-burn transfer along an ellipse touching both circular orbits.

    Works inward or outward. If r1 == r2 this degrades to a single plane-change
    burn (and needs a non-zero ``plane_change``).
    """
    body = get_body(body)
    _check_radii(body, r1, r2)
    mu = body.mu
    plane_change = abs(plane_change)

    if math.isclose(r1, r2, rel_tol=1e-9):
        if plane_change == 0:
            raise ValueError("Orbits are identical: nothing to transfer")
        v = circular_speed(mu, r1)
        burn = Burn("Plane change", r1, v, v, plane_change)
        return Transfer("Direct plane change", body, r1, r2, plane_change, (burn,), (), None, "departure")

    a = 0.5 * (r1 + r2)
    vc1, vc2 = circular_speed(mu, r1), circular_speed(mu, r2)
    vt1, vt2 = vis_viva(mu, r1, a), vis_viva(mu, r2, a)
    pairs = [(vc1, vt1), (vt2, vc2)]
    ang = _angles(plane_change_mode, _HOHMANN_MODES, pairs, plane_change)

    burns = (Burn("Departure burn", r1, vc1, vt1, ang[0]),
             Burn("Arrival burn", r2, vt2, vc2, ang[1]))
    return Transfer("Hohmann", body, r1, r2, plane_change, burns,
                    (half_period(mu, a),), None, plane_change_mode)


def bielliptic(body, r1: float, r2: float, rb: float, plane_change: float = 0.0,
               plane_change_mode: str = "optimal") -> Transfer:
    """Three-burn transfer that swings out to an intermediate apoapsis ``rb``.

    ``rb`` has to be larger than both r1 and r2. The bigger it is, the lower the
    delta-v (once r2/r1 is big enough) and the longer the trip takes.
    """
    body = get_body(body)
    _check_radii(body, r1, r2, rb)
    if rb <= max(r1, r2):
        raise ValueError("rb must be larger than both r1 and r2")
    mu = body.mu
    plane_change = abs(plane_change)

    aA, aB = 0.5 * (r1 + rb), 0.5 * (r2 + rb)      # first and second ellipse
    vc1, vc2 = circular_speed(mu, r1), circular_speed(mu, r2)
    vA1, vA_b = vis_viva(mu, r1, aA), vis_viva(mu, rb, aA)
    vB_b, vB2 = vis_viva(mu, rb, aB), vis_viva(mu, r2, aB)
    pairs = [(vc1, vA1), (vA_b, vB_b), (vB2, vc2)]
    ang = _angles(plane_change_mode, _BIELLIPTIC_MODES, pairs, plane_change, apoapsis_index=1)

    burns = (Burn("Departure burn", r1, vc1, vA1, ang[0]),
             Burn("Apoapsis burn", rb, vA_b, vB_b, ang[1]),
             Burn("Arrival burn", r2, vB2, vc2, ang[2]))
    return Transfer("Bi-elliptic", body, r1, r2, plane_change, burns,
                    (half_period(mu, aA), half_period(mu, aB)), rb, plane_change_mode)


# --------------------------------------------------------------------------- rb search

def default_rb_limit(body, r1: float, r2: float) -> float:
    """Biggest intermediate apoapsis I'm willing to consider by default.

    50x the outer orbit, but never more than half the body's sphere of influence
    (past that, other bodies' gravity wrecks the two-body assumption).
    """
    body = get_body(body)
    limit = 50.0 * max(r1, r2)
    if body.soi is not None:
        limit = min(limit, 0.5 * body.soi)
    return limit


def optimal_bielliptic_rb(body, r1: float, r2: float, plane_change: float = 0.0,
                          plane_change_mode: str = "optimal", rb_max: float | None = None,
                          rb_min_factor: float = 1.05) -> float:
    """Find the intermediate apoapsis that minimises total delta-v.

    Scans a log-spaced grid, then polishes the best spot with a golden-section
    search. If delta-v keeps dropping all the way to ``rb_max`` (which happens
    whenever bi-elliptic is clearly the better option) you simply get ``rb_max`` back.
    """
    body = get_body(body)
    outer = max(r1, r2)
    lo = outer * rb_min_factor
    hi = rb_max if rb_max is not None else default_rb_limit(body, r1, r2)
    if hi <= lo:
        raise ValueError(
            "No room for a bi-elliptic transfer: the rb limit is too close to the outer orbit "
            "(is the orbit near the sphere of influence?)")

    def dv(rb):
        return bielliptic(body, r1, r2, rb, plane_change, plane_change_mode).total_delta_v

    grid = np.geomspace(lo, hi, 60)
    vals = [dv(x) for x in grid]
    k = int(np.argmin(vals))
    if k == len(grid) - 1:
        return float(hi)
    a, b = grid[max(k - 1, 0)], grid[k + 1]

    phi = (math.sqrt(5) - 1) / 2
    c, d = b - phi * (b - a), a + phi * (b - a)
    fc, fd = dv(c), dv(d)
    for _ in range(40):
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - phi * (b - a)
            fc = dv(c)
        else:
            a, c, fc = c, d, fd
            d = a + phi * (b - a)
            fd = dv(d)
    return float(0.5 * (a + b))
