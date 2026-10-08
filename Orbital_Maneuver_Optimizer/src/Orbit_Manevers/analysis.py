"""Dimensionless analysis: how do the methods compare as the orbit ratio changes?

For planar circular-to-circular transfers the answer depends only on
R = r2/r1 (and B = rb/r1 for bi-elliptic), so everything here is in units of
the initial circular speed v1. That's what the classic "Hohmann vs bi-elliptic"
chart plots, and it's what produces the famous 11.94 / 15.58 numbers.
"""

from __future__ import annotations

import math

import numpy as np

from .bodies import get_body
from .orbits import circular_speed, vis_viva
from .plane_change import burn_delta_v


# --------------------------------------------------------------- planar, normalised

def hohmann_dv_norm(R):
    """Hohmann delta-v / v1 for r2/r1 = R."""
    R = np.asarray(R, dtype=float)
    return (np.sqrt(2 * R / (1 + R)) - 1) + (1 / np.sqrt(R)) * (1 - np.sqrt(2 / (1 + R)))


def bielliptic_dv_norm(R, B):
    """Bi-elliptic delta-v / v1 for r2/r1 = R and rb/r1 = B (needs B > max(1, R))."""
    R = np.asarray(R, dtype=float)
    B = np.asarray(B, dtype=float)
    d1 = np.sqrt(2 * B / (1 + B)) - 1
    d2 = np.sqrt(2 * R / (B * (R + B))) - np.sqrt(2 / (B * (1 + B)))
    d3 = np.sqrt(2 * B / (R * (R + B))) - 1 / np.sqrt(R)
    return d1 + np.abs(d2) + np.abs(d3)


def bi_parabolic_dv_norm(R):
    """Limit of the bi-elliptic as rb -> infinity (the trip takes forever)."""
    R = np.asarray(R, dtype=float)
    return (math.sqrt(2) - 1) * (1 + 1 / np.sqrt(R))


def find_crossover_ratios() -> tuple[float, float]:
    """The two magic numbers for planar transfers.

    Returns ``(r_low, r_high)``:
      * below r_low  (~11.94): Hohmann always wins
      * r_low..r_high (~15.58): bi-elliptic wins only if rb is big enough
      * above r_high: *any* bi-elliptic with rb > r2 beats Hohmann
    """
    # r_low: where the infinitely-far bi-parabolic limit just ties Hohmann
    lo, hi = 5.0, 20.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if hohmann_dv_norm(mid) > bi_parabolic_dv_norm(mid):
            hi = mid
        else:
            lo = mid
    r_low = 0.5 * (lo + hi)

    # r_high: at rb = r2 the bi-elliptic *is* the Hohmann transfer. Nudge rb just
    # past r2 and see whether delta-v goes down (bi-elliptic helps) or up.
    def nudge_helps(R: float) -> bool:
        return float(bielliptic_dv_norm(R, R * (1 + 1e-5))) < float(hohmann_dv_norm(R))

    lo, hi = r_low, 40.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if nudge_helps(mid):
            hi = mid
        else:
            lo = mid
    return r_low, 0.5 * (lo + hi)


# --------------------------------------------------------------- plane change studies

def plane_change_split_curve(body, r1: float, r2: float, plane_change: float, n: int = 241):
    """Hohmann delta-v as a function of how much plane change is done at departure.

    Returns ``(alpha, dv)``: alpha is the share (rad) done in the first burn, the
    rest happens at arrival. Handy for seeing *why* the optimum sits where it does.
    """
    body = get_body(body)
    mu = body.mu
    a = 0.5 * (r1 + r2)
    vc1, vc2 = circular_speed(mu, r1), circular_speed(mu, r2)
    vt1, vt2 = vis_viva(mu, r1, a), vis_viva(mu, r2, a)
    alpha = np.linspace(0.0, plane_change, n)
    dv = burn_delta_v(vc1, vt1, alpha) + burn_delta_v(vt2, vc2, plane_change - alpha)
    return alpha, dv


def plane_change_sweep(body, r1: float, r2: float, angles_deg, rb_max: float | None = None):
    """Best Hohmann vs best bi-elliptic delta-v across a range of plane changes.

    Returns a dict of numpy arrays: ``angle_deg``, ``hohmann``, ``bielliptic``, ``rb``
    (delta-v in m/s, rb in m). For each angle the bi-elliptic uses the best rb
    (up to ``rb_max``) and the best split of the plane change.
    """
    from .transfers import bielliptic, hohmann, optimal_bielliptic_rb  # avoid import cycle at module load

    body = get_body(body)
    out = {"angle_deg": [], "hohmann": [], "bielliptic": [], "rb": []}
    for ang in angles_deg:
        rad = math.radians(ang)
        h = hohmann(body, r1, r2, rad)
        rb = optimal_bielliptic_rb(body, r1, r2, rad, rb_max=rb_max)
        b = bielliptic(body, r1, r2, rb, rad)
        out["angle_deg"].append(ang)
        out["hohmann"].append(h.total_delta_v)
        out["bielliptic"].append(b.total_delta_v)
        out["rb"].append(rb)
    return {k: np.array(v) for k, v in out.items()}
