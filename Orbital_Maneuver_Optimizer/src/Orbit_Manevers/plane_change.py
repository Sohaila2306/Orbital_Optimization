"""Plane-change geometry and the "where should I do it?" optimiser.

The workhorse is :func:`burn_delta_v`. If a burn changes speed from ``v_before``
to ``v_after`` *and* swings the velocity vector through ``angle``, the law of
cosines gives the cost. With ``angle = 0`` it collapses to ``|v_after - v_before|``
and with ``v_before == v_after`` it gives the classic ``2 v sin(angle / 2)``.
"""

from __future__ import annotations

import math
from typing import Sequence

import numpy as np


def plane_angle(inc1: float, inc2: float, delta_raan: float = 0.0) -> float:
    """Angle between two orbital planes (radians).

    Same RAAN -> just |inc2 - inc1|. If the nodes differ you need the full
    spherical-trig version, which is what this does.
    """
    c = (math.cos(inc1) * math.cos(inc2)
         + math.sin(inc1) * math.sin(inc2) * math.cos(delta_raan))
    return math.acos(max(-1.0, min(1.0, c)))


def burn_delta_v(v_before, v_after, angle):
    """Delta-v of a burn that changes speed and rotates velocity by ``angle``.

    Works on scalars or numpy arrays.
    """
    v_before = np.asarray(v_before, dtype=float)
    v_after = np.asarray(v_after, dtype=float)
    val = v_before**2 + v_after**2 - 2.0 * v_before * v_after * np.cos(angle)
    return np.sqrt(np.maximum(val, 0.0))


def simple_plane_change(v: float, angle: float) -> float:
    """Pure plane change on a circular orbit (or at any point with speed v)."""
    return 2.0 * v * math.sin(angle / 2.0)


def optimal_split(speed_pairs: Sequence[tuple[float, float]], total_angle: float,
                  grid: int = 101, rounds: int = 6) -> list[float]:
    """Spread ``total_angle`` over several burns so the total delta-v is smallest.

    ``speed_pairs`` is one ``(v_before, v_after)`` tuple per burn. Supports one,
    two or three burns, which covers Hohmann and bi-elliptic transfers.

    The cost isn't convex (for equal speeds it's actually concave, which is why
    "do it all in one burn" sometimes wins), so a plain gradient method can get
    stuck. Instead this does a coarse grid search and then zooms in a few times.
    """
    n = len(speed_pairs)
    if n not in (1, 2, 3):
        raise ValueError("optimal_split supports 1 to 3 burns")
    if total_angle <= 0:
        return [0.0] * n
    if n == 1:
        return [float(total_angle)]

    T = float(total_angle)

    def cost(i, ang):
        vb, va = speed_pairs[i]
        return burn_delta_v(vb, va, ang)

    if n == 2:
        lo, hi = 0.0, T
        best = 0.0
        for _ in range(rounds):
            a = np.linspace(lo, hi, grid)
            c = cost(0, a) + cost(1, T - a)
            k = int(np.argmin(c))
            best = float(a[k])
            step = (hi - lo) / (grid - 1)
            lo, hi = max(0.0, best - step), min(T, best + step)
        return [best, T - best]

    # three burns: search over (a1, a3), the middle one gets whatever is left
    lo1, hi1, lo3, hi3 = 0.0, T, 0.0, T
    b1 = b3 = 0.0
    for _ in range(rounds):
        a1 = np.linspace(lo1, hi1, grid)[:, None]
        a3 = np.linspace(lo3, hi3, grid)[None, :]
        a2 = T - a1 - a3
        c = cost(0, a1) + cost(2, a3) + np.where(
            a2 >= -1e-12, cost(1, np.maximum(a2, 0.0)), np.inf)
        i, j = np.unravel_index(int(np.argmin(c)), c.shape)
        b1, b3 = float(a1[i, 0]), float(a3[0, j])
        s1 = (hi1 - lo1) / (grid - 1)
        s3 = (hi3 - lo3) / (grid - 1)
        lo1, hi1 = max(0.0, b1 - s1), min(T, b1 + s1)
        lo3, hi3 = max(0.0, b3 - s3), min(T, b3 + s3)
    return [b1, max(T - b1 - b3, 0.0), b3]
