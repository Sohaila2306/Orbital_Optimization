"""Numerical check of a transfer: actually fly it.

The analytic formulas in ``transfers.py`` assume perfect impulsive burns at
apsides. Here we integrate the two-body equations of motion with RK4 and apply
the burns as velocity changes to the *simulated* state. If the spacecraft ends
up on the intended circular orbit, in the intended plane, having spent the
predicted delta-v, the maths is doing what we think it is.

Geometry: the start orbit lies in the x-y plane, the spacecraft begins on the
+x axis (a node) moving toward +y. Plane changes rotate the velocity about the
x axis, so the target plane is the x-y plane tilted by the total plane change.
All burns happen on the x axis because apsides sit on the line of nodes.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .orbits import circular_speed, orbital_period


def _deriv(mu: float, y: np.ndarray) -> np.ndarray:
    r = y[:3]
    d = math.sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2])
    return np.concatenate((y[3:], -mu * r / d**3))


def _rk4(mu: float, y: np.ndarray, h: float) -> np.ndarray:
    k1 = _deriv(mu, y)
    k2 = _deriv(mu, y + 0.5 * h * k1)
    k3 = _deriv(mu, y + 0.5 * h * k2)
    k4 = _deriv(mu, y + h * k3)
    return y + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)


def propagate(mu: float, y0: np.ndarray, duration: float, eps: float = 0.02):
    """Integrate a coast arc. Returns (times, states), each including the start.

    The step size follows the local dynamical time sqrt(r^3/mu), so the
    integrator slows down near periapsis of eccentric arcs on its own.
    """
    y = np.array(y0, dtype=float)
    t = 0.0
    ts, ys = [0.0], [y.copy()]
    while t < duration - 1e-9:
        r = float(np.linalg.norm(y[:3]))
        h = min(eps * math.sqrt(r**3 / mu), duration - t)
        y = _rk4(mu, y, h)
        t += h
        ts.append(t)
        ys.append(y.copy())
    return np.array(ts), np.array(ys)


def _rotate_x(vec: np.ndarray, angle: float) -> np.ndarray:
    c, s = math.cos(angle), math.sin(angle)
    return np.array([vec[0], c * vec[1] - s * vec[2], s * vec[1] + c * vec[2]])


@dataclass
class SimulationResult:
    transfer: object                 # the Transfer that was flown
    times: np.ndarray                # s, shape (N,)
    states: np.ndarray               # [x, y, z, vx, vy, vz], shape (N, 6)
    burn_indices: list[int]          # index into times/states where each burn happens
    burn_vectors: list[np.ndarray]   # actual delta-v vectors applied (m/s)
    transfer_end: int                # index where the last burn happens (after it: final orbit)
    final_state: np.ndarray          # state right after the last burn

    @property
    def positions(self) -> np.ndarray:
        return self.states[:, :3]

    @property
    def simulated_delta_v(self) -> float:
        return float(sum(np.linalg.norm(v) for v in self.burn_vectors))

    def verification(self) -> dict:
        """How far did the simulation end up from what the formulas promised?"""
        t = self.transfer
        mu = t.body.mu
        r, v = self.final_state[:3], self.final_state[3:]
        h = np.cross(r, v)
        e_vec = np.cross(v, h) / mu - r / np.linalg.norm(r)
        tilt = math.acos(max(-1.0, min(1.0, h[2] / np.linalg.norm(h))))
        return {
            "delta_v_analytic": t.total_delta_v,
            "delta_v_simulated": self.simulated_delta_v,
            "delta_v_error": self.simulated_delta_v - t.total_delta_v,
            "final_radius_error_m": float(np.linalg.norm(r)) - t.r2,
            "final_eccentricity": float(np.linalg.norm(e_vec)),
            "plane_change_achieved_deg": math.degrees(tilt),
            "plane_change_error_deg": math.degrees(tilt - t.plane_change),
        }

    def verification_report(self) -> str:
        v = self.verification()
        return "\n".join([
            f"Simulation check for {self.transfer.label}",
            f"  delta-v: formula {v['delta_v_analytic']:.3f} m/s, simulated {v['delta_v_simulated']:.3f} m/s "
            f"(diff {v['delta_v_error']:+.3f} m/s)",
            f"  final radius off by {v['final_radius_error_m']:+.1f} m, "
            f"eccentricity {v['final_eccentricity']:.2e}",
            f"  plane change achieved {v['plane_change_achieved_deg']:.4f} deg "
            f"(error {v['plane_change_error_deg']:+.4f} deg)",
        ])


def simulate_transfer(transfer, eps: float = 0.02, final_orbit_fraction: float = 1.0) -> SimulationResult:
    """Fly a transfer numerically.

    ``eps`` controls the integrator step (smaller = more accurate and slower).
    ``final_orbit_fraction`` adds that fraction of an orbit in the target orbit
    after the last burn, which just makes plots look nicer.
    """
    mu = transfer.body.mu
    y = np.array([transfer.r1, 0.0, 0.0, 0.0, circular_speed(mu, transfer.r1), 0.0])

    times, states = [0.0], [y.copy()]
    burn_idx: list[int] = []
    burn_vecs: list[np.ndarray] = []
    clock = 0.0

    for i, burn in enumerate(transfer.burns):
        # --- apply the burn to the current simulated velocity
        v = y[3:]
        direction = _rotate_x(v / np.linalg.norm(v), burn.plane_change)
        v_new = burn.speed_after * direction
        burn_vecs.append(v_new - v)
        burn_idx.append(len(times) - 1)
        y = np.concatenate((y[:3], v_new))
        states[-1] = y.copy()        # the burn is instantaneous; overwrite the last sample

        # --- coast to the next burn (or around the final orbit)
        if i < len(transfer.coast_times):
            duration = transfer.coast_times[i]
        else:
            duration = final_orbit_fraction * orbital_period(mu, transfer.r2)
        if duration > 0:
            ts, ys = propagate(mu, y, duration, eps)
            times.extend(clock + ts[1:])
            states.extend(ys[1:])
            clock += duration
            y = ys[-1]

    states = np.array(states)
    final_state = states[burn_idx[-1]].copy()
    return SimulationResult(transfer, np.array(times), states, burn_idx, burn_vecs,
                            burn_idx[-1], final_state)
