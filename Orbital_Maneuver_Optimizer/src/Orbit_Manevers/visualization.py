"""Plots for mission planning.

Every function returns the matplotlib ``Figure`` and takes an optional
``save_path``. Nothing is shown automatically, so call ``plt.show()`` yourself
(or let the CLI do it with ``--plot``).
"""

from __future__ import annotations

import math
from typing import Sequence

import matplotlib.pyplot as plt
import numpy as np

from .analysis import (bi_parabolic_dv_norm, bielliptic_dv_norm, find_crossover_ratios,
                       hohmann_dv_norm, plane_change_split_curve)
from .bodies import get_body
from .simulation import SimulationResult, simulate_transfer
from .transfers import hohmann

PALETTE = ["#d1495b", "#2e86ab", "#edae49", "#66a182", "#8d6a9f", "#444444"]
_GRID = dict(alpha=0.25, linewidth=0.8)


def _finish(fig, save_path):
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=160, bbox_inches="tight")
    return fig


# ------------------------------------------------------------------ trajectories

def plot_trajectories(sims: Sequence[SimulationResult] | SimulationResult,
                      extent: float | None = None, title: str | None = None,
                      save_path: str | None = None):
    """Draw the simulated paths: start orbit, target orbit, and each transfer.

    Planar transfers get a top-down 2D view; if there's a plane change you get 3D.
    ``extent`` (in the body's display unit, e.g. km) crops the view, which is
    useful when a bi-elliptic apoapsis dwarfs everything else.
    """
    if isinstance(sims, SimulationResult):
        sims = [sims]
    ref = sims[0].transfer
    body, s = ref.body, ref.body.unit_scale
    three_d = ref.plane_change > 1e-9

    fig = plt.figure(figsize=(7.5, 7))
    ax = fig.add_subplot(projection="3d" if three_d else None)

    phi = np.linspace(0, 2 * np.pi, 400)
    tilt = ref.plane_change
    start = np.array([ref.r1 * np.cos(phi), ref.r1 * np.sin(phi), 0 * phi]) / s
    target = np.array([ref.r2 * np.cos(phi), ref.r2 * np.sin(phi) * math.cos(tilt),
                       ref.r2 * np.sin(phi) * math.sin(tilt)]) / s

    def line(arr, **kw):
        if three_d:
            ax.plot(arr[0], arr[1], arr[2], **kw)
        else:
            ax.plot(arr[0], arr[1], **kw)

    line(start, color="#888888", lw=1.4, ls="--", label="Start orbit")
    line(target, color="#222222", lw=1.4, ls="--", label="Target orbit")

    # the body itself
    if three_d:
        ax.scatter([0], [0], [0], color="#3b6ea5", s=60, label=body.name)
    else:
        ax.add_patch(plt.Circle((0, 0), body.radius / s, color="#3b6ea5", alpha=0.85, label=body.name))

    for i, sim in enumerate(sims):
        col = PALETTE[i % len(PALETTE)]
        end = sim.transfer_end + 1
        pos = sim.positions[:end] / s
        burn_pos = sim.positions[sim.burn_indices] / s
        if three_d:
            ax.plot(pos[:, 0], pos[:, 1], pos[:, 2], color=col, lw=2, label=sim.transfer.label)
            ax.scatter(burn_pos[:, 0], burn_pos[:, 1], burn_pos[:, 2], color=col, s=35, zorder=5)
        else:
            ax.plot(pos[:, 0], pos[:, 1], color=col, lw=2, label=sim.transfer.label)
            ax.scatter(burn_pos[:, 0], burn_pos[:, 1], color=col, s=35, zorder=5)

    if extent is None:
        extent = max(np.max(np.abs(sim.positions)) for sim in sims) / s * 1.08
    ax.set_xlim(-extent, extent)
    ax.set_ylim(-extent, extent)
    if three_d:
        ax.set_zlim(-extent, extent)
        ax.set_box_aspect((1, 1, 1))
        ax.set_zlabel(f"z ({body.unit_name})")
    else:
        ax.set_aspect("equal")
        ax.grid(**_GRID)
    ax.set_xlabel(f"x ({body.unit_name})")
    ax.set_ylabel(f"y ({body.unit_name})")
    ax.set_title(title or f"Transfer paths around {body.name}")
    ax.legend(loc="upper right", fontsize=9)
    return _finish(fig, save_path)


def plot_comparison_trajectories(comparison, extent: float | None = None,
                                 save_path: str | None = None):
    """Simulate every transfer in a Comparison and overlay the paths."""
    sims = [simulate_transfer(t, final_orbit_fraction=0.0) for t in comparison.transfers]
    title = None
    if comparison.plane_change:
        title = f"Transfers with a {math.degrees(comparison.plane_change):.1f} deg plane change"
    return plot_trajectories(sims, extent=extent, title=title, save_path=save_path)


# ------------------------------------------------------------------ delta-v and timing

def plot_delta_v_vs_ratio(rb_ratios: Sequence[float] = (15, 25, 50, 100), max_ratio: float = 40.0,
                          save_path: str | None = None):
    """The classic chart for planar transfers, in two panels.

    Top: total delta-v (in units of the starting circular speed v1) against the
    radius ratio r2/r1. Bottom: the same curves relative to Hohmann, so you can
    see where bi-elliptic starts to win (below zero) and by how much.
    The dotted lines mark the two crossover ratios.
    """
    R = np.linspace(2.0, max_ratio, 800)
    hoh = hohmann_dv_norm(R)
    lo, hi = find_crossover_ratios()

    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(8.5, 8), sharex=True,
                                  gridspec_kw={"height_ratios": [1.2, 1]})

    ax.plot(R, hoh, color="black", lw=2.4, label="Hohmann")
    ax2.axhline(0, color="black", lw=2.0)
    for i, B in enumerate(rb_ratios):
        mask = R < B
        col = PALETTE[i % len(PALETTE)]
        be = bielliptic_dv_norm(R[mask], B)
        ax.plot(R[mask], be, color=col, lw=1.8, label=f"Bi-elliptic, rb/r1 = {B:g}")
        ax2.plot(R[mask], 100 * (be / hoh[mask] - 1), color=col, lw=1.8)
    bp = bi_parabolic_dv_norm(R)
    ax.plot(R, bp, color="#777777", lw=1.6, ls="--", label="Bi-parabolic limit (rb -> inf)")
    ax2.plot(R, 100 * (bp / hoh - 1), color="#777777", lw=1.6, ls="--")

    for a in (ax, ax2):
        for x in (lo, hi):
            a.axvline(x, color="#999999", lw=1, ls=":")
        a.grid(**_GRID)
    ax.set_ylim(0.42, 0.62)
    ax2.set_ylim(-6, 6)
    for x in (lo, hi):
        ax2.text(x + 0.3, 5.3, f"{x:.2f}", color="#555555", fontsize=9)

    ax.set_ylabel("Total delta-v  /  v1")
    ax.set_title("Hohmann vs bi-elliptic (planar, circular to circular)")
    ax.legend(fontsize=9, loc="lower right")
    ax2.set_ylabel("Change vs Hohmann (%)\n(below 0 = bi-elliptic wins)")
    ax2.set_xlabel("Orbit radius ratio  r2 / r1")
    ax2.set_xlim(2, max_ratio)
    return _finish(fig, save_path)


def plot_plane_change_split(body, r1: float, r2: float, plane_change: float,
                            save_path: str | None = None):
    """Hohmann delta-v vs how the plane change is split between the two burns."""
    body = get_body(body)
    alpha, dv = plane_change_split_curve(body, r1, r2, plane_change)
    best = hohmann(body, r1, r2, plane_change, "optimal")
    a_best = best.burns[0].plane_change

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(np.degrees(alpha), dv / 1e3, color=PALETTE[1], lw=2.2)
    ax.scatter([0, math.degrees(plane_change)], [dv[0] / 1e3, dv[-1] / 1e3], color="#555555", zorder=5)
    ax.annotate("all at arrival", (0, dv[0] / 1e3), textcoords="offset points", xytext=(8, 8), fontsize=9)
    ax.annotate("all at departure", (math.degrees(plane_change), dv[-1] / 1e3), textcoords="offset points",
                xytext=(-85, 8), fontsize=9)
    ax.scatter([math.degrees(a_best)], [best.total_delta_v / 1e3], color=PALETTE[0], s=70, zorder=6)
    ax.annotate(f"optimum: {math.degrees(a_best):.2f} deg at departure,\n"
                f"{best.total_delta_v / 1e3:.3f} km/s",
                (math.degrees(a_best), best.total_delta_v / 1e3), textcoords="offset points",
                xytext=(20, 25), fontsize=9, arrowprops=dict(arrowstyle="-", color="#888888"))
    u, s = body.unit_name, body.unit_scale
    ax.set_xlabel("Plane change done in the departure burn (deg)")
    ax.set_ylabel("Total delta-v (km/s)")
    ax.set_title(f"Where to put the {math.degrees(plane_change):.1f} deg plane change  "
                 f"({r1 / s:,.0f} to {r2 / s:,.0f} {u})")
    ax.margins(y=0.1)
    ax.grid(**_GRID)
    return _finish(fig, save_path)


def plot_plane_change_sweep(sweep: dict, save_path: str | None = None):
    """Output of ``analysis.plane_change_sweep``: best delta-v of each method vs plane change."""
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(sweep["angle_deg"], sweep["hohmann"] / 1e3, color=PALETTE[0], lw=2.2,
            label="Hohmann (or direct plane change)")
    ax.plot(sweep["angle_deg"], sweep["bielliptic"] / 1e3, color=PALETTE[1], lw=2.2,
            label="Bi-elliptic (best rb)")
    ax.set_xlabel("Plane change (deg)")
    ax.set_ylabel("Total delta-v (km/s)")
    ax.set_title("Delta-v vs plane change")
    ax.grid(**_GRID)
    ax.legend()
    return _finish(fig, save_path)


def plot_delta_v_breakdown(comparison, save_path: str | None = None):
    """Stacked bars: how much of each method's delta-v is spent in each burn."""
    transfers = comparison.transfers
    fig, ax = plt.subplots(figsize=(7.5, 5))
    x = np.arange(len(transfers))
    bottom = np.zeros(len(transfers))
    for k, name in enumerate(["Departure burn", "Apoapsis burn", "Arrival burn", "Plane change"]):
        vals = np.array([sum(b.delta_v for b in t.burns if b.label == name) / 1e3 for t in transfers])
        if not vals.any():
            continue
        ax.bar(x, vals, bottom=bottom, color=PALETTE[(k + 1) % len(PALETTE)], width=0.55, label=name)
        bottom += vals
    for xi, tot in zip(x, bottom):
        ax.text(xi, tot + 0.03, f"{tot:.3f}", ha="center", fontsize=10, fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels([t.method for t in transfers])
    ax.set_ylabel("Delta-v (km/s)")
    ax.set_title("Delta-v by burn")
    ax.set_ylim(0, bottom.max() * 1.15)
    ax.grid(axis="y", **_GRID)
    ax.legend(fontsize=9)
    return _finish(fig, save_path)


def plot_tradeoff(comparison, save_path: str | None = None):
    """Delta-v vs trip time: the 'is saving fuel worth the wait?' view."""
    fig, ax = plt.subplots(figsize=(7.5, 5))
    for i, t in enumerate(comparison.transfers):
        ax.scatter(t.total_time / 3600, t.total_delta_v / 1e3, s=110, color=PALETTE[i % len(PALETTE)], zorder=5)
        ax.annotate(t.label, (t.total_time / 3600, t.total_delta_v / 1e3), textcoords="offset points",
                    xytext=(8, 8), fontsize=9)
    ax.set_xscale("log")
    ax.set_xlabel("Transfer time (hours, log scale)")
    ax.set_ylabel("Delta-v (km/s)")
    ax.set_title("Fuel vs time")
    ax.grid(**_GRID)
    return _finish(fig, save_path)
