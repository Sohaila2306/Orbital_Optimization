"""Command line interface.

    orbit-maneuvers compare --r1 6678 --r2 42164 --inc1 28.5 --mass 5000 --isp 320 --simulate
    orbit-maneuvers sweep-ratio --save-dir figs
    orbit-maneuvers sweep-inclination --r1 7000 --r2 7000 --max-angle 90
    orbit-maneuvers bodies

Distances are in km (or AU with --units au), angles in degrees.
"""

from __future__ import annotations

import argparse
import math
import os
import sys

from . import __version__
from .analysis import find_crossover_ratios, plane_change_sweep
from .bodies import BODIES, get_body
from .compare import compare_methods
from .constants import AU
from .plane_change import plane_angle
from .simulation import simulate_transfer


def _scale(units: str) -> float:
    return AU if units == "au" else 1e3


def _add_orbit_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--body", default="earth", help="central body (default: earth). See 'bodies'.")
    p.add_argument("--r1", type=float, required=True, help="start orbit radius (km, or AU with --units au)")
    p.add_argument("--r2", type=float, required=True, help="target orbit radius")
    p.add_argument("--altitude", action="store_true",
                   help="treat --r1/--r2/--rb as altitudes above the surface instead of radii from the centre")
    p.add_argument("--units", choices=["km", "au"], default="km", help="distance units for inputs (default: km)")


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="orbit-maneuvers",
        description="Compare Hohmann and bi-elliptic orbit transfers: delta-v, fuel and time.")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    # ---- compare
    c = sub.add_parser("compare", help="compare transfer methods for one mission")
    _add_orbit_args(c)
    c.add_argument("--inc1", type=float, default=0.0, help="start inclination, deg (default 0)")
    c.add_argument("--inc2", type=float, default=0.0, help="target inclination, deg (default 0)")
    c.add_argument("--draan", type=float, default=0.0, help="difference in RAAN between the planes, deg")
    c.add_argument("--plane-change", type=float, default=None,
                   help="total plane change in deg (overrides --inc1/--inc2/--draan)")
    c.add_argument("--rb", type=float, nargs="+", default=None,
                   help="bi-elliptic apoapsis radius(es). Omit to let the optimiser choose.")
    c.add_argument("--rb-max", type=float, default=None, help="upper limit for the rb optimiser")
    c.add_argument("--mode", default="optimal", choices=["optimal", "departure", "arrival", "apoapsis"],
                   help="where to do the plane change (default: optimal split)")
    c.add_argument("--mass", type=float, default=None, help="initial spacecraft mass, kg")
    c.add_argument("--isp", type=float, default=None, help="engine specific impulse, s")
    c.add_argument("--details", action="store_true", help="print a burn-by-burn breakdown")
    c.add_argument("--simulate", action="store_true", help="fly each transfer numerically and check the maths")
    c.add_argument("--plot", action="store_true", help="open plot windows")
    c.add_argument("--save-dir", default=None, help="save plots as PNGs into this folder")

    # ---- sweep-ratio
    s = sub.add_parser("sweep-ratio", help="classic Hohmann vs bi-elliptic chart + crossover ratios")
    s.add_argument("--max-ratio", type=float, default=40.0)
    s.add_argument("--plot", action="store_true")
    s.add_argument("--save-dir", default=None)

    # ---- sweep-inclination
    i = sub.add_parser("sweep-inclination", help="best delta-v of each method vs plane change angle")
    _add_orbit_args(i)
    i.add_argument("--max-angle", type=float, default=90.0)
    i.add_argument("--step", type=float, default=10.0)
    i.add_argument("--rb-max", type=float, default=None)
    i.add_argument("--plot", action="store_true")
    i.add_argument("--save-dir", default=None)

    sub.add_parser("bodies", help="list the built-in central bodies")
    return parser


def _to_radius(body, value: float, units: str, altitude: bool) -> float:
    metres = value * _scale(units)
    return body.altitude_to_radius(metres) if altitude else metres


def _wants_plots(args) -> bool:
    return bool(getattr(args, "plot", False) or getattr(args, "save_dir", None))


def _setup_matplotlib(args):
    """Pick a headless backend when we only need to save files."""
    import matplotlib
    if not getattr(args, "plot", False):
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def _save_path(args, name: str):
    d = getattr(args, "save_dir", None)
    if not d:
        return None
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, name)


def _cmd_compare(args) -> int:
    body = get_body(args.body)
    r1 = _to_radius(body, args.r1, args.units, args.altitude)
    r2 = _to_radius(body, args.r2, args.units, args.altitude)
    if args.plane_change is not None:
        dinc = math.radians(abs(args.plane_change))
    else:
        dinc = plane_angle(math.radians(args.inc1), math.radians(args.inc2), math.radians(args.draan))
    rb = None
    if args.rb:
        rb = [_to_radius(body, x, args.units, args.altitude) for x in args.rb]
    rb_max = None if args.rb_max is None else _to_radius(body, args.rb_max, args.units, args.altitude)

    result = compare_methods(body, r1, r2, dinc, rb=rb, plane_change_mode=args.mode, rb_max=rb_max,
                             initial_mass=args.mass, isp=args.isp)
    print(result.summary())

    if args.details:
        print()
        for t in result.transfers:
            print(t.burn_report())
            print()

    sims = None
    if args.simulate or _wants_plots(args):
        sims = [simulate_transfer(t, final_orbit_fraction=0.0 if _wants_plots(args) else 1.0)
                for t in result.transfers]
    if args.simulate:
        print()
        for sim in sims:
            print(sim.verification_report())
            print()

    if _wants_plots(args):
        plt = _setup_matplotlib(args)
        from . import visualization as viz
        viz.plot_trajectories(sims, save_path=_save_path(args, "trajectories.png"),
                              title=f"Transfer paths around {body.name}")
        viz.plot_delta_v_breakdown(result, save_path=_save_path(args, "delta_v_breakdown.png"))
        viz.plot_tradeoff(result, save_path=_save_path(args, "fuel_vs_time.png"))
        if dinc > 0 and not math.isclose(r1, r2, rel_tol=1e-9):
            viz.plot_plane_change_split(body, r1, r2, dinc, save_path=_save_path(args, "plane_change_split.png"))
        if args.save_dir:
            print(f"Plots saved to {args.save_dir}/")
        if args.plot:
            plt.show()
    return 0


def _cmd_sweep_ratio(args) -> int:
    lo, hi = find_crossover_ratios()
    print(f"Hohmann always wins below   r2/r1 = {lo:.2f}")
    print(f"Bi-elliptic always wins above r2/r1 = {hi:.2f}  (for any rb > r2)")
    print(f"Between the two it depends on how far out rb goes.")
    if _wants_plots(args):
        plt = _setup_matplotlib(args)
        from . import visualization as viz
        viz.plot_delta_v_vs_ratio(max_ratio=args.max_ratio, save_path=_save_path(args, "delta_v_vs_ratio.png"))
        if args.save_dir:
            print(f"Plot saved to {args.save_dir}/")
        if args.plot:
            plt.show()
    return 0


def _cmd_sweep_inclination(args) -> int:
    import numpy as np
    body = get_body(args.body)
    r1 = _to_radius(body, args.r1, args.units, args.altitude)
    r2 = _to_radius(body, args.r2, args.units, args.altitude)
    if math.isclose(r1, r2, rel_tol=1e-9):
        r2 = r1 * (1 + 1e-4)    # sweeps need two distinct radii; this is effectively a pure plane change
    rb_max = None if args.rb_max is None else _to_radius(body, args.rb_max, args.units, args.altitude)
    angles = np.arange(0.0, args.max_angle + 1e-9, args.step)
    sw = plane_change_sweep(body, r1, r2, angles, rb_max=rb_max)

    print(f"{'angle (deg)':>11}  {'Hohmann (km/s)':>15}  {'Bi-elliptic (km/s)':>19}  {'best rb (km)':>13}  winner")
    for a, h, b, rb in zip(sw["angle_deg"], sw["hohmann"], sw["bielliptic"], sw["rb"]):
        win = "bi-elliptic" if b < h - 1.0 else "Hohmann"
        print(f"{a:11.1f}  {h / 1e3:15.3f}  {b / 1e3:19.3f}  {rb / 1e3:13,.0f}  {win}")
    if _wants_plots(args):
        plt = _setup_matplotlib(args)
        from . import visualization as viz
        viz.plot_plane_change_sweep(sw, save_path=_save_path(args, "plane_change_sweep.png"))
        if args.save_dir:
            print(f"Plot saved to {args.save_dir}/")
        if args.plot:
            plt.show()
    return 0


def _cmd_bodies(_args) -> int:
    print(f"{'name':<9}{'mu (m^3/s^2)':>16}{'radius (km)':>14}{'SOI (km)':>16}")
    for key, b in BODIES.items():
        soi = f"{b.soi / 1e3:,.0f}" if b.soi else "-"
        print(f"{key:<9}{b.mu:>16.6e}{b.radius / 1e3:>14,.0f}{soi:>16}")
    return 0


def main(argv=None) -> int:
    args = _build_parser().parse_args(argv)
    handlers = {"compare": _cmd_compare, "sweep-ratio": _cmd_sweep_ratio,
                "sweep-inclination": _cmd_sweep_inclination, "bodies": _cmd_bodies}
    try:
        return handlers[args.command](args)
    except ValueError as err:
        print(f"error: {err}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
