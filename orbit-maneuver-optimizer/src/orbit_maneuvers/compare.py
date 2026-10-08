"""Run the methods side by side and get a readable answer."""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Sequence

from .bodies import CentralBody, get_body
from .fuel import FuelBudget, fuel_budget
from .transfers import (Transfer, bielliptic, default_rb_limit, hohmann,
                        optimal_bielliptic_rb)
from .utils import format_duration


@dataclass
class Comparison:
    body: CentralBody
    r1: float
    r2: float
    plane_change: float
    transfers: list[Transfer]
    budgets: list[FuelBudget] | None = None
    notes: list[str] = field(default_factory=list)

    def best_by_delta_v(self) -> Transfer:
        return min(self.transfers, key=lambda t: t.total_delta_v)

    def best_by_time(self) -> Transfer:
        return min(self.transfers, key=lambda t: t.total_time)

    def table(self) -> str:
        """Plain-text comparison table (ASCII only, so any terminal is happy)."""
        u, s = self.body.unit_name, self.body.unit_scale
        cheapest, fastest = self.best_by_delta_v(), self.best_by_time()

        headers = ["Method", f"rb ({u})", "Burns", "Delta-v (km/s)", "Duration"]
        if self.budgets:
            headers += ["Propellant (kg)", "Prop. %"]

        rows = []
        for i, t in enumerate(self.transfers):
            dv = f"{t.total_delta_v / 1e3:.3f}" + (" *" if t is cheapest else "")
            dur = format_duration(t.total_time) + (" +" if t is fastest else "")
            row = [t.method, self.body.format_distance(t.rb, with_unit=False) if t.rb else "-",
                   str(len(t.burns)), dv, dur]
            if self.budgets:
                b = self.budgets[i]
                row += [f"{b.propellant_total:,.1f}", f"{100 * b.propellant_fraction:.1f}"]
            rows.append(row)

        widths = [max(len(h), *(len(r[c]) for r in rows)) for c, h in enumerate(headers)]
        fmt = "  ".join(f"{{:<{w}}}" if c == 0 else f"{{:>{w}}}" for c, w in enumerate(widths))
        line = "  ".join("-" * w for w in widths)
        out = [fmt.format(*headers), line] + [fmt.format(*r) for r in rows]
        out.append("")
        out.append("* lowest delta-v    + shortest duration")
        return "\n".join(out)

    def summary(self) -> str:
        head = (f"{self.body.name}: r1 = {self.body.format_distance(self.r1)}  ->  "
                f"r2 = {self.body.format_distance(self.r2)}"
                f"   (ratio {self.r2 / self.r1:.2f}")
        if self.plane_change:
            head += f", plane change {math.degrees(self.plane_change):.2f} deg"
        head += ")"
        parts = [head, "", self.table()]
        if self.notes:
            parts += ["", "Notes:"] + [f"  - {n}" for n in self.notes]
        return "\n".join(parts)


def compare_methods(body, r1: float, r2: float, plane_change: float = 0.0,
                    rb: float | Sequence[float] | None = None,
                    plane_change_mode: str = "optimal",
                    rb_max: float | None = None,
                    initial_mass: float | None = None, isp: float | None = None) -> Comparison:
    """Compare Hohmann against bi-elliptic for one mission.

    Parameters
    ----------
    body : name ("earth", "mars", ...) or a CentralBody
    r1, r2 : initial / final circular orbit radii, metres
    plane_change : total inclination change between the orbits, radians
    rb : intermediate apoapsis radius for the bi-elliptic transfer (metres). Pass a
         number, a list of numbers to compare several designs, or leave it as None
         to let the optimiser pick the best one.
    plane_change_mode : "optimal" (default), "departure", "arrival", or
         "apoapsis" (bi-elliptic only). See ``transfers.py``.
    rb_max : upper limit on rb when optimising. Defaults to something sensible.
    initial_mass, isp : give both (kg, seconds) to also get propellant numbers.
    """
    body = get_body(body)
    transfers: list[Transfer] = []
    notes: list[str] = []

    # Hohmann can't use "apoapsis"; fall back to the nearest equivalent there
    h_mode = "arrival" if plane_change_mode == "apoapsis" else plane_change_mode
    hoh = hohmann(body, r1, r2, plane_change, h_mode)
    transfers.append(hoh)

    limit = rb_max if rb_max is not None else default_rb_limit(body, r1, r2)
    if rb is None:
        rb_list = [optimal_bielliptic_rb(body, r1, r2, plane_change, plane_change_mode, rb_max=limit)]
    elif isinstance(rb, (int, float)):
        rb_list = [float(rb)]
    else:
        rb_list = [float(x) for x in rb]

    bes = [bielliptic(body, r1, r2, x, plane_change, plane_change_mode) for x in rb_list]
    transfers.extend(bes)

    # a few honest remarks about what the numbers mean
    best_be = min(bes, key=lambda t: t.total_delta_v)
    saving = hoh.total_delta_v - best_be.total_delta_v
    ratio = max(r1, r2) / min(r1, r2)
    if saving > 1.0:
        if hoh.total_time > 0:
            when = f"but takes {best_be.total_time / hoh.total_time:,.1f}x as long as Hohmann"
        else:
            when = f"but needs {format_duration(best_be.total_time)} of coasting"
        notes.append(f"Bi-elliptic saves {saving:.1f} m/s ({100 * saving / hoh.total_delta_v:.1f}%) {when}.")
        if rb is None and math.isclose(best_be.rb, limit, rel_tol=1e-6):
            notes.append("The optimiser ran into the rb limit: a larger rb would save even more "
                         "(at the cost of a much longer trip). Raise rb_max if the dynamics allow it.")
    else:
        notes.append(f"Hohmann is cheaper here by {-saving:.1f} m/s.")
        if plane_change == 0 and ratio < 11.9:
            notes.append(f"Orbit ratio {ratio:.2f} is below ~11.94, so bi-elliptic can't beat "
                         "Hohmann for a planar transfer.")
        elif plane_change > 0:
            notes.append("Bi-elliptic tends to win for large plane changes or large orbit ratios; "
                         "this mission sits on the Hohmann side of that line.")

    budgets = None
    if initial_mass is not None and isp is not None:
        budgets = [fuel_budget(t, initial_mass, isp) for t in transfers]
    elif (initial_mass is None) != (isp is None):
        notes.append("Give both initial_mass and isp to get propellant numbers.")

    return Comparison(body, r1, r2, abs(plane_change), transfers, budgets, notes)
