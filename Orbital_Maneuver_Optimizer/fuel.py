"""Propellant bookkeeping with the Tsiolkovsky rocket equation."""

from __future__ import annotations

import math
from dataclasses import dataclass

from .constants import G0


def mass_ratio(delta_v: float, isp: float) -> float:
    """m_initial / m_final for a given delta-v (m/s) and specific impulse (s)."""
    return math.exp(delta_v / (isp * G0))


def propellant_mass(initial_mass: float, delta_v: float, isp: float) -> float:
    return initial_mass * (1.0 - math.exp(-delta_v / (isp * G0)))


def delta_v_from_masses(initial_mass: float, final_mass: float, isp: float) -> float:
    """Inverse of the above: how much delta-v does this much propellant buy?"""
    return isp * G0 * math.log(initial_mass / final_mass)


@dataclass(frozen=True)
class FuelBudget:
    initial_mass: float            # kg, before the first burn
    isp: float                     # s
    per_burn: tuple[float, ...]    # kg of propellant burned at each burn
    final_mass: float              # kg, after the last burn

    @property
    def propellant_total(self) -> float:
        return self.initial_mass - self.final_mass

    @property
    def propellant_fraction(self) -> float:
        return self.propellant_total / self.initial_mass


def fuel_budget(transfer, initial_mass: float, isp: float) -> FuelBudget:
    """Burn through the transfer's burns one at a time, tracking the shrinking mass.

    (Doing it burn by burn instead of in one lump gives the same total, but it
    also tells you how much each burn eats, which is handy when sizing tanks.)
    """
    if initial_mass <= 0 or isp <= 0:
        raise ValueError("initial_mass and isp must be positive")
    mass = initial_mass
    used = []
    for burn in transfer.burns:
        p = propellant_mass(mass, burn.delta_v, isp)
        used.append(p)
        mass -= p
    return FuelBudget(initial_mass, isp, tuple(used), mass)
