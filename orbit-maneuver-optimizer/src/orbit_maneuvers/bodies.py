"""A small catalogue of central bodies you can plan maneuvers around.

Numbers are the usual published values (gravitational parameter mu, mean
radius, and approximate sphere-of-influence radius). Feel free to add your own
by creating a ``CentralBody`` and passing it straight into any function.
"""

from __future__ import annotations

from dataclasses import dataclass

from .constants import AU


@dataclass(frozen=True)
class CentralBody:
    name: str
    mu: float                    # gravitational parameter, m^3/s^2
    radius: float                # mean equatorial radius, m
    soi: float | None = None     # sphere of influence radius, m (None = unbounded)
    unit_name: str = "km"        # distance unit used when printing / plotting
    unit_scale: float = 1e3      # metres per display unit
    unit_decimals: int = 0       # decimals to show when printing distances

    def format_distance(self, metres: float, with_unit: bool = True) -> str:
        """e.g. '42,164 km' or '1.52 AU', whatever suits this body."""
        text = f"{metres / self.unit_scale:,.{self.unit_decimals}f}"
        return f"{text} {self.unit_name}" if with_unit else text

    def altitude_to_radius(self, altitude: float) -> float:
        """Orbital radius (from the body's centre) for a given altitude, in metres."""
        return self.radius + altitude


BODIES: dict[str, CentralBody] = {
    "earth":   CentralBody("Earth",   3.986004418e14,  6.378137e6, 9.24e8),
    "moon":    CentralBody("Moon",    4.9048695e12,    1.7374e6,   6.61e7),
    "mars":    CentralBody("Mars",    4.282837e13,     3.3895e6,   5.77e8),
    "venus":   CentralBody("Venus",   3.24859e14,      6.0518e6,   6.16e8),
    "jupiter": CentralBody("Jupiter", 1.26686534e17,   7.1492e7,   4.82e10),
    # Heliocentric transfers (Earth -> Mars and friends). AU reads nicer than km here.
    "sun":     CentralBody("Sun",     1.32712440018e20, 6.957e8,   None, "AU", AU, 2),
}


def get_body(body: str | CentralBody) -> CentralBody:
    """Look a body up by name (case-insensitive), or pass a CentralBody straight through."""
    if isinstance(body, CentralBody):
        return body
    try:
        return BODIES[body.strip().lower()]
    except KeyError:
        known = ", ".join(sorted(BODIES))
        raise ValueError(f"Unknown body {body!r}. Available: {known}") from None
