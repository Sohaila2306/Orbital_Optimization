"""Earth to Mars around the Sun, with the timing side of things.

Delta-v is only half the story for an interplanetary trip: the planets have to
be in the right places. This prints the transfer, the phase angle you need at
departure and how often a launch window comes around.

Run:  python examples/02_earth_to_mars.py
"""

import math

from orbit_maneuvers import (get_body, hohmann, hohmann_phase_angle, synodic_period, wait_time)
from orbit_maneuvers.constants import AU

sun = get_body("sun")
r_earth, r_mars = 1.0 * AU, 1.523679 * AU

transfer = hohmann(sun, r_earth, r_mars)
print(transfer.burn_report())
print(f"Flight time: {transfer.total_time / 86400:.1f} days")

phase = math.degrees(hohmann_phase_angle(sun.mu, r_earth, r_mars))
print(f"Mars must be {phase:.1f} deg ahead of Earth at departure")
print(f"Launch windows repeat every {synodic_period(sun.mu, r_earth, r_mars) / 86400:.0f} days")

# Suppose Mars is currently 90 degrees ahead of Earth. How long until the window?
wait = wait_time(sun.mu, r_earth, r_mars, current_phase=math.radians(90))
print(f"If Mars is 90 deg ahead right now, the window opens in {wait / 86400:.0f} days")
