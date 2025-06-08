import numpy as np
import matplotlib.pyplot as plt

# Constants
EARTH_RADIUS = 6371e3        # Earth's radius in meters
MU = 3.986e14                # Earth's gravitational parameter (m^3/s^2)

# Starting orbit altitude (500 km above Earth's surface)
initial_orbit = EARTH_RADIUS + 500e3

# Apoapsis factors for bi-elliptic transfer (multiples of initial orbit radius)
apoapsis_factors = [10, 15, 20]

# Target altitudes for final orbit, from 1,000 km up to 1,000,000 km
target_altitudes_km = np.linspace(1000, 1_000_000, 300)
target_radii = EARTH_RADIUS + target_altitudes_km * 1e3

def hohmann_delta_v(r1, r2):
    """Calculate total delta-v for a Hohmann transfer between r1 and r2."""
    semi_major_axis = (r1 + r2) / 2
    delta_v1 = np.sqrt(MU/r1) * (np.sqrt(2*r2/(r1 + r2)) - 1)
    delta_v2 = np.sqrt(MU/r2) * (1 - np.sqrt(2*r1/(r1 + r2)))
    return abs(delta_v1) + abs(delta_v2)

def bielliptic_delta_v(r1, r2, rB):
    """Calculate total delta-v for a bi-elliptic transfer with apoapsis rB."""
    a1 = (r1 + rB) / 2
    a2 = (r2 + rB) / 2

    delta_v1 = np.sqrt(MU/r1) * (np.sqrt(2*rB/(r1 + rB)) - 1)
    delta_v2 = np.sqrt(MU/rB) * (np.sqrt(2*r2/(rB + r2)) - np.sqrt(2*r1/(rB + r1)))
    delta_v3 = np.sqrt(MU/r2) * (1 - np.sqrt(2*rB/(rB + r2)))
    return abs(delta_v1) + abs(delta_v2) + abs(delta_v3)

# Compute delta-v values for Hohmann transfer
hohmann_dv_values = [hohmann_delta_v(initial_orbit, r2) for r2 in target_radii]

# Compute delta-v values for Bi-elliptic transfer for each apoapsis factor
bielliptic_dv_values = {}
for factor in apoapsis_factors:
    rB = initial_orbit * factor
    bielliptic_dv_values[factor] = [bielliptic_delta_v(initial_orbit, r2, rB) for r2 in target_radii]

# Plot results
plt.figure(figsize=(10, 6))
plt.plot(target_altitudes_km / 1e3, np.array(hohmann_dv_values) / 1e3,
         label='Hohmann Transfer', linewidth=2)

for factor in apoapsis_factors:
    plt.plot(target_altitudes_km / 1e3, np.array(bielliptic_dv_values[factor]) / 1e3,
             linestyle='--', label=f'Bi-elliptic (rB = {factor} × r1)')

plt.xlabel("Target Altitude [×10³ km]")
plt.ylabel("Total Δv [km/s]")
plt.title("Delta-V Comparison: Hohmann vs. Bi-Elliptic Transfers")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
