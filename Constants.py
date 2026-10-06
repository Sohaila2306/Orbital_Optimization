"""
Constants.py
Shared physical constants, rocket parameters and input-validation helpers
used by every module in the project (SI units unless stated otherwise).
"""

# ---------------- Earth ----------------
G = 6.67430e-11             # gravitational constant [m^3 kg^-1 s^-2]
M_earth = 5.972e24          # mass of Earth [kg]
R_earth = 6371e3            # mean radius of Earth [m]
mu = 3.986004418e14         # Earth's standard gravitational parameter [m^3/s^2]
                            # (the accepted value; G*M_earth is slightly less accurate)
g0 = 9.80665                # standard gravity [m/s^2]

# ---------------- Orbit limits ----------------
MIN_ORBIT_ALTITUDE_KM = 160         # below this, atmospheric drag makes an orbit unstable
MAX_ORBIT_ALTITUDE_KM = 1_500_000   # rough limit of Earth's gravitational influence

# ---------------- Rocket / cost model ----------------
ISP_STAGE1 = 300            # specific impulse, first burn [s]
ISP_STAGE2 = 450            # specific impulse, later burns [s]
FUEL_COST_PER_KG = 5000     # [$ / kg]
INITIAL_MASS = 1000         # spacecraft mass before the first burn [kg]

# ---------------- Bi-elliptic search ----------------
# Intermediate apoapsis candidates, as multiples of the LARGER of the two orbit radii
# (this guarantees the apoapsis is always above both orbits).
BIELLIPTIC_FACTORS = (1.5, 2, 3, 5, 10, 15, 20)


# ---------------- Validation helpers ----------------
def check_altitude_km(altitude_km, label="Altitude"):
    """Raise ValueError if the altitude is outside the supported range."""
    if altitude_km < MIN_ORBIT_ALTITUDE_KM:
        raise ValueError(
            f"{label} must be at least {MIN_ORBIT_ALTITUDE_KM} km for a stable orbit.")
    if altitude_km > MAX_ORBIT_ALTITUDE_KM:
        raise ValueError(
            f"{label} must be below {MAX_ORBIT_ALTITUDE_KM:,} km to stay within Earth's gravity.")
    return altitude_km


def check_inclination_deg(inclination_deg, label="Inclination"):
    """Raise ValueError if the inclination is outside 0-180 degrees."""
    if not 0 <= inclination_deg <= 180:
        raise ValueError(f"{label} must be between 0 and 180 degrees.")
    return inclination_deg


def read_float(prompt):
    """Read a number from the user; accepts thousands separators like 35,786."""
    return float(input(prompt).replace(",", "").strip())
