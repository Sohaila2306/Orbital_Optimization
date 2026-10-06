"""
Orbital_Calculator.py
Convert between orbital period and altitude for circular Earth orbits.
"""
import numpy as np
from Constants import R_earth, mu, MIN_ORBIT_ALTITUDE_KM, MAX_ORBIT_ALTITUDE_KM

GEO_ALTITUDE_KM = 35786     # geostationary altitude
GEO_TOLERANCE_KM = 100


def orbital_period(alt_km):
    """Orbital period [hours] of a circular orbit at the given altitude [km]."""
    radius_m = R_earth + alt_km * 1000
    period_sec = 2 * np.pi * np.sqrt(radius_m**3 / mu)
    return period_sec / 3600


def altitude_from_period(period_hours):
    """Altitude [km] of the circular orbit with the given period [hours]."""
    period_sec = period_hours * 3600
    radius_m = (mu * (period_sec / (2 * np.pi))**2) ** (1 / 3)
    return (radius_m - R_earth) / 1000


MIN_PERIOD_HOURS = orbital_period(MIN_ORBIT_ALTITUDE_KM)


def classify_orbit(alt_km):
    if alt_km < MIN_ORBIT_ALTITUDE_KM:
        return "Below stable orbit (invalid)"
    if alt_km > MAX_ORBIT_ALTITUDE_KM:
        return "Beyond Earth's gravitational influence (invalid)"
    if abs(alt_km - GEO_ALTITUDE_KM) < GEO_TOLERANCE_KM:      # checked BEFORE MEO/HEO
        return "GEO (Geostationary Orbit)"
    if alt_km < 2000:
        return "LEO (Low Earth Orbit)"
    if alt_km < GEO_ALTITUDE_KM:
        return "MEO (Medium Earth Orbit)"
    return "HEO (High Earth Orbit)"


def parse_numeric_input(user_input):
    try:
        return float(user_input.replace(',', '').strip())
    except ValueError:
        return None


def get_choice():
    while True:
        choice = input(
            "Do you want to calculate orbital period or altitude? "
            "Type 'a' for altitude (km) or 't' for time (hours): "
        ).strip().lower()
        if choice in ('a', 't'):
            return choice
        print("Oops, that's not valid. Please enter 'a' or 't'.")


def get_float_input(prompt, allow_zero=False, is_altitude=True):
    kind = "altitude" if is_altitude else "time"
    while True:
        val = parse_numeric_input(input(prompt))
        if val is None or not np.isfinite(val):
            print("That's not a valid number. Please try again.")
            continue
        if val < 0:
            print(f"Negative values aren't allowed for {kind}. Try again.")
            continue
        if val == 0 and not allow_zero:
            print(f"{kind.capitalize()} must be greater than zero.")
            continue
        return val


def main():
    choice = get_choice()
    if choice == 't':      # altitude -> period
        while True:
            alt_km = get_float_input("Enter the altitude in kilometers: ", allow_zero=True)
            if alt_km < MIN_ORBIT_ALTITUDE_KM:
                print(f"Altitude's too low for a stable orbit. Needs to be at least "
                      f"{MIN_ORBIT_ALTITUDE_KM} km.")
            elif alt_km > MAX_ORBIT_ALTITUDE_KM:
                print(f"Altitude's too high to stay gravitationally bound to Earth. "
                      f"Keep it under {MAX_ORBIT_ALTITUDE_KM:,} km.")
            else:
                print(f"At {alt_km:,.2f} km altitude, the orbital period is roughly "
                      f"{orbital_period(alt_km):.2f} hours.")
                print(f"Orbit type: {classify_orbit(alt_km)}")
                break
    else:                  # period -> altitude
        while True:
            period_hr = get_float_input("Enter the orbital period in hours: ",
                                        allow_zero=False, is_altitude=False)
            if period_hr < MIN_PERIOD_HOURS:
                print(f"That period's too short for a valid orbit "
                      f"(minimum about {MIN_PERIOD_HOURS:.2f} hours). Please try again.")
                continue
            altitude = altitude_from_period(period_hr)
            if altitude > MAX_ORBIT_ALTITUDE_KM:
                print(f"Calculated altitude {altitude:,.2f} km exceeds Earth's gravitational "
                      f"limit (~{MAX_ORBIT_ALTITUDE_KM:,} km).")
            else:
                print(f"An orbital period of {period_hr:.2f} hours corresponds to an "
                      f"altitude of about {altitude:,.2f} km.")
                print(f"Orbit type: {classify_orbit(altitude)}")
                break


if __name__ == "__main__":
    main()
