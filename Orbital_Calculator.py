import numpy as np
from Constants import R_earth, mu

MIN_ORBIT_ALTITUDE_KM = 160       # Below this altitude, orbit isn't stable
MAX_ORBIT_ALTITUDE_KM = 1_500_000 # Rough upper limit of Earth's gravitational influence

def orbital_period(alt_km):
    radius_m = R_earth + alt_km * 1000
    period_sec = 2 * np.pi * np.sqrt(radius_m**3 / mu)
    return period_sec / 3600  # convert seconds to hours

def altitude_from_period(period_hours):
    period_sec = period_hours * 3600
    radius_m = (mu * (period_sec / (2 * np.pi))**2)**(1/3)
    alt_km = (radius_m - R_earth) / 1000
    return alt_km

def classify_orbit(alt_km):
    if alt_km < MIN_ORBIT_ALTITUDE_KM:
        return "Below stable orbit (invalid)"
    elif alt_km < 2000:
        return "LEO (Low Earth Orbit)"
    elif alt_km < 35786:
        return "MEO (Medium Earth Orbit)"
    elif abs(alt_km - 35786) < 100:
        return "GEO (Geostationary Orbit)"
    elif alt_km <= MAX_ORBIT_ALTITUDE_KM:
        return "HEO (High Earth Orbit)"
    else:
        return "Beyond Earth's gravitational influence (invalid)"

def parse_numeric_input(user_input):
    try:
        cleaned = user_input.replace(',', '').strip()
        return float(cleaned)
    except ValueError:
        return None

def get_choice():
    while True:
        choice = input(
            "Do you want to calculate orbital period or altitude? "
            "Type 'a' for altitude (km) or 't' for time (hours): "
        ).strip().lower()
        if choice in ['a', 't']:
            return choice
        print("Oops, that's not valid. Please enter 'a' or 't'.")

def get_float_input(prompt, allow_zero=False, is_altitude=True):
    while True:
        user_input = input(prompt)
        val = parse_numeric_input(user_input)
        if val is None:
            print("That's not a valid number. Please try again.")
            continue
        if val < 0:
            kind = "altitude" if is_altitude else "time"
            print(f"Negative values aren't allowed for {kind}. Try again.")
            continue
        if val == 0 and not allow_zero:
            kind = "altitude" if is_altitude else "time"
            print(f"{kind.capitalize()} must be greater than zero.")
            continue
        return val

def main():
    choice = get_choice()

    if choice == 't':
        while True:
            alt_km = get_float_input("Enter the altitude in kilometers: ", allow_zero=True, is_altitude=True)
            if alt_km < MIN_ORBIT_ALTITUDE_KM:
                print(f"Altitude's too low for a stable orbit. Needs to be at least {MIN_ORBIT_ALTITUDE_KM} km.")
            elif alt_km > MAX_ORBIT_ALTITUDE_KM:
                print(f"Altitude's too high to stay gravitationally bound to Earth. Keep it under {MAX_ORBIT_ALTITUDE_KM:,} km.")
            else:
                period = orbital_period(alt_km)
                classification = classify_orbit(alt_km)
                print(f"At {alt_km:,.2f} km altitude, the orbital period is roughly {period:.2f} hours.")
                print(f"Orbit type: {classification}")
                break

    else:  # choice == 'a'
        while True:
            period_hr = get_float_input("Enter the orbital period in hours: ", allow_zero=False, is_altitude=False)
            if period_hr < 1.4:
                print("That period's too short to be a valid orbit. Please try again.")
            else:
                altitude = altitude_from_period(period_hr)
                if altitude < MIN_ORBIT_ALTITUDE_KM:
                    print(f"Calculated altitude {altitude:,.2f} km is too low for stable orbit (minimum {MIN_ORBIT_ALTITUDE_KM} km).")
                elif altitude > MAX_ORBIT_ALTITUDE_KM:
                    print(f"Calculated altitude {altitude:,.2f} km exceeds Earth's gravitational field limit (~{MAX_ORBIT_ALTITUDE_KM:,} km).")
                else:
                    classification = classify_orbit(altitude)
                    print(f"An orbital period of {period_hr:.2f} hours corresponds to an altitude of about {altitude:,.2f} km.")
                    print(f"Orbit type: {classification}")
                    break

if __name__ == "__main__":
    main()
