import numpy as np
import matplotlib.pyplot as plt

# === Constants ===
R_earth = 6371e3         # Earth's radius in meters
mu = 3.986004418e14      # Earth's gravitational parameter, m^3/s^2
g0 = 9.80665             # Standard gravity, m/s^2
ISP_FIRST_STAGE = 300
ISP_SECOND_STAGE = 450
FUEL_COST_PER_KG = 5000

# === Hohmann transfer function ===
def hohmann_transfer(initial_alt_km, target_alt_km):
    r1 = R_earth + initial_alt_km * 1e3
    r2 = R_earth + target_alt_km * 1e3
    a_transfer = (r1 + r2) / 2

    v1 = np.sqrt(mu / r1)
    v2 = np.sqrt(mu / r2)
    v_transfer1 = np.sqrt(mu * (2 / r1 - 1 / a_transfer))
    v_transfer2 = np.sqrt(mu * (2 / r2 - 1 / a_transfer))

    delta_v1 = v_transfer1 - v1
    delta_v2 = v2 - v_transfer2
    total_delta_v = abs(delta_v1) + abs(delta_v2)

    m0 = 1000  # Initial mass (kg)
    mf1 = m0 / np.exp(abs(delta_v1) / (ISP_FIRST_STAGE * g0))
    mf2 = mf1 / np.exp(abs(delta_v2) / (ISP_SECOND_STAGE * g0))
    fuel1 = m0 - mf1
    fuel2 = mf1 - mf2
    total_fuel = fuel1 + fuel2
    fuel_cost = total_fuel * FUEL_COST_PER_KG

    T_transfer = np.pi * np.sqrt(a_transfer ** 3 / mu)

    return {
        "r1": r1, "r2": r2, "a_transfer": a_transfer,
        "delta_v1": delta_v1, "delta_v2": delta_v2,
        "total_delta_v": total_delta_v,
        "fuel1": fuel1, "fuel2": fuel2,
        "total_fuel": total_fuel, "fuel_cost": fuel_cost,
        "transfer_time_sec": T_transfer,
        "transfer_time_hr": T_transfer / 3600
    }

# === Plotting function (static) ===
def plot_transfer(r1, r2, a_transfer):
    c = a_transfer - r1
    b = np.sqrt(a_transfer**2 - c**2)

    theta = np.linspace(0, np.pi, 300)
    x_ellipse = a_transfer * np.cos(theta) - c
    y_ellipse = b * np.sin(theta)

    theta_full = np.linspace(0, 2 * np.pi, 500)
    x_circle1 = r1 * np.cos(theta_full)
    y_circle1 = r1 * np.sin(theta_full)
    x_circle2 = r2 * np.cos(theta_full)
    y_circle2 = r2 * np.sin(theta_full)

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.set_aspect('equal')
    ax.set_xlim(-1.2 * max(r1, r2), 1.2 * max(r1, r2))
    ax.set_ylim(-1.2 * max(r1, r2), 1.2 * max(r1, r2))
    ax.set_title("Hohmann Transfer Orbit")

    # Earth: more realistic color
    ax.add_patch(plt.Circle((0, 0), R_earth, color='#1f77b4', label='Earth'))  # Blueish

    ax.plot(x_circle1, y_circle1, 'r--', label='Initial Orbit')
    ax.plot(x_circle2, y_circle2, 'g--', label='Target Orbit')
    ax.plot(x_ellipse, y_ellipse, color='orange', label='Transfer Orbit')

    ax.legend(loc='upper left')
    plt.show()

# === Main execution ===
if __name__ == "__main__":
    try:
        MIN_ORBIT_ALTITUDE_KM = 160
        MAX_ORBIT_ALTITUDE_KM = 1_500_000

        initial_alt = float(input("Enter initial orbit altitude in km: ").replace(",", ""))
        target_alt = float(input("Enter target orbit altitude in km: ").replace(",", ""))

        if initial_alt < MIN_ORBIT_ALTITUDE_KM or target_alt < MIN_ORBIT_ALTITUDE_KM:
            raise ValueError(f"Altitudes must be at least {MIN_ORBIT_ALTITUDE_KM} km for a stable orbit.")

        if initial_alt > MAX_ORBIT_ALTITUDE_KM or target_alt > MAX_ORBIT_ALTITUDE_KM:
            raise ValueError(
                f"Altitudes must be below {MAX_ORBIT_ALTITUDE_KM:,} km to stay within Earth's gravity influence.")

        result = hohmann_transfer(initial_alt, target_alt)

        print("\n--- Hohmann Transfer Results ---")
        print(f"Initial orbit radius: {result['r1']/1e3:.1f} km")
        print(f"Target orbit radius: {result['r2']/1e3:.1f} km")
        print(f"Transfer semi-major axis: {result['a_transfer']/1e3:.1f} km")
        print(f"Delta-v1: {result['delta_v1']:.2f} m/s")
        print(f"Delta-v2: {result['delta_v2']:.2f} m/s")
        print(f"Total Delta-v: {result['total_delta_v']:.2f} m/s")
        print(f"Fuel used (1st burn): {result['fuel1']:.2f} kg")
        print(f"Fuel used (2nd burn): {result['fuel2']:.2f} kg")
        print(f"Total fuel: {result['total_fuel']:.2f} kg")
        print(f"Estimated fuel cost: ${result['fuel_cost']:.2f}")
        print(f"Transfer time: {result['transfer_time_hr']:.2f} hours")

        plot_transfer(result['r1'], result['r2'], result['a_transfer'])

    except ValueError as e:
        print(f"Input error: {e}")
