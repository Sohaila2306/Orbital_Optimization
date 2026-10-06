"""
hohmann_transfer.py
Two-impulse Hohmann transfer between two coplanar circular Earth orbits:
delta-v, fuel, cost, transfer time and a plot with the burn points.
"""
import numpy as np
import matplotlib.pyplot as plt

from Constants import (R_earth, mu, g0, ISP_STAGE1, ISP_STAGE2,
                       FUEL_COST_PER_KG, INITIAL_MASS, check_altitude_km, read_float)


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

    m0 = INITIAL_MASS
    mf1 = m0 / np.exp(abs(delta_v1) / (ISP_STAGE1 * g0))
    mf2 = mf1 / np.exp(abs(delta_v2) / (ISP_STAGE2 * g0))
    fuel1, fuel2 = m0 - mf1, mf1 - mf2
    total_fuel = fuel1 + fuel2

    T_transfer = np.pi * np.sqrt(a_transfer ** 3 / mu)

    return {
        "r1": r1, "r2": r2, "a_transfer": a_transfer,
        "delta_v1": delta_v1, "delta_v2": delta_v2, "total_delta_v": total_delta_v,
        "fuel1": fuel1, "fuel2": fuel2, "total_fuel": total_fuel,
        "fuel_cost": total_fuel * FUEL_COST_PER_KG,
        "transfer_time_sec": T_transfer, "transfer_time_hr": T_transfer / 3600,
    }


def plot_transfer(r1, r2, a_transfer):
    c = a_transfer - r1
    b = np.sqrt(r1 * r2)                      # semi-minor axis = sqrt(r_p * r_a)
    theta = np.linspace(0, np.pi, 300)
    x_ellipse = a_transfer * np.cos(theta) - c
    y_ellipse = b * np.sin(theta)
    t_full = np.linspace(0, 2 * np.pi, 500)

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.set_aspect('equal')
    lim = 1.2 * max(r1, r2)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_title("Hohmann Transfer Orbit")

    ax.add_patch(plt.Circle((0, 0), R_earth, color='#1f77b4', label='Earth'))
    ax.plot(r1 * np.cos(t_full), r1 * np.sin(t_full), 'r--', label='Initial Orbit')
    ax.plot(r2 * np.cos(t_full), r2 * np.sin(t_full), 'g--', label='Target Orbit')
    ax.plot(x_ellipse, y_ellipse, color='orange', label='Transfer Orbit')

    # maneuver points
    ax.plot([r1], [0], 'ko')
    ax.annotate(r"$\Delta v_1$", xy=(r1, 0), xytext=(r1 + 0.08 * lim, 0.08 * lim),
                arrowprops=dict(arrowstyle="->"))
    ax.plot([-r2], [0], 'ko')
    ax.annotate(r"$\Delta v_2$", xy=(-r2, 0), xytext=(-r2 - 0.08 * lim, 0.08 * lim),
                arrowprops=dict(arrowstyle="->"))

    ax.legend(loc='upper left')
    plt.show()


def main():
    try:
        initial_alt = check_altitude_km(read_float("Enter initial orbit altitude in km: "),
                                        "Initial altitude")
        target_alt = check_altitude_km(read_float("Enter target orbit altitude in km: "),
                                       "Target altitude")
    except ValueError as e:
        print(f"Input error: {e}")
        return

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


if __name__ == "__main__":
    main()
