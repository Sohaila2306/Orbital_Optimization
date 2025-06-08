import numpy as np
import matplotlib.pyplot as plt

# === Constants ===
R_earth = 6371e3         # Earth's radius in meters
mu = 3.986004418e14      # Earth's gravitational parameter, m^3/s^2
g0 = 9.80665             # Standard gravity, m/s^2
ISP_FIRST_STAGE = 300
ISP_SECOND_STAGE = 450
FUEL_COST_PER_KG = 5000

# === Bi-elliptic transfer function ===
def bielliptic_transfer(initial_alt_km, target_alt_km):
    r1 = R_earth + initial_alt_km * 1e3
    r3 = R_earth + target_alt_km * 1e3

    # Choose intermediate apoapsis as 10x target orbit by default
    r2 = 10 * r3

    # First burn (circular to transfer ellipse 1)
    v1 = np.sqrt(mu / r1)
    vA = np.sqrt(mu * (2 / r1 - 1 / ((r1 + r2) / 2)))
    delta_v1 = vA - v1

    # Second burn (apoapsis to second ellipse)
    vB = np.sqrt(mu * (2 / r2 - 1 / ((r1 + r2) / 2)))
    vC = np.sqrt(mu * (2 / r2 - 1 / ((r2 + r3) / 2)))
    delta_v2 = vC - vB

    # Third burn (insert into circular final orbit)
    vD = np.sqrt(mu * (2 / r3 - 1 / ((r2 + r3) / 2)))
    v2 = np.sqrt(mu / r3)
    delta_v3 = v2 - vD

    # Fuel estimation
    m0 = 1000
    mf1 = m0 / np.exp(abs(delta_v1) / (ISP_FIRST_STAGE * g0))
    mf2 = mf1 / np.exp(abs(delta_v2) / (ISP_SECOND_STAGE * g0))
    mf3 = mf2 / np.exp(abs(delta_v3) / (ISP_SECOND_STAGE * g0))
    fuel1 = m0 - mf1
    fuel2 = mf1 - mf2
    fuel3 = mf2 - mf3
    total_fuel = fuel1 + fuel2 + fuel3
    fuel_cost = total_fuel * FUEL_COST_PER_KG

    # Time: first + second half transfers
    T1 = np.pi * np.sqrt(((r1 + r2) / 2) ** 3 / mu)
    T2 = np.pi * np.sqrt(((r2 + r3) / 2) ** 3 / mu)
    total_time_sec = T1 + T2

    return {
        "r1": r1, "r2": r2, "r3": r3,
        "delta_v1": delta_v1, "delta_v2": delta_v2, "delta_v3": delta_v3,
        "fuel1": fuel1, "fuel2": fuel2, "fuel3": fuel3,
        "total_fuel": total_fuel, "fuel_cost": fuel_cost,
        "transfer_time_sec": total_time_sec, "transfer_time_hr": total_time_sec / 3600
    }

# === Plotting function ===
def plot_bielliptic(r1, r2, r3):
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_aspect('equal')
    lim = r2 * 1.2
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_title("Bi-elliptic Transfer Diagram")

    # Earth
    ax.add_patch(plt.Circle((0, 0), R_earth, color='blue', label='Earth'))

    # Orbits
    circle = lambda r: (r * np.cos(np.linspace(0, 2 * np.pi, 500)),
                        r * np.sin(np.linspace(0, 2 * np.pi, 500)))
    e1 = circle(r1)
    e3 = circle(r3)
    ax.plot(*e1, 'r--', label='Initial Orbit')
    ax.plot(*e3, 'g--', label='Final Orbit')

    # First transfer ellipse
    theta1 = np.linspace(0, np.pi, 300)
    a1 = (r1 + r2) / 2
    b1 = np.sqrt(r1 * r2)
    x1 = a1 * np.cos(theta1) - (a1 - r1)
    y1 = b1 * np.sin(theta1)
    ax.plot(x1, y1, 'orange', label='1st Transfer')

    # Second transfer ellipse
    theta2 = np.linspace(np.pi, 2 * np.pi, 300)
    a2 = (r2 + r3) / 2
    b2 = np.sqrt(r2 * r3)
    x2 = a2 * np.cos(theta2) - (a2 - r3)
    y2 = b2 * np.sin(theta2)
    ax.plot(x2, y2, 'purple', label='2nd Transfer')

    # Δv points
    ax.annotate(r"$\Delta v_1$", xy=(r1, 0), xytext=(r1 + 0.1*r1, 0.1*r1), arrowprops=dict(arrowstyle="->"))
    ax.annotate(r"$\Delta v_2$", xy=(0, r2), xytext=(0.1*r2, 1.1*r2), arrowprops=dict(arrowstyle="->"))
    ax.annotate(r"$\Delta v_3$", xy=(r3, 0), xytext=(r3 + 0.1*r3, -0.1*r3), arrowprops=dict(arrowstyle="->"))

    ax.legend(loc='upper left')
    plt.grid(True)
    plt.show()

# === Main Execution ===
if __name__ == '__main__':
    try:
        initial_alt = float(input("Enter initial orbit altitude (km): ").replace(",", ""))
        target_alt = float(input("Enter final orbit altitude (km): ").replace(",", ""))

        result = bielliptic_transfer(initial_alt, target_alt)

        print("\n--- Bi-elliptic Transfer Results ---")
        print(f"Initial radius: {result['r1']/1e3:.1f} km")
        print(f"Intermediate apoapsis: {result['r2']/1e3:.1f} km")
        print(f"Final radius: {result['r3']/1e3:.1f} km")
        print(f"Delta-v1: {result['delta_v1']:.2f} m/s")
        print(f"Delta-v2: {result['delta_v2']:.2f} m/s")
        print(f"Delta-v3: {result['delta_v3']:.2f} m/s")
        print(f"Total fuel: {result['total_fuel']:.2f} kg")
        print(f"Fuel cost: ${result['fuel_cost']:.2f}")
        print(f"Transfer time: {result['transfer_time_hr']:.2f} hours")

        plot_bielliptic(result['r1'], result['r2'], result['r3'])

    except ValueError as e:
        print(f"Input error: {e}")
