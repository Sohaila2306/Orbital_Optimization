import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Constants
R_EARTH = 6371e3  # Earth's radius in meters
MU = 3.986004418e14  # Earth's gravitational parameter (m^3/s^2)
G0 = 9.80665  # Standard gravity (m/s^2)

# Rocket parameters
ISP_STAGE1 = 300  # Specific impulse for stage 1 (seconds)
ISP_STAGE2 = 450  # Specific impulse for stage 2 (seconds)
FUEL_COST_PER_KG = 5000  # Cost per kg of fuel ($)
INITIAL_MASS = 1000  # Initial rocket mass (kg)


# === Delta-V Calculations ===

def hohmann_dv(r1, r2, di_rad):
    """
    Calculate total delta-v for a Hohmann transfer including plane change at first burn.

    Parameters:
        r1 (float): Initial orbit radius (m)
        r2 (float): Final orbit radius (m)
        di_rad (float): Inclination change in radians

    Returns:
        total_dv (float): Total delta-v required (m/s)
        dv1_total (float): First burn delta-v including plane change (m/s)
        dv2 (float): Second burn delta-v (m/s)
        dv_plane (float): Plane change delta-v component (m/s)
    """
    v1 = np.sqrt(MU / r1)  # Velocity in initial orbit
    v2 = np.sqrt(MU / r2)  # Velocity in final orbit

    # Delta-v for first burn (to enter transfer orbit)
    dv1 = v1 * (np.sqrt(2 * r2 / (r1 + r2)) - 1)
    # Delta-v for second burn (to circularize final orbit)
    dv2 = v2 * (1 - np.sqrt(2 * r1 / (r1 + r2)))

    # Plane change delta-v approximated at first burn velocity
    dv_plane = 2 * v1 * np.sin(di_rad / 2)

    # Combine first burn and plane change vectorially (Pythagorean)
    dv1_total = np.sqrt(dv1 ** 2 + dv_plane ** 2)

    # Sum total delta-v for the transfer
    total_dv = abs(dv1_total) + abs(dv2)
    return total_dv, abs(dv1_total), abs(dv2), dv_plane


def bielliptic_dv(r1, r2, rb, di_rad):
    """
    Calculate total delta-v for a bi-elliptic transfer including plane change at second burn.

    Parameters:
        r1 (float): Initial orbit radius (m)
        r2 (float): Final orbit radius (m)
        rb (float): Intermediate apoapsis radius (m)
        di_rad (float): Inclination change in radians

    Returns:
        total_dv (float): Total delta-v required (m/s)
        dv1 (float): First burn delta-v (m/s)
        dv2_total (float): Second burn delta-v including plane change (m/s)
        dv3 (float): Third burn delta-v (m/s)
        dv_plane (float): Plane change delta-v component at second burn (m/s)
    """
    v1 = np.sqrt(MU / r1)  # Initial orbit velocity

    # First burn: raise apoapsis to rb
    vA = np.sqrt(2 * MU * rb / (r1 * (r1 + rb)))
    dv1 = vA - v1

    # Second burn: change periapsis from rb to r2
    vB = np.sqrt(2 * MU * r1 / (rb * (r1 + rb)))
    vC = np.sqrt(2 * MU * r2 / (rb * (r2 + rb)))
    dv2 = vC - vB

    # Third burn: circularize at final orbit
    v2 = np.sqrt(MU / r2)
    vD = np.sqrt(2 * MU * rb / (r2 * (r2 + rb)))
    dv3 = v2 - vD

    # Plane change at second burn (at apoapsis orbit velocity)
    v_plane_burn = (vB + vC) / 2
    dv_plane = 2 * v_plane_burn * np.sin(di_rad / 2)

    # Combine second burn and plane change vectorially
    dv2_total = np.sqrt(dv2 ** 2 + dv_plane ** 2)

    # Sum total delta-v
    total_dv = abs(dv1) + abs(dv2_total) + abs(dv3)
    return total_dv, abs(dv1), abs(dv2_total), abs(dv3), dv_plane


# === Fuel and Cost Calculations ===

def compute_fuel_costs(dvs):
    """
    Calculate fuel consumption and cost based on delta-v sequence.

    Parameters:
        dvs (list of float): List of delta-v burns (m/s)

    Returns:
        fuel (list of float): Fuel used in each burn (kg)
        total_fuel (float): Total fuel consumed (kg)
        total_cost (float): Total cost of fuel ($)
    """
    m0 = INITIAL_MASS
    fuel = []
    for i, dv in enumerate(dvs):
        isp = ISP_STAGE1 if i == 0 else ISP_STAGE2
        mf = m0 / np.exp(dv / (isp * G0))
        fuel_used = m0 - mf
        fuel.append(fuel_used)
        m0 = mf
    total_fuel = sum(fuel)
    total_cost = total_fuel * FUEL_COST_PER_KG
    return fuel, total_fuel, total_cost


# === Plotting Functions ===

def plot_orbits(r1, r2, rb, inc1_deg, inc2_deg, method):
    """
    Plot Earth, initial and final orbits, and transfer orbit(s) in 3D.

    Parameters:
        r1 (float): Initial orbit radius (m)
        r2 (float): Final orbit radius (m)
        rb (float or None): Intermediate apoapsis radius for bi-elliptic (m)
        inc1_deg (float): Initial orbit inclination (degrees)
        inc2_deg (float): Final orbit inclination (degrees)
        method (str): "hohmann" or "bielliptic"
    """
    inc1_rad = np.radians(inc1_deg)
    inc2_rad = np.radians(inc2_deg)

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Draw Earth as a translucent sphere
    u, v = np.mgrid[0:2 * np.pi:100j, 0:np.pi:50j]
    x = R_EARTH * np.cos(u) * np.sin(v)
    y = R_EARTH * np.sin(u) * np.sin(v)
    z = R_EARTH * np.cos(v)
    ax.plot_surface(x, y, z, color='c', alpha=0.3)

    def plot_orbit(r, inc_rad, label, color):
        """Helper: plot a circular orbit with given radius and inclination."""
        theta = np.linspace(0, 2 * np.pi, 300)
        x_orbit = r * np.cos(theta)
        y_orbit = r * np.sin(theta)
        z_orbit = np.zeros_like(theta)

        # Rotate orbit plane about x-axis for inclination
        y_rot = y_orbit * np.cos(inc_rad) - z_orbit * np.sin(inc_rad)
        z_rot = y_orbit * np.sin(inc_rad) + z_orbit * np.cos(inc_rad)
        ax.plot(x_orbit, y_rot, z_rot, label=label, color=color, lw=2)

    # Plot initial and final circular orbits
    plot_orbit(r1, inc1_rad, "Initial Orbit", "blue")
    plot_orbit(r2, inc2_rad, "Final Orbit", "green")

    # Plot transfer orbits
    if method == "hohmann":
        a = (r1 + r2) / 2  # semi-major axis of transfer ellipse
        e = (r2 - r1) / (r2 + r1)  # eccentricity

        theta = np.linspace(0, np.pi, 300)  # half ellipse from r1 to r2
        r = a * (1 - e ** 2) / (1 + e * np.cos(theta))

        x_orbit = r * np.cos(theta)
        y_orbit = r * np.sin(theta)
        z_orbit = np.zeros_like(theta)

        y_rot = y_orbit * np.cos(inc1_rad) - z_orbit * np.sin(inc1_rad)
        z_rot = y_orbit * np.sin(inc1_rad) + z_orbit * np.cos(inc1_rad)

        ax.plot(x_orbit, y_rot, z_rot, label="Hohmann Transfer Orbit", color="orange", lw=2)

    else:  # Bi-elliptic transfer
        # First elliptical orbit r1 to rb
        a1 = (r1 + rb) / 2
        e1 = (rb - r1) / (rb + r1)
        theta1 = np.linspace(0, np.pi, 200)
        r_ellipse1 = a1 * (1 - e1 ** 2) / (1 + e1 * np.cos(theta1))

        x1 = r_ellipse1 * np.cos(theta1)
        y1 = r_ellipse1 * np.sin(theta1)
        y1_rot = y1 * np.cos(inc1_rad)
        z1_rot = y1 * np.sin(inc1_rad)
        ax.plot(x1, y1_rot, z1_rot, label="Bi-Elliptic Transfer Orbit 1", color="orange", lw=2)

        # Second elliptical orbit rb to r2
        a2 = (rb + r2) / 2
        e2 = (rb - r2) / (rb + r2)
        theta2 = np.linspace(np.pi, 2 * np.pi, 200)
        r_ellipse2 = a2 * (1 - e2 ** 2) / (1 + e2 * np.cos(theta2))

        x2 = r_ellipse2 * np.cos(theta2)
        y2 = r_ellipse2 * np.sin(theta2)
        y2_rot = y2 * np.cos(inc2_rad)
        z2_rot = y2 * np.sin(inc2_rad)
        ax.plot(x2, y2_rot, z2_rot, label="Bi-Elliptic Transfer Orbit 2", color="red", lw=2)

    # Set labels and title
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.set_zlabel('Z (m)')
    ax.set_title(f'Orbit Transfer: {method.capitalize()} Method')
    ax.legend()
    plt.show()


# === Main comparison function ===

def select_best_transfer(h1_km, h2_km, inc1_deg, inc2_deg):
    """
    Compare Hohmann and Bi-elliptic transfers for given orbit parameters,
    calculate fuel and cost, print results, and plot best transfer.
    """
    r1 = R_EARTH + h1_km * 1e3  # Convert altitude to radius in meters
    r2 = R_EARTH + h2_km * 1e3

    di_rad = np.abs(np.radians(inc2_deg) - np.radians(inc1_deg))  # Inclination change

    # Calculate Hohmann transfer delta-v and fuel
    dv_hohmann, dv1_hoh, dv2_hoh, dv_plane_hoh = hohmann_dv(r1, r2, di_rad)
    fuel_hoh, total_fuel_hoh, cost_hoh = compute_fuel_costs([dv1_hoh, dv2_hoh])

    # Calculate Bi-elliptic transfers for multiple rb values and pick minimum delta-v
    rb_candidates = [10 * r1, 15 * r1, 20 * r1]
    best_bielliptic = None

    for rb in rb_candidates:
        dv_bi, dv1_bi, dv2_bi, dv3_bi, dv_plane_bi = bielliptic_dv(r1, r2, rb, di_rad)
        fuel_bi, total_fuel_bi, cost_bi = compute_fuel_costs([dv1_bi, dv2_bi, dv3_bi])
        if best_bielliptic is None or dv_bi < best_bielliptic['dv']:
            best_bielliptic = {
                'dv': dv_bi,
                'rb': rb,
                'dv1': dv1_bi,
                'dv2': dv2_bi,
                'dv3': dv3_bi,
                'dv_plane': dv_plane_bi,
                'fuel': fuel_bi,
                'total_fuel': total_fuel_bi,
                'cost': cost_bi,
            }

    # Print results
    print("=== Transfer Orbit Comparison ===")
    print(f"Initial Orbit Altitude: {h1_km} km, Inclination: {inc1_deg}°")
    print(f"Final Orbit Altitude: {h2_km} km, Inclination: {inc2_deg}°")
    print("\n-- Hohmann Transfer --")
    print(f"Total Δv: {dv_hohmann:.2f} m/s")
    print(f"Fuel Used: {total_fuel_hoh:.2f} kg")
    print(f"Fuel Cost: ${cost_hoh:.2f}")

    print("\n-- Bi-Elliptic Transfer --")
    print(f"Intermediate Apoapsis Radius: {best_bielliptic['rb'] / 1e3:.0f} km")
    print(f"Total Δv: {best_bielliptic['dv']:.2f} m/s")
    print(f"Fuel Used: {best_bielliptic['total_fuel']:.2f} kg")
    print(f"Fuel Cost: ${best_bielliptic['cost']:.2f}")

    # Recommend best method based on total delta-v
    best_method = "Hohmann" if dv_hohmann < best_bielliptic['dv'] else "Bi-Elliptic"
    print(f"\nRecommended Transfer Method: {best_method}")

    # Plot best transfer orbit
    if best_method == "Hohmann":
        plot_orbits(r1, r2, None, inc1_deg, inc2_deg, "hohmann")
    else:
        plot_orbits(r1, r2, best_bielliptic['rb'], inc1_deg, inc2_deg, "bielliptic")


# === User Input ===

try:
    h1 = float(input("Enter initial orbit altitude (km): "))
    h2 = float(input("Enter final orbit altitude (km): "))
    inc1 = float(input("Enter initial orbit inclination (degrees): "))
    inc2 = float(input("Enter final orbit inclination (degrees): "))

    if h2 <= h1:
        print("Final orbit altitude must be greater than initial orbit altitude for transfer.")
    else:
        select_best_transfer(h1, h2, inc1, inc2)

except ValueError:
    print("Invalid input. Please enter numeric values.")
