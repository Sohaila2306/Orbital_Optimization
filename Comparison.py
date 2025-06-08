import numpy as np
import matplotlib.pyplot as plt

# === Constants ===
EARTH_RADIUS = 6371e3          # meters
MU = 3.986004418e14            # Earth's gravitational parameter (m^3/s^2)
G0 = 9.80665                  # Standard gravity (m/s^2)

# Fuel and rocket parameters
ISP_STAGE1 = 300               # Specific impulse stage 1 (s)
ISP_STAGE2 = 450               # Specific impulse stage 2 (s)
FUEL_COST_PER_KG = 5000        # Cost per kg of fuel ($)
INITIAL_MASS = 1000            # Initial mass of the spacecraft (kg)

# === Δv Calculations ===
def hohmann_delta_v(r1, r2):
    """Calculate total delta-v for a Hohmann transfer plus individual burns."""
    a = (r1 + r2) / 2
    dv1 = np.sqrt(MU / r1) * (np.sqrt(2 * r2 / (r1 + r2)) - 1)
    dv2 = np.sqrt(MU / r2) * (1 - np.sqrt(2 * r1 / (r1 + r2)))
    return abs(dv1) + abs(dv2), abs(dv1), abs(dv2)

def bielliptic_delta_v(r1, r2, rb):
    """Calculate total delta-v for a bi-elliptic transfer with apoapsis rb."""
    v1 = np.sqrt(MU / r1)
    vA = np.sqrt(2 * MU * rb / (r1 * (r1 + rb)))
    dv1 = vA - v1

    vB = np.sqrt(2 * MU * r1 / (rb * (r1 + rb)))
    vC = np.sqrt(2 * MU * r2 / (rb * (r2 + rb)))
    dv2 = vC - vB

    v2 = np.sqrt(MU / r2)
    vD = np.sqrt(2 * MU * rb / (r2 * (r2 + rb)))
    dv3 = v2 - vD

    return abs(dv1) + abs(dv2) + abs(dv3), abs(dv1), abs(dv2), abs(dv3)

# === Fuel cost calculation ===
def compute_fuel_costs(dv_list):
    """Calculate fuel consumption and cost given delta-v values per burn."""
    mass = INITIAL_MASS
    fuel_used_list = []
    for i, dv in enumerate(dv_list):
        isp = ISP_STAGE1 if i == 0 else ISP_STAGE2
        mf = mass / np.exp(dv / (isp * G0))   # final mass after burn
        fuel_used = mass - mf
        fuel_used_list.append(fuel_used)
        mass = mf   # update for next burn
    total_fuel = sum(fuel_used_list)
    total_cost = total_fuel * FUEL_COST_PER_KG
    return fuel_used_list, total_fuel, total_cost

# === Plotting comparison charts ===
def plot_comparison_subplots(data):
    labels = ['Hohmann', 'Bi-Elliptic']

    dv_values = [data['hohmann']['dv'] / 1e3, data['bi']['dv'] / 1e3]    # km/s
    fuel_values = [data['hohmann']['fuel'], data['bi']['fuel']]          # kg
    cost_values = [data['hohmann']['cost'], data['bi']['cost']]          # $

    fig, axs = plt.subplots(1, 3, figsize=(18, 5))

    # Δv subplot
    axs[0].barh(labels, dv_values, color=['steelblue', 'orange'])
    axs[0].set_title('Δv Comparison')
    axs[0].set_xlabel('Δv (km/s)')
    axs[0].grid(axis='x', linestyle='--', alpha=0.7)

    # Fuel subplot
    axs[1].barh(labels, fuel_values, color=['steelblue', 'orange'])
    axs[1].set_title('Fuel Consumption')
    axs[1].set_xlabel('Fuel (kg)')
    axs[1].grid(axis='x', linestyle='--', alpha=0.7)

    # Cost subplot
    axs[2].barh(labels, cost_values, color=['steelblue', 'orange'])
    axs[2].set_title('Cost Comparison')
    axs[2].set_xlabel('Cost ($)')
    axs[2].grid(axis='x', linestyle='--', alpha=0.7)

    plt.suptitle('Orbit Transfer Comparison by Metric', fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()

# === Transfer suggestion and summary ===
def suggest_transfer(r1, r2):
    hoh_total_dv, dv1, dv2 = hohmann_delta_v(r1, r2)
    fuel_h, total_fuel_h, cost_h = compute_fuel_costs([dv1, dv2])

    best_bi_dv = float('inf')
    best_factor = None
    best_fuel_b = None
    best_cost_b = None
    best_dvs = None

    # Check several apoapsis factors to find the most efficient bi-elliptic transfer
    for factor in [10, 15, 20]:
        rb = r1 * factor
        dv_total, dv1b, dv2b, dv3b = bielliptic_delta_v(r1, r2, rb)
        if dv_total < best_bi_dv:
            best_bi_dv = dv_total
            best_factor = factor
            best_dvs = [dv1b, dv2b, dv3b]
            fuel_b, total_fuel_b, cost_b = compute_fuel_costs(best_dvs)
            best_fuel_b = total_fuel_b
            best_cost_b = cost_b

    results = {
        "hohmann": {"dv": hoh_total_dv, "fuel": total_fuel_h, "cost": cost_h},
        "bi": {"dv": best_bi_dv, "fuel": best_fuel_b, "cost": best_cost_b},
    }

    # Print results and recommendation
    print(f"\nHohmann Δv: {hoh_total_dv:.2f} m/s, Fuel used: {total_fuel_h:.2f} kg, Estimated cost: ${cost_h:.2f}")
    print(f"Bi-Elliptic Δv (apoapsis = {best_factor}×r1): {best_bi_dv:.2f} m/s, Fuel used: {best_fuel_b:.2f} kg, Estimated cost: ${best_cost_b:.2f}")

    if hoh_total_dv < best_bi_dv:
        print("→ Suggestion: Hohmann transfer is more efficient ✅")
    elif best_bi_dv < hoh_total_dv:
        print("→ Suggestion: Bi-Elliptic transfer is more efficient ✅")
    else:
        print("→ Suggestion: Both methods are equally efficient ⚖️")

    plot_comparison_subplots(results)

# === Main program ===
if __name__ == "__main__":
    try:
        print("=== Orbit Transfer Analyzer ===")
        h1 = float(input("Enter initial orbit altitude (km): "))
        h2 = float(input("Enter target orbit altitude (km): "))

        if h2 <= h1:
            raise ValueError("Target altitude must be higher than the initial altitude.")

        r1 = EARTH_RADIUS + h1 * 1e3
        r2 = EARTH_RADIUS + h2 * 1e3

        suggest_transfer(r1, r2)

    except Exception as e:
        print(f"❌ Error: {e}")
