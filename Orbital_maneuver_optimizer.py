import Orbital_Calculator as oc
import hohmann_transfer as hoh
import Bi_Elptical_Transfer as bi
import Comparison as comp
import Non_Coplanar as nc

import matplotlib.pyplot as plt

# Earth's radius in km
R_earth_km = 6371

def main():
    try:
        # Get user input
        r1 = float(input("Enter initial orbit altitude (km): ")) + R_earth_km
        r2 = float(input("Enter final orbit altitude (km): ")) + R_earth_km
        inc1 = float(input("Enter initial orbit inclination (degrees): "))
        inc2 = float(input("Enter final orbit inclination (degrees): "))
    except ValueError:
        print("❌ Invalid input. Please enter numeric values.")
        return

    print(f"\nInitial Altitude: {r1 - R_earth_km:.1f} km, Final Altitude: {r2 - R_earth_km:.1f} km")
    print(f"Inclination Change: {abs(inc2 - inc1):.1f}°")

    # --- Planar transfer comparison ---
    print("\n=== Planar Transfer Comparison ===")
    comp.suggest_transfer(r1 * 1e3, r2 * 1e3)  # Convert km to meters

    # --- Non-Coplanar Transfer ---
    try:
        non_coplanar_result = nc.select_best_transfer(r1, r2, inc1_deg=inc1, inc2_deg=inc2)
        if not isinstance(non_coplanar_result, dict) or 'delta_v' not in non_coplanar_result:
            raise ValueError("select_best_transfer returned an invalid result.")

        print("\n=== Non-Coplanar Transfer ===")
        print(f"Total delta-v: {non_coplanar_result['delta_v'] / 1000:.3f} km/s")
        print(f"Fuel used: {non_coplanar_result['fuel']:.2f} kg")
        print(f"Cost estimate: ${non_coplanar_result['cost']:.2f}")
    except Exception as e:
        print(f"❌ Non-coplanar transfer computation failed: {e}")
        return

    # --- Delta-V Comparison Bar Chart ---
    try:
        hohmann_dv = comp.hohmann_delta_v(r1 * 1e3, r2 * 1e3)[0] / 1000  # km/s
        bi_factor = 10
        bi_dv = comp.bielliptic_delta_v(r1 * 1e3, r2 * 1e3, r1 * 1e3 * bi_factor)[0] / 1000
        noncoplanar_dv = non_coplanar_result['delta_v'] / 1000

        labels = ['Hohmann', 'Bi-Elliptic', 'Non-Coplanar']
        dv_values = [hohmann_dv, bi_dv, noncoplanar_dv]

        plt.figure(figsize=(8, 5))
        plt.bar(labels, dv_values, color=['blue', 'green', 'red'])
        plt.ylabel('Total Delta-V (km/s)')
        plt.title('Orbital Maneuver Delta-V Comparison')
        plt.grid(True)
        plt.tight_layout()
        plt.show()
    except Exception as e:
        print(f"⚠️ Failed to generate bar chart: {e}")

if __name__ == "__main__":
    main()
