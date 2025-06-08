Orbit Optimization Project

This repository contains the implementation and analysis of various orbital transfer strategies for optimizing spacecraft trajectories. The project focuses on comparing Hohmann and bi-elliptic transfers, including coplanar and non-coplanar cases, to minimize fuel consumption (delta-v) and transfer cost.

---

Project Structure and Subprojects

The project is organized into the following subprojects/modules, each responsible for a specific part of the orbit optimization workflow:

1- Constants.py
This module contains fundamental physical and astronomical constants used throughout the project. These constants include gravitational parameters, planetary radii, standard gravitational acceleration, and other fixed values essential for accurate orbital calculations and simulations.

2-Orbital_Calculator.py:
A utility module that provides conversions between orbital period and altitude (or orbital radius) for circular orbits. It allows users to calculate the orbital altitude given the orbital period, and vice versa, assuming a two-body problem with a central gravitational body (e.g., Earth).
Features:

- Calculate orbital altitude from a given orbital period
- Calculate orbital period from a given orbital altitude
- Based on standard gravitational parameter of the central body (e.g., Earth)
- Useful for quick estimations in mission planning and orbit design  

3-hohmann_transfer.py:
Implementation of the classical two-impulse Hohmann transfer maneuver between two coplanar circular orbits. 
This module not only computes the transfer orbit parameters and delta-v but also estimates fuel consumption and associated mission cost. It provides detailed trajectory data and visualization tools to plot the orbital transfer and burn points.

Features:

- Calculate transfer orbit parameters
- Compute total delta-v required for the transfer
- Estimate fuel consumption based on delta-v and spacecraft propulsion characteristics
- Analyze mission cost implications related to fuel usage
- Generate plots of the transfer trajectory, including maneuver points for clear visualization  

4-Bi_Elptical_Transfer.py:
Implementation of the three-impulse bi-elliptic transfer maneuver for large orbit ratio changes.
This module calculates delta-v for each maneuver, estimates total fuel consumption and cost based on engine efficiency, and visualizes the orbital trajectory. It uses a fixed intermediate apoapsis (10× the target orbit by default).

Features:

- Compute delta-v for all three burns:

	- First burn (initial orbit to transfer ellipse)
	- Second burn (at apoapsis to second ellipse)
	- Third burn (circularization at target orbit)

- Estimate fuel mass required and calculate total mission fuel cost
- Compute total transfer time for the maneuver
- Plot initial, transfer, and final orbits with labeled maneuver points
- Interactive input: accepts initial and target orbit altitudes from the user 

5-Hohmann-BiElleptical Compare.py:
This script compares the total delta-v requirements of classical Hohmann and bi-elliptic transfer maneuvers from a 500 km circular Earth orbit to a range of higher circular orbits. It computes and plots the delta-v needed for each transfer method, using different apoapsis distances for the bi-elliptic case. The resulting plot helps identify which method is more efficient at varying target altitudes.

6-Comparison.py:
This Python script compares two orbital transfer techniques—Hohmann transfer and bi-elliptic transfer—to move a spacecraft from one circular orbit to a higher one. The user provides the initial and target orbit altitudes, and the program:

- Calculates the required delta-v for each method.
- Estimates fuel consumption using the Tsiolkovsky rocket equation, based on a two-stage rocket.
- Computes the cost of fuel using a predefined cost per kilogram.
- Considers multiple apoapsis distances for the bi-elliptic transfer to find the most fuel-efficient option.
- Visualizes the comparison in a set of bar charts showing total delta-v, fuel use, and cost.
- Provides a recommendation on which transfer method is more efficient.

This tool is useful for preliminary mission planning, trade studies, and educational purposes in astrodynamics and space mission design.

7-Non_coplanar.py:
This code compares two orbital transfer methods—Hohmann and bi-elliptic—between user-specified initial and final circular orbits with inclination changes. It calculates the total velocity changes (delta-v) required for each transfer, including plane change maneuvers, then estimates the rocket fuel consumption and cost based on these delta-v values and given rocket parameters. The code also visualizes the Earth, the initial and final orbits, and the transfer trajectory in 3D plots. Finally, it selects and recommends the more efficient transfer method based on minimum total delta-v, printing the results and plotting the best transfer orbit.

Features:

- Defines constants for Earth and rocket specs (gravity, specific impulses, fuel cost, mass).
- Calculates delta-v for both transfers, incorporating orbital mechanics formulas and plane changes.
- Computes fuel mass needed and fuel cost using rocket equation and specific impulses for two rocket stages.
- Plots Earth as a sphere and orbits as circles or ellipses in 3D, accounting for orbital inclinations.
- Takes user inputs for initial/final orbit altitudes and inclinations, validates them, and runs the comparison. 

8- Orbital_maneuver_optimizer.py:
- The program prompts the user to input parameters defining two orbits, including their altitudes and inclinations.

- It then analyzes different orbital transfer strategies to move between these orbits, considering both planar and non-coplanar scenarios.

- The results include key performance metrics like velocity change needed, fuel usage, and estimated costs for each transfer option.

- Finally, it presents a visual comparison to help understand the efficiency of the transfer methods.

- Throughout, the program manages input errors and computational issues to ensure smooth execution.

 

 



