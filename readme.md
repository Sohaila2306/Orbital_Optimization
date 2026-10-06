# Orbital Mechanics Calculator Suite

A collection of Python tools for calculating circular Earth orbit characteristics, parameters, and two-impulse Hohmann transfers, complete with delta-v, fuel, cost estimates, and trajectory visualizations.

---

## Project Structure

This project consists of the following modular Python scripts:

*   **`Constants.py`**: Central repository for shared physical constants (Earth mass, gravitational constant, standard gravity), orbital boundaries (minimum stable altitude, gravitational limits), rocket parameters (specific impulse for stages, initial mass), and user input validation utilities.
*   **`Orbital_Calculator.py`**: Interactive tool for converting between circular orbit altitude and orbital period. It also classifies orbits into standard regimes (LEO, MEO, GEO, HEO).
*   **`hohmann_transfer.py`**: Computes the optimal two-impulse Hohmann transfer between two coplanar circular orbits, calculating $\Delta v$, fuel consumption, mission cost, transfer time, and generating a Matplotlib trajectory plot.

---

## Prerequisites & Dependencies

To run these scripts, you will need **Python 3.x** along with the following third-party libraries:

*   `numpy` (for vector and mathematical operations)
*   `matplotlib` (for generating transfer orbit visualizations)

You can install the required packages via pip:
```bash
pip install numpy matplotlib
```

---

## How to Run

1. Make sure all three files (`Constants.py`, `Orbital_Calculator.py`, and `hohmann_transfer.py`) are placed in the same directory.
2. Run the orbital calculator tool:
   ```bash
   python Orbital_Calculator.py
   ```
3. Run the Hohmann transfer calculator and visualizer:
   ```bash
   python hohmann_transfer.py
   ```

---

## How to Upload to GitHub

Follow these steps to upload your folder to a new GitHub repository:

### Step 1: Create a New Repository on GitHub
1. Log in to your [GitHub account](https://github.com/).
2. In the top-right corner of the page, click the **`+`** icon and select **New repository**.
3. Give your repository a name (e.g., `orbital-mechanics-calculator`).
4. Choose whether your repository should be **Public** or **Private**.
5. **Important:** Do *not* check the box to initialize with a README, .gitignore, or license, since you are uploading your own existing local folder.
6. Click **Create repository**.

### Step 2: Initialize Git in Your Local Folder
Open your terminal (or command prompt), navigate to the folder containing your `.py` files and this `README.md`, and run the following commands:

1. Initialize a local Git repository:
   ```bash
   git init
   ```
2. Stage all your files (including `Constants.py`, `Orbital_Calculator.py`, `hohmann_transfer.py`, and `README.md`):
   ```bash
   git add .
   ```
3. Commit the files with a descriptive message:
   ```bash
   git commit -m "Initial commit: Add orbital calculator and Hohmann transfer scripts"
   ```

### Step 3: Link and Push to GitHub
1. Rename your default branch to `main` (if it isn't already):
   ```bash
   git branch -M main
   ```
2. Connect your local folder to your remote GitHub repository (replace `your-username` and `your-repo-name` with your actual GitHub username and repository name):
   ```bash
   git remote add origin https://github.com/your-username/your-repo-name.git
   ```
3. Push your code up to GitHub:
   ```bash
   git push -u origin main
   ```

Once completed, refresh your GitHub repository page, and your code and README will be live!