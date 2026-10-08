import io
import math
from contextlib import redirect_stdout

from orbit_maneuvers import (compare_methods, delta_v_from_masses, fuel_budget, get_body, hohmann,
                             hohmann_phase_angle, mass_ratio, propellant_mass, synodic_period,
                             wait_time)
from orbit_maneuvers.cli import main
from helpers import close

AU = 1.495978707e11
MARS = 1.523679 * AU
MU_SUN = get_body("sun").mu


# ---- timing
def test_earth_mars_phase_angle_and_synodic_period():
    assert abs(math.degrees(hohmann_phase_angle(MU_SUN, AU, MARS)) - 44.3) < 0.1
    assert abs(synodic_period(MU_SUN, AU, MARS) / 86400 - 780) < 1


def test_inward_phase_angle_is_negative():
    assert hohmann_phase_angle(MU_SUN, MARS, AU) < 0


def test_wait_time_zero_at_the_window_and_bounded_by_synodic_period():
    req = hohmann_phase_angle(MU_SUN, AU, MARS)
    assert wait_time(MU_SUN, AU, MARS, req) < 1e-6
    for start in (0.0, 1.0, 2.5, -2.0):
        w = wait_time(MU_SUN, AU, MARS, start)
        assert 0 <= w <= synodic_period(MU_SUN, AU, MARS) + 1


def test_wait_time_actually_lands_on_the_window():
    n1, n2 = math.sqrt(MU_SUN / AU**3), math.sqrt(MU_SUN / MARS**3)
    start = 0.3
    w = wait_time(MU_SUN, AU, MARS, start)
    phase_then = start + (n2 - n1) * w
    req = hohmann_phase_angle(MU_SUN, AU, MARS)
    assert close(math.cos(phase_then), math.cos(req), rel=1e-9)
    assert close(math.sin(phase_then), math.sin(req), rel=1e-9, abs_=1e-9)


# ---- fuel
def test_rocket_equation_round_trip():
    m0, isp, dv = 5000.0, 320.0, 3893.0
    mf = m0 - propellant_mass(m0, dv, isp)
    assert close(delta_v_from_masses(m0, mf, isp), dv)
    assert close(m0 / mf, mass_ratio(dv, isp))


def test_fuel_budget_burn_by_burn_equals_lump_sum():
    t = hohmann("earth", 6678e3, 42164e3, math.radians(28.5))
    b = fuel_budget(t, 5000, 320)
    assert close(b.propellant_total, propellant_mass(5000, t.total_delta_v, 320), rel=1e-9)
    assert close(sum(b.per_burn), b.propellant_total, rel=1e-9)


# ---- compare
def test_compare_contains_both_methods_and_fuel():
    c = compare_methods("earth", 6678e3, 120000e3, initial_mass=2000, isp=450)
    assert [t.method for t in c.transfers] == ["Hohmann", "Bi-elliptic"]
    assert c.budgets is not None and len(c.budgets) == 2
    assert c.best_by_delta_v().method == "Bi-elliptic"   # ratio ~18 is past the crossover
    assert c.best_by_time().method == "Hohmann"
    assert "Hohmann" in c.summary()


def test_compare_accepts_several_rb_values():
    c = compare_methods("earth", 7000e3, 140000e3, rb=[300000e3, 400000e3])
    assert len(c.transfers) == 3


def test_compare_small_ratio_prefers_hohmann():
    c = compare_methods("earth", 6678e3, 42164e3)
    assert c.best_by_delta_v().method == "Hohmann"


# ---- cli
def run_cli(*argv):
    buf = io.StringIO()
    with redirect_stdout(buf):
        code = main(list(argv))
    return code, buf.getvalue()


def test_cli_compare_runs():
    code, out = run_cli("compare", "--r1", "300", "--r2", "35786", "--altitude", "--inc1", "28.5",
                        "--mass", "5000", "--isp", "320", "--simulate")
    assert code == 0
    assert "Hohmann" in out and "Simulation check" in out and "Propellant" in out


def test_cli_reports_bad_input_without_traceback():
    code, _ = run_cli("compare", "--r1", "100", "--r2", "200")
    assert code == 2


def test_cli_sweep_ratio_and_bodies():
    code, out = run_cli("sweep-ratio")
    assert code == 0 and "11.94" in out and "15.58" in out
    code, out = run_cli("bodies")
    assert code == 0 and "earth" in out
