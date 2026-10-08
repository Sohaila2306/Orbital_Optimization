import math

from orbit_maneuvers import bielliptic, get_body, hohmann, simulate_transfer

EARTH = get_body("earth")


def check(sim, dv_tol=0.5):
    v = sim.verification()
    assert abs(v["delta_v_error"]) < dv_tol
    assert abs(v["final_radius_error_m"]) < 1e-6 * sim.transfer.r2
    assert v["final_eccentricity"] < 1e-4
    assert abs(v["plane_change_error_deg"]) < 1e-3


def test_planar_hohmann_flies_as_predicted():
    check(simulate_transfer(hohmann(EARTH, 6678e3, 42164e3)))


def test_inward_hohmann():
    check(simulate_transfer(hohmann(EARTH, 42164e3, 7000e3)))


def test_non_coplanar_hohmann():
    sim = simulate_transfer(hohmann(EARTH, 6678e3, 42164e3, math.radians(28.5)))
    check(sim)
    assert abs(sim.verification()["plane_change_achieved_deg"] - 28.5) < 1e-3


def test_bielliptic_with_plane_change():
    t = bielliptic(EARTH, 7000e3, 7000e3 * 6, 7000e3 * 25, math.radians(45))
    check(simulate_transfer(t))


def test_direct_plane_change():
    check(simulate_transfer(hohmann(EARTH, 7000e3, 7000e3, math.radians(40))))


def test_heliocentric_earth_to_mars():
    sun = get_body("sun")
    au = 1.495978707e11
    sim = simulate_transfer(hohmann(sun, au, 1.523679 * au), final_orbit_fraction=0.0)
    check(sim, dv_tol=1.0)
    assert abs(sim.times[sim.transfer_end] / 86400 - 258.9) < 0.5


def test_burn_indices_line_up_with_burn_count():
    t = bielliptic(EARTH, 7000e3, 42164e3, 100000e3)
    sim = simulate_transfer(t, final_orbit_fraction=0.0)
    assert len(sim.burn_indices) == len(t.burns) == 3
    assert sim.burn_indices[0] == 0
