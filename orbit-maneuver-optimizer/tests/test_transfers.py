import math

from orbit_maneuvers import (bielliptic, bielliptic_dv_norm, find_crossover_ratios, get_body,
                             hohmann, hohmann_dv_norm, optimal_bielliptic_rb)
from helpers import close

EARTH = get_body("earth")
R_LEO, R_GEO = 6678e3, 42164e3


def raises(exc, fn, *a, **kw):
    try:
        fn(*a, **kw)
    except exc:
        return True
    return False


def test_hohmann_leo_to_geo_matches_textbook():
    t = hohmann(EARTH, R_LEO, R_GEO)
    assert abs(t.burns[0].delta_v - 2426) < 5
    assert abs(t.burns[1].delta_v - 1467) < 5
    assert abs(t.total_delta_v / 1e3 - 3.893) < 0.005
    assert abs(t.total_time / 3600 - 5.275) < 0.01


def test_hohmann_inward_costs_the_same_as_outward():
    out = hohmann(EARTH, R_LEO, R_GEO)
    back = hohmann(EARTH, R_GEO, R_LEO)
    assert close(out.total_delta_v, back.total_delta_v)
    assert close(out.total_time, back.total_time)


def test_bielliptic_with_rb_equal_to_r2_degenerates_to_hohmann():
    h = hohmann(EARTH, R_LEO, R_GEO)
    b = bielliptic(EARTH, R_LEO, R_GEO, R_GEO * (1 + 1e-9))
    assert close(h.total_delta_v, b.total_delta_v, rel=1e-6)


def test_bielliptic_matches_closed_form():
    r1, R, B = 7000e3, 18.0, 40.0
    v1 = math.sqrt(EARTH.mu / r1)
    t = bielliptic(EARTH, r1, R * r1, B * r1)
    assert close(t.total_delta_v / v1, float(bielliptic_dv_norm(R, B)), rel=1e-9)
    assert close(hohmann(EARTH, r1, R * r1).total_delta_v / v1, float(hohmann_dv_norm(R)), rel=1e-9)


def test_crossover_ratios_are_the_famous_numbers():
    lo, hi = find_crossover_ratios()
    assert abs(lo - 11.94) < 0.01
    assert abs(hi - 15.58) < 0.01


def test_who_wins_on_each_side_of_the_crossover():
    r1 = 7000e3
    # below 11.94: no rb helps
    for rb_factor in (9, 20, 100):
        assert bielliptic(EARTH, r1, 8 * r1, rb_factor * r1).total_delta_v > hohmann(EARTH, r1, 8 * r1).total_delta_v
    # above 15.58: any rb above r2 helps
    for rb_factor in (1.01, 1.5, 3):
        r2 = 20 * r1
        assert bielliptic(EARTH, r1, r2, rb_factor * r2).total_delta_v < hohmann(EARTH, r1, r2).total_delta_v


def test_bielliptic_takes_longer():
    r1 = 7000e3
    h = hohmann(EARTH, r1, 20 * r1)
    b = bielliptic(EARTH, r1, 20 * r1, 60 * r1)
    assert b.total_time > h.total_time


def test_direct_plane_change_when_radii_match():
    r = 7000e3
    t = hohmann(EARTH, r, r, math.radians(30))
    v = math.sqrt(EARTH.mu / r)
    assert close(t.total_delta_v, 2 * v * math.sin(math.radians(15)))
    assert t.total_time == 0.0


def test_plane_change_costs_extra():
    flat = hohmann(EARTH, R_LEO, R_GEO)
    tilted = hohmann(EARTH, R_LEO, R_GEO, math.radians(28.5))
    assert tilted.total_delta_v > flat.total_delta_v
    assert close(tilted.plane_change_penalty, tilted.total_delta_v - flat.total_delta_v, rel=1e-9)


def test_leo_to_geo_with_28_5_degrees_optimal_split():
    t = hohmann(EARTH, R_LEO, R_GEO, math.radians(28.5))
    assert abs(t.total_delta_v / 1e3 - 4.231) < 0.005
    # classic result: only a couple of degrees at departure, the rest at GEO
    assert 1.5 < math.degrees(t.burns[0].plane_change) < 3.0


def test_optimal_split_beats_fixed_choices():
    ang = math.radians(28.5)
    opt = hohmann(EARTH, R_LEO, R_GEO, ang, "optimal").total_delta_v
    dep = hohmann(EARTH, R_LEO, R_GEO, ang, "departure").total_delta_v
    arr = hohmann(EARTH, R_LEO, R_GEO, ang, "arrival").total_delta_v
    assert opt <= dep and opt <= arr
    assert arr < dep      # doing it at the far, slow end is better than in LEO


def test_textbook_bielliptic_plane_change_threshold():
    """With the whole plane change at apoapsis, bi-elliptic only pays off past ~38.9 deg."""
    r = 7000e3
    v = math.sqrt(EARTH.mu / r)
    def best(deg):
        rbs = [r * 1.005 * 1.03**k for k in range(0, 160)]
        return min(bielliptic(EARTH, r, r * (1 + 1e-7), rb, math.radians(deg), "apoapsis").total_delta_v
                   for rb in rbs)
    direct = lambda deg: 2 * v * math.sin(math.radians(deg) / 2)
    assert best(35) > direct(35)
    assert best(43) < direct(43)


def test_optimal_rb_runs_to_the_limit_for_large_ratios():
    r1 = 7000e3
    rb = optimal_bielliptic_rb(EARTH, r1, 25 * r1, rb_max=2000e6)
    assert close(rb, 2000e6)


def test_between_the_crossovers_only_a_far_rb_wins():
    r1, R = 7000e3, 13.5          # between 11.94 and 15.58
    h = hohmann(EARTH, r1, R * r1).total_delta_v
    assert bielliptic(EARTH, r1, R * r1, 20 * r1).total_delta_v > h     # not far enough
    assert bielliptic(EARTH, r1, R * r1, 100 * r1).total_delta_v < h    # far enough


def test_optimal_rb_is_interior_when_plane_change_is_moderate():
    r = 7000e3
    rb = optimal_bielliptic_rb(EARTH, r, r * 1.0001, math.radians(30))
    assert 1.05 < rb / r < 5          # found a real optimum, not a boundary


def test_bad_inputs_raise():
    assert raises(ValueError, hohmann, EARTH, 100e3, R_GEO)            # below the surface
    assert raises(ValueError, hohmann, EARTH, R_LEO, R_LEO)            # nothing to do
    assert raises(ValueError, bielliptic, EARTH, R_LEO, R_GEO, 30000e3)  # rb too small
    assert raises(ValueError, hohmann, EARTH, R_LEO, R_GEO, 0.1, "apoapsis")
    assert raises(ValueError, get_body, "krypton")
