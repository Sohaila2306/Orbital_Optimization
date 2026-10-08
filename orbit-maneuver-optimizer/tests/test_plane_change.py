import math

from orbit_maneuvers import burn_delta_v, optimal_split, plane_angle, simple_plane_change
from helpers import close


def test_burn_without_rotation_is_speed_difference():
    assert close(float(burn_delta_v(7000, 9000, 0.0)), 2000.0)


def test_burn_with_equal_speeds_is_simple_plane_change():
    v, ang = 7500.0, math.radians(20)
    assert close(float(burn_delta_v(v, v, ang)), simple_plane_change(v, ang))


def test_plane_angle_same_node_is_inclination_difference():
    assert close(plane_angle(math.radians(28.5), math.radians(51.6)), math.radians(51.6 - 28.5))


def test_plane_angle_with_node_offset():
    # two 90-degree-inclined planes whose nodes differ by 90 degrees are 90 degrees apart
    assert close(plane_angle(math.pi / 2, math.pi / 2, math.pi / 2), math.pi / 2)
    # equatorial and polar planes are 90 degrees apart regardless of node
    assert close(plane_angle(0.0, math.pi / 2, 1.2), math.pi / 2)


def test_optimal_split_sums_to_total_and_is_no_worse_than_extremes():
    pairs = [(7700.0, 10200.0), (1600.0, 3070.0)]
    total = math.radians(28.5)
    split = optimal_split(pairs, total)
    assert close(sum(split), total)
    cost = lambda a: sum(float(burn_delta_v(vb, va, th)) for (vb, va), th in zip(pairs, a))
    assert cost(split) <= cost([total, 0.0]) + 1e-9
    assert cost(split) <= cost([0.0, total]) + 1e-9


def test_optimal_split_three_burns_sums_to_total():
    pairs = [(7500.0, 7900.0), (6600.0, 6900.0), (7000.0, 7500.0)]
    total = math.radians(50)
    split = optimal_split(pairs, total)
    assert close(sum(split), total, rel=1e-9)
    assert all(a >= 0 for a in split)


def test_zero_plane_change_gives_zero_split():
    assert optimal_split([(1.0, 2.0), (2.0, 1.0)], 0.0) == [0.0, 0.0]
