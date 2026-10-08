# The maths behind the code

Short version of what each formula is and where it lives. Everything assumes impulsive burns between
circular orbits around one body. `mu` is the body's gravitational parameter.

## Basics (`orbits.py`)

Circular speed and the vis-viva equation:

```
v_circ = sqrt(mu / r)
v      = sqrt(mu * (2/r - 1/a))          # speed at radius r on an orbit with semi-major axis a
```

Time from one apsis to the other on an ellipse is half the period, `t = pi * sqrt(a^3 / mu)`.

## Hohmann transfer (`transfers.hohmann`)

One ellipse that touches both circles, `a = (r1 + r2) / 2`. Two burns:

```
dv1 = v_transfer(r1) - v_circ(r1)
dv2 = v_circ(r2)     - v_transfer(r2)
```

Normalised by `v1 = v_circ(r1)` and with `R = r2/r1`:

```
dv_H / v1 = (sqrt(2R/(1+R)) - 1) + (1/sqrt(R)) * (1 - sqrt(2/(1+R)))
```

## Bi-elliptic transfer (`transfers.bielliptic`)

Go out to an intermediate apoapsis `rb`, then come back in. Two ellipses (`a_A = (r1+rb)/2`,
`a_B = (r2+rb)/2`), three burns: raise apoapsis at `r1`, adjust periapsis at `rb`, circularise at `r2`.
The middle burn happens where the spacecraft is slowest, which is the whole point.

With `B = rb/r1`:

```
dv_BE / v1 = [sqrt(2B/(1+B)) - 1]
           + |sqrt(2R/(B(R+B))) - sqrt(2/(B(1+B)))|
           + |sqrt(2B/(R(R+B))) - 1/sqrt(R)|
```

As `B -> infinity` this tends to the bi-parabolic limit `(sqrt(2) - 1)(1 + 1/sqrt(R))`.

### The two crossover ratios (`analysis.find_crossover_ratios`)

* **R < 11.94**: Hohmann always wins. This is where the bi-parabolic limit just ties Hohmann.
* **11.94 < R < 15.58**: bi-elliptic wins only if `rb` is big enough.
* **R > 15.58**: any `rb > r2` beats Hohmann. At `rb = r2` the bi-elliptic *is* the Hohmann transfer, so
  15.58 is where nudging `rb` upward first starts to help.

The code finds both numerically (bisection) rather than hard-coding them, and the tests check them against
the textbook values.

## Plane changes (`plane_change.py`)

A burn that changes speed from `v_a` to `v_b` while turning the velocity vector by `theta` costs (law of cosines):

```
dv = sqrt(v_a^2 + v_b^2 - 2 v_a v_b cos(theta))
```

If the speed doesn't change this is the familiar `2 v sin(theta/2)`. Because combining the turn with a
speed change is cheaper than doing them separately, it pays to spread the plane change across burns.

For a Hohmann transfer with total plane change `delta_i` and `alpha` done at departure:

```
total(alpha) = dv(v_c1 -> v_t1, alpha) + dv(v_t2 -> v_c2, delta_i - alpha)
```

`plane_change.optimal_split` minimises this over `alpha` (and over two angles for a 3-burn bi-elliptic).
The cost is *not* convex (for equal speeds it is concave), so it uses a grid search that zooms in, instead
of a gradient method that could settle on the wrong answer.

For planes that differ in both inclination and node, `plane_angle` uses spherical trigonometry:

```
cos(angle) = cos(i1) cos(i2) + sin(i1) sin(i2) cos(delta_RAAN)
```

### Where the textbook threshold comes from

If the *entire* plane change is done at the apoapsis burn, a bi-elliptic transfer beats a direct plane change
only past about 38.9 degrees, with the bi-parabolic limit breaking even at 48.9 degrees. The 48.9 figure
falls straight out of `2 sin(theta/2) = 2(sqrt(2) - 1)`, and `tests/test_transfers.py` checks the 38.9 one
(bi-elliptic loses at 35 degrees and wins at 43, using `plane_change_mode="apoapsis"`). When the optimiser is
allowed to spread the plane change over all three burns, the break-even moves to lower angles.

## Timing (`timing.py`)

The target moves while you coast, so it must start `pi - n2 * t_transfer` ahead of you
(`n2` = target's angular rate). Windows repeat every synodic period `2*pi / |n1 - n2|`.

## Fuel (`fuel.py`)

Tsiolkovsky: `m_prop = m0 * (1 - exp(-dv / (Isp * g0)))`, applied burn by burn.

## Simulation (`simulation.py`)

Classical RK4 on `r'' = -mu r / |r|^3`, with a step size proportional to the local dynamical time
`sqrt(r^3/mu)` so periapsis passes of eccentric arcs get small steps automatically. Burns are applied to the
*simulated* velocity: rotate the current velocity direction about the line of nodes (the x axis) by the
burn's plane-change angle, scale to the planned speed, take the difference as the delta-v vector.
That means the check is independent of the formula for the burn cost.
