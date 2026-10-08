import math


def close(a, b, rel=1e-6, abs_=0.0):
    return math.isclose(a, b, rel_tol=rel, abs_tol=abs_)
