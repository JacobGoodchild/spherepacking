"""Erdős Problem #506: count distinct circles determined by planar point sets.

Exact integer arithmetic. A collinear triple determines a line, not a circle.
"""
import itertools
from math import gcd


def circle_key(p, q, r):
    """Canonical integer key (a,b,c,d) of a(x^2+y^2)+bx+cy+d=0 through p,q,r, or None if collinear."""
    rows = [(x * x + y * y, x, y, 1) for x, y in (p, q, r)]

    def det3(m):
        return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
                - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))

    minors = [det3([[row[c] for c in range(4) if c != k] for row in rows]) for k in range(4)]
    a, b, c, d = minors[0], -minors[1], minors[2], -minors[3]
    if a == 0:
        return None
    g = gcd(gcd(abs(a), abs(b)), gcd(abs(c), abs(d)))
    a, b, c, d = a // g, b // g, c // g, d // g
    if a < 0:
        a, b, c, d = -a, -b, -c, -d
    return (a, b, c, d)


def count_circles(points):
    keys = {circle_key(*t) for t in itertools.combinations(points, 3)}
    keys.discard(None)
    return len(keys)


def degenerate(points):
    """All on one line or all on one circle."""
    keys = {circle_key(*t) for t in itertools.combinations(points, 3)}
    return keys == {None} or (None not in keys and len(keys) == 1)
