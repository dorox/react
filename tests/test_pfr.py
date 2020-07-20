import unittest
import numpy as np
from chemreact import models, tools

rtol = 1e-2
plot = False
V = 10


def area(r):
    return np.trapz(r.solution["A"], r.solution["t"])


def tolearance(r, val):
    a = area(r)
    return abs((a - val) / a)


class Test_inputs(unittest.TestCase):
    def test_rect(self):
        r = models.PFR(V=V)
        r.inlet(A=tools.rect())
        _ = r.run(plot)


class Test_chemistry(unittest.TestCase):
    def test_simple(self):
        r = models.PFR(V=V)
        r.inlet(A=tools.gaussian())

        c = models.Chemistry()
        c.reaction("A=>B", k=0.01)
        r.chemistry = c
        _ = r.run(plot)
        1 + 1

