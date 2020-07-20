import unittest
import numpy as np
from chemreact import models, Flowsheet, tools

rtol = 1e-2
plot = False


def area(r):
    return np.trapz(r.solution["A"], r.solution["t"])


def tolearance(r, val):
    a = area(r)
    return abs((a - val) / a)


class Test_Flowsheets(unittest.TestCase):
    def test_simple(self):
        r1 = models.CSTR()
        r1.inlet(A=tools.rect())
        _ = r1.run(False)
        r2 = models.CSTR()
        r2.inlet(A=0)
        _ = r2.run(False)
        r3 = models.CSTR()
        r4 = models.CSTR()

        f = Flowsheet()
        f.connect(r1, r2, r3, r4)
        f.run(plot)

    def test_chemistry(self):
        r1 = models.CSTR()
        r2 = models.CSTR()
        r3 = models.CSTR()
        r4 = models.CSTR()

        c1 = models.Chemistry()
        c1.reaction("A=>0.5B")
        r2.chemistry = c1

        c2 = models.Chemistry()
        c2.reaction("B=>0.5C")
        r3.chemistry = c2

        f = Flowsheet()
        f.connect(r1, r2, r3, r4)
        r1.inlet(A=tools.step())

        f.run(plot)

    def test_inlet(self):

        r1 = models.CSTR(V=5)
        r2 = models.PFR(V=50)
        r3 = models.CSTR(V=10)

        f = Flowsheet()
        f.connect(r1, r2, r3)
        f.time_stop = 200

        r1.inlet(A=tools.rect())
        # in 3rd r - overstepping
        r3.solver_params(max_step=0.1)
        f.run(plot)
