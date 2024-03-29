import unittest
import numpy as np
import chemreact

rtol = 1e-2

models = chemreact.models
tools = chemreact.tools
plot = False
r = models.CSTR(1, 10)
r.time_start = 0
r.time_stop = 100


def area(r):
    return np.trapz(r.solution["A"], r.solution["t"])


def tolearance(r, val):
    a = area(r)
    return abs((a - val) / a)


class Test_Inputs(unittest.TestCase):
    def test_const(self):
        r.inlet(A=1)
        r.run(plot)
        self.assertLessEqual(tolearance(r, 100), rtol)

    def test_rect(self):
        r.inlet(A=tools.rect())
        r.run(plot)
        self.assertLessEqual(tolearance(r, 10), rtol)

    def test_tri(self):
        r.inlet(A=tools.triangle())
        r.run(plot)
        self.assertLessEqual(tolearance(r, 5), rtol)

    def test_ramp(self):
        r.inlet(A=tools.ramp())
        r.run(plot)
        self.assertLessEqual(tolearance(r, 75), rtol)

    def test_step(self):
        r.inlet(A=tools.step())
        r.run(plot)
        self.assertLessEqual(tolearance(r, 80), rtol)

    def test_gaus(self):
        r.inlet(A=tools.gaussian())
        r.run(plot)
        self.assertLessEqual(tolearance(r, 1), rtol)

    def test_exp(self):
        r.inlet(A=tools.exponential())
        r.run(plot)
        self.assertLessEqual(tolearance(r, 1), rtol)

    def test_gaus_hard(self):
        r2 = models.CSTR(0.4167, 2.2)
        r2.inlet(A=tools.gaussian(t1=20, sig=1.8))
        r2.run(plot)
        self.assertLessEqual(tolearance(r2, 1), rtol)

    def test_logistic(self):
        r.inlet(A=tools.logistic())
        r.run(plot)
        self.assertLessEqual(tolearance(r, 80), rtol)


class Test_Reactions(unittest.TestCase):
    def test_simple(self):
        c = models.Chemistry()
        c.reaction("A=>B")
        c.initial_concentrations(A=1)
        sol = c.run(plot=plot, output=True)
        self.assertTrue(sol.success)

        r.inlet(C=1)
        r.chemistry = c
        r.run(plot=plot)
        r.inlet(A=tools.step())
        r.inlet(D=1)
        sol = r.run(plot=plot, output=True)
        self.assertTrue(sol.success)


if __name__ == "__main__":
    unittest.main()

