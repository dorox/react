import numpy as np
import unittest
import timeit
import chemreact

plot = False

models = chemreact.models
tools = chemreact.tools
rtol = 1e-2


def area(s, r):
    return np.trapz(r.solution[s], r.solution["t"])


def tolearance(s, r, val):
    a = area(s, r)
    return abs((a - val) / a)


def tol(x, val):
    return abs(x - val) / x


class Test_Chem(unittest.TestCase):
    def test_source(self):
        c = models.Chemistry()
        c.reaction("X=>2X")
        sol = c.run(plot=plot, output=True)
        self.assertTrue(sol.success)
        c.reaction("X+Y=>Y")
        sol = c.run(plot=plot, output=True)
        self.assertTrue(sol.success)

    def test_timing(self):
        # Reversible reaction test
        c = models.Chemistry()
        c.reaction("A=>B")
        c.reaction("A+X<=>C+D")
        c.reaction("B+Y<=>C+E")
        c.reaction("E+D<=>F")
        c.initial_concentrations(A=1, X=1, B=1, Y=1)
        c.time_stop = 50

        def run():
            # c.solver_params(jac=None)
            c.run(plot=False)
            return

        n = 100
        t = timeit.timeit(run, number=n)
        t_n = t / n
        print(f"\nTime per simulation: {t_n:0.3f}s\n")
        self.assertLessEqual(t_n, 0.4)

    def test_brussaelator(self):
        # Brusselator test
        c = models.Chemistry()

        c.reaction("A=>X")
        c.reaction("2X+Y=>3X")
        c.reaction("B+X=>Y+D")
        c.reaction("X=>E")

        c.initial_concentrations(A=1, B=3, X=1, Y=1)

        # r(A) = 0
        c.stoichiometry[:, 0] = np.zeros(4)
        # r(B) = 0
        c.stoichiometry[:, 1] = np.zeros(4)

        c.time_stop = 100
        # c.solver_params(jac=None)
        sol = c.run(plot=plot, output=True)
        self.assertTrue(sol.success)
        self.assertLessEqual(tolearance("A", c, 100), rtol)
        self.assertLessEqual(tolearance("B", c, 300), rtol)
        self.assertLessEqual(tolearance("X", c, 97.603), rtol)
        self.assertLessEqual(tolearance("Y", c, 318.307), rtol)
        self.assertLessEqual(tolearance("E", c, 4784.089), rtol)
        self.assertLessEqual(tolearance("D", c, 14352.268), rtol)

    def test_reaction(self):
        # Reversible reaction test
        c = models.Chemistry()
        c.reaction("A=>B")
        c.reaction("A+X<=>C+D")
        c.reaction("B+Y<=>C+E")
        c.reaction("E+D<=>F")
        c.initial_concentrations(A=1, X=1, B=1, Y=1)
        c.time_stop = 50
        sol = c.run(plot=plot, output=True)
        self.assertTrue(sol.success)

    def test_import(self):
        c = models.Chemistry()
        c.reaction("A+X=>B+Y", k=1e-7)
        c.reaction("B+A=>C+Y", k=1e-7)
        c.initial_concentrations(A=1, X=80)
        c.import_data("data/Expt1.csv", plot)
        self.assertEqual(c.data["t"][0], c.time_start)
        self.assertEqual(c.data["t"][-1], c.time_stop)
        self.assertTrue(len(c.data) == 4)
        self.assertTrue(len(c.data["t"]) == 35)

    def test_fit(self):
        c = models.Chemistry()
        c.reaction("A+X=>B+Y", k=1e-7)
        c.reaction("B+A=>C+Y", k=1e-7)
        c.initial_concentrations(A=1, X=80)
        c.import_data("data/Expt1.csv", plot)
        c.fit(plot)
        self.assertLessEqual(tol(c.rate_constants[0], 7.436e-7), rtol)
        self.assertLessEqual(tol(c.rate_constants[1], 2.898e-7), rtol)

    if __name__ == "__main__":
        unittest.main()
