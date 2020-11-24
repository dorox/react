import unittest
from chemreact.models2 import Domain, Variable, Constant
import matplotlib.pyplot as plt


class Test_Domain(unittest.TestCase):
    def test_descriptors(self):
        v = Variable("A", 0)
        d = Domain("d")
        d.add_variable(v)
        d._y = [123]
        d.A
        d.A = 1
        self.assertTrue(d.A.initial_value == 1)
        d._is_running = True
        self.assertTrue(d.A == 123)

    def test_explicit_var_add(self):
        d = Domain("d")
        d.a = Variable("a", 1)
        self.assertTrue(d.a in d._vars)
        self.assertTrue(d.a.initial_value == 1)

    def test_adding_things(self):
        vars = []
        doms = dict()
        for i in range(4):
            doms[i] = Domain(name=f"d{i}")
            for j in range(3):
                vars.append(Variable(f"V{len(vars)}", len(vars)))
                doms[i].add_variable(vars[-1])
        doms[0].add_subdomain(doms[1])
        doms[1].add_subdomain(doms[2])
        doms[0].add_subdomain(doms[3])
        doms[0]._upd_idx()
        for i, v in enumerate(vars):
            self.assertTrue(v._idx == i)
        self.assertTrue(doms[0].y0 == list(range(len(vars))))

    def test_multiple_vars(self):
        d = Domain("d1")
        d.new_variables(["a", "b", "c"])
        d.new_constants(["D", "E", "F"])

    def test_run(self):
        d = Domain("d0")
        d.new_variables(("A", "B", "C"))
        d.A, d.B = (1, 1)

        def ode(t, y):
            dc = d.A * d.B
            return [-dc, -dc, dc]

        d.ode = ode
        sol = d.run()

    def test_lorentz(self):
        d = Domain("d0")
        d.new_variables(("x", "y", "z"))
        d.x, d.y, d.z = 1, 1, 1

        d.new_constants(("sig", "rho", "b"))
        d.sig, d.rho, d.b = 10, 28, 8 / 3

        def ode(t, y):
            dx = d.sig * (d.y - d.x)
            dy = d.x * (d.rho - d.z) - d.y
            dz = d.x * d.y - d.b * d.z
            return [dx, dy, dz]

        d.ode = ode
        sol = d.run()
        plt.plot(sol.y[0], sol.y[1])
        plt.plot(sol.y[1], sol.y[1])
        plt.plot(sol.y[2], sol.y[1])
        plt.show()

    def test_explicit(self):
        d = Domain("d")
        d.a.dt = lambda: d.a - d.b

        def dbdt():
            return d.a + d.b

        d.b.dt = dbdt
