import unittest
from chemreact.models2 import Domain, Variable, Constant
import matplotlib.pyplot as plt


class Test_Domain(unittest.TestCase):
    def test_descriptors(self):
        v = Variable("A", 0, index=0)
        d = Domain("d")
        d.add_variable(v)
        d._y = [123]
        d.A
        d.A = 1
        self.assertTrue(d.A.initial_value == 1)
        d._is_running = True
        self.assertTrue(d.A == 123)

    def test_adding_things(self):
        vars = []
        doms = dict()
        for i in range(4):
            doms[i] = Domain(name=f"d{i}")
            for j in range(3):
                vars.append(Variable(f"V{len(vars)}", len(vars), domain=doms[i]))
                doms[i].add_variable(vars[-1])
        doms[0].add_subdomain(doms[1].name, doms[1])
        doms[1].add_subdomain(doms[2].name, doms[2])
        doms[0].add_subdomain(doms[3].name, doms[3])
        doms[0]._upd_idx()
        for i, v in enumerate(vars):
            self.assertTrue(v._idx == i)
        self.assertTrue(doms[0].y0 == list(range(len(vars))))

    def test_run(self):
        d = Domain("d0")
        d.add_variable("A", initial_value=1)
        d.add_variable("B", initial_value=1)
        d.add_variable("C", initial_value=0)

        def ode(t, y):
            dc = d.A * d.B
            return [-dc, -dc, dc]

        d.ode = ode
        sol = d.run()

