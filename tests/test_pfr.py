import unittest
import numpy as np
from context import chemreact as react

rtol = 1e-2

models = react.models
tools = react.tools
plot = True
V=10
def area(r):
    return np.trapz(r.solution['A'], r.solution['t'])
def tolearance(r, val):
    a= area(r)
    return abs((a-val)/a)

class Test_inputs(unittest.TestCase):
    def test_rect(self):
        r = models.PFR(V=V)
        r.inlet(A=tools.rect())
        _=r.run()

class Test_chemistry(unittest.TestCase):
    def test_simple(self):
        r = models.PFR(V=V)
        r.inlet(A = tools.gaussian())

        c = models.Chemistry()
        c.reaction('A=>B', k=0.01)
        r.chemistry = c
        _=r.run()
        1+1