import unittest
import numpy as np
from context import react

rtol = 1e-2

plot = False

def area(r):
    return np.trapz(r.solution['A'], r.solution['t'])
def tolearance(r, val):
    a= area(r)
    return abs((a-val)/a)

class Test_Flowsheets(unittest.TestCase):

    def test_simple(self):
        r1= react.models.CSTR()
        r1.inlet(A=react.tools.rect())
        _=r1.run(False)
        r2 = react.models.CSTR()
        r2.inlet(A=0)
        _=r2.run(False)

        f = react.Flowsheet()
        f.connect(r1,r2)
        f.run()

    def test_chemistry(self):
        c = react.models.Chemistry()
        c.reaction('A=>B')
        c.reaction('B=>C')

        V=10
        r1= react.models.CSTR(V=V)
        r1.inlet(A=react.tools.step())
        r1.inlet(B=0)

        r2 = react.models.CSTR(V=V)
        r2.inlet(A=0)
        r2.add(c)

        r3 = react.models.CSTR(V=V)
        r3.inlet(A=0)

        f = react.Flowsheet()
        f.connect(r1,r2)
        f.connect(r2,r3)
        f.run()