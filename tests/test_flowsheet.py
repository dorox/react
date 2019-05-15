import unittest
import numpy as np
from context import chemreact as react

rtol = 1e-2

plot = True

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
        r3 = react.models.CSTR()
        r4 = react.models.CSTR()

        f = react.Flowsheet()
        f.connect(r1,r2,r3,r4)
        f.run(plot)

    def test_chemistry(self):
        r1 = react.models.CSTR()
        r2 = react.models.CSTR()
        r3 = react.models.CSTR()
        r4 = react.models.CSTR()
        
        c1 = react.models.Chemistry()
        c1.reaction('A=>0.5B')
        r2.chemistry = c1

        c2 = react.models.Chemistry()
        c2.reaction('B=>0.5C')
        r3.chemistry = c2

        f = react.Flowsheet()
        f.connect(r1,r2, r3, r4)
        r1.inlet(A = react.tools.step())

        f.run(plot)
    
