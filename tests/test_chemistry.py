import numpy as np
import unittest

from context import react
plot = False

models = react.models
tools = react.tools
rtol = 1e-2

def area(s, r):
    return np.trapz(r.solution[s], r.solution['t'])
def tolearance(s, r, val):
    a= area(s, r)
    return abs((a-val)/a)
def tol(x,val):
    return abs(x-val)/x

class Test_Chem(unittest.TestCase):

    def test_brussaelator(self):
        #Brusselator test
        chem = models.Chemistry()

        chem.reaction('A=>X')
        chem.reaction('2X+Y=>3X')
        chem.reaction('B+X=>Y+D')
        chem.reaction('X=>E')

        chem.initial_concentrations(A=1, B=3, X=1, Y= 1)

        chem.stoichiometry[:,3] =np.zeros(4)
        chem.stoichiometry[:,0] =np.zeros(4)

        chem.time_stop = 100
        sol = chem.run(plot)
        self.assertTrue(sol.success)
        self.assertLessEqual(tolearance('A',chem,100),rtol)
        self.assertLessEqual(tolearance('B',chem,300),rtol)
        self.assertLessEqual(tolearance('X',chem,97.603),rtol)
        self.assertLessEqual(tolearance('Y',chem,318.307),rtol)
        self.assertLessEqual(tolearance('E',chem,4784.089),rtol)
        self.assertLessEqual(tolearance('D',chem,14352.268),rtol)

    def test_reaction(self):
        #Reversible reaction test
        chem2 = models.Chemistry()
        chem2.reaction('A=>B')
        chem2.reaction('A+X<=>C+D')
        chem2.reaction('B+Y<=>C+E')
        chem2.reaction('E+D<=>F')
        chem2.initial_concentrations(A=1, X=1, B=1, Y=1)
        chem2.time_stop = 50
        sol = chem2.run(plot)
        self.assertTrue(sol.success)

    def test_import(self):
        c = models.Chemistry()
        c.reaction('A+X=>B+Y', k=1e-7)
        c.reaction('B+A=>C+Y', k=1e-7)
        c.initial_concentrations(A=1,X=80)
        c.import_data('data/Expt1.csv', plot)
        self.assertEqual(c.data['t'][0], c.time_start)
        self.assertEqual(c.data['t'][-1],c.time_stop)
        self.assertTrue(len(c.data)==4)
        self.assertTrue(len(c.data['t']==35))

    def test_fit(self):
        c = models.Chemistry()
        c.reaction('A+X=>B+Y', k=1e-7)
        c.reaction('B+A=>C+Y', k=1e-7)
        c.initial_concentrations(A=1,X=80)
        c.import_data('data/Expt1.csv', plot)
        c.fit(plot)
        self.assertLessEqual(tol(c.rate_constants[0],7.436e-7),rtol)
        self.assertLessEqual(tol(c.rate_constants[1],2.898e-7),rtol)

    if __name__=='__main__':
        unittest.main()