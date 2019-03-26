import numpy as np
from scipy.integrate import solve_ivp
import re #Regular expression module for reading reaction strings


class Reaction:
    def __init__(self, m_in = 1, k = 1):
        self.m_in = m_in #inlet mass flow
        self.time_start = 0
        self.time_stop = 10
        self.k = k
        self.c0 = np.array((1, 1, 0, 1))
        pass

    def _m_in(self, t):
        return self.m_in

    def _m_out(self, t):
        return self._m_in(t)

    def _c_in(self, t):
        pass

    def _mass_balance(self, t, m):
        dmdt = self._m_in(t) - self._m_out(t)
        return dmdt

    def _concentration_balance(self, t, c):
        #dcdt = self._m_in(t)*self._c_in(t) - self._m_out(t)*self._c_out(t) + self._c_gen(c)
        dcdt = self._c_gen(c)
        return dcdt
    
    def _c_gen(self, c):
        rates = self.rate_constants*(c[0]**2*c[1], c[1]*(c[3]**3))
        return np.dot(rates, self.stochiometry_matrix)

    def ode_solver(self):
        sol = solve_ivp(self._concentration_balance, [self.time_start,self.time_stop], self.c0, max_step = 0.1)
        return sol.t, sol.y

    def reaction(self, string):
        string = '''
        2 A + B -> 3 C
        B + 3 D -> 4 A
        stochiometry matrix:
               R1   R2
        A     -2    4
        B     -1   -1
        C      3    0
        D      0   -3
        '''
        self.stochiometry_matrix = [
            [-2, -1, 3, 0],
            [4, -1, 0, -3]
        ]
        k1 = 1
        k2 = 1
        self.rate_constants = np.array((1,1))

        #TODO:
        s='''
        A +B=C
        1B+ 4D> C
        '''
        for reaction in s.splitlines():
            print(re.split(r'\W+',reaction))
    
        re.findall(r'\w+', s)
        re.split(r'\W+',s)
