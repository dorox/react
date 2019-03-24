import numpy as np
from scipy.integrate import solve_ivp

class Reactor:
    def __init__(self, m_in = 1, k = 1):
        self.m_in = m_in #inlet mass flow
        self.time_start = 0
        self.time_stop = 10
        self.k = k
        pass

    def _m_in(t):
        return self.m_in

    def _m_out(t):
        return _m_in(t)

    def _c_in(t):

    def _mass_balance(self, t, m):
        dmdt = _m_in(t) - _m_out(t)
        return dmdt

    def _concentration_balance(self, t, c):
        dcdt = _m_in(t)*_c_in(t) - _m_out(t)*_c_out(t) + _c_gen(c)
        return dcdt
    
    def _c_gen(c):
        _dcdt = self.forces*c
        return np.sum(_dcdt, 0)

    def ode_solver(self):
        sol = solve_ivp(fun, [self.time_start,self.time_stop], [c0])
        return sol.t, sol.c

    def reaction(self):
        '''
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
        self.rate_constants = np.array([[k1], [k2]])
        self.forces = self.stochiometry_matrix*self.rate_constants