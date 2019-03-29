import numpy as np
from scipy.integrate import solve_ivp
import re #Regular expression module for reading reaction strings


class Reaction:
    def __init__(self, m_in = 1):
        self.stoichiometry = np.array([], dtype= int)
        self.rate_constants = np.array([])
        self.species = dict()
        self.m_in = m_in #inlet mass flow
        self.time_start = 0
        self.time_stop = 10
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
    
    def rate(self, c):
        rate = []
        for reaction in self.stoichiometry:
            rate.append(c[0]**self.stoichiometry[0])
    def _c_gen(self, c):
        rates = self.rate_constants*self.rate(c)(c[0]**2*c[1], c[1]*(c[3]**3))
        return np.dot(rates, self.stoichiometry)

    def ode_solver(self):
        sol = solve_ivp(
            self._concentration_balance,
            [self.time_start,self.time_stop],
            self.c0,
            max_step = 0.1
            )
        return sol.t, sol.y

    def reaction(self, string):
        string = string.replace(' ', '')

        irreversible = r'=>'
        reversible = r'<=>'
               
        if re.search(reversible, string):
            [reagents, products] = re.split(reversible, string)
            self._new_reaction(reagents, products)
            self._new_reaction(products, reagents)
        elif re.search(irreversible, string):
            [reagents, products] = re.split(irreversible, string)
            self._new_reaction(reagents, products)
        else:
            print('reactions not found')
        
    def _new_reaction(self, reagents, products):
        
        def coeff(s): return int(re.findall(r'\b\d+', s)[0]) or 1
        def species(s): return re.findall(r'[A-z]\w*', s)[0]

        self.species.update({species(i):-coeff(i) for i in re.findall(r'\w+', reagents)})
        self.species.update({species(i):coeff(i) for i in re.findall(r'\w+', products)})
        
        self.stoichiometry = np.append(
            self.stoichiometry, 
            [self.species.get(i) for i in self.species]
            )
        
        self.rate_constants = [1]
        self.c0 = [1, 1]
