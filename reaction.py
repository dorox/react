import numpy as np
from scipy.integrate import solve_ivp
import re #Regular expression module for reading reaction strings


class Reaction:
    def __init__(self, m_in = 1):
        self.stoichiometry = np.array([])
        self.rate_constants = np.array([])
        self.species = np.array([])
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
    
    def _c_gen(self, c):
        rates = self.rate_constants*(c[0]**2*c[1], c[1]*(c[3]**3))
        return np.dot(rates, self.stochiometry_matrix)

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
        species = re.findall(r'[A-z]\w*', string)
        
        irreversible = r'=>'
        reversible = r'<=>'
        
        species =  re.findall(r'[A-z]\w*', string)
        self.species = np.unique(
            np.append(self.species, species)
            )
        
        if re.search(reversible, string):
            [left, right] = re.split(reversible, string)
            self._new_reaction(species, left, right)
            self._new_reaction(species, right, left)
        elif re.search(irreversible, string):
            [left, right] = re.split(irreversible, string)
            self._new_reaction(species, left, right)
        else:
            print('reactions not found')
        
    def _new_reaction(self, species, left, right):
        
        coeff = [ -int(i) for i in re.findall(r'\b\d+', left)]
        species = re.findall(r'[A-z]\w*', left)
        stoichiometry = dict(zip(species, coeff))
        
        coeff = [ int(i) for i in re.findall(r'\b\d+', right)]
        species = re.findall(r'[A-z]\w*', right)
        stoichiometry.update(dict(zip(species, coeff)))
        
        self.stoichiometry = np.append(
            self.stoichiometry, 
            [stoichiometry.get(i, 0) for i in self.species]
            )
        
        print(stoichiometry)