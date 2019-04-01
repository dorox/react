import numpy as np
from scipy.integrate import solve_ivp
import re #Regular expression module for reading reaction strings


class Chemistry:
    def __init__(self):
        self.stoichiometry = np.array([], ndmin = 2)
        self.rate_constants = np.array([])
        self.orders = np.array([], ndmin = 2)
        self.species = dict()
        self.m_in = 1 #inlet mass flow
        self.time_start = 0
        self.time_stop = 10

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
        r = np.zeros(len(self.rate_constants))
        i=0
        for reaction in self.orders:
            if sum(abs(reaction))>0:
                r[i]=self.rate_constants[i]
                j=0
                for coeff in reaction:
                    r[i]*=c[j]**(coeff)
                    j+=1
                i+=1
        return r

    def _c_gen(self, c):
        rates = self.rate(c)
        return np.dot(rates, self.stoichiometry)

    def ode_solver(self):
        sol = solve_ivp(
            self._concentration_balance,
            [self.time_start,self.time_stop],
            self.c0,
            #max_step = 0.1
            )
        return sol.t, sol.y

    def reaction(self, string, k=1, k1=1, k2=1):

        for i in self.species:
            self.species[i] = 0

        string = string.replace(' ', '')

        irreversible = r'=>'
        reversible = r'<=>'
               
        if re.search(reversible, string):
            [reagents, products] = re.split(reversible, string)
            self._new_reaction(reagents, products, k1)
            self._new_reaction(products, reagents, k2)
        elif re.search(irreversible, string):
            [reagents, products] = re.split(irreversible, string)
            self._new_reaction(reagents, products, k)
        else:
            print('reactions not found')
        
    def _new_reaction(self, reagents, products, k):
        
        self.rate_constants = np.append(self.rate_constants, k)

        def coeff(s): 
            match = re.search(r'\b[\d\.]+', s)
            if match:
                return float(match.group())
            else:
                return 1

        def species(s): return re.findall(r'[A-z]\w*', s)[0]

        #---Creating stoichiometry matrix

        self.species.update({species(i):-coeff(i) for i in re.findall(r'[\.\w]+', reagents)})

        for i in re.findall(r'[\.\w]+', products):
            try:
                self.species[species(i)] += coeff(i)
            except KeyError:
                self.species.update({species(i):coeff(i)})

        new_s = [self.species.get(i) for i in self.species]
        s = self.stoichiometry

        if s.size<=1:
            self.stoichiometry = np.array(new_s, ndmin = 2)
        else:
            c = np.zeros((s.shape[0], len(new_s)))
            c[:s.shape[0], :s.shape[1]] = s
            c = np.vstack((c, new_s))
            self.stoichiometry = c
    
        #---Creating species rate order matrix

        for i in self.species:
            self.species[i] = 0
        self.species.update({species(i):abs(coeff(i)) for i in re.findall(r'[\.\w]+', reagents)})
        new_o = [self.species.get(i) for i in self.species]
        o = self.orders
        if o.size<=1:
            self.orders = np.array(new_o, ndmin = 2)
        else:
            c = np.zeros((o.shape[0], len(new_o)))
            c[:o.shape[0], :o.shape[1]] = o
            c = np.vstack((c, new_o))
            self.orders = c
