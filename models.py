import re  # Regular expression module for reading reaction strings
from collections import OrderedDict

import matplotlib.pyplot as p
import numpy as np
from scipy.integrate import solve_ivp


class Chemistry:
    def __init__(self):
        self.stoichiometry = np.array([], ndmin = 2)
        self.rate_constants = np.array([])
        self.orders = np.array([], ndmin = 2)
        self.species = OrderedDict()
        self.c0 = OrderedDict()
        self.time_start = 0
        self.time_stop = 100

    def _concentration_balance(self, t, c):
        dcdt = self._c_gen(c)
        return dcdt

    def rate(self, c):
        #TODO: rewrite as array ** array
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

    def initial_concentrations(self, **kwargs):
        #TODO error when kwards contains non-existing species
        for species, conc in kwargs.items():
            self.c0.update({species:conc})

    def run(self):
        self.ode_solver()
        self.plot(all = True)

    def plot(self, *args, all = False):
        if all:
            to_plot = self.species
        else:
            to_plot = args
        fig = p.figure()
        ax = fig.add_subplot(111)
        ax.set_ylabel('concentration')
        ax.set_xlabel('time')
        for species in to_plot:
            ax.plot(
                self.solution['t'], 
                self.solution[species],
                label = species,
                )
        ax.legend()
        p.show()
    
    def ode_solver(self):

        c0 = OrderedDict.fromkeys(self.species, 0)
        c0.update({i:self.c0[i] for i in self.c0})
        self.c0 = c0

        sol = solve_ivp(
            self._concentration_balance,
            [self.time_start,self.time_stop],
            list(self.c0.values())
            #max_step = 0.1
            )
        
        self.solution = OrderedDict.fromkeys(self.species)
        i = 0
        for k in self.solution:
            self.solution[k] = sol.y[i]
            i+=1
        self.solution['t'] = sol.t

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
        
        c0 = OrderedDict.fromkeys(self.species, 0)
        c0.update({i:self.c0[i] for i in self.c0})
        self.c0 = c0

    def _new_reaction(self, reagents, products, k):
        #TODO: rewrite from self.species to new dictionaries

        self.rate_constants = np.append(self.rate_constants, k)

        def get_coeff(s): 
            match = re.search(r'\b[\d\.]+', s)
            if match:
                return float(match.group())
            else:
                return 1

        def get_species(s): return re.findall(r'[A-z]\w*', s)[0]

        #---Creating stoichiometry matrix
        
        self.species.update({get_species(i):-get_coeff(i) for i in re.findall(r'[\.\w]+', reagents)})

        for i in re.findall(r'[\.\w]+', products):
            try:
                self.species[get_species(i)] += get_coeff(i)
            except KeyError:
                self.species.update({get_species(i):get_coeff(i)})

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
        
        self.species.update({get_species(i):abs(get_coeff(i)) for i in re.findall(r'[\.\w]+', reagents)})
        new_o = [self.species.get(i) for i in self.species]
        o = self.orders
        if o.size<=1:
            self.orders = np.array(new_o, ndmin = 2)
        else:
            c = np.zeros((o.shape[0], len(new_o)))
            c[:o.shape[0], :o.shape[1]] = o
            c = np.vstack((c, new_o))
            self.orders = c
        
        for i in self.species:
            self.species[i] = 0
