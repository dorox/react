import re  # Regular expression module for reading reaction strings
from collections import OrderedDict
from time import time
import matplotlib.pyplot as p
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
from tools import *

class Chemistry:
    def __init__(self):
        self.data = OrderedDict()
        self.stoichiometry = np.array([], ndmin = 2)
        self.rate_constants = np.array([])
        self.orders = np.array([], ndmin = 2)
        self.species = OrderedDict()
        self.c0 = OrderedDict()
        self.time_start = 0
        self.time_stop = 100

    def _obj_fun(self, rate_constants):
        self.rate_constants = rate_constants
        self.run(plot = False)
        obj = 0
        for key in self.data:
            if key != 'Time (s)':
                obj += np.sum((self.solution[key] - self.data[key])**2)
        return obj

    def fit(self):
        res = minimize(
            self._obj_fun, 
            self.rate_constants
            )
        print(res)
        self.rate_constants = res
        self.run()
        return res

    def import_data(self, file_name):
        data = get_data(file_name)
        self.data = data
        for k in data:
            if k != 'Time (s)':
                p.plot(data['Time (s)'],data[k], 'o', label = k)
        p.legend()
        p.show()

    def _concentration_balance(self, t, c):
        dcdt = self._c_gen(c)
        return dcdt

    def rate(self, c):
        r = np.zeros(len(self.rate_constants))
        i=0
        for reaction in self.orders:
            r[i]=self.rate_constants[i]

            # This is slower than a loop:
            #r[i]=np.prod(np.power(c,reaction))

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

    def run(self, plot = True):
        #If experimental data is imported:
        if self.data:
            time_stop = self.data['Time (s)'][-1]
            t_eval = self.data['Time (s)']
            
        #If no exp data: runnig default solver:
        else:
            t_eval = None
            time_stop = self.time_stop
        
        t0 = time()
        self.ode_solver(time_eval=t_eval, time_stop=time_stop)
        t1 = time()

        if plot:
            print(f'run time: {t1-t0:.3f}s')
            self.plot(all = True)

    def plot(self, *args, all = False):
        #TODO check if solution exists
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
    
    def ode_solver(self, time_eval = None, time_stop = None):

        c0 = OrderedDict.fromkeys(self.species, 0)
        c0.update({i:self.c0[i] for i in self.c0})
        self.c0 = c0
        if time_stop == None:
            time_stop = self.time_stop

        sol = solve_ivp(
            self._concentration_balance,
            [self.time_start, time_stop],
            list(self.c0.values()),
            t_eval = time_eval
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
