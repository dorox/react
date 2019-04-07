import re  # Regular expression module for reading reaction strings
from collections import OrderedDict
from time import time
import matplotlib.pyplot as p
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
from tools import get_data

class Domain:
    '''
    Basic functions of a modelling domain
    '''
    def __init__(self):
        self.data = OrderedDict() #keeps imported data
        self.solution = OrderedDict() #keeps simulation results
        self.variables = OrderedDict() #keeps variables(species) names
        self.initial_values = OrderedDict()
        self.time_start = 0 #Start of simulation
        self.time_stop = 100   #End of simulation
        self.time_eval = None   #Time steps for simulation

    def _ode(self, t, y):
        '''
        Empty ode function, 
        y has to be in the order of self.variables
        '''
        print('Empty ode function, y has to be in the order of self.variables')
        pass

    def import_data(self, file_name):
        '''
        Iporting .csv file with the format:

        time| variable1     | variable2     | ....
        ----------------------------------------
        0   | init.cond.1   | init.cond2.   |   ...
        ... |   ....        |   .....       |   ...

        '''
        data = get_data(file_name)
        self.data = data
        if self.data:
            self.time_start = self.data['t'][0]
            self.time_stop = self.data['t'][-1]
            self.time_eval = self.data['t']

        for k in data:
            if k != 't':
                p.plot(data['t'], data[k], 'o', label = k)
        p.legend()
        p.show()

    def run(self, plot = True):
        ''' 
        Running simulation of a modelling domain
        '''
        #TODO: add variable time-stepping - increasing time step with time
        t0 = time()
        sol = solve_ivp(
            self._ode,
            [self.time_start, self.time_stop],
            list(self.initial_values.values()),
            t_eval = self.time_eval,
            #method = 'BDF',
            )
        t1 = time()

        self.solution = OrderedDict.fromkeys(self.variables)
        
        i = 0
        for k in self.solution:
            self.solution[k] = sol.y[i]
            i+=1
        self.solution['t'] = sol.t

        if plot:
            print(f'run time: {t1-t0:.3f}s')
            self.plot(all = True)
        
        return sol

    def _obj_fun(self, parameters):
        self.rate_constants = parameters
        self.run(plot = False)
        obj = 0
        for key in self.data:
            if key != 't':
                obj += np.sum((self.solution[key] - self.data[key])**2)
        return obj

    def fit(self):
        '''
        Fitting the model into imported data
        '''
        res = minimize(
            self._obj_fun, 
            self.rate_constants,
            #bounds = tuple((0,None) for i in self.rate_constants),
            method = 'Nelder-Mead',
            )
        print(res)
        #self.rate_constants = res
        self.run()
        return res
    
    def plot(self, *args, all = False):
        '''
        Plotting selected variables as a function of time
        '''
        #TODO check if solution exists
        if all:
            to_plot = self.variables
        else:
            to_plot = args

        fig = p.figure()
        ax = fig.add_subplot(111)
        ax.set_ylabel('concentration')
        ax.set_xlabel('time')
        for variable in to_plot:
            ax.plot(
                self.solution['t'], 
                self.solution[variable],
                label = variable,
                )
            try:
                ax.plot(
                    self.data['t'],
                    self.data[variable],
                    'o',
                    label = 'exp.' + variable,
                )
            except KeyError:
                pass
        ax.legend()
        p.show()

class Chemistry(Domain):
    def __init__(self):
        super().__init__()
        self.stoichiometry = np.array([], ndmin = 2)
        self.rate_constants = np.array([])
        self.orders = np.array([], ndmin = 2)
        self.species = OrderedDict()
        self.c0 = OrderedDict()
        self.variables = self.species #will there be non-species variables?
        self.initial_values = self.c0 #temporary?
        self.parameters = self.rate_constants

    def _rate(self, c):
        '''
        Returns array of species generation rates
        '''
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

    def _ode(self, t, c):
        '''
        ODE to be solved:
        only species generation term - reaction rate
        '''
        dcdt = np.dot(self._rate(c), self.stoichiometry)
        return dcdt

    def initial_concentrations(self, **kwargs):
        '''
        Set the initial concentrations of species at t=0s
        '''
        #TODO error when kwards contains non-existing species
        for species, conc in kwargs.items():
            self.c0.update({species:conc})

    def reaction(self, string, k=1, k1=1, k2=1):
        '''
            Adding new chemical reaction to the system of reactions 
            in the current modelling domain.
            
            Parameters
            ----------
            string : str
                The chemical reaction in string form
            k : float
                rate constant in case of irreversible reaction
            k1 : float
                forward rate constant for reversible reaction
            k2 : float
                backward rate constant for reversible reaction
            
            Examples
            --------
            chem.reaction('A+B=>C')
            chem.reaction('A+B<=>C')
            chem.reaction('A+B=>C', k = 1e-6)
            chem.reaction('A+B<=>C', k1 = 3.5e-12, k2 = 4.5e-11)
        '''
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
        
        for s in self.species:
            try:
                self.c0.update({s:self.c0[s]})
            except KeyError:
                self.c0.update({s:0})

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
