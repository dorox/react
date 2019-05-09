import re  # Regular expression module for reading reaction strings
from collections import OrderedDict
from time import time
import matplotlib.pyplot as p
import numpy as np
from scipy.integrate import solve_ivp, OdeSolution
from scipy.optimize import minimize
from . import tools

class Domain:
    '''
    Basic functions of a modelling domain
    '''
    def __init__(self):
        self.data = OrderedDict() #keeps imported data
        self.solution = OrderedDict() #keeps simulation results
        self.variables = OrderedDict() #keeps variables(species) names and inlet functions
        self.initial_values = OrderedDict() # Initial conditions as const or f(t)
        self.time_start = 0 #Start of simulation
        self.time_stop = 100   #End of simulation
        self.time_eval = None   #Time steps for simulation
        self.events = []

    def _ode(self, t, y):
        '''
        Empty ode function, 
        y has to be in the order of self.variables
        '''
        print('Empty ode function, y has to be in the order of self.variables')
        pass

    def import_data(self, file_name, plot=True):
        '''
        Iporting .csv file with the format:

        time| variable1     | variable2     | ....
        ----------------------------------------
        0   | init.cond.1   | init.cond2.   |   ...
        ... |   ....        |   .....       |   ...

        '''
        data = tools.get_data(file_name)
        self.data = data
        if self.data:
            self.time_start = self.data['t'][0]
            self.time_stop = self.data['t'][-1]
            self.time_eval = self.data['t']
            for k in data:
                if k!='t':
                    self.initial_values[k] = data[k][0]

        if plot:
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
        #TODO: add sliders to get initial estimates for k values
        
        self.solution = OrderedDict.fromkeys(self.variables, [])
        self.solution['t'] = []
        self.dense = None
        t = self.time_start
        events = self.events.copy()

        iv = self.initial_values
        sol_initial_values = np.zeros(len(self.initial_values))
        for i, k in enumerate(iv):
            if callable(iv[k]):
                sol_initial_values[i] = iv[k](t)
            else:
                sol_initial_values[i] = iv[k]
    
        sol = None
        t0 = time()
        while True:
            sol = solve_ivp(
                self._ode,
                [t, self.time_stop],
                sol_initial_values,
                t_eval = self.time_eval,
                # max_step = 0.1,
                method = 'BDF',
                events = events,
                dense_output=True,
                # atol=1e-9,
                # rtol=1e-7,
                )

            for i,k in enumerate(self.solution):
                if not k=='t':
                    self.solution[k] = np.append(self.solution[k], sol.y[i])
                    sol_initial_values[i] = sol.y[i][-1]
            self.solution['t'] = np.append(self.solution['t'], sol.t)
            if not self.dense:
                self.dense = sol.sol
            else:
                self.dense.interpolants.extend(sol.sol.interpolants)
                self.dense.t_max = sol.sol.t_max
                self.dense.n_segments+=sol.sol.n_segments
                self.dense.ts = np.append(self.dense.ts, sol.sol.ts)
                self.dense.ts_sorted = np.append(self.dense.ts_sorted, sol.sol.ts_sorted)
            #If aborted by an event: continue
            if sol.status == 1:
                for i, e in enumerate(sol.t_events):
                    if e: 
                        events.pop(i)
                        t = e[0]
            #If aborted by t_stop: end
            elif sol.status == 0:
                break          
           
        t1 = time()
        print(f'run time: {t1-t0:.3f}s')
        if plot: self.plot(all = True)
        
        return sol

    def _obj_fun(self, parameters):
        #TODO: sensitivity matrix output
        self.rate_constants = parameters
        self.run(plot = False)
        obj = 0
        for key in self.data:
            if key != 't':
                obj += np.sum((self.solution[key] - self.data[key])**2)
        return obj

    def fit(self, plot=True):
        '''
        Fitting the model into imported data
        '''
        #TODO: Change to Gauss-Newton method with sensitivity estimation 
        #TODO: add checkboxes to choose components to optimise for
        res = minimize(
            self._obj_fun, 
            self.rate_constants,
            #bounds = tuple((0,None) for i in self.rate_constants),
            method = 'Nelder-Mead',
            )
        print(res)
        #self.rate_constants = res
        self.run(plot)
        return res
    
    def plot(self, *args, all = False):
        '''
        Plotting selected variables as a function of time
        '''
        #TODO check if solution exists before plotting
        #      add checkboxes for components to plot
        if all:
            to_plot = self.variables
        else:
            to_plot = args

        fig = p.figure()
        ax = fig.add_subplot(111)
        ax.set_ylabel('concentration')
        ax.set_xlabel('time')
        for variable in to_plot:
            l = ax.plot(
                self.solution['t'], 
                self.solution[variable],
                label = variable,
                )
            if callable(self.initial_values[variable]):
                l = ax.plot(
                    self.solution['t'], 
                    [self.initial_values[variable](t) for t in self.solution['t']],
                    label = 'inlet '+variable,
                    #drawstyle = 'steps-post'
                )
            try:
                ax.plot(
                    self.data['t'],
                    self.data[variable],
                    'o',
                    label = 'exp.' + variable,
                    color = l[-1].get_color()
                )
            except KeyError:
                pass
        ax.set_xlim(self.time_start,self.time_stop)
        ax.legend()
        p.show()
        #todo: return axis object
        return ax

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
        #TODO error when kwargs contains non-existing species
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

class CSTR(Domain):

    def __init__(self, q=1, V=1):
        super().__init__()
        self.source=None
        self.destination=None
        self.q = q
        self.V = V
        self._chemistry = OrderedDict()
    
    def _get_c_in(self, t):
        #used to get inlet variables values
        iv = self.initial_values
        c_in = np.zeros(len(iv))
        for i, k in enumerate(iv):
            if callable(iv[k]):
                c_in[i] = iv[k](t)
            else:
                c_in[i] = iv[k]
        return c_in

    def _ode(self, t, c):
        dmdt = self.q/self.V*(self._get_c_in(t)-c)
        if len(self._chemistry)>0:
            for k in self._chemistry:
                dmdt += self._chemistry[k]._ode(t,c)
        return dmdt
    
    def inlet(self, **kwargs):
        for k, v in kwargs.items():
            if type(v)==tuple:
                self.initial_values.update({k:v[0]})
                self.events = v[1]
            else:
                self.initial_values.update({k:v})
            self.variables.update({k:0})
        self._iv = [v for v in self.initial_values.values()]

    def add(self, domain):
        '''
        Adds new domain into the reactor
        '''
        if type(domain) == Chemistry:
            self._chemistry.update({len(self._chemistry):domain})
            for k in domain.species:
                try:
                    self.variables[k]
                except KeyError:
                    self.initial_values.update({k:0})
                    self.variables.update({k:0})
        pass