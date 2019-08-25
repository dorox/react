from collections import OrderedDict
from time import time
import matplotlib.pyplot as p
import numpy as np
from scipy.integrate import solve_ivp, OdeSolution
from scipy.sparse import bsr_matrix
from scipy.optimize import minimize
from . import tools, models

Chemistry = models.Chemistry
PFR = models.PFR

class Domain0d:
    '''
    Basic functions of a zero-dimension modelling domain
    '''
    def __init__(self):
        self.data = OrderedDict() #keeps imported data
        self.solution = OrderedDict() #keeps simulation results
        self.variables = OrderedDict() #keeps variables(species) names and inlet functions
        self.initial_values = OrderedDict() # Initial conditions as const or f(t)
        self.time_start = 0 #Start of simulation
        self.time_stop = 100   #End of simulation
        self.time_eval = None   #Time steps for simulation
        self._solver_params = None
        self.events = []
        self.stoichiometry = np.array([], ndmin = 2)
        self.orders = np.array([], ndmin = 2)

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
    
    def _sort_vars(self):
        '''
        Sort internal variables to simplify indexing during integration
        '''
        sorted_list = sorted(self.variables.keys())
        if type(self) is Chemistry:
            #sort stoichiometry
            for i, k in enumerate(self.variables.keys()):
                self.variables[k] = self.stoichiometry[:,i]
            #sorted stoichiometry
            self.stoichiometry = np.stack([self.variables[k] for k in sorted_list], 1)

            #sort orders
            for i, k in enumerate(self.variables.keys()):
                self.variables[k] = self.orders[:,i]
            #sorted orders
            self.orders = np.stack([self.variables[k] for k in sorted_list], 1)

        #sort initial values
        self.initial_values = OrderedDict({k:self.initial_values[k] for k in sorted_list})

        #sort species
        self.variables = OrderedDict({k:0 for k in sorted_list})

    def solver_params(self, **kwargs):
        self._solver_params = kwargs
        return
    
    def run(self, plot=True, output=False):
        ''' 
        Running simulation of a modelling domain
        '''
        #TODO: pre-initialise solver separately
        #TODO: optimise overhead
        
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
            params = {'t_eval' : self.time_eval,
                'method' : 'BDF',
                'events' : events,
                'dense_output': False,
                }
            if self._solver_params:
                params.update(self._solver_params)
            sol = solve_ivp(
                self._ode,
                [t, self.time_stop],
                sol_initial_values,
                **params
                )

            for i,k in enumerate(self.solution):
                if not k=='t':
                    self.solution[k] = np.append(self.solution[k], sol.y[i])
                    sol_initial_values[i] = sol.y[i][-1]
            self.solution['t'] = np.append(self.solution['t'], sol.t)
            
            #Dense output for solution interpolation
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
        
        if type(self) is PFR:
            print('fix for PFR')
            # self.solution['t'][1:]+=self.delay
        t1 = time()
        
        if plot: 
            self.plot(all = True)
            print(f'run time: {t1-t0:.3f}s')

        if output:
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
                t=np.linspace(self.time_start,self.time_stop,1000)
                l = ax.plot(
                    t, 
                    [self.initial_values[variable](t) for t in t],
                    label = 'inlet '+variable,
                    #drawstyle = 'steps-post'
                )
            if variable in self.data:
                ax.plot(
                    self.data['t'],
                    self.data[variable],
                    'o',
                    label = 'exp.' + variable,
                    color = l[-1].get_color()
                )
        #ax.set_xlim(self.time_start,self.time_stop)
        ax.legend()
        p.show()
        #todo: return axis object
        return ax

class Domain1d:
    '''
    Basic functions of a zero-dimension modelling domain
    '''
    def __init__(self):
        self.data = OrderedDict() #keeps imported data
        self.solution = OrderedDict() #keeps simulation results
        self.variables = OrderedDict() #keeps variables(species) names and inlet functions
        self.initial_values = OrderedDict() # Initial conditions as const or f(t)
        self.time_start = 0 #Start of simulation
        self.time_stop = 100   #End of simulation
        self.time_eval = None   #Time steps for simulation
        self._solver_params = None
        self.events = []

    def solver_params(self, **kwargs):
        self._solver_params = kwargs
        return
    
    def _ode(self, t, y):
        '''
        Empty ode function, 
        y has to be in the order of self.variables
        '''
        print('Empty ode function, y has to be in the order of self.variables')
        pass

    def run(self, plot=True, output=False):
        ''' 
        Running simulation of a modelling domain
        '''
        #TODO: pre-initialise solver separately
        #TODO: optimise overhead
        
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
            params = {'t_eval' : self.time_eval,
                'method' : 'BDF',
                'events' : events,
                'dense_output': False,
                }
            if self._solver_params:
                params.update(self._solver_params)
            sol = solve_ivp(
                self._ode,
                [t, self.time_stop],
                sol_initial_values,
                **params
                )

            for i,k in enumerate(self.solution):
                if not k=='t':
                    self.solution[k] = np.append(self.solution[k], sol.y[i])
                    sol_initial_values[i] = sol.y[i][-1]
            self.solution['t'] = np.append(self.solution['t'], sol.t)
            
            #Dense output for solution interpolation
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
        
        if type(self) is PFR:
            print('fix PFR')
            # self.solution['t'][1:]+=self.delay
        t1 = time()
        
        if plot: 
            self.plot(all = True)
            print(f'run time: {t1-t0:.3f}s')

        if output:
            return sol

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
                t=np.linspace(self.time_start,self.time_stop,1000)
                l = ax.plot(
                    t, 
                    [self.initial_values[variable](t) for t in t],
                    label = 'inlet '+variable,
                    #drawstyle = 'steps-post'
                )
            if variable in self.data:
                ax.plot(
                    self.data['t'],
                    self.data[variable],
                    'o',
                    label = 'exp.' + variable,
                    color = l[-1].get_color()
                )
        #ax.set_xlim(self.time_start,self.time_stop)
        ax.legend()
        p.show()
        #todo: return axis object
        return ax
