from collections import OrderedDict
from time import time
import numpy as np
from scipy.integrate import solve_ivp, OdeSolution
from scipy.optimize import minimize
from scipy.interpolate import interp1d

class Flowsheet(object):
    def __init__(self):
        self.units = OrderedDict()
        self.time_start = 0
        self.time_stop = 100
    
    def connect(self, source, destination):
        source.destination = destination
        destination.source = source
        if source not in self.units.values():
            self.units[len(self.units)] = source
        if destination not in self.units.values():
            self.units[len(self.units)] = destination

    def _ode(self, t, y):
        # Vectorized
        # y: m:n matrix, m-variables, n-units
        dydt = np.empty_like(y)
        for i, u in self.units.items():
            dydt[:,i] = u._ode(t,y[:,i])
            if u.destination:
                u.destination.initial_values = dydt[:,i]
        return dydt

    # def _initialise(self):
    #     for i, u in self.units.items():
    #         if source:
    #             u._get_c_in
    def run(self, plot=True):
        #Simple consecutive solution
        for i, u in self.units.items():
            if u.source:
                u.variables = u.source.variables
                u.events = u.source.events
                for k in u.source.variables.keys():
                    u.initial_values[k] = interp1d(
                        u.source.solution['t'],
                        u.source.solution[k],
                        kind = 'linear',
                        assume_sorted=False
                    )
            _=u.run(plot)

    def __run_simultaneous(self):
        ''' 
        Running simulation of a flowsheet
        '''
        self.solution = OrderedDict()
        self.solution['t'] = []
        t = self.time_start
        self.time_stop = 100
        fs_events = np.empty((2,0))
        fs_initial_values = np.empty((2,0))

        for i, u in self.units.items():
            fs_events = np.append(fs_events, u.events, 2)
            iv = u.initial_values
            initial_values = np.zeros(len(u.initial_values))
            for i, k in enumerate(iv):
                if callable(iv[k]):
                    initial_values[i] = iv[k](t)
                else:
                    initial_values[i] = iv[k]
            np.append(fs_initial_values, initial_values, 2)

        sol = None
        t0 = time()
        while True:
            sol = solve_ivp(
                self._ode,
                [t, self.time_stop],
                fs_initial_values,
                method = 'BDF',
                vectorized=True,
                events = fs_events,
                # max_step = 0.1,
                # atol=1e-9,
                # rtol=1e-7,
                )

            for i,k in enumerate(self.solution):
                if not k=='t':
                    self.solution[k] = np.append(self.solution[k], sol.y[i])
                    initial_values[i] = sol.y[i][-1]
            self.solution['t'] = np.append(self.solution['t'], sol.t)
            
            #If aborted by an event: continue
            if sol.status == 1:
                for i, e in enumerate(sol.t_events):
                    if e: 
                        self.events.pop(i)
                        t = e[0]
            #If aborted by t_stop: end
            elif sol.status == 0:
                break          
           
        t1 = time()
        print(f'run time: {t1-t0:.3f}s')
        if plot:
            
            self.plot(all = True)
        
        return sol
