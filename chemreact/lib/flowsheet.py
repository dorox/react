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
    
    def connect(self, *args):
        '''
        Create connections of units in series:

        connect(unit1, unit2, unit3, ...)
        '''
        if len(args)<=1:
            return
        
        for source, destination in ((args[i],args[i+1]) for i in range(len(args)-1)):
            source.destination = destination
            destination.source = source
            if source not in self.units.values():
                self.units[len(self.units)] = source
            if destination not in self.units.values():
                self.units[len(self.units)] = destination

    def run(self, plot=True):
        #Simple consecutive solution
        for u in self.units.values():
            u.time_start = self.time_start
            u.time_stop = self.time_stop
            if u.source:
                u.events = u.source.events
                u.q = u.source.q
                for k in u.source.variables.keys():
                    u.inlet(**{k:interp1d(
                        u.source.solution['t'],
                        u.source.solution[k],
                        kind = 'linear',
                        assume_sorted=False
                        )}
                    )
            _=u.run(plot)

