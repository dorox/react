from collections import OrderedDict
from time import time
import re
import csv
import matplotlib.pyplot as p
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import minimize

def rect(t0=10, t1=20, y0=0, y1=1):
    '''
    Rectangular function
    
    Parameters
    ----------
    t0 : number, optional
        Start time, by default 10
    t1 : number, optional
        Stop time, by default 20
    y0 : number, optional
        Initial value, by default 0
    y1 : number, optional
        Final value, by default 1
    
    Returns
    -------
    y(t):
        Function returnitng values at arbitrary time step
    [...]:
        Events array
    '''
    def _rect(t):
        if t<t0 or t>t1:
            return y0
        else:
            return y1
    
    def _event1(t,y): return t-t0
    _event1.terminal = True
    def _event2(t,y): return t-t1
    _event2.terminal = True
    return _rect, [_event1, _event2]

def triangle(t0=10, t1=15, t2=20, y0=0, y1=1):
    '''
    Triangular function
    Generates steady ramp up and steady ramp down
    
    Parameters
    ----------
    t0 : number, optional
        Start time of ramping up, by default 10
    t1 : number, optional
        Time of peak location, by default 15
    t2 : number, optional
        End time of ramp down, by default 20
    y0 : number, optional
        Initial value before ramp, by default 0
    y1 : number, optional
        Value at the peak, by default 1
    
    Returns
    -------
    y(t):
        Function returnitng values at arbitrary time step
    [...]:
        Events array
    '''
    dt1 = t1-t0
    dt2 = t2-t1
    def _triangle(t):
        if t0<t<=t1:
            return y0+y1*(t-t0)/dt1
        elif t1<t<=t2:
            return y0+y1-y1*(t-t1)/dt2
        else:
            return y0
    
    def _event1(t,y): return t-t0
    _event1.terminal = True
    def _event2(t,y): return t-t1
    _event2.terminal = True
    def _event3(t,y): return t-t2
    _event3.terminal = True

    return _triangle, [_event1, _event2, _event3]

def ramp(t0=10, t1=20, y0=0, y1=1):
    '''
    Linear rise until setpoint
    
    Parameters
    ----------
    t0 : number, optional
        Start time, by default 10
    t1 : number, optional
        End time, by default 20
    y0 : number, optional
        Intial value, by default 0
    y1 : number, optional
        End value, by default 1
    
    Returns
    -------
    y(t):
        Function returnitng values at arbitrary time step
    [...]:
        Events array
    '''
    dt = t1-t0
    def _ramp(t):
        if t0<t<t1:
            return y0+y1*(t-t0)/dt
        if t>=t1:
            return y1
        else:
            return y0

    def _event1(t,y): return t-t0
    _event1.terminal = True
    def _event2(t,y): return t-t1
    _event2.terminal = True
    return _ramp, [_event1, _event2]

def step(t1=10, y0=0, y1=1):
    '''
    Step function
    
    Parameters
    ----------
    t1 : number, optional
        Step time, by default 10
    y0 : number, optional
        Initial value, by default 0
    y1 : number, optional
        End value, by default 1
    
    Returns
    -------
    y(t):
        Function returnitng values at arbitrary time step
    [...]:
        Events array
    '''
    def _step(t):
        if t<t1:
            return y0
        else:
            return y1
    def _event1(t,y): return t-t1
    _event1.terminal = True
    return _step, [_event1]

def gaussian(t1=20, y_tot=1, sig=1):
    '''
    Gaussian pulse
    
    Parameters
    ----------
    t1 : number, optional
        Time of the peak, or expected value, by default 20
    y_tot : number, optional
        Area under the curve, by default 1
    sig : number, optional
        Variance, or peak width, by default 1
    
    Returns
    -------
    y(t):
        Function returnitng values at arbitrary time step
    '''
    y = y_tot /((2*np.pi)**0.5*sig)
    def _gaussian(t):
        return y*np.exp(-0.5*((t-t1)/sig)**2)
    return _gaussian

def exponential(t1=10, y_tot=1, c=1):
    '''
    Exponential growth until time t1
    
    Parameters
    ----------
    t1 : number, optional
        Time of the end of exponential curve, by default 10
    y_tot : number, optional
        Area under the curve, by default 1
    c : number, optional
        Rate of growth, by default 1
    
    Returns
    -------
    y(t):
        Function returnitng values at arbitrary time step
    [...]:
        Events array
    '''
    y = y_tot*c/np.exp(c*t1)
    def _exp(t):
        if t<t1:
            return y*np.exp(c*t)
        else:
            return 0
    def _event1(t,y): return t-t1
    _event1.terminal = True
    return _exp
    
def get_data(file_name):
    '''
    Reads data from file_name
    
    Parameters
    ----------
    file_name : string
        Filename relative to working directory
    
    Returns
    -------
    data : Odict
        {'t' : np.ndarray
        'Variable name' : np.ndarray
        '...' : np.ndarray
        }
    '''
    flen = sum(1 for line in open(file_name)) - 1
    with open(file_name) as csvfile:
        reader = csv.DictReader(csvfile)
        reader.fieldnames[0] = 't'
        data = OrderedDict().fromkeys(reader.fieldnames)
        i = 0
        for key in data:
            data[key] = np.zeros(flen)
        i = 0
        for row in reader:
            for key in row:
                data[key][i] = row[key]
            i+=1
    return data