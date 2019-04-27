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
    Triangle function
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
    Ramp function
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
    '''
    y = y_tot /((2*np.pi)**0.5*sig)
    def _gaussian(t):
        return y*np.exp(-0.5*((t-t1)/sig)**2)
    return _gaussian

def exponential(t1=10, y_tot=1, c=1):
    '''
    Exponential injection
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
    Imports data from a csv file
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