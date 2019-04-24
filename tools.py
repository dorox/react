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
    return np.vectorize(_rect)

def ramp(t1=10, y0=0, y1=1):
    '''
    Ramp function
    '''
    def _ramp(t):
        if t<t1:
            return y0+t*y1/t1
        else:
            return y1
    return np.vectorize(_ramp)

def step(t1=10, y0=0, y1=1):
    '''
    Step function
    '''
    def _step(t):
        if t<t1:
            return y0
        else:
            return y1
    return np.vectorize(_step)

def gaussian(t1=20, y_tot=1, sig=1):
    '''
    Gaussian pulse
    '''
    y = y_tot /((np.pi)**0.5*sig)
    def _gaussian(t):
        return y*np.exp(-((t-t1)/sig)**2)
    return np.vectorize(_gaussian)

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
    return np.vectorize(_exp)
    
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