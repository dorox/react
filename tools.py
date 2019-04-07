from collections import OrderedDict
from time import time
import re
import csv
import matplotlib.pyplot as p
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import minimize


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