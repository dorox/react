import numpy as np
import models
from tools import *

#Brusselator test
chem = models.Chemistry()

chem.reaction('A=>X')
chem.reaction('2X+Y=>3X')
chem.reaction('B+X=>Y+D')
chem.reaction('X=>E')

chem.initial_concentrations(A=1, B=3, X=1, Y= 1)

chem.stoichiometry[:,3] =np.zeros(4)
chem.stoichiometry[:,0] =np.zeros(4)

chem.time_stop = 100
chem.run()
chem.plot('X', 'Y')

#Reversible reaction test
chem2 = models.Chemistry()
chem2.reaction('A=>B')
chem2.reaction('A+X<=>C+D')
chem2.reaction('B+Y<=>C+E')
chem2.reaction('E+D<=>F')
chem2.initial_concentrations(A=1, X=1, B=1, Y=1)
chem2.time_stop = 50
chem2.run()

#Andrew test:
chem3 = models.Chemistry()
chem3.reaction('A+X=>B+Y')
chem3.reaction('B+A=>C+Y')
chem3.initial_concentrations(X=200)
# chem3.time_stop = 20
chem3.run()
chem3.import_data('data/Expt1.csv')
chem3.rate_constants = [1e-7,1e-7]
chem3.run()
chem3.fit()