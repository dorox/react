import models
import matplotlib.pyplot as p
#Brusselator test
chem = models.Chemistry()

chem.reaction('A=>X')
chem.reaction('2X+Y=>3X')
chem.reaction('B+X=>Y+D')
chem.reaction('X=>E')
print(chem.stoichiometry)
print(chem.species)

# chem.stoichiometry[:,3] =np.zeros(4)
# chem.stoichiometry[:,0] =np.zeros(4)

chem.time_stop = 100
chem.c0 = [100,1,1,300,0,0]
t, c = chem.ode_solver()

j=0
for i in c:
    p.plot(t, i, label = list(chem.species.keys())[j])
    j+=1
p.legend()
p.show()

