import models

chem = models.Chemistry()

chem.reaction('A+B=>C', k=2)
chem.reaction('2A+2B=>2C', k=3)
print(chem.stoichiometry)
print(chem.rate_constants)

