import reaction

r = reaction.Reaction()

r.reaction('4A+1B2O=>2C+3gd12')
print(r.species)
print(r.stoichiometry)
print(r.rate_constants)
