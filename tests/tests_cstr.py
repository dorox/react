import numpy as np
import models
import tools

plot = True

print('const')
print('------------')
r = models.CSTR(1,10)
r.inlet(A=1)
r.run(plot=plot)

print('rect')
print('------------')
r2 = models.CSTR(1,10)
r2.inlet(A=tools.rect(t0=10, t1=10.1))
s = r2.run(plot=plot)

print('tri')
print('------------')
r2 = models.CSTR(1,10)
r2.inlet(A=tools.triangle())
r2.run(plot=plot)

print('ramp')
print('------------')
r2 = models.CSTR(1,10)
r2.inlet(A=tools.ramp())
r2.run(plot=plot)

print('step')
print('------------')
r2 = models.CSTR(1,10)
r2.inlet(A=tools.step())
r2.run(plot=plot)

print('gaus')
print('------------')
r2 = models.CSTR(1,10)
r2.inlet(A=tools.gaussian())
r2.run(plot=plot)

print('exp')
print('------------')
r2 = models.CSTR(1,10)
r2.inlet(A=tools.exponential())
r2.run(plot=plot)


c = models.Chemistry()
c.reaction('A=>B+B')
c.reaction('B=>C')
c.reaction('C=>D')
c.reaction('D+B=>A')
c.initial_concentrations(A=1)
c.run(plot=plot)

r.add(c)
r.run(plot=plot)

r.inlet(A = tools.step())
r.run(plot=True)