import numpy as np
import models
import tools

print('const')
print('------------')
r = models.CSTR(1,10)
r.inlet(A=1)
r.run()

print('rect')
print('------------')
r2 = models.CSTR(1,10)
r2.inlet(A=tools.rect(t0=10, t1=10.1))
s = r2.run()

print('tri')
print('------------')
r2 = models.CSTR(1,10)
r2.inlet(A=tools.triangle())
r2.run()

print('ramp')
print('------------')
r2 = models.CSTR(1,10)
r2.inlet(A=tools.ramp())
r2.run()

print('step')
print('------------')
r2 = models.CSTR(1,10)
r2.inlet(A=tools.step())
r2.run()

print('gaus')
print('------------')
r2 = models.CSTR(1,10)
r2.inlet(A=tools.gaussian())
r2.run()

print('exp')
print('------------')
r2 = models.CSTR(1,10)
r2.inlet(A=tools.exponential())
r2.run()


c = models.Chemistry()
c.reaction('A=>B')
c.initial_concentrations(A=1)
c.run()

r.add(c)
r.run()

r.inlet(A = tools.rect())
r.run()