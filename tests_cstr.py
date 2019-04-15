import numpy as np
import models
import tools

r = models.CSTR(1,10)
r.inlet(A=1)
r.run()

c = models.Chemistry()
c.reaction('A=>B')
c.initial_concentrations(A=1)
c.run()

r.add(c)
r.run()

r.inlet(A = tools.rect())
r.run()