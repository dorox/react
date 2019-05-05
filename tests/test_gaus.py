import numpy as np
import models
import tools

print('gaus')
print('------------')
r2 = models.CSTR(0.4167,2.2)
r2.inlet(c_in=tools.rect())#t1=20, sig=1.8))
r2.import_data('data/cstr_pulse.csv')
r2.time_start = 0
r2.time_stop = 200
sol = r2.run()
ax = r2.plot(all=True)
ax.plot(range(0,100))