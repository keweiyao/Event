import numpy as np
import matplotlib.pyplot as plt
import sys
import event

# A static medium dictionary
medium = {'Temp': 0.4, 
		  'Vx'	: 0.0, 
		  'Vy'	: 0.0, 
		  'Vz'	: 0.0}

# define an event:
# Default arguments:
# mode="dynamic", hydrofile=None, mass=1.3, elastic=True, inelastic=False, table_folder='./tables'
# 1. mode = dynamic / static
# 2. hydrofile = user input filename
# 3. mass = M_c, M_b, ...
# 4. elastic = True / False
# 5. inelastic = True / False
# 6. table = where to put the tabulated cross-secitons and scattering rates.

# Static Meidum, only elastic
e1 = event.event(mode='static', inelastic=True, detailed_balance=False, mass=1.3)

# Static Meidum, elastic + inelastic
#e1 = event.event(mode='static', elastic=True, inelastic=True)

# Dynamic Meidum, only elastic
#e1 = event.event(mode='dynamic', hydrofile=sys.argv[1], inelastic=False)

# Dynamic Meidum, elastic + inelastic
#e1 = event.event(mode='dynamic', hydrofile=sys.argv[1], inelastic=True)

e1.initialize_HQ(NQ=3000)

plt.figure(figsize=(10, 10))
for i in range(50):
	print "t = ", e1.sys_time()
	status = e1.perform_hydro_step(StaticPropertyDictionary=medium)
	e1.HQ_hist()
	if not status:
		break
plt.show()


