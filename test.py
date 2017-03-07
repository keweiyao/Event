import numpy as np
import matplotlib.pyplot as plt
import sys
import event
import h5py

# A static medium dictionary
medium = {'Temp': 0.4, 
		  'Vx'	: 0.0, 
		  'Vy'	: 0.0, 
		  'Vz'	: 0.0}

def pT_weight(pt, y):
	return 1.0/(1.0+pt**4)/np.cosh(y)

# define an event:
# Default arguments:
# mode="dynamic", hydrofile=None, mass=1.3, elastic=True, inelastic=False, table_folder='./tables'
# 1. mode = dynamic / static
# 2. hydrofile = user input filename
# 3. mass = M_c, M_b, ...
# 4. elastic = True / False
# 5. inelastic = True / False
# 6. table = where to put the tabulated cross-secitons and scattering rates.

# Static Meidum
#e1 = event.event(mode='static', static_dt=0.1, elastic=True, inelastic=True, detailed_balance=True, mass=1.3)
#e1 = event.event(mode='static', static_dt=0.5, elastic=True, inelastic=True, detailed_balance=True, mass=1.3)


# Dynamic Meidum
Taa = np.loadtxt(sys.argv[1])
e1 = event.event(mode='dynamic', hydrofile=sys.argv[2], inelastic=True, detailed_balance=False)

f = h5py.File("particledata.hdf5", 'w')

e1.initialize_HQ(NQ=1000, Taa=Taa)

for i in range(200):
	print "t = %1.2f [fm/c]"%e1.sys_time()
	status = e1.perform_hydro_step(StaticPropertyDictionary=medium)
	dsp, dsx = e1.HQ_hist()
	f.create_dataset('p-%d'%i, data=dsp)
	if not status:
		break
f.close()


