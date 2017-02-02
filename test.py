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
e1 = event.event(mode='static', static_dt=0.2, elastic=True, inelastic=False, detailed_balance=False, mass=1.3)


# Dynamic Meidum
#e1 = event.event(mode='dynamic', hydrofile=sys.argv[1], inelastic=True, detailed_balance=False)


f = h5py.File("particledata.hdf5", 'w')

e1.initialize_HQ(NQ=100000, E0=2.)

for i in range(27):
	print "t = ", e1.sys_time()
	dsp, dsx = e1.HQ_hist()
	f.create_dataset("%d-p"%(i*2), data=dsp)
	f.create_dataset("%d-x"%(i*2), data=dsx)
	status = e1.perform_hydro_step(StaticPropertyDictionary=medium)
	dsp, dsx = e1.HQ_hist()
	f.create_dataset("%d-p"%(i*2+1), data=dsp)
	f.create_dataset("%d-x"%(i*2+1), data=dsx)
	e1.reset_HQ(E0=(3.+i))
	if not status:
		break
f.close()

