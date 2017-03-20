import numpy as np
#import matplotlib.pyplot as plt
import sys
import event
import h5py
from scipy.interpolate import interp1d

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
# A static medium dictionary
medium = {'Temp': 0.3, 'Vx'	: 0.0, 'Vy'	: 0.0, 'Vz'	: 0.0}
e1 = event.event(mode='static', static_dt=1., elastic=True, inelastic=True, detailed_balance=True, mass=1.3)
#e1 = event.event(mode='static', static_dt=1., elastic=True, inelastic=False, detailed_balance=False, mass=1.3)

# Dynamic Meidum
#e1 = event.event(mode='dynamic', hydrofile=sys.argv[2], inelastic=True, detailed_balance=True, mass=1.3)
#e1 = event.event(mode='dynamic', hydrofile=sys.argv[2], inelastic=False, detailed_balance=False, mass=1.3)
f = h5py.File("particledata.hdf5", 'w')

# Initialize HQ pT weight and xy-distribution

#TR = np.loadtxt(sys.argv[1]).T ## p=0.
#Taa = TR**2

#pt, dsigma = np.loadtxt(sys.argv[3]).T
#dfdpt2 = interp1d(pt, dsigma)

def pT_weight(pt, y):
	return dfdpt2(pt)

e1.initialize_HQ(NQ=10000)#, Taa=Taa, dxy=0.1, df_dpt2dy=pT_weight, oversample_power=-1.)

# Run Modgetel  
for i in range(100):
	print "t = %1.2f [fm/c]"%e1.sys_time()
	status = e1.perform_hydro_step(StaticPropertyDictionary=medium)
	if i%10 == 0:
		dsp, dsx, dsw = e1.HQ_hist()
		f.create_dataset('p-%d'%i, data=dsp)
		f.create_dataset('x-%d'%i, data=dsx)
		f.create_dataset('w-%d'%i, data=dsw)
	if not status:
		break
f.close()


