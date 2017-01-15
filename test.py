import numpy as np
import matplotlib.pyplot as plt
import sys
import event
import h5py

#-------------Thermal distribution function----------------------------------------
def dfdE(e, T, M):
	return np.exp(-e/T)*np.sqrt(e**2-M**2)*e

def dfdP(p, T, M):
	x = np.sqrt(p**2+M**2)/T
	return (x+1.)*np.exp(-x)

#-------------Corner plot function-------------------------------------------------
def corner(ds, ranges, bins=50):
	plt.clf()
	N = ds.shape[0]
	for i in range(N):
		for j in range(i+1):
			plt.subplot(N, N, i*N+j+1)
			if i==j:
				plt.hist(ds[i], bins=bins, range=[ranges[i,0], ranges[i,1]], histtype='step', normed=True)
				if i==0:
					e = np.linspace(1.3, 1.3*10, 100)
					de = e[1] - e[0]
					dfde = dfdE(e, 0.4, 1.3)
					dfde = dfde/np.sum(dfde)/de
					plt.plot(e, dfde, 'r-', linewidth=2.0)
				else:
					p = np.linspace(-1.3*5, 1.3*5, 100)
					dp = p[1] - p[0]
					dfdp = dfdP(p, 0.4, 1.3)
					dfdp = dfdp/np.sum(dfdp)/dp
					plt.plot(p, dfdp, 'r-', linewidth=2.0)
				plt.xlim(ranges[i,0], ranges[i,1])
			else:
				plt.hist2d(ds[j], ds[i], range=[[ranges[j,0], ranges[j,1]],[ranges[i,0], ranges[i,1]]], bins=bins)
				plt.xlim(ranges[j,0], ranges[j,1])
				plt.ylim(ranges[i,0], ranges[i,1])
	plt.pause(0.05)

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
e1 = event.event(mode='static', elastic=True, inelastic=True, detailed_balance=True, mass=1.3)

# Static Meidum, elastic + inelastic
#e1 = event.event(mode='static', elastic=True, inelastic=True)

# Dynamic Meidum, only elastic
#e1 = event.event(mode='dynamic', hydrofile=sys.argv[1], inelastic=False)

# Dynamic Meidum, elastic + inelastic
#e1 = event.event(mode='dynamic', hydrofile=sys.argv[1], inelastic=True)


e1.initialize_HQ(NQ=10000)

f = h5py.File("particledata.hdf5", 'w')
plt.figure(figsize=(10, 10))
for i in range(100):
	print "t = ", e1.sys_time()
	status = e1.perform_hydro_step(StaticPropertyDictionary=medium)
	ds = e1.HQ_hist()
	f.create_dataset("%d"%i, data=ds)
	if not status:
		break
f.close()



