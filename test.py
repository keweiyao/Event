import event
import sys
import h5py
import numpy as np
#import matplotlib.pyplot as plt

def dfdE(e, T, M):
	return np.exp(-e/T)*np.sqrt(e**2-M**2)*e
def dfdp(p, T, M):
	x = np.sqrt(p**2+M**2)/T
	return (x+1.)*np.exp(-x)

def corner(ds, ranges, bins=50):
	N = ds.shape[0]
	for i in range(N):
		for j in range(i+1):
			plt.subplot(N, N, i*N+j+1)
			if i==j:
				plt.hist(ds[i], bins=bins, range=[ranges[i,0], ranges[i,1]], histtype='step', normed=True)
				plt.xlim(ranges[i,0], ranges[i,1])
			else:
				plt.hist2d(ds[j], ds[i], range=[[ranges[j,0], ranges[j,1]],[ranges[i,0], ranges[i,1]]], bins=bins)
				plt.xlim(ranges[j,0], ranges[j,1])
				plt.ylim(ranges[i,0], ranges[i,1])


e1 = event.event(sys.argv[1])
e1.initialize_HQ()
for i in range(10):
	e1.perform_hydro_step()

"""
ds = h5py.File("history.hdf5", 'r')

T = 0.3
ranges = np.array([[0, 11], [-5, 5], [-5, 5], [-5, 10]])
e = np.linspace(1.3, 11, 1000)
de = e[1]-e[0]
dpde = dfdE(e, T, 1.3)
dpde = dpde/np.sum(dpde)/de
p = np.linspace(-5, 5, 1000)
dp = p[1]-p[0]
dpdp = dfdp(p, T, 1.3)
dpdp = dpdp/np.sum(dpdp)/dp
for i in range(0, 1000, 10):
	plt.clf()	
	corner(ds['t=%d'%i].value.T, ranges)
	plt.subplot(4,4,1)
	plt.plot(e, dpde, 'r-')
	plt.subplot(4,4,6)
	plt.plot(p, dpdp, 'r-')
	plt.subplot(4,4,11)
	plt.plot(p, dpdp, 'r-')
	plt.subplot(4,4,16)
	plt.plot(p, dpdp, 'r-')
	plt.pause(0.01)
plt.show()
"""
