import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
from scipy.interpolate import interp1d

x0, y0 = np.loadtxt(sys.argv[1]).T
ppref = interp1d(x0, y0)
pptot = (y0*x0).sum()*(x0[1]-x0[0])
x0, y1 = np.loadtxt(sys.argv[2]).T
AAtot = (y1*x0).sum()*(x0[1]-x0[0])
AAinit = interp1d(x0, y1)
ratio = AAtot/pptot
print(ratio)


flist = [h5py.File(fname, 'r') for fname in sys.argv[3:] ]

def get_Raa(pmu, pT0, AAtot):
	pT = np.sqrt(pmu[1]**2 + pmu[2]**2)
	weight = pT0*AAinit(pT0)
	H, be = np.histogram(pT, bins=25, range=[0.1, 50], normed=True, weights=weight)
	x = (be[1:] + be[:-1])*0.5
	ppy = ppref(x)
	return x, H/x/ppy*AAtot

def get_v2(pmu, pT0):
	N = 20
	n = 8
	pTbins = np.concatenate([np.linspace(0, 10, n)[:-1], np.linspace(10, 50, N+2-n)])
	pT = np.sqrt(pmu[1]**2 + pmu[2]**2)
	q2 = -(pmu[1]**2 - pmu[2]**2)/(pmu[1]**2 + pmu[2]**2)
	weight = pT0*AAinit(pT0)
	v2 = np.zeros(N)
	for i in range(N):
		index = (pT > pTbins[i]) & (pT < pTbins[i+1])
		v2[i] = np.average(q2[index], weights=weight[index])
	x = (pTbins[1:] + pTbins[:-1])*0.5
	return x, v2

for i in range(129):
	if i%10 != 0:
		continue
	plt.clf()
	
	for j, f in enumerate(flist):
		pmu = f['p-%d'%i].value.T
		pT0 = f['init_pT'].value
		
		plt.subplot(1,2,1)
		x, y = get_Raa(pmu, pT0, AAtot)
		plt.plot(x, y, label = sys.argv[3+j])
		
		plt.subplot(1,2,2)
		x, y = get_v2(pmu, pT0)
		plt.plot(x, y, 'o-', label = sys.argv[3+j])
		
	plt.subplot(1,2,1)
	plt.axis([0,50,0,1.4])
	plt.legend()
	
	plt.subplot(1,2,2)
	plt.plot([0,50],[0,0])
	plt.axis([0, 50, -0.2, 0.4])
	plt.legend()
	
	plt.pause(0.1)
plt.show()
