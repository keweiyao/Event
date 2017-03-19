import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
from scipy.interpolate import interp1d

x0, y0 = np.loadtxt(sys.argv[2]).T
ppref = interp1d(x0, y0)
pptot = (y0*x0).sum()*(x0[1]-x0[0])
x0, y1 = np.loadtxt(sys.argv[3]).T
AAtot = (y1*x0).sum()*(x0[1]-x0[0])
AAinit = interp1d(x0, y1)
ratio = AAtot/pptot
print(ratio)

f = h5py.File(sys.argv[1])
initial_pT = f['init_pT'].value
weight = initial_pT*AAinit(initial_pT)
H0, be = np.histogram(initial_pT, bins=40, range=[0.1, 70], normed=True, weights=weight)
x = (be[1:] + be[:-1])*0.5
plt.plot(x, H0/x, 'r-')
plt.plot(x, AAinit(x)/AAtot, 'b-')
ppy = ppref(x)
plt.semilogy()
plt.show()

for i in range(268):
	if i%10 != 0:
		continue
	ds = f['p-%d'%i].value.T
	pT = (ds[1]**2 + ds[2]**2)**0.5
	plt.clf()
	H1, be = np.histogram(pT, bins=40, range=[0.1, 70], normed=True, weights=weight)
	#plt.hist(pT, bins=100, range=[0.1, 70], normed=True, weights=weight, histtype='step')
	#plt.semilogy()
	dpt = be[1] - be[0]
	plt.plot(x, H0/x/(ppy)*AAtot, 'b-')
	plt.plot(x, H1/x/(ppy)*AAtot, 'r-')
	plt.axis([0,50,0,2])
	plt.pause(0.1)
plt.show()
