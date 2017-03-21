import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
from scipy.interpolate import interp1d

x0, y0 = np.loadtxt(sys.argv[4]).T
ppref = interp1d(x0, y0)
pptot = (y0*x0).sum()*(x0[1]-x0[0])
x0, y1 = np.loadtxt(sys.argv[5]).T
AAtot = (y1*x0).sum()*(x0[1]-x0[0])
AAinit = interp1d(x0, y1)
ratio = AAtot/pptot
print(ratio)

f1 = h5py.File(sys.argv[1])
initial_pT = f1['init_pT'].value
weight1 = initial_pT*AAinit(initial_pT)
f2 = h5py.File(sys.argv[2])
initial_pT = f2['init_pT'].value
weight2 = initial_pT*AAinit(initial_pT)
f3 = h5py.File(sys.argv[3])
initial_pT = f3['init_pT'].value
weight3 = initial_pT*AAinit(initial_pT)

for i in range(300):
	if i%10 != 0:
		continue
	ds = f1['p-%d'%i].value.T
	pT1 = (ds[1]**2 + ds[2]**2)**0.5
	ds = f2['p-%d'%i].value.T
	pT2 = (ds[1]**2 + ds[2]**2)**0.5
	ds = f3['p-%d'%i].value.T
	pT3 = (ds[1]**2 + ds[2]**2)**0.5
	plt.clf()
	H1, be = np.histogram(pT1, bins=100, range=[0.1, 50], normed=True, weights=weight1)
	H2, be = np.histogram(pT2, bins=100, range=[0.1, 50], normed=True, weights=weight2)
	H3, be = np.histogram(pT3, bins=100, range=[0.1, 50], normed=True, weights=weight2)
	dpt = be[1] - be[0]
	if i==0:
		x = (be[1:] + be[:-1])*0.5
		ppy = ppref(x)
	plt.plot(x, H1/x/(ppy)*AAtot, 'r-')
	plt.plot(x, H2/x/(ppy)*AAtot, 'g-')
	plt.plot(x, H3/x/(ppy)*AAtot, 'b-')
	plt.axis([0,50,0,1.4])
	plt.pause(0.1)
plt.show()
