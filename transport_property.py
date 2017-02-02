import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys

def calc(p0, p1, x0, x1):
	E0 = p0[0]
	pvec = p0[1:]
	q = (p1 - p0)[1:]
	dx = (x1 - x0)[1:]
	L = np.dot(dx, pvec)/np.dot(pvec, pvec)**0.5

	px, py, pz = pvec
	pt = np.sqrt(px**2+py**2)
	pabs = np.sqrt(pt**2+pz**2)
	R1 = np.array([[pz, 0., -pt],[0., pabs, 0.],[pt, 0., pz]])/pabs
	R2 = np.array([[px, py, 0.],[-py, px, 0.],[0., 0., pt]])/pt
	R = np.dot(R1, R2)
	#newip = np.dot(R, pvec)
	newiq = np.dot(R, q)
	return E0, L, p1[0]-p0[0], newiq[0], newiq[1], newiq[2], pabs
	
f = h5py.File(sys.argv[1])
T = 0.3 # GeV
M = 1.3 # GeV
DsL = []

E = np.linspace(2, 30, 27)
dEdx = []
print E
for k in range(0, 27):
	print k
	Lp0 = f['%d-p'%(2*k)].value
	Lp1 = f['%d-p'%(2*k+1)].value
	Lx0 = f['%d-x'%(2*k)].value
	Lx1 = f['%d-x'%(2*k+1)].value
	E0, L, dE, qx, qy, qz, pabs = np.array([calc(p0, p1, x0, x1) for p0, p1, x0, x1 in zip(Lp0, Lp1, Lx0, Lx1)]).T
	dEdx.append(-np.mean(dE/L, axis=0))
plt.plot(E, dEdx)
plt.show()
	
