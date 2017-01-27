import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys

def calc(p0, p1, x0, x1):
	E0 = p0[0]
	pvec = p0[1:]
	q = (p1 - p0)[1:]
	dx = (x1 - x0)[1:]
	L = (x1 - x0)[0]#np.sqrt(np.dot(dx, dx))

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
for i in range(0,25,2):
	Lp0 = f['%d-p'%(i)].value
	Lp1 = f['%d-p'%(i+1)].value
	Lx0 = f['%d-x'%(i)].value
	Lx1 = f['%d-x'%(i+1)].value
	E0, L, dE, qx, qy, qz, pabs = np.array([calc(p0, p1, x0, x1) for p0, p1, x0, x1 in zip(Lp0, Lp1, Lx0, Lx1)]).T
	sqL = np.sqrt(L)
	dEL = dE/L
	qL = np.array([qx/L/pabs*E0, qy/L/pabs*E0, qz/L/pabs*E0]).T
	qsqL = np.array([qx/sqL, qy/sqL, qz/sqL]).T

	index = ((E0>=1.3))
	print "Neff = ", np.sum(index)
	#print np.mean(dEL[index], axis=0)
	E0etaD = np.mean(qL[index], axis=0)[2]
	kappa = 1./3.*(np.cov(qsqL.T)[0,0] + np.cov(qsqL.T)[1,1] + np.cov(qsqL.T)[2,2])
	
	Ds = kappa/(2. * E0etaD**2)/0.197*2.*np.pi*T
	DsL.append(Ds)
	print Ds
print np.mean(DsL)
	
