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
N = 50
for k in range(0,4):
	Lp0 = np.concatenate([f['%d-p'%(i)].value for i in range(N*k,N*(k+1),2)], axis=0)
	Lp1 = np.concatenate([f['%d-p'%(i+1)].value for i in range(N*k,N*(k+1),2)], axis=0)
	Lx0 = np.concatenate([f['%d-x'%(i)].value for i in range(N*k,N*(k+1),2)], axis=0)
	Lx1 = np.concatenate([f['%d-x'%(i+1)].value for i in range(N*k,N*(k+1),2)], axis=0)
	E0, L, dE, qx, qy, qz, pabs = np.array([calc(p0, p1, x0, x1) for p0, p1, x0, x1 in zip(Lp0, Lp1, Lx0, Lx1)]).T
	sqL = np.sqrt(L)
	dEL = dE/L
	qL = np.array([qx/L/pabs*E0, qy/L/pabs*E0, qz/L/pabs*E0]).T
	qsqL = np.array([qx/sqL, qy/sqL, qz/sqL]).T

	index = ((E0>=1.3) & (E0<=1.31))
	print "Neff = ", np.sum(index)
	#print np.mean(dEL[index], axis=0)
	E0etaD = np.mean(qL[index], axis=0)[2]
	M = np.cov(qsqL.T)
	print M
	kappa = 1./3.*(M[0,0] + M[1,1] + M[2,2])
	print E0etaD
	Ds = kappa/(2. * E0etaD**2)/0.197*2.*np.pi*T
	DsL.append(Ds)
	print Ds
plt.plot(DsL)
plt.show()
print np.mean(DsL)
	
