import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys

def calc(p0, p1, x0, x1):
	E0 = p0[0]
	pvec = p0[1:]
	q = (p1 - p0)[1:]
	dx = (x1 - x0)[1:]
	dt = (x1 - x0)[0]
	L = np.sqrt(np.dot(dx, dx))

	px, py, pz = pvec
	pt = np.sqrt(px**2+py**2)
	pabs = np.sqrt(pt**2+pz**2)
	R1 = np.array([[pz, 0., -pt],[0., pabs, 0.],[pt, 0., pz]])/pabs
	R2 = np.array([[px, py, 0.],[-py, px, 0.],[0., 0., pt]])/pt
	R = np.dot(R1, R2)
	#newip = np.dot(R, pvec)
	newiq = np.dot(R, q)
	return E0, L, dt, p1[0]-p0[0], newiq[0], newiq[1], newiq[2]
	
f = h5py.File(sys.argv[1])

T = 0.15 # GeV

M = 1.3 # GeV
DsL = []

E = 1.4*np.exp(0.45*np.linspace(0, 9, 10))
p = np.sqrt(E**2-M**2)
dEdx = []
Kperp = []
Kpara = []
A = []
for k in range(0, 10):
	print k
	Lp0 = f['%d-p'%(2*k)].value
	Lp1 = f['%d-p'%(2*k+1)].value
	Lx0 = f['%d-x'%(2*k)].value
	Lx1 = f['%d-x'%(2*k+1)].value
	E0, dL, dt, dE, qx, qy, qz = np.array([calc(p0, p1, x0, x1) for p0, p1, x0, x1 in zip(Lp0, Lp1, Lx0, Lx1)]).T
	mean_dpz = np.mean(qz)
	mean_dE = np.mean(dE)
	mean_dt = np.mean(dt)
	mean_dL = np.mean(dL)
	dEdx.append( -mean_dE/mean_dL )
	K = np.cov([qx, qy, qz])
	Kperp.append( 0.5*(K[0,0] + K[1,1])/mean_dt )
	Kpara.append( K[2,2]/mean_dt )
	A.append( -mean_dpz/mean_dt )
	
dEdx = np.array(dEdx)
Kperp = np.array(Kperp)
Kpara = np.array(Kpara)
A = np.array(A)
eta = A/p

plt.subplot(3,1,1)
plt.plot(E, dEdx, 'r-', label=r'$dE/dx$')
plt.plot(E, A, 'g-', label=r'$dE/dt$')
plt.semilogx()
plt.legend(loc='best', fontsize=20., framealpha=0.)
plt.xlabel(r'$E$ [GeV]')
plt.ylabel(r'Energy Loss [GeV/fm]')

plt.subplot(3,1,2)
plt.plot(E, Kperp, 'b-', label=r'$\kappa_\perp$')
plt.plot(E, Kpara, 'b--', label=r'$\kappa_\|$')
plt.semilogx()
plt.legend(loc='best', fontsize=20., framealpha=0.)
plt.xlabel(r'$E$ [GeV]')
plt.ylabel(r'$\kappa$ [GeV${}^2$/fm]')

plt.subplot(3,1,3)
plt.plot(E, Kperp/E**2/eta**2 * 2.*np.pi*T/0.197, 'b-', label=r'$2\pi T D_\perp$')
plt.plot(E, Kpara/E**2/eta**2 * 2.*np.pi*T/0.197, 'b--', label=r'$2\pi T D_\|$')
plt.semilogx()
plt.legend(loc='best', fontsize=20., framealpha=0.)
plt.xlabel(r'$E$ [GeV]')
plt.ylabel(r'$2\pi T D_s$')
plt.show()
	
