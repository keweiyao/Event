import numpy as np
import matplotlib.pyplot as plt
import sys

#folder = sys.argv[1]
M = 1.3
M2 = M**2
GeVm2_to_mb = 0.3849
s12 = np.concatenate([np.linspace(1.01*M, 5.*M, 50), np.linspace(5.*M2, 30.*M, 50)])
s34 = np.linspace(1.01*M, 30.*M, 50)

s56 = np.linspace(1.01*M, 30.*M, 60)
t12 = np.linspace(0.1, 0.8, 32)
t34 = np.linspace(0.1, 0.8, 16)
t56 = np.linspace(0.1, 0.8, 8)
dtx = np.linspace(0.1, 5., 10.)

x1 = np.loadtxt('./tables/XQq2Qq.dat').reshape(100, 32).T*GeVm2_to_mb
x2 = np.loadtxt('./tables/XQg2Qg.dat').reshape(100, 32).T*GeVm2_to_mb
x3 = np.loadtxt('./tables/XQq2Qqg.dat').reshape(50, 16, 10).T*GeVm2_to_mb
x4 = np.loadtxt('./tables/XQg2Qgg.dat').reshape(50, 16, 10).T*GeVm2_to_mb
f5 = np.loadtxt('./tables/XQqg2Qq.dat').reshape(60, 8, 30, 30).T
f6 = np.loadtxt('./tables/XQgg2Qg.dat').reshape(60, 8, 30, 30).T


"""
for i, T in enumerate(t12):
	scale1 = 1.0/T**2
	scale2 = 1.0/(1.0 - M2/s12**2)**2/T**2
	plt.plot(s12, x1[i]/scale1, 'r-')
	plt.plot(s12, x2[i]/scale2, 'g-')
plt.show()

for it, dt in enumerate(dtx):
	plt.subplot(2, 5, it+1)
	for iT, T in enumerate(t34):
		scale = 1.0*dt**2/T**2
		plt.plot(s34, x3[it, iT]/scale, 'b-')
		plt.plot(s34, x4[it, iT]/scale, 'y-')
	#plt.ylim(1e-3, 1e0)
	#plt.semilogy()
plt.show()
"""


"""
N = 10
A1 = np.linspace(0.501, 0.999, 20)
A2 = np.linspace(-0.999, 0.999, 10)
for i1, a1 in enumerate(A1[::2]):
	for i2, a2 in enumerate(A2[::1]):
		plt.subplot(N, N, i1*N + i2 + 1)
		for iT, T in enumerate(t56):
			xk = 0.5*(a1*a2 + a1 - a2)
			x2 = 0.5*(-a1*a2 + a1 + a2)
			scale = s56**2/T/T/x2
			plt.plot(s56, f5[i2, 2*i1, iT]/scale, 'k-')
			plt.plot(s56, f6[i2, 2*i1, iT]/scale, 'c-')
		#plt.ylim(0, 200)
		#plt.semilogy()
plt.show()
"""


def myplot(folder):
	Er = np.linspace(1.01*M, 100.*M, 100)
	Er2 = np.linspace(1.01*M, 30.*M, 30)
	tr = np.linspace(0.13, 0.75, 8)
	dtr = np.linspace(0.1, 5., 10.)
	dtr2 = np.linspace(0.1, 5., 10.)

	r1 = np.loadtxt('./%s/RQq2Qq.dat'%folder).reshape(100, 16).T
	r2 = np.loadtxt('./%s/RQg2Qg.dat'%folder).reshape(100, 16).T
	r3 = np.loadtxt('./%s/RQq2Qqg.dat'%folder).reshape(100, 8, 10).T
	r4 = np.loadtxt('./%s/RQg2Qgg.dat'%folder).reshape(100, 8, 10).T
	r5 = np.loadtxt('./%s/RQqg2Qq.dat'%folder).reshape(30, 8, 10).T
	r6 = np.loadtxt('./%s/RQgg2Qg.dat'%folder).reshape(30, 8, 10).T

	
	for it, dt in enumerate(dtr):
		plt.subplot(2, 5, it+1)
		for iT, T in enumerate(tr):
			scale12 = 1.0#T**0.4
			plt.plot(Er, r1[iT]/scale12, 'r-', lw=1, label=r'$Qq\rightarrow Qq$' if iT==0 else '')
			plt.plot(Er, r2[iT]/scale12, 'g-', lw=1, label=r'$Qg\rightarrow Qg$' if iT==0 else '')
			scale34 = 1.0#dt**2/Er**0.75*T**2.6
			plt.plot(Er, r3[it, iT]/scale34, 'y-', lw=2, label=r'$Qq\rightarrow Qqg$' if iT==0 else '')
			plt.plot(Er, r4[it, iT]/scale34, 'b-', lw=2, label=r'$Qg\rightarrow Qgg$' if iT==0 else '')
			scale56 = 1.0#dtr2[it]**2/Er2**1.*T**3.5
			plt.plot(Er2, r5[it, iT]/scale56, 'y--', lw=2, label=r'$Qqg\rightarrow Qq$' if iT==0 else '')
			plt.plot(Er2, r6[it, iT]/scale56, 'b--', lw=2, label=r'$Qgg\rightarrow Qg$' if iT==0 else '')
		#plt.ylim(0, 1.5)
		#plt.semilogy()


plt.figure(figsize=(15, 6))

myplot("tables")


plt.legend(loc='best', fontsize=20, framealpha=0.0)
plt.ylabel(r"$\Gamma$ [GeV]", size=20)
plt.xlabel(r"$E_{HQ}$ [GeV]", size=20)
plt.show()

