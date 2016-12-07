import numpy as np
import matplotlib.pyplot as plt
import sys
#folder = sys.argv[1]
M = 1.3
M2 = M**2
"""
sx = np.concatenate([np.linspace(1.01*M2, 25.*M2, 60), np.linspace(25.*M2, 900.*M2, 60)])**0.5
tx = np.linspace(0.12, 0.8, 32)
x1 = np.loadtxt('./tables/XQq2Qq.dat').T
x2 = np.loadtxt('./tables/XQg2Qg.dat').T
x3 = np.loadtxt('./tables/XQq2Qqg.dat').T
x4 = np.loadtxt('./tables/XQg2Qgg.dat').T
for i, T in enumerate(tx):
	plt.plot(sx, x1[i], 'r-')
	plt.plot(sx, x2[i], 'g-')
	plt.plot(sx, x3[i], 'b-')
	plt.plot(sx, x4[i], 'y-')
plt.show()
"""
def myplot(folder, scale=1.0):
	Er = np.linspace(1.01*M, 100.*M, 200)
	tr = np.linspace(0.12, 0.8, 32)
	r1 = np.loadtxt('./%s/RQq2Qq.dat'%folder).T
	r2 = np.loadtxt('./%s/RQg2Qg.dat'%folder).T
	r3 = np.loadtxt('./%s/RQq2Qqg.dat'%folder).T/scale
	r4 = np.loadtxt('./%s/RQg2Qgg.dat'%folder).T/scale
	for i, T in enumerate(tr[::5]):
		plt.plot(Er, r1[i], 'r-', label=r'$Qq\rightarrow Qq$' if i==0 else '')
		plt.plot(Er, r2[i], 'g-', label=r'$Qg\rightarrow Qg$' if i==0 else '')
		plt.plot(Er, r3[i], 'b-', label=r'$Qq\rightarrow Qqg$ scaled by $\delta t^2$' if i==0 else '')
		plt.plot(Er, r4[i], 'y-', label=r'$Qg\rightarrow Qgg$ scaled by $\delta t^2$' if i==0 else '')
	
	plt.semilogy()

t = [0.1, 1.0, 10.0]
plt.figure(figsize=(15, 6))
for i in range(3):
	plt.subplot(1,3,i+1)
	myplot("tables-%1.1f"%t[i], t[i]**2)
	plt.ylim(1e-6, 1e1)
	plt.yticks([])
	plt.xticks([0, 20, 40, 60, 80, 100, 120])
	if i==2:
		plt.legend(loc='best', fontsize=20, framealpha=0.0)
	if i==0:
		plt.ylabel(r"$\Gamma$ [GeV]", size=20)
		plt.yticks([1e-6, 1e-4, 1e-2, 1e0])
	plt.xlabel(r"$E_{HQ}$ [GeV]", size=20)
	plt.title(r"$\delta t / \tau_k = %1.1f$"%t[i], size=20)
plt.subplots_adjust(wspace=0.0, bottom=0.15)
plt.show()

