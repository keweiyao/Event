import numpy as np
import matplotlib.pyplot as plt
import sys
#folder = sys.argv[1]
M = 1.3
M2 = M**2

sx = np.concatenate([np.linspace(1.01*M2, 25.*M2, 50), np.linspace(25.*M2, 900.*M2, 50)])**0.5
tr = np.linspace(0.12, 0.8, 16)
dtr = np.linspace(0.1, 10., 10.)
x1 = np.loadtxt('./tables/XQq2Qq.dat').reshape(100, 16).T
x2 = np.loadtxt('./tables/XQg2Qg.dat').reshape(100, 16).T
x3 = np.loadtxt('./tables/XQq2Qqg.dat').reshape(100, 16, 10).T
x4 = np.loadtxt('./tables/XQg2Qgg.dat').reshape(100, 16, 10).T
index = 9
for i, T in enumerate(tr):
	plt.plot(sx, x1[i]*T**2, 'r-')
	plt.plot(sx, x2[i]*T**2, 'g-')
	plt.plot(sx, x3[index, i]/dtr[index]**2, 'b-')
	plt.plot(sx, x4[index, i]/dtr[index]**2, 'y-')
plt.semilogy()
plt.show()


def myplot(folder, scale=1.0):
	Er = np.linspace(1.01*M, 100.*M, 100)
	tr = np.linspace(0.12, 0.8, 16)
	dtr = np.linspace(0.1, 10., 10.)
	r1 = np.loadtxt('./%s/RQq2Qq.dat'%folder).reshape(100, 16).T
	r2 = np.loadtxt('./%s/RQg2Qg.dat'%folder).reshape(100, 16).T
	r3 = np.loadtxt('./%s/RQq2Qqg.dat'%folder).reshape(100, 16, 10).T/scale
	r4 = np.loadtxt('./%s/RQg2Qgg.dat'%folder).reshape(100, 16, 10).T/scale
	step = 1
	index = 1
	print dtr
	for i, T in enumerate(tr[::step]):
		plt.plot(Er, r1[i*step], 'r-', label=r'$Qq\rightarrow Qq$' if i==0 else '')
		plt.plot(Er, r2[i*step], 'g-', label=r'$Qg\rightarrow Qg$' if i==0 else '')
		plt.plot(Er, r3[index, i*step], 'b-', label=r'$Qq\rightarrow Qqg$' if i==0 else '')
		plt.plot(Er, r4[index, i*step], 'y-', label=r'$Qg\rightarrow Qgg$' if i==0 else '')
	plt.semilogy()

plt.figure(figsize=(15, 6))

myplot("tables")
#plt.ylim(1e-6, 1e0)
plt.xticks([0, 20, 40, 60, 80, 100, 120])
plt.legend(loc='best', fontsize=20, framealpha=0.0)
plt.ylabel(r"$\Gamma$ [GeV]", size=20)
plt.xlabel(r"$E_{HQ}$ [GeV]", size=20)
plt.show()

