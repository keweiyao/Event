import numpy as np
import matplotlib.pyplot as plt
import sys
#folder = sys.argv[1]
M = 1.3
M2 = M**2
GeVm2_to_mb = 0.3849
s12 = np.concatenate([np.linspace(1.01*M, 5.*M, 50), np.linspace(5.*M2, 30.*M, 50)])
s34 = np.linspace(1.01*M, 30.*M, 50)
s56 = np.linspace(1.01*M, 10.*M, 10)
t12 = np.linspace(0.1, 0.8, 32)
t34 = np.linspace(0.1, 0.8, 16)
t56 = np.linspace(0.1, 0.8, 10)
x1 = np.loadtxt('./tables/XQq2Qq.dat').reshape(100, 32).T*GeVm2_to_mb
x2 = np.loadtxt('./tables/XQg2Qg.dat').reshape(100, 32).T*GeVm2_to_mb
x3 = np.loadtxt('./tables/XQq2Qqg.dat').reshape(50, 16, 10).T*GeVm2_to_mb
x4 = np.loadtxt('./tables/XQg2Qgg.dat').reshape(50, 16, 10).T*GeVm2_to_mb
f5 = np.loadtxt('./tables/XQqg2Qq.dat').reshape(10, 10, 10, 10).T
f6 = np.loadtxt('./tables/XQgg2Qg.dat').reshape(10, 10, 10, 10).T
index= 9
index1 = 0
index2 = 4
for i, T in enumerate(t12):
	plt.plot(s12, x1[i], 'r-')
	plt.plot(s12, x2[i], 'g-')

for i, T in enumerate(t34):
	plt.plot(s34, x3[index, i], 'b-')
	plt.plot(s34, x4[index, i], 'y-')
for i, T in enumerate(t56):
	plt.plot(s56, f5[index1, index2, i]/16./np.pi/(s56**2-M**2)**2, 'k-')
	plt.plot(s56, f6[index1, index2, i]/16./np.pi/(s56**2-M**2)**2, 'c-')

plt.semilogy()
#plt.semilogx()
plt.show()

def myplot(folder, scale=1.0):
	Er = np.linspace(1.01*M, 100.*M, 100)
	tr = np.linspace(0.12, 0.8, 8)
	dtr = np.linspace(0.1, 20., 10.)
	
	Er2 = np.linspace(1.01*M, 20.*M, 20)
	r1 = np.loadtxt('./%s/RQq2Qq.dat'%folder).reshape(100, 16).T
	r2 = np.loadtxt('./%s/RQg2Qg.dat'%folder).reshape(100, 16).T
	r3 = np.loadtxt('./%s/RQq2Qqg.dat'%folder).reshape(100, 16, 10).T/scale
	r4 = np.loadtxt('./%s/RQg2Qgg.dat'%folder).reshape(100, 16, 10).T/scale
	r5 = np.loadtxt('./%s/RQqg2Qq.dat'%folder).reshape(20, 8, 10).T/scale
	r6 = np.loadtxt('./%s/RQgg2Qg.dat'%folder).reshape(20, 8, 10).T/scale
	step = 1
	index = 5
	print dtr
	for i, T in enumerate(tr[::step]):
		plt.plot(Er, r1[i*step], 'r-', label=r'$Qq\rightarrow Qq$' if i==0 else '')
		plt.plot(Er, r2[i*step], 'g-', label=r'$Qg\rightarrow Qg$' if i==0 else '')
		plt.plot(Er, r3[index, i*step], 'b-', label=r'$Qq\rightarrow Qqg$' if i==0 else '')
		plt.plot(Er, r4[index, i*step], 'y-', label=r'$Qg\rightarrow Qgg$' if i==0 else '')
		plt.plot(Er2, r5[index, i*step], 'k-', label=r'$Qqg\rightarrow Qq$' if i==0 else '')
		plt.plot(Er2, r6[index, i*step], 'c-', label=r'$Qgg\rightarrow Qg$' if i==0 else '')
	#plt.semilogy()

plt.figure(figsize=(15, 6))

myplot("tables")
#plt.ylim(1e-6, 1e0)
plt.xticks([0, 20, 40, 60, 80, 100, 120])
plt.legend(loc='best', fontsize=20, framealpha=0.0)
plt.ylabel(r"$\Gamma$ [GeV]", size=20)
plt.xlabel(r"$E_{HQ}$ [GeV]", size=20)
plt.show()

