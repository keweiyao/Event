import numpy as np
import matplotlib.pyplot as plt

M = 1.3
M2 = M**2
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


Er = np.linspace(1.01*M, 100.*M, 200)
tr = np.linspace(0.12, 0.8, 32)
r1 = np.loadtxt('./tables/RQq2Qq.dat').T
r2 = np.loadtxt('./tables/RQg2Qg.dat').T
r3 = np.loadtxt('./tables/RQq2Qqg.dat').T
r4 = np.loadtxt('./tables/RQg2Qgg.dat').T
for i, T in enumerate(tx):
	plt.plot(Er, r1[i], 'r-')
	plt.plot(Er, r2[i], 'g-')
	plt.plot(Er, r3[i], 'b-')
	plt.plot(Er, r4[i], 'y-')
plt.show()
