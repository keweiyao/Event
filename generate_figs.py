import numpy as np
import matplotlib.pyplot as plt
import sys
import event
import h5py

#-------------Thermal distribution function----------------------------------------
def dfdE(e, T, M):
	return np.exp(-e/T)*np.sqrt(e**2-M**2)*e

def dfdP(p, T, M):
	x = np.sqrt(p**2+M**2)/T
	return (x+1.)*np.exp(-x)

#-------------Corner plot function-------------------------------------------------
def corner(ds, ranges, bins=50, title=''):
	plt.clf()
	N = ds.shape[0]
	for i in range(N):
		for j in range(i+1):
			plt.subplot(N, N, i*N+j+1)
			if i==j:
				plt.hist(ds[i], bins=bins, range=[ranges[i,0], ranges[i,1]], histtype='step', normed=True)
				if i==0:
					
					e = np.linspace(1.3, 1.3*10, 100)
					de = e[1] - e[0]
					dfde = dfdE(e, 0.4, 1.3)
					dfde = dfde/np.sum(dfde)/de
					plt.plot(e, dfde, 'r-', linewidth=2.0)
				else:
					p = np.linspace(-1.3*5, 1.3*5, 100)
					dp = p[1] - p[0]
					dfdp = dfdP(p, 0.4, 1.3)
					dfdp = dfdp/np.sum(dfdp)/dp
					plt.plot(p, dfdp, 'r-', linewidth=2.0)
				plt.xlim(ranges[i,0], ranges[i,1])
			else:
				plt.hist2d(ds[j], ds[i], range=[[ranges[j,0], ranges[j,1]],[ranges[i,0], ranges[i,1]]], bins=bins)
				plt.xlim(ranges[j,0], ranges[j,1])
				plt.ylim(ranges[i,0], ranges[i,1])
	plt.title(title)
	plt.pause(0.1)

E0 = (10.**2+1.3**2)**0.5
#-------------Corner plot function-------------------------------------------------
def compare(ds1, ds2, ds3, ranges, bins=50, title=''):
	position = [1, 2, 4, 3]
	plt.clf()
	for i in range(4):
		plt.subplot(2,2,position[i])
		plt.hist(ds1[i], bins=bins, range=[ranges[i,0], ranges[i,1]], histtype='step', normed=True, color='red', linewidth=1.5, label=r'$2\rightarrow 2$' if i==0 else ' ')
		plt.hist(ds2[i], bins=bins, range=[ranges[i,0], ranges[i,1]], histtype='step', normed=True, color='green', linewidth=1.5, label=r'$2\rightarrow 2 + 2\rightarrow 3$' if i==0 else ' ')
		plt.hist(ds3[i], bins=bins, range=[ranges[i,0], ranges[i,1]], histtype='step', normed=True, color='blue', linewidth=1.5, label=r'$2\rightarrow 2 + 2\rightarrow 3 + 3\rightarrow 2$' if i==0 else ' ')
		if i==0:
			plt.title(title, size=25)
			e = np.linspace(1.3, 1.3*10, 10000)
			de = e[1] - e[0]
			dfde = dfdE(e, 0.4, 1.3)
			dfde = dfde/np.sum(dfde)/de
			plt.plot(e, dfde, 'k-', linewidth=2.0, label='Boltzmann $T = 0.4$ GeV')
		else:
			p = np.linspace(-1.3*5, 1.3*5, 10000)
			dp = p[1] - p[0]
			dfdp = dfdP(p, 0.4, 1.3)
			dfdp = dfdp/np.sum(dfdp)/dp
			plt.plot(p, dfdp, 'k-', linewidth=2.0)
			

		plt.xlim(ranges[i,0], ranges[i,1])
		plt.semilogy()
		plt.ylim(5e-3, 2e0)
		if position[i]==1:
			plt.plot([E0, E0], [1e-5, 10.], 'k--', linewidth=2., label='Initial')
			plt.xlabel("$E$ [GeV]", size=25)
		if position[i]==2:
			plt.plot([0., 0.], [1e-5, 10.], 'k--', linewidth=2.)
			plt.xlabel("$p_x$ [GeV]", size=25)
		if position[i]==2 or position[i]==4:
			plt.yticks([])
		if position[i]==3:
			plt.plot([10., 10.], [1e-5, 10.], 'k--', linewidth=2.)
			plt.xlabel("$p_z$ [GeV]", size=25)
		if position[i]==4:
			plt.plot([0., 0.], [1e-5, 10.], 'k--', linewidth=2.)
			plt.xlabel("$p_y$ [GeV]", size=25)
		if i==0:
			plt.legend(loc = 'upper center', framealpha=0., fontsize=15)
	plt.subplots_adjust(wspace=0., hspace=0.2)
	plt.pause(0.1)

f1 = h5py.File("particledata-22.hdf5", "r")
f2 = h5py.File("particledata-22-23.hdf5", "r")
f3 = h5py.File("particledata-22-23-32.hdf5", "r")
plt.figure(figsize=(12,12))
for i in range(100):
	if i%1 != 0:
		continue
	ds1 = f1['%d'%i].value
	ds2 = f2['%d'%i].value
	ds3 = f3['%d'%i].value
	#corner(ds1.T, ranges=np.array([[0, 11], [-4, 4], [-4, 4], [-4, 11]]), title="%1.1f"%(i*0.5))
	compare(ds1.T, ds2.T, ds3.T, ranges=np.array([[1.3, 11.], [-4, 4], [-4, 4], [-3, 11]]), title=r"$t_{cell} = %1.1f$ fm/c"%(i*0.5))
plt.show()
