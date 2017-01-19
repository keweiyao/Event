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

#-------------Corner plot function-------------------------------------------------
def compare(ds1, ds2, ds3, ranges, bins=100, title=''):
	plt.clf()
	for i in range(4):
		plt.subplot(2,2,i+1)
		plt.hist(ds1[i], bins=bins, range=[ranges[i,0], ranges[i,1]], histtype='step', normed=True, color='red', linewidth=2.)
		plt.hist(ds2[i], bins=bins, range=[ranges[i,0], ranges[i,1]], histtype='step', normed=True, color='green', linewidth=2.)
		plt.hist(ds3[i], bins=bins, range=[ranges[i,0], ranges[i,1]], histtype='step', normed=True, color='blue', linewidth=2.)
		if i==0:
			e = np.linspace(1.3, 1.3*10, 1000)
			de = e[1] - e[0]
			dfde = dfdE(e, 0.4, 1.3)
			dfde = dfde/np.sum(dfde)/de
			plt.plot(e, dfde, 'k--', linewidth=1.0)
		else:
			p = np.linspace(-1.3*5, 1.3*5, 1000)
			dp = p[1] - p[0]
			dfdp = dfdP(p, 0.4, 1.3)
			dfdp = dfdp/np.sum(dfdp)/dp
			plt.plot(p, dfdp, 'k--', linewidth=1.0)
		plt.xlim(ranges[i,0], ranges[i,1])
	plt.title(title)
	plt.pause(0.1)

#f1 = h5py.File("particledata-22.hdf5", "r")
f2 = h5py.File("particledata-22-23.hdf5", "r")
#f3 = h5py.File("particledata-22-23-32.hdf5", "r")
for i in range(300):
	if i%3 != 0:
		continue
	#ds1 = f1['%d'%i].value
	ds2 = f2['%d'%i].value
	#ds3 = f3['%d'%i].value
	corner(ds2.T, ranges=np.array([[0, 11], [-4, 4], [-4, 4], [-4, 11]]), title="%1.1f"%(i*0.3))
	#compare(ds1.T, ds2.T, ds3.T, ranges=np.array([[0, 11], [-4, 4], [-4, 4], [-4, 11]]), title="%1.1f"%(i*0.5))
plt.show()
