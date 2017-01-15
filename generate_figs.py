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
				H,x= np.histogram(ds[i], bins=bins, range=[ranges[i,0], ranges[i,1]])
				#print np.sum(H)
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
	plt.pause(0.05)

f = h5py.File(sys.argv[1], "r")
for i in range(100):
	ds = f['%d'%i].value
	corner(ds.T, ranges=np.array([[0, 11], [-4, 4], [-4, 4], [-4, 11]]), title="%1.1f"%(i*0.5))
	
plt.show()
