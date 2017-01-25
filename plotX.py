import numpy as np
import matplotlib.pyplot as plt
import sys
import h5py

#folder = sys.argv[1]
M = 1.3
M2 = M**2
GeVm2_to_mb = 0.3849
"""
f1 = h5py.File('./tables/XQq2Qq.hdf5')
x1 = f1['Xsection-tab'].value
at = f1['Xsection-tab'].attrs
ss1 = np.concatenate([	np.linspace(at['sqrts_low'], at['sqrts_mid'], at['N_sqrt_half']),
						np.linspace(at['sqrts_mid'], at['sqrts_high'], at['N_sqrt_half'])])
T1 = np.linspace(at['T_low'], at['T_high'], at['N_T'])
f1.close()

f2 = h5py.File('./tables/XQg2Qg.hdf5')
x2 = f2['Xsection-tab'].value
at = f2['Xsection-tab'].attrs
ss2 = np.concatenate([	np.linspace(at['sqrts_low'], at['sqrts_mid'], at['N_sqrt_half']),
						np.linspace(at['sqrts_mid'], at['sqrts_high'], at['N_sqrt_half'])])
T2 = np.linspace(at['T_low'], at['T_high'], at['N_T'])
f2.close()

f3 = h5py.File('./tables/XQq2Qqg.hdf5')
x3 = f3['Xsection-tab'].value
at = f3['Xsection-tab'].attrs
ss3 = np.linspace(at['sqrts_low'], at['sqrts_high'], at['N_sqrt_half'])
T3 = np.linspace(at['T_low'], at['T_high'], at['N_T'])
dt3 = np.linspace(at['dt_low'], at['dt_high'], at['N_dt'])
f3.close()

f4 = h5py.File('./tables/XQg2Qgg.hdf5')
x4 = f4['Xsection-tab'].value
at = f4['Xsection-tab'].attrs
ss4 = np.linspace(at['sqrts_low'], at['sqrts_high'], at['N_sqrt_half'])
T4 = np.linspace(at['T_low'], at['T_high'], at['N_T'])
dt4 = np.linspace(at['dt_low'], at['dt_high'], at['N_dt'])
f4.close()

f5 = h5py.File('./tables/XQqg2Qq.hdf5')
x5 = f5['Xsection-tab'].value
at = f5['Xsection-tab'].attrs
ss5 = np.linspace(at['sqrts_low'], at['sqrts_high'], at['N_sqrt_half'])
T5 = np.linspace(at['T_low'], at['T_high'], at['N_T'])
a15 = np.linspace(at['a1_low'], at['a1_high'], at['N_a1'])
a25 = np.linspace(at['a2_low'], at['a2_high'], at['N_a2'])
f5.close()

f6 = h5py.File('./tables/XQgg2Qg.hdf5')
x6 = f6['Xsection-tab'].value
at = f6['Xsection-tab'].attrs
ss6 = np.linspace(at['sqrts_low'], at['sqrts_high'], at['N_sqrt_half'])
T6 = np.linspace(at['T_low'], at['T_high'], at['N_T'])
a65 = np.linspace(at['a1_low'], at['a1_high'], at['N_a1'])
a65 = np.linspace(at['a2_low'], at['a2_high'], at['N_a2'])
f6.close()

for i, T in enumerate(T1):
	scale1 = 1.0#/T**2
	scale2 = 1.0#/(1.0 - M2/ss1**2)**2/T**2
	plt.plot(ss1, x1[:,i]*scale1, 'r-')
	plt.plot(ss2, x2[:,i]*scale2, 'g-')
plt.show()

for it, dt in enumerate(dt3):
	plt.subplot(2, 5, it+1)
	for iT, T in enumerate(T3):
		scale = 1.0#*dt**2/T**2
		plt.plot(ss3, x3[:, iT, it]*scale, 'b-')
		plt.plot(ss4, x4[:, iT, it]*scale, 'y-')
plt.show()


N = 10
for i1, a1 in enumerate(a15[::1]):
	for i2, a2 in enumerate(a25[::1]):
		plt.subplot(N, N, i1*N + i2 + 1)
		for iT, T in enumerate(T5):
			x2 = 0.5*(-a1*a2 + a1 + a2)
			scale = 1.#*(ss5**2 - M2)/T**2/x2
			plt.plot(ss5, x5[:, iT, i1, i2]*scale, 'k-')
			plt.plot(ss6, x6[:, iT, i1, i2]*scale, 'c-')
plt.show()
"""

f1 = h5py.File('./tables/RQq2Qq.hdf5')
R1 = f1['Rates-tab'].value
at = f1['Rates-tab'].attrs
E1 = np.concatenate([	np.linspace(at['E1_low'], at['E1_high'], at['N_E1'])])
T1 = np.linspace(at['T_low'], at['T_high'], at['N_T'])
f1.close()

f2 = h5py.File('./tables/RQg2Qg.hdf5')
R2 = f2['Rates-tab'].value
at = f2['Rates-tab'].attrs
E2 = np.concatenate([	np.linspace(at['E1_low'], at['E1_high'], at['N_E1'])])
T2 = np.linspace(at['T_low'], at['T_high'], at['N_T'])
f2.close()

f3 = h5py.File('./tables/RQq2Qqg.hdf5')
R3 = f3['Rates-tab'].value
at = f3['Rates-tab'].attrs
E3 = np.concatenate([	np.linspace(at['E1_low'], at['E1_high'], at['N_E1'])])
T3 = np.linspace(at['T_low'], at['T_high'], at['N_T'])
dt3 = np.linspace(at['dt_low'], at['dt_high'], at['N_dt'])
f3.close()

f4 = h5py.File('./tables/RQg2Qgg.hdf5')
R4 = f4['Rates-tab'].value
at = f4['Rates-tab'].attrs
E4 = np.concatenate([	np.linspace(at['E1_low'], at['E1_high'], at['N_E1'])])
T4 = np.linspace(at['T_low'], at['T_high'], at['N_T'])
dt4 = np.linspace(at['dt_low'], at['dt_high'], at['N_dt'])
f4.close()

f5 = h5py.File('./tables/RQqg2Qq.hdf5')
R5 = f5['Rates-tab'].value
at = f5['Rates-tab'].attrs
E5 = np.concatenate([	np.linspace(at['E1_low'], at['E1_high'], at['N_E1'])])
T5 = np.linspace(at['T_low'], at['T_high'], at['N_T'])
dt5 = np.linspace(at['dt_low'], at['dt_high'], at['N_dt'])
f5.close()

f6 = h5py.File('./tables/RQgg2Qg.hdf5')
R6 = f6['Rates-tab'].value
at = f6['Rates-tab'].attrs
E6 = np.concatenate([	np.linspace(at['E1_low'], at['E1_high'], at['N_E1'])])
T6 = np.linspace(at['T_low'], at['T_high'], at['N_T'])
dt6 = np.linspace(at['dt_low'], at['dt_high'], at['N_dt'])
f6.close()

plt.figure(figsize=(15, 8))

for it, dt in enumerate(dt3):
	plt.subplot(2, 5, it+1)
	for iT, T in enumerate(T1):
		scale12 = 1.#T**0.4
		plt.plot(E1, R1[:, iT]/scale12, 'r-', lw=1, label=r'$Qq\rightarrow Qq$' if iT==0 and it==0 else '')
		plt.plot(E2, R2[:, iT]/scale12, 'g-', lw=1, label=r'$Qg\rightarrow Qg$' if iT==0 and it==0 else '')
	for iT, T in enumerate(T3):
		scale34 = 1.#dt**2/Er**0.75*T**2.6
		plt.plot(E3, R3[:, iT, it]/scale34, 'y-', lw=2, label=r'$Qq\rightarrow Qqg$' if iT==0 and it==1 else '')
		plt.plot(E4, R4[:, iT, it]/scale34, 'b-', lw=2, label=r'$Qg\rightarrow Qgg$' if iT==0 and it==1 else '')
		scale56 = 1.#dtr2[it]**2/Er2**1.4*T**3.5
		plt.plot(E5, R5[:, iT, it]/scale56, 'y--', lw=2, label=r'$Qqg\rightarrow Qq$' if iT==0 and it==2 else '')
		plt.plot(E6, R6[:, iT, it]/scale56, 'b--', lw=2, label=r'$Qgg\rightarrow Qg$' if iT==0 and it==2 else '')
	plt.ylim(1e-6, 2e0)
	if it == 0 or it == 5:
		plt.ylabel(r"$\Gamma$ [GeV]", size=20)
	if it >= 5:
		plt.xlabel(r"$E_{HQ}$ [GeV]", size=20)
	if it <5 and it > 0:
		plt.xticks([])
		plt.yticks([])
	if it <10 and it > 5:
		plt.yticks([])
	if it==0:
		plt.xticks([])
	plt.legend(loc='best', fontsize=20, framealpha=0.0)
	plt.title(r"$\langle \Delta t\rangle = %1.1f$ fm/c"%dt, size=15)
	plt.semilogy()
plt.subplots_adjust(wspace=0., hspace=0.15)
plt.show()


