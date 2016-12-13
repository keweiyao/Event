import sys
from libcpp cimport bool
from libcpp.vector cimport vector
from libc.stdlib cimport rand, RAND_MAX
from libc.math cimport *
from cython.operator cimport dereference as deref, preincrement as inc
import numpy as np
import matplotlib.pyplot as plt
sys.path.append('../HQ-Evo/')
sys.path.append('../Medium-Reader')
import HqEvo
import medium
from event cimport *

cdef double GeVm1_to_fmc = 0.197

#-------------Thermal distribution function----------------------------------------
def dfdE(e, T, M):
	return np.exp(-e/T)*np.sqrt(e**2-M**2)*e

def dfdp(p, T, M):
	x = np.sqrt(p**2+M**2)/T
	return (x+1.)*np.exp(-x)

#-------------Corner plot function-------------------------------------------------
def corner(ds, ranges, bins=50):
	N = ds.shape[0]
	for i in range(N):
		for j in range(i+1):
			plt.subplot(N, N, i*N+j+1)
			if i==j:
				plt.hist(ds[i], bins=bins, range=[ranges[i,0], ranges[i,1]], histtype='step', normed=True)
				plt.xlim(ranges[i,0], ranges[i,1])
			else:
				plt.hist2d(ds[j], ds[i], range=[[ranges[j,0], ranges[j,1]],[ranges[i,0], ranges[i,1]]], bins=bins)
				plt.xlim(ranges[j,0], ranges[j,1])
				plt.ylim(ranges[i,0], ranges[i,1])

#-------------C/Python Wrapper functions--------------------------------------------
cpdef double dot4(vector[double] & A, vector[double] & B):
	return product4(A, B)

cpdef print4(const vector[double] & A):
	print4vec(A)

cpdef rotate_ByEuler(vector[double] & A, double alpha, double beta, double gamma):
	cdef vector[double] Ap
	Ap.resize(4)
	rotate_Euler(A, Ap, alpha, beta, gamma)
	return Ap

cpdef rotate_ByAxis(vector[double] & A, double alpha, unsigned int dir):
	if not (dir in [1,2,3]):
		raise ValueError("Direction can only be 1(x), 2(y), 3(z)")
	cdef vector[double] Ap
	Ap.resize(4)
	rotate_axis(A, Ap, alpha, dir)
	return Ap

cpdef boost4_By3(vector[double] A, vector[double] v):
	cdef vector[double] Ap
	Ap.resize(4)
	boost_by3(A, Ap, v)
	return Ap

cpdef boost4_By4(vector[double] A, vector[double] u):
	cdef vector[double] Ap
	Ap.resize(4)
	boost_by4(A, Ap, u)
	return Ap

cpdef boost4_ByAxis(vector[double] A, double vd, unsigned int dir):
	if not (dir in [1,2,3]):
		raise ValueError("Direction can only be 1(x), 2(y), 3(z)")
	cdef vector[double] Ap
	Ap.resize(4)
	boost_axis(A, Ap, vd, dir)
	return Ap

cpdef boost4_ByCoM(vector[double] A, vector[double] B):
	cdef vector[double] Ap, Bp, Pcom
	Ap.resize(4)
	Bp.resize(4)
	Pcom.resize(4)
	cdef i
	for i in range(4):
		Pcom[i] = A[i] + B[i]
	go_to_CoM(Pcom, A, B, Ap, Bp)
	return Pcom, Ap, Bp

#-----------Particle data struct------------------------------------
cdef extern from "../src/utility.h":
	cdef struct particle:
		vector[double] p
		vector[double] x
		double t_last
		double weight
		bool freeze
		double T_dec
		vector[double] v_dec
		particle()
#-----------Event Class---------------------------------------------
cpdef freestream(vector[double] & x, vector[double] & p, double & dt):
	cdef vector[double] xnew = [x[0]+dt, x[1]+p[1]/p[0]*dt, x[2]+p[2]/p[0]*dt, x[3]+p[3]/p[0]*dt]
	return xnew

cdef class event:
	cdef object hydro_reader, hqsample, mode
	cdef double M
	cdef vector[particle] active_HQ
	cdef vector[particle] frzout_HQ
	cdef double tau0, dtau, tau
	X = []
	Y = []
	
	def __cinit__(self, mode="dynamic", hydrofile=None, mass=1.3, elastic=True, inelastic=False, table_folder='./tables'):
		self.mode = mode
		self.M = mass
		self.hydro_reader = medium.Medium(mode=mode, hydrofilename=hydrofile)
		self.hqsample = HqEvo.HqEvo(mass=mass, elastic=elastic, inelastic=inelastic, table_folder=table_folder)
		self.tau0 = self.hydro_reader.init_tau()
		self.dtau = self.hydro_reader.dtau()
		self.tau = self.tau0

	def sys_time(self):
		return self.tau

	cpdef initialize_HQ(self, NQ=100, XY_table=None, Pweight=None):
		print "Uniform sampling transverse momentum with in circle pt2 < Const"
		print "Zero longitudinal momentum"
		self.active_HQ.resize(NQ)
		cdef double pt, phipt, r, phir, E, pz, free_time
		cdef particle Q
		
		if XY_table==None and Pweight==None:
			print "Use default initialization pro"
		cdef vector[particle].iterator it = self.active_HQ.begin()
		while it != self.active_HQ.end():
			self.X.append([])
			self.Y.append([])
			pt = sqrt((9.*rand())/RAND_MAX)
			phipt = (2.*M_PI*rand())/RAND_MAX
			r = sqrt((4.*rand())/RAND_MAX)
			phir = (2.*M_PI*rand())/RAND_MAX
			E = sqrt(self.M**2 + pt**2)
			pz = 0.0
			deref(it).p = [E, pt*cos(phipt), pt*sin(phipt), pz]
			free_time = self.tau0/sqrt(1. - (pz/E)**2)
			deref(it).x = [0.0, r*cos(phir), r*sin(phir), 0.0]
			# free streaming to hydro starting time
			deref(it).x = freestream(deref(it).x, deref(it).p, free_time)
			inc(it)

	cpdef perform_hydro_step(self, StaticPropertyDictionary=None):
		status = True
		if self.mode == 'dynamic':
			status = self.hydro_reader.load_next()
		if self.mode == 'static':
			status = self.hydro_reader.load_next(StaticPropertyDictionary=StaticPropertyDictionary)
		self.tau += self.dtau
		cdef double t, x, y, z, tauQ
		cdef vector[double] pnew
		pnew.resize(4)
		cdef vector[particle].iterator it = self.active_HQ.begin()
		while it != self.active_HQ.end():
			if self.mode == 'static':
				while (deref(it).x[0] < self.tau):
					self.perform_HQ_step(it)
			if self.mode == 'dynamic':
				t, x, y, z = deref(it).x
				tauQ = sqrt(t**2 - z**2)
				while (tauQ < self.tau):
					self.perform_HQ_step(it)
					t, x, y, z = deref(it).x
					tauQ = sqrt(t**2 - z**2)
			inc(it)
		return status

	cdef perform_HQ_step(self, vector[particle].iterator it):
		cdef double t, x, y, z, t_elapse_lab, T, tauQ
		cdef int channel

		t, x, y, z = deref(it).x
		t_elapse_lab = t - deref(it).t_last
		
		tauQ = sqrt(t**2 - z**2)
		T, vx, vy, vz = self.hydro_reader.interpF(tauQ, [x, y, z, t], ['Temp', 'Vx', 'Vy', 'Vz'])		
		channel, dt, pnew = self.update_HQ(deref(it).p, [vx, vy, vz], T, t_elapse_lab)
		dt *= GeVm1_to_fmc
		if channel < 0:
			deref(it).x = freestream(deref(it).x, deref(it).p, dt)
		else:
			dt1 = (dt*rand())/RAND_MAX
			dt2 = dt - dt1
			deref(it).x = freestream(deref(it).x, deref(it).p, dt1)
			deref(it).x = freestream(deref(it).x, pnew, dt2)
			deref(it).p = pnew
			deref(it).t_last = t + dt1
		return

	cpdef update_HQ(self, vector[double] p1, vector[double] v3cell, double Temp, double t_elapse_lab):
		# Define local variables
		t_elapse_lab /= GeVm1_to_fmc
		cdef double E1_cell, alpha1_cell, beta1_cell, gamma1_cell, E2_cell, s
		cdef double t_elapse_cell, t_elapse_com
		cdef double alpha_com, beta_com, gamma_com
		cdef double p1z_cell_align, costheta2, sintheta2, cosphi2, sinphi2
		cdef double dt
		cdef int channel
		cdef vector[double] rv3cell = [-v3cell[0], -v3cell[1], -v3cell[2]]

		# Boost to cell frame and take down orientation of p1_cell
		cdef vector[double] p1_cell = boost4_By3(p1, v3cell)
		E1_cell = p1_cell[0]
		alpha1_cell = atan2(p1_cell[2], p1_cell[1]) + M_PI/2.
		beta1_cell = atan2(sqrt(p1_cell[1]**2+p1_cell[2]**2), p1_cell[3])
		gamma1_cell = 0.0

		# Boost time separation to cell frame
		cdef vector[double] dx4_lab = [t_elapse_lab, 			 t_elapse_lab*p1[1]/p1[0],
									   t_elapse_lab*p1[2]/p1[0], t_elapse_lab*p1[3]/p1[0]]
		cdef vector[double] dx4_cell = boost4_By3(dx4_lab, v3cell)
		t_elapse_cell = dx4_cell[0]
	
		# Sample channel in cell, constrained by dt_max
		channel, dt_cell = self.hqsample.sample_channel(E1_cell, Temp, t_elapse_cell)

		# Boost evolution time back to lab frame
		cdef double dt_lab = boost4_By3([dt_cell, 0., 0., 0.], rv3cell)[0]
		if channel < 0:
			return channel, dt_lab, p1
	
		# Sample E2_cell, s in cell
		E2_cell, s = self.hqsample.sample_initial(channel, E1_cell, Temp, t_elapse_cell)
	
		# Imagine rotate p1_cell to align with z-direction, construct p2_cell_align
		p1z_cell_align = sqrt(E1_cell**2 - self.M**2)
		costheta2 = (self.M**2 + 2.*E1_cell*E2_cell - s)/2./p1z_cell_align/E2_cell
		sintheta2 = sqrt(1. - costheta2**2)
		phi2 = (rand()*2.*M_PI)/RAND_MAX
		cosphi2 = cos(phi2)
		sinphi2 = sin(phi2) 
		cdef vector[double] p1_cell_align = [E1_cell, 0.0, 0.0, p1z_cell_align]
		cdef vector[double] p2_cell_align = [E2_cell, E2_cell*sintheta2*cosphi2, \
											 E2_cell*sintheta2*sinphi2, E2_cell*costheta2]


		# Center of mass frame of p1_cell_align and p2_cell_align, and take down orientation of p1_com
		cdef int i=0
		cdef vector[double] Pcom = [0., 0., 0., 0.], v3com = [0., 0., 0.]
		for i in range(4):	
			Pcom[i] = p1_cell_align[i] + p2_cell_align[i]
		for i in range(3):	
			v3com[i] = Pcom[i+1]/Pcom[0]
		cdef vector[double] p1_com = boost4_By3(p1_cell_align, v3com)
		alpha1_com = atan2(p1_com[2], p1_com[1]) + M_PI/2.
		beta1_com = atan2(sqrt(p1_com[1]**2+p1_com[2]**2), p1_com[3])
		gamma1_com = 0.0
	
		# Tranform t_elapse_cell to t_elapse_com
		t_elapse_com = boost4_By3(dx4_cell, v3com)[0]

		# Sample final state momentum in Com frame, with incoming paticles on z-axis
		cdef vector[double] p1_new_com_aligen = self.hqsample.sample_final(channel, s, Temp, t_elapse_com)[0]
	
		# Rotate final states back to original Com frame (not z-axis aligened)
		cdef vector[double] p1_new_com = rotate_ByEuler(p1_new_com_aligen, -gamma1_com, -beta1_com, -alpha1_com)
	
		# boost back to cell frame z-aligned
		cdef vector[double] rv3com = [-v3com[0], -v3com[1], -v3com[2]] 
		cdef vector[double] p1_new_cell_align = boost4_By3(p1_new_com, rv3com)
	
		# rotate back to original cell frame
		cdef vector[double] p1_new_cell = rotate_ByEuler(p1_new_cell_align, -gamma1_cell, -beta1_cell, -alpha1_cell)

		# boost back to lab frame
		cdef vector[double] p1_new = boost4_By3(p1_new_cell, rv3cell)
	
		# return updated momentum of heavy quark
		return channel, dt_lab, p1_new

	def HQ_hist(self):
		x, y, z = [], [], []
		E, px, py, pz = [], [], [], []
		plt.clf()
		cdef vector[particle].iterator it = self.active_HQ.begin()
		A = self.hydro_reader.get_current_frame('Temp')
		#plt.imshow(np.flipud(A), extent = [-15, 15, -15, 15])
		#cb = plt.colorbar()
		#cb.set_label(r'$T$ [GeV]')
		while it != self.active_HQ.end():
			E.append(deref(it).p[0])
			px.append(deref(it).p[1])
			py.append(deref(it).p[2])
			pz.append(deref(it).p[3])
			inc(it)
		E = np.array(E); px = np.array(px); py = np.array(py); pz = np.array(pz)
		corner(np.array([E, px, py, pz]), ranges=np.array([[0,6], [-4,4], [-4,4], [-4,4]]))
		plt.pause(0.02)

	def HQ_xy(self):
		x, y = [], []
		#plt.clf()
		cdef vector[particle].iterator it = self.active_HQ.begin()
		"""if self.mode == 'dynamic':
			A = self.hydro_reader.get_current_frame('Temp')
			XL, XH, YL, YH = self.hydro_reader.boundary()
			plt.imshow(np.flipud(A), extent = [XL, XH, YL, YH])
			cb = plt.colorbar()
			cb.set_label(r'$T$ [GeV]')"""
		while it != self.active_HQ.end():
			x.append(deref(it).x[1])
			y.append(deref(it).x[2])
			inc(it)
		plt.scatter(x, y, s=0.3, alpha=0.3)
		plt.axis([-15,15,-15,15])
		plt.pause(0.02)
















