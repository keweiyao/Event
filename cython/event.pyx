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

def dfdP(p, T, M):
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

#-------------C/Python Wrapper functions--------------------------------------------
cpdef boost4_ByCoM(vector[double] & A, vector[double] & B):
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
cdef freestream(vector[double] & x, vector[double] & p, double & dt):
	return [x[0]+dt, x[1]+p[1]/p[0]*dt, x[2]+p[2]/p[0]*dt, x[3]+p[3]/p[0]*dt]

cdef class event:
	cdef object hydro_reader, hqsample, mode
	cdef double M
	cdef vector[particle] active_HQ
	cdef vector[particle] frzout_HQ
	cdef double tau0, dtau, tau
	X = []
	Y = []
	C = []
	
	def __cinit__(self, mode="dynamic", hydrofile=None, mass=1.3, elastic=True, inelastic=False, detailed_balance=False, table_folder='./tables'):
		self.mode = mode
		self.M = mass
		self.hydro_reader = medium.Medium(mode=mode, hydrofilename=hydrofile)
		self.hqsample = HqEvo.HqEvo(mass=mass, elastic=elastic, inelastic=inelastic, detailed_balance=detailed_balance, table_folder=table_folder)
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
			pz = 30.
			E = sqrt(self.M**2 + pt**2 + pz**2)
			
			deref(it).p = [E, pt*cos(phipt), pt*sin(phipt), pz]
			deref(it).x = [0.0, r*cos(phir), r*sin(phir), 0.0]
			if self.mode == 'dynamic':
				free_time = self.tau0/sqrt(1. - (pz/E)**2)
				# free streaming to hydro starting time
				deref(it).x = freestream(deref(it).x, deref(it).p, free_time)
				deref(it).t_last = 0.
			else:
				# randomize t_last
				deref(it).t_last = -(5.*rand())/RAND_MAX
			inc(it)

	cpdef perform_hydro_step(self, StaticPropertyDictionary=None):

		status = True
		if self.mode == 'dynamic':
			status = self.hydro_reader.load_next()
		if self.mode == 'static':
			status = self.hydro_reader.load_next(StaticPropertyDictionary=StaticPropertyDictionary)
		self.tau += self.dtau
		cdef double t, x, y, z, tauQ2
		cdef vector[particle].iterator it = self.active_HQ.begin()
		while it != self.active_HQ.end():
			if self.mode == 'static':
				while (deref(it).x[0] < self.tau):
					self.perform_HQ_step(deref(it))
			if self.mode == 'dynamic':
				t, x, y, z = deref(it).x
				tauQ = sqrt(t**2 - z**2)
				while (tauQ < self.tau):
					self.perform_HQ_step(deref(it))
					t, x, y, z = deref(it).x
					tauQ = sqrt(t**2 - z**2)
			inc(it)
		return status

	cdef perform_HQ_step(self, particle & it):
		cdef double t, x, y, z, t_elapse_lab, tauQ
		cdef double T, vx, vy, vz
		cdef int channel

		t, x, y, z = it.x
		t_elapse_lab = t - it.t_last
		
		tauQ = sqrt(t**2 - z**2)
		T, vx, vy, vz = self.hydro_reader.interpF(tauQ, [x, y, z, t], ['Temp', 'Vx', 'Vy', 'Vz'])	
		# Note the change of units from GeV-1 (fm/c) to fm/c (GeV-1)
		t_elapse_lab *= GeVm1_to_fmc
		channel, dt, pnew = self.update_HQ(it.p, [vx, vy, vz], T, t_elapse_lab)
		dt *= GeVm1_to_fmc
		self.C.append(channel)
		if channel < 0:
			it.x = freestream(it.x, it.p, dt)
		else:
			dt1 = (dt*rand())/RAND_MAX
			dt2 = dt - dt1
			it.x = freestream(it.x, it.p, dt1)
			it.x = freestream(it.x, pnew, dt2)
			it.p = pnew
			it.t_last = t + dt1
		return

	cdef update_HQ(self, vector[double] & p1, vector[double] & v3cell, double & Temp, double & t_elapse_lab):
		# Define local variables
		cdef double E1_cell
		cdef double t_elapse_cell, t_elapse_com
		cdef double dt, dt_cell
		cdef int channel
		
		# Boost to cell frame and take down orientation of p1_cell
		cdef vector[double] p1_cell = boost4_By3(p1, v3cell)
		E1_cell = p1_cell[0]

		# Boost time separation to cell frame
		dx4_lab = [t_elapse_lab, t_elapse_lab*p1[1]/p1[0],
				   t_elapse_lab*p1[2]/p1[0], t_elapse_lab*p1[3]/p1[0]]
		cdef vector[double] dx4_cell = boost4_By3(dx4_lab, v3cell)
		t_elapse_cell = dx4_cell[0]
	
		# Sample channel in cell, constrained by dt_max
		channel, dt_cell = self.hqsample.sample_channel(E1_cell, Temp, t_elapse_cell)
		
		# Boost evolution time back to lab frame
		cdef double dt_lab = boost4_By3([dt_cell, 0., 0., 0.], [-v3cell[0], -v3cell[1], -v3cell[2]])[0]
		if channel < 0:
			return channel, dt_lab, p1
	
		# Sample initial state and return initial state particle four vectors 
		# Imagine rotate p1_cell to align with z-direction, construct p2_cell_align, ...
		cdef vector[ vector[double] ] Initial_States = self.hqsample.sample_initial(channel, E1_cell, Temp, t_elapse_cell)
		p1_cell_align = Initial_States[0]

		# Center of mass frame of p1_cell_align and other particles, and take down orientation of p1_com
		cdef size_t i=0
		cdef vector[double] Pcom = [0., 0., 0., 0.], v3com = [0., 0., 0.]
		for pp in Initial_States:
			for i in range(4):
				Pcom[i] += pp[i]
		cdef double s = product4(Pcom, Pcom)
		
		for i in range(3):	
			v3com[i] = Pcom[i+1]/Pcom[0]
		
		cdef vector[double] p1_com = boost4_By3(p1_cell_align, v3com)
	
		# Tranform t_elapse_cell to t_elapse_com
		t_elapse_com = boost4_By3(dx4_cell, v3com)[0]
		
		# Sample final state momentum in Com frame, with incoming paticles on z-axis
		cdef vector[double] p1_new_com_aligen = self.hqsample.sample_final(channel, s, Temp, t_elapse_com)[0]
		# Rotate final states back to original Com frame (not z-axis aligened)
		cdef vector[double] p1_new_com = rotate_back_from_D(p1_new_com_aligen, p1_com[1], p1_com[2], p1_com[3])
	
		# boost back to cell frame z-aligned
		cdef vector[double] p1_new_cell_align = boost4_By3(p1_new_com, [-v3com[0], -v3com[1], -v3com[2]] )
		
		# rotate back to original cell frame
		cdef vector[double] p1_new_cell = rotate_back_from_D(p1_new_cell_align, p1_cell[1], p1_cell[2], p1_cell[3])

		# boost back to lab frame
		cdef vector[double] p1_new = boost4_By3(p1_new_cell, [-v3cell[0], -v3cell[1], -v3cell[2]])
	
		# return updated momentum of heavy quark
		return channel, dt_lab, p1_new

	def HQ_hist(self):
		x, y, z = [], [], []
		E, px, py, pz = [], [], [], []
		plt.clf()
		cdef vector[particle].iterator it = self.active_HQ.begin()
		while it != self.active_HQ.end():
			E.append(deref(it).p[0])
			px.append(deref(it).p[1])
			py.append(deref(it).p[2])
			pz.append(deref(it).p[3])
			inc(it)
		E = np.array(E); px = np.array(px); py = np.array(py); pz = np.array(pz)
		corner(np.array([E, px, py, pz]), ranges=np.array([[0,32], [-4,4], [-4,4], [-4,32]]))
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
















