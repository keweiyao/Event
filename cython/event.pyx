import numpy as np
import sys
import matplotlib.pyplot as plt
from libcpp cimport bool
from libcpp.vector cimport vector
from libc.stdlib cimport rand, RAND_MAX
from libc.math cimport *
from cython.operator cimport dereference as deref, preincrement as inc
cimport numpy as np

sys.path.append('../HQ-Evo/')
sys.path.append('../Medium-Reader')
import HqEvo
import medium

from event cimport *

cdef double GeVm1_to_fmc = 0.197

cdef double default_df_dpt2dy(pT, y):
	return 1.0/(1.0+pT**4)*1.0/cosh(y)

cdef freestream(vector[double] & x, vector[double] & p, double dt):
	x[0] += dt
	x[1] += p[1]/p[0]*dt
	x[2] += p[2]/p[0]*dt
	x[3] += p[3]/p[0]*dt

#-----------Particle data struct------------------------------------
cdef extern from "../src/utility.h":
	cdef struct particle:
		vector[double] p
		vector[double] x
		float t_last,t_last2
		float Nc, Nc2
		float weight

#-----------Production vertex sampler class-------------------------
class XY_sampler:
	__slots__ = ['Taa']
	def __init__(self, Taa):
		self.Taa = Taa;
	def sample_xy(self):

#-----------Event Class---------------------------------------------
cdef class event:
	cdef object hydro_reader, hqsample, mode, C
	cdef double M
	cdef vector[particle] active_HQ
	cdef vector[particle] frzout_HQ
	cdef double tau0, dtau, tau
	
	def __cinit__(self, mode="dynamic", hydrofile=None, static_dt=0.5, mass=1.3, elastic=True, inelastic=False, detailed_balance=False, table_folder='./tables'):
		self.mode = mode
		self.M = mass
		self.hydro_reader = medium.Medium(mode=mode, hydrofilename=hydrofile, static_dt=static_dt)
		self.hqsample = HqEvo.HqEvo(mass=mass, elastic=elastic, inelastic=inelastic, detailed_balance=detailed_balance, table_folder=table_folder)
		self.tau0 = self.hydro_reader.init_tau()
		self.dtau = self.hydro_reader.dtau()
		self.tau = self.tau0

	def sys_time(self):
		return self.tau

	cpdef initialize_HQ(self, NQ=1000, df_dpt2dy=None):
		self.active_HQ.clear()
		self.active_HQ.resize(NQ)
		cdef double x, y, z
		cdef double pT, phipt, rapidity, mT, pTmax = 70., freetime
		cdef double p, cospz, sinpz, pmax=50.
		cdef vector[particle].iterator it

		if self.mode == 'dynamic':
			print "Initialize for dynamic medium"
			print "Uniform sampling HQ pT and rapidity"
			print "0. < pt < 70. [GeV], -1. < y < 1."
			print "Each heavy quark is weighted by df/dpt^2/dy"
			if df_dpt2dy==None:
				print "Notice: user define df/dpt^2/dy spectra not available, use defaule weight"
				df_dpt2dy = default_df_dpt2dy
			print "Heavy quark will be free streamed to the starting time of hydro"
			it = self.active_HQ.begin()
			while it != self.active_HQ.end():
				pT = sqrt(rand()*1./RAND_MAX)*pTmax
				mT = sqrt(pT**2 + self.M**2)
				phipt = rand()*2.*M_PI/RAND_MAX
				rapidity = rand()*2./RAND_MAX-1.
				deref(it).p.resize(4)
				deref(it).x.resize(4)
				deref(it).p = [mT*cosh(rapidity), pT*cos(phipt), pT*sin(phipt), mT*sinh(rapidity)]
				deref(it).t_last = 0.; deref(it).t_last2 = 0.
				deref(it).Nc = 0.; deref(it).Nc2 = 0.
				deref(it).weight = df_dpt2dy(pT, y)
		
				# free streaming to hydro starting time
				free_time = self.tau0/sqrt(1. - (deref(it).p[3]/deref(it).p[0])**2)
				x, y = 1., 1.
				z = 0.
				deref(it).x = [0.0, x, y, z]
				freestream(deref(it).x, deref(it).p, free_time)
				inc(it)

		if self.mode == 'static':
			print "Initialize for static medium"
			print "Uniform sampling HQ momentum"
			print "0. < |p| < 50. [GeV]"
			it = self.active_HQ.begin()
			while it != self.active_HQ.end():
				p = pow(rand()*1./RAND_MAX, 1./3.)*pmax
				phipt = rand()*2.*M_PI/RAND_MAX
				cospz = 0.999999*(rand()*2./RAND_MAX-1.)
				sinpz = sqrt(1.-cospz**2)
				deref(it).p.resize(4)
				deref(it).x.resize(4)
				deref(it).p = [sqrt(p**2+self.M**2), p*sinpz*cos(phipt), p*sinpz*sin(phipt), p*cospz]
				deref(it).t_last = 0.; deref(it).t_last2 = 0.
				deref(it).Nc = 0.; deref(it).Nc2 = 0.
				deref(it).weight = 1.0
				x = rand()*10./RAND_MAX - 5.
				y = rand()*10./RAND_MAX - 5.
				z = rand()*10./RAND_MAX - 5.
				deref(it).x = [0.1, x, y, z]
				inc(it)	
	
	cpdef perform_hydro_step(self, StaticPropertyDictionary=None):
		status = True
		if self.mode == 'dynamic':
			status = self.hydro_reader.load_next()
		if self.mode == 'static':
			if StaticPropertyDictionary==None:
				raise ValueError("static meidum property not defined")
			status = self.hydro_reader.load_next(StaticPropertyDictionary=StaticPropertyDictionary)
		self.tau += self.dtau
		cdef double t, x, y, z, tauQ2
		cdef vector[particle].iterator it
		if self.mode == 'static':
			it = self.active_HQ.begin()
			while it != self.active_HQ.end():
				while deref(it).x[0] <= self.tau:
					self.perform_HQ_step(deref(it))
				inc(it)
		if self.mode == 'dynamic':
			it = self.active_HQ.begin()
			while it != self.active_HQ.end():
				t, x, y, z = deref(it).x
				tauQ = sqrt(t**2 - z**2)
				while tauQ <= self.tau:
					self.perform_HQ_step(deref(it))
					t, x, y, z = deref(it).x
					tauQ = sqrt(t**2 - z**2)
				inc(it)
		return status

	cpdef perform_HQ_step(self, particle & it):
		cdef double t, x, y, z, t_elapse_lab, tauQ, dt
		cdef double T, vx, vy, vz
		cdef vector[double] pnew
		cdef int channel
		t, x, y, z = it.x
		t_elapse_lab = (t - it.t_last)/(it.Nc + 1.)
		t_elapse_lab2 = (t - it.t_last2)/(it.Nc2 + 1.)
		
		tauQ = sqrt(t**2 - z**2)
		T, vx, vy, vz	\
		= self.hydro_reader.interpF(tauQ, [x, y, z, t], 
					['Temp', 'Vx', 'Vy', 'Vz'])
		cdef double v2 = vx**2 + vy**2 + vz**2
		if v2 > 0.99999:
			vx /= v2*1.0001; vy /= v2*1.0001; vz /= v2*1.0001

		# Note the change of units from GeV-1 (fm/c) to fm/c (GeV-1)
		t_elapse_lab /= GeVm1_to_fmc
		t_elapse_lab2 /= GeVm1_to_fmc

		channel, dt, pnew= self.update_HQ(it.p, [vx, vy, vz], T, t_elapse_lab, t_elapse_lab2)		
		dt *= GeVm1_to_fmc
		if channel < 0:
			freestream(it.x, it.p, dt)
		else:
			freestream(it.x, it.p, dt)
			it.p = pnew
			if channel == 2 or channel == 3:	
				it.Nc = 0.
				it.t_last = t + dt
			if channel == 4 or channel == 5:	
				it.Nc2 = 0.	
				it.t_last2 = t + dt
			it.Nc += 1.
			it.Nc2 += 1.
		return

	cpdef update_HQ(self, vector[double] & p1_lab, vector[double] & v3cell, double & Temp, double & mean_dt23_lab, double & mean_dt32_lab):
		# Define local variables
		cdef double s, L1, L2, Lk, x2, xk, a1=0.6, a2=0.6
		cdef double dt_lab, dt_cell
		cdef int channel
		cdef vector[double] p1_cell, p1_cell_Z, p1_com, \
							p1_com_Z_new, p1_com_new, \
							p1_cell_Z_new, p1_cell_new,\
							p1_new,	fs,				\
							dx23_cell, dx32_cell, dx23_com, \
							v3com 
		cdef vector[ vector[double] ] IS, FS # initial states

		# Boost p1 (lab) to p1 (cell)
		p1_cell = boost4_By3(p1_lab, v3cell)

		# displacement vector within dx23_lab and dx32_lab seen from cell frame
		dx23_cell = [pi/p1_lab[0]*mean_dt23_lab for pi in p1_cell]
		dx32_cell = [pi/p1_lab[0]*mean_dt32_lab for pi in p1_cell]

		# Sample channel in cell, return channl index and evolution time seen from cell
		channel, dt_cell = self.hqsample.sample_channel(p1_cell[0], Temp, dx23_cell[0], dx32_cell[0])

		# Boost evolution time back to lab frame
		dt_lab = p1_lab[0]/p1_cell[0]*dt_cell

		# If not scattered, return channel=-1, evolution time in Lab frame and origin al p1(lab)
		if channel < 0:
			return channel, dt_lab, p1_lab
		# Sample initial state and return initial state particle four vectors 
		# Imagine rotate p1_cell to align with z-direction, construct p2_cell_align, ...
		IS = self.hqsample.sample_initial(channel, p1_cell[0], Temp, dx23_cell[0], dx32_cell[0])
		p1_cell_Z = IS[0]
		# Center of mass frame of p1_cell_align and other particles, and take down orientation of p1_com
		cdef size_t i=0
		Pcom = [0., 0., 0., 0.]
		for pp in IS:
			for i in range(4):
				Pcom[i] += pp[i]
		s = product4(Pcom, Pcom)
		
		v3com = [Pcom[i+1]/Pcom[0] for i in range(3)]

		p1_com = boost4_By3(p1_cell_Z, v3com)
		if channel >= 4: # for 3 -> 2 kinetics
			L1 = sqrt(p1_com[0]**2 - self.M**2)
			L2 = boost4_By3(IS[1], v3com)[0]
			Lk = boost4_By3(IS[2], v3com)[0]
			x2 = L2/(L1+L2+Lk)
			xk = Lk/(L1+L2+Lk)
			a1 = x2 + xk; 
			a2 = (x2 - xk)/(1. - a1)
		
		dx23_com = [pi/p1_cell_Z[0]*dx23_cell[0] for pi in p1_com]

		# Sample final state momentum in Com frame, with incoming paticles on z-axis
		FS = self.hqsample.sample_final(channel, s, Temp, dx23_com[0], a1, a2)

		p1_com_Z_new = FS[0]
		# Rotate final states back to original Com frame (not z-axis aligened)
		p1_com_new = rotate_back_from_D(p1_com_Z_new, p1_com[1], p1_com[2], p1_com[3])	
		# boost back to cell frame z-aligned
		p1_cell_Z_new = boost4_By3(p1_com_new, [-v3com[0], -v3com[1], -v3com[2]])	
		# rotate back to original cell frame
		p1_cell_new = rotate_back_from_D(p1_cell_Z_new, p1_cell[1], p1_cell[2], p1_cell[3])
		# boost back to lab frame
		p1_new = boost4_By3(p1_cell_new, [-v3cell[0], -v3cell[1], -v3cell[2]])
		# return updated momentum of heavy quark

		return channel, dt_lab, p1_new

	cpdef HQ_hist(self):
		cdef vector[particle].iterator it = self.active_HQ.begin()
		cdef vector[ vector[double] ] datap, datax
		datap.clear()
		datax.clear()
		while it != self.active_HQ.end():
			datap.push_back(deref(it).p)
			datax.push_back(deref(it).x)
			inc(it)
		return np.array(datap), np.array(datax)

	cpdef reset_HQ_energy(self, E0=10.):
		cdef vector[particle].iterator it = self.active_HQ.begin()
		cdef double pabs_new, pabs_old, px, py, pz
		while it != self.active_HQ.end():
			ratio = sqrt( (E0**2 - self.M**2)/(deref(it).p[0]**2 - self.M**2) )
			px = deref(it).p[1]*ratio
			py = deref(it).p[2]*ratio 
			pz = deref(it).p[3]*ratio 
			deref(it).p = [E0, px, py, pz]
			inc(it)















