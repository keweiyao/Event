import sys
sys.path.append('../HQ-Evo/')
sys.path.append('../Medium-Reader')
import HqEvo
import medium
import h5py
from libcpp cimport bool
from libcpp.vector cimport vector
from libc.stdlib cimport rand, RAND_MAX
from libc.math cimport *
import matplotlib.pyplot as plt
from cython.operator cimport dereference as deref, preincrement as inc
import numpy as np

cdef double GeVm1_to_fmc = 0.197

#------------Import c++ utility function for rotation and transformation-----------
cdef extern from "../src/utility.h":
	cdef double product4(const vector[double] & A, const vector[double] & B)
	cdef void print4vec(const vector[double] & A)
	cdef void rotate_Euler(const vector[double] & A, vector[double] & Ap, double alpha, double beta, double gamma)
	cdef void rotate_axis(const vector[double] & A, vector[double] & Ap, double alpha, unsigned int dir)
	cdef void boost_by3(const vector[double] & A, vector[double] & Ap, const vector[double] v)
	cdef void boost_by4(const vector[double] & A, vector[double] & Ap, const vector[double] u)
	cdef void boost_axis(const vector[double] & A, vector[double] & Ap, const double vd, unsigned int dir)
	cdef void go_to_CoM(const vector[double] & Pcom,
			   const vector[double] & A, const vector[double] & B,
			   vector[double] & Ap, vector[double] & Bp);

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
	cdef object hydro_reader, hqsample
	cdef double M
	cdef vector[particle] active_HQ
	cdef vector[particle] frzout_HQ
	cdef double tau0, dtau, tau
	X = []
	Y = []
	
	def __cinit__(self, hydrofile=None, mass=1.3, elastic=True, inelastic=False, table_folder='./tables'):
		self.M = mass
		self.hydro_reader = medium.Medium(hydrofile)
		self.hqsample = HqEvo.HqEvo(mass=mass, elastic=elastic, inelastic=inelastic, table_folder=table_folder)
		self.tau0 = self.hydro_reader.init_tau()
		self.dtau = self.hydro_reader.dtau()
		self.tau = self.tau0
	
	cpdef initialize_HQ(self, NQ=100, XY_table=None, Pweight=None):
		print "Uniform sampling transverse momentum with in circle pt2<"
		print "Zero longitudinal momentum"
		self.active_HQ.resize(NQ)
		cdef double pt, phipt, r, phir, E, pz, free_time
		cdef particle Q
		
		cdef vector[particle].iterator it = self.active_HQ.begin()
		while it != self.active_HQ.end():
			self.X.append([])
			self.Y.append([])
			pt = sqrt(1.0 + (8.*rand())/RAND_MAX)
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
	
	cpdef perform_hydro_step(self):
		self.hydro_reader.load_next()
		self.tau += self.dtau
		cdef double t, x, y, z, T, vx, vy, vz, dt, dt1, dt2
		cdef vector[double] pnew
		pnew.resize(4)
		cdef vector[particle].iterator it = self.active_HQ.begin()
		while it != self.active_HQ.end():
			t, x, y, z = deref(it).x
			tauQ = sqrt(t**2 - z**2)
			while (tauQ < self.tau):
				#Note that this vx and vy are at mid-rapidity and need to be boosted to obtain the solution at forward and backward rapidity
				T, vx, vy = self.hydro_reader.interpF(tauQ, [x, y], ['Temp', 'Vx', 'Vy'])
				vz = z/t
				# boost to get the solution at (t, x, y, z)
				gamma = 1.0/sqrt(1.0-vz*vz)
				vx = vx/gamma
				vy = vy/gamma
				
				dt, pnew = self.update_HQ(deref(it).p, [vx, vy, vz], T)
				dt *= GeVm1_to_fmc
				dt1 = (dt*rand())/RAND_MAX
				dt2 = dt - dt1
				deref(it).x = freestream(deref(it).x, deref(it).p, dt1)
				deref(it).x = freestream(deref(it).x, pnew, dt2)
				deref(it).p = pnew
				t, x, y, z = deref(it).x
				tauQ = sqrt(t**2 - z**2)
			inc(it)

	def plot_xy(self):
		x = []
		y = []
		plt.clf()
		cdef vector[particle].iterator it = self.active_HQ.begin()
		A = self.hydro_reader.get_current_frame('Temp')
		plt.imshow(np.flipud(A), extent = [-15, 15, -15, 15])
		cb = plt.colorbar()
		cb.set_label(r'$T$ [GeV]')
		while it != self.active_HQ.end():
			x.append(deref(it).x[1])
			y.append(deref(it).x[2])
			inc(it)
		
		x = np.array(x)
		y = np.array(y)
		plt.scatter(x, y)
		plt.axis([-15, 15, -15, 15])
		plt.xlabel(r"$x$ [fm]")
		plt.ylabel(r"$y$ [fm]")
		plt.savefig("figs/%1.3f.png"%self.tau)
		plt.pause(0.1)


	
	cpdef update_HQ(self, vector[double] p1, vector[double] v3cell, double Temp):
		# Define local variables
		cdef double E1_cell, alpha1_cell, beta1_cell, gamma1_cell, E2_cell, s
		cdef double alpha_com, beta_com, gamma_com
		cdef double p1z_cell_align, costheta2, sintheta2, cosphi2, sinphi2
		cdef double dt
		cdef int channel

		# Boost to cell frame and take down orientation of p1_cell
		cdef vector[double] p1_cell = boost4_By3(p1, v3cell)
		E1_cell = p1_cell[0]
		alpha1_cell = atan2(p1_cell[2], p1_cell[1]) + M_PI/2.
		beta1_cell = atan2(sqrt(p1_cell[1]**2+p1_cell[2]**2), p1_cell[3])
		gamma1_cell = 0.0

		# Sample channel in cell, constrained by dt_max
		channel, dt = self.hqsample.sample_channel(E1_cell, Temp)

		if channel < 0:
			return dt, p1
	
		# Sample E2_cell, s in cell
		E2_cell, s = self.hqsample.sample_initial(channel, E1_cell, Temp)
	
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

		# Sample final state momentum in Com frame, with incoming paticles on z-axis
		cdef vector[double] p1_new_com_aligen = self.hqsample.sample_final(channel, s, Temp)[0]
	
		# Rotate final states back to original Com frame (not z-axis aligened)
		cdef vector[double] p1_new_com = rotate_ByEuler(p1_new_com_aligen, -gamma1_com, -beta1_com, -alpha1_com)
	
		# boost back to cell frame z-aligned
		cdef vector[double] rv3com = [-v3com[0], -v3com[1], -v3com[2]] 
		cdef vector[double] p1_new_cell_align = boost4_By3(p1_new_com, rv3com)
	
		# rotate back to original cell frame
		cdef vector[double] p1_new_cell = rotate_ByEuler(p1_new_cell_align, -gamma1_cell, -beta1_cell, -alpha1_cell)

		# boost back to lab frame
		cdef vector[double] rv3cell = [-v3cell[0], -v3cell[1], -v3cell[2]]
		cdef vector[double] p1_new = boost4_By3(p1_new_cell, rv3cell)
	
		# return updated momentum of heavy quark
		return dt, p1_new
















