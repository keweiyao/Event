# cython: c_string_type=str, c_string_encoding=ascii
from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.string cimport string
from libc.stdlib cimport rand, RAND_MAX
from libc.math cimport *
from cython.operator cimport dereference as deref, preincrement as inc
from cpython.exc cimport PyErr_CheckSignals
import numpy as np
cimport numpy as np

import HqEvo
import HqLGV
import medium

cdef double GeVm1_to_fmc = 0.197
cdef double little_below_one = 1. - 1e-6
cdef double little_above_one = 1. + 1e-6		

#-----------Particle data struct------------------------------------
cdef extern from "../src/utility.h":
	cdef struct particle:
		bool freezeout
		vector[double] p
		vector[double] x
		double t_last_23, t_last_32
		vector[double] initp
		vector[double] vcell
		double Tf, s1, s2
		int pid

#-----------Production vertex sampler class-------------------------
cdef class XY_sampler:
	cdef np.ndarray Taa, IntTaa
	cdef double dxy, b
	cdef size_t Nx, Ny
	def __cinit__(self, Taa, dxy, b):
		self.b = b
		self.dxy = dxy
		self.Nx, self.Ny = Taa.shape
		self.Taa = Taa.reshape(-1)
		self.IntTaa = np.zeros_like(self.Taa, dtype=np.double)
		cdef double tot = self.Taa.sum()
		cdef int i
		cdef double dT
		for i, dT in enumerate(self.Taa):
			self.IntTaa[i] = self.IntTaa[i-1] + dT/tot
	cpdef sample_xy(self):
		cdef double r = np.random.rand()
		cdef int index = np.searchsorted(self.IntTaa, r)
		cdef double nx = np.floor((index-1.)/self.Ny)
		cdef double ny = index - 1 - nx*self.Ny
		cdef double x, y, s1, s2
		nx += np.random.rand()
		ny += np.random.rand()
		# to be examined
		x, y = (nx - self.Nx/2.)*self.dxy, (ny - self.Ny/2.)*self.dxy
		s1 = sqrt(y**2 + (x+self.b/2.)**2)
		s2 = sqrt(y**2 + (x-self.b/2.)**2)
		return x, y, s1, s2

#-----------Event Class---------------------------------------------
cdef void freestream(vector[particle].iterator it, double dt):
	cdef double da = dt/deref(it).p[0]
	deref(it).x[0] = deref(it).x[0] + dt
	deref(it).x[1] = deref(it).x[1] + deref(it).p[1]*da
	deref(it).x[2] = deref(it).x[2] + deref(it).p[2]*da
	deref(it).x[3] = deref(it).x[3] + deref(it).p[3]*da

cdef double proper_time(vector[double] x):
	return sqrt(x[0]**2 - x[3]**2)

cdef class event:
	cdef object hydro_reader, lbt, lgv, fh
	cdef str mode, transport
	cdef double M, Tc
	cdef vector[particle] HQ_list
	cdef vector[particle] frzout_HQ
	cdef double tau0, dtau, tau

	def __cinit__(self, medium_flags, physics_flags,
					table_folder='./tables', refresh_table=False, seed=None):
		if seed != None:
			np.random.seed(seed)
		# medium
		self.fh = open("Histroy.dat", 'w')
		self.mode = medium_flags['type']
		self.hydro_reader = medium.Medium(medium_flags=medium_flags)
		self.tau0 = self.hydro_reader.init_tau()
		self.dtau = self.hydro_reader.dtau()
		self.tau = self.tau0

		# transport
		self.transport = physics_flags['transport']['name']
		self.Tc = physics_flags['Tc']
		self.M = physics_flags['mass']

		if self.transport == "LBT":
			print("LBT mode")
			self.lbt = HqEvo.HqEvo(options=physics_flags,
									table_folder=table_folder,
									refresh_table=refresh_table)
		elif self.transport == "LGV":
			print("LGV mode")
			self.lgv = HqLGV.HqLGV(options=physics_flags)
		elif self.transport == "Hybrid":
			print("Hybrid mode")
			self.lgv = HqLGV.HqLGV(options=physics_flags)
			self.lbt = HqEvo.HqEvo(options=physics_flags,
									table_folder=table_folder,
									refresh_table=refresh_table)

	# Return the current time of the evolution.
	def sys_time(self) :
		return self.tau

	# Initilization
	cpdef initialize_HQ(self, NQ, init_flags):
		self.HQ_list.clear()
		self.HQ_list.resize(NQ)
		cdef double x, y, z, s1, s2
		# for A+B:
		cdef double pT, phipt, rapidity, mT, t0
		cdef double ymin, ymax, pTmax, pTmin
		# for box:
		cdef double p, cospz, sinpz
		cdef double pmax, L
		cdef vector[particle].iterator it

		if init_flags['type'] == 'A+B':
			print("Initialize for dynamic medium")
			HQ_xy_sampler = XY_sampler(init_flags['TAB'],
									   init_flags['dxy'],
									   init_flags['b'])
			pTmax = init_flags['pTmax']
			pTmin = init_flags['pTmin']
			ymax = init_flags['ymax']

			print("Heavy quarks are freestreamed to {} fm/c".format(self.tau0))
			it = self.HQ_list.begin()
			X = []
			Y = []
			while it != self.HQ_list.end():
				# Initialize momentum and spatial space
				pT = np.random.rand()*(pTmax-pTmin) + pTmin
				mT = sqrt(pT**2 + self.M**2)
				phipt = np.random.rand()*2.*M_PI
				rapidity = np.random.rand()*ymax*2. - ymax
				pcharm = [mT*cosh(rapidity), pT*cos(phipt), \
							   pT*sin(phipt), mT*sinh(rapidity)]
				panticharm = [mT*cosh(rapidity), -pT*cos(phipt), \
							   -pT*sin(phipt), -mT*sinh(rapidity)]
				x, y, s1, s2 = HQ_xy_sampler.sample_xy()
				r0 = [0.0, x, y, 0.0]
				t0 = self.tau0/sqrt(1. - (pcharm[3]/pcharm[0])**2)
				X.append(x)
				Y.append(y)
				# charm:
				# Initialize positional space at tau = 0+
				deref(it).p = pcharm
				deref(it).x = r0
				deref(it).s1 = s1
				deref(it).s2 = s2
				# free streaming to hydro starting time tau = tau0
				freestream(it, t0)
				# set last interaction vertex (assumed to be hydro start time)
				deref(it).t_last_23 = t0
				deref(it).t_last_32 = t0
				# initialize others
				deref(it).freezeout = False
				deref(it).initp = deref(it).p
				deref(it).vcell = [0., 0., 0.]
				deref(it).Tf = 0.
				deref(it).pid = 4
				inc(it)
				# anti-charm: only flip p^mu and pid
				# Initialize positional space at tau = 0+
				deref(it).p = panticharm
				deref(it).x = r0
				deref(it).s1 = s1
				deref(it).s2 = s2
				# free streaming to hydro starting time tau = tau0
				freestream(it, t0)
				# set last interaction vertex (assumed to be hydro start time)
				deref(it).t_last_23 = t0
				deref(it).t_last_32 = t0
				# initialize others
				deref(it).freezeout = False
				deref(it).initp = deref(it).p
				deref(it).vcell = [0., 0., 0.]
				deref(it).Tf = 0.
				deref(it).pid = 4
				inc(it)
			# check the variance of the sampling
			stdx, stdy = np.std(X), np.std(Y)
			print("std(x,y) = {:1.3f}, {:1.3f} [fm]".format(stdx, stdy) )

		if init_flags['type'] == 'box':
			print "Initialize for box simulation"
			pmax = init_flags['pmax']
			L = init_flags['L']
			it = self.HQ_list.begin()

			while it != self.HQ_list.end():
				p = np.random.uniform(0, pmax)
				phipt = np.random.uniform(0, 2.*np.pi)
				cospz = little_below_one*np.random.uniform(-1., 1.)
				sinpz = sqrt(1.-cospz**2)
				deref(it).p.resize(4)
				deref(it).x.resize(4)
				deref(it).p = [sqrt(p**2+self.M**2), p*sinpz*cos(phipt), \
								p*sinpz*sin(phipt), p*cospz]
				deref(it).t_last_23 = 0.0; deref(it).t_last_32 = 0.0
				deref(it).freezeout = False
				deref(it).initp = [sqrt(p**2+self.M**2), p*sinpz*cos(phipt), \
									p*sinpz*sin(phipt), p*cospz]
				deref(it).vcell = [0., 0., 0.]
				deref(it).Tf = 0.
				deref(it).pid = 4*np.random.choice([1, -1])
				x,y,z = np.random.uniform(-L,L,3)
				deref(it).x = [0.0, x, y, z]
				inc(it)

	cpdef bool perform_hydro_step(self, StaticProperty=None):
		PyErr_CheckSignals()
		if self.mode == 'dynamic':
			self.hydro_reader.load_next()
		elif self.mode == 'static':
			if StaticProperty==None:
				raise ValueError("static meidum property not defined")
			else:
				self.hydro_reader.load_next(StaticProperty=StaticProperty)
		else:
			raise ValueError("medium mode not defined")
		status = self.hydro_reader.hydro_status()

		self.tau += self.dtau
		cdef double t, x, y, z, tauQ
		cdef vector[particle].iterator it
		if self.mode == 'static':
			it = self.HQ_list.begin()
			while it != self.HQ_list.end():
				while deref(it).x[0] <= self.tau:
					self.perform_HQ_step(it)
				inc(it)
		elif self.mode == 'dynamic':
			it = self.HQ_list.begin()
			while it != self.HQ_list.end():
				if not deref(it).freezeout:
					tauQ = proper_time( deref(it).x )
					while tauQ <= self.tau and not deref(it).freezeout:
						self.perform_HQ_step(it)
						tauQ = proper_time( deref(it).x )
				inc(it)
		else:
			raise ValueError("Mode not implemented")
		return status

	cdef perform_HQ_step(self, vector[particle].iterator it):
		PyErr_CheckSignals()
		cdef double t, tauQ, T, vabs2, scale
		cdef double dt23_lab, dt32_lab, dt23_cell, dt32_cell
		cdef double dtHQ
		cdef vector[double] vcell
		cdef int channel

		###############################################################
		###############################################################
		# Get the cell temperature and velocity for this heavy quark, #
		# ensure |v| < 1. if T<Tc, label it as "freezout=True"        #
		###############################################################
		t = deref(it).x[0]
		tauQ = proper_time( deref(it).x )
		vcell.resize(3)
		T, vcell[0], vcell[1], vcell[2] = \
			self.hydro_reader.interpF(tauQ, deref(it).x,
									['Temp', 'Vx', 'Vy', 'Vz'])
		vabs2 = vcell[0]**2 + vcell[1]**2 + vcell[2]**2
		if vabs2 >= little_below_one**2:
			scale = 1.0/sqrt(vabs2)/little_above_one
			vcell[0] *= scale
			vcell[1] *= scale
			vcell[2] *= scale

		if T <= self.Tc:
			deref(it).freezeout = True
			deref(it).Tf = T
			deref(it).vcell = vcell
			return

		###############################################################
		###############################################################
		# update heavy quark status, returns which scatterig channel, #
		# how long has it should be evolved and the new momentum      #
		###############################################################
		# variables
		cdef vector[double] p1_cell, p1_cell_star, p1_cell_new, p1_new, p1,\
							k_gluon, k_gluon_cell = [0,0,0,0]
		p1 = deref(it).p # this is the initial momentum in lab frame

		# step<+0>: convert time[fm/c] to time[GeV^-1]
		# calculate the time elapse from the last emission / absorption
		dt23_lab = (t - deref(it).t_last_23)/GeVm1_to_fmc
		dt32_lab = (t - deref(it).t_last_32)/GeVm1_to_fmc

		# step<+1>: Boost from p1(lab) to p1(cell) by vcell
		boost4_By3(p1_cell, p1, vcell)
		dt23_cell = p1_cell[0]/p1[0]*dt23_lab
		dt32_cell = p1_cell[0]/p1[0]*dt32_lab

		# step<2>: evolution mode:
		if self.transport == "LGV": # Langevin mode
			dt_cell = 0.1 # GeV^-1, approxamtely 0.02 fm/c
			p1_cell_new = self.lgv_step(p1_cell, T, dt_cell)
		elif self.transport == "LBT": # Boltzmann mode
			channel, dt_cell = \
				self.lbt.sample_channel(p1_cell[0], T, dt23_cell, dt32_cell)
			if channel < 0:
				p1_cell_new, k_gluon_cell = p1_cell, [0,0,0,0]
			else:
				p1_cell_new, k_gluon_cell = self.lbt_step(channel, dt_cell,
											p1_cell, T, dt23_cell, dt32_cell)
		elif self.transport == "Hybrid": # Hybrid mode
			# Estimate a proper dt_cell using LBT scattering rates,
			channel, dt_cell = \
				self.lbt.sample_channel(p1_cell[0], T, dt23_cell, dt32_cell)
			# determine which steps comes first
			if rand()%2 == 0:	# lbt followed by lgv
				# since lbt comes first, we can use the channel sampled above
				if channel < 0:
					p1_cell_star, k_gluon_cell = p1_cell, [0,0,0,0]
				else:
					p1_cell_star, k_gluon_cell = self.lbt_step(channel, dt_cell,
											p1_cell, T, dt23_cell, dt32_cell)
				# then, perform lgv step
				p1_cell_new = self.lgv_step(p1_cell_star, T, dt_cell)
			else: # lbt followed by lgv
				# first, perform lgv step
				p1_cell_star = self.lgv_step(p1_cell, T, dt_cell)
				# since momentum is changed, we need to resample the channel,
				# but dt_cell is unchanged
				channel, _ = \
					self.lbt.sample_channel(p1_cell_star[0], T,
											dt23_cell, dt32_cell)
				if channel < 0:
					p1_cell_new, k_gluon_cell = p1_cell, [0,0,0,0] 
				else:
					p1_cell_new, k_gluon_cell = self.lbt_step(channel, dt_cell,
										p1_cell_star, T, dt23_cell, dt32_cell)
		else:
			raise ValueError("Transport mode has to be LGV, LBT or Hybrid")

		# step<-1>: Boost back from p1_new(cell) to p1_new by vcell
		boost4_By3_back(p1_new, p1_cell_new, vcell)
		if k_gluon_cell[0] > 1e-9:
			boost4_By3_back(k_gluon, k_gluon_cell, vcell)
		else:
			k_gluon = [0,0,0,0]
		dtHQ = p1[0]/p1_cell[0]*dt_cell

		# step<-0>: convert time[GeV^-1] to time[fm/c]
		dtHQ *= GeVm1_to_fmc

		#if channel >= 0:
		#self.fh.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(dtHQ, channel, *p1_new, *k_gluon) )

		# Update last emission/absorption time
		if channel == 2 or channel == 3:
			deref(it).t_last_23 = t + dtHQ
		elif channel == 4 or channel == 5:
			deref(it).t_last_32 = t + dtHQ
		else:
			pass

		# freestream and update momentum
		freestream(it, dtHQ)
		deref(it).p = p1_new

	cdef vector[double] lgv_step(self, vector[double] p1_cell,
												double T, double dt_lrf) :
		cdef vector[double] p1_cell_Z_new, p1_cell_new
		p1_cell_Z_new = self.lgv.update_by_Langevin(p1_cell[0], T, dt_lrf)
		rotate_back_from_D(p1_cell_new, p1_cell_Z_new,
							p1_cell[1], p1_cell[2], p1_cell[3])
		return p1_cell_new

	# this function needs to be cdef so that we can use the reference copy
	cdef (vector[double], vector[double]) lbt_step(self, int channel,
			double dt_cell, vector[double] p1_cell,	double T,
			double dt23_cell, double dt32_cell) :

		# Define local variables
		cdef double s, dt23_com
		cdef size_t i=0
		cdef vector[double] p1_cell_Z, Pcom, vcom, p1_com, pbuffer,\
			 				p1_com_Z_new, p1_com_new, \
			 				p1_cell_Z_new, p1_cell_new, k_gluon
		cdef double L1, L2, Lk, x2, xk, a1=0.6, a2=0.6

		# For different channel:
		if channel < 0:
			raise ValueError("Makesure lbt_step(...) gets a scattering channel")
		else:
			# Sample initial state particles' four vectors,
			# in the frame (C-Z) where p1_cell is rotated to align on z-axis
			self.lbt.sample_initial(channel, p1_cell[0],
									T, dt23_cell, dt32_cell)

			# Get the momentum of the initial state HQ in the (C-Z) frame
			p1_cell_Z = self.lbt.IS[0]

			# Calculate the Mendelstem variable s = (sum IS) dot (sum IS)
			Pcom.resize(4)
			for i in range(4):
				Pcom[i] = 0.
			for pbuffer in self.lbt.IS:
				for i in range(4):
					Pcom[i] += pbuffer[i]
			s = Pcom[0]**2 - Pcom[1]**2 - Pcom[2]**2 - Pcom[3]**2

			# Calculate center of mass velocity of initial states
			vcom.resize(3)
			for i in range(3):
				vcom[i] = Pcom[i+1]/Pcom[0];

			# <+0> Boost from (C-Z) to center of mass frame (COM)
			boost4_By3(p1_com, p1_cell_Z, vcom)
			dt23_com = p1_com[0]/p1_cell_Z[0]*dt23_cell

			# extra work for channel=4,5 for 3 -> 2 kinematics
			if channel in [4,5]:
				L1 = sqrt(p1_com[0]**2 - self.M**2)
				boost4_By3(pbuffer, self.lbt.IS[1], vcom)
				L2 = pbuffer[0]
				boost4_By3(pbuffer, self.lbt.IS[2], vcom)
				Lk = pbuffer[0]
				x2 = L2/(L1+L2+Lk)
				xk = Lk/(L1+L2+Lk)
				a1 = x2 + xk;
				a2 = (x2 - xk)/(1. - a1)

			# In CoM frame where hQ align on z-axis (COM-Z),
			# Sample final state momentum
			self.lbt.sample_final(channel, s, T, dt23_com, a1, a2)

			# Get the momentum of the final state HQ in the (COM-Z) frame
			p1_com_Z_new = self.lbt.FS[0]

			# Rotate back to from (Com-Z) frame to (COM) frame
			rotate_back_from_D(p1_com_new, p1_com_Z_new,
								p1_com[1], p1_com[2], p1_com[3])

			# <-0> Boost back from (COM) frame to (C-Z) frame
			boost4_By3_back(p1_cell_Z_new, p1_com_new, vcom)

			# Rotate back to from (C-Z) frame to (C) frame
			rotate_back_from_D(p1_cell_new, p1_cell_Z_new,
								p1_cell[1], p1_cell[2], p1_cell[3])

			### take down gluon absorption / emission (Not necessary)
			if channel in [4,5]:
				rotate_back_from_D(k_gluon, self.lbt.IS[2],
								p1_cell[1], p1_cell[2], p1_cell[3])
			elif channel in [2,3]:
				# Rotate back to from (Com-Z) frame to (COM) frame
				rotate_back_from_D(k_gluon, self.lbt.FS[2],
									p1_com[1], p1_com[2], p1_com[3])

				# <-0> Boost back from (COM) frame to (C-Z) frame
				boost4_By3_back(pbuffer, k_gluon, vcom)

				# Rotate back to from (C-Z) frame to (C) frame
				rotate_back_from_D(k_gluon, pbuffer,
									p1_cell[1], p1_cell[2], p1_cell[3])
			else:
				k_gluon = [0,0,0,0]
		# return updated momentum of heavy quark
		return p1_cell_new, k_gluon

	cpdef HQ_hist(self):
		cdef vector[particle].iterator it = self.HQ_list.begin()
		cdef vector[ vector[double] ] p, x
		p.clear()
		x.clear()
		while it != self.HQ_list.end():
			p.push_back(deref(it).p)
			x.push_back(deref(it).x)
			inc(it)
		return np.array(p), np.array(x)

	# Keep the direction of motion, rescale the energy of each heavy quark to E0
	cpdef reset_HQ_energy(self, E0):
		cdef vector[particle].iterator it = self.HQ_list.begin()
		cdef double pabs_new, pabs_old, px, py, pz
		while it != self.HQ_list.end():
			ratio = sqrt( (E0**2 - self.M**2)/(deref(it).p[0]**2 - self.M**2) )
			px = deref(it).p[1]*ratio
			py = deref(it).p[2]*ratio
			pz = deref(it).p[3]*ratio
			deref(it).p = [E0, px, py, pz]
			inc(it)

	# Reset heavy quark time and last interacting time to the current sys time
	cpdef reset_HQ_time(self):
		cdef vector[particle].iterator it = self.HQ_list.begin()
		while it != self.HQ_list.end():
			deref(it).x[0] = self.tau
			deref(it).t_last_23 = self.tau
			deref(it).t_last_32 = self.tau
			inc(it)

	# Return certain hydro field
	cpdef get_hydro_field(self, key):
		return self.hydro_reader.get_current_frame(key)



	cpdef output_oscar(self, filename):
		cdef vector[particle].iterator it = self.active_HQ.begin()
		cdef size_t i=0
		with open(filename, 'w') as f:
			line = ff.FortranRecordWriter('i10,19(2x,d12.6)')
			while it != self.active_HQ.end():
				f.write(line.write([deref(it).pid,
					deref(it).p[1],deref(it).p[2],
					deref(it).p[3],deref(it).p[0],
					self.M,
					deref(it).x[1],deref(it).x[2],
					deref(it).x[3],deref(it).x[0],
					deref(it).Tf,
					deref(it).vcell[0], deref(it).vcell[1], deref(it).vcell[2],
					deref(it).initp[1], deref(it).initp[2],
					deref(it).initp[3], deref(it).initp[0],
					deref(it).s1, deref(it).s2])+'\n')
				i += 1
				inc(it)
		return
