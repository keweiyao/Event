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
import sys
import fortranformat as ff


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
		double t_last, t_last2
		int Nc, Nc2
		int count22, count23, count32
		vector[double] initp
		vector[double] vcell
		double weight, Tf
		int pid

#-----------Production vertex sampler class-------------------------
cdef class XY_sampler:
	cdef np.ndarray Taa, IntTaa
	cdef double dxy
	cdef size_t Nx, Ny
	def __cinit__(self, Taa, dxy):
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
		nx += np.random.rand()
		ny += np.random.rand()
		# to be examined
		return (nx - self.Nx/2.)*self.dxy, (ny - self.Ny/2.)*self.dxy
		
#-----------Event Class---------------------------------------------
cdef void freestream(vector[particle].iterator it, double dt):
	cdef double da = dt/deref(it).p[0]
	deref(it).x[0] = deref(it).x[0] + dt
	deref(it).x[1] = deref(it).x[1] + deref(it).p[1]*da
	deref(it).x[2] = deref(it).x[2] + deref(it).p[2]*da
	deref(it).x[3] = deref(it).x[3] + deref(it).p[3]*da


cdef class event:
	cdef object hydro_reader, hqsample
	cdef str mode, transport
	cdef double M
	cdef vector[particle] active_HQ
	cdef vector[particle] frzout_HQ
	cdef double tau0, dtau, tau
	cdef double deltat_lrf
	cdef double Tc
	cdef double lambda_rescale, Kfactor

	def __cinit__(self, medium_flags, physics_flags, Tc=0.154, table_folder='./tables', refresh_table=False):
		self.mode = medium_flags['type']
		self.transport = physics_flags['physics']
		if self.transport == 'LBT':
			self.Kfactor = physics_flags['Kfactor']
			self.lambda_rescale = physics_flags['lambda_rescale']
		else:
			self.Kfactor = 1.
			self.lambda_rescale = 1.
		self.M = physics_flags['mass']
		self.hydro_reader = medium.Medium(medium_flags=medium_flags)
		self.tau0 = self.hydro_reader.init_tau()
		self.dtau = self.hydro_reader.dtau()
		self.tau = self.tau0
		self.Tc = Tc
		

		if self.transport == "LBT":
			self.hqsample = HqEvo.HqEvo(options=physics_flags,		
										table_folder=table_folder, refresh_table=refresh_table)
		elif self.transport == "LGV":
			self.deltat_lrf = physics_flags['dt_lrf']
			self.hqsample = HqLGV.HqLGV(options=physics_flags,
										table_folder=table_folder, refresh_table=refresh_table) 
	cpdef sys_time(self) :
		return self.tau

	cpdef initialize_HQ(self, NQ, init_flags):
		self.active_HQ.clear()
		self.active_HQ.resize(NQ)
		cdef double x, y, z
		# for A+B:
		cdef double pT, phipt, rapidity, mT, freetime
		cdef double ymin, ymax, pTmax, pTmin, alpha
		# for box:
		cdef double p, cospz, sinpz
		cdef double pmax, L
		cdef vector[particle].iterator it

		if init_flags['type'] == 'A+B':   
			print "Initialize for dynamic medium"
			HQ_xy_sampler = XY_sampler(init_flags['TAB'], init_flags['dxy'])
			pTmax = init_flags['pTmax']
			pTmin = init_flags['pTmin']
			ymax = init_flags['ymax']
			ymin = init_flags['ymin']
			alpha = init_flags['sample power']
			print "Heavy quark will be free streamed to the starting time of hydro"
			it = self.active_HQ.begin()
			while it != self.active_HQ.end():
				pT = -1.
				while pT < pTmin or pT > pTmax:
					if alpha > 0.:
						pT = (rand()*1./RAND_MAX)**alpha*(pTmax)
					if alpha < 0.:
						pT = (rand()*1./RAND_MAX)**alpha*(pTmin)
				mT = sqrt(pT**2 + self.M**2)
				phipt = rand()*2.*M_PI/RAND_MAX
				rapidity = rand()*(ymax-ymin)/RAND_MAX + ymin
				deref(it).p.resize(4)
				deref(it).x.resize(4)
				deref(it).p = [mT*cosh(rapidity), pT*cos(phipt), pT*sin(phipt), mT*sinh(rapidity)]
				deref(it).t_last = 0.; deref(it).t_last2 = 0.
				deref(it).Nc = 0; deref(it).Nc2 = 0
				deref(it).count22 = 0; deref(it).count23 = 0; deref(it).count32 = 0
				deref(it).freezeout = False
				deref(it).initp = deref(it).p
				deref(it).vcell = [0., 0., 0.]
				deref(it).Tf = 0.
				deref(it).pid = 4*np.random.choice([1, -1])
				# free streaming to hydro starting time
				freetime = self.tau0/sqrt(1. - (deref(it).p[3]/deref(it).p[0])**2)
				x, y = HQ_xy_sampler.sample_xy()
				z = 0.
				deref(it).x = [0.0, x, y, z]
				deref(it).t_last = freetime; deref(it).t_last2 = freetime
				freestream(it, freetime)
				inc(it)

		if init_flags['type'] == 'box':
			print "Initialize for box simulation"
			pmax = init_flags['pmax']
			L = init_flags['L']
			it = self.active_HQ.begin()
			while it != self.active_HQ.end():					  
				p = pow(rand()*1./RAND_MAX, 1./3.)*pmax
				phipt = rand()*2.*M_PI/RAND_MAX
				cospz = little_below_one*(rand()*2./RAND_MAX-1.)
				sinpz = sqrt(1.-cospz**2)
				deref(it).p.resize(4)
				deref(it).x.resize(4)
				deref(it).p = [sqrt(p**2+self.M**2), p*sinpz*cos(phipt), p*sinpz*sin(phipt), p*cospz]		   
				deref(it).t_last = 0.0; deref(it).t_last2 = 0.0
				deref(it).Nc = 0; deref(it).Nc2 = 0
				deref(it).count22 = 0; deref(it).count23 = 0; deref(it).count32 = 0
				deref(it).freezeout = False
				deref(it).initp = [sqrt(p**2+self.M**2), p*sinpz*cos(phipt), p*sinpz*sin(phipt), p*cospz]		   
				deref(it).vcell = [0., 0., 0.]
				deref(it).Tf = 0.
				deref(it).pid = 4*np.random.choice([1, -1])
				x = rand()*L*2./RAND_MAX - L
				y = rand()*L*2./RAND_MAX - L
				z = rand()*L*2./RAND_MAX - L
				deref(it).x = [0.0, x, y, z]
				inc(it) 

	cpdef bool perform_hydro_step(self, StaticPropertyDictionary=None) :
		PyErr_CheckSignals()
		if self.mode == 'dynamic':
			self.hydro_reader.load_next()
		elif self.mode == 'static':
			if StaticPropertyDictionary==None:
				raise ValueError("static meidum property not defined")
			else:
				self.hydro_reader.load_next(StaticPropertyDictionary=StaticPropertyDictionary)
		else:
			raise ValueError("medium mode not defined")
		status = self.hydro_reader.hydro_status()

		self.tau += self.dtau
		cdef double t, x, y, z, tauQ
		cdef vector[particle].iterator it
		if self.mode == 'static':
			it = self.active_HQ.begin()
			while it != self.active_HQ.end():
				while deref(it).x[0] <= self.tau:
					self.perform_HQ_step(it)
				inc(it)
		elif self.mode == 'dynamic':
			it = self.active_HQ.begin()
			while it != self.active_HQ.end():
				if deref(it).freezeout == False:
					t, x, y, z = deref(it).x
					tauQ = sqrt(t**2 - z**2)
					while tauQ <= self.tau:
						self.perform_HQ_step(it)
						if deref(it).freezeout:
							break
						t, x, y, z = deref(it).x
						tauQ = sqrt(t**2 - z**2)
				inc(it)
		else:
			raise ValueError("Mode not implemented")
		return status

	cdef perform_HQ_step(self, vector[particle].iterator it): 
		PyErr_CheckSignals()
		cdef double t, x, y, z, tauQ, T, vabs2, vabs
		cdef double t_elapse_lab, t_elapse_lab2
		cdef double dtHQ
		cdef vector[double] pnew, vcell
		cdef int channel
		
		t, x, y, z = deref(it).x
		tauQ = sqrt(t**2 - z**2)
		vcell.resize(3)
		T, vcell[0], vcell[1], vcell[2] = self.hydro_reader.interpF(tauQ, deref(it).x, ['Temp', 'Vx', 'Vy', 'Vz'])
		if T <= self.Tc:
			deref(it).freezeout = True
			deref(it).Tf = T
			deref(it).vcell = vcell
			return 
		
		vabs2 = vcell[0]**2 + vcell[1]**2 + vcell[2]**2
		if vabs2 >= little_below_one**2:
			scale = 1.0/sqrt(vabs2)/little_above_one
			vcell[0] *= scale
			vcell[1] *= scale
			vcell[2] *= scale
		
		if self.transport == 'LBT':
			t_elapse_lab = (t - deref(it).t_last)/(deref(it).Nc + 1.) / GeVm1_to_fmc	# convert to GeV-1
			t_elapse_lab2 = (t - deref(it).t_last2)/(deref(it).Nc2 + 1.) / GeVm1_to_fmc # convert to GeV-1
			channel, dtHQ, pnew = self.update_HQ_LBT(deref(it).p, vcell, T, t_elapse_lab, t_elapse_lab2)		      
			dtHQ *= GeVm1_to_fmc   # convert back to fm/c 
			if channel == 0 or channel == 1:
				deref(it).Nc += 1
				deref(it).Nc2 += 1
				deref(it).count22 += 1
			elif channel == 2 or channel == 3:
				deref(it).Nc = 0
				deref(it).Nc2 += 1
				deref(it).t_last = t + dtHQ
				deref(it).count23 += 1
			elif channel == 4 or channel == 5:
				deref(it).Nc2 = 0
				deref(it).Nc += 1
				deref(it).t_last2 = t + dtHQ
				deref(it).count32 += 1
			else:
				pass
			freestream(it, dtHQ)
			deref(it).p = pnew
		elif self.transport == 'LGV':
			channel, dtHQ, pnew = self.update_HQ_LGV(deref(it).p, vcell, T)
			freestream(it, dtHQ)
			deref(it).p = pnew
			# for Langevin transport, the time dtHQ is already in fm/c unit 
		else:
			pass
			#raise ValueError("Transport mode not recongized.")

	cdef (int, double, vector[double]) update_HQ_LGV(self, vector[double] p1_lab, vector[double] v3cell, double Temp) :
		cdef vector[double] p1_cell, p1_cell_Z_new, p1_cell_new 
		#Boost from p1(lab) to p1(cell)
		boost4_By3(p1_cell, p1_lab, v3cell)
		# there is no need to use p1_cell_Z, since we only need energy to do the Langevin transportation, and returns in Z
		self.hqsample.update_by_Langevin(p1_cell[0], Temp)
		p1_cell_Z_new = self.hqsample.post_result
		rotate_back_from_D(p1_cell_new, p1_cell_Z_new, p1_cell[1], p1_cell[2], p1_cell[3])
		cdef vector[double] pnew
		boost4_By3_back(pnew, p1_cell_new, v3cell)

		cdef double dt_cell = self.deltat_lrf
		cdef double dtHQ = p1_lab[0]/p1_cell[0] * dt_cell

		cdef int channel = 10  #? need to change this, reserve a spectial number for Langevin transport
		return channel, dtHQ, pnew


				
	# this function needs to be cdef so that we can use the reference copy
	cdef (int, double, vector[double]) update_HQ_LBT(self, vector[double] p1_lab, vector[double] v3cell, double Temp, double mean_dt23_lab, double mean_dt32_lab) :
		# Define local variables
		cdef double s, L1, L2, Lk, x2, xk, a1=0.6, a2=0.6
		cdef double dt_cell
		cdef size_t i=0
		cdef int channel
		cdef vector[double] p1_cell, p1_cell_Z, p1_com, \
			 p1_com_Z_new, p1_com_new, \
			 p1_cell_Z_new, p1_cell_new,\
			 pnew, fs, Pcom,	 \
			 dx23_cell, dx32_cell, dx23_com, \
			 v3com, pbuffer

		# Boost p1 (lab) to p1 (cell)
		boost4_By3(p1_cell, p1_lab, v3cell)

		# displacement vector within dx23_lab and dx32_lab seen from cell frame
		dx23_cell.resize(4)
		dx32_cell.resize(4)
		for i in range(4):
			dx23_cell[i] = p1_cell[i]/p1_lab[0]*mean_dt23_lab
			dx32_cell[i] = p1_cell[i]/p1_lab[0]*mean_dt32_lab

		# Sample channel in cell, return channl index and evolution time seen from cell
		channel, dt_cell = self.hqsample.sample_channel(p1_cell[0], Temp, dx23_cell[0]*self.lambda_rescale, dx32_cell[0]*self.lambda_rescale)
		# Boost evolution time back to lab frame
		cdef double dtHQ = p1_lab[0]/p1_cell[0]*dt_cell
		
		# If not scattered, return channel=-1, evolution time in Lab frame and origin al p1(lab)
		if channel < 0:
			pnew = p1_lab
		else:
			# Sample initial state and return initial state particle four vectors 
			# Imagine rotate p1_cell to align with z-direction, construct p2_cell_align, ...				
			self.hqsample.sample_initial(channel, p1_cell[0], Temp, dx23_cell[0]*self.lambda_rescale, dx32_cell[0]*self.lambda_rescale)
			p1_cell_Z = self.hqsample.IS[0]

			# Center of mass frame of p1_cell_align and other particles, and take down orientation of p1_com
			Pcom.resize(4)
			for i in range(4):
				Pcom[i] = 0.
			for pp in self.hqsample.IS:
				for i in range(4):
					Pcom[i] += pp[i]
			s = product4(Pcom, Pcom)
		
			v3com.resize(3); 
			for i in range(3):
				v3com[i] = Pcom[i+1]/Pcom[0];

			boost4_By3(p1_com, p1_cell_Z, v3com)
			if channel >= 4: # for 3 -> 2 kinetics
				L1 = sqrt(p1_com[0]**2 - self.M**2)
				boost4_By3(pbuffer, self.hqsample.IS[1], v3com)
				L2 = pbuffer[0]
				boost4_By3(pbuffer, self.hqsample.IS[2], v3com)
				Lk = pbuffer[0]
				x2 = L2/(L1+L2+Lk)
				xk = Lk/(L1+L2+Lk)
				a1 = x2 + xk; 
				a2 = (x2 - xk)/(1. - a1)
		
			dx23_com.resize(4)
			for i in range(4):
				dx23_com[i] = p1_com[i]/p1_cell_Z[0]*dx23_cell[0]

			# Sample final state momentum in Com frame, with incoming paticles on z-axis
			self.hqsample.sample_final(channel, s, Temp, dx23_com[0], a1, a2)

			p1_com_Z_new = self.hqsample.FS[0]
			# Rotate final states back to original Com frame (not z-axis aligened)
			rotate_back_from_D(p1_com_new, p1_com_Z_new, p1_com[1], p1_com[2], p1_com[3])
			# boost back to cell frame z-align
			boost4_By3_back(p1_cell_Z_new, p1_com_new, v3com)	   
			# rotate back to original cell frame
			rotate_back_from_D(p1_cell_new, p1_cell_Z_new, p1_cell[1], p1_cell[2], p1_cell[3])
			# boost back to lab frame
			boost4_By3_back(pnew, p1_cell_new, v3cell)
		# return updated momentum of heavy quark
			
		return channel, dtHQ, pnew

	cpdef HQ_hist(self):
		cdef vector[particle].iterator it = self.active_HQ.begin()
		cdef vector[ vector[double] ] p, x
		p.clear()
		x.clear()
		while it != self.active_HQ.end():
			p.push_back(deref(it).p)
			x.push_back(deref(it).x)
			inc(it)
		return np.array(p), np.array(x)

	cpdef Init_pT(self):
		cdef vector[particle].iterator it = self.active_HQ.begin()
		cdef vector[double] pT
		pT.clear()
		while it != self.active_HQ.end():
			pT.push_back(sqrt(deref(it).p[1]**2+deref(it).p[2]**2))
			inc(it)
		return np.array(pT)


	cpdef output_oscar(self, filename):
		cdef vector[particle].iterator it = self.active_HQ.begin()
		cdef size_t i=0
		with open(filename, 'w') as f:
			head3 = ff.FortranRecordWriter(
				'2(a8,2x),1x,i3,1x,i6,3x,i3,1x,i6,3x,a4,2x,e10.4,2x,i8')
			f.write('OSC1997A\n')
			f.write('final_id_p_x\n')
			f.write(head3.write(['lbt', '1.0alpha', 208, 82, 208, 82, 'aacm', 1380, 1])+'\n')
			eventhead = ff.FortranRecordWriter('i10,2x,i10,2x,f8.3,2x,f8.3,2x,i4,2x,i4,2X,i7')
			f.write(eventhead.write([1, self.active_HQ.size(), 0.001, 0.001, 1, 1, 1])+'\n')
			line = ff.FortranRecordWriter('i10,2x,i10,17(2x,d12.6)')
			while it != self.active_HQ.end():
				f.write(line.write([i, deref(it).pid, 
					deref(it).p[1],deref(it).p[2],deref(it).p[3],deref(it).p[0], 						self.M, 
					deref(it).x[1],deref(it).x[2],deref(it).x[3],deref(it).x[0], 						deref(it).Tf, 
					deref(it).vcell[0], deref(it).vcell[1], deref(it).vcell[2],
					 deref(it).initp[1],deref(it).initp[2],deref(it).initp[3], 
					1.])+'\n')
				i += 1
				inc(it)
		return

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

	cpdef get_hydro_field(self, key):
		return self.hydro_reader.get_current_frame(key)














