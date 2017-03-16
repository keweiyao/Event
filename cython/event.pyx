from libcpp cimport bool
from libcpp.vector cimport vector
from libc.stdlib cimport rand, RAND_MAX
from libc.math cimport *
from cython.operator cimport dereference as deref, preincrement as inc
import numpy as np
import sys
cimport numpy as np

sys.path.append('../HQ-Evo/')
sys.path.append('../Medium-Reader')
import HqEvo
import medium
import HqLGV

from event cimport *

cdef double GeVm1_to_fmc = 0.197
cdef double little_below_one = 1. - 1e-7
cdef double little_above_one = 1. + 1e-7

cdef double default_df_dpt2dy(double pT, double y):
        return 1.0/(1.0+pT**4)*1.0/cosh(y)

#-----------Particle data struct------------------------------------
cdef extern from "../src/utility.h":
        cdef struct particle:
                vector[double] p
                vector[double] x
                double t_last, t_last2
                int Nc, Nc2
                int count22, count23, count32
                double weight

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
                return (nx - self.Nx/2.)*self.dxy, (ny - self.Ny/2.)*self.dxy
#-----------Event Class---------------------------------------------

cdef bool freestream(vector[particle].iterator it, double dt):
        deref(it).x[0] = deref(it).x[0] + dt
        deref(it).x[1] = deref(it).x[1] + deref(it).p[1]/deref(it).p[0]*dt
        deref(it).x[2] = deref(it).x[2] + deref(it).p[2]/deref(it).p[0]*dt
        deref(it).x[3] = deref(it).x[3] + deref(it).p[3]/deref(it).p[0]*dt
        return True


cdef class event:
        cdef object hydro_reader, hqsample
        cdef str mode,  transport
        cdef double M
        cdef vector[particle] active_HQ
        cdef vector[particle] frzout_HQ
        cdef double tau0, dtau, tau
        cdef int channel
        cdef vector[double] pnew
        cdef double dtHQ
        cdef double deltat_lrf

        def __cinit__(self, mode="dynamic", transport="LBT", hydrofile=None, static_dt=0.5, mass=1.3, elastic=True, inelastic=False, detailed_balance=False, table_folder='./tables'):
                self.mode = mode
                self.transport = transport
                self.M = mass
                self.hydro_reader = medium.Medium(mode=mode, hydrofilename=hydrofile, static_dt=static_dt)
                self.tau0 = self.hydro_reader.init_tau()
                self.dtau = self.hydro_reader.dtau()
                self.tau = self.tau0
                self.deltat_lrf = min(self.dtau/10, 0.01)

                if transport == "LBT":
                        self.hqsample = HqEvo.HqEvo(mass=mass, elastic=elastic, inelastic=inelastic, detailed_balance=detailed_balance, table_folder=table_folder)
                elif transport == "LGV":
                        self.hqsample = HqLGV.HqLGV(mass=mass, elastic=elastic, EinR=False, deltat_lrf=self.deltat_lrf, table_folder=table_folder) 
        def sys_time(self):
                return self.tau

        cpdef initialize_HQ(self, NQ=1000, df_dpt2dy=None, Taa=None, dxy=0.1, oversample_power=1.):
                self.active_HQ.clear()
                self.active_HQ.resize(NQ)
                cdef double x, y, z
                cdef double pT, phipt, rapidity, mT, pTmax = 70., pTmin=0.1, freetime
                cdef double p, cospz, sinpz, pmax=50.
                cdef vector[particle].iterator it

                if self.mode == 'dynamic':      
                        print "Initialize for dynamic medium"
                        if Taa == None:
                                raise ValueError("Need Taa to sample x-y postion of heavy quark in dynamic mode")
                        else:
                                print "sample x-y position accoring to Taa"
                        HQ_xy_sampler = XY_sampler(Taa, dxy)
                        print "Uniform and independent sampling HQ df/dpT and df/dy"
                        print "0.1 < pt < 70. [GeV], -1. < y < 1."
                        print "Each heavy quark is then weighted by pt*df/dpt^2/dy"
                        if df_dpt2dy==None:
                                print "Notice: user define df/dpt^2/dy spectra not available, use defaule weight"
                                df_dpt2dy = default_df_dpt2dy
                        print "Heavy quark will be free streamed to the starting time of hydro"
                        it = self.active_HQ.begin()
                        while it != self.active_HQ.end():
                                pT = 0.
                                while pT < pTmin or pT > pTmax:
                                        if oversample_power > 0.:
                                                pT = (rand()*1./RAND_MAX)**oversample_power*(pTmax)
                                        if oversample_power < 0.:
                                                pT = (rand()*1./RAND_MAX)**oversample_power*(pTmin)
                                mT = sqrt(pT**2 + self.M**2)
                                phipt = rand()*2.*M_PI/RAND_MAX
                                rapidity = rand()*2./RAND_MAX-1.
                                deref(it).p.resize(4)
                                deref(it).x.resize(4)
                                deref(it).p = [mT*cosh(rapidity), pT*cos(phipt), pT*sin(phipt), mT*sinh(rapidity)]
                                deref(it).t_last = 0.; deref(it).t_last2 = 0.
                                deref(it).Nc = 0; deref(it).Nc2 = 0
                                deref(it).count22 = 0; deref(it).count23 = 0; deref(it).count32 = 0
                                deref(it).weight = pT**(2. - 1./oversample_power)*df_dpt2dy(pT, y)

                                # free streaming to hydro starting time
                                free_time = self.tau0/sqrt(1. - (deref(it).p[3]/deref(it).p[0])**2)
                                x, y = HQ_xy_sampler.sample_xy()
                                z = 0.
                                deref(it).x = [0.0, x, y, z]
                                deref(it).t_last = free_time; deref(it).t_last2 = free_time
                                freestream(it, free_time)
                                inc(it)

                if self.mode == 'static':
                        print "Initialize for static medium"
                        print "Uniform sampling HQ momentum"
                        print "0. < |p| < 50. [GeV]"
                        it = self.active_HQ.begin()
                        while it != self.active_HQ.end():                       
                                p = pow(rand()*1./RAND_MAX, 1./3.)*pmax
                                phipt = rand()*2.*M_PI/RAND_MAX
                                cospz = little_below_one*(rand()*2./RAND_MAX-1.)
                                sinpz = sqrt(1.-cospz**2)
                                deref(it).p.resize(4)
                                deref(it).x.resize(4)
                                deref(it).p = [sqrt(p**2+self.M**2), p*sinpz*cos(phipt), p*sinpz*sin(phipt), p*cospz]           
                                deref(it).p = [10., 0., 0., sqrt(100.-1.3**2)]
                                deref(it).t_last = 0.0; deref(it).t_last2 = 0.0
                                deref(it).Nc = 0; deref(it).Nc2 = 0
                                deref(it).count22 = 0; deref(it).count23 = 0; deref(it).count32 = 0
                                deref(it).weight = 1.0
                                x = rand()*10./RAND_MAX - 5.
                                y = rand()*10./RAND_MAX - 5.
                                z = rand()*10./RAND_MAX - 5.
                                deref(it).x = [0.0, x, y, z]
                                inc(it) 

        cpdef bool perform_hydro_step(self, StaticPropertyDictionary=None):
                cdef bool status = True
                if self.mode == 'dynamic':
                        status = self.hydro_reader.load_next()
                if self.mode == 'static':
                        if StaticPropertyDictionary==None:
                                raise ValueError("static meidum property not defined")
                        else:
                                status = self.hydro_reader.load_next(StaticPropertyDictionary=StaticPropertyDictionary)
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
                                count = 0
                                t, x, y, z = deref(it).x
                                tauQ = sqrt(t**2 - z**2)
                                while tauQ <= self.tau:
                                        self.perform_HQ_step(it)
                                        t, x, y, z = deref(it).x
                                        tauQ = sqrt(t**2 - z**2)
                                inc(it)
                else:
                        raise ValueError("Mode not implemented")
                return status

        cdef bool perform_HQ_step(self, vector[particle].iterator it):
                cdef double t, x, y, z, t_elapse_lab, t_elapse_lab2, tauQ
                cdef double T, vx, vy, vz
                t, x, y, z = deref(it).x
                t_elapse_lab = (t - deref(it).t_last)/(deref(it).Nc + 1.)
                t_elapse_lab2 = (t - deref(it).t_last2)/(deref(it).Nc2 + 1.)
                
                tauQ = sqrt(t**2 - z**2)
                T, vx, vy, vz = self.hydro_reader.interpF(tauQ, [x, y, z, t], ['Temp', 'Vx', 'Vy', 'Vz'])
                cdef double v2 = vx**2 + vy**2 + vz**2
                if v2 > little_below_one:
                        vx /= v2*little_above_one; vy /= v2*little_above_one; vz /= v2*little_above_one

                # Note the change of units from GeV-1 (fm/c) to fm/c (GeV-1)
                t_elapse_lab /= GeVm1_to_fmc
                t_elapse_lab2 /= GeVm1_to_fmc
                cdef vector[double] vcell = [vx, vy, vz]
                if self.transport == 'LBT':
                        self.update_HQ_LBT(deref(it).p, vcell, T, t_elapse_lab, t_elapse_lab2)              
                        self.dtHQ *= GeVm1_to_fmc
                elif self.transport == 'LGV':
                        self.update_HQ_LGV(deref(it).p, vcell, T)
                        # for Langevin transport, the time dtHQ is already in fm/c unit 
                
                if self.channel == 0 or self.channel == 1:
                        deref(it).Nc = deref(it).Nc + 1
                        deref(it).Nc2 = deref(it).Nc2 + 1
                if self.channel == 2 or self.channel == 3:
                        deref(it).Nc = 0
                        deref(it).Nc2 = deref(it).Nc2 + 1
                        deref(it).t_last = t + self.dtHQ
                if self.channel == 4 or self.channel == 5:
                        deref(it).Nc2 = 0
                        deref(it).Nc = deref(it).Nc + 1
                        deref(it).t_last2 = t + self.dtHQ
                freestream(it, self.dtHQ)
                deref(it).p = self.pnew
                return True


                
        cdef bool update_HQ_LGV(self, vector[double] p1_lab, vector[double] v3cell, double Temp):
                cdef vector[double] p1_cell, p1_cell_Z_new, p1_cell_new, p1_new 
                cdef double dt_lab, dt_cell
                #Boost from p1(lab) to p1(cell)
                p1_cell = boost4_By3(p1_lab, v3cell)
                # there is no need to use p1_cell_Z, since we only need energy to do the Langevin transportation, and returns in Z
                p1_cell_Z_new = self.hqsample.update_by_Langevin(p1_cell[0], Temp)
                p1_cell_new = rotate_back_from_D(p1_cell_Z_new, p1_cell[1], p1_cell[2], p1_cell[3])
                p1_new = boost4_By3(p1_cell_new, [-v3cell[0], -v3cell[1], -v3cell[2]])

                dt_cell = self.deltat_lrf
                dt_lab = p1_lab[0]/p1_cell[0] * dt_cell  #? what is this for??

                self.channel = 10  #? need to change this, reserve a spectial number for Langevin transport
                self.dtHQ = dt_lab
                self.pnew = p1_new
                return True


                

        cdef bool update_HQ_LBT(self, vector[double] p1_lab, vector[double] v3cell, double Temp, double mean_dt23_lab, double mean_dt32_lab):
                # Define local variables
                cdef double s, L1, L2, Lk, x2, xk, a1=0.6, a2=0.6
                cdef double dt_lab, dt_cell
                cdef int channel
                cdef size_t i=0
                cdef vector[double] p1_cell, p1_cell_Z, p1_com, \
                                                        p1_com_Z_new, p1_com_new, \
                                                        p1_cell_Z_new, p1_cell_new,\
                                                        p1_new, fs,                             \
                                                        dx23_cell, dx32_cell, dx23_com, \
                                                        v3com 

                # Boost p1 (lab) to p1 (cell)
                p1_cell = boost4_By3(p1_lab, v3cell)

                # displacement vector within dx23_lab and dx32_lab seen from cell frame
                dx23_cell = [pi/p1_lab[0]*mean_dt23_lab for pi in p1_cell]
                dx32_cell = [pi/p1_lab[0]*mean_dt32_lab for pi in p1_cell]

                # Sample channel in cell, return channl index and evolution time seen from cell
                channel, dt_cell = self.hqsample.sample_channel(p1_cell[0], Temp, 0.154, dx23_cell[0], dx32_cell[0])
                
                # Boost evolution time back to lab frame
                dt_lab = p1_lab[0]/p1_cell[0]*dt_cell

                # If not scattered, return channel=-1, evolution time in Lab frame and origin al p1(lab)
                if channel < 0:
                        self.channel = channel
                        self.dtHQ = dt_lab
                        self.pnew = p1_lab
                        return True
                else:
                # Sample initial state and return initial state particle four vectors 
                # Imagine rotate p1_cell to align with z-direction, construct p2_cell_align, ...                
                        self.hqsample.sample_initial(channel, p1_cell[0], Temp, dx23_cell[0], dx32_cell[0])
                        p1_cell_Z = self.hqsample.IS[0]

                        # Center of mass frame of p1_cell_align and other particles, and take down orientation of p1_com
                        i = 0
                        Pcom = [0., 0., 0., 0.]
                        for pp in self.hqsample.IS:
                                for i in range(4):
                                        Pcom[i] += pp[i]
                        s = product4(Pcom, Pcom)
                
                        v3com = [Pcom[i+1]/Pcom[0] for i in range(3)]
                
                        p1_com = boost4_By3(p1_cell_Z, v3com)
                        if channel >= 4: # for 3 -> 2 kinetics
                                L1 = sqrt(p1_com[0]**2 - self.M**2)
                                L2 = boost4_By3(self.hqsample.IS[1], v3com)[0]
                                Lk = boost4_By3(self.hqsample.IS[2], v3com)[0]
                                x2 = L2/(L1+L2+Lk)
                                xk = Lk/(L1+L2+Lk)
                                a1 = x2 + xk; 
                                a2 = (x2 - xk)/(1. - a1)
                
                        dx23_com = [pi/p1_cell_Z[0]*dx23_cell[0] for pi in p1_com]

                        # Sample final state momentum in Com frame, with incoming paticles on z-axis
                        self.hqsample.sample_final(channel, s, Temp, dx23_com[0], a1, a2)               

                        p1_com_Z_new = self.hqsample.FS[0]
                        # Rotate final states back to original Com frame (not z-axis aligened)
                        p1_com_new = rotate_back_from_D(p1_com_Z_new, p1_com[1], p1_com[2], p1_com[3])  
                        # boost back to cell frame z-aligned
                        p1_cell_Z_new = boost4_By3(p1_com_new, [-v3com[0], -v3com[1], -v3com[2]])       
                        # rotate back to original cell frame
                        p1_cell_new = rotate_back_from_D(p1_cell_Z_new, p1_cell[1], p1_cell[2], p1_cell[3])
                        # boost back to lab frame
                        p1_new = boost4_By3(p1_cell_new, [-v3cell[0], -v3cell[1], -v3cell[2]])
                        # return updated momentum of heavy quark

                        self.channel = channel
                        self.dtHQ = dt_lab
                        self.pnew = p1_new
                        return True

        cpdef HQ_hist(self):
                cdef vector[particle].iterator it = self.active_HQ.begin()
                cdef vector[ vector[double] ] p, x
                cdef vector[double] weight
                p.clear()
                x.clear()
                weight.clear()
                while it != self.active_HQ.end():
                        p.push_back(deref(it).p)
                        x.push_back(deref(it).x)
                        weight.push_back(deref(it).weight)
                        inc(it)
                return np.array(p), np.array(x), np.array(weight)

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
        def get_hydro_field(self, key):
                return self.hydro_reader.get_current_frame(key)













