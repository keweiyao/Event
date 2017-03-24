import numpy as np
import sys
import h5py
import matplotlib.pyplot as plt
# This is "main" script of LBT-LGV code
# First import the heavy quark event module
import event
# To declare an event module, we need a bunch of options
# They are split into two categories
# Event options: 
# 1. Medium options: options are organized in a dictionary by { 'keyward' : value }

# For example, a medium-option dictionary config a simulation with a static medium
static_config = {	'type'	  : 'static',	
					'static_dt' : 1.  }
# The following medium information is to be used later
box_info = {'Temp'  : 0.3, 
			'Vx'	: 0.0, 
			'Vy'	: 0.0, 
			'Vz'	: 0.0   }

# Another example, a medium-option dictionary config a dynamical medium
# with the path to the hydro-history file 
hydro_history_filepath = sys.argv[1]
dynamic_config = {  'type'	  : 'dynamic', 
					'hydrofile' : hydro_history_filepath	}

# Physics option
# 1. An example for linear-Boltzmann evolution
LBT_config = {  'physics'   : 'LBT',
                                '2->2'    : True,
                                '2->3'    : True,
                                '3->2'    : True,
                                'Nf'        : 3,
                                'mass'    : 1.3 }  
                                
# 1. An example for Langevin evolution
LGV_config = {  'physics'   : 'LGV',
                                'dt_lrf'        : 0.02,
                                'elastic'   : True,
                                'Einstein'  : True,
                                'Nf'        : 3,
                                'mass'    : 1.3 } 

# Initialization option
# Initizlize HQ in a static box [-L, L]^3
box_init = {	'type'  : 'box',
				'L'	 : 10.,
				'pmax'  : 10.   }

TAA = np.loadtxt(sys.argv[2])
realistic_init =  { 'type'		  : 'A+B',
					'sample power'  : 1.,
					'pTmin'		 : 0.1,
					'pTmax'		 : 70.,
					'ymin'		  : -1.,
					'ymax'		  : 1.,
					'TAB'		   : TAA,
					'dxy'		   : 0.1   }
sd = TAA**0.5
				
e1 = event.event(   medium_flags=dynamic_config , 
					physics_flags=LBT_config   )

e1.initialize_HQ(   NQ=40000,
					init_flags=realistic_init   )

# Run Model
f = h5py.File("particledata.hdf5", 'w')
Init_pT = e1.Init_pT()
f.create_dataset('init_pT', data=Init_pT)
"""
status = e1.perform_hydro_step()
print("t = %1.2f [fm/c]"%e1.sys_time() )

plt.clf()
plt.subplot(1,3,1,aspect='equal')
s = (e1.get_hydro_field('e') + e1.get_hydro_field('p') ) / (1e-9 + e1.get_hydro_field('Temp'))
plt.imshow(np.flipud(s.T), extent = [-15., 15., -15., 15.])

plt.subplot(1,3,2,aspect='equal')
plt.imshow(np.flipud(sd.T), extent = [-15., 15., -15., 15.])

plt.subplot(1,3,3,aspect='equal')
dsp, dsx = e1.HQ_hist()
plt.hist2d(dsx.T[1], dsx.T[2], bins=100, range=[[-15., 15.], [-15., 15.]])
plt.show()
"""

for i in range(500):
	print("t = %1.2f [fm/c]"%e1.sys_time() )
	status = e1.perform_hydro_step()#StaticPropertyDictionary=box_info)
	if i%10 == 0:
		dsp, dsx = e1.HQ_hist()
		f.create_dataset('p-%d'%i, data=dsp)
		f.create_dataset('x-%d'%i, data=dsx)
		"""
		plt.clf()
		plt.subplot(1,2,1,aspect='equal')
		s = (e1.get_hydro_field('e') + e1.get_hydro_field('p') ) / (1e-9 + e1.get_hydro_field('Temp'))
		plt.imshow(np.flipud(s.T), extent = [-15., 15., -15., 15.])
		
		plt.subplot(1,2,2,aspect='equal')
		plt.scatter(dsx.T[1, :1000], dsx.T[2, :1000],s=1., alpha=0.4)
		plt.axis([-15., 15., -15., 15.])
		#print ('x2 - y2', np.average((dsx.T[1]**2 - dsx.T[2]**2)/(dsx.T[1]**2 + dsx.T[2]**2)))
		#print ('px2 - py2',  np.average((dsp.T[1]**2 - dsp.T[2]**2)/(dsp.T[1]**2 + dsp.T[2]**2)))
		plt.pause(0.2)
		"""
	if not status:
		break

f.close()
