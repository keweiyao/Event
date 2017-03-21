import numpy as np
import sys
import event
import h5py
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

# Medium option
box_info = {'Temp'  : 0.3, 
			'Vx'	: 0.0, 
			'Vy'	: 0.0, 
			'Vz'	: 0.0   }
static_config = {	'type'	  : 'static',	
					 'static_dt' : 1.  }

hydro_history_filepath = sys.argv[1]
dynamic_config = {  'type'	  : 'dynamic', 
					'hydrofile' : hydro_history_filepath	}

# Physics option
LBT_config = {  'physics'   : 'LBT',
                                '2->2'    : True,
                                '2->3'    : True,
                                '3->2'    : False,
                                'Nf'        : 3,
                                'mass'    : 1.3 }  

LGV_config = {  'physics'   : 'LGV',
                                'dt_lrf'        : 0.02,
                                'elastic'   : True,
                                'Einstein'  : False,
                                'Nf'        : 3,
                                'mass'    : 1.3 } 

# Initialization option
box_init = {	'type'  : 'box',
				'L'	 : 10.,
				'pmax'  : 10.   }

TAA = np.loadtxt(sys.argv[2]).T ## p=0.
realistic_init =  { 'type'		  : 'A+B',
					'sample power'  : 1.,
					'pTmin'		 : 0.1,
					'pTmax'		 : 70.,
					'ymin'		  : -1.,
					'ymax'		  : 1.,
					'TAB'		   : TAA,
					'dxy'		   : 0.1   }

				
e1 = event.event(   medium_flags=dynamic_config , 
					physics_flags=LBT_config   )

e1.initialize_HQ(   NQ=10000,
					init_flags=realistic_init   )

# Run Model
f = h5py.File("particledata.hdf5", 'w')
Init_pT = e1.Init_pT()
f.create_dataset('init_pT', data=Init_pT)
for i in range(500):
	print("t = %1.2f [fm/c]"%e1.sys_time() )
	status = e1.perform_hydro_step()#StaticPropertyDictionary=box_info)
	if i%10 == 0:
		dsp, dsx = e1.HQ_hist()
		f.create_dataset('p-%d'%i, data=dsp)
		f.create_dataset('x-%d'%i, data=dsx)
	#plt.clf()
	#plt.scatter(dsx.T[1], dsx.T[2], s=0.4, color='orange', alpha=0.3)
	#T = e1.get_hydro_field('Temp')
	#plt.imshow(np.flipud(T.T), extent = [-13., 13., -13., 13.])
	#plt.pause(0.1)
	if not status:
		break
f.close()
