import numpy as np
import sys
import event
import h5py
from scipy.interpolate import interp1d

# FONLL pT spectra
"""
pt, dsigma = np.loadtxt(sys.argv[3]).T
dfdpt2 = interp1d(pt, dsigma)
def pT_weight(pt, y):
	return dfdpt2(pt)
"""

# Medium option
box_info = {'Temp'  : 0.3, 
			'Vx'	: 0.0, 
			'Vy'	: 0.0, 
			'Vz'	: 0.0   }
static_config = {	'type'	  : 'static',	
					 'static_dt' : 1.  }
"""		
hydro_history_filepath = sys.argv[1]
dynamic_config = {  'type'	  : 'dynamic', 
					'hydrofile' : hydro_history_filepath	}
"""
# Physics option
LBT_config = {  'physics'   : 'LBT',
				'2->2'	  : True,
				'2->3'	  : False,
				'3->2'	  : False,
				'Nf'		: 3,
				'mass'	  : 1.3 }  

LGV_config = {  'physics'   : 'LGV',
				'dt_lrf'	: 0.02,
				'elastic'   : True,
				'Einstein'  : True,
				'Nf'		: 3,
				'mass'	  : 1.3 } 

# Initialization option
box_init = {	'type'  : 'box',
				'L'	 : 10.,
				'pmax'  : 10.   }
"""
TR = np.loadtxt(sys.argv[2]).T ## p=0.
TAB = TR**2
realistic_init =  { 'type'		  : 'A+B',
					'sample power'  : 1.,
					'pTmin'		 : 0.1,
					'pTmax'		 : 70.,
					'ymin'		  : -1.,
					'ymax'		  : 1.,
					'TAB'		   : TAB,
					'dxy'		   : 0.1   }
"""
				
e1 = event.event(   medium_flags=static_config, 
					physics_flags=LBT_config   )

e1.initialize_HQ(   NQ=10000,
					init_flags=box_init   )

# Run Model
f = h5py.File("particledata.hdf5", 'w')
Init_pT = e1.Init_pT()
f.create_dataset('Init_PT', data=Init_pT)
for i in range(50):
	print("t = %1.2f [fm/c]"%e1.sys_time() )
	status = e1.perform_hydro_step(StaticPropertyDictionary=box_info)
	if i%1 == 0:
		dsp, dsx = e1.HQ_hist()
		f.create_dataset('p-%d'%i, data=dsp)
		f.create_dataset('x-%d'%i, data=dsx)
	if not status:
		break
f.close()
