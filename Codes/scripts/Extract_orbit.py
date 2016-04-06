#==================================================================================================
#Importing all libraries
#==================================================================================================
import matplotlib.pyplot as plt
import numpy as np
import numpy as np
import os
import h5py

from scripts.dynamical_friction import *

#Data folder
DataFolder = '/home/bustamsn/PhD/Data/hernquist_dynamical_friction/' 


#==================================================================================================
#Extracting
#==================================================================================================
#hernq_kick_BH_1e5 ----------------------------------------
#Potential centered
bh_s1 = black_hole_sim( simulation = "%s/hernq_kick_BH_1e5"%(DataFolder), 
                        snapbase = 'snapshot',
                        n_snap = 240,
                        center = [250,250,250] )
bh_s1.read_trajectory( pot_center = True )
np.savetxt( "../results/BH_hernq_kick_BH_1e5_p.dat", np.transpose( [bh_s1.t, bh_s1.r[:,0], bh_s1.r[:,1], bh_s1.r[:,2]] ) )

#Geometrically centered
bh_s1_g = black_hole_sim( simulation = "%s/hernq_kick_BH_1e5"%(DataFolder), 
                          snapbase = 'snapshot',
                          n_snap = 240,
                          center = [250,250,250] )
bh_s1_g.read_trajectory( )
np.savetxt( "../results/BH_hernq_kick_BH_1e5_g.dat", np.transpose( [bh_s1_g.t, bh_s1_g.r[:,0], bh_s1_g.r[:,1], bh_s1_g.r[:,2]] ) )


#hernq_kick_BH_1e6 ----------------------------------------
#Potential centered
bh_s2 = black_hole_sim( simulation = "%s/hernq_kick_BH_1e6"%(DataFolder), 
                        snapbase = 'snapshot',
                        n_snap = 240,
                        center = [250,250,250] )
bh_s2.read_trajectory( pot_center = True )
np.savetxt( "../results/BH_hernq_kick_BH_1e6_p.dat", np.transpose( [bh_s2.t, bh_s2.r[:,0], bh_s2.r[:,1], bh_s2.r[:,2]] ) )

#Geometrically centered
bh_s2_g = black_hole_sim( simulation = "%s/hernq_kick_BH_1e6"%(DataFolder), 
                          snapbase = 'snapshot',
                          n_snap = 240,
                          center = [250,250,250] )
bh_s2_g.read_trajectory( )
np.savetxt( "../results/BH_hernq_kick_BH_1e6_g.dat", np.transpose( [bh_s2_g.t, bh_s2_g.r[:,0], bh_s2_g.r[:,1], bh_s2_g.r[:,2]] ) )