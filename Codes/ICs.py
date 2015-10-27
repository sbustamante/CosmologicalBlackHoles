import numpy as np
import h5py

hdf5_file_name = '/home/bustamsn/PhD/CosmologicalBlackHoles/Data/MW_binary_BH.hdf5'

#Parameters of the BH binary system
d = 0.01                        # kpc
vc = 10                         # km/s


#Duplicating BHs information
with h5py.File(hdf5_file_name, 'r+') as f:
    for key in f['PartType5'].keys():
        X = list(f['PartType5'][key][:])
        X = np.array(X+X)
        #Deleting current information
        del f['PartType5'][key]
        dset = f['PartType5'].create_dataset(key, data=X)

    #Setting new coordinates and velocities (Circular orbit)
    f['PartType5']['Coordinates'][0] = np.array([500+d/2.0,500,500])
    f['PartType5']['Coordinates'][1] = np.array([500-d/2.0,500,500])
    f['PartType5']['Velocities'][0] = np.array([0,0+vc/2.0,0])
    f['PartType5']['Velocities'][1] = np.array([0,0-vc/2.0,0])
    f['PartType5']['ParticleIDs'][1] = 2900002




