import numpy as np
from pyscf import gto, scf
import rt_ghf

### initialize variables
timestep = 0.1
steps = 1000
total_steps = 40000

mag_z = 0.000085 # in au

### initialize static calculation
mol = gto.M(        
	verbose = 0,       
	atom='H 0 0 0',        
	basis='STO-3G',
    spin = 1)

#mol.incore_anyway(True)

### initial calculation to get t=0 hamiltonian
mf = scf.ghf.GHF(mol)
mf.scf()
# basis_func = [3.425250914, 0.1543289673] # pulled off of basis set exchange, not sure how to pull functions from pyscf

### initialize hamiltonian

# initialized with RHF
#arr_mag = np.zeros((2, 2))
#arr_h = np.zeros((2,2))
#mag = [0.5 * mag_z, -0.5 * mag_z]
#mag = np.array(mag)
#np.fill_diagonal(arr_mag, mag)
#h_core = mf.get_hcore()
#np.fill_diagonal(arr_h, h_core)
#h_mag = arr_h + arr_mag
#print(h_mag)

# initialized with GHF
arr = np.zeros((2,2))
mag = np.array([0.5 * mag_z, -0.5 * mag_z])
np.fill_diagonal(arr, mag)
h_mag = mf.get_hcore() + arr
print(h_mag)


### call dynamics 

#mf = scf.ghf.GHF(mol)
#mf.scf()

mf.get_hcore = lambda *args: h_mag

var = rt_ghf.GHF(mf, timestep, steps, total_steps)

var.dynamics()



































