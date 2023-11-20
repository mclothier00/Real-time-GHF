from pyscf import gto, scf
import rt_ghf

### initialize variables
timestep = 0.1
steps = 3
total_steps = 12

### initialize static calculation
mol = gto.M(        
	verbose = 0,       
	atom=(F'H 0 0 0'),        
	basis='STO-3G',
    spin = 1)

#mol.incore_anyway(True)

### initialize hamiltonian




mf = scf.ghf.GHF(mol)
mf.scf()

var = rt_ghf.GHF(mf, timestep, steps, total_steps)

var.dynamics()





































