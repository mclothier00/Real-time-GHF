import numpy as np
from scipy import integrate
from scipy.linalg import eigh, inv
from pyscf import gto, scf
import scipy

### initialize variables
timestep = 0.1
steps = 3
total_steps = 30

### initialize static calculation
distance = 1.1
mol = gto.M(        
	verbose = 0,       
	atom=(F'H 0 0 0; H 0 0 {distance}'),        
	basis='6-31g')
#        basis='cc-pvdz')

mf = scf.ghf.GHF(mol)
mf.scf()

### creating initial core hamiltonian 
#hcore = scf.hf.get_hcore(mol)
#hcore = scipy.linalg.block_diag(hcore, hcore)
hcore = mf.get_hcore()

den = mf.make_rdm1() # need to modify this in the future
ovlp = mf.get_ovlp()
x = scf.addons.canonical_orth_(ovlp)
fock = hcore + mf.get_veff(mol, den)
mo_old = mf.mo_coeff
nuc_e = mf.mol.energy_nuc()
mo_oth_old = []

####### DYNAMICS #######
for i in range(0, total_steps):
   ### transforming coefficients into an orthogonal matrix (step 10)
   mo_oth = np.dot(inv(x), mf.mo_coeff)

   ### create transformation matrix U from Fock matrix at time t (step 11)
   fock_oth = np.dot(x.T, np.dot(fock, x))
   u = scipy.linalg.expm(-1j*2*timestep*fock_oth)

   ### propagate MO coefficients (step 12)
   if i != 0:
       mo_oth_new = np.dot(u, mo_oth_old)
   else:
       mo_oth_new = np.dot(u, mo_oth)

   ### transform coefficients back into non-orthogonal basis and get density matrix
   mf.mo_coeff = np.dot(x, mo_oth_new)
   mo = mf.mo_coeff
   den = mf.make_rdm1(mo) #check what happens if dont specify the mos as input

   # calculate a new fock matrix
   fock = hcore + mf.get_veff(mol, den)

   # calculate energy and other observables
   if np.mod(i, steps)==0:
      elec_e = mf.energy_elec()
      print(F'Energy at step {i}: {elec_e + nuc_e}')

   mo_oth_old = mo_oth
