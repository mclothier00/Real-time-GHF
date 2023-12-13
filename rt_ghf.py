import numpy as np
from scipy.linalg import eigh, inv
from pyscf import gto, scf
import scipy
import matplotlib.pyplot as plt

# class needs mf, timestep, frequency, total_steps

class GHF:
    def __init__(self, mf, timestep, frequency, total_steps, orth=None):
        self.timestep = timestep
        self.frequency = frequency 
        self.total_steps = total_steps
        self._scf = mf
    
        if orth is None: self.orth = scf.addons.canonical_orth_(self._scf.get_ovlp())
        # orth is indexed not in terms of alpha and beta but in terms of eigenvalues; will not be block diagonal
 
    # may have to convert mo coefficients to complex array if they are real from pyscf 

    ####### DYNAMICS #######
    def dynamics(self):
        ### creating initial core hamiltonian
#        den = self._scf.make_rdm1() # need to modify this in the future
        fock = self._scf.get_fock()
        mag_x = []
        mag_y = []
        mag_z = []
        time = []
        energy = []

        ovlp = self._scf.get_ovlp()
        Nsp = int(ovlp.shape[0]/2)

        for i in range(0, self.total_steps):
            ### transforming coefficients into an orthogonal matrix (step 10)
            mo_oth = np.dot(inv(self.orth), self._scf.mo_coeff)

            ### create transformation matrix U from Fock matrix at time t (step 11)
            fock_oth = np.dot(self.orth.T, np.dot(fock, self.orth))
            u = scipy.linalg.expm(-1j*2*self.timestep*fock_oth)

            ### propagate MO coefficients (step 12)
            if i != 0:
                mo_oth_new = np.dot(u, mo_oth_old)
            else:
                mo_oth_new = np.dot(u, mo_oth)

            ### transform coefficients back into non-orthogonal basis and get density matrix
            self._scf.mo_coeff = np.dot(self.orth, mo_oth_new)
            den = self._scf.make_rdm1() 

            # calculate a new fock matrix
            fock = self._scf.get_fock()

            # calculate energy and other observables
            if np.mod(i, self.frequency)==0:
                ener_tot = self._scf.energy_tot()
                energy.append(ener_tot)
            #    print(F'Energy at step {i}: {ener_tot}')
                
                den = self._scf.make_rdm1()
                mag_x.append(np.sum(den[:Nsp,Nsp:] + den[Nsp:,:Nsp])) 
                mag_y.append(1j * (np.sum(den[:Nsp,Nsp:] - den[Nsp:,:Nsp])))
                mag_z.append(np.sum(den[:Nsp,:Nsp] - den[Nsp:,Nsp:]))

                t = (i * self.timestep) / 41341.374575751
                time.append(t)
                
            #    print(f'completed step {i}')                
            mo_oth_old = mo_oth

        ### graph magentization results (need to move out of class)
        mag_x = np.array(mag_x)
        mag_y = np.array(mag_y)
        mag_z = np.array(mag_z)
        print(mag_x)
        print(mag_y)
        print(mag_z)
        #print(energy)
        time = np.array(time)
        plt.plot(time, np.real(mag_x), 'r',time, np.real(mag_y), 'b', time, np.real(mag_z), 'g')
        plt.savefig('li_mag.png')
