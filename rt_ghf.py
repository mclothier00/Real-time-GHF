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
        #mo_oth_old = []
        shape = self.total_steps / self.frequency
        #mag_x = np.empty(shape = shape)
        #mag_y = np.empty(shape = shape)
        #mag_z = np.empty(shape = shape)
        #time = np.empty(shape = shape)
        mag_x = []
        mag_y = []
        mag_z = []
        time = []

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
             #   ener_tot = self._scf.energy_tot()
             #   print(F'Energy at step {i}: {ener_tot}')
                
                # 100% hard-coded for the h-atom example using values for BSE, need to change to iterate over matrix elements 
                den = self._scf.make_rdm1()
                mag_x.append((den[0,1] + den[1,0]) * 3.425250914) #* 0.1543289673)
                mag_y.append(1j * (den[0,1] - den[1,0]) * 3.425250914)# * 0.1543289673)
                mag_z.append((den[0,0] - den[1,1]) * 3.425250914)# * 0.1543289673)
                time.append(i)

            mo_oth_old = mo_oth

        ### graph magentization results (need to move out of class)
        mag_x = np.array(mag_x)
        mag_y = np.array(mag_y)
        mag_z = np.array(mag_z)
        print(mag_x)
        print(mag_y)
        print(mag_z)
        time = np.array(time)
        plt.plot(time, mag_x, 'r')#time, mag_y, 'b', time, mag_z, 'g')
#time, mag_x, 'r',
        plt.savefig('mag_z.png')
