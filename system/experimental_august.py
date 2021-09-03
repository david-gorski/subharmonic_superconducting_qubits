import numpy as np
from qutip import *

'''
The system as specified by Mingkang Xia through an email.
'''

# Parameters for Transmon

basis_size = 10

# transmon regime is where E_j / E_c >> 1
# E_j = 51.486
# E_c = 0.79168

c = destroy(basis_size)

# in units of MHz
# ? anharm should be negative?
anharm = - 0.20753 * 2 * np.pi # 207.53 MHz 
omega_not = 3.97117 * 2 * np.pi # 3.97117 GHz

omega = omega_not + anharm

# omega_j = (omega - anharm/2)*j + (anharm/2)*(j**2)

H_0 = omega*c.dag()*c + (anharm/2)*(c.dag()*c)*(c.dag()*c - 1) #order of operators here is probably off from Michael's
H0 = unperturbed_hamiltonian = H_0 # setting up aliases

basis_states = H0.eigenstates()[1]
energies = H0.eigenenergies()
first_energy_spacing = (energies[1] - energies[0])

H_1 = (c.dag() + c)
H1 = perturbation_hamiltonian = H_1 # setting up aliases