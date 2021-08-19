import numpy as np
from qutip import *

'''
The system I had used previously in intro_qubit_driving, but now reparameterized so that the hamiltonian is diagonal.
'''

# Parameters for Transmon

basis_size = 10

# transmon regime is where E_j / E_c >> 1
E_j = 51.486
E_c = 0.79168

c = destroy(basis_size)

anharm = -E_c # https://qiskit.org/textbook/ch-quantum-hardware/transmon-physics.html#4.-The-Quantized-Transmon-
omega_not = np.sqrt(8*E_c*E_j)

omega = omega_not + anharm

# omega_j = (omega - anharm/2)*j + (anharm/2)*(j**2)

H_0 = omega*c.dag()*c + (anharm/2)*(c.dag()*c)*(c.dag()*c - 1)
H0 = unperturbed_hamiltonian = H_0 # setting up aliases

basis_states = H0.eigenstates()[1]
energies = H0.eigenenergies()
first_energy_spacing = (energies[1] - energies[0])

H_1 = (c.dag() + c)
H1 = perturbation_hamiltonian = H_1 # setting up aliases