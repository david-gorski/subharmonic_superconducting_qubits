import numpy as np
from qutip import *

'''
The truncated hamiltonian with adjustable parameters
'''

def setup(_omega, _anharm, _basis_size=10):
    global omega, anharm, H0, energies, basis_states, first_energy_spacing, starting_state, target_state, H1, H2, basis_size, a
    basis_size = _basis_size
    omega = _omega
    anharm = _anharm
    a = destroy(basis_size)

    H0 = (omega)*a.dag()*a + (anharm/2)*(a.dag()*a)*(a.dag()*a - 1) # from definition of superconducting qubit

    basis_states = H0.eigenstates()[1]
    energies = H0.eigenenergies()
    first_energy_spacing = (energies[1] - energies[0])

    starting_state = basis_states[0]
    target_state = basis_states[1]

    H_1 = H1 = x_perturbation = perturbation_hamiltonian = (a.dag() + a)
    H2 = y_perturbation = H_2 = -(1j)*(a-a.dag())

# this data comes from Mingkang in September
setup(
    _omega =  3.97117 * 2 * np.pi, # 3.97117 GHz transition frequency between ground and excited 
    _anharm = - 0.20753 * 2 * np.pi, # 207.53 MHz 
)