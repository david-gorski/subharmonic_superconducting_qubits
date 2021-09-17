import numpy as np
from qutip import *

'''
The un-truncated hamiltonian with adjustable parameters
'''

def setup(_E_c, _E_j, _basis_size=10):
    global E_c, E_j, H0, energies, basis_states, first_energy_spacing, starting_state, target_state, H1, H2, basis_size, a
    basis_size = _basis_size
    E_c = _E_c
    E_j = _E_j
    a = destroy(basis_size)
    n_hat = Qobj((E_j / (32 * E_c)) ** (1 / 4) * (a.dag() - a))
    theta_hat = Qobj(
        ((2 * E_c) / E_j) ** (1 / 4) * (a.dag() + a)
    )
    H0 = -4 * E_c * n_hat ** 2 - E_j * theta_hat.cosm()
    basis_states = H0.eigenstates()[1]
    energies = H0.eigenenergies()
    first_energy_spacing = (energies[1] - energies[0])

    starting_state = basis_states[0]
    target_state = basis_states[1]

    H_1 = H1 = x_perturbation = perturbation_hamiltonian = (a.dag() + a)
    H2 = y_perturbation = H_2 = -(1j)*(a-a.dag())

# this data is from the find_params notebook
# this is pretuned to most closely match experimental_september
setup(1.15849102, 73.88873718)