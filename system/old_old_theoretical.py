import numpy as np
from qutip import *

'''
The system I had used previously in intro_qubit_driving.
'''

# Parameters for Transmon

basis_size = 10

# setup unperturbed hamiltonian

# transmon regime is where E_j / E_c >> 1
E_j = 51.486
E_c = 0.79168


def unperturbed_hamiltonian(E_c, E_j, basis_size=basis_size):
    a = destroy(N=basis_size)
    n_hat = Qobj((E_j / (32 * E_c)) ** (1 / 4) * (a.dag() - a))  # should have leading i
    theta_hat = Qobj(
        ((2 * E_c) / E_j) ** (1 / 4) * (a.dag() + a)
    )  # should have leading i

    UnperturbedHamiltonian = (
        -4 * E_c * n_hat ** 2 - E_j * theta_hat.cosm()
    )  # cosm is the cosine function
    return Qobj(UnperturbedHamiltonian)


H0 = unperturbed_hamiltonian(
    E_c, E_j, basis_size
)  # our hamiltonian with the pertubation
energies = H0.eigenenergies()


a = destroy(N=basis_size)
perturbation = Qobj(a.dag() + a)

basis_states = [(basis(basis_size, 0)), (basis(basis_size, 1)), (basis(basis_size, 2))]

H1 = perturbation
first_energy_spacing = energy_spacing = energies[1] - energies[0]