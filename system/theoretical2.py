import numpy as np
from qutip import *

'''
The system I as specified by Michael Hatridge in our first meeting, but now put into the form used in intro_qubit_driving.
'''

# Parameters for Transmon

basis_size = 10

a = destroy(basis_size)

# We convert from the anharm and omega parameters into E_j and E_c we are familiar with
# anharm = -E_c # https://qiskit.org/textbook/ch-quantum-hardware/transmon-physics.html#4.-The-Quantized-Transmon-
# omega_not = np.sqrt(8*E_c*E_j)
anharm = - 0.19 # -190 MHz
omega_not = 4.54 # 4.54 GHz
E_c = -1*anharm
E_j = (omega_not ** 2) / (8*E_c)

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
eigenvectors = H0.eigenstates()[1]

a = destroy(N=basis_size)
perturbation = Qobj(a.dag() + a)

basis_states = eigenvectors
fock_states = [(basis(basis_size, 0)), (basis(basis_size, 1)), (basis(basis_size, 2))]

H1 = perturbation
first_energy_spacing = energy_spacing = energies[1] - energies[0]