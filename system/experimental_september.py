import numpy as np
from qutip import *

'''
The system as specified by Mingkang Xia through an email.
'''

# Parameters for Transmon

basis_size = 10

a = destroy(basis_size)

anharm = - 0.20753 * 2 * np.pi # 207.53 MHz 
omega_0 = 3.97117 * 2 * np.pi # 3.97117 GHz

omega = omega_0 + anharm

# ! Hatlab form

H_hl = (omega - 2*anharm)*a.dag()*a + (anharm/6)*((a.dag()+a)**4) # from Mingkang slides 08/31/2021 (linearized)
hatlab_basis_states = H_hl.eigenstates()[1]
hatlab_energies = H_hl.eigenenergies()

# ! standard form
H_0 = omega*a.dag()*a + (anharm/2)*(a.dag()*a)*(a.dag()*a - 1) # from definition of superconducting qubit

H0 = unperturbed_hamiltonian = H_0 # setting up aliases

basis_states = H0.eigenstates()[1]
energies = H0.eigenenergies()
first_energy_spacing = (energies[1] - energies[0])

H_1 = (a.dag() + a)
H1 = perturbation_hamiltonian = H_1 # setting up aliases

# ! full cosine form

E_c = -2*anharm
E_j = (omega_0 + E_c)**2 / (8*E_c)

n_hat = Qobj((E_j / (32 * E_c)) ** (1 / 4) * (a.dag() - a))

theta_hat = Qobj(
    ((2 * E_c) / E_j) ** (1 / 4) * (a.dag() + a)
)

H_R = full_cosine_hamiltonian = -4 * E_c * n_hat ** 2 - E_j * theta_hat.cosm()

full_cosine_basis_states = full_cosine_hamiltonian.eigenstates()[1]
full_cosine_energies = full_cosine_hamiltonian.eigenenergies()