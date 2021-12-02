import qutip.ui.progressbar as progress_bar
from qutip import *
import utils.solver as solver

_atol = 1e-10
_rtol = 1e-10
_norm_tol = 1e-10

def time_evolve(
    H_d,
    H_c, # only single control is allowed
    pertubation_func,
    pertubation_ramp_up_time,
    pertubation_ramp_down_time,
    pertubation_period,
    tlist,
    starting_state,
    args={},
    pbar=False,
    store_states=True
):  
    ramp_up_tlist = []
    ramp_down_tlist = []
    for t in tlist:
        if t < pertubation_ramp_up_time:
            ramp_up_tlist.append(t)
        if t > pertubation_ramp_down_time:
            ramp_down_tlist.append(t)
    
    ramped_up_solution = solver.time_evolve(H_d, H_c, pertubation_func, ramp_up_tlist, starting_state, args, pbar, store_states)
    

    def steady_state_hamiltonian(t, args):
        return H_c*pertubation_func(t+pertubation_ramp_up_time)

    f_modes_0, f_energies = floquet_modes(steady_state_hamiltonian, pertubation_period, args)
    f_coeff = floquet_state_decomposition(f_modes_0, f_energies, ramped_up_solution.final_state) # initial state in terms of floquet modes

    f_modes_t = floquet_modes_t(f_modes_0, f_energies, ramp_down_tlist[-1], steady_state_hamiltonian, pertubation_period, args)

    psi_at_end_of_steady_state = floquet_wavefunction_t(f_modes_0, f_energies, f_coeff, ramp_down_tlist[-1], steady_state_hamiltonian, pertubation_period, args)

    ramped_down_solution = solver.time_evolve(H_d, H_c, pertubation_func, ramp_down_tlist, psi_at_end_of_steady_state, args, pbar, store_states)

    return ramped_down_solution

def set_tolerances(atol=1e-10, rtol=1e-10, norm_tol=1e-10):
    global _atol, _rtol, _norm_tol
    _atol = atol
    _rtol = rtol
    _norm_tol = norm_tol

    solver._atol = _atol
    solver._rtol = _rtol
    solver._norm_tol = _norm_tol
