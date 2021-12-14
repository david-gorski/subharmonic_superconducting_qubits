import qutip.ui.progressbar as progress_bar
from qutip import *
import utils.solver as solver
import numpy as np

_atol = 1e-10
_rtol = 1e-10
_norm_tol = 1e-10

def time_evolve(
    H_d,
    H_c, # only single control is allowed
    pulse,
    # perturbation_steady_state_func, # assuming that it will be the standard cosine form
    perturbation_ramp_up_time,
    perturbation_ramp_down_time,
    perturbation_period,
    tlist,
    starting_state,
    pbar=False
):  
    amp = pulse._amp
    freq = pulse._wd

    peaks_tlist = np.arange(0, tlist[-1], perturbation_period)
    ramp_up_tlist = []
    steady_state_tlist = []
    ramp_down_tlist = []
    for t in tlist:
        if t < perturbation_ramp_up_time:
            ramp_up_tlist.append(t)
        elif t > tlist[-1] - perturbation_ramp_down_time:
            ramp_down_tlist.append(t)
        else:
            steady_state_tlist.append(t)

    for t in peaks_tlist:
        if t > ramp_up_tlist[-1]:
            ramp_up_tlist.append(t)
            steady_state_tlist.insert(0, t)
            break
    ramp_up_tlist = np.linspace(ramp_up_tlist[0], ramp_up_tlist[-1], len(ramp_up_tlist)*2) # was getting errors because not enough points

    for i in range(len(peaks_tlist)):
        if peaks_tlist[i] > ramp_down_tlist[0]:
            ramp_down_tlist.insert(0, peaks_tlist[i-1])
            steady_state_tlist.append(peaks_tlist[i-1])
            break
    ramp_down_tlist = np.linspace(ramp_down_tlist[0], ramp_down_tlist[-1], len(ramp_down_tlist)*2) # was getting errors because not enough points

    ramp_up_pulse = []
    def ramp_up_pulse_func(t, args):
        answer = pulse.pulse_func(t)
        ramp_up_pulse.append(answer)
        return answer

    ramped_up_solution = solver.time_evolve(H_d, H_c, ramp_up_pulse_func, ramp_up_tlist, starting_state, {}, pbar, False)
        
    steady_state_tlist_start = steady_state_tlist[0]
    steady_state_tlist = np.array(steady_state_tlist) - steady_state_tlist_start

    def perturbation_steady_state_func(t, args):
        return amp * np.cos((t+steady_state_tlist_start) * freq)

    steady_state_hamiltonian = [H_d, [H_c, perturbation_steady_state_func]]
    # steady_state_solution = fsesolve(steady_state_hamiltonian, ramped_up_solution.final_state, steady_state_tlist, T=pertubation_period)
    # psi_at_end_of_steady_state = steady_state_solution.states[-1]

    f_modes_0, f_energies = floquet_modes(steady_state_hamiltonian, perturbation_period, {})
    f_coeff = floquet_state_decomposition(f_modes_0, f_energies, ramped_up_solution.final_state) # initial state in terms of floquet modes

    # f_modes_t = floquet_modes_t(f_modes_0, f_energies, ramp_down_tlist[0], steady_state_hamiltonian, pertubation_period, {})

    # only need to evolve the state for the time it is in steady state, irrespective of the absolute time of the pulse
    psi_at_end_of_steady_state = floquet_wavefunction_t(f_modes_0, f_energies, f_coeff, steady_state_tlist[-1], steady_state_hamiltonian, perturbation_period, {})

    ramp_down_tlist_start = ramp_down_tlist[0]
    ramp_down_tlist = np.array(ramp_down_tlist) - ramp_down_tlist_start
    ramp_down_pulse = []
    ramp_down_pulse_tlist = []
    def ramp_down_pulse_func(t, args):
        answer = pulse.pulse_func(t+ramp_down_tlist_start)
        ramp_down_pulse_tlist.append(t)
        ramp_down_pulse.append(answer)
        return answer

    ramped_down_solution = solver.time_evolve(H_d, H_c, ramp_down_pulse_func, ramp_down_tlist, psi_at_end_of_steady_state, {}, pbar, False)

    return ramped_down_solution

def set_tolerances(atol=1e-10, rtol=1e-10, norm_tol=1e-10):
    global _atol, _rtol, _norm_tol
    _atol = atol
    _rtol = rtol
    _norm_tol = norm_tol

    solver._atol = _atol
    solver._rtol = _rtol
    solver._norm_tol = _norm_tol
