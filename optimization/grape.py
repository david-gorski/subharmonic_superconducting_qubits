import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import copy
import importlib
import optimization.lowpass as lp
import utils.solver as solver
import utils.expectation_values as expect
import utils.graph as g
import utils.saver as saver
importlib.reload(lp)
importlib.reload(solver)
importlib.reload(expect)
importlib.reload(g)


'''
GRadient Ascent Pulse Engineering (GRAPE)

1. Create Discrete Pulse
2. Integrate through pulse duration
3. Optimize parameters for fidelity

'''

class _GRAPE_Config():
    pass

config = _GRAPE_Config()

def setup(system, pulse, offset, total_pulse_time, bin_size, smoothing_time, bin_lowpass_subdivisions=1000):
    global config
    config.pulse_t_slots = int( (total_pulse_time / bin_size) )
    config.padding_t_slots = int( ((smoothing_time*2) / bin_size) )
    config.num_t_slots = config.pulse_t_slots + config.padding_t_slots
    print("num_t_slots", config.num_t_slots)
    config.total_pulse_time = total_pulse_time
    config.bin_size = bin_size
    config.bin_lowpass_subdivisions = bin_lowpass_subdivisions
    config.amp_polling_rate = bin_size
    config.smoothing_time = smoothing_time
    config.time = np.linspace(0 , total_pulse_time+(smoothing_time*2), config.num_t_slots)
    config.system = system
    config.starting_state = system.starting_state
    config.target_state = system.target_state
    
    config.original_offset = offset

    config.original_unpadded_pulse = copy.deepcopy(pulse)
    config.original_pulse = np.array(list(np.zeros(int(config.padding_t_slots/2))) + list(copy.deepcopy(pulse)) + list(np.zeros(int(config.padding_t_slots/2))))
    if len(config.original_pulse) != config.num_t_slots:
        config.original_pulse = np.append(config.original_pulse, 0) # solve off by 1 error by adding padding
    smoothed_time, config.original_pulse_smoothed = lp.applyLowerPassFilter(config.original_pulse, config.time, config.bin_size, config.bin_lowpass_subdivisions, smoothing_time, debug = False, enforce_zero_at_Ends = False, extra_zero_padding = 0)
    config.smoothed_time = smoothed_time
    config.best_fidelity = 0

    config.last_offset = copy.deepcopy(config.original_offset)
    config.last_smooth_pulse = copy.deepcopy(config.original_pulse_smoothed)
    config.last_pulse = copy.deepcopy(config.original_pulse)

    config.iterations = 0

def cost(parameters, constants=None):
    global config, _current_smooth_pulse, _current_drive_frequency, _smooth_time
    [unsmoothed_pulse, other_params] = np.split(parameters, [config.num_t_slots])
    _smooth_time, _current_smooth_pulse = lp.applyLowerPassFilter(unsmoothed_pulse, config.time, config.bin_size, config.bin_lowpass_subdivisions, config.smoothing_time, debug = False, enforce_zero_at_Ends = False, extra_zero_padding = 0)

    offset = other_params[0]
    _current_drive_frequency = config.system.first_energy_spacing/3 - offset

    time = np.linspace(0, config.time[-1], int(config.time[-1]*10))
    solution = solver.time_evolve(H_d=config.system.H0, H_c=config.system.H1, pertubation_func=perturbation_func, tlist=time, starting_state = config.starting_state, args={}, pbar=False, store_states=False)
    
    fidelity = expect.expectation_value(solution.final_state, config.target_state)
    figure_of_merit = (1-fidelity)
    config.last_solution = solution
    config.last_offset = offset
    config.last_smooth_pulse = _current_smooth_pulse
    config.last_pulse = unsmoothed_pulse
    if fidelity > config.best_fidelity:
        config.best_fidelity = fidelity
        config.best_pulse_smoothed = _current_smooth_pulse
        config.best_pulse = unsmoothed_pulse
        config.best_drive_frequency = _current_drive_frequency
        config.best_offset = offset
        config.best_solution = solution
    config.iterations += 1
    print("%i : %f" % (config.iterations, fidelity))
    return figure_of_merit

def single_run():
    global config, _current_smooth_pulse, _current_drive_frequency, _smooth_time
    parameters = list(config.last_pulse)
    parameters.append(config.last_offset)
    parameters = np.array(parameters)
    constants = None
    [unsmoothed_pulse, other_params] = np.split(parameters, [config.num_t_slots])
    _smooth_time, _current_smooth_pulse = lp.applyLowerPassFilter(unsmoothed_pulse, config.time, config.bin_size, config.bin_lowpass_subdivisions, config.smoothing_time, debug = False, enforce_zero_at_Ends = False, extra_zero_padding = 0)

    offset = other_params[0]
    _current_drive_frequency = config.system.first_energy_spacing/3 - offset

    time = np.linspace(0, config.time[-1], int(config.time[-1]*10))
    solution = solver.time_evolve(H_d=config.system.H0, H_c=config.system.H1, pertubation_func=perturbation_func, tlist=time, starting_state = config.starting_state, args={}, pbar=False, store_states=False)
    config.last_solution = solution
    config.last_offset = offset
    config.last_smooth_pulse = _current_smooth_pulse
    config.last_pulse = unsmoothed_pulse

    return solution

_current_drive_frequency = 0
_current_smooth_pulse = []
_smooth_time = []
def perturbation_func(t, args=None):
    return (nearest(t, _smooth_time, _current_smooth_pulse) * np.cos(_current_drive_frequency * t))
    
def nearest(t, ts, vs):
    return vs[(np.abs(ts - t)).argmin()]

def optimize():
    global config
    parameters = list(config.last_pulse)
    parameters.append(config.last_offset)
    parameters = np.array(parameters)
    constants = None
    result1 = sp.optimize.minimize(cost, parameters, constants, method="Nelder-Mead", options={"disp":True})

    config.optimization_result = result1
    return result1

def graph_last_solution():
    global config
    expectation_values = expect.get_all_expectation_values(config.last_solution.states, config.system.basis_states)

    print("Fidelity = %f, Offset = %f" % (expectation_values[1][-1], config.last_offset))
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.set_xlabel("duration (ns)")
    ax.set_ylabel("fidelity")
    ax.set_title("Expectation Values")
    for i in range(0, len(expectation_values)):
        ax.plot(config.last_solution.times, expectation_values[i], label="$\psi_{%i}$" % i)
    plt.legend()
    plt.show()

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.set_xlabel("duration (ns)")
    ax.set_ylabel("fidelity")
    ax.set_title("Pulse")
    ax.plot(config.smoothed_time, config.last_smooth_pulse, label="smoothed")
    ax.plot(config.time, config.last_pulse, label="original")
    real_pulse = np.vectorize(perturbation_func)(config.smoothed_time)
    ax.plot(config.smoothed_time, real_pulse, label="actual")
    plt.legend()
    plt.show()
