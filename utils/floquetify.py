import copy
import numpy as np
import scipy as sp
import qutip
import qutip.floquet as floq

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    if np.amin(np.abs(array - value)) > 1e-12:
        print("error, duration does not exactly exist!", value, array[idx])
    return idx

def overlap(state1, state2):
    return (np.transpose(np.matrix(state1.full())) * np.matrix(state2.full())).item()

def floquetify(system, solution, pulse, tlist):

    freq = pulse._wd

    # get times on peaks of drive
    dur_points = np.arange(0, tlist[-1], (2*np.pi)/(freq))
    dur_points = dur_points[1:-1]
    durations = tlist

    # generate floquet kets
    kets = []
    energies = []
    for dur in dur_points:
        dur_amp = pulse.pulse_func(dur)
        def steady(t, args=None):
            return dur_amp * np.cos(freq * t)
        H = [system.H0, [system.H1, steady]]

        floquet_kets, quasi_energies = floq.floquet_modes(H, (2*np.pi)/freq)
        kets.append(floquet_kets)
        energies.append(quasi_energies)

    # sort floquet kets
    sorted_kets = [np.array(kets[0]).reshape((10,10))]
    sorted_energies = [np.array(energies[0])]
    for i in range(1, len(dur_points)):
        reshaped_ket = np.array(kets[i]).reshape((10,10))
        m = np.abs(np.dot( (sorted_kets[-1]), np.transpose(reshaped_ket)  ))**2
        row_ind, col_ind = sp.optimize.linear_sum_assignment(1-m)
        # print(col_ind)
        sorted_kets.append( reshaped_ket[col_ind] )
        sorted_energies.append( energies[i][col_ind] )

    # correctly order sorted kets
    expv_for_sorted_unordered_kets = []
    dur = dur_points[0]
    index_of_dur = find_nearest(durations, dur)
    ordered_energies = copy.deepcopy(sorted_energies)
    for k in sorted_kets[0]:
        expv_for_sorted_unordered_kets.append(qutip.expect(qutip.ket2dm(solution.states[index_of_dur]), qutip.Qobj(k)))
    for t in range(0, len(sorted_kets)):
        temp = []
        temp_energies = []
        for i in list(np.flip(np.argsort(expv_for_sorted_unordered_kets))):
            temp.append(sorted_kets[t][i])
            temp_energies.append(sorted_energies[t][i])
        sorted_kets[t] = temp
        ordered_energies[t] = temp_energies
    # sorted_kets = np.array(sorted_kets)[list(np.flip(np.argsort(expv_for_sorted_unordered_kets)))]


    # find overlaps of sorted kets
    overlaps = []
    expectation_values = []
    for i in range(len(dur_points)):
        ket = sorted_kets[i]
        dur = dur_points[i]
        index_of_dur = find_nearest(durations, dur)
        overlaps_at_this_dur = []
        expv_at_this_dur = [] 
        for k in ket:
            expv_at_this_dur.append(qutip.expect(qutip.ket2dm(solution.states[index_of_dur]), qutip.Qobj(k)))
            overlaps_at_this_dur.append(overlap(solution.states[index_of_dur], qutip.Qobj(k)))
        overlaps.append(np.array(overlaps_at_this_dur))
        expectation_values.append(np.array(expv_at_this_dur))

    overlaps = np.transpose(np.array(overlaps))
    expectation_values = np.transpose(np.array(expectation_values))

    return {
        "overlaps": overlaps,
        "expectation_values": expectation_values,
        "ordered_energies": ordered_energies
    }