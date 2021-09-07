import scipy.signal as sps
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as interp
import copy

def applyLowerPassFilter(
    pulse,  # pulse amplitudes at every time bin
    time,  # list of times for every bin
    time_per_bin,  # how many 10^-10 seconds are there for every bin?
    bin_subdivisions: float,  # how many times we want to subdivide pulse bins to get smooth result
    window_size: float,  # how many ns the window should be
    debug: bool = False,  # should print debug output?
    enforce_zero_at_Ends = False, # return function to zero at ends regardless of pulse shape
    extra_zero_padding = 0, #how many bins need to be padded as zero for the function to go to zero at ends
):
    # time scale is 10^-10 seconds per 1 time unit
    # thus how many bins occur per time will give us the appropriate width for our window

    if enforce_zero_at_Ends: #adds more zero padding
        for l in range(0, extra_zero_padding):
            pulse[l] = 0
            pulse[len(pulse) - l - 1] = 0

    og_pulse = copy.deepcopy(pulse)
    pulse = np.repeat(pulse, bin_subdivisions)
    og_time = copy.deepcopy(time)
    time = np.linspace(np.min(time), np.max(time), len(time) * bin_subdivisions)

    og_window_size = copy.deepcopy(window_size)
    window_size /= time_per_bin  # convert into bin units
    window_size *= bin_subdivisions  # scale with rest of system
    window_size = int(window_size)  # truncate to int

    window = sps.windows.hann(window_size)
    smoothed = sps.convolve(pulse, window, "same")  # zero padding partly included

    # normalize and convert back to original amplitude scale
    smoothed /= np.amax(smoothed)
    smoothed *= np.amax(pulse)

    if debug:
        fig, (ax_orig, ax_win, ax_filt) = plt.subplots(3, 1, sharex=True)
        ax_orig.step(time / 10, pulse)
        ax_orig.set_title("Original pulse")
        ax_orig.margins(0, 0.1)
        ax_win.step(time[0 : len(window)] / 10, window)
        ax_win.set_title("Filter impulse response")
        ax_win.margins(0, 0.1)
        ax_filt.step(time / 10, smoothed)
        ax_filt.set_title("Filtered signal")
        ax_filt.margins(0, 0.1)
        ax_filt.set_xlabel('time (ns)')
        fig.tight_layout()
        # fig.show()

    if enforce_zero_at_Ends and not checkIfPulseGoesToZero(smoothed):
        if debug:
            print("rezeroing pulse with more padding")
        return applyLowerPassFilter(og_pulse, og_time, time_per_bin, bin_subdivisions, og_window_size, debug, enforce_zero_at_Ends, extra_zero_padding + 1)

    return time, smoothed

_tolerance_for_pulse_zero_check = 1e-7
def checkIfPulseGoesToZero(pulse):
    x0 = np.abs(pulse[0])
    xf = np.abs(pulse[-1])
    if x0 > _tolerance_for_pulse_zero_check or xf > _tolerance_for_pulse_zero_check:
        return False
    return True