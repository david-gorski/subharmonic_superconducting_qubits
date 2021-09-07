import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sps

_width = 0
_amp = 0
_wd = 0
_k = 0
_b = 0

# based on Mingkang's mathematica notebook

def pulse_func(t, args=None):
    return _amp * ( np.tanh(_k*t-_b) - np.tanh(_k*(t-_width)+_b) ) / 2 * np.cos(_wd*t)

def setup(amplitude, drive_frequency, ramp_slope, cut_factor, tlist):
    global _k, _amp, _wd, _b, _width
    _k = ramp_slope
    _b = cut_factor
    _amp = amplitude
    _wd = drive_frequency
    _width = tlist[-1]

def get_pulse(tlist):
    return np.vectorize(pulse_func)(tlist)
