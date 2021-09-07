import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sps

_c = 0 #width
_d = 0 #sharpness
_t_max = 0
_amp = 0
_wd = 0
_normalization = 0

# ! Not fully tested

def pulse_func(t, args=None):
    return (
        _amp
        * _normalization
        * (
            np.cos(_wd * t)
            * np.abs(
                np.exp(-((t - _c / 2) ** 2) / (2 * (_c / _d) ** 2))
                - np.exp(-((_c / 2) ** 2) / (2 * (_c / _d) ** 2))
            )
        )
    )


def setup(amplitude, drive_frequency, sharpness, tlist):
    global _d, _c, _amp, _wd, _t_max, _width, _normalization
    _d = sharpness
    _wd = drive_frequency
    _t_max = tlist[-1]
    _c = _t_max
    _normalization = 1
    _amp = 1 #needs to be 1 for normalization run
    _normalization = 1/np.amax(get_pulse(tlist))
    _amp = amplitude

def get_pulse(tlist):
    return np.vectorize(pulse_func)(tlist)
