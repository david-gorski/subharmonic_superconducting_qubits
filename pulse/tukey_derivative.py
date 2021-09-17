import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sps

_width = 0
_ramp_coef = 0
_t_max = 0
_t_min = 0
_amp = 0
_wd = 0

# simple analytic derivative of tukey

def pulse_func(t, args=None):
    if t > _t_max:
        return 0
    elif t < _t_min:
        return 0
    elif t > _t_max-_width:
        # ramping down
        return (_amp * np.pi * 1/_ramp_coef * 1/_t_max) * (-np.sin(np.pi * (-2.0/_ramp_coef + 1 + 2.0*t/_ramp_coef/(_t_max)))) * np.cos(t*_wd)
    elif t > _width:
        # middle flat top
        return 0
    else:
        # ramping up
        return (_amp * np.pi * 1/_ramp_coef * 1/_t_max) * (-np.sin(np.pi * (-1 + 2*t/_ramp_coef/_t_max))) * np.cos(t*_wd)

def setup(amplitude, drive_frequency, ramp_coef, tlist):
    global _ramp_coef, _amp, _wd, _t_max, _width
    _ramp_coef = ramp_coef
    _amp = amplitude
    _wd = drive_frequency
    _t_max = tlist[-1]
    _t_min = tlist[0]
    _width = (_ramp_coef*(_t_max)/2.0)

def get_pulse(tlist):
    return np.vectorize(pulse_func)(tlist)
