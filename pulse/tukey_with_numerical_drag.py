import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sps

_width = 0
_ramp_coef = 0
_t_max = 0
_t_min = 0
_amp = 0
_wd = 0
_drag_coef = 1
_last_t = 0
_last_v = 0

def pulse_func(t, args=None):
    global _last_t, _last_v
    if t > _t_max:
        v = 0
    elif t < _t_min:
        v = 0
    elif t > _t_max-_width:
        # ramping down
        v = _amp * 0.5 * (1 + np.cos(np.pi * (-2.0/_ramp_coef + 1 + 2.0*t/_ramp_coef/(_t_max))))
    elif t > _width:
        # middle flat top
        v = _amp * 1
    else:
        # ramping up
        v = _amp * 0.5 * (1 + np.cos(np.pi * (-1 + 2.0*t/_ramp_coef/(_t_max))))

    # add drag term
    if((t - _last_t) == 0):
        two_point_derivative = 0
    else:
        two_point_derivative = (v - _last_v)/(t - _last_t)
    _last_t = t
    _last_v = v

    return (v+_drag_coef*two_point_derivative) * np.cos(t*_wd)


def setup(amplitude, drive_frequency, ramp_coef, drag_coef, tlist):
    global _ramp_coef, _amp, _wd, _t_max, _width, _drag_coef
    _ramp_coef = ramp_coef
    _drag_coef = drag_coef
    _amp = amplitude
    _wd = drive_frequency
    _t_max = tlist[-1]
    _t_min = tlist[0]
    _width = (_ramp_coef*(_t_max)/2.0)

def get_pulse(tlist):
    # pulse = []
    # for t in tlist:
    #     pulse.append(pulse_func(t))
    # return pulse
    return np.vectorize(pulse_func)(tlist)
