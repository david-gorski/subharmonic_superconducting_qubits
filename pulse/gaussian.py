import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sps

_c = 0 #width
_d = 0 #sharpness
_t_max = 0
_amp = 0
_wd = 0
_normalization = 0

# based on https://github.com/scipy/scipy/blob/v1.7.1/scipy/signal/windows/windows.py#L795-L875

# def pulse_func(t, args=None):
#     if t > _t_max:
#         return 0
#     elif t > _t_max-_width:
#         # ramping down
#         return _amp * 0.5 * (1 + np.cos(np.pi * (-2.0/_ramp_coef + 1 + 2.0*t/_ramp_coef/(_t_max)))) * np.cos(t*_wd)
#     elif t > _width:
#         # middle flat top
#         return _amp * 1 * np.cos(t*_wd)
#     else:
#         # ramping up
#         return _amp * 0.5 * (1 + np.cos(np.pi * (-1 + 2.0*t/_ramp_coef/(_t_max)))) * np.cos(t*_wd)


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
