import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sps

_pulse = []
_tlist = []
_bin_size = 0
_internal_method_function = _nearest

def pulse_func(t, args=None):
    return _internal_method_function(t, _tlist, _pulse, _bin_size)

def setup(amplitude, drive_frequency, ramp_coef, tlist, bin_subdivisions=10, method="nearest"):
    global _pulse, _tlist, _bin_size, _internal_method_function
    _tlist = np.linspace(np.min(tlist), np.max(tlist), len(tlist) * bin_subdivisions)
    _bin_size = _tlist[1] - _tlist[0]
    _pulse = sps.windows.tukey(len(_tlist), ramp_coef) * np.cos(drive_frequency * _tlist) * amplitude

    if method == "nearest":
        _internal_method_function = _nearest
    elif method == "interpolate":
        _internal_method_function = _interpolate
    else:
        print("error! invalid method")

def get_pulse(tlist):
    return np.vectorize(pulse_func)(tlist)

#faster
def _nearest(t, ts, vs, bin):
    return vs[(np.abs(np.array(ts) - t)).argmin()]

#slower
def _interpolate(t, ts, vs, bin):
    t_n_index = (np.abs(np.array(ts) - t)).argmin()
    t_n = ts[t_n_index]
    if (t_n_index == 0 or t_n_index == len(ts)-1):
        return vs[t_n_index]
    if t < t_n:
        t_plus_one = t_n
        v_plus_one = vs[t_n_index]
        t_minus_one = ts[t_n_index - 1]
        v_minus_one = vs[t_n_index - 1]
    else:
        t_plus_one = ts[t_n_index + 1]
        v_plus_one = vs[t_n_index + 1]
        t_minus_one = t_n
        v_minus_one = vs[t_n_index]
    d_plus_one = (t_plus_one - t) / bin
    d_minus_one = (t - t_minus_one) / bin
    return (d_plus_one*v_plus_one + d_minus_one*v_minus_one)

    