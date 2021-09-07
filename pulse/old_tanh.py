import numpy as np
import matplotlib.pyplot as plt

_args = {}

# todo this is horribly broken!
def pulse_func(t, args=None):
    return np.cos(_args["wd"] * t) * (
            _args["g"]
            * (
                (np.tanh((t - _args["s"]) * _args["cs"]) + 1) / 2
                - (np.tanh((t - _args["d"]) * _args["cd"]) + 1) / 2
            )
        ) * _args["n"] * _args["v"]

def setup(amplitude, drive_frequency, ramp_coef, ramp_up_midpoint, ramp_down_midpoint):
    global _args
    _args = {
        "wd": drive_frequency,
        "s" : ramp_coef,
        "d" : ramp_coef,
        "cs": ramp_up_midpoint,
        "cd": ramp_down_midpoint,
        "g" : 1,
        "v" : amplitude
    }

    # todo: setup proper normalization
    _args["n"] = 1

def plot(tlist):
    pulse_values = []
    for t in tlist:
        pulse_values.append(pulse_func(t))
    plt.plot(tlist, pulse_values)