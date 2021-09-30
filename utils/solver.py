import qutip.ui.progressbar as progress_bar
from qutip import *

_atol = 1e-10
_rtol = 1e-10
_norm_tol = 1e-10

def time_evolve(
    H_d,
    H_c, # only single control is allowed
    pertubation_func,
    tlist,
    starting_state,
    args={},
    pbar=False,
    store_states=True
):
    states_to_solve_for = []

    options = Options(
        store_states=store_states, store_final_state=True, atol=_atol, rtol=_rtol, norm_tol=_norm_tol
    )  # integrator will store results at each point

    # setup progress bar
    pbar = progress_bar.TextProgressBar()
    if pbar is not True:
        pbar = None

    return mesolve( # actually just reverts to shrodinger solver under the hood for our closed system
        [H_d, [H_c, pertubation_func]],
        starting_state,
        tlist,
        [],
        states_to_solve_for,
        args=args,
        options=options,
        progress_bar=pbar,
    )

def set_tolerances(atol=1e-10, rtol=1e-10, norm_tol=1e-10):
    global _atol, _rtol, _norm_tol
    _atol = atol
    _rtol = rtol
    _norm_tol = norm_tol
