import qutip.ui.progressbar as progress_bar
from qutip import *

def time_evolve(
    H_d,
    H_c, # only single control is allowed
    pertubation_func,
    tlist,
    starting_state,
    args={},
    pbar=False
):
    states_to_solve_for = []

    options = Options(
        store_states=True, store_final_state=True,
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
