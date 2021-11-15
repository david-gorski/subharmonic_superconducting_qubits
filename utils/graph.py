import matplotlib.pyplot as plt
import utils.expectation_values as expect
import utils.solver as solver
import qutip
import numpy as np
from tqdm import tqdm as tqdm


def graph_solution(solution, system):
    expectation_values = expect.get_all_expectation_values(solution.states, system.basis_states)
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.set_xlabel("duration (ns)")
    ax.set_ylabel("expectation value")
    for i in range(0, len(expectation_values)):
        ax.plot(solution.times, expectation_values[i], label="$\psi_{%i}$" % i)
    plt.legend()

def plot_on_bloch_sphere_from_solution(file_name, solution, skip=10):
    computational_subspace_kets = []
    for i in range(0, len(solution.states), skip):
        computational_subspace_ket = [[solution.states[i][0][0][0]], [solution.states[i][1][0][0]]]
        computational_subspace_kets.append(qutip.Qobj(computational_subspace_ket))

    xs = []
    ys = []
    zs = []
    vecs = []
    for ket in computational_subspace_kets:
        xs.append(qutip.expect(qutip.ket2dm(ket), qutip.sigmax()))
        ys.append(qutip.expect(qutip.ket2dm(ket), qutip.sigmay()))
        zs.append(qutip.expect(qutip.ket2dm(ket), qutip.sigmaz()))
        vecs.append([xs[-1], ys[-1], zs[-1]])

    from matplotlib import pyplot, animation
    from mpl_toolkits.mplot3d import Axes3D

    fig = pyplot.figure()
    ax = Axes3D(fig, azim=-40, elev=30)
    sphere = qutip.Bloch(axes=ax)

    def animate(i):
        sphere.clear()
        sphere.add_points([xs[:i+1], ys[:i+1], zs[:i+1]])
        sphere.make_sphere()
        return ax

    def init():
        sphere.vector_color = ['r']
        return ax

    ani = animation.FuncAnimation(fig, animate, np.arange(len(xs)),
                                init_func=init, blit=False, repeat=False)

    ani.save('%s_bloch_sphere.mp4' % file_name, fps=20)



def plot_on_bloch_sphere(file_name, kets, skip=1):
    computational_subspace_kets = []
    for i in range(0, len(kets), skip):
        ket = kets[i]
        computational_subspace_ket = [[ket[0][0][0]], [ket[1][0][0]]]
        computational_subspace_kets.append(qutip.Qobj(computational_subspace_ket))

    xs = []
    ys = []
    zs = []
    vecs = []
    for ket in computational_subspace_kets:
        xs.append(qutip.expect(qutip.ket2dm(ket), qutip.sigmax()))
        ys.append(qutip.expect(qutip.ket2dm(ket), qutip.sigmay()))
        zs.append(qutip.expect(qutip.ket2dm(ket), qutip.sigmaz()))
        vecs.append([xs[-1], ys[-1], zs[-1]])

    from matplotlib import pyplot, animation
    from mpl_toolkits.mplot3d import Axes3D

    fig = pyplot.figure()
    ax = Axes3D(fig, azim=-40, elev=30)
    sphere = qutip.Bloch(axes=ax)

    def animate(i):
        sphere.clear()
        sphere.add_points([xs[:i+1], ys[:i+1], zs[:i+1]])
        sphere.make_sphere()
        return ax

    def init():
        sphere.vector_color = ['r']
        return ax

    ani = animation.FuncAnimation(fig, animate, np.arange(len(xs)),
                                init_func=init, blit=False, repeat=False)

    ani.save('%s_bloch_sphere.mp4' % file_name, fps=20)


def generate_fidelity_landscape_tukey(system, pulse, freq_list, dur_list, t_multiplier=5):
    fids = [] # an array of columns in the landscape
    for freq in tqdm(freq_list):
        fid_column = []
        for dur in dur_list:
            tlist = np.linspace(0, dur, int(dur*t_multiplier))
            pulse.setup(amplitude=pulse._amp, drive_frequency=freq, ramp_coef=pulse._ramp_coef, tlist=tlist)
            s = solver.time_evolve(system.H0, system.H1, pulse.pulse_func, tlist, system.starting_state)
            fidelity = expect.expectation_value(s.states[-1], system.target_state)
            fid_column.append(fidelity)
        fids.append(fid_column)

    plt.pcolormesh(freq_list / (2*np.pi), dur_list, np.swapaxes(fids,0,1), shading='auto')
    
    return fids