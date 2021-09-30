import matplotlib.pyplot as plt
import utils.expectation_values as expect
import qutip
import numpy as np


def graph_solution(solution, system):
    expectation_values = expect.get_all_expectation_values(solution.states, system.basis_states)
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.set_xlabel("duration (ns)")
    ax.set_ylabel("expectation value")
    for i in range(0, len(expectation_values)):
        ax.plot(solution.times, expectation_values[i], label="$\psi_{%i}$" % i)
    plt.legend()

def plot_on_bloch_sphere(file_name, solution, skip=10):
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

