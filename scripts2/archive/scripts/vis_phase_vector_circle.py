import numpy as np
import matplotlib.pyplot as plt
import carpet
import carpet.visualize as vis
from carpet.various import  cexp

# def cexp(vec):
#     '''
#     Calculate complex exponent
#     '''
#     return np.exp(1j * vec)

def circle_canvas(facecolor='None', edgecolor='black', linestyle='-', linewidth=2, zorder=0, fig=None, ax=None,
                  **kwargs):
    if ax is None:
        ax = plt.gca()
    if fig is None:
        fig = plt.gcf()

    kwargs_circle = dict(facecolor=facecolor, edgecolor=edgecolor, linestyle=linestyle, linewidth=linewidth,
                         zorder=zorder)
    kwargs_circle.update(kwargs)

    circle = plt.Circle((0, 0), 1, **kwargs_circle)
    ax.add_artist(circle)
    ax.set_aspect('equal')
    ax.set_xlim([-1.05, 1.05])
    ax.set_ylim([-1.05, 1.05])
    ax.axis('off')
    return fig, ax


def circle_phase_scatter(phi, s=100, fig=None, ax=None, **kwargs):
    if ax is None:
        ax = plt.gca()
    if fig is None:
        fig = plt.gcf()

    ce = cexp(phi)

    kwargs_scatter = dict(s=s)
    kwargs_scatter.update(kwargs)
    ax.scatter(ce.real, ce.imag, **kwargs_scatter)
    return fig, ax


def circle_phase_arrows(phi, color='blue', length_includes_head=True, head_width=0.05,
                        fig=None, ax=None, **kwargs):
    if ax is None:
        ax = plt.gca()
    if fig is None:
        fig = plt.gcf()

    if type(color) is str or type(color) is tuple:
        colors = [color for _ in phi]
    else:
        colors = color

    kwargs_arrow = dict(length_includes_head=length_includes_head,head_width=head_width)
    kwargs_arrow.update(kwargs)
    for i, p in enumerate(phi):
        ce = cexp(p)
        plt.arrow(0, 0, ce.real, ce.imag, color=colors[i],  **kwargs_arrow)

    return fig, ax


N = 6
phi = carpet.get_phase_random(N)

fig, ax = plt.subplots(1,1)

circle_phase_scatter(phi, ax=ax,fig=fig)

circle_canvas(ax=ax)
plt.show()