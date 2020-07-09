"""
Fucntions for visualization
- Line styles: for better plots
- Create a figure
- Plot nodes (oscillators, possibly with phases), edges, nodes numbers
- phase_plot: plot a function of phi1,phi2, e.g. friction function
- Stability plots
"""
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
from matplotlib.colors import SymLogNorm, Normalize, LogNorm
from cycler import cycler

default_colormap = 'viridis'
default_midpoint_colormap = 'RdBu_r'
default_cyclic_colormap = 'hsv'

### Colors, linestyles, figures
def get_colors(num, cmap=default_colormap):
    '''
    :param num: How many colors to return
    :param cmap: e.g. 'jet' 'viridis' 'RdBu_r' 'hsv'
    :return:
    '''
    import matplotlib.cm as mcm
    cm = getattr(mcm, cmap)
    return cm(np.linspace(0, 1, num))


# See usage below
_lines = ["-", "--", "-.", ":"]
_cmap = plt.get_cmap("tab10")
_colors = [_cmap(i) for i in range(10)]

def get_colorcycler():  # first cycle colors, then line styles
    cyc = cycler('linestyle', _lines) * cycler('color', _colors)
    return cyc()


def get_linecycler():  # first cycle line styles, then colors
    cyc = cycler('color', _colors) * cycler('linestyle', _lines)
    return cyc()


def get_stylecycler():  # Only cycle linestyle; color will be changed by matplotlib automatically
    cyc = cycler('linestyle', _lines)
    return cyc()


## HOW TO USE STYLECYCLERS:
# (1)
# styles = vis.get_linecycler()
# for i in range(12):
#     plt.plot(np.linspace(0, 2), np.linspace(0, 2) ** 2 + i,**next(styles))
# plt.show()
# (2) OR Change globally
# rc('axes', prop_cycle = get_linecycler())

## More lines options
# lines = [(0,()), # solid
#          (0, (5, 5)), # dashed
#          (0, (3, 5, 1, 5)), # dashdotted
#          (0, (3, 5, 1, 5, 1, 5)), # dashdotdotted
#          (0, (1, 5)), # dotted
#          (0, (5, 1)), # densely dashed
#          (0, (3, 1, 1, 1)), # densely dashdotted
#          (0, (1, 1)), # densely dotted
#          (0, (5, 10)), # loosely dashed
#          (0, (3, 1, 1, 1, 1, 1))] # densely dashdotdotted?


def landscape_figure(height):
    fig = plt.figure(figsize=(1.618 * height, height))
    ax = plt.subplot(111)
    return fig, ax


### Plotting routines

def plot_nodes(coords, phi=None, color=(0.5, 0.5, 0.5), s=100,
                cmap = 'hsv',
                colorbar=True, vmin=0, vmax=2 * np.pi, zorder=2,
                fig=None, ax=None, **kwargs):
    '''
    :param color: if no phi is given, use this color to color all nodes
    :param s: point size - 1 number or a list
    :param zorder: zorder=2 to plot on top of edges
    :param kwargs: dictionary of keyword arguments for scatter
    '''
    # Get a figure and axis to draw on, unless they were already specified in input
    if ax is None:
        ax = plt.gca()
    if fig is None:
        fig = plt.gcf()

    if phi is not None:
        kwargs['c'] = phi
    else:
        kwargs['color'] = color
    # Plot nodes
    ax.set_aspect('equal')
    #  plt.gca().get_xaxis().set_visible(False)
    #  plt.gca().get_yaxis().set_visible(False)
    # plt.gcf().set_size_inches(8,6)

    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    # if isinstance(cmap, str):
    #     cmap = colors.Colormap(cmap)
    cax = ax.scatter(coords[:, 0], coords[:, 1], norm=norm, cmap=cmap, s=s, zorder=zorder,
                     **kwargs) # , markersize=24

    ## Colorbar
    if colorbar is True and phi is not None:
        cb = fig.colorbar(cax, cmap=cmap, norm=norm, fraction=0.046, pad=0.04) # Default fraction0.15, pad=0.05
        if vmin == 0 and vmax == 2 * np.pi:
            cb.set_ticks(ticks=[0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi])
            cb.set_ticklabels(['$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3 \pi}{2}$', r'$2 \pi$'])

    return fig, ax



def plot_edges(coords, T1, fig=None, ax=None, color='red', zorder=1, **kwargs):
    """
    :param zorder: 1 - to plot nodes on top
    :param kwargs: keyword arguments for plot
    """

    # Get a figure and axis to draw on, unless they were already specified in input
    if ax is None:
        ax = plt.gca()
    if fig is None:
        fig = plt.gcf()
    # Draw edges
    for (coord, translations) in zip(coords, T1):
        for translation in translations:
            points_to_plot = np.array([coord, coord + translation])
            ax.plot(points_to_plot[:, 0], points_to_plot[:, 1], '--', color=color, zorder=zorder, **kwargs)

    return fig, ax


def plot_node_numbers(coords, spacing, fig=None, ax=None, fontsize=12, shift=(-0.3,-0.1), zorder=3):
    '''
    :param shift: how much the numbers should be shifted with respect to the node center;
                  Will be shifted by shift * spacing
    :return:
    '''
    # Get a figure and axis to draw on, unless they were already specified in input
    if ax is None:
        ax = plt.gca()
    if fig is None:
        fig = plt.gcf()
    # Draw text
    for ix, coord in enumerate(coords):
        ax.text(coord[0] + shift[0] * spacing, coord[1] + shift[1] * spacing, str(ix), fontsize=fontsize, zorder=zorder)

    return fig, ax


## Plot local friction  gij(phi1,phi2)

# HELPER FUNCTIONS

class MidpointNormalize(colors.Normalize):
    '''
    Note:
         - check out colors.OffsetNorm() in newer versions
         - matplotlib 3.1 https://matplotlib.org/3.1.0/api/_as_gen/matplotlib.colors.DivergingNorm.html
    '''

    def __init__(self, vmin=None, vmax=None, midpoint=0, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        normalized_min = max(0, 1 / 2 * (1 - abs((self.midpoint - self.vmin) / (self.midpoint - self.vmax))))
        normalized_max = min(1, 1 / 2 * (1 + abs((self.vmax - self.midpoint) / (self.midpoint - self.vmin))))
        normalized_mid = 0.5
        x, y = [self.vmin, self.midpoint, self.vmax], [normalized_min, normalized_mid, normalized_max]
        return np.ma.masked_array(np.interp(value, x, y))


def _trajectory(ax, x, y, colors, *args, **kwargs):  # use 'quick.get_colors' as 'colors' argument
    _field(ax, x, y, list(range(len(x))), colors(len(x)), *args, **kwargs)


def _field(ax, x, y, z, colors, *args, **kwargs):
    mn = np.amin(z)
    mx = np.amax(z)
    for (x1, x2, y1, y2, zz) in zip(x, x[1:], y, y[1:], z):
        index = int((len(z) - 1) * (zz - mn) / (mx - mn))
        ax.plot([x1, x2], [y1, y2], *args,
                color=colors[index], **kwargs)


def _midpoint_imshow(vals, x1_min, x1_max, x2_min, x2_max, ax=None, colorbar=True,
                     title=None, xlabel=None, ylabel=None, midpoint=None, norm=None, cmap=None,
                     **kwargs):
    '''
    If =ax= exists, plot to ax, else plot to file or on screen if there is no filename given
	vals: [[f(x,y) for x in xs] for y in ys]
    '''
    if ax == None:
        ax = plt.gca()

    if cmap is None:
        if midpoint is not None:
            cmap = default_midpoint_colormap
        else:
            cmap = default_colormap

    if midpoint is not None:
        norm = MidpointNormalize(np.amin(vals), np.amax(vals), midpoint)

    mappable = ax.imshow(vals, extent=[x1_min, x1_max, x2_min, x2_max],
                         origin='lower', cmap=cmap, norm=norm,
                         **kwargs)  # origin lower - to specifiy how this matrix has to b

    ax.set_xlim((x1_min, x1_max))
    ax.set_ylim((x2_min, x2_max))
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    if title is not None:
        ax.set_title(title)
    ax.set_aspect(1)

    if colorbar:
        plt.colorbar(mappable)

    return ax


# END: HElPER FUNCTIONS

# Function to plot local friction: gij(phi1, phi2)


def phase_plot(vals, ax=None, colorbar=True,
               title=None, midpoint=None, norm=None, cmap=None,
               xlabel=r'$\varphi_1$', ylabel=r'$\varphi_2$',
               phase_colored_lines=True, num_phases=20,
               **kwargs):
    '''
    Makes a 2D colormap in the range from 0 to 2pi; add ticks proportional to pi/2;
    :param vals: 2D array
    :param phase_colored_lines: If True - add colors near the plot, representing each of phases, as defined by the default colormap
    :param num_phases: Used by phase_colored_lines;
    :param filename: If given - will save the plot to the file. If ax is given - it has to be a quick.Plot instance
    :param kwargs:
    :return:
    '''
    ax = _midpoint_imshow(vals, 0, 2 * np.pi, 0, 2 * np.pi,
                          ax=ax, colorbar=colorbar, title=title,
                          midpoint=midpoint, norm=norm, cmap=cmap,
                          xlabel=xlabel, ylabel=ylabel, **kwargs)
    # set x phase ticks
    ax.set_xticks(ticks=[0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi])
    ax.set_xticklabels(['$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3 \pi}{2}$', r'$2 \pi$'])
    # set y phase ticks
    ax.set_yticks(ticks=[0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi])
    ax.set_yticklabels(['$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3 \pi}{2}$', r'$2 \pi$'])

    if phase_colored_lines:
        color_func = lambda x: get_colors(x, cmap=default_cyclic_colormap)
        shift = 2 * np.pi / 200 * 2
        phi = np.array(range(num_phases)) * 2 * np.pi / (num_phases - 1)
        # color_coding
        magic_num = 10
        for i in range(magic_num):
            _trajectory(ax, phi, np.zeros(num_phases) - shift * (i + 1 / 2), color_func, lw=2)
            _trajectory(ax, np.zeros(num_phases) - shift * (i + 1 / 2), phi, color_func, lw=2)

        ax.set_xlim((-shift * (magic_num), np.pi * 2))
        ax.set_ylim((-shift * (magic_num), np.pi * 2))



if __name__ == '__main__':
    '''
    OK: plot nodes of a triangular lattice
    +-OK: plot node numbers
    '''
    import carpet.lattice.triangular as lattice

    a = 10
    nx = 4
    ny = 8

    coords, lattice_ids = lattice.get_nodes_and_ids(nx, ny, a)
    N1, T1 = lattice.get_neighbours_list(coords, nx, ny, a)


    plot_edges(coords, T1, color='blue')
    fig, ax = plot_nodes(coords, phi=np.linspace(0, 2 * np.pi, len(coords)), vmin=0, vmax=2 * np.pi)

    plot_node_numbers(coords, a, fig=fig, ax=ax)

    plt.savefig('1.png', bbox_inches='tight')
