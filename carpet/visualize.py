"""
Fucntions for visualization
- Line styles: for better plots
- Create a figure
- Plot nodes (oscillators, possibly with phases), edges, nodes numbers
- phase_plot: plot a function of phi1,phi2, e.g. friction function
- Stability plots
"""
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
from matplotlib.colors import SymLogNorm # , Normalize, LogNorm

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
    import matplotlib.cm as ccmm
    cm = getattr(ccmm, cmap)
    return cm(sp.linspace(0, 1, num))


_lines = ["-", "--", "-.", ":"]
_cmap = plt.get_cmap("tab10")
_colors = [_cmap(i) for i in range(10)]

def get_colorcycler():  # first cycle colors, then line styles
    from cycler import cycler
    cyc = cycler('linestyle', _lines) * cycler('color', _colors)
    return cyc()


def get_linecycler():  # first cycle line styles, then colors
    from cycler import cycler
    cyc = cycler('color', colors) * cycler('linestyle', _lines)
    return cyc()


def get_stylecycler():  # Only cycle linestyle; color will be changed by matplotlib automatically
    from cycler import cycler
    cyc = cycler('linestyle', _lines)
    return cyc()


## HOW TO USE STYLECYCLERS:
# (1)
# styles = get_linecycler()
# for i in range(12):
#     plt.plot(sp.linspace(0, 2), sp.linspace(0, 2) ** 2 + i,**next(styles))
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
                colorbar=True, vmin=0, vmax=2 * sp.pi, zorder=2,
                fig=None, ax=None):
    '''
    :param color: if no phi is given, use this color to color all nodes
    :param s: point size - 1 number or a list
    :param zorder: 2: to plot on top of edges
    '''
    # Get a figure and axis to draw on, unless they were already specified in input
    if ax is None:
        ax = plt.gca()
    if fig is None:
        fig = plt.gcf()

    if phi is not None:
        colors = phi
    else:
        colors = color
    # Plot nodes
    ax.set_aspect('equal')
    #  plt.gca().get_xaxis().set_visible(False)
    #  plt.gca().get_yaxis().set_visible(False)
    # plt.gcf().set_size_inches(8,6)

    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    # if isinstance(cmap, str):
    #     cmap = colors.Colormap(cmap)
    ax.scatter(coords[:, 0], coords[:, 1], c=colors, norm=norm, cmap=cmap, s=s, zorder=zorder) # , markersize=24

    if colorbar is True and phi is not None:
        legend = fig.add_axes([0.88, 0.25, 0.07, 0.5])  # [0.85, 0.25, 0.1, 0.5]
        cb = mpl.colorbar.ColorbarBase(legend, cmap=cmap, norm=norm, orientation='vertical')
        if vmin == 0 and vmax == 2 * sp.pi:
            cb.set_ticks(ticks=[0, sp.pi / 2, sp.pi, 3 * sp.pi / 2, 2 * sp.pi])
            cb.set_ticklabels(['$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3 \pi}{2}$', r'$2 \pi$'])

    return fig, ax



def plot_edges(coords, T1, fig=None, ax=None, zorder=1):
    """
    :param zorder: 1 - to plot nodes on top
    """

    # Get a figure and axis to draw on, unless they were already specified in input
    if ax is None:
        ax = plt.gca()
    if fig is None:
        fig = plt.gcf()
    # Draw edges
    for (coord, translations) in zip(coords, T1):
        for translation in translations:
            points_to_plot = sp.array([coord, coord + translation])
            ax.plot(points_to_plot[:, 0], points_to_plot[:, 1], 'r--', zorder=zorder)

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
    Note: check out colors.OffsetNorm() in newer versions
    '''

    def __init__(self, vmin=None, vmax=None, midpoint=0, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        normalized_min = max(0, 1 / 2 * (1 - abs((self.midpoint - self.vmin) / (self.midpoint - self.vmax))))
        normalized_max = min(1, 1 / 2 * (1 + abs((self.vmax - self.midpoint) / (self.midpoint - self.vmin))))
        normalized_mid = 0.5
        x, y = [self.vmin, self.midpoint, self.vmax], [normalized_min, normalized_mid, normalized_max]
        return sp.ma.masked_array(sp.interp(value, x, y))


def _trajectory(ax, x, y, colors, *args, **kwargs):  # use 'quick.get_colors' as 'colors' argument
    _field(ax, x, y, list(range(len(x))), colors(len(x)), *args, **kwargs)


def _field(ax, x, y, z, colors, *args, **kwargs):
    mn = sp.amin(z)
    mx = sp.amax(z)
    for (x1, x2, y1, y2, zz) in zip(x, x[1:], y, y[1:], z):
        index = int((len(z) - 1) * (zz - mn) / (mx - mn))
        ax.plot([x1, x2], [y1, y2], *args,
                color=colors[index], **kwargs)


def _midpoint_imshow(vals, x1_min, x1_max, x2_min, x2_max, ax=None, colorbar=True,
                     title=None, xlabel=None, ylabel=None, midpoint=None, norm=None, cmap=None,
                     **kwargs):
    '''
    If =ax= exists, plot to ax, else plot to file or on screen if there is no filename given
    '''
    if ax == None:
        ax = plt.gca()

    if cmap is None:
        if midpoint is not None:
            cmap = default_midpoint_colormap
        else:
            cmap = default_colormap

    if midpoint is not None:
        norm = MidpointNormalize(sp.amin(vals), sp.amax(vals), midpoint)

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
               xlabel=r'$\phi_1$', ylabel=r'$\phi_2$',
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
    ax = _midpoint_imshow(vals, 0, 2 * sp.pi, 0, 2 * sp.pi,
                          ax=ax, colorbar=colorbar, title=title,
                          midpoint=midpoint, norm=norm, cmap=cmap,
                          xlabel=xlabel, ylabel=ylabel, **kwargs)
    # set x phase ticks
    ax.set_xticks(ticks=[0, sp.pi / 2, sp.pi, 3 * sp.pi / 2, 2 * sp.pi])
    ax.set_xticklabels(['$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3 \pi}{2}$', r'$2 \pi$'])
    # set y phase ticks
    ax.set_yticks(ticks=[0, sp.pi / 2, sp.pi, 3 * sp.pi / 2, 2 * sp.pi])
    ax.set_yticklabels(['$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3 \pi}{2}$', r'$2 \pi$'])

    if phase_colored_lines:
        color_func = lambda x: get_colors(x, cmap=default_cyclic_colormap)
        shift = 2 * sp.pi / 200 * 2
        phi = sp.array(range(num_phases)) * 2 * sp.pi / (num_phases - 1)
        # color_coding
        magic_num = 10
        for i in range(magic_num):
            _trajectory(ax, phi, sp.zeros(num_phases) - shift * (i + 1 / 2), color_func, lw=2)
            _trajectory(ax, sp.zeros(num_phases) - shift * (i + 1 / 2), phi, color_func, lw=2)

        ax.set_xlim((-shift * (magic_num), sp.pi * 2))
        ax.set_ylim((-shift * (magic_num), sp.pi * 2))



## Visualize stability of fixpoints
## - For each fixpoint, plot in dual lattice space an eigenvalue with the biggest real part (exclude neutral perturbation)


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    import matplotlib.colors as colors
    if type(cmap) is str:
        cmap = plt.get_cmap(cmap)
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(sp.linspace(minval, maxval, n)))
    return new_cmap


def _collect_values(ev_dict, get_k, eps=10 ** -8):
    ks_stable = []  # vectors k
    vals_stable = []  # colors

    ks_unstable = []
    vals_unstable = []

    ks_saddle = []
    for m1 in range(nx):
        for m2 in range(ny):
            k = get_k(m1)
            if abs(ev_dict[m1, m2][-1]) < eps:
                ev = sp.real(ev_dict[m1, m2][-2])
            else:
                ev = sp.real(ev_dict[m1, m2][-1])

            # Stable
            if ev < 0:
                ks_stable.append(k)
                vals_stable.append(ev)

            # Unstable
            if ev > 0:
                ks_unstable.append(k)
                vals_unstable.append(ev)

                # Saddle nodes
            if ev > 0 and sp.real(ev_dict[m1, m2][0]) < 0:
                ks_saddle.append(k)

    ks_stable = sp.array(ks_stable)
    ks_unstable = sp.array(ks_unstable)
    ks_saddle = sp.array(ks_saddle)
    return (vals_stable, ks_stable), (vals_unstable, ks_unstable), ks_saddle


def plot_stability(ev_dict, get_k, range_stable=None, range_unstable=None):
    '''
    Must be followed by plt.show()

    - Eigenvalues to m-twist
    2019-08-05: include saddle nodes
    '''
    (vals_stable, ks_stable), (vals_unstable, ks_unstable), ks_saddle = _collect_values(ev_dict, get_k)

    # Plot Stable
    cmap = 'Greens_r'
    if range_stable is not None:
        cmap_norm = SymLogNorm(vmin=range_stable[0], vmax=range_stable[1], linthresh=10 ** -5)
    else:
        cmap_norm = SymLogNorm(vmin=sp.amin(vals_stable), vmax=sp.amax(vals_stable), linthresh=10 ** -5)
    plt.scatter(ks_stable[:, 0], ks_stable[:, 1], c=vals_stable, cmap=cmap, norm=cmap_norm)
    plt.colorbar()

    # Plot Unstable
    cmap = 'Reds'  # truncate_colormap('Reds', 0.6, 0.95)
    if range_stable is not None:
        cmap_norm = SymLogNorm(vmin=range_unstable[0], vmax=range_unstable[1], linthresh=10 ** -5)
    else:
        cmap_norm = SymLogNorm(vmin=sp.amin(vals_unstable), vmax=sp.amax(vals_unstable), linthresh=10 ** -5)
    plt.scatter(ks_unstable[:, 0], ks_unstable[:, 1], c=vals_unstable, cmap=cmap, norm=cmap_norm)
    plt.colorbar()

    # Highlight saddle nodes
    plt.scatter(ks_saddle[:, 0], ks_saddle[:, 1], facecolors='none', edgecolor='green', linewidths=2)

def plot_stability2(ev_dict, get_k, range_stable=None, range_unstable=None):
    '''
    Must be followed by plt.show()

    - Eigenvalues to m-twist
    2019-08-05: include saddle nodes
    '''
    (vals_stable, ks_stable), (vals_unstable, ks_unstable), ks_saddle = _collect_values(ev_dict, get_k)

    # Plot Stable
    cmap = 'Greens_r'
    if range_stable is not None:
        cmap_norm = SymLogNorm(vmin=range_stable[0], vmax=range_stable[1], linthresh=10 ** -5)
    else:
        cmap_norm = SymLogNorm(vmin=sp.amin(vals_stable), vmax=sp.amax(vals_stable), linthresh=10 ** -5)
    plt.scatter(ks_stable[:, 0], ks_stable[:, 1], c=vals_stable, cmap=cmap, norm=cmap_norm)
    plt.colorbar()

    # Plot Unstable
    cmap = 'Reds'  # truncate_colormap('Reds', 0.6, 0.95)
    if range_stable is not None:
        cmap_norm = SymLogNorm(vmin=range_unstable[0], vmax=range_unstable[1], linthresh=10 ** -5)
    else:
        cmap_norm = SymLogNorm(vmin=sp.amin(vals_unstable), vmax=sp.amax(vals_unstable), linthresh=10 ** -5)
    plt.scatter(ks_unstable[:, 0], ks_unstable[:, 1], c=vals_unstable, cmap=cmap, norm=cmap_norm)
    plt.colorbar()

    # Highlight saddle nodes
    plt.scatter(ks_saddle[:, 0], ks_saddle[:, 1], facecolors='none', edgecolor='green', linewidths=2)



## Plot eigenvalues vs k
## - evec2mtwist - in lattice
## - print decomposition

def plot_evals_k_vec(evals, evec2mtwist, get_k, range_stable=None, range_unstable=None):
    '''
    Must be followed by plt.show()
    Given a mapping from eigenvectors to mtwists, plot REAL parts of eigenvalues in dual lattice space.
    '''
    from matplotlib.colors import SymLogNorm
    ks_stable = []  # vectors k
    vals_stable = []  # colors

    ks_unstable = []
    vals_unstable = []

    for ievec, (m1, m2) in evec2mtwist.items():
        ev = evals[ievec].real  # Only consider real part
        (m1, m2) = evec2mtwist[ievec]
        k = get_k(m1, m2)
        # Stable
        if ev < 0:
            ks_stable.append(k)
            vals_stable.append(ev)

        # Unstable
        if ev > 0:
            ks_unstable.append(k)
            vals_unstable.append(ev)

    ks_stable = sp.array(ks_stable)
    ks_unstable = sp.array(ks_unstable)

    # Plot Stable
    if len(vals_stable) > 0:
        cmap = 'Greens_r'
        if range_stable is not None:
            cmap_norm = SymLogNorm(vmin=range_stable[0], vmax=range_stable[1], linthresh=10 ** -5)
        else:
            cmap_norm = SymLogNorm(vmin=sp.amin(vals_stable), vmax=sp.amax(vals_stable), linthresh=10 ** -5)
        plt.scatter(ks_stable[:, 0], ks_stable[:, 1], c=vals_stable, cmap=cmap, norm=cmap_norm)
        plt.colorbar()

    if len(vals_unstable) > 0:
        # Plot Unstable
        cmap = 'Reds'  # truncate_colormap('Reds', 0.6, 0.95)
        # norm = Normalize(vmin=sp.amin(vals_unstable), vmax=sp.amax(vals_unstable))
        if range_stable is not None:
            cmap_norm = SymLogNorm(vmin=range_unstable[0], vmax=range_unstable[1], linthresh=10 ** -5)
        else:
            cmap_norm = SymLogNorm(vmin=sp.amin(vals_unstable), vmax=sp.amax(vals_unstable), linthresh=10 ** -5)
        plt.scatter(ks_unstable[:, 0], ks_unstable[:, 1], c=vals_unstable, cmap=cmap, norm=cmap_norm)
        plt.colorbar()


if __name__ == '__main__':
    '''
    OK: plot nodes of a triangular lattice
    +-OK: plot node numbers
    '''
    import carpet.lattice_triangular as lattice

    a = 10
    nx = 5
    ny = 4

    coords, lattice_ids = lattice.get_nodes_and_ids(nx, ny, a)
    N1, T1 = lattice.get_neighbours_list(coords, nx, ny, a)


    plot_edges(coords, T1)
    fig, ax = plot_nodes(coords, phi=sp.linspace(0, 2 * sp.pi, len(coords)), vmin=0, vmax=2 * sp.pi)

    plot_node_numbers(coords, a, fig=fig, ax=ax)

    plt.savefig('1.png', bbox_inches='tight')
