import scipy as sp
import matplotlib.pyplot as plt
import colorsys
import matplotlib as mpl
import matplotlib.colors as colors

default_colormap = 'viridis'
default_midpoint_colormap = 'RdBu_r'
default_cyclic_colormap = 'hsv'


def get_colors(num, cmap=default_colormap):
    '''
    :param num: How many colors to return
    :param cmap: e.g. 'jet' 'viridis' 'RdBu_r' 'hsv'
    :return:
    '''
    import matplotlib.cm as ccmm
    cm = getattr(ccmm, cmap)
    return cm(sp.linspace(0, 1, num))


def simple_figure():
    fig = plt.figure(figsize=(1.618 * 3, 3))
    ax = plt.subplot(111)
    return fig, ax


def plot_nodes(coords, phi=None, color=(0.5, 0.5, 0.5), colorbar=True, fig=None, ax=None):
    '''
    TODO: use scatterplot instead of a loop
    '''
    # Get a figure and axis to draw on, unless they were already specified in input
    if ax is None:
        ax = plt.gca()
    if fig is None:
        fig = plt.gcf()

    # Plot nodes
    for ix, coord in enumerate(coords):
        if phi is None:
            rgb = color
        else:
            rgb = colorsys.hsv_to_rgb(phi[ix] / (2 * sp.pi), 1, 1)
        ax.plot(coord[0], coord[1], '.', color=rgb, markersize=24)
    ax.set_aspect('equal')
    #  plt.gca().get_xaxis().set_visible(False)
    #  plt.gca().get_yaxis().set_visible(False)
    # plt.gcf().set_size_inches(8,6)

    # Add colorbar
    if phi is not None and colorbar is True:
        norm = mpl.colors.Normalize(vmin=0, vmax=2 * sp.pi)
        cmap = 'hsv'
        # if isinstance(cmap, str):
        #     cmap = colors.Colormap(cmap)

        legend = fig.add_axes([0.88, 0.25, 0.07, 0.5])  # [0.85, 0.25, 0.1, 0.5]
        cb = mpl.colorbar.ColorbarBase(legend, cmap=cmap, norm=norm, orientation='vertical')
        cb.set_ticks(ticks=[0, sp.pi / 2, sp.pi, 3 * sp.pi / 2, 2 * sp.pi])
        cb.set_ticklabels(['$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3 \pi}{2}$', r'$2 \pi$'])

    return fig, ax


def plot_edges(coords, T1, fig=None, ax=None):
    # Get a figure and axis to draw on, unless they were already specified in input
    if ax is None:
        ax = plt.gca()
    if fig is None:
        fig = plt.gcf()
    # Draw edges
    for (coord, translations) in zip(coords, T1):
        for translation in translations:
            points_to_plot = sp.array([coord, coord + translation])
            ax.plot(points_to_plot[:, 0], points_to_plot[:, 1], 'r--')

    return fig, ax


def plot_node_numbers(coords, spacing, fig=None, ax=None):
    # Get a figure and axis to draw on, unless they were already specified in input
    if ax is None:
        ax = plt.gca()
    if fig is None:
        fig = plt.gcf()
    # Draw text
    for ix, coord in enumerate(coords):
        ax.text(coord[0] - 0.125 * spacing, coord[1] - 0.1 * spacing, str(ix), fontsize=12)

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


if __name__ == '__main__':
    '''
    OK: plot nodes of a triangular lattice
    +-OK: plot nod numbers
    '''
    import carpet2D

    e1 = sp.array([1, 0])
    #    e2 = sp.array([0,1])
    e2 = sp.array([1 / 2, 3 ** (1 / 2) / 2])
    a = 10
    l1 = 5
    l2 = 5

    coords, lattice_ids = carpet2D.get_nodes_and_ids(e1, e2, l1, l2, a)
    N1, T1 = carpet2D.get_neighbours_list(coords, e1, e2, l1, l2, a)

    #   plot_node_numbers(coords, a)
    plot_edges(coords, T1)
    fig, ax = plot_nodes(coords, phi=sp.linspace(0, 2 * sp.pi, len(coords)))
    plot_node_numbers(coords, a, fig=fig, ax=ax)

    plt.savefig('1.png', bbox_inches='tight')
