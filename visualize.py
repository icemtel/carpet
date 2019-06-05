import scipy as sp
import matplotlib.pyplot as plt
import colorsys
import matplotlib as mpl
import matplotlib.colors as colors

def get_colors(num, cmap='viridis'):
    '''
    :param num: How many colors to return
    :param cmap: e.g. 'jet' 'viridis' 'RdBu_r' 'hsv'
    :return:
    '''
    import matplotlib.cm as ccmm
    cm = getattr(ccmm, cmap)
    return cm(sp.linspace(0, 1, num))


def plot_nodes(coords, Phi=None, color=(0.5, 0.5, 0.5), colorbar=True):
    '''
    TODO: use scatterplot instead of a loop
    '''
    for ix, coord in enumerate(coords):
        if Phi is None:
            rgb = color
        else:
            rgb = colorsys.hsv_to_rgb(Phi[ix] / (2 * sp.pi), 1, 1)
        plt.plot(coord[0], coord[1], '.', color=rgb, markersize=24)
    plt.gca().set_aspect('equal')
    #  plt.gca().get_xaxis().set_visible(False)
    #  plt.gca().get_yaxis().set_visible(False)
    # plt.gcf().set_size_inches(6, 8)

    # Add colorbar
    if Phi is not None and colorbar is True:
        norm = mpl.colors.Normalize(vmin=0, vmax=2 * sp.pi)
        cmap = 'hsv'
        # if isinstance(cmap, str):
        #     cmap = colors.Colormap(cmap)

        legend = plt.gcf().add_axes([0.92, 0.25, 0.07, 0.5])  # [0.85, 0.25, 0.1, 0.5]
        cb = mpl.colorbar.ColorbarBase(legend, cmap=cmap, norm=norm, orientation='vertical')
        cb.set_ticks(ticks=[0, sp.pi / 2, sp.pi, 3 * sp.pi / 2, 2 * sp.pi])
        cb.set_ticklabels(['$0$', r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3 \pi}{2}$', r'$2 \pi$'])


def plot_edges(coords, T1):
    for (coord, translations) in zip(coords, T1):
        for translation in translations:
            points_to_plot = sp.array([coord, coord + translation])
            plt.plot(points_to_plot[:, 0], points_to_plot[:, 1], 'r--')


def plot_node_numbers(coords, spacing):
    for ix, coord in enumerate(coords):
        plt.text(coord[0] - 0.125 * spacing, coord[1] - 0.1 * spacing, str(ix), fontsize=12)


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
    plot_nodes(coords, Phi=sp.linspace(0, 2 * sp.pi, len(coords)))

    plt.savefig('1.png')
