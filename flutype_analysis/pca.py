"""
PCA analysis of data.
"""
from __future__ import print_function, absolute_import
from matplotlib import pyplot as plt
from flutype_analysis import utils


def pca_2dplot(output_pca, **kwargs):
    fig = plt.figure(**kwargs)
    plt.xlabel("princple component 1", size='large')
    plt.ylabel("princple component 2", size='large')
    plt.tick_params(labelsize="large")
    plt.scatter(output_pca[0], output_pca[1],
                c=output_pca["color"].values, cmap=utils.discrete_cmap(output_pca["color"].values.max(), 'jet'))

    color_mapping = output_pca[['color']].sort_values(by=["color"])

    cbar = plt.colorbar(ticks=sorted(output_pca["color"].unique()))
    cbar.set_ticklabels(color_mapping.index.droplevel(1).unique())
    cbar.ax.tick_params(labelsize="large")
    plt.clim(0.5, output_pca["color"].values.max() + 0.5)

    return fig
