"""
Utility and helper functions.
"""
from __future__ import print_function, absolute_import
from matplotlib import pyplot as plt


def map_strings_to_number(strings):
    """Transforms list of strings into numbers."""
    counter = 1
    d = {}
    numbers = []
    for s in strings:
        if s in d:
            pass
        else:
            d[s] = counter
            counter += 1
        numbers.append(d[s])
    return numbers


def discrete_cmap(N, base_cmap=None):
    """ Create an N-bin discrete colormap from the specified input map"""
    # Note that if base_cmap is a String or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for String, None, or a colormap instance:
    if N == 1:
        base = plt.cm.get_cmap(base_cmap)
        color_list = base(np.linspace(0, 1, N + 1))
        cmap_name = base.name + str(N)

    else:

        base = plt.cm.get_cmap(base_cmap)
        color_list = base(np.linspace(0, 1, N))
        cmap_name = base.name + str(N)

    return base.from_list(cmap_name, color_list, N)
