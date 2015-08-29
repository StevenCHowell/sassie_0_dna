import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np


def zoomout(ax, factor):

    xlim = ax.get_xlim()
    if 'log' == ax.get_yscale():
        xlim = np.log(xlim)
        xlim = (xlim[0] + xlim[1]) / 2 + np.array((-0.5, 0.5)) * \
            (xlim[1] - xlim[0]) * (1 + factor)
        xlim = np.exp(xlim)
    else:
        xlim = (xlim[0] + xlim[1]) / 2 + np.array((-0.5, 0.5)) * \
            (xlim[1] - xlim[0]) * (1 + factor)
    ax.set_xlim(xlim)

    ylim = ax.get_ylim()
    if 'log' == ax.get_yscale():
        ylim = np.log(ylim)
        ylim = (ylim[0] + ylim[1]) / 2 + np.array((-0.5, 0.5)) * \
            (ylim[1] - ylim[0]) * (1 + factor)
        ylim = np.exp(ylim)
    else:
        ylim = (ylim[0] + ylim[1]) / 2 + np.array((-0.5, 0.5)) * \
            (ylim[1] - ylim[0]) * (1 + factor)
    ax.set_ylim(ylim)


def xyplot(data, fmt='-', label=''):
    return plt.plot(data[:, 0], data[:, 1], fmt, label=label)


def xyerror(data, fmt='-', label=''):
    return plt.errorbar(data[:, 0], data[:, 1], data[:, 2], fmt=fmt, label=label)


def logx():
    return plt.xscale('log')


def logy():
    return plt.yscale('log')


def iqlabel():
    plt.xlabel(r'$Q (\AA^{-1})$')
    plt.ylabel(r'$I(Q)$')


def prlabel():
    plt.xlabel(r'$R (\AA)$')
    plt.xlabel(r'$P(R)$')


def symbol_order(i, l=False):
    symbols = ['s', 'o', '^', 'd', 'v', '*', '<', 'p', 'h', 'x']
    if l:
        return symbols[i % len(symbols)] + l
    else:
        return symbols[i % len(symbols)]


def color_order(i):
    colors = [[0,    0,    1],
              [0,  0.5,    0],
              [1,    0,    0],
              [0.75,    0, 0.75],
              [0, 0.75, 0.75],
              [0.75, 0.75,    0],
              # [   0,    1,    0],
              [0.25,    0, 0.25],
              # [   1,    1,    0],
              [1,  0.5,    0]]
    return colors[i % len(colors)]


def qual_color(i, style='set4'):
    '''
    http://www.personal.psu.edu/cab38/ColorBrewer/ColorBrewer.html
    9-class qualitative Set1,
    '''
    # 12-class paired
    pair = [[166, 206, 227],
            [31, 120, 180],
            [178, 223, 138],
            [51, 160,  44],
            [251, 154, 153],
            [227,  26,  28],
            [253, 191, 111],
            [255, 127,   0],
            [202, 178, 214],
            [106,  61, 154],
            [255, 255, 153],
            [177,  89,  40]]

    # 9-class Set1
    set1 = [[228,  26,  28],
            [55,  126, 184],
            [77,  175,  74],
            [152,  78, 163],
            [255, 127,   0],
            [255, 255,  51],
            [166,  86,  40],
            [247, 129, 191],
            [153, 153, 153]]

    # 8-class Dark2
    dark = [[27, 158, 119],
            [217,  95,   2],
            [117, 112, 179],
            [231,  41, 138],
            [102, 166,  30],
            [230, 171,   2],
            [166, 118,  29],
            [102, 102, 102]]


    #  http://tools.medialab.sciences-po.fr/iwanthue/
    set2 = [[206, 158, 154],
            [127, 206, 78],
            [201, 80, 202],
            [210, 89, 53],
            [78, 96, 56],
            [77, 59, 89],
            [138, 206, 163],
            [204, 178, 76],
            [199, 76, 125],
            [130, 173, 197],
            [117, 59, 44],
            [135, 115, 198]]

    set3 = [[24, 158, 174],
            [249, 117, 53],
            [212, 13, 158],
            [180, 192, 18],
            [79, 3, 41],
            [100, 101, 23],
            [31, 45, 57],
            [234, 126, 247],
            [29, 98, 167],
            [146, 94, 238]]

    set4 = [[56, 97, 138],
            [217, 61, 62],
            [75, 108, 35],
            [100, 68, 117],
            [228, 182, 48],
            [183, 92, 56],
            [107, 171, 215],
            [209, 84, 175],
            [177, 191, 57],
            [126, 116, 209]
            ]

    set5 = [[194, 141, 57],
            [173, 95, 211],
            [78, 156, 139],
            [108, 173, 68],
            [97, 77, 121],
            [166, 82, 84],
            [84, 94, 43],
            [202, 82, 147],
            [205, 73, 52],
            [128, 147, 203]]

    mpl_set = plt.cm.Set3(np.linspace(0, 1, 12))[:, :3] * 255.0

    styles = {'set1': set1, 'pair': pair, 'dark': dark, 'set2': set2,
              'set3': set3, 'set4': set4, 'set5': set5, 'mpl_set': mpl_set}
    colors = styles[style]
    return np.array(colors[i % len(colors)]) / 255.0


def seq_color(i, sort=True):
    '''
    http://geog.uoregon.edu/datagraphics/color_scales.htm
    Stepped-sequential scheme, 5 hues x 5 saturation/value levels
    '''
    colors = [[0.600, 0.060, 0.060],
              [0.700, 0.175, 0.175],
              [0.800, 0.320, 0.320],
              [0.900, 0.495, 0.495],
              [1.000, 0.700, 0.700],
              [0.600, 0.330, 0.060],
              [0.700, 0.438, 0.175],
              [0.800, 0.560, 0.320],
              [0.900, 0.697, 0.495],
              [1.000, 0.850, 0.700],
              [0.420, 0.600, 0.060],
              [0.525, 0.700, 0.175],
              [0.640, 0.800, 0.320],
              [0.765, 0.900, 0.495],
              [0.900, 1.000, 0.700],
              [0.060, 0.420, 0.600],
              [0.175, 0.525, 0.700],
              [0.320, 0.640, 0.800],
              [0.495, 0.765, 0.900],
              [0.700, 0.900, 1.000],
              [0.150, 0.060, 0.600],
              [0.262, 0.175, 0.700],
              [0.400, 0.320, 0.800],
              [0.562, 0.495, 0.900],
              [0.750, 0.700, 1.000]]
    if sort:
        j, k = divmod(i, 5)
        i = 5 * k + j

    return colors[i % len(colors)]


def make_colormap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    Example:
        rvb = make_colormap([c('red'), c('violet'), 0.33, c('violet'),
            c('blue'), 0.66, c('blue')])
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)


def diverge_map(low=qual_color(0), high=qual_color(1)):
    c = mcolors.ColorConverter().to_rgb
    if isinstance(low, basestring):
        low = c(low)
    if isinstance(high, basestring):
        high = c(high)
    return make_colormap([low, c('white'), 0.5, c('white'), high])

    # def diverge_map(high=(0.565, 0.392, 0.173), low=(0.094, 0.310, 0.635)):
    # c = mcolors.ColorConverter().to_rgb
    # if isinstance(low, basestring): low = c(low)
    # if isinstance(high, basestring): high = c(high)
    # return make_colormap([low, c('white'), 0.5, c('white'), high])

if __name__ == '__main__':
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    j = 11
    x = np.array(range(j))
    y = np.ones(j)
    fig = plt.figure()
    for i in x:
        y[:] = i
        s = symbol_order(i, '-')
        plt.plot(x, y, s, mec=color_order(i), c=color_order(i), ms=10,
                 mfc='none', label=str(color_order(i)), linewidth=5)
        # plt.plot(x, y, s, ms = 10)
        # plt.scatter(x, y, s, markeredgecolor=color_order(i), facecolors='none')
    dy = 0.1
    plt.ylim([-dy, x[-1] + dy])
    leg = plt.legend(scatterpoints=1, numpoints=1)
    plt.show()
    name = 'python_symbol_color'
    fig.savefig('/home/schowell/Dropbox/gw_phd/%s.eps' % name)
    fig.savefig('/home/schowell/Dropbox/gw_phd/%s.png' % name)
    fig.savefig('%s.eps' % name)
    fig.savefig('%s.png' % name)

    j = 12
    x = np.array(range(j))
    y = np.ones(j)
    fig = plt.figure()
    style = 'set1'
    style = 'dark'
    style = 'set2'
    style = 'set4'
    # style = 'mpl_set'
    for i in x:
        y[:] = i
        s = symbol_order(i, '-')
        plt.plot(x, y, s, mec=qual_color(i, style), c=qual_color(i, style), ms=15,
                 mfc='none', label=str(qual_color(i, style)), linewidth=10)
        # plt.plot(x, y, s, ms = 10)
        # plt.scatter(x, y, s, markeredgecolor=color_order(i, style), facecolors='none')
    dy = 0.1
    plt.ylim([-dy, x[-1] + dy])
    leg = plt.legend(scatterpoints=1, numpoints=1)
    plt.show()
    name = 'qual_color_' + style
    # fig.savefig('/home/schowell/Dropbox/gw_phd/%s.eps' % name)
    # fig.savefig('/home/schowell/Dropbox/gw_phd/%s.png' % name)
    fig.savefig('%s.eps' % name)
    fig.savefig('%s.png' % name)

    j = 25
    x = np.array(range(j))
    y = np.ones(j)
    fig = plt.figure()
    sort = False
    for i in x:
        y[:] = i
        s = symbol_order(i, '-')
        plt.plot(x, y, s, mec=seq_color(i, sort), c=seq_color(i, sort), ms=10,
                 mfc='none', label=str(seq_color(i, sort)), linewidth=2)
        # plt.plot(x, y, s, ms = 10)
        # plt.scatter(x, y, s, markeredgecolor=color_order(i), facecolors='none')
    dy = 0.1
    plt.ylim([-dy, x[-1] + dy])
    leg = plt.legend(scatterpoints=1, numpoints=1)
    plt.show()
    if sort:
        name = 'sorted_seq_color'
    else:
        name = 'seq_color'
    # fig.savefig('/home/schowell/Dropbox/gw_phd/%s.eps' % name)
    # fig.savefig('/home/schowell/Dropbox/gw_phd/%s.png' % name)
    fig.savefig('%s.eps' % name)
    fig.savefig('%s.png' % name)

    print '\m/ >.< \m/'
