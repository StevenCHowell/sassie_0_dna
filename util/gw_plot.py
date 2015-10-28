import matplotlib.pyplot as plt
# import matplotlib.colors as mcolors
import numpy as np


def zoomout(ax, factor, x=True, y=True):

    if x:
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

    if y:
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


def guinier(iq_data):
    '''
    return the Guinier data
    '''
    giq_data = np.zeros(iq_data.shape)
    giq_data[:, 0] = iq_data[:, 0] * iq_data[:, 0]

    # replace non-positive number by nan
    iq = np.copy(iq_data[:, 1])
    i_nonpositive = iq <= 0
    iq[i_nonpositive] = np.nan

    giq_data[:, 1] = np.log(iq)

    return giq_data


def guinier_plot(iq_data, fmt='-', label=''):
    '''
    plot the Guinier data
    '''
    giq_data = guinier(iq_data)
    return xyplot(giq_data, fmt=fmt, label=label)


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
    symbols = ['s', 'o', '^', 'd', 'v', '*', '<', 'p', 'h']
    if l:
        return symbols[i % len(symbols)] + l
    else:
        return symbols[i % len(symbols)]


def color_order(i):
    # Xiangyun's preferred colors
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
    # solarized
    solarized = [[ 38, 139, 210], # blue      #268bd2  4/4 blue      33 #0087ff
                 [220,  50,  47], # red       #dc322f  1/1 red      160 #d70000
                 [133, 153,   0], # green     #859900  2/2 green     64 #5f8700
                 [211,  54, 130], # magenta   #d33682  5/5 magenta  125 #af005f
                 [181, 137,   0], # yellow    #b58900  3/3 yellow   136 #af8700
                 [ 42, 161, 152], # cyan      #2aa198  6/6 cyan      37 #00afaf
                 [108, 113, 196], # violet    #6c71c4 13/5 brmagenta 61 #5f5faf
                 [203,  75,  22], # orange    #cb4b16  9/3 brred    166 #d75f00
                 [131, 148, 150], # base0     #839496 12/6 brblue   244 #808080
                 [  0,  43,  54]] # base03    #002b36  8/4 brblack  234 #1c1c1c

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

    styles = {'set1': set1, 'pair': pair, 'dark': dark, 'set2': set2, 'solarized': solarized,
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
    # matplotlib.use('Agg')
    import matplotlib.pyplot as plt


    title = 'Default Matplotlib Settings (from matplotlibrc file)'
    j = 11
    x = np.array(range(j))
    y = np.ones(j)
    fig = plt.figure()
    for i in x:
        y[:] = i
        s = symbol_order(i, '-')
        plt.plot(x, y, s, ms=10, label='default color: %d' %i, linewidth=5)
    dy = 0.1
    plt.ylim([-dy, x[-1] + dy])
    leg = plt.legend(scatterpoints=1, numpoints=1)
    name = 'default_colors'
    plt.title(title)
    # fig.savefig('/home/schowell/Dropbox/gw_phd/%s.eps' % name)
    # fig.savefig('/home/schowell/Dropbox/gw_phd/%s.png' % name)
    fig.savefig('%s.eps' % name)
    fig.savefig('%s.png' % name)
    plt.show()

    title = 'Duplicating MATLAB Color Settings'
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
    name = 'python_symbol_color'
    plt.title(title)
    # fig.savefig('/home/schowell/Dropbox/gw_phd/%s.eps' % name)
    # fig.savefig('/home/schowell/Dropbox/gw_phd/%s.png' % name)
    fig.savefig('%s.eps' % name)
    fig.savefig('%s.png' % name)
    plt.show()

    title = 'Qualitative Colors'
    j = 12
    x = np.array(range(j))
    y = np.ones(j)
    fig = plt.figure()
    style = 'set1'
    style = 'dark'
    style = 'set2'
    style = 'set4'
    # style = 'mpl_set'
    # style = 'solarized'
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
    name = 'qual_color_' + style
    plt.title(title)
    # fig.savefig('/home/schowell/Dropbox/gw_phd/%s.eps' % name)
    # fig.savefig('/home/schowell/Dropbox/gw_phd/%s.png' % name)
    fig.savefig('%s.eps' % name)
    fig.savefig('%s.png' % name)
    plt.show()

    title = 'Sequential Colors'
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
    if sort:
        name = 'sorted_seq_color'
    else:
        name = 'seq_color'
    plt.title(title)
    # fig.savefig('/home/schowell/Dropbox/gw_phd/%s.eps' % name)
    # fig.savefig('/home/schowell/Dropbox/gw_phd/%s.png' % name)
    fig.savefig('%s.eps' % name)
    fig.savefig('%s.png' % name)
    plt.show()

    print '\m/ >.< \m/'
