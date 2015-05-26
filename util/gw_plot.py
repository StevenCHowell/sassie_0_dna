import matplotlib.pyplot as plt    
import numpy as np

def zoomout(ax, factor):

    xlim = ax.get_xlim()
    if 'log' == ax.get_yscale():
        xlim = np.log(xlim)
        xlim = (xlim[0] + xlim[1])/2 + np.array((-0.5, 0.5)) * (xlim[1] - xlim[0]) * (1 + factor) 
        xlim = np.exp(xlim)
    else:
        xlim = (xlim[0] + xlim[1])/2 + np.array((-0.5, 0.5)) * (xlim[1] - xlim[0]) * (1 + factor) 
    ax.set_xlim(xlim)

    ylim = ax.get_ylim()
    if 'log' == ax.get_yscale():
        ylim = np.log(ylim)
        ylim = (ylim[0] + ylim[1])/2 + np.array((-0.5, 0.5)) * (ylim[1] - ylim[0]) * (1 + factor) 
        ylim = np.exp(ylim)
    else:
        ylim = (ylim[0] + ylim[1])/2 + np.array((-0.5, 0.5)) * (ylim[1] - ylim[0]) * (1 + factor) 
    ax.set_ylim(ylim)


def xyplot(data, fmt='-', label=''):
    return plt.plot(data[:,0], data[:,1], fmt, label=label)

def xyerror(data, fmt='-', label=''):
    return plt.errorbar(data[:,0], data[:,1], data[:,2], fmt=fmt, label=label)

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
        return symbols[i%len(symbols)] + l
    else:
        return symbols[i%len(symbols)]

def color_order(i):
    colors = [[   0,    0,    1],
              [   0,  0.5,    0],
              [   1,    0,    0],
              [0.75,    0, 0.75],
              [   0, 0.75, 0.75],
              [0.75, 0.75,    0],
              # [   0,    1,    0],
              [0.25,    0, 0.25],
              # [   1,    1,    0],
              [   1,  0.5,    0]]
    return colors[i%len(colors)]

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
        plt.plot(x, y, s, mec=color_order(i), c=color_order(i), ms = 10,
                 mfc='none', label=str(color_order(i)))
        # plt.plot(x, y, s, ms = 10)
        # plt.scatter(x, y, s, markeredgecolor=color_order(i), facecolors='none')
    dy = 0.1
    plt.ylim([-dy, x[-1]+dy])
    leg = plt.legend(scatterpoints=1, numpoints=1)
    plt.show()
    name = 'python_symbol_color'
    fig.savefig('/home/schowell/Dropbox/gw_phd/%s.eps' % name)
    fig.savefig('/home/schowell/Dropbox/gw_phd/%s.png' % name)
    
    print '\m/ >.< \m/'