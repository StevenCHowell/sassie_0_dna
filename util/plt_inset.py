# $Id$
import matplotlib.pyplot as plt
import numpy as np

def add_inset(ax,rect,axisbg='w'):
    '''
    taken from http://stackoverflow.com/a/17479417/3585557
    '''
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]  # <= Typo was here
    subax = fig.add_axes([x,y,width,height],axisbg=axisbg)
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax

def example1():
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111)
    rect = [0.2,0.2,0.7,0.7]
    ax1 = add_inset(ax,rect)
    ax2 = add_inset(ax1,rect)
    ax3 = add_inset(ax2,rect)
    plt.show()

def example2():
    fig = plt.figure(figsize=(10,10))
    axes = []
    subpos = [0.2,0.6,0.3,0.3]
    x = np.linspace(-np.pi,np.pi)
    for i in range(4):
        axes.append(fig.add_subplot(2,2,i))
    for axis in axes:
        axis.set_xlim(-np.pi,np.pi)
        axis.set_ylim(-1,3)
        axis.plot(x,np.sin(x))
        subax1 = add_inset(axis,subpos)
        subax2 = add_inset(subax1,subpos)
        subax1.plot(x,np.sin(x))
        subax2.plot(x,np.sin(x))

def example3():
    from matplotlib.gridspec import GridSpec
    fig = plt.figure(figsize=(11, 4))
    gs1 =  GridSpec(1, 2, left=0.075, right=0.75, wspace=0.1, hspace=0,
                    top=0.95)
    ax1 = plt.subplot(gs1[:, 0])
    ax_in1 = add_inset(ax1, [0.1, 0.2, 0.3, 0.3] , 'g')
    ax_in2 = add_inset(ax1, [0.5, 0.6, 0.3, 0.3] , 'b')

if __name__ == '__main__':
    example3()
    plt.show()

    example1()
    plt.show()

    example2()
    plt.show()

