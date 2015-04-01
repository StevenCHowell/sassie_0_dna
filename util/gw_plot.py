import matplotlib.pyplot as plt

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


    