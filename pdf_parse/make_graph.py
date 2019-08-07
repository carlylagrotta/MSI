import matplotlib.pyplot as plt
import numpy as np
import math
import random

# Generate graphs

def gplt(radii, mark='', lstyle='-', pt_num=1000, col='b', scatter=False):
    theta = np.linspace(0, 2*math.pi, num=pt_num);
    x = []
    y = []
    for i in theta:
        r = 0.8+random.random()*0.2
        x.append(math.cos(i)*r*radii)
        y.append(math.sin(i)*r*radii)
    
   
    if(scatter):
        plt.plot(x,y, lw=0, marker=mark, fillstyle='none', color=col,
                linestyle=lstyle)
    else:
        if(mark!=''):
            plt.plot(x, y, marker=mark,linestyle=lstyle, color=col)
        else:
            plt.plot(x, y,linestyle=lstyle, color=col)

    plt.show(block=True)


radii = 5
gplt(radii,lstyle='--', col='k', pt_num=200)
gplt(radii,mark='v', col='r', lstyle=None, scatter=True, pt_num=200)
gplt(radii,mark='v', col='g', lstyle=None, scatter=True, pt_num=200)
gplt(radii,mark='v', col='b', lstyle=None, scatter=True, pt_num=200)
gplt(radii,mark='v', col='k',lstyle=None, scatter=True, pt_num=200)
gplt(radii,mark='^', col='r', lstyle=None, scatter=True, pt_num=200)
gplt(radii,mark='^', col='g', lstyle=None, scatter=True, pt_num=200)
gplt(radii,mark='^', col='b', lstyle=None, scatter=True, pt_num=200)
gplt(radii,mark='^', col='k',lstyle=None, scatter=True, pt_num=200)









