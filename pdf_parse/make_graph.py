import matplotlib.pyplot as plt
import numpy as np
import math
import random

# Generate graphs
'''
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
'''

def coor_(radii, pt_num):
    theta = np.linspace(0, 2*math.pi, num=pt_num);    
    x = []
    y = []
    for i in theta:
        r = 0.8+random.random()*0.2
        x.append(math.cos(i)*r*radii)
        y.append(math.sin(i)*r*radii)
    return x, y

x1,y1 = coor_(5,100)
x2,y2 = coor_(4.8,100)
x3,y3 = coor_(5.3,100)
m1 = 'v'
m2 = 'o'
m3 = 's'
col1 = 'r'
col2 = 'b'
col3 = 'k'
plt.scatter(x1, y1, marker=m1, color=col1, facecolors = 'none')
plt.scatter(x2, y2, marker=m2, color=col2, facecolors = 'none')
#plt.scatter(x3, y3, marker=m3, color=col3, facecolors = 'none')
plt.show()
