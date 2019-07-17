import math
import skimage
from skimage import data, io

#!usr/bin/python3.7.3
### Graph_Parser Class ###
# Overview: Input graph, select points to mark axis, then interpolate data
#
#

#v1 - v2
def vsub(v1:tuple, v2:tuple):
    return (v1[0]-v2[0], v1[1]-v2[1])

def dot(v1:tuple, v2:tuple):
    return (v1[0]*v2[0]+v1[1]*v2[1])

def magn(v1:tuple):
    return math.sqrt(dot(v1, v1))

# project vector a onto b
def proj(a:tuple, b:tuple):
    factor = dot(a,b)/dot(b,b)
    return tuple([factor*x for x in b])

class Axis(object):
    # Graph axis object for graph parser
    # self.loc = ((x1,y1), (x2,y2))
    # self.val = ( val1  , val2   )
    def __init__(self, start:tuple, end:tuple, values:tuple=None):
        self.loc = (start, end);
        if(values):
            self.val = values

    def set_value(self, values:tuple):
        self.val = values;
    
    def pt_approx(self, pt:tuple):
        v1   = vsub(pt         , self.loc[0])
        v2   = vsub(self.loc[1], self.loc[0])
        vlen = magn(proj(v1, v2))

        return vlen/magn(v2)*(self.val[1]-self.val[0])+self.val[0]


class Graph_Parser(object):
    # vax=vertical axis, hax=horizontal Axis
    def __init__(self, vp1:tuple , vp2:tuple, 
                       hp1:tuple , hp2:tuple,
                       vval:tuple, hval:tuple):
        self.vax = Axis(vp1, vp2, vval)
        self.hax = Axis(hp1, hp2, hval)

