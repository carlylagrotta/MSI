import math
import skimage
from skimage import data, io
import numpy as np
import matplotlib.pyplot as plt

#!usr/bin/python3.7.3
### Graph_Parser Class ###
# Overview: Input graph, select points to mark axis, then interpolate data
#   Requirements: sci-kit image, matplotlib, numpy, math
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
            self.val = sorted(values)
            
    # acccessors        
    def get_ycoord(self):
        return sorted((self.loc[0][0], self.loc[1][0]))

    def get_xcoord(self):
        return sorted((self.loc[0][1], self.loc[1][1]))

    def get_range(self):
        return self.val

    # mutator
    def set_value(self, values:tuple):
        self.val = values;
    
    def pt_approx(self, pt:tuple):
        # Approximate value at given point relative to axis
        v1   = vsub(pt         , self.loc[0])
        v2   = vsub(self.loc[1], self.loc[0])
        vlen = magn(proj(v1, v2))
        return vlen/magn(v2)*(self.val[1]-self.val[0])+self.val[0]


class Graph_Parser(object):
    # vax=vertical axis, hax=horizontal Axis
    def __init__(self, vp1:tuple , vp2:tuple, 
                       hp1:tuple , hp2:tuple,
                       vval:tuple, hval:tuple,
                       pt:tuple = None,
                       path:str = None):
        self.vax      = Axis(vp2, vp1, vval)
        self.hax      = Axis(hp1, hp2, hval)
        
        # Optional parameters
        if(pt):self.color_pt = pt;        
        if(path):self.img = np.array(io.imread(path))
        if(pt and path): self.set_color()
        
    #mutators
    def set_image(self, path:str = ''):
        self.img = np.array(io.imread(path))

    def set_vax(self, vp1:tuple, vp2:tuple, vval:tuple):
        self.vax = Axis(vp1, vp2, vval)

    def set_hax(self, hp1:tuple, hp2:tuple, hval:tuple):
        self.hax = Axis(hp1, hp2, hval)

    def set_color(self, pt=None):
        if(pt):
            self.color = self.img[pt[0], pt[1]]
        else:
            self.color = self.img[self.color_pt[0], self.color_pt[1]]
    
    # Approximate time-series      
    def get_pts(self, step:int = 1, mode:str='color'):
        if(mode == 'color'):
            xax = self.hax.get_xcoord()
            yax = self.vax.get_ycoord()
            pts = []
            for t in range(xax[0]+2, xax[1]-2, step):
                vval = 0
                tval = 0
                c = 0
                
                for v in range(yax[0], yax[1]-3):
                    # If color matches
                    if(np.array_equal(self.img[v, t], self.color)):
                        vval += self.vax.pt_approx((v, t))
                        tval += self.hax.pt_approx((v, t))
                        c += 1
                        
                if(c != 0): pts.append((tval/c, vval/c)) # Identify colors/average values found and append
            
            return pts
