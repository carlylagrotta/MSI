import glob, sys, os
#import torch as tc
import numpy as np
import pandas as pd
import skimage as sk
from skimage import io
from skimage import color

# Written by Daniel Lee (2019)
# PGen class
# Purpose: Generate training data for the CNN to detect line graph
# patterns
#

if(len(sys.argv) > 3):
    load_dir = sys.argv[1]
    save_dir = sys.argv[2]
    lbl      = sys.argv[3]
    if(len(sys.argv) > 4):
        im_size = int(sys.argv[4])
    else:
        im_size = 30
        
else:
    print("Please run generate_pattern as follows: python3 generate_pattern.py <load-path> <save-path>")

class PGen(object):
    def __init__(self, ldir:str, sdir:str, labels:str, size:int=30):
        self.load_dir = ldir
        self.save_dir = sdir
        self.lbl      = labels
        self.win_size = size

    def get_file(self, fdir:str):
        ext = ['.png','.jpg','.pdf','.jpeg','.tiff']
        llist = os.listdir(fdir)
        flist = []

        for f in llist:
            for e in ext:
                if(f.endswith(e)):
                    flist.append(f)  
                    break

        return flist

    def window_data(self, fpath:str, c:int=0, step:int=1):
        # Extract image and convert to rgb
        img = io.imread(fpath)
        img = color.rgb2gray(img)
        img = np.array(img)
        
        if(img.shape[0] < self.win_size and img.shape[1] < self.win_size):
            # Check if size matches
            print(fpath+" is not large enough image")
        else:
            # Iterate through
            
            for row in range(0, img.shape[0]-self.win_size+1, step):
                for col in range(0, img.shape[0]-self.win_size+1, step):
                    tmp_im = img[row:row+self.win_size, col:col+self.win_size]
                    if(np.sum(np.sum(tmp_im))==self.win_size*self.win_size and
                            tmp_im.shape==(self.win_size,self.win_size)):
                        c += 1
                        filename  = os.path.join(self.save_dir, str(c)+'.png')
                        io.imsave(fname=filename, arr=tmp_im)
                             
        return c

    def generate_data(self, c:int=0, c_lim:int=300000):
        flist = self.get_file(self.load_dir)
        c_dif = 0;
        c_before= c;

        # Iterate through files
        for fname in flist:
            fpath = os.path.join(self.load_dir, fname)
            if(c_dif > c_lim):
                break
            c = self.window_data(fpath, c)
            c_dif = c-c_before;
            

def main():
    train_data = PGen(load_dir, save_dir, lbl, im_size)

    train_data.generate_data(c=0)

if __name__=='__main__':
    main()
