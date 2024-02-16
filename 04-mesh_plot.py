#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Andrea Arroyo
"""
#clear the console
try:
    from IPython import get_ipython
    get_ipython().magic('reset -sf')
    get_ipython().magic('clear')
except:
    pass



import numpy as np
import matplotlib.pylab as plt
# plt.close('all')
from pathlib import Path

import matplotlib


# matplotlib.use('Agg') #save figures in file without producing plots in spyder
# matplotlib.use('qt5agg') #plot in spyder
import glob
import argparse
import pdb

# matlab like structures
class structtype():
    pass


parser = argparse.ArgumentParser(description='Plot the mesh')
parser.add_argument('--stp', type=int, default=1, help='Step to plot the mesh')
parser.add_argument('--meshname', type=str, default='meshname', help='Names of meshes, separated by comma')
parser.add_argument('--monitorsfile', type=str, default='.dat', help='Names monitor points file, including extension')

args = parser.parse_args()

case_name = args.meshname.split(',')

monitorsfile = args.monitorsfile
monitors = False
if len(monitorsfile) > 4:
    monitors = True
    monitors_loc = pd.read_csv(monitorsfile,sep=' ', header='infer')


stp = args.stp

cols = ['','k','r','g','b','m','y','c','k','r','b','g','m','y','c']


for c in range(len(case_name)):
    
    Path('./meshes/'+case_name[c]).mkdir(parents=True, exist_ok=True)

    path = './meshes/' + case_name[c]

    z_r_files = glob.glob(path+'/z_r_grid_*')
    nb = len(z_r_files)
    grid = [structtype() for i in range(0,nb+1)]

    for ib in range(1,nb+1):
        print(ib)
        grid_path = path+'/z_r_grid_{0:d}.dat'.format(ib)
    
        file = open(grid_path, "r")
        data = file.readline().split()
        nx = int(data[0])
        ny = int(data[1])
            
        x = np.zeros([nx,ny])
        y = np.zeros([nx,ny])
        
        for j in range(ny):
            for i in range(nx):
                data = file.readline().split()
                x[i,j] = float(data[0])
                y[i,j] = float(data[1])
        
        grid[ib].x = x
        grid[ib].y = y
        grid[ib].ni = nx
        grid[ib].nj = ny
        #pdb.set_trace()
        
        file.close()
    
    
    plt.figure()
    plt.clf()
    
    # plt.title('every %d gridline '%stp+case_name[c])
    
    for ib in range(1,nb+1):
    #for ib in range(2,nb+1):
    #for ib in range(3,nb+1):
        plt.plot(grid[ib].x[ :, 0], grid[ib].y[ :, 0],'-', color=cols[ib],label='Block %d'%ib)
        plt.plot(grid[ib].x[ :,-1], grid[ib].y[ :,-1],'-', color=cols[ib])
        plt.plot(grid[ib].x[ 0, :], grid[ib].y[ 0, :],'-', color=cols[ib])
        plt.plot(grid[ib].x[-1, :], grid[ib].y[-1, :],'-', color=cols[ib])
        
        for i in range(0,grid[ib].ni,stp):
            plt.plot(grid[ib].x[i,:],grid[ib].y[i,:],'-', color=cols[ib])
        for i in range(0,grid[ib].nj,stp):
            plt.plot(grid[ib].x[:,i],grid[ib].y[:,i],'-', color=cols[ib])
        
        plt.axis('equal')
        plt.xlabel('X coordinate, m')
        plt.ylabel('Y coordinate, m')
        plt.gca().ticklabel_format(axis='both', style='plain')
        #plt.legend(loc='best')
        #plt.tight_layout()
        #plt.show()
    
    try:
        plt.plot(monitor_points[:,0],monitor_points[:,1],'ob')
    except:
        pass

    
plt.axis('equal')
plt.xlabel('$x/c$, non-dimensional')
plt.ylabel('$y/c$, non-dimensional')
plt.legend(loc='best')

plt.tight_layout()
plt.show()
