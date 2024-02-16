#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 15 16:48:52 2022

@author: a.arroyo-ramo
Viene de los codigos pruebadiscret_lib_msh_v3.py y 2-surface_distribution.py y 3-mesh_prepare.py
"""

#clear the console
try:
    from IPython import get_ipython
    get_ipython().magic('reset -sf')
    # get_ipython().magic('clear')
except:
    pass


from matplotlib import pyplot as plt
import matplotlib
# matplotlib.use('Agg') #save figures in file without producing plots in spyder
# matplotlib.use('qt5agg') #plot in spyder
# plt.close('all')
import pdb

import os
import sys

import numpy as np
import time
import copy
import math
import lib_hyperbolic_grid as hg
from scipy.interpolate import lagrange
from lib_msh import *
from lib_stretch_func import *
from hyperbolic_grid import *
from elliptic_grid import *
from pathlib import Path
from scipy.interpolate import interp1d
from scipy.interpolate import splprep, splev, splrep

try:
    os.makedirs('./meshes/')
except:
    pass


# matlab like structures
class structtype():
    pass


def rotate_grid_counterclockwise(block,angle):
    angle = -angle*np.pi/180
    newcoords = np.zeros(block.coord.shape)
    newcoords[:,:,0] = block.coord[:,:,0]*np.cos(angle) - block.coord[:,:,1]*np.sin(angle)  
    newcoords[:,:,1] = block.coord[:,:,0]*np.sin(angle) + block.coord[:,:,1]*np.cos(angle) 
    block.coord = newcoords
    return block


if(1):    
    
    foil_path = 'CD_airfoil2c1TE_1000'                      #airfoil case (in mesh_fine)
    case_name = 'COARSEfishtail'                              #name of the mesh
    
    AoA = 15                                                             #angle of attack
    
    bckg_n = 1                                                 #number of background blocks to use
    blks_n = bckg_n+2                                          #number of blocks to use (ogrid, wake and background)
    
    dy_rel = 0.5e-3                                          #first cell size in normal direction (size of first cell)
    #o-grid divided in several sections
    #section 1
    bl_exp = 0.02                                            #o-grid expansion ratio in normal direction (DNS 0.01-0.02)
    ogrid_ny = 5#80                                            #number of points in normal direction around ogrid
    nk = 96                                                     #number of points in spanwise direction (used to estimate the 3D mesh size)
    
    # #adding other expansion and number of layers to customise BL growth    
    # bl_exp_i = [3*bl_exp]
    # ogrid_ny_i = [60]
    #adding other expansion and number of layers to customise BL growth    
    bl_exp_i = [2*bl_exp,
                2*2*bl_exp,
                2*3*bl_exp]
    ogrid_ny_i = [5,
                  20,
                  35]  
    
    chord = 1                                                   #value to modify according to the airfoil. X scale of the domain
    offset = 0                                                  #offset to translate vertically the airfoil (used to center the airfoil with respect to background mesh)
    
    wake_stretch = 'geom' #'geom', 'tanh', 'erf' or 'sinh'
    # percent_chord_startwake = 0.1 #percentage of the chord to start the wake from trailing edge
    x_wake     = 1.0
    slope_wake_press = -75#-70                                          #angle wake pressure side
    slope_wake_suct = 76#70                                            #angle wake suction side
    x_b_corner = 0.18 							#the larger the value of n_bezier the wider the wake
    n_bezier_cp = 4 							#the larger the value of n_bezier the wider the wake
    ds_max_wake = 0.02   


    #------------------------- BACKGROUND BLOCK--------------------------
    ds_mesh_rel = 0.005                                        #background mesh cell size (make it match to o-grid outer layers)
    fact = 1.05                                                  #factor to compute background mesh cell with respect to minimum size of ogrid
    out_stretch = 'tanh' #'geom', 'tanh', 'erf' or 'sinh'     #stretching function towards the boundaries
    ds_max_out  = 0.1
    n_sponge = 0
    n_zcbc = 0
    
    #SIZE OF THE BOX
    blks_in = -2.15 * chord                                        #box information for background mesh (smaller region)
    blks_out = 2.5 * chord                                      #box information for background mesh (smaller region)
    blks_bottom = -0.9 * chord                                    #box information for background mesh (smaller region)
    blks_top = 1.15 * chord                                        #box information for background mesh (smaller region)
    
    stretch_inlet = False                                       #stretches mesh towards the inlet
    stretch_len_inlet = 5 * chord                             #domain extension
    expand_ratio_inlet = 1.015                                   #expansion ratio
    
    stretch_outlet = False                                      #stretches mesh towards the outlet 
    stretch_len_outlet = 4.5 * chord                           #domain extension
    expand_ratio_outlet = 1.015                                  #expansion ratio
    
    stretch_top = False                                         #stretches mesh towards the top
    stretch_len_top = 6 * chord                              #domain extension
    expand_ratio_top = 1.015                                     #expansion ratio
    
    stretch_bottom = False                                      #stretches mesh towards the bottom
    stretch_len_bottom = 6 * chord                           #domain extension
    expand_ratio_bottom = 1.015                                  #expansion ratio


#%% MESH SETUP

# name of the case
Path('./meshes/'+case_name).mkdir(parents=True, exist_ok=True)
mesh_info = open('meshes/'+case_name+'/'+case_name+'info.txt',"w+")

# camber file input
camb_n = 300                # number of points in the camber file (fixed value, help for HiPSTAR to set up)
scale = 1/chord

profs_n = 0

profs_n = profs_n + 1  
foil1 = np.loadtxt('./geom_fine/{}.dat'.format(foil_path))



grid = [structtype() for i in range(0,blks_n+1)]

mesh_info.write('Path\n')
mesh_info.write('foil_path = '+foil_path+'\n')
mesh_info.write('Constants\n')
mesh_info.write('chord = '+str(chord)+'\n')
mesh_info.write('scale = '+str(scale)+'\n')
mesh_info.write('offset = '+str(offset)+'\n')

mesh_info.write('\nMesh set up\n')
mesh_info.write('nk = '+str(nk)+'\n')
mesh_info.write('camb_n = '+str(camb_n)+'\n')

mesh_info.write('dy_rel = '+str(dy_rel)+'\n')
mesh_info.write('dy = '+str(dy_rel*scale)+'\n')
mesh_info.write('bl_exp = '+str(bl_exp)+'\n')
mesh_info.write('ogrid_ny = '+str(ogrid_ny)+'\n')
mesh_info.write('adding other expansion and number of layers to customise BL growth\n')
mesh_info.write('bl_exp_i = '+str(bl_exp_i)+'\n')
mesh_info.write('ogrid_ny_i = '+str(ogrid_ny_i)+'\n')

mesh_info.write('wake_stretch = '+wake_stretch+'\n')
mesh_info.write('ds_max_wake = '+str(ds_max_wake)+'\n')
mesh_info.write('x_b_corner = '+str(x_b_corner)+'\n')
mesh_info.write('slope_wake_press = '+str(slope_wake_press)+'\n')
mesh_info.write('slope_wake_suct = '+str(slope_wake_suct)+'\n')
mesh_info.write('n_bezier_cp = '+str(n_bezier_cp)+'\n')

mesh_info.write('\nBoundaries/Background\n')
mesh_info.write('fact = '+str(fact)+'\n')
mesh_info.write('stretch_len_outlet = '+str(stretch_len_outlet)+'\n')
mesh_info.write('stretch_len_inlet = '+str(stretch_len_inlet)+'\n')
mesh_info.write('stretch_top = '+str(stretch_len_top)+'\n')
mesh_info.write('stretch_len_bottom = '+str(stretch_len_bottom)+'\n')
mesh_info.write('out_stretch = '+out_stretch+'\n')

mesh_info.write('ds_max_out = '+str(ds_max_out)+'\n')
mesh_info.write('n_sponge = '+str(n_sponge)+'\n')
mesh_info.write('n_zcbc = '+str(n_zcbc)+'\n')

mesh_info.write('x_wake = '+str(x_wake)+'\n')












#%% GENERATE OGRID AIRFOIL MESH
print("-------- Ogrid block --------")

#----------------------------First o-grid section-----------------------------

    

x_new = foil1[:,0]
y_new = foil1[:,1]

# if not (foil1[0,:] == foil1[-1,:]).all():
#     x_new = np.append(x_new,x_new[0])
#     y_new = np.append(y_new,y_new[0])

x_new = x_new-x_new[0]
y_new = y_new-y_new[0] + offset        #add offset to center the airfoil

y_new = y_new + offset        #add offset to center the airfoil

x_new = x_new*scale
y_new = y_new*scale

foil1_x = x_new
foil1_y = y_new

# generate hyperbolic mesh (OGRID)
foil = np.column_stack([x_new,y_new])

dy = dy_rel*scale
val = lambda i : dy+dy*(i-6)*bl_exp if i>5 else dy
delta_s = [val(i) for i in range(ogrid_ny-1)] #list of cell thicknesses used to expand orthogonal grid

edge_0 = foil

#generate an orthogonal grid from an initial curve and a list of cell thicknesses
xy = hg.hyperbolic_grid(edge_0, delta_s, 'periodic',
                        alpha=2.0, eps_e=0.2, v_a=0.0,
                        alpha_mode='linear',
                        dissip_type='linear',
                        dist_sat=np.inf,
                        overlap=False)
# nij = np.shape(xy)

ds_new = np.sqrt(np.diff(x_new)**2 + np.diff(y_new)**2)
print('o-grid mesh info')
mesh_info.write('\no-grid mesh info\n')
print('num points, 1st dy size: %d, %8.6f' % (len(x_new),dy))
mesh_info.write('num points, 1st dy size: %d, %8.6f\n' % (len(x_new),dy))
print('max aspect ratio: %8.6f' % (np.max(np.diff(x_new)/dy)))
mesh_info.write('max aspect ratio: %8.6f' % (np.max(np.diff(x_new)/dy)))
print('min, max & max/min ds: %8.6f, %8.6f, %8.6f' % (ds_new.min(),ds_new.max(),ds_new.max()/ds_new.min()))
mesh_info.write('min, max & max/min ds: %8.6f, %8.6f, %8.6f\n\n' % (ds_new.min(),ds_new.max(),ds_new.max()/ds_new.min()))
print('')

#second/third/... ogrid sections to custom control BL growth
for i in range(len(bl_exp_i)):
    
    x_new2 = xy[:,-1,0]
    y_new2 = xy[:,-1,1]
    
    # generate hyperbolic mesh (OGRID)
    exterior = np.column_stack([x_new2,y_new2])
    
    dy2 = abs(min(xy[:,-1,1] - xy[:,-2,1])) 
    
    
    val2 = lambda j : dy2+dy2*j*bl_exp_i[i]
    delta_s2 = [val2(j) for j in range(ogrid_ny_i[i])] #list of cell thicknesses used to expand orthogonal grid
    
    edge_02 = exterior
    
    #generate an orthogonal grid from an initial curve and a list of cell thicknesses
    xy2 = hg.hyperbolic_grid(edge_02, delta_s2, 'periodic',
                            alpha=2.0, eps_e=0.2, v_a=0.0,
                            alpha_mode='linear',
                            dissip_type='linear',
                            dist_sat=np.inf,
                            overlap=False)
    
    xy = np.append(xy, xy2[:,1:,:], axis=1)



nij = np.shape(xy)

grid[-2].ni = nij[0]
grid[-2].nj = nij[1]
grid[-2].x = xy[:,:,0]
grid[-2].y = xy[:,:,1]


Block_ogrid = Block(np.stack((grid[-2].x,grid[-2].y),axis=2))

#%%
############################## Wake block ##############################

print("-------- Wake block --------")

# Split the O-grid in two blocks

#compute the slope of the airfoil normals
grid[-2].slopes = np.arctan2((grid[-2].y[:,-1]-grid[-2].y[:,0]),(grid[-2].x[:,-1]-grid[-2].x[:,0]))*180/np.pi



# slope_wake_press = -65-AoA        #in case that it is desired to rotate first the ogrid and then create the wake 
# slope_wake_suct = 75-AoA          #in case that it is desired to rotate first the ogrid and then create the wake 
for i1 in range(grid[-2].ni):
    if grid[-2].slopes[i1] <= slope_wake_press: break
for i2 in range(grid[-2].ni-1,-1,-1):
    if grid[-2].slopes[i2] >= slope_wake_suct: break


#create the portion of the wake that was part of the ogrid
ogrid_2 = structtype()
ogrid_2.x = np.concatenate((grid[-2].x[i2:,:],grid[-2].x[:i1,:]),axis=0)
ogrid_2.y = np.concatenate((grid[-2].y[i2:,:],grid[-2].y[:i1,:]),axis=0)
ogrid_2.ni = ogrid_2.x.shape[0]
ogrid_2.nj = ogrid_2.x.shape[1]

Block_ogrid_2 = Block(np.stack((ogrid_2.x,ogrid_2.y),axis=2))


#update the ogrid by removing the elements that are included in the wake
grid[-2].x = grid[-2].x[i1:i2,:]
grid[-2].y = grid[-2].y[i1:i2,:]
grid[-2].ni = grid[-2].x.shape[0]
grid[-2].nj = grid[-2].x.shape[1]
Block_ogrid = Block(np.stack((grid[-2].x,grid[-2].y),axis=2))



# Create the initial wake block

nj = 8
n_blend = 40


# x_bezier_cp = np.linspace(x_b_corner,x_back,n_bezier_cp)
x_bezier_cp = np.linspace(x_b_corner,x_wake,n_bezier_cp)


ds_max_wake = np.min([ds_max_wake, 4*np.max(np.sqrt((Block_ogrid[:,-1,0]-Block_ogrid[:,-2,0])**2 + (Block_ogrid[:,-1,1]-Block_ogrid[:,-2,1])**2))])
print('New ds_max_wake:',ds_max_wake)
mesh_info.write('ds_max_wake (recomputed)= '+str(ds_max_wake)+'\n')


xm = x_wake

for i in range(ogrid_2.ni):

    x_i = ogrid_2.x[i,-nj:]
    y_i = ogrid_2.y[i,-nj:]

    s = np.arange(1,nj+1)

    ###

    fx = lagrange(s,x_i)
    fy = lagrange(s,y_i)

    s = np.arange(nj,nj+n_blend)
    coord = np.empty((len(s),2))
    coord[:,0] = fx(s)
    coord[:,1] = fy(s)

    c_i = Curve(coord=coord)

    ###

    
    p0 = np.array([ogrid_2.x[i,-2],ogrid_2.y[i,-2]])
    p1 = np.array([ogrid_2.x[i,-1],ogrid_2.y[i,-1]])

    n = p1 - p0
    d = np.linalg.norm(n)
    n /= d

    a = (x_b_corner-p1[0])/n[0]
    y_cp = p1[1] + a*n[1]

    p_list = [p1]
    for x_cp in x_bezier_cp:
        p_list.append([x_cp,y_cp])
    cf = CurveFunc.bezier(p_list).natural()
    l = cf.length()

    if i == 0:
        s = stretch(wake_stretch,0,l,d_min=d,d_max=ds_max_wake)
        n_wake = len(s)
        b_coord = np.empty((ogrid_2.ni,n_wake,2))
    else:
        s1 = stretch(wake_stretch,0,l,d_min=d,n=n_wake)
        s2 = stretch(wake_stretch,0,l,n=n_wake,d_max=ds_max_wake)
        s = poly_blend(s1=s1,s2=s2)

    c = Curve(coord=cf(s))

    ###
    if n_blend > c.coord.shape[1]:
        n_blend = c.coord.shape[1]
        
    s = poly_blend(n_blend)

    for k in range(n_blend):
        c[k,:] = s[k]*c[k,:] + (1-s[k])*c_i[k,:]

    b_coord[i,:,:] = c.coord

Block_wake = Block(b_coord[:,1:,:])

grid[-1].ni = b_coord.shape[0]
grid[-1].nj = b_coord.shape[1]
grid[-1].x = b_coord[:,:,0]
grid[-1].y = b_coord[:,:,1]


# Improve orthogonality

coord_list = [b.coord for b in [Block_ogrid,Block_ogrid_2,Block_wake]]

intfs = [(0,'xm',1,'xp'),
         (0,'xp',1,'xm'),
         (1,'yp',2,'ym')]

(Block_new_ogrid,Block_new_ogrid_2,Block_wake) = Mesh(elliptic_grid(coord_list,intfs,n_iter=1)).blocks


for i in range(Block_ogrid.size()[0]):
    Block_ogrid[i,:,0] = poly_blend(s1=Block_ogrid[i,:,0],s2=Block_new_ogrid[i,:,0])
    Block_ogrid[i,:,1] = poly_blend(s1=Block_ogrid[i,:,1],s2=Block_new_ogrid[i,:,1])

for i in range(Block_ogrid_2.size()[0]):
    Block_ogrid_2[i,:,0] = poly_blend(s1=Block_ogrid_2[i,:,0],s2=Block_new_ogrid_2[i,:,0])
    Block_ogrid_2[i,:,1] = poly_blend(s1=Block_ogrid_2[i,:,1],s2=Block_new_ogrid_2[i,:,1])



# Remove some grid lines to avoid bad grid quality near the corner

n_crop = 4

Block_new_ogrid = np.concatenate(( Block_ogrid_2[-n_crop:,:-n_crop,:],
                             Block_ogrid[:,:-n_crop,:],
                             Block_ogrid_2[:n_crop,:-n_crop,:] ),axis=0)

Block_new_ogrid_2 = Block_ogrid_2[n_crop:-n_crop,:-n_crop,:]

Block_new_wake = np.concatenate(( Block_ogrid_2[n_crop:-n_crop,-n_crop:,:],
                            Block_wake[n_crop:-n_crop,:,:] ),axis=1)

Block_ogrid.coord = Block_new_ogrid
Block_ogrid_2.coord = Block_new_ogrid_2
Block_wake.coord = Block_new_wake


# Add overlapping points for overset grids method

n_ov = 8

Block_ogrid.coord = np.concatenate(( Block_ogrid_2[-n_ov:,:],
                               Block_ogrid[:,:,:],
                               Block_ogrid_2[:n_ov,:,:] ),axis=0)

Block_wake.coord = np.concatenate((Block_ogrid_2.coord,Block_wake.coord),axis=1)



############################## Rotate grids ############################
Block_ogrid = rotate_grid_counterclockwise(Block_ogrid,angle=AoA)
Block_wake = rotate_grid_counterclockwise(Block_wake,angle=AoA)

grid[-2].ni = Block_ogrid.coord.shape[0]
grid[-2].nj = Block_ogrid.coord.shape[1]
grid[-2].x = Block_ogrid.coord[:,:,0]
grid[-2].y = Block_ogrid.coord[:,:,1]

grid[-1].ni = Block_wake.coord.shape[0]
grid[-1].nj = Block_wake.coord.shape[1]
grid[-1].x = Block_wake.coord[:,:,0]
grid[-1].y = Block_wake.coord[:,:,1]


#%% GENERATE BACKGROUND MESH

print("-------- Background block --------")

xin  = blks_in      #limit of inlet x
xout = blks_out      #limit of outlet x
ytop = blks_top
ybottom = blks_bottom

ds_mesh = ds_mesh_rel*scale
ni = int((xout-xin)/ds_mesh)
nj = int((ytop-ybottom)/ds_mesh)
    
xs = np.linspace(xin, xin+(ni+1)*ds_mesh, ni)
ys = np.linspace(ybottom,ybottom+(nj+1)*ds_mesh,nj)

if stretch_inlet == True :
    dom_in = blks_in-stretch_len_inlet                                  #box information for background mesh boundaries
else :
    dom_in = blks_in
if stretch_outlet == True :
    dom_out = blks_out+stretch_len_outlet                               #box information for background mesh boundaries
else :
    dom_out = blks_out
if stretch_top == True :
    dom_top = blks_top+stretch_len_top                                  #box information for background mesh boundaries
else :
    dom_top = blks_top
if stretch_bottom == True :
    dom_bottom = blks_bottom-stretch_len_bottom                         #box information for background mesh boundaries
else :
    dom_bottom = blks_bottom

# outlet expansion
ni_ext = 0
if(stretch_outlet):
    ds_ext = xs[-1]-xs[-2]
    ds_ext2 = ds_ext
    ext = 0
    xs_ext = []
    while(ext<=stretch_len_outlet):
        ni_ext = ni_ext + 1
        ext = ext + ds_ext2
        xs_ext.append(xs[-1]+ext)
        if len(xs_ext) < 2 :
            ds_ext2 = ext*expand_ratio_outlet
        else :
            ds_ext2 = (xs_ext[-1]-xs_ext[-2])*expand_ratio_outlet
        
    blks_out = xs_ext[-1]
    ni = ni + ni_ext     
    xs_ext = np.array(xs_ext)
    xs = np.concatenate((xs,xs_ext),axis=None)

print('ni_ext outlet = %d\n' % (ni_ext))
mesh_info.write('ni_ext outlet = %d\n' % (ni_ext))


# inlet expansion
ni_ext = 0
if(stretch_inlet):
    ds_ext = xs[1]-xs[0]
    ds_ext2 = ds_ext
    ext = 0
    xs_ext = []
    while(ext<=stretch_len_inlet):
        ni_ext = ni_ext + 1
        ext = ext + ds_ext2
        xs_ext.append(xs[0]-ext)
        if len(xs_ext) < 2 :
            ds_ext2 = ext*expand_ratio_outlet
        else :
            ds_ext2 = (xs_ext[-2]-xs_ext[-1])*expand_ratio_inlet
        
    blks_in = xs_ext[-1]
    ni = ni + ni_ext     
    xs_ext = np.array(xs_ext)
    xs = np.concatenate((np.flip(xs_ext),xs),axis=None)

print('ni_ext inlet = %d\n' % (ni_ext))
mesh_info.write('ni_ext inlet = %d\n' % (ni_ext))



# top expansion
nj_ext = 0
if(stretch_top):
    ds_ext = ys[-1]-ys[-2]
    ds_ext2 = ds_ext
    ext = 0
    ys_ext = []
    while(ext<=stretch_len_top):
        nj_ext = nj_ext + 1
        ext = ext + ds_ext2
        ys_ext.append(ys[-1]+ext)
        if len(ys_ext) < 2 :
            ds_ext2 = ext*expand_ratio_top
        else :
            ds_ext2 = (ys_ext[-1]-ys_ext[-2])*expand_ratio_top
        
    blks_top = ys_ext[-1]
    nj = nj + nj_ext     
    ys_ext = np.array(ys_ext)
    ys = np.concatenate((ys,ys_ext),axis=None)

print('nj_ext top = %d\n' % (nj_ext))
mesh_info.write('nj_ext top = %d\n' % (nj_ext))


# bottom expansion
nj_ext = 0
if(stretch_bottom):
    ds_ext = ys[1]-ys[0]
    ds_ext2 = ds_ext
    ext = 0
    ys_ext = []
    while(ext<=stretch_len_bottom):
        nj_ext = nj_ext + 1
        ext = ext + ds_ext2
        ys_ext.append(ys[0]-ext)
        if len(ys_ext) < 2 :
            ds_ext2 = ext*expand_ratio_bottom
        else :
            ds_ext2 = (ys_ext[-2]-ys_ext[-1])*expand_ratio_bottom
        
    blks_bottom = ys_ext[-1]
    nj = nj + nj_ext     
    ys_ext = np.array(ys_ext)
    ys = np.concatenate((np.flip(ys_ext),ys),axis=None)

print('nj_ext bottom = %d\n' % (nj_ext))
mesh_info.write('nj_ext bottom = %d\n' % (nj_ext))

#add expansion until reaching the boundaries
   
if dom_out > xs[-1] :
    ni_ext = 0
    ds_ext = xs[-1]-xs[-2]
    ext = 0
    xs_ext2 = []
    while ext < dom_out :
        ni_ext += 1
        ext = xs[-1] + ds_ext*ni_ext
        xs_ext2.append(ext)
    xs_ext2 = np.array(xs_ext2)
    xs = np.concatenate((xs,xs_ext2),axis=None)
        
    
if dom_in < xs[0] :
    ni_ext = 0
    ds_ext = xs[1]-xs[0]
    ext = 0
    xs_ext2 = []
    while ext > dom_in :
        ni_ext += 1
        ext = xs[0] - ds_ext*ni_ext
        xs_ext2.append(ext)
    xs_ext = np.array(xs_ext2)
    xs = np.concatenate((np.flip(xs_ext2),xs),axis=None)
    


if dom_top > ys[-1] :
    nj_ext = 0
    ds_ext = ys[-1]-ys[-2]
    ext = 0
    ys_ext2 = []
    while ext < dom_top :
        nj_ext += 1
        ext = ys[-1] + ds_ext*nj_ext
        ys_ext2.append(ext)
    ys_ext2 = np.array(ys_ext2)
    ys = np.concatenate((ys,ys_ext2),axis=None)


if dom_bottom < ys[0] :
    nj_ext = 0
    ds_ext = ys[1]-ys[0]
    ext = 0
    ys_ext2 = []
    while ext > dom_bottom :
        nj_ext += 1
        ext = ys[0] - ds_ext*nj_ext
        ys_ext2.append(ext)
    ys_ext = np.array(ys_ext2)
    ys = np.concatenate((np.flip(ys_ext2),ys),axis=None)




# make mesh 2d
x, y = np.meshgrid(xs, ys, sparse=False, indexing='ij')

ni = x.shape[0]
nj = y.shape[1]

grid[1].ni = ni
grid[1].nj = nj
grid[1].x = x
grid[1].y = y

# Create block
Block_cartesian = Block(np.stack((grid[1].x,grid[1].y),axis=2))




#%% MESH SUMMARY    
alln_2d = 0
alln_3d = 0
print('-------- Mesh size summary --------')
# mesh_info.write('\nmesh size summary\n')
for nb in range(1,blks_n+1):
    print('block {}'.format(nb))
    mesh_info.write('block {}\n'.format(nb))
    print(grid[nb].ni,grid[nb].nj,grid[nb].ni*grid[nb].nj)
    mesh_info.write(str(grid[nb].ni)+'\t'+str(grid[nb].nj)+'\t'+str(grid[nb].ni*grid[nb].nj)+'\n')
    print(np.shape(grid[nb].x))
    mesh_info.write(str(np.shape(grid[nb].x))+'\n')
    alln_2d = alln_2d + grid[nb].ni*grid[nb].nj
    alln_3d = alln_3d + grid[nb].ni*grid[nb].nj*nk

print('')
mesh_info.write('')
print('2D mesh size: {} M'.format(alln_2d/1e6))
mesh_info.write('2D mesh size: {} M\n'.format(alln_2d/1e6))
print('3D mesh size: {} M'.format(alln_3d/1e6))
mesh_info.write('3D mesh size: {} M\n'.format(alln_3d/1e6))
print('')



#%% generate camb.dat

print('Generating camb.dat file')

blade_x = foil1_x
blade_y = foil1_y

LEind = np.argmin(blade_x)
TEind = 0

rot_ang = 0
x_new = blade_x
y_new = blade_y
while(abs(y_new[TEind])>1e-3):
    rt = 180*abs(y_new[TEind]/(x_new[LEind]-x_new[TEind]))/(2*np.pi)
    if(y_new[TEind]>0):
        rot_ang = rot_ang - rt
    else:
        rot_ang = rot_ang + rt
    theta = rot_ang*math.pi/180
    x_new = blade_x*np.cos(theta) - blade_y*np.sin(theta)
    y_new = blade_x*np.sin(theta) + blade_y*np.cos(theta)

PS_x = x_new[LEind:0:-1]
PS_y = y_new[LEind:0:-1]

SS_x = x_new[LEind:-1]
SS_y = y_new[LEind:-1]

camb_x = np.linspace(x_new[LEind],x_new[TEind],camb_n)

fps = interp1d(PS_x, PS_y, kind='linear', fill_value='extrapolate')
fss = interp1d(SS_x, SS_y, kind='linear', fill_value='extrapolate')

camb_y1 = np.array(fps(camb_x))
camb_y2 = np.array(fss(camb_x))
camb_y = 0.5*(camb_y1 + camb_y2)

theta = -rot_ang*math.pi/180
cambx1 = camb_x*np.cos(theta) - camb_y*np.sin(theta)
camby1 = camb_x*np.sin(theta) + camb_y*np.cos(theta)

cambx1 = cambx1[2:]
camby1 = camby1[2:]



#%% write camb.dat file

print('Writing camb.dat file')

mesh_info.close()

mesh_info = open('meshes/'+case_name+'/'+case_name+'info.txt',"a")
# block = get_mesh_CPUs(case_name,mesh_info,HPC,nk,cellsCPU=30000)
mesh_info.close()



cambx = cambx1
camby = camby1    
 
camb_n = len(cambx)
path = './meshes/' + case_name + '/camb.dat'
f = open(path,'w')
f.write('{} \n'.format(camb_n))
for i in range(camb_n):
    f.write('{} {} \n'.format(cambx[i],camby[i]))

f.close()

#%% save grids
# ############################## Save grids ##############################

# print("-------- Save grids --------")

# mesh = Mesh([c_block,ogrid,wake])
# mesh.saveas('z_r_grid')

print('Writing meshes to file')
for nb in range(1,blks_n+1):
    path = './meshes/' + case_name + '/z_r_grid_{0:d}.dat'.format(nb)
    
    f = open(path,'w')
    f.write('{} {} \n'.format(grid[nb].ni,grid[nb].nj))
    for j in range(grid[nb].nj):
        for i in range(grid[nb].ni):
            f.write('{} {} \n'.format(grid[nb].x[i,j],grid[nb].y[i,j]))

    f.close()


#%% PLOT MESHES

print('Plotting meshes')

stp = 1 #number of "jumps" to plot the mesh
cols = ['','k','r','b','g','m','y','c','k','r','b','g','m','y','c']

plt.figure(1)
plt.clf()

plt.title('Mesh every %d gridline'%stp)

for nb in range(1,blks_n+1):
#for nb in range(2,blks_n+1):
# for nb in range(1,blks_n):
    plt.plot(grid[nb].x[ :, 0], grid[nb].y[ :, 0],'-', color=cols[nb],label='Block %d'%nb)
    plt.plot(grid[nb].x[ :,-1], grid[nb].y[ :,-1],'-', color=cols[nb])
    plt.plot(grid[nb].x[ 0, :], grid[nb].y[ 0, :],'-', color=cols[nb])
    plt.plot(grid[nb].x[-1, :], grid[nb].y[-1, :],'-', color=cols[nb])
    
    for i in range(0,grid[nb].ni,stp):
        plt.plot(grid[nb].x[i,:],grid[nb].y[i,:],'-', color=cols[nb])
    for i in range(0,grid[nb].nj,stp):
        plt.plot(grid[nb].x[:,i],grid[nb].y[:,i],'-', color=cols[nb])

plt.axis('equal')
plt.xlabel('x coordinate, m')
plt.ylabel('y coordinate, m')
# plt.legend(loc='upper right')
plt.ticklabel_format(style='plain')
fig_path = './meshes/' + case_name + '/mesh.png'
plt.savefig(fig_path, dpi=300)
plt.show()


# # matplotlib.use('qt5agg')
# Mesh([Block_cartesian,Block_ogrid,Block_wake]).show()
# plt.title('Airfoil points: {}'.format(foil1.shape[0]))
# plt.xlabel('x coordinate, m')
# plt.ylabel('y coordinate, m')
# plt.tight_layout()


# Mesh([Block_cartesian,Block_wake,Block_ogrid]).show()
# plt.title('Airfoil points: {}'.format(foil1.shape[0]))
# plt.xlabel('x coordinate, m')
# plt.ylabel('y coordinate, m')
# plt.tight_layout()
