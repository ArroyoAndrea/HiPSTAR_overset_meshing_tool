#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Pawel Przytarski
Modifications: Andrea Arroyo
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
import matplotlib


# matplotlib.use('Agg') #no plots in spyder
# matplotlib.use('qt5agg') #plot in spyder

import lib_hyperbolic_grid as hg
import math

from pathlib import Path
from scipy.interpolate import interp1d


import os
try:
    os.makedirs('./meshes/')
except:
    pass


# matlab like structures
class structtype():
    pass


def rotate_grid_counterclockwise(meshdata,angle):
    angle = -angle*np.pi/180
    newcoords = np.zeros(meshdata.shape)
    newcoords[:,:,0] = meshdata[:,:,0]*np.cos(angle) - meshdata[:,:,1]*np.sin(angle)  
    newcoords[:,:,1] = meshdata*np.sin(angle) + meshdata*np.cos(angle) 
    return newcoords





#%% dev mesh settings REALLYCOARSE_v3
if(1):    
    
    foil_path = 'CD_airfoil2c1TE_1000'                      #airfoil case (in mesh_fine)
    case_name = 'COARSE4blks'                              #name of the mesh
    
    AoA = 8							#angle of attack to rotate
    
    bckg_n = 1                                                  #number of background blocks to use
    blks_n = bckg_n+1                                           #number of blocks to use (ogrid and background)
    
    dy_rel = 0.002                                            #first cell size in normal direction (size of first cell)
    #o-grid divided in several sections
    #section 1
    bl_exp = 0.02                                              #o-grid expansion ratio in normal direction (DNS 0.01-0.02)
    ogrid_ny = 50                                               #number of points in normal direction around ogrid
    nk = 32                                                     #number of points in spanwise direction (used to estimate the 3D mesh size)
    
    #adding other expansion and number of layers to customise BL growth    
    bl_exp_i = []
    ogrid_ny_i = []
    
    
    chord = 1#0.1#1.29 #0.08547                #value to modify according to the airfoil. X scale of the domain
    offset = 0                                                  #offset to translate vertically the airfoil (used to center the airfoil with respect to background mesh)
    
    #-------------------------FIRST BACKGROUND BLOCK--------------------------
    ds_mesh_rel = 0.009                                        #background mesh cell size (make it match to o-grid outer layers)
    
    #SIZE OF THE BOX
    blks_in = -2 * chord                                        #box information for background mesh (smaller region)
    blks_out = 2.5 * chord                                      #box information for background mesh (smaller region)
    blks_bottom = -1 * chord                                    #box information for background mesh (smaller region)
    blks_top = 1.5 * chord                                        #box information for background mesh (smaller region)
    
    stretch_inlet = False                                       #stretches mesh towards the inlet
    stretch_len_inlet = 0.5 * chord                             #domain extension
    expand_ratio_inlet = 1.05                                   #expansion ratio
    
    stretch_outlet = False                                      #stretches mesh towards the outlet 
    stretch_len_outlet = 0.75 * chord                           #domain extension
    expand_ratio_outlet = 1.05                                  #expansion ratio
    
    stretch_top = False                                         #stretches mesh towards the top
    stretch_len_top = 0.25 * chord                              #domain extension
    expand_ratio_top = 1.05                                     #expansion ratio
    
    stretch_bottom = False                                      #stretches mesh towards the bottom
    stretch_len_bottom = 0.25 * chord                           #domain extension
    expand_ratio_bottom = 1.05                                  #expansion ratio

    #-------------------------SECOND/THIRD/FOURTH,ETC BACKGROUND BLOCKS-------
    #refinement meshes
    ds_mesh_rel_i = [2e-3,
                     1e-3]                                         #background mesh cell size (make it match to o-grid outer layers)
    
    #SIZE OF THE BOX
    blks_in_i = [0.05,
                 -1.05]                                             #box information for background mesh (smaller region)
    blks_out_i = [0.65,
                  0.1]                                               #box information for background mesh (smaller region)
    blks_bottom_i = [-0.2,
                     -0.1]                                        #box information for background mesh (smaller region)
    blks_top_i = [0.1,
                  0.25]                                            #box information for background mesh (smaller region)


#%% dev mesh settings JET1000_2Blocks
if(1):    
    
    foil_path = 'CD_airfoil2c1TE_1000'                      #airfoil case (in mesh_fine)
    case_name = 'MEDIUM_2blocks'                              #name of the mesh
    
    AoA = 8							#angle of attack to rotate
    
    bckg_n = 1                                                 #number of background blocks to use
    blks_n = bckg_n+1                                           #number of blocks to use (ogrid and background)
    
    dy_rel = 0.0002                                            #first cell size in normal direction (size of first cell)
    #o-grid divided in several sections
    #section 1
    bl_exp = 0.02                                              #o-grid expansion ratio in normal direction (DNS 0.01-0.02)
    ogrid_ny = 60                                               #number of points in normal direction around ogrid
    nk = 32                                                     #number of points in spanwise direction (used to estimate the 3D mesh size)
    
    #adding other expansion and number of layers to customise BL growth    
    bl_exp_i = [3*bl_exp]
    ogrid_ny_i = [40]
    
    
    chord = 1#0.1#1.29 #0.08547                #value to modify according to the airfoil. X scale of the domain
    offset = 0                                                  #offset to translate vertically the airfoil (used to center the airfoil with respect to background mesh)
    
    #-------------------------FIRST BACKGROUND BLOCK--------------------------
    ds_mesh_rel = 0.0045                                        #background mesh cell size (make it match to o-grid outer layers)
    
    #SIZE OF THE BOX
    blks_in = -2 * chord                                        #box information for background mesh (smaller region)
    blks_out = 2.5 * chord                                      #box information for background mesh (smaller region)
    blks_bottom = -1 * chord                                    #box information for background mesh (smaller region)
    blks_top = 1.5 * chord                                        #box information for background mesh (smaller region)
    
    stretch_inlet = False                                       #stretches mesh towards the inlet
    stretch_len_inlet = 0.5 * chord                             #domain extension
    expand_ratio_inlet = 1.05                                   #expansion ratio
    
    stretch_outlet = False                                      #stretches mesh towards the outlet 
    stretch_len_outlet = 0.75 * chord                           #domain extension
    expand_ratio_outlet = 1.05                                  #expansion ratio
    
    stretch_top = False                                         #stretches mesh towards the top
    stretch_len_top = 0.25 * chord                              #domain extension
    expand_ratio_top = 1.05                                     #expansion ratio
    
    stretch_bottom = False                                      #stretches mesh towards the bottom
    stretch_len_bottom = 0.25 * chord                           #domain extension
    expand_ratio_bottom = 1.05                                  #expansion ratio

    #-------------------------SECOND/THIRD/FOURTH,ETC BACKGROUND BLOCKS-------
    #refinement meshes
    ds_mesh_rel_i = []                                         #background mesh cell size (make it match to o-grid outer layers)
    
    #SIZE OF THE BOX
    blks_in_i = []                                             #box information for background mesh (smaller region)
    blks_out_i = []                                               #box information for background mesh (smaller region)
    blks_bottom_i = []                                        #box information for background mesh (smaller region)
    blks_top_i = []                                            #box information for background mesh (smaller region)

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




grid = [structtype() for i in range(0,blks_n+1)]

mesh_info.write('Path\n')
mesh_info.write('foil_path = '+foil_path+'\n')
mesh_info.write('Constants\n')
mesh_info.write('chord = '+str(chord)+'\n')
# mesh_info.write('pitch = '+str(pitch)+'\n')
mesh_info.write('scale = '+str(scale)+'\n')
mesh_info.write('\nMesh set up\n')
mesh_info.write('dy_rel = '+str(dy_rel)+'\n')
mesh_info.write('dy = '+str(dy_rel*scale)+'\n')
mesh_info.write('bl_exp = '+str(bl_exp)+'\n')
mesh_info.write('ogrid_ny = '+str(ogrid_ny)+'\n')
mesh_info.write('adding other expansion and number of layers to customise BL growth\n')
mesh_info.write('bl_exp_i = '+str(bl_exp_i)+'\n')
mesh_info.write('ogrid_ny_i = '+str(ogrid_ny_i)+'\n')

mesh_info.write('nk = '+str(nk)+'\n')
mesh_info.write('ds_mesh_rel = '+str(ds_mesh_rel)+'\n')
mesh_info.write('ds_mesh = '+str(ds_mesh_rel*scale)+'\n')
mesh_info.write('camb_n = '+str(camb_n)+'\n')
mesh_info.write('stretch_outlet = '+str(stretch_outlet)+'\n')
mesh_info.write('stretch_len_outlet = '+str(stretch_len_outlet)+'\n')
mesh_info.write('expand_ratio_outlet = '+str(expand_ratio_outlet)+'\n')
mesh_info.write('stretch_inlet = '+str(stretch_inlet)+'\n')
mesh_info.write('stretch_len_inlet = '+str(stretch_len_inlet)+'\n')
mesh_info.write('expand_ratio_inlet = '+str(expand_ratio_inlet)+'\n')
mesh_info.write('stretch_top = '+str(stretch_len_top)+'\n')
mesh_info.write('expand_ratio_top = '+str(expand_ratio_top)+'\n')
mesh_info.write('stretch_bottom = '+str(stretch_bottom)+'\n')
mesh_info.write('stretch_len_bottom = '+str(stretch_len_bottom)+'\n')
mesh_info.write('expand_ratio_bottom = '+str(expand_ratio_bottom)+'\n')
mesh_info.write('blks_in = '+str(blks_in)+'\n')
mesh_info.write('blks_out = '+str(blks_out)+'\n')
mesh_info.write('blks_top = '+str(blks_top)+'\n')
mesh_info.write('blks_bottom = '+str(blks_bottom)+'\n')
mesh_info.write('offset = '+str(offset)+'\n')

mesh_info.write('SECOND/THIRD/FOURTH,ETC BACKGROUND BLOCKS\n')
mesh_info.write('ds_mesh_rel_i = '+str(ds_mesh_rel_i[:bckg_n-1])+'\n')
mesh_info.write('blks_in_i = '+str(blks_in_i[:bckg_n-1])+'\n')
mesh_info.write('blks_out_i = '+str(blks_out_i[:bckg_n-1])+'\n')
mesh_info.write('blks_bottom_i = '+str(blks_bottom_i[:bckg_n-1])+'\n')
mesh_info.write('blks_top_i = '+str(blks_top_i[:bckg_n-1])+'\n')


#%% GENERATE OGRID AIRFOIL MESH

#----------------------------First o-grid section-----------------------------
x_new = foil1[:,0]
y_new = foil1[:,1]

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

if AoA != 0:
    xy = rotate_grid_counterclockwise(xy,AoA)

nij = np.shape(xy)

grid[-1].ni = nij[0]
grid[-1].nj = nij[1]
grid[-1].x = xy[:,:,0]
grid[-1].y = xy[:,:,1]






#%% GENERATE BACKGROUND MESHES

#----------------------------First background-----------------------------
xin  = blks_in      #limit of inlet x
xout = blks_out      #limit of outlet x
ytop = blks_top
ybottom = blks_bottom

ds_mesh = ds_mesh_rel*scale
ni = int((xout-xin)/ds_mesh)
nj = int((ytop-ybottom)/ds_mesh)
    
xs = np.linspace(xin, xin+(ni+1)*ds_mesh, ni)
ys = np.linspace(ybottom,ybottom+(nj+1)*ds_mesh,nj)



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




#----------------------------Second/Third/Fourth background-------------------
#refinement meshes
if bckg_n > 1:
    for i in range(bckg_n-1):
        xin  = blks_in_i[i]      #limit of inlet x
        xout = blks_out_i[i]      #limit of outlet x
        ytop = blks_top_i[i]
        ybottom = blks_bottom_i[i]
        
        ds_mesh = ds_mesh_rel_i[i]*scale
        ni = int((xout-xin)/ds_mesh)
        nj = int((ytop-ybottom)/ds_mesh)
            
        xs = np.linspace(xin, xin+(ni+1)*ds_mesh, ni)
        ys = np.linspace(ybottom,ybottom+(nj+1)*ds_mesh,nj)
        
        
        # make mesh 2d
        x, y = np.meshgrid(xs, ys, sparse=False, indexing='ij')
        
        ni = x.shape[0]
        nj = y.shape[1]
        
        grid[2+i].ni = ni
        grid[2+i].nj = nj
        grid[2+i].x = x
        grid[2+i].y = y




#%% MESH SUMMARY    
alln_2d = 0
alln_3d = 0
print('mesh size summary')
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


#%% WRITE MESHES
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

plt.figure(20)
plt.clf()

plt.title('Mesh every %d gridline'%stp)

for nb in range(1,blks_n+1):

    plt.plot(grid[nb].x[ :, 0], grid[nb].y[ :, 0],'-', color=cols[nb],label='Block %d'%nb)
    plt.plot(grid[nb].x[ :,-1], grid[nb].y[ :,-1],'-', color=cols[nb])
    plt.plot(grid[nb].x[ 0, :], grid[nb].y[ 0, :],'-', color=cols[nb])
    plt.plot(grid[nb].x[-1, :], grid[nb].y[-1, :],'-', color=cols[nb])
    
    for i in range(0,grid[nb].ni,stp):
        plt.plot(grid[nb].x[i,:],grid[nb].y[i,:],'-', color=cols[nb])
    for i in range(0,grid[nb].nj,stp):
        plt.plot(grid[nb].x[:,i],grid[nb].y[:,i],'-', color=cols[nb])

plt.axis('equal')
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.legend(loc='upper right')
plt.show()

fig_path = './meshes/' + case_name + '/mesh.png'
plt.savefig(fig_path, dpi=300)






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
