#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Pawel Przytarski
Modifications: Andrea Arroyo

Distribute points along the blade based on its curvature 
Based on the method used in XFOIL developed by M. Drela
"""


#clear the console
try:
    from IPython import get_ipython
    get_ipython().magic('reset -sf')
    get_ipython().magic('clear')
except:
    pass


from scipy.interpolate import splprep, splev, splrep
import time
import numpy as np
import copy
import matplotlib.pylab as plt

import matplotlib
# matplotlib.use('Agg') #no plots in spyder
# matplotlib.use('qt5agg') #plot in spyder

plt.close('all')


import os
try:
    os.makedirs('./geom_fine/')
except:
    pass




#%% function definitions
def solve_tri(A, B, C, D, n):
    A2 = copy.copy(A)
    B2 = copy.copy(B)
    C2 = copy.copy(C)
    D2 = copy.copy(D)
    X = np.zeros(n,)
    for i in range(1,n):
        w = A2[i]/B2[i-1]
        B2[i] = B2[i] - w*C2[i-1]
        D2[i] = D2[i] - w*D2[i-1]
    X[-1] = D2[-1]/B2[-1]
    for i in range(n-2,-1,-1):
        X[i] = (D2[i] - C2[i]*X[i+1])/B2[i]
    return X

def rotate_airfoil(foil,angle) :
    print('Rotating airfoil {:d} deg'.format(angle))
    alpharad = angle*np.pi/180
    foilmodified = np.empty_like(foil)
    foilmodified[:,0] = foil[:,0]*np.cos(alpharad) + foil[:,1]*np.sin(alpharad)
    foilmodified[:,1] = -foil[:,0]*np.sin(alpharad) + foil[:,1]*np.cos(alpharad)
    return foilmodified

def scale_airfoil(foil,chord) :
    print('Scaling airfoil rate {:d}'.format(chord))
    foilmodified = chord*foil
    return foilmodified

def TE_origin(foil) :
    print('Locating origin of coordinates at TE')
    foilmodified = foil
    foilmodified = foil[:,:] - foil[0,:]
    return foilmodified

#%% default paneling parameters

curv_param = 3.5;
curv_te_ratio = 5.5;

ipfac = 10

# ----------------------- INPUT PARAMETERS --------------------------------
profile = 'CD_airfoil2'			#Name of airfoil in geom folder

n = 500                                 #number of points around the airfoil
smooth_iter = int(1e4)                  #in the range 1e4 - 1e6

chord = 1 				#chord length
TEorigin = True				#set the origin of coordinates at the trailing edge (True/False)

# -------------------------------------------------------------------------

#%% Loading airfoil


foil = np.loadtxt('./geom/{}.dat'.format(profile))

tic = time.time()

#%% Airfoil modification (scaling and location of origin of coordinates)
foil = scale_airfoil(foil,chord)
if TEorigin:
    foil = TE_origin(foil)


#%% parametrise the airfoil
Xb = foil[:,0]
Yb = foil[:,1]

DXb = np.diff(Xb)
DYb = np.diff(Yb)

DSb = np.sqrt(DXb**2 + DYb**2)
Sb = np.concatenate([[0],np.cumsum(DSb)])

nb = len(Sb);

sbref = 0.5*(Sb[-1]-Sb[0])      #reference surface of the airfoil (average value between cumulative curve pressure and suction sides)

tck, u = splprep([Xb, Yb], s=0, per=True)
Xbp, Ybp = splev(u, tck, der=1)
Xbpp, Ybpp = splev(u, tck, der=2)


#%% curvature
W5 = np.zeros(nb,)
W5 = abs((Xbp*Ybpp - Ybp*Xbpp)/(np.sqrt(Xbp**2 + Ybp**2))**3)*sbref


#%% locate leading edge
# (X-XTE,Y-YTE).(X',Y') = 0 at S = SLE
convg = 2*sbref*1e-5

xte = 0.5*(Xb[0]+Xb[-1])
yte = 0.5*(Yb[0]+Yb[-1])

le_cond = (Xb-xte)*Xbp + (Yb-yte)*Ybp
ile = np.where(le_cond<0)[0][0]

sle = u[ile]

for it in range(50):
    print(it)
    xle, yle = splev(sle, tck)
    xlep, ylep = splev(sle, tck, der=1)
    xlepp, ylepp = splev(sle, tck, der=2)

    res = (xle-xte)*xlep + (yle-yte)*ylep
    ress = xlep*xlep + ylep*ylep + (xle-xte)*xlepp + (yle-yte)*ylepp

    dsle = -res/ress
    dsle = max(dsle, -0.02*abs(xle-xte+yle-yte))
    dsle = min(dsle,  0.02*abs(xle-xte+yle-yte))
    sle = sle + dsle
    if(abs(dsle)<convg):
        break

sble = sle*Sb[-1]

curv_le = abs((xlep*ylepp - ylep*xlepp)/(np.sqrt(xlep**2 + ylep**2))**3)*sbref    


#%% average curvature around the leading edge
nk = 5
curv_sum = 0.0
for k in range(-nk,nk+1):
    frac = k/nk
    
    sbk = sble + frac*sbref/max(curv_le, 50.0)
    sk = sbk/Sb[-1]
    
    xk, yk = splev(sk, tck)
    xkp, ykp = splev(sk, tck, der=1)
    xkpp, ykpp = splev(sk, tck, der=2)
    
    curv_k = abs((xkp*ykpp - ykp*xkpp)/(np.sqrt(xkp**2 + ykp**2))**3)*sbref
    curv_sum = curv_sum + curv_k

curv_avg = curv_sum/(2*nk+1)

curv_coef = 6.0*curv_param

curv_te = curv_avg*curv_te_ratio
W5[0]  = curv_te
W5[-1] = curv_te


smth_len = max(1/max(curv_avg, 50.0), 0.25/(n/2))
smth = (smth_len*sbref)**2

#%% smooth curvature array for smoother point distribution
W1 = np.zeros(nb,)
W2 = np.zeros(nb,)
W3 = np.zeros(nb,)
W4 = np.zeros(nb,)

W2[0] = 1
W3[-1] = 1
for i in range(1,nb-1):
    ds1 = DSb[i-1]
    ds2 = 0.5*(DSb[i] + DSb[i-1])
    ds3 = DSb[i]
    W1[i] = smth*(-1/ds1)/ds2
    W2[i] = smth*( 1/ds1 + 1/ds3)/ds2 + 1
    W3[i] = smth*(-1/ds3)/ds2
W1[-1] = 0
W2[-1] = 1

W5 = solve_tri(W1, W2, W3, W5, nb)

#%%fix curvature at leading edge by modifying equations adjacent to leading edge
ile = np.where(sble<=Sb)[0][0]

if(sble == Sb[i]):
    W1[ile] = 0
    W2[ile] = 1
    W3[ile] = curv_param*curv_le
else:
    # modify equation at node just before leading edge
    ds1 = Sb[ile-1] - Sb[ile-2]
    ds2 = 0.5*(sble - Sb[ile-2])
    ds3 = sble - Sb[ile-1]
    
    W1[ile-1] = smth*(-1/ds1)/ds2
    W2[ile-1] = smth*( 1/ds1 + 1/ds3)/ds2 + 1
    W3[ile-1] = 0
    W5[ile-1] = W5[ile-1] + smth*curv_param*curv_le/(ds2*ds3)
 
    # modify equation at node just after leading edge
    ds1 = Sb[ile] - sble
    ds2 = 0.5*(Sb[ile+1] - sble)
    ds3 = Sb[ile+1] - Sb[ile]
 
    W1[ile] = 0
    W2[ile] = smth*( 1/ds1 + 1/ds3)/ds2 + 1
    W3[ile] = smth*(-1/ds3)/ds2
    W5[ile] = W5[ile] + smth*curv_param*curv_le/(ds1*ds2)
    

#%% sove for smoothed curvature array
W5 = solve_tri(W1, W2, W3, W5, nb)

curv_max = max(abs(W5))

W5 = W5/curv_max


#%% spline curvature array
spl = splrep(Sb, W5, s=0)
W6 = splev(Sb, spl, der=1)



#%% Set initial guess for node positions uniform in s.
# More nodes than specified (by factor of IPFAC) are
# temporarily used  for more reliable convergence.
nn = ipfac*n
Snew = np.zeros(nn,)

dsavg = (Sb[-1] - Sb[0])/(nn-1)
Snew[0] = Sb[0]
for i in range(1,nn): 
    Snew[i] = Sb[0] + dsavg*i

Snew[nn-1] = Sb[nb-1]


#%% Newton iteration loop for new node positions
W1 = np.zeros(nn,)
W2 = np.zeros(nn,)
W3 = np.zeros(nn,)
W4 = np.zeros(nn,)

itermax = 10000
for it in range(1,itermax):
    
    for i in range(1,nn-1):
        cv1  = splev(Snew[i-1], spl)
        cv2  = splev(Snew[i]  , spl)
        cv3  = splev(Snew[i+1], spl)
        cvs1 = splev(Snew[i-1], spl, der=1)
        cvs2 = splev(Snew[i]  , spl, der=1)
        cvs3 = splev(Snew[i+1], spl, der=1)

        cav1 = np.sqrt(cv1**2 + cv2**2)
        cav2 = np.sqrt(cv2**2 + cv3**2)
        
        cav1_s1 = cvs1*cv1/cav1
        cav1_s2 = cvs2*cv2/cav1

        cav2_s2 = cvs2*cv2/cav2
        cav2_s3 = cvs3*cv3/cav2
        
        ds1 = Snew[i] - Snew[i-1]
        ds2 = Snew[i] - Snew[i+1]
        
        f1 = curv_coef*cav1 + 1
        f2 = curv_coef*cav2 + 1
               
        # lower, main and upper diagonals
        W1[i] = -f1 + curv_coef*ds1*cav1_s1
        W2[i] =  f1 + f2 + curv_coef*(ds1*cav1_s2 + ds2*cav2_s2)
        W3[i] = -f2 + curv_coef*ds2*cav2_s3
        
        # residual requiring that (1 + C*curv)*deltaS is equal on both sides of node i
        res = ds1*f1 + ds2*f2
        W4[i] = -res

    # fix endpoints at TE
    W2[0] = 1
    W3[0] = 0
    W4[0] = 0
    W1[-1] = 0
    W2[-1] = 1
    W4[-1] = 0

    W4 = solve_tri(W1, W2, W3, W4, nn)
    
    # find under-relaxation to keep nodes from changing order
    rlx = 1
    for i in range(nn-1):
        ds = Snew[i+1] - Snew[i]
        dds = W4[i+1] - W4[i]
        ds_rat = 1.0 + rlx*dds/ds
        if(ds_rat>4):
            rlx = 3*ds/dds
        if(ds_rat<0.2):
            rlx = -0.8*ds/dds
    dmax = max(abs(W4))
    
    # update node positions
    for i in range(1,nn-1):
        Snew[i] = Snew[i] + rlx*W4[i]
    
    print(it,dmax)
    if(abs(dmax) < 1e-12):
        break
    
    if(0):
        plt.figure(51)
        plt.subplot(211)
        plt.plot(W4)
        plt.subplot(212)
        plt.plot(Snew)

    
#%% Diffuse the points to smooth out 2nd derivatives
alpha = 1.0

T0 = copy.copy(Snew)
T1 = copy.copy(Snew)

for it in range(1,smooth_iter+1):
    if(np.mod(it,1e4)==0):
        print('it: {}'.format(it))
    T1[1:-1] = T0[1:-1] + alpha*0.25*(T0[2:] - 2*T0[1:-1] + T0[:-2])
    T0 = copy.copy(T1)
Snew = copy.copy(T0)


#%% Set new panel mode coordinates
Sf = np.zeros(n,)
Xf = np.zeros(n,)
Yf = np.zeros(n,)
for i in range(n):
    ind = ipfac*i 
    Sf[i] = Snew[ind]
    u_temp = Snew[ind]/Snew[-1]
    Xf[i], Yf[i] = splev(u_temp, tck);
    
    
toc = time.time()
print('elapsed time:',toc-tic)
#%% Write distribution to file

if TEorigin:
    TEoriginstr = 'TE0'
else:
    TEoriginstr = ''
chordstr = str(chord)
chordstr = chordstr.replace('.','_')

newprofile = '{:s}c{:s}_{:s}'.format(profile,chordstr,TEoriginstr)   #profile + 'c' + TEoriginstr
path = './geom_fine/{:s}_{:d}.dat'.format(newprofile,n)
    
f = open(path,'w')
for i in range(n):
    f.write('{} {} \n'.format(Xf[i],Yf[i]))
f.close()

print('Airfoil perimeter:', Snew[-1]-Snew[0])
print('Chord length:', np.sqrt((xle-xte)**2+(yle-yte)**2))

surface_info = open('./geom_fine/{:s}_{:d}info.txt'.format(newprofile,n),"w+")
surface_info.write('curv_param: '+str(curv_param)+'\n')
surface_info.write('curv_te_ratio: '+str(curv_te_ratio)+'\n')
surface_info.write('ipfac: '+str(ipfac)+'\n')
surface_info.write('Chord length: '+str(np.sqrt((xle-xte)**2+(yle-yte)**2))+'\n')
surface_info.write('Airfoil perimeter: '+str(Snew[-1]-Snew[0])+'\n')
surface_info.write('smooth_iter: '+str(smooth_iter)+'\n')
surface_info.write('dS/dSref max: '+str(max(np.diff(Sf)/sbref))+'\n')
surface_info.write('dX/dSref max: '+str(max(np.diff(Xf)/sbref))+'\n')
surface_info.write('dY/dSref max: '+str(max(np.diff(Yf)/sbref))+'\n')
surface_info.write('ddS/dSref max: '+str(max(np.diff(np.diff(Sf))/sbref))+'\n')
surface_info.write('ddX/dSref max: '+str(max(np.diff(np.diff(Xf))/sbref))+'\n')
surface_info.write('ddY/dSref max: '+str(max(np.diff(np.diff(Yf))/sbref))+'\n')
surface_info.write('Elapsed time: '+str(toc-tic)+'\n')

surface_info.close()
#%%

plt.figure(5)
# plt.clf()
# plt.plot(Xb,Yb,'k-x')
# plt.plot(Xf,Yf,'r-x')
# plt.plot(Xb,Yb,'k-x',label='Reference')
# plt.plot(Xf,Yf,'-x',label='%d points, %d'%(n,smooth_iter))
plt.plot(Xf,Yf,'-x',label=newprofile.replace('_','\_'))
plt.legend(loc='upper right')
plt.xlabel('x')
plt.ylabel('y')
plt.axis('equal')



plt.figure(6)
#plt.clf()
# absolute 1st derivative
plt.subplot(231)
plt.title('dS')
plt.plot(np.diff(Sf),label='%d points, %d'%(n,smooth_iter))
plt.subplot(232)
plt.title('dX')
plt.plot(np.diff(Xf))
plt.subplot(233)
plt.title('dY')
plt.plot(np.diff(Yf))

# relative 1st derivative
plt.subplot(234)
plt.title('dS / Sref')
plt.plot(np.diff(Sf)/sbref,label='%d points, %d'%(n,smooth_iter))
plt.subplot(235)
plt.title('dX / Sref')
plt.plot(np.diff(Xf)/sbref)
plt.subplot(236)
plt.title('dY / Sref')
plt.plot(np.diff(Yf)/sbref)


plt.figure(7)
#plt.clf()
# absolute 2nd derivative
plt.subplot(231)
plt.title('ddS')
plt.plot(np.diff(np.diff(Sf)),label='%d points, %d'%(n,smooth_iter))
plt.subplot(232)
plt.title('ddX')
plt.plot(np.diff(np.diff(Xf)))
plt.subplot(233)
plt.title('ddY')
plt.plot(np.diff(np.diff(Yf)))

# relative 2nd derivative
plt.subplot(234)
plt.title('ddS / Sref')
plt.plot(np.diff(np.diff(Sf))/sbref,label='%d points, %d'%(n,smooth_iter))
plt.subplot(235)
plt.title('ddX / Sref')
plt.plot(np.diff(np.diff(Xf))/sbref)
plt.subplot(236)
plt.title('ddY / Sref')
plt.plot(np.diff(np.diff(Yf))/sbref)
