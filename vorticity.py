# -*- coding: utf-8 -*-
"""
To calculate enstrophy behind bubble and spike,
mixing layer defination:
https://search.proquest.com/docview/194682903?pq-origsite=gscholar
@author: Xin
"""
import pylab
import matplotlib.pyplot as plt
import h5py
import numpy as np
import os
import os.path

curfilePath = os.path.abspath(__file__)
curDir = os.path.abspath(os.path.join(curfilePath,os.pardir))
parentDir = os.path.abspath(os.path.join(curDir,os.pardir)) 

def high_order_gradient(fx,dx,order):
    length=len(fx)
    fxgrad = np.zeros(length)
    if order==4: 
        for i in range(2):
            fxgrad[i]=(-25*fx[i]+48*fx[i+1]-36*fx[i+2]+16*fx[i+3]-3*fx[i+4])/(12*dx)
        for i in range(2,length-2):
            fxgrad[i]=(-fx[i+2]+fx[i+1]*8-fx[i-1]*8+fx[i-2])/(12*dx)
        for i in range(length-2,length):
            fxgrad[i]=(25*fx[i]-48*fx[i-1]+36*fx[i-2]-16*fx[i-3]+3*fx[i-4])/(12*dx)
    if order==6:
        for i in range(3):
            fxgrad[i]=(-49/20*fx[i]+6*fx[i+1]-15/2*fx[i+2]+20/3*fx[i+3]-15/4*fx[i+4]+6/5*fx[i+5]-1/6*fx[i+6])/(dx)
        for i in range(3,length-3):
            fxgrad[i]=(fx[i+3]-9*fx[i+2]+45*fx[i+1]-45*fx[i-1]+9*fx[i-2]-fx[i-3])/(60*dx)
        for i in range(length-3,length):
            fxgrad[i]=(49/20*fx[i]-6*fx[i-1]+15/2*fx[i-2]-20/3*fx[i-3]+15/4*fx[i-4]-6/5*fx[i-5]+1/6*fx[i-6])/(dx)

        
    return fxgrad


#specify inout parameters here
g=1.0
inFile="tests_single_new.h5"
Lz=3.2
##maximum average height (grid point)
winPoint =  40
rhoH =1.0
rhoL =0.8
rho_inc = (rhoH + rhoL)/2.0

waveLen = 0.4
#input done

mylist = [parentDir,'/',inFile]
delimiter = ''
filepath = delimiter.join(mylist)
#nz enlarged only 
variable = ['PVx','PVy','PVz','PPress', 'Prho']
h5file = h5py.File(filepath,'r+')
#read dataset dimensions
mylist = ['Fields/','Prho','/','002000']
filepath = delimiter.join(mylist)
databk = h5file.get(filepath)
m1 = np.array(databk)
nz=m1.shape[0]
ny=m1.shape[1]
nx=m1.shape[2]

delimiter = ''
dz=dy=dx=Lz/nz
specout = 1000
skip = 20
seq = 0
step = []
for i in range(40):
    step.append(str((i+1)*specout).zfill(6))
#initialize time dependent mixing layer width, KE, PE, enstropy
vorx = np.zeros((nz, ny, nx))
vory = np.zeros((nz, ny, nx))
vorz = np.zeros((nz, ny, nx))
enstropy = np.zeros(len(step))  
sum_x = np.zeros(len(step))  


bub_loc_all = np.zeros(len(step))
bub_loc_all_ori = np.zeros(len(step))
sp_loc_all = np.zeros(len(step))
bub_velo_all = np.zeros(len(step))
bub_velo_all_aver = np.zeros(len(step))
sp_velo_all = np.zeros(len(step))
bub_velo_all_ori = np.zeros(len(step))


ii=0
#calculate bubble and spike location
for istep in step:
    mylist = ['Fields/', variable[4], '/', istep]
    filepath = delimiter.join(mylist)
    databk = h5file.get(filepath)
    rho_data = np.array(databk)
    if nx == 1:
	    m1 = (rho_data[:, ny/2-1, 0] + rho_data[:, ny/2, 0] )/2
    else:
	    m1 = (rho_data[:, ny/2-1, nx/2-1] + rho_data[:, ny/2, nx/2] 
	  + rho_data[:, ny/2-1, nx/2] + rho_data[:, ny/2, nx/2-1])/4.0
    m2 = rho_data[:, 0, 0]
    m1_filter=m1.copy();    
    m2_filter=m2.copy();    
    for jstep in range(2,nz-3):
        m1_filter[jstep]=(m1[jstep-2]+m1[jstep-1]+m1[jstep]+m1[jstep+1]+m1[jstep+2])/5;
        m2_filter[jstep]=(m2[jstep-2]+m2[jstep-1]+m2[jstep]+m2[jstep+1]+m2[jstep+2])/5;

    m1_grad = high_order_gradient(m1_filter,dx,6)
    m2_grad = high_order_gradient(m2_filter,dx,6)

    sp_loc = np.argmax(m1_grad)
    bub_loc = np.argmax(m2_grad)


    sp_loc_all[ii] = sp_loc
    bub_loc_all[ii] = bub_loc

    ii = ii + 1




enstrophy = 0.0
bubRegion = bub_loc_all[i] - winPoint
spkRegion = sp_loc_all[i] + winPoint
for j in range(ny):
    
    
    for k in reversed(range(nz)):
        if rho_data[k, j, 0] == rho_inc:
            bub_inc_loc = k
        break
    
    if bub_inc_loc >= bubRegion:
        ensbub[seq] = ensbub[seq] + np.sum(vorx[bubRegion:bub_inc_loc, j, :]**2
                 + vory[bubRegion:bub_inc_loc, j, :]**2
                 + vorz[bubRegion:bub_inc_loc, j, :]**2)*dx*dy*dz 
              
    for k in range(nz):
        if rho_data[k, j, 0] == rho_inc:
            spk_inc_loc = k
        break
    
    if spk_inc_loc <= spkRegion:
        enspk[seq] = enspk[seq] + np.sum(vorx[spk_inc_loc:spkRegion, j, :]**2
                 + vory[spk_inc_loc:spkRegion, j, :]**2
                 + vorz[spk_inc_loc:spkRegion, j, :]**2)*dx*dy*dz 
    
    
    
    #done

    
    max(bubRegion)
        
    enstrophy[seq] = enstrophy[seq]+ np.sum(vorx[:, ny/2-1, :]**2
                 +vory[:, ny/2-1, :]**2+vorz[:, ny/2-1, :]**2)*dx*dy*dz 
"""
	mylist = ['Fields/',variable[4],'/',istep]
	filepath = delimiter.join(mylist)
	databk = h5file.get(filepath)
	rho = np.array(databk)
	
	
	mylist = ['Fields/',variable[0],'/',istep]
	filepath = delimiter.join(mylist)
	databk = h5file.get(filepath)
	vx = np.array(databk)
	mylist = ['Fields/',variable[1],'/',istep]
	filepath = delimiter.join(mylist)
	databk = h5file.get(filepath)
	vy = np.array(databk)
	mylist = ['Fields/',variable[2],'/',istep]
	filepath = delimiter.join(mylist)
	databk = h5file.get(filepath)
	vz = np.array(databk)
	ke[seq] = 0.5 * np.sum((vx**2 + vy**2 + vz**2) * rho)*dx*dy*dz
	#internal energy
	mylist = ['Fields/',variable[3],'/',istep]
	filepath = delimiter.join(mylist)
	databk = h5file.get(filepath)
	press = np.array(databk)

	#enstrophy
	
	if nx == 1:
		vorx = np.gradient(vz, dz, axis=1) - np.gradient(vy, dz, axis=0)
	else:
		vorx = np.gradient(vz, dz, axis=1) - np.gradient(vy, dz, axis=0)
		vory = np.gradient(vx, dz, axis=0) - np.gradient(vz, dz, axis=2)
		vorz = np.gradient(vy, dz, axis=2) - np.gradient(vx, dz, axis=1)


	enstropy[seq] = 0.5*np.sum(vorx**2+vory**2+vorz**2)*dx*dy*dz 
	
		seq += 1

"""

h5file.close()


#f = open('output.d','w')
#for zz_ref in range(nz):
# f.write("%4s\t%10s\n" % (zz_ref, np.mean(m1[zz_ref,:,:])))
#f.close()
    
    
    
