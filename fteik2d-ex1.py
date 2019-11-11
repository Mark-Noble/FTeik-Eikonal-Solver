#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 1 13:10:39 2019

@author: mnoble
"""
import numpy as np
import eik2d
import time

# set up a few  constants and size of model
eps=10 ; n_sweep=4
nz=300 ; nx=300
dz=2. ; dx=2.
zsrc=0. ; xsrc=0.

# create velocity - homogeneous = 2000m/s and add high and low velocity anomaly
vel = np.zeros((nz,nx))
vel[:,:]=2000. ; vel[100:150,50:150]=1000. ; vel[100:150,150:250]=4900.
# convert to lownes
slow=1./vel

# call 2d Eokonal and print elapsed time
start=time.time()
tt = eik2d.fteik2d(slow,zsrc,xsrc,dz,dx,eps,n_sweep,nz,nx)
end=time.time()
print(end-start)

# write velocity and traveltimes to disk ( binary format)
def writebin(inp,flnam):
    # Write binary fila on disk (32 bits)
    with open(flnam,"wb") as fl:
        inp.T.astype('float32').tofile(fl)
        
writebin(vel,'vel.bin')
writebin(tt,'time.bin')

# That's all