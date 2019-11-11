# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

# Convert cm to inchs
def cm2inch(value):
    return value/2.54

# Default font size
plt.rcParams['font.size'] = 12

# Read in binary velocity model and traveltime
vfile = "./vel.bin"
tfile = "./time.bin"

# dimension 3D
nz, nx = 300, 300
dz , dx = 2., 2.
# origin of model
fz, fx = 0, 0

dtype = "float32"
# read in binaries (3D)
vel = np.fromfile(vfile, dtype = dtype).reshape((nz,nx), order = "F")
tt = np.fromfile(tfile, dtype = dtype).reshape((nz,nx), order = "F")

# Define range
z = np.arange(0+fz, dz*nz+fz, dz)
x = np.arange(0+fx, dx*nx+fx, dx)
# Define color range and contour lines
cmi, cma, dc = 1200, 4800, 50
clevels=np.arange(cmi,cma,((cma-cmi)/dc) )
llevels=np.arange(0.,100,0.01)
# Set up Plot alaong X axis


#### PLOT1 ###################################################################
# Get a 2D vertical slice for a constant Y
fig = plt.figure(figsize=(cm2inch(20),cm2inch(20)))
left,width=0.1,0.8
bottom,height=0.1,0.8
rect = [left, bottom, width, height]
ax1 = plt.axes(rect)
X, Z = np.meshgrid(x,z)
# Color plot, contour, and label contours
cc=ax1.contourf(X,Z,vel,clevels, cmap='rainbow',extend='both')
cl1=ax1.contour(X,Z,tt,llevels, colors='k',linewidths=(1,0.5,0.5,0.5,0.5))
plt.clabel(cl1,cl1.levels[::5], fontsize=4, colors='k', fmt='%4.2f', inline_spacing=1)
# axis
plt.ylim(np.max(z),np.min(z)) # flip axis for depth
plt.xlabel('X (m)')
plt.ylabel('Depth (m)')
plt.title('Vp model - P traveltime')
#Color bar
cb=plt.colorbar(cc,shrink=0.75, aspect=10, pad=0.03)
cb.ax.invert_yaxis()
cb.ax.set_title('m/s')
#fig.tight_layout()
plt.show()
#fig.savefig('test1.png', dpi=600)
fig.savefig('fig1.svg')


















