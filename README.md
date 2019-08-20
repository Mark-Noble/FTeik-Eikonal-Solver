 
    FTeik PACKAGE VERSION 1.0
    Copyright (c) 2019 Mark NOBLE, MINES ParisTech, France
    Email: mark.noble@mines-paristech.fr
    VERSION 1.0: 2019-08-01 , First public release 
 
    The FTeik package is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
    
    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
    
    See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program; if not, write to: Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

# FTeik-Eikonal-Solver

FTeik 2D and 3D Eikonal solver to compute first arrival traveltimes in a heterogeneous isotropic velocity model, with the possibility to use different grid spacing in all directions.

## Reference paper
Detailed implementation of local operator and global propagation scheme implemented in these subroutines come form the paper of : M. Noble, A. Gesret and N. Belayouni, 2014, Accurate 3-D finite difference computation of traveltimes in strongly heterogeneous media, Geophys.J.Int.,199,(3),1572-158.

**If you use this program for an academic project or a commercial project, citing this paper would be appreciated.**

## Compiling and running the code

The Eikonal solver subroutine includes 2 files:
- In 2D fteik2d.f90 and Include_FTeik2d.f

The subroutines are written in Fortran 90

The package comes along with a simple example of a main program to show how to call the Eikonal solver
- In 2D mainFTeik2d.f90

With the Gnu compiler "gfortran", it is recommended to use the options. Also included an example of a Makefile
  -O3 -ffree-form
  
  **Compile and execute code**
  
 ```console
 make all
 ./mainFTeik2d.exe
 ./mainFTeik3d.exe
 ```
 The program writes on disk the traveltimes in the file called ttmap.

## Arguments required to call the Eikonal Solver subroutine

 - In 2D: call fteik2d(slow, tt, nz, nx, zsrc, xsrc, dz, dx, eps, n_sweep)
 
 - In 3D: call fteik3d(slow, tt, nz, nx, ny, zsrc, xsrc, ysrc, dz, dx, dy, eps, n_sweep)
 
 **NOTE 1**: TravelTime field array and slowness field array have the same
           dimension. Slownesses are defined at center of cell, whereas times
           are computed on the corners. In practice the last row and last
           column of slowness field are not used in the computation. This is
           the same as in Podvin and Lecomte algorithm.

 **NOTE 2**: In order to get accurate traveltimes, all real numbers (scalars and arrays)
           must be decalred in double precision.

- integer*4 - nz,nx,ny : Dimensions of the time field array tt
                      in 2D tt(nz,nx) or in 3D tt(nz,nx,ny)
                      No dimension may be lower than 3.

 - real*8    - dz,dx,dy : Mesh spacing along the 3 axis

 - real*8    - tt       : Travel time field array: tt((nz,nx) or tt(nz,nx,ny)

 - real*8    - slow     : Slowness field array: slow(nz,nx) or slow(nz,nx,ny)

 - real*8    - zs,xs, : Point source coordinates referred expressed in meters
                    Licit ranges: [0.0,(nz-1.)*dzin][0.0,(nx-1.)*dxin]

 - integer*4 - epsin : radius in number of grid points arround source where then
                   spherical approximation will be used (for most applications
                   5 to 10 is enough.

 - integer*4 - nsweep : Number of sweeps over model. 1 is in general enough
