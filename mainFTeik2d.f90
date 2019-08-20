! 
! 
!                      FTeik PACKAGE VERSION 1.0
!                      ---------------------------
!  Copyright (c) 2019 Mark NOBLE, MINES ParisTech, France
!  Email: mark.noble@mines-paristech.fr
! 
!  2D and 3D Eikonal solver to compute first arrival traveltimes in a heterogeneous
!  isotropic velocity model, with the possibility to use different grid spacing in all
!  directions.
! 
!  RELEASE
!  -------
!  VERSION 1.0: 2019-08-01 , First public release 
! 
!  LEGAL STATEMENT
!  ---------------
!  Copyright (c) 2019 Mark NOBLE, MINES ParisTech, France
! 
!  The FTeik package is free software; you can redistribute it and/or modify it
!  under the terms of the GNU General Public License as published by the Free Software
!  Foundation; either version 2 of the License, or (at your option) any later version.
!  
!  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!  See the GNU General Public License for more details.
! 
!  You should have received a copy of the GNU General Public License along with this program;
!  if not, write to:
!  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
! 
!  BIBLIOGRAPHY
!  ------------
!  Detailed implementation of local operators and global propagation scheme implemented in
!  this sunroutine are inspired form this paper:
!  M. Noble,  A. Gesret and N. Belayouni, 2014, Accurate 3-D finite difference computation of
!  traveltimes in strongly heterogeneous media Geophys. J. Int.,  199 (3): 1572-158.
! 
!  If you find this algorithm useful, citing this paper would be apreciated
!______________________________________________________________________________
!            Arguments required to call the Eikonal Solver subroutine
!
! in 2D: call fteik2d(slow, tt, nz, nx, zsrc, xsrc, dz, dx, eps, n_sweep)
! in 3D: call fteik3d(slow, tt, nz, nx, ny, zsrc, xsrc, ysrc, dz, dx, dy, eps, n_sweep)
! 
! WARNING : TravelTime field array and slowness field array have the same
!           dimension. Slownesses are defined at center of cell, whereas times
!           are computed on the corners. In practice the last row and last
!           column of slowness field are not used in the computation. This is
!           the same as in Podvin and Lecomte algorithm.
!
! WARNING : In order to get accurate traveltimes, all reals (scalars and arays)
!           mist be decalred in double precision.
!
! integer*4 - nz,nx,ny : Dimensions of the time field array tt
!                      in 2D tt(nz,nx) or in 3D tt(nz,nx,ny)
!                      No dimension may be lower than 3.
!
! real*8    - dz,dx,dy : Mesh spacing along the 3 axis
!
! real*8    - tt       : Travel time field array: tt((nz,nx) or tt(nz,nx,ny)
!
! real*8    - slow     : Slowness field array: slow(nz,nx) or slow(nz,nx,ny)
!
! real*8    - zs,xs, : Point source coordinates referred expressed in meters
!                    Licit ranges: [0.0,(nz-1.)*dzin][0.0,(nx-1.)*dxin]
!
! integer*4 - eps : radius in number of grid points arround source where then
!                   spherical approximation will be used (for most applications
!                   5 to 10 is enough.
!
! integer*4 - n_sweep : Number of sweeps over model. 2 is in general enough
!________________________________________________________________________
!
! main program to call new Eikonal

program mainFTeik2d

implicit none

! Slowness model has same size as traveltime maps, last row and column not used
! slow(nz,nx)
real(8),allocatable,dimension(:,:) :: slow

! First arrival times(nz,nx)
real(8),allocatable,dimension(:,:) :: tt

! Number of grid points of traveltime array and velocity model
integer(4) :: nz,nx

! Number of points around source where pertubation operators used
! In general 10 is good enough
integer(4) :: eps

! Grid spacing in Z and X; and source position - ALL in meters
real(8) :: dz,dx,xs,zs

! Number of global sweeps
integer(4) :: n_sweep

! Handle error, when allocating arrays
integer(4) :: ierr

! In this example we are creating a simple 2 layer model with a
! slow perturbation in the middle
nz=1000 ; nx=1000
dz=10.d0 ; dx=10.d0
zs=150.d0 ; xs=150.d0

! set number of sweeps
n_sweep=1

! number of points to use spherical approximation (usually between 5 and 10)
eps=10

! Allocate velocity model
allocate(slow(nz,nx),stat=ierr)
if (ierr .ne. 0) stop "Error in allocating slow array"

! Allocate  first arrival times
allocate(tt(nz,nx),stat=ierr)
if (ierr .ne. 0) stop "Error in allocating tt array"

! Create a 2 layer 1D model,
slow(1:100,:)=1.d0 / 2000.d0
slow(101:nz,:)=1.d0 / 4000.d0
! insert velocity anomaly
slow(400:600,400:600)=1.d0 / 1000.d0

! call 2D eikonal solver
 call fteik2d(slow,tt,nz,nx,zs,xs,dz,dx,eps,n_sweep)

!write results in file in binary format
open(10,file='ttmap',access='direct',recl=4*(nx)*(nz))
write(10,rec=1) sngl(tt)
close(10)

! That's all
deallocate(tt)
deallocate(slow)

call exit
end program mainFTeik2d
