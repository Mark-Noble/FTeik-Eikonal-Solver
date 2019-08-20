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
! real*8    - zs,xs,ys : Point source coordinates referred expressed in meters
!                    Licit ranges: [0.0,(nz-1.)*dzin][0.0,(nx-1.)*dxin]
!
! integer*4 - eps : radius in number of grid points arround source where then
!                   spherical approximation will be used (for most applications
!                   5 to 10 is enough.
!
! integer*4 - n_sweep : Number of sweeps over model. 2 is in general enough
!________________________________________________________________________
!
  subroutine fteik3d(slow, tt, nz, nx, ny, zsrc, xsrc, ysrc, dz, dx, dy, eps, n_sweep)

  real(kind = 8), parameter :: zerr = 1.d-4
  real(kind = 8), parameter :: Big = 99999.d0
  real(kind = 8) :: eps

    real(kind = 8), dimension(nz,nx,ny), intent(in) :: slow
    real(kind = 8), dimension(nz,nx,ny), intent(out) :: tt
    integer(kind = 4), intent(in) :: nz, nx, ny, n_sweep
    real(kind = 8), intent(in) :: zsrc, xsrc, ysrc, dz, dx, dy

    integer(kind = 4) :: i, j, k, kk, iflag
    integer(kind = 4) :: zsi, xsi, ysi
    integer(kind = 4) :: sgntz, sgntx, sgnty, sgnvz, sgnvx, sgnvy

    real(kind = 8) :: zsa, xsa, ysa
    real(kind = 8) :: tzc, txc, tyc
    real(kind = 8) :: vzero

    ! Check inputs
    if ( nz .lt. 3 .or. nx .lt. 3 .or. ny .lt. 3 ) stop "Error: grid size nz, nx, ny too small"
    if ( dz .le. 0.d0 .or. dx .le. 0.d0 .or. dy .le. 0.d0 ) stop "Error: grid spacing dz, dx, dy too small"
    if ( n_sweep .lt. 0 ) stop "Error: wrong sweep number"
    if ( minval(slow) .le. 0.d0) stop "Error: slownesses are strange"
    if ( zsrc .lt. 0.d0 .or. zsrc .gt. dfloat(nz-1) * dz &
         .or. xsrc .lt. 0.d0 .or. xsrc .gt. dfloat(nx-1) * dx &
         .or. ysrc .lt. 0.d0 .or. ysrc .gt. dfloat(ny-1) * dy ) &
      stop "Error: source out of bounds"

    ! Convert src to grid position and try and take into account machine precision
    zsa = zsrc / dz + 1.d0
    xsa = xsrc / dx + 1.d0
    ysa = ysrc / dy + 1.d0

    ! Try to handle edges simply for source due to precision
    if ( zsa .ge. dfloat(nz) ) zsa = zsa - zerr
    if ( xsa .ge. dfloat(nx) ) xsa = xsa - zerr
    if ( ysa .ge. dfloat(ny) ) ysa = ysa - zerr

    ! Grid points to initialize source
    zsi = int(zsa)
    xsi = int(xsa)
    ysi = int(ysa)
    vzero = slow(zsi,xsi,ysi)

    ! Allocate work array for traveltimes
    tt = Big


    ! Initialize points around source
    tt(zsi,xsi,ysi) = t_ana(zsi, xsi, ysi, dz, dx, dy, zsa, xsa, ysa, vzero)
    tt(zsi+1,xsi,ysi) = t_ana(zsi+1, xsi, ysi, dz, dx, dy, zsa, xsa, ysa, vzero)
    tt(zsi,xsi+1,ysi) = t_ana(zsi, xsi+1, ysi, dz, dx, dy, zsa, xsa, ysa, vzero)
    tt(zsi,xsi,ysi+1) = t_ana(zsi, xsi, ysi+1, dz, dx, dy, zsa, xsa, ysa, vzero)
    tt(zsi+1,xsi+1,ysi) = t_ana(zsi+1, xsi+1, ysi, dz, dx, dy, zsa, xsa, ysa, vzero)
    tt(zsi+1,xsi,ysi+1) = t_ana(zsi+1, xsi, ysi+1, dz, dx, dy, zsa, xsa, ysa, vzero)
    tt(zsi,xsi+1,ysi+1) = t_ana(zsi, xsi+1, ysi+1, dz, dx, dy, zsa, xsa, ysa, vzero)
    tt(zsi+1,xsi+1,ysi+1) = t_ana(zsi+1, xsi+1, ysi+1, dz, dx, dy, zsa, xsa, ysa, vzero)

    ! Full sweeps
      call sweep3dinit(slow, tt, nz, nx, ny, dz, dx, dy, &
                   zsi, xsi, ysi, zsa, xsa, ysa, vzero)
    do kk = 1, n_sweep
      call sweep3d(slow, tt, nz, nx, ny, dz, dx, dy, &
                   zsi, xsi, ysi, zsa, xsa, ysa, vzero)
    end do

    return
  contains

    ! Function to calculate analytical times in homogeneous model
    real(kind = 8) function t_ana(i, j, k, dz, dx, dy, zsa, xsa, ysa, vzero)
      integer(kind = 4), intent(in) :: i, j, k
      real(kind = 8), intent(in) :: dz, dx, dy, zsa, xsa, ysa, vzero

      t_ana = vzero * ( ( ( dfloat(i) - zsa ) * dz )**2.d0 &
                      + ( ( dfloat(j) - xsa ) * dx )**2.d0 &
                      + ( ( dfloat(k) - ysa ) * dy )**2.d0 )**0.5d0
      return
    end function t_ana

    ! Function to calculate analytical times in homogeneous model + derivatives of times
    real(kind = 8) function t_anad(tzc, txc, tyc, i, j, k, dz, dx, dy, zsa, xsa, ysa, vzero)
      integer(kind = 4), intent(in) :: i, j, k
      real(kind = 8), intent(in) :: dz, dx, dy, zsa, xsa, ysa, vzero
      real(kind = 8) :: d0
      real(kind = 8), intent(out) :: tzc, txc, tyc

      d0 = ( ( dfloat(i) - zsa ) * dz )**2.d0 &
           + ( ( dfloat(j) - xsa ) * dx )**2.d0 &
           + ( ( dfloat(k) - ysa ) * dy )**2.d0
      t_anad = vzero * (d0**0.5d0)
      if ( d0 .gt. 0.d0 ) then
        tzc = ( d0**(-0.5d0) ) * ( dfloat(i) - zsa ) * dz * vzero
        txc = ( d0**(-0.5d0) ) * ( dfloat(j) - xsa ) * dx * vzero
        tyc = ( d0**(-0.5d0) ) * ( dfloat(k) - ysa ) * dy * vzero
      else
        tzc = 0.d0
        txc = 0.d0
        tyc = 0.d0
      end if
      return
    end function t_anad

    ! Function to perform full sweep
    subroutine sweep3d(slow, tt, nz, nx, ny, dz, dx, dy, &
                       zsi, xsi, ysi, zsa, xsa, ysa, vzero)
      real(kind = 8), intent(in) :: slow(nz,nx,ny)
      real(kind = 8), intent(inout) :: tt(nz,nx,ny)
      integer(kind = 4), intent(in) :: nz, nx, ny, zsi, xsi, ysi
      real(kind = 8), intent(in) :: dz, dx, dy, zsa, xsa, ysa, vzero
      integer(kind = 4) :: i, j, k, sgntz, sgntx, sgnty, sgnvz, sgnvx, sgnvy
      integer(kind = 4) :: i1, j1, k1, imin
      real(kind = 8) :: dzi, dxi, dyi, dz2i, dx2i, dy2i, dz2dx2, dz2dy2, dx2dy2, dsum
      real(kind = 8) :: vref, time_sol(8)
      real(kind = 8) :: tv, te, tn, ten, tnv, tev, tnve
      real(kind = 8) :: t1d, t2d, t3d, t1, t2, t3, ta, tb, tc
      real(kind = 8) :: t1d1, t1d2, t1d3, t2d1, t2d2, t2d3
      real(kind = 8) :: apoly, bpoly, cpoly, dpoly

      ! Precalculate constants
      dzi = 1.d0 / dz
      dxi = 1.d0 / dx
      dyi = 1.d0 / dy
      dz2i = 1.d0 / (dz*dz)
      dx2i = 1.d0 / (dx*dx)
      dy2i = 1.d0 / (dy*dy)
      dz2dx2 = dz2i * dx2i
      dz2dy2 = dz2i * dy2i
      dx2dy2 = dx2i * dy2i
      dsum = dz2i + dx2i + dy2i

        ! First sweeping: Top->Bottom ; West->East ; South->North
        sgntz = 1 ; sgntx = 1 ; sgnty = 1
        sgnvz = 1 ; sgnvx = 1 ; sgnvy = 1
        do k = 2, ny
          do j = 2, nx
            do i = 2, nz
              include "Include_FTeik3d.f"
            end do
          end do
        end do

        ! Second sweeping: Top->Bottom ; East->West ; South->North
        sgntz = 1 ; sgntx = -1 ; sgnty = 1
        sgnvz = 1 ; sgnvx = 0 ; sgnvy = 1
        do k = 2, ny
          do j = nx-1, 1, -1
            do i = 2, nz
              include "Include_FTeik3d.f"
            end do
          end do
        end do

        ! Third sweeping: Top->Bottom ; West->East ; North->South
        sgntz = 1 ; sgntx = 1 ; sgnty = -1
        sgnvz = 1 ; sgnvx = 1 ; sgnvy = 0
        do k = ny-1, 1, -1
          do j = 2, nx
            do i = 2, nz
              include "Include_FTeik3d.f"
            end do
          end do
        end do

        ! Fouth sweeping: Top->Bottom ; East->West ; North->South
        sgntz = 1 ; sgntx = -1 ; sgnty = -1
        sgnvz = 1 ; sgnvx = 0 ; sgnvy = 0
        do k = ny-1, 1, -1
          do j = nx-1, 1, -1
            do i = 2, nz
              include "Include_FTeik3d.f"
            end do
          end do
        end do

        ! Fifth sweeping: Bottom->Top ; West->East ; South->North
        sgntz = -1 ; sgntx = 1 ; sgnty = 1
        sgnvz = 0 ; sgnvx = 1 ; sgnvy = 1
        do k = 2, ny
          do j = 2, nx
            do i = nz-1, 1, -1
              include "Include_FTeik3d.f"
            end do
          end do
        end do

        ! Sixth sweeping: Bottom->Top ; East->West ; South->North
        sgntz = -1 ; sgntx = -1 ; sgnty = 1
        sgnvz = 0 ; sgnvx = 0 ; sgnvy = 1
        do k = 2, ny
          do j = nx-1, 1, -1
            do i = nz-1, 1, -1
              include "Include_FTeik3d.f"
            end do
          end do
        end do

        ! Seventh sweeping: Bottom->Top ; West->East ; North->South
        sgntz = -1 ; sgntx = 1 ; sgnty = -1
        sgnvz = 0 ; sgnvx = 1 ; sgnvy = 0
        do k = ny-1, 1, -1
          do j = 2, nx
            do i = nz-1, 1, -1
              include "Include_FTeik3d.f"
            end do
          end do
        end do

        ! Eighth sweeping: Bottom->Top ; East->West ; North->South
        sgntz = -1 ; sgntx = -1 ; sgnty = -1
        sgnvz = 0 ; sgnvx = 0 ; sgnvy = 0
        do k = ny-1, 1, -1
          do j = nx-1, 1, -1
            do i = nz-1, 1, -1
              include "Include_FTeik3d.f"
            end do
          end do
        end do

      return
    end subroutine sweep3d

    ! Function to perform full sweep
    subroutine sweep3dinit(slow, tt, nz, nx, ny, dz, dx, dy, &
                       zsi, xsi, ysi, zsa, xsa, ysa, vzero)
      real(kind = 8), intent(in) :: slow(nz,nx,ny)
      real(kind = 8), intent(inout) :: tt(nz,nx,ny)
      integer(kind = 4), intent(in) :: nz, nx, ny, zsi, xsi, ysi
      real(kind = 8), intent(in) :: dz, dx, dy, zsa, xsa, ysa, vzero
      integer(kind = 4) :: i, j, k, sgntz, sgntx, sgnty, sgnvz, sgnvx, sgnvy
      integer(kind = 4) :: i1, j1, k1, imin
      real(kind = 8) :: dzi, dxi, dyi, dz2i, dx2i, dy2i, dz2dx2, dz2dy2, dx2dy2, dsum
      real(kind = 8) :: vref, time_sol(8)
      real(kind = 8) :: tv, te, tn, ten, tnv, tev, tnve
      real(kind = 8) :: t1d, t2d, t3d, t1, t2, t3, ta, tb, tc
      real(kind = 8) :: t1d1, t1d2, t1d3, t2d1, t2d2, t2d3
      real(kind = 8) :: apoly, bpoly, cpoly, dpoly

      ! Precalculate constants
      dzi = 1.d0 / dz
      dxi = 1.d0 / dx
      dyi = 1.d0 / dy
      dz2i = 1.d0 / (dz*dz)
      dx2i = 1.d0 / (dx*dx)
      dy2i = 1.d0 / (dy*dy)
      dz2dx2 = dz2i * dx2i
      dz2dy2 = dz2i * dy2i
      dx2dy2 = dx2i * dy2i
      dsum = dz2i + dx2i + dy2i

        ! First sweeping: Top->Bottom ; West->East ; South->North
        sgntz = 1 ; sgntx = 1 ; sgnty = 1
        sgnvz = 1 ; sgnvx = 1 ; sgnvy = 1
        do k = max(2,ysi),ny
          do j = max(2,xsi),nx
            do i = max(2,zsi),nz
              include "Include_FTeik3d.f"
            end do
          end do
        end do

        ! Second sweeping: Top->Bottom ; East->West ; South->North
        sgntz = 1 ; sgntx = -1 ; sgnty = 1
        sgnvz = 1 ; sgnvx = 0 ; sgnvy = 1
        do k = max(2,ysi),ny
          do j = xsi+1,1,-1
            do i = max(2,zsi),nz
              include "Include_FTeik3d.f"
            end do
          end do
        end do

        ! Third sweeping: Top->Bottom ; West->East ; North->South
        sgntz = 1 ; sgntx = 1 ; sgnty = -1
        sgnvz = 1 ; sgnvx = 1 ; sgnvy = 0
        do k = ysi+1,1,-1
          do j = max(2,xsi),nx
            do i = max(2,zsi),nz
              include "Include_FTeik3d.f"
            end do
          end do
        end do

        ! Fouth sweeping: Top->Bottom ; East->West ; North->South
        sgntz = 1 ; sgntx = -1 ; sgnty = -1
        sgnvz = 1 ; sgnvx = 0 ; sgnvy = 0
        do k = ysi+1,1,-1
          do j = xsi+1,1,-1
            do i = max(2,zsi),nz
              include "Include_FTeik3d.f"
            end do
          end do
        end do

        ! Fifth sweeping: Bottom->Top ; West->East ; South->North
        sgntz = -1 ; sgntx = 1 ; sgnty = 1
        sgnvz = 0 ; sgnvx = 1 ; sgnvy = 1
        do k = max(2,ysi),ny
          do j = max(2,xsi),nx
            do i = zsi+1,1,-1
              include "Include_FTeik3d.f"
            end do
          end do
        end do

        ! Sixth sweeping: Bottom->Top ; East->West ; South->North
        sgntz = -1 ; sgntx = -1 ; sgnty = 1
        sgnvz = 0 ; sgnvx = 0 ; sgnvy = 1
        do k = max(2,ysi),ny
          do j = xsi+1,1,-1
            do i = zsi+1,1,-1
              include "Include_FTeik3d.f"
            end do
          end do
        end do

        ! Seventh sweeping: Bottom->Top ; West->East ; North->South
        sgntz = -1 ; sgntx = 1 ; sgnty = -1
        sgnvz = 0 ; sgnvx = 1 ; sgnvy = 0


       do k = ysi+1,1,-1
         do j = max(2,xsi),nx
           do i = zsi+1,1,-1
              include "Include_FTeik3d.f"
            end do
          end do
        end do

        ! Eighth sweeping: Bottom->Top ; East->West ; North->South
        sgntz = -1 ; sgntx = -1 ; sgnty = -1
        sgnvz = 0 ; sgnvx = 0 ; sgnvy = 0
        do k = ysi+1,1,-1
          do j = xsi+1,1,-1
            do i = zsi+1,1,-1
              include "Include_FTeik3d.f"
            end do
          end do
        end do

      return
    end subroutine sweep3dinit

  end subroutine fteik3d

