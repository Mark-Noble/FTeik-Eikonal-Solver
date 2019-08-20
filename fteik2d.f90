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
! real*8    - zs,xs,ys : Point source coordinates referred expressed in meters
!                    Licit ranges: [0.0,(nz-1.)*dzin][0.0,(nx-1.)*dxin]
!
! integer*4 - eps : radius in number of grid points arround source where then
!                   spherical approximation will be used (for most applications
!                   5 to 10 is enough.
!
! integer*4 - n_sweep : Number of sweeps over model. 1 is in general enough
!________________________________________________________________________
!
subroutine fteik2d(slow, tt, nz, nx, zsrc, xsrc, dz, dx, eps, n_sweep)

  implicit none
  
  integer(kind = 4), parameter :: nmax = 15000
  real(kind = 8), parameter :: zerr = 1.d-4
  real(kind = 8), parameter :: Big = 99999.d0
  
    real(kind = 8), dimension(nz,nx), intent(in) :: slow
    real(kind = 8), dimension(nz,nx), intent(out) :: tt
    integer(kind = 4), intent(in) :: nz, nx, n_sweep, eps
    real(kind = 8), intent(in) :: zsrc, xsrc, dz, dx

    integer(kind = 4) :: imin, i, j, kk, i1, j1, iflag
    integer(kind = 4) :: zsi, xsi
    integer(kind = 4) :: sgntz, sgntx, sgnvz, sgnvx

    real(kind = 8) :: zsa, xsa
    real(kind = 8) :: tzc, txc
    real(kind = 8) :: vzero, vref
    real(kind = 8) :: dzu, dzd, dxw, dxe
    real(kind = 8) :: tv, te, tev
    real(kind = 8) :: t0c, td(nmax), ta, tb, tauv, taue, tauev
    real(kind = 8) :: dzi, dxi, dz2i, dx2i
    real(kind = 8) :: sgnrz, sgnrx
    real(kind = 8) :: apoly, bpoly, cpoly, dpoly

    ! Check inputs
    if ( nz .lt. 3 .or. nx .lt. 3 ) stop "Error: grid size nz, nx too small"
    if ( max(nz, nx) .gt. nmax ) stop "Error: must increase size of NMAX"
    if ( dz .le. 0.d0 .or. dx .le. 0.d0 ) stop "Error: grid spacing dz, dx too small"
    if ( n_sweep .lt. 1 ) stop "Error: wrong sweep number"
    if ( minval(slow) .le. 0.d0 .or. maxval(slow) .ge. 1.d0 ) stop "Error: slownesses are strange"
    if ( zsrc .lt. 0.d0 .or. zsrc .gt. dfloat(nz-1) * dz &
         .or. xsrc .lt. 0.d0 .or. xsrc .gt. dfloat(nx-1) * dx ) &
      stop "Error: source out of bounds"

    ! Convert src to grid position and try and take into account machine precision
    zsa = zsrc / dz + 1.d0
    xsa = xsrc / dx + 1.d0

    ! Try to handle edges simply for source due to precision
    if ( zsa .ge. dfloat(nz) ) zsa = zsa - zerr
    if ( xsa .ge. dfloat(nx) ) xsa = xsa - zerr

    ! Grid points to initialize source
    zsi = int(zsa)
    xsi = int(xsa)
    vzero = slow(zsi,xsi)

    ! Allocate work array for traveltimes
    tt = Big

    ! Do our best to initialize source
    dzu = dabs( zsa - dfloat(zsi) )
    dzd = dabs( dfloat(zsi+1) - zsa )
    dxw = dabs( xsa - dfloat(xsi) )
    dxe = dabs( dfloat(xsi+1) - xsa )
    iflag = 0

    ! Source seems close enough to a grid point in X and Y direction
    if ( min(dzu, dzd) .lt. zerr .and. min(dxw, dxe) .lt. zerr ) then
      zsa = dnint(zsa)
      xsa = dnint(xsa)
      iflag = 1
    end if

    ! At least one of coordinates not close to any grid point in X and Y direction
    if ( min(dzu, dzd) .gt. zerr .or. min(dxw, dxe) .gt. zerr ) then
      if ( min(dzu, dzd) .lt. zerr) zsa = dnint(zsa)
      if ( min(dxw, dxe) .lt. zerr) xsa = dnint(xsa)
      iflag = 2
    end if

    ! Oops we are lost, not sure this happens - fix src to nearest grid point
    if ( iflag .ne. 1 .and. iflag .ne. 2 ) then
      zsa = dnint(zsa)
      xsa = dnint(xsa)
      iflag = 3
    end if

    ! We know where src is - start first propagation
    select case(iflag)
    case(1)
      tt(nint(zsa),nint(xsa)) = 0.d0
    case(3)
      tt(nint(zsa),nint(xsa)) = 0.d0
    case(2)
      dzu = dabs( zsa - dfloat(zsi) )
      dzd = dabs( dfloat(zsi+1) - zsa )
      dxw = dabs( xsa - dfloat(xsi) )
      dxe = dabs( dfloat(xsi+1) - xsa )

      ! First initialize 4 points around source
      tt(zsi,xsi) = t_ana(zsi, xsi, dz, dx, zsa, xsa, vzero)
      tt(zsi+1,xsi) = t_ana(zsi+1, xsi, dz, dx, zsa, xsa, vzero)
      tt(zsi,xsi+1) = t_ana(zsi, xsi+1, dz, dx, zsa, xsa, vzero)
      tt(zsi+1,xsi+1) = t_ana(zsi+1, xsi+1, dz, dx, zsa, xsa, vzero)

      td = Big
      td(xsi+1) = vzero * dxe * dx
      dx2i = 1.d0 / (dx*dx)
      do j = xsi+2, nx
        vref = slow(zsi,j-1)
        td(j) = td(j-1) + dx * vref
        tauv = td(j) - vzero * dabs( dfloat(j) - xsa ) * dx
        tauev = td(j-1) - vzero * dabs( dfloat(j-1) - xsa ) * dx

        sgntz = 1
        sgntx = 1
        taue = tt(zsi+1,j-1) - t_ana(zsi+1, j-1, dz, dx, zsa, xsa, vzero)
        t0c = t_anad(tzc, txc, zsi+1, j, dz, dx, zsa, xsa, vzero)
        sgnrz = dfloat(sgntz)
        sgnrx = dfloat(sgntx)
        ta = tauev + taue - tauv
        tb = tauev - taue + tauv
        dz2i = 1.d0 / (dzd*dzd) * dz
        apoly = dz2i + dx2i
        bpoly = 4.d0 * ( sgnrx * txc / dx + sgnrz * tzc / dzd ) - 2.d0 * ( ta * dx2i + tb * dz2i )
        cpoly = ( ta*ta * dx2i ) + ( tb*tb * dz2i ) &
                - 4.d0 * ( sgnrx * txc / dx * ta + sgnrz * tzc / dzd * tb ) &
                + 4.d0 * ( vzero*vzero - vref*vref )
        dpoly = bpoly*bpoly - 4.d0 * apoly * cpoly
        if ( dpoly .ge. 0.d0 ) tt(zsi+1,j) = 0.5d0 * ( sqrt(dpoly) - bpoly ) / apoly + t0c

        sgntz = -1
        sgntx = 1
        taue = tt(zsi,j-1) - t_ana(zsi, j-1, dz, dx, zsa, xsa, vzero)
        t0c = t_anad(tzc, txc, zsi, j, dz, dx, zsa, xsa, vzero)
        sgnrz = dfloat(sgntz)
        sgnrx = dfloat(sgntx)
        ta = tauev + taue - tauv
        tb = tauev - taue + tauv
        dz2i = 1.d0 / (dzu*dzu) * dz
        apoly = dz2i + dx2i
        bpoly = 4.d0 * ( sgnrx * txc / dx + sgnrz * tzc / dzu ) - 2.d0 * ( ta * dx2i + tb * dz2i )
        cpoly = ( ta*ta * dx2i ) + ( tb*tb * dz2i ) &
                - 4.d0 * ( sgnrx * txc / dx * ta + sgnrz * tzc / dzu * tb ) &
                + 4.d0 * ( vzero*vzero - vref*vref )
        dpoly = bpoly*bpoly - 4.d0 * apoly * cpoly
        if ( dpoly .ge. 0.d0 ) tt(zsi,j) = 0.5d0 * ( sqrt(dpoly) - bpoly ) / apoly + t0c
      end do

      td(xsi) = vzero * dxw * dx
      do j = xsi-1, 1, -1
        vref = slow(zsi,j)
        td(j) = td(j+1) + dx * vref
        tauv = td(j) - vzero * dabs( dfloat(j) - xsa ) * dx
        tauev = td(j+1) - vzero * dabs( dfloat(j+1) - xsa ) * dx

        sgntz = 1
        sgntx = -1
        taue = tt(zsi+1,j+1) - t_ana(zsi+1, j+1, dz, dx, zsa, xsa, vzero)
        t0c = t_anad(tzc, txc, zsi+1, j, dz, dx, zsa, xsa, vzero)
        sgnrz = dfloat(sgntz)
        sgnrx = dfloat(sgntx)
        ta = tauev + taue - tauv
        tb = tauev - taue + tauv
        dz2i = 1.d0 / (dzd*dzd) * dz
        apoly = dz2i + dx2i
        bpoly = 4.d0 * ( sgnrx * txc / dx + sgnrz * tzc / dzd ) - 2.d0 * ( ta * dx2i + tb * dz2i )
        cpoly = ( ta*ta * dx2i ) + ( tb*tb * dz2i ) &
                - 4.d0 * ( sgnrx * txc / dx * ta + sgnrz * tzc / dzd * tb ) &
                + 4.d0 * ( vzero*vzero - vref*vref )
        dpoly = bpoly*bpoly - 4.d0 * apoly * cpoly
        if ( dpoly .ge. 0.d0 ) tt(zsi+1,j) = 0.5d0 * ( sqrt(dpoly) - bpoly ) / apoly + t0c

        sgntz = -1
        sgntx = -1
        taue = tt(zsi+1,j+1) - t_ana(zsi+1, j+1, dz, dx, zsa, xsa, vzero)
        t0c = t_anad(tzc, txc, zsi, j, dz, dx, zsa, xsa, vzero)
        sgnrz = dfloat(sgntz)
        sgnrx = dfloat(sgntx)
        ta = tauev + taue - tauv
        tb = tauev - taue + tauv
        dz2i = 1.d0 / (dzu*dzu) * dz
        apoly = dz2i + dx2i
        bpoly = 4.d0 * ( sgnrx * txc / dx + sgnrz * tzc / dzu ) - 2.d0 * ( ta * dx2i + tb * dz2i )
        cpoly = ( ta*ta * dx2i ) + ( tb*tb * dz2i ) &
                - 4.d0 * ( sgnrx * txc / dx * ta + sgnrz * tzc / dzu * tb ) &
                + 4.d0 * ( vzero*vzero - vref*vref )
        dpoly = bpoly*bpoly - 4.d0 * apoly * cpoly
        if ( dpoly .ge. 0.d0 ) tt(zsi,j) = 0.5d0 * ( sqrt(dpoly) - bpoly ) / apoly + t0c
      end do

      td = Big
      td(zsi+1) = vzero * dzd * dz
      dz2i = 1.d0 / (dz*dz)
      do i = zsi+2, nz
        vref = slow(i-1,xsi)
        td(i) = td(i-1) + dz * vref
        taue = td(i) - vzero * dabs(dfloat(i) - zsa) * dz
        tauev = td(i-1) - vzero * dabs(dfloat(i-1) - zsa) * dz

        sgntz = 1
        sgntx = 1
        tauv = tt(i-1,xsi+1) - t_ana(i-1, xsi+1, dz, dx, zsa, xsa, vzero)
        t0c = t_anad(tzc, txc, i, xsi+1, dz, dx, zsa, xsa, vzero)
        sgnrz = dfloat(sgntz)
        sgnrx = dfloat(sgntx)
        ta = tauev + taue - tauv
        tb = tauev - taue + tauv
        dx2i = 1.d0 / (dxe*dxe) * dx
        apoly = dz2i + dx2i
        bpoly = 4.d0 * ( sgnrx * txc / dxe + sgnrz * tzc / dz ) - 2.d0 * ( ta * dx2i + tb * dz2i )
        cpoly = ( ta*ta * dx2i ) + ( tb*tb * dz2i ) &
                - 4.d0 * ( sgnrx * txc / dxe * ta + sgnrz * tzc / dz * tb ) &
                + 4.d0 * ( vzero*vzero - vref*vref )
        dpoly = bpoly*bpoly - 4.d0 * apoly * cpoly
        if ( dpoly .ge. 0.d0 ) tt(i,xsi+1) = 0.5d0 * ( sqrt(dpoly) - bpoly ) / apoly + t0c

        sgntz = 1
        sgntx = -1
        tauv = tt(i-1,xsi) - t_ana(i-1, xsi, dz, dx, zsa, xsa, vzero)
        t0c = t_anad(tzc, txc, i, xsi, dz, dx, zsa, xsa, vzero)
        sgnrz = dfloat(sgntz)
        sgnrx = dfloat(sgntx)
        ta = tauev + taue - tauv
        tb = tauev - taue + tauv
        dx2i = 1.d0 / (dxw*dxw) * dx
        apoly = dz2i + dx2i
        bpoly = 4.d0 * ( sgnrx * txc / dxw + sgnrz * tzc / dz ) - 2.d0 * ( ta * dx2i + tb * dz2i )
        cpoly = ( ta*ta * dx2i ) + ( tb*tb * dz2i ) &
                - 4.d0 * ( sgnrx * txc / dxw * ta + sgnrz * tzc / dz * tb ) &
                + 4.d0 * ( vzero*vzero - vref*vref )
        dpoly = bpoly*bpoly - 4.d0 * apoly * cpoly
        if ( dpoly .ge. 0.d0 ) tt(i,xsi) = 0.5d0 * ( sqrt(dpoly) - bpoly ) / apoly + t0c
      end do

      td(zsi) = vzero * dzu * dz
      do i = zsi-1,1,-1
        vref = slow(i,xsi)
        td(i) = td(i+1) + dz * vref
        taue = td(i) - vzero * dabs( dfloat(i) - zsa ) * dz
        tauev = td(i+1) - vzero * dabs( dfloat(i+1) - zsa) * dz

        sgntz = -1
        sgntx = 1
        tauv = tt(i+1,xsi+1) - t_ana(i+1, xsi+1, dz, dx, zsa, xsa, vzero)
        t0c = t_anad(tzc, txc, i, xsi+1, dz, dx, zsa, xsa, vzero)
        sgnrz = dfloat(sgntz)
        sgnrx = dfloat(sgntx)
        ta = tauev + taue - tauv
        tb = tauev - taue + tauv
        dx2i = 1.d0 / (dxe*dxe) * dx
        apoly = dz2i + dx2i
        bpoly = 4.d0 * ( sgnrx * txc / dxe + sgnrz * tzc / dz ) - 2.d0 * ( ta * dx2i + tb * dz2i )
        cpoly = ( ta*ta * dx2i ) + ( tb*tb * dz2i ) &
                - 4.d0 * ( sgnrx * txc / dxe * ta + sgnrz * tzc / dz * tb ) &
                + 4.d0 * ( vzero*vzero - vref*vref )
        dpoly = bpoly*bpoly - 4.d0 * apoly * cpoly
        if ( dpoly .ge. 0.d0 ) tt(i,xsi+1) = 0.5d0 * ( sqrt(dpoly) - bpoly ) / apoly + t0c

        sgntz = -1
        sgntx = -1
        tauv = tt(i+1,xsi) - t_ana(i+1, xsi, dz, dx, zsa, xsa, vzero)
        t0c = t_anad(tzc, txc, i, xsi, dz, dx, zsa, xsa, vzero)
        sgnrz = dfloat(sgntz)
        sgnrx = dfloat(sgntx)
        ta = tauev + taue - tauv
        tb = tauev - taue + tauv
        dx2i = 1.d0 / (dxw*dxw) * dx
        apoly = dz2i + dx2i
        bpoly = 4.d0 * ( sgnrx * txc / dxw + sgnrz * tzc / dz ) - 2.d0 * ( ta * dx2i + tb * dz2i )
        cpoly = ( ta*ta * dx2i ) + ( tb*tb * dz2i ) &
                - 4.d0 * ( sgnrx * txc / dxw * ta + sgnrz * tzc / dz * tb ) &
                + 4.d0 * ( vzero*vzero - vref*vref )
        dpoly = bpoly*bpoly - 4.d0 * apoly * cpoly
        if ( dpoly .ge. 0.d0 ) tt(i,xsi) = 0.5d0 * ( sqrt(dpoly) - bpoly ) / apoly + t0c
      end do
    end select

    ! Full sweeps
    do kk = 1, n_sweep
      call sweep2d(slow, tt, nz, nx, dz, dx, &
                   zsi, xsi, zsa, xsa, vzero, eps)
    end do

    return
    
  contains

    ! Function to calculate analytical times in homogeneous model
    real(kind = 8) function t_ana(i, j, dz, dx, zsa, xsa, vzero)
    implicit none
      integer(kind = 4), intent(in) :: i, j
      real(kind = 8), intent(in) :: dz, dx, zsa, xsa, vzero

      t_ana = vzero * ( ( ( dfloat(i) - zsa ) * dz )**2.d0 &
                      + ( ( dfloat(j) - xsa ) * dx )**2.d0 )**0.5d0
      return
    end function t_ana

    ! Function to calculate analytical times in homogeneous model + derivatives of times
    real(kind = 8) function t_anad(tzc, txc, i, j, dz, dx, zsa, xsa, vzero)
    implicit none
      integer(kind = 4), intent(in) :: i, j
      real(kind = 8), intent(in) :: dz, dx, zsa, xsa, vzero
      real(kind = 8) :: d0
      real(kind = 8), intent(out) :: tzc, txc

      d0 = ( ( dfloat(i) - zsa ) * dz )**2.d0 &
           + ( ( dfloat(j) - xsa ) * dx )**2.d0
      t_anad = vzero * (d0**0.5d0)
      if ( d0 .gt. 0.d0 ) then
        tzc = ( d0**(-0.5d0) ) * ( dfloat(i) - zsa ) * dz * vzero
        txc = ( d0**(-0.5d0) ) * ( dfloat(j) - xsa ) * dx * vzero
      else
        tzc = 0.d0
        txc = 0.d0
      end if
      return
    end function t_anad
    
        subroutine sweep2d(slow, tt, nz, nx, dz, dx, &
                       zsi, xsi, zsa, xsa, vzero, eps)
      implicit none
                       
      real(kind = 8), intent(in) :: slow(nz,nx)
      real(kind = 8), intent(inout) :: tt(nz,nx)
      integer(kind = 4), intent(in) :: nz, nx, zsi, xsi, eps
      real(kind = 8), intent(in) :: dz, dx, zsa, xsa, vzero
      real(kind = 8) :: tzc, txc
      integer(kind = 4) :: i, j, sgntz, sgntx, sgnvz, sgnvx
      integer(kind = 4) :: i1, j1, imin
      real(kind = 8) :: dzi, dxi, dz2i, dx2i
      real(kind = 8) :: vref, time_sol(4)
      real(kind = 8) :: tv, te, tev
      real(kind = 8) :: t2d, t1, t2, t3
      real(kind = 8) :: t1d1, t1d2
      real(kind = 8) :: apoly, bpoly, cpoly, dpoly
      real(kind = 8) :: t0c, ta, tb, tauv, taue, tauev

      ! Precalculate constants
      dzi = 1.d0 / dz
      dxi = 1.d0 / dx
      dz2i = 1.d0 / (dz*dz)
      dx2i = 1.d0 / (dx*dx)

        ! First sweep from Top to Bottom,
        ! Sweeping right and left
        sgntz = 1 ; sgnvz = 1
        do i = 2, nz
           sgntx = 1 ; sgnvx = 1
           do j = 2 , nx
                include "Include_FTeik2d.f"
           end do
           sgntx = -1 ; sgnvx = 0
           do j = nx-1, 1, -1
                include "Include_FTeik2d.f"
           end do
        end do
        
        ! Sweeping from bottom to top
        ! Sweeping right and left
        sgntz = -1 ; sgnvz = 0
        do i = nz-1, 1, -1
           sgntx = 1 ; sgnvx = 1
           do j = 2 , nx
                include "Include_FTeik2d.f"
           end do
           sgntx = -1 ; sgnvx = 0
           do j = nx-1, 1, -1
                include "Include_FTeik2d.f"
           end do
        end do
        
        ! First sweep from left to right
        ! Sweeping down and up
        sgntx = 1 ; sgnvx = 1
        do j = 2, nx
          sgntz = 1 ; sgnvz = 1
          do i = 2, nz
            include "Include_FTeik2d.f"
          end do

          sgntz= -1 ; sgnvz = 0
          do i = nz-1, 1, -1
            include "Include_FTeik2d.f"
          end do
        end do

        ! First sweep from right to left
        ! Sweeping down and up
        sgntx = -1 ; sgnvx = 0
        do j = nx-1, 1, -1
          sgntz = 1 ; sgnvz = 1
          do i = 2, nz
            include "Include_FTeik2d.f"
          end do

          sgntz = -1 ; sgnvz = 0
          do i = nz-1, 1, -1
            include "Include_FTeik2d.f"
          end do
        end do
        
      return
end subroutine sweep2d
    
  end subroutine fteik2d


