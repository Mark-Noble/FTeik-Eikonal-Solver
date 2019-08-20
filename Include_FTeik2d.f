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
 ! Index of velocity nodes
  i1 = i - sgnvz
  j1 = j - sgnvx

  ! Get local times of surrounding points
  tv = tt(i-sgntz,j)
  te = tt(i,j-sgntx)
  tev = tt(i-sgntz,j-sgntx)

  ! 1D operators (refracted times)
  vref = min( slow(i1,max(j-1,1)), slow(i1,min(j,nx-1)) )
  t1d1 = tv + dz * vref                   ! First dimension (Z axis)
  vref = min( slow(max(i-1,1),j1), slow(min(i,nz-1),j1) )
  t1d2 = te + dx * vref                   ! Second dimension (X axis)

  ! 2D operators
  t2d = Big
  t1 = Big
  t2 = Big
  t3 = Big
  vref = slow(i1,j1)

  ! Choose plane wave or spherical
  ! Test for plane wave
  if ( abs(i-zsi) .gt. eps .or. abs(j-xsi) .gt. eps ) then

    ! 4 points operator if possible, otherwise do three points
    if ( ( tv .le. te + dx*vref ) .and. ( te .le. tv + dz*vref ) &
      .and. ( te - tev .ge. 0.d0 ) .and. ( tv - tev .ge. 0.d0 ) ) then
      ta = tev + te - tv
      tb = tev - te + tv
      t1 = ( ( tb * dz2i + ta * dx2i ) + sqrt( 4.d0 * vref*vref * ( dz2i + dx2i ) &
           - dz2i * dx2i * ( ta - tb ) * ( ta - tb ) ) ) / ( dz2i + dx2i )

    ! Two 3 points operators
    else if ( ( te - tev ) .le. dz*dz * vref / sqrt( dx*dx + dz*dz ) &
      .and. ( te - tev ) .gt. 0.d0 ) then
      t2 = te + dx * sqrt( vref*vref - ( ( te - tev ) / dz )**2.d0 )

    else if ( ( tv - tev ) .le. dx*dx * vref / sqrt( dx*dx + dz*dz ) &
      .and. ( tv - tev ) .gt. 0.d0 ) then
      t3 = tv + dz * sqrt( vref*vref - ( ( tv - tev ) / dx )**2.d0 )
    end if

  ! Test for spherical
  else
    ! Do spherical operator if conditions ok
    if ( ( tv .lt. te + dx*vref ) .and. ( te .lt. tv + dz*vref ) &
      .and. ( te - tev .ge. 0.d0 ) .and. ( tv - tev .ge. 0.d0 ) ) then
      t0c = t_anad(tzc, txc, i, j, dz, dx, zsa, xsa, vzero)
      tauv = tv - t_ana(i-sgntz, j, dz, dx, zsa, xsa, vzero)
      taue = te - t_ana(i, j-sgntx, dz, dx, zsa, xsa, vzero)
      tauev = tev - t_ana(i-sgntz, j-sgntx, dz, dx, zsa, xsa, vzero)
      ta = tauev + taue - tauv
      tb = tauev - taue + tauv
      apoly = dz2i + dx2i
      bpoly = 4.d0 * ( dfloat(sgntx) * txc * dxi + dfloat(sgntz) * tzc * dzi ) &
              - 2.d0 * ( ta * dx2i + tb * dz2i )
      cpoly = ( ta*ta * dx2i ) + ( tb*tb * dz2i ) &
              - 4.d0 * ( dfloat(sgntx) * txc * dxi * ta + dfloat(sgntz) * tzc * dzi * tb ) &
              + 4.d0 * ( vzero*vzero - vref*vref )
      dpoly = bpoly*bpoly - 4.d0 * apoly * cpoly
      if ( dpoly .ge. 0.d0 ) t1 = 0.5d0 * ( sqrt(dpoly) - bpoly ) / apoly + t0c
      if ( t1-tv .lt. 0.d0 .or. t1-te .lt. 0.d0 ) t1 = Big
    end if

  end if

  t2d = min(t1, t2, t3)

  ! Select minimum time
  time_sol = [ tt(i,j), t1d1, t1d2, t2d ]
  imin = minloc(time_sol, dim = 1)
  tt(i,j) = time_sol(imin)
