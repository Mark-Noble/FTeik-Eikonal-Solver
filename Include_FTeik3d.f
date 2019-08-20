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
  k1 = k - sgnvy

  ! Get local times of surrounding points
  tv = tt(i-sgntz,j,k)
  te = tt(i,j-sgntx,k)
  tn = tt(i,j,k-sgnty)
  tev = tt(i-sgntz,j-sgntx,k)
  ten = tt(i,j-sgntx,k-sgnty)
  tnv = tt(i-sgntz,j,k-sgnty)
  tnve = tt(i-sgntz,j-sgntx,k-sgnty)

  ! 1D operators (refracted times)
  vref = min( slow(i1,max(j-1,1),max(k-1,1)), slow(i1,max(j-1,1),min(k,ny-1)), &
              slow(i1,min(j,nx-1),max(k-1,1)), slow(i1,min(j,nx-1),min(k,ny-1)) )
  t1d1 = tv + dz * vref                   ! First dimension (Z axis)
  vref = min( slow(max(i-1,1),j1,max(k-1,1)), slow(min(i,nz-1),j1,max(k-1,1)), &
              slow(max(i-1,1),j1,min(k,ny-1)), slow(min(i,nz-1),j1,min(k,ny-1)) )
  t1d2 = te + dx * vref                   ! Second dimension (X axis)
  vref = min( slow(max(i-1,1),max(j-1,1),k1), slow(max(i-1,1),min(j,nx-1),k1), &
              slow(min(i,nz-1),max(j-1,1),k1), slow(min(i,nz-1),min(j,nx-1),k1) )
  t1d3 = tn + dy * vref                   ! Third dimension (Y axis)
  t1d = min(t1d1, t1d2, t1d3)

  ! 2D operators
  ! ZX plane
  t2d1 = Big
  vref = min( slow(i1,j1,max(k-1,1)), slow(i1,j1,min(k,ny-1)) )
  if ( ( tv .lt. te + dx*vref ) .and. ( te .lt. tv + dz*vref ) ) then
    ta = tev + te - tv
    tb = tev - te + tv
    t2d1 = ( ( tb * dz2i + ta * dx2i ) + sqrt( 4.d0 * vref*vref * ( dz2i + dx2i ) &
           - dz2i * dx2i * ( ta - tb ) * ( ta - tb ) ) ) / ( dz2i + dx2i )
  end if

  ! ZY plane
  t2d2 = Big
  vref = min( slow(i1,max(j-1,1),k1), slow(i1,min(j,nx-1),k1) )
  if ( ( tv .lt. tn + dy*vref ) .and. ( tn .lt. tv + dz*vref ) ) then
    ta = tv - tn + tnv
    tb = tn - tv + tnv
    t2d2 = ( ( ta * dz2i + tb * dy2i ) + sqrt( 4.d0 * vref*vref * ( dz2i + dy2i ) &
           - dz2i * dy2i * ( ta - tb ) * ( ta - tb ) ) ) / ( dz2i + dy2i )
  end if

  ! XY plane
  t2d3 = Big
  vref = min( slow(max(i-1,1),j1,k1), slow(min(i,nz-1),j1,k1) )
  if ( ( te .lt. tn + dy*vref ) .and. ( tn .lt. te + dx*vref ) ) then
    ta = te - tn + ten
    tb = tn - te + ten
    t2d3 = ( ( ta * dx2i + tb * dy2i ) + sqrt( 4.d0 * vref*vref * ( dx2i + dy2i ) &
           - dx2i * dy2i * ( ta - tb ) * ( ta - tb ) ) ) / ( dx2i + dy2i )
  end if

  t2d = min(t2d1, t2d2, t2d3)

  ! 3D operator
  t3d = Big
  vref = slow(i1,j1,k1)
  ta = te - 0.5d0 * tn + 0.5d0 * ten - 0.5d0 * tv + 0.5d0 * tev - tnv + tnve
  tb = tv - 0.5d0 * tn + 0.5d0 * tnv - 0.5d0 * te + 0.5d0 * tev - ten + tnve
  tc = tn - 0.5d0 * te + 0.5d0 * ten - 0.5d0 * tv + 0.5d0 * tnv - tev + tnve
  if ( min(t1d, t2d) .gt. max(tv, te, tn) ) then
    t2 = vref*vref * dsum * 9.d0
    t3 = dz2dx2 * ( ta - tb ) * ( ta - tb ) &
       + dz2dy2 * ( tb - tc ) * ( tb - tc ) &
       + dx2dy2 * ( ta - tc ) * ( ta - tc )
    if ( t2 .ge. t3 ) then
      t1 = tb * dz2i + ta * dx2i + tc * dy2i
      t3d = ( t1 + sqrt( t2 - t3 ) ) / dsum
    end if
  end if

  ! Select minimum time
  time_sol = [ tt(i,j,k), t1d1, t1d2, t1d3, t2d1, t2d2, t2d3, t3d ]
  imin = minloc(time_sol, dim = 1)
  tt(i,j,k) = time_sol(imin)


