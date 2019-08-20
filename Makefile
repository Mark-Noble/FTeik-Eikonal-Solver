# 
#                      FTeik PACKAGE VERSION 1.0
#                      ---------------------------
#  Copyright (c) 2019 Mark NOBLE, MINES ParisTech, France
#  Email: mark.noble@mines-paristech.fr
# 
#  2D and 3D Eikonal solver to compute first arrival traveltimes in a heterogeneous
#  isotropic velocity model, with the possibility to use different grid spacing in all
#  directions.
# 
#  RELEASE
#  -------
#  VERSION 1.0: 2019-08-01 , First public release 
# 
#  LEGAL STATEMENT
#  ---------------
#  Copyright (c) 2019 Mark NOBLE, MINES ParisTech, France
# 
#  The FTeik package is free software; you can redistribute it and/or modify it
#  under the terms of the GNU General Public License as published by the Free Software
#  Foundation; either version 2 of the License, or (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
#  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  See the GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License along with this program;
#  if not, write to:
#  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
# 
#  BIBLIOGRAPHY
#  ------------
#  Detailed implementation of local operators and global propagation scheme implemented in
#  this sunroutine are inspired form this paper:
#  M. Noble,  A. Gesret and N. Belayouni, 2014, Accurate 3-D finite difference computation of
#  traveltimes in strongly heterogeneous media Geophys. J. Int.,  199 (3): 1572-158.
# 
#  If you find this algorithm useful, citing this paper would be apreciated
#______________________________________________________________________________
#
# Make changes where necessary to reflect the needs of your system
# Compilers       
FC= gfortran
# Options
FFLAGS= -O3 -ffree-form

#
#______________________________________________________________________________
#
# Place to put Binaries
BINDIR=.

B1	= mainFTeik2d.exe
F1	= mainFTeik2d.o fteik2d.o

B2	= mainFTeik3d.exe
F2	= mainFTeik3d.o fteik3d.o

ALL	= $(B1) $(B2) clean

all: $(ALL)

$(B1): $(F1)
	$(FC) $(FFLAGS) -o $@ $(F1)

$(B2): $(F2)
	$(FC) $(FFLAGS) -o $@ $(F2)


install:
	mv -f $(B1) $(BINDIR)
	mv -f $(B2) $(BINDIR)

clean:
	/bin/rm -f ./*.o ./*.mod *~
	
# Suffixes
SUFFIXES = .o .c .f90
.SUFFIXES: $(SUFFIXES)
# Fortran section
COMPILE.f90=$(FC) $(FFLAGS) -c
.f90.o:
	$(COMPILE.f90) $<
# C section
COMPILE.c=$(CC) $(CFLAGS) -c
.c.o:
	$(COMPILE.c) $<
