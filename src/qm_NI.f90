! MESRA software
! Molecular Electronic Structure Reorganization: Analysis
! Copyright (C) 2019 Thibaud Etienne
! More information at mesrasoftware.wordpress.com
! 
! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License v2
! as published by the Free Software Foundation.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to
! 
! Free Software Foundation, Inc. 
! 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

subroutine qm_NI

! Computation of the density-based descriptors from the numerical integration
! of detachment/attachment density functions on a grid of points.

use declare

! THR is a threshold used when checking that two cubes have consistent structure before
! performing operations on their elements.

THR = 1.0D-9

! Reads the file header of the first cube.

call rdcubehdr(file1,natom,x0,y0,z0,n1,n2,n3, &
& dx,dy,dz,dummy)

! Reads the file header of the second cube.

call rdcubehdr(file2,natom1,x01,y01,z01,n11,n21,n31 &
& ,dx1,dy1,dz1,dummy1)

! Checks the consistency of the two cube files.

call consistency(natom,natom1,n1,n11,n2,n21,n3,n31, &
& dx,dx1,dy,dy1,dz,dz1,x0,x01, &
& y0,y01,z0,z01,THR)

! Reads the density grids and performs numerical integrations.

call rddensities(x0,y0,z0,file1,file2,natom,n3,n2,n1, &
& dx,dy,dz,integ1,integ2, &
& integp, integm, xrp, yrp, zrp, &
& xrm, yrm, zrm, x1c, y1c, z1c, &
& x2c, y2c, z2c, rct, r12ct, phis, 1)

! Outputs the results.

call printingdadens(n1,integ1, integ2, integm, integp, &
& xrm, yrm, zrm, xrp, yrp, zrp, &
& x1c, y1c, z1c, x2c, y2c, z2c, &
& rct, r12ct, phis)

end
