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

subroutine LA_analysis(matrix_dLA,matrix_aLA,x_PA,y_PA)

! Performs the actual D/A population analysis.

use declare

real*8 :: x_PA,y_PA
real*8 :: matrix_dLA(nbs,nbs),matrix_aLA(nbs,nbs)

! Allocates the (S^x)P(S^y) arrays (P = D,A) 

allocate(SxDSy(nbs,nbs))
allocate(SxASy(nbs,nbs))

! Creates S^x, S^y and (S^x)P(S^y) (P = D,A).

call LA_mat(matrix_dLA,matrix_aLA,SxDSy,SxASy,x_PA,y_PA)

! Computes the descriptors from (S^x)P(S^y) (P = D,A).

call QuantumMetricsLA(SxDSy,SxASy)

! Deallocates the (S^x)P(S^y) matrices (P = D,A).

deallocate(SxDSy,SxASy)

end
