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

subroutine double_mat_prod(mat1,mat2,mat3,mat4,dimmat)

! Performs a double matrix product with one subroutine.

integer :: dimmat
real*8 :: mat1(dimmat,dimmat),mat2(dimmat,dimmat),mat3(dimmat,dimmat)
real*8 :: mat4(dimmat,dimmat),temp_mat(dimmat,dimmat)

temp_mat = matmul(mat1,mat2)
mat4 = matmul(temp_mat,mat3)


end
