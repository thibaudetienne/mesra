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

subroutine mo_to_ao(matmo,matao)

! Converts matrices from the MO space to the AO space.

use declare

real*8 :: matmo(norb,norb)
real*8 :: matao(nbs,nbs)

! If S is the atomic functions overlap matrix, C the LCAO coefficients matrix, Cdagger its adjoint, 
! and matmo an norb x norb matrix, then
! matao = C x matmo x Cdagger.

! tLK stands for "transformation from L into K", referring to the fact
! that we have L canonical orbitals written in the basis of K
! atomic functions.

! Allocates tLK, initializes tLK and matao.

allocate(tLK(norb,nbs))

tLK = 0.0d0
matao = 0.0d0

! Performs the matrix multiplication mentioned above, in two steps, and deallocates the intermediate, transformation matrix.

tLK=matmul(matmo,Cdagger)
matao=matmul(C,tLK)

deallocate(tLK)

end
