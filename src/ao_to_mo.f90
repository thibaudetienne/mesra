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

subroutine ao_to_mo(matao,matmo)

! Converts matrices from the AO space to the MO space.

use declare

real*8 :: matmo(norb,norb)
real*8 :: matao(nbs,nbs)

! If S is the atomic functions overlap matrix, C the LCAO coefficients matrix, Cdagger its adjoint, 
! and matao an nbs x nbs matrix, then
! matmo = (C^dagger) x S x matao x S x C.

! tLK stands for "transformation from L into K", referring to the fact
! that we have L canonical orbitals written in the basis of K
! atomic functions. This is also valid for tKL.

! Allocates tLK and tKL, initializes them and matmo.

allocate(tLK(norb,nbs))
allocate(tKL(nbs,norb))

tLK = 0.0d0
tKL = 0.0d0
matmo = 0.0d0

! Performs the matrix multiplication mentioned above, in four steps, 
! and deallocates the intermediate, transformation matrices.

tKL = matmul(S,C)
tKL = matmul(matao,tKL)
tLK = matmul(Cdagger,S)

matmo = matmul(tLK,tKL)

deallocate(tLK,tKL)

end
