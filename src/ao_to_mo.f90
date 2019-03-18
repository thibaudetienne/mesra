subroutine ao_to_mo(matao,matmo)

! Converts matrices from the AO space to the MO space.

use declare

real*8 :: matmo(norb,norb)
real*8 :: matao(nbs,nbs)

! If S is the atomic functions overlap matrix, C the LCAO coefficients matrix, Cdag its adjoint, 
! and matao an nbs x nbs matrix, then
! matmo = (C^dag) x S x matao x S x C.

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
tLK = matmul(Cdag,S)

matmo = matmul(tLK,tKL)

deallocate(tLK,tKL)

end
