subroutine mo_to_ao(matmo,matao)

! Converts matrices from the MO space to the AO space.

use declare

real*8 :: matmo(norb,norb)
real*8 :: matao(nbs,nbs)

! If S is the atomic functions overlap matrix, C the LCAO coefficients matrix, Cdag its adjoint, 
! and matmo an norb x norb matrix, then
! matao = C x matmo x Cdag.

! tLK stands for "transformation from L into K", referring to the fact
! that we have L canonical orbitals written in the basis of K
! atomic functions.

! Allocates tLK, initializes tLK and matao.

allocate(tLK(norb,nbs))

tLK = 0.0d0
matao = 0.0d0

! Performs the matrix multiplication mentioned above, in two steps, and deallocates the intermediate, transformation matrix.

tLK=matmul(matmo,Cdag)
matao=matmul(C,tLK)

deallocate(tLK)

end
