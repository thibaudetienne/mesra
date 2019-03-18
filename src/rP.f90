subroutine rP

! Reads the ground and excited state density matrices, and
! converts them from the AO space to the MO space.

use declare

! Allocates and initializes the matrices.

allocate(p(norb,norb))
allocate(px(norb,norb))
allocate(pK(nbs,nbs))
allocate(pxK(nbs,nbs))
allocate(pKS(nbs,nbs))
allocate(pxKS(nbs,nbs))

if (unr) then

allocate(palpha(norb,norb))
allocate(pbeta(norb,norb))
allocate(pKalpha(nbs,nbs))
allocate(pKalphaS(nbs,nbs))
allocate(pKbeta(nbs,nbs))
allocate(pKbetaS(nbs,nbs))
allocate(pxalpha(norb,norb))
allocate(pxbeta(norb,norb))
allocate(pxKalpha(nbs,nbs))
allocate(pxKalphaS(nbs,nbs))
allocate(pxKbeta(nbs,nbs))
allocate(pxKbetaS(nbs,nbs))

endif

! If the matrices originate from the Quantum Package:

if (soft .eq. 'qp') then

! Reads the density matrices and transforms them into the AO space.

call rP_qp
call rPx_qp

call mo_to_ao(p,pK)
call mo_to_ao(px,pxK)

! If they originate from Gaussian:

else if (soft .eq. 'g03' .or. soft .eq. 'g09' .or. soft .eq. 'g16') then

! Reads the density matrices and transforms them into the AO space.

call rP_g
call rPx_g

! If the computations involve the Z-vector, reads the appropriate matrices.

if (relax) then

 call rPrelax_g

 if (unr) then
  continue
 else
  allocate(pxrelaxed(norb,norb))
  call ao_to_mo(pxKrelaxed,pxrelaxed)
 endif

endif

! Transforms the matrices into the MO space.

call ao_to_mo(pK,p)
call ao_to_mo(pxK,px)

! If the excited state calculation was an unrestricted SCF, reads the appropriate 
! excited-state density matrices.

if (unr) then
 
 allocate(pxKspin(nbs,nbs))

! If the mesra calculation involves the Z-vector, and the excited-state calculation 
! was an unrestricted SCF, reads the appropriate matrices.

if (relax) then

 allocate(pxKspinRelaxed(nbs,nbs))
 allocate(pxKalphaRelaxed(nbs,nbs))
 allocate(pxKbetaRelaxed(nbs,nbs))
 allocate(pxalpharelaxed(norb,norb))
 allocate(pxbetarelaxed(norb,norb))

! countrlx is a key number in the rPxspin_g subroutines.
! Using countrlx allows to read the alpha and beta density matrices
! in the same subroutine as the density matrix of a restricted calculation.

countrlx = 2
 call rPxspin_g
pxKspinRelaxed = pxKspin
pxKspin = 0.0d0

countrlx = 1
 call rPxspin_g

! With the total and spin density matrices, one can deduce the 
! alpha and beta density matrices.

 pxKalphaRelaxed = (pxKrelaxed + pxKspinRelaxed)/2.0d0
 pxKbetaRelaxed = (pxKrelaxed - pxKspinRelaxed)/2.0d0 
 pxKalpha = (pxK + pxKspin)/2.0d0
 pxKbeta = (pxK - pxKspin)/2.0d0 

else

! If the excited-state calculation is not an unrestricted SCF,
! simply reads the excited-state spin density matrix and computes
! the alpha and beta density matrices from it and the total excited-state
! density matrix.

call rPxspin_g 

 pxKalpha = (pxK + pxKspin)/2.0d0
 pxKbeta = (pxK - pxKspin)/2.0d0 
endif

 deallocate(pxKspin)

! Transforms the ground- and excited-state alpha density matrix in the MO space.

call ao_to_mo(pKalpha,palpha)
call ao_to_mo(pxKalpha,pxalpha)

if (relax) call ao_to_mo(pxKalphaRelaxed,pxalphaRelaxed)

! Replacing the content of C by Cbeta (and Cdag by Cbdag) before calling the
! ao_to_mo subroutine allows to use a unique subroutine for converting AO matrices
! into the MO space.

C = Cb
Cdag = Cbdag

call ao_to_mo(pKbeta,pbeta)
call ao_to_mo(pxKbeta,pxbeta)

if (relax) call ao_to_mo(pxKbetaRelaxed,pxbetaRelaxed)

! Places back the Ca and Cadag values into the default C and Cdag matrices.

C = Ca
Cdag = Cadag

endif

endif

! If the difference density matrix is not built from an adiabatic connection of 
! the Z-vector, the density matrices are contracted with the atomic functions
! overlap matrix, and some verifications are performed relatively to the
! trace of the relevant matrices to check if the number of orbitals and
! electrons is correct.

pKS = matmul(pK,S)
pxKS = matmul(pxK,S)

call trace_mat(p,'p',norb)
call trace_mat(px,'px',norb)

call trace_mat(S,'S',nbs)

call trace_mat(pKS,'pKS',nbs)
call trace_mat(pxKS,'pxKS',nbs)

if (unr) then

pKalphaS = matmul(pKalpha,S)
pxKalphaS = matmul(pxKalpha,S)

pKbetaS = matmul(pKbeta,S)
pxKbetaS = matmul(pxKbeta,S)

call trace_mat(palpha,'palpha',norb)
call trace_mat(pbeta,'pbeta',norb)

call trace_mat(pKalphaS,'pKalphaS',nbs)
call trace_mat(pKbetaS,'pKbetaS',nbs)

call trace_mat(pxalpha,'pxalpha',norb)
call trace_mat(pxbeta,'pxbeta',norb)

call trace_mat(pxKalphaS,'pxKalphaS',nbs)
call trace_mat(pxKbetaS,'pxKbetaS',nbs)

if (relax) then

allocate(pxKalphaRelaxedS(nbs,nbs))
allocate(pxKbetaRelaxedS(nbs,nbs))

pxKalphaRelaxedS = matmul(pxKalphaRelaxed,S)
pxKbetaRelaxedS = matmul(pxKbetaRelaxed,S)

call trace_mat(pxKalphaRelaxedS,'pxKalphaRelaxedS',nbs)
call trace_mat(pxKbetaRelaxedS,'pxKbetaRelaxedS',nbs)

deallocate(pxKalphaRelaxed,pxKbetaRelaxed)

endif

endif

end
