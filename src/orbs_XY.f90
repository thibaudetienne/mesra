subroutine orbs_XY

! Computes the transition orbitals for XY-type excited-state calculations.
! Three types of transition orbitals can be calculated.

! For the moment, mesra only accepts calculations where the number of virtual orbitals
! is greater than the number of occupied orbitals, which is the case of the majority
! of calculations.

use declare

if (unr) then

if (nea .gt. norb-nea) then
write(6,*) 'Alpha occupied space bigger than virtual space'
write(6,*) 'Exit'
stop
endif

if (neb .gt. norb-neb) then
write(6,*) 'Beta occupied space bigger than virtual space'
write(6,*) 'Exit'
stop
endif

else

if (nea .gt. norb-nea) then
write(6,*) 'Occupied space bigger than virtual space'
write(6,*) 'Exit'
stop
endif

endif

! Allocates and reads the transition density matrix in the AO space.

allocate(tK(nbs,nbs))

call rTK

if (unr) then

write(6,*) '# Part A - Alpha orbitals'
write(6,*)

endif

! nel is a key number used because it allows us to use the same subroutines for alpha, beta 
! or closed-shell transformations. So does the lcao_coeff_mat array.

 if (nea .eq. 0) then
  continue
 else
  nel = nea
  countunr = 1
  allocate(lcao_coeff_mat(nbs,norb))
  lcao_coeff_mat = C
  call svd_and_rotation_XY
  deallocate(lcao_coeff_mat,t)
 endif

! Since t_1, t_2, and t_3 are allocated when using svd_and_rotation_XY, if we use 
! this subroutine twice (once for alpha and once for beta electrons), these matrices
! should be deallocated between the two procedures.

if (unr) then
 if (nea .eq. 0) then
  continue
 else
  deallocate(t_1,t_2,t_3)
 endif

 if (neb .eq. 0) then
write(6,*)
 else

! For open-shell molecules, repeats these operations for beta electrons.

write(6,*)
write(6,*) '# Part B - Beta orbitals'
write(6,*)

 nel = neb
 tK = tKb 
 countunr = 2

 allocate(lcao_coeff_mat(nbs,norb))
  lcao_coeff_mat = Cb
  call svd_and_rotation_XY
 deallocate(lcao_coeff_mat)

 endif

endif

! Deallocates the appropriate matrices.

deallocate(tK)

if (unr) then
 if (neb .eq. 0) then
  continue
 else
deallocate(t_1,t_2,t_3)
deallocate(tKb)
 endif
else
 deallocate(t_1,t_2,t_3)
endif

end
