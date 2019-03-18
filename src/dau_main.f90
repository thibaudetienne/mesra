subroutine dau_main

! Creates the detachment/attachment density matrices.

use declare

! Reads the density matrices.

call rP

! Computes the detachment/attachment density matrices
! in the canonical space.

! Allocates the difference density matrix, and its eigenvectors and eigenvalues matrices.

allocate(gD(norb,norb))
allocate(Uvec(norb,norb))
allocate(lvec(norb))

! The difference density matrix is the difference between the excited and ground
! state density matrices.

if (unr) then

write(6,*) '# Part A - Alpha density matrices'
write(50,*) '# Part A - Alpha density matrices'

! countunr is a key number useful for using the same subroutine for alpha and beta mesra computations.
  
 if (nea .eq. 0) then
continue
 else
 countunr = 1
 
 gD = pxalpha - palpha

! Prints the unrelaxed difference density matrix in a file.
 
 call print_mat_mo_to_ao_fchk(gD,'DeltaUalpha.fchk')
 
! Calls the detachment/attachment subroutine.

 call det_at(gD,'Ualpha',norb,Uvec,lvec)
 
 endif

 if (neb .eq. 0) then

continue

 else

! Repeats the operation for beta matrices.

write(6,*) '# Part B - Beta density matrices'
write(50,*) '# Part B - Beta density matrices'

 countunr = 2

! For transforming the matrices from one space to another, and since we are
! working now on the beta matrices, we change the values in the default C and Cdag
! since these two matrices are called in further subroutines.
 
 C = Cb
 Cdag = Cbdag

! Computes the total difference density matrix and prints it.

 gD = gD + pxbeta - pbeta

 call print_mat_mo_to_ao_fchk(gD,'DeltaU.fchk')
 
! Computes the beta difference density matrix

 gD = pxbeta - pbeta

Uvec = 0.0d0
lvec = 0.0d0

 call print_mat_mo_to_ao_fchk(gD,'DeltaUbeta.fchk')
 call det_at(gD,'Ubeta',norb,Uvec,lvec)

  endif 

else

! If the molecule is closed-shell, we do these operations only once.

 gD = px - p
 call det_at(gD,'U',norb,Uvec,lvec)
 call print_mat_mo_to_ao_fchk(gD,'DeltaU.fchk')

endif

! Deallocates the working matrices.

deallocate(gD,Uvec,lvec)

end
