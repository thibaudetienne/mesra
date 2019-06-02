subroutine no_adiab_z

! Actual computation of the detachment/attachment density matrices from the fully relaxed
! difference density matrix.

use declare

allocate(gD(norb,norb))
allocate(Uvec(norb,norb))
allocate(lvec(norb))

! Constructs the relaxed difference density matrix.

if (countunr .eq. 1) gD = pxalpharelaxed - palpha
if (countunr .eq. 2) gD = pxbetarelaxed - pbeta
if (countunr .eq. 3) gD = pxrelaxed - p

! Constructs the detachment/attachment density matrices from the relaxed difference density matrix.

if (countunr .eq. 1) call det_at(gD,'Ralpha',norb,Uvec,lvec)
if (countunr .eq. 2) call det_at(gD,'Rbeta',norb,Uvec,lvec)
if (countunr .eq. 3) call det_at(gD,'R',norb,Uvec,lvec)

! Prints the relaxed difference density matrix in an .fchk file.

if (countunr .eq. 1) call print_mat_mo_to_ao_fchk(gD,'DeltaRalpha.fchk')
if (countunr .eq. 2) call print_mat_mo_to_ao_fchk(gD,'DeltaRbeta.fchk')
if (countunr .eq. 3) call print_mat_mo_to_ao_fchk(gD,'DeltaR.fchk')

! Deallocates the working arrays.

deallocate(gD,Uvec,lvec)

end
