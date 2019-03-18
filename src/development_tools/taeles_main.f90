program taeles_main

use declare

! Reads input

call rinp

! Reads fchk file

call rfchk

! Extracts MO coefficients

call rC

! Reads the AO overlap matrix

call rS

if (jobtype .eq. 'da') then

call rP

! Computes the detachment/attachment density matrices
! in the canonical space

allocate(U(norb,norb))

call da(U,p,px,norb)

deallocate(U)

! Transforms them into the AO space

allocate(gamma_d_ao(nbs,nbs))
allocate(gamma_a_ao(nbs,nbs))

call daK(gamma_d,gamma_d_ao)
call daK(gamma_a,gamma_a_ao)

! Writes them into two fchk files

call triangledaK

call wfchk(triangled,'detachment.fchk',ntr)
call wfchk(trianglea,'attachment.fchk',ntr) 

endif

if (LA) then
 if (scanLA) then

  do i=0,100

xLA = i*0.01d0
yLA = 1.0d0 - xLA 

write(6,'(2f10.5)') xLA,yLA

call LA_analysis

  enddo

 else

xLA = xLA*0.01d0
yLA = 1.0d0 - xLA

call LA_analysis

 endif

else
 continue
endif

if (relax) call relax_D

if (jobtype .eq. 'orbs') then

call ntos

endif

! Deallocates the matrices

call deal

end
