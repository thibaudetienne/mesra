subroutine daXY

! Creates the detachment/attachment density matrices for XY-type excited-state calculations.

use declare

! Allocates and reads the AO-space transition density matrix.

allocate(tK(nbs,nbs))
call rTK

write(50,*)

if (unr) then

write(6,*) '# Part A - Alpha density matrices'
write(50,*) '# Part A - Alpha density matrices'
write(50,*)

if (nea .eq. 0) then
continue
else

! nel and countunr are key numbers allowing us to use a single subroutine for both alpha, beta
! electrons of open-shell molecules, as well as for closed-shell molecules.

countunr = 1
nel = nea

! Calls the actual detachment/attachment subroutine.

call da_XY

endif

! For open-shell molecules, repeats the operations for the beta electrons.

if (neb .eq. 0) then
continue
else

countunr = 2

write(6,*) '# Part B - Beta density matrices'
write(50,*) '# Part B - Beta density matrices'
write(50,*)

nel = neb
tK = tKb

call da_XY

endif

else

! For closed-shell molecules, simply launches the detachment/attachment subroutine.

nel = nea
call da_XY

endif

end
