subroutine rdagf

! Reads the name of the files containing the unrelaxed detachment/attachment density matrices
! and the detachment/attachment density matrix constructed by diagonalizing the Z-vector.

use declare

! shell_statement = 0 (1) means (un)restricted calculation.

read(10,*) shell_statement

if (shell_statement .eq. 0) then
read(10,*) detachmentU
read(10,*) attachmentU
read(10,*) detachmentR
read(10,*) attachmentR
else if (shell_statement .eq. 1) then
read(10,*) detachmentUalpha
read(10,*) attachmentUalpha
read(10,*) detachmentRalpha
read(10,*) attachmentRalpha
read(10,*) detachmentUbeta
read(10,*) attachmentUbeta
read(10,*) detachmentRbeta
read(10,*) attachmentRbeta
endif

end
