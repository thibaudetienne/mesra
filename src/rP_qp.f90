subroutine rP_qp

! Reads the ground state density matrix from a Quantum Package file.

use declare

p = 0.0d0
nL = 0

! Opens the file containing the ground state density matrix entries, reads it, and closes it.

open(10,file=fdens,status='old')

read(10,*) nL

do i=1,nL
 read(10,*) p(1:nL,i)
enddo

close(10)

end
