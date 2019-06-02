subroutine rPx_qp

! Reads the excited state density matrix from a Quantum Package file.

use declare

! Opens the file containing the excited state density matrix entries, reads it, and closes it.

open(10,file=fdensx,status='old')

read(10,*) nL

do i=1,nL
 read(10,*) px(1:nL,i)
enddo

close(10)

end
