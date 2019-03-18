subroutine weigen(eigenfile,eigenvalues,eigendim)

! Writes the eigenvalues into a file.

use declare

integer :: eigendim
character*(*) :: eigenfile
real*8 :: eigenvalues(eigendim)

! Opens the file.

open(10,file=eigenfile,form='formatted')

! Writes the eigenvalues.

do i=1,eigendim
 write(10,'(i5,f10.6)') i,eigenvalues(i)
enddo

! Closes the file.

close(10)

end
