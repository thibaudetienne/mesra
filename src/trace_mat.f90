subroutine trace_mat(mat,whichmat,matsize)

! Computes and outputs the trace of a matrix.

use declare

integer :: matsize
real*8 :: mat(matsize,matsize)
character*(*) :: whichmat

! Initializes the counter x and computes the trace of the matrix mat.

x = 0.0d0

do i=1,matsize
 x = x + mat(i,i)
enddo

! Outputs the result, and reinitializes the counter x.

write(50,*) 'Trace of ',whichmat

if (dabs(x) .lt. 0.00001d0) then
 write(50,'(f12.5)') 0.00000
else
 write(50,'(f12.5)') x
endif

x = 0.0d0

end
