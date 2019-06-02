subroutine triangleK(matrix,matriangle,dtmp)

! Converts a square matrix into its upper triangle, in the form of a vector.

use declare

integer :: dtmp
real*8 :: matrix(dtmp,dtmp),matriangle(dtmp)

! Initializes the counter for the components of the vector.

k = 0

! Transforms the square matrix elements into vector components.

do i=1,nbs
 do j=1,i
  k = k + 1
  matriangle(k) = matrix(i,j)
 enddo
enddo

end
