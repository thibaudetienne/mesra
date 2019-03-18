subroutine triangledaK

! Converts two square matrices into their upper triangle, in the form of two vectors.

use declare

allocate(triangled(ntr))
allocate(trianglea(ntr))

! Initializes the counter for the vector components.

k = 0

! Transforms the square matrix elements into vectors components.

do i=1,nbs
 do j=1,i
  k = k + 1
  triangled(k) = gamma_d_ao(i,j)
 enddo
enddo

! Re-initializes the counter for the vector components.

k = 0

! Transforms the square matrix elements into vectors components.

do i=1,nbs
 do j=1,i
  k = k + 1
  trianglea(k) = gamma_a_ao(i,j)
 enddo
enddo

end
