subroutine rS

! Reads the basis functions overlap matrix S.

use declare

! Opens the "overlap" file.

open(10,file=fov,status='old')

! Allocates and initializes the "triangle" array for reading the S elements before dispatching them into an nbs*nbs matrix.

allocate(triangle(ntr))
triangle = 0.0d0

! Reads the S elements and places them into a vector.

do
 read(10,'(a80)',iostat=error) line
  if (error .ne. 0) exit
  if (line(1:19) .eq. ' Dump of file   514') then
   read(10,'(1x,5e20.8)') (triangle(i),i=1,ntr)
   exit
  endif
enddo

! Allocates and initializes S.

allocate(S(nbs,nbs))
S = 0.0d0

! Dispatches the "triangle" vector components into the square S matrix.

k = 0

do i = 1,nbs
 do j = 1,i
  k = k + 1
  S(i,j) = triangle(k)
  S(j,i) = S(i,j)
 enddo
enddo

! Deallocates triangle and closes the overlap file.

deallocate(triangle)

close(10)

end
