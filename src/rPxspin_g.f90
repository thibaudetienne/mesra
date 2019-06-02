subroutine rPxspin_g

! Reads the relaxed excited state density matrix.

use declare

! Opens the .fchk file where the excited state density matrix is stored.

open(10,file=fchk,status='old')

! Allocates and initializes the working arrays.

allocate(triangle(ntr))

! Locates the starting point where the matrix is written, and reads it and first
! stores it into a vector.

if (relax) then
 if (countrlx .eq. 2) then
  do
   read(10,'(a80)',iostat=error) line
   if (line(1:15) .eq. 'Spin CI Density') then
    read(10,'(5e16.8)') (triangle(i), i=1,ntr)
    exit
   endif
  enddo
 else if (countrlx .eq. 1) then
  do
   read(10,'(a80)',iostat=error) line
   if (line(1:22) .eq. 'Spin CI Rho(1) Density') then
    read(10,'(5e16.8)') (triangle(i), i=1,ntr)
    exit
   endif
  enddo
 endif
else
 do
  read(10,'(a80)',iostat=error) line
  if (line(1:22) .eq. 'Spin CI Rho(1) Density') then
   read(10,'(5e16.8)') (triangle(i), i=1,ntr)
   exit
  endif
 enddo
endif

! Recasts the vector components into a symmetric matrix.

k = 0

do i = 1,nbs
 do j = 1,i
  k = k + 1
  pxKspin(i,j) = triangle(k)
  pxKspin(j,i) = pxKspin(i,j)
 enddo
enddo

deallocate(triangle)

close(10)

end
