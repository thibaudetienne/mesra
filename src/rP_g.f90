subroutine rP_g

! Reads the ground state density matrix from a Gaussian file.

use declare

! Opens the fchk file, allocates and initializes the triangle vector.
! (triangle stands for upper triangle of the symmetric matrix to read)

open(10,file=fchk,status='old')

allocate(triangle(ntr))

! Reads the AO-space ground state density matrix and stores its components
! into the triangle vector.

do
 read(10,'(a80)',iostat=error) line
  if (line(1:17) .eq. 'Total SCF Density') then
   read(10,'(5e16.8)') (triangle(i), i=1,ntr)
   exit
  endif
enddo

! Brings the triangle vector components into matrix elements.

k = 0

do i = 1,nbs
 do j = 1,i
  k = k + 1
  pK(i,j) = triangle(k)
  pK(j,i) = pK(i,j)
 enddo
enddo

if (unr) then

allocate(trianglespin(ntr))
allocate(pKspin(nbs,nbs))

do
 read(10,'(a80)',iostat=error) line
  if (line(1:16) .eq. 'Spin SCF Density') then
   read(10,'(5e16.8)') (trianglespin(i), i=1,ntr)
   exit
  endif
enddo

! Brings the triangle vector components into matrix elements.

k = 0

do i = 1,nbs
 do j = 1,i
  k = k + 1
  pKspin(i,j) = trianglespin(k)
  pKspin(j,i) = pKspin(i,j)
 enddo
enddo

pKalpha = (pK + pKspin)/2.0d0
pKbeta = (pK - pKspin)/2.0d0

endif

! Deallocates the triangle vectors and closes the fchk file.

deallocate(triangle)

if (unr) deallocate(trianglespin,pKspin)

close(10)

end
