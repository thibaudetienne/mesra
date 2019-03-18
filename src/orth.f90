subroutine orth(mat,whichmat,nlines,ncols)

use declare

integer :: nlines,ncols
real*8 :: mat(nlines,ncols),matdag(ncols,nlines)
character*(*) :: whichmat
character*128 :: orth_filename

orth_filename = whichmat//'.orth'

allocate(tLK(ncols,nlines))
allocate(tLL(ncols,ncols))

matdag = transpose(mat)

tLK = matmul(matdag,S)
tLL = matmul(tLK,mat)

open(15,file=orth_filename,form='formatted')

x = 0.0d0

do i=1,ncols
 write(15,*) i,tLL(i,i)
 x = x + tLL(i,i)
 if (dabs(1.0d0 - tLL(i,i)) .gt. 0.01d0) then
  write(6,*)
  write(6,*) 'Caution!'
  write(6,*) 'Deviation on the orthonormality of an orbital is superior to one percent'
  write(6,*)
  write(6,'(i5,f12.8)') i,tLL(i,i)
endif
enddo

write(50,*) 'Orthonormality test on matrix ', whichmat
write(50,'(f12.5)') x

deallocate(tLK,tLL)

close(15)

end

