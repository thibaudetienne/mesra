subroutine inv(mat,matm1,dimmat)

use declare

integer :: dimmat
real*8 :: mat(dimmat,dimmat),matm1(dimmat,dimmat)

allocate (tp(dimmat,dimmat))

tp = mat

call dpotrf('U',dimmat,tp,dimmat,error)

if (error .ne. 0) then
 write(6,*) 'Cholesky error'
 stop
endif

call dpotri('U',dimmat,tp,dimmat,error)

if (error .ne. 0) then
 write(6,*) 'Inversion error'
 stop
endif

k = 0

do i=1,dimmat
 do j=1,i
  k = k + 1
  matm1(i,j) = tp(j,i)
  matm1(j,i) = matm1(i,j)
 enddo
enddo

deallocate(tp)

end
