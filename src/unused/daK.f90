subroutine daK(damatmo,damatao)

use declare

real*8 :: damataoS(nbs,nbs)
real*8 :: damatmo(norb,norb)
real*8 :: damatao(nbs,nbs)

allocate (tp(norb,nbs))

tp = matmul(damatmo,Cdag)
damatao  = matmul(C,tp)
damataoS = matmul(damatao,S)

x = 0.0d0

do i=1,nbs
 x = x + damataoS(i,i)
enddo

write(6,*) 'Trace ', x

x = 0.0d0

deallocate(tp)

end
