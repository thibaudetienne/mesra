subroutine common_XY_unrestricted

use declare

!allocate(t(norb,norb))
!allocate(ttdag(norb,norb))
!allocate(ttrsp(norb,norb))

!call ao_to_mo(tK,t)

allocate(SC(nbs,norb))
allocate(Op(nbs,nel))
allocate(Opt(nel,nbs))
allocate(Vp(nbs,norb-nel))
allocate(Vpt(norb-nel,nbs))

SC = matmul(S,lcao_coeff_mat)

do i=1,nel
 do j=1,nbs
  Op(j,i) = SC(j,i)
 enddo
enddo

do i=1,norb-nel
 do j=1,nbs
  Vp(j,i) = SC(j,i+nel)
 enddo
enddo

Opt = transpose(Op)
Vpt = transpose(Vp)

!x = 0.0d0
!
!do i=1,norb
! x = x + t(i,i)
!enddo
!
!write(6,*) x
!
!ttrsp = transpose(t)
!
!ttdag = matmul(t,ttrsp)
!call trace_mat(ttdag,'TTdag',norb)
!
!xy_norm = 0.0d0
!
!do i=1,norb
! xy_norm = xy_norm + ttdag(i,i)
!enddo
!
!write(6,*) 'Transition orbitals normalization factor'
!write(6,*) xy_norm

allocate(xy_X(nel,norb-nel))
allocate(xy_Y(nel,norb-nel))
allocate(xy_Xt(norb-nel,nel))
allocate(xy_Yt(norb-nel,nel))

xy_X = 0.0d0
xy_Xt = 0.0d0
xy_Y = 0.0d0
xy_Yt = 0.0d0

allocate(tLK(nel,nbs))
tLK = 0.0d0
tLK = matmul(Opt,tK)
xy_X = matmul(tLK,Vp)
deallocate(tLK)

allocate(tLK(norb-nel,nbs))
tLK = 0.0d0
tLK = matmul(Vpt,tK)
xy_Yt = matmul(tLK,Op)
deallocate(tLK)

deallocate(SC,Op,Opt,Vp,Vpt)


xy_Xt = transpose(xy_X)

xy_Y = transpose(xy_Yt)


xy_X = xy_X*(dsqrt(2.0d0))/2.0d0
xy_Xt = xy_Xt*(dsqrt(2.0d0))/2.0d0
xy_Y = xy_Y*(dsqrt(2.0d0))/2.0d0
xy_Y = xy_Y*(dsqrt(2.0d0))/2.0d0

write(6,*) 'X matrix'

do i=1,nel
 do j=nel+1,norb
k = j-nel
write(6,'(2i5,f10.5)') i,k+nel,xy_X(i,k)
 enddo
enddo
write(6,*) 'Y matrix'

do i=1,nel
 do j=nel+1,norb
k = j-nel
write(6,'(2i5,f10.5)') i,k+nel,xy_Y(i,k)
 enddo
enddo

xy_norm = 0.0d0
x2y2 = 0.0d0

do i=1,nel
 do j=1,norb-nel
 xy_norm = xy_norm + xy_X(i,j)**2 + xy_Y(i,j)**2
 x2y2 = x2y2 + xy_X(i,j)**2 - xy_Y(i,j)**2
 enddo
enddo

write(6,*) 'Trace of XXt + YYt', xy_norm
write(6,*) 'Trace of XXt - YYt', x2y2

!deallocate(ttrsp,ttdag)

end
