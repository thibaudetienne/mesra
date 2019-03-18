subroutine split(input1)

implicit none
real*8, allocatable :: t1(:,:,:),tp(:),tn(:)
integer :: i,j,k,error,tt,n,n1,n2,n3,n4
integer :: i1,i2,i3,w
character*128 :: input1,input2,input3,lines,dummy
real*8 :: norm,dx,dy,dz,V,sum,sum1

! Reads the name of the cube input file; opens and reads it.

open(10,file=input1,form='formatted')

read(10,*)
read(10,*)
read(10,'(i5)') n
read(10,*) n1,dx,dummy,dummy
read(10,*) n2,dummy,dy,dummy
read(10,*) n3,dummy,dummy,dz

! The elementary volume is V:

V = dx*dy*dz

! The total number of grid points is tt; the dimension of the header is n4.

tt=n1*n2*n3
n4=6+n

! rewind(10)

do i=1,n
read(10,*) dummy
enddo

! Reads grid entries (t1).

allocate(t1(n3,n2,n1))
allocate(tp(tt))
allocate(tn(tt))

do i1 = 1,n1
do i2 = 1,n2
read(10,'(6E13.5)') (t1(i3,i2,i1) , i3=1,n3)
enddo
enddo

! Splits the entries according to their sign.

i = 0
do i1=1,n1
do i2=1,n2
do i3=1,n3
 i = i + 1
  if (t1(i3,i2,i1) .le. 0.0d0) then
   tn(i) = t1(i3,i2,i1)
   tp(i) = 0.0d0 
  else
   tp(i) = t1(i3,i2,i1)
   tn(i) = 0.0d0
  endif
enddo
enddo
enddo

! Integrates the two resulting functions.

do i=1,tt
sum = sum + tn(i)*V
sum1 = sum1 + tp(i)*V
enddo

! Reinitializes k and i, deallocates t1, tn and tp, closes
! the cube file and outputs the results.

k = 0
i = 0

deallocate(t1)
deallocate(tn,tp)

close(10)

write(6,*) 'Sum of the negative entries'
write(6,'(f10.4)') sum
write(6,*) 'Sum of the positive entries'
write(6,'(f10.4)') sum1
write(6,*) 'Integral over all the cube'
write(6,'(f10.4)') sum + sum1
write(6,*)

end
