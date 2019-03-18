subroutine cubeop(input1,input2,input3,input4,input5)

implicit none
real*8, allocatable :: t1(:,:,:),t2(:,:,:),t3(:,:,:)
integer :: i,j,k,error,tt,n,n1,n2,n3,n4
integer :: i1,i2,i3,w
character*128 :: input1,input2,input3,input4,input5,lines,dummy
real*8 :: norm,dx,dy,dz,V,sum,sum1,sum2,a,b

! Reads the input parameters; opens and reads the input files. For more
! details, see the Cubegen utility manual.

if (input4 .eq. 'p') a=1.0d0
if (input4 .eq. 'm') a=-1.0d0

if (input5 .eq. 'p') b=1.0d0
if (input5 .eq. 'm') b=-1.0d0

 open(10,file=input1,form='formatted')
 open(11,file=input2,form='formatted')
 open(12,file=input3,form='formatted')

 read(10,*)
 read(10,*)
 read(10,'(i5)') n
 read(10,*) n1,dx,dummy,dummy
 read(10,*) n2,dummy,dy,dummy
 read(10,*) n3,dummy,dummy,dz

  V = dx*dy*dz

! Total number of grid points (tt); dimension of the header (n4).

 tt=n1*n2*n3
 n4=6+n

  rewind(10)

! Copies the header in the output file.

 do k=1,n4
 read(10,'(A80)') lines
 read(11,'(A80)') lines
 write(12,'(A80)') lines
 enddo

! Reads grid entries and writes the sum.

 allocate(t1(n3,n2,n1))
 allocate(t2(n3,n2,n1))
 allocate(t3(n3,n2,n1))

do i1 = 1,n1
do i2 = 1,n2
read(10,'(6E13.5)') (t1(i3,i2,i1) , i3=1,n3)
read(11,'(6E13.5)') (t2(i3,i2,i1) , i3=1,n3)
enddo
enddo

  do i1=1,n1
  do i2=1,n2
  do i3=1,n3
  sum = sum + t1(i3,i2,i1)*V
  sum1 = sum1 + t2(i3,i2,i1)*V
  enddo
  enddo
  enddo

   t3 = a*t1 + b*t2

  sum2 = 0.0d0

  do i1=1,n1
  do i2=1,n2
  do i3=1,n3
  sum2 = sum2 + t3(i3,i2,i1)*V
  enddo
  enddo
  enddo

  do i1=1,n1
  do i2=1,n2
  write(12,'(6E13.5)') (t3(i3,i2,i1),i3=1,n3)
  enddo
  enddo

 deallocate(t1,t2,t3)

 close(10)
 close(11)
 close(12)

write(6,*) 'Integral of cube 1'
write(6,'(f10.4)') sum
write(6,*) 'Integral of cube 2'
write(6,'(f10.4)') sum1
write(6,*) 'Integral of the resulting cube'
write(6,'(f10.4)') sum2
write(6,*)


 end
