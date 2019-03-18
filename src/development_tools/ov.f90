program ov
implicit none
real*8,allocatable :: oGS(:,:),oES(:,:),overlap(:,:)
real*8,allocatable :: ToGS(:,:),ToES(:,:),norm_o(:),temp_mat(:,:)
character*128 :: input1,input2,input3,input4
integer :: i,j,nblocks,targ,il,ih
real*8 :: sum_o,highest,lowest
real*8,allocatable :: diff(:,:),diff_targ(:)

! reads the inputs arguments

! input 1 is the GS orbitals file

call getarg(1,input1)
input1=trim(input1)
open(10,file=input1,form='formatted')

! input 2 is the ES orbitals file

call getarg(2,input2)
input2=trim(input2)
open(11,file=input2,form='formatted')

! input 3 is the number of orbitals

call getarg(3,input3)
input3=trim(input3)

read(input3,'(i5)') nblocks

! input 4 is the target ES orbital to be projected on the GS ones

call getarg(4,input4)
input4=trim(input4)

read(input4,'(i5)') targ

! skips the one-line header of GS & ES orbitals files

read(10,*)
read(11,*)

! allocates the arrays

allocate(oGS(nblocks,nblocks))
allocate(oES(nblocks,nblocks))
allocate(overlap(nblocks,nblocks))
allocate(ToGS(nblocks,nblocks))
allocate(ToES(nblocks,nblocks))
allocate(norm_o(nblocks))
allocate(temp_mat(nblocks,nblocks))
allocate(diff(nblocks,nblocks))
allocate(diff_targ(nblocks))

! initializes the arrays

oGS = 0.0d0
oES = 0.0d0
overlap = 0.0d0
ToGS = 0.0d0
ToES = 0.0d0
norm_o = 0.0d0
temp_mat = 0.0d0
diff = 0.0d0
diff_targ = 0.0d0

! reads the GS & ES orbitals

do i=1,nblocks
 read(10,*)
 read(11,*)
 read(10,'(5D15.8)') (oGS(j,i),j=1,nblocks)
 read(11,'(5D15.8)') (oES(j,i),j=1,nblocks)
enddo

! computes the norm of GS orbitals

do j=1,nblocks
 do i=1,nblocks
 norm_o(j) = norm_o(j) + (oGS(i,j))**2
 enddo
enddo

norm_o = dsqrt(norm_o)

! normalizes the GS orbitals

do i=1,nblocks
 do j=1,nblocks
 oGS(i,j) = oGS(i,j)/norm_o(j)
 enddo
enddo 

! checks the orthonormality of the GS orbitals

ToGS = transpose(oGS)

temp_mat = matmul(ToGS,oGS)

sum_o = 0.0d0

do i=1,nblocks
  sum_o = sum_o + temp_mat(i,i)
enddo

write(6,'(A30,F10.4)') 'Trace of normalized GS UdagU', sum_o

! computes the norm of ES orbitals

norm_o = 0.0d0

do j=1,nblocks
 do i=1,nblocks
 norm_o(j) = norm_o(j) + (oES(i,j))**2
 enddo
enddo

norm_o = dsqrt(norm_o)

! normalizes the ES orbitals

do i=1,nblocks
 do j=1,nblocks
 oES(i,j) = oES(i,j)/norm_o(j)
 enddo
enddo 

! checks the orthonormality of the GS orbitals

ToES = transpose(oES)

temp_mat = matmul(ToES,oES)

sum_o = 0.0d0

do i=1,nblocks
  sum_o = sum_o + temp_mat(i,i)
enddo

write(6,'(A30,F10.4)') 'Trace of normalized ES UdagU', sum_o

! computes the overlap matrix

! ToGS=transpose(oGS)

! overlap=matmul(ToGS,oES)

! do i=1,nblocks
!  write(6,'(2i5,F10.4)') i,targ,overlap(i,targ)
! enddo

! compute deviation of GS orbitals coefficients wrt targ ES orbital ones

 do i=1,nblocks
  do j=1,nblocks
  diff(i,j) = oES(i,targ) - oGS(i,j)
  diff(i,j) = dabs(diff(i,j))
  enddo
 enddo

 ! computes the cumulated deviations

 do i=1,nblocks
  do j=1,nblocks
 diff_targ(i) = diff_targ(i) + diff(j,i)
  enddo
 enddo

 ! finds the orbital corresponding to the highest deviation

 highest = 0.0d0

 do i=1,nblocks
 if (diff_targ(i) .gt. highest) then
         highest = diff_targ(i)
         ih = i
 endif
 write(6,'(2i5,f10.4)') i,targ,diff_targ(i)
 enddo

 write(6,*) 'Highest deviation'
 write(6,'(F10.4,i5)') highest,ih

 ! finds the orbital corresponding to the lowest deviation

 lowest = diff_targ(1)

 do i=2,nblocks
 if (diff_targ(i) .lt. lowest) then
         lowest = diff_targ(i)
         il = i
 endif
 enddo

 write(6,*) 'Lowest deviation'
 write(6,'(F10.4,i5)') lowest,il

 ! deallocates arrays and closes files

deallocate(oGS,oES,overlap,ToGS,ToES,temp_mat,norm_o,diff)
deallocate(diff_targ)

close(10)
close(11)

end

