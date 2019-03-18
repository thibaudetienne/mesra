subroutine alpha_ddag(ddensU,adensU,ddensR,adensR,overF,inpX)

! This subroutine is self-contained.

implicit none

character*128 :: ddensU,ddensR,adensU,adensR,overF,line
integer :: i,j,k,error,nbasis,ntriang
real*8 :: sum1,sum2,sum3,inpX
real*8, allocatable :: ddU(:,:),adU(:,:),s(:,:)
real*8, allocatable :: ddR(:,:),adR(:,:)
real*8, allocatable :: temp(:),u(:,:),l(:),tempM(:,:),ut(:,:)
real*8, allocatable :: temp1(:),temp2(:),temp3(:)
real*8, allocatable :: sx(:,:),sy(:,:)
real*8, allocatable :: diffsxy(:)
real*8 :: x,y,trDSXYU,trASXYU,trDSXYZ,trASXYZ
real*8 :: qmXY,qpXY,phXY
real*8 :: trDSXY
real*8 :: trASXY
real*8,allocatable :: dsxyDU(:,:),dsxyAU(:,:),dsxyDR(:,:),dsxyAR(:,:)
real*8,allocatable :: delta_Dsxy(:,:),delta_Asxy(:,:)
real*8 :: tr_dsxyDU, tr_dsxyAU, tr_dsxyDR, tr_dsxyAR
real*8 :: qmXY_D, qpXY_D , qmXY_A, qpXY_A
real*8 :: alphaD, alphaA

! Sets the x and y exponents for the genralized population analysis.

x = inpX
y = 1.0d0 - x

! Opens the fchk and overlap files and reads the matrices.

open(10,file=ddensU,status='old')
open(11,file=adensU,status='old')
open(12,file=ddensR,status='old')
open(13,file=adensR,status='old')
open(14,file=overF,status='old')

do
 read(10,'(A80)',iostat=error) line
 if (error .ne. 0) exit
 if (line(1:25) .eq. 'Number of basis functions') then
  read(line,'(55x,i6)') nbasis
 endif
 if (line(1:17) .eq. 'Total SCF Density') exit
enddo

ntriang = nbasis*(nbasis+1)/2

do
 read(11,'(A80)',iostat=error) line
 if (error .ne. 0) exit
 if (line(1:17) .eq. 'Total SCF Density') exit
enddo

do
 read(12,'(A80)',iostat=error) line
 if (error .ne. 0) exit
 if (line(1:17) .eq. 'Total SCF Density') exit
enddo

do
 read(13,'(A80)',iostat=error) line
 if (error .ne. 0) exit
 if (line(1:17) .eq. 'Total SCF Density') exit
enddo

! Allocates the matrices before reading/constructing them.

allocate(ddU(nbasis,nbasis))
allocate(adU(nbasis,nbasis))
allocate(ddR(nbasis,nbasis))
allocate(adR(nbasis,nbasis))

allocate(s(nbasis,nbasis))

allocate(temp(ntriang))

allocate(tempM(nbasis,nbasis))

allocate(ut(nbasis,nbasis))

allocate(sx(nbasis,nbasis))
allocate(sy(nbasis,nbasis))


allocate(diffsxy(nbasis))

allocate(temp1(nbasis))
allocate(temp2(nbasis))
allocate(temp3(nbasis))

allocate(l(nbasis))
allocate(u(nbasis,nbasis))

allocate(dsxyDU(nbasis,nbasis))
allocate(dsxyDR(nbasis,nbasis))
allocate(dsxyAU(nbasis,nbasis))
allocate(dsxyAR(nbasis,nbasis))

allocate(delta_Dsxy(nbasis,nbasis))
allocate(delta_Asxy(nbasis,nbasis))

! Initializes reals and arrays.

temp1 = 0.0d0
temp2 = 0.0d0
temp3 = 0.0d0
temp = 0.0d0

sum1 = 0.0d0
sum2 = 0.0d0
sum3 = 0.0d0

ddU = 0.0d0
adU = 0.0d0
ddR = 0.0d0
adR = 0.0d0

s  = 0.0d0

tempM = 0.0d0

ut = 0.0d0

sx = 0.0d0
sy = 0.0d0

diffsxy = 0.0d0

qmXY = 0.0d0
qpXY = 0.0d0


write(50,*) 'Number of basis functions '
write(50,'(i8)') nbasis

! Reads the unrelaxed detachment density matrix in the AO space

read(10,'(5E16.8)') (temp(i), i=1,ntriang)

k = 0

do i=1,nbasis
 do j=1,i
 k = k + 1
 ddU(i,j) = temp(k)
 ddU(j,i) = ddU(i,j)
 enddo
enddo

temp = 0.0d0

! Reads the unrelaxed attachment density matrix in the AO space

read(11,'(5E16.8)') (temp(i), i=1,ntriang)

k = 0

do i=1,nbasis
 do j=1,i
 k = k + 1
 adU(i,j) = temp(k)
 adU(j,i) = adU(i,j)
 enddo
enddo

! Reads the relaxed detachment density matrix.

temp = 0.0d0

read(12,'(5E16.8)') (temp(i), i=1,ntriang)

k = 0

do i=1,nbasis
 do j=1,i
 k = k + 1
 ddR(i,j) = temp(k)
 ddR(j,i) = ddR(i,j)
 enddo
enddo

! Reads the relaxed attachment density matrix.

temp = 0.0d0

read(13,'(5E16.8)') (temp(i), i=1,ntriang)

k = 0

do i=1,nbasis
 do j=1,i
 k = k + 1
 adR(i,j) = temp(k)
 adR(j,i) = adR(i,j)
 enddo
enddo

! Reads the AOs overlap matrix.

temp = 0.0d0

do
 read(14,'(A80)',iostat=error) line
 if (error .ne. 0) exit
 if (line(1:19) .eq. ' Dump of file   514') then
  read(14,'(1x,5e20.8)') (temp(i), i = 1,ntriang)
  exit
 endif
enddo

k = 0

do i=1,nbasis
 do j=1,i
 k = k + 1
 s(i,j) = temp(k)
 s(j,i) = s(i,j)
 enddo
enddo

! Diagonalizes S, and checks the sum of its eigenvalues.

call diag(s,nbasis,l,u)

do i=1,nbasis
sum1 = sum1 + l(i)
enddo

write(50,*) 'Sum of S eigenvalues '
write(50,'(f10.4)') sum1

write(6,*)
write(50,*)

! Though by definition S is positive definite, checks that no numerical issue was 
! encountered during the diagonalization by checking that eigenvalues are positive.
! Constructs S^x and S^y.

do i=1,nbasis

if (l(i) .lt. 0.0d0) then
 l(i) = 0.0d0
endif

sx(i,i)=(l(i))**x
sy(i,i)=(l(i))**y
enddo

ut=transpose(u)

tempM=matmul(u,sx)
sx=matmul(tempM,ut)

tempM=matmul(u,sy)
sy=matmul(tempM,ut)

! Constructs (S^x)P(S^y) with P being the (un)relaxed detachment/attachment density matrices.

tempM = 0.0d0

tempM = matmul(sx,ddU)

dsxyDU = matmul(tempM,sy)

tempM = 0.0d0

tempM = matmul(sx,adU)

dsxyAU = matmul(tempM,sy)

tempM = 0.0d0

tempM = matmul(sx,ddR)

dsxyDR = matmul(tempM,sy)

tempM = 0.0d0

tempM = matmul(sx,adR)

dsxyAR = matmul(tempM,sy)

! Checks the trace of the constructed matrices.

tr_dsxyDU = 0.0d0
tr_dsxyAU = 0.0d0
tr_dsxyDR = 0.0d0
tr_dsxyAR = 0.0d0

do i=1,nbasis

tr_dsxyDU = tr_dsxyDU + dsxyDU(i,i)
tr_dsxyAU = tr_dsxyAU + dsxyAU(i,i)
tr_dsxyDR = tr_dsxyDR + dsxyDR(i,i)
tr_dsxyAR = tr_dsxyAR + dsxyAR(i,i)

enddo

! Divides each element of a matrix by its current trace.

dsxyDU =  ((tr_dsxyDU**(-1))*dsxyDU)
dsxyAU =  ((tr_dsxyAU**(-1))*dsxyAU)
dsxyDR =  ((tr_dsxyDR**(-1))*dsxyDR)
dsxyAR =  ((tr_dsxyAR**(-1))*dsxyAR)

do i=1,nbasis
if (dsxyDU(i,i) .lt. 0.0d0) dsxyDU(i,i) = 0.0d0
if (dsxyAU(i,i) .lt. 0.0d0) dsxyAU(i,i) = 0.0d0
if (dsxyDR(i,i) .lt. 0.0d0) dsxyDR(i,i) = 0.0d0
if (dsxyAR(i,i) .lt. 0.0d0) dsxyAR(i,i) = 0.0d0
enddo

! Write intermediate difference matrices.

delta_Dsxy = dsxyDR - dsxyDU
delta_Asxy = dsxyAR - dsxyAU

deallocate(temp)
allocate(temp(nbasis))

! Computes the detachment/attachment alpha^ddag.

sum1 = 0.0d0
sum2 = 0.0d0

do i=1,nbasis
diffsxy(i) = delta_Dsxy(i,i)
if (diffsxy(i) .lt. 0.0d0) temp1(i) = -1.0d0*diffsxy(i)
if (diffsxy(i) .gt. 0.0d0) temp2(i) = +1.0d0*diffsxy(i)
sum1 = sum1 + temp1(i)
sum2 = sum2 + temp2(i)
enddo

qmXY_D = sum1
qpXY_D = sum2

alphaD = 0.5d0*(tr_dsxyDR - tr_dsxyDU)*(qmXY_D+qpXY_D)

sum1 = 0.0d0
sum2 = 0.0d0

do i=1,nbasis
diffsxy(i) = delta_Asxy(i,i)
if (diffsxy(i) .lt. 0.0d0) temp1(i) = -1.0d0*diffsxy(i)
if (diffsxy(i) .gt. 0.0d0) temp2(i) = +1.0d0*diffsxy(i)
sum1 = sum1 + temp1(i)
sum2 = sum2 + temp2(i)
enddo

qmXY_A = sum1
qpXY_A = sum2

alphaA = 0.5d0*(tr_dsxyAR - tr_dsxyAU)*(qmXY_A+qpXY_A)

! Outputs the data.

write(6,*) 'lambda^Z '
write(6,'(f10.4)') 0.5d0*(tr_dsxyDR+tr_dsxyAR) - 0.5d0*(tr_dsxyDU+tr_dsxyAU)
write(6,*) 'Trace of Sx*DU*Sy '
write(6,'(f10.4)') tr_dsxyDU
write(6,*) 'Trace of Sx*AU*Sy'
write(6,'(f10.4)')tr_dsxyAU
write(6,*) 'Trace of Sx*DR*Sy '
write(6,'(f10.4)') tr_dsxyDR
write(6,*) 'Trace of Sx*AR*Sy '
write(6,'(f10.4)') tr_dsxyAR
write(6,*) 'alpha^ddag(D)'
write(6,'(f10.4)') alphaD 
write(6,*) 'alpha^ddag(A)'
write(6,'(f10.4)')  alphaA
write(6,*) 'alpha^ddag'
write(6,'(f10.4)') 0.5d0*(alphaD + alphaA)
write(6,*)
write(50,*) 'lambda^Z '
write(50,'(f10.4)') 0.5d0*(tr_dsxyDR+tr_dsxyAR) - 0.5d0*(tr_dsxyDU+tr_dsxyAU)
write(50,*) 'Trace of Sx*DU*Sy '
write(50,'(f10.4)') tr_dsxyDU
write(50,*) 'Trace of Sx*AU*Sy'
write(50,'(f10.4)')tr_dsxyAU
write(50,*) 'Trace of Sx*DR*Sy '
write(50,'(f10.4)') tr_dsxyDR
write(50,*) 'Trace of Sx*AR*Sy '
write(50,'(f10.4)') tr_dsxyAR
write(50,*) 'alpha^ddag(D)'
write(50,'(f10.4)') alphaD 
write(50,*) 'alpha^ddag(A)'
write(50,'(f10.4)')  alphaA
write(50,*) 'alpha^ddag'
write(50,'(f10.4)') 0.5d0*(alphaD + alphaA)
write(50,*)

! Deallocates the matrices and closes the files.

deallocate(ddU,adU,s,temp,temp1,temp2,temp3,tempM)
deallocate(l,u,ut)
deallocate(sx,sy,diffsxy)
deallocate(ddR,adR)
deallocate(dsxyDU,dsxyAU,dsxyDR,dsxyAR)
deallocate(delta_Dsxy,delta_Asxy)

close(10)
close(11)
close(12)
close(13)
close(14)

end
