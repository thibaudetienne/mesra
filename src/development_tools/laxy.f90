program laxy
implicit none
character*128 :: ddensF,adensF,overF,line,inpX
integer :: i,j,k,error,nbasis,ntriang
real*8 :: trDS,trAS,trDS12,trAS12,sum1,sum2,sum3,pi,psiM,psiL,psiXY
real*8, allocatable :: dd(:,:),ad(:,:),s(:,:),s12(:,:)
real*8, allocatable :: as(:,:),ds(:,:),ds12(:,:),as12(:,:)
real*8, allocatable :: temp(:),u(:,:),l(:),tempM(:,:),ut(:,:),l12(:)
real*8, allocatable :: diffs12(:),temp1(:),temp2(:),temp3(:),diffs(:)
real*8, allocatable :: sx(:,:),sy(:,:),ddsxy(:,:),adsxy(:,:),diffsxy(:)
real*8 :: x,y,trDSXY,trASXY,qmM,qpM,qmL,qpL,qmXY,qpXY,phM,phL,phXY

pi=dacos(-1.0d0)

write(6,*)

call getarg(1,ddensF)
ddensF=trim(ddensF)
call getarg(2,adensF)
adensF=trim(adensF)
call getarg(3,overF)
overF=trim(overF)
call getarg(4,inpX)

read(inpX,'(f10.4)') x

y = 1 - x

write(6,'(a30,2f10.2)') 'x and y ',x,y
write(6,*)

open(10,file=ddensF,status='old')
open(11,file=adensF,status='old')
open(12,file=overF,status='old')


do
 read(10,'(A80)',iostat=error) line
 if (error .ne. 0) exit
 if (line(1:25) .eq. 'Number of basis functions') then
         read(line,'(55x,i6)') nbasis
 endif
 if (line(1:17) .eq. 'Total SCF Density') exit
enddo

ntriang=nbasis*(nbasis+1)/2

do
 read(11,'(A80)',iostat=error) line
 if (error .ne. 0) exit
 if (line(1:17) .eq. 'Total SCF Density') exit
enddo

allocate(dd(nbasis,nbasis))
allocate(ad(nbasis,nbasis))
allocate(s(nbasis,nbasis))
allocate(s12(nbasis,nbasis))
allocate(temp(ntriang))
allocate(ds(nbasis,nbasis))
allocate(as(nbasis,nbasis))
allocate(ds12(nbasis,nbasis))
allocate(as12(nbasis,nbasis))
allocate(tempM(nbasis,nbasis))
allocate(l12(nbasis))
allocate(ut(nbasis,nbasis))
allocate(sx(nbasis,nbasis))
allocate(sy(nbasis,nbasis))
allocate(ddsxy(nbasis,nbasis))
allocate(adsxy(nbasis,nbasis))
allocate(diffsxy(nbasis))
allocate(diffs12(nbasis))
allocate(temp1(nbasis))
allocate(temp2(nbasis))
allocate(temp3(nbasis))
allocate(l(nbasis))
allocate(u(nbasis,nbasis))
allocate(diffs(nbasis))

dd = 0.0d0
ad = 0.0d0
s  = 0.0d0
s12 = 0.0d0
temp = 0.0d0
ds = 0.0d0
as = 0.0d0
ds12 = 0.0d0
as12 = 0.0d0
tempM = 0.0d0
l12 = 0.0d0
ut = 0.0d0
sx = 0.0d0
sy = 0.0d0
ddsxy = 0.0d0
adsxy = 0.0d0
diffsxy = 0.0d0
diffs = 0.0d0

qmM = 0.0d0
qpM = 0.0d0
qmL = 0.0d0
qpL = 0.0d0
qmXY = 0.0d0
qpXY = 0.0d0

write(6,'(a30,i10)') 'Number of basis functions ',nbasis
write(6,*)

read(10,'(5E16.8)') (temp(i), i=1,ntriang)

k = 0

do i=1,nbasis
 do j=1,i
 k = k + 1
 dd(i,j) = temp(k)
 dd(j,i) = dd(i,j)
 enddo
enddo

temp = 0.0d0

read(11,'(5E16.8)') (temp(i), i=1,ntriang)

k = 0

do i=1,nbasis
 do j=1,i
 k = k + 1
 ad(i,j) = temp(k)
 ad(j,i) = ad(i,j)
 enddo
enddo

temp = 0.0d0

do
 read(12,'(A80)',iostat=error) line
 if (error .ne. 0) exit
 if (line(1:19) .eq. ' Dump of file   514') then
  read(12,'(1x,5e20.8)') (temp(i), i = 1,ntriang)
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

ds = matmul(dd,s)
as = matmul(ad,s)

trDS = 0.0d0
trAS = 0.0d0

do i=1,nbasis
trDS = trDS + ds(i,i)
trAS = trAS + as(i,i)
enddo

diffs = 0.0d0
temp1 = 0.0d0
temp2 = 0.0d0
temp3 = 0.0d0
temp = 0.0d0

sum1 = 0.0d0
sum2 = 0.0d0
sum3 = 0.0d0

do i=1,nbasis
diffs(i) = as(i,i) - ds(i,i)
if (diffs(i) .lt. 0.0d0) temp1(i) = -1.0d0*diffs(i)
if (diffs(i) .gt. 0.0d0) temp2(i) = +1.0d0*diffs(i)
sum1 = sum1 + temp1(i)
sum2 = sum2 + temp2(i)
temp3(i) = (as(i,i))*(ds(i,i))
if (temp3(i) .lt. 0.0d0) temp3(i) = -0.0d0*temp3(i)
temp3(i) = (dsqrt(temp3(i)))/trDS
sum3 = sum3 + temp3(i)
enddo

qmM = sum1
qpM = sum2
phM = sum3

psiM = 2*(pi**(-1))*atan(phM*trDS/qpM)

call diag(s,nbasis,l,u)

sum1 = 0.0d0

do i=1,nbasis
sum1 = sum1 + l(i)
enddo

write(6,'(a30,f10.4)') 'Sum of S eigenvalues ', sum1
write(6,*)

do i=1,nbasis
if (l(i) .lt. 0.0d0) then
        l(i) = 0.0d0
endif
l12(i) = dsqrt(l(i))
s12(i,i)=l12(i)
sx(i,i)=(l(i))**x
sy(i,i)=(l(i))**y
enddo

ut=transpose(u)

tempM=matmul(u,s12)
s12=matmul(tempM,ut)

tempM=matmul(u,sx)
sx=matmul(tempM,ut)

tempM=matmul(u,sy)
sy=matmul(tempM,ut)

tempM = 0.0d0

tempM = matmul(s12,dd)

ds12 = matmul(tempM,s12)

tempM = 0.0d0

tempM = matmul(s12,ad)

as12 = matmul(tempM,s12)

tempM = 0.0d0

tempM = matmul(sx,dd)

ddsxy = matmul(tempM,sy)

tempM = 0.0d0

tempM = matmul(sx,ad)

adsxy = matmul(tempM,sy)

write(6,*) '(Mulliken)'
write(6,*)
write(6,'(a30,f10.4)') 'Trace of D*S ',trDS
write(6,*)
write(6,'(a30,f10.4)') 'Trace of A*S ',trAS
write(6,*)
write(6,'(a30,f10.4)') 'qCT- (0,1) ',qmM
write(6,*)
write(6,'(a30,f10.4)') 'qCT+ (0,1) ',qpM
write(6,*)
write(6,'(a30,f10.4)') 'phiS (0,1) ',phM
write(6,*)
write(6,'(a30,f10.4)') 'psi (0,1) ',psiM
write(6,*)

trDS12 = 0.0d0
trAS12 = 0.0d0
trDSXY = 0.0d0
trASXY = 0.0d0

do i=1,nbasis
trDS12 = trDS12 + ds12(i,i)
trAS12 = trAS12 + as12(i,i)
trDSXY = trDSXY + ddsxy(i,i)
trASXY = trASXY + adsxy(i,i)
enddo


write(6,*) '(Lowdin)'
write(6,*)
write(6,'(a30,f10.4)') 'Trace of S12*D*S12 ',trDS12
write(6,*)
write(6,'(a30,f10.4)') 'Trace of S12*A*S12 ',trAS12
write(6,*)

deallocate(temp)
allocate(temp(nbasis))

diffs12 = 0.0d0
temp1 = 0.0d0
temp2 = 0.0d0
temp3 = 0.0d0
temp = 0.0d0

sum1 = 0.0d0
sum2 = 0.0d0
sum3 = 0.0d0

do i=1,nbasis
diffs12(i) = as12(i,i) - ds12(i,i)
if (diffs12(i) .lt. 0.0d0) temp1(i) = -1.0d0*diffs12(i)
if (diffs12(i) .gt. 0.0d0) temp2(i) = +1.0d0*diffs12(i)
sum1 = sum1 + temp1(i)
sum2 = sum2 + temp2(i)
temp3(i) = (as12(i,i))*(ds12(i,i))
if (temp3(i) .lt. 0.0d0) temp3(i) = -0.0d0*temp3(i)
temp3(i) = (dsqrt(temp3(i)))/trDS12
sum3 = sum3 + temp3(i)
enddo

qmL = sum1
qpL = sum2
phL = sum3

psiL = 2*(pi**(-1))*atan(phL*trDS12/qpL)

write(6,'(a30,f10.4)') 'qCT- (0.5,0.5) ',qmL
write(6,*)
write(6,'(a30,f10.4)') 'qCT+ (0.5,0.5) ',qpL
write(6,*)
write(6,'(a30,f10.4)') 'phiS (0.5,0.5) ',phL
write(6,*)
write(6,'(a30,f10.4)') 'psi (0.5,0.5) ',psiL
write(6,*)

diffs12 = 0.0d0
temp1 = 0.0d0
temp2 = 0.0d0
temp3 = 0.0d0
temp = 0.0d0

sum1 = 0.0d0
sum2 = 0.0d0
sum3 = 0.0d0

do i=1,nbasis
diffsxy(i) = adsxy(i,i) - ddsxy(i,i)
if (diffsxy(i) .lt. 0.0d0) temp1(i) = -1.0d0*diffsxy(i)
if (diffsxy(i) .gt. 0.0d0) temp2(i) = +1.0d0*diffsxy(i)
sum1 = sum1 + temp1(i)
sum2 = sum2 + temp2(i)
temp3(i) = (adsxy(i,i))*(ddsxy(i,i))
if (temp3(i) .lt. 0.0d0) temp3(i) = -0.0d0*temp3(i)
temp3(i) = (dsqrt(temp3(i)))/trDSXY
sum3 = sum3 + temp3(i)
enddo

qmXY = sum1
qpXY = sum2
phXY = sum3

psiXY = 2*(pi**(-1))*atan(phXY*trDSXY/qpXY)

write(6,*) '(x,y)'
write(6,*)
write(6,'(a30,f10.4)') 'Trace of Sx*D*Sy ',trDSXY
write(6,*)
write(6,'(a30,f10.4)') 'Trace of Sx*A*Sy ',trASXY
write(6,*)
write(6,'(a30,f10.4)') 'qCT- (x,y) ',qmXY
write(6,*)
write(6,'(a30,f10.4)') 'qCT+ (x,y) ',qpXY
write(6,*)
write(6,'(a30,f10.4)') 'phiS (x,y) ',phXY
write(6,*)
write(6,'(a30,f10.4)') 'psi (x,y) ',psiXY
write(6,*)

deallocate(dd,ad,s,s12,temp,ds,as,ds12,as12,tempM)
deallocate(l,u,l12,ut,diffs12,temp1,temp2,temp3)
deallocate(sx,sy,ddsxy,adsxy,diffsxy,diffs)

close(10)
close(11)
close(12)

end
