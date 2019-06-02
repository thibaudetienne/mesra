subroutine rTK

! Reads the AO-space transition density matrix.

use declare

! First locates the starter in the 'density' file, according to the version of 
! the Gaussian software used and the difference between the number of alpha and
! beta electrons.

if (unr) then
 n0 = nt+1+nt*2*ntr+2*(ns-1)*nbs**2+1
else
 if (soft .eq. 'g03') then
  n0 = nt+1+nt*ntr+(ns-1)*nbs**2+1
 else if (soft .eq. 'g09' .or. soft .eq. 'g16') then
  n0 = nt+1+nt*ntr+2*(ns-1)*nbs**2+1
 endif
endif

! Opens the density file and reads the matrix.

open(10,file=fdens,status='old')

tK = 0.0d0

l = 0

! l is the number of elements to be read.

do
 read(10,'(a80)',iostat=error) line
 if (error .ne. 0) exit
 if (line(1:26) .eq. ' Dump of file   633 length') then
  read(line,'(31x,i9)') l
  allocate(lvec(l))

! If the number of alpha and beta electrons is different, there are nt+1 entries
! in the density file that might be different from a number and, as such, cannot
! be properly read by the usual routine. Instead, we count the number of lines
! in which this might occur and we transform the potential NaN entries into
! zeros.

  if (unr) then

! First, we count the number of lines in which there are NaN appearing.

  k = nint((nt+1)/5.0d0)

! We then allocate a dummy array to read them...

  allocate(dummy_array(k*5))
   read(10,'(1x5A20)') (dummy_array(i), i=1,5*k)

! ... and convert them into zeros.

   do i=1,5*k
    if (trim(adjustl(dummy_array(i))) .eq. 'NaN') then
     lvec(i) = 0.0d0
    else
     read(dummy_array(i),'(e20.8)') lvec(i)
    endif
   enddo

! Though the first zero entry is generally not NaN, its format is not proper and
! we should replace it simply by zero.

   lvec(1) = 0.0d0
   deallocate(dummy_array)

! Finally we read the other entries normally.

   read(10,'(1x5e20.8)') (lvec(i), i = 5*k+1,l-5*k)
  else
   read(10,'(1x5e20.8)') (lvec(i), i = 1,l)
  endif
   exit
  endif
enddo

close(10)

! Converts the lvec vector into the tK square matrix.

k = 0

do i=1,nbs
 do j=1,nbs
  tK(i,j) = lvec(n0 + k)
  k = k + 1
 enddo
enddo

! If the molecule is open-shell, repeats the operation for the beta transition density matrix.

if (unr) then
 if (neb .eq. 0) then
 continue
 else
 n0 = n0 + k
 k = 0
 allocate (tKb(nbs,nbs))
 do i = 1,nbs
  do j = 1,nbs
   tKb(i,j) = lvec(n0 + k)
   k = k + 1
  enddo
 enddo
 endif
endif

deallocate(lvec)

! Performs some tests on the AO-space transition density matrix.
! The trace of TS should be zero.

 if (nea .eq. 0) then
  continue
 else

allocate(TS(nbs,nbs))
TS = matmul(tK,S)
call trace_mat(TS,'TS',nbs)
deallocate(TS)

 endif

! If the molecule is open-shell, repeats the operation for the beta transition density matrix.

if (unr) then
 if (neb .eq. 0) then
 continue
 else
 allocate(TSb(nbs,nbs))
 TSb = matmul(tKb,S)
 call trace_mat(TSb,'TS beta',nbs)
 deallocate(TSb)
 endif
endif

end
