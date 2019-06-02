subroutine rPx_g

! Reads the excited state density matrix from a Gaussian file.

use declare

! Opens the "density" file, with the excited state density matrix in the AO space.

open(10,file=fdens,status='old')

! Initializes the relevant matrices and integers:
! pxK is the excited state density matrix in the AO space;
! l is the number of matrix elements stored in the "density" file;
! n0 is the starter.


pxK = 0.0d0
l = 0
n0 = 0

! Finds the value of l, reads the data from the "density file" and stores the entries into the lvec vector.

do
 read(10,'(a80)',iostat=error) line
  if (error .ne. 0) exit
  if (line(1:26) .eq. ' Dump of file   633 length') then
   read(line,'(31x,i9)') l
   allocate(lvec(l))

! If the number of alpha and beta electrons is different, there are nt+1 entries in the density file that might be different from a number and, as such, cannot be properly read by the usual routine. Instead, we count the number of lines in which this might occur and we transform the potential NaN entries into zeros.

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

! Though the first zero entry is generally not NaN, its format is not proper and we should replace it simply by zero.

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

! Selects the right starter n0 according to the odd/even number of electrons.
! Reads the excited state AO-density matrix elements from the lvec vector.
! If the number of alpha electrons (nea) is different from the number of beta electrons,
! The excited state AO-density matrix (pxK) is the sum of alpha (tKK1) and beta contributions (tKK2).

if (unr) then

n0 = nt+1+(ns-1)*2*ntr

allocate(tKK1(nbs,nbs))

tKK1 = 0.0d0

k = 0

do i = 1, nbs
 do j = 1,i
  k = k + 1
  tKK1(i,j) = lvec(n0+k)
  tKK1(j,i) = tKK1(i,j)
 enddo
enddo

allocate(tKK2(nbs,nbs))

tKK2 = 0.0d0

do i = 1, nbs
 do j = 1,i
  k = k + 1
  tKK2(i,j) = lvec(n0+k)
  tKK2(j,i) = tKK2(i,j)
 enddo
enddo

pxK = tKK1 + tKK2

deallocate(tKK1,tKK2)

else

n0 = nt+1+(ns-1)*ntr

k = 0

do i = 1, nbs
 do j = 1,i
  k = k + 1
  pxK(i,j) = lvec(n0+k)*2
  pxK(j,i) = pxK(i,j)
 enddo
enddo

endif

! Deallocates the lvec vector and closes the "density" file.

deallocate(lvec)

close(10)

end
