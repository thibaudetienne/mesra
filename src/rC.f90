subroutine rC

! Reads the LCAO coefficients matrix.

use declare

! Opens the file.

open(10,file=fchk,status='old')

! Allocates the vector (for reading the elements) and the matrix (to dispatch the elements after their first reading).

allocate(lvec(norb*nbs))
lvec = 0.0d0

allocate(C(nbs,norb))
C = 0.0d0

! Reads the fchk file until the orbital energies.

do
 read(10,'(A80)', iostat=error) line
  if (error.ne.0) exit
  if (line(1:22) .eq. 'Alpha Orbital Energies') then
   marker_Ea = line
   exit
  endif
enddo

if (unr) then
 do
  read(10,'(A80)', iostat=error) line
   if (error.ne.0) exit
   if (line(1:21) .eq. 'Beta Orbital Energies') then
    marker_Eb = line
    exit
   endif
 enddo
endif

! Reads the elements in the order given by Gaussian, and stores them into a vector.

do
 read(10,'(a80)',iostat=error) line
  if (error .ne. 0) exit
   if (line(1:21) .eq. 'Alpha MO coefficients') then
    marker_Ca = line
    read(10,'(5e16.8)') (lvec(i), i=1,norb*nbs)
    exit
   endif
enddo

! Dispatches the elements into an nbs*norb matrix.

k = 0

do i = 1,norb
 do j = 1,nbs
  k = k + 1
  C(j,i) = lvec(k)
 enddo
enddo

! Computes the adjoint of C, Cdag.

write(50,*)

if (unr) then
call orth(C,'Calpha',nbs,norb)
else
call orth(C,'C',nbs,norb)
endif

allocate(Cdag(norb,nbs))

Cdag = transpose(C)

! If the number of alpha and beta electrons is different, reads also the beta C matrix, Cb.

if (unr) then

! By default, Ca = C when the coefficients are read.

! Allocation of the Ca array.

allocate(Ca(nbs,norb))

Ca = C

! Allocation and initialization of Cb, re-intialization of lvec to read the Cb elements and store them
! in the form of a vector first.

allocate(Cb(nbs,norb))

Cb = 0.0d0
lvec = 0.0d0

! Reads the elements in the order given by Gaussian, and stores them into a vector.

do
 read(10,'(a80)',iostat=error) line
 if (error .ne. 0) exit
 if (line(1:20) .eq. 'Beta MO coefficients') then
  marker_Cb = line
  read(10,'(5e16.8)') (lvec(i), i=1,norb*nbs)
  exit
 endif
enddo

! Dispatches the Cb elements from the lvec vector into an nbs*norb Cb matrix.

k = 0

do i = 1,norb
 do j = 1,nbs
  k = k + 1
  Cb(j,i) = lvec(k)
 enddo
enddo

! Allocates and generates the adjoint of Ca and Cb.

allocate(Cadag(norb,nbs))
allocate(Cbdag(norb,nbs))

Cadag = transpose(Ca)
Cbdag = transpose(Cb)

call orth(Cb,'Cbeta',nbs,norb)

write(50,*)

endif

! Deallocates the intermediate lvec vector, and closes the .fchk file.

deallocate(lvec)

close(10)

end
