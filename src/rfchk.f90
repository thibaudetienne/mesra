subroutine rfchk

! Reads the .fchk file

use declare

! Initializes the parameters and quantities

nea = 0
neb = 0
nbs = 0
norb = 0
ntr = 0

unr = .false.

! Opens the .fchk file and reads nea, neb, nbs, norb.

open(10,file=fchk,status='old')

do
read(10,'(a80)',iostat=error) line
if (error .ne. 0) exit

if (line(1:25) .eq. 'Number of alpha electrons') then
read(line,'(56x,i5)') nea
endif

if (line(1:24) .eq. 'Number of beta electrons') then
read(line,'(56x,i5)') neb
endif

if (line(1:25) .eq. 'Number of basis functions') then
read(line,'(55x,i6)') nbs
endif

if (line(1:31) .eq. 'Number of independent functions') then
read(line,'(55x,i6)') norb
endif

enddo

! Outputs some information in the log file and on the user screen.

write(50,*) 'Number of alpha and beta electrons'
write(50,*) nea, neb
write(50,*) 'Number of atomic and molecular (spin)orbitals'
write(50,*) nbs, norb
write(6,*) 'Number of alpha and beta electrons'
write(6,*) nea, neb
write(6,*) 'Number of atomic and molecular (spin)orbitals'
write(6,*) nbs, norb

! If the number of alpha and beta electrons is different, we speak in terms of unrestricted SCF calculation.

if (nea .ne. neb) unr = .true.

! Rewinds and closes the file.

close(10)

! Calculates the number of elements in the upper triangle of an nbs*nbs square matrix.

ntr = (nbs*(nbs+1))/2

end
