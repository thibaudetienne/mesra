subroutine rinp_new

use declare

call getarg(1,input)
input=trim(input)

open(10,file=input,status='old')

read(10,*) soft

if (soft .eq. 'g03' .or. soft .eq. 'g09') then
 continue
else if (soft .eq. 'g16') then
continue
else if (soft .eq. 'qp') then
continue
else
 write(6,'(a128)') 'Error: Only Gaussian (g03,g09,g16) and Quantum Package (qp) softwares are&
                   & supported.'
 stop
endif

read(10,*) jobtype

if (jobtype(1:8) .eq. 'orbitals' .or. jobtype(1:9) .eq. 'densities') then
 continue
else if (jobtype(1:3) .eq. 'all') then
 continue
else
 write(6,'(a128)') 'Error: only "orbitals", "densities" and "all" job types are&
                   & supported"'
endif

if (jobtype(1:8) .eq. 'orbitals') then

 if (jobtype .eq. 'orbitals(1,1)') then
jobN = 111
  LA = .true.
  scanLA = .true.
  xy = .true.
  read(10,*) nLA
 else if (jobtype .eq. 'orbitals(1)(0)') then
jobN = 110
  LA = .true.
  scanLA = .true.
  xy = .false.
  read(10,*) nLA
 else if (jobtype .eq. 'orbitals(2)(1)') then
jobN = 121
  LA = .true.
  scanLA = .false.
  xy = .true.
  read(10,*) xLA
 else if (jobtype .eq. 'orbitals(2)(0)') then
jobN = 120
  LA = .true.
  scanLA = .false.
  xy = .false.
  read(10,*) xLA
 else if (jobtype .eq. 'orbitals(0)(1)') then
jobN = 101
  LA = .false.
  scanLA = .false.
  xy = .true.
 else if (jobtype .eq. 'orbitals(0)(0)') then
jobN = 100
  LA = .false.
  scanLA = .false.
  xy = .false.
 else
write(6,*) 'Error: Wrong options!'
 endif

else if (jobtype(1:9) .eq. 'densities') then

 if (jobtype .eq. 'densities(1,1)') then
jobN = 211
  LA   = .true.
  scanLA = .true.
  relax  = .true.
  adiab  = .true.
  read(10,*) nLA,nAdiab
 else if (jobtype .eq. 'densities(1)(2)') then
jobN = 212
  LA   = .true.
  scanLA = .true.
  relax  = .true.
  adiab  = .false.
  read(10,*) nLA
 else if (jobtype .eq. 'densities(1)(0)') then
jobN = 210
  LA   = .true.
  scanLA = .true.
  relax  = .false.
  adiab  = .false.
  read(10,*) nLA
 else if (jobtype .eq. 'densities(2)(1)') then
jobN = 221
  LA   = .true.
  scanLA = .false.
  relax  = .true.
  adiab  = .true.
  read(10,*) xLA,nAdiab
 else if (jobtype .eq. 'densities(2)(2)') then
jobN = 222
  LA   = .true.
  scanLA = .false.
  relax  = .true.
  adiab  = .false.
  read(10,*) xLA
 else if (jobtype .eq. 'densities(2)(0)') then
jobN = 220
  LA   = .true.
  scanLA = .false.
  relax  = .false.
  adiab  = .false.
  read(10,*) xLA
 else if (jobtype .eq. 'densities(0)(1)') then
jobN = 201
  LA   = .false.
  scanLA = .false.
  relax  = .true.
  adiab  = .true.
  read(10,*) nAdiab
 else if (jobtype .eq. 'densities(0)(2)') then
jobN = 202
  LA   = .false.
  scanLA = .false.
  relax  = .true.
  adiab  = .false.
 else if (jobtype .eq. 'densities(0)(0)') then
jobN = 200
  LA   = .false.
  scanLA = .false.
  relax  = .false.
  adiab  = .false.
 else
write(6,*) 'Error: Wrong options!'
 endif

else if (jobtype(1:3) .eq. 'all') then

 if (jobtype .eq. 'all(1,1,1)') then
jobN = 3111
  LA   = .true. 
  scanLA = .true.
  relax  = .true.
  adiab  = .true.
  xy     = .true.
  read(10,*) nLA,nAdiab
 else if (jobtype .eq. 'all(1)(1)(0)') then
jobN = 3110
  LA   = .true.
  scanLA = .true.
  relax  = .true.
  adiab  = .true.
  xy     = .false.
  read(10,*) nLa,nAdiab
 else if (jobtype .eq. 'all(1)(2)(1)') then
jobN = 3121
  LA   = .true.
  scanLA = .true.
  relax  = .true.
  adiab  = .false.
  xy     = .true.
  read(10,*) nLA
 else if (jobtype .eq. 'all(1)(2)(0)') then
jobN = 3120
  LA   = .true.
  scanLA = .true.
  relax  = .true.
  adiab  = .false.
  xy     = .false.
  read(10,*) nLA
 else if (jobtype .eq. 'all(1)(0)(1)') then
jobN = 3101
  LA   = .true.
  scanLA = .true.
  relax  = .false.
  adiab  = .false.
  xy     = .true.
  read(10,*) nLA
 else if (jobtype .eq. 'all(1)(0)(0)') then
jobN = 3100
  LA   = .true.
  scanLA = .true.
  relax  = .false.
  adiab  = .false.
  xy     = .false.
  read(10,*) nLA
 else if (jobtype .eq. 'all(2)(1)(1)') then
jobN = 3211
  LA   = .true.
  scanLA = .false.
  relax  = .true.
  adiab  = .true.
  xy     = .true.
  read(10,*) xLA,nAdiab
 else if (jobtype .eq. 'all(2)(1)(0)') then
jobN = 3210
  LA   = .true.
  scanLA = .false.
  relax  = .true.
  adiab  = .true.
  xy     = .false.
  read(10,*) xLA,nAdiab
 else if (jobtype .eq. 'all(2)(2)(1)') then
jobN = 3221
  LA   = .true.
  scanLA = .false.
  relax  = .true.
  adiab  = .false.
  xy     = .true.
  read(10,*) xLA
 else if (jobtype .eq. 'all(2)(2)(0)') then
jobN = 3220
  LA   = .true.
  scanLA = .false.
  relax  = .true.
  adiab  = .false.
  xy     = .false.
  read(10,*) xLA
 else if (jobtype .eq. 'all(2)(0)(1)') then
jobN = 3201
  LA   = .true.
  scanLA = .false.
  relax  = .false.
  adiab  = .false.
  xy     = .true.
  read(10,*) xLA
 else if (jobtype .eq. 'all(2)(0)(0)') then
jobN = 3200
  LA   = .true.
  scanLA = .false.
  relax  = .false.
  adiab  = .false.
  xy     = .false.
  read(10,*) xLA
 else if (jobtype .eq. 'all(0)(1)(1)') then
jobN = 3011
  LA   = .false.
  scanLA = .false.
  relax  = .true.
  adiab  = .true.
  xy     = .true.
  read(10,*) nAdiab
 else if (jobtype .eq. 'all(0)(1)(0)') then
jobN = 3010
  LA   = .false.
  scanLA = .false.
  relax  = .true.
  adiab  = .true.
  xy     = .false.
  read(10,*) nAdiab
 else if (jobtype .eq. 'all(0)(2)(1)') then
jobN = 3021
  LA   = .false.
  scanLA = .false.
  relax  = .true.
  adiab  = .false.
  xy     = .true.
 else if (jobtype .eq. 'all(0)(2)(0)') then
jobN = 3020
  LA   = .false.
  scanLA = .false.
  relax  = .true.
  adiab  = .false.
  xy     = .false.
 else if (jobtype .eq. 'all(0)(0)(1)') then
jobN = 3001
  LA   = .false.
  scanLA = .false.
  relax  = .false.
  adiab  = .false.
  xy     = .true.
 else if (jobtype .eq. 'all(0)(0)(0)') then
jobN = 3000
  LA   = .false.
  scanLA = .false.
  relax  = .false.
  adiab  = .false.
  xy     = .false.
 else
write(6,*) 'Error: Wrong options!'
 endif

endif


fov='overlap'
fdens='density'

read(10,*) fchk

if (soft .eq. 'qp') then
 read(10,*) fdensx
endif

read(10,*) nt,ns

close(10)

end
