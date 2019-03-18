subroutine rLA_status

! Reads the status to be decided relatively to the population analysis to perform.

use declare

! Reads the status from the input.

read(10,*) LA_status

! LA = population analysis (Linear Algebra) => which exponents (x)?
! noLA = no population analysis;
! scan = a scan of the x exponents from 0 to 1.

 if (LA_status .eq. 'LA') then
  LA = .true.
  scanLA = .false.
  read(10,*) xLA
 else if (LA_status .eq. 'noLA') then
  LA = .false.
  scanLA = .false.
 else if (LA_status .eq. 'scan') then
  LA = .true.
  scanLA = .true.
 endif

end
