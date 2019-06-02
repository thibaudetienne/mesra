subroutine printingdadens(sn1,sinteg1, sinteg2, sintegm, sintegp, &
& sxrm, syrm, szrm, sxrp, syrp, szrp, sx1c, sy1c, sz1c, sx2c, sy2c, &
& sz2c, srct, sr12ct, sphis)

! Outputs the results from the numerical integration of detachment/attachment densities.

use declare

integer :: sn1
real*8 :: sinteg1, sinteg2, sintegm, sintegp, sxrm, syrm, szrm, &
& sxrp, syrp, szrp, sx1c, sy1c, sz1c, &
& sx2c, sy2c, sz2c, srct, sr12ct, sphis, ata, sss12
character*128 :: string,string1,string2

! Bohr radius.

ata=0.52917725

! Averages the integral of detachment and attachment densities.

sss12=0.5d0*(sinteg1+sinteg2)

write(6,*) '1. Integration of detachment and attachment densities'
write(6,*) " " 
write(6,*) 'Integral of detachment density '
write(6,'(f13.4)') sinteg1
write(6,*) " " 
write(6,*) 'Integral of attachment density '
write(6,'(f13.4)') sinteg2
write(6,*) " " 
write(6,*) '2. Charge from n-/n+ densities'
write(6,*) " " 
write(6,*) 'chi- from negative density variation is '
write(6,'(f13.4)') (-1)*sintegm
write(6,*) " " 
write(6,*) 'chi+ from positive density variation is '
write(6,'(f13.4)') sintegp
write(6,*) " " 
write(6,*) 'chi and phi integrals'
write(6,'(f13.4,f13.4)') 0.5d0*(sintegp-sintegm),&
&0.5d0*(sintegp-sintegm)/sss12
write(6,*) " " 
write(6,*) '3. Centroid coordinates (bohr) from n-/n+ densities, zeta(+/-) '
write(6,*) " " 
write(6,*) 'Negative centroid coordinates'
write(6,'(3f13.4)') sxrm, syrm, szrm
write(6,*) " " 
write(6,*) 'Positive centroid coordinates'
write(6,'(3f13.4)')  sxrp, syrp, szrp
write(6,*) " "
write(6,*) 'CT distance from p-/p+ in Bohr and Angstrom'
write(6,'(2f13.4)') srct,srct*ata
write(6,*) " " 
write(6,*) '4. Centroid coordinates (bohr) from D/A, zeta'
write(6,*) " " 
write(6,*) 'Detachment density centroid coordinates'
write(6,'(3f13.4)') sx1c, sy1c, sz1c
write(6,*) " " 
write(6,*) 'Attachment density centroid coordinates'
write(6,'(3f13.4)') sx2c, sy2c, sz2c
write(6,*) " " 
write(6,*) 'CT distance from D/A in Bohr and Angstrom'
write(6,'(2f13.4)')  sr12ct,sr12ct*ata
write(6,*) " " 
write(6,*) '5. phiS (D/A overlap) and psi indices calculation'
write(6,*) " " 
write(6,*) 'phiS integral value'
write(6,'(f13.4)') sphis
write(6,*) ' ' 
write(6,*) 'psi descriptor value'
write(6,'(f13.4)') (2.0d0/dacos(-1.0d0))*datan(sss12*sphis/(0.5d0*(sintegp-sintegm)))
write(6,*) ' ' 
write(6,*) 'NB: Centroid coordinates are given with respect '
write(6,*) 'to the cube used for the numerical integration'
write(6,*) ' ' 

 end
