subroutine rlxy_LA

! Relaxes the descriptors.

use declare

! Launches rlxy_LAops a different number of times according to whether the 'scanpa' keyword
! has been given in the input.

! LA = 'Linear Algebra' = population analysis
! scanLA = x is scanned from 0 to 1.

if (scanLA) then

do iteration=0,100

 LA = .true. 
 scanLA = .false.
   write(6,'(a11,f5.2)') 'x value: ', iteration*0.01d0
   write(6,*)
   write(50,*)
   write(50,'(a11,f5.2)') 'x value: ', iteration*0.01d0
   write(50,*)
 xLA = iteration*0.01d0

 call rlxy_LAops

deallocate(p,px,pK,pxK,pKS,pxKS,pxKrelaxed,pxKrelaxedS,pxrelaxed)

enddo

else

call rlxy_LAops

endif

end
