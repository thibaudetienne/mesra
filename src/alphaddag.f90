subroutine alphaddag

! This routine launches the subroutines according to the difference between the number of alpha
! and beta electrons.

use declare

! shell_statement = 0 (1) for closed-(open-)shell molecules.

! fov = overlap file name.

if (shell_statement .eq. 0) then

if (scanLA) then
 do iteration=0,100
   write(6,'(a11,f5.2)') 'x value: ', iteration*0.01d0
   write(6,*)
   write(50,*)
   write(50,'(a11,f5.2)') 'x value: ', iteration*0.01d0
   write(50,*)
  xLA = iteration*0.01d0
  call alpha_ddag(detachmentU,attachmentU,detachmentR,attachmentR,fov,xLA)
 enddo
else
 xLA = 0.01d0*xlA
 call alpha_ddag(detachmentU,attachmentU,detachmentR,attachmentR,fov,xLA)
endif

else if (shell_statement .eq. 1) then

write(6,*) '# Part A - Alpha density matrices'
write(6,*)
write(50,*)
write(50,*) '# Part A - Alpha density matrices'
write(50,*)

if (scanLA) then
 do iteration=0,100
   write(6,'(a11,f5.2)') 'x value: ', iteration*0.01d0
   write(6,*)
   write(50,*)
   write(50,'(a11,f5.2)') 'x value: ', iteration*0.01d0
   write(50,*)
  xLA = iteration*0.01d0
  call alpha_ddag(detachmentUalpha,attachmentUalpha,detachmentRalpha,attachmentRalpha,fov,xLA)
 enddo
else
 xLA = 0.01d0*xlA
 call alpha_ddag(detachmentUalpha,attachmentUalpha,detachmentRalpha,attachmentRalpha,fov,xLA)
endif

write(6,*) '# Part B - Beta density matrices'
write(6,*)
write(50,*) '# Part B - Beta density matrices'
write(50,*)

if (scanLA) then
 do iteration=0,100
   write(6,'(a11,f5.2)') 'x value: ', iteration*0.01d0
   write(6,*)
   write(50,*)
   write(50,'(a11,f5.2)') 'x value: ', iteration*0.01d0
   write(50,*)
  xLA = iteration*0.01d0
  call alpha_ddag(detachmentUbeta,attachmentUbeta,detachmentRbeta,attachmentRbeta,fov,xLA)
 enddo
else
 xLA = 0.01d0*xlA
 call alpha_ddag(detachmentUbeta,attachmentUbeta,detachmentRbeta,attachmentRbeta,fov,xLA)
endif

endif

end
