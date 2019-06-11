! MESRA software
! Molecular Electronic Structure Reorganization: Analysis
! Copyright (C) 2019 Thibaud Etienne
! More information at mesrasoftware.wordpress.com
! 
! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License v2
! as published by the Free Software Foundation.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to
! 
! Free Software Foundation, Inc. 
! 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

subroutine alphaddagger

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
  call alpha_ddagger(detachmentU,attachmentU,detachmentR,attachmentR,fov,xLA)
 enddo
else
 xLA = 0.01d0*xlA
 call alpha_ddagger(detachmentU,attachmentU,detachmentR,attachmentR,fov,xLA)
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
  call alpha_ddagger(detachmentUalpha,attachmentUalpha,detachmentRalpha,attachmentRalpha,fov,xLA)
 enddo
else
 xLA = 0.01d0*xlA
 call alpha_ddagger(detachmentUalpha,attachmentUalpha,detachmentRalpha,attachmentRalpha,fov,xLA)
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
  call alpha_ddagger(detachmentUbeta,attachmentUbeta,detachmentRbeta,attachmentRbeta,fov,xLA)
 enddo
else
 xLA = 0.01d0*xlA
 call alpha_ddagger(detachmentUbeta,attachmentUbeta,detachmentRbeta,attachmentRbeta,fov,xLA)
endif

endif

end
