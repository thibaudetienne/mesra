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

subroutine rlxy_PA

! Relaxes the descriptors.

use declare

! Launches rlxy_PAops a different number of times according to whether the 'scanpa' keyword
! has been given in the input.

! PA = population analysis
! scanPA = x is scanned from 0 to 1.

if (scanPA) then

do iteration=0,100

 PA = .true. 
 scanPA = .false.
   write(6,'(a11,f5.2)') 'x value: ', iteration*0.01d0
   write(6,*)
   write(50,*)
   write(50,'(a11,f5.2)') 'x value: ', iteration*0.01d0
   write(50,*)
 xPA = iteration*0.01d0

 call rlxy_PAops

deallocate(p,px,pK,pxK,pKS,pxKS,pxKrelaxed,pxKrelaxedS,pxrelaxed)

enddo

else

call rlxy_PAops

endif

end
