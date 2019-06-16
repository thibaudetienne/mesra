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

subroutine orth(mat,whichmat,nlines,ncols)

use declare

integer :: nlines,ncols
real*8 :: mat(nlines,ncols),matdagger(ncols,nlines)
character*(*) :: whichmat
character*128 :: orth_filename

orth_filename = whichmat//'.orth'

allocate(tLK(ncols,nlines))
allocate(tLL(ncols,ncols))

matdagger = transpose(mat)

tLK = matmul(matdagger,S)
tLL = matmul(tLK,mat)

open(15,file=orth_filename,form='formatted')

x = 0.0d0

do i=1,ncols
 write(15,*) i,tLL(i,i)
 x = x + tLL(i,i)
 if (dabs(1.0d0 - tLL(i,i)) .gt. 0.01d0) then
  write(6,*)
  write(6,*) 'Caution!'
  write(6,*) 'Deviation on the orthonormality of an orbital is superior to one percent'
  write(6,*)
  write(6,'(i5,f12.8)') i,tLL(i,i)
  write(50,*) 'Caution!'
  write(50,*) 'Deviation on the orthonormality of an orbital is superior to one percent'
  write(50,*)
  write(50,'(i5,f12.8)') i,tLL(i,i)
endif
enddo

write(50,*) 'Trace of (',whichmat,'^dag)S',whichmat
write(50,'(f12.5)') x
write(50,*)

deallocate(tLK,tLL)

close(15)

end

