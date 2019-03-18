        subroutine phi_S
       
C       Produced at Théorie-Modélisation-Simulation
C       SRSMC - Université de Lorraine
C       Copyright 2014 by:
C       X. Assfeld, A. Monari, T. Very, T. Etienne

C   This file is part of NANCY_EX.

C   NANCY_EX is free software: you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation, either version 3 of the License, or
C   (at your option) any later version.

C   NANCY_EX is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
C   GNU General Public License for more details.

C   You should have received a copy of the GNU General Public License
C   along with NANCY_EX. If not, see <http://www.gnu.org/licenses/>.
  
C   Computes phi_S index based on the overlap
C   between the attachment and detachment densities
c   or the NTOs to quantify the locality of the charge-transfer.
      
       use declare

       implicit none

       THR = 1.D-09

c      outputs the header and reads the input

       call header

       read(5,*) input1

c      calls the proper subroutine, depending on which derivation of
c      phi_S selected in the input

       if (input1 .eq. 'dadens') then
        call phisdadens
       elseif (input1 .eq. 'ntoscouple') then
        call phisntos
       elseif (input1 .eq. 'ntosdensity') then
        call phisntostot
       elseif (input1 .eq. 'ctoscouple') then
        call phisctos
       elseif (input1 .eq. 'ctosdensity') then
        call phisctostot
       else
        write(6,'(a91)') "First argument must be 'dadens', 'ntoscouple'
     $, 'ntosdensity', 'ctoscouple' or 'ctosdensity'"
        write(6,*)
       endif

       end 
