 subroutine rddensities(sx0,sy0,sz0,sfile1,sfile2,snatom,sn3,sn2,&
&            sn1, sdx,sdy,sdz,sinteg1,sinteg2,sintegp,sintegm,&
&            sxrp, syrp, szrp, sxrm, syrm, szrm, sx1c, sy1c, sz1c, &
&            sx2c, sy2c, sz2c, srct,sr12ct,sphis, ss1)

! Reads the density grids and performs operations.

use declare

character*128 :: sfile1,sfile2
integer :: snatom, sn3, sn2, sn1 , ss1, sn1p
real*8 :: sdx, sdy, sdz, sinteg1, sinteg2, sintegp, sintegm, &
& sxrp, syrp, szrp, sxrm, syrm, szrm, &
& sx1c, sy1c, sz1c, sx2c, sy2c, sz2c, &
& srct,sr12ct, sphis, sx0, sy0, sz0

! Ensures that the number of steps in the x direction is a positive integer.

if (sn1 .lt. 0) then
 sn1p = (-1)*sn1
else
 sn1p = sn1
endif

! Opens the two cube files.

open(10,file=sfile1,form='formatted')
open(11,file=sfile2,form='formatted')

! Skips the header when reading the cube files.

do i=1,snatom+6
 read(10,*)
 read(11,*)
enddo

! Initializes the number of distance steps.

sinteg1 = 0.0D1
sinteg2 = 0.0D1

! Allocates the two density functions p1 and p2.

allocate (p1(sn3,sn2,sn1p))
allocate (p2(sn3,sn2,sn1p))

! Reads densities from cubes 1 and 2.

do i1 = 1,sn1p
  do i2 = 1,sn2
   read(10,'(6E13.5)') (p1(i3,i2,i1) , i3=1,sn3)
   read(11,'(6E13.5)') (p2(i3,i2,i1) , i3=1,sn3)
  enddo
enddo

! Performs the numerical integration of cubes 1 and 2 densities
! over all the space span by the cube.

do i1=1,sn1p
  do i2 = 1,sn2
    do i3 = 1,sn3
    sinteg1 = sinteg1 + dabs(p1(i3,i2,i1))*sdx*sdy*sdz
    sinteg2 = sinteg2 + dabs(p2(i3,i2,i1))*sdx*sdy*sdz
    enddo
  enddo
enddo

! Closes the two cube files.

close(10)
close(11)

allocate (diff(sn3,sn2,sn1p))
allocate (pp(sn3,sn2,sn1p)) ! "p" stands for "plus", i.e., positive entries.
allocate (pm(sn3,sn2,sn1p)) ! "m" stands for "minus", i.e., negative entries.

pp = 0.0d01
pm = 0.0d01

! Computes the difference between the densities read from cubes 1 and 2
! and splits the entries according to their sign.

do i1=1,sn1p
 do i2=1,sn2
  do i3=1,sn3
diff(i3,i2,i1) = p2(i3,i2,i1) - p1(i3,i2,i1)
   if (diff(i3,i2,i1) .gt. 0.0d1) then
    pp(i3,i2,i1) = diff(i3,i2,i1)
    pm(i3,i2,i1) = 0.0d1
   else
    pm(i3,i2,i1) = diff(i3,i2,i1)
    pp(i3,i2,i1) = 0.0d1
   endif
  enddo
 enddo
enddo

! Starts the integration of the densities, the computation of the
! charge transfer integral, locates the densities centroids and
! computes the inter-centroid distance together with the phiS index.

allocate (xNI(sn3,sn2,sn1p))
allocate (yNI(sn3,sn2,sn1p))
allocate (zNI(sn3,sn2,sn1p))

! Initializes phiS and some intermediate quantities used for computing the integrals.

sphis = 0.0d1

sintegp = 0.0d1
sintegm = 0.0d1

sxrp = 0.0d1
syrp = 0.0d1
szrp = 0.0d1
sxrm = 0.0d1
syrm = 0.0d1
szrm = 0.0d1
sx1c = 0.0d1
sy1c = 0.0d1
sz1c = 0.0d1
sx2c = 0.0d1
sy2c = 0.0d1
sz2c = 0.0d1

do i1=1,sn1p
 do i2=1,sn2
  do i3=1,sn3

! Computes the integral of pm and pp, and computes their centroid.

xNI(i3,i2,i1)=sx0+(i1-1)*sdx 
yNI(i3,i2,i1)=sy0+(i2-1)*sdy
zNI(i3,i2,i1)=sz0+(i3-1)*sdz

      sintegp = sintegp + pp(i3,i2,i1)*sdx*sdy*sdz
      sintegm = sintegm + pm(i3,i2,i1)*sdx*sdy*sdz

      sxrp = sxrp + &
& (xNI(i3,i2,i1)*pp(i3,i2,i1)*sdx*sdy*sdz)
      sxrm = sxrm + &
& (xNI(i3,i2,i1)*pm(i3,i2,i1)*sdx*sdy*sdz)

      syrp = syrp + &
& (yNI(i3,i2,i1)*pp(i3,i2,i1)*sdx*sdy*sdz)
      syrm = syrm + &
& (yNI(i3,i2,i1)*pm(i3,i2,i1)*sdx*sdy*sdz)

      szrp = szrp + &
& (zNI(i3,i2,i1)*pp(i3,i2,i1)*sdx*sdy*sdz)
      szrm = szrm + &
& (zNI(i3,i2,i1)*pm(i3,i2,i1)*sdx*sdy*sdz)

    enddo 
  enddo
enddo

! Divides the expectation values of centroid positions
! by the pm and pp integrals.

sxrm = sxrm/sintegm
syrm = syrm/sintegm
szrm = szrm/sintegm
sxrp = sxrp/sintegp
syrp = syrp/sintegp
szrp = szrp/sintegp

! Computes the inter-centroid distance.

srct = dsqrt((sxrp - sxrm)**2 + (syrp - syrm)**2 &
& + (szrp - szrm)**2)

! Computes the detachment/attachment inter-centroid distance, and
! the detachment/attachment overlap integral.

if (ss1 .eq. 1) then

do i1 = 1,sn1p
  do i2 = 1,sn2
    do i3 = 1,sn3

      sx1c = sx1c + &
& (xNI(i3,i2,i1)*p1(i3,i2,i1)*sdx*sdy*sdz)
      sx2c = sx2c + &
& (xNI(i3,i2,i1)*p2(i3,i2,i1)*sdx*sdy*sdz)

      sy1c = sy1c + &
& (yNI(i3,i2,i1)*p1(i3,i2,i1)*sdx*sdy*sdz)
      sy2c = sy2c + &
& (yNI(i3,i2,i1)*p2(i3,i2,i1)*sdx*sdy*sdz)

      sz1c = sz1c + &
& (zNI(i3,i2,i1)*p1(i3,i2,i1)*sdx*sdy*sdz)
      sz2c = sz2c + &
& (zNI(i3,i2,i1)*p2(i3,i2,i1)*sdx*sdy*sdz)

sphis = sphis + &
& dsqrt(dabs(p1(i3,i2,i1))*dabs(p2(i3,i2,i1))) &
& *sdx*sdy*sdz/((sinteg1 + sinteg2)/2)

    enddo
  enddo
enddo

sx1c = sx1c/sinteg1
sy1c = sy1c/sinteg1
sz1c = sz1c/sinteg1

sx2c = sx2c/sinteg2
sy2c = sy2c/sinteg2
sz2c = sz2c/sinteg2

sr12ct = dsqrt((sx1c - sx2c)**2 + (sy1c - sy2c)**2 &
& + (sz1c - sz2c)**2)

! Deallocates the arrays.

deallocate(p1)
deallocate(p2)
deallocate(diff)
deallocate(pp)
deallocate(pm)
deallocate(xNI)
deallocate(yNI)
deallocate(zNI)

else

deallocate(p2)
deallocate(diff)
deallocate(pp)
deallocate(pm)
deallocate(xNI)
deallocate(yNI)
deallocate(zNI)

endif

return 

end 
