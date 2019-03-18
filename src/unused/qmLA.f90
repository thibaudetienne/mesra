subroutine qmLA

use declare

allocate(fd(nbs))
allocate(fa(nbs))
allocate(fdta(nbs))
allocate(famd(nbs))
allocate(famdm(nbs))
allocate(famdp(nbs))

fd       = 0.0d0 
fa = 0.0d0
fdta = 0.0d0
famd = 0.0d0
famdm = 0.0d0
famdp = 0.0d0
chimLA = 0.0d0
chipLA = 0.0d0
phiSLA = 0.0d0
phiLA = 0.0d0
psiLA = 0.0d0
theta = 0.0d0
pi = 0.0d0

theta = 0.5d0*(trD + trA)

do i=1,nbs
fd(i) = SxDSy(i,i)
fa(i) = SxASy(i,i)
fdta(i) = fd(i)*fa(i)
if (fd(i) .gt. 0.0d0 .and. fa(i) .gt. 0.0d0) phiSLA = phiSLA + (dsqrt(fdta(i)))/theta
famd(i) = fa(i) - fd(i)
if (famd(i) .lt. 0.0d0) then
famdm(i) = -1.0d0*famd(i)
else
famdp(i) = famd(i)
endif
chimLA = chimLA + famdm(i)
chipLA = chipLA + famdp(i)
enddo

phiLA = 0.5d0*(chimLA + chipLA)/theta

pi = dacos(-1.0d0)

psiLA = 2*(datan(phiSLA/phiLA))/(pi)

write(6,*) 'phiSLA,phiLA,psiLA ',phiSLA,phiLA,psiLA

deallocate(fd,fa,fdta,famd,famdm,famdp)

end
