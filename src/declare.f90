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

module declare

implicit none

! General

integer :: error,i,j,k,nPA,nAdiab,jobN,countunr,countrlx,shell_statement,iteration
real*8 :: x,y,z,xi,xPA,yPA,trmat1,trmat2,PopX
character*128 :: PA_status,eigname

! Input reading

character*128 :: jobname,jobtype,qmPA_obj,subjobtype
character*128 :: soft,fchk,fov,fdens,input,fdensx
integer :: nt,ns
logical :: relax,adiab,PA,PAr,xy,scanPA,qmPA,t_PA,t_scanPA
character*128,allocatable :: dummy_array(:)

! Output
character*128 :: logfile

! Data related to the calculation

character*128 :: line
integer :: nea,neb,nel,nbs,norb,ntr
integer :: n,l,nL
logical :: unr

! LCAO matrices

real*8,allocatable :: triangle(:),trianglespin(:),lvec(:),triangled(:),trianglea(:)
real*8,allocatable :: C(:,:),Ca(:,:),Cb(:,:),Cdagger(:,:),Cadagger(:,:),Cbdagger(:,:)
character*128 :: marker_Ca, marker_Cb, marker_Ea, marker_Eb

! State density matrices

real*8,allocatable :: p(:,:),pK(:,:),pKS(:,:)
real*8,allocatable :: px(:,:),pxK(:,:),pxKS(:,:)
integer :: n0
real*8,allocatable :: pKspin(:,:),pKalpha(:,:),pKbeta(:,:)
real*8,allocatable :: pKalphaS(:,:),pKbetaS(:,:)
real*8,allocatable :: palpha(:,:),pbeta(:,:)
real*8,allocatable :: pxalpha(:,:),pxKalpha(:,:),pxbeta(:,:),pxKbeta(:,:)
real*8,allocatable :: pxKspin(:,:),pxKalphaS(:,:),pxKbetaS(:,:)
real*8,allocatable :: pxKspinRelaxed(:,:),pxKalphaRelaxed(:,:),pxKbetaRelaxed(:,:),pxalphaRelaxed(:,:),pxbetaRelaxed(:,:)
real*8,allocatable :: pxKalphaRelaxedS(:,:),pxKbetaRelaxedS(:,:)

! Overlap matrix

real*8,allocatable :: S(:,:)

! Temporary matrices

real*8,allocatable :: tp(:,:),tLK(:,:),tKL(:,:),tLL(:,:),tKK(:,:)
real*8,allocatable :: tKK1(:,:),tKK2(:,:)

! Difference density matrix, detachment/attachment

real*8,allocatable :: gD(:,:),gamma_d(:,:),gamma_a(:,:)
real*8,allocatable :: gamma_d1(:,:),gamma_a1(:,:),gamma_d_ao(:,:),gamma_a_ao(:,:)
real*8,allocatable :: triangle_gamma_d_ao(:),triangle_gamma_a_ao(:)
real*8,allocatable :: Mmg(:,:),Mmd(:,:)
real*8,allocatable :: SxDSy(:,:),SxASy(:,:),Uvec(:,:)
real*8 :: trD,trA,trDZ,trAZ,phiSPA,chimPA,chipPA,phiPA,psiPA,theta,theta0,thetaZ,alter_thetaZ
real*8 :: pi,phiSPA0,phiPA0
! Implemented by GB
real*8 :: phiSPAnew0,phiPAnew0
! End of GB implementation
real*8,allocatable :: fd(:),fa(:),fdta(:),famd(:),famdp(:),famdm(:)
real*8,allocatable :: gamma_d_XY_ao(:,:),gamma_a_XY_ao(:,:)

! Relaxed difference density matrix

real*8,allocatable :: pxKrelaxed(:,:),pxKrelaxedS(:,:),pxrelaxed(:,:)
real*8,allocatable :: zvec(:,:),zzd_zdz(:,:),zvecxi(:,:)
real*8,allocatable :: U0(:,:),U0t(:,:),U0tU(:,:)
real*8 :: Rcoef
! Implemented by GB
!real*8,allocatable :: newzzd_zdz(:,:)
! End of GB implementation

! Orbitals

real*8,allocatable :: tK(:,:),t(:,:),TS(:,:),diff_mat_unrelaxed(:,:)
real*8,allocatable :: tKb(:,:),tb(:,:),TSb(:,:)
real*8,allocatable :: ttrsp(:,:),ttdagger(:,:)
real*8,allocatable :: xy_X(:,:),xy_Y(:,:),xy_Xt(:,:),xy_Yt(:,:),XYt(:,:),YXt(:,:)
real*8,allocatable :: xy_XXt(:,:),xy_XtX(:,:),xy_YYt(:,:),xy_YtY(:,:)
real*8,allocatable :: t_0p(:,:),t_0m(:,:),t_02m(:,:),t_02p(:,:),t_0pt(:,:),t_0pt_0pt(:,:)
real*8,allocatable :: t_02mt(:,:),t_02mt_02mt(:,:)
real*8,allocatable :: t_1(:,:),t_2(:,:),t_3(:,:)
real*8 :: xy_norm,xy_residue,t3_norm,x2y2_norm
real*8,allocatable :: gamma_d_XY(:,:),gamma_a_XY(:,:)

character*128 :: orbs_filename

! SVD

real*8,allocatable :: lcao_coeff_mat(:,:)
real*8,allocatable :: lambda(:),left_eig(:,:),right_eig_t(:,:),mat_to_svd(:,:),right_eig(:,:)
real*8,allocatable :: O_lcao(:,:),V_lcao(:,:),rotated_O_lcao(:,:),rotated_V_lcao(:,:)
real*8,allocatable :: OdaggerSO(:,:),VdaggerSV(:,:)

! AO matrices

real*8,allocatable :: gD_ao(:,:),T_ao(:,:),zvec_ao(:,:)
real*8,allocatable :: triangle_gD_ao(:),triangle_T_ao(:),triangle_zvec_ao(:)

! Relaxed descriptors (rlxy_PA)

real*8 :: thetaUalpha, thetaUbeta
real*8 :: phiSPAUalpha, phiSPAUbeta
real*8 :: phiPAUalpha, phiPAUbeta
real*8 :: thetaZalpha, thetaZbeta

real*8 :: phiSPAU,phiPAU,psiPAU,thetaU,phiSPArlx,phiPArlx,psiPArlx
real*8 :: phiSPArlxalpha, phiSPArlxbeta
real*8 :: phiPArlxalpha, phiPArlxbeta
real*8 :: psiPAlxalpha, psiPArlxbeta

! Implemented by GB
real*8 :: zcoef,zcoefalpha,zcoefbeta
! Implemented by GB


! Numerical Integration

character*128 :: file1,file2
character*128 :: input1,input2,input3,filename
integer :: natom,n1,n2,n3,i1,i2,i3
integer :: natom1,n11,n21,n31,nci
real*8 :: dx,dy,dz,int1,integ1,integ2,phis,dummy,dummy1
real*8 :: integp, integm, integrp, integrm
real*8 ::dx1,dy1,dz1,THR,occ,virt
real*8 :: integ1r,integ2r,x1c,y1c,z1c,x2c,y2c,z2c,r12ct
real*8 :: xrp, yrp, zrp, xrm, yrm, zrm, rct
real*8, allocatable :: xNI(:,:,:), yNI(:,:,:), zNI(:,:,:)
real*8, allocatable :: nto1(:,:,:), nto2(:,:,:)
real*8 :: x0,y0,z0,x01,y01,z01
real*8, allocatable :: p1(:,:,:),p2(:,:,:),r(:,:,:)
real*8, allocatable :: pp(:,:,:),pm(:,:,:),diff(:,:,:)

! Split and cubeop

character*128 :: cubefile,cubefile1,cubefile2,cubefile3,op1,op2

! dagger

character*128 :: detachmentU,attachmentU,detachmentR,attachmentR
character*128 :: detachmentUalpha,attachmentUalpha,detachmentRalpha,attachmentRalpha
character*128 :: detachmentUbeta,attachmentUbeta,detachmentRbeta,attachmentRbeta

! qmNIrlx

real*8 :: theta_unrelaxed,phiS_unrelaxed,phi_unrelaxed,lambda_dagger,theta_z
real*8 :: phiS_relaxed,phi_relaxed,psi_relaxed,eta
real*8 :: phiS_relaxedPA0,phi_relaxedPA0
real*8 :: phiS_relaxedPA,phi_relaxedPA,psi_relaxedPA

end
