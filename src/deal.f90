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

subroutine deal

use declare

! Though many of these matrices are either not called during the calculations, or are already
! deallocated, we check all of them, and deallocate them if necessary.

if (allocated(dummy_array)) deallocate(dummy_array)
if (allocated(triangle)) deallocate(triangle)
if (allocated(trianglespin)) deallocate(trianglespin)
if (allocated(lvec)) deallocate(lvec)
if (allocated(triangled)) deallocate(triangled)
if (allocated(trianglea)) deallocate(trianglea)
if (allocated(c)) deallocate(c)
if (allocated(ca)) deallocate(ca)
if (allocated(cb)) deallocate(cb)
if (allocated(cdagger)) deallocate(cdagger)
if (allocated(cadagger)) deallocate(cadagger)
if (allocated(cbdagger)) deallocate(cbdagger)
if (allocated(p)) deallocate(p)
if (allocated(pk)) deallocate(pk)
if (allocated(pks)) deallocate(pks)
if (allocated(px)) deallocate(px)
if (allocated(pxk)) deallocate(pxk)
if (allocated(pxks)) deallocate(pxks)
if (allocated(pkspin)) deallocate(pkspin)
if (allocated(pkalpha)) deallocate(pkalpha)
if (allocated(pkbeta)) deallocate(pkbeta)
if (allocated(pkalphas)) deallocate(pkalphas)
if (allocated(pkbetas)) deallocate(pkbetas)
if (allocated(palpha)) deallocate(palpha)
if (allocated(pbeta)) deallocate(pbeta)
if (allocated(pxalpha)) deallocate(pxalpha)
if (allocated(pxkalpha)) deallocate(pxkalpha)
if (allocated(pxbeta)) deallocate(pxbeta)
if (allocated(pxkbeta)) deallocate(pxkbeta)
if (allocated(pxkspin)) deallocate(pxkspin)
if (allocated(pxkalphas)) deallocate(pxkalphas)
if (allocated(pxkbetas)) deallocate(pxkbetas)
if (allocated(pxkspinrelaxed)) deallocate(pxkspinrelaxed)
if (allocated(pxkalpharelaxed)) deallocate(pxkalpharelaxed)
if (allocated(pxkbetarelaxed)) deallocate(pxkbetarelaxed)
if (allocated(pxalpharelaxed)) deallocate(pxalpharelaxed)
if (allocated(pxbetarelaxed)) deallocate(pxbetarelaxed)
if (allocated(pxkalpharelaxeds)) deallocate(pxkalpharelaxeds)
if (allocated(pxkbetarelaxeds)) deallocate(pxkbetarelaxeds)
if (allocated(s)) deallocate(s)
if (allocated(tp)) deallocate(tp)
if (allocated(tlk)) deallocate(tlk)
if (allocated(tkl)) deallocate(tkl)
if (allocated(tll)) deallocate(tll)
if (allocated(tkk)) deallocate(tkk)
if (allocated(tkk1)) deallocate(tkk1)
if (allocated(tkk2)) deallocate(tkk2)
if (allocated(gd)) deallocate(gd)
if (allocated(gamma_d)) deallocate(gamma_d)
if (allocated(gamma_a)) deallocate(gamma_a)
if (allocated(gamma_d1)) deallocate(gamma_d1)
if (allocated(gamma_a1)) deallocate(gamma_a1)
if (allocated(gamma_d_ao)) deallocate(gamma_d_ao)
if (allocated(gamma_a_ao)) deallocate(gamma_a_ao)
if (allocated(triangle_gamma_d_ao)) deallocate(triangle_gamma_d_ao)
if (allocated(triangle_gamma_a_ao)) deallocate(triangle_gamma_a_ao)
if (allocated(mmg)) deallocate(mmg)
if (allocated(mmd)) deallocate(mmd)
if (allocated(sxdsy)) deallocate(sxdsy)
if (allocated(sxasy)) deallocate(sxasy)
if (allocated(uvec)) deallocate(uvec)
if (allocated(fd)) deallocate(fd)
if (allocated(fa)) deallocate(fa)
if (allocated(fdta)) deallocate(fdta)
if (allocated(famd)) deallocate(famd)
if (allocated(famdp)) deallocate(famdp)
if (allocated(gamma_d_xy_ao)) deallocate(gamma_d_xy_ao)
if (allocated(gamma_a_xy_ao)) deallocate(gamma_a_xy_ao)
if (allocated(pxkrelaxed)) deallocate(pxkrelaxed)
if (allocated(pxkrelaxeds)) deallocate(pxkrelaxeds)
if (allocated(pxrelaxed)) deallocate(pxrelaxed)
if (allocated(zvec)) deallocate(zvec)
if (allocated(zzd_zdz)) deallocate(zzd_zdz)
if (allocated(zvecxi)) deallocate(zvecxi)
if (allocated(u0)) deallocate(u0)
if (allocated(u0t)) deallocate(u0t)
if (allocated(u0tu)) deallocate(u0tu)
if (allocated(tk)) deallocate(tk)
if (allocated(t)) deallocate(t)
if (allocated(ts)) deallocate(ts)
if (allocated(diff_mat_unrelaxed)) deallocate(diff_mat_unrelaxed)
if (allocated(tkb)) deallocate(tkb)
if (allocated(tb)) deallocate(tb)
if (allocated(tsb)) deallocate(tsb)
if (allocated(ttrsp)) deallocate(ttrsp)
if (allocated(ttdagger)) deallocate(ttdagger)
if (allocated(xy_X)) deallocate(xy_X)
if (allocated(xy_Y)) deallocate(xy_Y)
if (allocated(xy_Xt)) deallocate(xy_Xt)
if (allocated(xy_yt)) deallocate(xy_yt)
if (allocated(xyt)) deallocate(xyt)
if (allocated(yxt)) deallocate(yxt)
if (allocated(xy_xxt)) deallocate(xy_xxt)
if (allocated(xy_xtx)) deallocate(xy_xtx)
if (allocated(xy_yyt)) deallocate(xy_yyt)
if (allocated(xy_yty)) deallocate(xy_yty)
if (allocated(t_0p)) deallocate(t_0p)
if (allocated(t_0m)) deallocate(t_0m)
if (allocated(t_02m)) deallocate(t_02m)
if (allocated(t_02p)) deallocate(t_02p)
if (allocated(t_0pt)) deallocate(t_0pt)
if (allocated(t_0pt_0pt)) deallocate(t_0pt_0pt)
if (allocated(t_02mt)) deallocate(t_02mt)
if (allocated(t_02mt_02mt)) deallocate(t_02mt_02mt)
if (allocated(t_1)) deallocate(t_1)
if (allocated(t_2)) deallocate(t_2)
if (allocated(t_3)) deallocate(t_3)
if (allocated(gamma_d_xy)) deallocate(gamma_d_xy)
if (allocated(gamma_a_xy)) deallocate(gamma_a_xy)
if (allocated(lcao_coeff_mat)) deallocate(lcao_coeff_mat)
if (allocated(lambda)) deallocate(lambda)
if (allocated(left_eig)) deallocate(left_eig)
if (allocated(right_eig_t)) deallocate(right_eig_t)
if (allocated(mat_to_svd)) deallocate(mat_to_svd)
if (allocated(right_eig)) deallocate(right_eig)
if (allocated(o_lcao)) deallocate(o_lcao)
if (allocated(v_lcao)) deallocate(v_lcao)
if (allocated(rotated_o_lcao)) deallocate(rotated_o_lcao)
if (allocated(rotated_v_lcao)) deallocate(rotated_v_lcao)
if (allocated(odaggerso)) deallocate(odaggerso)
if (allocated(vdaggersv)) deallocate(vdaggersv)
if (allocated(gd_ao)) deallocate(gd_ao)
if (allocated(t_ao)) deallocate(t_ao)
if (allocated(zvec_ao)) deallocate(zvec_ao)
if (allocated(triangle_gd_ao)) deallocate(triangle_gd_ao)
if (allocated(triangle_t_ao)) deallocate(triangle_t_ao)
if (allocated(triangle_zvec_ao)) deallocate(triangle_zvec_ao)
if (allocated(xni)) deallocate(xni)
if (allocated(yni)) deallocate(yni)
if (allocated(zni)) deallocate(zni)
if (allocated(nto1)) deallocate(nto1)
if (allocated(nto2)) deallocate(nto2)
if (allocated(p1)) deallocate(p1)
if (allocated(p2)) deallocate(p2)
if (allocated(r)) deallocate(r)
if (allocated(pp)) deallocate(pp)
if (allocated(pm)) deallocate(pm)
if (allocated(diff)) deallocate(diff)

end
