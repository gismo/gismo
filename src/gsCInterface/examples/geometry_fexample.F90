!** @file f_xml_io.c
!
!   @brief Tests Fortran interface
!
!   This file is part of the G+Smo library.
!
!   This Source Code Form is subject to the terms of the Mozilla Public
!   License, v. 2.0. If a copy of the MPL was not distributed with this
!   file, You can obtain one at http://mozilla.org/MPL/2.0/.
!
!   Author(s): E. Vollebregt
! 

program geometry_fexample
   use, intrinsic  :: iso_c_binding
   use Fgismo
   implicit none
   character(len=80, kind=C_CHAR) :: some_file
   type(t_gsgeometry) :: g
   integer            :: nu,  nv

   ! TODO: input option for xml-file?

   write(*,'(2(a,f5.1),a,i3)') 'reading XML for tensor B-spline'
   some_file = 'sw_tp.xml' // C_NULL_CHAR
   call f_gsCReadFile(some_file, g)

   write(*,'(a,i3)') 'done, g.dim=', f_gsFunctionSet_domainDim(g)

   call f_gsFunctionSet_print(g)

   if (.true.) then
      call show_basic_usage( g )
   endif

   if (.true.) then
      call show_recover_points( g )
   endif

   call f_gsfunctionset_delete(g)
   write(*,*) 'done.'

end program geometry_fexample

!-----------------------------------------------------------------------------------------------------------

subroutine show_basic_usage( g )
!--purpose: evaluate positions (x,y,z) and derivatives d[xyz]/d[uv] at some arbitrary (u,v) \in [0,1)^2
   use, intrinsic  :: iso_c_binding
   use Fgismo
   implicit none
!--subroutine arguments
   type(t_gsgeometry)           :: g
!--local variables
   integer(C_INT)               :: nRows, nCols, out_rows, out_cols, irow, icol, icoor, ipar
   type(t_gsmatrix)             :: uvm, xyzm, dxyzm
   real(C_DOUBLE), dimension(:,:), allocatable :: uv
   character(len=1), parameter  :: c_param(2)  = (/ 'u', 'v' /)
   character(len=1), parameter  :: c_coor(3)   = (/ 'x', 'y', 'z' /)
   character(len=5), parameter  :: c_deriv(6)  = (/ 'dx/du', 'dx/dv', 'dy/du', 'dy/dv', 'dz/du', 'dz/dv' /)

   nRows = 2
   nCols = 7
   allocate(uv(nrows,ncols))
   uv(1, 1:nCols) = (/ 0.0, 0.1, 0.501, 0.500, 0.500, 0.5, 1.0 /)
   uv(2, 1:nCols) = (/ 0.0, 0.0, 0.200, 0.200, 0.201, 0.9, 1.0 /)

   write(*,*) '------------------------------ show_basic_usage ------------------------------'
   write(*,'(2(a,i3))') 'Input #rows =', nRows, ', #cols =', nCols
   do irow = 1, nRows
      write(*,'(3a,10f10.3)') '      ',c_param(irow),': ', (uv(irow,icol), icol=1,nCols)
   enddo

   ! evaluate positions (x,y,z) at given parameter values

   uvm   = f_gsmatrix_create_rcd(nRows, nCols, uv)
   xyzm  = f_gsmatrix_create()
   dxyzm = f_gsmatrix_create()
   call f_gsFunctionSet_eval_into(G, uvm, xyzm)
   call f_gsFunctionSet_deriv_into(G, uvm, dxyzm)
   ! call f_gsmatrix_print(xyzm)

   ! show output data

   out_rows = f_gsmatrix_rows(xyzm)
   out_cols = f_gsmatrix_cols(xyzm)

   write(*,'(3(a,i3))') 'Values: #rows =', out_rows, ', #cols =', out_cols
   do irow = 1, out_rows
      write(*,'(3a,10f10.3)') '      ',c_coor(irow),': ', (xyzm%data(irow,icol), icol=1,out_cols)
   enddo

   ! show derivatives data

   out_rows = f_gsmatrix_rows(dxyzm)
   out_cols = f_gsmatrix_cols(dxyzm)

   write(*,'(3(a,i3))') 'Derivatives: #rows =', out_rows, ', #cols =', out_cols
   do irow = 1, out_rows
      write(*,'(3a,10f10.3)') '  ',c_deriv(irow),': ', (dxyzm%data(irow,icol), icol=1,out_cols)
   enddo

   call f_gsmatrix_delete(uvm)
   call f_gsmatrix_delete(xyzm)
   call f_gsmatrix_delete(dxyzm)
   deallocate(uv)

end subroutine show_basic_usage

!-----------------------------------------------------------------------------------------------------------

subroutine show_recover_points( g )
!--purpose: for some positions (x,y), determine z on the surface and corresponding (u,v)
   use, intrinsic  :: iso_c_binding
   use Fgismo
   implicit none
!--subroutine arguments
   type(t_gsgeometry)           :: g
!--local variables
   integer(C_INT), parameter    :: XDIR = 0, YDIR = 1, ZDIR = 2
   integer(C_INT)               :: nCols, irow, icol, out_rows, out_cols
   real(C_DOUBLE)               :: eps
   real(C_DOUBLE), dimension(:,:), allocatable :: xyz
   type(t_gsmatrix)             :: uvm, xyzm
   character(len=1), parameter  :: c_param(2) = (/ 'u', 'v' /)
   character(len=1), parameter  :: c_coor(3)  = (/ 'x', 'y', 'z' /)

   write(*,*) '----------------------------- show_recover_points -----------------------------'
   nCols = 4
   allocate(xyz(3,ncols))

   xyz(1, 1:ncols) = (/ 2451.0, 210.001, 708.0,   210.0 /)
   xyz(2, 1:ncols) = (/  122.3, -38.957,  18.568, -13.9 /)
   xyz(3, 1:ncols) = (/    0.0,   0.0,     0.0,     0.0 /)

   write(*,'(a,i3,a)') 'Input #cols =', nCols,', (x,y,z) =' 
   do irow = 1, 3
      write(*,'(3a,10f10.3)') '  ',c_coor(irow),': ', (xyz(irow,icol), icol=1,nCols)
   enddo

   ! evaluate z-positions and (u,v)-values for given (x,y)-positions

   ! evaluate positions (x,y,z) at given parameter values

   xyzm = f_gsmatrix_create_rcd(3, ncols, xyz)
   uvm  = f_gsmatrix_create()

   eps = 1d-6
   call f_gsGeometry_recoverPoints(G, uvm, xyzm, ZDIR, eps)

   ! print output data

   out_rows = f_gsmatrix_rows(uvm)
   out_cols = f_gsmatrix_cols(uvm)

   write(*,'(a)') 'Output (u,v) =' 
   do irow = 1, 2
      write(*,'(3a,10f10.3)') '  ',c_param(irow),': ', (uvm%data(irow,icol), icol=1, nCols)
   enddo
   write(*,'(a)') 'Output (x,y,z) =' 
   do irow = 1, 3
      write(*,'(3a,10f10.3)') '  ',c_coor(irow),': ', (xyzm%data(irow,icol), icol=1, nCols)
   enddo

   ! clean up input data, matrices used

   call f_gsmatrix_delete(xyzm)
   call f_gsmatrix_delete(uvm)
   deallocate(xyz)

end subroutine show_recover_points

!-----------------------------------------------------------------------------------------------------------
