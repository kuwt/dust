!./\\\\\\\\\\\...../\\\......./\\\..../\\\\\\\\\..../\\\\\\\\\\\\\. 
!.\/\\\///////\\\..\/\\\......\/\\\../\\\///////\\\.\//////\\\////..
!..\/\\\.....\//\\\.\/\\\......\/\\\.\//\\\....\///.......\/\\\......
!...\/\\\......\/\\\.\/\\\......\/\\\..\////\\.............\/\\\......
!....\/\\\......\/\\\.\/\\\......\/\\\.....\///\\...........\/\\\......
!.....\/\\\......\/\\\.\/\\\......\/\\\.......\///\\\........\/\\\......
!......\/\\\....../\\\..\//\\\...../\\\../\\\....\//\\\.......\/\\\......
!.......\/\\\\\\\\\\\/....\///\\\\\\\\/..\///\\\\\\\\\/........\/\\\......
!........\///////////........\////////......\/////////..........\///.......
!!=========================================================================
!!
!! Copyright (C) 2018-2019 Davide   Montagnani, 
!!                         Matteo   Tugnoli, 
!!                         Federico Fonte
!!
!! This file is part of DUST, an aerodynamic solver for complex
!! configurations.
!! 
!! Permission is hereby granted, free of charge, to any person
!! obtaining a copy of this software and associated documentation
!! files (the "Software"), to deal in the Software without
!! restriction, including without limitation the rights to use,
!! copy, modify, merge, publish, distribute, sublicense, and/or sell
!! copies of the Software, and to permit persons to whom the
!! Software is furnished to do so, subject to the following
!! conditions:
!! 
!! The above copyright notice and this permission notice shall be
!! included in all copies or substantial portions of the Software.
!! 
!! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
!! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
!! OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
!! NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
!! HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
!! WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
!! FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
!! OTHER DEALINGS IN THE SOFTWARE.
!! 
!! Authors: 
!!          Federico Fonte             <federico.fonte@outlook.com>
!!          Davide Montagnani       <davide.montagnani@gmail.com>
!!          Matteo Tugnoli                <tugnoli.teo@gmail.com>
!!=========================================================================

module mod_vtk_utils

use mod_param, only: &
  wp, max_char_len, nl

use mod_handling, only: &
  error, internal_error, warning, info, printout, new_file_unit


!---------------------------------------------------------------------
implicit none

!public ::

private

type :: t_output_var
  character(max_char_len) :: var_name
  logical :: vector
  real(wp), allocatable :: var(:,:)
end type t_output_var

interface add_output_var
  module procedure add_output_var_scal
  module procedure add_output_var_vec
end interface add_output_var

integer, parameter :: vtk_isize = 4
integer, parameter :: vtk_fsize = 4
character(len=*), parameter :: &
  this_mod_name = 'mod_vtk_utils'
!---------------------------------------------------------------------

contains

!---------------------------------------------------------------------

subroutine add_output_var_scal(output_vars, var, var_name)
 type(t_output_var), intent(inout), allocatable :: output_vars(:)
 real(wp), intent(in) :: var(:)
 character(len=*), intent(in) :: var_name


end subroutine add_outpu_var_scal

!---------------------------------------------------------------------

subroutine extend_vars(output_vars)
 type(t_output_var), intent(inout), allocatable :: output_vars(:)

 type(t_output_var), allocatable :: vars_tmp(:)
 integer :: varlen,i
 character(len=*), parameter :: this_sub_name='extend_vars'

  if(.not. allocated(output_vars)) &
    call error(this_sub_name, this_mod_name, 'Internal error: not &
    &allocated output_vars. This should have never happened, please &
    &report this error so that  a team of professionals could travel &
    &back in time and fix this')

end subroutine extend_vars

!---------------------------------------------------------------------

subroutine vtk_print_piece_header(offset, npoints, ncells, nquad, ntria)
 integer, intent(inout) :: offset
 integer, intent(in) :: npoints
 integer, intent(in) :: ncells

 character(len=200) :: buffer
 character(len=20) :: ostr, str1, str2

  write(str1,'(I0)') npoints
  write(str2,'(I0)') ncells
  buffer = '  <Piece NumberOfPoints="'//trim(str1)//'" NumberOfCells="'//&
                   trim(str2)//'">'//lf; write(fu) trim(buffer)
  !Points
  buffer =  '   <Points>'//lf; write(fu) trim(buffer)
  write(ostr,'(I0)') offset
  buffer='    <DataArray type="Float32" Name="Coordinates" &
                  &NumberOfComponents="3" format="appended" &
         &offset="'//trim(ostr)//'"/>'//lf;write(fu) trim(Buffer)

  buffer =  '   </Points>'//lf; write(fu) trim(buffer)
  offset = offset + vtk_isize + vtk_fsize*size(rr,1)*size(rr,2)

  !Cells
  buffer =  '   <Cells>'//lf; write(fu) trim(buffer)

  write(ostr,'(I0)') offset
  buffer = '    <DataArray type="Int32" Name="connectivity" &
           &Format="appended" offset="'//trim(ostr)//'"/>'//lf; 
  write(fu) trim(buffer)
  offset = offset + vtk_isize + vtk_isize*3*ntria + vtk_isize*4*nquad

  write(ostr,'(I0)') offset
  buffer =  '    <DataArray type="Int32" Name="offsets" &
    &Format="appended" offset="'//trim(ostr)//'"/>'//lf;
  write(fu) trim(buffer)
  offset = offset + vtk_isize + vtk_isize*ncells

  write(ostr,'(I0)') offset
  buffer = '    <DataArray type="Int32" Name="types" &
    &Format="appended" offset="'//trim(ostr)//'"/>'//lf;
  write(fu) trim(buffer)
  offset = offset + vtk_isize + vtk_isize*ncells

  buffer =  '   </Cells>'//lf; write(fu) trim(buffer)

  !Data
  if (nvars .gt. 0) then
    buffer =  '   <CellData Scalars="scalars">'//lf; 
    write(fu) trim(buffer)
    do i_v = 1,nvars
      write(ostr,'(I0)') offset
      buffer = '    <DataArray type="Float32" Name="'//trim(var_names(i_v))//'" &
        &Format="appended" offset="'//trim(ostr)//'"/>'//lf
      write(fu) trim(buffer)
      offset = offset + vtk_isize + vtk_fsize*ne
    enddo
    buffer = '   </CellData>'//lf; write(fu) trim(buffer)
  endif

  buffer = '  </Piece>'//lf; write(fu) trim(buffer)

end subroutine vtk_print_piece_header

end module mod_vtk_utils

