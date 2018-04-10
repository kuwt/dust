!!=====================================================================
!!
!! Copyright (C) 2018 Politecnico di Milano
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
!!          Federico Fonte             <federico.fonte@polimi.it>
!!          Davide Montagnani       <davide.montagnani@polimi.it>
!!          Matteo Tugnoli             <matteo.tugnoli@polimi.it>
!!=====================================================================

module mod_vtk_out

use mod_param, only: &
  wp, max_char_len, nl

use mod_handling, only: &
  error, warning, info, printout, new_file_unit


!---------------------------------------------------------------------
implicit none

public :: vtk_out_bin

private

integer, parameter :: vtk_isize = 4
integer, parameter :: vtk_fsize = 4
character(len=*), parameter :: &
  this_mod_name = 'mod_vtk_output'
!---------------------------------------------------------------------

contains

!---------------------------------------------------------------------

!> Output the processed data into a binary xml  .vtu file
!!
!! This is a very preliminary version, it is set to take stupid ee and rr
!! vectors, padded with zeros
subroutine vtk_out_bin (rr, ee, vort, w_rr, w_ee, w_vort, w_vel, out_filename)
 real(wp), intent(in) :: rr(:,:)
 integer, intent(in) :: ee(:,:)
 real(wp), intent(in) :: vort(:)
 real(wp), intent(in) :: w_rr(:,:)
 integer, intent(in) :: w_ee(:,:)
 real(wp), intent(in) :: w_vort(:)
 real(wp), intent(in), optional :: w_vel(:,:)
 character(len=*), intent(in) :: out_filename 

 integer :: fu, ierr, i, j, i_shift, ivar
 integer :: npoints, ncells, ne
 integer :: offset, nbytes
 character(len=200) :: buffer
 character(len=20) :: ostr, str1, str2
 character(len=1)  :: lf

 integer :: ie, nquad, ntria, etype
 integer npoints_w, nw

  lf = char(10) !line feed char
  ne = size(ee,2) !number of elements
  nw = size(w_ee,2) !number of wake elements
  npoints = size(rr,2)
  npoints_w = size(w_rr,2)
  ncells = ne

  ! First cycle the elements to get the number of quads and trias
  ! (this is really ugly, make it more flexible)
  nquad = 0; ntria = 0
  do ie = 1,ne
    if(ee(4,ie) .eq. 0) then
      ntria = ntria+1
    else
      nquad = nquad+1
    endif
  enddo



  call new_file_unit(fu,ierr)
  open(fu,file=trim(out_filename), &
        status='replace',access='stream',iostat=ierr)
  buffer = '<?xml version="1.0"?>'//lf; write(fu) trim(buffer)
  buffer = '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'//lf
  write(fu) trim(buffer)
  buffer = ' <UnstructuredGrid>'//lf; write(fu) trim(buffer)
  write(str1,'(I0)') npoints
  write(str2,'(I0)') ncells
  buffer = '  <Piece NumberOfPoints="'//trim(str1)//'" NumberOfCells="'//&
                   trim(str2)//'">'//lf; write(fu) trim(buffer)
  offset = 0
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
  offset = offset + vtk_isize + vtk_isize*ne

  write(ostr,'(I0)') offset
  buffer = '    <DataArray type="Int32" Name="types" &
    &Format="appended" offset="'//trim(ostr)//'"/>'//lf;
  write(fu) trim(buffer)
  offset = offset + vtk_isize + vtk_isize*ne

  buffer =  '   </Cells>'//lf; write(fu) trim(buffer)

  !Data
  buffer =  '   <CellData Scalars="scalars">'//lf; 
  write(fu) trim(buffer)
  write(ostr,'(I0)') offset
  buffer = '    <DataArray type="Float32" Name="Vort_intensity" &
    &Format="appended" offset="'//trim(ostr)//'"/>'//lf
  write(fu) trim(buffer)
  offset = offset + vtk_isize + vtk_fsize*ne

  buffer = '   </CellData>'//lf; write(fu) trim(buffer)

  !!Data
  !buffer =  '   <PointData Scalars="scalars">'//lf; 
  !write(fu) trim(buffer)
  !do ivar=1,pv_data%n_vars
  !  if(pv_data%vars(ivar)%output) then
  !    if(pv_data%vars(ivar)%vector) then
  !      write(ostr,'(I0)') offset
  !      buffer =  '    <DataArray type="Float32" Name="'// &
  !        trim(pv_data%vars(ivar)%v_name) // &
  !        '" NumberOfComponents="3" Format="appended" offset="'&
  !                                              //trim(ostr)//'"/>'//lf
  !      write(fu) trim(buffer)
  !      offset = offset + vtk_isize + vtk_fsize*3*n_v4e*ne
  !    else !not vector
  !      write(ostr,'(I0)') offset
  !      buffer = '    <DataArray type="Float32" Name="'// &
  !        trim(pv_data%vars(ivar)%v_name) //'" Format="appended" &
  !        &offset="'//trim(ostr)//'"/>'//lf
  !      write(fu) trim(buffer)
  !      offset = offset + vtk_isize + vtk_fsize*n_v4e*ne
  !    endif !vector
  !  endif !output
  !enddo
  !buffer = '   </PointData>'//lf; write(fu) trim(buffer)

  buffer = '  </Piece>'//lf; write(fu) trim(buffer)




  !WAKE ===
  write(str1,'(I0)') npoints_w
  write(str2,'(I0)') nw
  buffer = '  <Piece NumberOfPoints="'//trim(str1)//'" NumberOfCells="'//&
                   trim(str2)//'">'//lf; write(fu) trim(buffer)
  !Points
  buffer =  '   <Points>'//lf; write(fu) trim(buffer)
  write(ostr,'(I0)') offset
  buffer='    <DataArray type="Float32" Name="Coordinates" &
                  &NumberOfComponents="3" format="appended" &
         &offset="'//trim(ostr)//'"/>'//lf;write(fu) trim(Buffer)

  buffer =  '   </Points>'//lf; write(fu) trim(buffer)
  offset = offset + vtk_isize + vtk_fsize*size(w_rr,1)*size(w_rr,2)

  !Cells
  buffer =  '   <Cells>'//lf; write(fu) trim(buffer)

  write(ostr,'(I0)') offset
  buffer = '    <DataArray type="Int32" Name="connectivity" &
           &Format="appended" offset="'//trim(ostr)//'"/>'//lf; 
  write(fu) trim(buffer)
  offset = offset + vtk_isize  + vtk_isize*4*nw

  write(ostr,'(I0)') offset
  buffer =  '    <DataArray type="Int32" Name="offsets" &
    &Format="appended" offset="'//trim(ostr)//'"/>'//lf;
  write(fu) trim(buffer)
  offset = offset + vtk_isize + vtk_isize*nw

  write(ostr,'(I0)') offset
  buffer = '    <DataArray type="Int32" Name="types" &
    &Format="appended" offset="'//trim(ostr)//'"/>'//lf;
  write(fu) trim(buffer)
  offset = offset + vtk_isize + vtk_isize*nw

  buffer =  '   </Cells>'//lf; write(fu) trim(buffer)

  !Data
  buffer =  '   <CellData Scalars="scalars">'//lf; 
  write(fu) trim(buffer)
  write(ostr,'(I0)') offset
  buffer = '    <DataArray type="Float32" Name="Vort_intensity" &
    &Format="appended" offset="'//trim(ostr)//'"/>'//lf
  write(fu) trim(buffer)
  offset = offset + vtk_isize + vtk_fsize*nw

  buffer = '   </CellData>'//lf; write(fu) trim(buffer)

  if(present(w_vel)) then
    buffer =  '   <PointData Vectors="wake_vel">'//lf; 
    write(fu) trim(buffer)
    write(ostr,'(I0)') offset
    buffer =  '    <DataArray type="Float32" Name="'// &
            'wake_vel' // &
            '" NumberOfComponents="3" Format="appended" offset="'&
                                                  //trim(ostr)//'"/>'//lf
    write(fu) trim(buffer)
    offset = offset + vtk_isize + vtk_fsize*3*nw
    buffer = '   </PointData>'//lf; write(fu) trim(buffer)
  endif


  buffer = '  </Piece>'//lf; write(fu) trim(buffer)




  buffer = ' </UnstructuredGrid>'//lf; write(fu) trim(buffer)

  !All the appended data
  buffer = '  <AppendedData encoding="raw">'//lf; write(fu) trim(buffer)
  buffer = '_'; write(fu) trim(buffer) !mark the beginning of the data 

  !Points
  nbytes = vtk_fsize *  &
               size(rr,1)*size(rr,2)
  write(fu) nbytes
  do i=1,size(rr,2)
   write(fu) real(rr(:,i),vtk_fsize)
  enddo

  !Connectivity
  nbytes = vtk_isize*3*ntria+vtk_isize*4*nquad; write(fu) nbytes
  do i=1,size(ee,2)
    if (ee(4,i) .eq. 0) then
      write(fu) ee(1:3,i)-1
    else
      write(fu) ee(1:4,i)-1
    endif
  enddo

  !Offset
  nbytes =  vtk_isize*ne; write(fu) nbytes
  i_shift = 0
  do i=1,size(ee,2)
    if (ee(4,i) .eq. 0) then
      i_shift = i_shift + 3
    else
      i_shift = i_shift + 4
    endif
    write(fu) i_shift
  enddo

  !Cell types
  nbytes = vtk_isize*ne; write(fu) nbytes
  do i=1,size(ee,2)
    if (ee(4,i) .eq. 0) then
      etype = 5
    else
      etype = 9
    endif
    write(fu) etype
  enddo

  !Doublet data
  nbytes =  vtk_fsize*ne; write(fu) nbytes
  do i=1,size(vort,1)
    write(fu) real(vort(i), vtk_fsize)
  enddo

!=== Wake
  !Points
  nbytes = vtk_fsize *  &
               size(w_rr,1)*size(w_rr,2)
  write(fu) nbytes
  do i=1,size(w_rr,2)
   write(fu) real(w_rr(:,i),vtk_fsize)
  enddo

  !Connectivity
  nbytes = vtk_isize*4*nw; write(fu) nbytes
  do i=1,size(w_ee,2)
    write(fu) w_ee(1:4,i)-1
  enddo

  !Offset
  nbytes =  vtk_isize*nw; write(fu) nbytes
  i_shift = 0
  do i=1,size(w_ee,2)
    i_shift = i_shift + 4
    write(fu) i_shift
  enddo

  !Cell types
  nbytes = vtk_isize*nw; write(fu) nbytes
  do i=1,size(w_ee,2)
    etype = 9
    write(fu) etype
  enddo

  !Doublet data
  nbytes =  vtk_fsize*nw; write(fu) nbytes
  do i=1,size(w_vort,1)
    write(fu) real(w_vort(i), vtk_fsize)
  enddo
  
  !Wake points velocity
  if(present(w_vel)) then
    nbytes =  vtk_fsize*3*nw; write(fu) nbytes
    do i=1,size(w_vort,1)
      write(fu) real(w_vel(:,i), vtk_fsize)
    enddo
  endif
  
  !!All the variables data
  !do ivar=1,pv_data%n_vars
  !  if(pv_data%vars(ivar)%output) then
  !    if(pv_data%vars(ivar)%vector) then
  !      nbytes =  vtk_fsize*3*n_v4e*ne; write(fu) nbytes
  !      do j=1,size(pv_data%Xi,3)
  !        do i=1,size(pv_data%Xi,2)
  !          write(fu) real(pv_data%vars(ivar)%n_var(:,i,j), vtk_fsize)
  !        enddo
  !      enddo
  !    else !not vector
  !      nbytes =  vtk_fsize*n_v4e*ne; write(fu) nbytes
  !      do j=1,size(pv_data%Xi,3)
  !        do i=1,size(pv_data%Xi,2)
  !          write(fu) real(pv_data%vars(ivar)%n_var(1,i,j), vtk_fsize)
  !        enddo
  !      enddo
  !    endif !vector
  !  endif !output
  !enddo

  buffer = '  </AppendedData>'//lf; write(fu) trim(buffer)
  buffer = '</VTKFile>'//lf; write(fu) trim(buffer)



  !!Write all the variables
  !write(fu,'(A)') '   <PointData Scalars="scalars">'




  !!close the file
  !write(fu,'(A)') '   </PointData>'
  !write(fu,'(A)') '  </Piece>'
  !write(fu,'(A)') ' </UnstructuredGrid>'
  !write(fu,'(A)') '</VTKFile>'



  close(fu,iostat=ierr)
 
end subroutine vtk_out_bin

end module mod_vtk_out

