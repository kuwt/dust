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


!> Module to treat the most simple input-output from ascii formatted data
!! files
module mod_cgns_io

use mod_param, only: &
  wp, max_char_len, nl

use mod_handling, only: &
  !error, 
  warning, info, printout, new_file_unit


!----------------------------------------------------------------------

implicit none

include 'cgnslib_f.h'

public :: read_mesh_cgns

private

character(len=*), parameter :: this_mod_name = 'mod_cgns_io'

!----------------------------------------------------------------------

contains

!----------------------------------------------------------------------

!>Read a mesh in a cgns format employing the cgns API 
!!
!! WARNING: this is still experimental and based on an incomplete reverse
!! engineering of cgns mesh files, it shall be extended to be more reliable
!! and general
subroutine read_mesh_cgns(mesh_file,ee,rr)
 
 character(len=*), intent(in) :: mesh_file
 integer  , allocatable, intent(out) :: ee(:,:) 
 real(wp) , allocatable, intent(out) :: rr(:,:) 
 

 integer :: index_file , ier  , nbase , ibase
 character(len=32) :: basename , zonename , sectionname
 character(len=72) :: eltype
 integer :: icelldim , ndim , nzone , izone , id , isec , nsec
 integer :: ieltype , nstart , nend , nbelem , nelem , nnod
 integer :: isize(1,3), iprec(1,3) , ipar(1,3)
 integer :: nelem_zone
 
 integer, pointer     :: pdata(:), nmixed(:)
 integer, allocatable :: elemcg(:,:)
 
 integer, allocatable :: estart(:) , eend(:)
! integer, allocatable :: ee_cgns(:)
 
 integer :: i1 , i , ielem

 character(len=*), parameter :: this_sub_name = 'read_mesh_cgns'
 
 
  !write(*,'(A,A)')' Opening cgns file: ', trim(fname)
  call CG_OPEN_F(trim(mesh_file), MODE_READ, index_file, ier)
  if ( ier .ne. ALL_OK ) then
    write(*,*) ' ** error when reading cgns file'
    call CG_ERROR_EXIT_F()
  end if

! Number of bases

  call CG_NBASES_F(INDEX_FILE, nbase, ier)
  if ( ier /= ALL_OK ) call CG_ERROR_EXIT_F()
  write(*,'(A,I2)') ' Number of bases in the cgns file: ', nbase
  if (nbase /= 1) then
    write(*,*) ' ** error: Only one base supported'
    write(*,*) ' ** Program stopped'
    stop
  end if
  ibase = 1

! Base data
 
  call CG_BASE_READ_F(INDEX_FILE, ibase, basename, icelldim, ndim, ier)
  if ( ier /= ALL_OK ) call CG_ERROR_EXIT_F()

  !write(*,*) ' icelldim = ' , icelldim ! todo
  !write(*,*) ' ndim     = ' , ndim    
 
!  Special case for 2d: if icelldim=2 we assume a pure 2d mesh with
!  coordinates in x and y only (not a surface mesh in 3d)
 
  if ((icelldim == 2) .and. (ndim > icelldim)) ndim = icelldim

!  Number of zones. Assume a zone is a region. Loop over the zones
 
  call CG_NZONES_F(INDEX_FILE, ibase, nzone, ier)
  if ( ier /= ALL_OK ) call CG_ERROR_EXIT_F()
  write(*,'(A,I2)') ' Number of zones in the cgns base: ', nzone

  if ( nzone .ne. 1 ) then 
    write(*,*)  ' up to now, only one zone is allowed ! ' 
    stop  
  end if

  do izone = 1,nzone

! *** edge format ***
!   call APPEND_CHILD_N(pmsh,'region',overwrite=.false.)
!   preg => NTH_CHILD(pmsh,izone,'region')
!   call APPEND_CHILD_L0(preg,'region_name','volume_elements')
 
!  Read zone data
 
    call CG_ZONE_READ_F(INDEX_FILE, ibase, izone, zonename, isize, ier)
    if ( ier /= ALL_OK ) call CG_ERROR_EXIT_F()

!  Check if unstructured zone. todo
 
    call CG_ZONE_TYPE_F(INDEX_FILE, ibase, izone, id, ier)
    if ( ier /= ALL_OK ) call CG_ERROR_EXIT_F()
    if ( id /= unstructured ) then
      write(*,*) ' ** error: Unstructured grid required. Found ', &
                  trim(zonetypename(id))
      write(*,*) ' ** Program stopped'
      stop
    end if
!

    write(*,*)
    write(*,'(A,I9)') ' Number of nodes:      ', isize(1,1)
    write(*,'(A,I9)') ' Number of cells:      ', isize(1,2)
    write(*,'(A,I2)') ' Coordinate dimension: ', ndim
    write(*,*)

!  Coordinates, allocate temporary storage coo

! *** edge *** 
!   pcoo => DATA_SET_NEW('coordinates', REAL_KIND_CHAR()//'F',&
!                        ndim, isize(1,1), rcflag=.true.)

    allocate( rr(ndim,isize(1,1)) )

! *** edge *** todo 
!   call APPEND_CHILD(preg,pcoo)
!   allocate(coo(isize(1,1)), stat = iastat)
!   if (iastat /= 0) then
!     write(*,*) ' ** error: Allocating memory'
!     write(*,*) ' ** Program stopped'
!     stop
!   end if
 
!  Read coordinates from file
 
! Define working precision. todo
    iprec = realdouble


    !TODO: passing rr(1,:) is not a contiguous memory location and
    !needs a temporary, and rises a warning. Consider either using 
    !an explicit temporary or a transposition of rr
    call CG_COORD_READ_F(INDEX_FILE, ibase, izone, 'CoordinateX', &
                       iprec, 1, isize(1,1), rr(1,:), ier)
    if ( ier /= ALL_OK ) call CG_ERROR_EXIT_F()
!   pcoo%rd(:,1) = coo
    call CG_COORD_READ_F(INDEX_FILE, ibase, izone, 'CoordinateY', &
                       iprec, 1, isize(1,1), rr(2,:), ier)
    if ( ier /= ALL_OK ) call CG_ERROR_EXIT_F()
!   pcoo%rd(:,2) = coo
    if (ndim == 3) then
      call CG_COORD_READ_F(INDEX_FILE, ibase, izone, 'CoordinateZ', &
                           iprec, 1, isize(1,1), rr(3,:), ier)
      if ( ier /= ALL_OK ) call CG_ERROR_EXIT_F()
!     pcoo%rd(:,3) = coo
    end if

! Check ----
!  do i1 = 1 , 10
!    write(*,*) rr(:,i1)
!  end do

   open(unit=21,file='./rr.dat')
   do i1 = 1 , size(rr,2)
     write(21,*) rr(:,i1)
   end do
   close(21)
! Check ----

!  Find number of sections
!
    call CG_NSECTIONS_F(INDEX_FILE, ibase, izone, nsec, ier)
    if ( ier /= ALL_OK ) call CG_ERROR_EXIT_F()
    write(*,'(A,I3)') ' Number of sections in the zone: ', nsec

! ++++++++++++++++++++++++++++++++++++++++++
! find the total number of elements -----
    allocate(estart(nsec)) ; estart = 0
    allocate(eend(nsec))   ; eend   = 0
    do isec = 1 , nsec   
 
      call CG_SECTION_READ_F(INDEX_FILE, ibase, izone, isec, sectionname, &
                             ieltype, nstart, nend, nbelem, ipar, ier)
      if ( ier /= ALL_OK ) call CG_ERROR_EXIT_F()
      nelem = nend - nstart + 1

      eend(isec)   = nend  
      estart(isec) = nstart

      write(*,'(A,I3,2A)') ' Section number:', isec,&
                           ', section name: ', sectionname
      write(*,'(A,I9)') ' Number of elements: ', nelem
      write(*,'(A,A)')  ' Element type: ', trim(elementtypename(ieltype))
      write(*,'(A,I9)') ' First parent data: ', ipar(1,1)

    end do
!   write(*,*) ' estart = ' , estart
!   write(*,*) ' eend   = ' , eend  
    nelem_zone = nend 
    allocate(ee(4,nelem_zone)) ; ee = 0
! ++++++++++++++++++++++++++++++++++++++++++


    do isec = 1, nsec
      call CG_SECTION_READ_F(INDEX_FILE, ibase, izone, isec, sectionname, &
                             ieltype, nstart, nend, nbelem, ipar, ier)
      if ( ier /= ALL_OK ) call CG_ERROR_EXIT_F()
      nelem = nend - nstart + 1
 
! *** edge ***  Set ffa element name.
! Set dust element name.
 
      select case (ieltype)
      case(node)
        eltype = 'node'
        nnod = 1
      case(BAR_2)
        eltype = 'bar2'
        nnod = 2
      case(TRI_3)
        eltype = 'tria3'
        nnod = 3
      case(QUAD_4)
        eltype = 'quad4'
        nnod = 4
      case(mixed)
        eltype = 'mixed'
        nnod = 9
      case default
        write(*,*) ' error: unsupported element type: ', ieltype, &
                   trim(elementtypename(ieltype))
        stop
      end select


!  Read cgns element data

      if (eltype == 'mixed') then
        allocate(elemcg(nnod*nelem,1))
      else
        allocate(elemcg(nnod, nelem))
      end if

      if (ipar(1,1)==1) then
        allocate(pdata(nelem*4))
      else
        allocate(pdata(3))
      end if

      call CG_ELEMENTS_READ_F(INDEX_FILE, ibase, izone, isec, &
                              elemcg, pdata, ier)
      if ( ier /= ALL_OK ) call CG_ERROR_EXIT_F()

!  Store data, mixed elements get special treatment

      if (eltype == 'mixed') then
        allocate(nmixed(mixed))
        nmixed = 0
        ielem = 0
        I = 1
        do while (ielem < nelem)
          ielem = ielem + 1
          select case (elemcg(I,1))
          case(node)
            nnod = 1
          case(BAR_2)
            nnod = 2
          case(TRI_3)
            nnod = 3
          case(QUAD_4)
            nnod = 4
          case(TETRA_4)
            nnod = 4
          case(PYRA_5)
            nnod = 5
          case(PENTA_6)
            nnod = 6
          case(HEXA_8)
            nnod = 8
          case default
            write(*,*) ' error: unsupported element type: ', elemcg(I,1), &
                       trim(elementtypename(elemcg(I,1)))
            stop
          end select
          ee(1:nnod,ielem+estart(isec)-1) = elemcg(I+1:I+nnod,1)
          nmixed(elemcg(I,1)) = nmixed(elemcg(I,1)) + 1
          I = I + nnod + 1
        end do

        deallocate(nmixed)

      else
        do id = estart(isec),eend(isec)
          ee(:,id) = elemcg(id,:)
        end do
      end if

      
      deallocate(elemcg)

    end do


  end do




  
end subroutine read_mesh_cgns

!----------------------------------------------------------------------

end module mod_cgns_io
