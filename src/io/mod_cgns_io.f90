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

use mod_stringtools, only: &
  IsInList, StripSpaces

 use cgns

!----------------------------------------------------------------------

implicit none

!include 'cgnslib_f.h'

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
subroutine read_mesh_cgns(mesh_file, sectionNamesUsed, ee, rr)
 
 character(len=*), intent(in) :: mesh_file
 character(len=*), intent(in) :: sectionNamesUsed(:)
 integer  , allocatable, intent(out) :: ee(:,:) 
 real(wp) , allocatable, intent(out) :: rr(:,:) 
 
 

 character(len=max_char_len), allocatable :: sectionNames(:)
 logical, allocatable :: selectedSection(:)
 integer :: index_file , ier  , nbase , ibase
 character(len=32) :: basename , zonename , sectionname, coordname(3)
 character(len=72) :: eltype
 integer :: icelldim , ndim , nzone , izone , id , isec , nsec
 integer :: ieltype , nstart , nend , nbelem , nelem , nnod
 integer :: nNodes, nNodesUsed, posNode
 integer :: isize(3), iprec , parent_flag
 integer :: nelem_zone
 integer :: nSectionsUsed, posSection
 logical :: sectionFound

 real(wp) , allocatable :: coordinateList(:) 
 integer, allocatable :: nodeMap(:), nodeList(:)
 
 integer, pointer     :: pdata(:), nmixed(:)
 integer, allocatable :: elemcg(:,:)
 
 integer, allocatable :: estart(:) , eend(:)
! integer, allocatable :: ee_cgns(:)
 
 integer :: i1 , i , ielem, ielem_section, iNode

 character(len=*), parameter :: this_sub_name = 'read_mesh_cgns'
 
  ! Name of the arrays containint the grid coordinates
  coordname(1) = 'CoordinateX'
  coordname(2) = 'CoordinateY'
  coordname(3) = 'CoordinateZ'

  nSectionsUsed = size(sectionNamesUsed)
 
  write(*,*) ' '
  write(*,*) 'Opening cgns file: ', trim(mesh_file)
  write(*,*) ' '
  write(*,*) ' '
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
                  'NONE' ! trim(zonetypename(id))
      write(*,*) ' ** Program stopped'
      stop
    end if

    nNodes = isize(1)

    write(*,*)
    write(*,'(A,I9)') ' Number of nodes:      ', isize(1)
    write(*,'(A,I9)') ' Number of cells:      ', isize(2)
    write(*,'(A,I2)') ' Coordinate dimension: ', ndim
    write(*,*)



!  Read coordinates from file
 
! Define working precision. todo
    iprec = realdouble


!  Find number of sections
!
    call CG_NSECTIONS_F(INDEX_FILE, ibase, izone, nsec, ier)
    if ( ier /= ALL_OK ) call CG_ERROR_EXIT_F()
    write(*,'(A,I3)') ' Number of sections in the zone: ', nsec

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! find the total number of elements and section names
    allocate(estart(nsec)) ; estart = 0
    allocate(eend(nsec))   ; eend   = 0
    allocate(sectionNames(nsec));
    allocate(selectedSection(nsec)); selectedSection = .false.
    do isec = 1 , nsec   
 
      call CG_SECTION_READ_F(INDEX_FILE, ibase, izone, isec, sectionname, &
                             ieltype, nstart, nend, nbelem, parent_flag, ier)
      if ( ier /= ALL_OK ) call CG_ERROR_EXIT_F()
      nelem = nend - nstart + 1

      eend(isec)   = nend  
      estart(isec) = nstart

      call StripSpaces(sectionname)
      sectionNames(isec) = sectionname

      write(*,'(A,I3,A,I10,2A)') ' Section number:', isec,&
                           '   Number of elements: ', nelem, &
                           '   Name: ', sectionname
      ! TODO: why is no longer working with CGNS 3.3?
      ! write(*,'(A,A)')  ' Element type: ', trim(ElementTypeName(ieltype))
      ! write(*,'(A,I9)') ' First parent data: ', parent_flag

    end do

    write(*,*) ' '

    ! Get number of elements considering the selection of sections
    if (nSectionsUsed > 0) then
      nelem_zone = 0
      write(*,'(A,I3)') ' Number of sections used: ', nSectionsUsed

      do isec = 1, nSectionsUsed
        write(*,*) ' Section Name: ', trim(sectionNamesUsed(isec))

        sectionFound = IsInList(sectionNamesUsed(isec), sectionNames, posSection)
        if (sectionFound) then
          selectedSection(posSection) = .true.
          nelem_zone = nelem_zone + eend(posSection) - estart(posSection) + 1
        else
          write(*,*) ' ** error: section with name ', trim(sectionNamesUsed(isec)), ' not found'
          write(*,*) ' **        see above for the list of available sections'
          write(*,*) ' ** Program stopped'
          stop
        endif
      enddo
    else
      write(*,*)  ' All sections will be used' 
      ! All sections will be selected
      selectedSection = .true.
      nelem_zone = nend 
    endif

    write(*,*) ' '
    write(*,*) 'Total number of elements ', nelem_zone
    write(*,*) ' '

    allocate(ee(4,nelem_zone)) ; ee = 0
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Save connectivity
    ielem = 0
    do isec = 1, nsec
      if (selectedSection(isec)) then
        call CG_SECTION_READ_F(INDEX_FILE, ibase, izone, isec, sectionname, &
                               ieltype, nstart, nend, nbelem, parent_flag, ier)
        if ( ier /= ALL_OK ) call CG_ERROR_EXIT_F()
        nelem = nend - nstart + 1
 
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
                     'NONE' ! trim(ElementTypeName(ieltype)) TODO: not working in CGNS 3.3
          stop
        end select


        !  Read cgns element data
        if (eltype == 'mixed') then
          allocate(elemcg(nnod*nelem,1))
        else
          allocate(elemcg(nnod, nelem))
        end if

        if (parent_flag==1) then
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
          ielem_section = 0
          I = 1
          do while (ielem_section < nelem)
            ielem = ielem + 1
            ielem_section = ielem_section + 1
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
                         'NONE' ! trim(ElementTypeName(elemcg(I,1))) TODO: not working in CGNS 3.3
              stop
            end select
            ee(1:nnod,ielem) = elemcg(I+1:I+nnod,1)
            nmixed(elemcg(I,1)) = nmixed(elemcg(I,1)) + 1
            I = I + nnod + 1
          end do

          deallocate(nmixed)

        else
          do id = 1, eend(isec) - estart(isec) + 1
            ielem = ielem + 1
            ee(1:nnod,ielem) = elemcg(:,id)
          end do
        end if

        
        deallocate(elemcg)

      endif
    end do ! Loop on sections


    ! If a selection of elements is performed performs also a selection on nodes
    if (nSectionsUsed > 0) then
      nNodesUsed = 0
      allocate(nodeMap(nNodes)); nodeMap = 0
      allocate(nodeList(nNodes)); nodeList = 0

      ! Select nodes used
      do ielem = 1, nelem_zone
        do i1 = 1, 4
          if (ee(i1,ielem) > 0) then
            posNode = ee(i1,ielem)
            if (nodeMap(posNode) .eq. 0) then
              nNodesUsed = nNodesUsed + 1;
              nodeMap(posNode) = nNodesUsed
              nodeList(nNodesUsed) = posNode
            endif
            ee(i1,ielem) = nodeMap(posNode)
          endif
        enddo
      enddo

      ! Get nodes
      allocate(coordinateList(nNodes))
      allocate(rr(ndim,nNodesUsed))
      do i = 1, ndim
        coordinateList = 0.0_wp;
        call cg_coord_read_f(INDEX_FILE, ibase, izone, coordname(i), &
                             iprec, 1, isize(1), coordinateList, ier);
        if ( ier /= ALL_OK ) call CG_ERROR_EXIT_F()

        do iNode = 1, nNodesUsed
          rr(i,iNode) = coordinateList(nodeList(iNode))
        enddo
        
      enddo

    ! No selection of sections performed, read all nodes
    else 
      allocate(coordinateList(nNodes))
      allocate(rr(ndim,nNodes))
      do i = 1, ndim
        coordinateList = 0.0_wp;

        !TODO: passing rr(1,:) is not a contiguous memory location and
        !needs a temporary, and rises a warning. Consider either using 
        !an explicit temporary or a transposition of rr
        call cg_coord_read_f(INDEX_FILE, ibase, izone, coordname(i), &
                             iprec, 1, isize(1), coordinateList, ier);
        if ( ier /= ALL_OK ) call CG_ERROR_EXIT_F()

        rr(i,:) = coordinateList
        
      enddo

    endif

  end do ! Loop on zones


  
end subroutine read_mesh_cgns

!----------------------------------------------------------------------

end module mod_cgns_io
