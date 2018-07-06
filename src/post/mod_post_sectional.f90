module mod_post_sectional

use mod_param, only: &
  wp, nl, max_char_len, extended_char_len , pi

use mod_handling, only: &
  error, warning, info, printout, dust_time, t_realtime, new_file_unit

use mod_parse, only: &
  t_parse, &
  getrealarray, getreal, getint, getlogical, &
  getsuboption,  countoption
 
use mod_hdf5_io, only: &
   initialize_hdf5, destroy_hdf5, &
   h5loc, &
   new_hdf5_file, &
   open_hdf5_file, &
   close_hdf5_file, &
   new_hdf5_group, &
   open_hdf5_group, &
   close_hdf5_group, &
   write_hdf5, &
   read_hdf5, &
   read_hdf5_al, &
   check_dset_hdf5

use mod_geo_postpro, only: &
  load_components_postpro, update_points_postpro , prepare_geometry_postpro , &
  prepare_wake_postpro  ! expand_actdisk_postpro, 

use mod_geometry, only: &
  t_geo, t_geo_component

use mod_post_load, only: &
  load_refs, load_res

use mod_dat_out, only: &
   dat_out_sectional

use mod_tecplot_out, only: &
  tec_out_sectional

use mod_math, only: &
  cross

implicit none

! type rof box sectional loads
type t_box_secloads
  integer :: nelems
  integer  , allocatable :: elems(:)
  real(wp) , allocatable :: fracs(:)
  real(wp) , allocatable :: cen(:,:)
end type t_box_secloads

public :: post_sectional

private

character(len=max_char_len), parameter :: this_mod_name = 'mod_post_sectional'

contains

! ---------------------------------------------------------------------- 
!TODO:
! - deal with multiple components (DONE)
! - parametric components aligned with y-axis
! - tol_y_cen is hard-coded. write tol_y_cen as an input
! - ...
subroutine post_sectional ( sbprms , bxprms , basename , data_basename , an_name , ia , &
                            out_frmt , comps , components_names , all_comp , &
                            an_start , an_end , an_step )
type(t_parse), pointer :: sbprms
type(t_parse), pointer :: bxprms
character(len=*) , intent(in) :: basename
character(len=*) , intent(in) :: data_basename
character(len=*) , intent(in) :: an_name
integer          , intent(in) :: ia
character(len=*) , intent(in) :: out_frmt
type(t_geo_component), allocatable , intent(inout) :: comps(:)
character(len=max_char_len), allocatable , intent(inout) :: components_names(:)
logical , intent(in) :: all_comp
integer , intent(in) :: an_start , an_end , an_step

character(len=max_char_len) :: cname !, msg
integer(h5loc) :: floc, gloc, cloc
real(wp), allocatable :: refs_R(:,:,:), refs_off(:,:)
real(wp), allocatable :: refs_G(:,:,:), refs_f(:,:)
real(wp), allocatable :: vort(:), cp(:)
real(wp), allocatable :: points(:,:)
integer , allocatable :: elems(:,:)
integer :: nelem
integer :: n_comp , n_comp_tot , i_comp , id_comp , ax_coor , ref_id
character(len=max_char_len), allocatable :: all_components_names(:)
character(len=max_char_len), allocatable :: components_names_tmp(:)
integer :: nelem_span , nelem_chor , n_sect , n_time

real(wp), allocatable :: time(:)
real(wp) :: t
integer :: ires

! sectional loads: paramteric components -------------------
real(wp) :: axis_dir(3) , axis_nod(3) 
character(len=max_char_len) :: comp_name_stripped
real(wp) , allocatable :: sec_loads(:,:,:)
integer :: is 
integer :: n_loads = 4   ! F and moment around an axis
real(wp) , allocatable :: ref_mat(:,:) , off_mat(:,:)
real(wp) , allocatable :: y_cen(:) , y_span(:)
real(wp) , parameter   :: tol_y_cen = 1.0e-3_wp
real(wp) , allocatable :: r_axis(:,:) , r_axis_bas(:,:)

real(wp) :: F_bas(3) , F_bas1(3)
real(wp) :: M_bas(3) , M_bas1(3)
integer :: ie , ic , it , i1

! sectional loads: box -------------------------------------
character(len=max_char_len) :: comp_input
integer :: n_box
logical :: reshape_box
!TODO: clean declarations
real(wp) :: nCoordMaxBox
real(wp) , allocatable :: box_coord(:,:) , box_coord_tmp(:,:) , box_coord_rotated(:,:)
real(wp) :: vVec(3) , bVec(3) , wVec(3) , baseLen(2) , heigLen(2) , spanLen
real(wp) :: nVec(3) , ref_node(3) ! , axis_mom(3)
real(wp) :: nCoordMax , nCoordMin
real(wp) :: nCoord , dnCoord
real(wp) , allocatable :: nCoordSec(:)
real(wp) , allocatable :: nCoordCen(:,:)
real(wp) :: normalLateralFaces(3,4)
real(wp) :: refPoiLateralFaces(3,4)
real(wp) :: b1Vec(3) , b2Vec(3)
type(t_box_secloads) , allocatable  :: box_secloads(:)
real(wp) :: distance(4)
real(wp) , parameter :: nCoord_relOff = 0.05_wp
! "extra-param". TODO: check if they are called before #### box
real(wp) :: nCoordVert(4)   ! TRI o QUAD elements
integer  :: secVert(4)   ! TRI o QUAD elements
real(wp) :: interSectPoints(3,2)
real(wp) :: interSectAreas(2) , interSectCen(3,2)
integer  :: nInterSect , index2
real(wp) :: node1(3) , node2(3) , box_secloads_cen(3)
integer  :: iSec1 , iSec2 , sec1_nVer , sec2_nVer
integer  :: iSecWork , iSecOth
integer  :: nNodeInt(2,2)
integer  :: nver , iv 

character(len=max_char_len) :: filename

character(len=max_char_len), parameter :: this_sub_name = 'post_sectional'

  write(*,*) nl//' Analysis:',ia,' post_sectional() +++++++++ '//nl

  ! DONE DONE DONE DONE DONE DONE DONE DONE DONE DONE DONE DONE DONE DONE DONE !
  !TODO: deal with multiple components.                                        !
  ! A good strategy could be to avoid dealing with multiple components:        !
  ! --> check if load_components_postpro() subroutine is ok for this strategy  !
  !     OR stop everything before, reading only the names of the actual        !
  !        components (around l.90).                                           !
  !        Look for: ! *** Check if the component really exists *** !          !
  ! DONE DONE DONE DONE DONE DONE DONE DONE DONE DONE DONE DONE DONE DONE DONE !
  
  ! Some warnings and errors -------------------------------
  ! WARNING: sectional loads are computed for one components only
  !  at time. If no component is specified --> return
  n_comp = countoption(sbprms,'Component')
  if ( n_comp .le. 0 ) then
    call warning(trim(this_mod_name),trim(this_sub_name), &
           'No component specified for ''sectional_loads''&
         & analysis. Skipped analysis.')
    return
  else if ( n_comp .ge. 2 ) then
    call warning(trim(this_mod_name),trim(this_sub_name) , &
          'WARNING. More than one component specified &
         &for ''sectional_loads'' analysis: just the first one is considered. &
         &Please run another ''sectional_analysis'' if you need it on more than &
         &component.')
  end if
  
  ! check if components_names is allocated (it should alwasy be allocated)
  if ( .not. allocated(components_names) ) then
    call warning(trim(this_mod_name),trim(this_sub_name), &
           ' components_name .not. a, intent(inout) llocated. Something unexpected&
           & happened. Skipped analy, intent(inout) sis.')
    return
  end if 
  
  ! ######### *** Check if the component really exists *** ######### !
  call open_hdf5_file(trim(data_basename)//'_geo.h5', floc)
  call open_hdf5_group(floc,'Components',gloc)
  call read_hdf5(n_comp_tot,'NComponents',gloc)
  
  allocate(all_components_names(n_comp_tot))
  do i_comp = 1 , n_comp_tot
    write(cname,'(A,I3.3)') 'Comp',i_comp
    call open_hdf5_group(gloc,trim(cname),cloc)
    call read_hdf5(all_components_names(i_comp),'CompName',cloc)
    call close_hdf5_group(cloc)
  end do
  call close_hdf5_group(gloc)
  
  id_comp = -333
  do i_comp = 1 , n_comp_tot
     if ( trim(components_names(1)) .eq. trim(all_components_names(i_comp)) ) then
       id_comp = i_comp
     end if
  end do
  if ( id_comp .eq. -333 ) then ! component not found
    write(*,*) ' The available components are: '
    do i_comp = 1 , n_comp_tot 
      write(*,'(A,I0,A)') ' comp. ' , i_comp , &
                          ':  '//trim(all_components_names(i_comp))
    end do
    call error('dust_post','','No valid component defined. STOP')
  end if
  
  write(*,*) ' Component : ' , trim(components_names(1))//nl
  ! ######### *** Check if the component really exists *** # DONE ## !
  
  ! load only the first component just once -----------------
  call open_hdf5_file(trim(data_basename)//'_geo.h5', floc)
  
  allocate(components_names_tmp(1))
  components_names_tmp(1) = trim(components_names(1))
  
  call load_components_postpro(comps, points, nelem, floc, & 
                               components_names_tmp,  all_comp)
  call close_hdf5_file(floc)
  
  ! Prepare_geometry_postpro
  call prepare_geometry_postpro(comps)
  
  ! Read the axis for the computation of the sectional loads
  axis_dir = getrealarray(sbprms,'AxisDir',3)
  axis_nod = getrealarray(sbprms,'AxisNod',3)

  ! if a box is defined ---> use the box
  comp_input = trim(comps(1)%comp_input )
  n_box = countoption(sbprms,'BoxSect')
  !debug
  write(*,*) ' n_box : ' , n_box
  if ( n_box .eq. 1 ) comp_input = 'cgns' 
  
  select case( trim(comp_input) )
  
  case ( 'parametric' )
  
    !TODO: move to a subroutine
    
    ! Some assumptions ---------------
    id_comp = 1   ! 1. only one component is loaded
    ax_coor = 2   ! 2. the parametric elements is defined
                  !    along the y-axis
    
    nelem_span = comps(id_comp)%parametric_nelems_span
    nelem_chor = comps(id_comp)%parametric_nelems_chor
    n_sect = nelem_span
    
    
    ! ######################################################################
    ! Find the coordinates of the reference points in the local reference frame
    
    ! Find the coordinate along the axis of the sections -------------------
    ! the <ax_coor> coordinate y_cen of each section is built with the first
    ! panel in chord. Then the distance between the <ax_coor> of the centre
    ! of the other panels and y_cen is checked
    !TODO: find y-coord is general enough for all the 'parametric' components
    allocate(y_cen(n_sect),y_span(n_sect)) 
    ie = 0  
    do is = 1 , n_sect
      ie = ie + 1
      y_cen(is) = sum(comps(id_comp)%loc_points(ax_coor,comps(id_comp)%el(ie)%i_ver)) / &
                  dble(comps(id_comp)%el(ie)%n_ver)
      y_span(is) = abs ( comps(id_comp)%loc_points(ax_coor,comps(id_comp)%el(ie)%i_ver(1) ) - &
                         comps(id_comp)%loc_points(ax_coor,comps(id_comp)%el(ie)%i_ver(2) ) )
      do ic = 2 , nelem_chor ! check 
        ie = ie + 1
        if ( abs( y_cen(is) - & !comps(id_comp)%el(ie)%cen(2) ) .gt. tol_y_cen ) then
             sum(comps(id_comp)%loc_points(ax_coor,comps(id_comp)%el(ie)%i_ver)) / &
             dble(comps(id_comp)%el(ie)%n_ver) ) .gt. tol_y_cen ) then
          call error('dust_post','','Wrong section definition for component '//&
             trim(comps(id_comp)%comp_name)//' for ''sectional_loads'' analysis.&
             & STOP')
        end if
      end do
    end do
    
    ! Find the coordinate of the reference points on the axis --------------
    !  ( with coord. y_cen )
    if ( abs(axis_dir(2)) .lt. 1e-6 ) then
      call error('dust_post','','Wrong definition of the axis in&
           & sectional_loads analysis: abs(axis_dir(2)) .lt. 1e-6.&
           & STOP')
    end if
    allocate(r_axis(3,n_sect),r_axis_bas(3,n_sect))
    do is = 1 , n_sect
      r_axis(:,is) = axis_nod + &
              ( y_cen(is) - axis_nod(ax_coor) )/axis_dir(ax_coor) * axis_dir
    end do
    ! only initialisation here:
    ! - r_axis: coordinates in the local ref.frame
    ! - r_axis_bas: coordinates in the base ref.frame
    r_axis_bas = r_axis
    
    ! Find the ref_id or the reference frame where loads are projected -----
    ref_id = comps(id_comp)%ref_id
    !check
    write(*,'(A,I0,A)') ' ref. for force and moment projection (', &
         ref_id , ')  '//trim(comps(id_comp)%ref_tag)
    
    ! allocate tmp array to store the results --------------
    n_time = (an_end-an_start)/an_step + 1 ! int general eger division
    allocate( sec_loads(n_time,n_sect,n_loads) ) ; sec_loads = -333.0_wp ! initialisation for DEBUG
    allocate( time(n_time) ) ; time = -333.0_wp
    allocate( ref_mat(n_time,9) , off_mat(n_time,3) ) 
    ires = 0
    do it = an_start, an_end, an_step
      ires = ires + 1
      
      ! Open the file:
      write(filename,'(A,I4.4,A)') trim(data_basename)//'_res_',it,'.h5'
      call open_hdf5_file(trim(filename),floc)
    
      ! Load the references and move the points ---
      call load_refs(floc,refs_R,refs_off)
      ! Move the points ---------------------------
      call update_points_postpro(comps, points, refs_R, refs_off)
      ! Load the results --------------------------
      call load_res(floc, comps, vort, cp, t)
    
      call close_hdf5_file(floc)
    
      ! from local coordinates to coordinates in the base ref.frame
      do is = 1 , n_sect
        r_axis_bas(:,is) = matmul( refs_R(:,:,ref_id) , r_axis(:,is) ) + &
                           refs_off(:,ref_id)
      end do
    
      ! compute sectional loads. Loop over the panels of each section
      ie = 0
      do is = 1 , n_sect ! loop over sections
        ! force
        F_bas = 0.0_wp ; M_bas = 0.0_wp 
        do ic = 1 , nelem_chor ! loop over chord 
          ie = ie + 1
          F_bas1 = comps(id_comp)%el(ie)%dforce
          F_bas  = F_bas + F_bas1
          M_bas  = M_bas + cross( comps(id_comp)%el(ie)%cen &
                         - r_axis_bas(:,is) , F_bas1 )
        end do ! loop over chord
    
        ! From global to local coordinates of forces and moments 
        sec_loads(ires,is,1:3) = matmul( &
             transpose( refs_R(:,:, ref_id) ) , F_bas ) / y_span(is)
        ! moment ( only the component around the <ax_coord> )
        M_bas = matmul( &
             transpose( refs_R(:,:, ref_id) ) , M_bas )
        sec_loads(ires,is,4) = M_bas(ax_coor) / y_span(is)
    
      end do ! loop over sections
    
      ref_mat(ires,:) = reshape(refs_R(:,:,ref_id),(/ 9 /))
      off_mat(ires,:) = refs_off(:,ref_id)
      time(ires) = t
    
    end do
    
    select case(trim(out_frmt))
     case('dat')
      write(filename,'(A)') trim(basename)//'_'//trim(an_name) 
      call dat_out_sectional ( filename, y_cen, y_span, time, sec_loads, &
                              ref_mat , off_mat ) 
     case('tecplot')
      write(filename,'(A)') trim(basename)//'_'//trim(an_name)//'.plt' 
      call tec_out_sectional ( filename, time, sec_loads, y_cen, y_span ) 
    end select

    write(*,*) nl//nl//' end of sectional loads'//nl//nl
    
    deallocate(r_axis,r_axis_bas)
    deallocate(y_cen,y_span)
    deallocate(sec_loads,time,ref_mat,off_mat)
  
  
  case ( 'cgns' )

    !TODO: move to a subroutine
    
    ! Some assumptions ---------------
    id_comp = 1   ! 1. only one component is loaded
  
    ! Read input for box and build the 8 nodes of the box ----------------------------- 
    call getsuboption(sbprms,'BoxSect',bxprms)
    allocate(box_coord(8,3),box_coord_tmp(8,3))
    ref_node = getrealarray(bxprms,'refNode',3)
    vVec = getrealarray(bxprms,'faceVec',3)
    baseLen = getrealarray(bxprms,'faceBas',2) 
    heigLen = getrealarray(bxprms,'faceHei',2) 
    bVec = getrealarray(bxprms,'spanVec',3)
    spanLen = getreal(bxprms,'spanLen') 
    n_sect = getint(bxprms,'numSect')
    reshape_box = getlogical(bxprms,'reshapeBox') 
    nVec = getrealarray(sbprms,'AxisMom',3) 

    ! Normalise
    vVec = vVec / norm2(vVec)
    bVec = bVec / norm2(bVec)
    nVec = nVec / norm2(nVec)

    wVec = cross(vVec,nVec) ; wVec = wVec / norm2(wVec)

    !check
    write(*,*) ' n : ' , nVec
    write(*,*) ' b : ' , bVec
    write(*,*) ' v : ' , vVec
    write(*,*) ' w : ' , wVec

    box_coord(1,:) = ref_node
    box_coord(2,:) = box_coord(1,:) + baseLen(1) * vVec  
    box_coord(3,:) = box_coord(2,:) + heigLen(1) * wVec  
    box_coord(4,:) = box_coord(1,:) + heigLen(1) * wVec  

    box_coord(5,:) = box_coord(1,:) + spanLen * bVec
    box_coord(6,:) = box_coord(5,:) + baseLen(2) * vVec  
    box_coord(7,:) = box_coord(6,:) + heigLen(2) * wVec  
    box_coord(8,:) = box_coord(5,:) + heigLen(2) * wVec  

    !check
    write(*,*) nl//' box_coord : ' 
    do i1 = 1 , 8 
      write(*,*) box_coord(i1,:)
    end do
    write(*,*)

    ! Find if the whole component is contained in the box and/or reshape
    !  the box  ( in the direction perpendicolar to the nVec = axis_mom )
    ! nCoord: distance of a point from the first reference plane
    ie = 1
    nCoordMax = sum ( ( comps(id_comp)%loc_points( :, ie ) &
     - box_coord(1,:) ) * nVec )
    nCoordMin = nCoordMax !0.0_wp
    do ie = 2 , size( comps(id_comp)%loc_points , 2 ) ! size(comps(id_comp)%el)
      nCoord = sum ( ( comps(id_comp)%loc_points( :, ie ) &
        - box_coord(1,:) ) * nVec )
!debug !write(*,*) ' nCoord : ' , nCoord
      if ( nCoord .gt. nCoordMax ) nCoordMax = nCoord
      if ( nCoord .lt. nCoordMin ) nCoordMin = nCoord 
    end do
    !check
    write(*,*) nl//' nCoordMin,Max : ' , nCoordMin , nCoordMax , nl

    ! Add small offset to the min and max nCoord ---------------------
    dnCoord = ( nCoordMax - nCoordMin ) / dble(n_sect) 
    nCoordMin = nCoordMin - dnCoord * nCoord_relOff
    nCoordMax = nCoordMax + dnCoord * nCoord_relOff
    dnCoord = ( nCoordMax - nCoordMin ) / dble(n_sect) 

    !check. TODO: add the offset to min,max() in the following lines
    nCoordMaxBox = sum(bVec*nVec) * spanLen 
    nCoordMin = max(nCoordMin,0.0_wp)
    nCoordMax = min(nCoordMax,nCoordMaxBox)

    write(*,*) 'after the offset '
    write(*,*) ' nCoordMin,Max : ' , nCoordMin , nCoordMax , nl

    ! Reshape the box, if it is "too large" --------------------------
    if ( reshape_box ) then

      write(*,*) ' Reshape box: ' , reshape_box , nl

!     ! Update box corners: update only if the box is too large 
!     allocate(box_coord_tmp(8,3))
      box_coord_tmp = 0.0_wp
      ! check
      write(*,*) ' nCoordMaxBox : ' , nCoordMaxBox 
      box_coord_tmp(1,:) = box_coord(1,:) + ( box_coord(5,:) - box_coord(1,:) ) * nCoordMin / spanLen
      box_coord_tmp(2,:) = box_coord(2,:) + ( box_coord(6,:) - box_coord(2,:) ) * nCoordMin / spanLen
      box_coord_tmp(3,:) = box_coord(3,:) + ( box_coord(7,:) - box_coord(3,:) ) * nCoordMin / spanLen
      box_coord_tmp(4,:) = box_coord(4,:) + ( box_coord(8,:) - box_coord(4,:) ) * nCoordMin / spanLen
      box_coord_tmp(5,:) = box_coord(1,:) + ( box_coord(5,:) - box_coord(1,:) ) * nCoordMax / spanLen
      box_coord_tmp(6,:) = box_coord(2,:) + ( box_coord(6,:) - box_coord(2,:) ) * nCoordMax / spanLen
      box_coord_tmp(7,:) = box_coord(3,:) + ( box_coord(7,:) - box_coord(3,:) ) * nCoordMax / spanLen
      box_coord_tmp(8,:) = box_coord(4,:) + ( box_coord(8,:) - box_coord(4,:) ) * nCoordMax / spanLen

      box_coord = box_coord_tmp

!     deallocate(box_coord_tmp)
      !check
      write(*,*) ' Updated box_coord : ' 
      do i1 = 1 ,8 
        write(*,*) box_coord(i1,:)
      end do

    end if ! end reshape

    !check
    write(*,*) nl//' before sectional division: '
    write(*,*) ' nCoordMin,Max : ' , nCoordMin , nCoordMax

    ! Find the n-coord of the sections and the reference points at the centre of the sections
    ! obs: so far, the nCoord has been defined w.r.t. the original first face of the box.
    !      from now on, nCoordSec, nCoordCen are defined with respect to the first section,
    !      that can be moved
    allocate(nCoordSec(n_sect+1))
    allocate(nCoordCen(3,n_sect)) ; nCoordCen = 0.0_wp
    
    nCoordSec(1) = 0.0_wp ! + nCoordMin
    write(*,*)  nl//' nCoordSec(  ',1,') : ' , '                           ' , nCoordSec(1)
    do is = 1 , n_sect
      nCoordSec(is+1) = (nCoordMax-nCoordMin) * dble(is) / dble(n_sect) ! + nCoordMin
      nCoordCen(2,is) = (nCoordMax-nCoordMin) * ( dble(is) - 0.5_wp)  / dble(n_sect) ! + nCoordMin
      write(*,*)  ' nCoordCen(:,',is  ,') : ' , nCoordCen(:,is)
      write(*,*)  ' nCoordSec(  ',is+1,') : ' , '                            ', nCoordSec(is+1)
    end do

    !TODO: use generalised version of this formula.
    ! this is valid only if nVec = yVec  
    allocate(y_cen( n_sect)) ; y_cen = nCoordCen(2,:) ! TODO: ask if this offset is required: + ref_node(2)
    allocate(y_span(n_sect)) ; y_span= nCoordSec(2:) - nCoordSec(1:n_Sect) 

    allocate(r_axis(3,n_sect),r_axis_bas(3,n_sect))
    do is = 1 , n_sect
      r_axis(:,is) = axis_nod + &
              ( y_cen(is) - axis_nod(2) )/axis_dir(2) * axis_dir
    end do
    !debug
    write(*,*) ' r_axis : '
    do is = 1 , n_sect
      write(*,*) r_axis(:,is)
    end do
    ! only initialisation here:
    ! - r_axis: coordinates in the local ref.frame
    ! - r_axis_bas: coordinates in the base ref.frame
    r_axis_bas = r_axis

    ! Define the unit normal vectors to the 4 planar lateral faces of the boxes
    !    and the reference point, used as the origin for measuring distance from the plane
    refPoiLateralFaces(:,1) = box_coord(1,:)
    refPoiLateralFaces(:,2) = box_coord(1,:)
    refPoiLateralFaces(:,3) = box_coord(3,:)
    refPoiLateralFaces(:,4) = box_coord(4,:)

    b1Vec = box_coord(6,:) - box_coord(2,:) ; b1Vec = b1Vec / norm2(b1Vec) 
    b2Vec = box_coord(8,:) - box_coord(4,:) ; b2Vec = b2Vec / norm2(b2Vec) 
    normalLateralFaces(:,1) =  cross( bVec,wVec) / norm2(cross( bVec,wVec))
    normalLateralFaces(:,2) =  cross( vVec,bVec) / norm2(cross( vVec,bVec))
    normalLateralFaces(:,3) = -cross(b1Vec,wVec) / norm2(cross(b1Vec,wVec))
    normalLateralFaces(:,4) =  cross(b2Vec,vVec) / norm2(cross(b2Vec,vVec))

    !check
    write(*,*) nl//' normalLateralFaces '
    do i1 = 1 , 4
      write(*,*) normalLateralFaces(:,i1)
    end do

! ######################################################################    
! #      TO BE CHECKED       ###########################################    
! ######################################################################    
    ! Find the elements belonging to each section ---------------
    ! allocate and initialise
    allocate(box_secloads(0:n_sect+1))
    do is = 0 , n_sect+1 ! 0,n_sect+1 to add dummy ghost extrem sections
      box_secloads(is)%nelems = 0 
      allocate(box_secloads(is)%elems(size(comps(id_comp)%el))) ; box_secloads(is)%elems = -333
      allocate(box_secloads(is)%fracs(size(comps(id_comp)%el))) ; box_secloads(is)%fracs = -333.3_wp
      allocate(box_secloads(is)%cen(3,size(comps(id_comp)%el))) ; box_secloads(is)%cen   = -333.3_wp
    end do


    do ie = 1 , size(comps(id_comp)%el) ! loop over elems

      ! Reset some useful variables ----------------------------- 
      nCoordVert = 0.0_wp
      secVert = 0
      interSectPoints = 0.0_wp
      interSectAreas  = 0.0_wp
      nInterSect = 0
      nNodeInt = 0
      nver = comps(id_comp)%el(ie)%n_ver 

      ! Compute cen and area of the elements --------------------
      ! since no prepare or load geom has been called yet.
      ! TODO: move this part in the time loop, first timestep
      ! compute the centre and the area of the element
      comps(id_comp)%el(ie)%cen =  &
        sum(comps(id_comp)%loc_points(:,comps(id_comp)%el(ie)%i_ver-comps(id_comp)%i_points(1)+1),2) / &
                  dble(comps(id_comp)%el(ie)%n_ver)
      comps(id_comp)%el(ie)%area =  0.5_wp * norm2( &
                  cross ( comps(id_comp)%loc_points(:,comps(id_comp)%el(ie)%i_ver(3)) - & 
                          comps(id_comp)%loc_points(:,comps(id_comp)%el(ie)%i_ver(1)) , & 
                          comps(id_comp)%loc_points(:,comps(id_comp)%el(ie)%i_ver(2)) - & 
                          comps(id_comp)%loc_points(:,comps(id_comp)%el(ie)%i_ver(nver)) ) )

      ! Compute the distance of the centre of the elem from -----
      ! lateral faces of the box
      do i1 = 1 , 4
        distance(i1) = sum ( ( comps(id_comp)%el(ie)%cen - refPoiLateralFaces(:,i1) ) * normalLateralFaces(:,i1) )
!       !debug
!       write(*,*) distance(i1)
      end do      

      ! if the centre of the element belongs to the box -> compute the slice contributions
      if ( all( distance .gt. 0.0_wp ) ) then

        ! nCoord of the cen of the elem (useless??)
        nCoord = sum ( ( comps(id_comp)%el(ie)%cen - refPoiLateralFaces(:,1) ) * nVec )

        ! nCoord of the vertices of the elem ( (r-refPoi(1))\cdot nVec )
        do iv = 1 , nver
          nCoordVert(iv) = &
             sum ( ( comps(id_comp)%loc_points(:, &
                                      comps(id_comp)%el(ie)%i_ver(iv) ) &
                     - refPoiLateralFaces(:,1) ) * nVec )
        end do
 
        ! elems with at least one point s.t. nCoor \in (nCoorMix,nCoorMax) 
        if ( ( maxval(nCoordVert) .ge. 0.0_wp ) .and. &
             ( minval(nCoordVert) .le. nCoordMax-nCoordMin ) ) then

           ! find the sections where the ver of the elem belong -------------
           ! !!! if the one vert is outside the box, secVert keep 0 value !!!
           do iv = 1 , nver
             if ( nCoordVert(iv) - nCoordSec(n_sect+1) .gt. 0.0_wp  ) then
                 secVert(iv) = n_sect + 1
             elseif ( nCoordVert(iv) - nCoordSec(1) .le. 0.0_wp  ) then
                 secVert(iv) = 0
             else
               do is = 1 , n_sect !TODO: improve this loop; exit when found
                 if ( ( nCoordVert(iv) - nCoordSec(is)   .gt. 0.0_wp ) .and. & 
                      ( nCoordVert(iv) - nCoordSec(is+1) .le. 0.0_wp ) ) then
                   secVert(iv) = is
                 end if 
               end do
             end if 
           end do
           write(*,*) ' secVert : ' , secVert(1:nver)

           ! If the elements belong to: -------------------------------------
           ! (a) .gt. 2 sections ----> error.
           ! (b) .eq. 1 section  ----> easy  
           ! (c) .eq. 2 section  ----> find intersections and partial contributions
           !                          to the sections

           ! (a) check if the elements belong to less than three sections --> otherwise ERROR
           if ( maxval(secVert(1:nver) ) & 
              - minval(secVert(1:nver) ) .gt. 1 ) then
             write(*,*) ' comps(id_comp)%comp_name : ' , comps(id_comp)%comp_name 
             write(*,*) ' ie : ' , ie  
             write(*,*) ' secVert : ' , secVert(1:nver)
             call warning('dust_post','','The element above belongs to more than two&
                  & sections in sectional_loads analysis. STOP')
             return
           end if
           ! (b)
           if ( all(  secVert .eq. secVert(1) ) ) then
             box_secloads(secVert(1))%nelems = box_secloads(secVert(1))%nelems + 1
             box_secloads(secVert(1))%elems(box_secloads(secVert(1))%nelems) = ie
             box_secloads(secVert(1))%fracs(box_secloads(secVert(1))%nelems) = 1.0_wp
             box_secloads(secVert(1))%cen(:,box_secloads(secVert(1))%nelems) = & 
                                                            comps(id_comp)%el(ie)%cen  
           else !(c)

             !TODO: treat elements partially belonging to the box.  if ( ... )

             ! Find intersections with the plane delimiting the sections ----
             do iv = 1 , nver

               index2 = mod(iv,comps(id_comp)%el(ie)%n_ver)+1 ! following node

               ! if two consecutive nodes belong to different sections
               !  ---> find intersection 
               if ( secVert( index2 ) .ne. secVert(iv) ) then
                 nInterSect = nInterSect + 1
                 node1 = comps(id_comp)%loc_points(:, comps(id_comp)%el(ie)%i_ver(iv) )
                 node2 = comps(id_comp)%loc_points(:, comps(id_comp)%el(ie)%i_ver(index2) ) 
                 interSectPoints(:,nInterSect) = node1 + (node2-node1) * &
                      ( nCoordSec( max( secVert(iv), secVert(index2) ) ) - nCoordVert(iv) ) / &
                      ( nCoordVert(index2) - nCoordVert(iv) )  ! ^---- max(.,.): id. of the section

                 if ( secVert(iv) .eq. secVert(1) ) then ! sec.1
                   nNodeInt(1,nInterSect) = iv     ; nNodeInt(2,nInterSect) = index2
                 else ! sec.2
                   nNodeInt(1,nInterSect) = index2 ; nNodeInt(2,nInterSect) = iv
                 end if

               end if
             end do

             ! Count the nodes belonging to the splitted elem ---------------
             iSec1 = secVert(1) ; sec1_nVer = 1 ; sec2_nVer = 0
             do iv = 2 , comps(id_comp)%el(ie)%n_ver
               if ( secVert(iv) .ne. iSec1 ) then
                 iSec2 = secVert(iv) ; sec2_nVer = sec2_nVer + 1
               else
                 sec1_nVer = sec1_nVer + 1
               end if
             end do

             ! Compute the area of the subelem with the minimum n.of points -
             !  (since it could be only TRI or QUAD) and compute the area of
             !  the othe elem as a difference.
             if ( sec1_nVer .le. sec2_nVer ) then
               if ( sec1_nVer .gt. 2 ) then
                  write(*,*) ' error in sectional loads . stop ' ; stop
               end if

               interSectAreas(1) = 0.5_wp * norm2( &
                    cross( interSectPoints(:,2) - comps(id_comp)%loc_points(:, &
                             comps(id_comp)%el(ie)%i_ver( nNodeInt(1,1) ) ) , &
                           interSectPoints(:,1) - comps(id_comp)%loc_points(:, &
                             comps(id_comp)%el(ie)%i_ver( nNodeInt(1,2) ) ) ) ) ! nNodes(2,iSec1) ) )
               interSectAreas(2) = comps(id_comp)%el(ie)%area - interSectAreas(1)

               if ( sec1_nVer .eq. 1 ) then
                 interSectCen(:,1) = ( interSectPoints(:,1) + interSectPoints(:,2) + & 
                      comps(id_comp)%loc_points(:, &
                            comps(id_comp)%el(ie)%i_ver(nNodeInt(1,1)) ) ) / 3.0_wp
               elseif ( sec1_nVer .eq. 2 ) then
                 interSectCen(:,1) = ( interSectPoints(:,1) + interSectPoints(:,2) + & 
                      comps(id_comp)%loc_points(:, &
                            comps(id_comp)%el(ie)%i_ver(nNodeInt(1,1)) ) + &
                      comps(id_comp)%loc_points(:, &
                            comps(id_comp)%el(ie)%i_ver(nNodeInt(1,2)) ) ) / 4.0_wp
               else
                 write(*,*) 
               end if
               interSectCen(:,2) = ( comps(id_comp)%el(ie)%area * comps(id_comp)%el(ie)%cen - &
                       interSectAreas(1) * interSectCen(:,1) ) / interSectAreas(2)

             else
               if ( sec2_nVer .gt. 2 ) then
                  write(*,*) ' error in sectional loads . stop ' ; stop
               end if

               interSectAreas(2) = 0.5_wp * norm2( &
                    cross( interSectPoints(:,2) - comps(id_comp)%loc_points(:, &
                             comps(id_comp)%el(ie)%i_ver( nNodeInt(2,1) ) ) , &
                           interSectPoints(:,1) - comps(id_comp)%loc_points(:, &
                             comps(id_comp)%el(ie)%i_ver( nNodeInt(2,2) ) ) ) ) ! nNodes(2,iSec1) ) )
               interSectAreas(1) = comps(id_comp)%el(ie)%area - interSectAreas(2)

               if ( sec2_nVer .eq. 1 ) then
                 interSectCen(:,2) = ( interSectPoints(:,1) + interSectPoints(:,2) + & 
                      comps(id_comp)%loc_points(:, &
                            comps(id_comp)%el(ie)%i_ver(nNodeInt(2,1)) ) ) / 3.0_wp
               elseif ( sec2_nVer .eq. 2 ) then
                 interSectCen(:,2) = ( interSectPoints(:,1) + interSectPoints(:,2) + & 
                      comps(id_comp)%loc_points(:, &
                            comps(id_comp)%el(ie)%i_ver(nNodeInt(2,1)) ) + &
                      comps(id_comp)%loc_points(:, &
                            comps(id_comp)%el(ie)%i_ver(nNodeInt(2,2)) ) ) / 4.0_wp
               else
                 write(*,*) 
               end if
               interSectCen(:,1) = ( comps(id_comp)%el(ie)%area * comps(id_comp)%el(ie)%cen - &
                       interSectAreas(2) * interSectCen(:,2) ) / interSectAreas(1)

             end if

             ! Update structures
             box_secloads(isec1)%nelems = box_secloads(isec1)%nelems + 1
             box_secloads(isec1)%elems(box_secloads(isec1)%nelems) = ie
             box_secloads(isec1)%fracs(box_secloads(isec1)%nelems) = interSectAreas(1) / comps(id_comp)%el(ie)%area 
             box_secloads(isec1)%cen(:,box_secloads(isec1)%nelems) = interSectCen(:,1) 
             box_secloads(isec2)%nelems = box_secloads(isec2)%nelems + 1
             box_secloads(isec2)%elems(box_secloads(isec2)%nelems) = ie
             box_secloads(isec2)%fracs(box_secloads(isec2)%nelems) = interSectAreas(2) / comps(id_comp)%el(ie)%area 
             box_secloads(isec2)%cen(:,box_secloads(isec2)%nelems) = interSectCen(:,2) 


           end if

        end if ! elems with at least one point s.t. nCoor \in (nCoorMix,nCoorMax)

      end if ! centre of the elems between the lateral faces of the box

    end do ! loop over elems

    ! check ----
    do is = 1 , n_sect
      write(*,*) ' box_secloads(',is,')%nelems = ' , box_secloads(is)%nelems
    end do

    ! Allocations and time loop ++++++++++++++++++++++++++++++++++++++++++++++++++

    ! Find the id of the reference where the loads must be projected -------------
    write(filename,'(A,I4.4,A)') trim(data_basename)//'_res_',an_start,'.h5'
    call open_hdf5_file(trim(filename),floc)
    call load_refs(floc,refs_R,refs_off) ! ,refs_G,refs_f,refs_tag)
    call close_hdf5_file(floc) 

!   ref_id = -333
!   do it = lbound(refs_tag,1) , ubound(refs_tag,1)
!     if ( stricmp(refs_tag(it),  comps(id_comp)%ref_tag) ) ref_id = it
!   end do
!   if ( ref_id .eq. -333 ) then 
!     write(*,*)
!     write(*,*) ' Available references systems: '
!     do it = lbound(refs_tag,1) , ubound(refs_tag,1)
!       write(*,*) ' ref_id : ' , it , ' ref_tag ' , trim(refs_tag(it))
!     end do
!     call error('dust_post','','Unknown ref.sys. defined for loads output.&
!          & Your input in dust_post.in is '//trim(ref_tag)//'. All the&
!          & available ref.sys. are listed above.')
!   end if
!   write(*,*) ' Available references systems: '
!   do it = lbound(refs_tag,1) , ubound(refs_tag,1)
!     write(*,*) ' ref_id : ' , it , ' ref_tag ' , trim(refs_tag(it))
!   end do
!   write(*,*) ' ref_id for force and moment projection : ' , ref_id

    ref_id = comps(id_comp)%ref_id
    !check
    write(*,'(A,I0,A)') ' ref. for force and moment projection (', &
         ref_id , ')  '//trim(comps(id_comp)%ref_tag)

    ! allocate tmp array to store the results --------------
    n_time = (an_end-an_start)/an_step + 1 ! int general eger division
    allocate( sec_loads(n_time,n_sect,n_loads) ) ; sec_loads = -333.0_wp ! initialisation for DEBUG
    allocate( time(n_time) ) ; time = -333.0_wp
    allocate( ref_mat(n_time,9) , off_mat(n_time,3) ) 
    ires = 0
    do it = an_start, an_end, an_step

      write(*,*) ' it : ' , it
      ires = ires + 1
      
      ! Open the file:
      write(filename,'(A,I4.4,A)') trim(data_basename)//'_res_',it,'.h5'
      call open_hdf5_file(trim(filename),floc)

      ! Load the references and move the points ---
      call load_refs(floc,refs_R,refs_off)
      ! Move the points ---------------------------
      call update_points_postpro(comps, points, refs_R, refs_off)
      ! Load the results --------------------------
      call load_res(floc, comps, vort, cp, t)

      call close_hdf5_file(floc)

      !TODO: projection step to local frame
      ! ...
      do is = 1 , n_sect

        ! from local coordinates to coordinates in the base ref.frame
        r_axis_bas(:,is) = matmul( refs_R(:,:,ref_id) , r_axis(:,is) ) + &
                           refs_off(:,ref_id)

        ! force
        F_bas = 0.0_wp ; M_bas = 0.0_wp 
        do ie = 1 , box_secloads(is)%nelems

          !from local coordinates to coordinates in the base ref.frame
          box_secloads_cen = matmul( refs_R(:,:,ref_id) , &
                                     box_secloads(is)%cen(:,ie) ) 
          
          F_bas1 = comps(id_comp)%el( box_secloads(is)%elems(ie) )%dforce * &
                                      box_secloads(is)%fracs(ie)
          F_bas  = F_bas + F_bas1
          !TODO: add moment (it requires axis computation)
          M_bas  = M_bas + cross( box_secloads_cen - &
                                  r_axis_bas(:,is) , F_bas1 )
        end do

        ! From global to local coordinates of forces and moments 
        sec_loads(ires,is,1:3) = matmul( &
             transpose( refs_R(:,:, ref_id) ) , F_bas ) / y_span(is)
        ! TODO: moment
        ! moment ( only the component around the <ax_coord> )
        M_bas = matmul( &
             transpose( refs_R(:,:, ref_id) ) , M_bas )
        sec_loads(ires,is,4) = M_bas(2) / y_span(is)
      end do

      ref_mat(ires,:) = reshape(refs_R(:,:,ref_id),(/ 9 /))
      off_mat(ires,:) = refs_off(:,ref_id)
      time(ires) = t

    end do

    select case(trim(out_frmt))
     case('dat')
      write(filename,'(A)') trim(basename)//'_'//trim(an_name) 
      call dat_out_sectional ( filename, y_cen, y_span, time, sec_loads, &
                              ref_mat , off_mat ) 
     case('tecplot')
      write(filename,'(A)') trim(basename)//'_'//trim(an_name)//'.plt' 
      call tec_out_sectional ( filename, time, sec_loads, y_cen, y_span ) 
    end select

    write(*,*) nl//nl//' end of sectional loads'//nl//nl
! ######################################################################    
! #      TO BE CHECKED       ###########################################    
! ######################################################################    


    ! destroy box_secloads structure
    do is = lbound(box_secloads,1) , ubound(box_secloads,1)
      deallocate(box_secloads(is)%elems)
      deallocate(box_secloads(is)%fracs)
    end do

    deallocate(r_axis,r_axis_bas)
    deallocate(nCoordSec,nCoordCen)
    deallocate(box_coord_tmp,box_coord)
    deallocate(y_cen,y_span) 
  
  case default
  !  call error()
    write(*,*) ' error. ' ; stop
  
  
  end select
   
  
  !TODO: move deallocate(comps) outside this routine,
  !      because it is common to all the analyses
  deallocate(comps,components_names)

  write(*,*) nl//' post_sectional done.'//nl

end subroutine post_sectional

! ---------------------------------------------------------------------- 




! ---------------------------------------------------------------------- 

end module mod_post_sectional
