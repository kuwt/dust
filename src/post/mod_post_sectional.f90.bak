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

use mod_math, only: &
  cross

implicit none

! type rof box sectional loads
type t_box_secloads
  integer :: nelems
  integer  , allocatable :: elems(:)
! real(wp) , allocatable :: element_fraction(:)
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
real(wp) , allocatable :: y_cen(:)
real(wp) , parameter   :: tol_y_cen = 1.0e-3_wp
real(wp) , allocatable :: r_axis(:,:) , r_axis_bas(:,:)

real(wp) :: F_bas(3) , F_bas1(3)
real(wp) :: M_bas(3)
integer :: ie , ic , it , i1

! sectional loads: box -------------------------------------
character(len=max_char_len) :: comp_input
integer :: n_box
logical :: reshape_box
!TODO: clean declarations
real(wp) :: nCoordMaxBox
real(wp) , allocatable :: box_coord(:,:) , box_coord_tmp(:,:) , box_coord_rotated(:,:)
real(wp) :: vVec(3) , bVec(3) , wVec(3) , baseLen(2) , heigLen(2) , spanLen
real(wp) :: nVec(3) ! , axis_mom(3)
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
    allocate(y_cen(n_sect)) 
    ie = 0  
    do is = 1 , n_sect
      ie = ie + 1
      y_cen(is) = sum(comps(id_comp)%loc_points(ax_coor,comps(id_comp)%el(ie)%i_ver)) / &
                  dble(comps(id_comp)%el(ie)%n_ver)
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
             transpose( refs_R(:,:, ref_id) ) , F_bas )
        ! moment ( only the component around the <ax_coord> )
        M_bas = matmul( &
             transpose( refs_R(:,:, ref_id) ) , M_bas )
        sec_loads(ires,is,4) = M_bas(ax_coor)
    
      end do ! loop over sections
    
      ref_mat(ires,:) = reshape(refs_R(:,:,ref_id),(/ 9 /))
      off_mat(ires,:) = refs_off(:,ref_id)
      time(ires) = t
    
    end do
    
    write(filename,'(A)') trim(basename)//'_'//trim(an_name) ! and then appen Fx,Fy,Fz,M
    write(*,*) ' +++++++ size(y_cen) : ' , size(y_cen)
    call dat_out_sectional ( filename , y_cen , time , sec_loads , &
                              ref_mat , off_mat ) 
    
    write(*,*) nl//nl//' end of sectional loads'//nl//nl
    
    deallocate(r_axis,r_axis_bas)
    deallocate(y_cen)
    deallocate(sec_loads,time,ref_mat,off_mat)
  
  
  case ( 'cgns' )

    !TODO: move to a subroutine
    
    ! Some assumptions ---------------
    id_comp = 1   ! 1. only one component is loaded
  
    ! Read input for box and build the 8 nodes of the box ----------------------------- 
    call getsuboption(sbprms,'BoxSect',bxprms)
    allocate(box_coord(8,3))
    box_coord(1,:) = getrealarray(bxprms,'refNode',3)
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

      ! Update box corners: update only if the box is too large 
      allocate(box_coord_tmp(8,3))
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

      deallocate(box_coord_tmp)
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
    ! Find the elements belonging to each section ---------------------------------------
    ! allocate and initialise
    allocate(box_secloads(n_sect))
    do is = 1 , n_sect
      box_secloads(is)%nelems = 0 
      allocate(box_secloads(is)%elems(size(comps(id_comp)%el))) ; box_secloads(is)%elems = -333
    end do

    do ie = 1 , size(comps(id_comp)%el)

      ! compute the centre of the element
      comps(id_comp)%el(ie)%cen =  &
        sum(comps(id_comp)%loc_points(:,comps(id_comp)%el(ie)%i_ver-comps(id_comp)%i_points(1)+1),2) / &
                  dble(comps(id_comp)%el(ie)%n_ver)
!     write(*,*) nl , ' comps(id_comp)%el(ie)%cen : ' , comps(id_comp)%el(ie)%cen , nl
      !TODO: check the vertices of the elements. 
      ! Up to now, only rough check on the centres of the elements
      do i1 = 1 , 4
        distance(i1) = sum ( ( comps(id_comp)%el(ie)%cen - refPoiLateralFaces(:,i1) ) * normalLateralFaces(:,i1) )
!       write(*,*) distance(i1)
      end do 

      if ( all( distance .gt. 0.0_wp ) ) then

        nCoord = sum ( ( comps(id_comp)%el(ie)%cen - refPoiLateralFaces(:,1) ) * nVec )

        if ( ( nCoord .ge. 0.0_wp ) .and. ( nCoord .le. nCoordMax-nCoordMin ) ) then
!         write(*,*) ' nCoord : ' , nCoord
          ! Find the slice and update the structure
          do is = 1 , n_sect
            if ( ( nCoord - nCoordSec(is)   .gt. 0.0_wp ) .and. & 
                 ( nCoord - nCoordSec(is+1) .lt. 0.0_wp ) ) then
              box_secloads(is)%nelems = box_secloads(is)%nelems + 1
              box_secloads(is)%elems(box_secloads(is)%nelems) = ie
            end if
          end do
        end if
      end if
    end do 

    ! check ----
    do is = 1 , n_sect
      write(*,*) ' box_secloads(',is,')%nelems = ' , box_secloads(is)%nelems
    end do

! ######################################################################    
! #      TO BE CHECKED       ###########################################    
! ######################################################################    


    ! destroy box_secloads structure
    do is = 1 , size(box_secloads)
      deallocate(box_secloads(is)%elems)
    end do

    deallocate(nCoordSec,nCoordCen)
    deallocate(box_coord)
  
  
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
