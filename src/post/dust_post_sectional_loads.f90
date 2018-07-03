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

program dust_post

use mod_param, only: &
  wp, nl, max_char_len, extended_char_len , pi

use mod_handling, only: &
  error, warning, info, printout, dust_time, t_realtime, new_file_unit

use mod_geometry, only: &
  t_geo, t_geo_component

use mod_basic_io, only: &
  read_mesh_basic, write_basic

use mod_aero_elements, only: &
  c_elem, t_elem_p!, t_vp

!this is for the parsing
use mod_parse, only: &
  t_parse, &
  getstr, getlogical, getreal, getint, &
  getrealarray, getintarray , &
  ignoredParameters, finalizeParameters, &
  countoption, getsuboption, t_link, check_opt_consistency, &
  print_parse_debug

use mod_build_geo, only: &
  build_geometry

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

use mod_stringtools, only: &
  LowCase, isInList, stricmp , strip_mult_appendix

use mod_geo_postpro, only: &
  load_components_postpro, update_points_postpro , prepare_geometry_postpro, &
  expand_actdisk_postpro, prepare_wake_postpro

use mod_wake_pan, only: &
  t_wake_panels

use mod_wake_ring, only: &
  t_wake_rings

use mod_tecplot_out, only: &
  tec_out_viz, tec_out_probes, tec_out_box, tec_out_loads

use mod_vtk_out, only: &
  vtk_out_viz , vtr_write

use mod_dat_out, only: & 
  dat_out_probes_header, & 
  dat_out_loads_header, &
  dat_out_sectional

use mod_math, only: &
  cross

use mod_actuatordisk, only: &
  t_actdisk

implicit none

type t_box_secloads
  integer :: nelems
  integer  , allocatable :: elems(:)
  real(wp) , allocatable :: fracs(:)
end type t_box_secloads


!Input
character(len=*), parameter :: input_file_name_def = 'dust_post.in' 
character(len=max_char_len) :: input_file_name

!Geometry parameters
type(t_parse) :: prms
type(t_parse), pointer :: sbprms , bxprms

integer :: n_analyses, ia

character(len=max_char_len) :: basename, data_basename
character(len=max_char_len) :: an_name, an_type
integer :: an_start, an_end, an_step, nstep
integer :: n_comp, i_comp
character(len=max_char_len), allocatable :: components_names(:)
integer :: n_var, i_var
character(len=max_char_len), allocatable :: var_names(:)
character(len=max_char_len) :: lowstr
character(len=max_char_len) :: filename
character(len=max_char_len) :: out_frmt
logical :: all_comp
logical :: out_vort, out_vel, out_cp, out_press
logical :: out_wake

integer(h5loc) :: floc, geo_floc, gloc1, gloc2 , ploc

real(wp), allocatable :: points(:,:), points_exp(:,:)
integer, allocatable :: elems(:,:)
type(t_geo_component), allocatable :: comps(:)
integer :: nelem, nelem_out

character(len=max_char_len) , allocatable :: refs_tag(:)
real(wp), allocatable :: refs_R(:,:,:), refs_off(:,:)
real(wp), allocatable :: refs_G(:,:,:), refs_f(:,:)
real(wp), allocatable :: vort(:), cp(:)
real(wp), allocatable :: wvort(:), wvort_pan(:,:), wvort_rin(:,:)
real(wp), allocatable :: wpoints(:,:), wpoints_pan(:,:,:), wpoints_rin(:,:,:)
!real(wp), allocatable :: wcen(:,:,:)
integer,  allocatable :: wconn(:)

integer,  allocatable :: welems(:,:)
integer, allocatable  :: wstart(:,:)
integer :: nelem_w
real(wp) :: t

real(wp), allocatable :: print_vars(:,:)
character(len=max_char_len), allocatable :: print_var_names(:)
real(wp), allocatable :: print_vars_w(:,:)
character(len=max_char_len), allocatable :: print_var_names_w(:)
integer :: nprint, ivar

integer, allocatable :: print_elems(:,:)

! wake ------------
type(t_wake_panels) :: wake_pan
type(t_wake_rings)  :: wake_rin
type(t_elem_p), allocatable :: wake_elems(:)

! probe output ----
real(wp), allocatable :: probe_vars(:,:,:)
real(wp), allocatable :: time(:)
character(len=max_char_len), allocatable :: probe_var_names(:)
character(len=max_char_len), allocatable :: probe_loc_names(:)
character(len=max_char_len) :: in_type , str_a , filename_in , var_name
integer :: n_probes , n_vars , n_vars_int
real(wp), allocatable :: rr_probes(:,:)
logical :: probe_vel , probe_p , probe_vort
character(len=max_char_len) :: vars_str
real(wp) :: u_inf(3)
real(wp) :: P_inf , rho
real(wp) :: vel_probe(3) = 0.0_wp , vort_probe(3) = 0.0_wp 
real(wp) :: v(3) = 0.0_wp , w(3) = 0.0_wp
real(wp), allocatable , target :: sol(:) 
real(wp) :: pres_probe
integer :: fid_out , ip , ie, ic
! flow field ------
integer :: nxyz(3)
real(wp):: minxyz(3) , maxxyz(3)
real(wp), allocatable :: xbox(:) , ybox(:) , zbox(:)
real(wp) :: dxbox , dybox , dzbox
real(wp), allocatable :: box_vel(:,:) , box_p(:) , box_vort(:,:)
integer :: ix , iy , iz
integer , allocatable :: vars_n(:)
real(wp), allocatable :: vars(:,:) 
integer :: i_vars , i_var_v , i_var_p , i_var_w
! loads -----------
integer :: n_comps_meas
integer ,allocatable :: i_comps_meas(:)
character(len=max_char_len), allocatable :: comps_meas(:)
character(len=max_char_len) :: ref_tag
integer                     :: ref_id
real(wp) :: F_loc(3) , F_ref(3) , F_bas(3) , F_bas1(3)
real(wp) :: M_loc(3) , M_ref(3) , M_bas(3)
real(wp), allocatable :: force(:,:), moment(:,:)
integer :: ic2
real(wp), allocatable , target :: sol_p(:) 
! sectional loads ---
real(wp) :: axis_dir(3) , axis_nod(3) 
character(len=max_char_len), allocatable :: components_names_tmp(:)
character(len=max_char_len) :: comp_name_stripped
integer :: id_comp , nelem_span , nelem_chor , n_sect , n_time
real(wp) , allocatable :: sec_loads(:,:,:)
integer :: is 
integer :: n_loads = 4   ! F and moment around an axis
real(wp) , allocatable :: ref_mat(:,:) , off_mat(:,:)
real(wp) , allocatable :: y_cen(:)
real(wp) , parameter :: tol_y_cen = 1.0e-3_wp
real(wp) , allocatable :: r_axis(:,:) , r_axis_bas(:,:)
! sectional loads: box ---
! old character(len=max_char_len) :: box_filen_str
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
real(wp) , parameter :: nCoord_relOff = 0.01_wp ! 0.05_wp
real(wp) :: nCoordVert(4)   ! TRI o QUAD elements
integer  :: secVert(4)   ! TRI o QUAD elements
real(wp) :: interSectPoints(3,2)
real(wp) :: interSectAreas(2)
integer  :: nInterSect , index2
real(wp) :: node1(3) , node2(3)
integer  :: iSec1 , iSec2 , sec1_nVer , sec2_nVer
integer  :: iSecWork , iSecOth
integer  :: nNodeInt(2,2)


integer :: it , i1, ires , iv
integer :: ierr


call printout(nl//'>>>>>> DUST POSTPROCESSOR beginning >>>>>>'//nl)
call initialize_hdf5()

!------ Input reading ------

if(command_argument_count().gt.0) then                                         
  call get_command_argument(1,value=input_file_name)                           
else                                                                           
  input_file_name = input_file_name_def                                        
endif   

call prms%CreateStringOption('basename','Base name of the processed data')
call prms%CreateStringOption('data_basename','Base name of the data to be &
                              &processed')
call prms%CreateSubOption('Analysis','Definition of the motion of a frame', &
                          sbprms, multiple=.true.)
call sbprms%CreateStringOption('Type','type of analysis')
call sbprms%CreateStringOption('Name','specification of the analysis')
call sbprms%CreateIntOption('StartRes', 'Starting result of the analysis')
call sbprms%CreateIntOption('EndRes', 'Final result of the analysis')
call sbprms%CreateIntOption('StepRes', 'Result stride of the analysis')
call sbprms%CreateLogicalOption('Wake', 'Output also the wake for &
                                &visualization','T')
call sbprms%CreateStringOption('Format','Output format')
call sbprms%CreateStringOption('Component','Component to analyse', &
                               multiple=.true.)
call sbprms%CreateStringOption('Variable','Variables to be saved: velocity, pressure or&
                              & vorticity', multiple=.true.)

! probe output -------------
call sbprms%CreateStringOption('InputType','How to specify probe coordinates',&
                              multiple=.true.)
call sbprms%CreateRealArrayOption('Point','Point coordinates in dust_post.in',&
                              multiple=.true.)
call sbprms%CreateStringOption('File','File containing the coordinates of the probes',&
                              multiple=.true.)
! flow field output --------
call sbprms%CreateIntArrayOption( 'Nxyz','number of points per coordinate',&
                              multiple=.true.)
call sbprms%CreateRealArrayOption('Minxyz','lower bounds of the box',&
                              multiple=.true.)
call sbprms%CreateRealArrayOption('Maxxyz','upper bounds of the box',&
                              multiple=.true.)

! loads --------------------
call sbprms%CreateStringOption('CompName','Components where loads are computed',&
                              multiple=.true.)
call sbprms%CreateStringOption('Reference_Tag','Reference frame where loads&
                            & are computed',multiple=.true.)
! sectional loads ----------
call sbprms%CreateRealArrayOption('AxisDir','Direction of the axis defined the reference&
                            & points for sectional loads analisys', multiple=.true.)
call sbprms%CreateRealArrayOption('AxisNod','Node belonging to the axis used for sectional&
                            & loads analisys', multiple=.true.)
! sectional loads: box definition ---------
call sbprms%CreateSubOption('BoxSect','Definition of the box for sectional loads', &
                          bxprms)
call bxprms%CreateRealArrayOption('refNode','reference node to build the box')
call bxprms%CreateRealArrayOption('faceVec','vector identifying the direction of the base side &
                           &of the sections')
call bxprms%CreateRealArrayOption('faceBas','dimension along faceVec of the first and last sections')
call bxprms%CreateRealArrayOption('faceHei','dimension orthogonal to faceVec of the first and last sections')
call bxprms%CreateRealArrayOption('spanVec','vector defining the out-of plane direction of the box')
call bxprms%CreateRealOption('spanLen','dimension along the spanVec direction')
call bxprms%CreateIntOption('numSect','number of sections')
call sbprms%CreateRealArrayOption('AxisMom','axis for the computation of the moment. Perpendicular to sections')
! BoxSect = {
!  refNode = (/ -0.5 , 0.0 , -0.3 /)
!  faceVec = (/ 1.0 , 0.0 , 0.0 /)
!  faceBas = (/ 2.0 , 1.0 /) 
!  faceHei = (/ 1.0 , 1.0 /)
!  spanVec = (/ 0.0 , 1.0 , 0.0 /)
!  spanLen = 3.0 
! }
! AxisMom  = (/ 0.0 , 1.0 , 0.0 /) 
! old format
! call sbprms%CreateStringOption('BoxSect','Boxes delimiting the component for sectional&
!                             & loads analisys',multiple=.true.)
! call sbprms%CreateRealArrayOption('AxisMom','Axis for the computation of the moment&
!                             & in sectional loads analysis',multiple=.true.)


sbprms=>null()

call prms%read_options(input_file_name, printout_val=.false.)

basename = getstr(prms,'basename')
data_basename = getstr(prms,'data_basename')
n_analyses = countoption(prms,'Analysis')

!Cycle on all the analyses
do ia = 1,n_analyses
  
  !Get some of the options
  call getsuboption(prms,'Analysis',sbprms)
  an_type  = getstr(sbprms,'Type')

  !OUTPUT
  write(*,*) ' Analysis      : ' , ia , ' / ' , n_analyses
  write(*,*) ' Analysis Type : ' , trim(an_type) 

  call LowCase(an_type)
  an_name  = getstr(sbprms,'Name')
  out_frmt = getstr(sbprms,'Format') 
  call LowCase(out_frmt)
  an_start = getint(sbprms,'StartRes')
  an_end   = getint(sbprms,'EndRes')
  an_step  = getint(sbprms,'StepRes')
 
  !Check if we are analysing all the components or just some
  all_comp = .false.
  n_comp = countoption(sbprms, 'Component')
  if (n_comp .eq. 0)  then
    all_comp = .true.
  else
    allocate(components_names(n_comp))
    do i_comp = 1, n_comp
      components_names(i_comp) = getstr(sbprms, 'Component')
    enddo
    call LowCase(components_names(1),lowstr)    ! char 
    if(trim(lowstr) .eq. 'all') then
      all_comp = .true.
    endif
  endif

  !Fork the different kind of analyses
  select case(trim(an_type))

   !/////////////// Sectional Loads \\\\\\\\\\\\\\\\\
   case('sectional_loads')

    ! Some warnings and errors -------------------------------
    !WARNING: up to now, sectional loads are computed for one components only at time.
    if ( n_comp .le. 0 ) then
      call error('dust_post','','No component specified for ''sectional_loads''&
           & analysis. STOP')
    else if ( n_comp .ge. 2 ) then
      call warning('dust_post','','WARNING. More than one component specified &
           &for ''sectional_loads'' analysis: just the first one is considered. &
           &Please run another ''sectional_analysis'' if you need it on more than &
           &component.')
    end if

    ! load the geo components just once -----------------------
    call open_hdf5_file(trim(data_basename)//'_geo.h5', floc)
    !TODO: here get the run id

    !Strip the appendix of multiple components, to load all the multiple
    !components at once
    call strip_mult_appendix(components_names(1), comp_name_stripped, '__') 
    allocate(components_names_tmp(1)) ; components_names_tmp = trim(comp_name_stripped) ! components_names(1:1)
    call load_components_postpro(comps, points, nelem, floc, & 
                                 components_names_tmp,  all_comp)
    call close_hdf5_file(floc)
 
    !TODO: Check if the components in dust_post.in really exist
    !DEBUG 
    write(*,*) nl//' comps(:)%comp_name : '
    do i_comp = 1 , size(comps)
      write(*,*) trim(comps(i_comp)%comp_name)
    end do
    write(*,*) nl//' components_names   : '
    do i_comp = 1 , size(components_names)
      write(*,*) trim(components_names(i_comp))
    end do
    ! Prepare_geometry_postpro
    call prepare_geometry_postpro(comps)

    ! Read the axis for the computation of the sectional loads
    axis_dir = getrealarray(sbprms,'AxisDir',3)
    axis_nod = getrealarray(sbprms,'AxisNod',3)

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
    !TODO: work on postpro for rotors to allow separated output for all the blades
    !here: some workaround to cope with rotors.
    !      Get the actual Id of the first defined component, for sectional loads
    if ( .not. allocated(components_names) ) then
      call error('dust_post','','components_names NOT allocated. Something went&
           & wrong. STOP')
    end if
    id_comp = -333
    do i_comp = 1 , size(comps)
      if ( trim(comps(i_comp)%comp_name) .eq. trim(components_names(1)) ) then
        id_comp = i_comp
      end if
    end do

    ! check if the component really exists
    if ( id_comp .eq. -333 ) then
      do i_comp = 1 , size(comps)
        write(*,*) ' comps(i_comp)%comp_name ' , trim(comps(i_comp)%comp_name)
      end do
      call error('dust_post','','No valid component defined. STOP')
    end if
    
!   !DEBUG
!   write(*,*) ' size(comps) : ' , size(comps) , ' , id_comp '
!   write(*,*) ' trim(comps(id_comp)%comp_name) : ' , trim(comps(id_comp)%comp_name)
!   write(*,*) ' trim(comp_name_stripped)     ) : ' , trim(comp_name_stripped)
    if ( ( size(comps) .gt. 1 ) .and. &
       ( trim(comps(id_comp)%comp_name) .eq. trim(comp_name_stripped) ) ) then
      comps(id_comp)%ref_tag = trim(comps(id_comp)%ref_tag)//'__01'
      components_names(1) = trim(components_names(1))//'__01' 
!     write(*,*) trim(components_names(1))
    end if


    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

    select case( trim(comps(1)%comp_input) )

      !TODO: then move it to 'default', or add some logic ('Box') ??
      case( 'parametric' )

! ! sectional loads: box ---
! character(len=max_char_len) :: box_filen_str
! real(wp) , allocatable :: box_coord(8,3)
! real(wp) , allocatable :: axis_mom(3)

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
        nVec = getrealarray(sbprms,'AxisMom',3) 

        ! Normalise
        vVec = vVec / norm2(vVec)
        bVec = bVec / norm2(bVec)
        nVec = nVec / norm2(nVec)

        wVec = cross(vVec,nVec) ; wVec = wVec / norm2(wVec)
    
        box_coord(2,:) = box_coord(1,:) + baseLen(1) * vVec  
        box_coord(3,:) = box_coord(2,:) + heigLen(1) * wVec  
        box_coord(4,:) = box_coord(1,:) + heigLen(1) * wVec  

        box_coord(5,:) = box_coord(1,:) + spanLen * bVec
        box_coord(6,:) = box_coord(5,:) + baseLen(2) * vVec  
        box_coord(7,:) = box_coord(6,:) + heigLen(2) * wVec  
        box_coord(8,:) = box_coord(5,:) + heigLen(2) * wVec  


        ! Find if the whole component is contained in the box and/or reshape the box -----
        ! ( in the direction perpendicolar to the nVec = axis_mom )
        ie = 1
        nCoordMax = sum ( ( comps(id_comp)%loc_points( :, ie ) &
         - box_coord(1,:) ) * nVec )
        nCoordMin = nCoordMax !0.0_wp
        do ie = 2 , size( comps(id_comp)%loc_points , 2 ) ! size(comps(id_comp)%el)
          nCoord = sum ( ( comps(id_comp)%loc_points( :, ie ) &
            - box_coord(1,:) ) * nVec )
!         write(*,*) ' nCoord : ' , nCoord
          if ( nCoord .gt. nCoordMax ) nCoordMax = nCoord
          if ( nCoord .lt. nCoordMin ) nCoordMin = nCoord 
        end do

        !CHECK
        write(*,*) nl//' Before the offset'//nl//' nCoordMin : ' , nCoordMin
        write(*,*) 'nCoordMax : ' , nCoordMax , nl

        ! add the Offset
        dnCoord = ( nCoordMax - nCoordMin ) / dble(n_sect) 
        nCoordMin = nCoordMin - dnCoord * nCoord_relOff
        nCoordMax = nCoordMax + dnCoord * nCoord_relOff
        dnCoord = ( nCoordMax - nCoordMin ) / dble(n_sect) 

        ! Update box corners: TODO: add an IF statement: update only if the box is too large 
        allocate(box_coord_tmp(8,3))
        box_coord_tmp(1,:) = box_coord(1,:) + ( box_coord(5,:) - box_coord(1,:) ) * nCoordMin / spanLen
        box_coord_tmp(2,:) = box_coord(2,:) + ( box_coord(6,:) - box_coord(2,:) ) * nCoordMin / spanLen
        box_coord_tmp(3,:) = box_coord(3,:) + ( box_coord(7,:) - box_coord(3,:) ) * nCoordMin / spanLen
        box_coord_tmp(4,:) = box_coord(4,:) + ( box_coord(8,:) - box_coord(4,:) ) * nCoordMin / spanLen
        box_coord_tmp(5,:) = box_coord(1,:) + ( box_coord(5,:) - box_coord(1,:) ) * nCoordMax / spanLen
        box_coord_tmp(6,:) = box_coord(2,:) + ( box_coord(6,:) - box_coord(2,:) ) * nCoordMax / spanLen
        box_coord_tmp(7,:) = box_coord(3,:) + ( box_coord(7,:) - box_coord(3,:) ) * nCoordMax / spanLen
        box_coord_tmp(8,:) = box_coord(4,:) + ( box_coord(8,:) - box_coord(4,:) ) * nCoordMax / spanLen

        box_coord = box_coord_tmp
        !CHECK
        write(*,*) ' Updated box_coord : ' 
        do i1 = 1 ,8 
          write(*,*) box_coord(i1,:)
        end do

        ! Find the n-coord of the sections and the reference points at the centre of the sections
        allocate(nCoordSec(n_sect+1))
        allocate(nCoordCen(3,n_sect)) ; nCoordCen = 0.0_wp
        
        
        nCoordSec(1) = 0.0_wp ! + nCoordMin
!       write(*,*)  ' nCoordSec(',1,') : ' , nCoordSec(1)
        do is = 1 , n_sect
          nCoordSec(is+1) = (nCoordMax-nCoordMin) * dble(is) / dble(n_sect) ! + nCoordMin
          nCoordCen(2,is) = (nCoordMax-nCoordMin) * ( dble(is) - 0.5_wp)  / dble(n_sect) ! + nCoordMin
!         write(*,*)  '   nCoordCen(:,',is  ,') : ' , nCoordCen(:,is)
!         write(*,*)  ' nCoordSec(',is+1,') : ' , nCoordSec(is+1)
        end do

        !TODO: check next formula
        allocate(y_cen(n_sect)) ; y_cen = nCoordCen(2,:)

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

!       !CHECK
!       write(*,*) ' normalLateralFaces : '
!       do i1 = 1 , 4
!         write(*,*) normalLateralFaces(:,i1)
!       end do
        
        ! Find the elements belonging to each section ---------------------------------------
        ! allocate and initialise
        allocate(box_secloads(n_sect))
        do is = 1 , n_sect
          box_secloads(is)%nelems = 0 
          allocate(box_secloads(is)%elems(size(comps(id_comp)%el))) ; box_secloads(is)%elems = -333
          allocate(box_secloads(is)%fracs(size(comps(id_comp)%el))) ; box_secloads(is)%fracs = -333.3_wp
        end do

        do ie = 1 , size(comps(id_comp)%el)

          ! Reset some useful variables 
          nCoordVert = 0.0_wp
          secVert = 0
          interSectPoints = 0.0_wp
          interSectAreas  = 0.0_wp
          nInterSect = 0
          nNodeInt = 0

          ! compute the centre and the area of the element
          comps(id_comp)%el(ie)%cen =  &
            sum(comps(id_comp)%loc_points(:,comps(id_comp)%el(ie)%i_ver-comps(id_comp)%i_points(1)+1),2) / &
                      dble(comps(id_comp)%el(ie)%n_ver)
          comps(id_comp)%el(ie)%area =  0.5_wp * norm2( &
                      cross ( comps(id_comp)%loc_points(:,comps(id_comp)%el(ie)%i_ver(3)-comps(id_comp)%i_points(1)+1) - & 
                              comps(id_comp)%loc_points(:,comps(id_comp)%el(ie)%i_ver(1)-comps(id_comp)%i_points(1)+1) , & 
                              comps(id_comp)%loc_points(:,comps(id_comp)%el(ie)%i_ver(2)-comps(id_comp)%i_points(1)+1) - & 
                              comps(id_comp)%loc_points(:,comps(id_comp)%el(ie)%i_ver( & 
                                comps(id_comp)%el(ie)%n_ver )-comps(id_comp)%i_points(1)+1) ) )
!         write(*,*) nl , ' comps(id_comp)%el(ie)%cen : ' , comps(id_comp)%el(ie)%cen , nl

          ! check if the centre of the element belongs to the box 
          do i1 = 1 , 4
            distance(i1) = sum ( ( comps(id_comp)%el(ie)%cen - refPoiLateralFaces(:,i1) ) * normalLateralFaces(:,i1) )
!           write(*,*) distance(i1)
          end do 

          ! if the centre of the element belongs to the box -> compute the slice contributions
          if ( all( distance .gt. 0.0_wp ) ) then

            ! nCoord of the centre 
            nCoord = sum ( ( comps(id_comp)%el(ie)%cen - refPoiLateralFaces(:,1) ) * nVec )

            ! nCoord of the vertices
            do iv = 1 , comps(id_comp)%el(ie)%n_ver 
              nCoordVert(iv) = &
                 sum ( ( comps(id_comp)%loc_points(:, comps(id_comp)%el(ie)%i_ver(iv)      &
                                                     -comps(id_comp)%i_points(1)+1       ) &
                 - refPoiLateralFaces(:,1) ) * nVec )
            end do


            ! TODO: change this loop in order to deal with extreme sections of the box.
            ! if statements looks fine. check the contents
            if ( ( maxval(nCoordVert) .ge. 0.0_wp ) .and. &
                 ( minval(nCoordVert) .le. nCoordMax-nCoordMin ) ) then
               do iv = 1 , comps(id_comp)%el(ie)%n_ver   
                 do is = 1 , n_sect
                   if ( ( nCoordVert(iv) - nCoordSec(is)   .gt. 0.0_wp ) .and. & 
                        ( nCoordVert(iv) - nCoordSec(is+1) .lt. 0.0_wp ) ) then
                     secVert(iv) = is
                   end if 
                 end do
               end do

               ! check if the elements belong to less than three sections --> otherwise ERROR
               if ( maxval(secVert(1:comps(id_comp)%el(ie)%n_ver) ) & 
                  - minval(secVert(1:comps(id_comp)%el(ie)%n_ver) ) .gt. 1 ) then
                 write(*,*) ' comps(id_comp)%comp_name : ' , comps(id_comp)%comp_name 
                 write(*,*) ' ie : ' , ie  
                 write(*,*) ' secVert : ' , secVert(1:comps(id_comp)%el(ie)%n_ver)
                 call error('dust_post','','The element above belongs to more than two&
                      & sections in sectional_loads analysis. STOP')
               end if
               if ( all(  secVert .eq. secVert(1) ) ) then
                 box_secloads(secVert(1))%nelems = box_secloads(secVert(1))%nelems + 1
                 box_secloads(secVert(1))%elems(box_secloads(secVert(1))%nelems) = ie
                 box_secloads(secVert(1))%fracs(box_secloads(secVert(1))%nelems) = 1.0_wp
               else
!                write(*,*) ' ie : ' , ie  
!                write(*,*) ' secVert : ' , secVert(1:comps(id_comp)%el(ie)%n_ver)
!                write(*,*)


                 ! Find intersections with the plane delimiting the sections --------------
                 do iv = 1 , comps(id_comp)%el(ie)%n_ver
                   !TODO: treat elements partially belonging to the box.  if ( ... )
                   ! ... 
                   index2 = mod(iv,comps(id_comp)%el(ie)%n_ver)+1
                   if ( secVert( index2 ) .ne. secVert(iv) ) then
                     nInterSect = nInterSect + 1
                     node1 = comps(id_comp)%loc_points(:, comps(id_comp)%el(ie)%i_ver(iv) &
                                                         -comps(id_comp)%i_points(1)+1       ) 
                     node2 = comps(id_comp)%loc_points(:, comps(id_comp)%el(ie)%i_ver(index2) &
                                                         -comps(id_comp)%i_points(1)+1       ) 
                     interSectPoints(:,nInterSect) = node1 + (node2-node1) * &
                          ( nCoordSec( max( secVert(iv), secVert(index2) ) ) - nCoordVert(iv) ) / &
                          ( nCoordVert(index2) - nCoordVert(iv) )

!                    !DEBUG
!                    do i1 = 1 , comps(id_comp)%el(ie)%n_ver
!                      write(*,*) ' i1 : ' , i1 ,' . ' , &
!                       comps(id_comp)%loc_points(:, comps(id_comp)%el(ie)%i_ver(i1) &
!                                                          -comps(id_comp)%i_points(1)+1       ) 
!                    end do
!                    write(*,*) ' i , node1 : ' , iv , node1
!                    write(*,*) ' i , interS: ' , nInterSect , interSectPoints(:,nInterSect) 
!                    write(*,*) ' i , node2 : ' , index2 , node2
                     if ( secVert(iv) .eq. secVert(1) ) then ! sec.1
                       nNodeInt(1,nInterSect) = iv
                       nNodeInt(2,nInterSect) = index2
                     else ! sec.2
                       nNodeInt(1,nInterSect) = index2
                       nNodeInt(2,nInterSect) = iv
                     end if
                   end if
                 end do
!                write(*,*) nNodeInt(1,:)
!                write(*,*) nNodeInt(2,:)
!                write(*,*)

                 iSec1 = secVert(1) ; sec1_nVer = 1 ; sec2_nVer = 0
                 do iv = 2 , comps(id_comp)%el(ie)%n_ver
                   if ( secVert(iv) .ne. iSec1 ) then
                     iSec2 = secVert(iv) ; sec2_nVer = sec2_nVer + 1
                   else
                     sec1_nVer = sec1_nVer + 1
                   end if
                 end do

                 if ( sec1_nVer .le. sec2_nVer ) then
                   if ( sec1_nVer .gt. 2 ) then
                      write(*,*) ' error in sectional loads . stop ' ; stop
                   end if

                   interSectAreas(1) = 0.5_wp * norm2( &
                        cross( interSectPoints(:,2) - comps(id_comp)%loc_points(:, &
                                 comps(id_comp)%el(ie)%i_ver( nNodeInt(1,1) )-comps(id_comp)%i_points(1)+1 ) , &
                               interSectPoints(:,1) - comps(id_comp)%loc_points(:, &
                                 comps(id_comp)%el(ie)%i_ver( nNodeInt(1,2) )-comps(id_comp)%i_points(1)+1 ) ) ) ! nNodes(2,iSec1) ) )
                   interSectAreas(2) = comps(id_comp)%el(ie)%area - interSectAreas(1)
!                  write(*,*) ' el%area : ' , comps(id_comp)%el(ie)%area 
!                  write(*,*) ' interSectArea([1,2]) : ' , interSectAreas(1) , interSectAreas(2)
                 else
                   if ( sec2_nVer .gt. 2 ) then
                      write(*,*) ' error in sectional loads . stop ' ; stop
                   end if

                   interSectAreas(2) = 0.5_wp * norm2( &
                        cross( interSectPoints(:,2) - comps(id_comp)%loc_points(:, &
                                 comps(id_comp)%el(ie)%i_ver( nNodeInt(2,1) )-comps(id_comp)%i_points(1)+1 ) , &
                               interSectPoints(:,1) - comps(id_comp)%loc_points(:, &
                                 comps(id_comp)%el(ie)%i_ver( nNodeInt(2,2) )-comps(id_comp)%i_points(1)+1 ) ) ) ! nNodes(2,iSec1) ) )
                   interSectAreas(1) = comps(id_comp)%el(ie)%area - interSectAreas(2)
!                  write(*,*) ' el%area : ' , comps(id_comp)%el(ie)%area 
!                  write(*,*) ' interSectArea([2,1]) : ' , interSectAreas(2) , interSectAreas(1)

                 end if

                 ! Update structures
                 box_secloads(isec1)%nelems = box_secloads(isec1)%nelems + 1
                 box_secloads(isec1)%elems(box_secloads(isec1)%nelems) = ie
                 box_secloads(isec1)%fracs(box_secloads(isec1)%nelems) = interSectAreas(1) / comps(id_comp)%el(ie)%area 
                 box_secloads(isec2)%nelems = box_secloads(isec2)%nelems + 1
                 box_secloads(isec2)%elems(box_secloads(isec2)%nelems) = ie
                 box_secloads(isec2)%fracs(box_secloads(isec2)%nelems) = interSectAreas(2) / comps(id_comp)%el(ie)%area 


               end if

            end if

!           if ( ( nCoord .ge. 0.0_wp ) .and. ( nCoord .le. nCoordMax-nCoordMin ) ) then
!             ! Find the slice and update the structure
!             do is = 1 , n_sect
!               if ( ( nCoord - nCoordSec(is)   .gt. 0.0_wp ) .and. & 
!                    ( nCoord - nCoordSec(is+1) .lt. 0.0_wp ) ) then
!                 box_secloads(is)%nelems = box_secloads(is)%nelems + 1
!                 box_secloads(is)%elems(box_secloads(is)%nelems) = ie
!               end if
!             end do
!           end if
          end if
        end do 

        ! check ----
        do is = 1 , n_sect
          write(*,*) ' box_secloads(',is,')%nelems = ' , box_secloads(is)%nelems
          do ie = 1 , box_secloads(is)%nelems
            write(*,*) ' ie , frac : ' , box_secloads(is)%elems(ie) , box_secloads(is)%fracs(ie) 
          end do
        end do

        ! Allocations and time loop ++++++++++++++++++++++++++++++++++++++++++++++++++

        ! Find the id of the reference where the loads must be projected -------------
        write(filename,'(A,I4.4,A)') trim(data_basename)//'_res_',an_start,'.h5'
        call open_hdf5_file(trim(filename),floc)
        call load_refs(floc,refs_R,refs_off,refs_G,refs_f,refs_tag)
        call close_hdf5_file(floc) 

        ref_id = -333
        do it = lbound(refs_tag,1) , ubound(refs_tag,1)
          if ( stricmp(refs_tag(it),  comps(id_comp)%ref_tag) ) ref_id = it
        end do
        if ( ref_id .eq. -333 ) then 
          write(*,*)
          write(*,*) ' Available references systems: '
          do it = lbound(refs_tag,1) , ubound(refs_tag,1)
            write(*,*) ' ref_id : ' , it , ' ref_tag ' , trim(refs_tag(it))
          end do
          call error('dust_post','','Unknown ref.sys. defined for loads output.&
               & Your input in dust_post.in is '//trim(ref_tag)//'. All the&
               & available ref.sys. are listed above.')
        end if
        write(*,*) ' Available references systems: '
        do it = lbound(refs_tag,1) , ubound(refs_tag,1)
          write(*,*) ' ref_id : ' , it , ' ref_tag ' , trim(refs_tag(it))
        end do
        write(*,*) ' ref_id for force and moment projection : ' , ref_id

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
            ! force
            F_bas = 0.0_wp ; M_bas = 0.0_wp 
            do ie = 1 , box_secloads(is)%nelems
              F_bas1 = comps(id_comp)%el( box_secloads(is)%elems(ie) )%dforce * &
                                          box_secloads(is)%fracs(ie)
              F_bas  = F_bas + F_bas1
              !TODO: add moment (it requires axis computation)
              ! ...
            end do
!           !CHECK
!           write(*,*) ' is , F_bas ' , is , F_bas
            ! force
            sec_loads(ires,is,1:3) = matmul( &
                 transpose( refs_R(:,:, ref_id) ) , F_bas )
            ! TODO: moment
            sec_loads(ires,is,4) = 0.0_wp
          end do

          ref_mat(ires,:) = reshape(refs_R(:,:,ref_id),(/ 9 /))
          off_mat(ires,:) = refs_off(:,ref_id)
          time(ires) = t

        end do

        write(filename,'(A)') trim(basename)//'_'//trim(an_name) ! and then appen Fx,Fy,Fz,M
        write(*,*) ' +++++++ size(y_cen) : ' , size(y_cen)
        call dat_out_sectional ( filename , y_cen , time , sec_loads , &
                                  ref_mat , off_mat ) 

        write(*,*) nl//nl//' end of sectional loads'//nl//nl


        deallocate(y_cen)
        ! destroy box_secloads structure
        do is = 1 , size(box_secloads)
          deallocate(box_secloads(is)%elems)
          deallocate(box_secloads(is)%fracs)
        end do

        deallocate(nCoordSec,nCoordCen)

        deallocate(box_coord_tmp)
        deallocate(box_coord)


!     case( 'parametric' )
      case default

        ! find the sections in the structured mesh --------------------------------------- 
        nelem_span = comps(id_comp)%parametric_nelems_span
        nelem_chor = comps(id_comp)%parametric_nelems_chor
        
        n_sect = nelem_span ! n_span = n. sections
        
        !TODO: find y-coord is general enough for all the 'parametric' components
        ! Find the y-coord of the centre of the slices -----------------------------------
        allocate(y_cen(n_sect)) 
        write(*,*) ' +++++++ n_sect : ' , n_sect
        ie = 0  
        do is = 1 , n_sect
          ie = ie + 1
          y_cen(is) = sum(comps(id_comp)%loc_points(2,comps(id_comp)%el(ie)%i_ver-comps(id_comp)%i_points(1)+1)) / &
                      dble(comps(id_comp)%el(ie)%n_ver)
          do ic = 2 , nelem_chor
            ie = ie + 1
            if ( abs( y_cen(is) - & !comps(id_comp)%el(ie)%cen(2) ) .gt. tol_y_cen ) then
                   sum(comps(id_comp)%loc_points(2,comps(id_comp)%el(ie)%i_ver-comps(id_comp)%i_points(1)+1)) / &
                   dble(comps(id_comp)%el(ie)%n_ver) ) .gt. tol_y_cen ) then
              call error('dust_post','','Wrong section definition for component '//&
                   trim(comps(id_comp)%comp_name)//' for ''sectional_loads'' analysis. STOP')
            end if
          end do
        end do
        ! Find the coordinate of the reference points on the axis ( with coord. y_cen ) ---
        if ( abs(axis_dir(2)) .lt. 1e-6 ) then
          call error('dust_post','','Wrong definition of the axis in sectional_loads analysis:&
               & abs(axis_dir(2)) .lt. 1e-6 . STOP')
        end if
        allocate(r_axis(3,n_sect),r_axis_bas(3,n_sect))
        do is = 1 , n_sect
          r_axis(:,is) = axis_nod + ( y_cen(is) - axis_nod(2) )/axis_dir(2) * axis_dir
        end do
        r_axis_bas = r_axis

        ! Find the id of the reference where the loads must be projected -------------
        write(filename,'(A,I4.4,A)') trim(data_basename)//'_res_',an_start,'.h5'
        call open_hdf5_file(trim(filename),floc)
        call load_refs(floc,refs_R,refs_off,refs_G,refs_f,refs_tag)
        call close_hdf5_file(floc) 

        ref_id = -333
        do it = lbound(refs_tag,1) , ubound(refs_tag,1)
          if ( stricmp(refs_tag(it),  comps(id_comp)%ref_tag) ) ref_id = it
        end do
        if ( ref_id .eq. -333 ) then 
          write(*,*)
          write(*,*) ' Available references systems: '
          do it = lbound(refs_tag,1) , ubound(refs_tag,1)
            write(*,*) ' ref_id : ' , it , ' ref_tag ' , trim(refs_tag(it))
          end do
          call error('dust_post','','Unknown ref.sys. defined for loads output.&
               & Your input in dust_post.in is '//trim(ref_tag)//'. All the&
               & available ref.sys. are listed above.')
        end if
        write(*,*) ' Available references systems: '
        do it = lbound(refs_tag,1) , ubound(refs_tag,1)
          write(*,*) ' ref_id : ' , it , ' ref_tag ' , trim(refs_tag(it))
        end do
        write(*,*) ' ref_id for force and moment projection : ' , ref_id

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

          do is = 1 , n_sect
            r_axis_bas(:,is) = matmul( refs_R(:,:,ref_id) , r_axis(:,is) ) + &
                               refs_off(:,ref_id)
          end do

!         !DEBUG
!         write(*,*) ' r_axis , r_axis_bas : '
!         do is = 1 , n_sect
!           write(*,*) r_axis(:,is) , '    ' , r_axis_bas(:,is)
!         end do

          ! 
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
            sec_loads(ires,is,1:3) = matmul( &
                 transpose( refs_R(:,:, ref_id) ) , F_bas )
            ! moment
            M_bas = matmul( &
                 transpose( refs_R(:,:, ref_id) ) , M_bas )
            sec_loads(ires,is,4) = M_bas(2)

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

!     case default

!       call error('dust_post','','''sectional_loads'' analysis not implemented yet &
!                 &for non-parametric components')
    end select

    !deallocate
    deallocate(components_names_tmp)
    deallocate(comps)

   !//////////////////    Loads     \\\\\\\\\\\\\\\\\
   case('integral_loads')

    ! load the geo components just once just once
    call open_hdf5_file(trim(data_basename)//'_geo.h5', floc)
    !TODO: here get the run id
    call load_components_postpro(comps, points, nelem, floc, & 
                                 components_names,  all_comp)

    call close_hdf5_file(floc)

    ! Prepare_geometry_postpro
    call prepare_geometry_postpro(comps)

    !n_comps_meas = countoption(sbprms,'CompName')
    !allocate( comps_meas(n_comps_meas) )
    !do ic = 1 , n_comps_meas 
    !  comps_meas(ic) = getstr(sbprms,'CompName') 
    !end do
    !!TODO:loads 
    !! from string to id.s of the ic
    !allocate( i_comps_meas(n_comps_meas) ) ; i_comps_meas = 0
    !do ic = 1 , n_comps_meas 
    !  ! loop over the ref.sys 
    !  do ic2 = 1 , size(comps)
    !    if ( trim(comps_meas(ic)) .eq. trim(comps(ic2)%comp_name) ) then
    !      i_comps_meas(ic) = ic2 
    !      exit
    !    end if
    !  end do
    !end do
    if(allocated(components_names)) deallocate(components_names) 
    allocate(components_names(size(comps)))
    do ic = 1,size(comps)
      components_names(ic) = trim(comps(ic)%comp_name)
    enddo

    ! Allocate sol_p: ...%cp is not a pointer
    !  comps(ic)%el(ie)%cp = sol_p(ip) in the loop
    allocate(sol_p(nelem)) ; sol_p = 0.0_wp

    ref_tag = getstr(sbprms,'Reference_Tag')

    nstep = (an_end-an_start)/an_step + 1
    select case(trim(out_frmt))

     case('dat')
      ! Open output .dat file
      call new_file_unit(fid_out, ierr)
      write(filename,'(A)') trim(basename)//'_'//trim(an_name)//'.dat'
      open(unit=fid_out,file=trim(filename))
      call dat_out_loads_header( fid_out , components_names , ref_tag )

     case('tecplot')
      allocate(force(3,nstep), moment(3,nstep), time(nstep))

     case default
      call error('dust_post','','Unknown format '//trim(out_frmt)//&
                 ' for loads output')
    end select
    

    ! Find the id of the reference where the loads must be projected
    write(filename,'(A,I4.4,A)') trim(data_basename)//'_res_',an_start,'.h5'
    call open_hdf5_file(trim(filename),floc)
    call load_refs(floc,refs_R,refs_off,refs_G,refs_f,refs_tag)
    call close_hdf5_file(floc) 
    !DEBUG
    ref_id = -333
    do it = lbound(refs_tag,1) , ubound(refs_tag,1)
      if ( stricmp(refs_tag(it),  ref_tag) ) ref_id = it
    end do
    if ( ref_id .eq. -333 ) then 
      write(*,*)
      write(*,*) ' Available references systems: '
      do it = lbound(refs_tag,1) , ubound(refs_tag,1)
        write(*,*) ' ref_id : ' , it , ' ref_tag ' , trim(refs_tag(it))
      end do
      call error('dust_post','','Unknown ref.sys. defined for loads output.&
           & Your input in dust_post.in is '//trim(ref_tag)//'. All the&
           & available ref.sys. are listed above.')
    end if

    ires = 0
    do it=an_start, an_end, an_step
      ires = ires+1

      ! Open the file:
      write(filename,'(A,I4.4,A)') trim(data_basename)//'_res_',it,'.h5'
      call open_hdf5_file(trim(filename),floc)

      ! Load u_inf --------------------------------
      call open_hdf5_group(floc,'Parameters',ploc)
      call read_hdf5(u_inf,'u_inf',ploc)
      call read_hdf5(P_inf,'P_inf',ploc)
      call read_hdf5(rho,'rho_inf',ploc)
      call close_hdf5_group(ploc)

      ! Load the references and move the points ---
      call load_refs(floc,refs_R,refs_off)
      ! Move the points ---------------------------
      call update_points_postpro(comps, points, refs_R, refs_off)
      ! Load the results --------------------------
      call load_res(floc, comps, vort, cp, t)
      sol_p = cp
      !
      ip = 0
      do ic = 1 , size(comps)
       do ie = 1 , size(comps(ic)%el)
        ip = ip + 1
        comps(ic)%el(ie)%cp = sol_p(ip) 
       end do
      end do

      call close_hdf5_file(floc)

      ! Initialise integral loads in the desired ref.frame
      F_ref = 0.0_wp ; M_ref = 0.0_wp 

      ! Update the overall load with the comtribution from all the components
      do ic = 1 , size(comps)

        ! Initialise integral loads in the local ref.frame
        F_bas = 0.0_wp ; M_bas = 0.0_wp 
      
        ! Loads from the ic-th component in the base ref.frame
        do ie = 1 , size(comps(ic)%el)
          F_bas1 = comps(ic)%el(ie)%dforce

          F_bas = F_bas + F_bas1

          !TODO: check the radius distance in the cross product!!
          M_bas = M_bas + cross( comps(ic)%el(ie)%cen &
                         -refs_off(:,ref_id) , F_bas1 )
        end do !ie

        ! From the base ref.sys to the chosen ref.sys (offset and rotation)
        F_ref = F_ref + matmul( &
             transpose( refs_R(:,:, ref_id) ) , F_bas )
        M_ref = M_ref + matmul( &
             transpose( refs_R(:,:, ref_id) ) , M_bas )

      end do !ic
      
      select case(trim(out_frmt))

       case ('dat')
        ! Update output file
        write(fid_out,'(F12.6)'  ,advance='no') t 
        write(fid_out,'(3F16.6)' ,advance='no') F_ref
        write(fid_out,'(3F16.6)' ,advance='no') M_ref
        write(fid_out,'(9F16.10)',advance='no') refs_R(:,:, ref_id)
        write(fid_out,'(3F16.10)',advance='no') refs_off(:, ref_id)
        write(fid_out,*) ' '

       case('tecplot')
        time(ires) = t
        force(:,ires) = F_ref
        moment(:,ires) = M_ref

      end select


    end do !it

    
    select case(trim(out_frmt))

     case('dat')
      close(fid_out)
     
     case('tecplot')
      write(filename,'(A)') trim(basename)//'_'//trim(an_name)//'.plt'
      call tec_out_loads(filename, time, force, moment)
      deallocate(time, force, moment)

    end select

    deallocate(comps,components_names)
    deallocate(sol_p)

   !//////////////////Visualizations\\\\\\\\\\\\\\\\\
   case('viz') 

    !Check which variables to analyse
    out_vort = .false.; out_vel = .false.; out_press =.false.; out_cp = .false.
    n_var = countoption(sbprms, 'Variable')
    allocate(var_names(n_var))
    do i_var = 1, n_var 
      var_names(i_var) = getstr(sbprms, 'Variable')
!     !DEBUG
!     write(*,*) ' Variable : ' , var_names(i_var)
    enddo
    out_vort = isInList('vort',var_names)
    out_vel = isInList('vel',var_names)
    out_press = isInList('press',var_names)
    out_cp = isInList('cp',var_names)

    ! load the geo components just once just once
    call open_hdf5_file(trim(data_basename)//'_geo.h5', floc)
    !TODO: here get the run id
    !DEBUG
    write(*,*) ' allocated(components_names) : ' , allocated(components_names)
    write(*,*) ' size(components_names) : ' , size(components_names)
    write(*,*) ' n_comp                 : ' , n_comp 
    if ( allocated(components_names) .and. (size(components_names) .gt. 0 ) ) then
    do i_comp = 1 , size(components_names)
      write(*,*) trim(components_names(i_comp)) 
    end do
    end if
    call load_components_postpro(comps, points, nelem, floc, & 
                                 components_names,  all_comp)
    write(*,*) ' ++++ '
    call close_hdf5_file(floc)
    out_wake = getlogical(sbprms,'Wake')

    ! Prepare_geometry_postpro
    call prepare_geometry_postpro(comps)


    !time history
    do it =an_start, an_end, an_step

      ! Open the file:
      write(filename,'(A,I4.4,A)') trim(data_basename)//'_res_',it,'.h5'
      call open_hdf5_file(trim(filename),floc)

      ! Load u_inf
      call open_hdf5_group(floc,'Parameters',ploc)
      call read_hdf5(u_inf,'u_inf',ploc)
      call read_hdf5(P_inf,'P_inf',ploc)
      call read_hdf5(rho,'rho_inf',ploc)
      call close_hdf5_group(ploc)

      ! Load the references
      call load_refs(floc,refs_R,refs_off)

      ! Move the points
      call update_points_postpro(comps, points, refs_R, refs_off)
      !expand the actuator disks
      call expand_actdisk_postpro(comps, points, points_exp, elems)

      !Load the results
      call load_res(floc, comps, vort, cp, t)

      !Prepare the variable for output
      nelem_out = size(vort)
      nprint = 0
      if(out_vort) nprint = nprint+1
      if(out_cp)   nprint = nprint+1
      allocate(print_var_names(nprint), print_vars(nelem_out, nprint))
      
      ivar = 1
      if(out_vort) then
        print_vars(:,ivar) = vort
        print_var_names(ivar) = 'Vorticity'
        ivar = ivar +1
      endif
      if(out_cp) then
        print_vars(:,ivar) = cp
        print_var_names(ivar) = 'Cp'
        ivar = ivar +1
      endif

      write(filename,'(A,I4.4)') trim(basename)//'_'//trim(an_name)//&
                                                            '_',it
      
      if (out_wake) then
        
        call load_wake_viz(floc, wpoints, welems, wvort)
        nelem_w = size(welems,2)

        nprint = 0
        if(out_vort) nprint = nprint+1
        allocate(print_var_names_w(nprint), print_vars_w(nelem_w, nprint))
        
        ivar = 1
        if(out_vort) then
          !print_vars_w(:,ivar) = reshape(wvort,(/nelem_w/))
          print_vars_w(:,ivar) = wvort
          print_var_names_w(ivar) = 'Vorticity'
          ivar = ivar +1
        endif

        !Output the results (with wake)
        select case (trim(out_frmt))
         case ('tecplot')
          filename = trim(filename)//'.plt'
          call  tec_out_viz(filename, t, &
                       points_exp, elems, print_vars, print_var_names, &
                       w_rr=wpoints, w_ee=welems, w_vars=print_vars_w, &
                       w_var_names = print_var_names_w)
         case ('vtk')
          filename = trim(filename)//'.vtu'
          call  vtk_out_viz(filename, &
                       points_exp, elems, print_vars, print_var_names, &
                       w_rr=wpoints, w_ee=welems, w_vars=print_vars_w, &
                       w_var_names = print_var_names_w)
         case default
           call error('dust_post','','Unknown format '//trim(out_frmt)//&
                      ' for visualization output')
         end select
      
        deallocate (wpoints, welems,  wvort)
        deallocate(print_var_names_w, print_vars_w)

      else
        
        !Output the results (without wake)
        select case (trim(out_frmt))
         case ('tecplot')
          filename = trim(filename)//'.plt'
          call  tec_out_viz(filename, t, &
                       points_exp, elems, print_vars, print_var_names)
         case ('vtk')
          filename = trim(filename)//'.vtu'
          call  vtk_out_viz(filename, &
                       points_exp, elems, print_vars, print_var_names)
         case default
           call error('dust_post','','Unknown format '//trim(out_frmt)//&
                      ' for visualization output')
         end select

      endif

      call close_hdf5_file(floc)
      deallocate(refs_R, refs_off, vort, cp)
      deallocate(print_var_names, print_vars)

    enddo !Time history

    deallocate(comps, points)


   !//////////////////Domain probes \\\\\\\\\\\\\\\\\
   case('probes')

    ! Read probe coordinates: point_list or from_file
    in_type =  getstr(sbprms,'InputType')
    select case ( trim(in_type) )
     case('point_list')
      n_probes = countoption(sbprms,'Point')
      allocate(rr_probes(3,n_probes))
      do ip = 1 , n_probes
        rr_probes(:,ip) = getrealarray(sbprms,'Point',3)
      end do
     case('from_file')
      filename_in = getstr(sbprms,'File')   ! N probes and then their coordinates
      fid_out = 21
      open(unit=fid_out,file=trim(filename_in))
      read(fid_out,*) n_probes
      allocate(rr_probes(3,n_probes))
      do ip = 1 , n_probes
        read(fid_out,*) rr_probes(:,ip)
      end do
      close(fid_out)
     case default
      write(str_a,*) ia 
      call error('dust_post','','Unknown InputType: '//trim(in_type)//&
                 ' for analysis n.'//trim(str_a)//'.'//nl//&
                  'It must either "point_list" or "from_file".')
    end select

    ! Read variables to save : velocity | pressure | vorticity
    probe_vel = .false. ; probe_p = .false. ; probe_vort = .false.
    n_vars = countoption(sbprms,'Variable')

    if ( n_vars .eq. 0 ) then ! default: velocity | pressure | vorticity
      probe_vel = .true. ; probe_p = .true. ; probe_vort = .true.
    else
     do it = 1 , n_vars
      var_name = getstr(sbprms,'Variable')
      select case(trim(var_name))
       case ( 'velocity' ) ; probe_vel = .true.
       case ( 'pressure' ) ; probe_p   = .true.
       case ( 'vorticity') ; probe_vort= .true.
       case ( 'all') ; probe_vel = .true. ; probe_p   = .true. ; probe_vort= .true.
       case default
       call error('dust_post','','Unknown Variable: '//trim(var_name)//&
                  ' for analysis n.'//trim(str_a)//'.'//nl//&
                   'Choose "velocity", "pressure", "vorticity".')
      end select
     end do
    end if
     
    nprint = 0
    if(probe_vel) nprint = nprint+3
    if(probe_p)   nprint = nprint+1
    if(probe_vort) nprint = nprint+3
    allocate(probe_var_names(nprint))
    allocate(probe_loc_names(n_probes))
    probe_var_names = ''
    ivar = 1
    if(probe_vel) then
      probe_var_names(ivar) = 'ux'
      probe_var_names(ivar + 1) = 'uy'
      probe_var_names(ivar + 2) = 'uz'
      ivar = ivar + 3
    endif
    if(probe_p) then
      probe_var_names(ivar) = 'p'
      ivar = ivar +1
    endif
    if(probe_vort) then
      probe_var_names(ivar) = 'omx'
      probe_var_names(ivar + 1) = 'omy'
      probe_var_names(ivar + 2) = 'omz'
      ivar = ivar + 3
    endif
    vars_str = ''
    do ivar = 1,size(probe_var_names)
      vars_str = trim(vars_str)//'  '//trim(probe_var_names(ivar))
    enddo

    do ip = 1,n_probes
      write(probe_loc_names(ip),'(A,F10.5,A,F10.5,A,F10.5)') &
           'x=',rr_probes(1,ip), &
           'y=',rr_probes(2,ip), &
           'z=',rr_probes(3,ip)

    enddo

    !Allocate where the solution will be stored
    nstep = (an_end-an_start)/an_step + 1
    allocate(probe_vars(nprint, nstep, n_probes))
    allocate(time(nstep))


    ! load the geo components just once
    call open_hdf5_file(trim(data_basename)//'_geo.h5', floc)
    !TODO: here get the run id
    call load_components_postpro(comps, points, nelem,  floc, & 
                                 components_names,  all_comp)

    call close_hdf5_file(floc)

    ! Prepare_geometry_postpro
    call prepare_geometry_postpro(comps)

    ! Allocate and point to sol
    allocate(sol(nelem)) ; sol = 0.0_wp
    ip = 0
    do ic = 1 , size(comps)
     do ie = 1 , size(comps(ic)%el)
      ip = ip + 1
      comps(ic)%el(ie)%idou => sol(ip) 
     end do
    end do

    !time history
    ires = 0
    do it =an_start, an_end, an_step
      ires = ires+1

     ! Open the result file ----------------------
     write(filename,'(A,I4.4,A)') trim(data_basename)//'_res_',it,'.h5'
     call open_hdf5_file(trim(filename),floc)

     ! Load u_inf --------------------------------
     call open_hdf5_group(floc,'Parameters',ploc)
     call read_hdf5(u_inf,'u_inf',ploc)
     call read_hdf5(P_inf,'P_inf',ploc)
     call read_hdf5(rho,'rho_inf',ploc)
     call close_hdf5_group(ploc)

     ! Load the references and move the points ---
     call load_refs(floc,refs_R,refs_off)
     call update_points_postpro(comps, points, refs_R, refs_off)
     ! Load the results --------------------------
     call load_res(floc, comps, vort, cp, t)
     !sol = vort

     ! Load the wake -----------------------------
     call load_wake_pan(floc, wpoints_pan, wstart, wvort_pan)
     call load_wake_ring(floc, wpoints_rin, wconn, wvort_rin)
   
     call close_hdf5_file(floc)
     
     call prepare_wake_postpro( wpoints_pan, wpoints_rin, wstart, wconn, &
                 wvort_pan,  wvort_rin, wake_pan, wake_rin, wake_elems )

     time(ires) = t

     ! Compute velocity --------------------------
     do ip = 1 , n_probes ! probes

      vel_probe = 0.0_wp ; pres_probe = 0.0_wp ; vort_probe = 0.0_wp
      ivar = 1
      if ( probe_vel .or. probe_p ) then 

        ! compute velocity
        do ic = 1,size(comps)
         do ie = 1 , size( comps(ic)%el )

          call comps(ic)%el(ie)%compute_vel( rr_probes(:,ip) , u_inf , v )
          vel_probe = vel_probe + v/(4*pi) 
         
         end do
        end do

        do ie = 1, size(wake_elems)
          call wake_elems(ie)%p%compute_vel( rr_probes(:,ip) , u_inf , v )
          vel_probe = vel_probe + v/(4*pi) 
        enddo

        ! + u_inf
        vel_probe = vel_probe + u_inf
      end if

      if(probe_vel) then
        probe_vars(ivar:ivar+2, ires, ip) = vel_probe
        ivar = ivar+3
      endif

      ! compute pressure
      if ( probe_p ) then
        ! Bernoulli equation
        ! rho * dphi/dt + P + 0.5*rho*V^2 = P_infty + 0.5*rho*V_infty^2
        !TODO: add:
        ! - add the unsteady term: -rho*dphi/dt
        pres_probe = P_inf + 0.5_wp*rho*norm2(u_inf)**2 - 0.5_wp*rho*norm2(vel_probe)**2

        probe_vars(ivar, ires, ip) = pres_probe
        ivar = ivar+1
      end if
      
      ! compute vorticity
      if ( probe_vort ) then
        vort_probe = 2.0_wp
        !write(fid_out,'(3F12.6)',advance='no') vort_probe
        probe_vars(ivar:ivar+2, ires, ip) = vort_probe
        ivar = ivar+3
      end if

     end do  ! probes


    end do ! Time history

    !Output the results in the correct format
    select case (trim(out_frmt))

     case ('dat')
      call new_file_unit(fid_out, ierr)
      write(filename,'(A)') trim(basename)//'_'//trim(an_name)//'.dat'
      open(unit=fid_out,file=trim(filename)) 
      call dat_out_probes_header( fid_out , rr_probes , vars_str )

      do ires = 1, size(time)
        write(fid_out,'(F12.6)',advance='no') time(ires)
        do ip = 1, n_probes
          do ivar = 1, size(probe_vars,1)
            write(fid_out,'(3F12.6)',advance='no') probe_vars(ivar,ires,ip)
            write(fid_out,'(A)',advance='no') ' '
          enddo
        enddo
        write(fid_out,*) ' '
      enddo
      close(fid_out)

     case ('tecplot')
      write(filename,'(A)') trim(basename)//'_'//trim(an_name)//'.plt'
      call tec_out_probes(filename, time, probe_vars, probe_var_names, probe_loc_names)
    case default
      call error('dust_post','','Unknown format '//trim(out_frmt)//&
                 ' for probe output')
    end select

    deallocate(comps)
    deallocate(rr_probes,sol)
    deallocate(probe_vars, probe_var_names, time)


   !//////////////////Flow Field \\\\\\\\\\\\\\\\\
   case('flow_field')

    allocate(var_names(3)) ; var_names = ' '
    allocate(vars_n   (3)) ; vars_n = 0

    ! Read variables to save : velocity | pressure | vorticity
    probe_vel = .false. ; probe_p = .false. ; probe_vort = .false.
    n_vars = countoption(sbprms,'Variable')

    if ( n_vars .eq. 0 ) then ! default: velocity | pressure | vorticity
      probe_vel = .true. ; probe_p = .true. ; probe_vort = .true.
    else
     do it = 1 , n_vars
      var_name = getstr(sbprms,'Variable')
      select case(trim(var_name))
       case ( 'velocity' ) ; probe_vel = .true.
       case ( 'pressure' ) ; probe_p   = .true.
       case ( 'vorticity') ; probe_vort= .true.
       case ( 'all') ; probe_vel = .true. ; probe_p   = .true. ; probe_vort= .true.
       case default
       call error('dust_post','','Unknown Variable: '//trim(var_name)//&
                  ' for analysis n.'//trim(str_a)//'.'//nl//&
                   'Choose "velocity", "pressure", "vorticity".')
      end select
     end do
    end if

    ! load the geo components just once just once
    call open_hdf5_file(trim(data_basename)//'_geo.h5', floc)
    !TODO: here get the run id    !todo????
    call load_components_postpro(comps, points, nelem, floc, & 
                                 components_names,  all_comp)

    call close_hdf5_file(floc)

    ! Prepare_geometry_postpro
    call prepare_geometry_postpro(comps)

    ! Allocate and point to sol
    allocate(sol(nelem)) ; sol = 0.0_wp
    ip = 0
    do ic = 1 , size(comps)
     do ie = 1 , size(comps(ic)%el)
      ip = ip + 1
      comps(ic)%el(ie)%idou => sol(ip) 
     end do
    end do

    ! Read box dimensions
    nxyz   = getintarray( sbprms,'Nxyz'  ,3)
    minxyz = getrealarray(sbprms,'Minxyz',3)
    maxxyz = getrealarray(sbprms,'Maxxyz',3)

!   allocate( xbox(nxyz(1)), ybox(nxyz(2)) , zbox(nxyz(3)) ) 
    if ( nxyz(1) .gt. 1 ) then
      allocate( xbox(nxyz(1)) )
      dxbox = ( maxxyz(1) - minxyz(1) ) / dble( nxyz(1) - 1 ) 
      xbox = (/ ( minxyz(1) + dble(i1-1) * dxbox , i1 = 1 , nxyz(1) )/) 
    else
      allocate( xbox(1) )
      xbox(1) = minxyz(1)
    end if
    if ( nxyz(2) .gt. 1 ) then
      allocate( ybox(nxyz(2)) )
      dybox = ( maxxyz(2) - minxyz(2) ) / dble( nxyz(2) - 1 ) 
      ybox = (/ ( minxyz(2) + dble(i1-1) * dybox , i1 = 1 , nxyz(2) )/) 
    else
      allocate( ybox(1) )
      ybox(1) = minxyz(2)
    end if
    if ( nxyz(3) .gt. 1 ) then
      allocate( zbox(nxyz(3)) )
      dzbox = ( maxxyz(3) - minxyz(3) ) / dble( nxyz(3) - 1 ) 
      zbox = (/ ( minxyz(3) + dble(i1-1) * dzbox , i1 = 1 , nxyz(3) )/) 
    else
      allocate( zbox(1) )
      zbox(1) = minxyz(3)
    end if

    i_var = 0
    i_var_v = 0 ; i_var_p = 0 ; i_var_w = 0
    if ( probe_vel ) then
      allocate(box_vel (product(nxyz),3))
      i_var = i_var + 1
      var_names(i_var) = 'velocity'
      vars_n(i_var) = 3
      i_var_v = 3
    end if
    if ( probe_p   ) then
      allocate(box_p   (product(nxyz)  ))
      i_var = i_var + 1
      var_names(i_var) = 'pressure'
      vars_n(i_var) = 1
      i_var_p = i_var_v + 1
    end if 
    if ( probe_vort) then
      allocate(box_vort(product(nxyz),3))
      i_var = i_var + 1
      var_names(i_var) = 'vorticity'
      vars_n(i_var) = 3
      i_var_w = i_var_p + 3
    end if
    allocate(vars(sum(vars_n),product(nxyz))) ; vars = 1.0_wp 

    do it = an_start, an_end, an_step ! Time history

     ! Open the result file ----------------------
     write(filename,'(A,I4.4,A)') trim(data_basename)//'_res_',it,'.h5'
     call open_hdf5_file(trim(filename),floc)

     ! Load u_inf --------------------------------
     call open_hdf5_group(floc,'Parameters',ploc)
     call read_hdf5(u_inf,'u_inf',ploc)
     call read_hdf5(P_inf,'P_inf',ploc)
     call read_hdf5(rho,'rho_inf',ploc)
     call close_hdf5_group(ploc)

     ! Load the references and move the points ---
     call load_refs(floc,refs_R,refs_off,refs_G,refs_f)

!    !DEBUG
!    write(*,*) refs_G
!    write(*,*)
!    write(*,*) refs_f

!    !DEBUG
!    write(*,*) ' before update_points_postpro '
!    write(*,*) ' refs_G : ' , refs_G
     call update_points_postpro(comps, points, refs_R, refs_off, refs_G, refs_f)
!    write(*,*) '  after update_points_postpro '


     ! Load the results --------------------------
     call load_res(floc, comps, vort, cp, t)
     !sol = vort

     ! Load the wake -----------------------------
     call load_wake_pan(floc, wpoints_pan, wstart, wvort_pan)
     call load_wake_ring(floc, wpoints_rin, wconn, wvort_rin)
   
     call close_hdf5_file(floc)
     
     call prepare_wake_postpro( wpoints_pan, wpoints_rin, wstart, wconn, &
                 wvort_pan,  wvort_rin, wake_pan, wake_rin, wake_elems )

     ! Compute velocity --------------------------
     ip = 0
     ! Loop over the nodes of the box
     do iz = 1 , size(zbox)
      do iy = 1 , size(ybox)
       do ix = 1 , size(xbox)
        ip = ip + 1

        if ( probe_vel .or. probe_p ) then 
          ! compute velocity
          vel_probe = 0.0_wp ; pres_probe = 0.0_wp ; vort_probe = 0.0_wp
          do ic = 1,size(comps)
           do ie = 1 , size( comps(ic)%el )

!           call comps(ic)%el(ie)%compute_vel( rr_probes(:,ip) , u_inf , v )
            call comps(ic)%el(ie)%compute_vel( (/ xbox(ix) , ybox(iy) , zbox(iz) /) , & 
                                                u_inf , v )
            vel_probe = vel_probe + v/(4*pi) 
           
           end do
          end do
          ! wake contribution
          !do ic = 1 , size(wake%wake_panels,1)
          ! do ie = 1 , size(wake%wake_panels,2)
          !       
          !  call wake%wake_panels(ic,ie)%compute_vel( &
          !       (/ xbox(ix) , ybox(iy) , zbox(iz) /) , u_inf , v )
          !  vel_probe = vel_probe + v/(4*pi) 
          ! 
          ! end do
          !end do
  
          do ie = 1, size(wake_elems)
            call wake_elems(ie)%p%compute_vel( &
                     (/ xbox(ix) , ybox(iy) , zbox(iz) /) , u_inf , v )
            vel_probe = vel_probe + v/(4*pi) 
          enddo
         
          ! + u_inf
          vel_probe = vel_probe + u_inf
          
        end if
        if ( probe_vel ) then
!         vars(1:3,ix+(iy-1)*size(xbox)+(iz-1)*size(xbox)*size(ybox)) = vel_probe
          vars(1:3,ip) = vel_probe
        end if

        if ( probe_p ) then
          ! Bernoulli equation
          ! rho * dphi/dt + P + 0.5*rho*V^2 = P_infty + 0.5*rho*V_infty^2
          !TODO: add:
          ! - add the unsteady term: -rho*dphi/dt
          pres_probe = P_inf + 0.5_wp*rho*norm2(u_inf)**2 - 0.5_wp*rho*norm2(vel_probe)**2
          vars(i_var_v+1,ip) = pres_probe 
   !       
        end if

        if ( probe_vort ) then

        end if

       end do
      end do
     end do


    select case (trim(out_frmt))

    case ('vtk')
     write(filename,'(A,I4.4,A)') trim(basename)//'_'//trim(an_name)//&
                                                            '_',it,'.vtr'
     call vtr_write ( filename , xbox , ybox , zbox , &
                      vars_n(1:i_var) , var_names(1:i_var) , vars ) 

    case('tecplot')
     write(filename,'(A,I4.4,A)') trim(basename)//'_'//trim(an_name)//&
                                                            '_',it,'.plt'

     i_var = 0
     deallocate(var_names)
     allocate(var_names(7))
     if(probe_vel) then
       var_names(i_var + 1) = 'ux'
       var_names(i_var + 2) = 'uy'
       var_names(i_var + 3) = 'uz'
       i_var = i_var + 3
     endif
     if(probe_p) then
       var_names(i_var + 1) = 'p'
       i_var = i_var + 1
     endif
     if(probe_vort) then
       var_names(i_var + 1) = 'omx'
       var_names(i_var + 2) = 'omy'
       var_names(i_var + 3) = 'omz'
       i_var = i_var + 3
     endif

     call tec_out_box(filename, t, xbox, ybox, zbox, vars, var_names(1:i_var))

    case default
      call error('dust_post','','Unknown format '//trim(out_frmt)//&
                 ' for flowfield output')
    end select

!    !CHECK
!    close(fid_out)

    end do ! Time history

    if (allocated(box_vel )) deallocate(box_vel )
    if (allocated(box_p   )) deallocate(box_p   )
    if (allocated(box_vort)) deallocate(box_vort)
    deallocate(xbox,ybox,zbox)
    deallocate(comps,sol)

    deallocate(var_names,vars_n,vars)

   case default
    call error('dust_post','','Unknown type of analysis: '//trim(an_type))

  end select


  
  if(allocated(var_names)) deallocate(var_names)
  if(allocated(components_names)) deallocate(components_names)

  sbprms => null()
enddo !analysis

call destroy_hdf5()
call printout(nl//'<<<<<< DUST POSTPROCESSOR end       <<<<<<'//nl)

!----------------------------------------------------------------------
!----------------------------------------------------------------------
contains
!----------------------------------------------------------------------
!----------------------------------------------------------------------

subroutine load_refs(floc, refs_R, refs_off, refs_G, refs_f, refs_tag)
 integer(h5loc), intent(in) :: floc 
 real(wp), allocatable, intent(out) :: refs_R(:,:,:)
 real(wp), allocatable, intent(out) :: refs_off(:,:)
 real(wp), allocatable, intent(out) , optional :: refs_G(:,:,:)
 real(wp), allocatable, intent(out) , optional :: refs_f(:,:)
 character(len=max_char_len) , allocatable , intent(out) , optional :: refs_tag(:)

 integer(h5loc) :: gloc1, gloc2
 integer :: nrefs, iref
 character(len=max_char_len) :: rname

  call open_hdf5_group(floc,'References',gloc1)
  call read_hdf5(nrefs,'NReferences',gloc1)

  allocate(refs_R(3,3,0:nrefs-1), refs_off(3,0:nrefs-1))
  if (present(refs_G)  ) allocate(refs_G(3,3,0:nrefs-1))
  if (present(refs_f)  ) allocate(refs_f(3,0:nrefs-1))
  if (present(refs_tag)) allocate(refs_tag(0:nrefs-1))
  do iref = 0,nrefs-1
    write(rname,'(A,I3.3)') 'Ref',iref
    call open_hdf5_group(gloc1,trim(rname),gloc2)
   
    call read_hdf5(refs_R(:,:,iref),'R',gloc2)
    call read_hdf5(refs_off(:,iref),'Offset',gloc2)
    if (present(refs_tag)) call read_hdf5(refs_tag(  iref),'Tag',gloc2)
    if (present(refs_G  )) call read_hdf5(refs_G(:,:,iref),'RotVel',gloc2)
    if (present(refs_f  )) call read_hdf5(refs_f(:,iref),'Vel',gloc2)
    call read_hdf5(refs_off(:,iref),'Offset',gloc2)

    call close_hdf5_group(gloc2)
  enddo

  call close_hdf5_group(gloc1)

end subroutine load_refs

!----------------------------------------------------------------------

subroutine load_res(floc, comps, vort, cp, t)
 integer(h5loc), intent(in) :: floc 
 type(t_geo_component), intent(inout) :: comps(:)
 real(wp), allocatable, intent(out) :: vort(:)
 real(wp), allocatable, intent(out) :: cp(:)
 real(wp), intent(out) :: t

 integer :: ncomps, icomp, ie
 integer :: nelems, offset, nelems_comp
 integer(h5loc) :: gloc1, gloc2, gloc3
 character(len=max_char_len) :: cname
 real(wp), allocatable :: vort_read(:), cp_read(:)
 real(wp), allocatable :: pres_read(:), dforce_read(:,:)

  ncomps = size(comps)
  nelems = 0
  do icomp = 1, ncomps
    select type(el=>comps(icomp)%el)
     class default
      nelems = nelems + comps(icomp)%nelems 
     type is(t_actdisk)
      nelems = nelems + size(comps(icomp)%loc_points,2)
    end select
  enddo
  
  call read_hdf5(t,'time',floc)
  allocate(vort(nelems), cp(nelems))
  call open_hdf5_group(floc,'Components',gloc1)

  offset = 0
  do icomp = 1, ncomps
    
    nelems_comp = comps(icomp)%nelems
    write(cname,'(A,I3.3)') 'Comp',comps(icomp)%comp_id
    call open_hdf5_group(gloc1,trim(cname),gloc2)
    call open_hdf5_group(gloc2,'Solution',gloc3)

    !call read_hdf5(vort(offset+1:offset+nelems_comp),'Vort',gloc3)
    !call read_hdf5(cp(offset+1:offset+nelems_comp),'Cp',gloc3)
    call read_hdf5_al(vort_read,'Vort',gloc3)
    call read_hdf5_al(cp_read,'Cp',gloc3)
    call read_hdf5_al(pres_read,'Pres',gloc3)
    call read_hdf5_al(dforce_read,'dF',gloc3)

!   TODO: check if it is general enough *******
!   TODO: check if something is broken after changing intent(in to inout) for comps
    do ie = 1 , nelems_comp
      comps(icomp)%el(ie)%pres = pres_read(ie) 
      if( .not. allocated(comps(icomp)%el(ie)%dforce) ) then
        allocate(comps(icomp)%el(ie)%dforce(3))
      end if
      comps(icomp)%el(ie)%dforce = dforce_read(ie,:) 
    end do

    call close_hdf5_group(gloc3)
    call close_hdf5_group(gloc2)

    select type(el =>comps(icomp)%el)
     class default
      vort(offset+1:offset+nelems_comp) = vort_read
      cp(offset+1:offset+nelems_comp) = cp_read
      offset = offset + nelems_comp
      do ie = 1,nelems_comp
        if(associated(comps(icomp)%el(ie)%idou)) &
                        comps(icomp)%el(ie)%idou = vort_read(ie)
      enddo
     type is(t_actdisk)
      do ie = 1,nelems_comp
        vort(offset+1:offset+el(ie)%n_ver) = vort_read(ie)
        cp(offset+1:offset+el(ie)%n_ver) = cp_read(ie)
        offset = offset + el(ie)%n_ver
        if(associated(comps(icomp)%el(ie)%idou)) then
                        comps(icomp)%el(ie)%idou = vort_read(ie)
        endif
      enddo
    end select

    deallocate(vort_read, cp_read)

  enddo

  call close_hdf5_group(gloc1)

end subroutine load_res

!----------------------------------------------------------------------

subroutine load_wake_pan(floc, wpoints, wstart, wvort)
 integer(h5loc), intent(in) :: floc 
 real(wp), allocatable, intent(out) :: wpoints(:,:,:)
 integer, allocatable, intent(out) :: wstart(:,:)
 real(wp), allocatable, intent(out) :: wvort(:,:)

 integer(h5loc) :: gloc
  
  call open_hdf5_group(floc,'PanelWake',gloc)
  
  call read_hdf5_al(wpoints,'WakePoints',gloc)
  call read_hdf5_al(wstart,'StartPoints',gloc)
  call read_hdf5_al(wvort,'WakeVort',gloc)

  call close_hdf5_group(gloc)

end subroutine load_wake_pan

!----------------------------------------------------------------------

subroutine load_wake_ring(floc, wpoints, wconn, wvort)
 integer(h5loc), intent(in) :: floc 
 real(wp), allocatable, intent(out) :: wpoints(:,:,:)
 integer, allocatable, intent(out) :: wconn(:)
 !real(wp), allocatable, intent(out) :: wcen(:,:,:)
 real(wp), allocatable, intent(out) :: wvort(:,:)

 integer(h5loc) :: gloc
  
  call open_hdf5_group(floc,'RingWake',gloc)
  
  call read_hdf5_al(wpoints,'WakePoints',gloc)
  call read_hdf5_al(wconn,'Conn_pe',gloc)
!  call read_hdf5_al(wcen,'WakeCenters',gloc)
  call read_hdf5_al(wvort,'WakeVort',gloc)

  call close_hdf5_group(gloc)

end subroutine load_wake_ring

!----------------------------------------------------------------------

subroutine load_wake_viz(floc, wpoints, welems, wvort)
 integer(h5loc), intent(in) :: floc 
 real(wp), allocatable, intent(out) :: wpoints(:,:)
 integer, allocatable, intent(out)  :: welems(:,:)
 real(wp), allocatable, intent(out) :: wvort(:)

 integer(h5loc) :: gloc
 logical :: got_dset 
 real(wp), allocatable :: wpoints_read(:,:,:)
 real(wp), allocatable :: wpoints_pan(:,:), wpoints_rin(:,:)
 integer, allocatable  :: wstart(:,:), wconn(:)
 real(wp), allocatable :: wcen(:,:,:)
 real(wp), allocatable :: wvort_read(:,:)
 real(wp), allocatable :: wvort_pan(:), wvort_rin(:)
 integer, allocatable  :: welems_pan(:,:), welems_rin(:,:)
 integer :: nstripes, npoints_row, nrows, ndisks, nelem_w
 integer :: iew, ir, is, ip
 integer :: first_elem, act_disk, next_elem


 !get the panel wake
 got_dset = check_dset_hdf5('PanelWake',floc)
 if(got_dset) then

  call open_hdf5_group(floc,'PanelWake',gloc)
  call read_hdf5_al(wpoints_read,'WakePoints',gloc)
  call read_hdf5_al(wstart,'StartPoints',gloc)
  call read_hdf5_al(wvort_read,'WakeVort',gloc)


  nstripes = size(wvort_read,1); nrows = size(wvort_read,2);
  npoints_row = size(wpoints_read,2)
  nelem_w = nstripes*nrows

  allocate(wpoints_pan(3,npoints_row*(nrows+1)))
  allocate(welems_pan(4,nelem_w))
  allocate(wvort_pan(nelem_w))

  wpoints_pan = reshape(wpoints_read, (/3,npoints_row*(nrows+1)/))
  wvort_pan = reshape(wvort_read, (/nelem_w/))
  iew = 0
  do ir = 1,nrows
    do is = 1,nstripes
      iew = iew+1
      welems_pan(:,iew) = (/wstart(1,is)+npoints_row*(ir-1), &
                            wstart(2,is)+npoints_row*(ir-1), &
                            wstart(2,is)+npoints_row*(ir), &
                            wstart(1,is)+npoints_row*(ir)/)
    enddo
  enddo
 
  deallocate(wpoints_read, wstart, wvort_read)
  call close_hdf5_group(gloc)
 else
   !panel wake not present, allocate stuff at zero size
   allocate(wvort_pan(0), wpoints_pan(3,0), welems_pan(4,0))
 endif

 !get the ring wake
 got_dset = check_dset_hdf5('RingWake',floc)
 if(got_dset) then

  call open_hdf5_group(floc,'RingWake',gloc)
  call read_hdf5_al(wpoints_read,'WakePoints',gloc)
  call read_hdf5_al(wconn,'Conn_pe',gloc)
  call read_hdf5_al(wcen,'WakeCenters',gloc)
  call read_hdf5_al(wvort_read,'WakeVort',gloc)
  
  ndisks = size(wvort_read,1); nrows = size(wvort_read,2)
  npoints_row = size(wpoints_read,2)

  nelem_w = ndisks*nrows

  allocate(wpoints_rin(3,npoints_row*nrows+nelem_w))
  allocate(welems_rin(4,npoints_row*nrows))
  allocate(wvort_rin(npoints_row*nrows))

  wpoints_rin(:,1:npoints_row*nrows) = reshape(wpoints_read, &
                                          (/3,npoints_row*nrows/))
  wpoints_rin(:,npoints_row*nrows+1:size(wpoints_rin,2)) = &
         reshape(wcen, (/3,nelem_w/))

  iew = 0; act_disk = 0
  do ir = 1,nrows
    do ip = 1, npoints_row
      iew = iew+1

      if(ip .eq. 1) then
        first_elem = iew
      else
        if(wconn(ip-1) .ne. wconn(ip)) first_elem = iew
      endif

      if(ip.lt.npoints_row) then
        if(wconn(ip+1) .ne. wconn(ip)) then
          next_elem = first_elem
        else
          next_elem = iew + 1
        endif
      else
        next_elem = first_elem
      endif

      welems_rin(1,iew) = iew
      welems_rin(2,iew) = next_elem
      !welems_rin(3,iew) = npoints_row*nrows+(ir-1)*nelem_w+wconn(ip)
      welems_rin(3,iew) = npoints_row*nrows+(ir-1)*ndisks+wconn(ip)
      welems_rin(4,iew) = 0

      wvort_rin(iew) = wvort_read(wconn(ip),ir)
    enddo
  enddo

  call close_hdf5_group(gloc)
  deallocate(wpoints_read, wconn, wcen, wvort_read) 
 else
   allocate(wvort_rin(0), wpoints_rin(3,0), welems_rin(4,0))
 endif

 !Stitch together the two wakes
 allocate(wpoints(3,size(wpoints_pan,2)+size(wpoints_rin,2)))
 wpoints(:,1:size(wpoints_pan,2)) = wpoints_pan
 wpoints(:,size(wpoints_pan,2)+1:size(wpoints,2)) = wpoints_rin
 deallocate(wpoints_pan, wpoints_rin)

 allocate(welems(4,size(welems_pan,2)+size(welems_rin,2)))
 welems(:,1:size(welems_pan,2)) = welems_pan
 welems(1:3,size(welems_pan,2)+1:size(welems,2)) = welems_rin(1:3,:) + &
                                                   size(wpoints_pan,2)
 welems(4,size(welems_pan,2)+1:size(welems,2)) = 0
 deallocate(welems_pan, welems_rin)


 allocate(wvort(size(wvort_pan)+size(wvort_rin)))
 wvort(1:size(wvort_pan)) = wvort_pan
 wvort(size(wvort_pan)+1:size(wvort)) = wvort_rin
 deallocate(wvort_pan, wvort_rin)

end subroutine load_wake_viz

!----------------------------------------------------------------------

end program dust_post
