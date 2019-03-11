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

!> Module to treat the geometry of the solid bodies

module mod_geometry

use mod_param, only: &
  wp, max_char_len, nl, prev_tri, next_tri, prev_qua, next_qua

use mod_sim_param, only: &
  t_sim_param, sim_param

use mod_parse, only: &
  t_parse, getstr, getint, getreal, getrealarray, getlogical, countoption &
  , finalizeparameters

use mod_handling, only: &
  error, warning, info, printout, dust_time, t_realtime, check_preproc, &
  check_file_exists

use mod_basic_io, only: &
  read_mesh_basic, write_basic

use mod_cgns_io, only: &
  read_mesh_cgns

use mod_parametric_io, only: &
  read_mesh_parametric

!use mod_aero_elements, only: &
!  c_elem, t_elem_p

use mod_aeroel, only: &
  c_elem, c_pot_elem, c_vort_elem, c_impl_elem, c_expl_elem, &
  t_elem_p, t_pot_elem_p, t_vort_elem_p, t_impl_elem_p, t_expl_elem_p


use mod_surfpan, only: &
  t_surfpan

use mod_vortlatt, only: &
  t_vortlatt

use mod_liftlin, only: &
  t_liftlin

use mod_actuatordisk, only: &
  t_actdisk

use mod_c81, only: &
  t_aero_tab , read_c81_table , interp_aero_coeff

use mod_math, only: &
  cross

use mod_reference, only: &
  t_ref, build_references, update_all_references, destroy_references   ! , &
! update_relative_initial_conditions

use mod_hdf5_io, only: &
   h5loc, &
   new_hdf5_file, &
   open_hdf5_file, &
   close_hdf5_file, &
   new_hdf5_group, &
   new_hdf5_group, &
   open_hdf5_group, &
   close_hdf5_group, &
   write_hdf5, &
   write_hdf5_attr, &
   read_hdf5, &
   read_hdf5_al, &
   check_dset_hdf5

!----------------------------------------------------------------------

implicit none

public :: t_geo, t_geo_component, t_tedge, &
          create_geometry, update_geometry, destroy_geometry, &
          calc_geo_data_pan, calc_geo_data_ll, calc_geo_data_ad,&
          calc_node_vel , calc_geo_vel, destroy_elements

private

!----------------------------------------------------------------------

!> Geometry component
!!
!! A geometry component is a single geometrical component of the model,
!! such as a wing, an aileron, a rotor blade etc.
!!
!! It is used to load the mesh and link it to the proper refernce frame.
!!
!! It also contains the vector of all the aerodynamic elements generated
!! by the grid
!!
!! It also contains some temporary arrays employed during geometry generation
type :: t_geo_component

 !> Element type
 character(len=max_char_len) :: comp_el_type

 !> Name of the component
 character(len=max_char_len) :: comp_name

 !> Type of input: basic, cgns, parametric
 character(len=max_char_len) :: comp_input

 !> Id of the component (warning: not always defined)
 integer :: comp_id

 !> Number of elements of the component
 integer :: nelems

 !> Elements of the group
 !! The main memory allocation happens here, they will be pointed from
 !! somewhere else
 class(c_pot_elem), allocatable :: el(:)

 !> Global indexes of the points in the component
 integer, allocatable :: i_points(:)

 !> Reference frame tag
 character(len=max_char_len) :: ref_tag

 ! Reference frame id
 integer :: ref_id

 !> Points in local reference frame
 real(wp), allocatable :: loc_points(:,:)

 !> Number of surface panels in the component
 integer :: nSurfPan
 !> Number of vortex rings in the component
 integer :: nVortRin
 !> Number of lifting line elements in the component
 integer :: nLiftLin
 !> Number of lifting line ele
 integer :: nActDisk

 !> Is the component moving?
 logical :: moving

 !> Arrays for lifting line components
 character(len=max_char_len), allocatable :: airfoil_list(:)
 !> Number of ll elements in each region. Size(nelem_span_list) = nRegions
 integer, allocatable :: nelem_span_list(:)
 !> Normalised coordinate of sections in each wing region.
 !! Size(...) = nRegions+1
 real(wp),allocatable :: normalised_coord_e(:,:)
 !> Id of the airfoil elements (index in airfoil_list char array)
 integer ,allocatable :: i_airfoil_e(:,:)

 !> Dimensions of parametric elements only
 integer :: parametric_nelems_span , parametric_nelems_chor 

end type  t_geo_component

!-----------------------------------


!> Geometry type
!!
!! The type contains the whole geometry, with all the components and reference
!! frames and all the points of the mesh
type :: t_geo

 !> Total number of implicit elements
 integer :: nelem_impl

 !> Number of lifting line elements
 !integer :: nll

 !> Number of actuator disk elements
 !integer :: nad

 !> Number of explicit elements
 integer :: nelem_expl

 !> Number of surface panel elements
 integer :: nSurfPan
 !> Number of vortex ring elements
 integer :: nVortLatt
 !> Number of lifting line elements
 integer :: nLiftLin
 !> Number of Actuator disk elements
 integer :: nActDisk

 ! Only for implicit elements (for slicing)
 !> Global id of surface panel elements
 integer, allocatable :: idSurfPan(:)
 !> Global id of vortex ring elements
 integer, allocatable :: idVortLatt(:)
!!> Global id of lifting line elements
!integer, allocatable :: idLiftLin(:)
!!> Global id of actuator disk elements
!integer, allocatable :: idActDisk(:)
 !> Surfpan id of global surface panel elements (G2L: global to local)
 integer, allocatable :: idSurfPanG2L(:)

 !> Number of statical implicit elements
 integer :: nstatic_impl
 !> Number of moving implicit elements
 integer :: nmoving_impl

 !> Number of statical surfpan elements (for slicing)
 integer :: nstatic_SurfPan

 !> Number of static or moving lifting lines
 integer :: nstatic_expl, nmoving_expl

 !> All the components of the geometry
 type(t_geo_component), allocatable :: components(:)

 !> Points (element vertexes)
 real(wp), allocatable :: points(:,:)

 !> All the reference frames of the geometry
 type(t_ref), allocatable :: refs(:)

end type t_geo

!-----------------------------------

!> Trailing edge type
!!
!! The type contains the whole trailing edge info

type t_tedge

 !> Global id of the elements at the TE
 !integer , allocatable :: e(:,:)
 type(t_pot_elem_p), allocatable :: e(:,:)

 !> Global id of the nodes on the TE
 integer , allocatable :: i(:,:)

 !> Coordinates of the nodes of the TE
 !real(wp), allocatable :: rr(:,:)

 !> TE id of the nodes of the TE elements
 integer , allocatable :: ii(:,:)

 !> TE id of neighboring TE elements
 integer , allocatable :: neigh(:,:)

 !> Relative orientation of the neighboring TE elements
 integer , allocatable :: o(:,:)

 !> Unit vector at TE nodes
 real(wp), allocatable :: t(:,:)

 !> Reference frame of the TE nodes
 integer , allocatable :: ref(:)


end type t_tedge

!----------------------------------------------------------------------

interface destroy_elements

  module procedure destroy_elements_geo
  module procedure destroy_elements_comps

end interface 


!----------------------------------------------------------------------

character(len=*), parameter :: this_mod_name = 'mod_geo'

!----------------------------------------------------------------------

contains

!----------------------------------------------------------------------

!> Create all the goemetry
!!
!! The geometry creation proceeds in this way:
!! -# The references frames are read from the file and built
!! -# The components are built by loading them from the geometry file
!! -# The number of elements of each kind are counted and the geometrical
!!    quantities of each element are calculated
!! -# The pointer arrays are sorted with the static elements before and the
!!    moving elements after, and in the total elements surface panels
!!    and vortex rings before, then lifting lines and finally actuator disks
subroutine create_geometry(geo_file_name, ref_file_name, in_file_name,  geo, &
                      te, elems_impl, elems_expl, elems_ad, elems_ll, &
                    elems_tot, airfoil_data, sim_param, target_file, run_id)
 character(len=*), intent(in) :: geo_file_name
 character(len=*), intent(inout) :: ref_file_name
 character(len=*), intent(in) :: in_file_name
 type(t_geo), intent(out), target :: geo
 type(t_impl_elem_p), allocatable, intent(out) :: elems_impl(:)
 type(t_expl_elem_p), allocatable, intent(out) :: elems_expl(:)
 type(t_expl_elem_p), allocatable, intent(out) :: elems_ad(:)
 type(t_expl_elem_p), allocatable, intent(out) :: elems_ll(:)
 type(t_pot_elem_p),  allocatable, intent(out) :: elems_tot(:)
 type(t_tedge), intent(out) :: te
 type(t_aero_tab) , allocatable, intent(out) :: airfoil_data(:)
 type(t_sim_param) , intent(inout) :: sim_param
 character(len=*), intent(in) :: target_file
 integer, intent(in)              :: run_id(10)

 real(wp)                     :: tstart

 integer :: i, j, is, im,  i_comp, i_ll, i_ad, i_tot, i_expl, i_impl
 type(t_impl_elem_p), allocatable :: temp_static_i(:), temp_moving_i(:)
 type(t_expl_elem_p), allocatable :: temp_static_e(:), temp_moving_e(:)
 integer(h5loc) :: floc

 integer :: iSurfPan, iVortLatt
 
 character(len=max_char_len) :: msg
 character(len=*), parameter :: this_sub_name='create_geometry'

  tstart = sim_param%t0

  !build the reference frames
  !read which is the reference frame file, if default the main input file is
  !employed
  !reference_file = getstr(prms, 'ReferenceFile')
  if (trim(ref_file_name) .eq. 'no_set') then
    ref_file_name = trim(in_file_name)
  endif
  call check_file_exists(trim(ref_file_name), this_sub_name, this_mod_name)
  call build_references(geo%refs, trim(ref_file_name), sim_param)

! -----------------------------------------------------------------------------------------
! The following lines were needed when the motion was defined w.r.t. the initial
! time of the simulation, tstart, and not w.r.t. 0, as they are defined now
! 
! ! If rstart_from_file, set the right initial conditions!
! ! -> Correction needed for motion defined through its velocity
! if ( sim_param%restart_from_file .and. sim_param%reset_time ) then
!   call update_relative_initial_conditions(sim_param%restart_file, sim_param%ReferenceFile, geo%refs) 
! end if 
! -----------------------------------------------------------------------------------------

  !Create the output geometry file (if it is not restarting with the same name)
  if(trim(geo_file_name).ne.trim(target_file)) then
    call new_hdf5_file(trim(target_file), floc)
    call write_hdf5_attr(run_id, 'run_id', floc)
    call close_hdf5_file(floc)
  endif

  !Load the components from the file created by the preprocessor
  call check_preproc(geo_file_name)
  call load_components(geo, trim(geo_file_name), trim(target_file), &
                       sim_param, te)

  call import_aero_tab(geo,airfoil_data)

  ! Initialisation
  geo%nelem_impl   = 0
  geo%nelem_expl   = 0
  geo%nstatic_impl = 0
  geo%nmoving_impl = 0
  geo%nstatic_expl = 0
  geo%nmoving_expl = 0
  geo%nstatic_SurfPan = 0
  geo%nSurfPan   = 0
  geo%nVortLatt  = 0
  geo%nLiftLin   = 0
  geo%nActDisk   = 0

  ! count the elements
  do i_comp = 1,size(geo%components)

    select case(trim(geo%components(i_comp)%comp_el_type))
     case('p')
      if(geo%components(i_comp)%moving) then
        geo%nmoving_impl = geo%nmoving_impl + geo%components(i_comp)%nelems
      else
        geo%nstatic_impl = geo%nstatic_impl + geo%components(i_comp)%nelems
        geo%nstatic_SurfPan = geo%nstatic_SurfPan + geo%components(i_comp)%nelems
      endif
      geo%nSurfPan = geo%nSurfPan + geo%components(i_comp)%nelems
      geo%nelem_impl = geo%nelem_impl + geo%components(i_comp)%nelems

     case('v')
      if(geo%components(i_comp)%moving) then
        geo%nmoving_impl = geo%nmoving_impl + geo%components(i_comp)%nelems
      else
        geo%nstatic_impl = geo%nstatic_impl + geo%components(i_comp)%nelems
      endif
      geo%nVortLatt = geo%nVortLatt + geo%components(i_comp)%nelems
      geo%nelem_impl = geo%nelem_impl + geo%components(i_comp)%nelems

     case('l')
      if(geo%components(i_comp)%moving) then
        geo%nmoving_expl = geo%nmoving_expl + geo%components(i_comp)%nelems
      else
        geo%nstatic_expl = geo%nstatic_expl + geo%components(i_comp)%nelems
      endif
      geo%nLiftLin = geo%nLiftLin + geo%components(i_comp)%nelems
      geo%nelem_expl = geo%nelem_expl + geo%components(i_comp)%nelems
     
     case('a')
      if(geo%components(i_comp)%moving) then
        geo%nmoving_expl = geo%nmoving_expl + geo%components(i_comp)%nelems
      else
        geo%nstatic_expl = geo%nstatic_expl + geo%components(i_comp)%nelems
      endif
      geo%nActDisk = geo%nActDisk + geo%components(i_comp)%nelems
      geo%nelem_expl = geo%nelem_expl + geo%components(i_comp)%nelems

     case default
      call error (this_sub_name, this_mod_name, 'Unknow type of element: '&
                  //trim(geo%components(i_comp)%comp_el_type))
    end select

  enddo

  ! calculate the geometric quantities
  ! already update the geometry for the first time to get the right
  ! starting geometrical condition
  call prepare_geometry(geo)
  call update_geometry(geo, tstart, update_static=.true.)

  if(sim_param%debug_level .ge. 3) then
    call printout(nl//' Geometry details:' )
    write(msg,'(A,I9)') '  number of total elements:    ', &
                                                geo%nelem_impl + geo%nelem_expl
    call printout(msg)
    write(msg,'(A,I9)') '  number of implicit elements: ' ,geo%nelem_impl
    call printout(msg)
    write(msg,'(A,I9)') '  number of explicit elements: ' ,geo%nelem_expl
    call printout(msg)
    write(msg,'(A,I9)') '  number of static elements:   ' , &
                                            geo%nstatic_impl + geo%nstatic_expl
    call printout(msg)
    write(msg,'(A,I9)') '  number of moving elements:   ' , &
                                            geo%nmoving_impl + geo%nmoving_expl
    call printout(msg)
    write(msg,'(A,I9)') '  number of surface panels:    ' ,geo%nSurfpan
    call printout(msg)
    write(msg,'(A,I9)') '  number of vortex lattices:   ' ,geo%nVortLatt
    call printout(msg)
    write(msg,'(A,I9)') '  number of lifting lines:     ' ,geo%nLiftLin
    call printout(msg)
    write(msg,'(A,I9)') '  number of actuator disks:    ' ,geo%nActDisk
    call printout(msg)

  endif

  !Create the vector of pointers to all the elements
  allocate(elems_impl(geo%nelem_impl), elems_expl(geo%nelem_expl),  &
           elems_tot(geo%nelem_impl+geo%nelem_expl))
  allocate(elems_ad(geo%nactdisk), elems_ll(geo%nLiftLin))
  i_impl=0; i_expl=0; i_ad=0; i_ll=0; i_tot=0
  do i_comp = 1,size(geo%components)

    if (trim(geo%components(i_comp)%comp_el_type) .eq. 'p' .or. &
        trim(geo%components(i_comp)%comp_el_type) .eq. 'v') then

      do j = 1,size(geo%components(i_comp)%el)
        i_impl = i_impl+1
        i_tot = i_tot+1
        elems_tot(i_tot)%p => geo%components(i_comp)%el(j)
        select type(el => geo%components(i_comp)%el(j)); class is(c_impl_elem)
        elems_impl(i_impl)%p => el
        end select
      enddo

    elseif (trim(geo%components(i_comp)%comp_el_type) .eq. 'l') then

      do j = 1,size(geo%components(i_comp)%el)
        i_ll = i_ll+1
        i_expl = i_expl+1
        i_tot = i_tot+1
        elems_tot(i_tot)%p => geo%components(i_comp)%el(j)
        select type(el => geo%components(i_comp)%el(j)); class is(c_expl_elem)
          elems_ll(i_ll)%p => el
          elems_expl(i_expl)%p => el
        end select
      enddo

    elseif (trim(geo%components(i_comp)%comp_el_type) .eq. 'a') then

      do j = 1,size(geo%components(i_comp)%el)
        i_ad = i_ad+1
        i_expl = i_expl+1
        i_tot = i_tot+1
        elems_tot(i_tot)%p => geo%components(i_comp)%el(j)
        select type(el => geo%components(i_comp)%el(j)); class is(c_expl_elem)
          elems_ad(i_ad)%p => el
          elems_expl(i_expl)%p => el
        end select
      enddo

    endif
  enddo

  ! Sort implicit elements: first static, then moving -------
  !fill in the two temporaries
  allocate(temp_static_i(geo%nstatic_impl), temp_moving_i(geo%nmoving_impl))
  is = 0; im = 0;
  do i = 1,geo%nelem_impl
    if(elems_impl(i)%p%moving) then
      im = im+1
      temp_moving_i(im) = elems_impl(i)
    else
      is = is+1
      temp_static_i(is) = elems_impl(i)
    endif
  enddo

  !Now might be more bombproof to deallocate and allocate, but for the moment..
  elems_impl(1:geo%nstatic_impl) = temp_static_i
  elems_impl(geo%nstatic_impl+1:geo%nelem_impl) = temp_moving_i

  deallocate(temp_static_i, temp_moving_i)

  !Update the indexing since we re-ordered the vector
  do i = 1,geo%nelem_impl
    elems_impl(i)%p%id = i
  end do

  ! Save global indices of implicit elements for slicing
  allocate(geo%idSurfPan(geo%nSurfPan),geo%idVortLatt(geo%nVortLatt))
  allocate(geo%idSurfPanG2L(geo%nelem_impl)) ; geo%idSurfPanG2L = 0
  iSurfPan = 0 ; iVortLatt = 0
  do i = 1,geo%nelem_impl
    select type( el => elems_impl(i)%p ); class is(t_surfpan)
      iSurfPan = iSurfPan + 1
      geo%idSurfPan(iSurfPan) = el%id
      geo%idSurfPanG2L(i    ) = iSurfPan
    end select
    select type( el => elems_impl(i)%p ); class is(t_vortlatt)
      iVortLatt= iVortLatt + 1
      geo%idVortLatt(iVortLatt) = el%id
    end select
  end do

  !Now re-order the explicit elements
  allocate(temp_static_e(geo%nstatic_expl), temp_moving_e(geo%nmoving_expl))
  is = 0; im = 0;
  do i = 1,geo%nelem_expl
    if(elems_expl(i)%p%moving) then
      im = im+1
      temp_moving_e(im) = elems_expl(i)
    else
      is = is+1
      temp_static_e(is) = elems_expl(i)
    endif
  enddo
  if(geo%nelem_expl .gt. 0) then
    elems_expl(1:geo%nstatic_expl) = temp_static_e
    elems_expl(geo%nstatic_expl+1:geo%nelem_expl) = temp_moving_e
  endif
  deallocate(temp_static_e, temp_moving_e)

  !Update the indexing since we re-ordered the vector
  do i = 1,geo%nelem_expl
    elems_expl(i)%p%id = i+geo%nelem_impl
  end do

  call create_strip_connectivity(geo)

  ! Initialisation of the field surf_vel to zero (for surfpan only)
  do i=1,geo%nelem_impl
    select type(el=>elems_impl(i)%p) ; class is(t_surfpan)
      el%surf_vel = 0.0_wp
      call el%create_local_velocity_stencil( &    ! elems_tot, &
              geo%refs(geo%components(el%comp_id)%ref_id)%R_g )
      ! chtls stencil
      call el%create_chtls_stencil( &             ! elems_tot, &
              geo%refs(geo%components(el%comp_id)%ref_id)%R_g )
    end select
  end do
  
end subroutine create_geometry

!----------------------------------------------------------------------

!> Load the components mesh
!!
!! Here all the components from the geometry file are loaded.
!! The components linked to multiple reference frames are multiplied and
!! attached to their relative reference frames.
!! The points and the elements are then genereated.
!! Finally the trailing edge is created
subroutine load_components(geo, in_file, out_file, sim_param, te)
 type(t_geo), intent(inout),target :: geo
 character(len=*), intent(in) :: in_file
 character(len=*), intent(in) :: out_file
 type(t_sim_param) :: sim_param
 type(t_tedge), intent(out) :: te


 integer :: i1, i2, i3
 integer, allocatable :: ee(:,:)
 real(wp), allocatable :: rr(:,:)
 character(len=max_char_len) :: comp_el_type, comp_name, comp_input
 integer :: points_offset, n_vert , elems_offset
 real(wp), allocatable :: points_tmp(:,:)
 character(len=max_char_len) :: ref_tag, ref_tag_m
 integer :: ref_id, iref
 character(len=max_char_len) :: msg, cname, cname_write
 integer(h5loc) :: floc, gloc, cloc , geo_loc , te_loc, cloc2
 integer(h5loc) :: floc_out, gloc_out
 integer :: n_comp, i_comp, n_comp_input, i_comp_input, n_comp_write
 integer :: n_mult, i_mult
 logical :: mult

 ! Connectivity and te structures
 integer , allocatable :: neigh(:,:)
 ! Lifting Line elements
 real(wp), allocatable :: normalised_coord_e(:,:)
 integer                 , allocatable :: i_airfoil_e(:,:)
 character(max_char_len) , allocatable :: airfoil_list(:)
 integer                 , allocatable :: nelem_span_list(:)
 ! Parametric elements
 integer :: par_nelems_span , par_nelems_chor
 ! trailing edge ------
 integer , allocatable :: e_te(:,:) , i_te(:,:) , ii_te(:,:)
 integer , allocatable :: neigh_te(:,:) , o_te(:,:)
 real(wp), allocatable :: t_te(:,:)
 integer :: ne_te , nn_te
 ! tmp arrays --------
 type(t_pot_elem_p) , allocatable :: e_te_tmp(:,:)
 integer, allocatable  :: i_te_tmp(:,:) , ii_te_tmp(:,:)
 integer , allocatable :: neigh_te_tmp(:,:) , o_te_tmp(:,:)
 real(wp), allocatable :: t_te_tmp(:,:)
 integer , allocatable :: ref_te_tmp(:)
 ! # n. elements and nodes at TE ( of the prev. comps)
 integer :: ne_te_prev , nn_te_prev
 real(wp) :: trac, rad
 logical :: rewrite_geo

 character(len=*), parameter :: this_sub_name = 'load_components'

  !Re-write the geometry only if it is not restarting from a previous run
  !with the same name
  rewrite_geo = .true.
  if(trim(in_file).eq.trim(out_file)) rewrite_geo = .false.

  call printout('Reading geometry from file "'//trim(in_file)//'"')
  call open_hdf5_file(trim(in_file),floc)
  call open_hdf5_group(floc,'Components',gloc)
  call read_hdf5(n_comp,'NComponents',gloc)

  ! First it is necessary to count how many actual components will 
  ! be present, including the multiples. Due to pointers pointing to
  ! allocatables the components array cannot be expanded later
  n_comp_input = n_comp

  do i_comp_input = 1,n_comp_input
    write(cname,'(A,I3.3)') 'Comp',i_comp_input
    call open_hdf5_group(gloc,trim(cname),cloc)
    call read_hdf5(ref_tag,'RefTag',cloc)
    !Look for the reference frame of the component
    ref_id = -1
    do iref = 0,ubound(geo%refs,1)
      if (trim(geo%refs(iref)%tag) .eq. trim(ref_tag)) then
        !set id
        ref_id = iref
      endif
    enddo
    !if not found the reference
    if (ref_id .lt. 0) then
      write(msg,'(A,I2,A,A,A)') 'For component ',i_comp, &
                   ' a reference with tag ',trim(ref_tag),' was not found'
      call error(this_sub_name, this_mod_name, msg)
    endif
    mult = .false.
    n_mult = 1
    mult = geo%refs(ref_id)%multiple
    if(mult) then
      n_mult = geo%refs(ref_id)%n_mult
      n_comp = n_comp+n_mult-1 !(one was already counted)
    endif
    call close_hdf5_group(cloc)
  enddo
  
  ! Allocate the components of the right full size
  allocate(geo%components(n_comp))

  !Open and prepare the output file (if not restarting with the same name)
  if(rewrite_geo) then
    call open_hdf5_file(trim(out_file),floc_out)
    call new_hdf5_group(floc_out,'Components',gloc_out)
    call write_hdf5(n_comp,'NComponents',gloc_out)
  endif

  elems_offset = 0
  
  !TODO check this
  n_comp_write = n_comp
  

  i_comp = 1
  do i_comp_input = 1,n_comp_input

    write(cname,'(A,I3.3)') 'Comp',i_comp_input
    call open_hdf5_group(gloc,trim(cname),cloc)


    ! ====== REFERENCES ======

    call read_hdf5(ref_tag,'RefTag',cloc)

    !Look for the reference frame of the component
    ref_id = -1
    do iref = 0,ubound(geo%refs,1)
      if (trim(geo%refs(iref)%tag) .eq. trim(ref_tag)) then
        !set id
        ref_id = iref
      endif
    enddo
    !if not found the reference
    if (ref_id .lt. 0) then
      write(msg,'(A,I2,A,A,A)') 'For component ',i_comp, &
                   ' a reference with tag ',trim(ref_tag),' was not found'
      call error(this_sub_name, this_mod_name, msg)
    endif


    !Check if it is a multiple reference
    mult = .false.
    n_mult = 1
    mult = geo%refs(ref_id)%multiple

    ! ====== READING =====
    call read_hdf5(comp_el_type,'ElType',cloc)
    call read_hdf5(comp_name,'CompName',cloc)
    call read_hdf5(comp_input,'CompInput',cloc)

    if(mult) then
      n_mult = geo%refs(ref_id)%n_mult
    endif

    !cycle on all the possible multiple references, if it is not multiple
    !it will be passed just once
    !TODO: consider wrapping in another subroutine
    do i_mult = 1,n_mult
      ref_tag_m = ref_tag
      !set the component id
      geo%components(i_comp)%comp_id = i_comp
      if(mult)  then
        !look again for the multiple reference
        write(ref_tag_m,'(A,I2.2)') trim(ref_tag)//'__',i_mult
        ref_id = -1
        do iref = 0,ubound(geo%refs,1)
          if (trim(geo%refs(iref)%tag) .eq. trim(ref_tag_m)) then
            !set id
            ref_id = iref
          endif
        enddo
        !if not found the reference
        if (ref_id .lt. 0) then
          write(msg,'(A,I2,A,A,A)') 'For component ',i_comp, &
                      ' a reference with tag ',trim(ref_tag_m),' was not found'
          call error(this_sub_name, this_mod_name, msg)
        endif
      endif


      geo%components(i_comp)%ref_id  = ref_id
      geo%components(i_comp)%ref_tag = trim(ref_tag_m)
      geo%components(i_comp)%moving  = geo%refs(ref_id)%moving


!     ! ====== READING =====
      geo%components(i_comp)%comp_el_type = trim(comp_el_type)
      geo%components(i_comp)%comp_input   = trim(comp_input  )
 
      if(mult) then
        write(geo%components(i_comp)%comp_name,'(A,I2.2)') trim(comp_name)&
                                                            &//'__',i_mult
      else
        geo%components(i_comp)%comp_name = trim(comp_name)
      endif

      ! Geometry --------------------------
      call open_hdf5_group(cloc,'Geometry',geo_loc)
      call read_hdf5_al(ee   ,'ee'   ,geo_loc)
      call read_hdf5_al(rr   ,'rr'   ,geo_loc)
      call read_hdf5_al(neigh,'neigh',geo_loc)
      ! element-specific reads
      if ( comp_el_type(1:1) .eq. 'l' ) then
        call read_hdf5_al(airfoil_list      ,'airfoil_list'      ,geo_loc)!(:)
        call read_hdf5_al(nelem_span_list   ,'nelem_span_list'   ,geo_loc)!(:)
        call read_hdf5_al(i_airfoil_e       ,'i_airfoil_e'       ,geo_loc)!(:,:)
        call read_hdf5_al(normalised_coord_e,'normalised_coord_e',geo_loc)!(:,:)
        allocate(geo%components(i_comp)%airfoil_list(size(airfoil_list)))
        geo%components(i_comp)%airfoil_list = airfoil_list
        allocate(geo%components(i_comp)%nelem_span_list(size(nelem_span_list)))
        geo%components(i_comp)%nelem_span_list = nelem_span_list
        allocate(geo%components(i_comp)%i_airfoil_e( &
              size(i_airfoil_e,1),size(i_airfoil_e,2)) )
        geo%components(i_comp)%i_airfoil_e = i_airfoil_e
        allocate(geo%components(i_comp)%normalised_coord_e( &
              size(normalised_coord_e,1),size(normalised_coord_e,2)))
        geo%components(i_comp)%normalised_coord_e = normalised_coord_e
      else if (comp_el_type(1:1) .eq. 'a') then
        call read_hdf5(trac,'Traction',cloc)
        call read_hdf5(rad,'Radius',cloc)
      end if

      ! for PARAMETRIC elements only:
      ! parametric_nelems_span , parametric_nelems_chor 
      if ( trim(comp_input) .eq. 'parametric' ) then
        call read_hdf5(par_nelems_span,'parametric_nelems_span',geo_loc)
        call read_hdf5(par_nelems_chor,'parametric_nelems_chor',geo_loc)
      end if
 
      call close_hdf5_group(geo_loc)


      ! Trailing Edge (not all elements build the trailing edge)
      if( comp_el_type(1:1) .eq. 'p' .or. &
          comp_el_type(1:1) .eq. 'v' .or. &
          comp_el_type(1:1) .eq. 'l') then
        call open_hdf5_group(cloc,'Trailing_Edge',te_loc)
        call read_hdf5_al(    e_te,    'e_te',te_loc)
        call read_hdf5_al(    i_te,    'i_te',te_loc)
        !call read_hdf5_al(   rr_te,   'rr_te',te_loc)
        call read_hdf5_al(   ii_te,   'ii_te',te_loc)
        call read_hdf5_al(neigh_te,'neigh_te',te_loc)
        call read_hdf5_al(    o_te,    'o_te',te_loc)
        call read_hdf5_al(    t_te,    't_te',te_loc)
        !call read_hdf5_al(  ref_te,  'ref_te',te_loc)
        call close_hdf5_group(te_loc)
      else
        allocate(e_te(0,0), i_te(2,0), ii_te(2,0), neigh_te(2,0), o_te(2,0),&
                 t_te(3,0))
      endif
      
      !Re-write all that was read in the new output file (cannot just copy 
      !the read file since the components order will be changed when emplying
      !the multiple components). If it is restarting with the same name
      !avoid write

      write(cname_write,'(A,I3.3)') 'Comp',geo%components(i_comp)%comp_id
      if(rewrite_geo) then
        call new_hdf5_group(gloc_out,trim(cname_write),cloc2)
        call write_hdf5(ref_id,'RefId',cloc2)
        call write_hdf5(trim(ref_tag_m),'RefTag',cloc2)
        call write_hdf5(trim(comp_el_type),'ElType',cloc2)
        call write_hdf5(trim(geo%components(i_comp)%comp_name), &
                                                            'CompName',cloc2)
        call write_hdf5(trim(geo%components(i_comp)%comp_input),&
                                                           'CompInput',cloc2)
        call new_hdf5_group(cloc2,'Geometry',geo_loc)

        call write_hdf5(ee   ,'ee'   ,geo_loc)
        call write_hdf5(rr   ,'rr'   ,geo_loc)
        call write_hdf5(neigh,'neigh',geo_loc)
        if ( comp_el_type(1:1) .eq. 'l' ) then
          call write_hdf5(airfoil_list      ,'airfoil_list'      ,geo_loc)
          call write_hdf5(nelem_span_list   ,'nelem_span_list'   ,geo_loc)
          call write_hdf5(i_airfoil_e       ,'i_airfoil_e'       ,geo_loc)
          call write_hdf5(normalised_coord_e,'normalised_coord_e',geo_loc)

        else if (comp_el_type(1:1) .eq. 'a') then
          call write_hdf5(trac,'Traction',cloc2)
          call write_hdf5(rad,'Radius',cloc2)

        endif

        if ( trim(comp_input) .eq. 'parametric' ) then
          !write HDF5 fields
          call write_hdf5(par_nelems_span,'parametric_nelems_span',geo_loc)
          call write_hdf5(par_nelems_chor,'parametric_nelems_chor',geo_loc)
        end if

        call close_hdf5_group(geo_loc)

        if( comp_el_type(1:1) .eq. 'p' .or. &
            comp_el_type(1:1) .eq. 'v' .or. &
            comp_el_type(1:1) .eq. 'l') then
          call new_hdf5_group(cloc2,'Trailing_Edge',te_loc)
          call write_hdf5(    e_te,    'e_te',te_loc)
          call write_hdf5(    i_te,    'i_te',te_loc)
          call write_hdf5(   ii_te,   'ii_te',te_loc)
          call write_hdf5(neigh_te,'neigh_te',te_loc)
          call write_hdf5(    o_te,    'o_te',te_loc)
          call write_hdf5(    t_te,    't_te',te_loc)
          call close_hdf5_group(te_loc)
        endif

        call close_hdf5_group(cloc2)
      endif

      ! ======= CREATING ELEMENTS ======

      ! --- treat the points ---
      if(allocated(geo%points)) then
        points_offset = size(geo%points,2)
      else
        points_offset = 0
      endif

      !store the read points into the local points
      allocate(geo%components(i_comp)%loc_points(3,size(rr,2)))
      geo%components(i_comp)%loc_points = rr

      !Now for the moments the points are stored here without moving them,
      !will be moved later, consider not storing them here at all
      allocate(points_tmp(3,size(rr,2)+points_offset))
      if (points_offset .gt. 0) points_tmp(:,1:points_offset) = geo%points
      points_tmp(:,points_offset+1:points_offset+size(rr,2)) = rr
      call move_alloc(points_tmp, geo%points)
      allocate(geo%components(i_comp)%i_points(size(rr,2)))
      geo%components(i_comp)%i_points = &
                         (/((i3),i3=points_offset+1,points_offset+size(rr,2))/)


      ! --- treat the elements ---

      !allocate the elements of the component of the right kind
      geo%components(i_comp)%nelems = size(ee,2)
      select case(trim(geo%components(i_comp)%comp_el_type))
       case('p')
        allocate(t_surfpan::geo%components(i_comp)%el(size(ee,2)))
       case('v')
        allocate(t_vortlatt::geo%components(i_comp)%el(size(ee,2)))
       case('l')
        allocate(t_liftlin::geo%components(i_comp)%el(size(ee,2)))
       case('a')
        allocate(t_actdisk::geo%components(i_comp)%el(size(ee,2)))
       case default
        call error(this_sub_name, this_mod_name, &
          'Unknown type of element: '//geo%components(i_comp)%comp_el_type)
      end select

      !fill (some) of the elements fields
      do i2=1,size(ee,2)

        !vertices
        n_vert = count(ee(:,i2).ne.0)
        allocate(geo%components(i_comp)%el(i2)%i_ver(n_vert))
        allocate(geo%components(i_comp)%el(i2)%neigh(n_vert))
        geo%components(i_comp)%el(i2)%n_ver = n_vert
        geo%components(i_comp)%el(i2)%i_ver(1:n_vert) = &
                                              ee(1:n_vert,i2) + points_offset
        do i3=1,n_vert
          if ( neigh(i3,i2) .ne. 0 ) then
            geo%components(i_comp)%el(i2)%neigh(i3)%p => &
                                    geo%components(i_comp)%el(neigh(i3,i2))
          else
            ! do nothing, keep the neighbour pointer not associated
            geo%components(i_comp)%el(i2)%neigh(i3)%p => null()
          end if
        end do

        !motion
        geo%components(i_comp)%el(i2)%moving = geo%components(i_comp)%moving

      enddo

      geo%components(i_comp)%nSurfPan = 0; geo%components(i_comp)%nVortRin = 0;
      if(comp_el_type(1:1) .eq. 'p') &
                                  geo%components(i_comp)%nSurfPan = size(ee,2)
      if(comp_el_type(1:1) .eq. 'v') &
                                  geo%components(i_comp)%nVortRin = size(ee,2)
      ! TODO: add nLiftLin field ???

      !If it is an actuator disk read the traction
      if(geo%components(i_comp)%comp_el_type(1:1) .eq. 'a') then
        select type (el=>geo%components(i_comp)%el)
        type is(t_actdisk)
          do i2 = 1,size(el)
            el(i2)%traction = trac
            el(i2)%radius = rad
          enddo
        end select
      endif



      ! Trailing Edge ------------
      ne_te = 0; nn_te=0
      if(allocated(e_te)) ne_te = size(e_te,2)
      if(allocated(i_te)) nn_te = size(i_te,2)
      if (sim_param%debug_level .ge. 3) then
        write(msg,'(A,I0,A,I0)') ' Trailing edge: elements: ', &
                                                       ne_te, ' nodes ', nn_te
        call printout(msg)
      endif
      if (.not.allocated(te%e) ) then ! it should be enough
        allocate(te%e    (2,ne_te) )
        do i1 = 1,ne_te
          te%e(1,i1)%p => null()
          te%e(2,i1)%p => null()
          te%e(1,i1)%p  => geo%components(i_comp)%el(e_te(1,i1))
          if(e_te(2,i1) .gt. 0) &
            te%e(2,i1)%p  => geo%components(i_comp)%el(e_te(2,i1))
        enddo
        allocate(te%i    (2,nn_te) ) ; te%i     =     i_te
        !allocate(te%rr   (3,nn_te) ) ; te%rr    =    rr_te
        allocate(te%ii   (2,ne_te) ) ; te%ii    =    ii_te
        allocate(te%neigh(2,ne_te) ) ; te%neigh = neigh_te
        allocate(te%o    (2,ne_te) ) ; te%o     =     o_te
        allocate(te%t    (3,nn_te) ) ; te%t     =     t_te
        allocate(te%ref  (  nn_te) ) ; te%ref   = geo%components(i_comp)%ref_id

      elseif (ne_te .gt. 0) then
        nn_te_prev = size(te%i,2)
        ne_te_prev = size(te%e,2)
        allocate(e_te_tmp(2,size(te%e,2)+ne_te))
        e_te_tmp(:,             1:size(te%e,2)    ) = te%e
        do i1 = 1,ne_te
          e_te_tmp(1,size(te%e,2)+i1)%p => null()
          e_te_tmp(2,size(te%e,2)+i1)%p => null()
          e_te_tmp(1,size(te%e,2)+i1)%p  => &
                                       geo%components(i_comp)%el(e_te(1,i1))
          if(e_te(2,i1) .gt. 0) then
            e_te_tmp(2,size(te%e,2)+i1)%p  => &
                                       geo%components(i_comp)%el(e_te(2,i1))
          endif

        enddo
        call move_alloc(e_te_tmp,te%e)
        allocate(i_te_tmp(2,size(te%i,2)+nn_te))
        i_te_tmp(:,             1:size(te%i,2)    ) = te%i
        i_te_tmp(:,size(te%i,2)+1:size(i_te_tmp,2)) = i_te + points_offset
        call move_alloc(i_te_tmp,te%i)
        !allocate(rr_te_tmp(3,size(te%rr,2)+nn_te))
        !rr_te_tmp(:,              1:size(te%rr,2)    ) = te%rr
        !rr_te_tmp(:,size(te%rr,2)+1:size(rr_te_tmp,2)) = rr_te
        !call move_alloc(rr_te_tmp,te%rr)
        allocate(ii_te_tmp(2,size(te%ii,2)+ne_te))
        ii_te_tmp(:,              1:size(te%ii,2)    ) = te%ii
        ii_te_tmp(:,size(te%ii,2)+1:size(ii_te_tmp,2)) = ii_te + nn_te_prev
        call move_alloc(ii_te_tmp,te%ii)
        allocate(neigh_te_tmp(2,size(te%neigh,2)+ne_te))
        neigh_te_tmp(:,                 1:size(te%neigh,2)    ) = te%neigh
        where ( neigh_te .ne. 0 ) neigh_te = neigh_te + ne_te_prev
        neigh_te_tmp(:,size(te%neigh,2)+1:size(neigh_te_tmp,2)) = neigh_te
        call move_alloc(neigh_te_tmp,te%neigh)
        allocate( o_te_tmp(2,size(te%o ,2)+ne_te))
        o_te_tmp(:,              1:size(te%o ,2)    ) = te%o
        o_te_tmp(:,size(te%o ,2)+1:size( o_te_tmp,2)) =  o_te
        call move_alloc(o_te_tmp,te%o )
        allocate( t_te_tmp(3,size(te%t ,2)+nn_te))
        t_te_tmp(:,              1:size(te%t ,2)    ) = te%t
        t_te_tmp(:,size(te%t ,2)+1:size( t_te_tmp,2)) =  t_te
        call move_alloc(t_te_tmp,te%t )
        allocate(ref_te_tmp(size(te%ref   )+nn_te))
        ref_te_tmp(                 1:size(te%ref   )   ) =  te%ref
        ref_te_tmp(  size(te%ref  )+1:size(ref_te_tmp  )) = &
                                                  geo%components(i_comp)%ref_id

        call move_alloc(ref_te_tmp,te%ref)



      end if


      ! 2018-07-5. Deassociate neighboring elements through the TE
      if(trim(geo%components(i_comp)%comp_el_type) .eq. 'p') then
      do i1 = 1,ne_te
        do i2 = 1,geo%components(i_comp)%el(e_te(1,i1))%n_ver
          if ( neigh(i2,e_te(1,i1)) .eq. e_te(2,i1) ) then
            geo%components(i_comp)%el(e_te(1,i1))%neigh(i2)%p => null() 
          end if
        end do
        do i2 = 1,geo%components(i_comp)%el(e_te(2,i1))%n_ver
          if ( neigh(i2,e_te(2,i1)) .eq. e_te(1,i1) ) then
            geo%components(i_comp)%el(e_te(2,i1))%neigh(i2)%p => null() 
          end if
        end do
      end do
      endif

      deallocate(e_te, i_te, ii_te, neigh_te, o_te, t_te)
      !:::::::::::::::::::::::::::::::::::::::::::::::::::


      ! Update elems_offset for the next component
      elems_offset = elems_offset + size(ee,2)

      !cleanup
      deallocate(ee,rr,neigh)
      i_comp = i_comp + 1
    enddo !i_mult

    call close_hdf5_group(cloc)

  enddo !i_comp
  !update the total number of components
  !call write_hdf5(i_comp-1,'NComponents',gloc)
  call close_hdf5_group(gloc)
  call close_hdf5_file(floc)
  if(rewrite_geo) then
    call close_hdf5_group(gloc_out)
    call close_hdf5_file(floc_out)
  endif

  ! define comp_id for all the elems
  do i_comp = 1 , size(geo%components)
    do i1 = 1 , size(geo%components(i_comp)%el)
      geo%components(i_comp)%el(i1)%comp_id = i_comp
    end do
  end do

end subroutine load_components

!----------------------------------------------------------------------

subroutine import_aero_tab(geo,coeff)
 type(t_geo), intent(inout), target :: geo
 type(t_aero_tab) , allocatable , intent(inout) :: coeff(:)

 integer :: n_tmp , n_tmp2
 character(len=max_char_len) , allocatable :: list_tmp(:)
 character(len=max_char_len) , allocatable :: list_tmp_tmp(:)
 integer , allocatable :: i_airfoil_e_tmp (:,:)

 integer :: i_c, n_c, i_a, n_a, i_l

 n_tmp = 30
 allocate(list_tmp(n_tmp))

 ! Count # of different airfoil
 n_c = size(geo%components)
 n_a = 0
 do i_c = 1 ,  n_c

   if ( geo%components(i_c)%comp_el_type .eq. 'l' ) then

     allocate( i_airfoil_e_tmp( &
          size(geo%components(i_c)%i_airfoil_e,1), &
          size(geo%components(i_c)%i_airfoil_e,2) ) )
     i_airfoil_e_tmp = 0

     do i_a = 1 , size(geo%components(i_c)%airfoil_list)

       if (all( geo%components(i_c)%airfoil_list(i_a) .ne. &
                                                       list_tmp(1:n_a))) then
         ! new airfoil ----
         n_a = n_a + 1

         where (geo%components(i_c)%i_airfoil_e .eq. i_a) i_airfoil_e_tmp = n_a

         ! if n_a > n_tmp --> movalloc
         if ( n_a .gt. n_tmp ) then
           n_tmp2 = n_tmp + n_tmp
           allocate(list_tmp_tmp(n_tmp2))
           list_tmp_tmp(1:n_tmp) = list_tmp
           deallocate(list_tmp)
           call move_alloc(list_tmp_tmp,list_tmp)
           n_tmp = n_tmp2
         end if

         list_tmp(n_a) = geo%components(i_c)%airfoil_list(i_a)

       else
         ! airfoil already used: find the element and replace the global index
         do i_l = 1 , n_a
           if ( geo%components(i_c)%airfoil_list(i_a) .eq. list_tmp(i_l) ) exit
         end do

         where (geo%components(i_c)%i_airfoil_e .eq. i_a) i_airfoil_e_tmp = i_l

       end if


     end do

   geo%components(i_c)%i_airfoil_e = i_airfoil_e_tmp

   deallocate(i_airfoil_e_tmp)

   end if

 end do

 ! Read tables and fill coeff structure
 allocate(coeff(n_a))
 do i_a = 1 , n_a
   call read_c81_table( list_tmp(i_a) , coeff(i_a) )
 end do

end subroutine import_aero_tab

!----------------------------------------------------------------------

!> Prepare the geometry allocating all the relevant fields for each
!! kind of element
subroutine prepare_geometry(geo)
 type(t_geo), intent(inout), target :: geo

 integer :: i_comp, ie
 integer :: nsides
 class(c_pot_elem), pointer :: elem
 character(len=*), parameter :: this_sub_name = 'prepare_geometry'


 do i_comp = 1,size(geo%components)

   do ie = 1,size(geo%components(i_comp)%el)
     elem => geo%components(i_comp)%el(ie)

     nsides = size(elem%i_ver)

     !Fields common to each element
     allocate(elem%ver(3,nsides))
     allocate(elem%edge_vec(3,nsides))
     allocate(elem%edge_len(nsides))
     allocate(elem%edge_uni(3,nsides))
     
     !type-specific fields
     select type(elem)
      !Surface panel
      class is(t_surfpan)
       allocate(elem%verp(3,nsides))
       allocate(elem%cosTi(nsides))
       allocate(elem%sinTi(nsides))

      class is(t_liftlin)

       allocate(elem%tang_cen(3))
       allocate(elem%bnorm_cen(3))

       elem%csi_cen = 0.5_wp * &
                      sum(geo%components(i_comp)%normalised_coord_e(:,ie))
       elem%i_airfoil =  geo%components(i_comp)%i_airfoil_e(:,ie)

      class is(t_vortlatt)

      class is(t_actdisk)

      class default
       call error(this_sub_name, this_mod_name, 'Unknown element type')
     end select

   enddo
 enddo

end subroutine prepare_geometry

!----------------------------------------------------------------------

!> Calculate the geometrical quantities of an element
!!
!! The subroutine calculates all the relevant geometrical quantities of a
!! panel element (vortex ring or surface panel)
subroutine calc_geo_data_pan(elem,vert)
 class(c_pot_elem), intent(inout) :: elem
 real(wp), intent(in) :: vert(:,:)

 integer :: nsides, is
 real(wp):: nor(3), tanl(3)

  nsides = size(vert,2)

  elem%n_ver = nsides

  ! vertices
  elem%ver = vert

  ! center
  elem%cen =  sum ( vert,2 ) / real(nsides,wp)

  ! unit normal and area
  if ( nsides .eq. 4 ) then
    nor = cross( vert(:,3) - vert(:,1) , &
                 vert(:,4) - vert(:,2)     )
  else if ( nSides .eq. 3 ) then
    nor = cross( vert(:,3) - vert(:,2) , &
                 vert(:,1) - vert(:,2)     )
  end if

  elem%area = 0.5_wp * norm2(nor)
  elem%nor = nor / norm2(nor)

  ! local tangent unit vector as in PANAIR
  tanl = 0.5_wp * ( vert(:,nsides) + vert(:,1) ) - elem%cen

  elem%tang(:,1) = tanl / norm2(tanl)
  elem%tang(:,2) = cross( elem%nor, elem%tang(:,1)  )

! >>>>> USELESS LOOP: did we forget something?? Realised on 2018.11.14
!  select type(elem)
!  type is(t_surfpan)
!  do is = 1 , nsides
!  end do
!  end select
! <<<<< USELESS LOOP: did we forget something??

  ! vector connecting two consecutive vertices:
  ! edge_vec(:,1) =  ver(:,2) - ver(:,1)
  if ( nsides .eq. 3 ) then
    do is = 1 , nsides
      elem%edge_vec(:,is) = vert(:,next_tri(is)) - vert(:,is)
    end do
  else if ( nsides .eq. 4 ) then
    do is = 1 , nsides
      elem%edge_vec(:,is) = vert(:,next_qua(is)) - vert(:,is)
    end do
  end if

  ! edge: edge_len(:)
  do is = 1 , nsides
    elem%edge_len(is) = norm2(elem%edge_vec(:,is))
  end do

  ! unit vector
  do is = 1 , nSides
    elem%edge_uni(:,is) = elem%edge_vec(:,is) / elem%edge_len(is)
  end do

  !surface panels own fields
  select type(elem)
  type is(t_surfpan)
    do is = 1 , nsides
      ! cosTi , sinTi
      elem%cosTi(is) = sum( elem%edge_uni(:,is) * elem%tang(:,1) )
      elem%sinTi(is) = sum( elem%edge_uni(:,is) * elem%tang(:,2) )
      ! projection of the vertices on the mean plane
      elem%verp(:,is) = vert(:,is) - elem%nor * &
                      sum( (vert(:,is) - elem%cen ) * elem%nor )
    end do
  end select

  ! initialise %dforce
  elem%dforce = 0.0_wp

  ! initialise %dmom  
  elem%dmom   = 0.0_wp



end subroutine calc_geo_data_pan

!----------------------------------------------------------------------

!> Calculate the geometrical quantities of a lifting line element
!!
!! The subroutine calculates all the relevant geometrical quantities of a
!! lifting line element
subroutine calc_geo_data_ll(elem,vert)
 class(c_pot_elem), intent(inout) :: elem
 real(wp), intent(in) :: vert(:,:)

 integer :: is, nsides
 real(wp):: nor(3), tanl(3)

  nsides = size(vert,2)

  elem%n_ver = nsides

  ! vertices
  elem%ver = vert

  ! center, for the lifting line is the mid-point
  elem%cen =  sum ( vert(:,1:2),2 ) / 2.0_wp

  ! unit normal and area
  if ( nsides .eq. 4 ) then
    nor = cross( vert(:,3) - vert(:,1) , &
                 vert(:,4) - vert(:,2)     )
  else if ( nSides .eq. 3 ) then
    nor = cross( vert(:,3) - vert(:,2) , &
                 vert(:,1) - vert(:,2)     )
  end if

  elem%area = 0.5_wp * norm2(nor)
  elem%nor = nor / norm2(nor)

  ! local tangent unit vector as in PANAIR
  tanl = 0.5_wp * ( vert(:,nsides) + vert(:,1) ) - elem%cen

  elem%tang(:,1) = tanl / norm2(tanl)
  elem%tang(:,2) = cross( elem%nor, elem%tang(:,1)  )


  ! vector connecting two consecutive vertices:
  ! edge_vec(:,1) =  ver(:,2) - ver(:,1)
  if ( nsides .eq. 3 ) then
    do is = 1 , nsides
      elem%edge_vec(:,is) = vert(:,next_tri(is)) - vert(:,is)
    end do
  else if ( nsides .eq. 4 ) then
    do is = 1 , nsides
      elem%edge_vec(:,is) = vert(:,next_qua(is)) - vert(:,is)
    end do
  end if

  ! edge: edge_len(:)
  do is = 1 , nsides
    elem%edge_len(is) = norm2(elem%edge_vec(:,is))
  end do

  ! unit vector
  do is = 1 , nSides
    elem%edge_uni(:,is) = elem%edge_vec(:,is) / elem%edge_len(is)
  end do

  ! ll-specific fields
  select type(elem)
  type is(t_liftlin)
  elem%tang_cen = elem%edge_uni(:,2) - elem%edge_uni(:,4)
  elem%tang_cen = elem%tang_cen / norm2(elem%tang_cen)

  elem%bnorm_cen = cross(elem%tang_cen, elem%nor)
  elem%bnorm_cen = elem%bnorm_cen / norm2(elem%bnorm_cen)
  elem%chord = sum(elem%edge_len((/2,4/)))*0.5_wp
  end select

  ! initialise %dforce
  elem%dforce = 0.0_wp

  ! initialise %dmom  
  elem%dmom   = 0.0_wp

end subroutine calc_geo_data_ll

!----------------------------------------------------------------------

!> Calculate the geometrical quantities of an actuator disk
!!
!! The subroutine calculates all the relevant geometrical quantities of an
!! actuator disk
subroutine calc_geo_data_ad(elem,vert)
 class(c_pot_elem), intent(inout) :: elem
 real(wp), intent(in) :: vert(:,:)

 integer :: nsides, is
 real(wp):: nor(3), tanl(3)
 integer :: nxt

  nsides = size(vert,2)

  elem%n_ver = nsides

  ! vertices
  elem%ver = vert

  ! center, for the lifting line is the mid-point
  elem%cen =  sum ( vert,2 ) / real(nsides,wp)

  elem%area = 0.0_wp; elem%nor = 0.0_wp
  do is = 1, nsides
    nxt = 1+mod(is,nsides)
    nor = cross(vert(:,is) - elem%cen,&
                vert(:,nxt) - elem%cen )
    elem%area = elem%area + 0.5_wp * norm2(nor)
    elem%nor = elem%nor + nor/norm2(nor)
  enddo
    elem%nor = elem%nor/real(nsides,wp)

  ! local tangent unit vector: aligned with first node, normal to n
  tanl = (vert(:,1)-elem%cen)-sum((vert(:,1)-elem%cen)*elem%nor)*elem%nor

  elem%tang(:,1) = tanl / norm2(tanl)
  elem%tang(:,2) = cross( elem%nor, elem%tang(:,1)  )

  ! vector connecting two consecutive vertices:
  do is = 1 , nsides
    nxt = 1+mod(is,nsides)
    elem%edge_vec(:,is) = vert(:,nxt) - vert(:,is)
  end do

  ! edge: edge_len(:)
  do is = 1 , nsides
    elem%edge_len(is) = norm2(elem%edge_vec(:,is))
  end do

  ! unit vector
  do is = 1 , nsides
    elem%edge_uni(:,is) = elem%edge_vec(:,is) / elem%edge_len(is)
  end do

  ! initialise %dforce
  elem%dforce = 0.0_wp

  ! initialise %dmom  
  elem%dmom   = 0.0_wp

end subroutine calc_geo_data_ad

!----------------------------------------------------------------------

!> Calculate the local velocity on the panels to then enforce the
!! boundary condition
!!
subroutine calc_geo_vel(elem, G, f)
 class(c_pot_elem), intent(inout) :: elem
 real(wp), intent(in) :: f(3), G(3,3)

  elem%ub = f + matmul(G,elem%cen)

end subroutine calc_geo_vel

!----------------------------------------------------------------------

!> Calculate velocity of a point whose coordinate is rr
!! boundary condition
!!
subroutine calc_node_vel( r, G, f, v)
 real(wp), intent(in) :: r(3)
 real(wp), intent(in) :: f(3), G(3,3)
 real(wp), intent(out):: v(3)

 v = f + matmul(G,r)

end subroutine calc_node_vel

!----------------------------------------------------------------------

! Updatad node-element connectivity ( ee array )
! It is assumed that nodes and elements are sorted as follows
!
!  1-----5-----9-----13----17----21   <--- LE
!  |  1  |  4  |  7  |  10 |  13 |
!  2-----6-----10----14----18----22
!  |  2  |  5  |  8  |  11 |  14 |
!  3-----7-----11----15----19----23
!  |  3  |  6  |  9  |  12 |  15 |
!  4-----8-----12----16----20----24   <--- TE
! and the ee array is built as  5, 1, 2, 6
!                               6, 2, 3, 7
!                               7, 3, 4, 8
!                               9, 5, 6,10
!                                , ...
!                              23,19,20,24
subroutine create_strip_connectivity(geo)
 type(t_geo),  target, intent(inout) :: geo

 real(wp) :: h2

 integer :: io_te , io_tip
 integer :: n_el , ie_ind
 integer :: i_comp , i_el

 integer :: i_c , i_s , n_s , n_c , i_c2
 character(len=*), parameter :: this_sub_name = 'create_strip_connectivity'

 do i_comp = 1 , size(geo%components)
  associate(comp => geo%components(i_comp))

  ! connectivity built for lifting surface only
  if ( ( geo%components(i_comp)%comp_el_type(1:1) .eq. 'v' .or. &
         geo%components(i_comp)%comp_el_type(1:1) .eq. 'l' )  ) then

    n_el = size(comp%el)

    ! Some checks added: the first element must be at the LE,
    ! the last at the TE
    if ( associated(comp%el(  1 )%neigh(1)%p)) then
      call error(this_sub_name, this_mod_name, &
                'First element must be at the LE')
    end if
    if ( associated(comp%el(n_el)%neigh(3)%p) ) then
      call error(this_sub_name, this_mod_name, &
                 'Last element must be at the TE')
    end if

    ! initialisation
    io_tip = 1 ; io_te = 0 ; n_s = 0
    do i_el = 1 , n_el

     ie_ind = comp%el(i_el)%id

     if ( .not. associated(comp%el(i_el)%neigh(1)%p) ) then
       comp%el(i_el)%stripe_1%p => null()
       n_s = n_s + 1
     else
       comp%el(i_el)%stripe_1%p => comp%el(i_el)%neigh(1)%p
     end if

     comp%el(i_el)%dy = sum( cross( comp%el(i_el)%edge_uni(:,2) ,    &
                                    comp%el(i_el)%edge_vec(:,3)  ) * &
                             comp%el(i_el)%nor )

     h2 = sum( cross( comp%el(i_el)%edge_uni(:,2) ,                        &
                      comp%el(i_el)%ver(:,1) - comp%el(i_el)%ver(:,3) ) * &
               comp%el(i_el)%nor )

     if ( abs( h2 - comp%el(i_el)%dy ) .gt. 1e-4 .and. &
          sim_param%debug_level .ge. 7) then
       call warning(this_sub_name, this_mod_name, 'non parallel stripe edges (rough check)')
     end if

    end do

    ! allocate and fill comp%strip_elem array
    if ( mod(n_el,n_s) .ne. 0 ) then
      call error(this_sub_name, this_mod_name, ' The number of elements of&
         & a parametric element is not an exact multiple of the number of&
         & spanwise stripes. There is something wrong in the geometry input&
         & file')
    end if
    n_c = n_el / n_s  ! integer division to find number of chord panels

    do i_s = 1 , n_s
      do i_c = 1 , n_c
        allocate( comp%el(i_c+(i_s-1)*n_c)%stripe_elem(i_c) )
        do i_c2 = 1 , i_c
          comp%el(i_c+(i_s-1)*n_c)%stripe_elem(i_c2)%p => &
                                        comp%el(i_c2+(i_s-1)*n_c)
        end do
      end do
    end do

  end if

  end associate
 end do


end subroutine create_strip_connectivity

!----------------------------------------------------------------------

function move_points(pp, R, of)  result(rot_pp)
 real(wp), intent(in) :: pp(:,:)
 real(wp), intent(in)    :: R(:,:)
 real(wp), intent(in)    :: of(:)
 real(wp) :: rot_pp(size(pp,1),size(pp,2))

  rot_pp = matmul(R,pp)
  rot_pp(1,:) =rot_pp(1,:) + of(1)
  rot_pp(2,:) =rot_pp(2,:) + of(2)
  rot_pp(3,:) =rot_pp(3,:) + of(3)

end function move_points

!----------------------------------------------------------------------

!> Update all the geometry
!!
!! The geometry is updated at the beginning of the simulation, or after
!! a movement of the geometry. At the beginning the geometrical data of
!! all the elements is calculated, after only the moving elements are
!! updated.
!!
!! Also the velocity of the centerpoint of the elements is calculated
subroutine update_geometry(geo, t, update_static)
 type(t_geo), intent(inout) :: geo
 real(wp), intent(in) :: t
 logical, intent(in) :: update_static

 integer :: i_comp, ie

 !update all the references
 call update_all_references(geo%refs,t)

 do i_comp = 1,size(geo%components)
  associate(comp => geo%components(i_comp))

    if (comp%moving .or. update_static) then
      geo%points(:,comp%i_points) = move_points(comp%loc_points, &
                           geo%refs(comp%ref_id)%R_g, &
                           geo%refs(comp%ref_id)%of_g)

      do ie = 1,size(comp%el)
        !comp%el(ie)%ver = move_points(comp%loc_points(:,comp%el(ie)%i_ver), &
        !                   geo%refs(comp%ref_id)%R_g, &
        !                   geo%refs(comp%ref_id)%of_g)
        call comp%el(ie)%calc_geo_data(geo%points(:,comp%el(ie)%i_ver))
      enddo

! 2018.11.15 Moved outside, at the end of create_geometry
! !     !in the first pass compute also the velocity stencil for surfpans
! !     if(update_static) then
!         select type(els=>comp%el); class is(t_surfpan)
!         do ie = 1,size(els)
!           call els(ie)%create_local_velocity_stencil()
!         enddo
!         end select
! !     endif


      do ie = 1,size(comp%el)
        !Calculate the velocity of the centers to impose
        !the boundary condition
        call calc_geo_vel(comp%el(ie), geo%refs(comp%ref_id)%G_g, &
                                geo%refs(comp%ref_id)%f_g)
      enddo

    end if

  end associate
 enddo

end subroutine update_geometry

!----------------------------------------------------------------------

!> Destroy the pointers array of the components
!!
!! This subroutine is needed just for compatibility with gcc 4.8
!! Newer versions of the compiler and ifort deallocate also the 
!! vector of pointers alongside the geometry in destroy_geometry,
!! however gcc 4.8 gives a SIGSEGV and needs and explicit deallocation
subroutine destroy_elements_geo(geo)
 type(t_geo), intent(inout) :: geo

 integer :: ic

 do ic=1,size(geo%components)
   if(allocated(geo%components(ic)%el)) deallocate(geo%components(ic)%el)
 enddo

end subroutine destroy_elements_geo

!----------------------------------------------------------------------

!> Destroy the pointers array of the components
!!
!! This subroutine is needed just for compatibility with gcc 4.8
!! Newer versions of the compiler and ifort deallocate also the 
!! vector of pointers alongside the geometry in destroy_geometry,
!! however gcc 4.8 gives a SIGSEGV and needs and explicit deallocation
subroutine destroy_elements_comps(components)
 type(t_geo_component), intent(inout) :: components(:)

 integer :: ic

 do ic=1,size(components)
   if(allocated(components(ic)%el)) deallocate(components(ic)%el)
 enddo

end subroutine destroy_elements_comps

!----------------------------------------------------------------------

!> Destroy all the contents of a geometry type
!!
subroutine destroy_geometry(geo, elems)
 type(t_geo), intent(out) :: geo
 type(t_pot_elem_p), allocatable, intent(out) :: elems(:)

 integer :: i
 call destroy_references(geo%refs)

 !dummy to avoid warnings
 i = size(elems)


end subroutine destroy_geometry

!----------------------------------------------------------------------




!----------------------------------------------------------------------

end module mod_geometry
