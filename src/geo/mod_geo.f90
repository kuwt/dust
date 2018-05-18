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
  t_sim_param

use mod_parse, only: &
  t_parse, getstr, getint, getreal, getrealarray, getlogical, countoption &
  , finalizeparameters

use mod_handling, only: &
  error, warning, info, printout, dust_time, t_realtime, check_preproc

use mod_basic_io, only: &
  read_mesh_basic, write_basic

use mod_cgns_io, only: &
  read_mesh_cgns

use mod_parametric_io, only: &
  read_mesh_parametric

use mod_aero_elements, only: &
  c_elem, t_elem_p

use mod_surfpan, only: &
  t_surfpan

use mod_vortring, only: &
  t_vortring

use mod_liftlin, only: &
  t_liftlin

use mod_actuatordisk, only: &
  t_actdisk

use mod_c81, only: &
  t_aero_tab , read_c81_table , interp_aero_coeff

!use read_naca00xx, only: &
!  read_naca_arrays

use mod_math, only: &
  cross

use mod_reference, only: &
  t_ref, build_references, update_all_references, destroy_references
  
use mod_hdf5_io, only: &
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

!----------------------------------------------------------------------

implicit none

public :: t_geo, t_geo_component, t_tedge, set_parameters_geo, &
          create_geometry, update_geometry, destroy_geometry, &
          calc_geo_data_pan, calc_geo_data_ll

private

!----------------------------------------------------------------------

!> Element type
!!
!! Employed to build the geometry
type :: t_e
 
 !> Vertices indexes
 integer, allocatable :: vert(:)

 !> element type
 character :: etype

end type t_e

!-----------------------------------

!> Geometry component
!!
!! A geometry component is a single geometrical component of the model,
!! such as a wing, an aileron, a rotor blade etc.
!! 
!! It is used to load the mesh, build the geometry and in the future to
!! move the geometry.
!!
!! It also contains some temporary arrays employed during geometry generation
type :: t_geo_component

 !> Element type
 character(len=max_char_len) :: comp_el_type

 !> Name of the component
 character(len=max_char_len) :: comp_name

 !> Id of the component (warning: not always defined)
 integer :: comp_id

 !> Number of elements of the component
 integer :: nelems

 !> Elements of the group
 !! The main memory allocation happens here, they will be pointed from
 !! somewhere else
 class(c_elem), allocatable :: el(:)

 !> Global indexes of the points in the component
 integer, allocatable :: i_points(:)

 !> Reference frame tag 
 !integer :: ref_tag
 character(len=max_char_len) :: ref_tag

 ! Reference frame id
 integer :: ref_id
 
 !> Points in local reference frame
 real(wp), allocatable :: loc_points(:,:)

 !> Temporary to build the component
 !type(t_e), allocatable :: temp_elems(:)
 
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

! old idea ---> now each element knwos the indices of the previous elements
!               on the same strip
!!> Strip for vortex ring (structured) elements
!integer , allocatable :: strip_elem(:,:)

 !> Arrays for lifting line components
 character(len=max_char_len), allocatable :: airfoil_list(:)
 !> Number of ll elements in each region. Size(nelem_span_list) = nRegions
 integer, allocatable :: nelem_span_list(:)
 !> Normalised coordinate of sections in each wing region. Size(...) = nRegions+1
 real(wp),allocatable :: normalised_coord_e(:,:)
 !> Id of the airfoil elements (index in airfoil_list char array)
 integer ,allocatable :: i_airfoil_e(:,:)
 
!character(len=max_char_len), allocatable :: airfoil_table_p(:,:)

end type  t_geo_component

!-----------------------------------

!> Geometry group
!!
!! A geometry group is a set of elemenst/component logically grouped.
!!
!! Will be employed to calculate loads
type :: t_geo_group

 type(t_elem_p), allocatable :: e(:)

 real(wp) :: ref_base(3,3)

end type  t_geo_group

!-----------------------------------


!> Geometry type
!!
!! The type contains the whole geometry, with geometrical data and (probably)
!! also all the solution
type :: t_geo

 !> Total number of implicit elements
 integer :: nelem

 !> Number of lifting line elements
 integer :: nll

 !> Number of actuator disk elements
 integer :: nad

 !> Number of surface panel elements
 integer :: nSurfPan

 !> Number of vortex ring elements
 integer :: nVortRin

 !> Number of lifting line elements
 integer :: nLiftLin

 !> Number of Actuator disk elements
 integer :: nActDisk

 !> Number of statical elements
 integer :: nstatic
 !> Number of moving elements
 integer :: nmoving

 !> Number of static or moving lifting lines
 integer :: nstatic_ll, nmoving_ll

 !> Number of static or moving lifting lines
 integer :: nstatic_ad, nmoving_ad

 !> All the components of the geometry
 type(t_geo_component), allocatable :: components(:)

 !> All the groups of the geometry
 type(t_geo_group), allocatable :: groups(:)

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
 type(t_elem_p), allocatable :: e(:,:)

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

!-----------------------------------

character(len=*), parameter :: this_mod_name = 'mod_geo'
!----------------------------------------------------------------------

contains

!----------------------------------------------------------------------

!> Set the main parser parameters regarding the geometry input
subroutine set_parameters_geo(prms)
 type(t_parse), intent(inout) :: prms

  call prms%CreateStringOption('GeometryFile','Main geometry definition file')
  call prms%CreateStringOption('ReferenceFile','Reference frames file','no_set')

end subroutine set_parameters_geo

!----------------------------------------------------------------------

!> Create all the goemetry
!! 
!! The geometry creation proceeds in this way:
!! -# The main geometry file is read, and the components therein contained are
!!    created (reading the mesh)
!! -# If some additional geometry files are present, each one is read and all
!!    its component are created
!! -# The main geometry arrays, SurfPan and VortRin are allocated and filled
!!    in the prepare_geometry routine
!! -# The element pointer array used to build/solve the linear system is
!!    created, pointed at each element and then re-ordered with the static
!!    elements first, and the dynamic elements at the end
subroutine create_geometry(geo_file_name, ref_file_name, in_file_name,  geo, &
                           te, elems, elems_ll, elems_ad, elems_tot, &
                           airfoil_data, sim_param)
 character(len=*), intent(in) :: geo_file_name
 character(len=*), intent(inout) :: ref_file_name
 character(len=*), intent(in) :: in_file_name
 type(t_geo), intent(out), target :: geo
 type(t_elem_p), allocatable, intent(out) :: elems(:)
 type(t_elem_p), allocatable, intent(out) :: elems_ll(:)
 type(t_elem_p), allocatable, intent(out) :: elems_ad(:)
 type(t_elem_p), allocatable, intent(out) :: elems_tot(:)
 type(t_tedge), intent(out) :: te
 type(t_aero_tab) , allocatable, intent(out) :: airfoil_data(:)
 type(t_sim_param) , intent(inout) :: sim_param
 real(wp)                     :: tstart

 !character(len=max_char_len) :: reference_file
 !character(len=max_char_len) :: geo_file_name

 integer :: i, j, is, im,  i_comp, i_ll, i_ad, i_tot
 type(t_elem_p), allocatable :: temp_static(:), temp_moving(:)

 !integer , allocatable :: el_id_old(:), el_id_old_static(:), el_id_old_moving(:)

 character(len=max_char_len) :: msg

  tstart = sim_param%t0
 
  !build the reference frames
  !read which is the reference frame file, if default the main input file is 
  !employed
  !reference_file = getstr(prms, 'ReferenceFile')
  if (trim(ref_file_name) .eq. 'no_set') then
    ref_file_name = trim(in_file_name)
  endif

  call build_references(geo%refs, trim(ref_file_name), sim_param)

  !Load the components from the file created by the preprocessor
  !geo_file_name = getstr(prms, 'GeometryFile')
  call check_preproc(geo_file_name)
  call load_components(geo, trim(geo_file_name), sim_param, te)

  call import_aero_tab(geo,airfoil_data)

  ! Initialisation
  geo%nelem      = 0
  geo%nll        = 0
  geo%nad        = 0
  geo%nstatic    = 0
  geo%nmoving    = 0
  geo%nstatic_ll = 0
  geo%nmoving_ll = 0
  geo%nstatic_ad = 0
  geo%nmoving_ad = 0
  geo%nSurfPan   = 0
  geo%nVortRin   = 0
  geo%nLiftLin   = 0
  geo%nActDisk   = 0

  ! count the elements
  do i_comp = 1,size(geo%components)

    if (trim(geo%components(i_comp)%comp_el_type) .eq. 'p') then
      if(geo%components(i_comp)%moving) then
        geo%nmoving = geo%nmoving + geo%components(i_comp)%nelems
      else
        geo%nstatic = geo%nstatic + geo%components(i_comp)%nelems
      endif
      geo%nSurfPan = geo%nSurfPan + geo%components(i_comp)%nelems
      geo%nelem = geo%nelem + geo%components(i_comp)%nelems

    elseif (trim(geo%components(i_comp)%comp_el_type) .eq. 'v') then
      if(geo%components(i_comp)%moving) then
        geo%nmoving = geo%nmoving + geo%components(i_comp)%nelems
      else
        geo%nstatic = geo%nstatic + geo%components(i_comp)%nelems
      endif
      geo%nVortRin = geo%nVortRin + geo%components(i_comp)%nelems
      geo%nelem = geo%nelem + geo%components(i_comp)%nelems

    elseif (trim(geo%components(i_comp)%comp_el_type) .eq. 'l') then
      if(geo%components(i_comp)%moving) then
        geo%nmoving_ll = geo%nmoving_ll + geo%components(i_comp)%nelems
      else
        geo%nstatic_ll = geo%nstatic_ll + geo%components(i_comp)%nelems
      endif
      geo%nLiftLin = geo%nLiftLin + geo%components(i_comp)%nelems
      geo%nll = geo%nll + geo%components(i_comp)%nelems
    elseif (trim(geo%components(i_comp)%comp_el_type) .eq. 'a') then
      if(geo%components(i_comp)%moving) then
        geo%nmoving_ad = geo%nmoving_ad + geo%components(i_comp)%nelems
      else
        geo%nstatic_ad = geo%nstatic_ad + geo%components(i_comp)%nelems
      endif
      geo%nActDisk = geo%nActDisk + geo%components(i_comp)%nelems
      geo%nad = geo%nad + geo%components(i_comp)%nelems
    endif
    
  enddo

  ! calculate the geometric quantities
  ! already update the geometry for the first time to get the right 
  ! starting geometrical condition
  call prepare_geometry(geo)
  call update_geometry(geo, tstart, update_static=.true.)

  if(sim_param%debug_level .ge. 3) then
    call printout(nl//' Geometry details:' ) 
    write(msg,'(A,I9)') '  number of elements:        ' ,geo%nelem
    call printout(msg)
    write(msg,'(A,I9)') '  number of static elements: ' ,geo%nstatic
    call printout(msg)
    write(msg,'(A,I9)') '  number of moving elements: ' ,geo%nmoving
    call printout(msg)
    write(msg,'(A,I9)') '  number of surface panels:  ' ,geo%nsurfpan
    call printout(msg)
    write(msg,'(A,I9)') '  number of vortex rings:    ' ,geo%nvortrin
    call printout(msg)
    write(msg,'(A,I9)') '  number of lifting lines:   ' ,geo%nvortrin
    call printout(msg)
    write(msg,'(A,I9)') '  number of actuator disks:  ' ,geo%nvortrin
    call printout(msg)
  endif

  !Create the vector of pointers to all the elements
  allocate(elems(geo%nelem), elems_ll(geo%nll), elems_ad(geo%nad), &
           elems_tot(geo%nelem+geo%nll+geo%nad)) 
  i=0; i_ll=0; i_ad=0; i_tot=0
  do i_comp = 1,size(geo%components)
    
    if (trim(geo%components(i_comp)%comp_el_type) .eq. 'p' .or. &
        trim(geo%components(i_comp)%comp_el_type) .eq. 'v') then

      do j = 1,size(geo%components(i_comp)%el)
        i = i+1
        i_tot = i_tot+1
        elems(i)%p => geo%components(i_comp)%el(j)
      enddo

    elseif (trim(geo%components(i_comp)%comp_el_type) .eq. 'l') then

      do j = 1,size(geo%components(i_comp)%el)
        i_ll = i_ll+1
        i_tot = i_tot+1
        elems_ll(i_ll)%p => geo%components(i_comp)%el(j)
      enddo

    elseif (trim(geo%components(i_comp)%comp_el_type) .eq. 'a') then

      do j = 1,size(geo%components(i_comp)%el)
        i_ad = i_ad+1
        i_tot = i_tot+1
        elems_ad(i_ad)%p => geo%components(i_comp)%el(j)
      enddo

    endif
  enddo

  ! Sort elements: first static, then moving ------- 
  !fill in the two temporaries
  allocate(temp_static(geo%nstatic), temp_moving(geo%nmoving))
  !allocate(el_id_old(geo%nelem))          ; el_id_old = 0
  !allocate(el_id_old_static(geo%nstatic)) ; el_id_old_static = 0
  !allocate(el_id_old_moving(geo%nmoving)) ; el_id_old_moving = 0
  is = 0; im = 0;
  do i = 1,geo%nelem
    if(elems(i)%p%moving) then
      im = im+1
      temp_moving(im) = elems(i)
      !el_id_old_moving(im) = i
    else
      is = is+1
      temp_static(is) = elems(i)
      !el_id_old_static(is) = i
    endif
  enddo

  !Now might be more bombproof to deallocate and allocate, but for the moment..
  elems(1:geo%nstatic) = temp_static
  elems(geo%nstatic+1:geo%nelem) = temp_moving

  !el_id_old(1:geo%nstatic) = el_id_old_static
  !el_id_old(geo%nstatic+1:geo%nelem) = el_id_old_moving
 
  !Update the indexing since we re-ordered the vector
  do i = 1,geo%nelem
    elems(i)%p%id = i
  end do
  !Update elem-elem connectivity after re-ordering, for the moment only for
  !implicit panel elements
  !do i = 1,geo%nelem
  !  do j = 1,elems(i)%p%n_ver
  !    if ( elems(i)%p%i_neigh(j) .ne. 0 ) then
  !      elems(i)%p%i_neigh(j) = el_id_old( elems(i)%p%i_neigh(j) )
  !    else
  !      elems(i)%p%i_neigh(j) = 0 
  !    end if     
  !  end do 
  !end do 
  !Update te structure
  !do i = 1,size(te%e,2)
  !  do j = 1,2
  !    if ( te%e(j,i) .ne. 0 ) then
  !      te%e(j,i) = el_id_old( te%e(j,i) )
  !    else
  !      te%e(j,i) = 0 
  !    end if     
  !  end do 
  !end do

  !Now re-order the lifting line elements
  deallocate(temp_static, temp_moving)
  allocate(temp_static(geo%nstatic_ll), temp_moving(geo%nmoving_ll))
  is = 0; im = 0;
  do i = 1,geo%nll
    if(elems_ll(i)%p%moving) then
      im = im+1
      temp_moving(im) = elems_ll(i)
    else
      is = is+1
      temp_static(is) = elems_ll(i)
    endif
  enddo
  if(geo%nll .gt. 0) then
    elems_ll(1:geo%nstatic_ll) = temp_static
    elems_ll(geo%nstatic_ll+1:geo%nll) = temp_moving
  endif

  !Now re-order the actuator disks
  deallocate(temp_static, temp_moving)
  allocate(temp_static(geo%nstatic_ad), temp_moving(geo%nmoving_ad))
  is = 0; im = 0;
  do i = 1,geo%nad
    if(elems_ad(i)%p%moving) then
      im = im+1
      temp_moving(im) = elems_ad(i)
    else
      is = is+1
      temp_static(is) = elems_ad(i)
    endif
  enddo
  if(geo%nad .gt. 0) then
    elems_ad(1:geo%nstatic_ad) = temp_static
    elems_ad(geo%nstatic_ad+1:geo%nad) = temp_moving
  endif

  deallocate(temp_static, temp_moving)
  !deallocate(el_id_old, el_id_old_static, el_id_old_moving)

  !Patch together everything in elems_tot
  elems_tot(1:geo%nelem) = elems
  elems_tot(geo%nelem+1:geo%nelem+geo%nll) =elems_ll
  elems_tot(geo%nelem+geo%nll+1:geo%nelem+geo%nll+geo%nad) = elems_ad

  call create_local_velocity_stencil(geo,elems)    ! for surfpan only (3dP)

  call create_strip_connectivity(geo)

end subroutine create_geometry

!----------------------------------------------------------------------

subroutine load_components(geo, in_file, sim_param, te)
 type(t_geo), intent(inout),target :: geo
 character(len=*), intent(in) :: in_file
 type(t_sim_param) :: sim_param
 type(t_tedge), intent(out) :: te


 type(t_geo_component), allocatable :: comp_temp(:)
 integer :: i1, i2, i3
 integer, allocatable :: ee(:,:)
 real(wp), allocatable :: rr(:,:)
 character(len=max_char_len) :: comp_el_type, comp_name
 integer :: points_offset, n_vert , elems_offset
 real(wp), allocatable :: points_tmp(:,:)
 character(len=max_char_len) :: ref_tag, ref_tag_m
 integer :: ref_id, iref
 character(len=max_char_len) :: msg, cname
 integer(h5loc) :: floc, gloc, cloc , geo_loc , te_loc
 integer :: n_comp, i_comp, n_comp_input, i_comp_input
 integer :: n_mult, i_mult
 logical :: mult

 ! Connectivity and te structures 
 integer , allocatable :: neigh(:,:)
 ! Lifting Line elements
 real(wp), allocatable :: normalised_coord_e(:,:)
 integer                 , allocatable :: i_airfoil_e(:,:)
 character(max_char_len) , allocatable :: airfoil_list(:)
 integer                 , allocatable :: nelem_span_list(:)

 ! trailing edge ------
 integer , allocatable :: e_te(:,:) , i_te(:,:) , ii_te(:,:)
 integer , allocatable :: neigh_te(:,:) , o_te(:,:)
 real(wp), allocatable :: rr_te(:,:) , t_te(:,:)
 integer :: ne_te , nn_te
 ! tmp arrays --------
 type(t_elem_p) , allocatable :: e_te_tmp(:,:)
 integer, allocatable  ::i_te_tmp(:,:) , ii_te_tmp(:,:) 
 integer , allocatable :: neigh_te_tmp(:,:) , o_te_tmp(:,:)
 real(wp), allocatable ::rr_te_tmp(:,:) , t_te_tmp(:,:)
 integer , allocatable :: ref_te_tmp(:)
 integer :: ne_te_prev , nn_te_prev ! # n. elements and nodes at TE ( of the prev. comps) 
 real(wp) :: trac

 character(len=*), parameter :: this_sub_name = 'load_components'

  call open_hdf5_file(trim(in_file),floc)
  call open_hdf5_group(floc,'Components',gloc)
  call read_hdf5(n_comp,'NComponents',gloc)

  allocate(geo%components(n_comp))

  elems_offset = 0
  
  n_comp_input = n_comp
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
    if(mult) then
      n_mult = geo%refs(ref_id)%n_mult
      n_comp = n_comp+n_mult-1 !(one was already counted)
      allocate(comp_temp(n_comp))
      comp_temp(1:i_comp) = geo%components(1:i_comp)
      deallocate(geo%components)
      call move_alloc(comp_temp, geo%components)
    endif
    
    !cycle on all the possible multiple references, if it is not multiple
    !it will be passed just once 
    !TODO: consider wrapping in another subroutine
    do i_mult = 1,n_mult
      ref_tag_m = ref_tag
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

      !Re-write the reference tag to the component
      call write_hdf5(ref_id,'RefId',cloc)
        !TODO: fix this WARNING: it is not multiple references safe

    
 

      ! ====== READING =====

      call read_hdf5(comp_el_type,'ElType',cloc)
      geo%components(i_comp)%comp_el_type = trim(comp_el_type)

      call read_hdf5(comp_name,'CompName',cloc)
      geo%components(i_comp)%comp_name = trim(comp_name)


      ! Geometry --------------------------
      call open_hdf5_group(cloc,'Geometry',geo_loc)
      call read_hdf5_al(ee   ,'ee'   ,geo_loc)
      call read_hdf5_al(rr   ,'rr'   ,geo_loc)
      call read_hdf5_al(neigh,'neigh',geo_loc)
      ! element-specific reads 
      if ( comp_el_type(1:1) .eq. 'l' ) then
        call read_hdf5_al(airfoil_list      ,'airfoil_list'      ,geo_loc) ! (:)
        call read_hdf5_al(nelem_span_list   ,'nelem_span_list'   ,geo_loc) ! (:)
        call read_hdf5_al(i_airfoil_e       ,'i_airfoil_e'       ,geo_loc) ! (:,:)
        call read_hdf5_al(normalised_coord_e,'normalised_coord_e',geo_loc) ! (:,:)
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
        allocate(t_vortring::geo%components(i_comp)%el(size(ee,2)))
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
        
        !Component id
        geo%components(i_comp)%el(i2)%comp_id = i_comp
        
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
            !geo%components(i_comp)%el(i2)%i_neigh(i3) = &
            !                        neigh(i3,i2)
          end if
        end do
        
        !motion
        geo%components(i_comp)%el(i2)%moving = geo%components(i_comp)%moving
        allocate(geo%components(i_comp)%el(i2)%vel(3))

      enddo
      

      geo%components(i_comp)%nSurfPan = 0; geo%components(i_comp)%nVortRin = 0;
      if(comp_el_type(1:1) .eq. 'p') geo%components(i_comp)%nSurfPan = size(ee,2)
      if(comp_el_type(1:1) .eq. 'v') geo%components(i_comp)%nVortRin = size(ee,2)
! TODO: add nLiftLin field ???
!     if(comp_el_type(1:1) .eq. 'l') geo%components(i_comp)%nLiftLin = size(ee,2)

      !If it is an actuator disk read the traction
      if(geo%components(i_comp)%comp_el_type(1:1) .eq. 'a') then
        call read_hdf5(trac,'Traction',cloc)
        select type (el=>geo%components(i_comp)%el)
        type is(t_actdisk)
          do i2 = 1,size(el)
            el(i2)%traction = trac
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
      if (.not.allocated(te%e) .and. ne_te .gt. 0) then ! it should be enough
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
        allocate(te%t    (2,nn_te) ) ; te%t     =     t_te
        allocate(te%ref  (  nn_te) ) ; te%ref   =   geo%components(i_comp)%ref_id
        deallocate(e_te, i_te, ii_te, neigh_te, o_te, t_te)

      elseif (ne_te .gt. 0) then
        nn_te_prev = size(te%i,2)
        ne_te_prev = size(te%e,2)
        allocate(e_te_tmp(2,size(te%e,2)+ne_te)) 
        e_te_tmp(:,             1:size(te%e,2)    ) = te%e
        do i1 = 1,ne_te
          e_te_tmp(1,size(te%e,2)+i1)%p => null()
          e_te_tmp(2,size(te%e,2)+i1)%p => null()
          e_te_tmp(1,size(te%e,2)+i1)%p  => geo%components(i_comp)%el(e_te(1,i1)) 
          if(e_te(2,i1) .gt. 0) &
            e_te_tmp(2,size(te%e,2)+i1)%p  => geo%components(i_comp)%el(e_te(2,i1)) 
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
        deallocate(e_te, i_te, ii_te, neigh_te, o_te, t_te)
      end if 

      !:::::::::::::::::::::::::::::::::::::::::::::::::::
 
      ! Update elems_offset for the next component
      elems_offset = elems_offset + size(ee,2)

      !cleanup
      deallocate(ee,rr,neigh)
      i_comp = i_comp + 1
    enddo !i_mult

    call close_hdf5_group(cloc)

  enddo !i_comp
  call close_hdf5_group(gloc)
  call close_hdf5_file(floc)

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
          size(geo%components(i_c)%i_airfoil_e,1), size(geo%components(i_c)%i_airfoil_e,2) ) ) 
     i_airfoil_e_tmp = 0

     do i_a = 1 , size(geo%components(i_c)%airfoil_list)

       if ( all( geo%components(i_c)%airfoil_list(i_a) .ne. list_tmp(1:n_a) ) ) then
         ! new airfoil ----
         n_a = n_a + 1

         where ( geo%components(i_c)%i_airfoil_e .eq. i_a ) i_airfoil_e_tmp = n_a 
        
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

         where ( geo%components(i_c)%i_airfoil_e .eq. i_a ) i_airfoil_e_tmp = i_l 

       end if


     end do

   geo%components(i_c)%i_airfoil_e = i_airfoil_e_tmp

   deallocate(i_airfoil_e_tmp)
     
   ! check ----
   write(*,*) ' mod_geo.f89/import_aero_tab().  Component: ' , i_c
   write(*,*) ' size(geo%components(',i_c,')%i_airfoil_e : ' , shape(geo%components(i_c)%i_airfoil_e)
   do i_l = 1 , size(geo%components(i_c)%i_airfoil_e,2)
     write(*,*) geo%components(i_c)%i_airfoil_e(:,i_l)
   end do
   write(*,*)
   ! check ----

   end if

 end do

 ! Read tables and fill coeff structure
 allocate(coeff(n_a))
 write(*,*) ' Number of different airfoils : ' , n_a
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
 class(c_elem), pointer :: elem
 character(len=*), parameter :: this_sub_name = 'prepare_geometry'

 do i_comp = 1,size(geo%components)
   do ie = 1,size(geo%components(i_comp)%el)
     !associate(elem => geo%compnents(i_comp)%el(ie))
     elem => geo%components(i_comp)%el(ie)

     nsides = size(elem%i_ver)

     !Fields common to each element
     allocate(elem%ver(3,nsides))
     allocate(elem%cen(3))
     allocate(elem%nor(3))

     select type(elem)
      class is(t_surfpan)
       allocate(elem%tang(3,2))
       allocate(elem%verp(3,nsides))
       allocate(elem%edge_vec(3,nsides))
       allocate(elem%edge_len(nsides))
       allocate(elem%edge_uni(3,nsides))
       allocate(elem%cosTi(nsides))
       allocate(elem%sinTi(nsides))

      class is(t_vortring)
       allocate(elem%tang(3,2))
       allocate(elem%verp(3,nsides))
       allocate(elem%edge_vec(3,nsides))
       allocate(elem%edge_len(nsides))
       allocate(elem%edge_uni(3,nsides))
       allocate(elem%cosTi(nsides))
       allocate(elem%sinTi(nsides))

      class is(t_liftlin)
       allocate(elem%tang(3,2))
       allocate(elem%tang_cen(3))
       allocate(elem%bnorm_cen(3))
       allocate(elem%verp(3,nsides))
       allocate(elem%edge_vec(3,nsides))
       allocate(elem%edge_len(nsides))
       allocate(elem%edge_uni(3,nsides))
       allocate(elem%cosTi(nsides))
       allocate(elem%sinTi(nsides))

       
       
       elem%csi_cen = 0.5_wp * sum(geo%components(i_comp)%normalised_coord_e(:,ie))
       elem%i_airfoil =  geo%components(i_comp)%i_airfoil_e(:,ie)
       
       
       ! elem%chord    
       ! elem%csi_c       !! <- for interpolation of aerodynamic coefficients
       ! elem%airfoil(2)  !! 

       ! elem%theta   <-- ??? 

      class is(t_actdisk)
       allocate(elem%tang(3,2))
       allocate(elem%verp(3,nsides))
       allocate(elem%edge_vec(3,nsides))
       allocate(elem%edge_len(nsides))
       allocate(elem%edge_uni(3,nsides))
       allocate(elem%cosTi(nsides))
       allocate(elem%sinTi(nsides))

      class default
       call error(this_sub_name, this_mod_name, 'Unknown element type')
     end select

     !end associate
   enddo
 enddo

end subroutine prepare_geometry

!----------------------------------------------------------------------

!> Calculate the geometrical quantities of an element
!!
!! The subroutine calculates all the relevant geometrical quantities of a
!! panel element (vortex ring or surface panel) 
subroutine calc_geo_data_pan(elem,vert)
 class(c_elem), intent(inout) :: elem  
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

  ! projection of the vertices on the mean plane
  do is = 1 , nsides
    elem%verp(:,is) = vert(:,is) - elem%nor * &
                      sum( (vert(:,is) - elem%cen ) * elem%nor )
  end do
  
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

  ! cosTi , sinTi
  do is = 1 , nsides
    elem%cosTi(is) = sum( elem%edge_uni(:,is) * elem%tang(:,1) ) 
    elem%sinTi(is) = sum( elem%edge_uni(:,is) * elem%tang(:,2) ) 
  end do


end subroutine calc_geo_data_pan

!----------------------------------------------------------------------

!> Calculate the geometrical quantities of a lifting line element
!!
!! The subroutine calculates all the relevant geometrical quantities of a
!! lifting line element
subroutine calc_geo_data_ll(elem,vert)
 class(c_elem), intent(inout) :: elem  
 real(wp), intent(in) :: vert(:,:)

 integer :: nsides, is
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

  ! projection of the vertices on the mean plane
  do is = 1 , nsides
    elem%verp(:,is) = vert(:,is) - elem%nor * &
                      sum( (vert(:,is) - elem%cen ) * elem%nor )
  end do
  
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

  ! cosTi , sinTi
  do is = 1 , nsides
    elem%cosTi(is) = sum( elem%edge_uni(:,is) * elem%tang(:,1) ) 
    elem%sinTi(is) = sum( elem%edge_uni(:,is) * elem%tang(:,2) ) 
  end do


end subroutine calc_geo_data_ll

!----------------------------------------------------------------------

!> Calculate the geometrical quantities of an actuator disk
!!
!! The subroutine calculates all the relevant geometrical quantities of an
!! actuator disk
subroutine calc_geo_data_ad(elem,vert)
 class(c_elem), intent(inout) :: elem  
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

  ! projection of the vertices on the mean plane
  do is = 1 , nsides
    elem%verp(:,is) = vert(:,is) - elem%nor * &
                      sum( (vert(:,is) - elem%cen ) * elem%nor )
  end do
  
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
  

  ! cosTi , sinTi
  do is = 1 , nsides
    elem%cosTi(is) = sum( elem%edge_uni(:,is) * elem%tang(:,1) ) 
    elem%sinTi(is) = sum( elem%edge_uni(:,is) * elem%tang(:,2) ) 
  end do


end subroutine calc_geo_data_ad

!----------------------------------------------------------------------

!> Calculate the local velocity on the panels to then enforce the 
!! boundary condition
!!
subroutine calc_geo_vel(elem, G, f)
 class(c_elem), intent(inout) :: elem  
 real(wp), intent(in) :: f(3), G(3,3)

  if(.not.allocated(elem%ub)) allocate(elem%ub(3))
  elem%ub = f + matmul(G,elem%cen)

end subroutine calc_geo_vel

!----------------------------------------------------------------------

!> Subroutine used to double the mesh by reflecting it along a simmetry 
!! plane
!!
!! Given a plane defined by a center point and a normal vector, the mesh 
!! is doubled: all the points are reflected and new mirrored elements 
!! introduced. The elements and points arrays are doubled. 
!subroutine reflect_mesh(ee, rr, cent, norm)
! integer, allocatable, intent(inout) :: ee(:,:)
! real(wp), allocatable, intent(inout) :: rr(:,:)
! real(wp), intent(in) :: cent(3), norm(3)
!
! real(wp) :: n(3), d, l
! integer, allocatable :: ee_temp(:,:)
! real(wp), allocatable :: rr_temp(:,:)
! integer :: ip, np, ne
! integer :: ie, iv, nv
!
!  ne = size(ee,2); np = size(rr,2)
!  
!  ! enlarge size
!  allocate(ee_temp(size(ee,1),2*ne))
!  allocate(rr_temp(size(rr,1),2*np))
! 
!  !first part equal
!  ee_temp(:,1:ne) = ee
!  rr_temp(:,1:np) = rr
!  
!  !second part of the elements: index incremented, need to rearrange the 
!  !connectivity to preserve the normal
!  ee_temp(:,ne+1:2*ne) = 0
!  do ie = 1,ne
!    nv = count(ee(:,ie).ne.0)
!    ee_temp(1,ne+ie) = np+ee(1,ie)
!    do iv = 2,nv
!      ee_temp(iv,ne+ie) = np+ee(nv-iv+2,ie)
!    enddo
!  enddo
!  !ee_temp(:,ne+1:2*ne) = ee+np
! 
!  !calculate normal unit vector and distance from origin
!  n = norm/norm2(norm) 
!  l = sum(cent * n)
! 
!  !now reflect the points
!  do ip=1,np
!    d = sum( rr(:,ip) * n) - l 
!    rr_temp(:,np+ip) = rr_temp(:,ip) - 2*d*n
!  enddo
!
!  !move alloc back to the original vectors
!  call move_alloc(rr_temp, rr)
!  call move_alloc(ee_temp, ee)
!
!end subroutine reflect_mesh

!----------------------------------------------------------------------

!subroutine build_connectivity(elems)
! type(t_elem_p), intent(inout) :: elems(:)
!
! integer :: nelem, nvert
! integer :: ie1, iv1, ie2, iv2
! integer :: vert1, vert2
!
! character(len=max_char_len) :: msg
! real(t_realtime) :: t0, t1
!
! character(len=*), parameter :: this_sub_name = 'build_connectivity'
!
!  t0 = dust_time()
!
!  nelem = size(elems)
! 
!  do ie1 = 1,nelem
!    nvert = elems(ie1)%p%n_ver
!    allocate(elems(ie1)%p%i_neigh(nvert))
!    elems(ie1)%p%i_neigh = 0
!  enddo
!
!  do ie1 = 1,nelem
!    nvert = elems(ie1)%p%n_ver
!    do iv1 = 1,nvert
!      if(elems(ie1)%p%i_neigh(iv1) .ne. 0) cycle
!      vert1 = elems(ie1)%p%i_ver(iv1)
!      vert2 = elems(ie1)%p%i_ver(mod(iv1,nvert)+1)
!      do ie2 = ie1+1,nelem
!        if (any(elems(ie2)%p%i_ver .eq. vert1) .and. &
!                              any(elems(ie2)%p%i_ver .eq. vert2)) then
!          !assign first element neighbour
!          elems(ie1)%p%i_neigh(iv1) = elems(ie2)%p%id
!          ! assign the second in the correct position
!          ! the element HAS to have the nodes in the correct position,
!          ! with vert1 after vert2, indicating that the neighbouring 
!          ! element has the same normal orientation
!          do iv2 = 1,elems(ie2)%p%n_ver
!            if (elems(ie2)%p%i_ver(iv2) .eq. vert2) then
!              if (elems(ie2)%p%i_ver(mod(iv2, elems(ie2)%p%n_ver)+1) & 
!                                                            .eq.vert1) then
!                elems(ie2)%p%i_neigh(iv2) = elems(ie1)%p%id
!              else
!                ! debug ----
!                write(*,*) ' elements ie1 , ie2 , elems(ie1)%p%id , elems(ie2)%p%id : ' 
!                write(*,*)            ie1 , ie2 , elems(ie1)%p%id , elems(ie2)%p%id
!                write(*,*) ' elems. i1 ' , elems(ie1)%p%i_ver
!                write(*,*) ' elems. i2 ' , elems(ie2)%p%i_ver
!                write(*,*) ' vertices el.i1 ' , vert1 , vert2
!                write(*,*) ' vertices el.i2 ' , elems(ie2)%p%i_ver(iv1) , elems(ie2)%p%i_ver(iv2)
!                ! debug ----
!                call error(this_sub_name, this_mod_name, &
!                  'Neighbouring elements with opposed normal orientation,&
!                  &this is not allowed. Stop')
!              endif
!            endif
!          enddo
!        endif
!      enddo
!    enddo
!  enddo
!  t1 = dust_time()
!  !write(msg,'(A,F9.3,A)') 'Connectivity built in ', t1-t0, ' s.'
!  call printout(msg)
!
!end subroutine build_connectivity

!----------------------------------------------------------------------

!! TODO : 
!! - add a reference velocity to identify the trailing edge.
!!  Up to now, the relative velocity is assumed to be "mainly oriented" as
!!  the x-axis in the local reference frame. 
!! - identify the trailing edge for vortex ring elements
!! TODO: REMOVE THE SUBROUTINE it is not used anymore
!subroutine build_te(geo,elems,te)
! type(t_geo)   , intent(in), target :: geo
! type(t_elem_p), intent(in)  :: elems(:)
! type(t_tedge) , intent(out) :: te
! 
! real(wp) , parameter :: tol_sewing = 3e-3_wp   ! 3e-6_wp  ! 1e-0 for nasa-crm , 3e-3_wp for naca0012
!
! type(t_elem_p) , allocatable :: elems_m(:)
!
! integer :: i_comp , i_e
!
! write(*,*) nl//' Start TE connectivity '//nl  
!
! 
!
!
! do i_comp = 1,size(geo%components) 
!!  associate(comp => geo%components(i_comp))
!
!   select type(elements => geo%components(i_comp)%el)
!    type is(t_surfpan)
!   !if ( geo%components(i_comp)%comp_el_type(1:1) .eq. 'p' ) then
!
!
!     ! Merge nodes to obtain "closed-TE" connectivity ----- 
!     ! merge_nodes(rri,eei,tol, rr_m,ee_m,i_m)
!    
!     !call merge_nodes( geo%points , geo%components(i_comp)%i_points ,  &
!     !                  geo%components(i_comp)%el ,                      &
!     !                  tol_sewing , elems_m  )
!
!     call merge_nodes( geo%points , geo%components(i_comp)%i_points ,  &
!                       elements ,                      &
!                       tol_sewing , elems_m  )
!
!     ! Find "closed-TE" connectivity ---------------------- 
!     call build_connectivity( elems_m )
!
! 
!     ! Find TE -------------------------------------------- 
!     call find_te( geo , geo%components(i_comp) , elems , elems_m , te )
!
!     ! WARNING: deallocate elems_m !
!     do i_e = 1,size(elems_m)
!       deallocate(elems_m(i_e)%p) 
!     end do
!     deallocate(elems_m)
!
!    type is(t_vortring)
!   !elseif ( geo%components(i_comp)%comp_el_type(1:1) .eq. 'v' ) then
!  
!
!     ! Find TE -------------------------------------------- 
!     call find_te_v( geo , geo%components(i_comp) , elems , te ) ! No "closed-TE" connectivity is required
!
!    class default
!   !else
!     write(*,*) ' In build_te(). Component ', i_comp 
!     write(*,*) '  Wrong type. Stop '  ! 
!   
!   end select
!   !end if  
! 
!!  end associate
! end do
!
! write(*,*) nl//' End of TE connectivity '//nl  
!
!
!! -------------
!contains
!
!! 01. Inside build_te : merge_nodes ------
!subroutine merge_nodes( ri , indi , elemsi , tol ,   elems  )
!!type(t_geo)   , intent(in), target :: geo
! real(wp)      , intent(in) :: ri(:,:)
! integer       , intent(in) :: indi(:)
! !type(t_elem_p), intent(in) :: elemsi(:)
! class(c_elem), intent(in) :: elemsi(:)
! !type(t_surfpan), intent(in) :: elemsi(:)
! real(wp) :: tol
! type(t_elem_p), allocatable , intent(out) :: elems(:)
!! type(t_surfpan),allocatable , target      :: c_elems(:)
! integer       , allocatable  :: im(:,:)
! real(wp)      , allocatable  :: rr(:,:)
! integer  , allocatable :: im_tmp(:,:)
!
! integer :: n_nodes , n_merge
!
! integer :: in1 , in2
! integer :: i1 , i2 , i_e , i_v
!
! n_nodes = size(indi,1)
!
!
! allocate(elems(size(elemsi))) 
! do i1 = 1 , size(elemsi)
!   allocate(t_surfpan::elems(i1)%p)
!   allocate(elems(i1)%p%i_ver(elemsi(i1)%n_ver))
!   elems(i1)%p%id    = elemsi(i1)%id
!   elems(i1)%p%n_ver = elemsi(i1)%n_ver
!   elems(i1)%p%i_ver = elemsi(i1)%i_ver
!
! end do
!
!
! allocate(im_tmp(2,n_nodes) ) ; im_tmp = 0
! allocate(    rr(3,n_nodes) ) ; rr = ri(:,indi)
!
! n_merge = 0
! do i1 = 1 , n_nodes
!
!  in1 = indi(i1)
!
!  do i2 = i1 + 1 , n_nodes
!
!   in2 = indi(i2)
! 
!   if ( norm2(ri(:,in1)-ri(:,in2)) .lt. tol ) then
!    
!    n_merge = n_merge + 1 
! 
!    rr(:,i1) = 0.5_wp * (ri(:,in1) + ri(:,in2))
!    rr(:,i2) = rr(:,i1) !   0.5_wp * (ri(:,i1) + ri(:,i2))
!
!
!    do i_e = 1 , size(elems) 
!     do i_v = 1 , size(elemsi(i_e)%i_ver)
!
!      if ( elems(i_e)%p%i_ver(i_v) .eq. in2 ) then
!
!        if ( all(elems(i_e)%p%i_ver .ne. in1 ) ) then
!
!          elems(i_e)%p%i_ver(i_v) = in1 
!  
!        end if
!      end if 
!     end do
!    end do
! 
!    im_tmp(1,n_merge) = in1
!    im_tmp(2,n_merge) = in2
! 
!   end if
! 
!  end do
!
! end do
! write(*,*) ' n_merge : ' , n_merge
!
! allocate(im(2,n_merge)) ; im = im_tmp(:,1:n_merge)
!
! deallocate(im_tmp)
! deallocate(rr)
!
!end subroutine merge_nodes
!
!! 02. Inside build_te : merge_nodes ------
!
!! 03. Inside build_te : find_te ----------
!! WARNING !!!! subroutine not completed
!! Look for TODO ++++
!subroutine find_te( geo , comp , elems , elems_m , te )
! type(t_geo)   , intent(in), target :: geo
! type(t_geo_component), intent(in) :: comp
! type(t_elem_p)       , intent(in) :: elems(:)
! type(t_elem_p)       , intent(in) :: elems_m(:) ! read "actual" neighbors
! type(t_tedge)        , intent(inout) :: te
!
! ! tmp arrays --------
! integer , allocatable               :: e_te_tmp(:,:) , i_te_tmp(:,:) , ii_te_tmp(:,:) 
! integer , allocatable               :: neigh_te_tmp(:,:) , o_te_tmp(:,:)
! real(wp), allocatable               ::rr_te_tmp(:,:) , t_te_tmp(:,:)
! integer , allocatable               :: ref_te_tmp(:)
! integer , allocatable               :: i_el_nodes_tmp(:,:)
! ! actual arrays -----
! integer , allocatable               :: e_te(:,:) , i_te(:,:) , ii_te(:,:) 
! integer , allocatable               :: neigh_te(:,:) , o_te(:,:)
! real(wp), allocatable               :: rr_te(:,:) , t_te(:,:) 
! integer , allocatable               :: ref_te(:)
!
! integer :: n_el
! integer :: ne_te , nn_te           ! # n. elements and nodes at TE (for the actual comp)
! integer :: ne_te_prev , nn_te_prev ! # n. elements and nodes at TE ( of the prev. comps) 
!
! integer :: neigh_ib_ie , ie_ind
! integer  :: ind1 , ind2 , i_node1 , i_node2 , i_node1_max , i_node2_max
! real(wp) ::  mi1 , mi2 
!
! real(wp) , parameter :: inner_prod_thresh = - 0.5d0
!
! integer :: i_e , i_b , i_n
!
! integer :: nSides1 , nSides2
! integer :: i1 , i2 , e1 , e2
!
! real(wp) , dimension(3) :: vec1
!
! real(wp) , dimension(3) , parameter  :: u_rel = (/ 1.0_wp , 0.0_wp , 0.0_wp /) ! hard-coded parameter ... 
! real(wp) , dimension(3) :: v_rel = u_rel / norm2(u_rel)
! real(wp) , dimension(3) , parameter :: side_dir = (/ 0.0_wp , 1.0_wp , 0.0_wp /) ! hard-coded parameter ... 
! ! TODO: - read as an input of the component
! !       - Projection of the tangent unit vector with no side slip
!
! n_el = size(comp%el)
!
! ! allocate tmp structures --
! allocate( e_te_tmp(2,(n_el+1)/2) )                ;  e_te_tmp = 0 
! allocate( i_el_nodes_tmp(2,size(comp%i_points)) ) ;  i_el_nodes_tmp = 0
! allocate( i_te_tmp      (2,size(comp%i_points)) ) ;  i_te_tmp = 0
! allocate(ii_te_tmp      (2,n_el               ) ) ; ii_te_tmp = 0
! allocate(rr_te_tmp      (3,size(comp%i_points)) ) ; rr_te_tmp = 0.0d0
!
!
! ! initialise counters ------
! ne_te = 0 ; nn_te = 0
!
! do i_e = 1 , n_el
!
!
!   do i_b = 1 , comp%el(i_e)%n_ver
!
!    if ( ( elems_m(i_e)%p%i_neigh(i_b) .ne. 0 ) .and. &
!         ( all(e_te_tmp(1,1:ne_te) .ne. elems_m(i_e)%p%i_neigh(i_b)) ) ) then ! should be enough 
!
!     ! Check normals ----
!     ! 1.
!     if ( sum( comp%el(i_e)%nor * elems( elems_m(i_e)%p%i_neigh(i_b) )%p%nor ) &
!           < inner_prod_thresh ) then
!       ne_te = ne_te + 1
!
!       ! 2. ... other criteria to find te elements --------
!     
!       ! surface elements at the trailing edge ------------
!       neigh_ib_ie = elems_m(i_e)%p%i_neigh(i_b)
!       ie_ind      = elems_m(i_e)%p%id
!       e_te_tmp(1,ne_te) = ie_ind
!       e_te_tmp(2,ne_te) = neigh_ib_ie 
!
!!      ! nodes on the trailing edge ----
!!      i_el_nodes_tmp(1,nete) = ee( ib,ie )
!!      i_el_nodes_tmp(2,nete) = ee( mod(ib,size(neigh,1)) + 1 , ie ) ! <--- size(neigh) ****
!       i_el_nodes_tmp(1,ne_te) = elems(ie_ind)%p%i_ver(i_b)  
!       i_el_nodes_tmp(2,ne_te) = elems(ie_ind)%p%i_ver( mod(i_b,elems(ie_ind)%p%n_ver)+1 )   
!  
!!      ! Find the corresponding node on the other side of the traling edge (for open te)
!!      ind1 = 1 ; mi1 = norm2( rr(:,i_el_nodes_tmp(1,nete)) - rr(:,ee(1,neigh(ib,ie))) )
!!      ind2 = 1 ; mi2 = norm2( rr(:,i_el_nodes_tmp(2,nete)) - rr(:,ee(1,neigh(ib,ie))) )
!       ind1 = 1 ; mi1 = norm2( geo%points(:,i_el_nodes_tmp(1,ne_te)) - &
!                               geo%points(:, elems(neigh_ib_ie)%p%i_ver(1) ) ) 
!       ind2 = 1 ; mi2 = norm2( geo%points(:,i_el_nodes_tmp(2,ne_te)) - & 
!                               geo%points(:, elems(neigh_ib_ie)%p%i_ver(1) ) ) 
!
!!      do i1 = 2 , elems(ie_ind)%p%n_ver       ! <--- WRONG
!       do i1 = 2 , elems(neigh_ib_ie)%p%n_ver  ! <--- mod 2018-03-38 
!!        write(*,*) ' mi1        ' , mi1
!!        write(*,*) ' norm2(...) ' , norm2( rr(:,i_el_nodes_tmp(1,nete)) - rr(:,ee(i1,neigh(ib,ie) )))
!!        write(*,*) ' mi2        ' , mi1
!!        write(*,*) ' norm2(...) ' , norm2( rr(:,i_el_nodes_tmp(2,nete)) - rr(:,ee(i1,neigh(ib,ie) )))
!         if (     norm2( geo%points(:,i_el_nodes_tmp(1,ne_te)) - &  
!                         geo%points(:, elems(neigh_ib_ie)%p%i_ver(i1) ) ) .lt. mi1 ) then
!           ind1 = i1
!           mi1  = norm2( geo%points(:,i_el_nodes_tmp(1,ne_te)) - &  
!                         geo%points(:, elems(neigh_ib_ie)%p%i_ver(i1) ) )
!         end if
!         if (     norm2( geo%points(:,i_el_nodes_tmp(2,ne_te)) - &
!                         geo%points(:, elems(neigh_ib_ie)%p%i_ver(i1) ) ) .lt. mi2 ) then
!           ind2 = i1
!           mi2  = norm2( geo%points(:,i_el_nodes_tmp(2,ne_te)) - &
!                         geo%points(:, elems(neigh_ib_ie)%p%i_ver(i1) ) )  
!         end if
!       end do
! 
!       ! Select the minimum value of the corresponding nodes at the te, to avoid double nodes ...
!!      i_node1 = min(i_el_nodes_tmp(1,ne_te),ee(ind1,neigh(ib,ie)) )
!!      i_node2 = min(i_el_nodes_tmp(2,ne_te),ee(ind2,neigh(ib,ie)) )
!       i_node1     = min(i_el_nodes_tmp(1,ne_te), elems(neigh_ib_ie)%p%i_ver(ind1) )
!       i_node2     = min(i_el_nodes_tmp(2,ne_te), elems(neigh_ib_ie)%p%i_ver(ind2) )
!       i_node1_max = max(i_el_nodes_tmp(1,ne_te), elems(neigh_ib_ie)%p%i_ver(ind1) )
!       i_node2_max = max(i_el_nodes_tmp(2,ne_te), elems(neigh_ib_ie)%p%i_ver(ind2) )
!       if ( all( i_te_tmp(1,1:nn_te) .ne. i_node2 ) ) then
!         nn_te = nn_te + 1
!         i_te_tmp(1,nn_te) = i_node2
!         i_te_tmp(2,nn_te) = i_node2_max
!!        rr_te_tmp(:,nn_te) = 0.5d0 * ( & 
!!                      geo%points(:,i_el_nodes_tmp(2,ne_te)) +     &
!!                      geo%points(:,elems(neigh_ib_ie)%p%i_ver(ind2))     )
!         rr_te_tmp(:,nn_te) = 0.5d0 * ( & 
!                       geo%points(:,i_node2) + geo%points(:,i_node2_max) )
!         ii_te_tmp(1,ne_te) = nn_te
!       else 
!         do i1 = 1 , nn_te   ! find ...
!           if ( i_te_tmp(1,i1) .eq. i_node2 ) then 
!             ii_te_tmp(1,ne_te) = i1 ; exit
!           end if
!         end do
!       end if 
!       if ( all( i_te_tmp(1,1:nn_te) .ne. i_node1 ) ) then
!         nn_te = nn_te + 1
!         i_te_tmp(1,nn_te) = i_node1
!         i_te_tmp(2,nn_te) = i_node1_max
!!        rr_te_tmp(:,nn_te) = 0.5d0 * ( & 
!!                      geo%points(:,i_el_nodes_tmp(1,ne_te)) +     &
!!                      geo%points(:,elems(neigh_ib_ie)%p%i_ver(ind1))     )
!         rr_te_tmp(:,nn_te) = 0.5d0 * ( & 
!                       geo%points(:,i_node1) + geo%points(:,i_node1_max) )
!         ii_te_tmp(2,ne_te) = nn_te
!       else
!         do i1 = 1 , nn_te   ! find ...
!           if ( i_te_tmp(1,i1) .eq. i_node1 ) then 
!             ii_te_tmp(2,ne_te) = i1 ; exit
!           end if
!         end do
!       end if 
!      
!     end if
!    
!   end if
!     
!  end do
!
! end do
!
! ! From tmp to actual e_te array --------------------
! if (allocated(e_te)) deallocate(e_te)
! allocate(e_te(2,ne_te)) ; e_te = e_te_tmp(:,1:ne_te)
! 
! ! From tmp to actual e_te array --------------------
! if (allocated(i_te)) deallocate(i_te)
! allocate(i_te(2,nn_te)) ; i_te = i_te_tmp(:,1:nn_te)
! 
! ! From tmp to actual r_te array --------------------
! if (allocated(rr_te))  deallocate(rr_te)
! allocate(rr_te(3,nn_te)) ; rr_te = rr_te_tmp(:,1:nn_te)
! 
! ! From tmp to actual r_te array --------------------
! if (allocated(ii_te))  deallocate(ii_te)
! allocate(ii_te(2,ne_te)) ; ii_te = ii_te_tmp(:,1:ne_te)
!
! 
! write(*,*) ' n.edges at the te : ' , ne_te
! write(*,*) ' n.nodes at the te : ' , nn_te
!
! deallocate( e_te_tmp, i_te_tmp, rr_te_tmp, ii_te_tmp )
!
! ! Trailing edge elements connectivity --------------
! 
! allocate(neigh_te(2,ne_te) ) ; neigh_te = 0
! allocate(    o_te(2,ne_te) ) ;     o_te = 0
! 
! do e1 = 1 , ne_te
!  nSides1 = 2
!  do i1 = 1 , nSides1
!   do e2 = e1+1 , ne_te
!    nSides2 = 2
!    do i2 = 1 , nSides2
!     if ( ii_te(i1,e1) .eq. ii_te(i2,e2) ) then
!      neigh_te(i1,e1) = e2
!      neigh_te(i2,e2) = e1
!      if ( i1 .ne. i2 ) then
!       o_te(i1,e1) = 1
!       o_te(i2,e2) = 1
!      else
!       o_te(i1,e1) = -1
!       o_te(i2,e2) = -1
!      end if 
!     end if
!    end do    
!   end do
!  end do
! end do
!
! ! Compute unit vector at the te nodes ----------------
! !  for the first implicit wake panel  ----------------
! allocate(t_te  (3,nn_te)) ; t_te = 0.0_wp  
! allocate(ref_te(  nn_te)) ; ref_te = comp%ref_id 
! do i_n = 1 , nn_te
!   vec1 = 0.0_wp
!   do i_e = 1 , ne_te 
!     if ( any( i_n .eq. ii_te(:,i_e)) ) then
!       t_te(:,i_n) =  t_te(:,i_n) + elems(e_te(1,i_e))%p%nor &
!                                  + elems(e_te(2,i_e))%p%nor
!
!       ! TODO: refine the definition of the vec1
!!      vec1 = cross(nor(:,e_te(1,i_e)) , nor(:,e_te(2,i_e)) )
!       vec1 = vec1 + cross(elems(e_te(1,i_e))%p%nor , elems(e_te(2,i_e))%p%nor)
!
!     end if
!   end do
!
!   vec1 = vec1 - sum(vec1*v_rel) * v_rel
!   vec1 = vec1 / norm2(vec1)
!
!   ! projection ****
!!  write(*,*) ' vec1 ' , vec1
!   t_te(:,i_n) = t_te(:,i_n) - sum(t_te(:,i_n)*vec1) * vec1
!
!   ! normalisation
!   t_te(:,i_n) = t_te(:,i_n) / norm2(t_te(:,i_n))
!
!!  transpose transformation to obtain the "transformed" unit vecotr at the te
!   t_te(:,i_n) = matmul(transpose(geo%refs(comp%ref_id)%R_g),t_te(:,i_n))
!
! end do
!
!
! ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ! Update overall structure --------
! if (.not.allocated(te%e)) then ! it should be enough
!   allocate(te%e    (2,ne_te) ) ; te%e     =     e_te
!   allocate(te%i    (2,nn_te) ) ; te%i     =     i_te
!   allocate(te%rr   (3,nn_te) ) ; te%rr    =    rr_te
!   allocate(te%ii   (2,ne_te) ) ; te%ii    =    ii_te
!   allocate(te%neigh(2,ne_te) ) ; te%neigh = neigh_te
!   allocate(te%o    (2,ne_te) ) ; te%o     =     o_te
!   allocate(te%t    (2,nn_te) ) ; te%t     =     t_te
!   allocate(te%ref  (  nn_te) ) ; te%ref   =   ref_te
! else
!   nn_te_prev = size(te%i,2)
!   ne_te_prev = size(te%e,2)
!   allocate(e_te_tmp(2,size(te%e,2)+ne_te)) 
!   e_te_tmp(:,             1:size(te%e,2)    ) = te%e
!   e_te_tmp(:,size(te%e,2)+1:size(e_te_tmp,2)) = e_te
!   call move_alloc(e_te_tmp,te%e) 
!   allocate(i_te_tmp(2,size(te%i,2)+nn_te)) 
!   i_te_tmp(:,             1:size(te%i,2)    ) = te%i
!   i_te_tmp(:,size(te%i,2)+1:size(i_te_tmp,2)) = i_te
!   call move_alloc(i_te_tmp,te%i) 
!   allocate(rr_te_tmp(3,size(te%rr,2)+nn_te)) 
!   rr_te_tmp(:,              1:size(te%rr,2)    ) = te%rr
!   rr_te_tmp(:,size(te%rr,2)+1:size(rr_te_tmp,2)) = rr_te
!   call move_alloc(rr_te_tmp,te%rr) 
!   allocate(ii_te_tmp(2,size(te%ii,2)+ne_te)) 
!   ii_te_tmp(:,              1:size(te%ii,2)    ) = te%ii
!   ii_te_tmp(:,size(te%ii,2)+1:size(ii_te_tmp,2)) = ii_te + nn_te_prev 
!   call move_alloc(ii_te_tmp,te%ii) 
!   allocate(neigh_te_tmp(2,size(te%neigh,2)+ne_te)) 
!   neigh_te_tmp(:,                 1:size(te%neigh,2)    ) = te%neigh
!   where ( neigh_te .ne. 0 ) neigh_te = neigh_te + ne_te_prev
!   neigh_te_tmp(:,size(te%neigh,2)+1:size(neigh_te_tmp,2)) = neigh_te
!   call move_alloc(neigh_te_tmp,te%neigh) 
!   allocate( o_te_tmp(2,size(te%o ,2)+ne_te)) 
!   o_te_tmp(:,              1:size(te%o ,2)    ) = te%o 
!   o_te_tmp(:,size(te%o ,2)+1:size( o_te_tmp,2)) =  o_te 
!   call move_alloc(o_te_tmp,te%o ) 
!   allocate( t_te_tmp(3,size(te%t ,2)+nn_te)) 
!   t_te_tmp(:,              1:size(te%t ,2)    ) = te%t 
!   t_te_tmp(:,size(te%t ,2)+1:size( t_te_tmp,2)) =  t_te 
!   call move_alloc(t_te_tmp,te%t ) 
!   allocate(ref_te_tmp(size(te%ref   )+nn_te)) 
!   ref_te_tmp(                 1:size(te%ref   )   ) =  te%ref 
!   ref_te_tmp(  size(te%ref  )+1:size(ref_te_tmp  )) = ref_te 
!   call move_alloc(ref_te_tmp,te%ref) 
! end if
!
!
!end subroutine find_te
!
!! 04. Inside build_te : find_te_v --------
!! It is assumed that nodes and elements are sorted as follows
!!
!!  1-----5-----9-----13----17----21   <--- LE
!!  |  1  |  4  |  7  |  10 |  13 |
!!  2-----6-----10----14----18----22
!!  |  2  |  5  |  8  |  11 |  14 |
!!  3-----7-----11----15----19----23
!!  |  3  |  6  |  9  |  12 |  15 |
!!  4-----8-----12----16----20----24   <--- TE
!
!subroutine find_te_v( geo , comp , elems , te )
! type(t_geo)   , intent(in), target :: geo
! type(t_geo_component), intent(in) :: comp
! type(t_elem_p)       , intent(in) :: elems(:)
!!type(t_elem_p)       , intent(in) :: elems_m(:) ! read "actual" neighbors
! type(t_tedge)        , intent(inout) :: te
!
! ! tmp arrays --------
! integer , allocatable               :: e_te_tmp(:,:) , i_te_tmp(:,:) , ii_te_tmp(:,:) 
! integer , allocatable               :: neigh_te_tmp(:,:) , o_te_tmp(:,:)
! real(wp), allocatable               ::rr_te_tmp(:,:) , t_te_tmp(:,:)
! integer , allocatable               :: ref_te_tmp(:)
! integer , allocatable               :: i_el_nodes_tmp(:,:)
! ! actual arrays -----
! integer , allocatable               :: e_te(:,:) , i_te(:,:) , ii_te(:,:) 
! integer , allocatable               :: neigh_te(:,:) , o_te(:,:)
! real(wp), allocatable               ::rr_te(:,:) , t_te(:,:)
! integer , allocatable               :: ref_te(:)
!
! integer :: n_el
! integer :: ne_te , nn_te           ! # n. elements and nodes at TE (for the actual comp)
! integer :: ne_te_prev , nn_te_prev ! # n. elements and nodes at TE ( of the prev. comps) 
!
! integer :: io_te , io_tip
! integer :: i_el , n_free , n_free_le_te , i_v , ie_ind
! integer :: i_node1 , i_node2 , ind_te 
! integer :: nSides1 , nSides2 
! integer :: i_node3 , i_node4
!
! integer :: i1, i2, e1, e2
!
!
! n_el = size(comp%el)
!
! ! allocate tmp structures --
! allocate( e_te_tmp(2,(n_el+1)/2) )                ;  e_te_tmp = 0 
! allocate( i_el_nodes_tmp(2,size(comp%i_points)) ) ;  i_el_nodes_tmp = 0
! allocate( i_te_tmp      (2,size(comp%i_points)) ) ;  i_te_tmp = 0
! allocate(ii_te_tmp      (2,n_el               ) ) ; ii_te_tmp = 0
! allocate(rr_te_tmp      (3,size(comp%i_points)) ) ; rr_te_tmp = 0.0d0
! allocate( t_te_tmp      (3,size(comp%i_points)) ) ; rr_te_tmp = 0.0d0
!
!
! ! initialise counters ------
! ne_te = 0 ; nn_te = 0
!
! ! Find elements on the TE --------
! ! "logical" to know if you are on a tip element an on the TE/LE elements
! io_tip = 1  ;  io_te  = 0 ! initialisation 
! do i_el = 1 , n_el
!
!   ie_ind = comp%el(i_el)%id
!
!   n_free = 0 
!   do i_v = 1 , comp%el(i_el)%n_ver
!     if ( comp%el(i_el)%i_neigh(i_v) .eq. 0 ) n_free = n_free + 1
!   end do
!
!   ! Back on the tip ----
!   if ( n_free .eq. 2 ) then
!      n_free_le_te = 2 ! on the tip
!      if ( io_tip .eq. 0 ) io_tip = 1
!   end if
!   ! n. free edges on the te,le: tip -> 2 ; inner elements -> 1
!   if ( io_tip .eq. 1 ) then
!     n_free_le_te = 2
!   elseif ( io_tip .eq. 0 ) then
!     n_free_le_te = 1
!   end if 
!
!   if ( ( n_free .eq. n_free_le_te ) .and. ( io_te .eq. 1 ) ) then ! Element at TE
!
!     ! element number: first component only of e_te ----------   
!     ne_te = ne_te + 1
!     e_te_tmp(1,ne_te) = ie_ind 
!
!     ! find i_node1 , i_node2 --------------------------------
!     if     ( n_free .eq. 2 ) then     ! on wing tips --------
!       ! two sides have no neighbors: TE must be identified
!  
!       do i1 = 1 , comp%el(i_el)%n_ver
!         if ( ( comp%el(i_el)%i_neigh(i1)   .eq. 0 ) .and. &
!              ( comp%el(i_el-1)%i_neigh(i1) .ne. 0 ) ) then
!             ind_te = i1 ; exit 
!         end if
!       end do
!
!!      i_node1 = elems(ie_ind)%p%i_ver(ind_te)  
!!      i_node2 = elems(ie_ind)%p%i_ver( mod(ind_te,elems(ie_ind)%p%n_ver)+1 )   
!
!     elseif ( n_free .eq. 1 ) then     ! inner elems ---------
!
!       do i1 = 1 , comp%el(i_el)%n_ver
!         if ( comp%el(i_el)%i_neigh(i1) .eq. 0 ) then 
!           ind_te = i1 ; exit 
!         end if
!       end do
!
!     end if
!
!     i_node1 = elems(ie_ind)%p%i_ver(ind_te)  
!     i_node2 = elems(ie_ind)%p%i_ver( mod(ind_te  ,elems(ie_ind)%p%n_ver)+1 )   
!     i_node3 = elems(ie_ind)%p%i_ver( mod(ind_te+1,elems(ie_ind)%p%n_ver)+1 )   
!     i_node4 = elems(ie_ind)%p%i_ver( mod(ind_te+2,elems(ie_ind)%p%n_ver)+1 )   
!
!     ! build i_te , ii_te , rr_te ----------------------------
!     if ( all ( i_te_tmp(:,1:nn_te) .ne. i_node2 ) ) then ! new node
!       nn_te = nn_te + 1
!       i_te_tmp(:,nn_te) = i_node2
!       rr_te_tmp(:,nn_te) = geo%points(:,i_node2)
!       ii_te_tmp(1,ne_te) = nn_te
!       t_te_tmp(:,nn_te)  = geo%points(:,i_node2) - geo%points(:,i_node3)
!       t_te_tmp(:,nn_te)  = t_te_tmp(:,nn_te) / norm2(t_te_tmp(:,nn_te))
!
!!      transpose transformation to obtain the "transformed" unit vecotr at the te
!       t_te_tmp(:,nn_te) = matmul(transpose(geo%refs(comp%ref_id)%R_g), &
!                                                         t_te_tmp(:,nn_te))
!
!     else ! this node has been already found
!       do i1 = 1 , nn_te   ! find ...
!         if ( i_te_tmp(1,i1) .eq. i_node2 ) then 
!           ii_te_tmp(1,ne_te) = i1 ; exit
!         end if
!       end do
!     end if
!
!     if ( all ( i_te_tmp(:,1:nn_te) .ne. i_node1 ) ) then ! new node
!       nn_te = nn_te + 1
!       i_te_tmp(:,nn_te) = i_node1
!       rr_te_tmp(:,nn_te) = geo%points(:,i_node1)
!       ii_te_tmp(2,ne_te) = nn_te
!       t_te_tmp(:,nn_te)  = geo%points(:,i_node1) - geo%points(:,i_node4)
!       t_te_tmp(:,nn_te)  = t_te_tmp(:,nn_te) / norm2(t_te_tmp(:,nn_te))
!
!!      transpose transformation to obtain the "transformed" unit vecotr at the te
!       t_te_tmp(:,nn_te) = matmul(transpose(geo%refs(comp%ref_id)%R_g), &
!                                                         t_te_tmp(:,nn_te))
!
!     else ! this node has been already found
!       do i1 = 1 , nn_te   ! find ...
!         if ( i_te_tmp(1,i1) .eq. i_node1 ) then 
!           ii_te_tmp(2,ne_te) = i1 ; exit
!         end if
!       end do
!     end if
!
!
!   end if
!
!   ! Next element ----
!   if     ( ( n_free .eq. n_free_le_te ) .and. ( io_te .eq. 0 ) ) then
!      io_te = 1 
!   elseif ( ( n_free .eq. n_free_le_te ) .and. ( io_te .eq. 1 ) ) then
!      io_te = 0 ; io_tip = 0
!   end if
!
! end do
!
!  
! ! From tmp to actual e_te array --------------------
! if (allocated(e_te)) deallocate(e_te)
! allocate(e_te(2,ne_te)) ; e_te = e_te_tmp(:,1:ne_te)
! 
! ! From tmp to actual i_te array --------------------
! if (allocated(i_te)) deallocate(i_te)
! allocate(i_te(2,nn_te)) ; i_te = i_te_tmp(:,1:nn_te)
! 
! ! From tmp to actual rr_te array -------------------
! if (allocated(rr_te))  deallocate(rr_te)
! allocate(rr_te(3,nn_te)) ; rr_te = rr_te_tmp(:,1:nn_te)
! 
! ! From tmp to actual ii_te array -------------------
! if (allocated(ii_te))  deallocate(ii_te)
! allocate(ii_te(2,ne_te)) ; ii_te = ii_te_tmp(:,1:ne_te)
!
! ! From tmp to actual rr_te array -------------------
! if (allocated( t_te))  deallocate( t_te)
! allocate( t_te(3,nn_te)) ;  t_te =  t_te_tmp(:,1:nn_te)
!
! 
! write(*,*) ' n.edges at the te : ' , ne_te
! write(*,*) ' n.nodes at the te : ' , nn_te
!
! deallocate( e_te_tmp, i_te_tmp, rr_te_tmp, ii_te_tmp , t_te_tmp )
!
!
! ! Trailing edge elements connectivity --------------
! 
! allocate(neigh_te(2,ne_te) ) ; neigh_te = 0
! allocate(    o_te(2,ne_te) ) ;     o_te = 0
! 
! do e1 = 1 , ne_te
!  nSides1 = 2
!  do i1 = 1 , nSides1
!   do e2 = e1+1 , ne_te
!    nSides2 = 2
!    do i2 = 1 , nSides2
!     if ( ii_te(i1,e1) .eq. ii_te(i2,e2) ) then
!      neigh_te(i1,e1) = e2
!      neigh_te(i2,e2) = e1
!      if ( i1 .ne. i2 ) then
!       o_te(i1,e1) = 1
!       o_te(i2,e2) = 1
!      else
!       o_te(i1,e1) = -1
!       o_te(i2,e2) = -1
!      end if 
!     end if
!    end do    
!   end do
!  end do
! end do
!
! ! Compute unit vector at the te nodes ----------------
! !  for the first implicit wake panel  ----------------
! allocate(ref_te(  nn_te)) ; ref_te = comp%ref_id 
! 
!
! ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ! Update overall structures and TE connectivity as in find_te()
! if (.not.allocated(te%e)) then ! it should be enough
!   allocate(te%e    (2,ne_te) ) ; te%e     =     e_te
!   allocate(te%i    (2,nn_te) ) ; te%i     =     i_te
!   allocate(te%rr   (3,nn_te) ) ; te%rr    =    rr_te
!   allocate(te%ii   (2,ne_te) ) ; te%ii    =    ii_te
!   allocate(te%neigh(2,ne_te) ) ; te%neigh = neigh_te
!   allocate(te%o    (2,ne_te) ) ; te%o     =     o_te
!   allocate(te%t    (2,nn_te) ) ; te%t     =     t_te
!   allocate(te%ref  (  nn_te) ) ; te%ref   =   ref_te
! else
!   nn_te_prev = size(te%i,2)
!   ne_te_prev = size(te%e,2)
!   allocate(e_te_tmp(2,size(te%e,2)+ne_te)) 
!   e_te_tmp(:,             1:size(te%e,2)    ) = te%e
!   e_te_tmp(:,size(te%e,2)+1:size(e_te_tmp,2)) = e_te
!   call move_alloc(e_te_tmp,te%e) 
!   allocate(i_te_tmp(2,size(te%i,2)+nn_te)) 
!   i_te_tmp(:,             1:size(te%i,2)    ) = te%i
!   i_te_tmp(:,size(te%i,2)+1:size(i_te_tmp,2)) = i_te
!   call move_alloc(i_te_tmp,te%i) 
!   allocate(rr_te_tmp(3,size(te%rr,2)+nn_te)) 
!   rr_te_tmp(:,              1:size(te%rr,2)    ) = te%rr
!   rr_te_tmp(:,size(te%rr,2)+1:size(rr_te_tmp,2)) = rr_te
!   call move_alloc(rr_te_tmp,te%rr) 
!   allocate(ii_te_tmp(2,size(te%ii,2)+ne_te)) 
!   ii_te_tmp(:,              1:size(te%ii,2)    ) = te%ii
!   ii_te_tmp(:,size(te%ii,2)+1:size(ii_te_tmp,2)) = ii_te + nn_te_prev 
!   call move_alloc(ii_te_tmp,te%ii) 
!   allocate(neigh_te_tmp(2,size(te%neigh,2)+ne_te)) 
!   neigh_te_tmp(:,                 1:size(te%neigh,2)    ) = te%neigh
!   where ( neigh_te .ne. 0 ) neigh_te = neigh_te + ne_te_prev
!   neigh_te_tmp(:,size(te%neigh,2)+1:size(neigh_te_tmp,2)) = neigh_te
!   call move_alloc(neigh_te_tmp,te%neigh) 
!   allocate( o_te_tmp(2,size(te%o ,2)+ne_te)) 
!   o_te_tmp(:,              1:size(te%o ,2)    ) = te%o 
!   o_te_tmp(:,size(te%o ,2)+1:size( o_te_tmp,2)) =  o_te 
!   call move_alloc(o_te_tmp,te%o ) 
!   allocate( t_te_tmp(3,size(te%t ,2)+nn_te)) 
!   t_te_tmp(:,              1:size(te%t ,2)    ) = te%t 
!   t_te_tmp(:,size(te%t ,2)+1:size( t_te_tmp,2)) =  t_te 
!   call move_alloc(t_te_tmp,te%t ) 
!   allocate(ref_te_tmp(size(te%ref   )+nn_te)) 
!   ref_te_tmp(                 1:size(te%ref   )   ) =  te%ref 
!   ref_te_tmp(  size(te%ref  )+1:size(ref_te_tmp  )) = ref_te 
!   call move_alloc(ref_te_tmp,te%ref) 
! end if
!
!
!end subroutine find_te_v
!
!
!end subroutine build_te

!----------------------------------------------------------------------

!> Compute the coefficients (pot_vel_stencil, for surfpan elements) for
!!  computing the velocity from the velocity potential (phi = -mu).
!! On-body analysis for 3dPanels (surfpan). This coefficients are constant
!!  in the local frame, associated with the component. In order to obtain
!!  the components of the velcoity in the base frame, the global rotation
!!  matrix is needed. 
subroutine create_local_velocity_stencil (geo, elems)
 type(t_geo), intent(inout) :: geo
 type(t_elem_p), intent(in)  :: elems(:)

 real(wp) :: surf_bubble
 integer  :: i_comp , i_el , i_v

 do i_comp = 1 , size(geo%components)

  ! Field pot_vel_stencil belongs to surfpan elements only!
  if ( geo%components(i_comp)%comp_el_type(1:1) .eq. 'p' ) then

   do i_el = 1 , size(geo%components(i_comp)%el) 

    if ( allocated(geo%components(i_comp)%el(i_el)%pot_vel_stencil) ) then
      write(*,*) ' WARNING. Already allocated pot_vel_stencil array for '
      write(*,*) ' component , element ' , i_comp , i_el
      deallocate(geo%components(i_comp)%el(i_el)%pot_vel_stencil)
    end if
    allocate(geo%components(i_comp)%el(i_el)%pot_vel_stencil & 
             (3,geo%components(i_comp)%el(i_el)%n_ver) ) 

    surf_bubble = geo%components(i_comp)%el(i_el)%area

    do i_v = 1 , geo%components(i_comp)%el(i_el)%n_ver

      ! Update surf_bubble
      surf_bubble = surf_bubble + & 
           elems(geo%components(i_comp)%el(i_el)%id)%p%area / &
           dble(elems(geo%components(i_comp)%el(i_el)%id)%p%n_ver)

      geo%components(i_comp)%el(i_el)%pot_vel_stencil(:,i_v) = &
               cross( geo%components(i_comp)%el(i_el)%edge_vec(:,i_v) , &
                      geo%components(i_comp)%el(i_el)%nor )  


    end do

    geo%components(i_comp)%el(i_el)%pot_vel_stencil = &
      geo%components(i_comp)%el(i_el)%pot_vel_stencil / surf_bubble

   end do 

! else if ( geo%components(i_comp)%comp_el_type(1:1) .eq. 'v' ) then
!   no stencil for the velocity ...

  end if

 end do
 

end subroutine create_local_velocity_stencil

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
  if ( geo%components(i_comp)%comp_el_type(1:1) .eq. 'v' ) then

    n_el = size(comp%el)

    ! Some checks added: the first element must be at the LE, the last at the TE
    if ( associated(comp%el(  1 )%neigh(1)%p)) then
      call error(this_sub_name, this_mod_name, 'First element must be at the LE')
    end if
    if ( associated(comp%el(n_el)%neigh(3)%p) ) then
      call error(this_sub_name, this_mod_name, 'Last element must be at the TE')
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
    
     if ( abs( h2 - comp%el(i_el)%dy ) .gt. 1e-4 ) then
       write(*,*) ' WARNING: non parallel stripe edges (rough check) '
     end if
    
    end do

    ! allocate and fill comp%strip_elem array
    if ( mod(n_el,n_s) .ne. 0 ) then
      call error(this_sub_name, this_mod_name, ' mod(n_elem,n_s) .ne. 0. &
         & Wrong input file for vortring component.')
    end if
    n_c = n_el / n_s  ! integer division to find number of chord panels

    do i_s = 1 , n_s
      do i_c = 1 , n_c
!       write(*,*) ' comp%el(',i_c+(i_s-1)*n_c,')%stripe_elem : '
        allocate( comp%el(i_c+(i_s-1)*n_c)%stripe_elem(i_c) )
        do i_c2 = 1 , i_c
          comp%el(i_c+(i_s-1)*n_c)%stripe_elem(i_c2)%p => &
                                        comp%el(i_c2+(i_s-1)*n_c)
!         write(*,*) comp%el(i_c+(i_s-1)*n_c)%stripe_elem(i_c2)
        end do
      end do
    end do

! elseif ( el_type .eq. 'v' ) do_nothing
 
  end if 

  end associate
 end do


end subroutine create_strip_connectivity


!! OLD node-element connectivity ( ee array )
!! ! It is assumed that nodes and elements are sorted as follows
!! !
!! !  1-----5-----9-----13----17----21   <--- LE
!! !  |  1  |  4  |  7  |  10 |  13 |
!! !  2-----6-----10----14----18----22
!! !  |  2  |  5  |  8  |  11 |  14 |
!! !  3-----7-----11----15----19----23
!! !  |  3  |  6  |  9  |  12 |  15 |
!! !  4-----8-----12----16----20----24   <--- TE
!! ! and the ee array is built is  1, 2, 6, 5
!! !                               2, 3, 7, 6
!! !                               3, 4, 8, 7
!! !                               5, 6,10, 9
!! !                              ...
!! !                              19,20,24,23
!! subroutine create_strip_connectivity(geo)
!!  type(t_geo),  intent(inout) :: geo
!! 
!!  real(wp) :: h2
!! 
!!  integer :: io_te , io_tip
!!  integer :: n_el , ie_ind
!!  integer :: i_comp , i_el 
!! 
!! 
!!  do i_comp = 1 , size(geo%components)
!!   associate(comp => geo%components(i_comp))
!! 
!!   ! connectivity built for lifting surface only
!!   if ( geo%components(i_comp)%comp_el_type(1:1) .eq. 'v' ) then
!! 
!!     n_el = size(comp%el)
!!     
!!     io_tip = 1 ; io_te = 0 ! initialisation
!!     do i_el = 1 , n_el     
!!     
!!      ie_ind = comp%el(i_el)%id
!!     
!!      if ( comp%el(i_el)%i_neigh(4) .eq. 0 ) then
!!       comp%el(i_el)%stripe_1 = 0
!!      else 
!!       comp%el(i_el)%stripe_1 = comp%el(i_el)%i_neigh(4)
!!      end if
!!     
!!      comp%el(i_el)%dy = sum( cross( comp%el(i_el)%edge_uni(:,1) ,    &
!!                                     comp%el(i_el)%edge_vec(:,2)  ) * &
!!                              comp%el(i_el)%nor )
!!     
!!      h2 = sum( cross( comp%el(i_el)%edge_uni(:,1) ,                        &
!!                       comp%el(i_el)%ver(:,4) - comp%el(i_el)%ver(:,2) ) * &  
!!                comp%el(i_el)%nor )
!!     
!!      if ( abs( h2 - comp%el(i_el)%dy ) .gt. 1e-4 ) then
!!        write(*,*) ' WARNING: non parallel stripe edges (rough check) '
!!      end if
!!     
!!     end do
!!  
!!   end if 
!! 
!!   write(*,*) ' In create_strip_conectivity. comp.' , i_comp  
!!   end associate
!!  end do
!! 
!! 
!! end subroutine create_strip_connectivity

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

subroutine update_geometry(geo, t, update_static)
 type(t_geo), intent(inout) :: geo
 real(wp), intent(in) :: t
 logical, intent(in) :: update_static

 integer :: i_comp, ie
 !type(t_elem_p) :: el
 
 !update all the references
 call update_all_references(geo%refs,t)

 do i_comp = 1,size(geo%components)
  associate(comp => geo%components(i_comp))

    if (comp%moving .or. update_static) then
      geo%points(:,comp%i_points) = move_points(comp%loc_points, &
                           geo%refs(comp%ref_id)%R_g, &
                           geo%refs(comp%ref_id)%of_g)
      select case(trim(comp%comp_el_type))
       case('p','v')
        do ie = 1,size(comp%el)
          call calc_geo_data_pan(comp%el(ie),geo%points(:,comp%el(ie)%i_ver))

          !Calculate the velocity of the centers to impose the boundary condition
          call calc_geo_vel(comp%el(ie), geo%refs(comp%ref_id)%G_g, &
                                  geo%refs(comp%ref_id)%f_g)
        enddo
       case('l')
        do ie = 1,size(comp%el)
          call calc_geo_data_ll(comp%el(ie),geo%points(:,comp%el(ie)%i_ver))

          !Calculate the velocity of the centers to impose the boundary condition
          call calc_geo_vel(comp%el(ie), geo%refs(comp%ref_id)%G_g, &
                                  geo%refs(comp%ref_id)%f_g)
        enddo
       case('a')
        do ie = 1,size(comp%el)
          call calc_geo_data_ad(comp%el(ie),geo%points(:,comp%el(ie)%i_ver))

          !Calculate the velocity of the centers to impose the boundary condition
          call calc_geo_vel(comp%el(ie), geo%refs(comp%ref_id)%G_g, &
                                  geo%refs(comp%ref_id)%f_g)
        enddo
      end select
    end if

  end associate
 enddo


end subroutine update_geometry

!----------------------------------------------------------------------

!> Destroy all the contents of a geometry type
!!
!! TODO: check that passing it with intent out is enough to deallocate it
subroutine destroy_geometry(geo, elems)
 type(t_geo), intent(out) :: geo
 type(t_elem_p), allocatable, intent(out) :: elems(:)
 
 integer :: i
 !CHECK: passing the elements with intent out should be enough to deallocate
 !everything, but should be checked
 !deallocate(elems)
 call destroy_references(geo%refs)

 !dummy to avoid warnings
 i = size(elems)
 

end subroutine destroy_geometry

!----------------------------------------------------------------------




!----------------------------------------------------------------------

end module mod_geometry
