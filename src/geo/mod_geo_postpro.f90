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


!> This is a module dedicated to the handling of the geometry during the 
!! postprocessing
!!
!! It uses the same data structures of the postprocessing as well as
!! performing similar tasks, but much less cumbersome 

module mod_geo_postpro

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

use mod_geometry, only: &
  t_geo, t_geo_component , calc_geo_data_pan , calc_geo_vel

use mod_wake_pan, only: &
  t_wake_panels

use mod_wake_ring, only: &
  t_wake_rings

use mod_stringtools, only: &
  LowCase, IsInList, strip_mult_appendix

!----------------------------------------------------------------------

implicit none

public :: load_components_postpro, update_points_postpro , prepare_geometry_postpro, expand_actdisk_postpro, &
          prepare_wake_postpro

private

!----------------------------------------------------------------------

character(len=*), parameter :: this_mod_name = 'mod_geo_postpro'

!----------------------------------------------------------------------

contains

!----------------------------------------------------------------------

!----------------------------------------------------------------------

!!!> Create all the goemetry
!!!! 
!!!! The geometry creation proceeds in this way:
!!!! -# The main geometry file is read, and the components therein contained are
!!!!    created (reading the mesh)
!!!! -# If some additional geometry files are present, each one is read and all
!!!!    its component are created
!!!! -# The main geometry arrays, SurfPan and VortRin are allocated and filled
!!!!    in the prepare_geometry routine
!!!! -# The element pointer array used to build/solve the linear system is
!!!!    created, pointed at each element and then re-ordered with the static
!!!!    elements first, and the dynamic elements at the end
!!subroutine create_geometry(geo_file_name, ref_file_name, in_file_name,  geo, &
!!                           te, elems, elems_ll, elems_tot, airfoil_data, &
!!                           sim_param)
!! character(len=*), intent(in) :: geo_file_name
!! character(len=*), intent(inout) :: ref_file_name
!! character(len=*), intent(in) :: in_file_name
!! type(t_geo), intent(out), target :: geo
!! type(t_elem_p), allocatable, intent(out) :: elems(:)
!! type(t_elem_p), allocatable, intent(out) :: elems_ll(:)
!! type(t_elem_p), allocatable, intent(out) :: elems_tot(:)
!! type(t_tedge), intent(out) :: te
!! type(t_aero_tab) , allocatable, intent(out) :: airfoil_data(:)
!! type(t_sim_param) , intent(inout) :: sim_param
!! real(wp)                     :: tstart
!!
!! !character(len=max_char_len) :: reference_file
!! !character(len=max_char_len) :: geo_file_name
!!
!! integer :: i, j, is, im,  i_comp, i_ll, i_tot
!! type(t_elem_p), allocatable :: temp_static(:), temp_moving(:)
!!
!! integer , allocatable :: el_id_old(:), el_id_old_static(:), el_id_old_moving(:)
!!
!! character(len=max_char_len) :: msg
!!
!!  tstart = sim_param%t0
!! 
!!  !build the reference frames
!!  !read which is the reference frame file, if default the main input file is 
!!  !employed
!!  !reference_file = getstr(prms, 'ReferenceFile')
!!  if (trim(ref_file_name) .eq. 'no_set') then
!!    ref_file_name = trim(in_file_name)
!!  endif
!!
!!  call build_references(geo%refs, trim(ref_file_name), sim_param)
!!
!!  !Load the components from the file created by the preprocessor
!!  !geo_file_name = getstr(prms, 'GeometryFile')
!!  call check_preproc(geo_file_name)
!!  call load_components(geo, trim(geo_file_name), sim_param, te)
!!
!!  call import_aero_tab(geo,airfoil_data)
!!
!!  ! Initialisation
!!  geo%nelem      = 0
!!  geo%nll        = 0
!!  geo%nstatic    = 0
!!  geo%nmoving    = 0
!!  geo%nstatic_ll = 0
!!  geo%nmoving_ll = 0
!!  geo%nSurfPan   = 0
!!  geo%nVortRin   = 0
!!  geo%nLiftLin   = 0
!!
!!  ! count the elements
!!  do i_comp = 1,size(geo%components)
!!
!!    if (trim(geo%components(i_comp)%comp_el_type) .eq. 'p') then
!!      if(geo%components(i_comp)%moving) then
!!        geo%nmoving = geo%nmoving + geo%components(i_comp)%nelems
!!      else
!!        geo%nstatic = geo%nstatic + geo%components(i_comp)%nelems
!!      endif
!!      geo%nSurfPan = geo%nSurfPan + geo%components(i_comp)%nelems
!!      geo%nelem = geo%nelem + geo%components(i_comp)%nelems
!!
!!    elseif (trim(geo%components(i_comp)%comp_el_type) .eq. 'v') then
!!      if(geo%components(i_comp)%moving) then
!!        geo%nmoving = geo%nmoving + geo%components(i_comp)%nelems
!!      else
!!        geo%nstatic = geo%nstatic + geo%components(i_comp)%nelems
!!      endif
!!      geo%nVortRin = geo%nVortRin + geo%components(i_comp)%nelems
!!      geo%nelem = geo%nelem + geo%components(i_comp)%nelems
!!
!!    elseif (trim(geo%components(i_comp)%comp_el_type) .eq. 'l') then
!!      if(geo%components(i_comp)%moving) then
!!        geo%nmoving_ll = geo%nmoving_ll + geo%components(i_comp)%nelems
!!      else
!!        geo%nstatic_ll = geo%nstatic_ll + geo%components(i_comp)%nelems
!!      endif
!!      geo%nLiftLin = geo%nLiftLin + geo%components(i_comp)%nelems
!!      geo%nll = geo%nll + geo%components(i_comp)%nelems
!!    endif
!!    
!!  enddo
!!
!!  ! calculate the geometric quantities
!!  ! already update the geometry for the first time to get the right 
!!  ! starting geometrical condition
!!  call prepare_geometry(geo)
!!  call update_geometry(geo, tstart, update_static=.true.)
!!
!!  if(sim_param%debug_level .ge. 3) then
!!    call printout(nl//' Geometry details:' ) 
!!    write(msg,'(A,I9)') '  number of elements:        ' ,geo%nelem
!!    call printout(msg)
!!    write(msg,'(A,I9)') '  number of static elements: ' ,geo%nstatic
!!    call printout(msg)
!!    write(msg,'(A,I9)') '  number of moving elements: ' ,geo%nmoving
!!    call printout(msg)
!!    write(msg,'(A,I9)') '  number of surface panels:  ' ,geo%nsurfpan
!!    call printout(msg)
!!    write(msg,'(A,I9)') '  number of vortex rings:    ' ,geo%nvortrin
!!    call printout(msg)
!!  endif
!!
!!  !Create the vector of pointers to all the elements
!!  allocate(elems(geo%nelem), elems_ll(geo%nll), elems_tot(geo%nelem+geo%nll)) 
!!  i=0; i_ll=0; i_tot=0
!!  do i_comp = 1,size(geo%components)
!!    
!!    if (trim(geo%components(i_comp)%comp_el_type) .eq. 'p' .or. &
!!        trim(geo%components(i_comp)%comp_el_type) .eq. 'v') then
!!
!!      do j = 1,size(geo%components(i_comp)%el)
!!        i = i+1
!!        i_tot = i_tot+1
!!        elems(i)%p => geo%components(i_comp)%el(j)
!!  !      elems_tot(i_tot)%p => geo%components(i_comp)%el(j)
!!      enddo
!!
!!    elseif (trim(geo%components(i_comp)%comp_el_type) .eq. 'l') then
!!
!!      do j = 1,size(geo%components(i_comp)%el)
!!        i_ll = i_ll+1
!!        i_tot = i_tot+1
!!        elems_ll(i_ll)%p => geo%components(i_comp)%el(j)
!!  !      elems_tot(i_tot)%p => geo%components(i_comp)%el(j)
!!      enddo
!!
!!    endif
!!  enddo
!!
!!  ! Sort elements: first static, then moving ------- 
!!  !fill in the two temporaries
!!  allocate(temp_static(geo%nstatic), temp_moving(geo%nmoving))
!!  allocate(el_id_old(geo%nelem))          ; el_id_old = 0
!!  allocate(el_id_old_static(geo%nstatic)) ; el_id_old_static = 0
!!  allocate(el_id_old_moving(geo%nmoving)) ; el_id_old_moving = 0
!!  is = 0; im = 0;
!!  do i = 1,geo%nelem
!!    if(elems(i)%p%moving) then
!!      im = im+1
!!      temp_moving(im) = elems(i)
!!      el_id_old_moving(im) = i
!!    else
!!      is = is+1
!!      temp_static(is) = elems(i)
!!      el_id_old_static(is) = i
!!    endif
!!  enddo
!!
!!  !Now might be more bombproof to deallocate and allocate, but for the moment..
!!  elems(1:geo%nstatic) = temp_static
!!  elems(geo%nstatic+1:geo%nelem) = temp_moving
!!
!!  el_id_old(1:geo%nstatic) = el_id_old_static
!!  el_id_old(geo%nstatic+1:geo%nelem) = el_id_old_moving
!! 
!!  !Update the indexing since we re-ordered the vector
!!  do i = 1,geo%nelem
!!    elems(i)%p%id = i
!!  end do
!!  !Update elem-elem connectivity after re-ordering, for the moment only for
!!  !implicit panel elements
!!  !do i = 1,geo%nelem
!!  !  do j = 1,elems(i)%p%n_ver
!!  !    if ( elems(i)%p%i_neigh(j) .ne. 0 ) then
!!  !      elems(i)%p%i_neigh(j) = el_id_old( elems(i)%p%i_neigh(j) )
!!  !    else
!!  !      elems(i)%p%i_neigh(j) = 0 
!!  !    end if     
!!  !  end do 
!!  !end do 
!!  !Update te structure
!!  !do i = 1,size(te%e,2)
!!  !  do j = 1,2
!!  !    if ( te%e(j,i) .ne. 0 ) then
!!  !      te%e(j,i) = el_id_old( te%e(j,i) )
!!  !    else
!!  !      te%e(j,i) = 0 
!!  !    end if     
!!  !  end do 
!!  !end do
!!
!!  !Now re-order the lifting line elements
!!  deallocate(temp_static, temp_moving)
!!  allocate(temp_static(geo%nstatic_ll), temp_moving(geo%nmoving_ll))
!!  is = 0; im = 0;
!!  do i = 1,geo%nll
!!    if(elems_ll(i)%p%moving) then
!!      im = im+1
!!      temp_moving(im) = elems_ll(i)
!!    else
!!      is = is+1
!!      temp_static(is) = elems_ll(i)
!!    endif
!!  enddo
!!  if(geo%nll .gt. 0) then
!!    elems_ll(1:geo%nstatic_ll) = temp_static
!!    elems_ll(geo%nstatic_ll+1:geo%nll) = temp_moving
!!  endif
!!
!!  !Update te%neigh NOT NEEDED, because in te numbering
!!
!!  deallocate(temp_static, temp_moving)
!!  deallocate(el_id_old, el_id_old_static, el_id_old_moving)
!!
!!  !Patch together everything in elems_tot
!!  elems_tot(1:geo%nelem) = elems
!!  elems_tot(geo%nelem+1:geo%nelem+geo%nll) =elems_ll
!!
!!  call create_local_velocity_stencil(geo,elems)    ! for surfpan only (3dP)
!!
!!  call create_strip_connectivity(geo)
!!
!!end subroutine create_geometry

!----------------------------------------------------------------------

subroutine load_components_postpro(comps, points, nelem, floc, &
             components_names, all_comp)
 type(t_geo_component), allocatable, intent(inout) :: comps(:)
 real(wp), allocatable, intent(out) :: points(:,:)
 integer, intent(out)               :: nelem
 !integer, allocatable, intent(out)  :: ep_conn(:,:)
 !character(len=*), intent(in) :: in_file
 integer(h5loc), intent(in) :: floc
 character(len=*), allocatable, intent(inout) :: components_names(:)
 logical, intent(in) :: all_comp

 type(t_geo_component), allocatable :: comp_temp(:)
 integer :: i1 , i2, i3
 integer, allocatable :: ee(:,:)
 real(wp), allocatable :: rr(:,:)
 character(len=max_char_len) :: comp_el_type, comp_name, comp_name_stripped
 integer :: points_offset, n_vert! , elems_offset
 real(wp), allocatable :: points_tmp(:,:)
 character(len=max_char_len) :: ref_tag, ref_tag_m
 integer :: ref_id
 character(len=max_char_len) :: cname !, msg
 integer(h5loc) :: gloc, cloc , geo_loc
 integer :: n_comp, i_comp, n_comp_tot, i_comp_tot , i_comp_tmp

 character(len=max_char_len), allocatable :: components(:) , components_tmp(:)
 character(len=max_char_len) :: component , component_stripped

 !integer :: ie_t, ie
 !integer :: n_mult, i_mult
 !logical :: mult

 ! Connectivity and te structures 
 !integer , allocatable :: neigh(:,:)
 ! Lifting Line elements
 !real(wp), allocatable :: normalised_coord_e(:,:)
 !integer                 , allocatable :: i_airfoil_e(:,:)
 !character(max_char_len) , allocatable :: airfoil_list(:)
 !integer                 , allocatable :: nelem_span_list(:)

 ! trailing edge ------
 !integer , allocatable :: e_te(:,:) , i_te(:,:) , ii_te(:,:)
 !integer , allocatable :: neigh_te(:,:) , o_te(:,:)
 !real(wp), allocatable :: rr_te(:,:) , t_te(:,:)
 !integer :: ne_te , nn_te
 ! tmp arrays --------
 !type(t_elem_p) , allocatable :: e_te_tmp(:,:)
 !integer, allocatable  ::i_te_tmp(:,:) , ii_te_tmp(:,:) 
 !integer , allocatable :: neigh_te_tmp(:,:) , o_te_tmp(:,:)
 !real(wp), allocatable ::rr_te_tmp(:,:) , t_te_tmp(:,:)
 !integer , allocatable :: ref_te_tmp(:)
 !integer :: ne_te_prev , nn_te_prev ! # n. elements and nodes at TE ( of the prev. comps) 

 character(len=*), parameter :: this_sub_name = 'load_components_postpro'

  ! Read all the components
  call open_hdf5_group(floc,'Components',gloc)
  call read_hdf5(n_comp_tot,'NComponents',gloc)

! !DEBUG
! write(*,*) nl//' All the components: '
  allocate(components(n_comp_tot))
  do i_comp_tot = 1 , n_comp_tot
    write(cname,'(A,I3.3)') 'Comp',i_comp_tot
    call open_hdf5_group(gloc,trim(cname),cloc)
    call read_hdf5(components(i_comp_tot),'CompName',cloc)
    call close_hdf5_group(cloc)
!   !DEBUG
!   write(*,*) trim(components(i_comp_tot))
  end do
  
! if ( allocated(components_names) ) then
!   write(*,*) ' components_names : '
!   do i_comp_tot = 1 , size(components_names)
!     write(*,*) trim(components_names(i_comp_tot))
!   end do
! end if

! RE-BUILD components_names() to host multiple components ++++++++++++++++++
! components_names: input for post-processing
! comp_name: all the components of the model (blades: Hub01__01, ..., Hub01__Nb)
! comp_name_stripped: Hub01__01, ..., Hub01__Nb  ---> Hub01
! multiple components: if components_names(i1) is <Hub>, then add all the blades to the output
 
  allocate(components_tmp(n_comp_tot))
  i_comp_tmp = 0 
  if ( allocated(components_names) ) then
    do i1 = 1 , size(components_names)

      do i2 = 1 , n_comp_tot
        
        call strip_mult_appendix(components(i2), component_stripped, '__') 
        
        if ( trim(components_names(i1)) .eq. trim(component_stripped) ) then
          i_comp_tmp = i_comp_tmp + 1
          components_tmp(i_comp_tmp) = trim(components(i2))
          write(*,*) ' a. +1 '
        end if

        if ( trim(components_names(i1)) .eq. trim(components(i2)) ) then
          i_comp_tmp = i_comp_tmp + 1
          components_tmp(i_comp_tmp) = trim(components(i2))
          write(*,*) ' b. +1 '
        end if

      end do
 
    end do

    deallocate(components_names)
    write(*,*) ' i_comp_tmp : ' , i_comp_tmp
    allocate(components_names(i_comp_tmp)) 
    components_names = components_tmp(1:i_comp_tmp)
 
  ! If components_names was not allocated, all the components are required
  else ! ( .not. allocated(components_names) ) then

    allocate(components_names(n_comp_tot))
    do i_comp_tot = 1 , n_comp_tot
      write(cname,'(A,I3.3)') 'Comp',i_comp_tot
      call open_hdf5_group(gloc,trim(cname),cloc)
      call read_hdf5(components_names(i_comp_tot),'CompName',cloc)
      call close_hdf5_group(cloc)
    end do

  end if

! !DEBUG
! i_comp_tot = size(components_names)
! write(*,*) ' i_comp_tot : ' , i_comp_tot
! write(*,*) ' In load_components_postro. Components_names: ' 
! do i_comp_tot = 1 , size(components_names)
!   write(*,*) i_comp_tot , ' : ' , trim(components_names(i_comp_tot))
! end do

!  allocate(comps(n_comp))
  nelem = 0
  i_comp = 0; n_comp = 0 
  do i_comp_tot = 1,n_comp_tot
    
    write(cname,'(A,I3.3)') 'Comp',i_comp_tot
    call open_hdf5_group(gloc,trim(cname),cloc)
    
    call read_hdf5(comp_name,'CompName',cloc)

    !Strip the appendix of multiple components, to load all the multiple
    !components at once
!   write(*,*) ' comp_name : ' , trim(comp_name)
    call strip_mult_appendix(comp_name, comp_name_stripped, '__') 
!   write(*,*) ' comp_name_stripped : ' , trim(comp_name_stripped)
    
!   if(IsInList(comp_name_stripped, components_names) .or. all_comp) then
!   write(*,*) ' all_comp  :' , all_comp
!   write(*,*) ' comp_name :' , trim(comp_name)
    if(IsInList(comp_name, components_names) .or. all_comp) then

      i_comp = i_comp+1; n_comp = n_comp+1
      allocate(comp_temp(n_comp))
      if(allocated(comps)) comp_temp(1:n_comp-1) = comps 
      call move_alloc(comp_temp, comps)

      comps(i_comp)%comp_id = i_comp_tot
      call read_hdf5(ref_id,'RefId',cloc)
      call read_hdf5(ref_tag,'RefTag',cloc)

      comps(i_comp)%ref_id  = ref_id
      comps(i_comp)%ref_tag = trim(ref_tag)
      !geo%components(i_comp)%moving  = geo%refs(ref_id)%moving

      ! ====== READING =====
      call read_hdf5(comp_el_type,'ElType',cloc)
      comps(i_comp)%comp_el_type = trim(comp_el_type)

      comps(i_comp)%comp_name = trim(comp_name)

      ! Geometry --------------------------
      call open_hdf5_group(cloc,'Geometry',geo_loc)
      call read_hdf5_al(ee   ,'ee'   ,geo_loc)
      call read_hdf5_al(rr   ,'rr'   ,geo_loc)
      !!call read_hdf5_al(neigh,'neigh',geo_loc)
      !! !element-specific reads 
      !!if ( comp_el_type(1:1) .eq. 'l' ) then
      !!  call read_hdf5_al(airfoil_list      ,'airfoil_list'      ,geo_loc) ! (:)
      !!  call read_hdf5_al(nelem_span_list   ,'nelem_span_list'   ,geo_loc) ! (:)
      !!  call read_hdf5_al(i_airfoil_e       ,'i_airfoil_e'       ,geo_loc) ! (:,:)
      !!  call read_hdf5_al(normalised_coord_e,'normalised_coord_e',geo_loc) ! (:,:)
      !!  allocate(geo%components(i_comp)%airfoil_list(size(airfoil_list))) 
      !!  geo%components(i_comp)%airfoil_list = airfoil_list 
      !!  allocate(geo%components(i_comp)%nelem_span_list(size(nelem_span_list))) 
      !!  geo%components(i_comp)%nelem_span_list = nelem_span_list 
      !!  allocate(geo%components(i_comp)%i_airfoil_e( &
      !!        size(i_airfoil_e,1),size(i_airfoil_e,2)) ) 
      !!  geo%components(i_comp)%i_airfoil_e = i_airfoil_e 
      !!  allocate(geo%components(i_comp)%normalised_coord_e( &
      !!        size(normalised_coord_e,1),size(normalised_coord_e,2))) 
      !!  geo%components(i_comp)%normalised_coord_e = normalised_coord_e
      !!end if
      call close_hdf5_group(geo_loc)

      ! ======= CREATING ELEMENTS ======

      ! --- treat the points ---
      if(allocated(points)) then
        points_offset = size(points,2) 
      else
        points_offset = 0
      endif

      !store the read points into the local points
      allocate(comps(i_comp)%loc_points(3,size(rr,2)))
      comps(i_comp)%loc_points = rr
      
      !Now for the moments the points are stored here without moving them, 
      !will be moved later, consider not storing them here at all
      allocate(points_tmp(3,size(rr,2)+points_offset))
      if (points_offset .gt. 0) points_tmp(:,1:points_offset) = points
      points_tmp(:,points_offset+1:points_offset+size(rr,2)) = rr
      call move_alloc(points_tmp, points)
      allocate(comps(i_comp)%i_points(size(rr,2)))
      comps(i_comp)%i_points = &
                         (/((i3),i3=points_offset+1,points_offset+size(rr,2))/)


      ! --- treat the elements ---

      !allocate the elements of the component of the right kind
      comps(i_comp)%nelems = size(ee,2)
      nelem = nelem + comps(i_comp)%nelems
      select case(trim(comps(i_comp)%comp_el_type))
       case('p')
        allocate(t_surfpan::comps(i_comp)%el(size(ee,2)))
       case('v')
        allocate(t_vortlatt::comps(i_comp)%el(size(ee,2)))
       case('l')
        allocate(t_liftlin::comps(i_comp)%el(size(ee,2)))
       case('a')
        allocate(t_actdisk::comps(i_comp)%el(size(ee,2)))
       case default
        call error(this_sub_name, this_mod_name, &
                 'Unknown type of element: '//comps(i_comp)%comp_el_type)
      end select
        
      !fill (some) of the elements fields
      do i2=1,size(ee,2)
        
        !Component id
        !comps(i_comp)%el(i2)%comp_id = i_comp
        
        !vertices
        n_vert = count(ee(:,i2).ne.0)
        allocate(comps(i_comp)%el(i2)%i_ver(n_vert))
        !allocate(geo%components(i_comp)%el(i2)%neigh(n_vert))
        comps(i_comp)%el(i2)%n_ver = n_vert
        comps(i_comp)%el(i2)%i_ver(1:n_vert) = &
                                              ee(1:n_vert,i2) + points_offset
        !do i3=1,n_vert
        !  if ( neigh(i3,i2) .ne. 0 ) then
        !    geo%components(i_comp)%el(i2)%neigh(i3)%p => &
        !                            geo%components(i_comp)%el(neigh(i3,i2))
        !  else
        !    ! do nothing, keep the neighbour pointer not associated
        !    geo%components(i_comp)%el(i2)%neigh(i3)%p => null()
        !  end if
        !end do
        
        !motion
        !geo%components(i_comp)%el(i2)%moving = geo%components(i_comp)%moving

      enddo
 
      ! Update elems_offset for the next component
      !elems_offset = elems_offset + size(ee,2)

      !cleanup
      deallocate(ee,rr)

    endif !load the element because in list

    call close_hdf5_group(cloc)

  enddo !i_comp
  call close_hdf5_group(gloc)

  !generate the "global" connectivity
  !allocate(ep_conn(4,nelem))
  !ep_conn = 0
  !ie_t = 0
  !do i_comp = 1,n_comp
  !  do ie = 1,comps(i_comp)%nelems
  !  ie_t = ie_t + 1
  !  ep_conn(1:comps(i_comp)%el(ie)%n_ver,ie_t) = comps(i_comp)%el(ie)%i_ver

  !  enddo
  !enddo



end subroutine load_components_postpro

!----------------------------------------------------------------------

!!subroutine import_aero_tab(geo,coeff)
!! type(t_geo), intent(inout), target :: geo
!! type(t_aero_tab) , allocatable , intent(inout) :: coeff(:)
!!
!! integer :: n_tmp , n_tmp2
!! character(len=max_char_len) , allocatable :: list_tmp(:) 
!! character(len=max_char_len) , allocatable :: list_tmp_tmp(:) 
!! integer , allocatable :: i_airfoil_e_tmp (:,:)
!!
!! integer :: i_c, n_c, i_a, n_a, i_l
!!
!! n_tmp = 30 
!! allocate(list_tmp(n_tmp)) 
!!
!! ! Count # of different airfoil
!! n_c = size(geo%components) 
!! n_a = 0
!! do i_c = 1 ,  n_c
!!
!!   if ( geo%components(i_c)%comp_el_type .eq. 'l' ) then
!!    
!!     allocate( i_airfoil_e_tmp( &
!!          size(geo%components(i_c)%i_airfoil_e,1), size(geo%components(i_c)%i_airfoil_e,2) ) ) 
!!     i_airfoil_e_tmp = 0
!!
!!     do i_a = 1 , size(geo%components(i_c)%airfoil_list)
!!
!!       if ( all( geo%components(i_c)%airfoil_list(i_a) .ne. list_tmp(1:n_a) ) ) then
!!         ! new airfoil ----
!!         n_a = n_a + 1
!!
!!         where ( geo%components(i_c)%i_airfoil_e .eq. i_a ) i_airfoil_e_tmp = n_a 
!!        
!!         ! if n_a > n_tmp --> movalloc
!!         if ( n_a .gt. n_tmp ) then
!!           n_tmp2 = n_tmp + n_tmp
!!           allocate(list_tmp_tmp(n_tmp2))
!!           list_tmp_tmp(1:n_tmp) = list_tmp
!!           deallocate(list_tmp)
!!           call move_alloc(list_tmp_tmp,list_tmp)
!!           n_tmp = n_tmp2
!!         end if
!! 
!!         list_tmp(n_a) = geo%components(i_c)%airfoil_list(i_a)
!!
!!       else
!!         ! airfoil already used: find the element and replace the global index
!!         do i_l = 1 , n_a
!!           if ( geo%components(i_c)%airfoil_list(i_a) .eq. list_tmp(i_l) ) exit
!!         end do
!!
!!         where ( geo%components(i_c)%i_airfoil_e .eq. i_a ) i_airfoil_e_tmp = i_l 
!!
!!       end if
!!
!!
!!     end do
!!
!!   geo%components(i_c)%i_airfoil_e = i_airfoil_e_tmp
!!
!!   deallocate(i_airfoil_e_tmp)
!!     
!!   ! check ----
!!   write(*,*) ' mod_geo.f89/import_aero_tab().  Component: ' , i_c
!!   write(*,*) ' size(geo%components(',i_c,')%i_airfoil_e : ' , shape(geo%components(i_c)%i_airfoil_e)
!!   do i_l = 1 , size(geo%components(i_c)%i_airfoil_e,2)
!!     write(*,*) geo%components(i_c)%i_airfoil_e(:,i_l)
!!   end do
!!   write(*,*)
!!   ! check ----
!!
!!   end if
!!
!! end do
!!
!! ! Read tables and fill coeff structure
!! allocate(coeff(n_a))
!! write(*,*) ' Number of different airfoils : ' , n_a
!! do i_a = 1 , n_a
!!   call read_c81_table( list_tmp(i_a) , coeff(i_a) )
!! end do
!!
!!end subroutine import_aero_tab

!----------------------------------------------------------------------

!> Prepare the geometry allocating all the relevant fields for each 
!! kind of element
subroutine prepare_geometry_postpro(comps)
 type(t_geo_component), intent(inout), target :: comps(:)

 integer :: i_comp, ie
 integer :: nsides
 class(c_pot_elem), pointer :: elem
 character(len=*), parameter :: this_sub_name = 'prepare_geometry_postpro'

 do i_comp = 1,size(comps)
   do ie = 1,size(comps(i_comp)%el)
     elem => comps(i_comp)%el(ie)

     nsides = size(elem%i_ver)

     !Fields common to each element
     allocate(elem%ver(3,nsides))

     select type(elem)
      class is(t_surfpan)
       allocate(elem%verp(3,nsides))
       allocate(elem%edge_vec(3,nsides))
       allocate(elem%edge_len(nsides))
       allocate(elem%edge_uni(3,nsides))
       allocate(elem%cosTi(nsides))
       allocate(elem%sinTi(nsides))

      class is(t_vortlatt)
       allocate(elem%edge_vec(3,nsides))
       allocate(elem%edge_len(nsides))
       allocate(elem%edge_uni(3,nsides))

      class is(t_liftlin)
       allocate(elem%edge_vec(3,nsides))
       allocate(elem%edge_len(nsides))
       allocate(elem%edge_uni(3,nsides))
     !!  
     !!  elem%csi_cen = 0.5_wp * sum(geo%components(i_comp)%normalised_coord_e(:,ie))
     !!  elem%i_airfoil =  geo%components(i_comp)%i_airfoil_e(:,ie)
     !!  
     !!  
     !!  ! elem%chord    
     !!  ! elem%csi_c       !! <- for interpolation of aerodynamic coefficients
     !!  ! elem%airfoil(2)  !! 

     !!  ! elem%theta   <-- ??? 
      class is(t_actdisk)
       allocate(elem%edge_vec(3,nsides))
       allocate(elem%edge_len(nsides))
       allocate(elem%edge_uni(3,nsides))

      class default
       call error(this_sub_name, this_mod_name, 'Unknown element type')
     end select

     !end associate
   enddo
 enddo

end subroutine prepare_geometry_postpro

!----------------------------------------------------------------------

subroutine prepare_wake_postpro( wpoints_pan, wpoints_rin, wstart, wconn, &
                  wvort_pan, wvort_rin,  wake_pan, wake_rin, wake_elems )
 real(wp), intent(in) :: wpoints_pan(:,:,:)
 real(wp), intent(in) :: wpoints_rin(:,:,:)
 integer , intent(in) :: wstart(:,:)
 integer , intent(in) :: wconn(:)
 real(wp), intent(in) :: wvort_pan(:,:)
 real(wp), intent(in) :: wvort_rin(:,:)
 type(t_wake_panels), target, intent(out) :: wake_pan
 type(t_wake_rings), target,  intent(out) :: wake_rin
 type(t_elem_p), allocatable, intent(out) :: wake_elems(:)


 integer :: n_wake_stripes , npan, ndisks, nrows
 integer :: nsides
 integer :: p1 , p2 
 integer :: ip , iw, id, ir, iconn, ie
 integer :: npt_disk
 integer, allocatable :: disk_pts(:)

  !First get all the panel wake
  n_wake_stripes = size(wstart ,2)
  npan           = size(wvort_pan  ,2) 

  !TODO: check if the following dimensions are right
  wake_pan%npan = npan
  wake_pan%n_wake_stripes = n_wake_stripes
  wake_pan%wake_len = npan
! wake%n_wake_points = ...
  allocate(wake_pan%i_start_points(2,wake_pan%n_wake_stripes))
  allocate(wake_pan%w_points(3,wake_pan%n_wake_points,npan+1))
  allocate(wake_pan%wake_panels(wake_pan%n_wake_stripes,npan))
  allocate(wake_pan%idou(wake_pan%n_wake_stripes,npan))

  wake_pan%i_start_points = wstart
  wake_pan%w_points       = wpoints_pan

  nsides = 4 
  do ip = 1,npan
    do iw=1,wake_pan%n_wake_stripes
     wake_pan%wake_panels(iw,ip)%mag => wake_pan%idou(iw,ip)
     allocate(wake_pan%wake_panels(iw,ip)%ver(3,nsides))
     allocate(wake_pan%wake_panels(iw,ip)%edge_vec(3,nsides))
     allocate(wake_pan%wake_panels(iw,ip)%edge_len(nsides))
     allocate(wake_pan%wake_panels(iw,ip)%edge_uni(3,nsides))
    enddo
  enddo

  do ip = 1,wake_pan%wake_len
   do iw = 1,wake_pan%n_wake_stripes
       p1 = wake_pan%i_start_points(1,iw)
       p2 = wake_pan%i_start_points(2,iw)
       call calc_geo_data_postpro(wake_pan%wake_panels(iw,ip), &
       reshape((/wake_pan%w_points(:,p1,ip),   wake_pan%w_points(:,p2,ip), &
                 wake_pan%w_points(:,p2,ip+1), wake_pan%w_points(:,p1,ip+1)/),&
                                                                   (/3,4/)))
       wake_pan%wake_panels(iw,ip)%mag = wvort_pan(iw,ip) 
   end do
  end do

  !Then get all the ring wake
  ndisks = size(wvort_rin,1)
  nrows = size(wvort_rin,2)
  wake_rin%ndisks = ndisks; wake_rin%wake_len = nrows
  allocate(wake_rin%wake_rings(wake_rin%ndisks,wake_rin%wake_len))
  allocate(wake_rin%idou(wake_rin%ndisks,wake_rin%wake_len))

  do id = 1,wake_rin%ndisks
    !reverse the connectivity, from pts2disk to disk2pts
    npt_disk = count(wconn .eq. id)
    allocate(disk_pts(npt_disk))
    iconn = 1
    do ip = 1,size(wconn)
      if(wconn(ip) .eq. id) then
        disk_pts(iconn) = ip
        iconn = iconn+1
      endif
    enddo

    nsides = npt_disk
    do ir = 1,wake_rin%wake_len
      wake_rin%wake_rings(id,ir)%mag => wake_rin%idou(id,ir)
      wake_rin%wake_rings(id,ir)%n_ver = nsides
      allocate(wake_rin%wake_rings(id,ir)%ver(3,nsides))
      allocate(wake_rin%wake_rings(id,ir)%edge_vec(3,nsides))
      allocate(wake_rin%wake_rings(id,ir)%edge_len(nsides))
      allocate(wake_rin%wake_rings(id,ir)%edge_uni(3,nsides))
      call calc_geo_data_postpro(wake_rin%wake_rings(id,ir), &
                  wpoints_rin(:,disk_pts,ir))
      wake_rin%wake_rings(id,ir)%mag = wvort_rin(id,ir)
      
    enddo
    deallocate(disk_pts)
  enddo

  !Stitch together everything
  allocate(wake_elems(wake_pan%n_wake_stripes*wake_pan%wake_len + &
                      wake_rin%ndisks*wake_rin%wake_len))
  ie = 1
  do ip = 1,wake_pan%wake_len
   do iw = 1,wake_pan%n_wake_stripes
     wake_elems(ie)%p => wake_pan%wake_panels(iw,ip)
     ie = ie + 1
   end do
  end do

  do ir = 1,wake_rin%wake_len
    do id = 1,wake_rin%ndisks
      wake_elems(ie)%p => wake_rin%wake_rings(id,ir)
      ie = ie+1
    enddo
  enddo


end subroutine prepare_wake_postpro

!----------------------------------------------------------------------

!> Calculate the geometrical quantities of an element
!!
!! The subroutine calculates all the relevant geometrical quantities of a
!! panel element (vortex ring or surface panel) 
subroutine calc_geo_data_postpro(elem,vert)
 class(c_pot_elem), intent(inout) :: elem  
 real(wp), intent(in) :: vert(:,:)

 integer :: nsides, is
 real(wp):: nor(3), tanl(3)
 integer :: nxt

  nsides = size(vert,2)

  elem%n_ver = nsides
  
  ! vertices
  elem%ver = vert

  ! center
  elem%cen =  sum ( vert,2 ) / real(nsides,wp)

  select type(elem)

   class default

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

   class is(t_actdisk)
  
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

  end select

  ! vector connecting two consecutive vertices: 
  ! edge_vec(:,1) =  ver(:,2) - ver(:,1)
  do is = 1 , nsides
    nxt = 1+mod(is,nsides)
    elem%edge_vec(:,is) = vert(:,nxt) - vert(:,is)
  end do

  ! edge: edge_len(:) 
  do is = 1 , nsides
    elem%edge_len(is) = norm2(elem%edge_vec(:,is)) 
  end do

  ! unit vector 
  do is = 1 , nSides
    elem%edge_uni(:,is) = elem%edge_vec(:,is) / elem%edge_len(is)
  end do

  ! variables for surfpan elements only
  select type(elem)
    class is(t_surfpan)
      do is = 1 , nsides
        elem%verp(:,is) = vert(:,is) - elem%nor * &
                          sum( (vert(:,is) - elem%cen ) * elem%nor )
        elem%cosTi(is) = sum( elem%edge_uni(:,is) * elem%tang(:,1) ) 
        elem%sinTi(is) = sum( elem%edge_uni(:,is) * elem%tang(:,2) ) 
      end do

    class default

  end select


end subroutine calc_geo_data_postpro

! in geo/mod_geo.f90 ---------------------------------------------------
!
!!> Calculate the local velocity on the panels to then enforce the 
!!! boundary condition
!!!
!subroutine calc_geo_vel(elem, G, f)
! class(c_elem), intent(inout) :: elem  
! real(wp), intent(in) :: f(3), G(3,3)
!
!  if(.not.allocated(elem%ub)) allocate(elem%ub(3))
!  elem%ub = f + matmul(G,elem%cen)
!
!end subroutine calc_geo_vel

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

subroutine update_points_postpro(comps, points, refs_R, refs_off, refs_G , refs_f)
 type(t_geo_component), intent(inout) :: comps(:)
 real(wp), intent(inout) :: points(:,:)
 real(wp), intent(in)    :: refs_R(:,:,0:)
 real(wp), intent(in)    :: refs_off(:,0:)
 real(wp), optional , intent(in)    :: refs_G(:,:,0:)
 real(wp), optional , intent(in)    :: refs_f(:,0:)

 integer :: i_comp, ie

 do i_comp = 1,size(comps)
  associate(comp => comps(i_comp))
  points(:,comp%i_points) = move_points(comp%loc_points, &
                           refs_R(:,:,comp%ref_id), &
                           refs_off(:,comp%ref_id))

  do ie = 1,size(comp%el)
    call calc_geo_data_postpro(comp%el(ie),points(:,comp%el(ie)%i_ver))

    if ( present(refs_G) .and. present(refs_f) ) then
      !Calculate the velocity of the centers to impose the boundary condition
      call calc_geo_vel(comp%el(ie), refs_G(:,:,comp%ref_id) , &
                                       refs_f(:,comp%ref_id) )
    end if
  enddo

  end associate
 enddo

!DEBUG
!write(*,*) ' ***** '


end subroutine update_points_postpro

!----------------------------------------------------------------------

subroutine expand_actdisk_postpro(comps, points, points_exp, elems)
 type(t_geo_component), intent(in) :: comps(:)
 real(wp), intent(in) :: points(:,:)
 real(wp), allocatable, intent(out) :: points_exp(:,:)
 integer, allocatable, intent(out)  :: elems(:,:)


 real(wp), allocatable :: pt_tmp(:,:)
 integer, allocatable  :: ee_tmp(:,:)
 integer :: i_comp, ie, extra_offset, iv, ipt, ipt1!, ie_t, next
 integer :: start_pts, start_cen

  extra_offset = 0
  allocate(points_exp(3,0), elems(4,0))
  do i_comp = 1,size(comps)
    associate(cmp=>comps(i_comp))
    select type(el => cmp%el)
     type is(t_actdisk) 
      !make space also for the centers
      allocate(pt_tmp(3,size(points_exp,2) + &
               size(cmp%loc_points,2) + cmp%nelems))
      pt_tmp(:,1:size(points_exp,2)) = points_exp
      start_pts = size(points_exp,2)
      start_cen = size(points_exp,2) + size(cmp%i_points)
      pt_tmp(:,start_pts+1 : start_cen) = points(:,cmp%i_points)
      do ie = 1,cmp%nelems
        pt_tmp(:,start_cen+ie) = el(ie)%cen
      enddo
      call move_alloc(pt_tmp, points_exp)
      extra_offset = extra_offset + cmp%nelems

      allocate(ee_tmp(4,size(elems,2)+size(cmp%i_points)))
      ee_tmp(:,1:size(elems,2)) = elems
      ee_tmp(:,size(elems,2)+1:size(ee_tmp,2)) = 0

      ipt = 1
      do ie = 1,cmp%nelems
        ipt1 = ipt
        do iv = 1,el(ie)%n_ver-1
          ee_tmp(1,size(elems,2)+ipt) = start_pts+ipt
          ee_tmp(2,size(elems,2)+ipt) = start_pts+ipt+1
          ee_tmp(3,size(elems,2)+ipt) = start_cen+ie
          ipt = ipt+1
        enddo
        !last element
          ee_tmp(1,size(elems,2)+ipt) = start_pts+ipt
          ee_tmp(2,size(elems,2)+ipt) = start_pts+ipt1
          ee_tmp(3,size(elems,2)+ipt) = start_cen+ie
          ipt = ipt+1
      enddo
      call move_alloc(ee_tmp, elems)

     class default
      allocate(pt_tmp(3,size(points_exp,2)+size(cmp%loc_points,2)))
      pt_tmp(:,1:size(points_exp,2)) = points_exp
      pt_tmp(:,size(points_exp,2)+1:size(pt_tmp,2)) = points(:,cmp%i_points)
      call move_alloc(pt_tmp, points_exp)

      allocate(ee_tmp(4,size(elems,2)+cmp%nelems))
      ee_tmp(:,1:size(elems,2)) = elems
      ee_tmp(:,size(elems,2)+1:size(ee_tmp,2)) = 0
      do ie = 1,cmp%nelems
        ee_tmp(1:el(ie)%n_ver,size(elems,2)+ie) =  &
                                                el(ie)%i_ver+extra_offset 
      enddo
      call move_alloc(ee_tmp, elems)
      


    end select
    end associate
  enddo
end subroutine
!----------------------------------------------------------------------

end module mod_geo_postpro
