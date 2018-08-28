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


!> Module to generate the geometry from different kinds of inputs, from mesh
!! files or parametric input

module mod_build_geo

use mod_param, only: &
  wp, max_char_len, nl, prev_tri, next_tri, prev_qua, next_qua

use mod_parse, only: &
  t_parse, getstr, getint, getreal, getrealarray, getlogical, countoption &
  , finalizeparameters

use mod_handling, only: &
  error, warning, info, printout, dust_time, t_realtime

use mod_basic_io, only: &
  read_mesh_basic, write_basic

use mod_cgns_io, only: &
  read_mesh_cgns

use mod_parametric_io, only: &
  read_mesh_parametric, read_actuatordisk_parametric

use mod_ll_io, only: &
  read_mesh_ll

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

use mod_math, only: &
   cross  

implicit none

public :: build_geometry

private

character(len=*), parameter :: this_mod_name = 'mod_build_geo'

contains 

subroutine build_geometry(geo_files, ref_tags, comp_names, output_file, & 
                          global_tol_sew , global_inner_prod_te )
 character(len=*), intent(in) :: geo_files(:)
 character(len=*), intent(in) :: ref_tags(:)
 character(len=*), intent(in) :: comp_names(:)
 character(len=*), intent(in) :: output_file
 real(wp),intent(in) :: global_tol_sew , global_inner_prod_te

 integer :: n_geo, i_geo
 integer(h5loc) :: file_loc, group_loc

  n_geo = size(geo_files)
  
  call new_hdf5_file(output_file, file_loc)
  call new_hdf5_group(file_loc, 'Components', group_loc)
  call write_hdf5(n_geo,'NComponents',group_loc)
  do i_geo = 1,n_geo
    call build_component(group_loc, trim(geo_files(i_geo)), &
                         trim(ref_tags(i_geo)), trim(comp_names(i_geo)), i_geo , &
                         global_tol_sew , global_inner_prod_te )
  enddo

  call close_hdf5_group(group_loc)
  call close_hdf5_file(file_loc)

end subroutine build_geometry

!-----------------------------------------------------------------------

subroutine build_component(gloc, geo_file, ref_tag, comp_tag, comp_id, &
                           global_tol_sew , global_inner_prod_te)
 integer(h5loc), intent(in)   :: gloc
 character(len=*), intent(in) :: geo_file
 character(len=*), intent(in) :: ref_tag
 character(len=*), intent(in) :: comp_tag
 integer, intent(in)          :: comp_id
 real(wp),intent(in) :: global_tol_sew , global_inner_prod_te

 type(t_parse) :: geo_prs
 character(len=max_char_len) :: mesh_file
 integer, allocatable :: ee(:,:)
 real(wp), allocatable :: rr(:,:)
 real(wp), allocatable :: theta_p(:) , chord_p(:)
 character(len=max_char_len) :: comp_el_type
 character :: ElType
 character(len=max_char_len) :: mesh_file_type
 logical :: mesh_reflection
 real(wp) :: reflection_point(3), reflection_normal(3)
 character(len=max_char_len) :: comp_name
 integer(h5loc) :: comp_loc , geo_loc , te_loc

 character(len=max_char_len), allocatable :: airfoil_list(:)
 integer , allocatable                    :: nelem_span_list(:)
 integer , allocatable                    :: i_airfoil_e(:,:)
 real(wp), allocatable                    :: normalised_coord_e(:,:)
 real(wp) :: trac, radius

 ! Section names for CGNS
 integer :: nSections, iSection
 character(len=max_char_len), allocatable :: sectionNamesCGNS(:)

 integer :: npoints_chord_tot, nelems_span, nelems_span_tot
 ! Connectivity and te structures 
 integer , allocatable :: neigh(:,:)

 ! Parameters for gap sewing and te identification
 real(wp) :: tol_sewing , inner_product_threshold
 integer :: i_count
 ! Projection of the unit tangent exiting from te
 logical  :: te_proj_logical
 character(len=max_char_len) :: te_proj_dir
 real(wp) , allocatable :: te_proj_vec(:)

 ! trailing edge ------
 integer , allocatable :: e_te(:,:) , i_te(:,:) , ii_te(:,:)
 integer , allocatable :: neigh_te(:,:) , o_te(:,:)
 real(wp), allocatable :: rr_te(:,:) , t_te(:,:)

 character(len=*), parameter :: this_sub_name = 'build_component'


  !Prepare all the parameters to be read in the file
  !mesh file
  call geo_prs%CreateStringOption('MeshFile','Mesh file definition')
  call geo_prs%CreateStringOption('MeshFileType','Mesh file type')
  !element types
  call geo_prs%CreateStringOption('ElType', &
              'element type (temporary) p:panel, v:vortex ring, l:lifting line')
  !reference frame
  !call geo_prs%CreateStringOption('Reference_Tag',&
  !                                 'reference frame tag of the component','0')
  !reflections
  call geo_prs%CreateLogicalOption('mesh_reflection',&
               'Has all the file a custom reference frame', 'F')
  call geo_prs%CreateRealArrayOption('reflection_point', &
               'Center point of reflection plane, (x, y, z)', &
               '(/0.0, 0.0, 0.0/)')
  call geo_prs%CreateRealArrayOption('reflection_normal', &
               'Normal of reflection plane, (xn, yn, zn)', &
               '(/0.0, 0.0, 1.0/)')
  
  call geo_prs%CreateRealOption('traction', &
               'Traction of the rotor')
  call geo_prs%CreateRealOption('Radius', &
               'Radius of the rotor')

  ! Parameters for gap sewing and edge identification
  call geo_prs%CreateRealOption('TolSewing', &
               'Dimension of the gaps that will be filled, to close TE', &
               '0.001')  ! very rough default value
  call geo_prs%CreateRealOption('InnerProductTe', &
               'Threshold of the inner product of adiacent elements &
               &across the TE','-0.5')  ! very rough default value

  ! Parameters for te unit tangent vector generation
  call geo_prs%CreateLogicalOption('ProjTe','Remove some component &
               &from te vectors.')
  call geo_prs%CreateStringOption('ProjTeDir','Parallel or normal to &
               &ProjTeVector direction.','parallel') 
  call geo_prs%CreateRealArrayOption('ProjTeVector','Vector used for &
               &the te projection.')

  ! Section name from CGNS file
  call geo_prs%CreateStringOption('SectionName', &
               'Section name from CGNS file', multiple=.true.)
  
  !read the parameters
  call geo_prs%read_options(geo_file,printout_val=.false.)

  mesh_file_type = getstr(geo_prs,'MeshFileType')
  !ref_tag        = getstr(geo_prs,'Reference_Tag')

  mesh_reflection   = getlogical(geo_prs, 'mesh_reflection')
  reflection_point  = getrealarray(geo_prs, 'reflection_point',3)
  reflection_normal = getrealarray(geo_prs, 'reflection_normal',3)

  comp_el_type = getstr(geo_prs,'ElType')
  ElType = comp_el_type(1:1)
  if ( ( trim(ElType) .ne. 'p' ) .and. ( trim(ElType) .ne. 'v') .and. & 
       ( trim(ElType) .ne. 'l' ) .and. ( trim(ElType) .ne. 'a') ) then
    call error (this_sub_name, this_mod_name, 'Wrong ElType ('// &
         trim(ElType)//') for component'//trim(comp_tag)// &
         ', defined in file: '//trim(geo_file)//'. ElType must be:'// &
         ' p,v,l,a.')
  end if

  ! Parameters for gap sewing and edge identification ------------------
  ! TolSewing --------------------------------------
  i_count = countoption(geo_prs,'TolSewing') 
  if ( i_count .eq. 1 ) then
    tol_sewing = getreal(geo_prs,'TolSewing')
  else if ( i_count .eq. 0 ) then
    ! Inherit tol_sewing from general parameter
    tol_sewing = global_tol_sew 
  else
    tol_sewing = getreal(geo_prs,'TolSewing')

    call warning(this_sub_name, this_mod_name, 'More than one &
       &TolSewing parameter defined for component'//trim(comp_tag)// &
       ', defined in file: '//trim(geo_file)//'.')
  end if
 
  ! InnerProductTe ---------------------------------
  i_count = countoption(geo_prs,'InnerProductTe') 
  if ( i_count .eq. 1 ) then
    inner_product_threshold = getreal(geo_prs,'InnerProductTe')
  else if ( i_count .eq. 0 ) then
    ! Inherit tol_sewing from general parameter
    inner_product_threshold = global_inner_prod_te
  else
    inner_product_threshold = getreal(geo_prs,'InnerProductTe')

    call warning(this_sub_name, this_mod_name, 'More than one &
       &InnerProductTe parameter defined for component'//trim(comp_tag)// &
       ', defined in file: '//trim(geo_file)//'. First value used.')
  end if 

  ! te projection ----------------------------------
  i_count = countoption(geo_prs,'ProjTe')
  te_proj_logical = .false.
  if ( i_count .eq. 1 ) then
    te_proj_logical = getlogical(geo_prs,'ProjTe') 
  end if
   
  ! ------
  if ( te_proj_logical ) then
    i_count = countoption(geo_prs,'ProjTeDir' )
    if ( i_count .eq. 1 ) then
      te_proj_dir  = getstr(geo_prs, 'ProjTeDir')
    elseif ( i_count .lt. 1 ) then
      te_proj_dir = 'parallel'
      call warning (this_sub_name, this_mod_name, 'ProjTe = T, but &
         & no ProjTeVDir defined for component'//trim(comp_tag)// &
         ', defined in file: '//trim(geo_file)//" Default: 'parallel'.")
    else
      te_proj_dir = getstr(geo_prs, 'ProjTeDir') 
      call warning(this_sub_name, this_mod_name, 'ProjTe = T, and &
         & more than one ProjTeDir defined for component'//trim(comp_tag)// &
         ', defined in file: '//trim(geo_file)//'. First value used.')
    end if
    if ( ( trim(te_proj_dir) .ne. 'parallel' ) .and. &
         ( trim(te_proj_dir) .ne. 'normal'   )         ) then 
      call error (this_sub_name, this_mod_name, "Wrong 'ProjTeDir' &
         & input for component"//trim(comp_tag)// &
         ', defined in file: '//trim(geo_file)//'.')
    end if 
  end if
  ! ------
  if ( te_proj_logical ) then
    i_count = countoption(geo_prs,'ProjTeVector' )
    if ( i_count .eq. 1 ) then
      te_proj_vec  = getrealarray(geo_prs, 'ProjTeVector',3)
    elseif ( i_count .lt. 1 ) then
      call error (this_sub_name, this_mod_name, 'ProjTe = T, but &
         & no ProjTeVector defined for component'//trim(comp_tag)// &
         ', defined in file: '//trim(geo_file)//'.')
    else 
      te_proj_vec  = getrealarray(geo_prs, 'ProjTeVector',3)
      call warning(this_sub_name, this_mod_name, 'ProjTe = T, but &
         & more than one ProjTeVector defined for component'//trim(comp_tag)// &
         ', defined in file: '//trim(geo_file)//'. First value used.')
    end if
  end if


  !Build the group -----------------------------------------------------
  write(comp_name,'(A,I3.3)')'Comp',comp_id
  call new_hdf5_group(gloc, trim(comp_name), comp_loc)
  call write_hdf5(comp_tag,'CompName',comp_loc)
  call write_hdf5(mesh_file_type,'CompInput',comp_loc)
  call write_hdf5(ref_tag,'RefTag',comp_loc)

! call new_hdf5_group(file_loc, 'Components', group_loc)
  call write_hdf5(trim(comp_el_type),'ElType',comp_loc)

  call new_hdf5_group(comp_loc, 'Geometry', geo_loc)

  ! read the files
  select case (trim(mesh_file_type))

   case('basic')
    mesh_file = getstr(geo_prs,'MeshFile')
    call read_mesh_basic(trim(mesh_file),ee, rr)

   case('cgns')
    mesh_file = getstr(geo_prs,'MeshFile')

    ! Check the selection of mesh sections
    nSections = countoption(geo_prs, 'SectionName')

    allocate(sectionNamesCGNS(nSections))

    do iSection = 1, nSections
      sectionNamesCGNS(iSection) = trim(getstr(geo_prs, 'SectionName'))
    enddo

    call read_mesh_cgns(trim(mesh_file), sectionNamesCGNS,  ee, rr)

   case('parametric')
    mesh_file = geo_file
    if ( (ElType .eq. 'v') .or. (ElType .eq. 'p')  ) then
      ! TODO : actually it is possible to define the parameters in the GeoFile 
      !directly, find a good way to do this
      call read_mesh_parametric(trim(mesh_file),ee, rr, &
                                npoints_chord_tot, nelems_span )
      ! nelems_span_tot will be overwritten if symmetry is required (around l.220)
      nelems_span_tot =   nelems_span

    elseif(ElType .eq. 'l') then ! LIFTING LINE element
      call read_mesh_ll(trim(mesh_file),ee,rr, &
                        airfoil_list   , nelem_span_list   , &
                        i_airfoil_e    , normalised_coord_e, &
                                npoints_chord_tot, nelems_span, &
                                chord_p,theta_p )
      ! nelems_span_tot will be overwritten if symmetry is required (around l.220)
      nelems_span_tot =   nelems_span

      call write_hdf5(airfoil_list   ,'airfoil_list'   ,geo_loc)
      call write_hdf5(nelem_span_list,'nelem_span_list',geo_loc)
      call write_hdf5(theta_p,'theta_p',geo_loc)
      call write_hdf5(chord_p,'chord_p',geo_loc)
      call write_hdf5(i_airfoil_e,'i_airfoil_e',geo_loc)
      call write_hdf5(normalised_coord_e,'normalised_coord_e',geo_loc)

    elseif(ElType .eq. 'a') then ! ACTUATOR DISK
      call read_actuatordisk_parametric(trim(mesh_file),ee,rr)
      trac = getreal(geo_prs,'traction')
      call write_hdf5(trac,'Traction',comp_loc)
      radius = getreal(geo_prs,'Radius')
      call write_hdf5(radius,'Radius', comp_loc)

    end if
   case default
    call error(this_sub_name, this_mod_name, 'Unknown mesh file type')

  end select

  ! reflect the mesh (if requested)
  write(*,*) 'mesh_reflection' , mesh_reflection
  if ( mesh_reflection ) then
    select case (trim(mesh_file_type))
      case('cgns', 'basic' )  ! TODO: check basic
        call reflect_mesh(ee, rr, reflection_point, reflection_normal)
      case('parametric')
!! !       if ( ElType .ne. 'l' ) then
          call reflect_mesh_structured(ee, rr,  &
                                       npoints_chord_tot , nelems_span , &
                                       reflection_point, reflection_normal)
          nelems_span_tot = 2*nelems_span

!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !!
!! Same routines used for all parametric elements: p,v,l        !!
!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !!
!! !       else ! LIFTING LINE element
!! !         call reflect_mesh_structured(ee, rr,  &
!! !                                      npoints_chord_tot , nelems_span , &
!! !                                      reflection_point, reflection_normal)
!! !         nelems_span_tot = 2*nelems_span
!! !       end if
      case default
       call error(this_sub_name, this_mod_name,&
             'Symmetry routines not implemented for this kind of input.')
    end select
  end if

  if ( trim(mesh_file_type) .eq. 'parametric' ) then
    !DEBUG
    write(*,*) ' +++ nelems_span_tot : ' , nelems_span_tot
    write(*,*) ' +++ nelems_chor_tot : ' , npoints_chord_tot - 1
    !write HDF5 fields
    call write_hdf5(  nelems_span_tot  ,'parametric_nelems_span',geo_loc)
    call write_hdf5(npoints_chord_tot-1,'parametric_nelems_chor',geo_loc)

  end if


  call write_hdf5(rr,'rr',geo_loc)
  call write_hdf5(ee,'ee',geo_loc)

  !:::::::::::::::::::::::::::::::::::::::::::::::::::

  ! selectcase('cgns','parametric','basic') ---> build structures
  ! ---- local numbering ----
  !  -> build_connectivity                     || <-- these routines depend 
  !  -> build_te                               ||     on the way a component 
  !  -> create_local_velocity_stencil ( 3dP )  ||     is defined
  !  -> create_strip_connecivity      ( vr )   ||
  selectcase(trim(mesh_file_type))
    case( 'basic' )

! TODO: add WARNING and CHECK conventions if ElType = 'v'

      call build_connectivity_general( ee , neigh )

      if ( ElType .eq. 'v' ) then
        write(*,*) nl//' WARNING: component with id.', comp_id
        write(*,*) '  defined as ''basic'' with ''vortex ring'' elements:'
        write(*,*) '  be sure to your input files comply with the conventions '
        write(*,*) '  of parameteric component structures.'//nl
        call build_te_parametric( ee , rr , ElType ,  &
           npoints_chord_tot , nelems_span_tot , &
           e_te, i_te, rr_te, ii_te, neigh_te, o_te, t_te ) !te as an output
      else
        call build_te_general ( ee , rr , ElType , &
                  tol_sewing , inner_product_threshold , &
                  te_proj_logical , te_proj_dir , te_proj_vec , &
                  e_te, i_te, rr_te, ii_te, neigh_te, o_te, t_te ) 
                                                            !te as an output
      end if

    case( 'cgns' )

      call build_connectivity_general( ee , neigh )

      call build_te_general ( ee , rr , ElType , &
                tol_sewing , inner_product_threshold , &
                te_proj_logical , te_proj_dir , te_proj_vec , &
                e_te, i_te, rr_te, ii_te, neigh_te, o_te, t_te ) 
                                                          !te as an output

    case( 'parametric' )
     if ( ElType .eq. 'l' .or. ElType .eq. 'v' .or. ElType .eq. 'p' ) then
        call build_connectivity_parametric( ee ,     &
                       npoints_chord_tot , nelems_span_tot , &
                       neigh )
        call build_te_parametric( ee , rr , ElType ,  &
           npoints_chord_tot , nelems_span_tot , &
           e_te, i_te, rr_te, ii_te, neigh_te, o_te, t_te ) !te as an output
     elseif(ElType .eq. 'a') then
       allocate(neigh(size(ee,1), size(ee,2)))
       neigh = 0 !All actuator disk are isolated
     endif


!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !!
!! Same routines used for all parametric elements: p,v,l        !!
!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !!
!!       else
!!         call build_connectivity_parametric( trim(mesh_file) , ee ,     &
!!                       ElType , npoints_chord_tot , nelems_span_tot , &
!!                       mesh_reflection , neigh )
!!         call build_te_parametric( ee , rr , ElType ,  &
!!            npoints_chord_tot , nelems_span_tot , &
!!            e_te, i_te, rr_te, ii_te, neigh_te, o_te, t_te ) !te as an output
!! 
!! !       call create_local_velocity_stencil_parametric()
!! !       call create_strip_connectivity_parametric()
!! 
!!       end if
    case default
       call error(this_sub_name, this_mod_name,&
             'Unknown option for building connectivity.')
  end select    

  call write_hdf5(neigh,'neigh',geo_loc)
  call close_hdf5_group(geo_loc)

  if (ElType .ne. 'a') then
    call new_hdf5_group(comp_loc, 'Trailing_Edge', te_loc)
    call write_hdf5(    e_te,    'e_te',te_loc)
    call write_hdf5(    i_te,    'i_te',te_loc)
    call write_hdf5(   rr_te,   'rr_te',te_loc)
    call write_hdf5(   ii_te,   'ii_te',te_loc)
    call write_hdf5(neigh_te,'neigh_te',te_loc)
    call write_hdf5(    o_te,    'o_te',te_loc)
    call write_hdf5(    t_te,    't_te',te_loc)
    call close_hdf5_group(te_loc)
  endif

 
  call close_hdf5_group(comp_loc)
  !cleanup
  deallocate(ee,rr)
  call finalizeparameters(geo_prs)

end subroutine build_component

!----------------------------------------------------------------------

!> Subroutine used to double the mesh by reflecting it along a simmetry 
!! plane
!!
!! Given a plane defined by a center point and a normal vector, the mesh 
!! is doubled: all the points are reflected and new mirrored elements 
!! introduced. The elements and points arrays are doubled. 
subroutine reflect_mesh(ee, rr, cent, norm)
 integer, allocatable, intent(inout) :: ee(:,:)
 real(wp), allocatable, intent(inout) :: rr(:,:)
 real(wp), intent(in) :: cent(3), norm(3)

 real(wp) :: n(3), d, l
 integer, allocatable :: ee_temp(:,:)
 real(wp), allocatable :: rr_temp(:,:)
 integer :: ip, np, ne
 integer :: ie, iv, nv

  ne = size(ee,2); np = size(rr,2)
  
  ! enlarge size
  allocate(ee_temp(size(ee,1),2*ne))
  allocate(rr_temp(size(rr,1),2*np))
 
  !first part equal
  ee_temp(:,1:ne) = ee
  rr_temp(:,1:np) = rr
  
  !second part of the elements: index incremented, need to rearrange the 
  !connectivity to preserve the normal
  ee_temp(:,ne+1:2*ne) = 0
  do ie = 1,ne
    nv = count(ee(:,ie).ne.0)
    ee_temp(1,ne+ie) = np+ee(1,ie)
    do iv = 2,nv
      ee_temp(iv,ne+ie) = np+ee(nv-iv+2,ie)
    enddo
  enddo
 
  !calculate normal unit vector and distance from origin
  n = norm/norm2(norm) 
  l = sum(cent * n)
 
  !now reflect the points
  do ip=1,np
    d = sum( rr(:,ip) * n) - l 
    rr_temp(:,np+ip) = rr_temp(:,ip) - 2*d*n
  enddo

  !move alloc back to the original vectors
  call move_alloc(rr_temp, rr)
  call move_alloc(ee_temp, ee)

end subroutine reflect_mesh

!----------------------------------------------------------------------

!> Subroutine used to double the mesh by reflecting it along a simmetry 
!! plane
!!
!! Given a plane defined by a center point and a normal vector, the mesh 
!! is doubled: all the points are reflected and new mirrored elements 
!! introduced. The elements and points arrays are doubled.
subroutine reflect_mesh_structured( ee, rr, &
              npoints_chord_tot , nelems_span , cent, norm )
 integer, allocatable, intent(inout) :: ee(:,:)
 real(wp), allocatable, intent(inout) :: rr(:,:)
 real(wp), intent(in) :: cent(3), norm(3)
 integer , intent(in) :: npoints_chord_tot , nelems_span ! <- unused input

 integer :: nelems_chord
 real(wp) :: n(3), d, l
 integer , allocatable :: ee_temp(:,:) , ee_sort(:,:)
 real(wp), allocatable :: rr_temp(:,:)
 integer :: ip, np, ne
 integer :: ie, nv
 ! some checks  and warnings -----
 real(wp) :: mabs  , m
 real(wp) :: minmaxPn 
 real(wp), parameter :: eps = 1e-6_wp ! TODO: move it as an input
 integer :: imabs, i1, nsew

 character(len=*), parameter :: this_sub_name = 'reflect_mesh_structured'

 ! some checks ---------------------------------
 mabs = 0.0_wp
 do i1  = 1 , size(rr,2)
   if (      abs( sum((rr(:,i1)-cent)*norm) ) .gt. mabs ) then
     mabs  = abs( sum((rr(:,i1)-cent)*norm) ) ; imabs = i1
   end if
 end do
 m = sum( (rr(:,imabs) - cent) * norm )
 if ( m .gt. 0.0_wp ) then
   minmaxPn = sum( (rr(:,1)-cent)*norm )
   do i1  = 2 , size(rr,2)
     if (     sum( (rr(:,i1)-cent)*norm ) .lt. minmaxPn ) then 
       minmaxPn = sum( (rr(:,1)-cent)*norm )
     end if
   end do
   if ( minmaxPn .gt. eps ) then
!    write(*,*) ' warning: discontinuous component '
     call error(this_sub_name, this_mod_name, 'Discontinuous component.')
   else if ( minmaxPn .lt. -eps ) then
     call error(this_sub_name, this_mod_name, 'Body compenetration.')
   else
     write(*,*) ' sewing ' 
   end if
 else if ( m .lt. 0.0_wp ) then
   minmaxPn = sum( (rr(:,1)-cent)*norm )
   do i1  = 2 , size(rr,2)
     if (     sum( (rr(:,i1)-cent)*norm ) .gt. minmaxPn ) then 
       minmaxPn = sum( (rr(:,1)-cent)*norm )
     end if
   end do
   if ( minmaxPn .lt. -eps ) then
!    write(*,*) ' warning: discontinuous component '
     call error(this_sub_name, this_mod_name, 'Discontinuous component.')
   else if ( minmaxPn .gt. eps ) then
     call error(this_sub_name, this_mod_name, 'Body compenetration.')
   else
     write(*,*) ' sewing ' 
   end if
 end if

 ! a whole section must be sewed ------
 nsew = 0
 do i1 = 1 , size(rr,2)
   if ( abs(sum( (rr(:,i1)-cent)*norm )) .lt. eps ) then
     nsew = nsew + 1
   end if
 end do

 if ( nsew .ne. npoints_chord_tot ) then
   write(*,*) ' rr. size : ' , size(rr,1) , size(rr,2)
   write(*,*) ' npoints_chord_tot ' , npoints_chord_tot
   write(*,*) ' nsew = ' , nsew
   call error(this_sub_name, this_mod_name, 'Wrong sewing of the first section')
 else
   write(*,*) ' nsew = ' , nsew
 end if
 ! some checks ---------------------------------


  ne = size(ee,2); np = size(rr,2)
  
  ! enlarge size
  allocate(ee_temp(size(ee,1),2*ne))
  allocate(rr_temp(size(rr,1),2*np-npoints_chord_tot))

  !first part equal
  ee_temp(:,1:ne) = ee
  rr_temp(:,1:np) = rr
  
  !second part of the elements: index incremented, need to rearrange the 
  !connectivity to preserve the normal
  ee_temp(:,ne+1:2*ne) = 0
  do ie = 1,ne
    nv = count(ee(:,ie).ne.0)
    ! Quad elements only for structured meshes
    ee_temp(1,ne+ie) = np-npoints_chord_tot+ee(2,ie)
    ee_temp(2,ne+ie) = np-npoints_chord_tot+ee(1,ie)
    ee_temp(3,ne+ie) = np-npoints_chord_tot+ee(4,ie)
    ee_temp(4,ne+ie) = np-npoints_chord_tot+ee(3,ie)
  enddo
  ! correct the first elements ----- 
  ee_temp(1,ne+(/(i1,i1=1,npoints_chord_tot-1)/)) = (/(i1,i1=1,npoints_chord_tot-1)/)
  ee_temp(4,ne+(/(i1,i1=1,npoints_chord_tot-1)/)) = (/(i1,i1=2,npoints_chord_tot  )/)
 
  !calculate normal unit vector and distance from origin
  n = norm/norm2(norm) 
  l = sum(cent * n)
 
  !now reflect the points
  do ip=1,np-npoints_chord_tot  
    d = sum( rr(:,ip+npoints_chord_tot) * n) - l 
    rr_temp(:,np+ip) = rr_temp(:,npoints_chord_tot+ip) - 2*d*n
  enddo

  ! sort ee array  
  nelems_chord = npoints_chord_tot-1
  allocate(ee_sort(size(ee_temp,1),size(ee_temp,2)) ) ; ee_sort = 0
  do i1 = 1 , nelems_span
    ee_sort(:,1+(i1-1)*nelems_chord:i1*nelems_chord) = &
       ee_temp(:,2*ne-i1*nelems_chord+1:2*ne-(i1-1)*nelems_chord)
  end do
  do i1 = 1 , nelems_span
    ee_sort(:,1+(i1-1)*nelems_chord+ne:i1*nelems_chord+ne) = &
       ee_temp(:,1+(i1-1)*nelems_chord:i1*nelems_chord)
  end do

! !move alloc back to the original vectors
  call move_alloc(rr_temp, rr)
  call move_alloc(ee_sort, ee)


  if ( allocated(ee_temp) )  deallocate(ee_temp)

end subroutine reflect_mesh_structured

!----------------------------------------------------------------------
!----------------------------------------------------------------------
! routines : build_connectivity            , 
!            build_te                      ,
!            create_local_velocity_stencil ,
!            create_strip_connectivity
! for ('cgns','basic') and ('parametric')
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! TODO : connectivity is lost at the symmetry plane ---> fix it
subroutine build_connectivity_general ( ee , neigh )

 integer, allocatable, intent(in)  :: ee(:,:)
 integer, allocatable, intent(out) :: neigh(:,:)
 
 integer :: nelems , nverts
 
 integer :: nSides1 , nSides2
 integer :: vert1 , vert2
 integer :: ie1 , ie2 , iv1 , iv2

 character(len=*), parameter :: this_sub_name = 'build_connectivity_general'

 nelems = size( ee , 2 )
 nverts = 4
 
 if ( size(ee,1) .ne. nverts ) then
   write(*,*) ' In build_connectivity_general(), wrong size(ee,1).'
   write(*,*) '  size(ee,1)=',size(ee,1),' .ne. 4. Stop.'
   stop
 end if
 
 allocate( neigh(size(ee,1),size(ee,2)) ) ; neigh = 0
 
 do ie1 = 1 , nelems
   nSides1 = nverts - count( ee(:,ie1) .eq. 0 )
 
   do iv1 = 1 , nSides1
 
     ! if ( neigh(iv1,ie1) .ne. 0 ) cycle   ! <---- is it useful ?
     vert1 = ee(iv1,ie1)
     vert2 = ee(mod(iv1,nSides1)+1,ie1)
     
     do ie2 = ie1+1  , nelems
         
       if ( any(ee(:,ie2) .eq. vert1 ) .and. &
            any(ee(:,ie2) .eq. vert2 ) ) then
       
         nSides2 = nverts - count( ee(:,ie2) .eq. 0 )
         
         do iv2 = 1 , nSides2
           
           if ( ee(iv2,ie2) .eq. vert2 ) then
 
             if ( ee(mod(iv2,nSides2)+1,ie2) .eq. vert1 ) then
               neigh(iv1,ie1) = ie2
               neigh(iv2,ie2) = ie1
             else
               ! debug -----         
               write(*,*) ' elements ie1 , ie2 ' , ie1 , ie2 
               write(*,*) ' vertices ie1 :     ' , ee(:,ie1) 
               write(*,*) ' vertices ie2 :     ' , ee(:,ie2) 
               write(*,*) ' vert1 , vert2:     ' , vert1 , vert2
               ! debug -----          
               call error(this_sub_name, this_mod_name, &
                 'Neighbouring elements with opposed normal orientation,&
                 &this is not allowed. Stop')
             
             end if
 
           end if
         end do
       end if
     end do
   end do
 end do

 write(*,*) this_mod_name , ':' ,this_sub_name ,' ... ' , 'done.'


end subroutine build_connectivity_general

!----------------------------------------------------------------------
!> Parameteric definition of the component.
!! Quad elements allowed only

! WARNING: no fairing at the wing tip is allowed (up to now)
subroutine build_connectivity_parametric ( ee , & 
                 npoints_chord_tot , nelems_span , neigh )

integer , allocatable , intent(in) :: ee(:,:)
integer , intent(in) :: npoints_chord_tot , nelems_span
integer , allocatable , intent(out):: neigh(:,:)

integer :: nelems_chord ! , nelems_span
integer :: i1 , i2 , iel
character(len=*), parameter :: this_sub_name = 'build_component_parametric'

if ( allocated(neigh) )   deallocate(neigh)
allocate(neigh(4,size(ee,2))) ; neigh = 0

nelems_chord = npoints_chord_tot - 1
  ! inner strips -----
  do i1 = 1 , nelems_span
    do i2 = 1 , nelems_chord
      iel = i2+(i1-1) * nelems_chord
      neigh(1,iel) = iel - 1
      neigh(2,iel) = iel - nelems_chord
      neigh(3,iel) = iel + 1
      neigh(4,iel) = iel + nelems_chord
    end do
  end do
  ! correct tip   elements
  neigh(2,(/ (    i1                          , i1 = 1,nelems_chord) /)) = 0
  neigh(4,(/ (i1+(nelems_span-1)*nelems_chord , i1 = 1,nelems_chord) /)) = 0
  ! correct le-te elements
  neigh(1,(/ ( 1+(i1-1)*nelems_chord , i1 = 1,nelems_span) /)) = 0
  neigh(3,(/ (   (i1  )*nelems_chord , i1 = 1,nelems_span) /)) = 0


end subroutine build_connectivity_parametric

!----------------------------------------------------------------------

subroutine build_te_general ( ee , rr , ElType , &
                 tol_sewing , inner_prod_thresh ,  &
                 te_proj_logical , te_proj_dir , te_proj_vec , &
                 e_te, i_te, rr_te, ii_te, neigh_te, o_te, t_te ) 
                                                          !te as an output
 integer   , intent(in) :: ee(:,:)
 real(wp)  , intent(in) :: rr(:,:)
 character , intent(in) :: ElType
 real(wp)  , intent(in) :: tol_sewing
 real(wp)  , intent(in) :: inner_prod_thresh
 logical   , intent(in) :: te_proj_logical
 character(len=*), intent(in) :: te_proj_dir
 real(wp)  , intent(in) :: te_proj_vec(:)

 ! te structures
 integer , allocatable :: e_te(:,:) , i_te(:,:) , ii_te(:,:)
 integer , allocatable :: neigh_te(:,:) , o_te(:,:)
 real(wp), allocatable :: rr_te(:,:) , t_te(:,:)

 ! merge ------
!! 1e-0 for nasa-crm , 3e-3_wp for naca0012      !!! old, hardcoded params
!real(wp),   parameter :: tol_sewing = 1e-3_wp   !!! now tol_sewing is an input 
!!real(wp),   parameter :: tol_sewing = 1e-0_wp  !!!
 real(wp), allocatable :: rr_m(:,:)
 integer , allocatable :: ee_m(:,:) , i_m(:,:)
 ! 'closed-te' connectivity -----
 integer , allocatable :: neigh_m(:,:)
 character(len=*), parameter :: this_sub_name = 'build_te_general'

 if ( ElType .ne. 'p' ) then
   call error(this_sub_name, this_mod_name, &
       'element type for a cgns file can only be panel. ElType = ''p'' ' )
 end if

 write(*,*) ' In ' , this_sub_name

 call merge_nodes_general ( rr , ee , tol_sewing , rr_m , ee_m , i_m  ) 

 write(*,*) '   merge_nodes_general ... done.'

 call build_connectivity_general ( ee_m , neigh_m )

 write(*,*) '   build_connectivity_general ... done.'


 call find_te_general ( rr , ee , neigh_m , inner_prod_thresh , &  
                te_proj_logical , te_proj_dir , te_proj_vec , &
                e_te, i_te, rr_te, ii_te, neigh_te, o_te, t_te ) 
                                                         !te as an output


contains

! -------------

subroutine merge_nodes_general ( ri , ei , tol , rr , ee , im )
 real(wp), intent(in) :: ri(:,:)
 integer , intent(in) :: ei(:,:)
 real(wp), intent(in) :: tol
 real(wp), allocatable , intent(out) :: rr(:,:)
 integer , allocatable , intent(out) :: ee(:,:) , im(:,:)
 
 integer , allocatable :: im_tmp(:,:)
 integer :: nSides1
 integer :: n_nodes , n_elems
 integer :: n_merge
 
 integer :: i1 , i2 , i_e , i_v 
 
 character(len=*), parameter :: this_sub_name = 'merge_nodes_general'

  n_nodes = size(ri,2)
  n_elems = size(ei,2)
  
  allocate(rr(size(ri,1),size(ri,2))) ; rr = ri
  allocate(ee(size(ei,1),size(ei,2))) ; ee = ei
  
  allocate(im_tmp(2,n_nodes)) ; im_tmp = 0
  
  n_merge = 0
  do i1 = 1 , n_nodes
  
    do i2 = i1 + 1 , n_nodes 
    
      if ( norm2(ri(:,i2)-ri(:,i1)) .lt. tol ) then
  
        n_merge = n_merge + 1
        rr(:,i1) = 0.5_wp * (ri(:,i1)+ri(:,i2))
        rr(:,i2) = rr(:,i1)
      
        do i_e = 1 , size(ei,2)
          
          ! 0 assumed to be possible only in ei(4,:)
          nSides1 = count( ei(:,i_e) .ne. 0 ) 
        
          do i_v = 1 , nSides1
            ! avoid merging nodes belonging to the same element
            if ( ei(i_v,i_e) .eq. i2 ) then
              if ( all(ee(:,i_e) .ne. i1 ) ) then
        
                ee(i_v,i_e) = i1
  
                im_tmp(:,n_merge) = (/ i1 , i2 /)
        
              end if
            end if
          end do
        end do
      end if
    end do
  end do
  write(*,*) ' n_merge : ' , n_merge
  
  allocate(im(2,n_merge)) ; im = im_tmp(:,1:n_merge)
  
  deallocate(im_tmp)

end subroutine merge_nodes_general

! -------------

subroutine find_te_general ( rr , ee , neigh_m , inner_prod_thresh , & 
                te_proj_logical , te_proj_dir , te_proj_vec ,  & 
                e_te, i_te, rr_te, ii_te, neigh_te, o_te, t_te ) 
                                                         !te as an output

 real(wp), intent(in) :: rr(:,:)
 integer , intent(in) :: ee(:,:) , neigh_m(:,:)
 real(wp), intent(in) :: inner_prod_thresh
 logical   , intent(in) :: te_proj_logical
 character(len=*), intent(in) :: te_proj_dir
 real(wp)  , intent(in) :: te_proj_vec(:)
 ! actual arrays -----
 integer , allocatable , intent(out) :: e_te(:,:) , i_te(:,:) , ii_te(:,:) 
 integer , allocatable , intent(out) :: neigh_te(:,:) , o_te(:,:)
 real(wp), allocatable , intent(out) :: rr_te(:,:) , t_te(:,:) 
 
 ! tmp arrays --------
 integer , allocatable               :: e_te_tmp(:,:) , i_te_tmp(:,:) , ii_te_tmp(:,:) 
 real(wp), allocatable               ::rr_te_tmp(:,:)
 integer , allocatable               :: i_el_nodes_tmp(:,:)

 real(wp), allocatable :: nor(:,:) , cen(:,:) , area(:)

 integer :: nSides , nSides_b , nSides1 , nSides2
 integer :: i_e  , i_b , i_n , e1 , e2
 integer :: n_el , n_p , ne_te , nn_te

 integer  :: ind1 , ind2 , i_node1 , i_node2 , i_node1_max , i_node2_max
 real(wp) ::  mi1 ,  mi2

 real(wp) :: t_te_len
 integer  :: t_te_nelem

!real(wp), parameter :: inner_prod_thresh = - 0.5d0 ! old hardcoded, now an input

!real(wp), dimension(3) , parameter :: u_rel = (/ 1.0_wp , 0.0_wp , 0.0_wp /) 
!                                               ! hard-coded parameter ... 
!real(wp), dimension(3) , parameter :: side_dir = (/ 0.0_wp , 1.0_wp , 0.0_wp /) 
!                                                ! hard-coded parameter ...
 real(wp) , dimension(3) :: side_dir 
 ! TODO: read as an input of the component

 integer :: i1 , i2 
 character(len=*), parameter :: this_sub_name = 'find_te_general'

 n_el = size(ee,2)
 n_p  = size(rr,2)

 ! allocate tmp structures --
 allocate( e_te_tmp(2,(n_el+1)/2) )                ;  e_te_tmp = 0 
 allocate( i_el_nodes_tmp(2,n_p ) ) ;  i_el_nodes_tmp = 0
 allocate( i_te_tmp      (2,n_p ) ) ;  i_te_tmp = 0
 allocate(ii_te_tmp      (2,n_el) ) ; ii_te_tmp = 0
 allocate(rr_te_tmp      (3,n_p ) ) ; rr_te_tmp = 0.0d0

 ! initialise counters ------
 ne_te = 0 ; nn_te = 0

 ! compute normals ---------- TODO: move out of this routine ---
 allocate(nor(3,n_el)) ; nor = 0.0_wp
 allocate(cen(3,n_el)) ; cen = 0.0_wp
 allocate(area( n_el)) ; area= 0.0_wp
 do i_e = 1 , n_el

   nSides = count( ee(:,i_e) .ne. 0 )

   nor(:,i_e) = 0.5_wp * cross( rr(:,ee(     3,i_e)) - rr(:,ee(1,i_e)) , & 
                                rr(:,ee(nSides,i_e)) - rr(:,ee(2,i_e))     )
   area( i_e) = norm2(nor(:,i_e))
   nor(:,i_e) = nor(:,i_e) / area(i_e) 
  
   do i1 = 1 , nSides
     cen(:,i_e) = cen(:,i_e) + rr(:,ee(i1,i_e))
   end do
   cen(:,i_e) = cen(:,i_e) / dble(nSides)

 end do
 
 do i_e = 1 , n_el 

   nSides = count( ee(:,i_e) .ne. 0 )

   do i_b = 1 , nSides
   
     if ( ( neigh_m(i_b,i_e) .ne. 0 ) .and. &
          ( all( e_te_tmp(1,1:ne_te) .ne. neigh_m(i_b,i_e) ) ) ) then
     
       ! Check normals ----
       ! 1. 
       if ( sum( nor(:,i_e) * nor(:,neigh_m(i_b,i_e)) ) .le. &
                                                       inner_prod_thresh ) then

         ne_te = ne_te + 1
         
         ! 2. ... other criteria to find te elements -------
         
         ! surface elements at the trailing edge -----------
         e_te_tmp(1,ne_te) = i_e
         e_te_tmp(2,ne_te) = neigh_m(i_b,i_e)

         ! nodes on the trailing edge ----
         i_el_nodes_tmp(1,ne_te) = ee( i_b , i_e )
         i_el_nodes_tmp(2,ne_te) = ee( mod(i_b,nSides) + 1 , i_e ) 

         ! Find the corresponding node on the other side 
         ! of the traling edge (for open te)
         ind1 = 1 ; mi1 = norm2( rr(:,i_el_nodes_tmp(1,ne_te)) -  &
                                                 rr(:,ee(1,neigh_m(i_b,i_e))) )
         ind2 = 1 ; mi2 = norm2( rr(:,i_el_nodes_tmp(2,ne_te)) - &
                                                 rr(:,ee(1,neigh_m(i_b,i_e))) )

         nSides_b = count( ee(:,neigh_m(i_b,i_e)) .ne. 0 )

         do i1 = 2 , nSides_b
           if (     norm2( rr(:,i_el_nodes_tmp(1,ne_te)) - &
                               rr(:,ee(i1,neigh_m(i_b,i_e)) ) ) .lt. mi1 ) then
             ind1 = i1
             mi1  = norm2( rr(:,i_el_nodes_tmp(1,ne_te)) - &
                                               rr(:,ee(i1,neigh_m(i_b,i_e)) ) )
           end if
           if (     norm2( rr(:,i_el_nodes_tmp(2,ne_te)) - &
                               rr(:,ee(i1,neigh_m(i_b,i_e)) ) ) .lt. mi2 ) then
             ind2 = i1
             mi2  = norm2( rr(:,i_el_nodes_tmp(2,ne_te)) - &
                                               rr(:,ee(i1,neigh_m(i_b,i_e)) ) )
           end if
         end do
         
         i_node1     = min(i_el_nodes_tmp(1,ne_te), ee(ind1,neigh_m(i_b,i_e)) )  
         i_node2     = min(i_el_nodes_tmp(2,ne_te), ee(ind2,neigh_m(i_b,i_e)) )  
         i_node1_max = max(i_el_nodes_tmp(1,ne_te), ee(ind1,neigh_m(i_b,i_e)) )  
         i_node2_max = max(i_el_nodes_tmp(2,ne_te), ee(ind2,neigh_m(i_b,i_e)) )  

         if ( all( i_te_tmp(1,1:nn_te) .ne. i_node2 ) ) then
           nn_te = nn_te + 1
           i_te_tmp(1,nn_te) = i_node2
           i_te_tmp(2,nn_te) = i_node2_max
           rr_te_tmp(:,nn_te) = 0.5d0 * ( & 
                                 rr(:,i_node2) + rr(:,i_node2_max) )
           ii_te_tmp(1,ne_te) = nn_te
         else 
           do i1 = 1 , nn_te   ! find ...
             if ( i_te_tmp(1,i1) .eq. i_node2 ) then 
               ii_te_tmp(1,ne_te) = i1 ; exit
             end if
           end do
         end if 
         if ( all( i_te_tmp(1,1:nn_te) .ne. i_node1 ) ) then
           nn_te = nn_te + 1
           i_te_tmp(1,nn_te) = i_node1
           i_te_tmp(2,nn_te) = i_node1_max
           rr_te_tmp(:,nn_te) = 0.5d0 * ( & 
                                 rr(:,i_node1) + rr(:,i_node1_max) )
           ii_te_tmp(2,ne_te) = nn_te
         else
           do i1 = 1 , nn_te   ! find ...
             if ( i_te_tmp(1,i1) .eq. i_node1 ) then 
               ii_te_tmp(2,ne_te) = i1 ; exit
             end if
           end do
         end if 
 
       end if

     end if  
 
   end do

 end do

 ! From tmp to actual e_te array --------------------
 if (allocated(e_te)) deallocate(e_te)
 allocate(e_te(2,ne_te)) ; e_te = e_te_tmp(:,1:ne_te)
 
 ! From tmp to actual e_te array --------------------
 if (allocated(i_te)) deallocate(i_te)
 allocate(i_te(2,nn_te)) ; i_te = i_te_tmp(:,1:nn_te)
 
 ! From tmp to actual r_te array --------------------
 if (allocated(rr_te))  deallocate(rr_te)
 allocate(rr_te(3,nn_te)) ; rr_te = rr_te_tmp(:,1:nn_te)
 
 ! From tmp to actual r_te array --------------------
 if (allocated(ii_te))  deallocate(ii_te)
 allocate(ii_te(2,ne_te)) ; ii_te = ii_te_tmp(:,1:ne_te)

 
 write(*,*) ' n.edges at the te : ' , ne_te
 write(*,*) ' n.nodes at the te : ' , nn_te

 deallocate( e_te_tmp, i_te_tmp, rr_te_tmp, ii_te_tmp )

 ! Trailing edge elements connectivity --------------
 
 allocate(neigh_te(2,ne_te) ) ; neigh_te = 0
 allocate(    o_te(2,ne_te) ) ;     o_te = 0
 
 do e1 = 1 , ne_te
  nSides1 = 2
  do i1 = 1 , nSides1
   do e2 = e1+1 , ne_te
    nSides2 = 2
    do i2 = 1 , nSides2
     if ( ii_te(i1,e1) .eq. ii_te(i2,e2) ) then
      neigh_te(i1,e1) = e2
      neigh_te(i2,e2) = e1
      if ( i1 .ne. i2 ) then
       o_te(i1,e1) = 1
       o_te(i2,e2) = 1
      else
       o_te(i1,e1) = -1
       o_te(i2,e2) = -1
      end if 
     end if
    end do    
   end do
  end do
 end do

 ! Compute unit vector at the te nodes ----------------
 !  for the first implicit wake panel  ----------------
 allocate(t_te  (3,nn_te)) ; t_te = 0.0_wp  
 do i_n = 1 , nn_te
!  vec1 = 0.0_wp
   t_te_len = 0.0_wp
   t_te_nelem = 0
   do i_e = 1 , ne_te 
     if ( any( i_n .eq. ii_te(:,i_e)) ) then
       t_te(:,i_n) =  t_te(:,i_n) + nor(:,e_te(1,i_e)) &
                                  + nor(:,e_te(2,i_e))

! **** mod 2018-07-04: avoid projection on v_rel direction.
! **** - avoid using hard-coded or params of the solver
! **** - for simulating fixed elements w/o free-stream velocity 

! **** ! <<<<<<<<<<
! **** ! TODO: refine the definition of the vec1
! **** vec1 = vec1 + cross(nor(:,e_te(1,i_e)) , nor(:,e_te(2,i_e)) )

! **** ! >>>>>>>>>>
       t_te_len = t_te_len + sqrt(area(e_te(1,i_e))) + sqrt(area(e_te(2,i_e)))
       t_te_nelem = t_te_nelem + 2
! **** ! >>>>>>>>>>

     end if
   end do

! tests: projection tests
!  ! 1. projection along the 'relative velocity'
!  t_te(:,i_n) = sum(t_te(:,i_n)*v_rel) * v_rel
!  ! 2. projection perpendicular to the side_slip direction
!  t_te(:,i_n) = t_te(:,i_n) - sum(t_te(:,i_n)*side_dir) * side_dir

! **** ! <<<<<<<<<<
! **** vec1 = vec1 - sum(vec1*v_rel) * v_rel
! **** vec1 = vec1 / norm2(vec1)
! ****
! **** ! projection ****
! **** write(*,*) ' vec1 ' , vec1
! **** t_te(:,i_n) = t_te(:,i_n) - sum(t_te(:,i_n)*vec1) * vec1
! **** ! <<<<<<<<<<

   ! normalisation
   t_te(:,i_n) = t_te(:,i_n) / norm2(t_te(:,i_n)) ! TODO: check if it is right
! **** ! >>>>>>>>>>
   t_te_len = 1.0_wp    ! 10.0_wp TODO: clean
   t_te(:,i_n) = t_te(:,i_n) * t_te_len ! / dble(t_te_nelem)
! **** ! >>>>>>>>>>

   if ( trim(te_proj_dir) .eq. 'normal' ) then
   
     side_dir = te_proj_vec
     if ( te_proj_logical ) then
       t_te(:,i_n) = t_te(:,i_n) - sum(t_te(:,i_n)*side_dir) * side_dir
     end if

     !debug
     write(*,*) ' proj_te normal : ' , t_te(:,i_n)

   elseif ( trim(te_proj_dir) .eq. 'parallel' ) then
   
     side_dir = te_proj_vec
     if ( te_proj_logical ) then
       t_te(:,i_n) = sum(t_te(:,i_n)*side_dir) * side_dir
     end if

   end if

 end do


! WARNING: in mod_geo.f90 the inverse transformation was introduced in the 
!  definition of the streamwise tangent unit vector at the te



end subroutine find_te_general

! -------------

end subroutine build_te_general

!----------------------------------------------------------------------

subroutine build_te_parametric ( ee , rr , ElType , &
                npoints_chord_tot , nelems_span , &
                e_te, i_te, rr_te, ii_te, neigh_te, o_te, t_te ) !te as an output
 
 integer , intent(in) :: ee(:,:)
 real(wp), intent(in) :: rr(:,:)
 character,intent(in) :: ElType
 integer , intent(in) :: npoints_chord_tot , nelems_span
 
 ! te structures
 integer , allocatable :: e_te(:,:) , i_te(:,:) , ii_te(:,:)
 integer , allocatable :: neigh_te(:,:) , o_te(:,:)
 real(wp), allocatable :: rr_te(:,:) , t_te(:,:)
 
 integer :: nelems_chord
 integer :: i1

 
 nelems_chord = npoints_chord_tot - 1
 
 ! e_te ----------
 allocate(e_te(2,nelems_span)) ; e_te = 0
 if ( ElType .eq. 'p' ) then
! <<<<<<<<<
!  e_te(1,:) = (/ ( 1+(i1-1)*nelems_chord , i1=1,nelems_span) /)   
!  e_te(2,:) = (/ (   (i1  )*nelems_chord , i1=1,nelems_span) /)
! >>>>>>>>>
   e_te(1,:) = (/ (   (i1  )*nelems_chord , i1=1,nelems_span) /)
   e_te(2,:) = (/ ( 1+(i1-1)*nelems_chord , i1=1,nelems_span) /)   
! >>>>>>>>>
 else
   e_te(1,:) = (/ (   (i1  )*nelems_chord , i1=1,nelems_span) /)
 ! e_te(2,:) = 0 
 end if
 
 ! i_te ----------
 allocate(i_te(2,nelems_span+1)) ; i_te = 0
 if ( ElType .eq. 'p' ) then
   i_te(:,1) = (/ ee(2,1) , ee(3,nelems_chord)/)
   i_te(1,2:nelems_span+1) = (/ ( ee(1,1+(i1-1)*nelems_chord) ,  &
                                                           i1=1,nelems_span) /)
   i_te(2,2:nelems_span+1) = (/ ( ee(4,  (i1  )*nelems_chord) , &
                                                           i1=1,nelems_span) /)
 else
   i_te(:,1) = (/ ee(3,nelems_chord) , ee(3,nelems_chord)/)
   i_te(1,2:nelems_span+1) = (/ ( ee(4,  (i1  )*nelems_chord) , &
                                                           i1=1,nelems_span) /)
   i_te(2,2:nelems_span+1) = (/ ( ee(4,  (i1  )*nelems_chord) , &
                                                           i1=1,nelems_span) /)
 end if
 
 ! rr_te ---------
 allocate(rr_te(3,nelems_span+1)) ; rr_te = 0.0_wp
 rr_te = 0.5_wp * ( rr(:,i_te(1,:)) + rr(:,i_te(2,:)) )
 
 ! ii_te ---------
 allocate(ii_te(2,nelems_span)) ; ii_te = 0
 if ( ElType .eq. 'p' ) then 
! <<<<<<<<<
!  ii_te(1,:) = (/ ( i1 , i1 = 1,nelems_span   ) /)
!  ii_te(2,:) = (/ ( i1 , i1 = 2,nelems_span+1 ) /)
! >>>>>>>>>
   ii_te(1,:) = (/ ( i1 , i1 = 2,nelems_span+1 ) /)
   ii_te(2,:) = (/ ( i1 , i1 = 1,nelems_span   ) /)
! >>>>>>>>>
 elseif ( ElType .eq. 'v' .or. ElType .eq. 'l') then 
   ii_te(1,:) = (/ ( i1 , i1 = 2,nelems_span+1 ) /)
   ii_te(2,:) = (/ ( i1 , i1 = 1,nelems_span   ) /)
 end if
 
 ! neigh_te ------
 allocate(neigh_te(2,nelems_span)) ; neigh_te = 0
 neigh_te(1,:) = (/ ( i1+1 , i1 = 1,nelems_span ) /)
 neigh_te(2,:) = (/ ( i1-1 , i1 = 1,nelems_span ) /) 
 neigh_te(2,1) = 0 
 neigh_te(1,nelems_span) = 0
 
 ! o_te ----------
 allocate(o_te(2,nelems_span)) ; o_te = 1  
 o_te(2,1) = 0
 o_te(1,nelems_span) = 0
 
 
 ! t_te ----------
 allocate(t_te(3,nelems_span+1)) ; t_te = 0.0_wp
 if ( ElType .eq. 'p' ) then
! <<<<<<<<<
!  t_te(:,1) = 0.5 * ( rr(:,ee(2,e_te(1,1))) - rr(:,ee(3,e_te(1,1))) + & 
!                      rr(:,ee(3,e_te(2,1))) - rr(:,ee(2,e_te(2,1)))     )
! >>>>>>>>>
   t_te(:,1) = 0.5 * ( rr(:,ee(3,e_te(1,1))) - rr(:,ee(2,e_te(1,1))) + & 
                       rr(:,ee(2,e_te(2,1))) - rr(:,ee(3,e_te(2,1)))     )
   t_te(:,1) = t_te(:,1) / norm2(t_te(:,1))
! >>>>>>>>>
   do i1 = 1 , nelems_span 
! <<<<<<<<<
!    t_te(:,i1+1) = 0.5 * ( rr(:,ee(1,e_te(1,i1))) - rr(:,ee(4,e_te(1,i1))) + & 
!                           rr(:,ee(4,e_te(2,i1))) - rr(:,ee(1,e_te(2,i1)))   )
! >>>>>>>>>
     t_te(:,i1+1) = 0.5 * ( rr(:,ee(4,e_te(1,i1))) - rr(:,ee(1,e_te(1,i1))) + & 
                            rr(:,ee(1,e_te(2,i1))) - rr(:,ee(4,e_te(2,i1)))  )
     t_te(:,i1+1) = t_te(:,i1+1) / norm2(t_te(:,i1+1))
   end do
 else
   t_te(:,1) = 0.5 * (  rr(:,ee(3,e_te(1,1))) - rr(:,ee(2,e_te(1,1))) )
   t_te(:,1) = t_te(:,1) / norm2(t_te(:,1))
   do i1 = 1 , nelems_span 
     t_te(:,i1+1) = 0.5 * ( rr(:,ee(4,e_te(1,i1))) - rr(:,ee(1,e_te(1,i1))) )
     t_te(:,i1+1) = t_te(:,i1+1) / norm2(t_te(:,i1+1))
   end do
 
 end if


end subroutine build_te_parametric

!----------------------------------------------------------------------
! subroutines :
!   create_local_velocity_stencil_...      ... = general, parametric
!   create_strip_connectivity_...





!----------------------------------------------------------------------



end module mod_build_geo
