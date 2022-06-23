!./\\\\\\\\\\\...../\\\......./\\\..../\\\\\\\\\..../\\\\\\\\\\\\\.
!.\/\\\///////\\\..\/\\\......\/\\\../\\\///////\\\.\//////\\\////..
!..\/\\\.....\//\\\.\/\\\......\/\\\.\//\\\....\///.......\/\\\......
!...\/\\\......\/\\\.\/\\\......\/\\\..\////\\.............\/\\\......
!....\/\\\......\/\\\.\/\\\......\/\\\.....\///\\...........\/\\\......
!.....\/\\\......\/\\\.\/\\\......\/\\\.......\///\\\........\/\\\......
!......\/\\\....../\\\..\//\\\...../\\\../\\\....\//\\\.......\/\\\......
!.......\/\\\\\\\\\\\/....\///\\\\\\\\/..\///\\\\\\\\\/........\/\\\......
!........\///////////........\////////......\/////////..........\///.......
!!=========================================================================
!!
!! Copyright (C) 2018-2022 Politecnico di Milano,
!!                           with support from A^3 from Airbus
!!                    and  Davide   Montagnani,
!!                         Matteo   Tugnoli,
!!                         Federico Fonte
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
!!          Federico Fonte
!!          Davide Montagnani
!!          Matteo Tugnoli
!!=========================================================================

!> This is the main file of the DUST postprocessor
program dust_post

use mod_param, only: &
  wp, nl, max_char_len, extended_char_len , pi

use mod_sim_param, only: &
  sim_param

use mod_handling, only: &
  error, warning, info, printout, dust_time, t_realtime, new_file_unit, &
  check_basename, check_file_exists

use mod_geometry, only: &
  t_geo, t_geo_component

use mod_basic_io, only: &
  read_mesh_basic, write_basic

use mod_aeroel, only: &
  c_elem, c_pot_elem, c_vort_elem, c_impl_elem, c_expl_elem, &
  t_elem_p, t_pot_elem_p, t_vort_elem_p, t_impl_elem_p, t_expl_elem_p

use mod_doublet, only: &
  initialize_doublet

use mod_surfpan, only: &
  initialize_surfpan

use mod_vortline, only: &
  initialize_vortline

use mod_vortpart, only: &
  initialize_vortpart

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
  LowCase, isInList, stricmp

use mod_geo_postpro, only: &
  load_components_postpro, update_points_postpro , prepare_geometry_postpro, &
  expand_actdisk_postpro

use mod_tecplot_out, only: &
  tec_out_viz, tec_out_probes, tec_out_box, tec_out_loads

use mod_vtk_out, only: &
  vtk_out_viz , vtr_write

use mod_dat_out, only: &
  dat_out_probes_header, &
  dat_out_loads_header

use mod_math, only: &
  cross

use mod_actuatordisk, only: &
  t_actdisk

! "unravelled" analysis
use mod_post_load, only: &
  load_refs , load_res , load_wake_pan , load_wake_ring

use mod_post_viz, only: &
  post_viz

use mod_post_probes, only: &
  post_probes

use mod_post_flowfield, only: &
  post_flowfield

use mod_post_integral, only: &
  post_integral, post_hinge_loads

use mod_post_sectional, only: &
  post_sectional

use mod_post_chordwise, only: & 
  post_chordwise 

use mod_post_aa, only: &
  post_aeroacoustics

implicit none

!Input
character(len=*), parameter               :: input_file_name_def = 'dust_post.in'
character(len=max_char_len)               :: input_file_name

!Geometry parameters
type(t_parse)                             :: prms
type(t_parse), pointer                    :: sbprms , bxprms
integer                                   :: n_analyses, ia             
character(len=max_char_len)               :: basename, data_basename
character(len=max_char_len)               :: an_name, an_type
integer                                   :: an_start, an_end, an_step
integer                                   :: n_comp, i_comp 
integer                                   :: n_hinge, i_hinge
character(len=max_char_len), allocatable  :: hinge_tag(:)
character(len=max_char_len), allocatable  :: components_names(:)
character(len=max_char_len), allocatable  :: var_names(:)
character(len=max_char_len)               :: lowstr
character(len=max_char_len)               :: out_frmt
logical                                   :: all_comp
logical                                   :: all_hinge
logical                                   :: average

call printout(nl//'>>>>>> DUST POSTPROCESSOR beginning >>>>>>'//nl)
call initialize_hdf5()

!------ Input reading ------

if(command_argument_count().gt.0) then
  call get_command_argument(1,value=input_file_name)
else
  input_file_name = input_file_name_def
endif

call printout(nl//'Reading input parameters from file "'//&
                trim(input_file_name)//'"'//nl)

call prms%CreateStringOption('basename','Base name of the processed data')
call prms%CreateStringOption('data_basename','Base name of the data to be &
                              &processed')

call prms%CreateRealOption( 'far_field_ratio_doublet', &
      "Multiplier for far field threshold computation on doublet", '10.0')
call prms%CreateRealOption( 'far_field_ratio_source', &
      "Multiplier for far field threshold computation on sources", '10.0')
call prms%CreateRealOption( 'doublet_threshold', &
      "Thresold for considering the point in plane in doublets", '1.0e-6')
call prms%CreateRealOption( 'rankine_rad', &
      "Radius of Rankine correction for vortex induction near core", '0.1')
call prms%CreateRealOption( 'vortex_rad', &
      "Radius of vortex core, for particles", '0.1')
call prms%CreateRealOption( 'cutoff_rad', &
      "Radius of complete cutoff  for vortex induction near core", '0.001')

call prms%CreateSubOption('analysis','Definition of the motion of a frame', &
                          sbprms, multiple=.true.)
call sbprms%CreateStringOption('type','type of analysis')
call sbprms%CreateStringOption('name','specification of the analysis')
call sbprms%CreateIntOption('start_res', 'Starting result of the analysis')
call sbprms%CreateIntOption('end_res', 'Final result of the analysis')
call sbprms%CreateIntOption('step_res', 'Result stride of the analysis')
call sbprms%CreateLogicalOption('average', 'Perform time averaging','F')
call sbprms%CreateLogicalOption('wake', 'Output also the wake for &
                                &visualization','T')
call sbprms%CreateLogicalOption('separate_wake', 'Output the wake in a separate &
                                &way','F')
call sbprms%CreateStringOption('format','Output format')
call sbprms%CreateStringOption('component','Component to analyse', &
                              multiple=.true.)
call sbprms%CreateStringOption('hinge_tag','Hinge to analyse', &
                              multiple=.true.)                              
call sbprms%CreateStringOption('variable','Variables to be saved: velocity, pressure or&
                              & vorticity', multiple=.true.)

! probe output -------------
call sbprms%CreateStringOption('input_type','How to specify probe coordinates',&
                              multiple=.true.)
call sbprms%CreateRealArrayOption('point','Point coordinates in dust_post.in',&
                              multiple=.true.)
call sbprms%CreateStringOption('file','File containing the coordinates of the probes',&
                              multiple=.true.)
! flow field output --------
call sbprms%CreateIntArrayOption( 'n_xyz','number of points per coordinate',&
                              multiple=.true.)
call sbprms%CreateRealArrayOption('min_xyz','lower bounds of the box',&
                              multiple=.true.)
call sbprms%CreateRealArrayOption('max_xyz','upper bounds of the box',&
                              multiple=.true.)
! loads --------------------
call sbprms%CreateStringOption('comp_name','Components where loads are computed',&
                              multiple=.true.)
call sbprms%CreateStringOption('reference_tag','Reference frame where loads&
                            & are computed',multiple=.true.)
! sectional loads ----------
call sbprms%CreateRealArrayOption('axis_dir','Direction of the axis defined the reference&
                            & points for sectional loads analisys', multiple=.true.)
call sbprms%CreateRealArrayOption('axis_nod','Node belonging to the axis used for sectional&
                            & loads analisys', multiple=.true.)
call sbprms%CreateLogicalOption('lifting_line_data', 'Output lifting line data&
                            & alongside sectional loads','F')
call sbprms%CreateLogicalOption('vortex_lattice_data', 'Output of corrected vortex lattice data&
                            & alongside sectional loads','F')
! chordwise loads
call sbprms%CreateIntOption('n_station','Number of stations where the loads are extracted' &
                            , multiple=.true.)
call sbprms%CreateRealArrayOption('span_station','Spanwise coordinates in the component reference frame&
                            & where they are extracted', multiple=.true.)


! sectional loads: box -----
call sbprms%CreateSubOption('box_sect','Definition of the box for sectional loads', &
                          bxprms)
call bxprms%CreateRealArrayOption('ref_node','reference node to build the box')
call bxprms%CreateRealArrayOption('face_vec','vector identifying the direction of the base side &
                          &of the sections')
call bxprms%CreateRealArrayOption('face_bas','dimension along faceVec of the first and last sections')
call bxprms%CreateRealArrayOption('face_hei','dimension orthogonal to faceVec of the first and last sections')
call bxprms%CreateRealArrayOption('span_vec','vector defining the out-of plane direction of the box')
call bxprms%CreateRealOption('span_len','dimension along the spanVec direction')
call bxprms%CreateIntOption('num_sect','number of sections')
call bxprms%CreateLogicalOption('reshape_box','logical input to reshape the box if &
                          &it is "too large"')
call sbprms%CreateRealArrayOption('axis_mom','axis for the computation of the moment. Perpendicular to sections')

sbprms=>null()

call check_file_exists(input_file_name, 'dust postprocessor')
call prms%read_options(input_file_name, printout_val=.false.)

basename = getstr(prms,'basename')
data_basename = getstr(prms,'data_basename')

sim_param%FarFieldRatioDoublet  = getreal(prms, 'far_field_ratio_doublet')
sim_param%FarFieldRatioSource  = getreal(prms, 'far_field_ratio_source')
sim_param%DoubletThreshold   = getreal(prms, 'doublet_threshold')
sim_param%RankineRad = getreal(prms, 'rankine_rad')
sim_param%VortexRad = getreal(prms, 'vortex_rad')
sim_param%CutoffRad  = getreal(prms, 'cutoff_rad')

call initialize_doublet()
call initialize_surfpan()
call initialize_vortline()

n_analyses = countoption(prms,'analysis')

!> Check that the basenames are valid
call check_basename(trim(basename),'dust postprocessor')
call check_basename(trim(data_basename),'dust postprocessor')

!> Cycle on all the analyses
do ia = 1, n_analyses

  !> Get general parameter for the analysis
  call getsuboption(prms,'Analysis',sbprms)
  an_type  = getstr(sbprms,'type')
  
  call LowCase(an_type)
  an_name  = getstr(sbprms,'name')
  out_frmt = getstr(sbprms,'format')
  
  call LowCase(out_frmt)
  an_start = getint(sbprms,'start_res')
  an_end   = getint(sbprms,'end_res')
  an_step  = getint(sbprms,'step_res')
  average  = getlogical(sbprms, 'average')

  !> Check if we are analysing all the components or just some
  all_comp = .false.
  n_comp = countoption(sbprms, 'component')

  if (n_comp .eq. 0)  then
    all_comp = .true.
  else
    allocate(components_names(n_comp))
    do i_comp = 1, n_comp
      components_names(i_comp) = getstr(sbprms, 'component')
    enddo
    call LowCase(components_names(1),lowstr)    ! char
    if(trim(lowstr) .eq. 'all') then
      all_comp = .true.
      n_comp = 0
      deallocate(components_names)
    endif
  endif

  if (trim(an_type) .eq. 'hinge_loads') then 
    all_hinge = .false. 
    n_hinge = countoption(sbprms, 'hinge_tag')
    if (n_hinge .eq. 0)  then
      all_hinge = .true.
    else
      allocate(hinge_tag(n_hinge))
      do i_hinge = 1, n_hinge
        hinge_tag(i_hinge) = getstr(sbprms, 'hinge_tag')
      enddo
      call LowCase(hinge_tag(1), lowstr)    ! char
      if(trim(lowstr) .eq. 'all') then
        all_hinge = .true.
        n_hinge = 0
        deallocate(hinge_tag)
      endif
    endif
  endif 

  !> Fork the different kind of analyses
  select case( trim(an_type) )

    !> Integral Loads
    case('integral_loads')
      call post_integral( sbprms , basename , data_basename , an_name , ia , &
                          out_frmt , components_names , all_comp , &
                          an_start , an_end , an_step, average )

    !> Hinge Moment
    case('hinge_loads')
      call post_hinge_loads( sbprms , basename , data_basename , an_name , ia , &
                          out_frmt , components_names,  all_comp, hinge_tag, all_hinge,   &
                          an_start , an_end , an_step, average )

    !> Visualizations
    case('viz')
      call post_viz( sbprms , basename , data_basename , an_name , ia , &
                    out_frmt , components_names , all_comp , &
                    an_start , an_end , an_step, average )

    !> Domain probes
    case('probes')
      call post_probes( sbprms , basename , data_basename , an_name , ia , &
                      out_frmt , components_names , all_comp , &
                      an_start , an_end , an_step )

    !> Flow Field
    case('flow_field')
      call post_flowfield ( sbprms , basename , data_basename , an_name , ia , &
                          out_frmt , components_names , all_comp , &
                          an_start , an_end , an_step, average)

    !> Sectional Loads
    case('sectional_loads')
      call post_sectional ( sbprms , bxprms , basename , data_basename , an_name , ia , &
                          out_frmt , components_names , all_comp , &
                          an_start , an_end , an_step, average )

    !> Chordwise Loads 
    case('chordwise_loads') 
      call post_chordwise (sbprms, basename, data_basename, an_name, ia, &
                            out_frmt, components_names, all_comp, &
                            an_start, an_end, an_step, average)

    !> Aeroacoustics
    case('aeroacoustics')
      call post_aeroacoustics ( sbprms, basename, data_basename, &
                                an_name, ia, out_frmt, components_names, &
                                all_comp, an_start, an_end, an_step, average )

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

end program dust_post
