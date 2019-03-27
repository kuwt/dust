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
!! Copyright (C) 2018-2019 Davide   Montagnani, 
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
!!          Federico Fonte             <federico.fonte@outlook.com>
!!          Davide Montagnani       <davide.montagnani@gmail.com>
!!          Matteo Tugnoli                <tugnoli.teo@gmail.com>
!!=========================================================================

!> This is the main file of the DUST postprocessor
program dust_post

use mod_param, only: &
  wp, nl, max_char_len, extended_char_len , pi

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
  post_integral

use mod_post_sectional, only: &
  post_sectional

implicit none

!Input
character(len=*), parameter :: input_file_name_def = 'dust_post.in' 
character(len=max_char_len) :: input_file_name

!doublet parameters
real(wp) :: ff_ratio_dou, ff_ratio_sou, eps_dou, r_Rankine, r_Vortex, r_cutoff

!Geometry parameters
type(t_parse) :: prms
type(t_parse), pointer :: sbprms , bxprms

integer :: n_analyses, ia
character(len=max_char_len) :: basename, data_basename
character(len=max_char_len) :: an_name, an_type
integer :: an_start, an_end, an_step
integer :: n_comp, i_comp
character(len=max_char_len), allocatable :: components_names(:)
character(len=max_char_len), allocatable :: var_names(:)
character(len=max_char_len) :: lowstr
character(len=max_char_len) :: out_frmt
logical :: all_comp
logical :: average

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

call prms%CreateRealOption( 'FarFieldRatioDoublet', &
      "Multiplier for far field threshold computation on doublet", '10.0')
call prms%CreateRealOption( 'FarFieldRatioSource', &
      "Multiplier for far field threshold computation on sources", '10.0')
call prms%CreateRealOption( 'DoubletThreshold', &
      "Thresold for considering the point in plane in doublets", '1.0e-6')
call prms%CreateRealOption( 'RankineRad', &
      "Radius of Rankine correction for vortex induction near core", '0.1')
call prms%CreateRealOption( 'VortexRad', &
      "Radius of vortex core, for particles", '0.1')
call prms%CreateRealOption( 'CutoffRad', &
      "Radius of complete cutoff  for vortex induction near core", '0.001')

call prms%CreateSubOption('Analysis','Definition of the motion of a frame', &
                          sbprms, multiple=.true.)
call sbprms%CreateStringOption('Type','type of analysis')
call sbprms%CreateStringOption('Name','specification of the analysis')
call sbprms%CreateIntOption('StartRes', 'Starting result of the analysis')
call sbprms%CreateIntOption('EndRes', 'Final result of the analysis')
call sbprms%CreateIntOption('StepRes', 'Result stride of the analysis')
call sbprms%CreateLogicalOption('Average', 'Perform time averaging','F')
call sbprms%CreateLogicalOption('Wake', 'Output also the wake for &
                                &visualization','T')
call sbprms%CreateLogicalOption('SeparateWake', 'Output the wake in a separate &
                                &way','F')
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
! sectional loads: box -----
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
call bxprms%CreateLogicalOption('reshapeBox','logical input to reshape the box if &
                           &it is "too large"')
call sbprms%CreateRealArrayOption('AxisMom','axis for the computation of the moment. Perpendicular to sections')

sbprms=>null()

call check_file_exists(input_file_name, 'dust postprocessor')
call prms%read_options(input_file_name, printout_val=.false.)

basename = getstr(prms,'basename')
data_basename = getstr(prms,'data_basename')

ff_ratio_dou  = getreal(prms, 'FarFieldRatioDoublet')
ff_ratio_sou  = getreal(prms, 'FarFieldRatioSource')
eps_dou   = getreal(prms, 'DoubletThreshold')
r_Rankine = getreal(prms, 'RankineRad')
r_Vortex  = getreal(prms, 'VortexRad')
r_cutoff  = getreal(prms, 'CutoffRad')

call initialize_doublet(ff_ratio_dou, eps_dou, r_Rankine, r_cutoff);
call initialize_surfpan(ff_ratio_sou)
call initialize_vortline(r_Rankine, r_cutoff)
call initialize_vortpart(r_Vortex, r_cutoff)

n_analyses = countoption(prms,'Analysis')

!Check that the basenames are valid
call check_basename(trim(basename),'dust postprocessor')
call check_basename(trim(data_basename),'dust postprocessor')

!Cycle on all the analyses
do ia = 1,n_analyses
  
  !Get general parameter for the analysis
  call getsuboption(prms,'Analysis',sbprms)
  an_type  = getstr(sbprms,'Type')
  call LowCase(an_type)
  an_name  = getstr(sbprms,'Name')
  out_frmt = getstr(sbprms,'Format')
  call LowCase(out_frmt)
  an_start = getint(sbprms,'StartRes')
  an_end   = getint(sbprms,'EndRes')
  an_step  = getint(sbprms,'StepRes')
  average  = getlogical(sbprms, 'Average')
 
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
      n_comp = 0
      deallocate(components_names)
    endif
  endif

  !Fork the different kind of analyses
  select case(trim(an_type))

   !//////////////// Integral Loads  \\\\\\\\\\\\\\\\
   case('integral_loads')

    ! look in mod_post_integral
    call post_integral( sbprms , basename , data_basename , an_name , ia , &
                        out_frmt , components_names , all_comp , &
                        an_start , an_end , an_step, average)

   !//////////////////Visualizations\\\\\\\\\\\\\\\\\
   case('viz') 

    call post_viz( sbprms , basename , data_basename , an_name , ia , &
                   out_frmt , components_names , all_comp , &
                   an_start , an_end , an_step, average )

   !//////////////////Domain probes \\\\\\\\\\\\\\\\\
   case('probes')

    call post_probes( sbprms , basename , data_basename , an_name , ia , &
                      out_frmt , components_names , all_comp , &
                      an_start , an_end , an_step )

   !////////////////// Flow Field  \\\\\\\\\\\\\\\\\\
   case('flow_field')

    call post_flowfield ( sbprms , basename , data_basename , an_name , ia , &
                          out_frmt , components_names , all_comp , &
                          an_start , an_end , an_step, average)

   !///////////////  Sectional Loads  \\\\\\\\\\\\\\\
   case('sectional_loads')
    
    call post_sectional ( sbprms , bxprms , basename , data_basename , an_name , ia , &
                          out_frmt , components_names , all_comp , &
                          an_start , an_end , an_step, average )

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
! contains
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!
! moved to ./post/mod_post_load.f90 
! 
! subroutine load_refs(floc, refs_R, refs_off, refs_G, refs_f, refs_tag)
! subroutine load_res(floc, comps, vort, press, t)
! subroutine load_wake_pan(floc, wpoints, wstart, wvort)
! subroutine load_wake_ring(floc, wpoints, wconn, wvort)
! subroutine load_wake_viz(floc, wpoints, welems, wvort)
! 
!----------------------------------------------------------------------

end program dust_post
