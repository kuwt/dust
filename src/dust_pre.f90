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
!! Copyright (C) 2018-2020 Davide   Montagnani,
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

!> This is the main file of the DUST preprocessor
program dust_pre

use mod_param, only: &
  wp, nl, max_char_len, extended_char_len

use mod_handling, only: &
  error, warning, info, printout, dust_time, t_realtime, check_file_exists

use mod_geometry, only: &
  t_geo, create_geometry, update_geometry, &
  t_tedge,  destroy_geometry

use mod_basic_io, only: &
  read_mesh_basic, write_basic

!this is for the parsing
use mod_parse, only: &
  t_parse, &
  getstr, getlogical, getreal, getint, &
  ignoredParameters, finalizeParameters, &
  countoption, getsuboption, t_link, check_opt_consistency

use mod_build_geo, only: &
  build_geometry

use mod_hdf5_io, only: &
initialize_hdf5, destroy_hdf5

implicit none

!Input
character(len=*), parameter :: input_file_name_def = 'dust_pre.in'
character(len=max_char_len) :: input_file_name
character(len=max_char_len) :: output_file_name
character(len=max_char_len) :: output_file_name_read
logical :: cmd_set_filename
character(len=max_char_len), allocatable :: comp_names(:)
character(len=max_char_len), allocatable :: geo_files(:)
character(len=max_char_len), allocatable :: ref_tags(:)
!character(len=extended_char_len) :: message

!Geometry parameters
type(t_parse) :: prms
type(t_link), pointer :: lnk

!General parameters: tol_sew,inner_prod_thr
real(wp) :: tol_sew
real(wp) :: inner_prod_thr

integer :: n_geo, n_tag, n_name, i


call printout(nl//'>>>>>> DUST PREPROCESSOR beginning >>>>>>'//nl)
call initialize_hdf5()

!------ Input reading ------

if(command_argument_count().gt.0) then
  call get_command_argument(1,value=input_file_name)
else
  input_file_name = input_file_name_def
endif


if(command_argument_count().gt.1) then
  call get_command_argument(2,value=output_file_name)
  cmd_set_filename = .true.
else
  cmd_set_filename = .false.
endif


!call prms%CreateRealOption( 'tstart', "Starting time", '0.0')

call prms%CreateStringOption('CompName','Component Name', multiple=.true.)
call prms%CreateStringOption('GeoFile','Geometry definition files', multiple=.true.)
call prms%CreateStringOption('RefTag','Reference Tag of the component', multiple=.true.)

call prms%CreateRealOption('TolSewing','Global parameter for closing gaps','0.001')
call prms%CreateRealOption('InnerProductTe','Global parameter for edge identification','-0.5')

call prms%CreateStringOption('FileName','Preprocessor output file')

call check_file_exists(input_file_name,'dust preprocessor')
call prms%read_options(input_file_name, printout_val=.false.)

n_name = countoption(prms,'CompName')
n_geo  = countoption(prms,'GeoFile')
n_tag  = countoption(prms,'RefTag')

tol_sew = getreal(prms,'TolSewing')
inner_prod_thr = getreal(prms,'InnerProductTe')

if(n_geo .ne. n_name)  call error('dust_pre','','Different number of components &
  & and components names in input file "'//trim(input_file_name)//'"')
if(n_geo .ne. n_tag)  call error('dust_pre','','Different number of components &
  & and references tags in input file "'//trim(input_file_name)//'"')

output_file_name_read = getstr(prms, 'FileName')

allocate(comp_names(n_geo), geo_files(n_geo), ref_tags(n_geo))
do i=1,n_geo
  comp_names(i) = getstr(prms,'CompName')
  geo_files(i) = getstr(prms,'GeoFile',olink=lnk)
  call check_opt_consistency(lnk,prev=.true.,prev_opt='CompName')
  call check_opt_consistency(lnk,next=.true.,next_opt='RefTag')
  ref_tags(i)  = getstr(prms,'RefTag')
enddo

if (.not. cmd_set_filename) output_file_name = output_file_name_read

call build_geometry(geo_files,ref_tags,comp_names,trim(output_file_name), &
                 tol_sew , inner_prod_thr )


deallocate(comp_names, geo_files, ref_tags)
call destroy_hdf5()
call printout(nl//'<<<<<< DUST PREPROCESSOR end       <<<<<<'//nl)

end program dust_pre
