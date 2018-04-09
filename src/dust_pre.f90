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

program dust_pre

use mod_param, only: &
  wp, nl, max_char_len, extended_char_len

use mod_handling, only: &
  error, warning, info, printout, dust_time, t_realtime

use mod_geometry, only: &
  t_geo, set_parameters_geo, create_geometry, update_geometry, &
  t_tedge,  destroy_geometry

use mod_basic_io, only: &
  read_mesh_basic, write_basic

use mod_aero_elements, only: &
  c_elem, t_elem_p!, t_vp

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
character(len=max_char_len), allocatable :: geo_files(:)
character(len=extended_char_len) :: message

!Time parameters
real(wp) :: tstart, tend, dt, time

!Geometry parameters
type(t_elem_p), allocatable :: elems(:)
type(t_parse) :: prms
type(t_parse), pointer :: sbprms, ssbprms
integer :: aa, bb, cc, dd, ee, ff
type(t_geo) :: geo
type(t_tedge) :: te


integer :: n_geo, i

!debug
integer :: n_refs,n,nmov, i_ref, intparam
character(len=max_char_len) :: Ref, mov_type
logical :: moving
type(t_link), pointer :: o_link

!write(*,*) 'DUST beginning'
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

call prms%CreateStringOption('GeoFile','Geometry definition files', multiple=.true.)
call prms%CreateStringOption('FileName','Preprocessor output file')

call prms%CreateStringOption('Ref','  ',multiple=.true.)
call prms%CreateLogicalOption('Moving','  ',multiple=.true.)

call prms%CreateSubOption('Movement',' ',sbprms,multiple=.true.)
call sbprms%CreateStringOption('Movement_Type','  ')
call sbprms%CreateIntOption('first_option','   ')
call sbprms%CreateIntOption('second_option','   ')

sbprms => null()

call prms%read_options(input_file_name, printout_val=.false.)

n_geo = countoption(prms,'GeoFile')
n_refs = countoption(prms,'Ref')

n = countoption(prms,'Moving')
if (n .ne. n_refs) call error('','','Wrong number of Moving inserted')
nmov = 0


do i_ref = 1,n_refs
  Ref = getstr(prms,'Ref')
    write(*,*) 'Reference: ',trim(Ref)
  moving = getlogical(prms,'Moving',olink=o_link)
    write(*,*) 'Is it moving: ',moving
  if(moving) then
    nmov = nmov + 1
    call check_opt_consistency(o_link, next=.true., next_opt='Movement')
    call getsuboption(prms,'Movement',sbprms) 

    mov_type =getstr(sbprms,'Movement_Type')
    select case(trim(mov_type))
    case('first')
      intparam = getint(sbprms,'first_option')
    case('second')
      intparam = getint(sbprms,'second_option')
    case default
      call error('','','Wrong kinda movement')
    end select
    write(*,*) 'Movement: ',trim(mov_type)
    write(*,*) 'Parameter: ',intparam
  endif
enddo

n = countoption(prms,'Movement')

if(n .ne. nmov) call error('','','Wrong number of movement with respect to Moving = T')

output_file_name_read = getstr(prms, 'FileName')

allocate(geo_files(n_geo))
do i=1,n_geo
  geo_files(i) = getstr(prms,'GeoFile')
enddo


n = 3
do i = 1,n

  n = n+1
  write(*,*) 'i',i
enddo 
if (.not. cmd_set_filename) output_file_name = output_file_name_read

call build_geometry(geo_files,trim(output_file_name))

! Read the file list from the main file
!call set_parameters_geo(prms)



! get the parameters and print them out
!call printout(new_line('a')//'====== Input parameters: ======')


!call printout(nl//'====== Geometry Creation ======')
!call create_geometry(prms, input_file_name, geo, te, elems, tstart)



call destroy_hdf5()
call printout(nl//'<<<<<< DUST PREPROCESSOR end       <<<<<<'//nl)

end program dust_pre
