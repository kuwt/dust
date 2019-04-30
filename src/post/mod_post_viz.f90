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

!> Module containing the subroutines to perform visualizations during
!! postprocessing
module mod_post_viz

use mod_param, only: &
  wp, nl, max_char_len, extended_char_len , pi

use mod_handling, only: &
  error, warning, info, printout, dust_time, t_realtime, new_file_unit

use mod_parse, only: &
  t_parse, &
  getstr, getlogical, &
  countoption

use mod_hdf5_io, only: &
!  initialize_hdf5, destroy_hdf5, &
   h5loc, &
!  new_hdf5_file, &
   open_hdf5_file, &
   close_hdf5_file, & ! , &
!  new_hdf5_group, &
   open_hdf5_group, &
   close_hdf5_group, &
!  write_hdf5, &
   read_hdf5 
!  read_hdf5_al, &
!  check_dset_hdf5

use mod_stringtools, only: &
  LowCase, isInList
! LowCase, isInList, stricmp

use mod_geometry, only: &
  t_geo, t_geo_component, destroy_elements

use mod_geo_postpro, only: &
  load_components_postpro, update_points_postpro, prepare_geometry_postpro, &
  expand_actdisk_postpro

use mod_tecplot_out, only: &
  tec_out_viz ! , tec_out_probes, tec_out_box, tec_out_loads

use mod_vtk_out, only: &
  vtk_out_viz ! , vtr_write

use mod_post_load, only: &
  load_refs, load_res, load_wake_viz , check_if_components_exist

use mod_vtk_utils, only: &
  t_output_var, add_output_var, copy_output_vars, clear_output_vars

implicit none

public :: post_viz

private

character(len=*), parameter :: this_mod_name='mod_post_viz'
character(len=max_char_len) :: msg

contains

! ---------------------------------------------------------------------- 

subroutine post_viz( sbprms , basename , data_basename , an_name , ia , &
                     out_frmt , components_names , all_comp , &
                     an_start , an_end , an_step, average )
 type(t_parse), pointer :: sbprms
 character(len=*) , intent(in) :: basename
 character(len=*) , intent(in) :: data_basename
 character(len=*) , intent(in) :: an_name
 integer          , intent(in) :: ia
 character(len=*) , intent(in) :: out_frmt
 character(len=max_char_len), allocatable , intent(inout) :: components_names(:)
 logical , intent(in) :: all_comp
 integer , intent(in) :: an_start , an_end , an_step
 logical, intent(in) :: average
 
 type(t_geo_component), allocatable :: comps(:)
 character(len=max_char_len) :: filename
 integer(h5loc) :: floc , ploc
 logical :: out_vort, out_vel, out_cp, out_press , out_wake, separate_wake
 integer :: n_var , i_var
 character(len=max_char_len), allocatable :: var_names(:)
 real(wp), allocatable :: points(:,:), points_exp(:,:) , wpoints(:,:)
 real(wp), allocatable :: vppoints(:,:), vpvort(:)
 integer , allocatable :: elems(:,:) , welems(:,:)
 integer :: nelem , nelem_w, nelem_vp

 real(wp), allocatable :: points_ave(:,:), ave_vars(:,:)
 
 real(wp) :: u_inf(3)
 real(wp) :: P_inf , rho

 real(wp), allocatable :: refs_R(:,:,:), refs_off(:,:)
 real(wp), allocatable :: vort(:), cp(:), vel(:), press(:), wvort(:)

 real(wp), allocatable :: print_vars(:,:)
 character(len=max_char_len), allocatable :: print_var_names(:)
 character(len=max_char_len), allocatable :: ave_var_names(:)
 
 real(wp), allocatable :: print_vars_w(:,:)
 real(wp), allocatable :: print_vars_vp(:,:)
 character(len=max_char_len), allocatable :: print_var_names_w(:)

 type(t_output_var), allocatable :: out_vars(:), ave_out_vars(:)
 type(t_output_var), allocatable :: out_vars_w(:), out_vars_vp(:)
 integer :: nprint , nprint_w, nelem_out
 
 integer :: it, ires
 real(wp) :: t
 character(len=*), parameter :: this_sub_name='post_viz'

  write(msg,'(A,I0,A)') nl//'++++++++++ Analysis: ',ia,' visualization'//nl
  call printout(trim(msg))
  
  ! Print the wake or not 
  out_wake = getlogical(sbprms,'Wake')
  separate_wake = getlogical(sbprms,'SeparateWake')

  !Check which variables to analyse
  out_vort = .false.; out_vel = .false.; out_press =.false.; out_cp = .false.
  n_var = countoption(sbprms, 'Variable')
  allocate(var_names(n_var))
  do i_var = 1, n_var 
    var_names(i_var) = getstr(sbprms, 'Variable') ; call LowCase(var_names(i_var))
  enddo
  out_vort = isInList('vorticity',var_names) ! Always lower case string in the code !
  out_vel  = isInList('velocity' ,var_names)
  out_press= isInList('pressure' ,var_names)
  out_cp   = isInList('cp'       ,var_names)
  nprint = 0; nprint_w = 0
  if(out_vort)  nprint = nprint+1
  if(out_cp)    nprint = nprint+1
  if(out_vel)   nprint = nprint+1  !<--- *** TODO ***
  if(out_press) nprint = nprint+1  !<--- *** TODO ***
  allocate(out_vars(nprint))
  if(average) allocate(ave_out_vars(nprint))
  !for the wake 
  if(out_wake) then
    allocate(out_vars_w(nprint))
    allocate(out_vars_vp(nprint))
  endif
  
  
  ! Load the components (just once)
  call open_hdf5_file(trim(data_basename)//'_geo.h5', floc)
  
  
  call load_components_postpro(comps, points, nelem, floc, &
                               components_names, all_comp)
  
  call close_hdf5_file(floc)
  
  
  if(out_wake .and. average) call error(this_sub_name, this_mod_name, &
  'Cannot output an averaged wake visualization. Remove the wake or avoid &
  &averaging')

  ! Prepare_geometry_postpro
  call prepare_geometry_postpro(comps)
  
  ! Time loop
  ires = 0
  do it = an_start, an_end, an_step
    ires = ires+1
  
    ! Open the file 
    write(filename,'(A,I4.4,A)') trim(data_basename)//'_res_',it,'.h5'
    call open_hdf5_file(trim(filename),floc)
  
    ! Load free-stream parameters
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
    if(average) then
      if( .not. allocated(points_ave)) then
        allocate(points_ave(size(points_exp,1),size(points_exp,2)))
        points_ave = 0.0_wp
      endif
      points_ave = points_ave*(real(ires-1,wp)/real(ires,wp)) + &
                   points_exp/real(ires,wp)
    endif
  
    !Load the results ! TODO: check this routine and the content of the files to be read
    ! TODO : compute the missing quantities
    call load_res(floc, comps, vort, press, t)
  
    !Prepare the variable for output
    nelem_out = size(vort)
    ! TODO: compute/or read pressure and velocity field. Now set equal to zero
    allocate(  vel(size( vort,1)) ) ; vel = 0.0_wp
    allocate(   cp(size(press,1)) ) ;  cp = 0.0_wp
  
    !allocate(print_var_names(nprint), print_vars(nelem_out, nprint))
    
    i_var = 1
    if(out_vort) then
      !print_vars(:,i_var) = vort
      !print_var_names(i_var) = 'Vorticity'
      call add_output_var(out_vars(i_var), vort, 'Singularity_Intensity', &
                          .false.)
      i_var = i_var +1
    endif
    if(out_cp) then
      !print_vars(:,i_var) = cp
      !print_var_names(i_var) = 'Cp'
      call add_output_var(out_vars(i_var), cp, 'Cp',.false.)
      i_var = i_var +1
    endif
    if(out_vel) then
      !print_vars(:,i_var) = vel
      !print_var_names(i_var) = 'Velocity'
      call add_output_var(out_vars(i_var), vel, 'Velocity',.false.)
      i_var = i_var +1
    endif
    if(out_press) then
      !print_vars(:,i_var) = press
      !print_var_names(i_var) = 'Pressure'
      call add_output_var(out_vars(i_var), press, 'Pressure',.false.)
      i_var = i_var +1
    endif
  
    if(average) then
      if( ires .eq. 1) then
        call copy_output_vars(out_vars, ave_out_vars, .true.)
      endif
      !if( .not. allocated(ave_vars)) then
      !  allocate(ave_vars(size(print_vars,1),size(print_vars,2)))
      !  ave_vars = 0.0_wp
      !endif
      !if( .not. allocated(ave_var_names)) then
      !  allocate(ave_var_names(size(print_var_names,1)))
      !  ave_var_names = print_var_names
      !endif
      do i_var = 1,size(ave_out_vars)
      ave_out_vars(i_var)%var = ave_out_vars(i_var)%var * &
                                (real(ires-1,wp)/real(ires,wp)) + &
                                out_vars(i_var)%var/real(ires,wp)
      enddo
    endif
    
    if(.not. average) then
      ! Output filename
      write(filename,'(A,I4.4)') trim(basename)//'_'//trim(an_name)//'-',it
      
      if (out_wake) then
        
        call load_wake_viz(floc, wpoints, welems, wvort, vppoints, vpvort)
        nelem_w = size(welems,2)
        nelem_vp = size(vppoints,2)
        
        i_var = 1
        if(out_vort) then
          !print_vars_w(:,ivar) = reshape(wvort,(/nelem_w/))
          !print_vars_w(:,i_var) = wvort
          !print_vars_vp(:,i_var) = vpvort
          !print_var_names_w(i_var) = 'Vorticity'
          call add_output_var(out_vars_w(i_var), wvort, &
                 'Singularity_Intensity',.false.)
          call add_output_var(out_vars_vp(i_var), vpvort, &
                 'Singularity_Intensity',.false.)
          i_var = i_var +1
        endif
        if(out_cp) then
          call add_output_var(out_vars_w(i_var), cp, &
                 'Cp',.true.)
          call add_output_var(out_vars_vp(i_var), cp, &
                 'Cp',.true.)
          i_var = i_var +1
        endif
        if(out_vel) then
          call add_output_var(out_vars_w(i_var), vel, &
                 'Velocity',.true.)
          call add_output_var(out_vars_vp(i_var), vel, &
                 'Velocity',.true.)
          i_var = i_var +1
        endif
        if(out_press) then
          call add_output_var(out_vars_w(i_var), press, &
                 'Pressure',.true.)
          call add_output_var(out_vars_vp(i_var), press, &
                 'Pressure',.true.)
          i_var = i_var +1
        endif
  
        !Output the results (with wake)
        select case (trim(out_frmt))
         case ('tecplot')
          filename = trim(filename)//'.plt'
          call  tec_out_viz(filename, t, &
                       points_exp, elems, print_vars, print_var_names, &
                       w_rr=wpoints, w_ee=welems, w_vars=print_vars_w, &
                       w_var_names = print_var_names_w, &
                       vp_rr=vppoints, vp_vars=print_vars_vp, &
                       vp_var_names = print_var_names_w)
         case ('vtk')
          filename = trim(filename)//'.vtu'
          !call  vtk_out_viz(filename, &
          !             points_exp, elems, print_vars, print_var_names, &
          !             w_rr=wpoints, w_ee=welems, w_vars=print_vars_w, &
          !             w_var_names = print_var_names_w, &
          !             vp_rr=vppoints, vp_vars=print_vars_vp, &
          !             vp_var_names = print_var_names_w, &
          !             separate_wake = separate_wake)
          call  vtk_out_viz(filename, &
                       points_exp, elems, out_vars, &
                       w_rr=wpoints, w_ee=welems, w_vars=out_vars_w, &
                       vp_rr=vppoints, vp_vars=out_vars_vp, &
                       separate_wake = separate_wake)
         case default
           call error('dust_post','','Unknown format '//trim(out_frmt)//&
                      ' for visualization output')
         end select
      
        deallocate (wpoints, welems,  wvort)
        call clear_output_vars(out_vars_w)
        call clear_output_vars(out_vars_vp)
        !deallocate(print_var_names_w, print_vars_w)
        !deallocate(print_vars_vp)
  
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
                       points_exp, elems, out_vars)
         case default
           call error('dust_post','','Unknown format '//trim(out_frmt)//&
                      ' for visualization output')
         end select
  
      endif !output wake

    endif !not average
  
    call close_hdf5_file(floc)
  
    deallocate(refs_R, refs_off)
    !deallocate(print_var_names, print_vars)
    call clear_output_vars(out_vars)
  
    if (allocated(vort ) ) deallocate(vort )
    if (allocated(press) ) deallocate(press)
    if (allocated(vel  ) ) deallocate(vel  )
    if (allocated(cp   ) ) deallocate(cp   )
  
  
  end do ! Time loop

  !Print the average
  if(average) then
    write(filename,'(A)') trim(basename)//'_'//trim(an_name)//'_ave'
    !Output the results (without wake)
    select case (trim(out_frmt))
     case ('tecplot')
      filename = trim(filename)//'.plt'
      call  tec_out_viz(filename, t, &
                   points_ave, elems, ave_vars, ave_var_names)
     case ('vtk')
      filename = trim(filename)//'.vtu'
      call  vtk_out_viz(filename, &
                   points_ave, elems, ave_out_vars)
     case default
       call error('dust_post','','Unknown format '//trim(out_frmt)//&
                  ' for visualization output')
    end select
    call clear_output_vars(ave_out_vars)
    deallocate(ave_out_vars)
  endif
  
  deallocate(points)
  call destroy_elements(comps)
  deallocate(comps)
  !deallocate(var_names)
  deallocate(out_vars)
  if(out_wake) deallocate(out_vars_w, out_vars_vp)
  
  
  write(msg,'(A,I0,A)') nl//'++++++++++ Visualization done'//nl
  call printout(trim(msg))

end subroutine post_viz

! ---------------------------------------------------------------------- 

end module mod_post_viz
