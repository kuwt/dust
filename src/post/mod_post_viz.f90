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
  t_geo, t_geo_component

use mod_geo_postpro, only: &
  load_components_postpro, update_points_postpro, prepare_geometry_postpro, &
  expand_actdisk_postpro !, prepare_wake_postpro

use mod_tecplot_out, only: &
  tec_out_viz ! , tec_out_probes, tec_out_box, tec_out_loads

use mod_vtk_out, only: &
  vtk_out_viz ! , vtr_write

use mod_post_load, only: &
  load_refs, load_res, load_wake_viz

implicit none

public :: post_viz

private

contains

! ---------------------------------------------------------------------- 

subroutine post_viz( sbprms , basename , data_basename , an_name , ia , &
                     out_frmt , comps , components_names , all_comp , &
                     an_start , an_end , an_step )
type(t_parse), pointer :: sbprms
character(len=*) , intent(in) :: basename
character(len=*) , intent(in) :: data_basename
character(len=*) , intent(in) :: an_name
integer          , intent(in) :: ia
character(len=*) , intent(in) :: out_frmt
type(t_geo_component), allocatable , intent(inout) :: comps(:)
character(len=max_char_len), allocatable , intent(inout) :: components_names(:)
logical , intent(in) :: all_comp
integer , intent(in) :: an_start , an_end , an_step

character(len=max_char_len) :: filename
integer(h5loc) :: floc , ploc
logical :: out_vort, out_vel, out_cp, out_press , out_wake
integer :: n_var , i_var
character(len=max_char_len), allocatable :: var_names(:)
real(wp), allocatable :: points(:,:), points_exp(:,:) , wpoints(:,:)
integer , allocatable :: elems(:,:) , welems(:,:)
integer :: nelem , nelem_w

real(wp) :: u_inf(3)
real(wp) :: P_inf , rho

real(wp), allocatable :: refs_R(:,:,:), refs_off(:,:)
real(wp), allocatable :: vort(:), cp(:), vel(:), press(:), wvort(:)

real(wp), allocatable :: print_vars(:,:)
character(len=max_char_len), allocatable :: print_var_names(:)
real(wp), allocatable :: print_vars_w(:,:)
character(len=max_char_len), allocatable :: print_var_names_w(:)
integer :: nprint , nelem_out

integer :: it
real(wp) :: t

    write(*,*) nl//' Analysis:',ia,' post_viz() +++++++++++++++ '//nl

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

! Load the components (just once)
call open_hdf5_file(trim(data_basename)//'_geo.h5', floc)

!DEBUG
write(*,*) nl//' Before load_components ... '
if ( allocated(components_names) ) then
  write(*,*) ' size(components_names) ' , size(components_names)
  do i_var = 1 , size(components_names)
    write(*,*) ' components_names(',i_var,'): ' , trim(components_names(i_var))
  end do
else 
  write(*,*) ' components_names .not. allocated '  
end if

call load_components_postpro(comps, points, nelem, floc, &
                             components_names, all_comp)
!DEBUG
write(*,*) nl//' After load_components_postpro. size(comps) : ' , size(comps)
do i_var = 1 , size(comps)
  write(*,*) ' comps(',i_var,')%comp_name : ' , trim(comps(i_var)%comp_name)
! write(*,*) ' size(comps(',i_var,')%el)  : ' , size(comps(i_var)%el)
end do

call close_hdf5_file(floc)

! Print the wake or not 
out_wake = getlogical(sbprms,'Wake')

! Prepare_geometry_postpro
call prepare_geometry_postpro(comps)

! Time loop
do it = an_start, an_end, an_step

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

      !Load the results ! TODO: check this routine and the content of the files to be read
      ! TODO : compute the missing quantities
      call load_res(floc, comps, vort, cp, t)

      !Prepare the variable for output
      nelem_out = size(vort)
      nprint = 0
      if(out_vort)  nprint = nprint+1
      if(out_cp)    nprint = nprint+1
      if(out_vel)   nprint = nprint+1  !<--- *** TODO ***
      if(out_press) nprint = nprint+1  !<--- *** TODO ***
      ! TODO: compute/or read pressure and velocity field. Now set equal to zero
      allocate(  vel(size(vort,1)) ) ; vel   = 0.0_wp
      allocate(press(size(  cp,1)) ) ; press = 0.0_wp

      allocate(print_var_names(nprint), print_vars(nelem_out, nprint))
      
      i_var = 1
      if(out_vort) then
        print_vars(:,i_var) = vort
        print_var_names(i_var) = 'Vorticity'
        i_var = i_var +1
      endif
      if(out_cp) then
        print_vars(:,i_var) = cp
        print_var_names(i_var) = 'Cp'
        i_var = i_var +1
      endif
      if(out_vel) then
        print_vars(:,i_var) = vel
        print_var_names(i_var) = 'Velocity'
        i_var = i_var +1
      endif
      if(out_press) then
        print_vars(:,i_var) = press
        print_var_names(i_var) = 'Pressure'
        i_var = i_var +1
      endif

      ! Output filename
      write(filename,'(A,I4.4)') trim(basename)//'_'//trim(an_name)//'_',it
      
      if (out_wake) then
        
        call load_wake_viz(floc, wpoints, welems, wvort)
        nelem_w = size(welems,2)

        nprint = 0
        if(out_vort) nprint = nprint+1
        allocate(print_var_names_w(nprint), print_vars_w(nelem_w, nprint))
        
        i_var = 1
        if(out_vort) then
          !print_vars_w(:,ivar) = reshape(wvort,(/nelem_w/))
          print_vars_w(:,i_var) = wvort
          print_var_names_w(i_var) = 'Vorticity'
          i_var = i_var +1
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

      deallocate(refs_R, refs_off, vort, cp, vel, press)
      deallocate(print_var_names, print_vars)



end do ! Time loop

deallocate(comps, points,components_names)
deallocate(var_names)


    write(*,*) nl//' post_viz done.'//nl
! subroutine post_viz( )

end subroutine post_viz

! ---------------------------------------------------------------------- 

end module mod_post_viz
