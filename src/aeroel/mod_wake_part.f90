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


!> Module to treat the vortex particles wake
module mod_wake_part

use mod_param, only: &
  wp, nl, pi

use mod_sim_param, only: &
  t_sim_param

use mod_handling, only: &
  error, warning, info, printout, dust_time, t_realtime

use mod_geometry, only: &
  t_geo, t_tedge, calc_geo_data_ad

!use mod_aero_elements, only: &
!  c_elem, t_elem_p

use mod_aeroel, only: &
  c_elem, c_pot_elem, c_vort_elem, c_impl_elem, c_expl_elem, &
  t_elem_p, t_pot_elem_p, t_vort_elem_p, t_impl_elem_p, t_expl_elem_p

use mod_vortlatt, only: &
  t_vortlatt

use mod_actuatordisk, only: &
  t_actdisk

use mod_doublet, only: &
  velocity_calc_doublet

use mod_hdf5_io, only: &
   h5loc, &
   new_hdf5_file, &
   open_hdf5_file, &
   close_hdf5_file, &
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

public :: t_wake_parts, initialize_wake_parts, update_wake_parts, &
          load_wake_parts, destroy_wake_parts

private


type t_wake_parts

 !> Maximum number of parts in the part wake
 integer :: nparts

 !> Actual number of parts in the part wake
 integer :: wake_len

 !> Number of generating disks
 integer :: ndisks

 !> Number of points for each "row"
 integer :: np_row

 !> Generating actuator disk elements
 type(t_pot_elem_p), allocatable :: gen_elems(:)


 !> Ring elements
 type(t_actdisk), allocatable :: wake_parts(:,:)

 !> Wake points (unstructured, 3xnpoints)
 !real(wp), allocatable :: w_points(:,:)

 !> vortex intensities
 real(wp), allocatable :: idou(:,:)

 !> pointer to the wake elements to be passed to the linsys
 !! solver
 type(t_pot_elem_p), allocatable :: pan_p(:)

end type t_wake_parts

character(len=*), parameter :: this_mod_name = 'mod_wake_part'
!----------------------------------------------------------------------
contains
!----------------------------------------------------------------------

!> Initialize the panel wake
subroutine initialize_wake_parts(wake, geo, nparts)
 type(t_wake_parts), intent(out),target :: wake
 type(t_geo), intent(in), target :: geo
 integer, intent(in) :: nparts

 integer :: ic, nad, ie, npt, id, ir
 integer :: nsides

  !count the number of actuator disks and points
  nad = 0; npt = 0;
  do ic = 1, size(geo%components)
    if(geo%components(ic)%comp_el_type(1:1) .eq. 'a') then
      nad = nad + geo%components(ic)%nelems
      do ie = 1, size(geo%components(ic)%el)
        npt = npt + geo%components(ic)%el(ie)%n_ver
      enddo
    endif
  enddo


  !set and allocate all the relevant variables
  wake%nparts = nparts
  wake%ndisks = nad
  wake%np_row = npt
  allocate(wake%gen_elems(wake%ndisks))
  !allocate(wake%w_points(3,wake%n_wake_points,npan+1))
  allocate(wake%wake_parts(wake%ndisks,wake%nparts))
  allocate(wake%idou(wake%ndisks,wake%nparts))
  !allocate(wake%pan_p(0))

  !Associate
  id = 1
  do ic = 1, size(geo%components)
    if(geo%components(ic)%comp_el_type(1:1) .eq. 'a') then
      do ie = 1,size(geo%components(ic)%el)
        wake%gen_elems(id)%p => geo%components(ic)%el(ie)
        id = id+1
      enddo
    endif
  enddo

  do id = 1,wake%ndisks
    do ir = 1,wake%nparts
      wake%wake_parts(id,ir)%mag => wake%idou(id,ir)
      nsides = wake%gen_elems(id)%p%n_ver
      wake%wake_parts(id,ir)%n_ver = nsides
      allocate(wake%wake_parts(id,ir)%ver(3,nsides))
      !allocate(wake%wake_parts(id,ir)%cen(3))
      !allocate(wake%wake_parts(id,ir)%nor(3))
      !allocate(wake%wake_parts(id,ir)%tang(3,2))
      !TODO: check if the following are really used
      allocate(wake%wake_parts(id,ir)%edge_vec(3,nsides))
      allocate(wake%wake_parts(id,ir)%edge_len(nsides))
      allocate(wake%wake_parts(id,ir)%edge_uni(3,nsides))
    enddo
  enddo


  wake%idou = 0.0_wp

  !Starting length of the wake is
  wake%wake_len = 0


  allocate(wake%pan_p(0))
  do id = 1,wake%ndisks
    wake%wake_parts(id,:)%moving = wake%gen_elems(id)%p%moving
  enddo

end subroutine initialize_wake_parts

!----------------------------------------------------------------------

!> Destroy a wake panels type by simply passing it as intent(out)
subroutine destroy_wake_parts(wake)
 type(t_wake_parts), intent(out) :: wake

 !dummy to avoid compiler warnings
 wake%nparts = -1

end subroutine

!----------------------------------------------------------------------

!> Load the wake parts solution from a previous result
subroutine load_wake_parts(filename, wake)
 character(len=*), intent(in) :: filename
 type(t_wake_parts), intent(inout), target :: wake

 integer(h5loc) :: floc, gloc
 real(wp), allocatable :: wpoints(:,:,:), wvort(:,:)
 integer, allocatable :: conn_pe(:)
 integer :: id, ir, ip, np, ipt
 character(len=*), parameter :: this_sub_name = 'load_rin_wake'

  !Read the past results
  call open_hdf5_file(filename, floc)
  call open_hdf5_group(floc, 'RingWake', gloc)

  call read_hdf5_al(wpoints,'WakePoints',gloc)
  call read_hdf5_al(conn_pe,'Conn_pe',gloc)
  call read_hdf5_al(wvort,'WakeVort',gloc)

  call close_hdf5_group(gloc)
  call close_hdf5_file(floc)

  !check the consistency
  ip = 1
  do id = 1,wake%ndisks
    np = wake%gen_elems(id)%p%n_ver
    if (.not. all(conn_pe(ip:ip+np-1) .eq. id)) call error( &
    this_sub_name, this_mod_name, 'Different wake diske connectivity&
    & between the loded and built geometry')
    ip = ip+np
  enddo

  !the wake length is the loaded one, or less if imposed so
  wake%wake_len = min(wake%nparts, size(wvort,2))

  !Load the old results
  ip = 1
  do id = 1,wake%ndisks
    np = wake%gen_elems(id)%p%n_ver
    do ir = 1,wake%wake_len
      wake%wake_parts(id,ir)%ver(:,:) = wpoints(:,ip:ip+np-1,ir)
    enddo
    ip = ip+np
  enddo
  wake%idou(:,1:wake%wake_len) = wvort(:,1:wake%wake_len)

  deallocate(wake%pan_p); allocate(wake%pan_p(wake%ndisks*wake%wake_len))

  !Update the geometrical quantities
  ipt = 1
  do ir = 1,wake%wake_len
    do id = 1,wake%ndisks
      !call calc_geo_data_ad(wake%wake_parts(id,ir), &
      !              wake%wake_parts(id,ir)%ver)
      call wake%wake_parts(id,ir)%calc_geo_data( &
                    wake%wake_parts(id,ir)%ver)
      wake%pan_p(ipt)%p => wake%wake_parts(id,ir)
      ipt = ipt + 1
    enddo
  enddo

end subroutine load_wake_parts

!----------------------------------------------------------------------

!> Update the position and the intensities of the wake panels
!!
!! Note: at this subroutine is passed the whole array of elements,
!! comprising both the implicit panels and the explicit (ll)
!! elements
subroutine update_wake_parts(wake, elems, wake_pan_p, sim_param)
 type(t_wake_parts), intent(inout), target :: wake
 type(t_pot_elem_p), intent(in) :: elems(:)
 type(t_pot_elem_p), intent(in) :: wake_pan_p(:)
 type(t_sim_param), intent(in) :: sim_param

 integer :: ip,id,ir, ie, np
 real(wp) :: pos_p(3), vel_p(3), v(3)
 type(t_pot_elem_p), allocatable :: pan_p_temp(:)
 real(wp), allocatable :: points(:,:,:)
 logical :: increase_wake
 integer :: size_old

  !==> 1) Update wake points position ==

  !Save the old positions for the integration
  increase_wake = .false.
  if(wake%wake_len .lt. wake%nparts) then
    wake%wake_len = wake%wake_len+1
    increase_wake = .true.
  endif
  allocate(points(3,wake%np_row,wake%wake_len))

  !Store at the beginning the disk points
  ip=1
  do id = 1,wake%ndisks
    np = wake%gen_elems(id)%p%n_ver
    points(:,ip:ip+np-1,1) = wake%gen_elems(id)%p%ver
    ip = ip+np
  enddo

  !Then store the old points of the rest of the wake (the shift forward
  ! of the points in the array is happening here)
  do ir = 1,wake%wake_len-1
    ip=1
    do id = 1,wake%ndisks
      np = wake%wake_parts(id,ir)%n_ver
      points(:,ip:ip+np-1,ir+1) = wake%wake_parts(id,ir)%ver(:,:)
      ip = ip+np
    enddo
  enddo


  !calculate the velocities at the old positions of the points and then
  !update the positions

!$omp parallel do private(pos_p, vel_p, ip, v)
  do ip = 1,size(points,2)
    do ir = 1,size(points,3)
      pos_p = points(:,ip,ir)
      vel_p = 0.0_wp

      !calculate the influence of the solid bodies
      do ie=1,size(elems)
        v = 0.0_wp
        call elems(ie)%p%compute_vel(pos_p, sim_param%u_inf, v)
        vel_p = vel_p + v/(4*pi)
      enddo

      ! calculate the influence of the panel wake
      do ie=1,size(wake_pan_p)
        v = 0.0_wp
        call wake_pan_p(ie)%p%compute_vel(pos_p, sim_param%u_inf, v)
        vel_p = vel_p + v/(4*pi)
      enddo

      ! calculate the influence of the part wake
      do ie=1,size(wake%pan_p)
        v = 0.0_wp
        call wake%pan_p(ie)%p%compute_vel(pos_p, sim_param%u_inf, v)
        vel_p = vel_p + v/(4*pi)
      enddo

      !calculate the influence of particles

      vel_p    = vel_p   +sim_param%u_inf

      !update the position in time
      points(:,ip,ir) = points(:,ip,ir) + vel_p*sim_param%dt
    enddo !ir
  enddo !ip
!$omp end parallel do

  !redistribute the points
  do ir = 1,wake%wake_len
    ip = 1
    do id = 1,wake%ndisks
      np = wake%wake_parts(id,ir)%n_ver
      wake%wake_parts(id,ir)%ver = points(:,ip:ip+np-1,ir)
      ip = ip + np
    enddo
  enddo

  deallocate(points)


  !==> 2) Increase the length of the wake, if it is necessary
  if (increase_wake) then
    allocate(pan_p_temp(wake%ndisks*wake%wake_len))
    size_old = 0; if(allocated(wake%pan_p)) size_old = size(wake%pan_p)

    do ip = 1,size_old
      pan_p_temp(ip) = wake%pan_p(ip)
    enddo
    do id = 1,wake%ndisks
      pan_p_temp(size_old+id)%p => wake%wake_parts(id,wake%wake_len)
    enddo

    !manual move alloc
    if(allocated(wake%pan_p)) deallocate(wake%pan_p)
    allocate(wake%pan_p(size(pan_p_temp)))
    do ip = 1,size(wake%pan_p)
      wake%pan_p(ip) = pan_p_temp(ip)
    enddo
    deallocate(pan_p_temp)
  endif

  !==> 3) Update the intensities of the panels
  !       From the back, all the vortex intensities come from the 
  !       previous panel
  do ir = wake%wake_len,2,-1
    do id = 1,wake%ndisks
      wake%wake_parts(id,ir)%mag = wake%wake_parts(id,ir-1)%mag
    enddo
  enddo
  do id = 1,wake%ndisks
    wake%wake_parts(id,1)%mag  = wake%gen_elems(id)%p%mag
  enddo

  !Update the geometrical quantities
  do ir = 1,wake%wake_len
    do id = 1,wake%ndisks
      !call calc_geo_data_ad(wake%wake_parts(id,ir), &
      !              wake%wake_parts(id,ir)%ver)
      call wake%wake_parts(id,ir)%calc_geo_data( &
                    wake%wake_parts(id,ir)%ver)
    enddo
  enddo



end subroutine update_wake_parts

!----------------------------------------------------------------------

end module mod_wake_part
