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


!> Module to treat the vortex ring wake of actuator disks
module mod_wake_ring

use mod_param, only: &
  wp, nl, pi

use mod_handling, only: &
  error, warning, info, printout, dust_time, t_realtime

use mod_geometry, only: &
  t_geo, t_tedge, calc_geo_data_ad

use mod_aero_elements, only: &
  c_elem, t_elem_p

use mod_vortring, only: &
  t_vortring

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

public :: t_wake_rings, initialize_wake_rings, update_wake_rings, &
          load_wake_rings, destroy_wake_rings

private


type t_wake_rings

 !> Maximum number of rings in the ring wake
 integer :: nrings

 !> Actual number of rings in the ring wake
 integer :: wake_len

 !> Number of generating disks
 integer :: ndisks

 !> Number of points for each "row"
 integer :: np_row

 !> Generating actuator disk elements
 type(t_elem_p), allocatable :: gen_elems(:)


 !> Ring elements
 type(t_actdisk), allocatable :: wake_rings(:,:)

 !> Wake points (unstructured, 3xnpoints)
 !real(wp), allocatable :: w_points(:,:)

 !> vortex intensities
 real(wp), allocatable :: ivort(:,:)

 !> pointer to the wake elements to be passed to the linsys
 !! solver
 type(t_elem_p), allocatable :: pan_p(:)

end type t_wake_rings

character(len=*), parameter :: this_mod_name = 'mod_wake_ring'
!----------------------------------------------------------------------
contains
!----------------------------------------------------------------------

!> Initialize the panel wake
subroutine initialize_wake_rings(wake, geo, nrings)
 type(t_wake_rings), intent(out),target :: wake
 type(t_geo), intent(in), target :: geo
 integer, intent(in) :: nrings

 integer :: ic, nad, ie, npt, id, ir
 integer :: nsides

  !count the number of actuator disks and points
  nad = 0; npt = 0;
  do ic = 1, size(geo%components)
    if(geo%components(ic)%comp_el_type(1:1) .eq. 'a') then
      nad = nad + geo%components(ic)%nelems
      do ie = 1, nad
        npt = npt + geo%components(ic)%el(ie)%n_ver
      enddo
    endif
  enddo


  !set and allocate all the relevant variables
  wake%nrings = nrings
  wake%ndisks = nad
  wake%np_row = npt
  allocate(wake%gen_elems(wake%ndisks))
  !allocate(wake%w_points(3,wake%n_wake_points,npan+1))
  allocate(wake%wake_rings(wake%ndisks,wake%nrings))
  allocate(wake%ivort(wake%ndisks,wake%nrings))
  !allocate(wake%pan_p(0))
  
  !Associate
  id = 1
  do ic = 1, size(geo%components)
    if(geo%components(ic)%comp_el_type(1:1) .eq. 'a') then
      do ie = 1, nad
        wake%gen_elems(id)%p => geo%components(ic)%el(ie)
        id = id+1
      enddo
    endif
  enddo

  do id = 1,wake%ndisks
    do ir = 1,wake%nrings
      wake%wake_rings(id,ir)%idou => wake%ivort(id,ir)
      nsides = wake%gen_elems(id)%p%n_ver
      wake%wake_rings(id,ir)%n_ver = nsides
      allocate(wake%wake_rings(id,ir)%ver(3,nsides))
      allocate(wake%wake_rings(id,ir)%cen(3))
      allocate(wake%wake_rings(id,ir)%nor(3))
      allocate(wake%wake_rings(id,ir)%tang(3,2))
      !TODO: check if the following are really used
      allocate(wake%wake_rings(id,ir)%verp(3,nsides))
      allocate(wake%wake_rings(id,ir)%edge_vec(3,nsides))
      allocate(wake%wake_rings(id,ir)%edge_len(nsides))
      allocate(wake%wake_rings(id,ir)%edge_uni(3,nsides))
      allocate(wake%wake_rings(id,ir)%cosTi(nsides))
      allocate(wake%wake_rings(id,ir)%sinTi(nsides))
    enddo
  enddo


  wake%ivort = 0.0_wp

  !the first line of points is calculated from the mesh points
  !do id = 1,wake%ndisks
  !  wake%wake_rings(id,1)%ver = wake%gen_elems(id)%p%ver
  !enddo
  
  !Starting length of the wake is 
  wake%wake_len = 0


  !TODO : initialize first row of wake here
  allocate(wake%pan_p(0))
  do id = 1,wake%ndisks
    wake%wake_rings(id,:)%moving = wake%gen_elems(id)%p%moving
  !  wake%pan_p(id)%p => wake%wake_panels(id,1)
  enddo

end subroutine initialize_wake_rings

!----------------------------------------------------------------------

!> Destroy a wake panels type by simply passing it as intent(out)
subroutine destroy_wake_rings(wake)
 type(t_wake_rings), intent(out) :: wake

 !dummy to avoid compiler warnings
 wake%nrings = -1

end subroutine

!----------------------------------------------------------------------

!> Prepare the first row of panels to be inserted inside the linear system
!!
!subroutine prepare_wake_panels(wake, geo, dt, uinf)
! type(t_wake_panels), intent(inout) :: wake
! type(t_geo), intent(in) :: geo
! real(wp), intent(in) :: dt
! real(wp), intent(in) :: uinf(3)
! 
! integer :: p1, p2
! integer :: ip, iw, ipan
! real(wp) :: dist(3)
!
!  !Update the first row of panels: set points positions
!
!  !first row  of new points comes from geometry
!  wake%w_start_points = 0.5_wp * (geo%points(:,wake%gen_points(1,:)) + &
!                                  geo%points(:,wake%gen_points(2,:)))
!  wake%w_points(:,:,1) = wake%w_start_points
!
!  !Second row of points: first row + 0.3*|uinf|*t with t = R*t0
!  do ip=1,wake%n_wake_points
!    dist = matmul(geo%iefs(wake%gen_ref(ip))%R_g,wake%gen_dir(:,ip))
!    wake%w_points(:,ip,2) = wake%w_points(:,ip,1) + dist*0.3_wp*norm2(uinf)*dt
!  enddo
!
!  ! Update the panels geometrical quantities of all the panels, the 
!  ! first two row of points have just been updated, the other rows of points
!  ! were updated at the end of the last iteration
!  do ipan = 1,wake%wake_len
!    do iw = 1,wake%n_wake_stripes
!      p1 = wake%i_start_points(1,iw)
!      p2 = wake%i_start_points(2,iw)
!      call calc_geo_data_pan(wake%wake_panels(iw,ipan), &
!           reshape((/wake%w_points(:,p1,ipan),   wake%w_points(:,p2,ipan), &
!                     wake%w_points(:,p2,ipan+1), wake%w_points(:,p1,ipan+1)/),&
!                                                                     (/3,4/)))
!    enddo
!  enddo
!
!end subroutine prepare_wake_panels

!----------------------------------------------------------------------

subroutine load_wake_rings(filename, wake)
 character(len=*), intent(in) :: filename
 type(t_wake_rings), intent(inout), target :: wake

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
  wake%wake_len = min(wake%nrings, size(wvort,2))

  !Load the old results
  ip = 1
  do id = 1,wake%ndisks
    np = wake%gen_elems(id)%p%n_ver
    do ir = 1,wake%wake_len
      !DEBUG:
      write(*,*) 'size ver', shape(wake%wake_rings(id,ir)%ver(:,:))
      write(*,*) 'size results', shape(wpoints(:,ip:ip+np-1,ir))
      wake%wake_rings(id,ir)%ver(:,:) = wpoints(:,ip:ip+np-1,ir)
    enddo
    ip = ip+np
  enddo
  wake%ivort(:,1:wake%wake_len) = wvort(:,1:wake%wake_len)

  deallocate(wake%pan_p); allocate(wake%pan_p(wake%ndisks*wake%wake_len))

  !Update the geometrical quantities
  ipt = 1
  do ir = 1,wake%wake_len
    do id = 1,wake%ndisks
      call calc_geo_data_ad(wake%wake_rings(id,ir), &
                    wake%wake_rings(id,ir)%ver)
      wake%pan_p(ipt)%p => wake%wake_rings(id,ir)
      ipt = ipt + 1
    enddo
  enddo

end subroutine load_wake_rings

!----------------------------------------------------------------------

!> Update the position and the intensities of the wake panels
!!
!! Note: at this subroutine is passed the whole array of elements,
!! comprising both the implicit panels and the explicit (ll) 
!! elements
subroutine update_wake_rings(wake, elems, wake_pan_p, dt, uinf)
 type(t_wake_rings), intent(inout), target :: wake
 type(t_elem_p), intent(in) :: elems(:)
 type(t_elem_p), intent(in) :: wake_pan_p(:)
 real(wp), intent(in) :: dt
 real(wp), intent(in) :: uinf(3)

 integer :: ip,id,ir, ie, np
 real(wp) :: pos_p(3), vel_p(3), v(3)
 type(t_elem_p), allocatable :: pan_p_temp(:)
 real(wp), allocatable :: points(:,:,:)
 logical :: increase_wake
 integer :: size_old

  !==> 1) Update wake points position ==
  
  !Save the old positions for the integration
  increase_wake = .false.
  if(wake%wake_len .lt. wake%nrings) then
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
      np = wake%wake_rings(id,ir)%n_ver
      points(:,ip:ip+np-1,ir+1) = wake%wake_rings(id,ir)%ver(:,:)
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
        call elems(ie)%p%compute_vel(pos_p, uinf, v)
        vel_p = vel_p + v/(4*pi)
      enddo

      ! calculate the influence of the panel wake
      do ie=1,size(wake_pan_p)
        v = 0.0_wp
        call wake_pan_p(ie)%p%compute_vel(pos_p, uinf, v)
        vel_p = vel_p + v/(4*pi)
      enddo

      ! calculate the influence of the ring wake
      do ie=1,size(wake%pan_p)
        v = 0.0_wp
        call wake%pan_p(ie)%p%compute_vel(pos_p, uinf, v)
        vel_p = vel_p + v/(4*pi)
      enddo

      !calculate the influence of particles

      vel_p    = vel_p   +uinf   

      !update the position in time
      points(:,ip,ir) = points(:,ip,ir) + vel_p*dt
    enddo !ir
  enddo !ip
!$omp end parallel do
 
  !redistribute the points
  do ir = 1,wake%wake_len
    ip = 1
    do id = 1,wake%ndisks
      np = wake%wake_rings(id,ir)%n_ver
      wake%wake_rings(id,ir)%ver = points(:,ip:ip+np-1,ir)
      ip = ip + np
    enddo
  enddo

  deallocate(points)

 
  !==> 3) Increase the length of the wake, if it is necessary
  if (increase_wake) then
    allocate(pan_p_temp(wake%ndisks*wake%wake_len))
    size_old = 0; if(allocated(wake%pan_p)) size_old = size(wake%pan_p)

    do ip = 1,size_old
      pan_p_temp(ip) = wake%pan_p(ip)
    enddo
    do id = 1,wake%ndisks
      pan_p_temp(size_old+id)%p => wake%wake_rings(id,wake%wake_len)
    enddo

    !manual move alloc
    if(allocated(wake%pan_p)) deallocate(wake%pan_p)
    allocate(wake%pan_p(size(pan_p_temp)))
    do ip = 1,size(wake%pan_p)
      wake%pan_p(ip) = pan_p_temp(ip)
    enddo
    deallocate(pan_p_temp)
  endif

  !==> 4) Update the intensities of the panels
  !       From the back, all the vortex intensities come from the previous panel
  !==> 1) Update the first row of vortex intensities: 
  !      it was already calculated (implicitly) in the linear system
  do ir = wake%wake_len,2,-1
    do id = 1,wake%ndisks
      wake%wake_rings(id,ir)%idou = wake%wake_rings(id,ir-1)%idou
    enddo
  enddo
  do id = 1,wake%ndisks
    wake%wake_rings(id,1)%idou  = wake%gen_elems(id)%p%idou
  enddo

  !Update the geometrical quantities
  do ir = 1,wake%wake_len
    do id = 1,wake%ndisks
      call calc_geo_data_ad(wake%wake_rings(id,ir), &
                    wake%wake_rings(id,ir)%ver)
    enddo
  enddo



end subroutine update_wake_rings

!----------------------------------------------------------------------

end module mod_wake_ring
