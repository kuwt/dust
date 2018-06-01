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


!> Module to treat the vortex panel wake
module mod_wake_pan

use mod_param, only: &
  wp, nl, pi

use mod_handling, only: &
  error, warning, info, printout, dust_time, t_realtime

use mod_geometry, only: &
  t_geo, t_tedge, calc_geo_data_pan, calc_node_vel

use mod_aero_elements, only: &
  c_elem, t_elem_p

use mod_vortring, only: &
  t_vortring

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

public :: t_wake_panels, initialize_wake_panels, update_wake_panels, &
          prepare_wake_panels, load_wake_panels, destroy_wake_panels

private

!> Type containing wake panels information
type :: t_wake_panels
 
 !> Number of maximum streamwise panels
 integer :: npan

 !> Number of actual streamwise panels
 integer :: wake_len

 !> Number of wake stripes ("spanwise" panels)
 integer :: n_wake_stripes

 !> Number of wake points in the "spanwise" direction
 integer :: n_wake_points

 !> Index of the 2 generating elements of the wake
 !! (2 x n_wake_stripes)
 !integer, allocatable :: gen_elems(:,:)
 type(t_elem_p), allocatable :: gen_elems(:,:)
 integer, allocatable :: gen_elems_id(:,:)

 !> Index of the 2 generating points of each wake point
 !! (2 x n_wake_points)
 integer, allocatable :: gen_points(:,:)

 !> Direction of the wake at the trailing edge
 !! (3xn_wake_points)
 real(wp), allocatable :: gen_dir(:,:)

 !> Reference frame of the generating points (id)
 !! (n_wake_points)
 integer, allocatable :: gen_ref(:)

 !> Wake starting points: calculated from the 2 (or 1) starting points of the 
 !! geometry possibly moved 
 real(wp), allocatable :: w_start_points(:,:)

 !> Index of the 2 wake starting points for each wake stripe
 !! (2 x n_wake_stripes)
 integer, allocatable :: i_start_points(:,:)

 !> Points of the wake, in a structured way
 !! (3 x n_wake_points x npan+1)
 real(wp), allocatable :: w_points(:,:,:)
 
 !> Relative velocity ( u_inf - ub ) of the nodes at the TE
 !! (3 x n_wake_points) 
 !! used to determine the first prescribed panel of the wake
 real(wp), allocatable :: w_vel_te(:,:)

 !> Velocity at the nodes of the wake. For output only
 !! (3 x n_wake_points x npan+1)
 real(wp), allocatable :: w_vel(:,:,:)

 !> elements of the panels
 type(t_vortring), allocatable :: wake_panels(:,:)

 !> vortex intensities
 real(wp), allocatable :: ivort(:,:)

 !> pointer to the wake elements to be passed to the linsys
 !! solver
 type(t_elem_p), allocatable :: pan_p(:)
 
end type 

character(len=*), parameter :: this_mod_name='mod_wake_pan'

!----------------------------------------------------------------------
contains
!----------------------------------------------------------------------

!> Initialize the panel wake
subroutine initialize_wake_panels(wake, geo, te,  npan)
 type(t_wake_panels), intent(out),target :: wake
 type(t_geo), intent(in) :: geo
 type(t_tedge), intent(in) :: te
 integer, intent(in) :: npan

 integer :: iw, ip, nsides


  !set and allocate all the relevant variables
  wake%npan = npan
  wake%n_wake_stripes = size(te%e,2)
  wake%n_wake_points  = size(te%i,2)
  allocate(wake%gen_elems(2,wake%n_wake_stripes))
  allocate(wake%gen_elems_id(2,wake%n_wake_stripes))
  allocate(wake%gen_points(2,wake%n_wake_points))
  allocate(wake%gen_dir(3,wake%n_wake_points))
  allocate(wake%gen_ref(wake%n_wake_points))
  allocate(wake%w_start_points(3,wake%n_wake_points))
  allocate(wake%i_start_points(2,wake%n_wake_stripes))
  allocate(wake%w_points(3,wake%n_wake_points,npan+1))
  allocate(wake%w_vel(3,wake%n_wake_points,npan+1))
  allocate(wake%wake_panels(wake%n_wake_stripes,npan))
  allocate(wake%ivort(wake%n_wake_stripes,npan))
  !allocate(wake%pan_p(0))

  !Associate for all the panels the relevant intensity and allocate all the
  !relevant fields
  nsides = 4
  do ip = 1,npan
    do iw=1,wake%n_wake_stripes
     wake%wake_panels(iw,ip)%idou => wake%ivort(iw,ip)
     allocate(wake%wake_panels(iw,ip)%ver(3,nsides))
     allocate(wake%wake_panels(iw,ip)%cen(3))
     allocate(wake%wake_panels(iw,ip)%nor(3))
     allocate(wake%wake_panels(iw,ip)%tang(3,2))
     allocate(wake%wake_panels(iw,ip)%verp(3,nsides))
     allocate(wake%wake_panels(iw,ip)%edge_vec(3,nsides))
     allocate(wake%wake_panels(iw,ip)%edge_len(nsides))
     allocate(wake%wake_panels(iw,ip)%edge_uni(3,nsides))
     allocate(wake%wake_panels(iw,ip)%cosTi(nsides))
     allocate(wake%wake_panels(iw,ip)%sinTi(nsides))
    enddo
  enddo

  !Set the generating elements and points from the trailing edge
  wake%gen_elems  = te%e
  wake%gen_points = te%i
  wake%gen_dir = te%t
  wake%gen_ref = te%ref
  wake%i_start_points = te%ii

  do iw=1,wake%n_wake_stripes
    wake%gen_elems_id(1,iw) = wake%gen_elems(1,iw)%p%id
    if(associated(wake%gen_elems(2,iw)%p)) then
      wake%gen_elems_id(2,iw) = wake%gen_elems(2,iw)%p%id
    else
      wake%gen_elems_id(2,iw) = 0
    endif
  enddo

  wake%ivort = 0.0_wp


  !the first line of points is calculated from the mesh points
  wake%w_start_points = 0.5_wp * (geo%points(:,wake%gen_points(1,:)) + &
                                  geo%points(:,wake%gen_points(2,:)))
  wake%w_points(:,:,1) = wake%w_start_points
  
  !Starting length of the wake is 
  wake%wake_len = 1


  !TODO : initialize first row of wake here
  allocate(wake%pan_p(wake%n_wake_stripes))
  do iw = 1,wake%n_wake_stripes
    wake%wake_panels(iw,:)%moving = geo%refs(wake%gen_ref(iw))%moving
    wake%pan_p(iw)%p => wake%wake_panels(iw,1)
  enddo

end subroutine initialize_wake_panels

!----------------------------------------------------------------------

!> Destroy a wake panels type by simply passing it as intent(out)
subroutine destroy_wake_panels(wake)
 type(t_wake_panels), intent(out) :: wake

 !dummy to avoid compiler warnings
 wake%npan = -1

end subroutine

!----------------------------------------------------------------------

!> Prepare the first row of panels to be inserted inside the linear system
!!
subroutine prepare_wake_panels(wake, geo, dt, uinf)
 type(t_wake_panels), intent(inout) :: wake
 type(t_geo), intent(in) :: geo
 real(wp), intent(in) :: dt
 real(wp), intent(in) :: uinf(3)
 
 integer :: p1, p2
 integer :: ip, iw, ipan
 real(wp) :: dist(3) , vel_te(3)

  !Update the first row of panels: set points positions

  !first row  of new points comes from geometry
  wake%w_start_points = 0.5_wp * (geo%points(:,wake%gen_points(1,:)) + &
                                  geo%points(:,wake%gen_points(2,:)))
  wake%w_points(:,:,1) = wake%w_start_points

  !Second row of points: first row + 0.3*|uinf|*t with t = R*t0
  do ip=1,wake%n_wake_points
    dist = matmul(geo%refs(wake%gen_ref(ip))%R_g,wake%gen_dir(:,ip))
    call calc_node_vel( wake%w_start_points(:,ip), &
            geo%refs(wake%gen_ref(ip))%G_g, &
            geo%refs(wake%gen_ref(ip))%f_g, &
            vel_te )
!   !DEBUG
!   write(*,*) 'ip , wake%gen_ref(ip) : ' , ip , wake%gen_ref(ip)
!   write(*,*) ' vel_te , dt : ' ,vel_te ,dt
    wake%w_points(:,ip,2) = wake%w_points(:,ip,1) + dist*0.3_wp*norm2(uinf-vel_te)*dt
  enddo

  ! Update the panels geometrical quantities of all the panels, the 
  ! first two row of points have just been updated, the other rows of points
  ! were updated at the end of the last iteration
  do ipan = 1,wake%wake_len
    do iw = 1,wake%n_wake_stripes
      p1 = wake%i_start_points(1,iw)
      p2 = wake%i_start_points(2,iw)
      call calc_geo_data_pan(wake%wake_panels(iw,ipan), &
           reshape((/wake%w_points(:,p1,ipan),   wake%w_points(:,p2,ipan), &
                     wake%w_points(:,p2,ipan+1), wake%w_points(:,p1,ipan+1)/),&
                                                                     (/3,4/)))
    enddo
  enddo

end subroutine prepare_wake_panels

!----------------------------------------------------------------------

subroutine load_wake_panels(filename, wake)
 character(len=*), intent(in) :: filename
 type(t_wake_panels), intent(inout), target :: wake

 integer(h5loc) :: floc, gloc
 real(wp), allocatable :: wpoints(:,:,:), wvort(:,:)
 integer, allocatable :: start_points(:,:)
 integer :: ipan, iw, p1, p2, ipt
 character(len=*), parameter :: this_sub_name = 'load_pan_wake'

  !Read the past results
  call open_hdf5_file(filename, floc)
  call open_hdf5_group(floc, 'PanelWake', gloc)

  call read_hdf5_al(wpoints,'WakePoints',gloc)
  call read_hdf5_al(start_points,'StartPoints',gloc)
  call read_hdf5_al(wvort,'WakeVort',gloc)

  call close_hdf5_group(gloc)
  call close_hdf5_file(floc)

  !Perform a few checks to be sure that the correct size solution was loaded
  if(.not. all(start_points .eq. wake%i_start_points)) call error( &
    this_sub_name, this_mod_name, 'Different wake trailing edge connectivity&
    & between the loded and built geometry')
  
  !the wake length is the loaded one, or less if imposed so
  wake%wake_len = min(wake%npan, size(wvort,2))

  !store points position and doublets intensity
  wake%w_points(:,:,1:wake%wake_len+1) = wpoints(:,:,1:wake%wake_len+1)
  wake%ivort(:,1:wake%wake_len) = wvort(:,1:wake%wake_len)

  deallocate(wake%pan_p)
  allocate(wake%pan_p(wake%n_wake_stripes*wake%wake_len))
  ! Update the panels geometrical quantities of all the panels
  ipt = 1
  do ipan = 1,wake%wake_len
    do iw = 1,wake%n_wake_stripes
      p1 = wake%i_start_points(1,iw)
      p2 = wake%i_start_points(2,iw)
      call calc_geo_data_pan(wake%wake_panels(iw,ipan), &
           reshape((/wake%w_points(:,p1,ipan),   wake%w_points(:,p2,ipan), &
                    wake%w_points(:,p2,ipan+1), wake%w_points(:,p1,ipan+1)/),&
                                                                     (/3,4/)))
      wake%pan_p(ipt)%p => wake%wake_panels(iw,ipan)
      ipt = ipt + 1 
    enddo
  enddo


  deallocate(wpoints, wvort, start_points)
end subroutine load_wake_panels


!----------------------------------------------------------------------

!> Update the position and the intensities of the wake panels
!!
!! Note: at this subroutine is passed the whole array of elements,
!! comprising both the implicit panels and the explicit (ll) 
!! elements
subroutine update_wake_panels(wake, elems, wake_rings_p, dt, uinf)
 type(t_wake_panels), intent(inout), target :: wake
 type(t_elem_p), intent(in) :: elems(:)
 type(t_elem_p), intent(in) :: wake_rings_p(:)
 real(wp), intent(in) :: dt
 real(wp), intent(in) :: uinf(3)

 integer :: iw, ipan, ie, ip, np
 real(wp) :: pos_p(3), vel_p(3), v(3)
 type(t_elem_p), allocatable :: pan_p_temp(:)
 real(wp), allocatable :: point_old(:,:,:)

  wake%w_vel = 0.0_wp
  
  !==> 1) Update the first row of vortex intensities: 
  !      it was already calculated (implicitly) in the linear system
  do iw = 1,wake%n_wake_stripes
    ! 
    if      ( associated(wake%gen_elems(2,iw)%p) ) then
      !wake%wake_panels(iw,1)%idou  = elems(wake%gen_elems(1,iw))%p%idou - &
      !                               elems(wake%gen_elems(2,iw))%p%idou
      wake%wake_panels(iw,1)%idou  = wake%gen_elems(1,iw)%p%idou - &
                                     wake%gen_elems(2,iw)%p%idou
    else if ( .not. associated(wake%gen_elems(2,iw)%p) ) then
      !wake%wake_panels(iw,1)%idou  = elems(wake%gen_elems(1,iw))%p%idou
      wake%wake_panels(iw,1)%idou  = wake%gen_elems(1,iw)%p%idou
    end if
  enddo

  !==> 2) Update wake points position ==
  
  !Save the old positions for the integration
  allocate(point_old(size(wake%w_points,1),size(wake%w_points,2), &
                                                        size(wake%w_points,3)))
  point_old = wake%w_points


  !calculate the velocities at the old positions of the points and then
  !update the positions (from the third row of points: the first is the 
  !trailing edge, the second is extrapolated from the trailing edge)
  np = wake%wake_len+1
  if(wake%wake_len .lt. wake%npan) np = np + 1
  
!$omp parallel do collapse(2) private(pos_p, vel_p, ie, ipan, iw, v)
  do ipan = 3,np
    do iw = 1,wake%n_wake_points
      pos_p = point_old(:,iw,ipan-1) 
      vel_p = 0.0_wp

      !calculate the influence of the solid bodies 
      do ie=1,size(elems)
        v = 0.0_wp
        call elems(ie)%p%compute_vel(pos_p, uinf, v)
        vel_p = vel_p + v/(4*pi)
      enddo

      ! calculate the influence of the wake panels
      do ie=1,size(wake%pan_p)
        v = 0.0_wp
        call wake%pan_p(ie)%p%compute_vel(pos_p, uinf, v)
        vel_p = vel_p + v/(4*pi)
      enddo

      ! calculate the influence of the wake rings
      do ie=1,size(wake_rings_p)
        v = 0.0_wp
        call wake_rings_p(ie)%p%compute_vel(pos_p, uinf, v)
        vel_p = vel_p + v/(4*pi)
      enddo
      !calculate the influence of particles

      ! for OUTPUT only -----
      wake%w_vel(:,iw,ipan-1) = vel_p

!     vel_p(1) = vel_p(1)+uinf(1)
      vel_p    = vel_p   +uinf   

      !update the position in time
      wake%w_points(:,iw,ipan) = point_old(:,iw,ipan-1) + vel_p*dt
! Check ----
      ! Prescribed rigid wake ----
!      wake%w_points(:,iw,ipan) = point_old(:,iw,ipan-1) + (/1.0_wp, 0.0_wp, 0.0_wp/)*dt
! Check ----
      
    enddo
  enddo
!$omp end parallel do

  deallocate(point_old)

 
  !==> 3) Increase the length of the wake, if it is necessary
  if (wake%wake_len .lt. wake%npan) then
      wake%wake_len = wake%wake_len + 1
      allocate(pan_p_temp(wake%n_wake_stripes*wake%wake_len))

      !debug: doing it all explicitly
      !if (wake%wake_len .gt. 2) then
      !  pan_p_temp(1:(wake%n_wake_stripes*(wake%wake_len-1))) = wake%pan_p
      !endif
      !do iw = 1,wake%n_wake_stripes
      !  pan_p_temp(wake%n_wake_stripes*(wake%wake_len-1)+iw)%p &
      !                                  => wake%wake_panels(iw,wake%wake_len)
      !enddo
      !call move_alloc(pan_p_temp,wake%pan_p)

      do ip = 1,size(wake%pan_p)
        pan_p_temp(ip) = wake%pan_p(ip)
      enddo
      do iw = 1,wake%n_wake_stripes
        pan_p_temp(wake%n_wake_stripes*(wake%wake_len-1)+iw)%p &
                                        => wake%wake_panels(iw,wake%wake_len)
      enddo
      if(allocated(wake%pan_p)) deallocate(wake%pan_p)
      allocate(wake%pan_p(size(pan_p_temp)))
      do ip = 1,size(wake%pan_p)
        wake%pan_p(ip) = pan_p_temp(ip)
      enddo
      deallocate(pan_p_temp)

  endif

  !==> 4) Update the intensities of the panels
  !       From the back, all the vortex intensities come from the previous panel
  do ipan = wake%wake_len,2,-1
    do iw = 1,wake%n_wake_stripes
      wake%wake_panels(iw,ipan)%idou = wake%wake_panels(iw,ipan-1)%idou
    enddo
  enddo

  ! The geometrical quantities of the panels will be all updated at the 
  ! beginning of the next iteration in prepare_wake


end subroutine update_wake_panels

!----------------------------------------------------------------------

end module mod_wake_pan
