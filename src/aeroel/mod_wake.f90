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


!> Module to treat the whole wake
module mod_wake

use mod_param, only: &
  wp, nl, pi, max_char_len

use mod_sim_param, only: &
  t_sim_param

use mod_handling, only: &
  error, warning, info, printout, dust_time, t_realtime

use mod_geometry, only: &
  t_geo, t_tedge, calc_geo_data_pan, calc_node_vel

!use mod_aero_elements, only: &
!  c_elem, t_elem_p

use mod_aeroel, only: &
  c_elem, c_pot_elem, c_vort_elem, c_impl_elem, c_expl_elem, &
  t_elem_p, t_pot_elem_p, t_vort_elem_p, t_impl_elem_p, t_expl_elem_p

use mod_surfpan, only: &
  t_surfpan

use mod_vortlatt, only: &
  t_vortlatt

use mod_vortline, only: &
  t_vortline

use mod_vortpart, only: &
  t_vortpart, t_vortpart_p

use mod_actuatordisk, only: &
  t_actdisk

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

use mod_octree, only: &
  t_octree, sort_particles, calculate_multipole, apply_multipole
!----------------------------------------------------------------------

implicit none

public :: t_wake, initialize_wake, update_wake, &
          prepare_wake, load_wake, destroy_wake

private


!> Type containing wake panels information
type :: t_wake

 !> Number of maximum streamwise panels
 integer :: nmax_pan

 !> Number of maximum "rows" of disks
 integer :: nmax_rin

 !> Number of actual streamwise panels
 integer :: pan_wake_len

 !> Actual number of rings in the ring wake
 integer :: rin_wake_len

 !> Number of wake stripes ("spanwise" panels)
 integer :: n_pan_stripes

 !> Number of wake points in the "spanwise" direction
 integer :: n_pan_points

 !> Number of generating disks
 integer :: ndisks

 !> Number of points for each "row"
 integer :: np_row

 !> Pointer and index of the 2 generating elements of the panel wake
 !! (2 x n_pan_stripes)
 type(t_pot_elem_p), allocatable :: pan_gen_elems(:,:)
 integer, allocatable :: pan_gen_elems_id(:,:)

 !> Generating actuator disk elements
 type(t_pot_elem_p), allocatable :: rin_gen_elems(:)

 !> Index of the 2 generating points of each wake point
 !! (2 x n_pan_points)
 integer, allocatable :: pan_gen_points(:,:)

 !> Direction of the wake at the trailing edge
 !! (3xn_pan_points)
 real(wp), allocatable :: pan_gen_dir(:,:)

 !> Reference frame of the generating points (id)
 !! (n_pan_points)
 integer, allocatable :: pan_gen_ref(:)

 !> Panels neighbours in wake numbering
 integer, allocatable :: pan_neigh(:,:)

 !> Relative orientation of neighbours
 integer, allocatable :: pan_neigh_o(:,:)

 !> Wake starting points: calculated from the 2 (or 1) starting points of the
 !! geometry possibly moved
 real(wp), allocatable :: w_start_points(:,:)

 !> Index of the 2 wake starting points for each wake stripe
 !! (2 x n_pan_stripes)
 integer, allocatable :: i_start_points(:,:)

 !> Points of the wake, in a structured way
 !! (3 x n_pan_points x npan+1)
 real(wp), allocatable :: pan_w_points(:,:,:)
 !> Velocities of the wake panels
 !! (3 x n_pan_points x npan+1) !!!! now used for rotational effects on pressure !!!!
 real(wp), allocatable :: pan_w_vel(:,:,:)

 !> Relative velocity ( u_inf - ub ) of the nodes at the TE
 !! (3 x n_pan_points)
 !! used to determine the first prescribed panel of the wake
 real(wp), allocatable :: w_vel_te(:,:)

 !> Velocity at the nodes of the wake. For output only
 !! (3 x n_pan_points x npan+1)
 real(wp), allocatable :: w_vel(:,:,:)

 !> elements of the panels
 type(t_vortlatt), allocatable :: wake_panels(:,:)

 !> Ring elements
 type(t_actdisk), allocatable :: wake_rings(:,:)

 !> end vortices
 type(t_vortline), allocatable :: end_vorts(:)

 !> doublets intensities
 real(wp), allocatable :: pan_idou(:,:)

 !> vortex intensities
 real(wp), allocatable :: rin_idou(:,:)

 !> pointer to the wake elements to be passed to the linsys
 !! solver
 type(t_pot_elem_p), allocatable :: pan_p(:)

 !> pointer to the wake elements to be passed to the linsys
 !! solver
 type(t_pot_elem_p), allocatable :: rin_p(:)


 !! Particles data

 !> Maximum number of particles
 integer :: nmax_prt

 !> Actual number of particles
 integer :: n_prt

 !> Wake particles
 type(t_vortpart), allocatable :: wake_parts(:)

 !> Magnitude of particles vorticity
 real(wp), allocatable :: prt_ivort(:)

 !> Velocity of the particles !!!! now used for rotational effects on pressure !!!!
 !real(wp), allocatable :: prt_vel(:,:)

 !> Wake particles pointer
 type(t_vortpart_p), allocatable :: part_p(:)

 !> Bounding box
 real(wp) :: part_box_min(3), part_box_max(3)

 type(t_vort_elem_p), allocatable :: vort_p(:)

 !> Last vortex intensity from removed panels
 real(wp), allocatable :: last_pan_idou(:)

 !> Are the panels full? (and so need to produce particles...)
 logical :: full_panels=.false.

 !> Are the rings full? (and so need to produce particles...)
 logical :: full_rings=.false.

end type

!module variables to share among the different subroutines
 real(wp), allocatable :: points_end(:,:)
 real(wp), allocatable :: points_end_ring(:,:)

!> Class to change methods from different wake implementations
type, abstract :: c_wake_mov
  contains
  procedure(i_get_vel), deferred, pass(this) :: get_vel
end type

abstract interface 
 subroutine i_get_vel(this, elems, wake, pos, sim_param, vel)
   import :: c_wake_mov, wp, t_pot_elem_p, t_wake, t_sim_param
   class(c_wake_mov) :: this
   type(t_pot_elem_p), intent(in) :: elems(:)
   type(t_wake), intent(in) :: wake
   real(wp), intent(in) :: pos(3)
   type(t_sim_param), intent(in) :: sim_param
   real(wp), intent(out) :: vel(3)
 end subroutine
end interface

type, extends(c_wake_mov) :: t_free_wake
contains
  procedure, pass(this) :: get_vel => get_vel_free
end type

type, extends(c_wake_mov) :: t_rigid_wake
contains
  procedure, pass(this) :: get_vel => get_vel_rigid
end type

class(c_wake_mov), allocatable :: wake_movement
character(len=max_char_len) :: msg
real(t_realtime) :: t1 , t0
character(len=*), parameter :: this_mod_name='mod_wake'

!----------------------------------------------------------------------
contains
!----------------------------------------------------------------------

!> Initialize the panel wake
subroutine initialize_wake(wake, geo, te,  npan, nrings, &
                           nparts, part_box_min, part_box_max, sim_param)
 type(t_wake), intent(out),target :: wake
 type(t_geo), intent(in), target :: geo
 type(t_tedge), intent(in) :: te
 integer, intent(in) :: npan
 integer, intent(in) :: nrings
 integer, intent(in) :: nparts
 real(wp), intent(in) :: part_box_min(3), part_box_max(3)
 type(t_sim_param), intent(in) :: sim_param

 integer :: iw, ip, nsides
 integer :: ic, nad, ie, npt, id, ir
 integer :: p1, p2
 real(wp) :: dist(3) , vel_te(3)

! real(wp) , parameter :: te_min_v = 1.0_wp ! hard-coded
! replaced with sim_param%min_vel_at_te

  if (sim_param%rigid_wake) then
    allocate(t_rigid_wake::wake_movement)
  else
    allocate(t_free_wake::wake_movement)
  endif
  !panel wake: set and allocate all the relevant variables
  wake%nmax_pan = npan
  wake%n_pan_stripes = size(te%e,2)
  wake%n_pan_points  = size(te%i,2)
  allocate(wake%pan_gen_elems(2,wake%n_pan_stripes))
  allocate(wake%pan_gen_elems_id(2,wake%n_pan_stripes))
  allocate(wake%pan_gen_points(2,wake%n_pan_points))
  allocate(wake%pan_gen_dir(3,wake%n_pan_points))
  allocate(wake%pan_gen_ref(wake%n_pan_points))
  allocate(wake%w_start_points(3,wake%n_pan_points))
  allocate(wake%i_start_points(2,wake%n_pan_stripes))
  allocate(wake%pan_neigh(2,wake%n_pan_stripes))
  allocate(wake%pan_neigh_o(2,wake%n_pan_stripes))
  allocate(wake%pan_w_points(3,wake%n_pan_points,npan+1))
  allocate(wake%pan_w_vel(   3,wake%n_pan_points,npan+1))
  allocate(wake%w_vel(3,wake%n_pan_points,npan+1))
  allocate(wake%wake_panels(wake%n_pan_stripes,npan))
  allocate(wake%end_vorts(wake%n_pan_stripes))
  allocate(wake%pan_idou(wake%n_pan_stripes,npan))

  !ring wake: count the number of actuator disks and points
  nad = 0; npt = 0;
  do ic = 1, size(geo%components)
    if(geo%components(ic)%comp_el_type(1:1) .eq. 'a') then
      nad = nad + geo%components(ic)%nelems
      do ie = 1, size(geo%components(ic)%el)
        npt = npt + geo%components(ic)%el(ie)%n_ver
      enddo
    endif
  enddo
  !ring wake: set and allocate all the relevant variables
  wake%nmax_rin = nrings
  wake%ndisks = nad
  wake%np_row = npt
  allocate(wake%rin_gen_elems(wake%ndisks))
  allocate(wake%wake_rings(wake%ndisks,wake%nmax_rin))
  allocate(wake%rin_idou(wake%ndisks,wake%nmax_rin))

  !panel wake: Associate for each panel the relevant intensity 
  ! and allocate all the relevant fields
  nsides = 4
  do ip = 1,npan
    do iw=1,wake%n_pan_stripes
     wake%wake_panels(iw,ip)%mag => wake%pan_idou(iw,ip)
     wake%wake_panels(iw,ip)%n_ver = nsides
     allocate(wake%wake_panels(iw,ip)%ver(3,nsides))
     allocate(wake%wake_panels(iw,ip)%edge_vec(3,nsides))
     allocate(wake%wake_panels(iw,ip)%edge_len(nsides))
     allocate(wake%wake_panels(iw,ip)%edge_uni(3,nsides))
    enddo
  enddo

  !ring wake: associate the generating elements
  id = 1
  do ic = 1, size(geo%components)
    if(geo%components(ic)%comp_el_type(1:1) .eq. 'a') then
      do ie = 1,size(geo%components(ic)%el)
        wake%rin_gen_elems(id)%p => geo%components(ic)%el(ie)
        id = id+1
      enddo
    endif
  enddo
  
  !ring wake: allocate the relevant fields of each ring
  do id = 1,wake%ndisks
    do ir = 1,wake%nmax_rin
      wake%wake_rings(id,ir)%mag => wake%rin_idou(id,ir)
      nsides = wake%rin_gen_elems(id)%p%n_ver
      wake%wake_rings(id,ir)%n_ver = nsides
      allocate(wake%wake_rings(id,ir)%ver(3,nsides))
      allocate(wake%wake_rings(id,ir)%edge_vec(3,nsides))
      allocate(wake%wake_rings(id,ir)%edge_len(nsides))
      allocate(wake%wake_rings(id,ir)%edge_uni(3,nsides))
    enddo
    wake%wake_rings(id,:)%moving = wake%rin_gen_elems(id)%p%moving
  enddo
  !Starting length of the ring wake is
  wake%rin_wake_len = 0
  allocate(wake%rin_p(0))

  !panel wake: set the generating elements and points from the trailing edge
  wake%pan_gen_elems  = te%e
  !do iw=1,size(te%e,2)
  !!TODO: consider shifting the whole t.e. to c_pot_elem to avoid this nightmare
  !select type(el => te%e(1,iw)%p); class is(c_pot_elem)
  !  wake%gen_elems(1,iw)%p => el
  !end select
  !select type(el => te%e(2,iw)%p); class is(c_pot_elem)
  !  wake%gen_elems(2,iw)%p => el
  !end select
  !enddo
  wake%pan_gen_points = te%i
  wake%pan_gen_dir = te%t
  wake%pan_gen_ref = te%ref
  wake%i_start_points = te%ii
  wake%pan_neigh = te%neigh
  wake%pan_neigh_o = te%o
  do iw=1,wake%n_pan_stripes
    wake%pan_gen_elems_id(1,iw) = wake%pan_gen_elems(1,iw)%p%id
    if(associated(wake%pan_gen_elems(2,iw)%p)) then
      wake%pan_gen_elems_id(2,iw) = wake%pan_gen_elems(2,iw)%p%id
    else
      wake%pan_gen_elems_id(2,iw) = 0
    endif
  enddo

  wake%pan_idou = 0.0_wp
  wake%rin_idou = 0.0_wp
  
  !panel wake: The whole first line of (implicit) panels is initialized 
  !here to be used during linear system initialization

  !the first line of points is calculated from the mesh points
  wake%w_start_points = 0.5_wp * (geo%points(:,wake%pan_gen_points(1,:)) + &
                                  geo%points(:,wake%pan_gen_points(2,:)))
  wake%pan_w_points(:,:,1) = wake%w_start_points
  wake%pan_w_vel(   :,:,:) = 0.0_wp 

  !Second row of points: first row + 0.3*|uinf|*t with t = R*t0
  do ip=1,wake%n_pan_points
!   dist = matmul(geo%refs(wake%pan_gen_ref(ip))%R_g,wake%pan_gen_dir(:,ip))
    call calc_node_vel( wake%w_start_points(:,ip), &
            geo%refs(wake%pan_gen_ref(ip))%G_g, &
            geo%refs(wake%pan_gen_ref(ip))%f_g, &
            vel_te )
    if ( norm2(sim_param%u_inf-vel_te) .gt. sim_param%min_vel_at_te ) then

      dist = matmul(geo%refs(wake%pan_gen_ref(ip))%R_g,wake%pan_gen_dir(:,ip))

      wake%pan_w_points(:,ip,2) = wake%pan_w_points(:,ip,1) +  &
                  dist*sim_param%first_panel_scaling* &
                  norm2(sim_param%u_inf-vel_te)*sim_param%dt / norm2(dist)
  ! normalisation occurs here! --------------------------------------^

    else

      dist = matmul(geo%refs(wake%pan_gen_ref(ip))%R_g,wake%pan_gen_dir(:,ip))

      wake%pan_w_points(:,ip,2) = wake%pan_w_points(:,ip,1) +  &
                  dist*sim_param%first_panel_scaling * & ! next line may be commented
                  sim_param%min_vel_at_te*sim_param%dt
    end if

  enddo

  !Starting length of the panel wake is 1
  wake%pan_wake_len = 1
  
  !Associate the first implicit row
  allocate(wake%pan_p(wake%n_pan_stripes))
  do iw = 1,wake%n_pan_stripes
    wake%wake_panels(iw,:)%moving = geo%refs(wake%pan_gen_ref(iw))%moving
    wake%pan_p(iw)%p => wake%wake_panels(iw,1)
  enddo

  !set the end vortexes to null
  do iw = 1,wake%n_pan_stripes
    wake%end_vorts(iw)%mag => null()
  enddo
  
  ! Calculate geometrical quantities of first row
  do iw = 1,wake%n_pan_stripes
    p1 = wake%i_start_points(1,iw)
    p2 = wake%i_start_points(2,iw)
    call wake%wake_panels(iw,1)%calc_geo_data( &
         reshape((/wake%pan_w_points(:,p1,1),   wake%pan_w_points(:,p2,1), &
                   wake%pan_w_points(:,p2,1+1), wake%pan_w_points(:,p1,1+1)/),&
                                                                   (/3,4/)))
  enddo


  !Particles
  wake%nmax_prt = nparts
  allocate(wake%wake_parts(wake%nmax_prt))
  allocate(wake%prt_ivort(wake%nmax_prt))
  wake%n_prt = 0
  allocate(wake%part_p(0))
  do ip = 1,wake%nmax_prt
    wake%wake_parts(ip)%mag => wake%prt_ivort(ip)
  enddo
  wake%part_box_min = part_box_min
  wake%part_box_max = part_box_max

  allocate(wake%vort_p(0))
  allocate(wake%last_pan_idou(wake%n_pan_stripes))
  wake%last_pan_idou = 0.0_wp

  wake%full_panels = .false.

end subroutine initialize_wake

!----------------------------------------------------------------------

!> Destroy a wake panels type by simply passing it as intent(out)
subroutine destroy_wake(wake)
 type(t_wake), intent(out) :: wake

 !dummy to avoid compiler warnings
 wake%nmax_pan = -1
 if(allocated(points_end)) deallocate(points_end)
 if(allocated(points_end_ring)) deallocate(points_end_ring)

end subroutine

!----------------------------------------------------------------------

!> Prepare the first row of panels to be inserted inside the linear system
!!
subroutine prepare_wake(wake, geo, sim_param)
 type(t_wake), target, intent(inout) :: wake
 type(t_geo), intent(in) :: geo
 type(t_sim_param), intent(in) :: sim_param

 integer :: p1, p2
 integer :: ip, iw, ipan, id, is, nprev
 real(wp) :: dist(3) , vel_te(3), pos_p(3)
 real(wp) :: dir(3), partvec(3), ave
 integer :: ir, k

 ! flow separation variables
 integer :: i_comp , i_elem , n_elem
 integer :: n_end_vort   ! for structure update (around l.750)

 character(len=max_char_len) :: msg
 character(len=*), parameter :: this_sub_name='prepare_wake'

  !==> Panels: update the first rows of panels

  !first row  of new points comes from geometry
  wake%w_start_points = 0.5_wp * (geo%points(:,wake%pan_gen_points(1,:)) + &
                                  geo%points(:,wake%pan_gen_points(2,:)))
  wake%pan_w_points(:,:,1) = wake%w_start_points

  !Second row of points: first row + 0.3*|uinf|*t with t = R*t0
  do ip=1,wake%n_pan_points
    dist = matmul(geo%refs(wake%pan_gen_ref(ip))%R_g,wake%pan_gen_dir(:,ip))
    call calc_node_vel( wake%w_start_points(:,ip), &
            geo%refs(wake%pan_gen_ref(ip))%G_g, &
            geo%refs(wake%pan_gen_ref(ip))%f_g, &
            vel_te )
    if ( norm2(sim_param%u_inf-vel_te) .gt. sim_param%min_vel_at_te ) then
      wake%pan_w_points(:,ip,2) = wake%pan_w_points(:,ip,1) + &
                          dist*sim_param%first_panel_scaling* &
                          norm2(sim_param%u_inf-vel_te)*sim_param%dt / norm2(dist)
  ! normalisation occurs here! -------------------------------------------^
    else
      wake%pan_w_points(:,ip,2) = wake%pan_w_points(:,ip,1) +  &
                  dist*sim_param%first_panel_scaling * & ! next line may be commented
                  sim_param%min_vel_at_te*sim_param%dt / norm2(dist)
    end if
  enddo

  !==> Panels:  update the panels geometrical quantities of all the panels, 
  !      the first two row of points have just been updated, the other 
  !      rows of points were updated at the end of the last iteration
  do ipan = 1,wake%pan_wake_len
    do iw = 1,wake%n_pan_stripes
      p1 = wake%i_start_points(1,iw)
      p2 = wake%i_start_points(2,iw)
      call wake%wake_panels(iw,ipan)%calc_geo_data( &
          reshape((/wake%pan_w_points(:,p1,ipan),   wake%pan_w_points(:,p2,ipan), &
                    wake%pan_w_points(:,p2,ipan+1), wake%pan_w_points(:,p1,ipan+1)/),&
                                                                     (/3,4/)))
    enddo
  enddo

  !==> Particles: if the panel wake is at the end, create a particle
  if(wake%full_panels) then
    k = 1
    do iw = 1,wake%n_pan_stripes
      p1 = wake%i_start_points(1,iw)
      p2 = wake%i_start_points(2,iw)
      partvec = 0.0_wp
      !Left side
      dir = wake%pan_w_points(:,p1,wake%nmax_pan+1)-points_end(:,p1)
      if (wake%pan_neigh(1,iw) .gt. 0) then
        ave = wake%wake_panels(iw,wake%pan_wake_len)%mag - &
              wake%pan_neigh_o(1,iw)* &
              wake%wake_panels(wake%pan_neigh(1,iw),wake%pan_wake_len)%mag
        ave = ave/2.0_wp
      else
        ave = wake%wake_panels(iw,wake%pan_wake_len)%mag
      endif
      partvec = partvec + dir*ave

      !Right side
      dir = -wake%pan_w_points(:,p2,wake%nmax_pan+1)+points_end(:,p2)
      if (wake%pan_neigh(2,iw) .gt. 0) then
        ave = wake%wake_panels(iw,wake%pan_wake_len)%mag - &
              wake%pan_neigh_o(2,iw)*wake%wake_panels(wake%pan_neigh(2,iw),wake%pan_wake_len)%mag
        ave = ave/2.0_wp
      else
        ave = wake%wake_panels(iw,wake%pan_wake_len)%mag
      endif
      partvec = partvec + dir*ave

      !End side
      dir = points_end(:,p1) - points_end(:,p2)
      ave = wake%wake_panels(iw,wake%pan_wake_len)%mag-wake%last_pan_idou(iw)
      wake%last_pan_idou(iw) = wake%wake_panels(iw,wake%pan_wake_len)%mag
      partvec = partvec + dir*ave

      !Calculate the center
      pos_p = (points_end(:,p1)+points_end(:,p2)+ &
              wake%pan_w_points(:,p1,wake%nmax_pan+1) + &
              wake%pan_w_points(:,p2,wake%nmax_pan+1) )/4.0_wp

      !pos_p = (1.5_wp*points_end(:,p1)+1.5_wp*points_end(:,p2)+ &
      !        wake%pan_w_points(:,p1,wake%nmax_pan+1) + &
      !        wake%pan_w_points(:,p2,wake%nmax_pan+1) )/5.0_wp

      !pos_p = (points_end(:,p1)+points_end(:,p2))/2.0_wp

      !Add the particle (if it is in the box)
      if(all(pos_p .ge. wake%part_box_min) .and. &
         all(pos_p .le. wake%part_box_max)) then
        do ip = k, size(wake%wake_parts)
          if (wake%wake_parts(ip)%free) then
            
            wake%wake_parts(ip)%free = .false.
            k = ip+1
            wake%n_prt = wake%n_prt+1
            wake%wake_parts(ip)%mag = norm2(partvec)
            if(wake%wake_parts(ip)%mag .gt. 1.0e-13_wp) then
              wake%wake_parts(ip)%dir = partvec/wake%wake_parts(ip)%mag
            else
              wake%wake_parts(ip)%dir = partvec
            endif
            wake%wake_parts(ip)%cen = pos_p

            wake%wake_parts(ip)%vel = 0.5_wp * &
                             ( wake%pan_w_vel(:,p1,wake%nmax_pan+1) + &
                               wake%pan_w_vel(:,p2,wake%nmax_pan+1) )
            exit
          endif
        enddo
        if (ip .gt. wake%nmax_prt) then
          write(msg,'(A,I0,A)') 'Exceeding the maximum number of ', &
            wake%nmax_prt, ' wake particles introduced. Stopping. Consider &
            &restarting with a higher number of maximum wake particles'
        call error(this_sub_name, this_mod_name, trim(msg))
        endif !max number of particles
      endif !inside the box
      
    enddo !iw
  endif !full wake 

  if(wake%full_rings) then
    k = 1
    nprev = 0
    do id = 1,wake%ndisks
      do is = 1,wake%wake_rings(id,wake%rin_wake_len)%n_ver
        p1 = is + nprev
        p2 = nprev + 1 + mod(is,wake%wake_rings(id,wake%rin_wake_len)%n_ver)
        partvec = 0.0_wp

        dir = points_end_ring(:,p2)-points_end_ring(:,p1)
        ave = wake%wake_rings(id,wake%rin_wake_len)%mag
        partvec =dir*ave

        !Calculate the center
        pos_p = (points_end_ring(:,p1)+points_end_ring(:,p2))/2.0_wp

        !Add the particle
        if(all(pos_p .ge. wake%part_box_min) .and. &
           all(pos_p .le. wake%part_box_max)) then
          do ip = k, size(wake%wake_parts)
            if (wake%wake_parts(ip)%free) then
              
              wake%wake_parts(ip)%free = .false.
              k = ip+1
              wake%n_prt = wake%n_prt+1
              wake%wake_parts(ip)%mag = norm2(partvec)
              if(wake%wake_parts(ip)%mag .gt. 1.0e-13_wp) then
                wake%wake_parts(ip)%dir = partvec/wake%wake_parts(ip)%mag
              else
                wake%wake_parts(ip)%dir = partvec
              endif
              wake%wake_parts(ip)%cen = pos_p
              exit
            endif
          enddo  !ip
          if (ip .gt. wake%nmax_prt) then
            write(msg,'(A,I0,A)') 'Exceeding the maximum number of ', &
              wake%nmax_prt, ' wake particles introduced. Stopping. Consider &
              &restarting with a higher number of maximum wake particles'
          call error(this_sub_name, this_mod_name, trim(msg))
          endif
        endif
      enddo !is
      nprev = nprev + wake%wake_rings(id,wake%rin_wake_len)%n_ver
    enddo !id
  endif


  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Treat flow separations 
  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if(sim_param%use_ve) then
    do i_comp = 1 , size( geo%components)! ***** loop #1 over components *****

     if ( k .lt. 1 ) k = 1

     ! flow separation allowed only for surfpan elements -----
      if ( trim( geo%components(i_comp)%comp_el_type ) .eq. 'p' ) then

        n_elem = size( geo%components(i_comp)%el )

        do i_elem = 1 , n_elem     ! ***** loop #2 over elements   *****

          select type( el => geo%components(i_comp)%el(i_elem) ) 
          type is (t_surfpan)

          ! flow separation
          if ( el % al_free .gt. 0.0_wp ) then 
            pos_p = el%cen + el%nor * el%h_bl + & 
                       el % surf_vel * sim_param%dt

            if(all(pos_p .ge. wake%part_box_min) .and. &
               all(pos_p .le. wake%part_box_max)) then
              !Add the particle
              do ip = k, size(wake%wake_parts)

                if (wake%wake_parts(ip)%free) then
                  
                  wake%wake_parts(ip)%free = .false.
                  k = ip+1
                  wake%n_prt = wake%n_prt+1
                  wake%wake_parts(ip)%mag = norm2(el%free_vort)
                  if(wake%wake_parts(ip)%mag .gt. 1.0e-13_wp) then
                    wake%wake_parts(ip)%dir = el%free_vort/wake%wake_parts(ip)%mag
                  else
                    wake%wake_parts(ip)%dir = el%free_vort
                  endif
                  !wake%wake_parts(ip)%cen = el%cen + el%nor * el%h_bl + & 
                  !       el % surf_vel * sim_param%dt
                  wake%wake_parts(ip)%cen = pos_p
                  wake%wake_parts(ip)%vel = el % surf_vel
                  exit 
                endif
              enddo
              if (ip .gt. wake%nmax_prt) then
                write(msg,'(A,I0,A)') 'Exceeding the maximum number of ', &
                  wake%nmax_prt, ' wake particles introduced. Stopping. Consider &
                  &restarting with a higher number of maximum wake particles'
              call error(this_sub_name, this_mod_name, trim(msg))
              endif
            endif !particle in box
          end if !generate the particle
          end select!select panels
        end do ! ***** loop #2 over elements   *****
      end if !if is a 3dp
    end do ! ***** loop #1 over components *****
  endif !calculate the viscous particles detachment

  ! Recreate sturctures and pointers, if full
  if(wake%full_panels .or. wake%full_rings .or. (wake%n_prt.gt.0) ) then
    !Recreate the pointer vector
    if(allocated(wake%part_p)) deallocate(wake%part_p)
    allocate(wake%part_p(wake%n_prt))
    deallocate(wake%vort_p)
    ! to add or not line vortices at the (only when ring or panel wakes are full )
    n_end_vort = 0
    if ( wake%full_panels .or. wake%full_rings ) then
      n_end_vort = wake%n_pan_stripes
    end if
    allocate(wake%vort_p(wake%n_prt + n_end_vort))
    !TODO: consider inverting these two cycles
    k = 1
    do ip = 1, wake%n_prt
      do ir=k,wake%nmax_prt
        if(.not. wake%wake_parts(ir)%free) then
          k = ir+1
          wake%part_p(ip)%p => wake%wake_parts(ir)
          wake%vort_p(ip)%p => wake%wake_parts(ir)
          exit
        endif
      enddo
    enddo
    !Add the end vortex to the votical elements pointer
!   if ( wake%full_panels .or. wake%full_rings ) then ! useless if ( n_end_vort may be 0 )
    do iw = 1, n_end_vort
      wake%vort_p(wake%n_prt+iw)%p => wake%end_vorts(iw)
    enddo
!   end if
  endif

  !If the wake is full, attach the end vortex
  if (wake%full_panels) then
    do iw = 1,wake%n_pan_stripes
      wake%end_vorts(iw)%mag => wake%wake_panels(iw,wake%pan_wake_len)%mag
      p1 = wake%i_start_points(1,iw)
      p2 = wake%i_start_points(2,iw)
      call wake%end_vorts(iw)%calc_geo_data( &
          reshape((/wake%pan_w_points(:,p1,wake%pan_wake_len+1),  &
                    wake%pan_w_points(:,p2,wake%pan_wake_len+1)/), (/3,2/)))
      wake%end_vorts(iw)%ver_vel(:,1) = wake%pan_w_vel(:,p1,wake%pan_wake_len+1)
      wake%end_vorts(iw)%ver_vel(:,2) = wake%pan_w_vel(:,p2,wake%pan_wake_len+1) 
    enddo
  endif

end subroutine prepare_wake

!----------------------------------------------------------------------

!> Load the wake panels solution from a previous result
subroutine load_wake(filename, wake)
 character(len=*), intent(in) :: filename
 type(t_wake), intent(inout), target :: wake

 integer(h5loc) :: floc, gloc
 real(wp), allocatable :: wpoints(:,:,:), wvels(:,:,:), wvort(:,:)
 real(wp), allocatable :: vppoints(:,:), vpvort(:,:) , vpvels(:,:)
 integer, allocatable :: start_points(:,:)
 integer, allocatable :: conn_pe(:)
 integer :: ipan, iw, p1, p2, ipt
 integer :: id, ir, ip, np
 character(len=*), parameter :: this_sub_name = 'load_wake'
   
  call open_hdf5_file(filename, floc)

  !=== Panels ===
  !Read the past results
  call open_hdf5_group(floc, 'PanelWake', gloc)

  call read_hdf5_al(wpoints,'WakePoints',gloc)
  call read_hdf5_al(wvels  ,'WakeVels'  ,gloc)
  call read_hdf5_al(start_points,'StartPoints',gloc)
  call read_hdf5_al(wvort,'WakeVort',gloc)

  call close_hdf5_group(gloc)
  !call close_hdf5_file(floc)

  !Perform a few checks to be sure that the correct size solution was loaded
  if(.not. all(start_points .eq. wake%i_start_points)) call error( &
    this_sub_name, this_mod_name, 'Different wake trailing edge connectivity&
    & between the loded and built geometry')

  !the wake length is the loaded one, or less if imposed so
  wake%pan_wake_len = min(wake%nmax_pan, size(wvort,2))

  !store points position and doublets intensity
  wake%pan_w_points(:,:,1:wake%pan_wake_len+1) = wpoints(:,:,1:wake%pan_wake_len+1)
  wake%pan_w_vel(   :,:,1:wake%pan_wake_len+1) = wvels(  :,:,1:wake%pan_wake_len+1)
  wake%pan_idou(:,1:wake%pan_wake_len) = wvort(:,1:wake%pan_wake_len)

  deallocate(wake%pan_p)
  allocate(wake%pan_p(wake%n_pan_stripes*wake%pan_wake_len))
  ! Update the panels geometrical quantities of all the panels
  ipt = 1
  do ipan = 1,wake%pan_wake_len
    do iw = 1,wake%n_pan_stripes
      p1 = wake%i_start_points(1,iw)
      p2 = wake%i_start_points(2,iw)
      call wake%wake_panels(iw,ipan)%calc_geo_data( &
           reshape((/wake%pan_w_points(:,p1,ipan),   wake%pan_w_points(:,p2,ipan), &
                    wake%pan_w_points(:,p2,ipan+1), wake%pan_w_points(:,p1,ipan+1)/),&
                                                                     (/3,4/)))
      wake%pan_p(ipt)%p => wake%wake_panels(iw,ipan)
      ipt = ipt + 1
    enddo
  enddo

  deallocate(wpoints, wvels, wvort, start_points)

  !=== Rings ===
  !Read the past results
  !call open_hdf5_file(filename, floc)
  call open_hdf5_group(floc, 'RingWake', gloc)

  call read_hdf5_al(wpoints,'WakePoints',gloc)
  call read_hdf5_al(conn_pe,'Conn_pe',gloc)
  call read_hdf5_al(wvort,'WakeVort',gloc)

  call close_hdf5_group(gloc)
  !call close_hdf5_file(floc)

  !check the consistency
  ip = 1
  do id = 1,wake%ndisks
    np = wake%rin_gen_elems(id)%p%n_ver
    if (.not. all(conn_pe(ip:ip+np-1) .eq. id)) call error( &
    this_sub_name, this_mod_name, 'Different wake disk connectivity&
    & between the loded and built geometry')
    ip = ip+np
  enddo

  !the wake length is the loaded one, or less if imposed so
  wake%rin_wake_len = min(wake%nmax_rin, size(wvort,2))

  !Load the old results
  ip = 1
  do id = 1,wake%ndisks
    np = wake%rin_gen_elems(id)%p%n_ver
    do ir = 1,wake%rin_wake_len
      wake%wake_rings(id,ir)%ver(:,:) = wpoints(:,ip:ip+np-1,ir)
    enddo
    ip = ip+np
  enddo
  wake%rin_idou(:,1:wake%rin_wake_len) = wvort(:,1:wake%rin_wake_len)

  deallocate(wake%rin_p); allocate(wake%rin_p(wake%ndisks*wake%rin_wake_len))

  !Update the geometrical quantities
  ipt = 1
  do ir = 1,wake%rin_wake_len
    do id = 1,wake%ndisks
      call wake%wake_rings(id,ir)%calc_geo_data( &
                    wake%wake_rings(id,ir)%ver)
      wake%rin_p(ipt)%p => wake%wake_rings(id,ir)
      ipt = ipt + 1
    enddo
  enddo

  !=== Particles ===
  call open_hdf5_group(floc, 'ParticleWake', gloc)
  call read_hdf5_al(vppoints,'WakePoints',gloc)
  if ( check_dset_hdf5('WakeVels',gloc) ) then
    call read_hdf5_al(vpvels,'WakeVels',gloc) ! <<<< restart with Bernoulli integral equation
  end if
  call read_hdf5_al(vpvort,'WakeVort',gloc)
  call read_hdf5(wake%last_pan_idou,'LastPanIdou',gloc)
  call close_hdf5_group(gloc)
  
  if(size(vpvort,2) .gt. wake%nmax_prt) call error(this_sub_name, &
    this_mod_name, 'Loading a number of particles &
    & greater than the maximum allowed for this run: cannot truncate the &
    & particles in a meaningful way. Consider restarting the run with a &
    & higher amount of maximum particles')

  wake%n_prt = size(vpvort,2)

  deallocate(wake%part_p)
  allocate(wake%part_p(wake%n_prt))
  if(wake%n_prt .gt. 0) then
    deallocate(wake%vort_p)
    allocate(wake%vort_p(wake%n_prt+wake%n_pan_stripes))
    do iw = 1, wake%n_pan_stripes 
      wake%vort_p(wake%n_prt+iw)%p => wake%end_vorts(iw)
    enddo
  endif

  do ip = 1,wake%n_prt
    wake%wake_parts(ip)%cen = vppoints(:,ip)
    wake%wake_parts(ip)%vel = vpvels(:,ip)
    wake%wake_parts(ip)%mag = norm2(vpvort(:,ip))
    if(wake%wake_parts(ip)%mag .gt. 1.0e-13_wp) then
      wake%wake_parts(ip)%dir = vpvort(:,ip)/wake%wake_parts(ip)%mag
    else
      wake%wake_parts(ip)%dir = vpvort(:,ip)
    endif
    wake%wake_parts(ip)%free = .false.
    wake%part_p(ip)%p => wake%wake_parts(ip)
    wake%vort_p(ip)%p => wake%wake_parts(ip)
  enddo
  
  deallocate(vppoints, vpvort, vpvels)
  
  call close_hdf5_file(floc)


  if(wake%pan_wake_len .eq. wake%nmax_pan) wake%full_panels=.true.
  !If the wake is full, attach the end vortex
  if (wake%full_panels) then
    do iw = 1,wake%n_pan_stripes
      wake%end_vorts(iw)%mag => wake%wake_panels(iw,wake%pan_wake_len)%mag
      p1 = wake%i_start_points(1,iw)
      p2 = wake%i_start_points(2,iw)
      call wake%end_vorts(iw)%calc_geo_data( &
          reshape((/wake%pan_w_points(:,p1,wake%pan_wake_len+1),  &
                    wake%pan_w_points(:,p2,wake%pan_wake_len+1)/), (/3,2/)))
    enddo
  endif
end subroutine load_wake


!----------------------------------------------------------------------

!> Update the position and the intensities of the wake panels
!!
!! Note: at this subroutine is passed the whole array of elements,
!! comprising both the implicit panels and the explicit (ll)
!! elements
subroutine update_wake(wake, elems, octree, sim_param)
 type(t_wake), intent(inout), target :: wake
 type(t_pot_elem_p), intent(in) :: elems(:)
 type(t_octree), intent(inout) :: octree
 type(t_sim_param), intent(in) :: sim_param

 integer :: iw, ipan, ie, ip, np, iq
 integer :: id, ir, n_part
 real(wp) :: pos_p(3), vel_p(3), alpha_p(3)
 real(wp) :: str(3), stretch(3)
 real(wp) :: df(3), diff(3)
 type(t_pot_elem_p), allocatable :: pan_p_temp(:)
 real(wp), allocatable :: point_old(:,:,:)
 real(wp), allocatable :: points(:,:,:)
 real(wp), allocatable, target :: vortevol_prt(:,:)
 logical :: increase_wake
 integer :: size_old
 character(len=*), parameter :: this_sub_name='update_wake'

  wake%w_vel = 0.0_wp

  if(wake%pan_wake_len .eq. wake%nmax_pan) wake%full_panels=.true.

  !==> Panels:  Update the first row of vortex intensities:
  !      it was already calculated (implicitly) in the linear system
  do iw = 1,wake%n_pan_stripes
    !
    if      ( associated(wake%pan_gen_elems(2,iw)%p) ) then
      wake%wake_panels(iw,1)%mag  = wake%pan_gen_elems(1,iw)%p%mag - &
                                     wake%pan_gen_elems(2,iw)%p%mag
    else if ( .not. associated(wake%pan_gen_elems(2,iw)%p) ) then
      wake%wake_panels(iw,1)%mag  = wake%pan_gen_elems(1,iw)%p%mag
    end if

  enddo

  !==> Panels: Update wake points position ==

  !Save the old positions for the integration
  allocate(point_old(size(wake%pan_w_points,1),size(wake%pan_w_points,2), &
                                               size(wake%pan_w_points,3)))
  point_old = wake%pan_w_points

  !calculate the velocities at the old positions of the points and then
  !update the positions (from the third row of points: the first is the
  !trailing edge, the second is extrapolated from the trailing edge)
  np = wake%pan_wake_len+1
  !if(wake%pan_wake_len .lt. wake%nmax_pan) np = np + 1
  if(.not.wake%full_panels) np = np + 1

!$omp parallel do collapse(2) private(pos_p, vel_p, ie, ipan, iw) schedule(dynamic)
  do ipan = 3,np
    do iw = 1,wake%n_pan_points
      pos_p = point_old(:,iw,ipan-1)
      vel_p = 0.0_wp

      !call compute_vel_from_all(elems, wake, pos_p, sim_param, vel_p)

      !! for OUTPUT only -----
      !wake%w_vel(:,iw,ipan-1) = vel_p

      !vel_p    = vel_p   + sim_param%u_inf
      call wake_movement%get_vel(elems, wake, pos_p, sim_param, vel_p)

      !update the position in time
      wake%pan_w_vel(   :,iw,ipan) = vel_p
      wake%pan_w_points(:,iw,ipan) = point_old(:,iw,ipan-1) + vel_p*sim_param%dt
    enddo
  enddo
!$omp end parallel do
  
  !if the wake is full, calculate another row of points to generate the 
  !particles
  if(wake%full_panels) then
    if(.not.allocated(points_end)) allocate(points_end(3,wake%n_pan_points)) 

    ! create another row of points
!$omp parallel do private(iw, pos_p, vel_p) schedule(dynamic)
    do iw = 1,wake%n_pan_points
      pos_p = point_old(:,iw,wake%pan_wake_len+1)
      vel_p = 0.0_wp

      !call compute_vel_from_all(elems, wake, pos_p, sim_param, vel_p)

      !vel_p    = vel_p   + sim_param%u_inf
      call wake_movement%get_vel(elems, wake, pos_p, sim_param, vel_p)

      !update the position in time
      points_end(:,iw) = pos_p + vel_p*sim_param%dt
    enddo
!$omp end parallel do
  endif
  

  deallocate(point_old)

  !==> Rings: Update wake points position ==

  !Save the old positions for the integration
  increase_wake = .false.
  if(wake%rin_wake_len .lt. wake%nmax_rin) then
    wake%rin_wake_len = wake%rin_wake_len+1
    increase_wake = .true.
  else
    wake%full_rings = .true.
  endif
  allocate(points(3,wake%np_row,wake%rin_wake_len))

  !Store at the beginning the disk points
  ip=1
  do id = 1,wake%ndisks
    np = wake%rin_gen_elems(id)%p%n_ver
    points(:,ip:ip+np-1,1) = wake%rin_gen_elems(id)%p%ver
    ip = ip+np
  enddo

  !Then store the old points of the rest of the wake (the shift forward
  ! of the points in the array is happening here)
  do ir = 1,wake%rin_wake_len-1
    ip=1
    do id = 1,wake%ndisks
      np = wake%wake_rings(id,ir)%n_ver
      points(:,ip:ip+np-1,ir+1) = wake%wake_rings(id,ir)%ver(:,:)
      ip = ip+np
    enddo
  enddo

  !calculate the velocities at the old positions of the points and then
  !update the positions

!$omp parallel do private(pos_p, vel_p, ip)
  do ip = 1,size(points,2)
    do ir = 1,size(points,3)
      pos_p = points(:,ip,ir)
      vel_p = 0.0_wp

      call wake_movement%get_vel(elems, wake, pos_p, sim_param, vel_p)

      !update the position in time
      points(:,ip,ir) = points(:,ip,ir) + vel_p*sim_param%dt
    enddo !ir
  enddo !ip
!$omp end parallel do

  !if the wake is full, calculate another row of points to generate the 
  !particles
  if(wake%full_rings) then
    if(.not.allocated(points_end_ring)) allocate(points_end_ring(3,wake%np_row)) 
    !Store the end points
    ip=1
    do id = 1,wake%ndisks
      np = wake%wake_rings(id,ir)%n_ver
      points_end_ring(:,ip:ip+np-1) = wake%wake_rings(id,ir)%ver(:,:)
      ip = ip+np
    enddo
    ! calculate velocity and evolve them
    do ip = 1,wake%np_row
      pos_p = points_end_ring(:,ip)
      vel_p = 0.0_wp

      !call compute_vel_from_all(elems, wake, pos_p, sim_param, vel_p)

      !vel_p    = vel_p   + sim_param%u_inf
      call wake_movement%get_vel(elems, wake, pos_p, sim_param, vel_p)

      !update the position in time
      points_end_ring(:,ip) = pos_p + vel_p*sim_param%dt
    enddo
  endif

  !redistribute the points
  do ir = 1,wake%rin_wake_len
    ip = 1
    do id = 1,wake%ndisks
      np = wake%wake_rings(id,ir)%n_ver
      wake%wake_rings(id,ir)%ver = points(:,ip:ip+np-1,ir)
      ip = ip + np
    enddo
  enddo

  deallocate(points)


  !==>    Particles: evolve the position in time

  allocate(vortevol_prt(3,wake%n_prt))
  vortevol_prt = 0.0_wp

  !if ( allocated(wake%prt_vel) ) deallocate(wake%prt_vel)
  !allocate(wake%prt_vel(3,wake%n_prt))


  !calculate the velocities at the points
!$omp parallel do private(pos_p, vel_p, ip, iq,  stretch, diff, df, str)
  do ip = 1, wake%n_prt

    if (sim_param%use_vs .or. sim_param%use_vd) then
      wake%part_p(ip)%p%stretch => vortevol_prt(:,ip)
    endif
    
    !If not using the fast multipole, update particles position now
    if (.not.sim_param%use_fmm) then
      pos_p = wake%part_p(ip)%p%cen

      call wake_movement%get_vel(elems, wake, pos_p, sim_param, vel_p)

      !wake%prt_vel(:,ip) = vel_p                            ! *****(old) vel_prt(:,ip) = vel_p
      !wake%part_p(ip)%p%vel     =  wake%prt_vel(:,ip)       ! *****(old) vel_prt(:,ip)
      wake%part_p(ip)%p%vel     =  vel_p

      !if using vortex stretching, calculate it now
      if(sim_param%use_vs) then
        stretch = 0.0_wp
        do iq = 1, wake%n_prt
        if (ip.ne.iq) then
          call wake%part_p(iq)%p%compute_stretch(wake%part_p(ip)%p%cen, &
               wake%part_p(ip)%p%dir*wake%part_p(ip)%p%mag, str)
          !stretch = stretch + str/(4.0_wp*pi)
           stretch = stretch +(str - &
           sum(str*wake%part_p(ip)%p%dir)*wake%part_p(ip)%p%dir)/(4.0_wp*pi)
          !removed the parallel component
        endif 
        enddo
        !do ie=1,size(wake%end_vorts)
        !  call wake%end_vorts(ie)%compute_stretch(wake%part_p(ip)%p%cen, &
        !             wake%part_p(ip)%p%dir*wake%part_p(ip)%p%mag, str)
        !  stretch = stretch + str/(4.0_wp*pi)
        !enddo
        vortevol_prt(:,ip) = vortevol_prt(:,ip) + stretch
      endif !use_vs

      !if using the vortex diffusion, calculate it now
      if(sim_param%use_vd) then
        diff = 0.0_wp
        do iq = 1, wake%n_prt
        if (ip.ne.iq) then
          call wake%part_p(iq)%p%compute_diffusion(wake%part_p(ip)%p%cen, &
               wake%part_p(ip)%p%dir*wake%part_p(ip)%p%mag, df)
          diff = diff + df*sim_param%nu_inf
        endif 
        enddo !iq
        vortevol_prt(:,ip) = vortevol_prt(:,ip) + diff
      endif !use_vd
    end if !use_fmm

    
  enddo
!$omp end parallel do
  
  if (sim_param%use_fmm) then
    t0 = dust_time()
    call sort_particles(wake%part_p, wake%n_prt, octree, sim_param)
    call calculate_multipole(wake%part_p, octree)
    call apply_multipole(wake%part_p, octree, elems, wake%pan_p, wake%rin_p, &
                         wake%end_vorts, sim_param)
    t1 = dust_time()
    write(msg,'(A,F9.3,A)') 'Multipoles calculation: ' , t1 - t0,' s.'
    call printout(msg)
    write(msg,'(A,I0)') 'Number of particles: ' , wake%n_prt
    call printout(msg)
  endif

  !Check the difference
  !err = norm2(points_prt-points_prt_fmm)/norm2(points_prt)
  !write(*,*) 'error',err

  !Assign the moved points, if they get outside the bounding box free the 
  !particles
  n_part = wake%n_prt
  do ip = 1, n_part
    if(sim_param%use_pa) call avoid_collision(elems, wake, &
                        wake%part_p(ip)%p, sim_param, wake%part_p(ip)%p%vel)
                           !wake%part_p(ip)%p, sim_param, wake%prt_vel(:,ip))
    if(.not. wake%part_p(ip)%p%free) then
      !pos_p = wake%part_p(ip)%p%cen + wake%prt_vel(:,ip)*sim_param%dt
      pos_p = wake%part_p(ip)%p%cen + wake%part_p(ip)%p%vel*sim_param%dt
      if(all(pos_p .ge. wake%part_box_min) .and. &
         all(pos_p .le. wake%part_box_max)) then
        !wake%part_p(ip)%p%cen = points_prt(:,ip)
        wake%part_p(ip)%p%cen = pos_p
        if(sim_param%use_vs .or. sim_param%use_vd) then
          alpha_p = wake%part_p(ip)%p%dir*wake%part_p(ip)%p%mag + &
                          vortevol_prt(:,ip)*sim_param%dt
          wake%part_p(ip)%p%mag = norm2(alpha_p)
          if(wake%part_p(ip)%p%mag .ne. 0.0_wp) &
             wake%part_p(ip)%p%dir = alpha_p/wake%part_p(ip)%p%mag
        endif
      else
        wake%part_p(ip)%p%free = .true.
        wake%n_prt = wake%n_prt -1
      endif
    endif
    !nullify(wake%part_p(ip)%p%npos)
    !nullify(wake%part_p(ip)%p%vel)
    if(sim_param%use_vs) nullify(wake%part_p(ip)%p%stretch)
  enddo

  !==> Panels:  Increase the length of the wake, if it is necessary
  !if (wake%pan_wake_len .lt. wake%nmax_pan) then
  if (.not.wake%full_panels) then
      wake%pan_wake_len = wake%pan_wake_len + 1
      allocate(pan_p_temp(wake%n_pan_stripes*wake%pan_wake_len))

      do ip = 1,size(wake%pan_p)
        pan_p_temp(ip) = wake%pan_p(ip)
      enddo
      do iw = 1,wake%n_pan_stripes
        pan_p_temp(wake%n_pan_stripes*(wake%pan_wake_len-1)+iw)%p &
                                        => wake%wake_panels(iw,wake%pan_wake_len)
      enddo
      if(allocated(wake%pan_p)) deallocate(wake%pan_p)
      allocate(wake%pan_p(size(pan_p_temp)))
      do ip = 1,size(wake%pan_p)
        wake%pan_p(ip) = pan_p_temp(ip)
      enddo
      deallocate(pan_p_temp)

  endif


  !==> Rings:  Increase the length of the wake, if it is necessary
  if (increase_wake) then
    allocate(pan_p_temp(wake%ndisks*wake%rin_wake_len))
    size_old = 0; if(allocated(wake%rin_p)) size_old = size(wake%rin_p)

    do ip = 1,size_old
      pan_p_temp(ip) = wake%rin_p(ip)
    enddo
    do id = 1,wake%ndisks
      pan_p_temp(size_old+id)%p => wake%wake_rings(id,wake%rin_wake_len)
    enddo

    !manual move alloc
    if(allocated(wake%rin_p)) deallocate(wake%rin_p)
    allocate(wake%rin_p(size(pan_p_temp)))
    do ip = 1,size(wake%rin_p)
      wake%rin_p(ip) = pan_p_temp(ip)
    enddo
    deallocate(pan_p_temp)
  endif

  !==> Panels:  Update the intensities of the panels
  !       From the back, all the vortex intensities come from
  !       the previous panel
  do ipan = wake%pan_wake_len,2,-1
    do iw = 1,wake%n_pan_stripes
      wake%wake_panels(iw,ipan)%mag = wake%wake_panels(iw,ipan-1)%mag
    enddo
  enddo

  !==> End vortices: If the wake is full, attach the end vortex
  !if (wake%pan_wake_len .eq. wake%nmax_pan) then
  if (wake%full_panels) then
    do iw = 1,wake%n_pan_stripes
      wake%end_vorts(iw)%mag => wake%wake_panels(iw,wake%pan_wake_len)%mag
    enddo
  endif

  !==> Rings: Update the intensities of the rings
  !       From the back, all the vortex intensities come from the 
  !       previous panel
  do ir = wake%rin_wake_len,2,-1
    do id = 1,wake%ndisks
      wake%wake_rings(id,ir)%mag = wake%wake_rings(id,ir-1)%mag
    enddo
  enddo
  do id = 1,wake%ndisks
    wake%wake_rings(id,1)%mag  = wake%rin_gen_elems(id)%p%mag
  enddo

  !==> Rings: Update rings geometrical quantities
  do ir = 1,wake%rin_wake_len
    do id = 1,wake%ndisks
      call wake%wake_rings(id,ir)%calc_geo_data( &
                    wake%wake_rings(id,ir)%ver)
    enddo
  enddo

  ! The geometrical quantities of the panels will be all updated in prepare
  ! wake after the geometry have been updated


end subroutine update_wake

!----------------------------------------------------------------------

subroutine compute_vel_from_all(elems, wake, pos, sim_param, vel)
 type(t_pot_elem_p), intent(in) :: elems(:)
 type(t_wake), intent(in) :: wake
 real(wp), intent(in) :: pos(3)
 type(t_sim_param), intent(in) :: sim_param
 real(wp), intent(out) :: vel(3)

 integer :: ie
 real(wp) :: v(3)
  
  vel = 0.0_wp

  !calculate the influence of the solid bodies
  do ie=1,size(elems)
    call elems(ie)%p%compute_vel(pos, sim_param%u_inf, v)
    vel = vel + v/(4*pi)
  enddo

  ! calculate the influence of the wake panels
  do ie=1,size(wake%pan_p)
    call wake%pan_p(ie)%p%compute_vel(pos, sim_param%u_inf, v)
    vel = vel + v/(4*pi)
  enddo

  ! calculate the influence of the wake rings
  do ie=1,size(wake%rin_p)
    call wake%rin_p(ie)%p%compute_vel(pos, sim_param%u_inf, v)
    vel = vel+ v/(4*pi)
  enddo

  !calculate the influence of the end vortex
  !TODO: check what happens when it is not active
  do ie=1,size(wake%end_vorts)
    call wake%end_vorts(ie)%compute_vel(pos, sim_param%u_inf, v)
    vel = vel+ v/(4*pi)
  enddo

  !calculate the influence of particles
  do ie=1,size(wake%part_p)
    call wake%part_p(ie)%p%compute_vel(pos, sim_param%u_inf, v)
    vel = vel+ v/(4*pi)
  enddo



end subroutine compute_vel_from_all

!----------------------------------------------------------------------

subroutine get_vel_free(this, elems, wake, pos, sim_param, vel)
 class(t_free_wake) :: this
 type(t_pot_elem_p), intent(in) :: elems(:)
 type(t_wake), intent(in) :: wake
 real(wp), intent(in) :: pos(3)
 type(t_sim_param), intent(in) :: sim_param
 real(wp), intent(out) :: vel(3)

  call compute_vel_from_all(elems, wake, pos, sim_param, vel)
  
  vel = vel + sim_param%u_inf

end subroutine get_vel_free

!----------------------------------------------------------------------

subroutine get_vel_rigid(this, elems, wake, pos, sim_param, vel)
 class(t_rigid_wake) :: this
 type(t_pot_elem_p), intent(in) :: elems(:)
 type(t_wake), intent(in) :: wake
 real(wp), intent(in) :: pos(3)
 type(t_sim_param), intent(in) :: sim_param
 real(wp), intent(out) :: vel(3)

! vel = sim_param%u_inf
 vel = sim_param%rigid_wake_vel

end subroutine get_vel_rigid
!----------------------------------------------------------------------

subroutine avoid_collision(elems, wake, part, sim_param, vel)
 type(t_pot_elem_p), intent(in) :: elems(:)
 type(t_wake), intent(inout) :: wake
 type(t_vortpart), intent(inout) :: part
 type(t_sim_param), intent(in) :: sim_param
 real(wp), intent(inout) :: vel(3)

 integer :: ie
 real(wp) :: dist(3), n(3)
 real(wp) :: pos(3)
 real(wp) :: distn, distnor, normvel
 real(wp) :: damp_radius, cont, rad_mult, k

 !damp_radius = 0.3
 cont = 0.95_wp
 !cont = 1.0
 rad_mult = 1.3_wp
 pos = part%cen
 k = 0.85_wp
 

  !calculate the influence of the solid bodies
  do ie=1,size(elems)
    dist = pos-elems(ie)%p%cen
    damp_radius = sim_param%dt*sim_param%u_ref*rad_mult + &
                                          maxval(elems(ie)%p%edge_len)/2.0_wp
    distn = norm2(dist)
    !if it is in the check radius perform calculations
    if ((distn .lt. damp_radius)) then
      n = elems(ie)%p%nor
      normvel = sum(vel*n) 
      distnor = sum(dist * elems(ie)%p%nor)
      if(normvel .lt. -k*distnor/sim_param%dt) then
        !normvel = -distnor/sim_param%dt*(1 - &
        !                         1/(-normvel*sim_param%dt/8.0_wp/distnor+1)**8)
        normvel = -distnor/sim_param%dt*(1.0-k)*(1.0 - &
                 1.0/(-(normvel+k*distnor/sim_param%dt)* &
                 sim_param%dt/8.0_wp/distnor/(1-k) + 1.0)**8)
        !normvel = max(normvel,-distnor/sim_param%dt)
      endif
      if(normvel .lt. -cont*distnor/sim_param%dt) then
        part%free = .true.
        wake%n_prt = wake%n_prt -1
        return
      endif
      vel = vel - (sum(vel*n) - normvel) * n
      !vel = vel - (sum(vel*n) + cont*distnor/sim_param%dt) * n
      !vel = vel - sum(vel*n) * (1-damp) * n
      !distnor = sum(dist * elems(ie)%p%nor)
      !if (abs(distnor) .lt. 0.03_wp*minval(elems(ie)%p%edge_len)) then
      !  part%free = .true.
      !  wake%n_prt = wake%n_prt -1
      !  return
      !endif
    endif
  enddo



end subroutine avoid_collision

end module mod_wake
