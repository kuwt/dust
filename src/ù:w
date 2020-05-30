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

!> Module for introducing hinges and rotating parts in a component
module mod_hinges

use mod_param, only: &
  wp, max_char_len, pi

use mod_math, only: &
  cross

use mod_parse, only: &
  t_parse, getstr, getint, getreal, getrealarray, getlogical, getsuboption, &
  countoption, finalizeparameters

implicit none

public :: t_hinge, t_hinge_input, build_hinges, hinge_input_parser, &
          initialize_hinge_config

private

! ---------------------------------------------------------------
!> Hinge input data type
type :: t_hinge_input
  character(len=max_char_len) :: tag
  character(len=max_char_len) :: nodes_input
  real(wp) :: node1(3), node2(3)
  integer  :: n_nodes
  real(wp), allocatable :: rr(:,:)
  character(len=max_char_len) :: node_file
  real(wp) :: ref_dir(3)
  real(wp) :: offset
  character(len=max_char_len) :: rotation_input
  real(wp) :: rotation_amplitude
  real(wp) :: rotation_omega
  real(wp) :: rotation_phase
end type t_hinge_input

! ---------------------------------------------------------------
!> Hinge connectivity, meant for rigid rotation and bleding regions
type :: t_hinge_conn
  !> Local index of the nodes performing the desired motion
  integer, allocatable :: node_id(:)
  !> Surface node performing motion vs. hinge node connectivity
  integer , allocatable :: ind(:,:)
  !> Surface node performing motion vs. hinge node connectivity, weights
  ! for the weighted average of the motion
  real(wp), allocatable :: wei(:,:)
  !> Node to hinge connectivity array of objs
  type(t_n2h_conn), allocatable :: n2h(:)
end type t_hinge_conn

! ---------------------------------------------------------------
!> Hinge node to surface node connectivity. The different number of surface
! nodes per hinge nodes requires an array of obj, containing indices and
! weights
type :: t_n2h_conn
  !> Indices of the surface nodes
  integer, allocatable :: p2h(:)
  !> Weights
  real(wp),allocatable :: w2h(:)
end type t_n2h_conn

! ---------------------------------------------------------------
!> Hinge node configurations: to be used in defining reference
! and actual configurations
type :: t_hinge_config

  !> Position of the first and last point of the hinge in the local
  ! reference frame of the component
  real(wp) :: rr0(3), rr1(3)
  !> Local coordinates of the hinge nodes
  real(wp), allocatable :: rr(:,:)

  !> Unit vectors of the hinge node reference frames
  ! User input in local reference frame, for reference configuration: h, v
  !    h,
  !    n = cross(v,h) , n = n/norm2(n)
  !    v = cross(h,n)
  !> h: rotation axis
  real(wp), allocatable :: h(:,:)
  !> v: zero direction (user inputs are overwritten, in order to 
  ! build ortonormal reference frameshinge%rot  
  real(wp), allocatable :: v(:,:)
  !> n: normal direction
  real(wp), allocatable :: n(:,:)

end type t_hinge_config

!> Hinge type
type :: t_hinge

  !> Number of points for transferring the motion from the hinge to the
  ! surface, through a weighted averageamplitude
  ! *** to do *** hardcoded, so far. Read as an input, with a default value,
  ! equal to 2 (or 1?)
  integer :: n_wei = 2

  !> Order of the norm used for computing distance-based weights
  ! *** to do *** hardcoded, so far. Read as an input, with a default value,
  ! equal to 1 (or 2?)
  integer :: w_order = 1

  !> Type: parametric, from_file
  character(len=max_char_len) :: nodes_input

  !> N.nodes of the hinge
  integer :: n_nodes

  !> Offset for avoiding irregular behavior
  real(wp) :: offset

  !> Offset for avoiding irregular behavior
  real(wp) :: ref_dir(3)

  !> Type: constant, function, coupling
  character(len=max_char_len) :: input_type

  !> function:const, :sin, :cos amplitude, angular velocity and phase
  real(wp) :: f_ampl
  real(wp) :: f_omega
  real(wp) :: f_phase

  !> Array of the rotation angle (read as an input, or from coupling nodes)
  real(wp), allocatable :: theta(:)
  real(wp), allocatable :: theta_old(:)

  !> Reference configuration
  type(t_hinge_config) :: ref
  !> Actual configuration
  type(t_hinge_config) :: act

  !> Connectivity for the rigid rotation motion of a ``control surface''
  type(t_hinge_conn) :: rot
  !> Connectivity of the blending region to avoid irregular behavior 
  ! with a ``control surface'' large rotation (blending region extends
  ! from -offset to offset qualitatively in the ref_dir of the hinge)
  type(t_hinge_conn) :: blen

  contains
 
  procedure, pass(this) :: build_connectivity
  procedure, pass(this) :: init_theta
  procedure, pass(this) :: from_reference_to_actual_config ! empty (useless?)
  procedure, pass(this) :: update_hinge_nodes
  procedure, pass(this) :: update_theta
  procedure, pass(this) :: hinge_deflection

end type t_hinge

! ---------------------------------------------------------------

contains

! ---------------------------------------------------------------
!> Build hinge connectivity from component nodes rr, expressed in the
! local reference frame, and the hinge nodes.
subroutine build_connectivity(this, loc_points)
  class(t_hinge), intent(inout) :: this
  real(wp),       intent(in)    :: loc_points(:,:)

  real(wp) :: hinge_width
  integer  :: nb, nh, ib, ih, iw
  real(wp), allocatable :: rrb(:,:), rrh(:,:)
  real(wp) :: Rot(3,3) = 0.0_wp

  real(wp), allocatable :: dist_all(:), wei_v(:)
  integer , allocatable ::              ind_v(:)

  integer , allocatable :: rot_node_id(:), ble_node_id(:)
  integer , allocatable :: rot_ind(:,:)  , ble_ind(:,:)
  real(wp), allocatable :: rot_wei(:,:)  , ble_wei(:,:)

  integer , allocatable :: rot_p2h(:,:)   ! *** to do *** improve the actual inefficient
  real(wp), allocatable :: rot_w2h(:,:)   ! and memory intensive implementation
  integer , allocatable :: rot_i2h(:)     ! Look for ***** in this routine
  integer , allocatable :: ble_p2h(:,:)   ! *** to do *** improve the actual inefficient
  real(wp), allocatable :: ble_w2h(:,:)   ! and memory intensive implementation
  integer , allocatable :: ble_i2h(:)     ! Look for ***** in this routine

  integer :: nrot, nble

  !> N. of surfae and hinge nodes
  nb = size(loc_points,2)
  nh = this%n_nodes
 
  !> Coordinates in the hinge reference frame
  ! Rotation matrix, build with the local ortonormal ref.frame of 
  ! the first hinge node
  Rot(1,:) = this % ref % v(:,1)
  Rot(2,:) = this % ref % h(:,1)
  Rot(3,:) = this % ref % n(:,1)

  allocate( rrb(3,nb) );  allocate( rrh(3,nh) )
  do ib = 1, nb
    rrb(:,ib) = matmul( Rot, loc_points(:,ib) - this%ref%rr(:,1) )
  end do
  do ib = 1, nh
    rrh(:,ib) = matmul( Rot, this%ref%rr(:,ib) - this%ref%rr(:,1) )
  end do

  ! hinge width, measured in the hinge direction
  hinge_width = rrh(2,nh) - rrh(2,1) 

  !> Compute connectivity and weights
  ! Allocate auxiliary node_id(:), ind(:,:), wei(:,:) arrays
  allocate(rot_node_id(        nb)); rot_node_id = 0
  allocate(ble_node_id(        nb)); ble_node_id = 0
  allocate(rot_ind(this%n_wei, nb));     rot_ind = 0
  allocate(ble_ind(this%n_wei, nb));     ble_ind = 0
  allocate(rot_wei(this%n_wei, nb));     rot_wei = 0.0_wp
  allocate(ble_wei(this%n_wei, nb));     ble_wei = 0.0_wp

  ! diff_all and dist_all auxiliaary arrays
  allocate(dist_all(   nb)); dist_all = 0.0_wp

  ! auxiliary matrix for hinge to surf connectivity 
  ! ***** to do ***** improve the implementation
  allocate(rot_p2h(nh,nb)); rot_p2h = 0      ;  allocate(ble_p2h(nh,nb)); ble_p2h = 0
  allocate(rot_w2h(nh,nb)); rot_w2h = 0.0_wp ;  allocate(ble_w2h(nh,nb)); ble_w2h = 0.0_wp
  allocate(rot_i2h(nh   )); rot_i2h = 0      ;  allocate(ble_i2h(nh   )); ble_i2h = 0

  nrot = 0; nble = 0
  ! Loop over all the surface points
  do ib = 1, nb

    if ( ( rrb(2,ib) .gt. 0.0_wp ) .and. ( rrb(2,ib) .lt. hinge_width ) ) then

      if ( rrb(1,ib) .lt. -this%offset ) then
        ! do nothing
      elseif ( rrb(1,ib) .gt. this%offset ) then ! rigid rotation

        nrot = nrot + 1
        rot_node_id(nrot) = ib

        do ih = 1, nh
          dist_all(ih) = abs( rrb(2,ib) - rrh(2,ih) )
        end do

        call sort_vector_real( dist_all, this%n_wei, wei_v, ind_v )
        
        wei_v = 1.0_wp / wei_v**this%w_order
        wei_v = wei_v / sum(wei_v)

        rot_wei(:,nrot) = wei_v
        rot_ind(:,nrot) = ind_v

        ! *****
        do iw = 1, this%n_wei
          rot_i2h(ind_v(iw)) = rot_i2h(ind_v(iw)) + 1
          rot_p2h(ind_v(iw),   rot_i2h(ind_v(iw))) = ib
          rot_w2h(ind_v(iw),   rot_i2h(ind_v(iw))) = wei_v(iw)
        end do
        ! *****

      else ! blending region
        
        nble = nble + 1
        ble_node_id(nble) = ib

        do ih = 1, nh
          dist_all(ih) = abs( rrb(2,ib) - rrh(2,ih) )
        end do

        call sort_vector_real( dist_all, this%n_wei, wei_v, ind_v )
        
        wei_v = 1.0_wp / wei_v**this%w_order
        wei_v = wei_v / sum(wei_v)

        ble_wei(:,nble) = wei_v
        ble_ind(:,nble) = ind_v

        ! *****
        do iw = 1, this%n_wei
          ble_i2h(ind_v(iw)) = ble_i2h(ind_v(iw)) + 1
          ble_p2h(ind_v(iw),   ble_i2h(ind_v(iw))) = ib
          ble_w2h(ind_v(iw),   ble_i2h(ind_v(iw))) = wei_v(iw)
        end do
        ! *****

      end if

    end if

  end do

  !> Fill hinge object, with the connectivity and weight arrays
  allocate(this%rot %node_id(        nrot)); this%rot %node_id = rot_node_id(1:nrot)
  allocate(this%rot %ind(this%n_wei, nrot)); this%rot %ind     = rot_ind( :, 1:nrot)
  allocate(this%rot %wei(this%n_wei, nrot)); this%rot %wei     = rot_wei( :, 1:nrot)
  allocate(this%blen%node_id(        nble)); this%blen%node_id = ble_node_id(1:nble)
  allocate(this%blen%ind(this%n_wei, nble)); this%blen%ind     = ble_ind( :, 1:nble)
  allocate(this%blen%wei(this%n_wei, nble)); this%blen%wei     = ble_wei( :, 1:nble)
  allocate(this%rot %n2h(nh))
  allocate(this%blen%n2h(nh))
  do ih = 1, nh
    allocate(this%rot %n2h(ih)%p2h( rot_i2h(ih) )) ; this%rot %n2h(ih)%p2h = rot_p2h( ih, 1:rot_i2h(ih) )
    allocate(this%rot %n2h(ih)%w2h( rot_i2h(ih) )) ; this%rot %n2h(ih)%w2h = rot_w2h( ih, 1:rot_i2h(ih) )
    allocate(this%blen%n2h(ih)%p2h( ble_i2h(ih) )) ; this%blen%n2h(ih)%p2h = ble_p2h( ih, 1:ble_i2h(ih) )
    allocate(this%blen%n2h(ih)%w2h( ble_i2h(ih) )) ; this%blen%n2h(ih)%w2h = ble_w2h( ih, 1:ble_i2h(ih) )
  end do


  !> Explicit deallocations
  deallocate(rrb, rrh, dist_all, wei_v, ind_v)
  deallocate(rot_node_id, rot_ind, rot_wei)
  deallocate(ble_node_id, ble_ind, ble_wei)
  deallocate(rot_i2h, rot_p2h, rot_w2h) ! *****
  deallocate(ble_i2h, ble_p2h, ble_w2h) ! *****


end subroutine build_connectivity

! ---------------------------------------------------------------
!> Naif sort
subroutine sort_vector_real( vec, nel, sor, ind )
  real(wp), intent(inout) :: vec(:)
  integer , intent(in) :: nel
  real(wp), allocatable, intent(out):: sor(:)
  integer , allocatable, intent(out):: ind(:)

  real(wp):: minv
  integer :: i

  allocate(sor(nel)); sor = 0.0_wp
  allocate(ind(nel)); ind = 0

  minv = minval( vec )
  do i = 1, nel
    sor(i) = maxval( vec, 1 )
    ind(i) = maxloc( vec, 1 )
    vec( ind(i) ) = minv - 0.1_wp ! naif
  end do


end subroutine sort_vector_real

! ---------------------------------------------------------------
!> Compute actual configuration of the hinge nodes
! *** to do *** unpredictable (or better, wrong) behavior when
! restarting a simulation
subroutine init_theta(this, t)
  class(t_hinge), intent(inout) :: this
  real(wp)      , intent(in)    :: t

  if ( t .ne. 0.0_wp ) then
    write(*,*) ' Error in t_hinge % init_theta: t .ne. 0.0_wp. &
               &This argument is meant for future restart capabilities. &
               &So far, must be passed to the function equal to 0.0_wp. Stop'
    stop
  end if

  call this % update_theta( t )

  this%theta_old = 0.0_wp


end subroutine init_theta

! ---------------------------------------------------------------
!> Compute actual configuration of the hinge nodes
subroutine from_reference_to_actual_config(this)
  class(t_hinge), intent(inout) :: this

  ! *** to do ***
  ! still usefull? anything else to do, that is missing in
  ! update_hinge_nodes() subroutine?
  

end subroutine from_reference_to_actual_config

! ---------------------------------------------------------------
!> Update hinge nodes, for non-coupled components
subroutine update_hinge_nodes( this, R, of )
  class(t_hinge), intent(inout) :: this
  real(wp)      , intent(in)    ::  R(:,:)
  real(wp)      , intent(in)    :: of(:)

  !> Actual configuration: node position
  this % act % rr = matmul( R, this % ref % rr )
  this % act % rr(1,:) = this % act % rr(1,:) + of(1) 
  this % act % rr(2,:) = this % act % rr(2,:) + of(2)
  this % act % rr(3,:) = this % act % rr(3,:) + of(3)

  !> Actual configuration: node orientation
  this % act % h = matmul( R, this % ref % h )
  this % act % v = matmul( R, this % ref % v )
  this % act % n = matmul( R, this % ref % n )


end subroutine update_hinge_nodes

! ---------------------------------------------------------------
!> Update hinge nodes, for non-coupled components
subroutine update_theta( this, t )
  class(t_hinge), intent(inout) :: this
  real(wp)      , intent(in)    :: t

  if ( trim(this%input_type) .eq. 'function:const' ) then
    this%theta = this%f_ampl
  elseif ( trim(this%input_type) .eq. 'function:sin' ) then
    this%theta = this%f_ampl * sin( this%f_omega * t - this%f_phase )
  elseif ( trim(this%input_type) .eq. 'function:cos' ) then
    this%theta = this%f_ampl * cos( this%f_omega * t - this%f_phase )
  else
    write(*,*) ' Error in t_hinge % init_theta(): only working &
               &with function:const, :sin, :cos, so far. Stop'; stop
  end if


end subroutine update_theta

! ---------------------------------------------------------------
!> Update the coordinates rr of the points of the surface, after
! hinge deflection
subroutine hinge_deflection( this, rr, t, postpro )
  class(t_hinge), intent(inout) :: this
  real(wp)      , intent(inout) :: rr(:,:)
  real(wp)      , intent(in)    :: t
  logical, optional, intent(in) :: postpro
  logical            ::   local_postpro
  logical, parameter :: default_postpro = .false.

  real(wp) :: th, th1, thp, yc, xq, yq, xqp, yqp
  real(wp) :: nx(3,3), Rot_I(3,3), Rot(3,3)
  integer :: nrot, nble, ib, ih, ii

  !> Use the same routine for solver (incremental rotation) and
  ! postpro (non-incremental rotation)
  if ( present(postpro) ) then;  local_postpro = postpro
  else                        ;  local_postpro = default_postpro
  end if

  !> n.nodes in the rigid-rotation and in the blending regions
  nrot = size(this%rot %node_id)
  nble = size(this%blen%node_id)

  do ih = 1, this%n_nodes

    if ( .not. local_postpro ) then
      th = ( this % theta(ih) - this % theta_old(ih) ) * pi/180.0_wp
    else 
      th =   this % theta(ih)                          * pi/180.0_wp
    end if

    if ( th .ne. 0.0_wp ) then ! (equality check on real?)

      !Rotation matrix
      nx(1,:) = (/            0.0_wp, -this%act%h(3,ih),  this%act%h(2,ih) /)
      nx(2,:) = (/  this%act%h(3,ih),            0.0_wp, -this%act%h(1,ih) /)
      nx(3,:) = (/ -this%act%h(2,ih),  this%act%h(1,ih),            0.0_wp /)

      Rot_I = sin(th) * nx + ( 1.0_wp - cos(th) ) * matmul( nx, nx )
      Rot = Rot_I
      Rot(1,1) = Rot_I(1,1) + 1.0_wp
      Rot(2,2) = Rot_I(2,2) + 1.0_wp
      Rot(3,3) = Rot_I(3,3) + 1.0_wp

      !> Rigid rotation
      do ib = 1, size(this%rot%n2h(ih)%p2h)
        ii = this%rot%n2h(ih)%p2h(ib)
        rr(:,ii) = rr(:,ii) + &
                   this%rot%n2h(ih)%w2h(ib) * &
                   matmul( Rot_I, rr(:,ii)-this%act%rr(:,ih) )
        ! rr(:,ii) = this%rot%n2h(ih)%w2h(ib) * &
        !          ( this%act%rr(:,ih) + &
        !            matmul( Rot, rr(:,ii)-this%act%rr(:,ih) ) )
      end do
    
      !> Blending region
      do ib = 1, size(this%blen%n2h(ih)%p2h)

        ii = this%blen%n2h(ih)%p2h(ib)
        th = -th  ! dirty implementation *** to do ***

        !> coordinate of the centre of the circle used for blending, 
        ! in the n-direction
        yc = cos(th)/sin(th) * this%offset * ( 1.0_wp + cos(th) ) + &
                               this%offset * sin(th)
        !> Some auxiliary quantities
        xq = sum( ( rr(:,ii) - this%act%rr(:,ih) ) * this%act%v(:,ih) )
        yq = sum( ( rr(:,ii) - this%act%rr(:,ih) ) * this%act%n(:,ih) )
        thp = 0.5_wp * ( xq + this%offset ) / this%offset * th
        xqp = yc*sin(thp)          - this%offset - yq*sin(thp) - xq
        yqp = yc*(1.0_wp-cos(thp))               + yq*cos(thp) - yq

        !> Update coordinates
        rr(:,ii) = rr(:,ii) + &
                   this%blen%n2h(ih)%w2h(ib) * &
                 ( xqp * this%act%v(:,ih) + yqp * this%act%n(:,ih) )

        th = -th  ! dirty implementation *** to do ***

      end do

    end if

    !> Update theta_old
    this % theta_old(ih) = this % theta(ih)

  end do


end subroutine hinge_deflection

! ---------------------------------------------------------------
!> Build hinge configuration
! Used to build reference configuration in load_component(), and
! initialize actual configuration, as the reference configuration
subroutine initialize_hinge_config( h_config, hinge )
  type(t_hinge_config), intent(inout) :: h_config
  type(t_hinge)       , intent(inout) :: hinge

  real(wp) :: hv(3), vv(3), nv(3)
  integer :: i

  !> rr0, rr1 (useless?)
  h_config%rr0 = hinge%ref%rr(:,1)
  h_config%rr1 = hinge%ref%rr(:,hinge%n_nodes)
  hv = h_config%rr1 - h_config%rr0;  hv = hv / norm2(hv)
  vv = hinge%ref_dir / norm2(hinge%ref_dir)
  nv = cross( vv, hv ); nv = nv / norm2(nv)
  vv = cross( hv, nv ); vv = vv / norm2(vv)  ! useless normalization

  !> h, v, n
  allocate( h_config%h(3, hinge%n_nodes) )
  allocate( h_config%v(3, hinge%n_nodes) )
  allocate( h_config%n(3, hinge%n_nodes) )

  do i = 1, hinge%n_nodes

    h_config%h(:,i) = hv;  h_config%v(:,i) = vv;  h_config%n(:,i) = nv 

  end do

end subroutine initialize_hinge_config


! ---------------------------------------------------------------
!> build hinges_input, used in dust_pre for generating geometry input file
! for the solver
subroutine build_hinges( geo_prs, n_hinges, hinges )
  type(t_parse), intent(inout) :: geo_prs
  integer      , intent(in)    :: n_hinges
  type(t_hinge_input), allocatable, intent(inout) :: hinges(:)

  type(t_parse), pointer :: hinge_prs
  integer :: i, j

  if ( allocated(hinges) ) deallocate(hinges)
  allocate( hinges(n_hinges) )

  do i = 1, n_hinges

    !> De-associate, before reading next hinge
    hinge_prs => null()
    !> Open hinge sub-parser
    call getsuboption(geo_prs, 'Hinge', hinge_prs)

     hinges(i) % tag = getstr(hinge_prs, 'Hinge_Tag')
     hinges(i) % nodes_input = getstr(hinge_prs, 'Hinge_Nodes_Input')

     if ( trim(hinges(i) % nodes_input) .eq. 'parametric' ) then
       hinges(i) % node1 = getrealarray(hinge_prs, 'Node1', 3)
       hinges(i) % node2 = getrealarray(hinge_prs, 'Node2', 3)
       hinges(i) % n_nodes = getint(hinge_prs, 'N_Nodes')
       hinges(i) % node_file = 'hinge with parametric input. If you read this &
           &string, something probabily went wrong.'
       allocate( hinges(i)%rr( 3, hinges(i)%n_nodes ) )
       do j = 1, hinges(i)%n_nodes
         hinges(i) % rr(:,j) = hinges(i) % node1 + &
                             ( hinges(i) % node2 - hinges(i) % node1 ) * &
                             dble(j-1)/dble(hinges(i)%n_nodes-1)
       end do

     elseif ( trim(hinges(i) % nodes_input) .eq. 'from_file' ) then
       ! *** to do *** fill dummy node1, node2, n_nodes fields
       hinges(i) % node_file = getstr(hinge_prs,'Node_File')
       call read_hinge_nodes( hinges(i)%node_file, &
                              hinges(i)%n_nodes  , &
                              hinges(i)%rr )
     else
       ! *** to do *** use error message handling
       write(*,*) ' Wrong Hinge_Nodes_Input input. Stop '; stop
     end if

     !> Hinge reference direction (zero rotation direction) and offset to
     ! avoid irregular behavior during hinge rotation
     hinges(i) % ref_dir = getrealarray(hinge_prs, 'Hinge_Ref_Dir', 3)
     hinges(i) % offset  = getreal(hinge_prs, 'Hinge_Offset')

     !> Hinge input: function, amplitude
     !> Overwrite 'constant' rotation_input with 'function:const'
     if ( trim( hinges(i)%rotation_input ) .eq. 'constant' ) then
       hinges(i)%rotation_input = 'function:const'
     end if
     hinges(i) % rotation_input     = getstr(hinge_prs, 'Hinge_Rotation_Input')
     if ( ( trim(hinges(i)%rotation_input) .ne. 'function:const' ) .or. &
          ( trim(hinges(i)%rotation_input) .ne. 'function:sin'   ) .or. &
          ( trim(hinges(i)%rotation_input) .ne. 'function:cos'   ) .or. &
          ( trim(hinges(i)%rotation_input) .ne. 'from_file'      ) .or. &
          ( trim(hinges(i)%rotation_input) .ne. 'coupling'       ) ) then
     else
       if ( trim(hinges(i)%rotation_input) .eq. 'from_file' ) then
         write(*,*) ' Error in t_hinge%build_hinge(): rotation_input = from_file &
                    &not implemented yet.'
       end if
       if ( trim(hinges(i)%rotation_input) .eq. 'coupling' ) then
         write(*,*) ' Error in t_hinge%build_hinge(): rotation_input = coupling &
                    &not implemented yet.'
       end if
     end if
     hinges(i) % rotation_amplitude = getreal(hinge_prs,'Hinge_Rotation_Amplitude')
     hinges(i) % rotation_omega     = getreal(hinge_prs,'Hinge_Rotation_Omega')
     hinges(i) % rotation_phase     = getreal(hinge_prs,'Hinge_Rotation_Phase')
     

!    ! check ---
!    write(*,*) ' Hinge id:', i
!    write(*,*) ' _Tag        : ', trim(hinges(i)%tag)
!    write(*,*) ' _Nodes_Input: ', trim(hinges(i)%nodes_input)
!    ! check ---

  end do


end subroutine build_hinges

! ---------------------------------------------------------------
!> Read ascii file of hinge node coordinates, w/ input type from_file
!  Three-column ascii file is expected
subroutine read_hinge_nodes( filen, n_nodes, rr )
  character(len=max_char_len), intent(in)  :: filen
  integer                    , intent(out) :: n_nodes
  real(wp), allocatable      , intent(out) :: rr(:,:)

  integer :: i, io
  integer :: fid = 21

  if ( allocated(rr) ) deallocate(rr)

  n_nodes = 0;  io = 0

  !> Preliminary read for determining the n. of lines
  open(unit=fid, file=trim(filen))
  do while ( io .eq. 0 )
    read(fid,*,iostat=io) ! dummy
    n_nodes = n_nodes + 1
  end do
  close(fid) 

  !> Allocate and fill rr array
  allocate(rr(3,n_nodes))
  open(unit=fid, file=trim(filen))
  do i = 1, n_nodes
    read(fid,*) rr(:,i)
  end do
  close(fid) 


end subroutine read_hinge_nodes

! ---------------------------------------------------------------
!> Hinge input parser, called in mod_build_geo.f90 by dust_pre preprocessor
subroutine hinge_input_parser( geo_prs, hinge_prs )
  type(t_parse),          intent(inout) :: geo_prs
  type(t_parse), pointer, intent(inout) :: hinge_prs

  call geo_prs%CreateIntOption('n_hinges', &
              'N. of hinges and rotating parts (e.g. aileron) of the component', &
              '0') ! default: no hinges -> n_hinges = 0
  
  call geo_prs%CreateSubOption('Hinge', 'Parser for hinge input', &
               hinge_prs, multiple=.true. )

  call hinge_prs%CreateStringOption('Hinge_Tag','Name of the hinge')
  call hinge_prs%CreateStringOption('Hinge_Nodes_Input', &
         'Type of hinge nodes input: parametric or from_file.')
  call hinge_prs%CreateIntOption('N_Nodes','N.hinge nodes')
  call hinge_prs%CreateRealArrayOption('Node1', &
      'First node of the hinge. Components in the local ref.frame of the component')
  call hinge_prs%CreateRealArrayOption('Node2', &
      'Second node of the hinge. Components in the local ref.frame of the component')
  call hinge_prs%CreateRealArrayOption('Hinge_Ref_Dir', &
      'Reference direction of the hinges, indicating zero-deflection direction in &
      &the local ref.frame of the component')
  call hinge_prs%CreateRealOption('Hinge_Offset','Offset in the Ref_Dir needed for &
      &avoiding irregular behavior of the surface for large deflections')
  call hinge_prs%CreateStringOption('Hinge_Rotation_Input', &
      'Input type of the rotation: constant, function, from_file, coupling')
  call hinge_prs%CreateRealOption('Hinge_Rotation_Amplitude', &
      'Amplitude of the rotation, for constant, function:const, :sin, &
      &:cos Rotation_Input')
  call hinge_prs%CreateRealOption('Hinge_Rotation_Omega', &
      'Angular velocity of the rotation, for constant, function:const, &
      &:sin, :cos Rotation_Input', '0.0')
  call hinge_prs%CreateRealOption('Hinge_Rotation_Phase', &
      'Phase of the rotation, for constant, function:const, :sin, :cos &
      &Rotation_Input', '0.0')
  ! *** to do ***
  ! add all the fields required for all the input types


end subroutine hinge_input_parser

! ---------------------------------------------------------------

end module mod_hinges
