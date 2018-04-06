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


!> Module to store and update all the reference frames
!!
!! The reference frames are built in a hierarchical manner. The base reference 
!! frame is the fixed cartesian frame and all the subsequent frames are defined
!! upon that, and all the geometry inside dust is calculated with respect to 
!! that frame. The base frame cannot be user defined and has internal index, 
!! as well as tag equal to 0.
!!
!! All the following reference frames can be defined by the user in an input 
!! file (see \ref build_references for details). They are defined by:
!! - A parent tag referring to the reference frame upon which the frame is
!!   built
!! - An origin of the reference frame, expressed in coordinates in the parent
!!   reference frame
!! - An orientation matrix, defining the orientation of the frame axis in 
!!   components of the parent reference frame
!!
!! Additionaly it is possible to impose a movement of the reference frame with
!! respect to the parent reference frame. 
!! At the moment the only movement that can be defined is a constant rate 
!! rotation around a pole and an axis. 

module mod_reference


use mod_param, only: &
  wp, max_char_len, nl

use mod_parse, only: &
  t_parse, getstr, getint, getreal, getrealarray, getlogical, countoption, &
  finalizeparameters, t_link, check_opt_consistency

use mod_handling, only: &
  error, warning, info, printout, dust_time, t_realtime

!----------------------------------------------------------------------

implicit none

public :: t_ref, build_references, update_all_references, destroy_references

private

!----------------------------------------------------------------------

!> Reference Frame type
!!
!! Employed to define a reference frame relative to another one
type :: t_ref

  !> Reference id
  integer :: id

  !> Reference tag (can be non-consecutive)
  integer :: tag

  !> Reference Name
  character(len=max_char_len) :: name

  !> Parent reference id, used to index the references array
  integer :: parent_id

  !> Parent tag, can be non consecutive
  integer :: parent_tag

  !> Number of children references
  integer :: n_chil

  !> Id of children references (is possible to use pointers here)
  integer, allocatable :: chil_id(:)

  !> Origin in the parent
  real(wp), allocatable :: orig(:)

  !> Frame with respect to the parent
  real(wp), allocatable :: frame(:,:)

  !> Is the frame moving with respect to the parent?
  logical :: self_moving

  !> Is the reference frame moving ? (might be stationary with respect to a 
  !! moving parent)
  logical :: moving


  !> Moviment type
  character(len=max_char_len) :: mov_type

  !> Rotation pole
  real(wp), allocatable :: pole(:)

  !> Rotation axis
  real(wp), allocatable :: axis(:)

  !> Rotation rate around the axis
  real(wp) :: Omega

  !> Starting rotation angle
  real(wp) :: psi_0

  !> Total offset with respect to the base reference
  real(wp), allocatable :: of_g(:)

  !> Total rotation with respect to the base reference
  real(wp), allocatable :: R_g(:,:)

  
  !> Total frame velocity with respect to the base reference
  real(wp), allocatable :: f_g(:)

  !> Total frame rotation rate with respect to the base reference
  real(wp), allocatable :: G_g(:,:)


  contains

  procedure, pass(this) :: update_ref  => reference_update_ref

  procedure, pass(this) :: update_self => reference_update_self


end type t_ref

!-----------------------------------

character(len=*), parameter :: this_mod_name = 'mod_reference'

!----------------------------------------------------------------------

contains

!----------------------------------------------------------------------

!> Build all the set of reference frames, reading them from the reference
!! frames file
!!
!! The input file used is either the one present in the main input file, or 
!! if not specified the main input file is used.
!!
!! The base reference frame is the fixed one, and cannot be custom defined, 
!! and has tag nr 0.
!! All the other reference frames must hierarchically refer to the fixed 
!! reference frame. 
!!
!! First all the details of the reference are read from the file. Then
!! the tree of parent-children reference frames is built, and then the 
!! moving properties of eache frame is set traversing the frames tree
!!
!! Each reference frame can be both self_moving, if moving with respect to the
!! parent, or moving if either moving or fixed on a moving parent
subroutine build_references(refs, reference_file)
 type(t_ref), allocatable, intent(out)   :: refs(:)
 character(len=*), intent(in) :: reference_file

 type(t_parse) :: ref_prs
 integer :: n_refs
 integer :: iref, iref2
 integer :: n_child
 integer, allocatable :: temp_chil(:)
 character(len=max_char_len) :: msg
 type(t_link), pointer :: lnk

 character(len=*), parameter :: this_sub_name = 'build_references'

  !Define all the parameters to be read
  call ref_prs%CreateIntOption('Reference_Tag','Integer tag of reference frame',&
               multiple=.true.)
  call ref_prs%CreateIntOption('Parent_Tag','Tag of the parent reference', &
               multiple=.true.)
  call ref_prs%CreateLogicalOption('Moving','Is the reference moving', &
               multiple=.true.)
  call ref_prs%CreateRealArrayOption('Origin','Origin of reference frame with&
               &respect to the parent', multiple=.true.)
  call ref_prs%CreateRealArrayOption('Orientation','Orientation of reference &
               &frame with respect to the parent', multiple=.true.)
  call ref_prs%CreateStringOption('MovementType','Kind of moving imposed', &
               multiple=.true.)
  !For the moment allowing just rotation with a certain angular speed around
  !an axis
  call ref_prs%CreateRealArrayOption('Rot_Pole','Pole of rotation in parent &
               &reference frame', multiple=.true.)
  call ref_prs%CreateRealArrayOption('Rot_Axis','Axis of rotation in parent &
               &reference frame', multiple=.true.)
  call ref_prs%CreateRealOption('Rot_Rate','Rate of rotation around axis', &
                multiple=.true.)
  call ref_prs%CreateRealOption('Psi_0','Starting rotation angle', &
                multiple=.true.)

  
  !read the file
  call ref_prs%read_options(trim(reference_file),printout_val=.true.)

  n_refs = countoption(ref_prs,'Reference_Tag') 
  
  !TODO: here we should check that all the other options have the same number
  !of occurrencies

  ! IMPORTANT: the references are allocated starting from zero, zero is the 
  ! base reference, and then all the other references occupy the following
  ! position. Remember to correctly set the numbering starting from zero for
  ! all the other subroutines that employ the references
  allocate(refs(0:n_refs))


  !Setup the base reference
  refs(0)%id = 0
  refs(0)%tag = 0
  refs(0)%parent_id = -1
  refs(0)%parent_tag = -1
  refs(0)%n_chil = 0
  allocate(refs(0)%orig(3), refs(0)%frame(3,3))
  refs(0)%orig = (/0.0_wp, 0.0_wp, 0.0_wp/)
  refs(0)%frame = reshape((/1.0_wp, 0.0_wp, 0.0_wp, &
                            0.0_wp, 1.0_wp, 0.0_wp, &
                            0.0_wp, 0.0_wp, 1.0_wp/),(/3,3/))
  refs(0)%self_moving = .false.
  refs(0)%moving = .false.
  allocate(refs(0)%of_g(3), refs(0)%R_g(3,3))
  refs(0)%of_g = (/0.0_wp, 0.0_wp, 0.0_wp/)
  refs(0)%R_g = reshape((/1.0_wp, 0.0_wp, 0.0_wp, &
                          0.0_wp, 1.0_wp, 0.0_wp, & 
                          0.0_wp, 0.0_wp, 1.0_wp/),(/3,3/))
  !Setup the other references
  do iref = 1,n_refs
    refs(iref)%id = iref
    refs(iref)%tag = getint(ref_prs,'Reference_Tag')
    refs(iref)%parent_tag = getint(ref_prs,'Parent_Tag')
    refs(iref)%n_chil = 0

    allocate(refs(iref)%orig(3), refs(iref)%frame(3,3))
    refs(iref)%orig  = getrealarray(ref_prs,'Origin',3, olink=lnk)
    call check_opt_consistency(lnk,next=.true.,next_opt='Orientation')
    refs(iref)%frame = reshape(getrealarray(ref_prs,'Orientation',9),(/3,3/))
    !allocated here, will be set in update_all_refs
    allocate(refs(iref)%of_g(3), refs(iref)%R_g(3,3))

    refs(iref)%self_moving = getlogical(ref_prs,'Moving')
    refs(iref)%moving = .false. !standard, will be checked later
    if (refs(iref)%self_moving) then

      refs(iref)%mov_type = trim(getstr(ref_prs,'MovementType'))
      select case (trim(refs(iref)%mov_type))

       case('constant_rotation')
        allocate(refs(iref)%pole(3), refs(iref)%axis(3))
        refs(iref)%pole  = getrealarray(ref_prs,'Rot_Pole',3)
        refs(iref)%axis  = getrealarray(ref_prs,'Rot_Axis',3)
        refs(iref)%Omega = getreal(ref_prs,'Rot_Rate')
        refs(iref)%psi_0 = getreal(ref_prs,'Psi_0')

       case default
         call error(this_sub_name, this_mod_name, 'Unknown type of movement')

      end select
    endif
  enddo

  !Generate the parent-children references
  do iref = 1,n_refs
    refs(iref)%parent_id = -1
    
    !look for the parent in all the other references (master included)
    do iref2 = 0,n_refs
      !if found
      if (refs(iref2)%tag .eq. refs(iref)%parent_tag) then
        !set parent id
        refs(iref)%parent_id = iref2
        !Update parent's children list
        refs(iref2)%n_chil = refs(iref2)%n_chil + 1
        allocate(temp_chil(refs(iref2)%n_chil))
        if (allocated(refs(iref2)%chil_id)) &
              temp_chil(1:refs(iref2)%n_chil-1) = refs(iref2)%chil_id
        temp_chil(refs(iref2)%n_chil) = refs(iref)%id
        call move_alloc(temp_chil,refs(iref2)%chil_id)
        exit
      endif
    enddo
  
    !if not found a parent
    if (refs(iref)%parent_id .lt. 0) then
      write(msg,'(A,I2,A,I2,A)') 'For reference tag ',refs(iref)%tag, &
                   ' a parent with tag ',refs(iref)%parent_tag,' was not found'
      call error(this_sub_name, this_mod_name, msg)
    endif
    
  enddo

  !transverse the tree to set which things are moving and which are not
  call set_movement(refs, 0, .false.)

  !cleanup
  call finalizeparameters(ref_prs)

end subroutine build_references

!----------------------------------------------------------------------

subroutine destroy_references(refs)
 type(t_ref), allocatable, intent(out)   :: refs(:)

end subroutine destroy_references

!----------------------------------------------------------------------

recursive subroutine set_movement(refs, iref, parent_moving)
 type(t_ref), intent(inout) :: refs(0:)
 integer, intent(in) :: iref
 logical, intent(in) :: parent_moving

 integer :: i

  refs(iref)%moving = parent_moving .or. refs(iref)%self_moving
 
  do i=1,refs(iref)%n_chil
    call set_movement(refs, refs(iref)%chil_id(i), refs(iref)%moving)
  enddo

end subroutine set_movement

!----------------------------------------------------------------------

!> Chech that the inserted references are geometrically correct, correct 
!! where possible, abort if the input is wrong
subroutine check_references(refs)
 type(t_ref), intent(inout) :: refs(0:)

 !TODO TODO TODO
 !check that the orientation is orthogonal

 !chech that the orientation is unitarian, otherwise normalize

 !chech that the rotation axis is unitarian, otherwise normalize

end subroutine check_references

!----------------------------------------------------------------------

!> Update the single reference frame
!!
!! This function is called recursively traversing the reference frames tree
!!
!! First the local movement with respect to the parent is calculated. Then
!! it is combined with the movement coming from the parent with respect to 
!! the base reference, to obtain the transformation of the present frame with 
!! respect to the base frame. Then it is passed to all the children to 
!! update their data
recursive subroutine reference_update_ref(this,refs,t,R,of,G,f)
 class(t_ref), intent(inout) :: this
 type(t_ref), intent(inout) :: refs(0:)
 real(wp), intent(in) :: t
 real(wp), intent(in) :: R(:,:)
 real(wp), intent(in) :: of(:)
 real(wp), intent(in) :: G(:,:)
 real(wp), intent(in) :: f(:)

 real(wp) :: R_loc(3,3), of_loc(3)
 real(wp) :: G_loc(3,3), f_loc(3)
 integer :: i1

  !Update the local transformation with respect to the parent 
  
  call this%update_self(t,R, of,  R_loc, of_loc, G_loc, f_loc)

  !Update the transformation with respect to the base system
  this%R_g  = matmul( R , R_loc )
  this%of_g = matmul( R , of_loc ) + of

  this%G_g = G + matmul(R, G_loc)
  this%f_g = f + matmul(R, f_loc)
  
  
  !For all the reference children call the same updating subroutine, passing
  !this reference global matrices
  do i1 = 1 , this%n_chil
 
    call refs(this%chil_id(i1))%update_ref(refs, t, this%R_g, this%of_g, &
                                                       this%G_g, this%f_g) 
 
  end do 


end subroutine reference_update_ref

!----------------------------------------------------------------------

!> Update the local transformation with respect to the parent reference 
!! frame
!!
!! The local transformation is just the orientation of the local frame if the 
!! frame is not moving, or it is some pre-defined motion law. At the moment 
!! just a constant rotation around a pole and an axis is implemented
subroutine reference_update_self(this, t, R_par, of_par, R_loc, of_loc, &
                                                                 G_loc, f_loc)
 class(t_ref), intent(inout) :: this
 real(wp), intent(in)        :: t
 real(wp), intent(in)        :: R_par(:,:)
 real(wp), intent(in)        :: of_par(:)
 real(wp), intent(out)       :: R_loc(:,:)
 real(wp), intent(out)       :: of_loc(:)
 real(wp), intent(out)       :: G_loc(:,:)
 real(wp), intent(out)       :: f_loc(:)

 real(wp)  :: Psi
 real(wp)  :: R_01_t(3,3)   , R_10_0(3,3)
 real(wp)  :: Omega_vec(3,3)

 character(len=*), parameter :: this_sub_name = 'reference_update_self'

  if (.not.this%self_moving) then 
    
    !not moving with respect to the parent: orientation is the original
    ! one, and motion related matrices are zero
    R_loc = this%frame
    of_loc = this%orig
    G_loc = 0.0_wp
    f_loc = 0.0_wp

  else

    select case (trim(this%mov_type))

     case('constant_rotation')
      !Actual rotation angle
      Psi = this%psi_0 + this%Omega * t
      ! R01(t) -----
      call rot_mat_axis_angle( this%axis, Psi, R_01_t )
      R_loc = matmul(R_01_t, this%frame) 
      ! R10(0) -----
      R_10_0 = transpose(this%frame)
      of_loc = matmul( matmul(R_loc,R_10_0) , this%orig-this%pole) + this%pole

      !Matrix of the vector product of Omega: OmegaX_
      Omega_vec = reshape( &
                 (/0.0_wp, this%Omega*this%axis(3), -this%Omega*this%axis(2), &
                   -this%Omega*this%axis(3), 0.0_wp, this%Omega*this%axis(1), &
                   this%Omega*this%axis(2), -this%Omega*this%axis(1), 0.0_wp/)&
                  ,(/3,3/))
      G_loc = matmul(Omega_vec,transpose(R_par))
      f_loc = matmul(Omega_vec,(-matmul(transpose(R_par),of_par)-this%pole))

     case default

       call error(this_sub_name, this_mod_name, 'Unknown type of movement')
     
     end select

  endif

end subroutine reference_update_self

!----------------------------------------------------------------------

!> Trigger the update of all the reference frame by starting the recursive
!! tree traversing
subroutine update_all_references(refs, t)
 type(t_ref), intent(inout) :: refs(0:)
 real(wp), intent(in) :: t

 real(wp) :: of(3), R(3,3), G(3,3), f(3)
 
 !Give to the first reference identity/null data to start the traversing
 of = (/0.0_wp, 0.0_wp, 0.0_wp/)
 R  = reshape((/1.0_wp, 0.0_wp, 0.0_wp, &
                0.0_wp, 1.0_wp, 0.0_wp, &
                0.0_wp, 0.0_wp, 1.0_wp /),(/3,3/))
 f = (/0.0_wp, 0.0_wp, 0.0_wp/)
 G  = reshape((/0.0_wp, 0.0_wp, 0.0_wp, &
                0.0_wp, 0.0_wp, 0.0_wp, &
                0.0_wp, 0.0_wp, 0.0_wp /),(/3,3/))

 call refs(0)%update_ref(refs, t, R, of, G, f)

end subroutine update_all_references

!----------------------------------------------------------------------

!> Calculate the rotation matrix  generated by the rotation of an angle 
!! around an axis
subroutine rot_mat_axis_angle ( axis, angle, R )

 real(kind=8)     , intent(in)  :: axis(3)
 real(kind=8)     , intent(in)  :: angle 
 real(kind=8)     , intent(out) :: R(3,3)
 real(kind=8)                   :: r0(9) , r1(9) , r2(9)

 real(kind=8) :: nx , ny , nz

 nx = axis(1)
 ny = axis(2)
 nz = axis(3)

 r0 = (/ 1.0d0 , 0.0d0 , 0.0d0 , &
         0.0d0 , 1.0d0 , 0.0d0 , &
         0.0d0 , 0.0d0 , 1.0d0     /)
 r1 = (/ 0.0d0 , - nz  ,   ny  , &
           nz  , 0.0d0 , - nx  , &
         - ny  ,   nx  , 0.0d0     /)

 r2 = (/ nx*nx , nx*ny , nx*nz , &
         ny*nx , ny*ny , ny*nz , &
         nz*nx , nz*ny , nz*nz     /)

 R = reshape ( r0*cos(angle)+r1*sin(angle)+r2*(1-cos(angle)) ,  &
               (/3,3/) , order=(/2,1/) ) 

end subroutine rot_mat_axis_angle

!----------------------------------------------------------------------

end module mod_reference
