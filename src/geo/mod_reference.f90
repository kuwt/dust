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

use mod_math, only: &
  linear_interp, cross

use mod_param, only: &
  wp, eps, max_char_len, nl , pi

use mod_sim_param, only: &
  sim_param

use mod_parse, only: &
  t_parse, getstr, getint, getintarray, getreal, getrealarray, getlogical, countoption, &
  getsuboption, &
  finalizeparameters, t_link, check_opt_consistency
  
use mod_stringtools, only: &
  lowcase

use mod_handling, only: &
  error, warning, info, printout, dust_time, t_realtime

use mod_basic_io, only: &
  read_real_array_from_file

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

public :: t_ref, build_references, update_all_references, destroy_references !, &
!          update_relative_initial_conditions

private

!----------------------------------------------------------------------

!> Reference Frame type
!!
!! Employed to define a reference frame relative to another one
type :: t_ref

  !> Reference id
  integer :: id

  !> Reference tag (can be non-consecutive)
  !integer :: tag
  character(len=max_char_len) :: tag

  !> Reference Name
  !character(len=max_char_len) :: name

  !> Parent reference id, used to index the references array
  integer :: parent_id

  !> Parent tag, can be non consecutive
  !integer :: parent_tag
  character(len=max_char_len) :: parent_tag

  !> Number of children references
  integer :: n_chil

  !> Id of children references (is possible to use pointers here)
  integer, allocatable :: chil_id(:)

  !> Origin in the parent
  real(wp) :: orig(3)

  !> Frame with respect to the parent
  real(wp) :: frame(3,3)

  !> Is the frame moving with respect to the parent?
  logical :: self_moving

  !> Is the reference frame moving ? (might be stationary with respect to a 
  !! moving parent)
  logical :: moving

  !> Moviment type
  character(len=max_char_len) :: mov_type

  !> Is the reference frame the root of multiple reference frames?
  logical :: multiple

  !> Which kind of multiple reference frame is employed
  character(len=max_char_len) :: mult_type

  !> Number of multiple (final) reference frames
  integer :: n_mult

  !> Rotation pole
  real(wp) :: pole(3)

  !> Rotation axis
  real(wp) :: axis(3)

  !> Rotation rate around the axis
  real(wp) :: Omega

  !> Starting rotation angle
  real(wp) :: psi_0

  !> Total offset with respect to the base reference
  real(wp) :: of_g(3)

  !> Total rotation with respect to the base reference
  real(wp) :: R_g(3,3)
  
  !> Total frame velocity with respect to the base reference
  real(wp) :: f_g(3)

  !> Total frame rotation rate with respect to the base reference
  real(wp) :: G_g(3,3)

  !> Angular velocity in the root base reference 
  real(wp) :: angVel_g(3)

  !> General motion arrays
! type(t_motion) :: motion
  !> Position of the origin w.r.t. the position of the pole (at t = 0)
  real(wp), allocatable :: orig_pol_0(:)
  !> Position of the pole
  real(wp), allocatable :: pol_pos(:,:)
  !> Velocity of the pole
  real(wp), allocatable :: pol_vel(:,:)
  !> Time arrays containing the time instants when the motion of the pole is defined
  real(wp), allocatable :: pol_tim(:)
  !> Rotation around the axis
  real(wp), allocatable :: rot_pos(:)
  !> Angular Velocity around the axis
  real(wp), allocatable :: rot_vel(:)
  !> Time arrays containing the time instants when the rotation around the pole is defined
  real(wp), allocatable :: rot_tim(:)

  !> "Local" relative motion: position of the pole and rotation around the axis (used for restart)
  real(wp) :: relative_pol_pos(3)
  real(wp) :: relative_rot_pos


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

 type(t_ref), allocatable :: refs_temp(:)
 type(t_parse) :: ref_prs
 type(t_parse), pointer :: sbprms , sbprms_pol , sbprms_rot , sbprms_hin 
 integer :: n_refs, n_refs_input
 integer :: n_mult_refs
 integer :: iref, iref_input, iref2, i_mult_ref
 integer, allocatable :: temp_chil(:)
 character(len=max_char_len) :: msg
 integer :: prev_id
 type(t_link), pointer :: lnk
!real(wp), allocatable :: psi_0(:)     ! psi_0 *****
 real(wp) :: psi_0
 real(wp) :: hub_offset, norm(3), rot_axis(3), rot_rate
 integer :: n_in, n_mov, n_mult

 character(len=*), parameter :: this_sub_name = 'build_references'

 character(len=max_char_len) :: rot_filen , pol_filen 
 real(wp) , allocatable :: rot_mat(:,:) , pol_mat(:,:)
 integer , allocatable :: pol_fun_int(:)
 real(wp), allocatable :: pol_vec(:) , pol_ome(:) , pol_pha(:) 
 real(wp), allocatable :: pol_off(:) , pol_pos0(:)
 real(wp) :: pol_amp
 integer  :: rot_fun_int
 real(wp) :: rot_amp, rot_ome, rot_pha, rot_off, rot_pos0
 integer :: it , nt
 character(len=max_char_len) :: ref_tag_str
 ! complex multiple components -----
 integer :: i_mult_blades , n_mult_blades , count_dofs
 integer :: hub_id
 integer :: i_dof , n_dofs
 character(len=max_char_len), allocatable :: hinge_type(:)
 real(wp) , allocatable :: hinge_offs(:,:) , hinge_coll(:) , hinge_cyAm(:) , hinge_cyPh(:)

 !> Transient option
 logical                     :: transient_logical
 character(len=max_char_len) :: transient_fun
 real(wp)                    :: transient_time, init_rot_rate
 
 integer :: i1


  !Define all the parameters to be read
  call ref_prs%CreateStringOption('Reference_Tag','Integer tag of reference frame',&
               multiple=.true.)
  call ref_prs%CreateStringOption('Parent_Tag','Tag of the parent reference', &
               multiple=.true.)
  call ref_prs%CreateRealArrayOption('Origin','Origin of reference frame with&
               &respect to the parent', multiple=.true.)
  call ref_prs%CreateRealArrayOption('Orientation','Orientation of reference &
               &frame with respect to the parent', multiple=.true.)
  call ref_prs%CreateLogicalOption('Moving','Is the reference moving', &
               multiple=.true.)
  call ref_prs%CreateLogicalOption('Multiple','Is the reference multiple', &
               multiple=.true.)

  ! Motion sub-parser ---------------------------------------------
  call ref_prs%CreateSubOption('Motion','Definition of the motion of a frame',sbprms, &
               multiple=.true.)

  ! Pole motion sub-parser ----------------------------------------
  call sbprms%CreateSubOption('Pole','Definition of the motion of the pole', &
              sbprms_pol)
  call sbprms_pol%CreateStringOption('Input','Input: velocity or position')
  call sbprms_pol%CreateStringOption('Input_Type','from_file or &
                                                             &simple_function')
  ! TODO: add %CreateStrinArrayOption(...) to options/mod_parse.f90
  call sbprms_pol%CreateIntArrayOption('Function','fun definition. &
                                                        &0:constant,1:sin,...')
  call sbprms_pol%CreateStringOption('File','file .dat containing the motion')
  call sbprms_pol%CreateRealOption('Amplitude','Multiplicative factor &
                                            &for the motion, defult 1.0','1.0')
  call sbprms_pol%CreateRealArrayOption('Vector','Relative amplitude of the &
                        &three coordinates, default 1 1 1','(/1.0, 1.0, 1.0/)')
  call sbprms_pol%CreateRealArrayOption('Omega','Pulsation of the motion for &
                          &each coordinate, default 1 1 1','(/1.0, 1.0, 1.0/)')
  call sbprms_pol%CreateRealArrayOption('Phase','Phase of the motion for &
                          &each coordinate, default 0 0 0','(/0.0, 0.0, 0.0/)')
  call sbprms_pol%CreateRealArrayOption('Offset','Phase of the motion for &
                          &each coordinate, default 0 0 0','(/0.0, 0.0, 0.0/)')
  call sbprms_pol%CreateRealArrayOption('Position_0','Initial position &
                      &for each coordinate, default 0 0 0','(/0.0, 0.0, 0.0/)')
  ! End Pole motion sub-parser ----------------------------------------

  ! Rotation motion sub-parser ------------------------------------
  call sbprms%CreateSubOption('Rotation','Definition of the rotation of &
                              &the frame', sbprms_rot)
  call sbprms_rot%CreateStringOption('Input','Input: velocity or position')
  call sbprms_rot%CreateStringOption('Input_Type','from_file or &
                                                             &simple_function')
  ! TODO: add %CreateStrinArrayOption(...) to options/mod_parse.f90
  call sbprms_rot%CreateIntOption('Function','fun definition.&
                                                    &0:constant,1:sin,...')
  call sbprms_rot%CreateStringOption('File','file .dat containing the motion')
  call sbprms_rot%CreateRealArrayOption('Axis','axis of rotation')
  call sbprms_rot%CreateRealOption('Amplitude','Multiplicative factor for &
                                                  &the motion, default 1','1')
  call sbprms_rot%CreateRealOption('Omega','Pulsation of the motion for &
                                          &each coordinate, default 1','1.0')
  call sbprms_rot%CreateRealOption('Phase','Phase of the motion &
                            &each coordinate, default 0','0.0')
  call sbprms_rot%CreateRealOption('Offset','Offset of the motion &
                        &for each coordinate, default 0','0.0')
  call sbprms_rot%CreateRealOption('Psi_0','Initial position &
                        &for each coordinate, default 0','0.0')
  ! End Rotation motion sub-parser ------------------------------------

  ! End Motion sub-parser ---------------------------------------------
  sbprms => null()

  ! Multiple sub-parser -------------------------------------------
  call ref_prs%CreateSubOption('Multiplicity','Parameters for multiple frames',&
                sbprms, multiple=.true.)
  call sbprms%CreateStringOption('MultType','Kind of multiplicity')
  call sbprms%CreateIntOption('N_Frames', 'Number of reference frames')
  call sbprms%CreateIntOption('N_Blades', 'Number of reference repeated structures,&
               & blades or whatever')
  call sbprms%CreateRealArrayOption('Rot_Axis','Axis of rotation in parent &
               &reference frame')
  call sbprms%CreateRealOption('Hub_Offset','Offset from the pole')
  call sbprms%CreateRealOption('Rot_Rate','Rate of rotation around axis')
  call sbprms%CreateLogicalOption('Transient','Starting transient, .t. or .f.','F')
  call sbprms%CreateStringOption('Transient_Fun', &
                                 'Linear or cosine (1-cosine) transient','linear')
  call sbprms%CreateRealOption('Transient_Time','Time interval of the &
                               &transient','1.0')
  call sbprms%CreateRealOption('Init_Rot_Rate','Initial value of the angular velocity', &
                               '1.0')
  call sbprms%CreateRealOption('Psi_0','Starting rotation angle') 
  call sbprms%CreateIntOption('N_Dofs', 'Number of dofs for each blade')
  ! for complex rotor configurations ---------
  call sbprms%CreateSubOption('Dof','Definition of a hinge dof',&
                 sbprms_hin,multiple=.true.)
  call sbprms_hin%CreateStringOption('Hinge_Type','Char to define the rotation axis: &
                   & it can be Flap, Pitch, Lag')
  call sbprms_hin%CreateRealArrayOption('Hinge_Offset','Position of the hinge')
  call sbprms_hin%CreateRealOption('Collective','Amplitude of the collective motion')
  call sbprms_hin%CreateRealOption('Cyclic_Ampl','Amplitude of the cyclic motion')
  call sbprms_hin%CreateRealOption('Cyclic_Phas','Phase of the cyclic motion')

  ! End Multiple sub-parser -------------------------------------------
 
 
  !read the file
  call ref_prs%read_options(trim(reference_file),printout_val=.true.)

  !Get the number of reference frames and check that all the required params
  !are actually present
  n_refs = countoption(ref_prs,'Reference_Tag') 
  n_refs_input = n_refs  

  n_in = countoption(ref_prs,'Parent_Tag')
  if(n_in .ne. n_refs_input) call error(this_sub_name, this_mod_name, &
   'Inconsistent number of Reference_Tag and Parent_Tag inputs in reference& 
   & frames. Forgot a "Parent_tag = ..." ?')
  n_in = countoption(ref_prs,'Origin')
  if(n_in .ne. n_refs_input) call error(this_sub_name, this_mod_name, &
   'Inconsistent number of Reference_Tag and Origin  inputs in reference& 
   & frames. Forgot a "Origin = ..." ?')
  n_in = countoption(ref_prs,'Orientation')
  if(n_in .ne. n_refs_input) call error(this_sub_name, this_mod_name, &
   'Inconsistent number of Reference_Tag and Orientation inputs in reference& 
   & frames. Forgot a "Orientation = ..." ?')
  n_in = countoption(ref_prs,'Moving')
  if(n_in .ne. n_refs_input) call error(this_sub_name, this_mod_name, &
   'Inconsistent number of Reference_Tag and Moving inputs in reference& 
   & frames. Forgot a "Moving = ..." ?')
  n_in = countoption(ref_prs,'Multiple')
  if(n_in .ne. n_refs_input) call error(this_sub_name, this_mod_name, &
   'Inconsistent number of Reference_Tag and Multiple inputs in reference& 
   & frames. Forgot a "Multiple = ..." ?')

  ! IMPORTANT: the references are allocated starting from zero, zero is the 
  ! base reference, and then all the other references occupy the following
  ! position. Remember to correctly set the numbering starting from zero for
  ! all the other subroutines that employ the references
  allocate(refs(0:n_refs))


  !Setup the base reference
  refs(0)%id = 0
  refs(0)%tag = '0'
  refs(0)%parent_id = -1
  refs(0)%parent_tag = '-1'
  refs(0)%n_chil = 0
  refs(0)%orig = (/0.0_wp, 0.0_wp, 0.0_wp/)
  refs(0)%frame = reshape((/1.0_wp, 0.0_wp, 0.0_wp, &
                            0.0_wp, 1.0_wp, 0.0_wp, &
                            0.0_wp, 0.0_wp, 1.0_wp/),(/3,3/))
  refs(0)%self_moving = .false.
  refs(0)%multiple = .false.
  refs(0)%moving = .false.
  refs(0)%of_g = (/0.0_wp, 0.0_wp, 0.0_wp/)
  refs(0)%R_g = reshape((/1.0_wp, 0.0_wp, 0.0_wp, &
                          0.0_wp, 1.0_wp, 0.0_wp, & 
                          0.0_wp, 0.0_wp, 1.0_wp/),(/3,3/))
  refs(0)%relative_pol_pos = (/ 0.0_wp , 0.0_wp , 0.0_wp /)
  refs(0)%relative_rot_pos = 0.0_wp
  !Setup the other references
  iref = 0; n_mov = 0; n_mult = 0
  do iref_input = 1,n_refs_input

    iref = iref+1

    refs(iref)%id = iref
    refs(iref)%tag = getstr(ref_prs,'Reference_Tag')
    refs(iref)%parent_tag = getstr(ref_prs,'Parent_Tag')
    refs(iref)%n_chil = 0

    refs(iref)%orig  = getrealarray(ref_prs,'Origin',3, olink=lnk)
    call check_opt_consistency(lnk,next=.true.,next_opt='Orientation')
    refs(iref)%frame = reshape(getrealarray(ref_prs,'Orientation',9),(/3,3/))
    !allocated here, will be set in update_all_refs

    refs(iref)%self_moving = getlogical(ref_prs,'Moving', olink=lnk)
    refs(iref)%moving = .false. !standard, will be checked later

    
    if (refs(iref)%self_moving) then
     
      n_mov = n_mov + 1
      call check_opt_consistency(lnk, next=.true., next_opt='Motion')
      call getsuboption(ref_prs,'Motion',sbprms)
      call getsuboption(sbprms,'Pole',sbprms_pol)

      select case( trim(getstr(sbprms_pol,'Input')) )

        case('position')

          if ( countoption(sbprms_pol,'Input_Type') .eq. 0 ) then
            write(ref_tag_str,'(I5)') refs(iref)%tag
            call error(this_sub_name, this_mod_name, '"Input" field not defined &
                  &in Motion={Pole={ for Ref.Frame with Reference_Tag'//trim(ref_tag_str))
          end if

          ! Input type : from_file , simple_function
          select case( trim(getstr(sbprms_pol,'Input_Type')) )    

            case('from_file')
              if ( countoption(sbprms_pol,'File') .eq. 0 ) then
                write(ref_tag_str,'(I5)') refs(iref)%tag
                call error(this_sub_name, this_mod_name, '"File" field not defined &
                      &in Motion={Pole={ for Ref.Frame with Reference_Tag'//trim(ref_tag_str))
              end if 

              pol_filen = trim(getstr(sbprms_pol,'File'))
              call read_real_array_from_file ( 4 , trim(pol_filen) , pol_mat )
              nt = size(pol_mat,1)

              allocate(refs(iref)%pol_pos(3,nt)) 
              allocate(refs(iref)%pol_vel(3,nt)) 
              allocate(refs(iref)%pol_tim(  nt)) 

              ! Read time and position
              refs(iref)%pol_tim = pol_mat(:,1) 
              refs(iref)%pol_pos = transpose(pol_mat(:,2:4))

              ! Check that min(pol_tim) <= min(sim_param%time_vec) -----
              call check_input_from_file( refs(iref)%tag , 'Pole' , &
                                          refs(iref)%pol_tim , sim_param%time_vec ) 

              ! Compute velocity with Finite Difference 
              do it = 1,nt
                if ( it .eq. 1) then
                  refs(iref)%pol_vel(:,it) = ( refs(iref)%pol_pos(:,2) - &
                                               refs(iref)%pol_pos(:,1) ) /  &
                   ( refs(iref)%pol_tim(2) - refs(iref)%pol_tim(1) )
                elseif ( it .lt. nt ) then
                  refs(iref)%pol_vel(:,it) = ( refs(iref)%pol_pos(:,it+1) - &
                                               refs(iref)%pol_pos(:,it-1) ) /  &
                  ( refs(iref)%pol_tim(it+1) - refs(iref)%pol_tim(it-1) )
                elseif ( it .eq. nt ) then
                  refs(iref)%pol_vel(:,nt) = ( refs(iref)%pol_pos(:,nt) - &
                                               refs(iref)%pol_pos(:,nt-1) ) /  &
                    ( refs(iref)%pol_tim(nt) - refs(iref)%pol_tim(nt-1) )
                end if
              end do
 
            case('simple_function')
              if ( countoption(sbprms_pol,'Function') .eq. 0 ) then
                write(ref_tag_str,'(A)') refs(iref)%tag
                call error(this_sub_name, this_mod_name, '"Function" field not defined &
                  &in Motion={Pole={ for Ref.Frame with Reference_Tag'//trim(ref_tag_str))
              end if 

              allocate(refs(iref)%pol_pos(3,sim_param%n_timesteps)) 
              allocate(refs(iref)%pol_vel(3,sim_param%n_timesteps)) 
              allocate(refs(iref)%pol_tim(  sim_param%n_timesteps)) 

              ! Read the integer id for the user defined functions 0:constant,1:sin              
              allocate(pol_fun_int(3), pol_vec(3), pol_ome(3), &
                       pol_pha(3), pol_off(3))
              pol_fun_int = getintarray(sbprms_pol,'Function',3)
              pol_amp = getreal(sbprms_pol,'Amplitude')
              pol_vec = getrealarray(sbprms_pol,'Vector',3)
              pol_ome = getrealarray(sbprms_pol,'Omega',3)
              pol_pha = getrealarray(sbprms_pol,'Phase',3)
              pol_off = getrealarray(sbprms_pol,'Offset',3)
               
              refs(iref)%pol_tim = sim_param%time_vec

              do i1 = 1 , 3
                if ( pol_fun_int(i1) .eq. 0 ) then ! constant function


                  do it = 1 , sim_param%n_timesteps
                    refs(iref)%pol_pos(i1,it) = pol_amp * pol_vec(i1) + pol_off(i1)
                    refs(iref)%pol_vel(i1,it) = 0.0_wp
                  end do

                elseif ( pol_fun_int(i1) .eq. 1 ) then ! sin function

                  do it = 1 , sim_param%n_timesteps 
                    refs(iref)%pol_pos(i1,it) = pol_amp * pol_vec(i1) * &
                       sin( pol_ome(i1) * refs(iref)%pol_tim(it) - pol_pha(i1) ) + &
                       pol_off(i1)
                    refs(iref)%pol_vel(i1,it) = pol_amp * pol_vec(i1) * pol_ome(i1) * &
                       cos( pol_ome(i1) * refs(iref)%pol_tim(it) - pol_pha(i1) )
                  end do

                else
                call error(this_sub_name, this_mod_name, 'Undefined "Function" field &
                  &in Motion={Pole={ for Ref.Frame with Reference_Tag'//trim(ref_tag_str)//&
                  &'. It must be 0:constant or 1:sin')
                end if
              end do
              deallocate(pol_fun_int, pol_vec, pol_ome, pol_pha, pol_off)

            case default
              write(ref_tag_str,'(I5)') refs(iref)%tag
              call error(this_sub_name, this_mod_name, 'Undefined "Input_Type" field &
                  &in Motion={Pole={ for Ref.Frame with Reference_Tag'//trim(ref_tag_str)//&
                  &'. It must be "from_file" or "simple_function".')
     
          
          end select

        case('velocity')

          if ( countoption(sbprms_pol,'Input_Type') .eq. 0 ) then
            write(ref_tag_str,'(I5)') refs(iref)%tag
            call error(this_sub_name, this_mod_name, '"Input" field not defined &
                  &in Motion={Pole={ for Ref.Frame with Reference_Tag'//trim(ref_tag_str))
          end if
          !if ( countoption(sbprms_pol,'Position_0') .eq. 0 ) then
          !  write(ref_tag_str,'(I5)') refs(iref)%tag
          !  call error(this_sub_name, this_mod_name, '"Position_0" field not defined &
          !        &in Motion={Pole={ for Ref.Frame with Reference_Tag'//trim(ref_tag_str))
          !end if

          allocate(pol_pos0(3))
          pol_pos0 = getrealarray(sbprms_pol,'Position_0',3)

          ! Read the integer id for the user defined functions 0:constant,1:sin              
          allocate(pol_fun_int(3), pol_vec(3), pol_ome(3), &
                   pol_pha(3), pol_off(3))
          pol_fun_int = getintarray(sbprms_pol,'Function',3)
          pol_amp = getreal(sbprms_pol,'Amplitude')
          pol_vec = getrealarray(sbprms_pol,'Vector',3)
          pol_ome = getrealarray(sbprms_pol,'Omega',3)
          pol_pha = getrealarray(sbprms_pol,'Phase',3)
          pol_off = getrealarray(sbprms_pol,'Offset',3)
               
          ! Input type : from_file , simple_function
          select case( trim(getstr(sbprms_pol,'Input_Type')) )    

            case('from_file')
              if ( countoption(sbprms_pol,'File') .eq. 0 ) then
                write(ref_tag_str,'(I5)') refs(iref)%tag
                call error(this_sub_name, this_mod_name, '"File" field not defined &
                      &in Motion={Pole={ for Ref.Frame with Reference_Tag'//trim(ref_tag_str))
              end if 

              pol_filen = trim(getstr(sbprms_pol,'File'))
              call read_real_array_from_file ( 4 , trim(pol_filen) , pol_mat )
              nt = size(pol_mat,1)

              allocate(refs(iref)%pol_pos(3,nt)) 
              allocate(refs(iref)%pol_vel(3,nt)) 
              allocate(refs(iref)%pol_tim(  nt)) 

              ! Read time and velocity
              refs(iref)%pol_tim = pol_mat(:,1) 
              refs(iref)%pol_vel = transpose(pol_mat(:,2:4))

              ! Check that min(pol_tim) <= min(sim_param%time_vec) -----
              call check_input_from_file( refs(iref)%tag , 'Pole' , &
                                          refs(iref)%pol_tim , sim_param%time_vec ) 

              ! Compute velocity with Esplicit Euler integration
              refs(iref)%pol_pos(:,1) = pol_pos0  ! Initial condition ...
              do it = 2 , nt
                refs(iref)%pol_pos(:,it) = refs(iref)%pol_pos(:,it-1) + &
                     refs(iref)%pol_vel(:,it-1) * &
                   ( refs(iref)%pol_tim(it) - refs(iref)%pol_tim(it-1) )
              end do
 
            case('simple_function')
              if ( countoption(sbprms_pol,'Function') .eq. 0 ) then
                write(ref_tag_str,'(I5)') refs(iref)%tag
                call error(this_sub_name, this_mod_name, '"Function" field not defined &
                  &in Motion={Pole={ for Ref.Frame with Reference_Tag'//trim(ref_tag_str))
              end if 

              allocate(refs(iref)%pol_pos(3,sim_param%n_timesteps)) 
              allocate(refs(iref)%pol_vel(3,sim_param%n_timesteps)) 
              allocate(refs(iref)%pol_tim(  sim_param%n_timesteps)) 

              refs(iref)%pol_tim = sim_param%time_vec

              do i1 = 1 , 3
                if ( pol_fun_int(i1) .eq. 0 ) then ! constant function

                  do it = 1 , sim_param%n_timesteps
                    refs(iref)%pol_vel(i1,it) = pol_amp * pol_vec(i1)
                    refs(iref)%pol_pos(i1,it) = pol_pos0(i1) + & 
                            pol_amp * pol_vec(i1) * &
                          ( refs(iref)%pol_tim(it) )
! old: "relative    refs(iref)%pol_pos(i1,it) = pol_pos0(i1) + & 
!       time"               pol_amp * pol_vec(i1) * &
!                         ( refs(iref)%pol_tim(it) - refs(iref)%pol_tim(1) )
                  end do

                elseif ( pol_fun_int(i1) .eq. 1 ) then ! sin function

                  do it = 1 , sim_param%n_timesteps 
                    refs(iref)%pol_vel(i1,it) = pol_amp * pol_vec(i1) * &
                       sin( pol_ome(i1) * refs(iref)%pol_tim(it) - pol_pha(i1) ) + &
                       pol_off(i1)
                    if ( pol_ome(i1) .ne. 0.0_wp ) then
                      refs(iref)%pol_pos(i1,it) = - pol_amp * pol_vec(i1) / pol_ome(i1) * &
                       ( cos( pol_ome(i1) * refs(iref)%pol_tim(it) - pol_pha(i1) ) - &
                         cos( - pol_pha(i1) ) ) + &
                         pol_off(i1) * ( refs(iref)%pol_tim(it) ) + &
                         pol_pos0(i1)
! old: "relative      refs(iref)%pol_pos(i1,it) = - pol_amp * pol_vec(i1) / pol_ome(i1) * &
!       time"          ( cos( pol_ome(i1) * refs(iref)%pol_tim(it) - pol_pha(i1) ) - &
!                        cos( pol_ome(i1) * refs(iref)%pol_tim( 1) - pol_pha(i1) ) ) + &
!                        pol_off(i1) * ( refs(iref)%pol_tim(it) - refs(iref)%pol_tim(1) ) + &
!                        pol_pos0(i1)
! old: wrong ?        refs(iref)%pol_pos(i1,it) = - pol_amp * pol_vec(i1) / pol_ome(i1) * &
!                        cos( pol_ome(i1) * refs(iref)%pol_tim(it) - pol_pha(i1) ) + &
!                        pol_off(i1) * ( refs(iref)%pol_tim(it) - refs(iref)%pol_tim(it-1) ) + &
!                        pol_pos0(i1)
                    else 
                      refs(iref)%pol_pos(i1,it) = &
                       ( pol_amp * pol_vec(i1) * sin( - pol_pha(i1) ) + pol_off(i1) ) * &
                       ( refs(iref)%pol_tim(it) ) + &
                         pol_pos0(i1)
! old: "relative      refs(iref)%pol_pos(i1,it) = &
!       time"          ( pol_amp * pol_vec(i1) * sin( - pol_pha(i1) ) + pol_off(i1) ) * &
!                      ( refs(iref)%pol_tim(it) - refs(iref)%pol_tim(1) ) + &
!                        pol_pos0(i1)
! old: wrong ?        refs(iref)%pol_pos(i1,it) = &
!                      ( pol_amp * pol_vec(i1) * sin( - pol_pha(i1) ) + pol_off(i1) ) * &
!                      ( refs(iref)%pol_tim(it) - refs(iref)%pol_tim(it-1) ) + &
!                        pol_pos0(i1)
                    end if
                  end do

                else
                call error(this_sub_name, this_mod_name, 'Undefined "Function" field &
                  &in Motion={Pole={ for Ref.Frame with Reference_Tag'//trim(ref_tag_str)//&
                  &'. It must be 0:constant or 1:sin')
                end if
              end do

            case default
              write(ref_tag_str,'(I5)') refs(iref)%tag
              call error(this_sub_name, this_mod_name, 'Undefined "Input_Type" field &
                  &in Motion={Pole={ for Ref.Frame with Reference_Tag'//trim(ref_tag_str)//&
                  &'. It must be "from_file" or "simple_function".')
     
          end select

          deallocate(pol_fun_int, pol_vec, pol_ome, pol_pha, pol_off, pol_pos0)

        case default
          write(ref_tag_str,'(I5)') refs(iref)%tag
          call error(this_sub_name, this_mod_name, 'Undefined "Input" field &
                &in Motion={Pole={ for Ref.Frame with Reference_Tag'//trim(ref_tag_str)//&
                &'. It must be "position" or "velocity".')

      end select


      call getsuboption(sbprms,'Rotation',sbprms_rot)

      select case(trim(getstr(sbprms_rot,'Input')) )

        case('velocity')

          if ( countoption(sbprms_rot,'Input_Type') .eq. 0 ) then
            write(ref_tag_str,'(A)') trim(refs(iref)%tag)
            call error(this_sub_name, this_mod_name, '"Input" field not defined &
                  &in Motion={Rotation={ for Ref.Frame with Reference_Tag'//trim(ref_tag_str))
          end if
          if ( countoption(sbprms_rot,'Axis') .eq. 0 ) then
            write(ref_tag_str,'(A)') trim(refs(iref)%tag)
            call error(this_sub_name, this_mod_name, '"Axis" field not defined &
                  &in Motion={Rotation={ for Ref.Frame with Reference_Tag'//trim(ref_tag_str))
          end if
          !if ( countoption(sbprms_rot,'Psi_0') .eq. 0 ) then
          !  write(ref_tag_str,'(A)') trim(refs(iref)%tag)
          !  call error(this_sub_name, this_mod_name, '"Psi_0" field not defined &
          !        &in Motion={Rotation={ for Ref.Frame with Reference_Tag'//trim(ref_tag_str))
          !end if

          refs(iref)%axis = getrealarray(sbprms_rot,'Axis',3)  ! %axis is not allocatable
          rot_pos0 = getreal(sbprms_rot,'Psi_0')

          ! Input type : from_file , simple_function
          select case( trim(getstr(sbprms_rot,'Input_Type')) )    

            case('from_file')
              if ( countoption(sbprms_rot,'File') .eq. 0 ) then
                write(ref_tag_str,'(I5)') refs(iref)%tag
                call error(this_sub_name, this_mod_name, '"File" field not defined &
                      &in Motion={Pole={ for Ref.Frame with Reference_Tag'//trim(ref_tag_str))
              end if 


              rot_filen = trim(getstr(sbprms_rot,'File'))
              call read_real_array_from_file ( 2 , trim(rot_filen) , rot_mat )
              nt = size(rot_mat,1)


              allocate(refs(iref)%rot_pos(nt)) 
              allocate(refs(iref)%rot_vel(nt)) 
              allocate(refs(iref)%rot_tim(nt)) 

              ! Read time and velocity
              refs(iref)%rot_tim = rot_mat(:,1) 
              refs(iref)%rot_vel = rot_mat(:,2)

              ! Check that min(pol_tim) <= min(sim_param%time_vec) -----
              call check_input_from_file( refs(iref)%tag , 'Rotation' , &
                                          refs(iref)%pol_tim , sim_param%time_vec ) 

              ! Compute velocity with Esplicit Euler integration
              refs(iref)%rot_pos(1) = rot_pos0  ! Initial condition ...
              do it = 2 , nt
                refs(iref)%rot_pos(it) = refs(iref)%rot_pos(it-1) + &
                     refs(iref)%rot_vel(it-1) * &
                   ( refs(iref)%rot_tim(it) - refs(iref)%rot_tim(it-1) )
              end do
 
            case('simple_function')
              if ( countoption(sbprms_rot,'Function') .eq. 0 ) then
                write(ref_tag_str,'(I5)') refs(iref)%tag
                call error(this_sub_name, this_mod_name, '"Function" field not defined &
                  &in Motion={Rotation={ for Ref.Frame with Reference_Tag'//trim(ref_tag_str))
              end if 

              allocate(refs(iref)%rot_pos(sim_param%n_timesteps)) 
              allocate(refs(iref)%rot_vel(sim_param%n_timesteps)) 
              allocate(refs(iref)%rot_tim(sim_param%n_timesteps)) 

              ! Read the integer id for the user defined functions 0:constant,1:sin              
              rot_fun_int = getint(sbprms_rot,'Function')
              rot_amp = getreal(sbprms_rot,'Amplitude')
              rot_ome = getreal(sbprms_rot,'Omega')
              rot_pha = getreal(sbprms_rot,'Phase')
              rot_off = getreal(sbprms_rot,'Offset')
              
              refs(iref)%rot_tim = sim_param%time_vec

              if ( rot_fun_int .eq. 0 ) then ! constant function

                do it = 1 , sim_param%n_timesteps
                  refs(iref)%rot_vel(it) = rot_amp 
                  refs(iref)%rot_pos(it) = rot_pos0 + rot_amp * &
                        ( refs(iref)%rot_tim(it) )
! old: "relative  refs(iref)%rot_pos(it) = rot_pos0 + rot_amp * &
!       time"           ( refs(iref)%rot_tim(it) - refs(iref)%rot_tim(1) )
                end do

              elseif ( rot_fun_int .eq. 1 ) then ! sin function
                refs(iref)%rot_vel(1) = rot_amp * &
                   sin( rot_ome * refs(iref)%rot_tim(1) - rot_pha ) + rot_off
                refs(iref)%rot_pos(1) = rot_pos0  ! Initial condition ...
                do it = 2 , sim_param%n_timesteps 
                  refs(iref)%rot_vel(it) = rot_amp * &
                     sin( rot_ome * refs(iref)%rot_tim(it) - rot_pha ) + &
                     rot_off
                  if ( rot_ome .ne. 0.0_wp ) then
                    refs(iref)%rot_pos(it) = - rot_amp / rot_ome * &
                      ( cos( rot_ome * refs(iref)%rot_tim(it) - rot_pha ) - &
                        cos( - rot_pha ) ) + &
                       rot_off * ( refs(iref)%rot_tim(it) ) + &
                       rot_pos0
! old: "relative    refs(iref)%rot_pos(it) = - rot_amp / rot_ome * &
!       time"         ( cos( rot_ome * refs(iref)%rot_tim(it) - rot_pha ) - &
!                       cos( rot_ome * refs(iref)%rot_tim( 1) - rot_pha ) ) + &
!                      rot_off * ( refs(iref)%rot_tim(it) - refs(iref)%rot_tim(1) ) + &
!                      rot_pos0
! old: wrong ?      refs(iref)%rot_pos(it) = - rot_amp / rot_ome * &
!                      (cos( rot_ome * refs(iref)%rot_tim(it) - rot_pha ) - &
!                       cos( - rot_pha )) + &
!                      rot_off * ( refs(iref)%rot_tim(it) - refs(iref)%rot_tim(1) ) + &
!                      rot_pos0
                  else 
                    refs(iref)%rot_pos(it) = &
                     ( rot_amp * sin( - rot_pha ) + rot_off ) * &
                     ( refs(iref)%rot_tim(it) ) + &
                       rot_pos0
! old: "relative    refs(iref)%rot_pos(it) = &
!       time"        ( rot_amp * sin( - rot_pha ) + rot_off ) * &
!                    ( refs(iref)%rot_tim(it) - refs(iref)%rot_tim(1) ) + &
!                      rot_pos0
! old: wrong ?      refs(iref)%rot_pos(it) = &
!                    ( rot_amp * sin( - rot_pha ) + rot_off ) * &
!                    ( refs(iref)%rot_tim(it) - refs(iref)%rot_tim(it-1) ) + &
!                      rot_pos0
                  end if
                end do
              else
              call error(this_sub_name, this_mod_name, 'Undefined "Function" field &
                &in Motion={RotatioRotationfor Ref.Frame with Reference_Tag'//trim(ref_tag_str)//&
                &'. It must be 0:constant or 1:sin')
              end if

            case default
              write(ref_tag_str,'(A)') refs(iref)%tag
              call error(this_sub_name, this_mod_name, 'Undefined "Input_Type" field &
                  &in Motion={Rotation={ for Ref.Frame with Reference_Tag'//trim(ref_tag_str)//&
                  &'. It must be "from_file" or "simple_function".')
     
          end select

        case('position')

          if ( countoption(sbprms_rot,'Input_Type') .eq. 0 ) then
            write(ref_tag_str,'(A)') trim(refs(iref)%tag)
            call error(this_sub_name, this_mod_name, '"Input" field not defined &
                  &in Motion={Rotation={ for Ref.Frame with Reference_Tag'//trim(ref_tag_str))
          end if
          if ( countoption(sbprms_rot,'Axis') .eq. 0 ) then
            write(ref_tag_str,'(A)') trim(refs(iref)%tag)
            call error(this_sub_name, this_mod_name, '"Axis" field not defined &
                  &in Motion={Rotation={ for Ref.Frame with Reference_Tag'//trim(ref_tag_str))
          end if
          !if ( countoption(sbprms_rot,'Psi_0') .eq. 0 ) then
          !  write(ref_tag_str,'(A)') trim(refs(iref)%tag)
          !  call error(this_sub_name, this_mod_name, '"Psi_0" field not defined &
          !        &in Motion={Rotation={ for Ref.Frame with Reference_Tag'//trim(ref_tag_str))
          !end if

          refs(iref)%axis = getrealarray(sbprms_rot,'Axis',3) ! %axis is not allocatable
          !rot_pos0 = getreal(sbprms_rot,'Psi_0')

          ! Input type : from_file , simple_function
          select case( trim(getstr(sbprms_rot,'Input_Type')) )    

            case('from_file')
              if ( countoption(sbprms_rot,'File') .eq. 0 ) then
                write(ref_tag_str,'(I5)') refs(iref)%tag
                call error(this_sub_name, this_mod_name, '"File" field not defined &
                      &in Motion={Pole={ for Ref.Frame with Reference_Tag'//trim(ref_tag_str))
              end if 


              rot_filen = trim(getstr(sbprms_rot,'File'))
              call read_real_array_from_file ( 2 , trim(rot_filen) , rot_mat )
              nt = size(rot_mat,1)


              allocate(refs(iref)%rot_pos(nt)) 
              allocate(refs(iref)%rot_vel(nt)) 
              allocate(refs(iref)%rot_tim(nt)) 

              ! Read time and velocity
              refs(iref)%rot_tim = rot_mat(:,1) 
              refs(iref)%rot_pos = rot_mat(:,2)

              ! Check that min(pol_tim) <= min(sim_param%time_vec) -----
              call check_input_from_file( refs(iref)%tag , 'Rotation' , &
                                          refs(iref)%pol_tim , sim_param%time_vec ) 

              ! Compute velocity with Finite Difference 
              do it = 1,nt
                if ( it .eq. 1) then
                  refs(iref)%rot_vel(it) = ( refs(iref)%rot_pos(2) - &
                                               refs(iref)%rot_pos(1) ) /  &
                   ( refs(iref)%rot_tim(2) - refs(iref)%rot_tim(1) )
                elseif ( it .lt. nt ) then
                  refs(iref)%rot_vel(it) = ( refs(iref)%rot_pos(it+1) - &
                                               refs(iref)%rot_pos(it-1) ) /  &
                  ( refs(iref)%rot_tim(it+1) - refs(iref)%rot_tim(it-1) )
                elseif ( it .eq. nt ) then
                  refs(iref)%rot_vel(nt) = ( refs(iref)%rot_pos(nt) - &
                                               refs(iref)%rot_pos(nt-1) ) /  &
                    ( refs(iref)%rot_tim(nt) - refs(iref)%rot_tim(nt-1) )
                end if
              end do
 
            case('simple_function')
              if ( countoption(sbprms_rot,'Function') .eq. 0 ) then
                write(ref_tag_str,'(A)') refs(iref)%tag
                call error(this_sub_name, this_mod_name, '"Function" field not defined &
                  &in Motion={Rotation={ for Ref.Frame with Reference_Tag'//trim(ref_tag_str))
              end if 


              allocate(refs(iref)%rot_pos(sim_param%n_timesteps)) 
              allocate(refs(iref)%rot_vel(sim_param%n_timesteps)) 
              allocate(refs(iref)%rot_tim(sim_param%n_timesteps)) 

              ! Read the integer id for the user defined functions 0:constant,1:sin              
              rot_fun_int = getint(sbprms_rot,'Function')
              rot_amp = getreal(sbprms_rot,'Amplitude')
              rot_ome = getreal(sbprms_rot,'Omega')
              rot_pha = getreal(sbprms_rot,'Phase')
              rot_off = getreal(sbprms_rot,'Offset')
               
              refs(iref)%rot_tim = sim_param%time_vec

              if ( rot_fun_int .eq. 0 ) then ! constant function

                do it = 1 , sim_param%n_timesteps
                  refs(iref)%rot_pos(it) = rot_amp + rot_off
                  refs(iref)%rot_vel(it) = 0.0_wp
                end do

              elseif ( rot_fun_int .eq. 1 ) then ! sin function
                refs(iref)%rot_vel(1) = rot_amp * &
                   sin( rot_ome * refs(iref)%rot_tim(1) - rot_pha ) + rot_off
                refs(iref)%rot_pos(1) = rot_pos0  ! Initial condition ...
                do it = 1 , sim_param%n_timesteps 
                  refs(iref)%rot_pos(it) = rot_amp * &
                     sin( rot_ome * refs(iref)%rot_tim(it) - rot_pha ) + &
                     rot_off
                  refs(iref)%rot_vel(it) = rot_amp * &
                     cos( rot_ome * refs(iref)%rot_tim(it) - rot_pha ) * rot_ome
                end do
              else
              call error(this_sub_name, this_mod_name, 'Undefined "Function" field &
                &in Motion={RotatioRotationfor Ref.Frame with Reference_Tag'//trim(ref_tag_str)//&
                &'. It must be 0:constant or 1:sin')
              end if

            case default
              write(ref_tag_str,'(A)') refs(iref)%tag
              call error(this_sub_name, this_mod_name, 'Undefined "Input_Type" field &
                  &in Motion={Rotation={ for Ref.Frame with Reference_Tag'//trim(ref_tag_str)//&
                  &'. It must be "from_file" or "simple_function".')
     
          end select

        case default
          write(ref_tag_str,'(A)') refs(iref)%tag
          call error(this_sub_name, this_mod_name, 'Undefined "Input" field &
                &in Motion={Pole={ for Ref.Frame with Reference_Tag'//trim(ref_tag_str)//&
                &'. It must be "velocity" or "position".')

      end select

      allocate(refs(iref)%orig_pol_0(3))
      ! CHECK -----> interpolate value at t = 0
      refs(iref)%orig_pol_0 = refs(iref)%orig - refs(iref)%pol_pos(:,1)

      ! Motion sub-parser ---------------------------------------------
      sbprms => null() 
      ! Motion sub-parser ---------------------------------------------

    endif

    ! Motion sub-parser ---------------------------------------------
    sbprms     => null() 
    sbprms_pol => null() 
    sbprms_rot => null() 
    ! Motion sub-parser ---------------------------------------------

    refs(iref)%multiple = getlogical(ref_prs,'Multiple')

    if (refs(iref)%multiple) then
      
      n_mult = n_mult + 1
      call getsuboption(ref_prs,'Multiplicity',sbprms)
      refs(iref)%mult_type = trim(getstr(sbprms,'MultType'))

      select case(trim(refs(iref)%mult_type))
       ! ++++++++++++++++++++++ Simple rotor ++++++++++++++++++++++++
       case('simple_rotor')

        n_mult_refs = getint(sbprms,'N_Frames')
        refs(iref)%n_mult = n_mult_refs

        !1) allocate a series of extra reference frames and move-alloc everything
        n_refs = n_refs+n_mult_refs
        !allocate(refs_temp(0:n_refs+n_mult_refs))
        allocate(refs_temp(0:n_refs))
        refs_temp(0:iref) = refs(0:iref)
        deallocate(refs)
        call move_alloc(refs_temp, refs)


        !2) Read the inputs
!       allocate( psi_0(n_mult_refs))     ! psi_0 *****
        rot_axis = getrealarray(sbprms,'Rot_Axis',3)
        rot_rate = getreal(sbprms,'Rot_Rate')
        psi_0 = getreal(sbprms,'Psi_0')
        hub_offset = getreal(sbprms,'Hub_Offset')


        !3) for each new reference insert all the parameters
        prev_id = refs(iref)%id
        do i_mult_ref=1,n_mult_refs
          iref = iref+1
          refs(iref)%id = iref
          refs(iref)%multiple = .false.
          write(msg,'(A,I2.2)') trim(refs(prev_id)%tag)//'__',i_mult_ref
          refs(iref)%tag = trim(msg)
          refs(iref)%parent_tag = trim(refs(prev_id)%tag)
          refs(iref)%n_chil = 0
          refs(iref)%axis  = rot_axis
          refs(iref)%pole  = (/0.0_wp, 0.0_wp, 0.0_wp/)
          refs(iref)%Omega = rot_rate
          refs(iref)%psi_0 = psi_0 - 2*pi*real(i_mult_ref-1,wp)/&
                                          real(n_mult_refs,wp)
          norm = -cross(refs(iref)%axis,(/0.0_wp, 1.0_wp, 0.0_wp/))
          if (norm2(norm) .le. eps) &
            norm = cross(refs(iref)%axis,(/1.0_wp, 0.0_wp, 0.0_wp/))
          refs(iref)%orig = hub_offset * norm/norm2(norm)
          refs(iref)%frame(1:3,3) = refs(iref)%axis/norm2(refs(iref)%axis)
          refs(iref)%frame(1:3,2) = norm/norm2(norm)
          refs(iref)%frame(1:3,1) = cross(refs(iref)%frame(1:3,2), &
                                          refs(iref)%frame(1:3,3))
          !allocated here, will be set in update_all_refs
          refs(iref)%self_moving = .true.
          refs(iref)%moving = .false. !standard, will be checked later

          allocate(refs(iref)%pol_pos(3,sim_param%n_timesteps)) 
          allocate(refs(iref)%pol_vel(3,sim_param%n_timesteps)) 
          allocate(refs(iref)%pol_tim(  sim_param%n_timesteps)) 
          allocate(refs(iref)%rot_pos(  sim_param%n_timesteps)) 
          allocate(refs(iref)%rot_vel(  sim_param%n_timesteps)) 
          allocate(refs(iref)%rot_tim(  sim_param%n_timesteps)) 
          do it = 1,sim_param%n_timesteps
            refs(iref)%pol_pos(:,it) = refs(iref)%pole
            refs(iref)%pol_vel(:,it) = (/ 0.0_wp , 0.0_wp , 0.0_wp /)
            refs(iref)%pol_tim(  it) = sim_param%time_vec(it) 
            refs(iref)%rot_pos(  it) = refs(iref)%psi_0 + &
                 refs(iref)%Omega * ( sim_param%time_vec(it) )    ! <---- CHECK !!!!
                 !refs(iref)%Omega * ( sim_param%time_vec(it) - sim_param%time_vec(1) )    ! <---- CHECK !!!!
            refs(iref)%rot_vel(  it) = refs(iref)%Omega
            refs(iref)%rot_tim(  it) = sim_param%time_vec(it)
          end do 
          allocate(refs(iref)%orig_pol_0(3))
          refs(iref)%orig_pol_0 = refs(iref)%orig - refs(iref)%pol_pos(:,1)

        enddo
!       deallocate(psi_0)     ! psi_0 *****
       ! ++++++++++++++++++++++     Rotor    ++++++++++++++++++++++++
       case('rotor')
      
         n_mult_blades = getint(sbprms,'N_Blades')
         refs(iref)%n_mult = n_mult_blades
         n_dofs   = getint(sbprms,'N_Dofs')

         if ( n_dofs .lt. 0 ) then
           call warning(this_sub_name, this_mod_name, ' Warning in the input file&
                  & N_Dofs is set to 0, since it was .lt. 0')
           n_dofs = 0
         end if
         ! check if countoption(Dof) .eq. N_Dofs
         count_dofs = countoption(sbprms,'Dof')
         if ( count_dofs .ne. n_dofs ) then
           write(*,*) " countoption (...,'Dof') = " , count_dofs
           write(*,*) " n_dofs                  = " , n_dofs
           call error(this_sub_name, this_mod_name, ' Error in the input file&
                  & of the Reference systems (References.in?): number of Dof&
                  & fields .ne. N_Dofs')
         end if


         !1) allocate a series of extra reference frames and move-alloc everything
         n_refs = n_refs+n_mult_blades*(n_dofs+1)
         !allocate(refs_temp(0:n_refs+n_mult_refs))
         allocate(refs_temp(0:n_refs))
         refs_temp(0:iref) = refs(0:iref)
         deallocate(refs)
         call move_alloc(refs_temp, refs)

         !2) Read the inputs
         rot_axis = getrealarray(sbprms,'Rot_Axis',3)
         rot_rate = getreal(sbprms,'Rot_Rate')
         transient_logical = getlogical(sbprms,'Transient')
         transient_fun     = getstr( sbprms,'Transient_Fun' ) ; call LowCase(transient_fun)
         transient_time    = getreal(sbprms,'Transient_Time')
         init_rot_rate     = getreal(sbprms,'Init_Rot_Rate')
         if ( transient_logical ) then
           if ( ( countoption(sbprms,'Transient_Time') .eq. 0 ) .or. &
                ( countoption(sbprms,'Transient_Fun ') .eq. 0 ) ) then
             call error(this_sub_name, this_mod_name, '"Transient" set as .TRUE.,&
                   & but no "Transient_Time or Transient_Fun is defined. Stop')
           else
             if ( ( trim(transient_fun) .ne. 'linear' ) .and. &
                  ( trim(transient_fun) .ne. 'cosine' ) ) then
               call error(this_sub_name, this_mod_name, '"Transient_Fun" &
                     & must be "Linear" or "Cosine", while input is "'// &
                     trim(transient_fun)//'". Stop')
             end if
           end if
         end if

         psi_0 = getreal(sbprms,'Psi_0')
         hub_offset = getreal(sbprms,'Hub_Offset')

         ! read and allocate some tmp arrays to describe the motion around the hinges
         if ( n_dofs .gt. 0 ) then
           allocate(hinge_type(n_dofs  )) 
           allocate(hinge_offs(n_dofs,3)) ; hinge_offs = 0.0_wp 
           allocate(hinge_coll(n_dofs  )) ; hinge_coll = 0.0_wp
           allocate(hinge_cyAm(n_dofs  )) ; hinge_cyAm = 0.0_wp
           allocate(hinge_cyPh(n_dofs  )) ; hinge_cyPh = 0.0_wp
           do i_dof = 1 , n_dofs
             call getsuboption(sbprms,'Dof',sbprms_hin)
             hinge_type(i_dof  ) = getstr(sbprms_hin,'Hinge_Type') 
             hinge_offs(i_dof,:) = getrealarray(sbprms_hin,'Hinge_Offset',3)
             hinge_coll(i_dof  ) = getreal(sbprms_hin,'Collective') 
             hinge_cyAm(i_dof  ) = getreal(sbprms_hin,'Cyclic_Ampl') 
             hinge_cyPh(i_dof  ) = getreal(sbprms_hin,'Cyclic_Phas') 

             ! from deg to rad
             hinge_coll(i_dof) = hinge_coll(i_dof) * pi / 180.0_wp
             hinge_cyAm(i_dof) = hinge_cyAm(i_dof) * pi / 180.0_wp
             hinge_cyPh(i_dof) = hinge_cyPh(i_dof) * pi / 180.0_wp

           end do
         end if

         !3) for each new reference insert all the parameters
         hub_id = refs(iref)%id
         do i_mult_blades = 1 , n_mult_blades ! loop over the blades
           ! i_dof=0, rotation around the axis
           ! i_dof=1:n_dofs, hinge dofs
           do i_dof = 0 , n_dofs ! loop over the dofs

             prev_id = refs(iref)%id

             iref = iref + 1
             refs(iref)%id = iref 
             refs(iref)%multiple = .false.

             ! tag of the ref.sys.----
             ! for the postpro, the "last" ref.sys. must be <hub_refsys>__<i_mult_blades>
             !TODO: rationalise tag logics (driven by postpro routines)
             if ( i_dof .eq. n_dofs ) then
!              ! double underscore for the last dof
!              write(msg,'(A,I2.2,A,I2.2)') trim(refs(hub_id)%tag)//'__Body',i_mult_blades, &
!                         '_Mov',i_dof
               write(msg,'(A,I2.2)') trim(refs(hub_id)%tag)//'__',i_mult_blades
             else
               write(msg,'(A,I2.2,A,I2.2)') trim(refs(hub_id)%tag)//'_Body',i_mult_blades, &
                          '_Mov',i_dof
             end if
             refs(iref)%tag = trim(msg)

             if ( i_dof .eq. 0 ) then
               refs(iref)%parent_tag = trim(refs(hub_id)%tag)
             else 
               refs(iref)%parent_tag = trim(refs(prev_id)%tag)
             end if

             ! n_chil: initialisation only
             refs(iref)%n_chil = 0 
!            if ( i_dof .eq. n_dofs ) refs(iref)%n_chil = 0 

             ! if ( i_dof .eq. 0 ) ---> rotation of the body around the rotor axis
             !   --> constant rotation
             ! if ( i_dof .ge. 1 ) ---> rotation of the body around the hinges 
             !   --> harmonic rotation
             if ( i_dof .eq. 0 ) then ! constant rotation around the rotor axis 

               refs(iref)%axis  = rot_axis
               refs(iref)%pole  = (/0.0_wp, 0.0_wp, 0.0_wp/)
               refs(iref)%Omega = rot_rate
               refs(iref)%psi_0 = psi_0 - 2*pi*real(i_mult_blades-1,wp) &
                                              /real(n_mult_blades,wp)
               !TODO: check the norm vector to define an origin of the psi angle 
               norm = -cross(refs(iref)%axis,(/0.0_wp, 1.0_wp, 0.0_wp/))
               if (norm2(norm) .le. eps) &
                 norm = cross(refs(iref)%axis,(/1.0_wp, 0.0_wp, 0.0_wp/))
               refs(iref)%orig = hub_offset * norm/norm2(norm)
               refs(iref)%frame(1:3,3) = refs(iref)%axis/norm2(refs(iref)%axis)
               refs(iref)%frame(1:3,2) = norm/norm2(norm)
               refs(iref)%frame(1:3,1) = cross(refs(iref)%frame(1:3,2), &
                                               refs(iref)%frame(1:3,3))
               refs(iref)%self_moving = .true.
               refs(iref)%moving = .false. !standard, will be checked later

               allocate(refs(iref)%pol_pos(3,sim_param%n_timesteps)) 
               allocate(refs(iref)%pol_vel(3,sim_param%n_timesteps)) 
               allocate(refs(iref)%pol_tim(  sim_param%n_timesteps)) 
               allocate(refs(iref)%rot_pos(  sim_param%n_timesteps)) 
               allocate(refs(iref)%rot_vel(  sim_param%n_timesteps)) 
               allocate(refs(iref)%rot_tim(  sim_param%n_timesteps))

               if ( .not. transient_logical ) then ! default: no transient
                 do it = 1,sim_param%n_timesteps
                   refs(iref)%pol_pos(:,it) = refs(iref)%pole
                   refs(iref)%pol_vel(:,it) = (/ 0.0_wp , 0.0_wp , 0.0_wp /)
                   refs(iref)%pol_tim(  it) = sim_param%time_vec(it) 
                   refs(iref)%rot_pos(  it) = refs(iref)%psi_0 + &
                        refs(iref)%Omega * ( sim_param%time_vec(it) )  ! <---- CHECK !!!!
                        !refs(iref)%Omega * ( sim_param%time_vec(it) - sim_param%time_vec(1) )  ! <---- CHECK !!!!
                   refs(iref)%rot_vel(  it) = refs(iref)%Omega
                   refs(iref)%rot_tim(  it) = sim_param%time_vec(it)
                 end do
               else
                 if ( trim(transient_fun) .eq. 'linear' ) then
!                  ! debug ---
!                  write(*,*) ' transient time: ' , transient_time
!                  ! debug ---
                   do it = 1,sim_param%n_timesteps
                     refs(iref)%pol_pos(:,it) = refs(iref)%pole
                     refs(iref)%pol_vel(:,it) = (/ 0.0_wp , 0.0_wp , 0.0_wp /)
                     refs(iref)%pol_tim(  it) = sim_param%time_vec(it) 
                     refs(iref)%rot_tim(  it) = sim_param%time_vec(it)
                     if ( sim_param%time_vec(it) .lt. transient_time ) then
                       refs(iref)%rot_vel(  it) = &
                          ( refs(iref)%Omega - init_rot_rate ) * &
                            sim_param%time_vec(it) / transient_time  &
                          + init_rot_rate
                       refs(iref)%rot_pos(  it) = refs(iref)%psi_0 &
                          + init_rot_rate * sim_param%time_vec(it) &  
                          + 0.5_wp * ( refs(iref)%Omega - init_rot_rate ) * &
                            sim_param%time_vec(it)**2.0_wp / transient_time
                     else
                       refs(iref)%rot_vel(  it) = refs(iref)%Omega
                       refs(iref)%rot_pos(  it) = refs(iref)%psi_0 &
                          + 0.5_wp * ( refs(iref)%Omega + init_rot_rate ) * &
                            transient_time &
                          + refs(iref)%Omega * ( sim_param%time_vec(it) - transient_time )
                     end if
!                    ! debug ---
!                    write(*,*) refs(iref)%rot_vel(it) , refs(iref)%rot_pos(it)
!                    ! debug ---
                   end do
                 elseif ( trim(transient_fun) .eq. 'cosine' ) then
                   do it = 1,sim_param%n_timesteps
                     refs(iref)%pol_pos(:,it) = refs(iref)%pole
                     refs(iref)%pol_vel(:,it) = (/ 0.0_wp , 0.0_wp , 0.0_wp /)
                     refs(iref)%pol_tim(  it) = sim_param%time_vec(it) 
                     refs(iref)%rot_tim(  it) = sim_param%time_vec(it)
                     if ( sim_param%time_vec(it) .lt. transient_time ) then
                       refs(iref)%rot_vel(  it) = &
                            0.5_wp * ( refs(iref)%Omega - init_rot_rate ) * &
                          ( 1.0_wp - cos( pi * sim_param%time_vec(it) / transient_time ) ) &
                          + init_rot_rate
                       refs(iref)%rot_pos(  it) = refs(iref)%psi_0 &
                          + init_rot_rate * sim_param%time_vec(it) &
                          + 0.5_wp * ( refs(iref)%Omega - init_rot_rate ) * &
                          ( sim_param%time_vec(it) - transient_time / pi * &
                            sin( pi*sim_param%time_vec(it)/transient_time) )
                     else
                       refs(iref)%rot_vel(  it) = refs(iref)%Omega
                       refs(iref)%rot_pos(  it) = refs(iref)%psi_0 &
                          + 0.5_wp * ( refs(iref)%Omega + init_rot_rate ) * &
                            transient_time &
                          + refs(iref)%Omega * ( sim_param%time_vec(it) - transient_time )
                     end if
                     ! debug ---
                     write(*,*) refs(iref)%rot_vel(it) , refs(iref)%rot_pos(it)
                     ! debug ---
                   end do
                 end if
               end if
 
               allocate(refs(iref)%orig_pol_0(3))
               refs(iref)%orig_pol_0 = refs(iref)%orig - refs(iref)%pol_pos(:,1)

             else ! harmonic rotation around the hinge axis 

               if ( trim(hinge_type(i_dof)) .eq. 'Flap' ) then
                 refs(iref)%axis  = (/ 1.0_wp , 0.0_wp , 0.0_wp /)
               else if ( trim(hinge_type(i_dof)) .eq. 'Pitch' ) then
                 refs(iref)%axis  = (/ 0.0_wp , 1.0_wp , 0.0_wp /)
               else if ( trim(hinge_type(i_dof)) .eq. 'Lag' ) then
                 refs(iref)%axis  = (/ 0.0_wp , 0.0_wp , 1.0_wp /)
               else
                 call error(this_sub_name, this_mod_name, 'Unknown hinge type:&
                      & it must be: Flap, Pitch ot Lag')
               end if
               refs(iref)%pole  = hinge_offs(i_dof,:) ! (/ 0.0_wp , 0.0_wp , 0.0_wp /) ! *******
               refs(iref)%Omega = rot_rate
               refs(iref)%psi_0 = psi_0 - 2*pi*real(i_mult_blades-1,wp)&
                                              /real(n_mult_blades,wp)

               ! position of the hinge. 
               ! Orientation is described by the motion only (frame = eye(3))
               refs(iref)%orig =  hinge_offs(i_dof,:) ! hinge_offs(i_dof,:) ! *******
               refs(iref)%frame(1:3,1) = (/ 1.0_wp , 0.0_wp , 0.0_wp /)
               refs(iref)%frame(1:3,2) = (/ 0.0_wp , 1.0_wp , 0.0_wp /)
               refs(iref)%frame(1:3,3) = (/ 0.0_wp , 0.0_wp , 1.0_wp /)

               refs(iref)%self_moving = .true.
               refs(iref)%moving = .false. !standard, will be checked later

               allocate(refs(iref)%pol_pos(3,sim_param%n_timesteps)) 
               allocate(refs(iref)%pol_vel(3,sim_param%n_timesteps)) 
               allocate(refs(iref)%pol_tim(  sim_param%n_timesteps)) 
               allocate(refs(iref)%rot_pos(  sim_param%n_timesteps)) 
               allocate(refs(iref)%rot_vel(  sim_param%n_timesteps)) 
               allocate(refs(iref)%rot_tim(  sim_param%n_timesteps)) 

               do it = 1,sim_param%n_timesteps
                 ! position of the pole ---------
                 refs(iref)%pol_pos(:,it) = refs(iref)%pole
                 refs(iref)%pol_vel(:,it) = (/ 0.0_wp , 0.0_wp , 0.0_wp /)
                 refs(iref)%pol_tim(  it) = sim_param%time_vec(it) 
                 ! rotation ---------------------
                 refs(iref)%rot_tim(  it) = sim_param%time_vec(it)
                 if ( refs(iref)%Omega .ne. 0 ) then
                   refs(iref)%rot_pos(  it) = hinge_coll(i_dof) + hinge_cyAm(i_dof) * &
                      cos( refs(iref)%Omega * refs(iref)%rot_tim(it) &
                         + refs(iref)%psi_0 - hinge_cyPh(i_dof) )
                   refs(iref)%rot_vel(  it) = -refs(iref)%Omega * hinge_cyAm(i_dof) * &
                      sin( refs(iref)%Omega * refs(iref)%rot_tim(it) &
                         + refs(iref)%psi_0 - hinge_cyPh(i_dof) )
                 else
                   refs(iref)%rot_pos(  it) = hinge_coll(i_dof) 
                   refs(iref)%rot_vel(  it) = 0.0_wp
                 end if
               end do

               allocate(refs(iref)%orig_pol_0(3))
               refs(iref)%orig_pol_0 = refs(iref)%orig - refs(iref)%pol_pos(:,1)

             end if

   
           end do ! loop over the dofs
           
           !
         end do ! loop over the blades

         if ( n_dofs .gt. 0 ) then
           deallocate(hinge_type) 
           deallocate(hinge_offs)
           deallocate(hinge_coll)
           deallocate(hinge_cyAm)
           deallocate(hinge_cyPh)
         end if


       ! ++++++++++++++++++++++    Default   ++++++++++++++++++++++++
       case default
         call error(this_sub_name, this_mod_name, 'Unknown type of movement')
      end select
      
      ! Motion sub-parser ---------------------------------------------
      sbprms => null() 
      ! Motion sub-parser ---------------------------------------------

    endif
    
  enddo

  n_in = countoption(ref_prs,'Motion')
  if(n_in .ne. n_mov) call error(this_sub_name, this_mod_name, &
   'Inconsistent number of Motion sections and Moving=T references')
  n_in = countoption(ref_prs,'Multiplicity')
  if(n_in .ne. n_mult) call error(this_sub_name, this_mod_name, &
   'Inconsistent number of Multiplicity sections and Multiple=T references')



  !Generate the parent-children references
  do iref = 1,n_refs
    refs(iref)%parent_id = -1
    
    !look for the parent in all the other references (master included)
    do iref2 = 0,n_refs
      !if found
      if (trim(refs(iref2)%tag) .eq. trim(refs(iref)%parent_tag)) then
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
      write(msg,'(A,A,A,A,A)') 'For reference tag ',trim(refs(iref)%tag), &
                   ' a parent with tag ',trim(refs(iref)%parent_tag),' was not found'
      call error(this_sub_name, this_mod_name, msg)
    endif
    
  enddo

  ! check reference frame input
  call check_references(refs)

  !transverse the tree to set which things are moving and which are not
  call set_movement(refs, 0, .false.)

  !cleanup
  call finalizeparameters(ref_prs)

end subroutine build_references

!----------------------------------------------------------------------

subroutine destroy_references(refs)
 type(t_ref), allocatable, intent(out)   :: refs(:)

 integer :: i

 i = size(refs) !dummy operation to avoid warnings

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

 real(wp) :: rot(3,3)
 integer :: i_ref , n_refs

 real(wp) :: det
 real(wp) :: eps = 1.0e-4_wp
 character(len=max_char_len) :: msg

 character(len=*), parameter :: this_sub_name = 'check_references'


 n_refs = ubound(refs,1)

 do i_ref = 0 , n_refs

   ! Orientation matrix
   rot = refs(i_ref)%frame

   ! === check that the orientation is orthogonal ===
   
   !> Unitary determinant
   det = rot(1,1)*rot(2,2)*rot(3,3) &
       + rot(1,2)*rot(2,3)*rot(3,1) &
       + rot(1,3)*rot(2,1)*rot(3,2) &
       - rot(1,1)*rot(2,3)*rot(3,2) &
       - rot(1,2)*rot(2,1)*rot(3,3) &
       - rot(1,3)*rot(2,2)*rot(3,1)
   !> Orthogonal columns
   if ( abs( det - 1.0_wp ) .gt. eps ) then
     write(msg,'(F12.6)') det
     call error(this_sub_name, this_mod_name, 'Input Orientation matrix &
                   &of reference frame "'//trim(refs(i_ref)%tag)//'" has determinant &
                   &equal to det = '//trim(msg)//' instead of (approximately) 1.0')
   end if

   !> Orthogonal columns
   if ( abs( sum( rot(:,1)* rot(:,2) ) ) .gt. eps ) then
     write(msg,'(F12.6)') sum( rot(:,1)* rot(:,2) )
     call error(this_sub_name, this_mod_name, 'Input Orientation matrix &
                   &of reference frame "'//trim(refs(i_ref)%tag)//'" is not &
                   &orthogonal: sum( rot(:,1) * rot(:,2) ) = '//trim(msg))
   end if

   if ( abs( sum( rot(:,1)* rot(:,3) ) ) .gt. eps ) then
     write(msg,'(F12.6)') sum( rot(:,1)* rot(:,3) )
     call error(this_sub_name, this_mod_name, 'Input Orientation matrix &
                   &of reference frame "'//trim(refs(i_ref)%tag)//'" is not &
                   &orthogonal: sum( rot(:,1) * rot(:,3) ) = '//trim(msg))
   end if

   if ( abs( sum( rot(:,2)* rot(:,3) ) ) .gt. eps ) then
     write(msg,'(F12.6)') sum( rot(:,2)* rot(:,3) )
     call error(this_sub_name, this_mod_name, 'Input Orientation matrix &
                   &of reference frame "'//trim(refs(i_ref)%tag)//'" is not &
                   &orthogonal: sum( rot(:,2) * rot(:,3) ) = '//trim(msg))
   end if

   !> Unitary columns
   if ( abs( sum( rot(:,1)* rot(:,1) ) - 1.0_wp ) .gt. eps ) then
     write(msg,'(F12.6)') sum( rot(:,1)* rot(:,1) )
     call error(this_sub_name, this_mod_name, 'Input Orientation matrix &
                   &of reference frame "'//trim(refs(i_ref)%tag)//'" does not have &
                   &unitary columns: sum( rot(:,1) * rot(:,1) ) = '//trim(msg))
   end if
   if ( abs( sum( rot(:,2)* rot(:,2) ) - 1.0_wp ) .gt. eps ) then
     write(msg,'(F12.6)') sum( rot(:,2)* rot(:,2) )
     call error(this_sub_name, this_mod_name, 'Input Orientation matrix &
                   &of reference frame "'//trim(refs(i_ref)%tag)//'" does not have &
                   &unitary columns: sum( rot(:,2) * rot(:,2) ) = '//trim(msg))
   end if
   if ( abs( sum( rot(:,3)* rot(:,3) ) - 1.0_wp ) .gt. eps ) then
     write(msg,'(F12.6)') sum( rot(:,3)* rot(:,3) )
     call error(this_sub_name, this_mod_name, 'Input Orientation matrix &
                   &of reference frame "'//trim(refs(i_ref)%tag)//'" does not have &
                   &unitary columns: sum( rot(:,3) * rot(:,3) ) = '//trim(msg))
   end if


   ! === check that the rotation axis is unitary, otherwise normalize ===
   if ( refs(i_ref)%moving ) then

     if ( norm2(refs(i_ref)%axis - 1.0_wp ) .gt. eps ) then

       write(msg,'(F12.6)') norm2(refs(i_ref)%axis)
       call warning(this_sub_name, this_mod_name, 'Input Rotational axis &
                     &of reference frame "'//trim(refs(i_ref)%tag)//'" had norm2() = &
                     &'//trim(msg)//'. Normalised.')

       refs(i_ref)%axis = refs(i_ref)%axis / norm2(refs(i_ref)%axis)
       
     end if
 
   end if   


  end do

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
recursive subroutine reference_update_ref(this,refs,t,R,of,G,f,angVel)
 class(t_ref), intent(inout) :: this
 type(t_ref), intent(inout) :: refs(0:)
 real(wp), intent(in) :: t
 real(wp), intent(in) :: R(:,:)
 real(wp), intent(in) :: of(:)
 real(wp), intent(in) :: G(:,:)
 real(wp), intent(in) :: f(:)
 real(wp), intent(in) :: angVel(:)

 real(wp) :: R_loc(3,3), of_loc(3)
 real(wp) :: G_loc(3,3), f_loc(3) , angVel_loc(3)
 integer :: i1
  

  !Update the local transformation with respect to the parent 

  call this%update_self(t,R, of,  R_loc, of_loc, G_loc, f_loc, angVel_loc)

  !Update the transformation with respect to the base system
  this%R_g  = matmul( R , R_loc )
  this%of_g = matmul( R , of_loc ) + of

  this%G_g = G + matmul(R, G_loc)
  this%f_g = f + matmul(R, f_loc)
  
  this%angVel_g = angVel + matmul(R, angVel_loc)

  
  !For all the reference children call the same updating subroutine, passing
  !this reference global matrices
  do i1 = 1 , this%n_chil

    call refs(this%chil_id(i1))%update_ref(refs, t, this%R_g, this%of_g, &
                                                    this%G_g, this%f_g, &
                                                    this%angVel_g ) 
 
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
                                                         G_loc,  f_loc, angVel_loc)
 class(t_ref), intent(inout) :: this
 real(wp), intent(in)        :: t
 real(wp), intent(in)        :: R_par(:,:)
 real(wp), intent(in)        :: of_par(:)
 real(wp), intent(out)       :: R_loc(:,:)
 real(wp), intent(out)       :: of_loc(:)
 real(wp), intent(out)       :: G_loc(:,:)
 real(wp), intent(out)       :: f_loc(:)
 real(wp), intent(out)       :: angVel_loc(:)

 real(wp) :: Psi , Omega
 real(wp), allocatable :: xPole(:) , vPole(:)
 real(wp) :: R_01_t(3,3)   , R_10_0(3,3)
 real(wp) :: Omega_vec(3,3)

 character(len=*), parameter :: this_sub_name = 'reference_update_self'

  if (.not.this%self_moving) then 
    
    !not moving with respect to the parent: orientation is the original
    ! one, and motion related matrices are zero
    R_loc  = this%frame
    of_loc = this%orig
    G_loc = 0.0_wp
    f_loc = 0.0_wp
    angVel_loc = 0.0_wp

  else

    ! Actual rotation angle
    call linear_interp( this%rot_pos, this%rot_tim , t , Psi )
    call linear_interp( this%rot_vel, this%rot_tim , t , Omega )
    call linear_interp( this%pol_pos, this%pol_tim , t , xPole )
    call linear_interp( this%pol_vel, this%pol_tim , t , vPole )

    this%relative_pol_pos = xPole
    this%relative_rot_pos = Psi

    ! R01(t) -----
    call rot_mat_axis_angle( this%axis, Psi, R_01_t )
    R_loc = matmul(R_01_t, this%frame) 
    ! R10(0) -----
    R_10_0 = transpose(this%frame)
    of_loc = matmul( matmul(R_loc,R_10_0) , this%orig_pol_0) + xPole

    ! Matrix of the vector product of Omega: OmegaX_
    Omega_vec = reshape( &
               (/             0.0_wp,  Omega*this%axis(3), -Omega*this%axis(2), &
                 -Omega*this%axis(3),              0.0_wp,  Omega*this%axis(1), &
                  Omega*this%axis(2), -Omega*this%axis(1),              0.0_wp/)&
                ,(/3,3/))
    G_loc = matmul(Omega_vec,transpose(R_par))
    f_loc = matmul(Omega_vec,(-matmul(transpose(R_par),of_par)- xPole))  + &
                      vPole

    AngVel_loc = Omega*this%axis 

  endif

end subroutine reference_update_self

!----------------------------------------------------------------------

!> Trigger the update of all the reference frame by starting the recursive
!! tree traversing
subroutine update_all_references(refs, t)
 type(t_ref), intent(inout) :: refs(0:)
 real(wp), intent(in) :: t

 real(wp) :: of(3), R(3,3), G(3,3), f(3) , angVel(3)
 
 !Give to the first reference identity/null data to start the traversing
 of = (/0.0_wp, 0.0_wp, 0.0_wp/)
 R  = reshape((/1.0_wp, 0.0_wp, 0.0_wp, &
                0.0_wp, 1.0_wp, 0.0_wp, &
                0.0_wp, 0.0_wp, 1.0_wp /),(/3,3/))
 f  = (/0.0_wp, 0.0_wp, 0.0_wp/)
 G  = reshape((/0.0_wp, 0.0_wp, 0.0_wp, &
                0.0_wp, 0.0_wp, 0.0_wp, &
                0.0_wp, 0.0_wp, 0.0_wp /),(/3,3/))
 angVel = (/ 0.0_wp , 0.0_wp , 0.0_wp /)

 call refs(0)%update_ref(refs, t, R, of, G, f,angVel)

end subroutine update_all_references

!----------------------------------------------------------------------

!> Calculate the rotation matrix  generated by the rotation of an angle 
!! around an axis
subroutine rot_mat_axis_angle ( axis, angle, R )

 real(wp)     , intent(in)  :: axis(3)
 real(wp)     , intent(in)  :: angle 
 real(wp)     , intent(out) :: R(3,3)
 real(wp)                   :: r0(9) , r1(9) , r2(9)

 real(wp) :: nx , ny , nz

 nx = axis(1)
 ny = axis(2)
 nz = axis(3)

 r0 = (/ 1.0_wp , 0.0_wp , 0.0_wp , &
         0.0_wp , 1.0_wp , 0.0_wp , &
         0.0_wp , 0.0_wp , 1.0_wp     /)
 r1 = (/ 0.0_wp , - nz  ,   ny  , &
           nz  , 0.0_wp , - nx  , &
         - ny  ,   nx  , 0.0_wp     /)

 r2 = (/ nx*nx , nx*ny , nx*nz , &
         ny*nx , ny*ny , ny*nz , &
         nz*nx , nz*ny , nz*nz     /)

 R = reshape ( r0*cos(angle)+r1*sin(angle)+r2*(1.0_wp-cos(angle)) ,  &
               (/3,3/) , order=(/2,1/) ) 

end subroutine rot_mat_axis_angle

!----------------------------------------------------------------------

subroutine check_input_from_file( ref_tag_str , pol_rot_str , time_from_file , sim_param_time )
 character(len=*)  , intent(in) :: ref_tag_str
 character(len=*)  , intent(in) :: pol_rot_str
 real(wp) , intent(in) :: time_from_file(:) , sim_param_time(:)

 character(len=*), parameter :: this_sub_name = 'check_input_from_file'

 if ( time_from_file(1) .gt. sim_param_time(1) ) then
    write(*,*) ' input time_vec(1)     : ' , time_from_file(1) 
    write(*,*) ' sim_param%time_vec(1) : ' , sim_param_time(1) 
    call error(this_sub_name, this_mod_name, 'Wrong input file &
                  &in Motion={'//trim(pol_rot_str)// &
                   '={ for Ref.Frame with Reference_Tag'//trim(ref_tag_str)//&
                  &'. Initial time value > sim_param%time_vec(1).')
 end if
 if ( time_from_file(size(time_from_file)) .gt. sim_param_time(size(sim_param_time)) ) then
    write(*,*) ' input time_vec(end)     : ' , time_from_file(size(time_from_file)) 
    write(*,*) ' sim_param%time_vec(end) : ' , sim_param_time(size(sim_param_time)) 
    call error(this_sub_name, this_mod_name, 'Wrong input file &
                  &in Motion={'//trim(pol_rot_str)// &
                   '={ for Ref.Frame with Reference_Tag'//trim(ref_tag_str)//&
                  &'. Final time value < sim_param%time_vec(end).')
 end if

end subroutine check_input_from_file

!----------------------------------------------------------------------

! subroutine update_relative_initial_conditions (restart_file, ref_file , refs )
!  character(len=max_char_len), intent(in) :: restart_file
!  character(len=max_char_len), intent(in) :: ref_file
!  type(t_ref), intent(inout) :: refs(0:)
! 
!  type(t_parse) :: ref_prs
!  type(t_parse), pointer :: sbprms , sbprms_pol , sbprms_rot 
! 
!  integer(h5loc) :: floc , refs_gloc , ref_loc !, gloc1, gloc2, gloc3
!  character(len=max_char_len) :: ref_title
! 
!  real(wp) :: relative_pos_0(3)
!  real(wp) :: relative_rot_0
!  real(wp) :: relative_rot_01
! 
!  logical :: multiple
! 
!  integer :: iref
!  integer :: nDofs , nBlades
!  character(len=max_char_len) :: multType
!  integer :: i , it , nref_ref_in , ib , idof
! 
!  character(len=*), parameter :: this_sub_name = 'update_relative_initial_conditions'
! 
!  !Define all the parameters to be read
!  call ref_prs%CreateStringOption('Reference_Tag','Integer tag of reference frame',&
!               multiple=.true.)
!  call ref_prs%CreateLogicalOption('Moving','Is the reference moving', &
!               multiple=.true.)
!  call ref_prs%CreateLogicalOption('Multiple','Is the reference multiple', &
!               multiple=.true.)
! 
!  ! Motion sub-parser ---------------------------------------------
!  call ref_prs%CreateSubOption('Motion','Definition of the motion of a frame',sbprms, &
!               multiple=.true.)
!  ! Pole motion sub-parser ----------------------------------------
!  call sbprms%CreateSubOption('Pole','Definition of the motion of the pole', &
!              sbprms_pol)
!  call sbprms_pol%CreateStringOption('Input','Input: velocity or position')
!  ! End Pole motion sub-parser ----------------------------------------
!  ! Rotation motion sub-parser ------------------------------------
!  call sbprms%CreateSubOption('Rotation','Definition of the rotation of &
!                              &the frame', sbprms_rot)
!  call sbprms_rot%CreateStringOption('Input','Input: velocity or position')
!  ! End Rotation motion sub-parser ------------------------------------
!  ! End Motion sub-parser ---------------------------------------------
!  sbprms => null()
! 
!  ! Multiple sub-parser -------------------------------------------
!  call ref_prs%CreateSubOption('Multiplicity','Parameters for multiple frames',&
!                sbprms, multiple=.true.)
!  call sbprms%CreateStringOption('MultType','Kind of multiplicity')
!  call sbprms%CreateIntOption('N_Frames', 'Number of reference frames')
!  call sbprms%CreateIntOption('N_Blades', 'Number of reference repeated structures,&
!               & blades or whatever')
!  call sbprms%CreateIntOption('N_Dofs', 'Number of dofs for each blade')
!  ! End Multiple sub-parser -------------------------------------------
! 
! 
!  !read the file
!  call ref_prs%read_options(trim(ref_file),printout_val=.false.)
! 
!  ! debug -----
!  write(*,*) ' shape(refs) : ' , shape(refs)
!  ! debug -----
!  
!  nref_ref_in = countoption(ref_prs,'Reference_Tag') 
! 
!  ! open restart file to read the relative initial conditions
!  call open_hdf5_file(trim(restart_file), floc)
!  call open_hdf5_group(floc, 'References', refs_gloc)
! 
!  iref = 0
!  
!  do i = 1 , nref_ref_in
!     
!    iref = iref + 1
!    ! debug ----
!    write(*,*) ' refs(',iref,')%tag : ' , trim(refs(iref)%tag)
!    ! debug ----
! 
!    write(ref_title,'(A,I3.3)')'Ref',iref
!    call open_hdf5_group(refs_gloc, trim(ref_title), ref_loc)
! 
!    if ( refs(iref)%self_moving ) then
! 
!      call getsuboption(ref_prs,'Motion',sbprms)
! 
!      ! === Pole ===
!      call getsuboption(sbprms,'Pole',sbprms_pol)
! 
!      select case( trim(getstr(sbprms_pol,'Input')) )
!         case('velocity')
! !         ! add the i.c. from restart file
! !         write(*,*) ' pole velocity, iref : ' , iref
! 
!           call read_hdf5(relative_pos_0,'RelativePolPos',ref_loc)
! 
!           do it = 1 , size(refs(iref)%pol_pos,2)
!             refs(iref)%pol_pos(:,it) = refs(iref)%pol_pos(:,it) + &
!                                        relative_pos_0
!           end do
! 
!         case default
!           ! do nothing
!       end select
! 
!       ! === Rotation ===
!       call getsuboption(sbprms,'Rotation',sbprms_rot)
! 
!       select case(trim(getstr(sbprms_rot,'Input')) )
!         case('velocity')
! !         ! add the i.c. from restart file
! !         write(*,*) ' rot  velocity, iref : ' , iref
! 
!           call read_hdf5(relative_rot_0,'RelativeRotPos',ref_loc)
! 
!           do it = 1 , size(refs(iref)%rot_pos)
!             refs(iref)%rot_pos(it) = refs(iref)%rot_pos(it) + &
!                                      relative_rot_0
!           end do
! 
!         case default
!           ! do nothing
!       end select
! 
!       sbprms_rot => null()
!       sbprms_pol => null()
!       sbprms     => null()
! 
!    end if
! 
!    call close_hdf5_group(ref_loc)
! 
! ! Multiplicity -----
!    multiple = getlogical(ref_prs,'Multiple')
! 
!    if ( multiple ) then ! loop over multiplicity
! 
!      call getsuboption(ref_prs,'Multiplicity',sbprms)
!      
!      multType = trim(getstr(sbprms,'MultType'))
!      nBlades  = getint(sbprms,'N_Blades') 
!      nDofs    = getint(sbprms,'N_dofs')
! 
!      ! debug ----
!      write(*,*) ' multType , nBlades , nDofs ' , trim(multType) , nBlades , nDofs
!      ! debug ----
!  
!      select case(trim(multType) )
!  
!        case('simple rotor')
!          ! to be implemented .... 
! 
!        case('rotor')
!          ! multiple ref.sys. defined by means of offset and rotation ---------------
!          ! read only the rotation offset of the first ref.sys. of the first blade
!          ! and add to the first ref.sys. of all the other blades
! 
!          write(ref_title,'(A,I3.3)')'Ref',iref+1
!          call open_hdf5_group(refs_gloc, trim(ref_title), ref_loc)
!          call read_hdf5(relative_rot_01,'RelativeRotPos',ref_loc)
!          call close_hdf5_group(ref_loc)
!          ! debug ----
!          write(*,*) trim(ref_title)
!          ! debug ----
! 
!          
! 
!          do ib = 1 , nBlades
!  
!            ! Update first ref.sys. of the blade (~ first hinge)
!            iref = iref + 1
!            write(ref_title,'(A,I3.3)')'Ref',iref
!            ! debug ----
!            write(*,*) ' Update : ' , trim(ref_title) , relative_rot_01 
!            ! debug ----
!            do it = 1 , size(refs(iref)%rot_pos)
!              refs(iref)%rot_pos(it) = refs(iref)%rot_pos(it) + &
!                                       relative_rot_01
!            end do
!  
!            ! Update the following ref.sys. of rotors 
!            do idof = 1 , nDofs
!    
!              iref = iref + 1
!              write(ref_title,'(A,I3.3)')'Ref',iref
! !            call open_hdf5_group(refs_gloc, trim(ref_title), ref_loc)
! !            call read_hdf5(relative_rot_0,'RelativeRotPos',ref_loc)
! !            call close_hdf5_group(ref_loc)
!              ! debug ----
!              write(*,*) ' Update : ' , trim(ref_title)  , relative_rot_0 
!              ! debug ----
! 
!              do it = 1 , size(refs(iref)%rot_pos)
!                refs(iref)%rot_pos(it) = refs(iref)%rot_pos(it) + &
!                                         relative_rot_0
!              end do
!    
!            end do 
!   
!          end do 
!  
!  
!      end select
! 
!      sbprms => null()
!  
!    end if
! 
! 
!  end do
! 
!  call close_hdf5_group(refs_gloc)
!  call close_hdf5_file(floc)
! 
! ! ! stop ----
! !  write(*,*) ' stop l.1714 update_relative_initial_conditions '
! !  stop
! ! ! stop ----
! 
! end subroutine update_relative_initial_conditions

!----------------------------------------------------------------------

end module mod_reference
