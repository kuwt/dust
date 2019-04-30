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


!> Module to handle the input/output of the solution of DUST in hdf5 format,
!! meant principally to exchange data with the post processor and restart
module mod_dust_io

use mod_param, only: &
  wp, max_char_len, nl

use mod_sim_param, only: &
  sim_param

use mod_handling, only: &
  error, warning, info, printout, dust_time, t_realtime, check_preproc

!use mod_aero_elements, only: &
!  c_elem, t_elem_p

use mod_aeroel, only: &
  c_elem, c_pot_elem, c_vort_elem, c_impl_elem, c_expl_elem, &
  t_elem_p, t_pot_elem_p, t_vort_elem_p, t_impl_elem_p, t_expl_elem_p

use mod_surfpan, only: &
  t_surfpan

use mod_vortlatt, only: &
  t_vortlatt

use mod_liftlin, only: &
  t_liftlin

use mod_reference, only: &
  t_ref
  
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

use mod_geometry, only: &
  t_geo, t_geo_component

use mod_wake, only: &
  t_wake

use mod_stringtools, only: &
  stricmp

use mod_version, only: &
  git_sha1, version


use mod_parse, only: &
  t_parse, getstr, getint, getintarray, getreal, getrealarray, getlogical, countoption, &
  getsuboption,&
  finalizeparameters, t_link, check_opt_consistency


!----------------------------------------------------------------------

implicit none

public :: save_status, load_solution, load_time

private


character(len=*), parameter :: this_mod_name = 'mod_dust_io'

!----------------------------------------------------------------------

contains

!----------------------------------------------------------------------
subroutine save_status(geo, wake, it, time, run_id)
 type(t_geo), intent(in)         :: geo
 type(t_wake), intent(in) :: wake
 integer, intent(in)             :: it
 real(wp), intent(in)            :: time
 integer, intent(in)             :: run_id(10)

 integer(h5loc) :: floc, gloc1, gloc2, gloc3, ploc
 character(len=max_char_len) :: comp_name, sit
 integer :: icomp, ncomp
 integer :: iref, nref
 character(len=max_char_len) :: ref_name
 integer :: ie, ne
 real(wp), allocatable :: vort(:), cp(:) , pres(:)
 real(wp), allocatable :: dforce(:,:), surf_vel(:,:)
 real(wp), allocatable :: points_w(:,:,:), cent(:,:,:) , vel_w(:,:,:)
 real(wp), allocatable :: vort_v(:,:)
 integer, allocatable :: conn_pe(:)
 integer :: ir, id, ip, np

  !create the output file
  write(sit,'(I4.4)') it
  call new_hdf5_file(trim(sim_param%basename)//'_res_'//trim(sit)//'.h5', &
                     floc)
  call write_hdf5_attr(run_id, 'run_id', floc)
  call write_hdf5_attr(git_sha1, 'git_sha1', floc)
  call write_hdf5_attr(version, 'version', floc)
  call sim_param%save_param(floc)
  call write_hdf5(time,'time',floc)

  call new_hdf5_group(floc, 'Parameters', ploc)
  call write_hdf5(sim_param%u_inf,'u_inf', ploc)
  call write_hdf5(sim_param%P_inf,'P_inf', ploc)
  call write_hdf5(sim_param%rho_inf,'rho_inf', ploc)
  call close_hdf5_group(ploc)
  

  ! 1) %%%%% Component solution: 
  ! mimic the structure of the components inside the geometry input file

  call new_hdf5_group(floc, 'Components', gloc1)
  ncomp = size(geo%components)
  call write_hdf5(ncomp, 'NComponents', gloc1)
  !Cycle the components
  do icomp = 1, ncomp
    !write(comp_name,'(A,I3.3)')'Comp',icomp
    write(comp_name,'(A,I3.3)')'Comp',geo%components(icomp)%comp_id
    call new_hdf5_group(gloc1, trim(comp_name), gloc2)

    call write_hdf5(trim(geo%components(icomp)%comp_name),'CompName',gloc2)
    call write_hdf5(trim(geo%components(icomp)%ref_tag),'RefTag',gloc2)
    call write_hdf5(geo%components(icomp)%ref_id,'RefId',gloc2)
    
    call new_hdf5_group(gloc2, 'Solution', gloc3)

    ne = size(geo%components(icomp)%el)
    allocate(vort(ne), cp(ne) , pres(ne) , dforce(3,ne) )
    do ie = 1,ne
     vort(ie) = geo%components(icomp)%el(ie)%mag
     !cp(ie) = geo%components(icomp)%el(ie)%cp
     pres(ie) = geo%components(icomp)%el(ie)%pres
     dforce(:,ie) = geo%components(icomp)%el(ie)%dforce
    enddo
    call write_hdf5(vort,'Vort',gloc3)
    !call write_hdf5(cp,'Cp',gloc3)
    call write_hdf5(pres,'Pres',gloc3)
    call write_hdf5(dforce,'dF',gloc3)
    deallocate(vort, cp, pres, dforce)
    
    !Output the surface velocity
    if ( trim( geo%components(icomp)%comp_el_type ) .eq. 'p' ) then
      allocate(surf_vel(3,ne))    
      do ie = 1,ne
        select type( el => geo%components(icomp)%el(ie) ) ; type is (t_surfpan)
          surf_vel(:,ie) = el%surf_vel
        end select
      end do 
      call write_hdf5(surf_vel,'surf_vel',gloc3)
      deallocate(surf_vel)
    endif

    call close_hdf5_group(gloc3)
    call close_hdf5_group(gloc2)
  enddo
  call close_hdf5_group(gloc1)
  
  ! 2) %%%% Wake:
  ! just print the whole points, the solution and the starting row
  ! connectivity to build connectivity after
  !=== Panels ===
  call new_hdf5_group(floc, 'PanelWake', gloc1)

  !call write_hdf5(wake_pan%w_points(:,:,1:wake_pan%wake_len+1),'WakePoints',gloc1 )
  !call write_hdf5(wake_pan%i_start_points,'StartPoints',gloc1)
  !call write_hdf5(wake_pan%idou(:,1:wake_pan%wake_len),'WakeVort',gloc1)
  call write_hdf5(wake%pan_w_points(:,:,1:wake%pan_wake_len+1),'WakePoints',gloc1 )
! ! debug ----
! write(*,*) ' shape(wake%pan_w_points(   :,:,:))                     : ' , &
!              shape(wake%pan_w_points(   :,:,:))
! write(*,*) ' shape(wake%pan_w_points(   :,:,1:wake%pan_wake_len+1)) : ' , &
!              shape(wake%pan_w_points(   :,:,1:wake%pan_wake_len+1))
! write(*,*) ' shape(wake%pan_w_vel(   :,:,:))                        : ' , &
!              shape(wake%pan_w_vel(   :,:,:))
! write(*,*) ' shape(wake%pan_w_vel(   :,:,1:wake%pan_wake_len+1))    : ' , &
!              shape(wake%pan_w_vel(   :,:,1:wake%pan_wake_len+1))
! ! debug ----
  call write_hdf5(wake%pan_w_vel(   :,:,1:wake%pan_wake_len+1),'WakeVels'  ,gloc1 ) ! <<<< restart with Bernoulli integral equation
  call write_hdf5(wake%i_start_points,'StartPoints',gloc1)
  call write_hdf5(wake%pan_idou(:,1:wake%pan_wake_len),'WakeVort',gloc1)

  call close_hdf5_group(gloc1)

  !=== Rings ===
  call new_hdf5_group(floc, 'RingWake', gloc1)
  allocate(points_w(3,wake%np_row,wake%rin_wake_len))
  allocate(conn_pe(wake%np_row))
  allocate(cent(3,wake%ndisks,wake%rin_wake_len))
  ip = 1
  do id = 1,wake%ndisks
    np = wake%rin_gen_elems(id)%p%n_ver
    do ir = 1,wake%rin_wake_len
      points_w(:,ip:ip+np-1,ir) = wake%wake_rings(id,ir)%ver(:,:)
      cent(:,id,ir) = wake%wake_rings(id,ir)%cen(:)
    enddo
    conn_pe(ip:ip+np-1) = id
    ip = ip+np
  enddo
  call write_hdf5(points_w,'WakePoints',gloc1)
  call write_hdf5(conn_pe,'Conn_pe',gloc1)
  call write_hdf5(cent,'WakeCenters',gloc1)
  call write_hdf5(wake%rin_idou(:,1:wake%rin_wake_len),'WakeVort',gloc1)
  call close_hdf5_group(gloc1)
  deallocate(points_w, conn_pe, cent)

  !=== Particles ===
  call new_hdf5_group(floc, 'ParticleWake', gloc1)
  allocate(points_w(3,wake%n_prt,1))
  allocate(vort_v(3,wake%n_prt))
  allocate(vel_w(3,wake%n_prt,1))
  do ip = 1, wake%n_prt
    points_w(:,ip,1) = wake%part_p(ip)%p%cen
    vel_w(:,ip,1) = wake%part_p(ip)%p%vel
    vort_v(:,ip) = wake%part_p(ip)%p%dir * wake%part_p(ip)%p%mag
  enddo
  call write_hdf5(points_w(:,:,1),'WakePoints',gloc1)
  call write_hdf5(   vel_w(:,:,1),'WakeVels'  ,gloc1)
  
  call write_hdf5(vort_v,'WakeVort',gloc1)
  call write_hdf5(wake%last_pan_idou,'LastPanIdou',gloc1)
  call close_hdf5_group(gloc1)
  deallocate(points_w, vort_v, vel_w)

  ! 3) %%%% References
  ! save the whole list of references
  call new_hdf5_group(floc, 'References', gloc1)
  nref = size(geo%refs)
  call write_hdf5(nref, 'NReferences', gloc1)
  do iref = 0,nref-1
    write(ref_name,'(A,I3.3)')'Ref',iref
    call new_hdf5_group(gloc1, trim(ref_name), gloc2)
    
    call write_hdf5(geo%refs(iref)%tag, 'Tag', gloc2)
    call write_hdf5(geo%refs(iref)%of_g, 'Offset',gloc2)
    call write_hdf5(geo%refs(iref)%R_g, 'R',gloc2)
    call write_hdf5(geo%refs(iref)%f_g, 'Vel',gloc2)
    call write_hdf5(geo%refs(iref)%G_g, 'RotVel',gloc2)
    call write_hdf5(geo%refs(iref)%relative_pol_pos, 'RelativePolPos',gloc2)
    call write_hdf5(geo%refs(iref)%relative_rot_pos, 'RelativeRotPos',gloc2)

    call close_hdf5_group(gloc2)
  enddo

  call close_hdf5_group(gloc1)

  !close the output file
  call close_hdf5_file(floc)


end subroutine save_status

!----------------------------------------------------------------------
!> Load the solution from a pre-existent solution file
!!
!! Loads just the bodies solution, which is stored into the components,
!! the loading of the wakes is left to other subroutines
subroutine load_solution(filename,comps,refs)
 character(len=max_char_len), intent(in) :: filename
 type(t_geo_component), intent(inout) :: comps(:)
 type(t_ref), intent(in) :: refs(0:)

 integer(h5loc) :: floc, gloc1, gloc2, gloc3
 integer :: ncomp, icomp, icomp2
 integer :: ne, ie
 character(len=max_char_len) :: comp_name_read, comp_id
 real(wp), allocatable :: idou(:), pres(:), dF(:,:)
 

 character(len=*), parameter :: this_sub_name = 'load_solution'

  call open_hdf5_file(filename, floc)

  !TODO: check the run id with the geo file
  call open_hdf5_group(floc, 'Components', gloc1)
  
  !Check if the number of components is consistent
  call read_hdf5(ncomp, 'NComponents', gloc1)
  if(ncomp .ne. size(comps)) call error(this_sub_name, this_mod_name, &
  'Mismatching number of components between geometry and solution files')

  !Cycle on all the components on the file and then on all the local
  !components: they come from different files and can be mixed up in the
  !order
  do icomp = 1, ncomp
    write(comp_id,'(A,I3.3)')'Comp',icomp
    call open_hdf5_group(gloc1, trim(comp_id), gloc2)

    call read_hdf5(comp_name_read,'CompName',gloc2)
    do icomp2 = 1,size(comps)
      if (stricmp(comp_name_read, comps(icomp2)%comp_name)) then
        !We found the correct local component, read all
        call check_ref(gloc2, floc, refs(comps(icomp2)%ref_id))
        call open_hdf5_group(gloc2, 'Solution', gloc3)
        call read_hdf5_al(idou,'Vort',gloc3)
        call read_hdf5_al(pres,'Pres',gloc3)
        call read_hdf5_al(dF,'dF',gloc3)
        call close_hdf5_group(gloc3)
        ne = size(comps(icomp2)%el)
        do ie =1,ne
          comps(icomp2)%el(ie)%mag   = idou(ie)
          comps(icomp2)%el(ie)%pres   = pres(ie)
          comps(icomp2)%el(ie)%dforce = dF(:,ie)
        enddo
        deallocate(idou, pres, dF)
        exit
      endif

    enddo

    !If the component was not found there is a problem
    if(icomp2 .gt. size(comps)) call error(this_sub_name, this_mod_name, &
    'Component loaded from geometry file not found in result file')

    call close_hdf5_group(gloc2)
  enddo

  call close_hdf5_group(gloc1)
  call close_hdf5_file(floc)

end subroutine load_solution

!----------------------------------------------------------------------

!> Chech that the component that we are about to load will be placed, in
!! the actual reference frames, in the same position in which it was saved
!!
!! This is done comparing offsets and rotations between the actual references
!! and the one saved
subroutine check_ref(gloc, floc, ref)
 integer(h5loc), intent(in) :: gloc
 integer(h5loc), intent(in) :: floc
 type(t_ref), intent(in)    :: ref

 integer :: iref
 integer(h5loc) :: refs_gloc, ref_loc
 character(len=max_char_len) :: ref_title
 real(wp) :: R(3,3), of(3)
 character(len=*), parameter :: this_sub_name='check_ref'
 
  call read_hdf5(iref, 'RefId', gloc)
  call open_hdf5_group(floc, 'References', refs_gloc)
  write(ref_title,'(A,I3.3)')'Ref',iref
  call open_hdf5_group(refs_gloc, trim(ref_title), ref_loc)

  call read_hdf5(R,'R',ref_loc)
  call read_hdf5(of,'Offset',ref_loc)

! ! debug ----
! write(*,*) ' offset : read vs computed '
! write(*,*) real(of) , '     ' , real(ref%of_g)
! ! debug ----

  if ((.not. all(real(R) .eq. real(ref%R_g))) .or. &
      (.not. all(real(of) .eq. real(ref%of_g))) ) then
    call warning(this_sub_name, this_mod_name, 'Mismatching initial position &
    & while loading reference '//trim(ref%tag)//' This may lead to &
    &unexpected behaviour')
  endif
  
  
  call close_hdf5_group(ref_loc)
  call close_hdf5_group(refs_gloc)


end subroutine check_ref

!----------------------------------------------------------------------

!> Load the time value from a result file
subroutine load_time(filename, time)
 character(len=*), intent(in) :: filename
 real(wp), intent(out) :: time

 integer(h5loc) :: floc

  call open_hdf5_file(filename, floc)

  call read_hdf5(time,'time',floc)

  call close_hdf5_file(floc)
end subroutine load_time

!----------------------------------------------------------------------

! moved to geo/mod_references.f90

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
! 
!  integer :: iref , nref , idref
!  integer :: i , it , nref_ref_in
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
!  nref_ref_in = countoption(ref_prs,'Reference_Tag') 
! 
!  ! open restart file to read the relative initial conditions
!  call open_hdf5_file(trim(restart_file), floc)
!  call open_hdf5_group(floc, 'References', refs_gloc)
! 
!  iref = 0
!  write(*,*) ' @@@@@@@@@@@@@@ '
!  do i = 1 , nref_ref_in
!     
!    iref = iref + 1
!    write(*,*) ' iref : ' , iref
! 
!    write(ref_title,'(A,I3.3)')'Ref',iref
!    call open_hdf5_group(refs_gloc, trim(ref_title), ref_loc)
! 
!    if ( refs(iref)%self_moving ) then
! 
! 
!      call getsuboption(ref_prs,'Motion',sbprms)
! 
!      ! === Pole ===
!      call getsuboption(sbprms,'Pole',sbprms_pol)
! 
!      select case( trim(getstr(sbprms_pol,'Input')) )
!         case('velocity')
!           ! add the i.c. from restart file
!           write(*,*) ' pole velocity, iref : ' , iref
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
!           ! add the i.c. from restart file
!           write(*,*) ' rot  velocity, iref : ' , iref
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
!    end if
! 
!    call close_hdf5_group(ref_loc)
! 
! ! Multiplicity -----
! !  if ( ) then ! loop over multiplicity
! !
! !    select case(trim(   mult_type))
! !
! !      case('simple rotor')
! !
! !      case('rotor')
! !
! !
! !    end select
! !
! !  end if
! 
! 
!  end do
! 
!  call close_hdf5_group(refs_gloc)
!  call close_hdf5_file(floc)
! 
! end subroutine update_relative_initial_conditions
! 
! !----------------------------------------------------------------------

end module mod_dust_io
