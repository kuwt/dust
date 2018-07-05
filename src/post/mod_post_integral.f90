module mod_post_integral

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
  LowCase, isInList, stricmp

use mod_geo_postpro, only: &
  load_components_postpro, update_points_postpro , prepare_geometry_postpro , &
  prepare_wake_postpro  ! expand_actdisk_postpro, 

use mod_geometry, only: &
  t_geo, t_geo_component

use mod_post_load, only: &
  load_refs, load_res ! , load_wake_viz , load_wake_pan , load_wake_ring

use mod_tecplot_out, only: &
  tec_out_loads

use mod_dat_out, only: & 
  dat_out_loads_header

use mod_math, only: &
  cross

implicit none

public :: post_integral

private

character(len=max_char_len), parameter :: this_mod_name = 'mod_post_integral'

contains

! ---------------------------------------------------------------------- 

subroutine post_integral( sbprms, basename, data_basename, an_name , ia , &
                          out_frmt, comps , components_names, all_comp , &
                          an_start, an_end, an_step )
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

integer(h5loc) :: floc , ploc
real(wp), allocatable :: points(:,:)
integer , allocatable :: elems(:,:)
integer :: nelem

character(len=max_char_len) :: filename
character(len=max_char_len) , allocatable :: refs_tag(:)
real(wp), allocatable :: refs_R(:,:,:), refs_off(:,:)
real(wp), allocatable :: refs_G(:,:,:), refs_f(:,:)
real(wp), allocatable :: vort(:), cp(:)
character(len=max_char_len) :: ref_tag
integer                     :: ref_id
real(wp) :: F_loc(3) , F_ref(3) , F_bas(3) , F_bas1(3)
real(wp) :: M_loc(3) , M_ref(3) , M_bas(3)
real(wp), allocatable :: force(:,:), moment(:,:)
real(wp) :: u_inf(3)
real(wp) :: P_inf , rho
integer :: ic2 , ic , it , ie , ierr , ires , fid_out , nstep
real(wp), allocatable :: time(:)
real(wp) :: t


character(len=max_char_len), parameter :: this_sub_name = 'post_integral'

    write(*,*) nl//' Analysis:',ia,' post_integral() ++++++++++ '//nl

!DEBUG
write(*,*) ' before load_components_postpro '
do ic = 1 , size(components_names)
  write(*,*) trim(components_names(ic))
end do

! load the geo components just once just once
call open_hdf5_file(trim(data_basename)//'_geo.h5', floc)
!TODO: here get the run id
call load_components_postpro(comps, points, nelem, floc, & 
                             components_names,  all_comp)
call close_hdf5_file(floc)

!DEBUG
write(*,*) '  after load_components_postpro '
do ic = 1 , size(components_names)
  write(*,*) trim(components_names(ic))
end do

! Prepare_geometry_postpro
call prepare_geometry_postpro(comps)

! Reference system where the loads are projected:
!  read the input and check if it exist
ref_tag = getstr(sbprms,'Reference_Tag')

write(filename,'(A,I4.4,A)') trim(data_basename)//'_res_',an_start,'.h5'
call open_hdf5_file(trim(filename),floc)
call load_refs(floc,refs_R,refs_off,refs_G,refs_f,refs_tag)
call close_hdf5_file(floc) 

ref_id = -333  ! check input
do it = lbound(refs_tag,1) , ubound(refs_tag,1)
  if ( stricmp(refs_tag(it),  ref_tag) ) ref_id = it
end do
if ( ref_id .eq. -333 ) then 
  write(*,*)
  write(*,*) ' Available references systems: '
  do it = lbound(refs_tag,1) , ubound(refs_tag,1)
    write(*,*) ' ref_id : ' , it , ' ref_tag ' , trim(refs_tag(it))
  end do
  call warning('dust_post','','Unknown ref.sys. defined for loads output.&
       & Your input in dust_post.in is '//trim(ref_tag)//'. All the&
       & available ref.sys. are listed above.')
  return ! jump this analysis if the reference frame is not available
end if

nstep = (an_end-an_start)/an_step + 1
! Output format
select case(trim(out_frmt))

 case('dat')
  ! Open output .dat file
  call new_file_unit(fid_out, ierr)
  write(filename,'(A)') trim(basename)//'_'//trim(an_name)//'.dat'
  open(unit=fid_out,file=trim(filename))
  call dat_out_loads_header( fid_out , components_names , ref_tag )

 case('tecplot')
  allocate(force(3,nstep), moment(3,nstep), time(nstep))

 case default
  call error('dust_post','','Unknown format '//trim(out_frmt)//&
             ' for integral_loads analysis. Choose dat or tecplot.')

end select

ires = 0
do it=an_start, an_end, an_step ! Time loop
  ires = ires+1

  ! Open the file:
  write(filename,'(A,I4.4,A)') trim(data_basename)//'_res_',it,'.h5'
  call open_hdf5_file(trim(filename),floc)

  ! Load u_inf --------------------------------
  call open_hdf5_group(floc,'Parameters',ploc)
  call read_hdf5(u_inf,'u_inf',ploc)
  call read_hdf5(P_inf,'P_inf',ploc)
  call read_hdf5(rho,'rho_inf',ploc)
  call close_hdf5_group(ploc)

  ! Load the references and move the points ---
  call load_refs(floc,refs_R,refs_off)
  ! Move the points ---------------------------
  call update_points_postpro(comps, points, refs_R, refs_off)
  ! Load the results --------------------------
  call load_res(floc, comps, vort, cp, t)

  call close_hdf5_file(floc)


  ! Initialise integral loads in the desired ref.frame
  F_ref = 0.0_wp ; M_ref = 0.0_wp 

  ! Update the overall load with the comtribution from all the components
  do ic = 1 , size(comps)

    ! Initialise integral loads in the local ref.frame
    F_bas = 0.0_wp ; M_bas = 0.0_wp 
  
    ! Loads from the ic-th component in the base ref.frame
    do ie = 1 , size(comps(ic)%el)
      F_bas1 = comps(ic)%el(ie)%dforce

      F_bas = F_bas + F_bas1

      M_bas = M_bas + cross( comps(ic)%el(ie)%cen &
                     -refs_off(:,ref_id) , F_bas1 )

    end do !ie

    ! From the base ref.sys to the chosen ref.sys (offset and rotation)
    F_ref = F_ref + matmul( &
         transpose( refs_R(:,:, ref_id) ) , F_bas )
    M_ref = M_ref + matmul( &
         transpose( refs_R(:,:, ref_id) ) , M_bas )

  end do !ic
  
  ! Update output dat file / update output arrays for tecplot
  select case(trim(out_frmt))

   case ('dat')
    write(fid_out,'(E12.3)'  ,advance='no') t 
    write(fid_out,'(3E12.3)' ,advance='no') F_ref
    write(fid_out,'(3E12.3)' ,advance='no') M_ref
    write(fid_out,'(9E12.3)',advance='no') refs_R(:,:, ref_id)
    write(fid_out,'(3E12.3)',advance='no') refs_off(:, ref_id)
    write(fid_out,*) ' '

   case('tecplot')
    time(ires) = t
    force(:,ires) = F_ref
    moment(:,ires) = M_ref

  end select

end do ! Time loop

! Close dat file / write tec file
select case(trim(out_frmt))

 case('dat')
  close(fid_out)
 
 case('tecplot')
  write(filename,'(A)') trim(basename)//'_'//trim(an_name)//'.plt'
  call tec_out_loads(filename, time, force, moment)
  deallocate(time, force, moment)

end select


!TODO: move deallocate(comps) outside this routine,
!      because it is common to all the analyses
deallocate(comps,components_names)

    write(*,*) nl//' post_integral done.'//nl

end subroutine post_integral

! ---------------------------------------------------------------------- 

end module mod_post_integral
