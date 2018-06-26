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

use mod_geometry, only: &
  t_geo, t_geo_component

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

integer(h5loc) :: floc ! , ploc
real(wp), allocatable :: points(:,:)
integer , allocatable :: elems(:,:)
integer :: nelem

integer :: ic

character(len=max_char_len), parameter :: this_sub_name = 'post_integral'

    write(*,*) nl//' Analysis:',ia,' post_viz() +++++++++++++++ '//nl

! load the geo components just once just once
call open_hdf5_file(trim(data_basename)//'_geo.h5', floc)
!TODO: here get the run id
call load_components_postpro(comps, points, nelem, floc, & 
                             components_names,  all_comp)
call close_hdf5_file(floc)

! Prepare_geometry_postpro
call prepare_geometry_postpro(comps)

!! ! Select all the components
!! if ( allocated(components_names) ) then
!!   call warning(trim(this_sub_name), trim(this_mod_name), &
!!      'All the components are used. <Components> input &
!!      &is ignored, and deallocated.' )
!!   deallocate(components_names)
!! end if 
!! allocate(components_names(size(comps)))
!! do ic = 1,size(comps)
!!   components_names(ic) = trim(comps(ic)%comp_name)
!!   write(*,*) ' components_names(',ic,') : ' , &
!!     trim(components_names(ic))
!! enddo
!! 
!! 
!! 
!TODO: move deallocate(comps) outside this routine,
!      because it is common to all the analyses
deallocate(comps)


end subroutine post_integral

! ---------------------------------------------------------------------- 

end module mod_post_integral
