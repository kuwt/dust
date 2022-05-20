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
!! Copyright (C) 2018-2022 Politecnico di Milano,
!!                           with support from A^3 from Airbus
!!                    and  Davide   Montagnani,
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
!!          Federico Fonte
!!          Davide Montagnani
!!          Matteo Tugnoli
!!=========================================================================

module mod_dat_out

use mod_param, only: &
  wp, nl, max_char_len, extended_char_len, ascii_real

use mod_handling, only: &
  error, warning, internal_error, new_file_unit

!---------------------------------------------------------------------
implicit none

public :: dat_out_probes_header, dat_out_loads_header, dat_out_hinge_header,  dat_out_aa_header, &
          dat_out_sectional, dat_out_sectional_ll, dat_out_sectional_vl, dat_out_aa

private

character(len=*), parameter :: &
  this_mod_name = 'mod_dat_output'
!---------------------------------------------------------------------

contains

!---------------------------------------------------------------------

subroutine dat_out_loads_header ( fid , comps_meas , ref_sys , average)
  integer , intent(in)          :: fid
  character(len=*), intent(in)  :: comps_meas(:)
  character(len=*), intent(in)  :: ref_sys
  logical, intent(in)           :: average

  character(len=max_char_len)   :: istr
  integer :: n_comps , ic

  n_comps = size(comps_meas)

  write(istr,'(I0)') n_comps
  write(fid,*) '# Integral loads: N.components: ' , trim(istr)
  write(fid,*) '#                 Ref.sys     : ' , trim(ref_sys)
  
  !> three-dimensional space
  write(fid,'(A)',advance='no') ' #                 Components  : '
  
  do ic = 1 , n_comps - 1
    write(fid,'(A)',advance='no') trim(comps_meas(ic))//' , '
  end do

  write(fid,'(A)') trim(comps_meas(n_comps))
  
  if (.not. average)  then
    write(fid,*) '#  t , Fx , Fy , Fz , Mx , My , Mz , ref_mat(9) , ref_off(3) '
  else
    write(fid,*) '# Fx_average , Fy_average , Fz_average ,&
                  & Mx_average , My_average , Mz_average ,&
                  & ref_mat(9) , ref_off(3) '
  endif

end subroutine dat_out_loads_header

subroutine dat_out_hinge_header ( fid , comps_meas , hinge_tag , average)
  integer , intent(in)          :: fid
  character(len=*), intent(in)  :: comps_meas
  character(len=*), intent(in)  :: hinge_tag
  logical, intent(in)           :: average

  write(fid,*) '# Hinge Moment: '
  !> three-dimensional space
  write(fid,'(A)',advance='no') ' #               Components  : '
  write(fid,'(A)') trim(comps_meas)

  !> three-dimensional space
  write(fid,'(A)',advance='no') ' #                    Hinge  : '
  write(fid,'(A)') trim(hinge_tag)

  if (.not. average)  then
    write(fid,*) '#  t , Fv , Fh , Fn , Mv , Mh , Mn , axis_mat(9) , node_hinge(3) '
  else
    write(fid,*) '# Fv_average , Fh_average , Fn_average ,&
                  & Mv_average , Mh_average , Mn_average ,&
                  & axis_mat(9) , node_hinge(3) '
  endif

end subroutine dat_out_hinge_header

!---------------------------------------------------------------------

subroutine dat_out_aa_header ( fid , t)
  integer , intent(in) :: fid
  real(wp), intent(in) :: t

  write(fid,'(A)') '# Aeroacoustic Data'

  write(fid,'(A)') '# Time'
  write(fid,'('//ascii_real//')') t

  write(fid,'(A)') '# Element Center (3), Element Normal (3), &
                &Element Center Velocity (3), Element Area (1), &
                &Force Acting on Element (3)'

end subroutine dat_out_aa_header

!---------------------------------------------------------------------

subroutine dat_out_aa ( fid , cen, n, vel, area, f)
  integer , intent(in) :: fid
  real(wp), intent(in) :: cen(3), n(3), vel(3), area, f(3)

  write(fid,'(13'//ascii_real//')') cen, n, vel, area, f

end subroutine dat_out_aa

!---------------------------------------------------------------------

subroutine dat_out_probes_header ( fid , rr_probes , vars_str )
  integer , intent(in) :: fid
  real(wp), intent(in) :: rr_probes(:,:)
  character(len=*), intent(in) :: vars_str

  character(len=max_char_len) :: istr
  integer :: n_probes , ic
  character(len=8) :: nnum

  n_probes = size(rr_probes,2)

  if ( size(rr_probes,1) .ne. 3 ) then
    call error(trim(this_mod_name),'','Wrong format of the rr_probes inputs.&
            & size(rr_probes,1) .ne. 3. Stop ')
  end if

  write(fid,*) '# N. of point probes:' , n_probes
  !> three-dimensional space
  do ic = 1 , 3
    write(nnum,'(I0)') size(rr_probes,2)
    write(fid,'('//trim(nnum)//ascii_real//')') rr_probes(ic,:)
  end do

  write(istr,'(I0)') n_probes
  write(fid,*) '#    t     '//trim(istr)//'( '//trim(vars_str)//' )'

end subroutine dat_out_probes_header

!---------------------------------------------------------------------

subroutine dat_out_sectional (basename, compname, y_cen, y_span, time, &
                              sec_loads, ref_mat, off_mat, average )
  character(len=*) , intent(in) :: basename
  character(len=*) , intent(in) :: compname
  real(wp) , intent(in) :: y_cen(:)
  real(wp) , intent(in) :: y_span(:)
  real(wp) , intent(in) :: time(:)
  real(wp) , intent(in) :: sec_loads(:,:,:)
  real(wp) , intent(in) :: ref_mat(:,:)
  real(wp) , intent(in) :: off_mat(:,:)
  logical,   intent(in) :: average

  character(len=2) :: load_str(4)
  character(len=8) :: nnum
  character(len=max_char_len) :: filename
  integer :: it , nt , fid , i1

  load_str = (/ 'Fx' , 'Fy' , 'Fz' , 'Mo' /)

  nt = size(time)

  ! Some checks --------
  if ( size(y_cen) .ne. size(sec_loads,2) ) then
    call internal_error(trim(this_mod_name),'','Inconsistent inputs.&
            & size(sec_loads,2) .ne. size(y_cen). Stop ')
  end if
  
  if ( size(sec_loads,1) .ne. nt ) then
    call internal_error(trim(this_mod_name),'','Inconsistent inputs.&
            & size(sec_loads,1) .ne. size(time). Stop ')
  end if
  
  if ( size(sec_loads,3) .ne. 4 ) then
    call internal_error(trim(this_mod_name),'','Inconsistent inputs.&
            & size(sec_loads,3) .ne. 4. Stop ')
  end if

  ! Print out .dat files
  fid = 21
  do i1 = 1 , 4
    if(average) then
      write(filename,'(A)') trim(basename)//'_'//trim(load_str(i1))//'_ave.dat'
    else
      write(filename,'(A)') trim(basename)//'_'//trim(load_str(i1))//'.dat'
    endif

    
    open(unit=fid,file=trim(filename))
    ! Header -----------
    write(fid,'(A)') '# Sectional load '//trim(load_str(i1))//&
                    &' of component: '//trim(compname)
    write(fid,'(A,I0,A,I0,A)') '# n_sec : ' , size(sec_loads,2) , ' ; n_time : ' , nt , '. Next lines: y_cen , y_span'
    write(nnum,'(I0)') size(y_cen)
    write(fid,'('//trim(nnum)//ascii_real//')') y_cen
    write(fid,'('//trim(nnum)//ascii_real//')') y_span

    if(average) then
      write(fid,'(A)') '#sec(n_sec)'
      write(nnum,'(I0)') size(y_cen)
      write(fid,'('//trim(nnum)//ascii_real//')') sec_loads(1,:,i1)
    else
      write(fid,'(A)') '# t , sec(n_sec) , ref_mat(9) , ref_off(3) '
      ! Dump data --------
      do it = 1 , nt
        write(nnum,'(I0)') 1+size(y_cen)+9+3
        write(fid,'('//trim(nnum)//ascii_real//')') time(it), &
                          sec_loads(it,:,i1) , ref_mat(it,:) , off_mat(it,:)
      end do
    endif
    close(fid)

  end do


end subroutine dat_out_sectional

!---------------------------------------------------------------------

subroutine dat_out_sectional_ll (basename, compname, y_cen, y_span, time, &
                                ll_sec, average )
  character(len=*) , intent(in) :: basename
  character(len=*) , intent(in) :: compname
  real(wp) , intent(in)         :: y_cen(:)
  real(wp) , intent(in)         :: y_span(:)
  real(wp) , intent(in)         :: time(:)
  real(wp) , intent(in)         :: ll_sec(:,:,:)
  logical,   intent(in)         :: average

  character(len=8)              :: nnum
  character(len=max_char_len)   :: filename
  integer                       :: it , nt , fid , ierr, il
  character(len=*), parameter   :: this_sub_name = 'dat_out_sectional_ll'
  character(len=21)             :: load_str(9)
  character(len=max_char_len)   :: description_str(9)

  load_str(1) = 'Cl' 
  load_str(2) = 'Cd' 
  load_str(3) = 'Cm'
  load_str(4) = 'alpha'
  load_str(5) = 'alpha_isolated'
  load_str(6) = 'vel_2d'
  load_str(7) = 'vel_2d_isolated'
  load_str(8) = 'vel_outplane'
  load_str(9) = 'vel_outplane_isolated'
  
  description_str(1) = 'lift coefficient'
  description_str(2) = 'drag coefficient'
  description_str(3) = 'moment coefficient'
  description_str(4) = 'angle of attack'
  description_str(5) = 'isolated angle of attack'
  description_str(6) = 'in section plane velocity'
  description_str(7) = 'in section plane velocity isolated'
  description_str(8) = 'out of section (spanwise) velocity'
  description_str(9) = 'out of section (spanwise) isolated velocity'

  nt = size(time)

  ! Some checks --------
  if ( size(y_cen) .ne. size(ll_sec,2) ) then
    call internal_error(trim(this_mod_name),'','Inconsistent inputs.&
            & different length of nodes and solution ')
  end if
  if ( size(ll_sec,1) .ne. nt ) then
    call internal_error(trim(this_mod_name),'','Inconsistent inputs.&
            & size(sec_loads,1) .ne. size(time). Stop ')
  end if

  !> Print out .dat files
  do il = 1 , size(load_str)
    call new_file_unit(fid, ierr)
    if(average) then
      write(filename,'(A)') trim(basename)//'_'//trim(load_str(il))//'_ave.dat'
    else
      write(filename,'(A)') trim(basename)//'_'//trim(load_str(il))//'.dat'
    endif
    open(unit=fid,file=trim(filename))
    !> Header 
    write(fid,'(A)') '# Sectional '//trim(description_str(il))//&
                    &' of component: '//trim(compname)
    write(fid,'(A,I0,A,I0,A)') '# n_sec : ' , size(y_cen) , ' ; n_time : ' , nt , '. Next lines: y_cen , y_span'
    write(nnum,'(I0)') size(y_cen)
    write(fid,'('//trim(nnum)//ascii_real//')') y_cen
    write(fid,'('//trim(nnum)//ascii_real//')') y_span

    if(average) then
      write(fid,'(A)') '# '//trim(load_str(il))//'(n_sec)'
      write(nnum,'(I0)') size(y_cen)
      write(fid,'('//trim(nnum)//ascii_real//')') ll_sec(1,:,il)
    else
      write(fid,'(A)') '# t , '//trim(load_str(il))//'(n_sec) '
      do it = 1 , nt
        write(nnum,'(I0)') 1+size(y_cen)
        write(fid,'('//trim(nnum)//ascii_real//')') time(it), ll_sec(it,:,il)
      end do
    endif
    close(fid)
  enddo


end subroutine dat_out_sectional_ll


subroutine dat_out_sectional_vl (basename, compname, y_cen, y_span, time, &
  vl_sec, average )

  character(len=*) , intent(in) :: basename
  character(len=*) , intent(in) :: compname
  real(wp) , intent(in)         :: y_cen(:)
  real(wp) , intent(in)         :: y_span(:)
  real(wp) , intent(in)         :: time(:)
  real(wp) , intent(in)         :: vl_sec(:,:,:)
  logical,   intent(in)         :: average

  character(len=8)              :: nnum
  character(len=max_char_len)   :: filename
  integer                       :: it , nt , fid , ierr, il
  character(len=*), parameter   :: this_sub_name = 'dat_out_sectional_vl'
  character(len=21)             :: load_str(9)
  character(len=max_char_len)   :: description_str(9)

  load_str(1) = 'Cl'; 
  load_str(2) = 'Cd'; 
  load_str(3) = 'Cm'
  load_str(4) = 'alpha'; 
  load_str(5) = 'alpha_isolated'
  load_str(6) = 'vel_2d'; 
  load_str(7) = 'vel_2d_isolated'
  load_str(8) = 'vel_outplane'; 
  load_str(9) = 'vel_outplane_isolated'

  description_str(1) = 'lift coefficient'
  description_str(2) = 'drag coefficient'
  description_str(3) = 'moment coefficient'
  description_str(4) = 'angle of attack'
  description_str(5) = 'isolated angle of attack'
  description_str(6) = 'in section plane velocity'
  description_str(7) = 'in section plane velocity isolated'
  description_str(8) = 'out of section (spanwise) velocity'
  description_str(9) = 'out of section (spanwise) isolated velocity'

  nt = size(time)

  ! Some checks --------
  if ( size(y_cen) .ne. size(vl_sec,2) ) then
    call internal_error(trim(this_mod_name),'','Inconsistent inputs.&
    & different length of nodes and solution ')
  end if
  if ( size(vl_sec,1) .ne. nt ) then
    call internal_error(trim(this_mod_name),'','Inconsistent inputs.&
    & size(sec_loads,1) .ne. size(time). Stop ')
  end if

  ! Print out .dat files
  do il = 1 , size(load_str)
    call new_file_unit(fid, ierr)
    if(average) then
      write(filename,'(A)') trim(basename)//'_'//trim(load_str(il))//'_ave.dat'
    else
      write(filename,'(A)') trim(basename)//'_'//trim(load_str(il))//'.dat'
    endif
    open(unit=fid,file=trim(filename))
    ! Header -----------
    write(fid,'(A)') '# Sectional '//trim(description_str(il))//&
    &' of component: '//trim(compname)
    write(fid,'(A,I0,A,I0,A)') '# n_sec : ' , size(y_cen) , ' ; n_time : ' , nt , '. Next lines: y_cen , y_span'
    write(nnum,'(I0)') size(y_cen)
    write(fid,'('//trim(nnum)//ascii_real//')') y_cen
    write(fid,'('//trim(nnum)//ascii_real//')') y_span

    if(average) then
      write(fid,'(A)') '# '//trim(load_str(il))//'(n_sec)'
      write(nnum,'(I0)') size(y_cen)
      write(fid,'('//trim(nnum)//ascii_real//')') vl_sec(1,:,il)
    else
      write(fid,'(A)') '# t , '//trim(load_str(il))//'(n_sec) '
      do it = 1 , nt
        write(nnum,'(I0)') 1+size(y_cen)
        write(fid,'('//trim(nnum)//ascii_real//')') time(it), vl_sec(it,:,il)
      end do
    endif
    close(fid)
  enddo

end subroutine dat_out_sectional_vl

!---------------------------------------------------------------------


end module mod_dat_out
