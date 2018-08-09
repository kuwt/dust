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

module mod_dat_out

use mod_param, only: &
  wp, nl, max_char_len, extended_char_len, ascii_real

use mod_handling, only: &
  error, warning ! , info, printout, dust_time, t_realtime

!---------------------------------------------------------------------
implicit none

public :: dat_out_probes_header , dat_out_loads_header , &
          dat_out_sectional

private

character(len=*), parameter :: &
  this_mod_name = 'mod_dat_output'
!---------------------------------------------------------------------

contains

!---------------------------------------------------------------------

subroutine dat_out_loads_header ( fid , comps_meas , ref_sys , average)
 integer , intent(in) :: fid
 character(len=*), intent(in) :: comps_meas(:)
 character(len=*), intent(in) :: ref_sys 
 logical, intent(in)          :: average

 character(len=max_char_len) :: istr
 integer :: n_comps , ic

 n_comps = size(comps_meas)

 write(istr,'(I0)') n_comps 
 write(fid,*) '# comments ... n.components: ' , trim(istr)
 write(fid,*) '# comments ... ref.sys     : ' , trim(ref_sys)
 ! three-dimensional space
 do ic = 1 , n_comps - 1
   write(fid,'(A)',advance='no') trim(comps_meas(ic))//' , '
 end do
 write(fid,'(A)') trim(comps_meas(n_comps))
 if (.not. average)  then
   write(fid,*) '#  t  Fx  Fy  Fz  Mx  My  Mz  ' 
 else
   write(fid,*) '# Fx_average  Fy_average  Fz_average &
                    & Mx_average  My_average  Mz_average ' 
 endif
   

end subroutine dat_out_loads_header

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
 
  write(fid,*) '# comments ...' , n_probes
  ! three-dimensional space
  do ic = 1 , 3
    write(nnum,'(I0)') size(rr_probes,2)
    write(fid,'('//trim(nnum)//ascii_real//')') rr_probes(ic,:)
  end do
 
  write(istr,'(I0)') n_probes 
  write(fid,*) '#    t     '//trim(istr)//'( '//trim(vars_str)//' )'

end subroutine dat_out_probes_header

!---------------------------------------------------------------------

subroutine dat_out_sectional (basename, y_cen ,y_span ,time ,sec_loads , &
              ref_mat , off_mat, average )
 character(len=*) , intent(in) :: basename
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
    write(*,*) ' size(sec_loads,2) : ' , size(sec_loads,2)
    write(*,*) ' size(y_cen)       : ' , size(y_cen)
    call error(trim(this_mod_name),'','Inconsistent inputs.& 
            & size(sec_loads,2) .ne. size(y_cen). Stop ')
  end if
  if ( size(sec_loads,1) .ne. nt ) then
    call error(trim(this_mod_name),'','Inconsistent inputs.& 
            & size(sec_loads,1) .ne. size(time). Stop ')
  end if
  if ( size(sec_loads,3) .ne. 4 ) then
    call error(trim(this_mod_name),'','Inconsistent inputs.& 
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
    write(fid,'(A)') '# TODO: component name and other beautiful stuff '
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

end module mod_dat_out
