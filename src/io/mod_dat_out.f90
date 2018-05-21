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
  wp, nl, max_char_len, extended_char_len

use mod_handling, only: &
  error, warning ! , info, printout, dust_time, t_realtime

!---------------------------------------------------------------------
implicit none

public :: dat_out_probes_header , dat_out_probes 

private

character(len=*), parameter :: &
  this_mod_name = 'mod_dat_output'
!---------------------------------------------------------------------

contains

!---------------------------------------------------------------------

subroutine dat_out_probes_header ( fid , rr_probes , vars_str )
 integer , intent(in) :: fid
 real(wp), intent(in) :: rr_probes(:,:)
 character(len=*), intent(in) :: vars_str

 character(len=max_char_len) :: istr
 integer :: n_probes , ic

 n_probes = size(rr_probes,2)

 if ( size(rr_probes,1) .ne. 3 ) then
   call error(trim(this_mod_name),'','Wrong format of the rr_probes inputs.& 
           & size(rr_probes,1) .ne. 3. Stop ')
 end if

 write(fid,*) '# comments ...' , n_probes
 ! three-dimensional space
 do ic = 1 , 3
   write(fid,*) rr_probes(ic,:) 
 end do

 write(istr,'(I0)') n_probes 
 write(fid,*) '#    t     '//istr//'( '//trim(vars_str)//')'

end subroutine dat_out_probes_header

!---------------------------------------------------------------------

subroutine dat_out_probes ( fid , vars )
 integer , intent(in) :: fid
 real(wp), intent(in) :: vars(:)

 

end subroutine dat_out_probes

!---------------------------------------------------------------------

end module mod_dat_out
