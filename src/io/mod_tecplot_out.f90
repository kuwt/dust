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

module mod_tecplot_out
 
use mod_param, only: &
  wp, max_char_len, nl

use mod_handling, only: &
  error, warning, info, printout, new_file_unit

!---------------------------------------------------------------------
implicit none

public :: tec_out_sol_bin, tec_out_viz, tec_out_probes, tec_out_box

private

integer, parameter :: s_size = 4
integer, parameter :: d_size = 8
character(len=*), parameter :: &
  this_mod_name = 'mod_tecplot_output'

!---------------------------------------------------------------------

contains

!---------------------------------------------------------------------

subroutine tec_out_sol_bin(rr, ee, vort, w_rr, w_ee, w_vort, t, out_filename)
 real(wp), intent(in) :: rr(:,:)
 integer, intent(in)  :: ee(:,:)
 real(wp), intent(in) :: vort(:)
 real(wp), intent(in) :: w_rr(:,:)
 integer, intent(in)  :: w_ee(:,:)
 real(wp), intent(in) :: w_vort(:)
 real(wp), intent(in) :: t
 character(len=*), intent(in) :: out_filename 

 character, parameter :: zc = char(0)
 integer :: npoints, ncells, ne
 integer :: fu, ierr, i1, i2
 real(kind=4)  , parameter :: zoneMarker = 299.0
 real(kind=4)  , parameter :: eohMarker  = 357.0
 character(len=max_char_len) :: buffer_char
 integer :: ie, nquad, ntria
 integer npoints_w, nw



  ne = size(ee,2) !number of elements
  nw = size(w_ee,2) !number of wake elements
  npoints = size(rr,2)
  npoints_w = size(w_rr,2)
  ncells = ne

  ! First cycle the elements to get the number of quads and trias
  ! (this is really ugly, make it more flexible)
  nquad = 0; ntria = 0
  do ie = 1,ne
    if(ee(4,ie) .eq. 0) then
      ntria = ntria+1
    else
      nquad = nquad+1
    endif
  enddo

  call new_file_unit(fu,ierr)
  open(unit=fu,file=trim(out_filename),status='replace',access='stream', &
       form='unformatted',iostat=ierr)

  !magic number
  buffer_char = "#!TDV112" ;  write(fu) trim(buffer_char) 

  !integer 1 to set the byte order
  write(fu) int(1,s_size)

  !integer 0 to mark a full file (1=solution only, 2=grid only)
  write(fu) int(0,s_size)

  !title
  call put_tec_string('DUST solution',fu)

  !number of variables
  write(fu) int(4,s_size) !x-y-z-vort

  !Variables names
  call put_tec_string('X',fu)
  call put_tec_string('Y',fu)
  call put_tec_string('Z',fu)
  call put_tec_string('Vort',fu)
 
  !BODY
  !Zone marker
  write(fu) zoneMarker

  !Zone name
  call put_tec_string('Body',fu)

  !Parent zone: -1 for no parent (not too clear)
  write(fu) int(-1,s_size)

  !Strand id: -2 for automatic generation
  write(fu) int(-2,s_size)

  !Solution time
  write(fu) real(t,d_size)

  !"not used, set to -1"
  write(fu) int(-1,s_size) 
  
  !Zone type: 0 ordered, 3 unstructured quadrilateral
  write(fu) int(3,s_size)

  !Toggle variable location specification
  write(fu) int(1,s_size)

  !x, y, z stand at nodes
  write(fu) int(0,s_size)
  write(fu) int(0,s_size)
  write(fu) int(0,s_size)
  ! Vort stands at cells centers
  write(fu) int(1,s_size)

  !are local neighbour supplied? No.
  write(fu) int(0,s_size)

  !Number of miscellaneous somethings
  write(fu) int(0,s_size)

  ! n.points, n.elems
  write(fu) int(size(rr,2),s_size)
  write(fu) int(size(ee,2),s_size)
  
  ! All zero, three fields, for future use
  write(fu) int(0,s_size) , int(0,s_size) , int(0,s_size)
  
  !"no more auxilliary name/value pairs"
  write(fu) int(0,s_size)


  !WAKE ZONE HERE

  !Zone marker
  write(fu) zoneMarker

  !Zone name
  call put_tec_string('Wake',fu)

  !Parent zone: -1 for no parent (not too clear)
  write(fu) int(-1,s_size)

  !Strand id: -2 for automatic generation
  write(fu) int(-2,s_size)

  !Solution time
  write(fu) real(t,d_size)

  !"not used, set to -1"
  write(fu) int(-1,s_size) 
  
  !Zone type: 0 ordered, 3 unstructured quadrilateral
  write(fu) int(3,s_size)

  !Toggle variable location specification
  write(fu) int(1,s_size)

  !x, y, z stand at nodes
  write(fu) int(0,s_size)
  write(fu) int(0,s_size)
  write(fu) int(0,s_size)
  ! Vort stands at cells centers
  write(fu) int(1,s_size)

  !are local neighbour supplied? No.
  write(fu) int(0,s_size)

  !Number of miscellaneous somethings
  write(fu) int(0,s_size)

  ! n.points, n.elems
  write(fu) int(size(w_rr,2),s_size)
  write(fu) int(size(w_ee,2),s_size)
  
  ! All zero, three fields, for future use
  write(fu) int(0,s_size) , int(0,s_size) , int(0,s_size)
  
  !"no more auxilliary name/value pairs"
  write(fu) int(0,s_size)


  !DATA
  write(fu) eohMarker
  
  !BODY
  !Zone marker
  write(fu) zoneMarker

  !Size of the variables: 2=double
  write(fu) int(2,s_size)
  write(fu) int(2,s_size)
  write(fu) int(2,s_size)
  write(fu) int(2,s_size)

  !Has passive variables?
  write(fu) int(0,s_size)

  !If has passive vars here must be the list of the passive vars

  !Has variable sharing
  write(fu) int(0,s_size)

  !sharing connectivity? -1 no sharing
  write(fu) int(-1,s_size)
  
  !Min and max of all the variables
  !x
  write(fu) dble(minval( rr(1,:) ))
  write(fu) dble(maxval( rr(1,:) ))
  !y
  write(fu) dble(minval( rr(2,:) ))
  write(fu) dble(maxval( rr(2,:) ))
  !z
  write(fu) dble(minval( rr(3,:) ))
  write(fu) dble(maxval( rr(3,:) ))
  !Vort
  write(fu) dble(minval( vort ))
  write(fu) dble(maxval( vort ))

  !The actual data
  ! x, y, z
  do i1 = 1 , 3
    do i2 = 1 , size(rr,2)
      write(fu) real(rr(i1,i2), d_size)
    end do
  end do
  ! Vort
  do i2 = 1 , size(vort)
    write(fu) real(vort(i2),d_size)
  end do

  !Connectivity
  do i1 = 1 , size(ee,2)
    if (ee(4, i1) .eq. 0) then
      !it is a triangle, write two times the last node
      write(fu) int(ee(1:3,i1)-1, s_size)
      write(fu) int(ee(3,i1)-1, s_size)
    else
      write(fu) int(ee(:,i1)-1, s_size)
    endif
  end do
 

  !WAKE
  !Zone marker
  write(fu) zoneMarker

  !Size of the variables: 2=double
  write(fu) int(2,s_size)
  write(fu) int(2,s_size)
  write(fu) int(2,s_size)
  write(fu) int(2,s_size)

  !Has passive variables?
  write(fu) int(0,s_size)

  !If has passive vars here must be the list of the passive vars

  !Has variable sharing
  write(fu) int(0,s_size)

  !sharing connectivity? -1 no sharing
  write(fu) int(-1,s_size)
  
  !Min and max of all the variables
  !x
  write(fu) dble(minval( w_rr(1,:) ))
  write(fu) dble(maxval( w_rr(1,:) ))
  !y
  write(fu) dble(minval( w_rr(2,:) ))
  write(fu) dble(maxval( w_rr(2,:) ))
  !z
  write(fu) dble(minval( w_rr(3,:) ))
  write(fu) dble(maxval( w_rr(3,:) ))
  !Vort
  write(fu) dble(minval( w_vort ))
  write(fu) dble(maxval( w_vort ))

  !The actual data
  ! x, y, z
  do i1 = 1 , 3
    do i2 = 1 , size(w_rr,2)
      write(fu) real(w_rr(i1,i2), d_size)
    end do
  end do
  ! Vort
  do i2 = 1 , size(w_vort)
    write(fu) real(w_vort(i2),d_size)
  end do

  !Connectivity
  do i1 = 1 , size(w_ee,2)
    if (w_ee(4, i1) .eq. 0) then
      !it is a triangle, write two times the last node
      write(fu) int(w_ee(1:3,i1)-1, s_size)
      write(fu) int(w_ee(3,i1)-1, s_size)
    else
      write(fu) int(w_ee(:,i1)-1, s_size)
    endif
  end do



  close(fu)
  

end subroutine tec_out_sol_bin

!----------------------------------------------------------------------

subroutine tec_out_viz(out_filename, t, &
                       rr, ee, vars, var_names,  &
                     w_rr, w_ee, w_vars, w_var_names)
 character(len=*), intent(in) :: out_filename 
 real(wp), intent(in) :: t
 real(wp), intent(in) :: rr(:,:)
 integer, intent(in)  :: ee(:,:)
 real(wp), intent(in) :: vars(:,:)
 character(len=*), intent(in) :: var_names(:)
 real(wp), intent(in), optional :: w_rr(:,:)
 integer, intent(in), optional  :: w_ee(:,:)
 real(wp), intent(in), optional :: w_vars(:,:)
 character(len=*), intent(in), optional :: w_var_names(:)

 character, parameter :: zc = char(0)
 integer :: npoints, ncells, ne
 integer :: fu, ierr, i1, i2
 real(kind=4)  , parameter :: zoneMarker = 299.0
 real(kind=4)  , parameter :: eohMarker  = 357.0
 character(len=max_char_len) :: buffer_char
 integer :: ie, nquad, ntria
 integer npoints_w, nw
 logical :: got_wake
 integer :: nvars, w_nvars, iv


  got_wake = present(w_var_names)

  ne = size(ee,2) !number of elements
  npoints = size(rr,2)
  ncells = ne
  nvars = size(var_names)
  
  if(got_wake) then
    nw = size(w_ee,2) !number of wake elements
    npoints_w = size(w_rr,2)
    w_nvars = size(w_var_names)
  endif

  !WARNING: for the moment we are assuming that the number and type of 
  ! wake variables is lower or equal to the body ones

  ! First cycle the elements to get the number of quads and trias
  ! (this is really ugly, make it more flexible)
  nquad = 0; ntria = 0
  do ie = 1,ne
    if(ee(4,ie) .eq. 0) then
      ntria = ntria+1
    else
      nquad = nquad+1
    endif
  enddo

  call new_file_unit(fu,ierr)
  open(unit=fu,file=trim(out_filename),status='replace',access='stream', &
       form='unformatted',iostat=ierr)

  !magic number
  buffer_char = "#!TDV112" ;  write(fu) trim(buffer_char) 

  !integer 1 to set the byte order
  write(fu) int(1,s_size)

  !integer 0 to mark a full file (1=solution only, 2=grid only)
  write(fu) int(0,s_size)

  !title
  call put_tec_string('DUST solution',fu)

  !number of variables
  write(fu) int(3+nvars,s_size) !x-y-z-vort

  !Variables names
  call put_tec_string('X',fu)
  call put_tec_string('Y',fu)
  call put_tec_string('Z',fu)
  do iv = 1,nvars
    call put_tec_string(trim(var_names(iv)),fu)
  enddo
 
  !BODY
  !Zone marker
  write(fu) zoneMarker

  !Zone name
  call put_tec_string('Body',fu)

  !Parent zone: -1 for no parent (not too clear)
  write(fu) int(-1,s_size)

  !Strand id: -2 for automatic generation
  write(fu) int(-2,s_size)

  !Solution time
  write(fu) real(t,d_size)

  !"not used, set to -1"
  write(fu) int(-1,s_size) 
  
  !Zone type: 0 ordered, 3 unstructured quadrilateral
  write(fu) int(3,s_size)

  !Toggle variable location specification
  write(fu) int(1,s_size)

  !x, y, z stand at nodes
  write(fu) int(0,s_size)
  write(fu) int(0,s_size)
  write(fu) int(0,s_size)
  ! Other variables stand at cells centers
  do iv = 1,nvars
    write(fu) int(1,s_size)
  enddo

  !are local neighbour supplied? No.
  write(fu) int(0,s_size)

  !Number of miscellaneous somethings
  write(fu) int(0,s_size)

  ! n.points, n.elems
  write(fu) int(size(rr,2),s_size)
  write(fu) int(size(ee,2),s_size)
  
  ! All zero, three fields, for future use
  write(fu) int(0,s_size) , int(0,s_size) , int(0,s_size)
  
  !"no more auxilliary name/value pairs"
  write(fu) int(0,s_size)


  !WAKE ZONE HERE
  if (got_wake) then
    !Zone marker
    write(fu) zoneMarker

    !Zone name
    call put_tec_string('Wake',fu)

    !Parent zone: -1 for no parent (not too clear)
    write(fu) int(-1,s_size)

    !Strand id: -2 for automatic generation
    write(fu) int(-2,s_size)

    !Solution time
    write(fu) real(t,d_size)

    !"not used, set to -1"
    write(fu) int(-1,s_size) 
    
    !Zone type: 0 ordered, 3 unstructured quadrilateral
    write(fu) int(3,s_size)

    !Toggle variable location specification
    write(fu) int(1,s_size)

    !x, y, z stand at nodes
    write(fu) int(0,s_size)
    write(fu) int(0,s_size)
    write(fu) int(0,s_size)
    ! Other variables stand at cells centers
    do iv = 1,nvars
      write(fu) int(1,s_size)
    enddo

    !are local neighbour supplied? No.
    write(fu) int(0,s_size)

    !Number of miscellaneous somethings
    write(fu) int(0,s_size)

    ! n.points, n.elems
    write(fu) int(size(w_rr,2),s_size)
    write(fu) int(size(w_ee,2),s_size)
    
    ! All zero, three fields, for future use
    write(fu) int(0,s_size) , int(0,s_size) , int(0,s_size)
    
    !"no more auxilliary name/value pairs"
    write(fu) int(0,s_size)
  endif


  !DATA
  write(fu) eohMarker
  
  !BODY
  !Zone marker
  write(fu) zoneMarker

  !Size of the variables: 2=double
  write(fu) int(2,s_size)
  write(fu) int(2,s_size)
  write(fu) int(2,s_size)
  do iv = 1,nvars
    write(fu) int(2,s_size)
  enddo

  !Has passive variables?
  write(fu) int(0,s_size)

  !If has passive vars here must be the list of the passive vars

  !Has variable sharing
  write(fu) int(0,s_size)

  !sharing connectivity? -1 no sharing
  write(fu) int(-1,s_size)
  
  !Min and max of all the variables
  !x
  write(fu) dble(minval( rr(1,:) ))
  write(fu) dble(maxval( rr(1,:) ))
  !y
  write(fu) dble(minval( rr(2,:) ))
  write(fu) dble(maxval( rr(2,:) ))
  !z
  write(fu) dble(minval( rr(3,:) ))
  write(fu) dble(maxval( rr(3,:) ))
  !Vars
  do iv = 1,nvars
    write(fu) dble(minval( vars(:,iv) ))
    write(fu) dble(maxval( vars(:,iv) ))
  enddo

  !The actual data
  ! x, y, z
  do i1 = 1 , 3
    do i2 = 1 , size(rr,2)
      write(fu) real(rr(i1,i2), d_size)
    end do
  end do
  ! Vars
  do iv = 1,nvars
    do i2 = 1 , size(vars,1)
      write(fu) real(vars(i2,iv),d_size)
    end do
  enddo

  !Connectivity
  do i1 = 1 , size(ee,2)
    if (ee(4, i1) .eq. 0) then
      !it is a triangle, write two times the last node
      write(fu) int(ee(1:3,i1)-1, s_size)
      write(fu) int(ee(3,i1)-1, s_size)
    else
      write(fu) int(ee(:,i1)-1, s_size)
    endif
  end do
 

  !WAKE
  if(got_wake) then
    !Zone marker
    write(fu) zoneMarker

    !Size of the variables: 2=double
    write(fu) int(2,s_size)
    write(fu) int(2,s_size)
    write(fu) int(2,s_size)
    do iv = 1,nvars
      write(fu) int(2,s_size)
    enddo
    
    if (w_nvars .lt. nvars) then
      !Has passive variables?
      write(fu) int(1,s_size)
      !x-y-z are not passive
      write(fu) int(0,s_size)
      write(fu) int(0,s_size)
      write(fu) int(0,s_size)
      !first variables are not passive
      do iv=1,w_nvars
        write(fu) int(0,s_size)
      enddo
      !last variables are passive
      do iv = 1,(nvars-w_nvars)
        write(fu) int(1,s_size)
      enddo
    else
      !Has passive variables?
      write(fu) int(0,s_size)
    endif

    !Has variable sharing
    write(fu) int(0,s_size)

    !sharing connectivity? -1 no sharing
    write(fu) int(-1,s_size)
    
    !Min and max of all the variables
    !x
    write(fu) dble(minval( w_rr(1,:) ))
    write(fu) dble(maxval( w_rr(1,:) ))
    !y
    write(fu) dble(minval( w_rr(2,:) ))
    write(fu) dble(maxval( w_rr(2,:) ))
    !z
    write(fu) dble(minval( w_rr(3,:) ))
    write(fu) dble(maxval( w_rr(3,:) ))
    !Vars
    do iv = 1,w_nvars
      write(fu) dble(minval( vars(:,iv) ))
      write(fu) dble(maxval( vars(:,iv) ))
    enddo

    !The actual data
    ! x, y, z
    do i1 = 1 , 3
      do i2 = 1 , size(w_rr,2)
        write(fu) real(w_rr(i1,i2), d_size)
      end do
    end do
    ! Vort
    do iv = 1,w_nvars
      do i2 = 1 , size(w_vars,1)
        write(fu) real(w_vars(i2,iv),d_size)
      end do
    enddo

    !Connectivity
    do i1 = 1 , size(w_ee,2)
      if (w_ee(4, i1) .eq. 0) then
        !it is a triangle, write two times the last node
        write(fu) int(w_ee(1:3,i1)-1, s_size)
        write(fu) int(w_ee(3,i1)-1, s_size)
      else
        write(fu) int(w_ee(:,i1)-1, s_size)
      endif
    end do

    close(fu)
  endif
  

end subroutine tec_out_viz

!----------------------------------------------------------------------

subroutine tec_out_box(out_filename, time, xpoints, ypoints, zpoints, vars, var_names)
 character(len=*), intent(in) :: out_filename 
 real(wp), intent(in) :: time
 real(wp), intent(in) :: xpoints(:)
 real(wp), intent(in) :: ypoints(:)
 real(wp), intent(in) :: zpoints(:)
 real(wp), intent(in) :: vars(:,:)
 character(len=*), intent(in) :: var_names(:)
 
 character, parameter :: zc = char(0)
 integer :: fu, ierr, i1, i2, i3
 real(kind=4)  , parameter :: zoneMarker = 299.0
 real(kind=4)  , parameter :: eohMarker  = 357.0
 character(len=max_char_len) :: buffer_char
 integer :: nvars, iv


  nvars = size(var_names)

  call new_file_unit(fu,ierr)
  open(unit=fu,file=trim(out_filename),status='replace',access='stream', &
       form='unformatted',iostat=ierr)

  !magic number
  buffer_char = "#!TDV112" ;  write(fu) trim(buffer_char) 

  !integer 1 to set the byte order
  write(fu) int(1,s_size)

  !integer 0 to mark a full file (1=solution only, 2=grid only)
  write(fu) int(0,s_size)

  !title
  call put_tec_string('DUST flowfield',fu)

  !number of variables
  write(fu) int(3+nvars,s_size) !x y z+vars

  !Variables names
  call put_tec_string('X',fu)
  call put_tec_string('Y',fu)
  call put_tec_string('Z',fu)
  do iv = 1,nvars
    call put_tec_string(trim(var_names(iv)),fu)
  enddo
  
  !Zone marker
  write(fu) zoneMarker

  !Zone name
  call put_tec_string('flowfield',fu)

  !Parent zone: -1 for no parent (not too clear)
  write(fu) int(-1,s_size)

  !Strand id: -2 for automatic generation
  write(fu) int(-2,s_size)

  !Solution time
  write(fu) real(time,d_size)

  !"not used, set to -1"
  write(fu) int(-1,s_size) 
  
  !Zone type: 0 ordered, 3 unstructured quadrilateral
  write(fu) int(0,s_size)

  !Toggle variable location specification: no, default on nodes
  write(fu) int(0,s_size)

  !are local neighbour supplied? No.
  write(fu) int(0,s_size)

  !Number of miscellaneous somethings
  write(fu) int(0,s_size)

  !Ordered zone: Imax, Jmax, Kmax
  write(fu) int(size(xpoints),s_size)
  write(fu) int(size(ypoints),s_size)
  write(fu) int(size(zpoints),s_size)

  !"no more auxilliary name/value pairs"
  write(fu) int(0,s_size)


  !DATA
  write(fu) eohMarker

  !Zone marker
  write(fu) zoneMarker

  !Size of the variables: 2=double
  write(fu) int(2,s_size) !x
  write(fu) int(2,s_size) !y
  write(fu) int(2,s_size) !z
  do iv = 1,nvars
    write(fu) int(2,s_size)
  enddo

  !Has passive variables?
  write(fu) int(0,s_size)

  !Has variable sharing?
  write(fu) int(0,s_size)

  !sharing connectivity? -1 no sharing
  write(fu) int(-1,s_size)

  !Min e max
  write(fu) dble(minval( xpoints ))
  write(fu) dble(maxval( xpoints ))
  write(fu) dble(minval( ypoints ))
  write(fu) dble(maxval( ypoints ))
  write(fu) dble(minval( zpoints ))
  write(fu) dble(maxval( zpoints ))
  do iv = 1,nvars
    write(fu) dble(minval( vars(iv,:) ))
    write(fu) dble(maxval( vars(iv,:) ))
  enddo

  !The actual data
  do i3 = 1,size(zpoints); do i2 = 1,size(ypoints); do i1 = 1, size(xpoints)
    write(fu) real(xpoints(i1), d_size)
    write(fu) real(ypoints(i2), d_size)
    write(fu) real(zpoints(i3), d_size)
  enddo; enddo; enddo

  do iv = 1,nvars
    do i1 = 1, size(vars,2)
      write(fu) real(vars(iv,i1),d_size)
    enddo
  enddo

end subroutine

!---------------------------------------------------------------------

!---------------------------------------------------------------------

subroutine tec_out_probes(out_filename, time, vars, var_names, zone_names)
 character(len=*), intent(in) :: out_filename 
 real(wp), intent(in) :: time(:)
 real(wp), intent(in) :: vars(:,:,:)
 character(len=*), intent(in) :: var_names(:)
 character(len=*), intent(in) :: zone_names(:)
 
 character, parameter :: zc = char(0)
 integer :: fu, ierr, i1
 real(kind=4)  , parameter :: zoneMarker = 299.0
 real(kind=4)  , parameter :: eohMarker  = 357.0
 character(len=max_char_len) :: buffer_char
 integer :: nvars, ip, nprobes, varlen, iv


  nvars = size(var_names)
  nprobes = size(zone_names)
  varlen = size(vars, 2)

  call new_file_unit(fu,ierr)
  open(unit=fu,file=trim(out_filename),status='replace',access='stream', &
       form='unformatted',iostat=ierr)

  !magic number
  buffer_char = "#!TDV112" ;  write(fu) trim(buffer_char) 

  !integer 1 to set the byte order
  write(fu) int(1,s_size)

  !integer 0 to mark a full file (1=solution only, 2=grid only)
  write(fu) int(0,s_size)

  !title
  call put_tec_string('DUST probes',fu)

  !number of variables
  write(fu) int(1+nvars,s_size) !t+vars

  !Variables names
  call put_tec_string('t',fu)
  do iv = 1,nvars
    call put_tec_string(trim(var_names(iv)),fu)
  enddo
  
  !A single zone for each probe
  do ip = 1,nprobes
    !Zone marker
    write(fu) zoneMarker

    !Zone name
    call put_tec_string(trim(zone_names(ip)),fu)

    !Parent zone: -1 for no parent (not too clear)
    write(fu) int(-1,s_size)

    !Strand id: -2 for automatic generation
    write(fu) int(-2,s_size)

    !Solution time, useless in this case
    write(fu) real(0.0,d_size)

    !"not used, set to -1"
    write(fu) int(-1,s_size) 
    
    !Zone type: 0 ordered, 3 unstructured quadrilateral
    write(fu) int(0,s_size)

    !Toggle variable location specification: no, default on nodes
    write(fu) int(0,s_size)

    !are local neighbour supplied? No.
    write(fu) int(0,s_size)

    !Number of miscellaneous somethings
    write(fu) int(0,s_size)

    !Ordered zone: Imax, Jmax, Kmax
    write(fu) int(varlen,s_size)
    write(fu) int(1,s_size)
    write(fu) int(1,s_size)

    !Number of elemenst
    !write(fu) int(varlen,s_size)

    ! All zero, three fields, for future use
    !write(fu) int(0,s_size) , int(0,s_size) , int(0,s_size)

    !"no more auxilliary name/value pairs"
    write(fu) int(0,s_size)

  enddo

  !DATA
  write(fu) eohMarker

  do ip=1,nprobes

    !Zone marker
    write(fu) zoneMarker

    !Size of the variables: 2=double
    write(fu) int(2,s_size) !time
    do iv = 1,nvars
      write(fu) int(2,s_size)
    enddo

    !Has passive variables?
    write(fu) int(0,s_size)

    !Has variable sharing?
    write(fu) int(0,s_size)

    !sharing connectivity? -1 no sharing
    write(fu) int(-1,s_size)

    !Min e max
    write(fu) dble(minval( time ))
    write(fu) dble(maxval( time ))
    do iv = 1,nvars
      write(fu) dble(minval( vars(iv,:,ip) ))
      write(fu) dble(maxval( vars(iv,:,ip) ))
    enddo

    !The actual data
    do i1 = 1, varlen
      write(fu) real(time(i1), d_size)
    enddo
    do iv = 1,nvars
      do i1 = 1, varlen
        write(fu) real(vars(iv,i1,ip),d_size)
      enddo
    enddo


  enddo

end subroutine

!---------------------------------------------------------------------

subroutine put_tec_string(str, fu)
 character(len=*), intent(in) :: str
 integer, intent(in)          :: fu

 integer :: sl, i
 character, parameter :: zc = char(0)
 
 sl = len(trim(str))

 do i = 1, sl
   write(fu) ichar(str(i:i))
 enddo 
 write(fu) ichar(zc)

end subroutine put_tec_string

end module mod_tecplot_out
