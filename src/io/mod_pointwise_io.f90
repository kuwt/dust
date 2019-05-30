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

! list of subroutines and functions
! - read_mesh_pointwise
! - read_points
! - read_lines
! - set_parser_pointwise

module mod_pointwise_io

use mod_param, only: &
  wp, max_char_len, nl, pi

use mod_handling, only: &
  error, warning, info, printout, new_file_unit, check_file_exists

use mod_parse, only: &
  t_parse, getstr, getint, getreal, getrealarray, getlogical, &
  countoption, getsuboption, getintarray

use mod_spline, only: &
  t_spline, hermite_spline , deallocate_spline

!----------------------------------------------------------------------

implicit none

public :: read_mesh_pointwise

private

character(len=*), parameter :: this_mod_name = 'mod_pointwise_io'

!----------------------------------------------------------------------
!> module variables


!----------------------------------------------------------------------
!> point type
type :: t_point
  integer                 :: id
  real(wp)                :: coord(3)
  character(max_char_len) :: airfoil
  real(wp)                :: chord   !
  real(wp)                :: theta   ! deg
  real(wp)                :: sec_nor(3)
end type t_point

!> line type
type :: t_line
  character(max_char_len) :: l_type
  integer                 :: end_points(2)
  integer                 :: nelems
  character(max_char_len) :: type_span  ! discretisation in span
  integer                 :: neigh_line(2) = 0
  real(wp), allocatable   :: t_vec1(:) , t_vec2(:)
end type t_line


!----------------------------------------------------------------------

contains

!----------------------------------------------------------------------
!> read pointwise input file and build ee, rr arrays of connectivity
!  and point coordinates
subroutine read_mesh_pointwise ( mesh_file , ee , rr , &
                     npoints_chord_tot , nelem_span_tot )
              
 character(len=*), intent(in) :: mesh_file
 integer  , allocatable, intent(out) :: ee(:,:) 
 real(wp) , allocatable, intent(out) :: rr(:,:) 
 integer  , intent(out), optional    :: npoints_chord_tot, nelem_span_tot

 ! parser and sub-parsers
 type(t_parse) :: pmesh_prs
 type(t_parse) , pointer :: point_prs , line_prs

 !>
 integer   :: nelem_chord
 character :: ElType
 real(wp)  :: ref_chord_fraction

 !> point and line structures
 type(t_point) , allocatable :: points(:)
 type(t_line ) , allocatable :: lines(:)

 !> ee, rr size
 integer ::  nelem_chord_tot, npoint_chord_tot, npoint_span_tot
 integer :: ee_size , rr_size
 integer :: i_ch , i_sp , i_point , i_elem

 real(wp) , allocatable :: ref_line_points(:,:)
 real(wp) , allocatable :: ref_line_normal(:,:)


 integer :: i 

 
 !> Prepare parser and sub-parsers
 call set_parser_pointwise( pmesh_prs , point_prs , line_prs )

 !> Enable option reading 
 call pmesh_prs%read_options(mesh_file,printout_val=.true.)

 nelem_chord = getint(pmesh_prs,'nelem_chord')
 ElType  = trim(getstr(pmesh_prs,'ElType'))
 ref_chord_fraction = getreal(pmesh_prs,'reference_chord_fraction')
 
 !> Read points and lines 
 call read_points ( pmesh_prs , point_prs , points )
 call read_lines  ( pmesh_prs ,  line_prs , lines  , nelem_span_tot )


 !> Set dimensions of ee, rr mat
 if (ElType.eq.'p') then ;  nelem_chord_tot = 2 * nelem_chord 
 else                    ;  nelem_chord_tot =     nelem_chord
 endif

 npoint_chord_tot = nelem_chord_tot + 1
 npoint_span_tot  = nelem_span_tot  + 1

 ee_size =  nelem_chord_tot *  nelem_span_tot
 rr_size = npoint_chord_tot * npoint_span_tot

! ! check ---
! write(*,*) ' ee_size , rr_size : ' , ee_size , rr_size

 ! === connectivity matrix, ee ===
 allocate( ee( 4 , ee_size ) ) ; ee = 0
 do i_ch = 1 , nelem_chord_tot
  do i_sp = 1 , nelem_span_tot
   i_elem  =  nelem_chord_tot*( i_sp - 1 ) + i_ch
   i_point = npoint_chord_tot*( i_sp - 1 ) + i_ch
   ee(1,i_elem) = i_point + npoint_chord_tot         ! iPoint
   ee(2,i_elem) = i_point                            ! iPoint + 1
   ee(3,i_elem) = i_point + 1                        ! iPoint + npoint_chord_tot + 1
   ee(4,i_elem) = i_point + npoint_chord_tot + 1     ! iPoint + npoint_chord_tot
  end do
 end do

 ! === define the reference line ===
 !> Check line connectivity, re-order lines and define tangent vectors
 call sort_lines( lines ) 

 !> check input consistency (todo: improve this subroutine)
 call check_point_line_inputs ( points , lines )

 !>  
 call fill_line_tan_vec( points , lines )

 !>
 call build_reference_line( npoint_span_tot , &
                            points , lines  , &
                            ref_line_points , &
                            ref_line_normal     )


  ! check ---
  do i = 1 , size(ref_line_points,1)
    write(*,*) ref_line_points(i,:) , ref_line_normal(i,:)
  end do




 write(*,*) ' stop in mod_pointwise_io ' ; stop


end subroutine read_mesh_pointwise

!----------------------------------------------------------------------
!> build reference line
subroutine build_reference_line( npoint_span_tot , &
                                 points, lines   , &
                                 ref_line_points , &
                                 ref_line_normal     )

 integer                     , intent(in)    :: npoint_span_tot
 type(t_point)               , intent(inout) :: points(:)
 type(t_line )               , intent(inout) :: lines( :)
 real(wp)      , allocatable , intent(out)   :: ref_line_points(:,:)
 real(wp)      , allocatable , intent(out)   :: ref_line_normal(:,:)

 type(t_spline) :: spl

 integer :: i , i1 , i2 , j , n

!write(*,*) ' n. point in spanwise direction: ' , npoint_span_tot
 allocate(ref_line_points(npoint_span_tot,3))
 allocate(ref_line_normal(npoint_span_tot,3))
 
 !> Starting point
 ref_line_points(1,:) = points( lines(1)%end_points(1) ) % coord
 i2 = 1 
 do i = 1 , size(lines)
   ! end points
   i1 = i2
   i2 = i1 + lines(i)%nelems
   ref_line_points(i2,:) = points( lines(i)%end_points(2) ) % coord

   if (      trim(lines(i)%l_type) .eq. 'Straight' ) then
     
     call straight_line( points( lines(i)%end_points(1) )%coord , &
                         points( lines(i)%end_points(2) )%coord , &
                         lines(i)%nelems                        , &
                         lines(i)%type_span                     , &
                         ref_line_points(i1:i2,:)               , &
                         ref_line_normal(i1:i2,:) )
 
   else if ( trim(lines(i)%l_type) .eq. 'Spline'   ) then
     
     !> build t_spline structure
     n = lines(i)%end_points(2) - lines(i)%end_points(1) + 1
     allocate(spl%rr( n , 3 ))
     allocate(spl%d0( 3 )) ; spl%d0 = lines(i)%t_vec1
     allocate(spl%d1( 3 )) ; spl%d1 = lines(i)%t_vec2
     spl%e_bc = (/ 'derivative' , 'derivative' /) 

     do j = 1 , n
       spl%rr( j , : ) = points( lines(i)%end_points(1)-1+j )%coord
     end do

     !> compute ref_line_points on the spline
     call hermite_spline( spl , lines(i)%nelems           , &
                                lines(i)%type_span        , &
                                ref_line_points(i1:i2,:)  , &
                                ref_line_normal(i1:i2,:) )

     call deallocate_spline( spl )
     
   else
     write(*,*) ' error in build_reference_line(): '
     write(*,*) ' lines(',i,')%l_type = ' , lines(i)%l_type
     write(*,*) ' but l_type must be either Spline or Straight.'
     write(*,*) ' Stop. '; stop
   end if

 end do 

end subroutine build_reference_line

!----------------------------------------------------------------------
!> straight_line subdivision
subroutine straight_line( r1 , r2 , nelems , type_span , rr , nor )
  real(wp)                , intent(in) :: r1(3) , r2(3)
  integer                 , intent(in) :: nelems
  character(max_char_len) , intent(in) :: type_span
  real(wp)                , intent(inout) :: rr(:,:)
  real(wp)                , intent(inout) ::nor(:,:)

  real(wp) :: nor_v(3)

  integer :: i

  ! allocate(rr(nelems+1,3))

  nor_v = r2 - r1 ; nor_v = nor_v / norm2(nor_v)

  do i = 1 , nelems+1

    !> rr
    if ( trim(type_span) .eq. 'uniform' ) then !> uniform spacing
      rr(i,:) = r1 * dble(nelems+1-i)/dble(nelems) + &
                r2 * dble(       i-1)/dble(nelems)
    else
      write(*,*) ' error in straight_line. Only uniform spacing '
      write(*,*) ' implemented so far. Stop ' ; stop
    end if

    !> nor
    nor(i,:) = nor_v

  end do

end subroutine straight_line

!----------------------------------------------------------------------
!> sort lines: first line has the first point with id = 1
!  and then chain last-to-first points until the last line
subroutine sort_lines( lines )
  type(t_line ) , intent(inout) :: lines( :)
  type(t_line ) , allocatable :: lines_tmp(:)
  integer :: n_lines 
  integer :: i , i1 , i2 , i2_old , il , n_1
 
  n_lines = size(lines)

  allocate(lines_tmp(n_lines))
  lines_tmp = lines

  !> some preliminary checks:
  ! -> only one line has end_points(1) = 1 -> first line
  n_1 = 0
  do i = 1 , n_lines
    if ( lines_tmp(i)%end_points(1) .eq. 1 ) then
      n_1 = n_1 + 1 ; i1 = i
    end if
  end do
  if ( n_1 .ne. 1 ) then
    write(*,*) ' error in sort_lines: only 1 line must have'
    write(*,*) ' end_points(1) = 1. Stop ' ; stop
  end if  
  lines(1) = lines_tmp(i1)
  ! -> no multiple start/end points
  ! ...

  !> find the chain of lines 
  i1 = lines(1)%end_points(1) ; i2 = lines(1)%end_points(2)
  do i = 2 , n_lines

    il = 0 ; i2_old = i2 ! re-initialisation

    do while ( ( i2 .eq. i2_old ) .and. ( il .lt. n_lines ) )
      il = il + 1
      if (      lines_tmp(il)%end_points(1) .eq. i2_old ) then
        i1 =    lines_tmp(il)%end_points(1) 
        i2 =    lines_tmp(il)%end_points(2)
      end if

    end do
    !> 
    if ( i2 .eq. i2_old ) then
      write(*,*) ' error in sort_lines: broken chain. stop ' ; stop
    end if

    lines(i) = lines_tmp(il)

  end do

  


end subroutine sort_lines

!----------------------------------------------------------------------
!> check input consistency
subroutine check_point_line_inputs( points , lines )
  type(t_point) , intent(inout) :: points(:)
  type(t_line ) , intent(inout) :: lines( :)

  integer :: i , n_lines

  n_lines = size(lines)

  !> all the points are present
  if ( size(points) .lt. lines(size(lines))%end_points(2) ) then
    write(*,*) ' n. points .lt. last line end point id. stop ' ; stop
  end if
  
  !> if first or last lines are Splines, the user must provide the 
  !  initial or final direction
  if ( ( trim(lines(1)%l_type) .eq. 'Spline' ) .and. &
       ( .not. allocated(lines(1)%t_vec1) ) ) then
    write(*,*) ' error in check_point_lines_input: the first line is a Spline and &
                &has no TangentVec1 provided as an input. Stop ' ; stop
  end if
  if ( ( trim(lines(size(lines))%l_type) .eq. 'Spline' ) .and. &
       ( .not. allocated(lines(size(lines))%t_vec2) ) ) then
    write(*,*) ' error in check_point_lines_input: the last line is a Spline and &
                &has no TangentVec2 provided as an input. Stop. ' ; stop
  end if

  !> no consecutive splines are allowed so far
  do i = 1 , n_lines-1
    if ( ( trim(lines(i  )%l_type) .eq. 'Spline' ) .and. &
         ( trim(lines(i+1)%l_type) .eq. 'Spline' ) ) then
      write(*,*) ' error in check_point_lines_input: two consecutive splines are not &
                  &allowed so far. '
      write(*,*) ' In future release, this coul not be true anymore: '
      write(*,*) ' two consecutive splines could be joined together if no'
      write(*,*) ' TangentVec is provided, or two consecutive spline with'
      write(*,*) ' TangentVec provided could be treated as well. So far, Stop.'
      stop
    end if
  end do

  

 
end subroutine

!----------------------------------------------------------------------
!> 
subroutine fill_line_tan_vec( points , lines )
 type(t_point) , intent(inout) :: points(:)
 type(t_line ) , intent(inout) :: lines( :)

 integer :: n_lines , i

 n_lines = size(lines)

 ! first straight lines
 do i = 1 , n_lines

   if ( trim(lines(i)%l_type) .eq. 'Straight' ) then

     allocate( lines(i)%t_vec1(3) , lines(i)%t_vec2(3) ) ;

     lines(i)%t_vec1 = points( lines(i)%end_points(2) ) %coord - &
                       points( lines(i)%end_points(1) ) %coord

     if ( norm2(lines(i)%t_vec1) .lt. 1.0e-9_wp ) then
       write(*,*) ' error in fill_line_tan_vec ... stop ' ; stop
     end if

     lines(i)%t_vec1 = lines(i)%t_vec1 / norm2(lines(i)%t_vec1)
     lines(i)%t_vec2 = lines(i)%t_vec1

   end if 

 end do

 !> then splines
 do i = 1 , n_lines

   if ( trim(lines(i)%l_type) .eq. 'Spline' ) then

     if ( .not. allocated(lines(i)%t_vec1) ) then
 
       if  ( i .eq. 1 ) then  ! check
         write(*,*) ' error in fill_line_tan_vec: '
         write(*,*) ' first line is a spline w/o tangentVec1 input '
         write(*,*) ' stop. ' ; stop
       end if

       allocate(lines(i)%t_vec1(3))
       lines(i)%t_vec1 = lines(i-1)%t_vec2
       
     end if

     if ( .not. allocated(lines(i)%t_vec2) ) then

       if  ( i .eq. n_lines ) then  ! check
         write(*,*) ' error in fill_line_tan_vec: '
         write(*,*) ' last line is a spline w/o tangentVec2 input '
         write(*,*) ' stop. ' ; stop

       end if

       allocate(lines(i)%t_vec2(3))
       lines(i)%t_vec2 = lines(i+1)%t_vec1

     end if

   end if

 end do



end subroutine fill_line_tan_vec


!----------------------------------------------------------------------
!> fill point structure
subroutine read_points ( pmesh_prs , point_prs , points )
 type(t_parse) ,               intent(inout) :: pmesh_prs
 type(t_parse) , pointer     , intent(inout) :: point_prs
 type(t_point) , allocatable , intent(out)   :: points(:)

 integer :: nPoints , i

 ! === Read point groups ===
 nPoints = countoption(pmesh_prs,'Point') ; allocate( points(nPoints) )

 ! loop over Point groups
 do i = 1 , nPoints
   call getsuboption( pmesh_prs , 'Point' , point_prs )
   points(i) % id       = getint(      point_prs, 'Id')
   points(i) % coord    = getrealarray(point_prs, 'Coordinates',3)
   points(i) % airfoil  = getstr(      point_prs, 'Airfoil')
   points(i) % chord    = getreal(     point_prs, 'Chord')
   points(i) % theta    = getreal(     point_prs, 'Twist')
   points(i) % sec_nor  = getrealarray(point_prs, 'SectionNormal',3)
 end do

end subroutine read_points

!----------------------------------------------------------------------
!> fill line structure
subroutine read_lines ( pmesh_prs , line_prs , lines  , nelems_span_tot)
 type(t_parse) ,               intent(inout) :: pmesh_prs
 type(t_parse) , pointer     , intent(inout) ::  line_prs
 type(t_line ) , allocatable , intent(out)   :: lines(:)
 integer                     , intent(out)   :: nelems_span_tot

 integer :: nLines , i

 ! === Read line groups ===
 nLines  = countoption(pmesh_prs,'Line') ; allocate( lines(nLines) ) 

 nelems_span_tot = 0
 do i = 1 , nLines

   call getsuboption( pmesh_prs , 'Line' , line_prs )

   lines(i) % l_type     = getstr(      line_prs , 'Type'   )
   lines(i) % end_points = getintarray( line_prs , 'EndPoints' , 2 )
   lines(i) % nelems     = getint(      line_prs , 'Nelems' )
   lines(i) % type_span  = getstr(      line_prs , 'Type_span')

   !> allocate tvec for spline if they are provided as a input
   if ( trim(lines(i)%l_type) .eq.'Spline' ) then 
 
     !> allocate tangent vec1 if specified as an input 
     if ( countoption( line_prs , 'TangentVec1' ) .eq. 0 ) then 
       ! do nothing: tan vec must be inherited 
     elseif ( countoption( line_prs , 'TangentVec1' ) .eq. 1 ) then
       allocate(lines(i)%t_vec1(3))
       lines(i)%t_vec1 = getrealarray(line_prs , 'TangentVec1' , 3 )
     else
       write(*,*) ' Error in read_lines(). Provided more than one &
                   &TangentVec1 as an input. Stop' ; stop
     end if

     !> allocate tangent vec2 if specified as an input 
     if ( countoption( line_prs , 'TangentVec2' ) .eq. 0 ) then
       ! do nothing: tan vec must be inherited 
     elseif ( countoption( line_prs , 'TangentVec2' ) .eq. 1 ) then
       allocate(lines(i)%t_vec2(3))
       lines(i)%t_vec2 = getrealarray(line_prs , 'TangentVec2' , 3 )
     else
       write(*,*) ' Error in read_lines(). Provided more than one &
                   &TangentVec2 as an input. Stop' ; stop
     end if
 
   end if

   nelems_span_tot = nelems_span_tot + lines(i)%nelems

 end do

end subroutine read_lines

!----------------------------------------------------------------------
!> prepare all the parser (and sub-parsers) options
subroutine set_parser_pointwise( pmesh_prs , point_prs , line_prs )
 type(t_parse)           , intent(out) :: pmesh_prs
 type(t_parse) , pointer , intent(out) :: point_prs , line_prs


 ! === Prepare fields to be read by the parser ===
 call pmesh_prs%CreateStringOption('ElType', &
               'element type (temporary) p panel v vortex ring', &
               multiple=.false.)
 call pmesh_prs%CreateIntOption('nelem_chord',  &
               'number of chord-wise elements', &
               multiple=.false.)
 call pmesh_prs%CreateStringOption('type_chord', &
               'type of chord-wise division: uniform, cosine, cosineLE, cosineTE',&
               'uniform', &
               multiple=.false.)
 ! no need for a starting point
 call pmesh_prs%CreateRealOption('reference_chord_fraction',&
              'Reference chord fraction', &
              '0.0',&
              multiple=.false.)
 
 ! === Point subparser === 
 call pmesh_prs%CreateSubOption('Point','Point group',point_prs, &
              multiple=.true.)

 call point_prs%CreateIntOption('Id', 'Point Id.' )
 call point_prs%CreateRealArrayOption('Coordinates', &
               'coordinates of the points used to define the comp' )
 call point_prs%CreateStringOption(  'airfoil', 'section airfoil' )
 call point_prs%CreateRealOption(      'chord', 'section chord' )
 call point_prs%CreateRealOption(      'twist', 'section twist angle' )
 call point_prs%CreateRealArrayOption('SectionNormal', &
               'normal vector of the plane section containing the airfoil &
               &points', &
               '(/ 0.0, 1.0, 0.0 /)' ) ! default y-axis

 ! === Line sub-parser ===
 call pmesh_prs%CreateSubOption('Line','Line group',line_prs, &
               multiple=.true.)
 ! --- line group ---
 call line_prs%CreateStringOption(       'Type', &
               'type of the line connecting points: Straight or Spline' )
 call line_prs%CreateIntArrayOption('EndPoints', &
               'list of point id.s belonging to the line' )
 call line_prs%CreateIntOption(        'Nelems', &
               'n. spanwise elems of the line section' )
 call line_prs%CreateStringOption('type_span', 'type of span-wise division: &
               &uniform, cosine, cosineIB, cosineOB', &
               'uniform' ) ! defualt
 !> TangentVec1,2: in the code:
 ! straight lines: useles input -> tan vec computed
 ! spline lines  : either assigned or inherited from neighbouring lines
 !                        
 call line_prs%CreateRealArrayOption('TangentVec1' , &
               'assign tangent vector to the line at its first point' )
 call line_prs%CreateRealArrayOption('TangentVec2' , &
               'assign tangent vector to the line at its first point' )
                



end subroutine set_parser_pointwise

!----------------------------------------------------------------------
! ! === checks ===
! ! --- points, lines structures ---
! write(*,*) ' nPoints , nLines : ' , size(points) , size(lines)
! do i = 1 , size(points)
!   write(*,*)      points(i) % id    
!   write(*,*)      points(i) % coord
!   write(*,*) trim(points(i) % airfoil)
!   write(*,*)      points(i) % chord
!   write(*,*)      points(i) % theta
! end do
! do i = 1 , size(lines)
!   write(*,*) trim(lines(i) % l_type)
!   write(*,*)      lines(i) % end_points
!   write(*,*)      lines(i) % nelems
! end do

!  !> Build reference line
!  ! test ---
!  allocate(spl%rr(4,3)) ; allocate(spl%d0(3)) ; allocate(spl%d1(3))
!  spl%rr(1,:) = (/ 0.0_wp , 0.0_wp , 0.0_wp /) 
!  spl%rr(2,:) = (/ 1.0_wp , 1.0_wp , 0.0_wp /) 
!  spl%rr(3,:) = (/ 1.0_wp , 2.0_wp , 0.0_wp /) 
!  spl%rr(4,:) = (/ 0.0_wp , 3.0_wp , 0.0_wp /) 
!  spl%d0 = (/ 1.0_wp , 0.0_wp , 0.0_wp /)
!  spl%d1 = (/-1.0_wp , 0.0_wp , 0.0_wp /)
!  spl%e_bc(1) = 'derivative'
!  spl%e_bc(2) = 'derivative'
! 
!  call hermite_spline( spl , rr_spl ) 
! 
! ! check ---
! do i = 1 , size(rr_spl,1)
!   write(*,*) rr_spl(i,:)
! end do


end module mod_pointwise_io
