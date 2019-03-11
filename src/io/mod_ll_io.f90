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


!> Module to treat input-output for lifting line components
!!
module mod_ll_io

use mod_param, only: &
  wp, max_char_len, nl, pi

use mod_handling, only: &
  error, warning, info, printout, new_file_unit, check_file_exists

use mod_parse, only: &
  t_parse, getstr, getint, getreal, getrealarray, getlogical, countoption

!----------------------------------------------------------------------

implicit none

public :: read_mesh_ll

private

character(len=*), parameter :: this_mod_name = 'mod_parametric_io'

!----------------------------------------------------------------------

contains

!----------------------------------------------------------------------

subroutine read_mesh_ll(mesh_file,ee,rr, &
                        airfoil_list_actual, nelem_span_list, &
                        i_airfoil_e , normalised_coord_e    , &
                        npoints_chord_tot , nelem_span_tot  , &
                                            chord_p,theta_p)

 character(len=*), intent(in) :: mesh_file
 integer  , allocatable, intent(out) :: ee(:,:) 
 real(wp) , allocatable, intent(out) :: rr(:,:) 
 character(len=max_char_len), allocatable , intent(out) :: airfoil_list_actual(:)
 character(len=max_char_len), allocatable               :: airfoil_list(:)
 integer  , allocatable, intent(out) :: nelem_span_list(:)
 integer  , allocatable, intent(out) :: i_airfoil_e(:,:)
 real(wp) , allocatable, intent(out) :: normalised_coord_e(:,:)
 integer  , intent(out), optional    :: npoints_chord_tot, nelem_span_tot
 real(wp) , allocatable, intent(out), optional :: chord_p(:),theta_p(:)

 type(t_parse) :: pmesh_prs

 integer :: nelem_chord, nelem_chord_tot ! , nelem_span_tot <--- moved as an output
 integer :: npoint_chord_tot, npoint_span_tot
 integer :: nRegions, nSections, nAirfoils, rr_size , ee_size 
 integer :: iRegion , iSection , iAirfoil , Span , iChord , iElement , iPoint 
 real(wp):: ref_chord_fraction
 real(wp), allocatable :: ref_point(:)
 ! data read from file
 ! sections ---
 real(wp)         , allocatable :: chord_list(:) , twist_list(:) 
!character(len=max_str_len), allocatable , intent(out) :: airfoil_list(:)
 ! regions  ---
!integer , allocatable , intent(out) :: nelem_span_list(:)
 real(wp), allocatable :: span_list(:) , sweep_list(:) , dihed_list(:)
 character(len=max_char_len), allocatable :: type_span_list(:)
 integer :: n_type_span
 character :: ElType
 logical :: symmetry
 ! Sections 1. 2.
 real(wp), allocatable :: rrSection1(:,:) , rrSection2(:,:) 
!real(wp) :: th1 , th2 , ch1 , ch2 
 real(wp) :: dx_ref , dy_ref , dz_ref
 integer :: ista , iend , ich

 integer :: i_aero1, i_aero2, iSpan, i

 character(len=*), parameter :: this_sub_name = 'read_mesh_parametric'

 
  !Prepare all the parameters to be read in the file
  ! Global parameters
  call pmesh_prs%CreateStringOption('ElType', &
                'element type (temporary) p panel v vortex ring')
  call pmesh_prs%CreateLogicalOption('mesh_symmetry', &
                'symmetry yes/no','F' )
! !!! USELESS !!! for LIFTING LINE components !!!!!!!!!!!!
! call pmesh_prs%CreateIntOption('nelem_chord',&
!               'number of chord-wise elements', &
!               multiple=.false.);
! call pmesh_prs%CreateStringOption('type_chord',&
!               'type of chord-wise division: uniform, cosine, cosineLE, cosineTE',&
!               'uniform', &
!               multiple=.false.);
! !!! USELESS !!! for LIFTING LINE components !!!!!!!!!!!!
  call pmesh_prs%CreateRealArrayOption('starting_point',&
               'Starting point (inboard TE), (x, y, z)', &
               '(/0.0, 0.0, 0.0/)',&
               multiple=.false.);
  call pmesh_prs%CreateRealOption('reference_chord_fraction',&
               'Reference chord fraction', &
               '0.0',&
               multiple=.false.);
  ! TODO: do we really need to allow the user to define the parameter above?
  ! the reference chord fraction is a number in the range [0,1] that defines
  ! - where the reference point is located in the root chord
  ! - the line related to the sweep angle
  ! - the point around which the sections are rotated to impose twist


  ! Section parameters
! already there (few lines above) as a characteristic of the whole component !!!
! call pmesh_prs%CreateRealOption('ref_point_chord', 'reference point of the section (as a fraction of the chord)',&
!               multiple=.true.);
  call pmesh_prs%CreateRealOption(     'chord', 'section chord', &
                multiple=.true.);
  call pmesh_prs%CreateRealOption(     'twist', 'section twist angle',&
                multiple=.true.);
  call pmesh_prs%CreateStringOption( 'airfoil', 'section airfoil',&
                multiple=.true.);

  ! Region parameters
  call pmesh_prs%CreateRealOption(    'span', 'region span',&
                multiple=.true.);
  call pmesh_prs%CreateRealOption(   'sweep', 'region sweep angle [degrees]',&
                multiple=.true.);
  call pmesh_prs%CreateRealOption(   'dihed', 'region dihedral angle [degrees]',&
                multiple=.true.);
  call pmesh_prs%CreateIntOption( 'nelem_span', 'number of span-wise elements in the region',&
                multiple=.true.);
  call pmesh_prs%CreateStringOption('type_span', 'type of span-wise division: uniform, cosine, cosineIB, cosineOB', &
                multiple=.true.);


  !read the parameters
  call pmesh_prs%read_options(mesh_file,printout_val=.true.)

  nelem_chord     = 1          !  getint(pmesh_prs,'nelem_chord')   !!! USELESS !!!
  nelem_chord_tot = 1          !  getint(pmesh_prs,'nelem_chord')   !!! USELESS !!!
  ElType  = getstr(pmesh_prs,'ElType')
  symmetry= getlogical(pmesh_prs,'mesh_symmetry')
   
  nSections = countoption(pmesh_prs,'chord')
  nRegions  = countoption(pmesh_prs,'span')

  !DEBUG
  write(*,*) ' nSections : ' , nSections
  write(*,*) ' nRegions  : ' , nRegions 
  
  ! Check that nSections = nRegion + 1
  if ( nSections .ne. nRegions + 1 ) then
    call error(this_sub_name, this_mod_name, 'Unconsistent input: &
         &nSections .ne. nRegions. Stop.')
  end if 

  ref_chord_fraction = getreal(pmesh_prs,'reference_chord_fraction')
  ref_point          = getrealarray(pmesh_prs,'starting_point',3)

  ! Get total number of elements and initialize arrays
  allocate(nelem_span_list(nRegions));   nelem_span_list = 0
  allocate(      span_list(nRegions));         span_list = 0.0_wp
  allocate(     sweep_list(nRegions));        sweep_list = 0.0_wp
  allocate(     dihed_list(nRegions));        dihed_list = 0.0_wp
  nelem_span_tot = 0
  do iRegion = 1,nRegions
    nelem_span_list(iRegion) = getint(pmesh_prs,'nelem_span')
    span_list(iRegion)  = getreal(pmesh_prs,'span' )
    sweep_list(iRegion) = getreal(pmesh_prs,'sweep')
    dihed_list(iRegion) = getreal(pmesh_prs,'dihed')
    nelem_span_tot = nelem_span_tot + nelem_span_list(iRegion)
  enddo

  ! type_span
  allocate( type_span_list(nRegions));!   type_span_list = ''
  n_type_span = countoption(pmesh_prs,'type_span')
  if ( n_type_span .eq. 0 ) then ! default is 'uniform'
    do iRegion = 1 , nRegions
      type_span_list(iRegion) = 'uniform'
    end do
  else if ( n_type_span .eq. nRegions ) then
    do iRegion = 1 , nRegions
      type_span_list(iRegion) = getstr(pmesh_prs,'type_span')
    end do
  else
   write(*,*) ' mesh_file   : ' , trim(mesh_file)
   write(*,*) ' n_type_span : ' , n_type_span
   write(*,*) ' nRegions    : ' , nRegions    
   call error(this_sub_name, this_mod_name, 'Unconsistent input: &
         &n_type_span .ne. nRegions. Stop.')
  end if

  allocate(chord_list  (nSections))  ; chord_list = 0.0d0
  allocate(twist_list  (nSections))  ; twist_list = 0.0d0
  allocate(airfoil_list(nSections)) 

  do iSection= 1,nSections
    chord_list(iSection)   = getreal(pmesh_prs,'chord')
    twist_list(iSection)   = getreal(pmesh_prs,'twist')
    airfoil_list(iSection) = getstr(pmesh_prs,'airfoil')
  enddo


  npoint_chord_tot = nelem_chord_tot + 1
  npoint_span_tot  = nelem_span_tot  + 1

  ee_size = nelem_span_tot
  rr_size = npoint_span_tot * npoint_chord_tot

  ! Directly build connectivity matrix
  allocate(ee(4,ee_size)) ; ee = 0
  do iChord = 1,nelem_chord_tot
    do iSpan = 1,nelem_span_tot
      iElement =  nelem_chord_tot*(iSpan-1) + iChord
      iPoint   = npoint_chord_tot*(iSpan-1) + iChord
      ee(1,iElement) = iPoint + npoint_chord_tot         ! iPoint
      ee(2,iElement) = iPoint                            ! iPoint + 1
      ee(3,iElement) = iPoint + 1                        ! iPoint + npoint_chord_tot + 1
      ee(4,iElement) = iPoint + npoint_chord_tot + 1     ! iPoint + npoint_chord_tot
    enddo
  enddo

  allocate(rr   (3,rr_size)) ; rr      = 0.0_wp
  allocate(chord_p(npoint_span_tot)) ; chord_p = 0.0_wp
  allocate(theta_p(npoint_span_tot)) ; theta_p = 0.0_wp

  allocate(i_airfoil_e       (2,nelem_span_tot)) ; i_airfoil_e        = 0
  allocate(normalised_coord_e(2,nelem_span_tot)) ; normalised_coord_e = 0.0_wp
! old
! allocate(airfoil_table_p (2,npoint_span_tot)) ;
! allocate(normalised_coord_p(npoint_span_tot)) ; normalised_coord_p = 0.0_wp

! ! Initialize the span division to the maximum dimension
  allocate(rrSection1(3,npoint_chord_tot)) ; rrSection1 = 0.0_wp
  allocate(rrSection2(3,npoint_chord_tot)) ; rrSection2 = 0.0_wp
! allocate(rrSection1(3,nRegions)) ; rrSection1 = 0.0_wp
! allocate(rrSection2(3,nRegions)) ; rrSection2 = 0.0_wp

  ! Initialise dr_ref for the definition of the actual reference_point 
  !  of each bay (by updating)
  dx_ref = 0.0_wp
  dy_ref = 0.0_wp
  dz_ref = 0.0_wp


  twist_list = twist_list * pi / 180
  ista = 1 ; iend = npoint_chord_tot 
  ich = 1

  do iRegion = 1,nRegions
  
! check ----
    write(*,*) ' Region ' , iRegion , ' / ' , nRegions
! check ----
 
    if ( iRegion .gt. 1 ) then  ! first section = last section of the previous region 
      rrSection1 = rrSection2
    else                        ! build points
      rrSection1(:,1) = ref_point
      rrSection1(:,2) = rrSection1(:,1) + chord_list(iRegion) * &
                (/ cos(twist_list(iRegion)) , 0.0_wp , -sin(twist_list(iRegion)) /)
      chord_p(ich) = chord_list(1)
      theta_p(ich) = twist_list(1)

! The following lines should not have existed even before 2019-02-06:
! fields overwritten at the end of the outer do loop
! === old-2019-02-06 === !
!      i_airfoil_e(1,ich) = iRegion 
!      i_airfoil_e(2,ich) = iRegion+1
! === old-2019-02-06 === !

      ! Update rr
      rr(:,ista:iend) = rrSection1

    end if
  
    dx_ref =  span_list(iRegion) * tan( sweep_list(iRegion)* pi / 180.0_wp ) ! + dx_ref 
    dy_ref =  span_list(iRegion)                                             ! + dy_ref 
    dz_ref =  span_list(iRegion) * tan( dihed_list(iRegion)* pi / 180.0_wp ) ! + dz_ref 
  
    rrSection2(:,1) = rrSection1(:,1) + (/ dx_ref , dy_ref , dz_ref /)
    rrSection2(:,2) = rrSection2(:,1) + chord_list(iRegion+1) * &
              (/ cos(twist_list(iRegion+1)) , 0.0_wp , -sin(twist_list(iRegion+1)) /)
 
    ! Interpolation of the nodes of the region i (between sections i and i+1)
    do iSpan = 1 , nelem_span_list(iRegion)
      ista = iend + 1 
      iend = iend + npoint_chord_tot
      ich = ich + 1

      
      if ( trim(type_span_list(iRegion)) .eq. 'uniform' ) then    
        ! uniform spacing in span
        rr(:,ista:iend) = rrSection1 + dble(iSpan) / dble(nelem_span_list(iRegion)) * &
                   ( rrSection2 - rrSection1 )
        chord_p(ich) = chord_list(iRegion) + dble(iSpan) / dble(nelem_span_list(iRegion)) * &
                   ( chord_list(iRegion+1) - chord_list(iRegion) )
        theta_p(ich) = twist_list(iRegion) + dble(iSpan) / dble(nelem_span_list(iRegion)) * &
                   ( twist_list(iRegion+1) - twist_list(iRegion) )

      else if ( trim(type_span_list(iRegion)) .eq. 'cosine' ) then    
        ! cosine  spacing in span
        rr(:,ista:iend) = 0.5_wp * ( rrSection1 + rrSection2 ) - &
                          0.5_wp * ( rrSection2 - rrSection1 ) * &
                                     cos(iSpan*pi/ dble(nelem_span_list(iRegion)) ) 
        chord_p(ich) = 0.5_wp * ( chord_list(iRegion+1) + chord_list(iRegion) ) - &
                       0.5_wp * ( chord_list(iRegion+1) - chord_list(iRegion) ) * &
                                     cos(iSpan*pi/ dble(nelem_span_list(iRegion)) ) 
        theta_p(ich) = 0.5_wp * ( twist_list(iRegion+1) + twist_list(iRegion) ) - &
                       0.5_wp * ( twist_list(iRegion+1) - twist_list(iRegion) ) * &
                                     cos(iSpan*pi/ dble(nelem_span_list(iRegion)) ) 
      else if ( trim(type_span_list(iRegion)) .eq. 'cosineOB' ) then    
        ! cosine  spacing in span: outboard refinement
        rr(:,ista:iend) = rrSection1 + &
                        ( rrSection2 - rrSection1 ) * &
                                     sin(0.5_wp*iSpan*pi/ dble(nelem_span_list(iRegion)) ) 
        chord_p(ich) =          (                         chord_list(iRegion) ) + &
                                ( chord_list(iRegion+1) - chord_list(iRegion) ) * &
                                     sin(0.5_wp*iSpan*pi/ dble(nelem_span_list(iRegion)) ) 
        theta_p(ich) =          (                         twist_list(iRegion) ) - &
                                ( twist_list(iRegion+1) - twist_list(iRegion) ) * &
                                     sin(0.5_wp*iSpan*pi/ dble(nelem_span_list(iRegion)) ) 
      else if ( trim(type_span_list(iRegion)) .eq. 'cosineIB' ) then    
        ! cosine  spacing in span: inboard refinement
        rr(:,ista:iend) = rrSection2 - &
                        ( rrSection2 - rrSection1 ) * &
                                     cos(0.5_wp*iSpan*pi/ dble(nelem_span_list(iRegion)) ) 
        chord_p(ich) =          (                         chord_list(iRegion+1))- &
                                ( chord_list(iRegion+1) - chord_list(iRegion) ) * &
                                     cos(0.5_wp*iSpan*pi/ dble(nelem_span_list(iRegion)) ) 
        theta_p(ich) =          (                         twist_list(iRegion+1))- &
                                ( twist_list(iRegion+1) - twist_list(iRegion) ) * &
                                     cos(0.5_wp*iSpan*pi/ dble(nelem_span_list(iRegion)) ) 
      else
        write(*,*) ' mesh_file   : ' , trim(mesh_file)
        write(*,*) ' type_span_list(',iRegion,') : ' , trim(type_span_list(iRegion)) 
        call error(this_sub_name, this_mod_name, 'Unconsistent input: &
              & type_span must be equal to uniform, cosine, cosineIB, cosineOB.')
      end if 

! ===
! Moved outside this loop (see below), to allow the specification of aerodynamic
! characteristic on some sections only
! ===
!       ! airfoil file and non-dimensional coordinate between two sections
!       i_airfoil_e( 1,ich-1) = iRegion 
!       i_airfoil_e( 2,ich-1) = iRegion+1
! !     normalised_coord_e(2,ich-1) = ( rr(2,ista) - rrSection1(2,iRegion) ) / &
! !                        ( rrSection2(2,iRegion) - rrSection1(2,iRegion) )
!       normalised_coord_e(2,ich-1) = ( rr(2,ista) - rrSection1(2,1) ) / &
!                          ( rrSection2(2,1) - rrSection1(2,1) )
!       if ( iSpan .ne. 1 ) then ! else = 0.0_wp
!         normalised_coord_e(1,ich-1) = normalised_coord_e(2,ich-2)
!       end if
! ===

    end do
  
  end do


! === new-2019-02-06 ===
! Save the fields related to the definition of actual airfoils
! airfoil_list(i) may be equal to 'interp' --> no def. of airfoil on the i-th section
!                                              --> inteprolation
  ! Check airfoil_list()
  if ( trim(airfoil_list(1)) .eq. 'interp' ) then
    call error(this_sub_name, this_mod_name, 'The first "airfoil"&
          & cannot be set as "interp".')
  end if
  if ( trim(airfoil_list(nSections)) .eq. 'interp' ) then
    call error(this_sub_name, this_mod_name, 'The last "airfoil"&
          & cannot be set as "interp".')
  end if

  nAirfoils = 0
  do i = 1 , nSections
    if ( trim(airfoil_list(i)) .ne. 'interp' ) nAirfoils = nAirfoils + 1 
  end do

  allocate(airfoil_list_actual(nAirfoils))
  iAirfoil = 0
  do i = 1 , nSections
    if ( trim(airfoil_list(i)) .ne. 'interp' ) then
      iAirfoil = iAirfoil + 1
      airfoil_list_actual(iAirfoil) = trim(airfoil_list(i))
      call check_file_exists(airfoil_list_actual(iAirfoil), this_sub_name, &
           this_mod_name)
    end if
  end do

  i_aero2 = 1 ; ista = 1 ; iAirfoil = 0
  do iAirfoil = 1 , nAirfoils - 1

    i_aero1 = i_aero2
    do i = i_aero1 + 1 , nSections
      if ( airfoil_list(i) .ne. 'interp' ) then
        i_aero2 = i
        exit
      end if
    end do

    iend = sum(nelem_span_list(1:i_aero2-1))

    do iSpan = ista , iend 
      i_airfoil_e(1, iSpan) = iAirfoil 
      i_airfoil_e(2, iSpan) = iAirfoil+1 

      normalised_coord_e(2,iSpan) = ( rr(2,(iSpan+1)*2) - rr(2,ista*2) ) / &
                                    ( rr(2,iend *2+1) - rr(2,ista*2) )
      if ( iSpan .ne. ista ) then ! else = 0.0_wp
        normalised_coord_e(1,iSpan) = normalised_coord_e(2,iSpan-1)
      end if

    end do

    ista = iend + 1

  end do

! === new-2019-02-06 ===


  ! optional output ----
  npoints_chord_tot = npoint_chord_tot
  ! optional output ----


end subroutine read_mesh_ll

!----------------------------------------------------------------------

end module mod_ll_io
