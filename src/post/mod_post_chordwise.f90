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
!!          Alessandro Cocco
!!          Alberto Savino
!!=========================================================================

!> Module containing the subroutines to perform chordwise loads
!! analysis during postprocessing
module mod_post_chordwise

  use mod_param, only: &
    wp, nl, max_char_len, extended_char_len , pi

  use mod_handling, only: &
    error, internal_error, warning, info, printout, dust_time, t_realtime, &
    new_file_unit

  use mod_parse, only: &
    t_parse, &
    getrealarray, getreal, getint, getlogical, &
    getsuboption,  countoption

  use mod_hdf5_io, only: &
    initialize_hdf5, destroy_hdf5, &
    h5loc, &
    new_hdf5_file, &
    open_hdf5_file, &
    close_hdf5_file, &
    new_hdf5_group, &
    open_hdf5_group, &
    close_hdf5_group, &
    write_hdf5, &
    read_hdf5, &
    read_hdf5_al, &
    check_dset_hdf5

  use mod_geo_postpro, only: &
    load_components_postpro, update_points_postpro , prepare_geometry_postpro

  use mod_geometry, only: &
    t_geo, t_geo_component, destroy_elements

  use mod_post_load, only: &
    load_refs, load_res, load_ll, load_vl

  !use mod_dat_out, only: &
  !  dat_out_chordwise

  use mod_tecplot_out, only: &
    tec_out_sectional

  use mod_math, only: &
    cross, linear_interp 

  implicit none

  public :: post_chordwise

  private

  character(len=*), parameter :: this_mod_name = 'mod_post_chordwise'
  character(len=max_char_len) :: msg

  contains


!TODO (to be tested):
! - parametric components aligned with y-axis
! - tol_y_cen is hard-coded. write tol_y_cen as an input
! - preCICE coupled components 
subroutine post_chordwise(sbprms, basename, data_basename, an_name, &
                          ia, out_frmt, components_names, all_comp, &
                          an_start, an_end, an_step, average )
  type(t_parse), pointer                                  :: sbprms
  character(len=*) , intent(in)                           :: basename
  character(len=*) , intent(in)                           :: data_basename
  character(len=*) , intent(in)                           :: an_name
  integer          , intent(in)                           :: ia
  character(len=*) , intent(in)                           :: out_frmt
  character(len=max_char_len), allocatable, intent(inout) :: components_names(:)
  logical , intent(in)                                    :: all_comp
  integer , intent(in)                                    :: an_start , an_end , an_step
  logical , intent(in)                                    :: average

  type(t_geo_component), allocatable                      :: comps(:)
  character(len=max_char_len)                             :: cname
  integer(h5loc)                                          :: floc, gloc, cloc, ploc 
  real(wp), allocatable                                   :: refs_R(:,:,:), refs_off(:,:)
  real(wp), allocatable                                   :: vort(:), cp(:)
  real(wp), allocatable                                   :: points(:,:)
  integer                                                 :: nelem
  integer                                                 :: n_comp , n_comp_tot  
  integer                                                 :: i_comp , id_comp , ax_coor , ref_id
  character(len=max_char_len), allocatable                :: all_components_names(:)
  character(len=max_char_len), allocatable                :: components_names_tmp(:)
  integer                                                 :: nelem_span , nelem_chor , n_sect , n_time

  real(wp), allocatable                                   :: time(:)
  real(wp)                                                :: t
  integer                                                 :: ires

  !> chordwise load: paramteric components
  real(wp)                                                :: axis_dir(3) , axis_nod(3)
  integer                                                 :: is                                 
  integer, parameter                                      :: n_loads = 3   ! 
  real(wp), allocatable                                   :: y_cen(:) , y_span(:)
  real(wp), parameter                                     :: tol_y_cen = 1.0e-3_wp
  
  integer                                                 :: ie , ic , it , i1, ista, idir
  integer                                                 :: ipan_minus, ipan_plus 
  integer                                                 :: n_station 
  real(wp), allocatable                                   :: span_station(:)
  real(wp), allocatable                                   :: y_cen_tras(:) 
  integer, allocatable                                    :: id_minus(:), id_plus(:) 
  real(wp), allocatable                                   :: nor(:,:), tang(:,:) 
  real(wp)                                                :: u_inf(3), P_inf, rho 
  !> field to interpolate  
  real(wp), allocatable                                   :: force_minus(:),  force_plus(:) 
  real(wp), allocatable                                   :: tang_minus(:),   tang_plus(:) 
  real(wp), allocatable                                   :: nor_minus(:),    nor_plus(:) 
  real(wp), allocatable                                   :: cen_minus(:),    cen_plus(:) 
  real(wp)                                                :: pres_minus,     pres_plus 
  real(wp)                                                :: cp_minus,       cp_plus 
  !> interpolated field
  real(wp), allocatable                                   :: force_int(:,:,:,:)   
  real(wp), allocatable                                   :: tang_int(:,:,:,:)    
  real(wp), allocatable                                   :: nor_int(:,:,:,:)    
  real(wp), allocatable                                   :: cen_int(:,:,:,:) 
  real(wp), allocatable                                   :: pres_int(:,:,:)      
  real(wp), allocatable                                   :: cp_int(:,:,:)      

  character(len=max_char_len)                             :: filename
  character(len=max_char_len)                             :: comp_input

  character(len=*), parameter :: this_sub_name = 'post_chordwise' 

  write(msg,'(A,I0,A)') nl//'++++++++++ Analysis: ',ia,' chordwise loads'//nl
  call printout(trim(msg))


  ! Some warnings and errors -------------------------------
  ! WARNING: sectional loads are computed for one components only
  !  at time. If no component is specified --> return
  n_comp = countoption(sbprms,'component')
  if ( n_comp .le. 0 ) then
    call warning(this_mod_name, this_sub_name, 'No component specified for &
                &chordwise_loads analysis. Skipped analysis.')
    return
  else if ( n_comp .ge. 2 ) then
    call warning(this_mod_name, this_sub_name, &
        'More than one component specified &
      &for ''chordwise_loads'' analysis: just the first one is considered. &
      &Please run another ''chordwise_analysis'' if you need it on more than &
      &component.')
  end if

  !> check if components_names is allocated (it should alwasy be allocated)
  if ( .not. allocated(components_names) ) then
    call internal_error(this_mod_name,this_sub_name, &
                        'components_name not allocated.')
  end if

  !------------ Check if the component really exists ---------- !
  call open_hdf5_file(trim(data_basename)//'_geo.h5', floc)
  call open_hdf5_group(floc,'Components',gloc)
  call read_hdf5(n_comp_tot,'NComponents',gloc)

  allocate(all_components_names(n_comp_tot))
  do i_comp = 1 , n_comp_tot
    write(cname,'(A,I3.3)') 'Comp',i_comp
    call open_hdf5_group(gloc,trim(cname),cloc)
    call read_hdf5(all_components_names(i_comp),'CompName',cloc)    
    call close_hdf5_group(cloc)
  end do
  call close_hdf5_group(gloc)

  id_comp = -333
  do i_comp = 1 , n_comp_tot
      if ( trim(components_names(1)) .eq. trim(all_components_names(i_comp)) ) then
        id_comp = i_comp
      end if
  end do
  if ( id_comp .eq. -333 ) then ! component not found
    call printout('  Component not found. The available components are: ')
    do i_comp = 1 , n_comp_tot
      write(msg,'(A,I0,A)') '  comp. ' , i_comp , &
                          ':  '//trim(all_components_names(i_comp))
      call printout(trim(msg))
    end do
    call error(this_sub_name,this_mod_name,'No valid component found.')
  end if

  write(msg,*) '  Analysing component : ' , trim(components_names(1))//nl
  call printout(trim(msg))
  
  !> load only the first component just once 
  call open_hdf5_file(trim(data_basename)//'_geo.h5', floc)

  allocate(components_names_tmp(1))
  components_names_tmp(1) = trim(components_names(1))

  call load_components_postpro(comps, points, nelem, floc, &
                              components_names_tmp,  all_comp)
  call close_hdf5_file(floc)

  ! Prepare_geometry_postpro
  call prepare_geometry_postpro(comps)

  ! Read the axis for the computation of the sectional loads
  axis_dir = getrealarray(sbprms,'axis_dir',3)
  axis_nod = getrealarray(sbprms,'axis_nod',3)

  n_station = getint(sbprms,'n_station')
  span_station = getrealarray(sbprms,'span_station', n_station)
  comp_input = trim(comps(1)%comp_input ) 
  select case( trim(comp_input) )

  case ( 'parametric' )
    ! Some assumptions ---------------
    id_comp = 1   ! 1. only one component is loaded
    ax_coor = 2   ! 2. the parametric elements is defined
                  !    along the y-axis 

    nelem_span = comps(id_comp)%parametric_nelems_span
    nelem_chor = comps(id_comp)%parametric_nelems_chor
    n_sect = nelem_span

    do i1 = 1, size(comps(id_comp)%loc_points, 2) 
      comps(id_comp)%loc_points(:, i1) = matmul(comps(id_comp)%coupling_node_rot,comps(id_comp)%loc_points(:,i1))
    end do
    
    allocate(y_cen(n_sect),y_span(n_sect))
    ie = 0
    do is = 1 , n_sect
      ie = ie + 1
      y_cen(is) = &
              sum(comps(id_comp)%loc_points(ax_coor,comps(id_comp)%el(ie)%i_ver)) / &
              real(comps(id_comp)%el(ie)%n_ver,wp)
      do ic = 2 , nelem_chor ! check
        ie = ie + 1
        if (abs( y_cen(is) - &
          sum(comps(id_comp)%loc_points(ax_coor,comps(id_comp)%el(ie)%i_ver))&
            / real(comps(id_comp)%el(ie)%n_ver,wp) ) &
            .gt. tol_y_cen ) then
          call error(this_sub_name,this_mod_name,'Wrong section definition &
          &for component '//trim(comps(id_comp)%comp_name)//' for sectional&
          & loads analysis.')
        end if
      end do
    end do

    allocate(id_minus(n_station))
    allocate(id_plus(n_station))
    allocate(y_cen_tras(n_sect)) 
    do ista = 1, n_station
      !> translate centers on reference station  
      y_cen_tras = y_cen - span_station(ista) 
      do is = 1, n_sect - 1
        !> get index of the stripes across the span stripe 
        if (y_cen_tras(is) .le. 0.0_wp .and. y_cen_tras(is + 1) .ge. 0.0_wp)  then 
          id_minus(ista) = is
          id_plus(ista) = is + 1
        endif
      enddo     
    enddo  


    ! Find the coordinate of the reference points on the axis --------------
    !  ( with coord. y_cen )
    if ( abs(axis_dir(2)) .lt. 1e-6_wp ) then
      call error('dust_post','','Wrong definition of the axis in&
            & sectional_loads analysis: abs(axis_dir(2)) .lt. 1e-6.&
            & STOP')
    end if

    ! Find the ref_id or the reference frame where loads are projected -----
    ref_id = comps(id_comp)%ref_id
    write(msg,'(A,I0,A)') '   Employing the local reference frame for moment &
      &and force projection.'//nl//'   Reference tag: '&
      //trim(comps(id_comp)%ref_tag)//'  (ID: ',ref_id,')'
    call printout(trim(msg))

    n_time = (an_end-an_start)/an_step + 1 ! int general eger division
    allocate( time(n_time) ) ; time = -333.0_wp
    allocate( nor(n_time,3), tang(n_time,3) )

    ires = 0

    allocate(force_minus(3)); force_minus = 0.0_wp
    allocate(tang_minus(3));  tang_minus = 0.0_wp
    allocate(nor_minus(3));   nor_minus = 0.0_wp
    allocate(cen_minus(3));   cen_minus = 0.0_wp 
    pres_minus = 0.0_wp
    cp_minus = 0.0_wp 

    allocate(force_plus(3)); force_plus = 0.0_wp
    allocate(tang_plus(3));  tang_plus = 0.0_wp
    allocate(nor_plus(3));   nor_plus = 0.0_wp
    allocate(cen_plus(3));   cen_plus = 0.0_wp 
    pres_plus = 0.0_wp
    cp_plus = 0.0_wp 

    allocate(force_int(n_time,n_station,nelem_chor,3)); force_int = 0.0_wp 
    allocate(tang_int(n_time,n_station,nelem_chor,3));  tang_int = 0.0_wp 
    allocate(nor_int(n_time,n_station,nelem_chor,3));   nor_int = 0.0_wp 
    allocate(cen_int(n_time,n_station,nelem_chor,3));   cen_int = 0.0_wp 
    allocate(pres_int(n_time,n_station,nelem_chor));    pres_int = 0.0_wp 
    allocate(cp_int(n_time,n_station,nelem_chor));      cp_int = 0.0_wp 
    
    do it = an_start, an_end, an_step
      ires = ires + 1

      !> Open the file:
      write(filename,'(A,I4.4,A)') trim(data_basename)//'_res_',it,'.h5'
      call open_hdf5_file(trim(filename),floc)
      !> Load free-stream parameters
      call open_hdf5_group(floc,'Parameters',ploc)
      call read_hdf5(u_inf,'u_inf',ploc)
      call read_hdf5(P_inf,'P_inf',ploc)
      call read_hdf5(rho,'rho_inf',ploc)
      call close_hdf5_group(ploc) 

      ! Load the references and move the points ---
      call load_refs(floc, refs_R, refs_off)
      ! Move the points ---------------------------
      call update_points_postpro(comps, points, refs_R, refs_off, &
                                  filen = trim(filename) )
      ! Load the results --------------------------
      call load_res(floc, comps, vort, cp, t)
      
      call close_hdf5_file(floc)

      
      ! compute chordwise loads. Loop over the panels of each section
      do ista = 1, n_station 
        ipan_minus = 0
        ipan_plus = 0
        ipan_minus = (id_minus(ista)- 1)*nelem_chor
        ipan_plus = (id_plus(ista)- 1)*nelem_chor 

        do ic = 1, nelem_chor ! loop over chord
          ipan_minus = ipan_minus + 1  
          !> 4d matrix: time, station, chord, [x, y, z] 
          force_minus = comps(id_comp)%el(ipan_minus)%dforce
          tang_minus = comps(id_comp)%el(ipan_minus)%tang(:,1)
          nor_minus = comps(id_comp)%el(ipan_minus)%nor
          cen_minus = comps(id_comp)%el(ipan_minus)%cen 
          !> 3d matrix: time, station, chord(pres) 
          pres_minus = comps(id_comp)%el(ipan_minus)%pres 
          cp_minus = (pres_minus - P_inf)/&
                                    (0.5_wp*rho*norm2(u_inf)**2) 

          ipan_plus = ipan_plus + 1 
          !> 4d matrix: time, station, chord, [x, y, z] 
          force_plus = comps(id_comp)%el(ipan_plus)%dforce
          tang_plus = comps(id_comp)%el(ipan_plus)%tang(:,1)
          nor_plus = comps(id_comp)%el(ipan_plus)%nor
          !> 3d matrix: time, station, chord(pres) 
          pres_plus = comps(id_comp)%el(ipan_plus)%pres 
          cp_plus  = (pres_plus - P_inf)/& 
                                    (0.5_wp*rho*norm2(u_inf)**2) 

          ! From global to local coordinates of forces and moments
          force_minus = matmul(transpose(refs_R(:,:, ref_id)), &
                                                  force_minus) & 
                                                  / y_span(id_minus(ista)) 
          tang_minus = matmul(transpose( refs_R(:,:, ref_id)), &
                                                  tang_minus)                                                   
          nor_minus = matmul(transpose( refs_R(:,:, ref_id)), &
                                                  nor_minus)                                                   
          cen_minus = matmul(transpose(refs_R(:,:, ref_id)), &
                                                  cen_minus) 

          force_plus = matmul(transpose(refs_R(:,:, ref_id)), &
                                                  force_plus) & 
                                                  / y_span(id_plus(ista)) 
          tang_plus = matmul(transpose( refs_R(:,:, ref_id)), &
                                                  tang_plus)                                                   
          nor_plus = matmul(transpose( refs_R(:,:, ref_id)), &
                                                  nor_plus)                                                   
          cen_plus = matmul(transpose(refs_R(:,:, ref_id)), &
                                                  cen_plus)           
          !> interpolate  
          do idir = 1,3
            call linear_interp((/force_minus(idir), force_plus(idir)/), & 
                                (/y_cen(id_minus(ista)), y_cen(id_plus(ista))/),& 
                                span_station(ista), force_int(ires,ista,ic,idir)) 
            call linear_interp((/tang_minus(idir), tang_plus(idir)/), & 
                                (/y_cen(id_minus(ista)), y_cen(id_plus(ista))/),& 
                                span_station(ista), tang_int(ires,ista,ic,idir)) 
            call linear_interp((/nor_minus(idir), nor_plus(idir)/), & 
                                (/y_cen(id_minus(ista)), y_cen(id_plus(ista))/),& 
                                span_station(ista), nor_int(ires,ista,ic,idir)) 
            call linear_interp((/cen_minus(idir), cen_plus(idir)/), & 
                                (/y_cen(id_minus(ista)), y_cen(id_plus(ista))/),& 
                                span_station(ista), cen_int(ires,ista,ic,idir)) 
          enddo 
          
          call linear_interp((/pres_minus, pres_plus/), & 
                            (/y_cen(id_minus(ista)), y_cen(id_plus(ista))/),& 
                            span_station(ista), pres_int(ires,ista,ic)) 
          call linear_interp((/cp_minus, cp_plus/), & 
                            (/y_cen(id_minus(ista)), y_cen(id_plus(ista))/),& 
                            span_station(ista), cp_int(ires,ista,ic)) 
          
          
        end do ! loop over chord
        
          
        

      end do ! loop over sections

      time(ires) = t

    end do 

    if(average) then 
      !> todo
    else 
      select case(trim(out_frmt))
      case('dat')
        write(filename,'(A)') trim(basename)//'_'//trim(an_name)
      case('tecplot')
        write(filename,'(A)') trim(basename)//'_'//trim(an_name)//'.plt' 
      end select 
    endif 


  case('pointwise')
    call error(this_sub_name,this_mod_name, 'case pointwise not implemented so far')
  case('cgns') 
    call error(this_sub_name,this_mod_name, 'case cgns not implemented so far')
  case default
    call error(this_sub_name,this_mod_name, 'type not known')   
  end select 

end subroutine 
end module 
