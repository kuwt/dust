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

!> Module containing the subroutines to perform flowfield visualizations
!! during postprocessing
module mod_post_flowfield

use mod_param, only: &
  wp, nl, max_char_len, extended_char_len , pi

use mod_handling, only: &
  error, warning

use mod_geometry, only: &
  t_geo, t_geo_component

use mod_parse, only: &
  t_parse, &
  getstr, getlogical, getintarray, getrealarray, &
  countoption

use mod_aeroel, only: &
  t_elem_p 


use mod_wake, only: &
  t_wake

use mod_hdf5_io, only: &
   h5loc, &
   open_hdf5_file, &
   close_hdf5_file, & 
   open_hdf5_group, &
   close_hdf5_group, &
   read_hdf5 

use mod_stringtools, only: &
  LowCase

use mod_geo_postpro, only: &
  load_components_postpro, update_points_postpro , &
  prepare_geometry_postpro  

use mod_tecplot_out, only: &
  tec_out_box 

use mod_vtk_out, only: &
  vtr_write  

use mod_post_load, only: &
  load_refs , load_res, load_wake_post

implicit none

public :: post_flowfield

private

character(len=max_char_len), parameter :: this_mod_name = 'mod_post_flowfield'

contains

! ---------------------------------------------------------------------- 

subroutine post_flowfield( sbprms , basename , data_basename , an_name , ia , &
                            out_frmt , components_names , all_comp , &
                            an_start , an_end , an_step )
type(t_parse), pointer :: sbprms
character(len=*) , intent(in) :: basename
character(len=*) , intent(in) :: data_basename
character(len=*) , intent(in) :: an_name
integer          , intent(in) :: ia
character(len=*) , intent(in) :: out_frmt
character(len=max_char_len), allocatable , intent(inout) :: components_names(:)
logical , intent(inout) :: all_comp
integer , intent(in) :: an_start , an_end , an_step

type(t_geo_component), allocatable  :: comps(:)
integer , parameter :: n_max_vars = 3 !vel,p,vort, ! TODO: 4 with cp
character(len=max_char_len), allocatable :: var_names(:)
integer , allocatable :: vars_n(:)
integer :: i_var , n_vars
logical :: probe_vel , probe_p , probe_vort 
real(wp) :: u_inf(3)
real(wp) :: P_inf , rho
real(wp) :: vel_probe(3) = 0.0_wp , vort_probe(3) = 0.0_wp
real(wp) :: pres_probe 
real(wp) :: v(3) = 0.0_wp

real(wp), allocatable , target :: sol(:) 
integer(h5loc) :: floc , ploc
real(wp), allocatable :: points(:,:)
integer :: nelem

real(wp), allocatable :: refs_R(:,:,:), refs_off(:,:)
real(wp), allocatable :: refs_G(:,:,:), refs_f(:,:)
real(wp), allocatable :: vort(:), cp(:)

type(t_wake)        :: wake
type(t_elem_p), allocatable :: wake_elems(:)

real(wp) :: t

integer :: nxyz(3)
real(wp):: minxyz(3) , maxxyz(3)
real(wp), allocatable :: xbox(:) , ybox(:) , zbox(:)
real(wp) :: dxbox , dybox , dzbox
real(wp), allocatable :: box_vel(:,:) , box_p(:) , box_vort(:,:)
integer :: ix , iy , iz
real(wp), allocatable :: vars(:,:) 
integer :: i_var_v , i_var_p , i_var_w

integer :: ip , ic , ie , i1 , it

character(len=max_char_len) :: str_a , var_name 
character(len=max_char_len) :: filename

character(len=max_char_len), parameter :: & 
   this_sub_name = 'post_flowfield'

    write(*,*) nl//' Analysis:',ia,' post_flowfield() +++++++++ '//nl

! Select all the components
!   components_names is allocated in load_components_postpro()
!   and deallocated in dust_post at the end of each analysis
if ( allocated(components_names) ) then
  call warning(trim(this_sub_name), trim(this_mod_name), &
     'All the components are used. <Components> input &
     &is ignored, and deallocated.' )
  deallocate(components_names)
end if 
all_comp = .true.

! Read variables to save : velocity | pressure | vorticity
! TODO: add Cp
probe_vel = .false. ; probe_p = .false. ; probe_vort = .false.
n_vars = countoption(sbprms,'Variable')

if ( n_vars .eq. 0 ) then ! default: velocity | pressure | vorticity
  probe_vel = .true. ; probe_p = .true. ; probe_vort = .true.
else
 do i_var = 1 , n_vars
  var_name = getstr(sbprms,'Variable')
  write(*,*) ' trim(var_name) : ' , trim(var_name) ; call LowCase(var_name)
  select case(trim(var_name))
   case ( 'velocity' ) ; probe_vel = .true.
   case ( 'pressure' ) ; probe_p   = .true.
   case ( 'vorticity') ; probe_vort= .true.
   case ( 'cp'       ) 
    write(str_a,*) ia 
    call error('dust_post','','Unknown Variable: '//trim(var_name)//&
               ' for analysis n.'//trim(str_a)//'.'//nl//&
                'Choose "velocity", "pressure", "vorticity".')
   case ( 'all') ; probe_vel = .true. ; probe_p   = .true. ; probe_vort= .true.
   case default
    write(str_a,*) ia 
    call error('dust_post','','Unknown Variable: '//trim(var_name)//&
               ' for analysis n.'//trim(str_a)//'.'//nl//&
                'Choose "velocity", "pressure", "vorticity".')
  end select
 end do
end if

! load the geo components just once
call open_hdf5_file(trim(data_basename)//'_geo.h5', floc)
!TODO: here get the run id    !todo????
call load_components_postpro(comps, points, nelem, floc, & 
                             components_names,  all_comp)
call close_hdf5_file(floc)

! Prepare_geometry_postpro
call prepare_geometry_postpro(comps)

! Allocate and point to sol
allocate(sol(nelem)) ; sol = 0.0_wp
ip = 0
do ic = 1 , size(comps)
 do ie = 1 , size(comps(ic)%el)
  ip = ip + 1
  comps(ic)%el(ie)%mag => sol(ip) 
 end do
end do

! Read box dimensions ...
nxyz   = getintarray( sbprms,'Nxyz'  ,3)
minxyz = getrealarray(sbprms,'Minxyz',3)
maxxyz = getrealarray(sbprms,'Maxxyz',3)

! ... allocate 'box' and deal with inconsistent input
if ( nxyz(1) .gt. 1 ) then ! x-coord
  allocate( xbox(nxyz(1)) )
  dxbox = ( maxxyz(1) - minxyz(1) ) / dble( nxyz(1) - 1 ) 
  xbox = (/ ( minxyz(1) + dble(i1-1) * dxbox , i1 = 1 , nxyz(1) )/) 
else
  allocate( xbox(1) )
  xbox(1) = minxyz(1)
end if
if ( nxyz(2) .gt. 1 ) then ! y-coord
  allocate( ybox(nxyz(2)) )
  dybox = ( maxxyz(2) - minxyz(2) ) / dble( nxyz(2) - 1 ) 
  ybox = (/ ( minxyz(2) + dble(i1-1) * dybox , i1 = 1 , nxyz(2) )/) 
else
  allocate( ybox(1) )
  ybox(1) = minxyz(2)
end if
if ( nxyz(3) .gt. 1 ) then ! z-coord
  allocate( zbox(nxyz(3)) )
  dzbox = ( maxxyz(3) - minxyz(3) ) / dble( nxyz(3) - 1 ) 
  zbox = (/ ( minxyz(3) + dble(i1-1) * dzbox , i1 = 1 , nxyz(3) )/) 
else
  allocate( zbox(1) )
  zbox(1) = minxyz(3)
end if

! Allocate box_vel, box_p, box_vort, ...
allocate(var_names(n_max_vars)) ; var_names = ' '
allocate(vars_n   (n_max_vars)) ; vars_n = 0
i_var = 0
i_var_v = 0 ; i_var_p = 0 ; i_var_w = 0
if ( probe_vel ) then
  allocate(box_vel (product(nxyz),3))
  i_var = i_var + 1
  var_names(i_var) = 'velocity'
  vars_n(i_var) = 3
  i_var_v = 3
end if
if ( probe_p   ) then
  allocate(box_p   (product(nxyz)  ))
  i_var = i_var + 1
  var_names(i_var) = 'pressure'
  vars_n(i_var) = 1
  i_var_p = i_var_v + 1
end if 
if ( probe_vort) then
  allocate(box_vort(product(nxyz),3))
  i_var = i_var + 1
  var_names(i_var) = 'vorticity'
  vars_n(i_var) = 3
  i_var_w = i_var_p + 3
end if
!TODO: cp
! if ( probe_cp  ) then
!   allocate(box_cp  (product(nxyz)  ))
!   i_var = i_var + 1
!   var_names(i_var) = 'cp'
!   vars_n(i_var) = 1
!   i_var_cp = i_var_w + 1
! end if

! Allocate and fill vars array, for output +++++++++++++++++
!  sum(vars_n): # of scalar fields to be plotted
!  product(nxyz) # of points where the vars are plotted 
allocate(vars(sum(vars_n),product(nxyz))) ; vars = 0.0_wp 

write(*,'(A,I0,A,I0,A,I0)') nl//' it_start,it_end,an_step : ' , &
  an_start , ' , ' , an_end , ' , ' , an_step
do it = an_start, an_end, an_step ! Time history

  ! Show timing, since this analysis is quite slow
  write(*,'(A,I0,A,I0)') ' it : ' , it , ' / ' , &
    ( an_end - an_start + 1 ) / an_step

  ! Open the result file ----------------------
  write(filename,'(A,I4.4,A)') trim(data_basename)// &
                                        '_res_',it,'.h5'
  call open_hdf5_file(trim(filename),floc)

  ! Load u_inf --------------------------------
  call open_hdf5_group(floc,'Parameters',ploc)
  call read_hdf5(u_inf,'u_inf',ploc)
  call read_hdf5(P_inf,'P_inf',ploc)
  call read_hdf5(rho,'rho_inf',ploc)
  call close_hdf5_group(ploc)

  ! Load the references and move the points ---
  call load_refs(floc,refs_R,refs_off,refs_G,refs_f)

  call update_points_postpro(comps, points, refs_R, refs_off, refs_G, refs_f)

  ! Load the results --------------------------
  call load_res(floc, comps, vort, cp, t)
  !sol = vort

  ! Load the wake -----------------------------
  call load_wake_post(floc, wake, wake_elems) 
  call close_hdf5_file(floc)

  ! Compute fields to be plotted +++++++++++++++++++++++++++++
  ip = 0
  ! Loop over the nodes of the box
  do iz = 1 , size(zbox)  ! z-coord
   do iy = 1 , size(ybox)  ! y-coord
    do ix = 1 , size(xbox)  ! x-coord
     ip = ip + 1

      if ( probe_vel .or. probe_p ) then 

        ! Compute velocity
        vel_probe = 0.0_wp ; pres_probe = 0.0_wp ; vort_probe = 0.0_wp

        ! body
        do ic = 1,size(comps) ! Loop on components
         do ie = 1 , size( comps(ic)%el ) ! Loop on elems of the comp
          call comps(ic)%el(ie)%compute_vel( (/ xbox(ix) , ybox(iy) , zbox(iz) /) , & 
                                              u_inf , v )
          vel_probe = vel_probe + v/(4*pi) 
         end do
        end do

        ! wake
        do ie = 1, size(wake_elems)
          call wake_elems(ie)%p%compute_vel( &
                   (/ xbox(ix) , ybox(iy) , zbox(iz) /) , &
                   u_inf , v )
          vel_probe = vel_probe + v/(4*pi) 
        enddo
       
        ! + u_inf
        vel_probe = vel_probe + u_inf
        
      end if

      if ( probe_vel ) then
        vars(1:3,ip) = vel_probe
      end if

      if ( probe_p ) then
        ! Bernoulli equation
        ! rho * dphi/dt + P + 0.5*rho*V^2 = P_infty + 0.5*rho*V_infty^2
        !TODO: add:
        ! - add the unsteady term: -rho*dphi/dt
        pres_probe = P_inf + 0.5_wp*rho*norm2(u_inf)**2 - 0.5_wp*rho*norm2(vel_probe)**2
        vars(i_var_v+1,ip) = pres_probe 
         
      end if

      if ( probe_vort ) then
        !TODO: add vorticity output
        ! ...
      end if

    end do  ! x-coord
   end do  ! y-coord
  end do  ! z-coord

  ! Output +++++++++++++++++++++++++++++++++++++++++++++++++

  select case (trim(out_frmt))

  case ('vtk')
    write(filename,'(A,I4.4,A)') trim(basename)//'_'//&
                             trim(an_name)//'_',it,'.vtr'
    call vtr_write ( filename , xbox , ybox , zbox , &
                     vars_n(1:i_var) , var_names(1:i_var) , &
                     vars ) 
  case('tecplot')
   write(filename,'(A,I4.4,A)') trim(basename)//'_'//&
                            trim(an_name)//'_',it,'.plt'
   i_var = 0
   deallocate(var_names)
   allocate(var_names(7))
   if(probe_vel) then
     var_names(i_var + 1) = 'ux'
     var_names(i_var + 2) = 'uy'
     var_names(i_var + 3) = 'uz'
     i_var = i_var + 3
   endif
   if(probe_p) then
     var_names(i_var + 1) = 'p'
     i_var = i_var + 1
   endif
   if(probe_vort) then
     var_names(i_var + 1) = 'omx'
     var_names(i_var + 2) = 'omy'
     var_names(i_var + 3) = 'omz'
     i_var = i_var + 3
   endif

   call tec_out_box(filename, t, xbox, ybox, zbox, &
                    vars, var_names(1:i_var))

  case default
    call error('dust_post','','Unknown format '//trim(out_frmt)//&
               ' for flowfield output. Choose: vtk or tecplot.')
  end select

end do  ! Time loop




! Nullify and deallocate
do ic = 1 , size(comps)
 do ie = 1 , size(comps(ic)%el)
  ip = ip + 1
  nullify(comps(ic)%el(ie)%mag) ! => null()  
 end do
end do
deallocate(sol) 

if ( allocated(box_vel ) ) deallocate(box_vel )
if ( allocated(box_p   ) ) deallocate(box_p   )
if ( allocated(box_vort) ) deallocate(box_vort)
!if ( allocated(box_cp  ) ) deallocate(box_cp  )
deallocate(xbox,ybox,zbox)
deallocate(var_names,vars_n)

!TODO: move deallocate(comps) outside this routine.
! Check if partial deallocation or nullification is needed.
deallocate(comps,components_names)

    write(*,*) nl//' post_flowfield done.'//nl

end subroutine post_flowfield

! ---------------------------------------------------------------------- 

end module mod_post_flowfield