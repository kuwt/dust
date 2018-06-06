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


!> Module to handle the input/output of the solution of DUST in hdf5 format,
!! meant principally to exchange data with the post processor and restart
module mod_dust_io

use mod_param, only: &
  wp, max_char_len, nl

use mod_sim_param, only: &
  t_sim_param

use mod_handling, only: &
  error, warning, info, printout, dust_time, t_realtime, check_preproc

use mod_aero_elements, only: &
  c_elem, t_elem_p

use mod_surfpan, only: &
  t_surfpan

use mod_vortring, only: &
  t_vortring

use mod_liftlin, only: &
  t_liftlin

use mod_reference, only: &
  t_ref
  
use mod_hdf5_io, only: &
   h5loc, &
   new_hdf5_file, &
   open_hdf5_file, &
   close_hdf5_file, &
   new_hdf5_group, &
   open_hdf5_group, &
   close_hdf5_group, &
   write_hdf5, &
   write_hdf5_attr, &
   read_hdf5, &
   read_hdf5_al, &
   check_dset_hdf5

use mod_geometry, only: &
  t_geo, t_geo_component

use mod_wake_pan, only: &
  t_wake_panels

use mod_wake_ring, only: &
  t_wake_rings

use mod_stringtools, only: &
  stricmp
!----------------------------------------------------------------------

implicit none

public :: save_status, load_solution, load_time

private


character(len=*), parameter :: this_mod_name = 'mod_dust_io'

!----------------------------------------------------------------------

contains

!----------------------------------------------------------------------
subroutine save_status(geo, wake_pan, wake_rin, sim_params, it, time, run_id)
 type(t_geo), intent(in)         :: geo
 type(t_wake_panels), intent(in) :: wake_pan
 type(t_wake_rings), intent(in) :: wake_rin
 type(t_sim_param), intent(in)  :: sim_params
 integer, intent(in)             :: it
 real(wp), intent(in)            :: time
 integer, intent(in)             :: run_id(10)

 integer(h5loc) :: floc, gloc1, gloc2, gloc3, ploc
 character(len=max_char_len) :: comp_name, sit
 integer :: icomp, ncomp
 integer :: iref, nref
 character(len=max_char_len) :: ref_name
 integer :: ie, ne
 real(wp), allocatable :: vort(:), cp(:) , pres(:)
 real(wp), allocatable :: dforce(:,:) 
 real(wp), allocatable :: points_w(:,:,:), cent(:,:,:)
 integer, allocatable :: conn_pe(:)
 integer :: ir, id, ip, np

  !create the output file
  write(sit,'(I4.4)') it
  call new_hdf5_file(trim(sim_params%basename)//'_res_'//trim(sit)//'.h5', &
                     floc)
  call write_hdf5_attr(run_id, 'run_id', floc)
  call write_hdf5(time,'time',floc)

  call new_hdf5_group(floc, 'Parameters', ploc)
  call write_hdf5(sim_params%u_inf,'u_inf', ploc)
  call write_hdf5(sim_params%P_inf,'P_inf', ploc)
  call write_hdf5(sim_params%rho_inf,'rho_inf', ploc)
  call close_hdf5_group(ploc)
  

  ! 1) %%%%% Component solution: 
  ! mimic the structure of the components inside the geometry input file

  call new_hdf5_group(floc, 'Components', gloc1)
  ncomp = size(geo%components)
  call write_hdf5(ncomp, 'NComponents', gloc1)
  !Cycle the components
  do icomp = 1, ncomp
    write(comp_name,'(A,I3.3)')'Comp',icomp
    call new_hdf5_group(gloc1, trim(comp_name), gloc2)

    call write_hdf5(trim(geo%components(icomp)%comp_name),'CompName',gloc2)
    call write_hdf5(trim(geo%components(icomp)%ref_tag),'RefTag',gloc2)
    call write_hdf5(geo%components(icomp)%ref_id,'RefId',gloc2)
    
    call new_hdf5_group(gloc2, 'Solution', gloc3)

    ne = size(geo%components(icomp)%el)
    allocate(vort(ne), cp(ne) , pres(ne) , dforce(3,ne) )
    do ie = 1,ne
     vort(ie) = geo%components(icomp)%el(ie)%idou
     cp(ie) = geo%components(icomp)%el(ie)%cp
     pres(ie) = geo%components(icomp)%el(ie)%pres
     dforce(:,ie) = geo%components(icomp)%el(ie)%dforce
    enddo
    call write_hdf5(vort,'Vort',gloc3)
    call write_hdf5(cp,'Cp',gloc3)
    call write_hdf5(pres,'Pres',gloc3)
    call write_hdf5(dforce,'dF',gloc3)
    deallocate(vort, cp, pres, dforce)

    call close_hdf5_group(gloc3)
    call close_hdf5_group(gloc2)
  enddo
  call close_hdf5_group(gloc1)
  
  ! 2) %%%% Wake:
  ! just print the whole points, the solution and the starting row
  ! connectivity to build connectivity after
  call new_hdf5_group(floc, 'PanelWake', gloc1)

  call write_hdf5(wake_pan%w_points(:,:,1:wake_pan%wake_len+1),'WakePoints',gloc1 )
  call write_hdf5(wake_pan%i_start_points,'StartPoints',gloc1)
  call write_hdf5(wake_pan%ivort(:,1:wake_pan%wake_len),'WakeVort',gloc1)

  call close_hdf5_group(gloc1)

  call new_hdf5_group(floc, 'RingWake', gloc1)
  allocate(points_w(3,wake_rin%np_row,wake_rin%wake_len))
  allocate(conn_pe(wake_rin%np_row))
  allocate(cent(3,wake_rin%ndisks,wake_rin%wake_len))
  ip = 1
  do id = 1,wake_rin%ndisks
    np = wake_rin%gen_elems(id)%p%n_ver
    do ir = 1,wake_rin%wake_len
      points_w(:,ip:ip+np-1,ir) = wake_rin%wake_rings(id,ir)%ver(:,:)
      cent(:,id,ir) = wake_rin%wake_rings(id,ir)%cen(:)
    enddo
    conn_pe(ip:ip+np-1) = id
    ip = ip+np
  enddo
  call write_hdf5(points_w,'WakePoints',gloc1)
  call write_hdf5(conn_pe,'Conn_pe',gloc1)
  call write_hdf5(cent,'WakeCenters',gloc1)
  call write_hdf5(wake_rin%ivort(:,1:wake_rin%wake_len),'WakeVort',gloc1)
  call close_hdf5_group(gloc1)
  deallocate(points_w, conn_pe, cent)

  ! 3) %%%% References
  ! save the whole list of references
  call new_hdf5_group(floc, 'References', gloc1)
  nref = size(geo%refs)
  call write_hdf5(nref, 'NReferences', gloc1)
  do iref = 0,nref-1
    write(ref_name,'(A,I3.3)')'Ref',iref
    call new_hdf5_group(gloc1, trim(ref_name), gloc2)
    
    call write_hdf5(geo%refs(iref)%tag, 'Tag', gloc2)
    call write_hdf5(geo%refs(iref)%of_g, 'Offset',gloc2)
    call write_hdf5(geo%refs(iref)%R_g, 'R',gloc2)
    call write_hdf5(geo%refs(iref)%f_g, 'Vel',gloc2)
    call write_hdf5(geo%refs(iref)%G_g, 'RotVel',gloc2)

    call close_hdf5_group(gloc2)
  enddo

  call close_hdf5_group(gloc1)

  !close the output file
  call close_hdf5_file(floc)


end subroutine save_status

!----------------------------------------------------------------------
!> Load the solution from a pre-existent solution file
!!
!! Loads just the bodies solution, which is stored into the components,
!! the loading of the wakes is left to other subroutines
subroutine load_solution(filename,comps)
 character(len=max_char_len), intent(in) :: filename
 type(t_geo_component), intent(inout) :: comps(:)

 integer(h5loc) :: floc, gloc1, gloc2, gloc3
 integer :: ncomp, icomp, icomp2
 integer :: ne, ie
 character(len=max_char_len) :: comp_name_read, comp_id
 real(wp), allocatable :: idou(:), pres(:), dF(:,:)
 

 character(len=*), parameter :: this_sub_name = 'load_solution'

  call open_hdf5_file(filename, floc)

  !TODO: check the run id with the geo file
  call open_hdf5_group(floc, 'Components', gloc1)
  
  !Check if the number of components is consistent
  call read_hdf5(ncomp, 'NComponents', gloc1)
  if(ncomp .ne. size(comps)) call error(this_sub_name, this_mod_name, &
  'Mismatching number of components between geometry and solution files')

  !Cycle on all the components on the file and then on all the local
  !components: they come from different files and can be mixed up in the
  !order
  do icomp = 1, ncomp
    write(comp_id,'(A,I3.3)')'Comp',icomp
    call open_hdf5_group(gloc1, trim(comp_id), gloc2)

    call read_hdf5(comp_name_read,'CompName',gloc2)
    do icomp2 = 1,size(comps)
      if (stricmp(comp_name_read, comps(icomp2)%comp_name)) then
        !We found the correct local component, read all
        call open_hdf5_group(gloc2, 'Solution', gloc3)
        call read_hdf5_al(idou,'Vort',gloc3)
        call read_hdf5_al(pres,'Pres',gloc3)
        call read_hdf5_al(dF,'dF',gloc3)
        call close_hdf5_group(gloc3)
        ne = size(comps(icomp2)%el)
        do ie =1,ne
          comps(icomp2)%el(ie)%idou   = idou(ie)
          comps(icomp2)%el(ie)%pres   = pres(ie)
          comps(icomp2)%el(ie)%dforce = dF(:,ie)
        enddo
        deallocate(idou, pres, dF)
        exit
      endif

    enddo

    !If the component was not found there is a problem
    if(icomp2 .gt. size(comps)) call error(this_sub_name, this_mod_name, &
    'Component loaded from geometry file not found in result file')

    call close_hdf5_group(gloc2)
  enddo

  call close_hdf5_group(gloc1)
  call close_hdf5_file(floc)

end subroutine load_solution

!----------------------------------------------------------------------
subroutine load_time(filename, time)
 character(len=*), intent(in) :: filename
 real(wp), intent(out) :: time

 integer(h5loc) :: floc

  call open_hdf5_file(filename, floc)

  call read_hdf5(time,'time',floc)

  call close_hdf5_file(floc)
end subroutine load_time

!----------------------------------------------------------------------

end module mod_dust_io
