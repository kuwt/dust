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
  t_geo

use mod_wake_pan, only: &
  t_wake_panels

use mod_wake_ring, only: &
  t_wake_rings
!----------------------------------------------------------------------

implicit none

public :: save_status

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
 real(wp), allocatable :: vort(:), cp(:)
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
    allocate(vort(ne), cp(ne))
    do ie = 1,ne
     vort(ie) = geo%components(icomp)%el(ie)%idou
     cp(ie) = geo%components(icomp)%el(ie)%cp
    enddo
    call write_hdf5(vort,'Vort',gloc3)
    call write_hdf5(cp,'Cp',gloc3)
    deallocate(vort, cp)

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

end module mod_dust_io
