module mod_post_load

use mod_param, only: &
  wp, nl, max_char_len, pi

use mod_handling, only: &
  error, warning, info, printout, dust_time, t_realtime, new_file_unit

use mod_geometry, only: &
  t_geo, t_geo_component

use mod_hdf5_io, only: &
   h5loc, &
   open_hdf5_file, &
   close_hdf5_file, &
   open_hdf5_group, &
   close_hdf5_group, &
   read_hdf5, &
   read_hdf5_al, &
   check_dset_hdf5

use mod_actuatordisk, only: &
  t_actdisk

implicit none

public :: load_refs , load_res , load_wake_viz , load_wake_pan , load_wake_ring

private

contains

!----------------------------------------------------------------------

subroutine load_refs(floc, refs_R, refs_off, refs_G, refs_f, refs_tag)
 integer(h5loc), intent(in) :: floc 
 real(wp), allocatable, intent(out) :: refs_R(:,:,:)
 real(wp), allocatable, intent(out) :: refs_off(:,:)
 real(wp), allocatable, intent(out) , optional :: refs_G(:,:,:)
 real(wp), allocatable, intent(out) , optional :: refs_f(:,:)
 character(len=max_char_len) , allocatable , intent(out) , optional :: refs_tag(:)

 integer(h5loc) :: gloc1, gloc2
 integer :: nrefs, iref
 character(len=max_char_len) :: rname

  call open_hdf5_group(floc,'References',gloc1)
  call read_hdf5(nrefs,'NReferences',gloc1)

  allocate(refs_R(3,3,0:nrefs-1), refs_off(3,0:nrefs-1))
  if (present(refs_G)  ) allocate(refs_G(3,3,0:nrefs-1))
  if (present(refs_f)  ) allocate(refs_f(3,0:nrefs-1))
  if (present(refs_tag)) allocate(refs_tag(0:nrefs-1))
  do iref = 0,nrefs-1
    write(rname,'(A,I3.3)') 'Ref',iref
    call open_hdf5_group(gloc1,trim(rname),gloc2)
   
    call read_hdf5(refs_R(:,:,iref),'R',gloc2)
    call read_hdf5(refs_off(:,iref),'Offset',gloc2)
    if (present(refs_tag)) call read_hdf5(refs_tag(  iref),'Tag',gloc2)
    if (present(refs_G  )) call read_hdf5(refs_G(:,:,iref),'RotVel',gloc2)
    if (present(refs_f  )) call read_hdf5(refs_f(:,iref),'Vel',gloc2)
    call read_hdf5(refs_off(:,iref),'Offset',gloc2)

    call close_hdf5_group(gloc2)
  enddo

  call close_hdf5_group(gloc1)

end subroutine load_refs

!----------------------------------------------------------------------

subroutine load_res(floc, comps, vort, press, t)
 integer(h5loc), intent(in) :: floc 
 type(t_geo_component), intent(inout) :: comps(:)
 real(wp), allocatable, intent(out) :: vort(:)
 real(wp), allocatable, intent(out) :: press(:)
 real(wp), intent(out) :: t

 integer :: ncomps, icomp, ie
 integer :: nelems, offset, nelems_comp
 integer(h5loc) :: gloc1, gloc2, gloc3
 character(len=max_char_len) :: cname
 real(wp), allocatable :: vort_read(:)!, cp_read(:)
 real(wp), allocatable :: pres_read(:), dforce_read(:,:)

  ncomps = size(comps)
  nelems = 0
  do icomp = 1, ncomps
    select type(el=>comps(icomp)%el)
     class default
      nelems = nelems + comps(icomp)%nelems 
     type is(t_actdisk)
      nelems = nelems + size(comps(icomp)%loc_points,2)
    end select
  enddo
  
  call read_hdf5(t,'time',floc)
  allocate(vort(nelems), press(nelems))
  call open_hdf5_group(floc,'Components',gloc1)

  offset = 0
  do icomp = 1, ncomps
    
    nelems_comp = comps(icomp)%nelems
    write(cname,'(A,I3.3)') 'Comp',comps(icomp)%comp_id
    call open_hdf5_group(gloc1,trim(cname),gloc2)
    call open_hdf5_group(gloc2,'Solution',gloc3)

    !call read_hdf5(vort(offset+1:offset+nelems_comp),'Vort',gloc3)
    !call read_hdf5(cp(offset+1:offset+nelems_comp),'Cp',gloc3)
    call read_hdf5_al(vort_read,'Vort',gloc3)
    !call read_hdf5_al(cp_read,'Cp',gloc3)
    call read_hdf5_al(pres_read,'Pres',gloc3)
    call read_hdf5_al(dforce_read,'dF',gloc3)

!   TODO: check if it is general enough *******
!   TODO: check if something is broken after changing intent(in to inout) for comps

    do ie = 1 , nelems_comp
      comps(icomp)%el(ie)%pres = pres_read(ie) 
      comps(icomp)%el(ie)%dforce = dforce_read(:,ie) 
    end do

    call close_hdf5_group(gloc3)
    call close_hdf5_group(gloc2)

    select type(el =>comps(icomp)%el)
     class default
      vort(offset+1:offset+nelems_comp) = vort_read
      press(offset+1:offset+nelems_comp) = pres_read
      offset = offset + nelems_comp
      do ie = 1,nelems_comp
        if(associated(comps(icomp)%el(ie)%mag)) &
                        comps(icomp)%el(ie)%mag = vort_read(ie)
      enddo
     type is(t_actdisk)
      do ie = 1,nelems_comp
        vort(offset+1:offset+el(ie)%n_ver) = vort_read(ie)
        press(offset+1:offset+el(ie)%n_ver) = pres_read(ie)
        offset = offset + el(ie)%n_ver
        if(associated(comps(icomp)%el(ie)%mag)) then
                        comps(icomp)%el(ie)%mag = vort_read(ie)
        endif
      enddo
    end select

    deallocate(vort_read)!, cp_read)

  enddo

  call close_hdf5_group(gloc1)

end subroutine load_res

!----------------------------------------------------------------------

subroutine load_wake_viz(floc, wpoints, welems, wvort, vppoints,  vpvort)
 integer(h5loc), intent(in) :: floc 
 real(wp), allocatable, intent(out) :: wpoints(:,:)
 integer, allocatable, intent(out)  :: welems(:,:)
 real(wp), allocatable, intent(out) :: wvort(:)
 real(wp), allocatable, intent(out) :: vppoints(:,:)
 real(wp), allocatable, intent(out) :: vpvort(:)

 integer(h5loc) :: gloc
 logical :: got_dset 
 real(wp), allocatable :: wpoints_read(:,:,:)
 real(wp), allocatable :: wpoints_pan(:,:), wpoints_rin(:,:)
 integer, allocatable  :: wstart(:,:), wconn(:)
 real(wp), allocatable :: wcen(:,:,:)
 real(wp), allocatable :: wvort_read(:,:)
 real(wp), allocatable :: wvort_pan(:), wvort_rin(:)
 integer, allocatable  :: welems_pan(:,:), welems_rin(:,:)
 integer :: nstripes, npoints_row, nrows, ndisks, nelem_w
 integer :: iew, ir, is, ip
 integer :: first_elem, act_disk, next_elem


  !get the panel wake
  got_dset = check_dset_hdf5('PanelWake',floc)
  if(got_dset) then
 
    call open_hdf5_group(floc,'PanelWake',gloc)
    call read_hdf5_al(wpoints_read,'WakePoints',gloc)
    call read_hdf5_al(wstart,'StartPoints',gloc)
    call read_hdf5_al(wvort_read,'WakeVort',gloc)
 
 
    nstripes = size(wvort_read,1); nrows = size(wvort_read,2);
    npoints_row = size(wpoints_read,2)
    nelem_w = nstripes*nrows
 
    allocate(wpoints_pan(3,npoints_row*(nrows+1)))
    allocate(welems_pan(4,nelem_w))
    allocate(wvort_pan(nelem_w))
 
    wpoints_pan = reshape(wpoints_read, (/3,npoints_row*(nrows+1)/))
    wvort_pan = reshape(wvort_read, (/nelem_w/))
    iew = 0
    do ir = 1,nrows
      do is = 1,nstripes
        iew = iew+1
        welems_pan(:,iew) = (/wstart(1,is)+npoints_row*(ir-1), &
                              wstart(2,is)+npoints_row*(ir-1), &
                              wstart(2,is)+npoints_row*(ir), &
                              wstart(1,is)+npoints_row*(ir)/)
      enddo
    enddo
   
    deallocate(wpoints_read, wstart, wvort_read)
    call close_hdf5_group(gloc)
  else
    !panel wake not present, allocate stuff at zero size
    allocate(wvort_pan(0), wpoints_pan(3,0), welems_pan(4,0))
  endif
 
  !get the ring wake
  got_dset = check_dset_hdf5('RingWake',floc)
  if(got_dset) then
 
    call open_hdf5_group(floc,'RingWake',gloc)
    call read_hdf5_al(wpoints_read,'WakePoints',gloc)
    call read_hdf5_al(wconn,'Conn_pe',gloc)
    call read_hdf5_al(wcen,'WakeCenters',gloc)
    call read_hdf5_al(wvort_read,'WakeVort',gloc)
    
    ndisks = size(wvort_read,1); nrows = size(wvort_read,2)
    npoints_row = size(wpoints_read,2)
 
    nelem_w = ndisks*nrows
 
    allocate(wpoints_rin(3,npoints_row*nrows+nelem_w))
    allocate(welems_rin(4,npoints_row*nrows))
    allocate(wvort_rin(npoints_row*nrows))
 
    wpoints_rin(:,1:npoints_row*nrows) = reshape(wpoints_read, &
                                            (/3,npoints_row*nrows/))
    wpoints_rin(:,npoints_row*nrows+1:size(wpoints_rin,2)) = &
           reshape(wcen, (/3,nelem_w/))
 
    iew = 0; act_disk = 0
    do ir = 1,nrows
      do ip = 1, npoints_row
        iew = iew+1
 
        if(ip .eq. 1) then
          first_elem = iew
        else
          if(wconn(ip-1) .ne. wconn(ip)) first_elem = iew
        endif
 
        if(ip.lt.npoints_row) then
          if(wconn(ip+1) .ne. wconn(ip)) then
            next_elem = first_elem
          else
            next_elem = iew + 1
          endif
        else
          next_elem = first_elem
        endif
 
        welems_rin(1,iew) = iew
        welems_rin(2,iew) = next_elem
        !welems_rin(3,iew) = npoints_row*nrows+(ir-1)*nelem_w+wconn(ip)
        welems_rin(3,iew) = npoints_row*nrows+(ir-1)*ndisks+wconn(ip)
        welems_rin(4,iew) = 0
 
        wvort_rin(iew) = wvort_read(wconn(ip),ir)
      enddo
    enddo
 
    call close_hdf5_group(gloc)
    deallocate(wpoints_read, wconn, wcen, wvort_read) 
  else
    allocate(wvort_rin(0), wpoints_rin(3,0), welems_rin(4,0))
  endif
 
  !Stitch together the two wakes
  allocate(wpoints(3,size(wpoints_pan,2)+size(wpoints_rin,2)))
  wpoints(:,1:size(wpoints_pan,2)) = wpoints_pan
  wpoints(:,size(wpoints_pan,2)+1:size(wpoints,2)) = wpoints_rin
  deallocate(wpoints_pan, wpoints_rin)
 
  allocate(welems(4,size(welems_pan,2)+size(welems_rin,2)))
  welems(:,1:size(welems_pan,2)) = welems_pan
  welems(1:3,size(welems_pan,2)+1:size(welems,2)) = welems_rin(1:3,:) + &
                                                    size(wpoints_pan,2)
  welems(4,size(welems_pan,2)+1:size(welems,2)) = 0
  deallocate(welems_pan, welems_rin)
 
 
  allocate(wvort(size(wvort_pan)+size(wvort_rin)))
  wvort(1:size(wvort_pan)) = wvort_pan
  wvort(size(wvort_pan)+1:size(wvort)) = wvort_rin
  deallocate(wvort_pan, wvort_rin)
  
  got_dset = check_dset_hdf5('ParticleWake',floc)
  if(got_dset) then
 
    call open_hdf5_group(floc,'ParticleWake',gloc)
    call read_hdf5_al(vppoints,'WakePoints',gloc)
    call read_hdf5_al(wvort_read,'WakeVort',gloc)
    allocate(vpvort(size(wvort_read,2)))
    do ip = 1,size(vpvort)
      vpvort(ip) = norm2(wvort_read(:,ip))
    enddo
    deallocate(wvort_read)
    call close_hdf5_group(gloc)
  endif

end subroutine load_wake_viz

!----------------------------------------------------------------------

subroutine load_wake_pan(floc, wpoints, wstart, wvort)
 integer(h5loc), intent(in) :: floc 
 real(wp), allocatable, intent(out) :: wpoints(:,:,:)
 integer, allocatable, intent(out) :: wstart(:,:)
 real(wp), allocatable, intent(out) :: wvort(:,:)

 integer(h5loc) :: gloc
  
  call open_hdf5_group(floc,'PanelWake',gloc)
  
  call read_hdf5_al(wpoints,'WakePoints',gloc)
  call read_hdf5_al(wstart,'StartPoints',gloc)
  call read_hdf5_al(wvort,'WakeVort',gloc)

  call close_hdf5_group(gloc)

end subroutine load_wake_pan

!----------------------------------------------------------------------

subroutine load_wake_ring(floc, wpoints, wconn, wvort)
 integer(h5loc), intent(in) :: floc 
 real(wp), allocatable, intent(out) :: wpoints(:,:,:)
 integer, allocatable, intent(out) :: wconn(:)
 !real(wp), allocatable, intent(out) :: wcen(:,:,:)
 real(wp), allocatable, intent(out) :: wvort(:,:)

 integer(h5loc) :: gloc
  
  call open_hdf5_group(floc,'RingWake',gloc)
  
  call read_hdf5_al(wpoints,'WakePoints',gloc)
  call read_hdf5_al(wconn,'Conn_pe',gloc)
!  call read_hdf5_al(wcen,'WakeCenters',gloc)
  call read_hdf5_al(wvort,'WakeVort',gloc)

  call close_hdf5_group(gloc)

end subroutine load_wake_ring

end module mod_post_load
