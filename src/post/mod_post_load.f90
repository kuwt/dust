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
   check_dset_hdf5, &
   get_dset_dimensions_hdf5

use mod_actuatordisk, only: &
  t_actdisk

use mod_wake, only: &
  t_wake

use mod_aeroel, only: &
  t_elem_p

implicit none

public :: load_refs , load_res ,  &
          load_wake_post, load_wake_viz , &
          load_wake_pan , load_wake_ring , &
          check_if_components_exist

private

character(len=*), parameter :: this_mod_name = 'mod_post_load'

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
 integer :: ncomps_sol

 character(len=*), parameter :: this_sub_name = 'load_res'

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
  call read_hdf5(ncomps_sol,'NComponents',gloc1)
  if(ncomps_sol .lt. ncomps) call error(this_sub_name, this_mod_name, &
  'Different number of components between solution and geometry')

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

    !check consistency of geometry and solution
    if((size(vort_read,1) .ne. nelems_comp) .or. &
       (size(pres_read,1) .ne. nelems_comp) .or. &
       (size(dforce_read,2) .ne. nelems_comp)) call error(this_mod_name, &
       this_sub_name, 'inconsistent number of elements between geometry and&
       & solution')   

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

!----------------------------------------------------------------------

!> Load the wake for the postprocessing. 
!! do everything 
subroutine load_wake_post(floc, wake, wake_p)
 integer(h5loc), intent(in) :: floc 
 type(t_wake), target, intent(out)  :: wake
 type(t_elem_p), allocatable, intent(out) :: wake_p(:)

 integer(h5loc) :: gloc
 real(wp), allocatable :: wpoints_pan(:,:,:)
 integer,  allocatable :: wstart_pan(:,:)
 real(wp), allocatable :: wvort_pan(:,:)
 real(wp), allocatable :: wpoints_rin(:,:,:)
 integer,  allocatable :: wconn_rin(:)
 real(wp), allocatable :: wvort_rin(:,:)
 real(wp), allocatable :: vppoints(:,:), vpvort(:,:)
 integer :: n_wake_stripes , npan, ndisks, nrows
 integer :: nsides
 integer :: p1 , p2 
 integer :: ip , iw, id, ir, iconn, i
 integer :: npt_disk
 integer, allocatable :: disk_pts(: )

  !=== Panels ===
  call open_hdf5_group(floc,'PanelWake',gloc)
  call read_hdf5_al(wpoints_pan,'WakePoints',gloc)
  call read_hdf5_al(wstart_pan,'StartPoints',gloc)
  call read_hdf5_al(wvort_pan,'WakeVort',gloc)
  call close_hdf5_group(gloc)

  n_wake_stripes = size(wstart_pan ,2)
  npan           = size(wvort_pan  ,2) 

  wake%nmax_pan = npan
  wake%n_pan_stripes = n_wake_stripes
  wake%pan_wake_len = npan
  wake%n_pan_points = size(wpoints_pan,2)
  allocate(wake%i_start_points(2,wake%n_pan_stripes))
  allocate(wake%pan_w_points(3,wake%n_pan_points,npan+1))
  allocate(wake%wake_panels(wake%n_pan_stripes,npan))
  allocate(wake%pan_idou(wake%n_pan_stripes,npan))

  wake%i_start_points     = wstart_pan
  wake%pan_w_points       = wpoints_pan

  nsides = 4 
  do ip = 1,npan
    do iw=1,wake%n_pan_stripes
     wake%wake_panels(iw,ip)%mag => wake%pan_idou(iw,ip)
     wake%wake_panels(iw,ip)%n_ver = nsides
     allocate(wake%wake_panels(iw,ip)%ver(3,nsides))
     allocate(wake%wake_panels(iw,ip)%edge_vec(3,nsides))
     allocate(wake%wake_panels(iw,ip)%edge_len(nsides))
     allocate(wake%wake_panels(iw,ip)%edge_uni(3,nsides))
    enddo
  enddo

  do ip = 1,wake%pan_wake_len
   do iw = 1,wake%n_pan_stripes
       p1 = wake%i_start_points(1,iw)
       p2 = wake%i_start_points(2,iw)
       call wake%wake_panels(iw,ip)%calc_geo_data( &
       reshape((/wake%pan_w_points(:,p1,ip),   wake%pan_w_points(:,p2,ip), &
                 wake%pan_w_points(:,p2,ip+1), wake%pan_w_points(:,p1,ip+1)/),&
                                                                   (/3,4/)))
       wake%wake_panels(iw,ip)%mag = wvort_pan(iw,ip) 
   end do
  end do

  !=== Rings ===
  call open_hdf5_group(floc,'RingWake',gloc)
  call read_hdf5_al(wpoints_rin,'WakePoints',gloc)
  call read_hdf5_al(wconn_rin,'Conn_pe',gloc)
  call read_hdf5_al(wvort_rin,'WakeVort',gloc)
  call close_hdf5_group(gloc)

  ndisks = size(wvort_rin,1)
  nrows = size(wvort_rin,2)
  wake%ndisks = ndisks; wake%rin_wake_len = nrows
  allocate(wake%wake_rings(wake%ndisks,wake%rin_wake_len))
  allocate(wake%rin_idou(wake%ndisks,wake%rin_wake_len))

  do id = 1,wake%ndisks
    !reverse the connectivity, from pts2disk to disk2pts
    npt_disk = count(wconn_rin .eq. id)
    allocate(disk_pts(npt_disk))
    iconn = 1
    do ip = 1,size(wconn_rin)
      if(wconn_rin(ip) .eq. id) then
        disk_pts(iconn) = ip
        iconn = iconn+1
      endif
    enddo

    nsides = npt_disk
    do ir = 1,wake%rin_wake_len
      wake%wake_rings(id,ir)%mag => wake%rin_idou(id,ir)
      wake%wake_rings(id,ir)%n_ver = nsides
      allocate(wake%wake_rings(id,ir)%ver(3,nsides))
      allocate(wake%wake_rings(id,ir)%edge_vec(3,nsides))
      allocate(wake%wake_rings(id,ir)%edge_len(nsides))
      allocate(wake%wake_rings(id,ir)%edge_uni(3,nsides))
      call wake%wake_rings(id,ir)%calc_geo_data(wpoints_rin(:,disk_pts,ir))
      wake%wake_rings(id,ir)%mag = wvort_rin(id,ir)
      
    enddo
    deallocate(disk_pts)
  enddo

  !=== Particles + Line vortex ===
  call open_hdf5_group(floc, 'ParticleWake', gloc)
  call read_hdf5_al(vppoints,'WakePoints',gloc)
  call read_hdf5_al(vpvort,'WakeVort',gloc)
  call close_hdf5_group(gloc)

  wake%n_prt = size(vpvort,2)
  wake%nmax_prt = size(vpvort,2)

  allocate(wake%wake_parts(wake%nmax_prt))
  allocate(wake%prt_ivort(wake%nmax_prt))
  allocate(wake%part_p(wake%n_prt))
  if(wake%n_prt .gt. 0) then
    allocate(wake%vort_p(wake%n_prt+wake%n_pan_stripes))
    allocate(wake%end_vorts(wake%n_pan_stripes))
    do iw = 1, wake%n_pan_stripes 
      wake%vort_p(wake%n_prt+iw)%p => wake%end_vorts(iw)
      wake%end_vorts(iw)%mag => wake%wake_panels(iw,wake%pan_wake_len)%mag
      p1 = wake%i_start_points(1,iw)
      p2 = wake%i_start_points(2,iw)
      call wake%end_vorts(iw)%calc_geo_data( &
          reshape((/wake%pan_w_points(:,p1,wake%pan_wake_len+1),  &
                    wake%pan_w_points(:,p2,wake%pan_wake_len+1)/), (/3,2/)))
    enddo
  else
    allocate(wake%vort_p(wake%n_prt))
  endif

  do ip = 1,wake%n_prt
    wake%wake_parts(ip)%cen = vppoints(:,ip)
    wake%wake_parts(ip)%mag => wake%prt_ivort(ip)
    wake%wake_parts(ip)%mag = norm2(vpvort(:,ip))
    if(wake%wake_parts(ip)%mag .gt. 1.0e-13_wp) then
      wake%wake_parts(ip)%dir = vpvort(:,ip)/wake%wake_parts(ip)%mag
    else
      wake%wake_parts(ip)%dir = vpvort(:,ip)
    endif
    wake%wake_parts(ip)%free = .false.
    wake%part_p(ip)%p => wake%wake_parts(ip)
    wake%vort_p(ip)%p => wake%wake_parts(ip)
  enddo
  
  deallocate(vppoints, vpvort)

  !Stitch everything together
  if(wake%n_prt .gt. 0) then
    allocate(wake_p(wake%n_pan_stripes*wake%pan_wake_len &
    + wake%ndisks*wake%rin_wake_len + &
      wake%n_pan_stripes+ wake%n_prt))
  else
    allocate(wake_p(wake%n_pan_stripes*wake%pan_wake_len &
    + wake%ndisks*wake%rin_wake_len))
  endif
 

  i=0
  !panels
  do ip = 1,wake%pan_wake_len
    do iw=1,wake%n_pan_stripes
    i = i+1
    wake_p(i)%p => wake%wake_panels(iw,ip)
    enddo
  enddo
  !rings
  do ir = 1,wake%rin_wake_len
    do id = 1,wake%ndisks
    i = i+1
    wake_p(i)%p => wake%wake_rings(id,ir)
    enddo
  enddo

  if(wake%n_prt .gt. 0) then
    !end vortices
    do iw=1,wake%n_pan_stripes
      i = i+1
      wake_p(i)%p => wake%end_vorts(iw)
    enddo
    !particles
    do ip = 1,wake%n_prt
      i = i+1
      wake_p(i)%p => wake%wake_parts(ip)
    enddo
  endif

end subroutine

!----------------------------------------------------------------------

!TODO: include the possibility of defining multiple components as an input
subroutine check_if_components_exist( components_list , filename )
 character(len=*) , intent(in) :: components_list(:)
 character(len=*) , intent(in) :: filename

 character(len=max_char_len) , allocatable :: components(:)
 character(len=max_char_len) :: cname !, msg
 integer(h5loc) :: floc , gloc , cloc
 integer :: n_comp_tot , n_comp_inp

 integer :: i1 , i2 , i_check , i_comp

 character(len=*), parameter :: this_sub_name = 'check_if_components_exist'

 n_comp_inp = size(components_list) ; write(*,*) ' n_comp_inp : ' , n_comp_inp
 !debug
  do i_comp = 1 , n_comp_inp
    write(*,*) 'comp',i_comp,':',components_list(i_comp)
  end do

 write(*,*) nl//' ++++ Check component input : '
 write(*,*) ' checking file: ', trim(filename)
 call open_hdf5_file(trim(filename),floc)
 call open_hdf5_group(floc,'Components',gloc)
 call read_hdf5(n_comp_tot,'NComponents',gloc)

 allocate(components(n_comp_tot))
 do i_comp = 1 , n_comp_tot
   write(cname,'(A,I3.3)') 'Comp',i_comp
   call open_hdf5_group(gloc,trim(cname),cloc)
   call read_hdf5(components(i_comp),'CompName',cloc)
   call close_hdf5_group(cloc)
 end do
 call close_hdf5_group(gloc)
 call close_hdf5_file(floc)

 do i1 = 1 , n_comp_inp
   i_check = 0 ; write(*,*) ' i_check : ' , i_check 
   do i2 = 1 , n_comp_tot
     if ( trim(components(i2)) .eq. trim(components_list(i1)) ) &
                                                        i_check = 1  
   end do
   if ( i_check .eq. 0 ) then
     write(*,*) ' All the available SINGLE components are: '
     do i2 = 1 , n_comp_tot
       write(*,*) i2 , ' : ' , trim(components(i2))
     end do
     call error(this_sub_name, this_mod_name, &
              'Component '//trim(components_list(i1))//' received as &
             &an input does not exist.') 
   end if 
 end do 

 write(*,*) nl//' ++++ Check component input : OK. '//nl 

end subroutine check_if_components_exist

!----------------------------------------------------------------------

end module mod_post_load