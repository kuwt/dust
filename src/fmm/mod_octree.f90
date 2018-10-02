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


!> Module to handle the octree grid
module mod_octree

use mod_param, only: &
  wp, nl, pi, max_char_len

use mod_sim_param, only: &
  t_sim_param

use mod_handling, only: &
  error, warning, info, printout, dust_time, t_realtime

use mod_vortpart, only: &
  t_vortpart, t_vortpart_p

use mod_multipole, only: &
  t_multipole, t_polyexp, t_ker_der

!----------------------------------------------------------------------

implicit none

public :: initialize_octree, sort_particles, t_octree, perform_multipole

private

!> Pointer to a cell
type :: t_cell_p
  type(t_cell), pointer :: p => null()
end type

!----------------------------------------------------------------------

!> Type containing all the data relative to a single octree cell
type :: t_cell

  logical :: leaf

  logical :: branch

  logical :: active

  type(t_cell_p) :: parent

  type(t_cell_p) :: children(8)

  type(t_cell_p) :: neighbours(-1:1,-1:1,-1:1)

  type(t_cell_p), allocatable :: interaction_list(:)

  integer :: cart_index(3)

  integer :: level

  integer :: npart

  type(t_vortpart_p), allocatable :: cell_parts(:)

  real(wp) :: cen(3)

  type(t_multipole) :: mp

end type

!----------------------------------------------------------------------

!> Type containing the pointers to a layer in the octree
type :: t_cell_layer

  real(wp) :: cell_size
  
  type(t_cell), allocatable :: lcells(:,:,:)

  type(t_ker_der), allocatable :: ker_der(:,:,:)

end type


!----------------------------------------------------------------------

!> Type containing the whole octree structure
type :: t_octree
 
 !> maximum and minimum of the coordinates of the whole octree
 real(wp) :: xmin(3), xmax(3)

 !> number of level 1 cubic boxes in each direction
 integer :: nbox(3)

 !> number of levels of octree division
 integer :: nlevels
 
 !> total number of cells (sum of cells at each level)
 integer :: ncells_tot

 !> all the actual cells, in a 1d array
 !type(t_cell), allocatable :: cells(:)
 
 !> cell pointers organized in layers, each layer corresponding to 
 !! a level of the octree
 type(t_cell_layer), allocatable :: layers(:)
 
 !> pointer array to all the leaves
 type(t_cell_p), allocatable :: leaves(:)

 !> number of leaves
 integer :: nleaves

 !> Degree of multipole expansion
 integer :: degree

 !> Polynomial expansion tools
 type(t_polyexp) :: pexp
 type(t_polyexp) :: pexp_der

 !> Regularized kernel radius
 real(wp) :: delta

end type

!----------------------------------------------------------------------
! Interfaces
interface push_ptr
  module procedure push_cell_ptr_scalar, push_cell_ptr_vec
  module procedure push_part_ptr_scalar, push_part_ptr_vec
end interface

!----------------------------------------------------------------------
character(len=*), parameter :: this_mod_name='mod_octree'
character(len=max_char_len) :: msg
real(t_realtime) :: t1 , t0

integer :: min_part_4_cell

!----------------------------------------------------------------------
contains
!----------------------------------------------------------------------

!> Initialize the octree
subroutine initialize_octree(box_length, nbox, origin, nlevels, min_part, &
                             degree, delta, octree)
 real(wp), intent(in) :: box_length
 integer,  intent(in) :: nbox(3)
 real(wp), intent(in) :: origin(3)
 integer,  intent(in) :: nlevels
 integer,  intent(in) :: min_part
 integer,  intent(in) :: degree
 real(wp), intent(in) :: delta
 type(t_octree), intent(out), target :: octree

 integer :: l, i,j,k, ic,jc,kc
 integer :: imax, jmax, kmax
 integer :: c, p, child
 integer :: indx(3)

  octree%xmin = origin
  octree%xmax = origin + real(nbox,wp)*box_length
  octree%nbox = nbox
  octree%nlevels = nlevels
  octree%delta = delta

  min_part_4_cell = min_part

  octree%degree = degree
  call octree%pexp%set_degree(degree)
  call octree%pexp_der%set_degree(2*degree)

  octree%ncells_tot = product(nbox)
  do l = 2,nlevels
    octree%ncells_tot = octree%ncells_tot + product(nbox)*8**(l-1)
  enddo

  allocate(octree%leaves(product(nbox)*8**(nlevels-1)))
  !allocate(octree%cells(octree%ncells_tot))
  !allocate(octree%leaves(0))

  !build the hierarchical tree
  allocate(octree%layers(nlevels))

  ! == Prepare the hierarchical tree
  !first level
  !c = 1 !start the counter for the global cells array
  octree%layers(1)%cell_size = box_length

  allocate(octree%layers(1)%lcells(nbox(1),nbox(2),nbox(3)))
  do k = 1,nbox(3)
    do j = 1,nbox(2)
      do i = 1,nbox(1)
        !octree%layers(1)%lcells(i,j,k)%p => octree%cells(c)
        !octree%cells(c)%cart_index = (/i,j,k/)
        !octree%cells(c)%level = 1
        octree%layers(1)%lcells(i,j,k)%cart_index = (/i,j,k/)
        octree%layers(1)%lcells(i,j,k)%level = 1


        octree%layers(1)%lcells(i,j,k)%cen = octree%xmin + &
                    (real((/i,j,k/),wp)-0.5_wp)*(octree%xmax-octree%xmin) &
                    /real(nbox,wp)
   !     c = c+1
      enddo
    enddo
  enddo
  !dive down the levels
  do l = 2,nlevels
    
    octree%layers(l)%cell_size = box_length/(2**(l-1))
    allocate(octree%layers(l)%lcells(nbox(1)*2**(l-1),nbox(2)*2**(l-1),nbox(3)*2**(l-1)))
    !cycle on the parents
    do k = 1,nbox(3)*2**(l-2)
    do j = 1,nbox(2)*2**(l-2)
    do i = 1,nbox(1)*2**(l-2)

      p = 1 !reset the counter for the parent to set the children
      do kc = 1,2
      do jc = 1,2
      do ic = 1,2
        !point the cell in the ordered layer
        !octree%layers(l)%lcells(2*(i-1)+ic,2*(j-1)+jc,2*(k-1)+kc)%p => octree%cells(c)
        !octree%cells(c)%cart_index = (/2*(i-1)+ic, 2*(j-1)+jc, 2*(k-1)+kc/)
        !octree%cells(c)%level = l
        octree%layers(l)%lcells(2*(i-1)+ic,2*(j-1)+jc,2*(k-1)+kc)%cart_index &
                             = (/2*(i-1)+ic, 2*(j-1)+jc, 2*(k-1)+kc/)        
        octree%layers(l)%lcells(2*(i-1)+ic,2*(j-1)+jc,2*(k-1)+kc)%level = l

        !set the parent
        octree%layers(l)%lcells(2*(i-1)+ic,2*(j-1)+jc,2*(k-1)+kc)%parent%p &
                                     => octree%layers(l-1)%lcells(i,j,k) 
        octree%layers(l)%lcells(i,j,k)%level = l
        !set the present cell as children of the parent
        octree%layers(l-1)%lcells(i,j,k)%children(p)%p => &
        octree%layers(l)%lcells(2*(i-1)+ic,2*(j-1)+jc,2*(k-1)+kc)

        octree%layers(l)%lcells(2*(i-1)+ic,2*(j-1)+jc,2*(k-1)+kc)%cen = &
               octree%xmin + &
               (real((/2*(i-1)+ic,2*(j-1)+jc,2*(k-1)+kc/),wp)-0.5_wp)* &
               (octree%xmax-octree%xmin)/real(nbox*2**(l-1),wp)

        !update the counters
        p = p+1
      
      enddo; enddo; enddo !children ic,jc,kc
    enddo; enddo; enddo !parents i,j,k
  enddo


  t0 = dust_time()
  ! == Set the neighbours and perform initializations
  do l = 1,nlevels

    imax=nbox(1)*2**(l-1); jmax=nbox(2)*2**(l-1); kmax=nbox(3)*2**(l-1);
    !cycle on the elements on the level
    do k = 1,kmax; do j = 1,jmax; do i = 1,imax
      !cycle on all the possible neighbours
      do kc = -1,1; do jc = -1,1; do ic = -1,1
        
        if( ((i+ic).ge.1 .and. (i+ic).le.imax) .and. &
            ((j+jc).ge.1 .and. (j+jc).le.jmax) .and. & 
            ((k+kc).ge.1 .and. (k+kc).le.kmax) .and. & 
            .not.(ic.eq.0 .and. jc.eq.0 .and. kc.eq.0) ) then
          octree%layers(l)%lcells(i,j,k)%neighbours(ic,jc,kc)%p =>  &
             octree%layers(l)%lcells(i+ic,j+jc,k+kc)
        endif

      enddo; enddo; enddo !neighbours ic,jc,kc

      !initialize few things
      allocate(octree%layers(l)%lcells(i,j,k)%cell_parts(0))
      octree%layers(l)%lcells(i,j,k)%npart = 0
      octree%layers(l)%lcells(i,j,k)%active = .true.
      octree%layers(l)%lcells(i,j,k)%leaf = .false.
      allocate(octree%layers(l)%lcells(i,j,k)%interaction_list(0) )
      call octree%layers(l)%lcells(i,j,k)%mp%init(octree%pexp)

    enddo; enddo; enddo !layer cells i,j,k
  enddo
  t1 = dust_time()
  write(msg,'(A,F9.3,A)') 'Initialized neighbours in: ' , t1 - t0,' s.'
  call printout(msg)



  ! == Set the interaction list
  !first level: interaction list is all the elements which are not neighbours
  imax=nbox(1); jmax=nbox(2); kmax=nbox(3);
  !cycle on the elements on the level
  do k = 1,kmax; do j = 1,jmax; do i = 1,imax
    !cycle again on all the cells on first level
    do kc = 1,kmax; do jc = 1,jmax; do ic = 1,imax
          indx = (/ic,jc,kc/)
          if(any(abs(indx - octree%layers(1)%lcells(i,j,k)%cart_index).gt.1)) then
            call push_ptr(octree%layers(1)%lcells(i,j,k)%interaction_list, &
            octree%layers(1)%lcells(ic,jc,kc))
          endif
    enddo; enddo; enddo !layer cells i,j,k
  enddo; enddo; enddo !layer cells i,j,k


  t0 = dust_time()
  !second level and downwards: the interaction list are all the children of
  !parent which are not a neighbour of the cell
  do l = 2,nlevels
    imax=nbox(1)*2**(l-1); jmax=nbox(2)*2**(l-1); kmax=nbox(3)*2**(l-1);
    !cycle on the elements on the level
    do k = 1,kmax; do j = 1,jmax; do i = 1,imax
      !cycle on all the neighbours of the parent
      do kc = -1,1; do jc = -1,1; do ic = -1,1
        !cycle on all the childs
        do child = 1,8
          if(associated(octree%layers(l)%lcells(i,j,k)%parent%p%neighbours(ic,jc,kc)%p)) then
            indx = octree%layers(l)%lcells(i,j,k)%parent%p%neighbours(ic,jc,kc)%p%children(child)%p%cart_index
            if(any(abs(indx - octree%layers(l)%lcells(i,j,k)%cart_index).gt.1)) then
              call push_ptr(octree%layers(l)%lcells(i,j,k)%interaction_list, &
            octree%layers(l)%lcells(i,j,k)%parent%p%neighbours(ic,jc,kc)%p%children(child)%p)
            endif
          endif

        enddo
      enddo; enddo; enddo !parents neighbours ic,jc,kc

    enddo; enddo; enddo !layer cells i,j,k
  enddo
  t1 = dust_time()
  write(msg,'(A,F9.3,A)') 'Initialized interaction lists in: ' , t1 - t0,' s.'
  call printout(msg)

  
  !pre-build all the kernel derivatives for cell-cell interactions
  do l=1,nlevels
    allocate(octree%layers(l)%ker_der(-3:3,-3:3,-3:3))
    do k=-3,3; do j=-3,3; do i=-3,3
      if(any(abs((/i,j,k/)) .ge. 2)) then
        call octree%layers(l)%ker_der(i,j,k)%&
        compute_der(octree%layers(l)%cell_size*real((/-i,-j,-k/),wp), &
        octree%delta, octree%pexp_der)
      endif
    enddo; enddo; enddo
  enddo

end subroutine initialize_octree

!----------------------------------------------------------------------

!> Destroy a wake panels type by simply passing it as intent(out)
subroutine destroy_octree(octree)
 type(t_octree), intent(out) :: octree


end subroutine

!----------------------------------------------------------------------

!> Sort particles inside the octree grid
subroutine sort_particles(part,octree)
 type(t_vortpart_p), intent(in), target :: part(:)
 type(t_octree), intent(inout), target :: octree

 integer :: ip
 integer :: idx(3)
 real(wp) :: csize
 integer :: l, i,j,k, child
 integer :: imax, jmax, kmax, ll, nl
 logical :: got_leaves
 
  ll = octree%nlevels
  csize = octree%layers(ll)%cell_size
  
  !set all the cells to empty
  do l = 1,ll
    imax=octree%nbox(1)*2**(l-1); 
    jmax=octree%nbox(2)*2**(l-1); 
    kmax=octree%nbox(3)*2**(l-1);
    !cycle on the elements on the level
    do k = 1,kmax; do j = 1,jmax; do i = 1,imax
      call reset_cell(octree%layers(l)%lcells(i,j,k))
    enddo; enddo; enddo !layer cells i,j,k
  enddo
  !deallocate(octree%leaves); allocate(octree%leaves(0))

  t0 = dust_time()
  !cycle on all the particles
  do ip=1,size(part)
    
    !check in which cell at the lowest level it is located
    idx = ceiling((part(ip)%p%cen-octree%xmin)/csize)
    !add the particle to the lowest level
    octree%layers(ll)%lcells(idx(1),idx(2),idx(3))%npart = &
                    octree%layers(ll)%lcells(idx(1),idx(2),idx(3))%npart + 1
    call push_ptr(octree%layers(ll)%lcells(idx(1),idx(2),idx(3))%cell_parts,part(ip)%p)

  enddo
  t1 = dust_time()
  write(msg,'(A,F9.3,A)') 'Sorted particles in: ' , t1 - t0,' s.'
  call printout(msg)



  nl = 0
  t0 = dust_time()
  !Bottom level: just check if are leaves
  imax=octree%nbox(1)*2**(ll-1); 
  jmax=octree%nbox(2)*2**(ll-1); 
  kmax=octree%nbox(3)*2**(ll-1);
  do k = 1,kmax; do j = 1,jmax; do i = 1,imax
    if( octree%layers(ll)%lcells(i,j,k)%npart .ge. min_part_4_cell) then
      octree%layers(ll)%lcells(i,j,k)%leaf = .true.
      !call push_ptr(octree%leaves, octree%layers(ll)%lcells(i,j,k))
      nl = nl+1
      octree%leaves(nl)%p => octree%layers(ll)%lcells(i,j,k)
    endif
  enddo; enddo; enddo !layer cells i,j,k


  !From the bottom level upwards
  do l = ll-1,1,-1
    t0 = dust_time()
    imax=octree%nbox(1)*2**(l-1); 
    jmax=octree%nbox(2)*2**(l-1); 
    kmax=octree%nbox(3)*2**(l-1);
    !cycle on the elements on the level
    do k = 1,kmax; do j = 1,jmax; do i = 1,imax

    !cycle on the children, gather the number of particles and if are 
    !leaves
    got_leaves = .false.
    do child = 1,8
      if(octree%layers(l)%lcells(i,j,k)%children(child)%p%leaf) &
                                                          got_leaves = .true.
      octree%layers(l)%lcells(i,j,k)%npart = &
            octree%layers(l)%lcells(i,j,k)%npart + &
            octree%layers(l)%lcells(i,j,k)%children(child)%p%npart
      !push all the particles pointers
      call push_ptr(octree%layers(l)%lcells(i,j,k)%cell_parts, &
      octree%layers(l)%lcells(i,j,k)%children(child)%p%cell_parts)
    enddo !child
    if (.not. octree%layers(l)%lcells(i,j,k)%branch) then
      !if it is not a branch analyse the situation of the children
      if(got_leaves) then
        !Some of the children are leaves: set all children as leaves and then
        !set the current cell and all the branch upwards as branch
        do child = 1,8
          if(.not. octree%layers(l)%lcells(i,j,k)%children(child)%p%leaf) then
            octree%layers(l)%lcells(i,j,k)%children(child)%p%leaf = .true.
            !call push_ptr(octree%leaves, octree%layers(l)%lcells(i,j,k)%children(child)%p)
            nl = nl+1
            octree%leaves(nl)%p => octree%layers(l)%lcells(i,j,k)%children(child)%p
          endif
        enddo
        call set_branch(octree%layers(l)%lcells(i,j,k))
      else
        !None of the children are leaves. Inherit the particles and check 
        !if the current cell is a leaf
        do child = 1,8
          octree%layers(l)%lcells(i,j,k)%children(child)%p%active = .false.
          !call push_ptr(octree%layers(l)%lcells(i,j,k)%cell_parts, &
          !octree%layers(l)%lcells(i,j,k)%children(child)%p%cell_parts)
        enddo
      endif
    else
      !if it is a branch we still need to check that one of the children
      !is not orphaned, i.e. not a leaf nor a branch. Since it is a sibling
      !of a branch it must become a leaf even if it has not enough particles
      do child = 1,8
        if( .not. octree%layers(l)%lcells(i,j,k)%children(child)%p%leaf &
      .and. .not. octree%layers(l)%lcells(i,j,k)%children(child)%p%branch) then
          octree%layers(l)%lcells(i,j,k)%children(child)%p%leaf = .true.
          !call push_ptr(octree%leaves, octree%layers(l)%lcells(i,j,k)%children(child)%p)
          nl = nl+1
          octree%leaves(nl)%p => octree%layers(l)%lcells(i,j,k)%children(child)%p
        endif
      enddo !child
    endif

    enddo; enddo; enddo !layer cells i,j,k
  enddo
  t1 = dust_time()
  write(msg,'(A,F9.3,A)') 'Checked leaves  in: ' , t1 - t0,' s.'
  call printout(msg)

  octree%nleaves = nl


  !!From the lowest level to the first, add all the particles
  !!TODO: this might not be 100% needed. When in the following step I will
  !! need to push the particles upward, I will need to know the particles in 
  !! the siblings of the parents. I could calculate them on the fly there, 
  !! but this way it is simpler
  !t0 = dust_time()
  !do l = ll-1,1,-1
  !  imax=octree%nbox(1)*2**(l-1); 
  !  jmax=octree%nbox(2)*2**(l-1); 
  !  kmax=octree%nbox(3)*2**(l-1);
  !  !cycle on the elements on the level
  !  do k = 1,kmax; do j = 1,jmax; do i = 1,imax
  !    !cycle on all the childs
  !    do child = 1,8
  !      !gather the particles from all the childs
  !      octree%layers(l)%lcells(i,j,k)%p%npart = &
  !                             octree%layers(l)%lcells(i,j,k)%p%npart + &
  !                octree%layers(l)%lcells(i,j,k)%p%children(child)%p%npart
  !                             
  !    enddo
  !  enddo; enddo; enddo !layer cells i,j,k
  !enddo
  !t1 = dust_time()
  !write(msg,'(A,F9.3,A)') 'Summed particles upward in: ' , t1 - t0,' s.'
  !call printout(msg)

  !!Only on the lowest level: set the leaves (potentially recursively)
  !t0 = dust_time()
  !deallocate(octree%leaves); allocate(octree%leaves(0))
  !imax=octree%nbox(1)*2**(ll-1); 
  !jmax=octree%nbox(2)*2**(ll-1); 
  !kmax=octree%nbox(3)*2**(ll-1);
  !!cycle on the elements on the level
  !do k = 1,kmax; do j = 1,jmax; do i = 1,imax
  !  call check_cell_content(octree%layers(ll)%lcells(i,j,k)%p, octree%leaves)
  !enddo; enddo; enddo !layer cells i,j,k
  !t1 = dust_time()
  !write(msg,'(A,F9.3,A)') 'Leaves set in: ' , t1 - t0,' s.'
  !call printout(msg)

end subroutine sort_particles


!----------------------------------------------------------------------

subroutine perform_multipole(part,octree)
 type(t_vortpart_p), intent(in), target :: part(:)
 type(t_octree), intent(inout) :: octree

 integer :: i, j, k, lv, l, child, il, ip, ine, ipp, m
 integer :: imax, jmax, kmax, ll
 integer :: idx_diff(3)
 real(wp) :: Rnorm2

  ll = octree%nlevels
  !reset everything (might be done elsewhere)
  do l = ll,1,-1
    imax=octree%nbox(1)*2**(l-1); 
    jmax=octree%nbox(2)*2**(l-1); 
    kmax=octree%nbox(3)*2**(l-1);
    !cycle on the elements on the level
    do k = 1,kmax; do j = 1,jmax; do i = 1,imax
         octree%layers(l)%lcells(i,j,k)%mp%a = 0.0_wp
         octree%layers(l)%lcells(i,j,k)%mp%b = 0.0_wp
    enddo; enddo; enddo !layer cells i,j,k

  enddo


  !For all the leaves, calculate the multipole
  t0 = dust_time()
  do lv = 1, octree%nleaves
    call octree%leaves(lv)%p%mp%leaf_M(&
               octree%leaves(lv)%p%cen, &
               octree%leaves(lv)%p%cell_parts, &
               octree%pexp  )
  enddo
  t1 = dust_time()
  write(msg,'(A,F9.3,A)') 'Calculated multipoles in leaves in: ' , t1 - t0,' s.'
  call printout(msg)


  !layer by layer (except the last, in which are only leaves or inactive)
  !perform multipole to multipole
  t0 = dust_time()
  do l = ll-1,1,-1
    imax=octree%nbox(1)*2**(l-1); 
    jmax=octree%nbox(2)*2**(l-1); 
    kmax=octree%nbox(3)*2**(l-1);

    !cycle on the elements on the level
    do k = 1,kmax; do j = 1,jmax; do i = 1,imax
      if(octree%layers(l)%lcells(i,j,k)%active .and. &
         octree%layers(l)%lcells(i,j,k)%branch) then
         !if it is active and is a branch perform the M2M from all the
         !children
        do child = 1,8
          call octree%layers(l)%lcells(i,j,k)%mp%M2M(&
               octree%layers(l)%lcells(i,j,k)%cen, &
               octree%layers(l)%lcells(i,j,k)%children(child)%p%mp, &
               octree%layers(l)%lcells(i,j,k)%children(child)%p%cen, &
               octree%pexp )
        enddo 
      endif
    enddo; enddo; enddo !layer cells i,j,k
  enddo
  t1 = dust_time()
  write(msg,'(A,F9.3,A)') 'Calculated M2M in: ' , t1 - t0,' s.'
  call printout(msg)

  
  !from the top, interact in the interaction list, then pass to the children
  !(if not already at the lowest level)
  do l = 1,ll
    imax=octree%nbox(1)*2**(l-1); 
    jmax=octree%nbox(2)*2**(l-1); 
    kmax=octree%nbox(3)*2**(l-1);

  t0 = dust_time()
    !cycle on the elements on the level
    do k = 1,kmax; do j = 1,jmax; do i = 1,imax
      if(octree%layers(l)%lcells(i,j,k)%active) then
         !if it is active perform M2L with all the interaction list
        do il = 1,size(octree%layers(l)%lcells(i,j,k)%interaction_list)
          idx_diff = octree%layers(l)%lcells(i,j,k)%interaction_list(il)%p%cart_index -&
                  (/i,j,k/)
          call octree%layers(l)%lcells(i,j,k)%mp%M2L(&
               octree%layers(l)%ker_der(idx_diff(1),idx_diff(2),idx_diff(3)), &
               octree%pexp, &
               octree%pexp_der, &
               octree%layers(l)%lcells(i,j,k)%interaction_list(il)%p%mp)
        enddo 
      endif
    enddo; enddo; enddo !layer cells i,j,k

  t1 = dust_time()
  write(msg,'(A,I0,A,F9.3,A)') 'Calculated M2L layer ',l,' in: ' , t1 - t0,' s.'
  call printout(msg)

  t0 = dust_time()
    do k = 1,kmax; do j = 1,jmax; do i = 1,imax
      if(octree%layers(l)%lcells(i,j,k)%active .and. l.lt.ll) then
        !if active and not in the bottom layer, bring down the local expansion
        ! to the children
        do child = 1,8
          call octree%layers(l)%lcells(i,j,k)%children(child)%p%mp%L2L(&
               octree%layers(l)%lcells(i,j,k)%children(child)%p%cen, &
               octree%layers(l)%lcells(i,j,k)%mp, &
               octree%layers(l)%lcells(i,j,k)%cen, &
               octree%pexp )
        enddo 
      endif
    enddo; enddo; enddo !layer cells i,j,k
  t1 = dust_time()
  write(msg,'(A,I0,A,F9.3,A)') 'Calculated L2L layer ',l,' in: ' , t1 - t0,' s.'
  call printout(msg)
  enddo

  !for all the leaves apply the local expansion and then local interactions 
  !WARNING: this is all temporary due to the fact that we are just
  !calculating the potential, will be changed when calculating the velocity,
  !maybe moving it to the multipole module
  t0 = dust_time()
  do lv = 1, octree%nleaves
    !I am on a leaf, cycle on all the particles inside the leaf
    do ip = 1,octree%leaves(lv)%p%npart

      !Reset the interaction LAZY: the potential is stored in the direction..
      octree%leaves(lv)%p%cell_parts(ip)%p%dir = 0.0_wp

      !first apply the local multipole expansion
      do m = 1,size(octree%leaves(lv)%p%mp%b,2)
        octree%leaves(lv)%p%cell_parts(ip)%p%dir(1) = &
          octree%leaves(lv)%p%cell_parts(ip)%p%dir(1) + &
          octree%leaves(lv)%p%mp%b(1,m)* &
          product((octree%leaves(lv)%p%cell_parts(ip)%p%cen- &
          octree%leaves(lv)%p%cen)**octree%pexp%pwr(:,m))
      enddo

      !then interact with all the neighbouring cell particles
      do k=-1,1; do j=-1,1; do i=-1,1
        !if(all((/i,j,k/).ne.0)) then
        if(associated(octree%leaves(lv)%p%neighbours(i,j,k)%p)) then
          do ipp = 1,octree%leaves(lv)%p%neighbours(i,j,k)%p%npart

            Rnorm2 = sum(( &
            octree%leaves(lv)%p%cell_parts(ip)%p%cen - &
            octree%leaves(lv)%p%neighbours(i,j,k)%p%cell_parts(ipp)%p%cen &
            )**2) + octree%delta**2 

            octree%leaves(lv)%p%cell_parts(ip)%p%dir(1) = &
              octree%leaves(lv)%p%cell_parts(ip)%p%dir(1) + &
              octree%leaves(lv)%p%neighbours(i,j,k)%p%cell_parts(ipp)%p%mag/(4.0_wp*pi*sqrt(Rnorm2))
          enddo
        endif
      enddo; enddo; enddo

      !finally interact with the particles inside the cell 
          do ipp = 1,octree%leaves(lv)%p%npart
            if (ipp .ne. ip) then
              Rnorm2 = sum(( &
              octree%leaves(lv)%p%cell_parts(ip)%p%cen - &
              octree%leaves(lv)%p%cell_parts(ipp)%p%cen &
              )**2) + octree%delta**2 

              octree%leaves(lv)%p%cell_parts(ip)%p%dir(1) = &
                octree%leaves(lv)%p%cell_parts(ip)%p%dir(1) + &
                octree%leaves(lv)%p%cell_parts(ipp)%p%mag/(4.0_wp*pi*sqrt(Rnorm2))
            endif
          enddo

    enddo

  enddo
  t1 = dust_time()
  write(msg,'(A,F9.3,A)') 'Calculated leaves interactions in: ' , t1 - t0,' s.'
  call printout(msg)

end subroutine perform_multipole


!----------------------------------------------------------------------

!> Check the quantity of particles in a cell
!!
!! If there are too few particles, these are bound together and applied to
!! the parent
recursive subroutine check_cell_content(cell, leaves)
 type(t_cell), intent(inout) :: cell
 type(t_cell_p), allocatable, intent(inout) :: leaves(:)
 
 integer :: sib
 logical :: enough
  

  if(cell%npart .ge. min_part_4_cell) then
    ! Enough particles, since we started from the bottom, it is a leaf, set
    ! it as leaf and cut out the descendant
    call set_leaf(cell,.true.)  
    call push_ptr(leaves, cell)

  else
    !Not enough particles in this cell
    enough = .false.
    !We have to check that all the siblings have not enough particles
    !oviously only if the cell has siblings, otherwise is top level and
    !it will be a leaf anyway
    if(associated(cell%parent%p)) then
      do sib = 1,8
        if(cell%parent%p%children(sib)%p%npart .ge. min_part_4_cell) then
          enough = .true. 
        endif
      enddo
    else
      !it does not have a parent, it is top level, force to be a leaf
      enough = .true.
    endif

    if (.not. enough) then
      !all the siblings have not enough particles,  all are going to be 
      ! set inactive, push the pointers to particles to the parent
      do sib = 1,8
        !TODO: for some reasons the generic does not work here
        call push_ptr(cell%parent%p%cell_parts, cell%parent%p%children(sib)%p%cell_parts)
      enddo
      !escalate this call upwards
      call check_cell_content(cell%parent%p, leaves)

    else
      !the present cell has not enough particles, however the siblings
      !do, so the present cell becomes a leaf anyway
      call set_leaf(cell,.true.)  
      call push_ptr(leaves, cell)

    endif


  endif

end subroutine check_cell_content

!----------------------------------------------------------------------

!> Set a whole branch inactive under a leaf. 
!!
!! Used to set a whole branch below a leaf inactive. When called with 
!! leaf=.true. the present cell is set as leaf, and all the branch 
!! underneath is set to inactive (children and all descendants)
recursive subroutine set_leaf(cell,leaf)
 type(t_cell), intent(inout) :: cell
 logical :: leaf

 integer :: i 
  
  if (.not.leaf) cell%active = .false.
 
  do i=1,8
    if(associated(cell%children(i)%p)) then
      !call set_inactive_tree(cell%children(i)%p, leaf=.false.)
      call set_leaf(cell%children(i)%p, .false.)
    endif
  enddo

end subroutine set_leaf

!----------------------------------------------------------------------

recursive subroutine set_branch(cell)
 type(t_cell), intent(inout) :: cell

 cell%branch = .true.

 if(associated(cell%parent%p)) call set_branch(cell%parent%p)

end subroutine set_branch

!----------------------------------------------------------------------

!> Add a certain number of particles upward in the tree
recursive subroutine add_part_upwards(cell, npart)
 type(t_cell), intent(inout) :: cell
 integer, intent(in) :: npart

  cell%npart = cell%npart + npart
  if(associated(cell%parent%p)) then
    call add_part_upwards(cell%parent%p, npart)
  endif

end subroutine add_part_upwards

!----------------------------------------------------------------------

!> Push another pointer to a list of element pointers
subroutine push_cell_ptr_scalar(pointers, element)
 type(t_cell_p), allocatable, target, intent(inout) :: pointers(:)
 type(t_cell), intent(in), target :: element

 type(t_cell_p) :: tmp_ptr(size(pointers)+1)
 integer :: i, l

  l = size(pointers)
  do i=1,l
    tmp_ptr(i)%p => pointers(i)%p
  enddo
  deallocate(pointers); allocate(pointers(l+1))
  do i=1,l
    pointers(i)%p => tmp_ptr(i)%p
  enddo
  pointers(l+1)%p => element

end subroutine push_cell_ptr_scalar

!----------------------------------------------------------------------

!> Push another pointer list to a list of element pointers
subroutine push_cell_ptr_vec(pointers, elements)
 type(t_cell_p), allocatable, target, intent(inout) :: pointers(:)
 type(t_cell_p), target, intent(in) :: elements(:)

 type(t_cell_p) :: tmp_ptr(size(pointers)+1)
 integer :: i, l, m

  l = size(pointers)
  m = size(elements)
  do i=1,l
    tmp_ptr(i)%p => pointers(i)%p
  enddo
  deallocate(pointers); allocate(pointers(l+m))
  do i=1,l
    pointers(i)%p => tmp_ptr(i)%p
  enddo
  do i = 1,m
    pointers(l+i)%p => elements(i)%p
  enddo

end subroutine push_cell_ptr_vec

!----------------------------------------------------------------------

!> Push another pointer to a list of particles pointers
subroutine push_part_ptr_scalar(pointers, particle)
 type(t_vortpart_p), allocatable, target, intent(inout) :: pointers(:)
 type(t_vortpart), intent(in), target :: particle

 type(t_vortpart_p) :: tmp_ptr(size(pointers)+1)
 integer :: i, l

  l = size(pointers)
  do i=1,l
    tmp_ptr(i)%p => pointers(i)%p
  enddo
  deallocate(pointers); allocate(pointers(l+1))
  do i=1,l
    pointers(i)%p => tmp_ptr(i)%p
  enddo
  pointers(l+1)%p => particle

end subroutine push_part_ptr_scalar

!----------------------------------------------------------------------

!> Push another pointer list to a list of particle pointers
subroutine push_part_ptr_vec(pointers, elements)
 type(t_vortpart_p), allocatable, target, intent(inout) :: pointers(:)
 type(t_vortpart_p), target, intent(in) :: elements(:)

 type(t_vortpart_p) :: tmp_ptr(size(pointers)+1)
 integer :: i, l, m

  l = size(pointers)
  m = size(elements)
  do i=1,l
    tmp_ptr(i)%p => pointers(i)%p
  enddo
  deallocate(pointers); allocate(pointers(l+m))
  do i=1,l
    pointers(i)%p => tmp_ptr(i)%p
  enddo
  do i = 1,m
    pointers(l+i)%p => elements(i)%p
  enddo

end subroutine push_part_ptr_vec

!----------------------------------------------------------------------

!> Reset the particles stored inside the cell
subroutine reset_cell(cell)
 type(t_cell), intent(inout) :: cell

 deallocate(cell%cell_parts)
 allocate(cell%cell_parts(0))
 cell%npart = 0
 cell%active = .true.
 cell%branch = .false.
 cell%leaf = .false.

end subroutine reset_cell

!----------------------------------------------------------------------

end module mod_octree
