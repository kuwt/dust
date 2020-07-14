module mod_precice_rbf

use mod_param, only: &
    wp

implicit none

private

public :: t_precice_rbf

!> Connectivity, indices and weights of the structural nodes driving
! the surface points
type :: t_rbf_conn

  !> Indices (surface to structure)
  integer,  allocatable :: ind(:,:) 
  !> Weights (surface to structure)
  real(wp), allocatable :: wei(:,:)

end type t_rbf_conn

!> RBF coupling structures
type :: t_precice_rbf

  !> Coupling nodes
  real(wp), allocatable :: nodes(:,:)

  !> Grid nodes connectivity
  type(t_rbf_conn) :: nod

  !> Elem centers connectivity
  type(t_rbf_conn) :: cen

  !> --- Parameters of rbf interpolation ---
  !> Number of points for transferring the motion from the structure to the
  ! surface, through a weighted average
  ! *** to do *** hardcoded, so far. Read as an input, with a default value,
  ! equal to 2 (or 1?)
  integer :: n_wei = 2

  !> Order of the norm used for computing distance-based weights
  ! *** to do *** hardcoded, so far. Read as an input, with a default value,
  ! equal to 1 (or 2?)
  integer :: w_order = 1

  contains

  procedure, pass(this) :: build_connectivity

end type t_precice_rbf

!> --------------------------------------------------------------
contains
!----------------------------------------------------------------

!----------------------------------------------------------------
!> Read local coordinates of surface nodes, rr, and build data for
! structure to surface interpolation of the motion
subroutine build_connectivity(this, rr, ee)
  class(t_precice_rbf), intent(inout) :: this
  real(wp),             intent(in)    :: rr(:,:)
  integer ,             intent(in)    :: ee(:,:)
  
  real(wp), allocatable :: dist_all(:), wei_v(:)
  integer , allocatable ::              ind_v(:)
  real(wp) :: cen(3)

  integer :: np, ns, ne, n
  integer :: ip, is, ie

  !> Number of surface points
  np = size(rr,2)
  !> Number of surface elems
  ne = size(ee,2)
  !> Number of coupling nodes of the structure
  ns = size(this%nodes,2)

  ! ! debug ---
  ! write(*,*) ' shape(rr): ', shape(rr)
  ! write(*,*) ' shape(ee): ', shape(ee)

  ! write(*,*); write(*,*) ' rr: '
  ! do ip = 1, np
  !   write(*,*) rr(:,ip)
  ! end do
  ! write(*,*); write(*,*) ' this%nodes: '
  ! do is = 1, ns
  !   write(*,*) this%nodes(:,is)
  ! end do
  ! ! debug ---

  !> === Surface nodes ===
  allocate(this%nod%ind(this%n_wei, np)); this%nod%ind = 0
  allocate(this%nod%wei(this%n_wei, np)); this%nod%wei = 0.0_wp

  allocate(dist_all(ns)); dist_all = 0.0_wp

  do ip = 1, np

    !> Distance of the surface nodes from the structural nodes
    do is = 1, ns
      dist_all(is) = norm2( rr(:,ip) - this%nodes(:,is) )
      ! write(*,*) this%nodes(:,is), '         ', dist_all(is)
    end do

    call sort_vector_real( dist_all, this%n_wei, wei_v, ind_v )

    wei_v = 1.0_wp / max( wei_v, 1e-9_wp )**this%w_order
    wei_v = wei_v / sum(wei_v)

    ! ! debug ---
    ! write(*,*) ' ip, rr(:,ip)', ip, rr(:,ip)
    ! write(*,*) ' wei_v: ' , wei_v
    ! write(*,*) ' ind_v: ' , ind_v
    ! ! debug ---

    this%nod%wei(:,ip) = wei_v
    this%nod%ind(:,ip) = ind_v

  enddo

  ! ! debug ---
  ! write(*,*) ' stop in mod_precice_rbf/build_connectivity().'; stop
  ! ! debug ---

 
  !> === Surface centers ===
  allocate(this%cen%ind(this%n_wei, np)); this%cen%ind = 0
  allocate(this%cen%wei(this%n_wei, np)); this%cen%wei = 0.0_wp

  deallocate(dist_all); allocate(dist_all(ns)); dist_all = 0.0_wp

  do ie = 1, ne

    !> Compute element center
    ! *** to do *** elem center for 'll'
    cen = 0.0_wp; n = 0
    do ip = 1, 4
      if ( ee(ip,ie) .ne. 0 ) then
        n = n + 1
        cen = cen + rr(:,ee(ip,ie))
      end if
    end do
    cen = cen / dble(n)

    !> Distance of the surface nodes from the structural nodes
    do is = 1, ns
      dist_all(is) = norm2( cen - this%nodes(:,is) )
    end do

    call sort_vector_real( dist_all, this%n_wei, wei_v, ind_v )

    !> Weight, inverse of the norm, avoid singularities
    wei_v = 1.0_wp / max( wei_v, 1e-9_wp ) **this%w_order
    wei_v = wei_v / sum(wei_v)

    this%cen%wei(:,ie) = wei_v
    this%cen%ind(:,ie) = ind_v

  enddo

  ! stop

! ! check ---
! write(*,*) 
! write(*,*) ' Check in t_precice_rbf % build_connectivity, %nod '
! do ip = 1, np
!   write(*,*) this%nod%ind(:,ip), this%nod%wei(:,ip)
! end do
! write(*,*) 
! write(*,*) ' Check in t_precice_rbf % build_connectivity, %cen '
! do ie = 1, ne
!   write(*,*) this%cen%ind(:,ie), this%cen%wei(:,ie)
! end do
! write(*,*) 
! write(*,*) ' Stop.'
! write(*,*) 
! stop
! ! check ---

  !> Deallocate and cleaning
  if ( allocated(dist_all) )  deallocate(dist_all)
  if ( allocated(wei_v   ) )  deallocate(wei_v   )
  if ( allocated(ind_v   ) )  deallocate(ind_v   )


end subroutine build_connectivity

! ---------------------------------------------------------------
!> Naif sort, copied from mod_hinges
! *** to do *** clean the implementation, moving sort_ routines
! into math module
subroutine sort_vector_real( vec, nel, sor, ind )
  real(wp), intent(inout) :: vec(:)
  integer , intent(in)    :: nel
  real(wp), allocatable, intent(out):: sor(:)
  integer , allocatable, intent(out):: ind(:)

  real(wp):: maxv
  integer :: i

  allocate(sor(nel)); sor = 0.0_wp
  allocate(ind(nel)); ind = 0

  maxv = maxval( vec )
  do i = 1, nel
    sor(i) = minval( vec, 1 )
    ind(i) = minloc( vec, 1 )
    vec( ind(i) ) = maxv + 0.1_wp ! naif
  end do


end subroutine sort_vector_real


end module mod_precice_rbf
