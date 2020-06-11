module mod_precice_rbf

use mod_param, only: &
    wp

implicit none

private

public :: t_precice_rbf


!> RBF coupling structures
type :: t_precice_rbf
  !> Coupling nodes
  real(wp), allocatable :: nodes(:,:)
  !> Indices (surface to structure)
  integer,  allocatable :: ind(:,:) 
  !> Weights (surface to structure)
  real(wp), allocatable :: wei(:,:)

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
subroutine build_connectivity(this, rr)
  class(t_precice_rbf), intent(inout) :: this
  real(wp),             intent(in)    :: rr(:,:)
  
  real(wp), allocatable :: dist_all(:), wei_v(:)
  integer , allocatable ::              ind_v(:)

  integer :: np, ns
  integer :: ip, is

  !> Number of surface points
  np = size(rr,2)
  !> Number of coupling nodes of the structure
  ns = size(this%nodes,2)

  allocate(this%ind(this%n_wei, np)); this%ind = 0
  allocate(this%wei(this%n_wei, np)); this%wei = 0.0_wp

  allocate(dist_all(np)); dist_all = 0.0_wp

  do ip = 1, np

    !> Distance of the surface nodes from the structural nodes
    do is = 1, ns
      dist_all(is) = norm2( rr(:,ip) - this%nodes(:,is) )
    end do

    call sort_vector_real( dist_all, this%n_wei, wei_v, ind_v )

    wei_v = 1.0_wp / wei_v**this%w_order
    wei_v = wei_v / sum(wei_v)

    this%wei(:,ip) = wei_v
    this%ind(:,ip) = ind_v

  enddo

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
  integer , intent(in) :: nel
  real(wp), allocatable, intent(out):: sor(:)
  integer , allocatable, intent(out):: ind(:)

  real(wp):: minv
  integer :: i

  allocate(sor(nel)); sor = 0.0_wp
  allocate(ind(nel)); ind = 0

  minv = minval( vec )
  do i = 1, nel
    sor(i) = maxval( vec, 1 )
    ind(i) = maxloc( vec, 1 )
    vec( ind(i) ) = minv - 0.1_wp ! naif
  end do


end subroutine sort_vector_real


end module mod_precice_rbf
