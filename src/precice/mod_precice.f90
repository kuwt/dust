module mod_precice

use mod_param, only: &
    wp

use mod_geometry, only: &
    t_geo

implicit none

private

public :: t_precice

!> Parameters
integer, parameter :: precice_mcl = 50 ! precice_max_char_len

!> PreCICE mesh -------------------------------------------------
type t_precice_mesh
  character(len=precice_mcl) :: mesh_name
  integer                    :: mesh_id
  integer , allocatable      :: node_ids(:)
  real(wp), allocatable      :: nodes(:,:)
  integer                    :: nnodes
  integer                    :: ndim
end type t_precice_mesh

!> PreCICE field ------------------------------------------------
type :: t_precice_field
  integer :: fid
  character(len=precice_mcl) :: fname
  character(len=precice_mcl) :: fio    ! 'read'/'write'
  character(len=precice_mcl) :: ftype  ! 'scalar'/'vector'
  real(wp), allocatable :: fdata(:,:)  ! (1,nnodes)/(nd,nnodes)
end type t_precice_field

!> PreCICE type -------------------------------------------------
type :: t_precice

  !> PreCICE configuration
  character(len=precice_mcl) :: config_file_name
  character(len=precice_mcl) :: solver_name
  character(len=precice_mcl) :: mesh_name     ! *** to do *** just one mesh?
  integer :: comm_rank, comm_size

  !> PreCICE mesh
  type(t_precice_mesh) :: mesh

  !> PreCICE fields
  integer :: n_fields
  type(t_precice_field), allocatable :: fields(:)

  !> PreCICE variables
  real(wp) :: dt_precice
  integer  :: is_ongoing
  character(len=precice_mcl) :: write_initial_data, &
                                read_it_checkp, write_it_checkp

  contains

  procedure, pass(this) :: initialize
  procedure, pass(this) :: initialize_mesh
  procedure, pass(this) :: initialize_fields

end type t_precice
!> --------------------------------------------------------------

contains

!----------------------------------------------------------------
!>
subroutine initialize(this)
  class(t_precice), intent(inout) :: this

  write(*,*) ' Using PreCICE '

  !> Default input for dust in a preCICE coupled simulation
  this % config_file_name = './../precice-config.xml'
  this % solver_name = 'dust'
  this %   mesh_name = 'dust_mesh'
  this % comm_rank = 0
  this % comm_size = 1

  !> Initialize PreCICE participant and mesh
  call precicef_create( this % solver_name, & 
                        this % config_file_name, &
                        this % comm_rank, &
                        this % comm_size )

  this % mesh % mesh_name = this % mesh_name
  call precicef_get_dims( this % mesh % ndim )
  call precicef_get_mesh_id( this % mesh % mesh_name, this % mesh % mesh_id )

  !> Initialize some preCICE variables
  this % write_initial_data(1:precice_mcl)='                                                  '
  this % read_it_checkp(    1:precice_mcl)='                                                  '
  this % write_it_checkp(   1:precice_mcl)='                                                  '

  call precicef_action_write_initial_data(this % write_initial_data)
  call precicef_action_read_iter_checkp(  this % read_it_checkp)
  call precicef_action_write_iter_checkp( this %write_it_checkp)

end subroutine initialize

!----------------------------------------------------------------
!> 
subroutine initialize_mesh( this, geo )
  class(t_precice), intent(inout) :: this
  type(t_geo), intent(in) :: geo

  integer :: i, i_comp, n_comp
  integer :: dnnodes, nnodes

  !> TODO: add component field, that describe if the component
  ! participates to the coupling. So far, all the components 
  ! participate.
  logical :: comp_coupling = .true.

  n_comp = size(geo%components)

  !> Count n. of nodes ==========================================
  nnodes = 0
  do i_comp = 1, n_comp
    if ( comp_coupling ) then
      select case( trim(geo%components(i_comp)%comp_el_type) )
        case('l')
          !> Set PreCICE nodes with LE only! ***to do***
          dnnodes = size(geo%components(i_comp)%i_points) / 2
          !> Increment number of nodes
          nnodes = nnodes + dnnodes
      ! case('v')
      ! case('p')
        case default
      end select
    end if
  end do

  !> Allocate participant%mesh fields ===========================
  allocate(this%mesh%node_ids(nnodes)); this%mesh%node_ids = 0
  allocate(this%mesh%nodes( this%mesh%ndim, nnodes ))
  nnodes = 0
  do i_comp = 1, n_comp
    if ( comp_coupling ) then
      select case( trim(geo%components(i_comp)%comp_el_type) )
        case('l')
          
          dnnodes = size(geo%components(i_comp)%i_points) / 2
          !> Here in the local reference frame ! ***to do***
          this%mesh%nodes(:,nnodes+1:nnodes+dnnodes) = &
            geo%components(i_comp)%loc_points(:,1:2*dnnodes:2) + 0.001_wp
          nnodes = nnodes + dnnodes

      ! case('v')
      ! case('p')
        case default
      end select
    end if
  end do

  !> check ---
  write(*,*) ' PreCICE, mesh % nodes: '
  do i = 1 , nnodes
    write(*,*) this%mesh%nodes(:,i)
  end do

  !> Initialize mesh ============================================
  !> Nodes
  call precicef_set_vertices( this%mesh%mesh_id, nnodes,  &
                              this%mesh%nodes, this%mesh%node_ids )
  this%mesh%nnodes = nnodes

  !> Connectivity
  !> Edges
  ! ...
  !> Elements
  ! ...


end subroutine initialize_mesh

!----------------------------------------------------------------
!> 
subroutine initialize_fields( this )
  class(t_precice), intent(inout) :: this

  !> Expected fields
  integer, parameter :: n_max_fields = 6
  character(len=precice_mcl) :: field_list(n_max_fields)
  character(len=precice_mcl) :: type_list(n_max_fields)
  character(len=precice_mcl) :: io_list(n_max_fields)

  integer :: i, nf, hasdata

  !> Fields that can be exchanged through PreCICE. Add a field if needed
  field_list(1) = 'Position'        ; type_list(1) = 'vector' ; io_list(1) = 'read' 
  field_list(2) = 'Velocity'        ; type_list(2) = 'vector' ; io_list(2) = 'read' 
  field_list(3) = 'Rotation'        ; type_list(3) = 'vector' ; io_list(3) = 'read' 
  field_list(4) = 'AngularVelocity' ; type_list(4) = 'vector' ; io_list(4) = 'read' 
  field_list(5) = 'Force'           ; type_list(5) = 'vector' ; io_list(5) = 'write'
  field_list(6) = 'Moment'          ; type_list(6) = 'vector' ; io_list(6) = 'write'

  !> Count n. of exchanged fields
  nf = 0
  do i = 1, n_max_fields
    call precicef_has_data( trim(field_list(i)), this%mesh%mesh_id, hasdata )
    if ( hasdata .eq. 1 ) nf = nf + 1
  end do

  allocate( this%fields( nf ) )
  nf = 0
  do i = 1, n_max_fields
    call precicef_has_data( trim(field_list(i)), this%mesh%mesh_id, hasdata )
    if ( hasdata .eq. 1 ) then
      nf = nf + 1
      this%fields(nf)%fname = trim(field_list(i))
      this%fields(nf)%fio   = trim(   io_list(i))
      this%fields(nf)%ftype = trim( type_list(i))
      call precicef_get_data_id( trim(this%fields(nf)%fname), &
                                 this%mesh%mesh_id, this%fields(nf)%fid )
      !> Allocate and initialize fields to zero *** to do *** where to write initial data?
      if ( trim(this%fields(nf)%ftype) .eq. 'scalar' ) then
        allocate( this%fields(nf)%fdata(1, this%mesh%nnodes) )
        this%fields(nf)%fdata = 0.0_wp
      else if ( trim(this%fields(nf)%ftype) .eq. 'vector' ) then
        allocate( this%fields(nf)%fdata( this%mesh%ndim, &
                                         this%mesh%nnodes) )
        this%fields(nf)%fdata = 0.0_wp
      end if
    end if
  end do


end subroutine initialize_fields

!----------------------------------------------------------------

end module mod_precice
