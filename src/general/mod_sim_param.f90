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


!> Module to define the structrues containing the simulation parameters
!!
module mod_sim_param

use mod_param, only: &
  wp, max_char_len

use mod_hdf5_io, only: &
   h5loc, &
   write_hdf5_attr

implicit none

public :: t_sim_param

private

type t_sim_param
  !Time:
  !> Start time
  real(wp) :: t0
  !> Time step
  real(wp) :: dt
  !> Final time
  real(wp) :: tend
  !> Number of timesteps
  integer  :: n_timesteps
  !> Vector of time instants
  real(wp) , allocatable :: time_vec(:)

  !Physical parameters:
  !> Free stream pressure
  real(wp) :: P_inf
  !> Free stream density
  real(wp) :: rho_inf
  !> Free stream velocity
  real(wp) :: u_inf(3)
  !> Reference velocity (magnitude of u_inf unless specified)
  real(wp) :: u_ref
  !> Free stream speed of sound
  real(wp) :: a_inf
  !> Free stream dynamic viscosity
  real(wp) :: mu_inf
  !> Free stream kinematic viscosity
  real(wp) :: nu_inf

  !Wake
  !> Scaling of the first implicit panel
  real(wp) :: first_panel_scaling
  !> Minimum velocity at the trailing edge
  real(wp) :: min_vel_at_te
  !> Is the wake rigid?
  logical :: rigid_wake
    !> Velocity of the rigid wake
    real(wp)  :: rigid_wake_vel(3)
  !> Number of wake panels
  integer :: n_wake_panels
  !> Number of wake particles
  integer :: n_wake_particles
  !> Minimum and maximum of the particles box
  real(wp) :: particles_box_min(3)
  real(wp) :: particles_box_max(3)

  !Method parameters
  !> Multiplier for far field threshold computation on doublet
  real(wp) :: FarFieldRatioDoublet
  !> Multiplier for far field threshold computation on sources
  real(wp) :: FarFieldRatioSource
  !> Thresold for considering the point in plane in doublets
  real(wp) :: DoubletThreshold
  !> Rankine Radius for vortices
  real(wp) :: RankineRad
  !> Vortex Radius for vortex particles
  real(wp) :: VortexRad
  !> Complete cutoff radius
  real(wp) :: CutoffRad
  !> use the vortex stretching or not
  logical :: use_vs
  !> use the vorticity diffusion or not
  logical :: use_vd
  !> use the penetration avoidance
  logical :: use_pa
  !> simulate viscosity effects or not
  logical :: use_ve

  !FMM parameters
  !> Employing the FMM method
  logical :: use_fmm
    !> Size of the Octree box
    real(wp) :: BoxLength
    !> Number of boxes in each direction
    integer :: NBox(3)
    !> Origin of the Octree system of boxes
    real(wp) :: OctreeOrigin(3)
    !> Number of Octree levels
    integer :: NOctreeLevels
    !> Minimum number of particles for each box
    integer :: MinOctreePart
    !> Multipole expansion degree
    integer :: MultipoleDegree
    !> Use dynamic levels
    logical :: use_dyn_layers
      !> Maximum number of octree levels
      integer :: NMaxOctreeLevels
      !> Time ratio that triggers the increase of levels
      real(wp) :: LeavesTimeRatio
    !> use particles redistribution
    logical :: use_pr


  !Handling parameters:
  !> Debug level
  integer :: debug_level
  !> Output interval 
  real(wp) :: dt_out
  !> Basename
  character(len=max_char_len) :: basename
  !> Geometry file
  character(len=max_char_len) :: GeometryFile
  !> References file
  character(len=max_char_len) :: ReferenceFile
  !> Restart from file
  logical :: restart_from_file
    !> Restart file
    character(len=max_char_len) :: restart_file
    !> Reset the time after restart
    logical :: reset_time

contains

  procedure, pass(this) :: save_param => save_sim_param

end type t_sim_param

!----------------------------------------------------------------------
contains
!----------------------------------------------------------------------

subroutine save_sim_param(this, loc)
 class(t_sim_param) :: this
 integer(h5loc), intent(in) :: loc

  call write_hdf5_attr(this%t0, 't0', loc)
  call write_hdf5_attr(this%dt, 'dt', loc)
  call write_hdf5_attr(this%tend, 'tend', loc)
  call write_hdf5_attr(this%P_inf, 'P_inf', loc)
  call write_hdf5_attr(this%rho_inf, 'rho_inf', loc)
  call write_hdf5_attr(this%u_inf, 'u_inf', loc)
  call write_hdf5_attr(this%u_ref, 'u_ref', loc)
  call write_hdf5_attr(this%a_inf, 'a_inf', loc)
  call write_hdf5_attr(this%mu_inf, 'mu_inf', loc)
  call write_hdf5_attr(this%first_panel_scaling, 'first_panel_scaling', loc)
  call write_hdf5_attr(this%min_vel_at_te, 'min_vel_at_te', loc)
  call write_hdf5_attr(this%rigid_wake, 'rigid_wake', loc)
  if(this%rigid_wake) &
    call write_hdf5_attr(this%rigid_wake_vel, 'rigid_wake_vel', loc)
  call write_hdf5_attr(this%n_wake_panels, 'n_wake_panels', loc)
  call write_hdf5_attr(this%n_wake_particles, 'n_wake_particles', loc)
  call write_hdf5_attr(this%particles_box_min, 'particles_box_min', loc)
  call write_hdf5_attr(this%particles_box_max, 'particles_box_max', loc)
  call write_hdf5_attr(this%FarFieldRatioDoublet, 'FarFieldRatioDoublet', loc)
  call write_hdf5_attr(this%FarFieldRatioSource, 'FarFieldRatioSource', loc)
  call write_hdf5_attr(this%DoubletThreshold, 'DoubletThreshold', loc)
  call write_hdf5_attr(this%RankineRad, 'RankineRad', loc)
  call write_hdf5_attr(this%VortexRad, 'VortexRad', loc)
  call write_hdf5_attr(this%CutoffRad, 'CutoffRad', loc)
  call write_hdf5_attr(this%use_vs, 'vortstretch', loc)
  call write_hdf5_attr(this%use_vd, 'vortdiff', loc)
  call write_hdf5_attr(this%use_pa, 'PenetrationAvoidance', loc)
  call write_hdf5_attr(this%use_ve, 'ViscosityEffects', loc)
  call write_hdf5_attr(this%use_fmm, 'use_fmm', loc)
  if(this%use_fmm) then
    call write_hdf5_attr(this%BoxLength, 'BoxLength', loc)
    call write_hdf5_attr(this%Nbox, 'Nbox', loc)
    call write_hdf5_attr(this%OctreeOrigin, 'OctreeOrigin', loc)
    call write_hdf5_attr(this%NOctreeLevels, 'NOctreeLevels', loc)
    call write_hdf5_attr(this%MinOctreePart, 'MinOctreePart', loc)
    call write_hdf5_attr(this%MultipoleDegree, 'MultipoleDegree', loc)
    call write_hdf5_attr(this%use_dyn_layers, 'use_dyn_layers', loc)
    if(this%use_dyn_layers) then
      call write_hdf5_attr(this%NMaxOctreeLevels, 'NMaxOctreeLevels', loc)
      call write_hdf5_attr(this%LeavesTimeRatio, 'LeavesTimeRatio', loc)
    endif
    call write_hdf5_attr(this%use_pr, 'ParticlesRedistribution', loc)
  endif
  call write_hdf5_attr(this%debug_level, 'debug_level', loc)
  call write_hdf5_attr(this%dt_out, 'dt_out', loc)
  call write_hdf5_attr(this%basename, 'basename', loc)
  call write_hdf5_attr(this%GeometryFile, 'GeometryFile', loc)
  call write_hdf5_attr(this%ReferenceFile, 'ReferenceFile', loc)
  call write_hdf5_attr(this%restart_from_file, 'restart_from_file', loc)
  if(this%restart_from_file) then
    call write_hdf5_attr(this%restart_file, 'restart_file', loc)
    call write_hdf5_attr(this%reset_time, 'reset_time', loc)
  endif

end subroutine save_sim_param

end module mod_sim_param
