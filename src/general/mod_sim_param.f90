!./\\\\\\\\\\\...../\\\......./\\\..../\\\\\\\\\..../\\\\\\\\\\\\\.
!.\/\\\///////\\\..\/\\\......\/\\\../\\\///////\\\.\//////\\\////..
!..\/\\\.....\//\\\.\/\\\......\/\\\.\//\\\....\///.......\/\\\......
!...\/\\\......\/\\\.\/\\\......\/\\\..\////\\.............\/\\\......
!....\/\\\......\/\\\.\/\\\......\/\\\.....\///\\...........\/\\\......
!.....\/\\\......\/\\\.\/\\\......\/\\\.......\///\\\........\/\\\......
!......\/\\\....../\\\..\//\\\...../\\\../\\\....\//\\\.......\/\\\......
!.......\/\\\\\\\\\\\/....\///\\\\\\\\/..\///\\\\\\\\\/........\/\\\......
!........\///////////........\////////......\/////////..........\///.......
!!=========================================================================
!!
!! Copyright (C) 2018-2022 Politecnico di Milano,
!!                           with support from A^3 from Airbus
!!                    and  Davide   Montagnani,
!!                         Matteo   Tugnoli,
!!                         Federico Fonte
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
!!          Federico Fonte
!!          Davide Montagnani
!!          Matteo Tugnoli
!!=========================================================================


!> Module to define the structrues containing the simulation parameters
!!
module mod_sim_param

use mod_param, only: &
  wp, max_char_len

use mod_hdf5_io, only: &
   h5loc, &
   write_hdf5_attr

use mod_parse, only: &
  t_parse, &
  countoption , &
  getstr, getlogical, getreal, getint, getrealarray, getintarray, &
  ignoredParameters, finalizeParameters

implicit none

public :: t_sim_param, sim_param

private

type t_sim_param

  !For debugging purpose
  character(max_char_len) :: basename_debug

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
  !> Actual time
  real(wp) :: time
  !> Previous time
  real(wp) :: time_old
  !> ndt between 2 wake updates
  integer :: ndt_update_wake

  !> Output detailed geometry each timestep
  logical :: output_detailed_geo

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
  !> Join close trailing edges
  logical :: join_te
    !> All trailing edges closer than join_te_factor will be joined
    real(wp) :: join_te_factor

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
    !> use the vortex stretching from elements
    logical :: vs_elems
    !> use the divergence filtering
    logical :: use_divfilt
      !> time scale of the divergence filter
      real(wp) :: filt_eta
  !> use the vorticity diffusion or not
  logical :: use_vd
  !> use turbulent viscosity or not
  logical :: use_tv
  !> use the penetration avoidance
  logical :: use_pa
    !> check radius for penetration avoidance
    real(wp) :: pa_rad_mult
    !> element impact radius for penetration avoidance
    real(wp) :: pa_elrad_mult
  !> simulate viscosity effects or not
  logical :: use_ve

  !Lifting Lines
  character(len=max_char_len) :: llSolver
  !> Reynolds corrections of .c81 tables
  logical  :: llReynoldsCorrections
  !> n factor for Reynolds corrections of .c81 tables: (Re/Re_T)^n
  real(wp) :: llReynoldsCorrectionsNfact
  !> Maximum number of iteration in LL algorithm
  integer  :: llMaxIter
  !> Tolerance for the relative error in fixed point iteration for LL
  real(wp) :: llTol
  !> Damping param in fixed point iteration for LL used to avoid oscillations
  real(wp) :: llDamp
  !> Avoid "unphysical" separations in inner sections of LL? :: llTol
  logical  :: llStallRegularisation
  !> Number of "unphysical" separations that can be removed
  integer  :: llStallRegularisationNelems
  !> Number of iterations between two regularisation processes
  integer  :: llStallRegularisationNiters
  !> Reference stall AOA for regularisation
  real(wp) :: llStallRegularisationAlphaStall
  !> Constant Artificial Viscosity for regularisation
  real(wp) :: llArtificialViscosity
  !> Adaptive Artificial Viscosity algorithm
  logical  :: llArtificialViscosityAdaptive
  !> Adaptive Artificial Viscosity algorithm, reference AOA
  real(wp) :: llArtificialViscosityAdaptive_Alpha
  !> Adaptive Artificial Viscosity algorithm, blending interval
  ! between full regularisation and no regularisation ( AOA )
  real(wp) :: llArtificialViscosityAdaptive_dAlpha
  !> Use AVL expression for inviscid load computation ( ~ VL )
  logical  :: llLoadsAVL

  !FMM parameters
  !> Employing the FMM method
  logical :: use_fmm
    !> Employing the FMM method also for panels
    logical :: use_fmm_pan
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
      !> Level at which is checked the presence of panels
      integer :: lvl_solid
      real(wp) :: part_redist_ratio

  !HCAS parameters
  !> Use hcas
  logical :: hcas
    !> Time of deployment of the hcas
    real(wp) :: hcas_time
    !> Velocity of the hcas
    real(wp) :: hcas_vel(3)


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

  !Variable wind
  !> Gust
  logical :: use_gust
  !> Gust type
  character(len=max_char_len) :: GustType
  !> Gust parameters
  real(wp) :: gust_origin(3)
  real(wp) :: gust_front_direction(3)
  real(wp) :: gust_front_speed
  real(wp) :: gust_u_des
  real(wp) :: gust_perturb_direction(3)
  real(wp) :: gust_gradient
  real(wp) :: gust_time

contains

  procedure, pass(this) :: save_param => save_sim_param
end type t_sim_param

type(t_sim_param) :: sim_param

!----------------------------------------------------------------------
contains
!----------------------------------------------------------------------


subroutine save_sim_param(this, loc)
 class(t_sim_param) :: this
 integer(h5loc), intent(in) :: loc

  call write_hdf5_attr(this%t0, 't0', loc)
  call write_hdf5_attr(this%dt, 'dt', loc)
  call write_hdf5_attr(this%tend, 'tend', loc)
  call write_hdf5_attr(this%output_detailed_geo, 'output_detailed_geo', loc)
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
  call write_hdf5_attr(this%join_te, 'join_te', loc)
  if(this%join_te) &
    call write_hdf5_attr(this%join_te_factor, 'join_te_factor', loc)


  call write_hdf5_attr(this%FarFieldRatioDoublet, 'FarFieldRatioDoublet', loc)
  call write_hdf5_attr(this%FarFieldRatioSource, 'FarFieldRatioSource', loc)
  call write_hdf5_attr(this%DoubletThreshold, 'DoubletThreshold', loc)
  call write_hdf5_attr(this%RankineRad, 'RankineRad', loc)
  call write_hdf5_attr(this%VortexRad, 'VortexRad', loc)
  call write_hdf5_attr(this%CutoffRad, 'CutoffRad', loc)
  call write_hdf5_attr(this%use_vs, 'Vortstretch', loc)
  if(this%use_vs) then
    call write_hdf5_attr(this%vs_elems, 'VortstretchFromElems', loc)
    call write_hdf5_attr(this%use_divfilt, 'DivergenceFiltering', loc)
    call write_hdf5_attr(1.0_wp/this%filt_eta*this%dt, 'FilterTimescale', loc)
  endif
  call write_hdf5_attr(this%use_vd, 'vortdiff', loc)
  call write_hdf5_attr(this%use_tv, 'turbvort', loc)
  call write_hdf5_attr(this%use_pa, 'PenetrationAvoidance', loc)
  if(this%use_pa) then
    call write_hdf5_attr(this%pa_rad_mult, 'PenetrationAvoidanceCheckRadius', loc)
    call write_hdf5_attr(this%pa_elrad_mult, 'PenetrationAvoidanceElementRadius', loc)
  endif
  call write_hdf5_attr(this%use_ve, 'ViscosityEffects', loc)
  call write_hdf5_attr(this%use_fmm, 'use_fmm', loc)
  if(this%use_fmm) then
    call write_hdf5_attr(this%use_fmm_pan, 'use_fmm_panels', loc)
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
    if(this%use_pr) then
      call write_hdf5_attr(this%lvl_solid, 'OctreeLevelSolid', loc)
      call write_hdf5_attr(this%part_redist_ratio, &
                                          'ParticlesRedistributionRatio', loc)
    endif
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
  call write_hdf5_attr(this%hcas, 'HCAS', loc)
  if(this%hcas) then
    call write_hdf5_attr(this%hcas_time, 'HCAS_time', loc)
    call write_hdf5_attr(this%hcas_vel, 'HCAS_velocity', loc)
  endif

end subroutine save_sim_param

end module mod_sim_param
