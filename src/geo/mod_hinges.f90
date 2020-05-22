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
!! Copyright (C) 2018-2020 Davide   Montagnani, 
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
!!          Federico Fonte             <federico.fonte@outlook.com>
!!          Davide Montagnani       <davide.montagnani@gmail.com>
!!          Matteo Tugnoli                <tugnoli.teo@gmail.com>
!!=========================================================================

!> Module for introducing hinges and rotating parts in a component
module mod_hinges

use mod_param, only: &
  wp, max_char_len

implicit none

public :: t_hinge

private

! ---------------------------------------------------------------
!> Hinge node configurations: to be used in defining reference
! and actual configurations
type :: t_hinge_config

  !> Position of the first and last point of the hinge in the local
  ! reference frame of the component
  real(wp) :: rr0(3), rr1(3)
  !> Local coordinates of the hinge nodes
  real(wp), allocatable :: rr(:,:)

  !> Unit vectors of the hinge node reference frames
  ! User input in local reference frame, for reference configuration: h, v
  !    h,
  !    n = cross(v,h) , n = n/norm2(n)
  !    v = cross(h,n)
  !> h: rotation axis
  real(wp), allocatable :: h(:,:)
  !> v: zero direction (user inputs are overwritten, in order to 
  ! build ortonormal reference frames
  real(wp), allocatable :: v(:,:)
  !> n: normal direction
  real(wp), allocatable :: n(:,:)

end type t_hinge_config

!> Hinge type
type :: t_hinge

  !> Type: constant, function, coupling
  character(len=max_char_len) :: input_type

  !> N.nodes of the hinge
  integer :: n_nodes

  !> Offset for avoiding irregular behavior
  real(wp) :: offset

  !> Array of the rotation angle (read as an input, or from coupling nodes)
  real(wp), allocatable :: theta(:)

  !> Reference configuration
  type(t_hinge_config) :: ref
  !> Actual configuration
  type(t_hinge_config) :: act

  contains
 
  procedure, pass(this) :: from_reference_to_actual_config

end type t_hinge
! ---------------------------------------------------------------

contains

subroutine from_reference_to_actual_config(this)
  class(t_hinge), intent(inout) :: this

end subroutine from_reference_to_actual_config


end module mod_hinges
