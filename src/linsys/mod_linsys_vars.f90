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



!> Module to expose only the linsys type
!!
!! needed to break circular dependencies

module mod_linsys_vars


use mod_param, only: &
  wp

!----------------------------------------------------------------------

implicit none

public :: t_linsys

private

!----------------------------------------------------------------------
!> Type containing all the data necessary to solve the linear system
!! for the aerodynamic elements doublets intensities
type :: t_linsys
 
 !> Rank of the linear system
 integer :: rank
 
 !> Linear system matrix
 real(wp), allocatable :: A(:,:)

 !> Linear system right hand side
 real(wp), allocatable :: b(:)
 
 !> Static part of the right hand side
 !!
 !! The right hand side contains contributions from all the surface panels.
 !! The contributions from surface panels that do not actually move can be 
 !! calculated just once, than multiplied for the freestream velocity and 
 !! finally the contribution of the moving panels is added 
 real(wp), allocatable :: b_static(:,:)

 !> Static part of the lifting line contribution
 !!
 !! The contribution of the lifting line is similar to wake contribution,
 !! however a part could be static on static and can be pre-computed
 real(wp), allocatable :: L_static(:,:)

 !> Result of the linear system solution (intensity of the panels doublets)
 real(wp), allocatable :: res(:)

 !> Number of static and moving panels
 integer :: nstatic, nmoving

 !> Number of static and moving lifting lines
 integer :: nstatic_ll, nmoving_ll, n_ll

end type t_linsys

!----------------------------------------------------------------------

end module mod_linsys_vars
