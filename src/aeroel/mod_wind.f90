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
!!
!!          Andrea Colli          <andrea.colli@polimi.it>
!!
!!=========================================================================

module mod_wind

use mod_param, only: &
  wp, pi
  
use mod_sim_param, only: &
  sim_param

implicit none

contains

function variable_wind(pos, time) result(wind) 

real(wp), intent(in) :: pos(3)
real(wp), intent(in) :: time

real(wp) :: wind(3)

real(wp) :: gust_origin(3), gust_direction(3), gust_u_ds, gust_gradient

real(wp) :: s

!TODO implement reference frames

wind = sim_param%u_inf

if (sim_param%use_gust) then
  select case(trim(sim_param%GustType))
    case('ACM')
      gust_origin = sim_param%gust_origin
      gust_direction = sim_param%gust_direction
      gust_u_ds = sim_param%gust_u_ds
      gust_gradient = sim_param%gust_gradient
      
      s = sum((pos-(gust_origin+sim_param%u_inf*time))*sim_param%u_inf)/norm2(sim_param%u_inf)
      
      wind = wind + gust_u_ds/2*(1-COS(pi*s/gust_gradient))
    
    case('linear')
      gust_origin = sim_param%gust_origin
      gust_direction = sim_param%gust_direction
      gust_u_ds = sim_param%gust_u_ds
      gust_gradient = sim_param%gust_gradient
      
      s = sum((pos-(gust_origin+sim_param%u_inf*time))*sim_param%u_inf)/norm2(sim_param%u_inf)
      
      wind = wind + s*gust_u_ds/gust_gradient*(/0.0, 0.0, 1.0/)
    case default
  end select
end if

end function variable_wind


end module mod_wind
