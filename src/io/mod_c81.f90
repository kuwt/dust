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


!> Module to treat the c81 tables
module mod_c81

use mod_param, only: &
  wp, max_char_len

implicit none

!-----------------------------------

! public and private

!-----------------------------------

type t_aero_par_val
  real(wp) , allocatable :: val(:)
end type t_aero_par_val

!-----------------------------------

type t_aero_par
  integer :: n_par
  type(t_aero_par_val) , allocatable :: par(:)
end type t_aero_par

!-----------------------------------

type t_aero_2d_tab

 !> value of the third parameter, Re
 real(wp) :: Re

 !> cl(a,M,Re) 3d-data (extrapolation on Re or no dependence from Re)
 real(wp) , allocatable :: cl(:,:)
 type(t_aero_par) :: cl_par

 !> cd(a,M,Re) 3d-data
 real(wp) , allocatable :: cd(:,:)
 type(t_aero_par) :: cd_par

 !> cm(a,M,Re) 3d-data
 real(wp) , allocatable :: cm(:,:)
 type(t_aero_par) :: cm_par

end type t_aero_2d_tab

!-----------------------------------

!> Tables conatining aerodynamic coefficients
!!
!! Tabulated data for lifting line elements 
type t_aero_tab

 !> airfoil file
 character(len=max_char_len) :: airfoil_file 

 !> airfoil data id
 integer :: id

 !> tables cl(a,M) for each Re number
 type(t_aero_2d_tab) , allocatable :: aero_coeff(:)

end type t_aero_tab

!-----------------------------------

contains

!-----------------------------------
! Vahana input file style '../vahana_input/vahana4.inp'
! HP: the aerodynamic ocefficients are given for the same values
! of al_i, M_i
subroutine read_c81_table ( filen , coeff )
  character(len=max_char_len) , intent(in) :: filen
  type(t_aero_tab) , intent(inout) :: coeff

  character(len=max_char_len) :: line , string
  integer :: cl1 , cl2 , cd1 , cd2 , cm1 , cm2 !  table_dim(2)
  real(wp) :: dummy
  integer :: dummy_int , iblnk
  integer :: fid , i1 , iRe , nRe
  real(wp) :: Re

  integer , parameter :: n_coeff = 3   ! cl , cd , cm

  fid = 21
  open(unit=fid,file=trim(adjustl(filen))) 
  ! First tree lines containing unintelligilbe parameters
  read(fid,*) nRe , dummy_int , dummy_int ! 2 0 0
  read(fid,*) ! 0 1
  read(fid,*) ! 0.158 0.158

  allocate( coeff%aero_coeff(nRe) )

  do iRe = 1 , nRe

    ! For each Reynolds number (1 Reynolds number, 3 coefficients)
    read(fid,*) ! COMMENT#1
    read(fid,*) Re  , dummy
    read(fid,'(A)') line
      iblnk = index(line,' ') 
      string = trim(adjustl(line(iblnk:)))
      read(string,'(6I2)') cl2 , cl1 , cd2 , cd1 , cm2 , cm1
    write(*,*) ' table_dim = ' , cl2 , cl1 , cd2 , cd1 , cm2 , cm1

    ! Reynolds number
    coeff%aero_coeff(iRe)%Re = Re
    ! cl
    coeff%aero_coeff(iRe)%cl_par%n_par = 2
    allocate(coeff%aero_coeff(iRe)%cl_par%par(2))
    allocate(coeff%aero_coeff(iRe)%cl_par%par(1)%val(cl1))
    allocate(coeff%aero_coeff(iRe)%cl_par%par(2)%val(cl2))
    allocate(coeff%aero_coeff(iRe)%cl(cl1,cl2))
    ! cd
    coeff%aero_coeff(iRe)%cd_par%n_par = 2
    allocate(coeff%aero_coeff(iRe)%cd_par%par(2))
    allocate(coeff%aero_coeff(iRe)%cd_par%par(1)%val(cd1))
    allocate(coeff%aero_coeff(iRe)%cd_par%par(2)%val(cd2))
    allocate(coeff%aero_coeff(iRe)%cd(cd1,cd2))
    ! cm
    coeff%aero_coeff(iRe)%cm_par%n_par = 2
    allocate(coeff%aero_coeff(iRe)%cm_par%par(2))
    allocate(coeff%aero_coeff(iRe)%cm_par%par(1)%val(cm1))
    allocate(coeff%aero_coeff(iRe)%cm_par%par(2)%val(cm2))
    allocate(coeff%aero_coeff(iRe)%cm(cm1,cm2))

    read(fid,*) coeff%aero_coeff(iRe)%cl_par%par(2)%val
!   write(*,*)  coeff%aero_coeff(iRe)%cl_par%par(2)%val
    do i1 = 1 , cl1
      read(fid,*) coeff%aero_coeff(iRe)%cl_par%par(1)%val(i1) , & 
                  coeff%aero_coeff(iRe)%cl(i1,:)
!     write(*,*)  coeff%aero_coeff(iRe)%cl(i1,:)
    end do

    read(fid,*) coeff%aero_coeff(iRe)%cd_par%par(2)%val
!   write(*,*)  coeff%aero_coeff(iRe)%cd_par%par(2)%val
    do i1 = 1 , cd1
      read(fid,*) coeff%aero_coeff(iRe)%cd_par%par(1)%val(i1) , & 
                  coeff%aero_coeff(iRe)%cd(i1,:)
!     write(*,*)  coeff%aero_coeff(iRe)%cd(i1,:)
    end do

    read(fid,*) coeff%aero_coeff(iRe)%cm_par%par(2)%val
!   write(*,*) coeff%aero_coeff(iRe)%cm_par%par(2)%val
    do i1 = 1 , cm1
      read(fid,*) coeff%aero_coeff(iRe)%cm_par%par(1)%val(i1) , & 
                  coeff%aero_coeff(iRe)%cm(i1,:)
!     write(*,*)  coeff%aero_coeff(iRe)%cm(i1,:)
    end do
!   write(*,*)

  end do

! write(*,*) line(iblnk:) 
! write(*,*) table_dim  
  
  close(fid)



end subroutine read_c81_table

!-----------------------------------

!-----------------------------------


end module mod_c81
