!  QuAC - Quantum Atomic Collisions
!  Copyright (C) 2009, 2010, 2011, 2012  Yuri Aoto
!
!  This program is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!  
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Module for calculating the 
!
!
Module ModSystem

  USE ModUtil, only : real_kind, error_msg, Error, pi, one, zero, units_len
  USE ModPot, only : Potential, pot

  IMPLICIT NONE

  private

  public :: TypeParticles, particles, lambda, Sec_Deriv, Deriv, Sec_Deriv_tilde, r_to_y_meshkov, y_to_r_meshkov

  TYPE TypeParticles
     REAL(kind = real_kind) :: mass1=-one, mass2=-one, reduced_mass=-one
     CHARACTER(LEN = 7) :: atom1, atom2
     CHARACTER(LEN = units_len) :: mass_uni = "me"
     CHARACTER(LEN = 15) :: name1,name2
     INTEGER :: Z1,Z2
     REAL(kind = real_kind) :: abundance1,abundance2
  END TYPE TypeParticles

  TYPE(TypeParticles), POINTER :: particles=>NULL()

  INTEGER :: lambda

CONTAINS
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  Real(kind = real_kind) Function Sec_Deriv(J,r,wave_vector,energy,already_calc_pot)

    IMPLICIT NONE

    Integer, Intent(in) :: J
    Real(kind = real_kind), Intent(in) :: r
    Real(kind = real_kind), intent(in), optional :: energy, wave_vector
    real(kind = real_kind), intent(in), optional :: already_calc_pot

    real(kind = real_kind) :: one_or_minus_one
    
    if(.not.present(energy).or.energy.ge.zero) then
       one_or_minus_one = -one ! Positive values of energy
    else
       one_or_minus_one = one
    end if

    if(present(energy).and.present(wave_vector))then
       if(present(already_calc_pot))then
          Sec_Deriv = (J - lambda**2)/(r**2) + already_calc_pot/abs(energy) + one_or_minus_one
       else
          Sec_Deriv = (J - lambda**2)/(r**2) + Potential(r/wave_vector)/abs(energy) + one_or_minus_one
       end if
    else
       if(present(already_calc_pot))then
          Sec_Deriv = 2*Particles%reduced_mass*already_calc_pot + (J - lambda**2)/r**2
       else
          Sec_Deriv = 2*Particles%reduced_mass*Potential(r) + (J - lambda**2)/r**2
       end if
    end if

  END Function Sec_Deriv
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  Subroutine Deriv(func,J,r,wave_vector,energy,diff,already_calc_pot)

    IMPLICIT NONE

    Real(kind = real_kind), Dimension(2), Intent(in) :: func
    Integer, Intent(in) :: J
    Real(kind = real_kind), Intent(in) :: r
    Real(kind = real_kind), intent(in), optional :: wave_vector, energy
    Real(kind = real_kind), Dimension(2), Intent(out) :: diff
    real(kind = real_kind), intent(in), optional :: already_calc_pot

    diff(1) = func(2)
    diff(2) = Sec_Deriv(J,r,wave_vector,energy,already_calc_pot)*func(1)

  END Subroutine Deriv
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  Real(kind = real_kind) Function Sec_Deriv_tilde(J,Y,r_bar,beta,already_calc_pot)

    IMPLICIT NONE

    Integer, Intent(in) :: J
    Real(kind = real_kind), Intent(in) :: Y, r_bar, beta
    real(kind = real_kind), intent(in), optional :: already_calc_pot

    Real(kind = real_kind) :: g

! calculating from r
!    g = pi*r_bar/(2*beta)*(1 + (beta*(r/r_bar - 1))**2)
!    Sec_Deriv_tilde = -(g**2)*Sec_Deriv(J,r) + pi**2/4

! calculating from y
    g = pi*r_bar/(2*beta)*(1 + tan(pi*Y/2)**2)
    Sec_Deriv_tilde = -(g**2)*Sec_Deriv(J,y_to_r_meshkov(Y,r_bar,beta),already_calc_pot=already_calc_pot) + pi**2/4

  End Function Sec_Deriv_tilde
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  Real(kind = real_kind) Function r_to_y_meshkov(r,r_bar,beta)

    IMPLICIT NONE

    Real(kind = real_kind), Intent(in) :: r, r_bar, beta

    R_to_y_meshkov = 2*atan(beta*(r/r_bar - 1))/pi

  End Function r_to_y_meshkov
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  Real(kind = real_kind) Function y_to_r_meshkov(y,r_bar,beta)

    IMPLICIT NONE

    Real(kind = real_kind), Intent(in) :: Y, r_bar, beta
    Real(kind = real_kind) :: Ymod

    if(Y.gt.one)then
       Ymod = one
    else
       Ymod = Y
    end if

    Y_to_r_meshkov = r_bar*(1 + tan(pi*Ymod/2)/beta)

  End Function Y_to_r_meshkov
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
END Module ModSystem
