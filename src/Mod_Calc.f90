!  QuAC - Quantum Atomic Collisions
!  Copyright (C) 2009, 2010, 2011, 2012  Yuri Aoto
!
!  This file is part of QuAC.
!
!  QuAC is free software: you can redistribute it and/or modify
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
! Module for the calculations
!
! Any unit check or conversion is carried in the module.
! The Integration subroutine integrate the radial Shrödinger equation, using the parameters passed in an TypeIntParameter
! variable. Its internal subprograms take car of integration steps and tests over the wave function. Calc_Cross_Sec calculate
! the (total) cross section by partial-wave method. It allow non "continuous" values of l. Calc_Partial_Cross_Sec calculate
! the partial cross section of a given defined angular momentum. Calc_Dif_Cross_Sec calculate the differential cross section,
! in the center of mass or laboratory reference system. Calc_Phase_Shift calculate the phase shift by linear regression over
! the given (final) wave function points. Calc_Scatt_Length calculate the scattering length, by direct linear regression over 
! the given (final) wave function points or by Meshkov procedure.
!
!
Module ModCalc

  USE ModUtil, only : real_kind, pi, zero, error_msg, uni_log, Error, quad_met, one, two, three, Sph_Bessel_J, Sph_Bessel_N
  USE ModPot, only: TypeIntegerLinkedList, Potential, Pot, formats
  USE ModSystem, only: particles, Sec_Deriv, Deriv, Sec_Deriv_tilde, r_to_y_meshkov, y_to_r_meshkov

  IMPLICIT NONE

  private

  public :: TypeIntParameters, Integration, Calc_Cross_Sec, Calc_Dif_Cross_Sec,&
       Calc_Phase_Shift, Calc_Scatt_Length, Analytical_LenJonN_Nm2_Scatt_Length, Analytical_Square_Phase_Shift,&
       Analytical_HardSphere_Phase_Shift

  Type TypeIntParameters
     CHARACTER(LEN = 15) :: method                              ! Numerical method - RK4, ABM, NUMEROV, JOHNSON
     LOGICAL :: save_last_deriv                                 ! save the derivative in the last point
     Real(kind = real_kind) :: step_size                        ! Step size
     Real(kind = real_kind) :: r_min                            ! Intitial point of the integration
     Real(kind = real_kind), Dimension(2) :: initial_wf         ! Wave function and its derivative at r_min
     Real(kind = real_kind) :: r_bar_mesh, beta_mesh            ! Meshkov parameters, to change from ordinary to log derivative variable

     ! Iterations to increase the step_size, to end this increasing, to start saving points, and to end the integration
     Integer :: inc_step, end_inc_step, start_save, end_integration

     Integer :: save_at                      ! Save points, after start_save, at save_at steps
     LOGICAL :: save_points                  ! Print points saved in log file
     Integer :: print                        ! Print steps at each print steps. Don't print if is negative
  END Type TypeIntParameters

CONTAINS
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  Subroutine Integration(Arg,l,energy,wave_vector,points,IntSq,no_knots)

    IMPLICIT NONE

    TYPE(TypeIntParameters), Intent(in) :: Arg
    Integer, Intent(in) :: l
    Real(kind = real_kind), Intent(in), optional :: energy, wave_vector
    Real(kind = real_kind), Allocatable, Dimension(:,:), Intent(out), optional :: points
    Real(kind = real_kind), Intent(out), optional :: IntSq
    Integer, Intent(out), optional :: no_knots

    Integer :: i_iter, i_points, J, increase_step, last_non_inter_pt, ini_simp_rule
    Real(kind = real_kind) :: step_size, r, r_at_uni, r1_at_uni, r2_at_uni, func_tmp, func_prev, pot_next_1, pot_next_2
    Real(kind = real_kind), Dimension(2) :: func,func_ant
    Real(kind = real_kind), Dimension(1:4,2) :: func_4_steps
    real(kind = real_kind), Dimension(-4:0) :: pot_vector

    if((.not.present(energy).and.present(wave_vector)).or.(present(energy).and..not.present(wave_vector)))then
       error_msg = "Passed just one of energy/wave_vector in integration subroutine"
       CALL Error(2)
    end if

    if(present(points))then
       If(Arg%method.eq."JOHNSON")then
          Allocate(points(1,2))
       else
          if(arg%save_last_deriv)then
             Allocate(points(2,2))
          else
             Allocate(points((Arg%end_integration-Arg%start_save)/Arg%save_at+1,2))
          end if
       end If
    end if

!!$    if(Arg%save_global_potential)then
!!$       save_Potential = .true.
!!$       select case(Arg%method)
!!$       case("RK4")
!!$          Allocate(Global_Pot_vector(0:2*Arg%end_integration))
!!$       case("ABM")
!!$          Allocate(Global_Pot_vector(0:Arg%end_integration+3))
!!$       case default
!!$          Allocate(Global_Pot_vector(0:Arg%end_integration))
!!$       end select
!!$    else
!!$       save_Potential = .false.
!!$    end if
!!$    if(Arg%use_global_potential)then
!!$       use_Pot_vector = .true.
!!$    else
!!$       use_Pot_vector = .false.
!!$    end if

    func_ant(1) = zero
    func_ant(2) = zero
    if(present(IntSq))then
       IntSq = zero
    end if
    if(arg%save_last_deriv)then
       last_non_inter_pt = Arg%end_integration
    else
       last_non_inter_pt = Arg%end_integration-2
    end if

    ini_simp_rule = mod(last_non_inter_pt,2)

    if(present(no_knots)) no_knots=0

    increase_step = Arg%inc_step
    step_size = Arg%step_size

    i_points = 1
    i_iter = 0
    J = l*(l + 1)

    func = Arg%initial_wf
    r = Arg%r_min

    if(Arg%method.eq."JOHNSON")then
       r_at_uni = y_to_r_meshkov(r,Arg%r_bar_mesh,Arg%beta_mesh)
       r1_at_uni = y_to_r_meshkov(r+step_size,Arg%r_bar_mesh,Arg%beta_mesh)
       r2_at_uni = y_to_r_meshkov((r+step_size)+step_size,Arg%r_bar_mesh,Arg%beta_mesh)
    else
       r_at_uni = r
       r1_at_uni = r+step_size
       r2_at_uni = (r+step_size)+step_size ! to ensure equivalence for value of r after two steps

       if(present(energy)) r_at_uni = r_at_uni/wave_vector
       if(present(energy)) r1_at_uni = r1_at_uni/wave_vector
       if(present(energy)) r2_at_uni = r2_at_uni/wave_vector

    end if

    pot_vector(0) = Potential(r_at_uni)
    pot_next_1 = Potential(r1_at_uni)
    pot_next_2 = Potential(r2_at_uni)

    ! Print first point
    If(Arg%print.GE.0)Then

       Write(uni_log,*) 'Numerical integration - ',Trim(Arg%method)
       Write(uni_log,*)
       
       SELECT CASE(Arg%method)
          
       CASE("RK4","ABM","NUMEROV")

          if(present(energy))then

             Write(uni_log,'("distance (dimension less)   distance (a0)   radial function   radial potential ",'//&
                  '"  effective potential")')
             Write(uni_log,FMT='('//TRIM(formats%distance)//'," ",'//TRIM(formats%distance)//'," ",'&
                  //TRIM(formats%wave_function)//'," ",'//TRIM(formats%potential)//'," ",'//TRIM(formats%potential)//')')&
                  r,r_at_uni,func(1),pot_vector(0),pot_vector(0)+J/r_at_uni**2
          else
             
             Write(uni_log,FMT='("distance   radial function   radial potential   effective potential")')
             Write(uni_log,FMT='('//TRIM(formats%distance)//'," ",'//TRIM(formats%wave_function)//'," ",'&
                  //TRIM(formats%potential)//'," ",'//TRIM(formats%potential)//')')&
                  r,func(1),pot_vector(0),pot_vector(0)+J/r**2
          end if
          
       CASE("JOHNSON")
          
          Write(uni_log,FMT='(" distance (modified, in [a,1])   logaritmic derivative   distance (original)   radial potential",'//&
               '"effective potential")')
          Write(uni_log,FMT='('//TRIM(formats%distance)//'," ",'//TRIM(formats%wave_function)//'," ",'&
               //TRIM(formats%potential)//'," ",'//TRIM(formats%potential)//'," ",'//TRIM(formats%potential)//')') &
               r,func(1),r_at_uni,pot_vector(0),pot_vector(0)+J/r_at_uni**2
          
       CASE DEFAULT
          
          error_msg = "Inconsistency of integration metod in printing first point: "//TRIM(Arg%method)
          CALL Error(2)
          
       END SELECT

    End If

    CALL Tests ()

    ! Integration
    Do i_iter = 1,Arg%end_integration,1

       r = r+step_size ! Step where the wave function will be calculated
       if(Arg%method.eq."JOHNSON")then
          r_at_uni = y_to_r_meshkov(r,Arg%r_bar_mesh,Arg%beta_mesh)
       else
          r_at_uni = r
          if(present(energy)) r_at_uni = r_at_uni/wave_vector
       end if

       pot_vector(-4) = pot_vector(-3)
       pot_vector(-3) = pot_vector(-2)
       pot_vector(-2) = pot_vector(-1)
       pot_vector(-1) = pot_vector(0)

       select case (i_iter)
       case(1)
          pot_vector(0) = pot_next_1
       case(2)
          pot_vector(0) = pot_next_2
       case default
          pot_vector(0) = Potential(r_at_uni)
       end select

       SELECT CASE(Arg%method)

          CASE ("NUMEROV")
             ! func(1) -> newest wave function
             ! func(2) -> wave function of the previous step

             func_tmp = func(1)

             if(i_iter.eq.1)then
                func(1) = Numerov_First_Step(Arg%r_min,Arg%initial_wf(1),Arg%initial_wf(2),(/pot_vector(0),pot_next_1,pot_next_2/))
             else
                func(1) = Numerov_Step(r,func(1),func(2),pot_vector(-2:0))
             end if
             
             func(2) = func_tmp
             
          CASE("JOHNSON")
             ! funcao(1) -> logaritmic derivative times integration step, recently calculated

             func(1) = Johnson_Step (r,func(1),pot_vector(0))

          CASE("ABM")
             ! funcao(1) -> newest wave function
             ! funcao(2) -> derivative of the newest wave function
             ! fun4pas(1,:) -> wave function in the previous step
             ! fun4pas(2,:) -> wave function in the step before of the previous step
             ! fun4pas(3,:) -> wave function in the step before of the step before of the previous step
             ! fun4pas(4,:) -> wave function in the step before of the step before of the step before of the previous step

             if(i_iter.le.3)then
                func_4_steps(i_iter,:) = RK4_Step(r,func,pot_vector(-1:0))
                func = func_4_steps(i_iter,:)
             else
                func = ABM_Step(r,func_4_steps,pot_vector(-4:0))
             end if
             
             func_4_steps(4,:) = func_4_steps(3,:)
             func_4_steps(3,:) = func_4_steps(2,:)
             func_4_steps(2,:) = func_4_steps(1,:)
             func_4_steps(1,:) = func

          CASE("RK4")
             ! funcao(1) -> newest wave function
             ! funcao(2) -> derivative of the newest wave function

             func = RK4_step(r,func,pot_vector(-1:0))

          CASE DEFAULT

             error_msg = "Inconsistency of integration method: "//TRIM(Arg%method)
             CALL Error(2)

       END SELECT

       CALL Tests ()

    End Do

    If(Arg%print.GT.0)Then
       Write(uni_log,*)
    EndIf

    ! Print points
    If(Arg%save_points.and.present(points))Then
       Write(uni_log,FMT='("Points saved for phase shift/scattering length calculation:")')
       Write(uni_log,*)
       
       Do i_points=1,size(points,1),1
          Write(uni_log,FMT='('//TRIM(formats%distance)//'," ",'//TRIM(formats%wave_function)//')')&
               points(i_points,1), points(i_points,2)
       End Do
          
       Write(uni_log,*)

    EndIf

  Contains
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
    Real(kind = real_kind) Function Numerov_First_Step(x_ini,y_ini,y_ini_der,al_calc_pot)
      
      Real(kind = real_kind), intent(in) :: x_ini,y_ini,y_ini_der
      real(kind = real_kind), intent(in), dimension(0:2), optional :: al_calc_pot

      Real(kind = real_kind) :: tmp, V0, V1, V2
      
      
      V0 = Sec_Deriv(J,x_ini                  ,wave_vector,energy,al_calc_pot(0))
      V1 = Sec_Deriv(J,x_ini +   step_size,wave_vector,energy,al_calc_pot(1))
      V2 = Sec_Deriv(J,(x_ini + step_size) + step_size,wave_vector,energy,al_calc_pot(2))

      tmp = y_ini*(1 - V2*step_size**2/24) + step_size*y_ini_der*(1 - V2*step_size**2/12)
      tmp = tmp + (step_size**2/24)*7*V0*Arg%initial_wf(1) - (V2*step_size**4/36)*V0*y_ini
      Numerov_First_Step = tmp/(1 - V1*step_size**2/4 + V1*V2*step_size**4/18)

    End Function Numerov_First_Step
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
    Real(kind = real_kind) Function Numerov_Step (x,y_current,y_prev,al_calc_pot)
      
      Real(kind = real_kind), intent(in) :: X,y_current,y_prev
      real(kind = real_kind), intent(in), dimension(-2:0), optional :: al_calc_pot

      Real(kind = real_kind) :: tmp

      tmp = 10*Sec_Deriv(J,X-1*step_size,wave_vector,energy,al_calc_pot(-1))*y_current + &
           &Sec_Deriv(J,X-2*step_size,wave_vector,energy,al_calc_pot(-2))*y_prev
      tmp = 2*y_current - y_prev + (step_size**2/12)*tmp

      Numerov_Step = tmp/(1 - (Sec_Deriv(J,X,wave_vector,energy,al_calc_pot(0))*step_size**2)/12)

    End Function Numerov_Step
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$x$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
    Real(kind = real_kind) Function Johnson_Step (x,y_current,al_calc_pot)
      
      IMPLICIT NONE
      
      Real(kind = real_kind), intent(in) :: x,y_current
      real(kind = real_kind), intent(in), optional :: al_calc_pot

      Real(kind = real_kind) :: Q_tilde, u, w
      
      Q_tilde=Sec_Deriv_tilde(J,x,Arg%r_bar_mesh,Arg%beta_mesh,al_calc_pot)
      if(MOD(i_iter,2).eq.0)then
         u=Q_tilde
         if(i_iter.eq.0.or.i_iter.eq.Arg%end_integration) then
            w=1
         else
            w=2
         end if
      else
         u=Q_tilde/(1+Q_tilde*step_size**2/6)
         w=4
      end if

      Johnson_Step = y_current/(1+y_current) - u*w*step_size**2/3

    End Function Johnson_Step
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$x$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
    Function ABM_Step (x,y_current_4_steps,al_calc_pot)
      
      IMPLICIT NONE
      
      Real(kind = real_kind), Intent(in) :: x
      Real(kind = real_kind), Dimension(4,2), intent(in) :: y_current_4_steps
      real(kind = real_kind), intent(in), dimension(-4:0), optional :: al_calc_pot

      Real(kind = real_kind), Dimension(2) :: ABM_Step
      Real(kind = real_kind), Dimension(4,2) :: diff
      Real(kind = real_kind), Dimension(2) :: func_try

      CALL Deriv(y_current_4_steps(1,:),J,x-1*step_size,wave_vector,energy,diff(1,:),al_calc_pot(-1))
      CALL Deriv(y_current_4_steps(2,:),J,x-2*step_size,wave_vector,energy,diff(2,:),al_calc_pot(-2))
      CALL Deriv(y_current_4_steps(3,:),J,x-3*step_size,wave_vector,energy,diff(3,:),al_calc_pot(-3))
      CALL Deriv(y_current_4_steps(4,:),J,x-4*step_size,wave_vector,energy,diff(4,:),al_calc_pot(-4))
      
      func_try = y_current_4_steps(1,:) + (step_size/24)*(55*diff(1,:) - 59*diff(2,:) + 37*diff(3,:) - 9*diff(4,:))
      
      CALL Deriv(func_try,J,X,wave_vector,energy,Diff(4,:),al_calc_pot(0))

      ABM_Step = y_current_4_steps(1,:) + (step_size/24)*(9*diff(4,:) + 19*diff(1,:) - 5*diff(2,:) + diff(3,:))
      
    End Function ABM_Step
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
    Function RK4_Step (x,y_current,al_calc_pot)
      
      IMPLICIT NONE

      Real(kind = real_kind), Intent(in) :: x
      Real(kind = real_kind), Dimension(2), Intent(in) :: y_current
      real(kind = real_kind), intent(in), dimension(-1:0), optional :: al_calc_pot

      real(kind = real_kind) :: pot_half_way
      Real(kind = real_kind), Dimension(2) :: RK4_Step
      Real(kind = real_kind), Dimension(2) :: K1, K2, K3, K4, diff

      if(present(energy).and.present(wave_vector))then
         pot_half_way = Potential((X-step_size/2)/wave_vector)
      else
         pot_half_way = Potential(X-step_size/2)
      end if

      CALL Deriv(y_current,J,X-step_size,wave_vector,energy,diff,al_calc_pot(-1))
      K1 = step_size*diff
      
      CALL Deriv(y_current+(K1/2),J,X-step_size/2,wave_vector,energy,diff,pot_half_way)
      K2 = step_size*diff
      
      CALL Deriv(y_current+(K2/2),J,X-step_size/2,wave_vector,energy,diff,pot_half_way)
      K3 = step_size*diff
      
      CALL Deriv(y_current+K3,J,X,wave_vector,energy,diff,al_calc_pot(0))
      K4 = step_size*diff

      RK4_Step = func + (K1 + 2*K2 + 2*K3 + K4)/6
      
    End Function RK4_Step
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
    Subroutine Tests()
      
      IMPLICIT NONE

      Real(kind = real_kind) :: coef_quad

      if(i_iter.gt.0)then
         
         If(Arg%print.GT.0.AND.MOD(i_iter,Arg%print).EQ.0)Then
            
            SELECT CASE(Arg%method)
               
            CASE("RK4","ABM","NUMEROV")
               if(present(energy))then
                  Write(uni_log,FMT='('//TRIM(formats%distance)//'," ",'//TRIM(formats%distance)//'," ",'&
                       //TRIM(formats%wave_function)//'," ",'//TRIM(formats%potential)//'," ",'//TRIM(formats%potential)//')')&
                       r,r_at_uni,func(1),pot_vector(0),pot_vector(0)+J/r_at_uni**2
                  
               else
                  Write(uni_log,FMT='('//TRIM(formats%distance)//'," ",'//TRIM(formats%wave_function)//'," ",'&
                       //TRIM(formats%potential)//'," ",'//TRIM(formats%potential)//')')&
                       r,func(1),pot_vector(0),pot_vector(0)+J/r**2
                  
               end if
               
            CASE("JOHNSON")
               
               Write(uni_log,FMT='('//TRIM(formats%distance)//'," ",'//TRIM(formats%wave_function)//'," ",'&
                    //TRIM(formats%potential)//'," ",'//TRIM(formats%potential)//'," ",'//TRIM(formats%potential)//')') &
                    r,func(1),r_at_uni,pot_vector(0),pot_vector(0)+J/r_at_uni**2
               
            CASE DEFAULT
               
               error_msg = "Inconsistency in printing in numerical integration."//TRIM(Arg%method)
               CALL Error(2)
               
            END SELECT
            
         Endif
         
         ! Count the number of knots !!! IT WILL COUNT TWICE THE KNOT IN THE LAST 5 POINTS OF THE INTEGRATION IN BOUND! CORRECT IT!
         If(present(no_knots))then
            if(i_iter.le.last_non_inter_pt) then !Arg%end_integration-2.OR.arg%save_last_deriv) then
               if(func_ant(1)*func(1).lt.zero.or.func(1).eq.zero) no_knots = no_knots+1
            end if
         end If

         func_ant = func
         
         ! Increase step size
         If(i_iter.ge.increase_step.and.i_iter.le.Arg%end_inc_step.and.&
              (Arg%method.eq."NUMEROV"))then
            
            if (i_iter.eq.increase_step) then
               func_prev=func(1)
            End if
            
            if (i_iter.eq.increase_step+2) then
               step_size=2*step_size
               func(2)=func_prev
               func_prev=func(1)
               increase_step=i_iter
            end if
            
         end If
         
      end if

      ! Acumulate quadrature sum
      If(present(IntSq).AND.i_iter.le.last_non_inter_pt)then

         select case(quad_met)

         case("TRAPEZOIDAL")
            if(i_iter.eq.0.or.i_iter.eq.last_non_inter_pt)then
               coef_quad = one/two
            else
               coef_quad = one
            end if

         CASE("SIMPSON")
            
            if(i_iter.eq.last_non_inter_pt.or.&
                 mod(last_non_inter_pt,2).eq.0.and.i_iter.eq.0.or.&
                 mod(last_non_inter_pt,2).eq.1.and.i_iter.eq.1)then
               coef_quad = one/three
               if(i_iter.EQ.1)then
                  coef_quad = coef_quad + one/two ! Second point: First interval by trapezoidal rule
               end if
            elseif(i_iter.GT.ini_simp_rule)then
               if(mod(i_iter+ini_simp_rule,2).eq.0)then
                  coef_quad = two/three
               else
                  coef_quad = two*two/three
               end if
            else
               coef_quad = one/two ! First point: First interval by trapezoidal rule
            end if
            
         CASE DEFAULT
            error_msg = trim(quad_met)//" is not an know quadrature name."
            CALL Error(2)

         end select

         IntSq = IntSq + coef_quad*(abs(step_size/wave_vector))*func(1)**2
         
      end if

      ! Check if save point
      If(present(points).and.ABS(i_iter).GE.Arg%start_save.AND.MOD(i_iter-Arg%start_save,Arg%save_at).EQ.0&
           &.AND.i_points.LE.size(points,1))Then
         if(arg%save_last_deriv)then
            points(1,1) = r
            points(2,1) = r

            points(1,2) = func(1)
            points(2,2) = func(2)
         else
            points(i_points,1) = r
            points(i_points,2) = func(1)
            If(arg%method.eq."JOHNSON")then
               points(i_points,2) = points(i_points,2)/step_size
            end If
         end if

         i_points = i_points + 1

      Endif

    End Subroutine Tests
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  End Subroutine Integration
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  Real(kind = real_kind) Function Calc_Cross_Sec(phase_shift,wave_vector,l_max,l_list,l_single)

    IMPLICIT NONE

    Real(kind = real_kind), Dimension(:), Intent(in) :: phase_shift
    Real(kind = real_kind), intent(in) :: wave_vector
    Integer, Intent(in), optional :: l_max
    TYPE(TypeIntegerLinkedList), intent(in), POINTER, optional :: l_list
    Integer, Intent(in), optional :: l_single

    Integer :: i_l
    TYPE(TypeIntegerLinkedList), POINTER :: current_int => NULL()

!
! Consistency between l_list/lmax and size of phase shift not verified
!

    if(count((/present(l_max),present(l_list),present(l_single)/)).ne.1)then
       error_msg = "in Calc_Cross_Sec one and only one l specification must be passed."
       CALL Error(2)
    end if

    Calc_Cross_Sec = zero

    if(present(l_max))then

       do i_l=0,l_max,1
          Calc_Cross_Sec = Calc_Cross_Sec + (2*i_l + 1)*sin(phase_shift(i_l+1))**2
       End do

    elseif(present(l_list))then

       current_int=>l_list
       i_l=0
       do while(associated(current_int))
          i_l=i_l+1
          Calc_Cross_Sec = Calc_Cross_Sec + (2*current_int%value + 1)*sin(phase_shift(i_l))**2
          current_int=>current_int%next
       end do
       
    else

       Calc_Cross_Sec = Calc_Cross_Sec + (2*l_single + 1)*sin(phase_shift(1))**2

    end if
    
    Calc_Cross_Sec = 4*pi*(1/wave_vector**2)*Calc_Cross_Sec

!!$    do i=0,size(l),1
!!$       Calc_Cross_Sec = Calc_Cross_Sec + Partial_Cross_Sec(phase_shift(i),l(i),wave_vector)
!!$    End do
    

  END Function Calc_Cross_Sec
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
!!$  Real(kind = real_kind) Function Calc_Partial_Cross_Sec(phase_shift,l,wave_vector)
!!$
!!$    IMPLICIT NONE
!!$
!!$    Real(kind = real_kind), Intent(in) :: phase_shift
!!$    Integer, intent(in) :: l
!!$    Real(kind = real_kind), intent(in) :: wave_vector
!!$
!!$    Calc_Partial_Cross_Sec = 4*pi*(1/wave_vector**2)*(2*l + 1)*sin(phase_shift)**2
!!$
!!$  END Function Calc_Partial_Cross_Sec
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  Real(kind = real_kind) Function Calc_Dif_Cross_Sec(phase_shift,angle,wave_vector,is_center_of_mass_ref)

    IMPLICIT NONE

    Real(kind = real_kind), Dimension(0:), Intent(in) :: phase_shift
    Real(kind = real_kind), Intent(in) :: angle
    Real(kind = real_kind), Intent(in) :: wave_vector
    Logical, Intent(in) :: is_center_of_mass_ref

    Real(kind = real_kind) :: cs_real, cs_imag
    Real(kind = real_kind) :: ang_rad, cos_ang
    Real(kind = real_kind) :: legendre, legendre_0, legendre_1
    Real(kind = real_kind) :: tau, A
    Integer :: l

    ang_rad = angle*pi/180
    A = cos(ang_rad)

    If(is_center_of_mass_ref)Then
       cos_ang = cos(ang_rad)
    Else
       tau = MIN(Particles%mass1/Particles%mass2,Particles%mass2/Particles%mass1)
       cos_ang = tau*(A**2-1) + A*sqrt(((A*tau)**2) - (tau**2) + 1)
    End If

    cs_real = 0
    cs_imag = 0

    legendre_0 = 1
    legendre_1 = cos_ang

    cs_real = cs_real + (2*0 + 1)*cos(phase_shift(0))*sin(phase_shift(0))*legendre_0
    cs_imag = cs_imag + (2*0 + 1)*sin(phase_shift(0))*sin(phase_shift(0))*legendre_0

    If(size(phase_shift).NE.0)Then

       cs_real = cs_real + (2*1 + 1)*cos(phase_shift(1))*sin(phase_shift(1))*legendre_1
       cs_imag = cs_imag + (2*1 + 1)*sin(phase_shift(1))*sin(phase_shift(1))*legendre_1
       
       If(size(phase_shift).NE.1)Then
          
          do l=2,size(phase_shift),1
             
             legendre = (cos_ang*(2*l - 1)*legendre_1 - (l - 1)*legendre_0)/l
             
             cs_real = cs_real + (2*l + 1)*cos(phase_shift(l))*sin(phase_shift(l))*legendre
             cs_imag = cs_imag + (2*l + 1)*sin(phase_shift(l))*sin(phase_shift(l))*legendre
             
             legendre_0 = legendre_1
             legendre_1 = legendre
             
          End do

       End If

    End If

    cs_real = cs_real/wave_vector
    cs_imag = cs_imag/wave_vector

    Calc_Dif_Cross_Sec = cs_real**2 + cs_imag**2 ! Esta é a seção de choque no sistema do centro de massas

    If(is_center_of_mass_ref)Then
       Calc_Dif_Cross_Sec = Calc_Dif_Cross_Sec
    Else
       Calc_Dif_Cross_Sec = Calc_Dif_Cross_Sec*sqrt((1 + (2*tau*cos_ang) + (tau**2))**3)/ABS(1 + tau*cos_ang)
    End If

  END Function Calc_Dif_Cross_Sec
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  Real(kind = real_kind) Function Calc_Phase_Shift(l,points,wv)   !energy,wv)

    IMPLICIT NONE

    Integer, Intent(in) :: l
    Real(kind = real_kind), Dimension(:,:), Intent(in) :: points
!    Real(kind = real_kind), Intent(in) :: energy
    Real(kind = real_kind), Intent(in) :: wv
!    Logical, Intent(in) :: Integracao

!    Type(PontosIgEsp) :: integrando                          ! Integrando

    Real(kind = real_kind), Allocatable, Dimension(:,:) :: points_mod
    Real(kind = real_kind) :: sum_x, sum_x_sqr
    Real(kind = real_kind) :: sum_x_y, sum_y

!    Real(kind = real_kind) :: ValMAX     ! Valor máximo da função no final da integração. Diferencas para encontrar ValMAX

    integer :: i                                             ! Contador

    If(wv.lt.0)then
       error_msg = "Não calcula-se deslocamento de fase para energia nula por esta maneira (vetor de onda negativo)."
       CALL Error(2)
    End If

!!$    If(Integracao)Then
!!$       
!!$       integrando%Esp = pts%Coord(2,1) - pts%Coord(1,1)
!!$       
!!$       integrando%n = (pts%n - 1)
!!$       Allocate(integrando%y(0:integrando%n))
!!$       
!!$       !       DerivANT = pts%Coord(pts%n,2) - pts%Coord(pts%n-1,2)
!!$       
!!$       ! Encontra valor máximo assumido pela função no final da integração
!!$       !       Do i = pts%n-1,1,-1
!!$       !
!!$       !          DerivATUAL = pts%Coord(i,2)-pts%Coord(i-1,2)
!!$       !          
!!$       !          If(DerivATUAL*DerivANT.LT.0)Then
!!$       !             EXIT
!!$       !          Else
!!$       !             DerivANT = DerivATUAL
!!$       !          End If
!!$       !          
!!$       !       End Do
!!$       
!!$       ! OBSOLETO! FORNECE VALORES NÃO CONTÍNUOS SE NÃO IR PARA R GRANDE
!!$       !       ValMAX = ABS(pts%Coord(i,2))
!!$       
!!$       !Derivada no final
!!$       ValMax=(pts%Coord(pts%n,2)-pts%Coord(pts%n-2,2))/(2*integrando%Esp)
!!$       !Coeficiente da função em R grande
!!$       ValMax=sqrt(pts%Coord(pts%n-1,2)**2 + ValMax**2)
!!$       
!!$       Do i = 1,pts%n,1
!!$          integrando%y(i-1) = BesselJ(l,pts%Coord(i,1))*Potential(pts%Coord(i,1)/wv)*pts%Coord(i,2)/ValMAX
!!$       End do
!!$       
!!$       delta_l = ATAN(-QuadIgEsp(integrando,1)/Ener)
!!$       !       delta_l = ATAN(-QuadIgEsp(integrando,5)*wv/Ener)
!!$       
!!$       Deallocate(Integrando%y)
!!$       
!!$    Else
       
       sum_x = 0
       sum_x_sqr = 0
       sum_x_y = 0
       sum_y = 0
       
       Allocate(points_mod(size(points,1),2))

       Do i = 1,size(points,1),1
          points_mod(i,2) = points(i,2)/Sph_Bessel_J(l,points(i,1))
          points_mod(i,1) = Sph_Bessel_N(l,points(i,1))/Sph_Bessel_J(l,points(i,1))
       End Do

       Do i = 1,size(points,1),1
          sum_x = sum_x + points_mod(i,1)
          sum_y = sum_y + points_mod(i,2)
          sum_x_sqr = sum_x_sqr + points_mod(i,1)**2
          sum_x_y = sum_x_y + points_mod(i,1)*points_mod(i,2)
       End Do

       Calc_Phase_Shift = atan(-(size(points,1)*sum_x_y - sum_x*sum_y)/(sum_x_sqr*sum_y - sum_x*sum_x_y))
       
!!    End if
    
     END Function Calc_Phase_Shift
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  Real(kind = real_kind) Function Calc_Scatt_Length(points,met,wave_vector,r_bar,beta)

    IMPLICIT NONE

    Real(kind = real_kind), Dimension(:,:), Intent(in) :: points
    Character(LEN = 20), intent(in) :: met
    Real(kind = real_kind), Intent(in), optional :: wave_vector, r_bar, beta


    Real(kind = real_kind) :: sum_x, sum_x_sqr
    Real(kind = real_kind) :: sum_x_y, sum_y

    integer :: i

    SELECT CASE(met)

    CASE("DIRECT")
       sum_x = zero
       sum_x_sqr = zero
       sum_x_y = zero
       sum_y = zero

       ! Convenção: se wave_vector<0 então o cálculo foi feito a energia 0 e pts%coord(i,1) já está em unidade de distância
       Do i = 1,size(points,1),1
          
          sum_y = sum_y + points(i,2)
          
          if(present(wave_vector))then
             sum_x = sum_x + points(i,1)/wave_vector
             sum_x_sqr = sum_x_sqr + (points(i,1)/wave_vector)**2
             sum_x_y = sum_x_y + (points(i,1)/wave_vector)*points(i,2)
          else
             sum_x = sum_x + points(i,1)
             sum_x_sqr = sum_x_sqr + points(i,1)**2
             sum_x_y = sum_x_y + points(i,1)*points(i,2)
          end if

       End Do
       
       Calc_Scatt_Length = -(sum_x_sqr*sum_y - sum_x*sum_x_y)/(size(points,1)*sum_x_y - sum_x*sum_y)

    CASE("MESHKOV")

       if(.not.present(r_bar).or..not.present(beta))then
          error_msg = "The Meshkov parameters should be passed for MESHKOV method of scattering length calculation."
          CALL Error(2)
       end if

       Calc_Scatt_Length = r_bar*(2*points(1,2)/(pi*beta) + one)

    CASE DEFAULT

       error_msg = "Inconsistency in method of scattering length calculation: "//TRIM(met)
       CALL Error(2)

    END SELECT
    
  END Function Calc_Scatt_Length
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  Real(kind = real_kind) Function Analytical_LenJonN_Nm2_Scatt_Length()

    IMPLICIT NONE

    Real(kind = real_kind) :: tmp1, tmp2, tmp3, tmp4, tmp5
    Real(kind = real_kind) :: x

    x=Pot%radius*sqrt(2*Pot%epsilon*Particles%reduced_mass)

    tmp1 = (Pot%n_lenjon-x-1)/(2*Pot%n_lenjon-4)
    tmp2 = (Pot%n_lenjon-x-3)/(2*Pot%n_lenjon-4)
    tmp3 = (Pot%n_lenjon-3)*one/(Pot%n_lenjon-2)
    tmp4 = (Pot%n_lenjon-1)*one/(Pot%n_lenjon-2)
    tmp5 = Pot%radius*(2*x/(Pot%n_lenjon-2))**(one/(Pot%n_lenjon-2))
  
    Analytical_LenJonN_Nm2_Scatt_Length=tmp5*gamma(tmp1)*gamma(tmp3)/(gamma(tmp2)*gamma(tmp4))

  End Function Analytical_LenJonN_Nm2_Scatt_Length
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  Real(kind = real_kind) Function Analytical_Square_Phase_Shift(l,energy,wave_vector)
    
    IMPLICIT NONE
    
    Integer, intent(in) :: l
    real(kind = real_kind), intent(in) :: energy

    real(kind = real_kind) :: wave_vector
    real(kind = real_kind) :: tmp1, tmp2, tmp3, a1, a2
    real(kind = real_kind) :: derbessJa1, derbessJa2, derbessNa1

    tmp1 = sqrt(1-Pot%epsilon/energy)
    a1 = wave_vector*Pot%radius
    a2 = a1*tmp1

    derbessJa1 = (l*(Sph_Bessel_J(l-1,a1)/a1)-(l+1)*(Sph_Bessel_J(l+1,a1)/a1))/(2*l+1)
    derbessJa2 = (l*(Sph_Bessel_J(l-1,a2)/a2)-(l+1)*(Sph_Bessel_J(l+1,a2)/a2))/(2*l+1)
    derbessNa1 = (l*(Sph_Bessel_N(l-1,a1)/a1)-(l+1)*(Sph_Bessel_N(l+1,a1)/a1))/(2*l+1)
    
    tmp2 = ((tmp1*derbessJa2*(Sph_Bessel_J(l,a1)/a1)/(Sph_Bessel_J(l,a2)/a2))-derbessJa1)
    tmp3 = ((tmp1*derbessJa2*(Sph_Bessel_N(l,a1)/a1)/(Sph_Bessel_J(l,a2)/a2))-derbessNa1)
    
    Analytical_Square_Phase_Shift = atan(tmp2/tmp3)

  End Function Analytical_Square_Phase_Shift
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  Real(kind = real_kind) Function Analytical_HardSphere_Phase_Shift(l,wave_vector)
    
    IMPLICIT NONE
    
    Integer, intent(in) :: l
    real(kind = real_kind), intent(in) :: wave_vector

    real(kind = real_kind) :: tmp

    tmp = Sph_Bessel_J(l,wave_vector*Pot%radius)/Sph_Bessel_N(l,wave_vector*Pot%radius)
    analytical_HardSphere_Phase_Shift = atan(tmp)

  End Function Analytical_HardSphere_Phase_Shift
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
END Module ModCalc
