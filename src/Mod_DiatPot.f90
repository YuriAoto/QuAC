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
! Module for potential definitions and operations
!
! All the potential information, including the masses of particle, are obteined by Define_Potential subroutine
! and are stored in the protected Pot variable. Additional information and conversion to atomic units are carried by
! Check_Potential subroutine. All information in Pot variable are stored in atomic units. The potential at a given 
! internuclar distance, in atomic units, a0, is calculated by Potential function, and retorned in hartree.
!
!
MODULE ModPot

  USE ModUtil, ONLY : real_kind, at_uni_dist, at_uni_ener, at_uni_dip_mom, status, error_msg, ios_var, ios_ok, uni_inp, uni_log, &
       units_len, zero, one, two, pi, eps_NR_default, Unit_Converter, Newton_Raphson, Ordered_Search, Error 
  USE, INTRINSIC :: iso_c_binding, ONLY: c_char, c_null_char

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: Pot, TypePrintFormats, formats, SLL_Size, RLL_Size,stop_error,&
       ILL_Size, Free_SLL, Free_RLL, Free_ILL, Set_uni_EE_input, Total_pot, curr_pot_knot, save_Potential, use_Pot_vector, &
       POTENTIAL_kw, DIP_MOM_kw, Define_Potential,  Select_Potential, Check_Potential, Potential, TypeIntegerLinkedList, &
       TypeRealLinkedList, TypeStringLinkedList, Non_Consistent_Info, Print_Potential_Data, Check_Unit, Free_Potential_Memory, &
       Classical_Turning_Points

  TYPE TypeIntegerLinkedList
     INTEGER :: VALUE
     TYPE(TypeIntegerLinkedList), POINTER :: next=>NULL()
  END TYPE TypeIntegerLinkedList

  TYPE TypeRealLinkedList
     REAL(kind = real_kind) :: VALUE
     TYPE(TypeRealLinkedList), POINTER :: next=>NULL()
  END TYPE TypeRealLinkedList

  TYPE TypeStringLinkedList
     CHARACTER(LEN=100) :: VALUE
     TYPE(TypeStringLinkedList), POINTER :: next=>NULL()
  END TYPE TypeStringLinkedList

  TYPE TypePot_DipMomCurve
     CHARACTER(LEN = 100) :: name = ""
     CHARACTER(LEN = 30) :: TYPE = ""
     LOGICAL :: is_dip_mom = .FALSE.

     ! Variables for each potential
     INTEGER :: rep_exp=0, n_lenjon=0, m_lenjon=0
     REAL(kind = real_kind) :: radius=zero, Beta=zero, epsilon=zero
     REAL(kind = real_kind) :: coef_6=zero, coef_8=zero, coef_10=zero
     REAL(kind = real_kind) :: r_infty=zero, ener_infty=zero
     REAL(kind = real_kind) :: TolLimInf = 1.0E-10_real_kind
     REAL(kind = real_kind), ALLOCATABLE, DIMENSION(:) :: coef
     REAL(kind = real_kind), ALLOCATABLE, DIMENSION(:) :: r
     REAL(kind = real_kind), ALLOCATABLE, DIMENSION(:) :: V
     REAL(kind = real_kind), ALLOCATABLE, DIMENSION(:) :: z

     INTEGER :: lambda=0

     ! Variables for direct calculation
     CHARACTER(LEN=100) :: ES_command="", rm_command=""
     CHARACTER(LEN=20) :: format_ener="", format_r="", keyword=""
     TYPE(TypeStringLinkedList), POINTER :: input=>NULL()
     CHARACTER(LEN=100,kind=c_char) :: get_ener_command_c=""
     INTEGER :: uni_EE_input=0

     ! Variables for general properties
     REAL(kind = real_kind) :: zero_energy_ret_p=zero
     REAL(kind = real_kind) :: lower_bound=zero
     REAL(kind = real_kind) :: upper_bound=zero
     REAL(kind = real_kind) :: eq_dist=zero
     REAL(kind = real_kind) :: abs_min_energy=zero
     REAL(kind = real_kind) :: diss_ener=zero
     REAL(kind = real_kind) :: force_const=zero

     ! Units used in potential definition
     CHARACTER(LEN = units_len) :: energy_uni = "hartree"
     CHARACTER(LEN = units_len) :: dist_uni = "a0"

     TYPE(TypePot_DipMomCurve), POINTER :: next=>NULL()
  END TYPE TypePot_DipMomCurve

  TYPE TypePrintFormats
     CHARACTER(LEN = 20) :: energy,               cross_section,               angle,               ang_momentum
     LOGICAL ::             energy_def = .false., cross_section_def = .false., angle_def = .false., ang_momentum_def = .false.

     CHARACTER(LEN = 20) :: distance,               wave_function,               potential,               energy_cm1
     LOGICAL ::             distance_def = .false., wave_function_def = .false., potential_def = .false., energy_cm1_def = .false.

     CHARACTER(LEN = 20) :: einstein_coef,               dip_mom,               phase_shift,               dunhan
     LOGICAL ::             einstein_coef_def = .false., dip_mom_def = .false., phase_shift_def = .false., dunhan_def = .false.

     CHARACTER(LEN = 20) :: vib_level
     LOGICAL ::             vib_level_def = .false.

     LOGICAL :: default_values = .false.
  END TYPE TypePrintFormats

  REAL(kind = real_kind) :: Energy_TMP_module

  INTEGER :: uni_EE_input_beg

  TYPE(TypePrintFormats), POINTER :: formats=>NULL()

  TYPE(TypePot_DipMomCurve), PROTECTED, SAVE, POINTER :: Pot_list=>NULL()

  TYPE(TypePot_DipMomCurve), POINTER :: Pot=>NULL()

  CHARACTER(len=40), PARAMETER :: free_particle_name = "FREE_PARTICLE"

  CHARACTER(len = 30), PARAMETER :: POTENTIAL_kw="POTENTIAL"
  CHARACTER(len = 30), PARAMETER :: DIP_MOM_kw="DIP_MOM"

  LOGICAL :: stop_error = .TRUE.

  Type(TypeRealLinkedList), POINTER :: Total_pot
  Type(TypeRealLinkedList), POINTER :: curr_pot_knot
  LOGICAL :: save_Potential = .FALSE.
  LOGICAL :: use_Pot_vector = .FALSE.

CONTAINS
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  INTEGER FUNCTION SLL_Size(list)
    
    IMPLICIT NONE

    TYPE(TypeStringLinkedList), INTENT(in), POINTER :: list
    TYPE(TypeStringLinkedList), POINTER :: current_string

    SLL_Size = 0
    current_string=>list
    DO WHILE(ASSOCIATED(current_string))
       current_string => current_string%next
       SLL_Size = SLL_Size + 1
    END DO

  END FUNCTION SLL_Size
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  INTEGER FUNCTION RLL_Size(list)
    
    IMPLICIT NONE

    TYPE(TypeRealLinkedList), INTENT(in), POINTER :: list
    TYPE(TypeRealLinkedList), POINTER :: current_real

    RLL_Size = 0
    current_real=>list
    DO WHILE(ASSOCIATED(current_real))
       current_real => current_real%next
       RLL_Size = RLL_Size + 1
    END DO

  END FUNCTION RLL_Size
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  INTEGER FUNCTION ILL_Size(list)
    
    IMPLICIT NONE

    TYPE(TypeIntegerLinkedList), INTENT(in), POINTER :: list
    TYPE(TypeIntegerLinkedList), POINTER :: current_int

    ILL_Size = 0
    current_int=>list
    DO WHILE(ASSOCIATED(current_int))
       current_int => current_int%next
       ILL_Size = ILL_Size + 1
    END DO

  END FUNCTION ILL_Size
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  SUBROUTINE Free_SLL(list)

    IMPLICIT NONE

    TYPE(TypeStringLinkedList), INTENT(inout), POINTER :: list
    TYPE(TypeStringLinkedList), POINTER :: next=>NULL()

    DO WHILE(ASSOCIATED(list))
       next => list%next
       DEALLOCATE(list)
       list => next
    END DO

  END SUBROUTINE Free_SLL
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  SUBROUTINE Free_RLL(list)

    IMPLICIT NONE

    TYPE(TypeRealLinkedList), INTENT(inout), POINTER :: list
    TYPE(TypeRealLinkedList), POINTER :: next=>NULL()

    DO WHILE(ASSOCIATED(list))
       next => list%next
       DEALLOCATE(list)
       list => next
    END DO

  END SUBROUTINE Free_RLL
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  SUBROUTINE Free_ILL(list)

    IMPLICIT NONE

    TYPE(TypeIntegerLinkedList), INTENT(inout), POINTER :: list
    TYPE(TypeIntegerLinkedList), POINTER :: next=>NULL()

    DO WHILE(ASSOCIATED(list))
       next => list%next
       DEALLOCATE(list)
       list => next
    END DO

  END SUBROUTINE Free_ILL
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  SUBROUTINE Set_uni_EE_input(unit)

    IMPLICIT NONE

    INTEGER :: unit

    uni_EE_input_beg = unit

  END SUBROUTINE Set_uni_EE_input
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  SUBROUTINE Define_Potential()

    IMPLICIT NONE

    CHARACTER(LEN = 130) :: info, info2

    CHARACTER(LEN=3) :: lin_log_print
    REAL(kind = real_kind) :: r_min_print, r_max_print, delta_print
    CHARACTER(LEN=15) :: dist_uni_print, ener_uni_print
    REAL(kind = real_kind) :: r_print, pot_print
    LOGICAL :: print_pot

    TYPE(TypePot_DipMomCurve), POINTER :: current_pot=>NULL()
    TYPE(TypeRealLinkedList), POINTER :: current_real=>NULL(), r_tmp=>NULL()
    TYPE(TypeRealLinkedList), POINTER :: current_real2=>NULL(), V_tmp=>NULL()
    TYPE(TypeRealLinkedList), POINTER :: current_real3=>NULL(), coef_tmp=>NULL()
    TYPE(TypeStringLinkedList), POINTER :: current_string=>NULL()

    CHARACTER(LEN = 2) :: Pot_or_Mom_dip
    CHARACTER(LEN = 100) :: name_candidate
    CHARACTER(LEN = 100) :: input_line
    INTEGER :: aux_int, aux_int2
    REAL(kind = real_kind) :: aux_real, aux_real2

    BACKSPACE(uni_inp)

    print_pot = .FALSE.

    ! CHANGE!! GET PRINT INFORMATION LATER
    READ(uni_inp,*,Iostat=ios_var) info
    if((TRIM(info).NE.TRIM(POTENTIAL_kw).AND.&
         TRIM(info).NE.TRIM(DIP_MOM_kw)).OR.&
         ios_var.ne.ios_ok)THEN
       Pot=>NULL()
       RETURN
    END IF

    IF(.NOT.ASSOCIATED(Pot))THEN
       ALLOCATE(Pot_list)
       Pot=>Pot_list
    ELSE
       ALLOCATE(Pot%next)
       Pot => Pot%next
    END IF

    BACKSPACE(uni_inp)
    READ(uni_inp,FMT='(A130)',Iostat=ios_var) info

    IF(INDEX(ADJUSTL(info),TRIM(DIP_MOM_kw)).EQ.1) Pot%is_dip_mom = .TRUE.

    uni_EE_input_beg = uni_EE_input_beg+1
    Pot%uni_EE_input = uni_EE_input_beg

    info = adjustl(info)

    if(Pot%is_dip_mom)then
       READ(info(LEN_TRIM(DIP_MOM_kw)+2:),FMT='(A100)',iostat=ios_var) name_candidate
    else
       READ(info(LEN_TRIM(POTENTIAL_kw)+2:),FMT='(A100)',iostat=ios_var) name_candidate
    end if

    IF(name_candidate.EQ.free_particle_name.AND..NOT.Pot%is_dip_mom)THEN
       error_msg = TRIM(name_candidate)//" is not an allowed potential name."
       CALL Error(1)
    END IF

    current_pot=>Pot_list

    DO WHILE(ASSOCIATED(current_pot))
       IF(TRIM(name_candidate).EQ.TRIM(current_pot%name))THEN
          error_msg = "The name "//TRIM(name_candidate)//" was already used."
          CALL Error(1)
       END IF
       current_pot=>current_pot%next
    END DO

    Pot%name=name_candidate

    READ(uni_inp,*,Iostat=ios_var) info2
    Backspace(uni_inp)
    ! Get print information
    IF(info2.EQ."Print")THEN
       
       READ(uni_inp,*,Iostat=ios_var) info2,lin_log_print,r_min_print,r_max_print,delta_print,dist_uni_print,ener_uni_print
       
       IF(ios_var.NE.0)THEN
          IF(Pot%is_dip_mom)THEN
             error_msg = "Error reading print information for dipole moment, in "//TRIM(Pot%name)
          ELSE
             error_msg = "Error reading print information for potential, in "//TRIM(Pot%name)
          END IF
          CALL Error(0)
       END IF
       
       IF(lin_log_print.NE."Lin".AND.lin_log_print.NE."Log") THEN
          error_msg = "Pass 'Lin' or 'Log' after 'Print'."
          CALL Error(0)
       END IF
       
       IF(lin_log_print.EQ."Lin".AND.(r_min_print.LT.0.OR.r_max_print.LT.0))THEN
          IF(Pot%is_dip_mom)THEN
             error_msg = "Found negative value for distance in print information of dipole moment."
          ELSE
             error_msg = "Found negative value for distance in print information of potential."
          END IF
          CALL Error(0)
       END IF
       
       IF((r_max_print-r_min_print)*delta_print.LT.0)THEN
          IF(Pot%is_dip_mom)THEN
             error_msg = "The variation of the distance should be consistent with the initial and final values, "//&
                  &"in print information of dipole moment."
          ELSE
             error_msg = "The variation of the distance should be consistent with the initial and final values, "//&
                  &"in print information of potential."             
          END IF
          CALL Error(0)
       END IF
       
       aux_real = Unit_Converter(dist_uni_print,at_uni_dist)
       CALL Check_Unit(dist_uni_print,"distance")
       
       IF(Pot%is_dip_mom)THEN
          aux_real = Unit_Converter(ener_uni_print,at_uni_dip_mom)
          CALL Check_Unit(ener_uni_print,"dipole moment")
       ELSE
          aux_real = Unit_Converter(ener_uni_print,at_uni_ener)
          CALL Check_Unit(ener_uni_print,"energy")
       END IF
       
       print_pot = .TRUE.
       
    END IF

    ! Read name and kind of potential
    READ(uni_inp,*,Iostat=ios_var) Pot%type

    ! Select expoents of a Lennard Jones type of potential
    IF(INDEX(Pot%type,"LenJon").EQ.1.AND.INDEX(Pot%type,"-").EQ.7.OR.&
         INDEX(Pot%type,"LSLenJon").EQ.1.AND.INDEX(Pot%type,"-").EQ.9&
         )THEN

       IF(INDEX(Pot%type,"LenJon").EQ.1)then
          READ(unit=Pot%type,FMT='(7X,I1,1X,I1)',Iostat=ios_var) Pot%n_lenjon,Pot%m_lenjon
       else
          READ(unit=Pot%type,FMT='(9X,I1,1X,I1)',Iostat=ios_var) Pot%n_lenjon,Pot%m_lenjon
       end IF

       IF(ios_var.NE.0)THEN
          IF(INDEX(Pot%type,"LenJon").EQ.1)then
             error_msg = "Error reading Lennard-Jones type of curve. Pass like ""LenJon-N-M"", with M and N just one digit."
             CALL Error(0)
          else
             error_msg = "Error reading Lennard-Jones type of curve. Pass like ""LSLenJon-N-M"", with M and N just one digit."
             CALL Error(0)
          end IF
       END IF

       CALL CheckGreaterThanZero("m-LenJon",int=Pot%m_lenjon)
       CALL CheckGreaterThanZero("n-LenJon",int=Pot%n_lenjon)

       IF(ios_var.NE.0)THEN
          IF(INDEX(Pot%type,"LenJon").EQ.1)then
             Pot%type = "LenJonNM" 
          else
             Pot%type = "LSLenJonNM"
          end IF
       END IF

    END IF

    ! Get units and check type of potential
    BACKSPACE(uni_inp)
    SELECT CASE(Pot%type)

    CASE("HardSphere")
       IF(Pot%is_dip_mom)THEN
          READ(uni_inp,*,Iostat=ios_var) info2,Pot%dist_uni
       ELSE
          READ(uni_inp,*,Iostat=ios_var) info2,Pot%lambda,Pot%dist_uni
       END IF
    CASE("SoftSphere","Square","LenJon","LenJonNM","Morse","MurSorbie","LinearSpl","NatCubSpl","CubSplr6",&
         &"CubSplr6r8","CubSplFitr6","CubSplFitr6r8","CubSplFitr6r8r10","LSLenJon","LSLenJonNM","AppMorse","Direct")
       IF(Pot%is_dip_mom)THEN
          READ(uni_inp,*,Iostat=ios_var) info2,Pot%dist_uni,Pot%energy_uni
       ELSE
          READ(uni_inp,*,Iostat=ios_var) info2,Pot%lambda,Pot%dist_uni,Pot%energy_uni
       END IF
    CASE DEFAULT
       error_msg = TRIM(Pot%type)//" is an unknow type of curve."
       CALL Error(0)
    END SELECT
    IF(ios_var.NE.0)THEN
       error_msg = "Error reading potential definition."
       CALL Error(0)
    END IF

    aux_real = Unit_Converter(Pot%dist_uni,at_uni_dist)
    CALL Check_Unit(Pot%dist_uni,"distance")

    IF(Pot%is_dip_mom)THEN
       aux_real = Unit_Converter(Pot%energy_uni,at_uni_dip_mom)
       CALL Check_Unit(Pot%energy_uni,"dipole moment")
    ELSE
       aux_real = Unit_Converter(Pot%energy_uni,at_uni_ener)
       CALL Check_Unit(Pot%energy_uni,"energy")
    END IF

    ! Get the parameters
    SELECT CASE(Pot%type)

    CASE("HardSphere")

       READ(uni_inp,*,Iostat=ios_var) Pot%radius
       CALL CheckIos()
       CALL CheckGreaterThanZero("radius",float=Pot%radius)

       READ(uni_inp,*,Iostat=ios_var) ! Read anything to set correctly the position

    CASE("SoftSphere")

       READ(uni_inp,*,Iostat=ios_var) Pot%radius, Pot%epsilon, Pot%rep_exp
       CALL CheckIos()
       CALL CheckGreaterThanZero("radius",float=Pot%radius)
       CALL CheckGreaterThanZero("epsilon",float=Pot%epsilon)
       CALL CheckGreaterThanZero("rep_exp",int=Pot%rep_exp)

       READ(uni_inp,*,Iostat=ios_var) ! Read anything to set correctly the position

    CASE("Square")

       READ(uni_inp,*,Iostat=ios_var) Pot%radius,Pot%epsilon
       CALL CheckIos()
       CALL CheckGreaterThanZero("radius",float=Pot%radius)

       READ(uni_inp,*,Iostat=ios_var) ! Read anything to set correctly the position

    CASE("LenJon","LenJonNM")

       READ(uni_inp,*,Iostat=ios_var) Pot%radius,Pot%epsilon
       CALL CheckIos()
       CALL CheckGreaterThanZero("radius",float=Pot%radius)
       CALL CheckGreaterThanZero("dissociation energy",float=Pot%epsilon)

       READ(uni_inp,*,Iostat=ios_var) ! Read anything to set correctly the position

    CASE("Morse")

       READ(uni_inp,*,Iostat=ios_var) Pot%radius,Pot%epsilon,Pot%Beta
       CALL CheckIos()
       CALL CheckGreaterThanZero("radius",float=Pot%radius)
       CALL CheckGreaterThanZero("dissociation energy",float=Pot%epsilon)
       CALL CheckGreaterThanZero("beta",float=Pot%beta)

       READ(uni_inp,*,Iostat=ios_var) ! Read anything to set correctly the position

    CASE("MurSorbie")

       READ(uni_inp,*,Iostat=ios_var) Pot%radius,Pot%epsilon
       CALL CheckIos()
       CALL CheckGreaterThanZero("radius",float=Pot%radius)
       CALL CheckGreaterThanZero("dissociation energy",float=Pot%epsilon)

       current_real => NULL()
       aux_int=0

       DO WHILE(.TRUE.)

          READ(uni_inp,*,Iostat=ios_var) aux_real
          IF(ios_var.NE.0) EXIT

          IF(.NOT.ASSOCIATED(current_real))THEN
             ALLOCATE(coef_tmp)
             current_real=>coef_tmp
          ELSE
             ALLOCATE(current_real%next)
             current_real=>current_real%next
          END IF

          aux_int = aux_int + 1
          current_real%value = aux_real

       END DO

       ALLOCATE(Pot%coef(aux_int))

       current_real=>coef_tmp
       DO aux_int2 = 1,aux_int,1
          Pot%coef(aux_int2) = current_real%value

          current_real2=>current_real%next
          DEALLOCATE(current_real)
          current_real=>current_real2
       END DO

    CASE ("LinearSpl","NatCubSpl","CubSplr6","CubSplr6r8","CubSplFitr6","CubSplFitr6r8","CubSplFitr6r8r10","LSLenJon",&
         "LSLenJonNM","AppMorse")

       current_real => NULL()
       current_real2 => NULL()
       aux_int=0

       IF(Pot%type.EQ."CubSplFitr6")THEN
          READ(uni_inp,*,Iostat=ios_var) Pot%coef_6
          CALL CheckIos()
       END IF
       
       IF(Pot%type.EQ."CubSplFitr6r8")THEN
          READ(uni_inp,*,Iostat=ios_var) Pot%coef_6, Pot%coef_8
          CALL CheckIos()
       END IF

       IF(Pot%type.EQ."CubSplFitr6r8r10")THEN
          READ(uni_inp,*,Iostat=ios_var) Pot%coef_6, Pot%coef_8, Pot%coef_10
          CALL CheckIos()
       END IF

       DO
          READ(uni_inp,*,Iostat=ios_var) aux_real, aux_real2
          IF(ios_var.NE.0) EXIT

          CALL CheckGreaterThanZero("distance",float=aux_real)
          IF(ASSOCIATED(current_real).AND.aux_real.LE.current_real%value)THEN
             error_msg = "Not sorted curve. It should be passed increasing in distance."
             CALL Error(0)
          END IF

          IF(.NOT.ASSOCIATED(current_real))THEN
             ALLOCATE(r_tmp)
             current_real=>r_tmp
             ALLOCATE(V_tmp)
             current_real2=>V_tmp
          ELSE
             ALLOCATE(current_real%next)
             current_real=>current_real%next
             ALLOCATE(current_real2%next)
             current_real2=>current_real2%next
          END IF

          aux_int = aux_int + 1
          current_real%value = aux_real
          current_real2%value = aux_real2

       END DO

       IF(aux_int.LT.2)THEN
          error_msg = "Less than two points found in curve definition."
          CALL Error(0)
       END IF

       IF(Pot%type.EQ."CubSplr6r8")THEN ! Last point just to take the limit energy
          aux_int = aux_int-1
       END IF

       ALLOCATE(Pot%r(aux_int),Pot%V(aux_int))

       current_real=>r_tmp
       current_real2=>V_tmp
       DO aux_int2 = 1,aux_int,1
          Pot%r(aux_int2) = current_real%value
          Pot%V(aux_int2) = current_real2%value

          current_real3=>current_real%next
          DEALLOCATE(current_real)
          current_real=>current_real3

          current_real3=>current_real2%next
          DEALLOCATE(current_real2)
          current_real2=>current_real3

       END DO

       IF(Pot%type.EQ."CubSplr6r8")THEN ! Last point just to take the limit energy
          Pot%ener_infty = current_real2%value
          DEALLOCATE(current_real)
          DEALLOCATE(current_real2)
       END IF

    CASE ("Direct")

       READ(uni_inp,FMT='(100A)',Iostat=ios_var) Pot%ES_command
       CALL CheckIos()
       READ(uni_inp,FMT='(100A)',Iostat=ios_var) Pot%get_ener_command_c
       CALL CheckIos()
       READ(uni_inp,FMT='(100A)',Iostat=ios_var) Pot%rm_command
       CALL CheckIos()
       READ(uni_inp,*,Iostat=ios_var) Pot%keyword,Pot%format_ener,Pot%format_r
       CALL CheckIos()
       READ(uni_inp,*,Iostat=ios_var) Pot%r_infty,Pot%zero_energy_ret_p,Pot%lower_bound,Pot%upper_bound,Pot%radius,Pot%force_const
       CALL CheckIos()

       Pot%format_r = '(A,'//TRIM(Pot%format_r)//')'
       Pot%format_ener = '('//TRIM(Pot%format_ener)//')'

       current_string => NULL()

       DO
          READ(uni_inp,FMT='(100A)',iostat=ios_var) input_line
          IF(ios_var.NE.0.OR.INDEX(ADJUSTL(input_line),"POTENTIAL").EQ.1) EXIT

          IF(.NOT.ASSOCIATED(current_string))THEN
             ALLOCATE(Pot%input)
             current_string=>Pot%input
          ELSE
             ALLOCATE(current_string%next)
             current_string=>current_string%next
          END IF

          current_string%value = input_line

       END DO
       Pot%get_ener_command_c=TRIM(Pot%get_ener_command_c)//c_null_char

       Pot%ener_infty = zero
       Pot%eq_dist = Pot%radius
       
       IF(.NOT.Pot%is_dip_mom) then
          Pot%ener_infty = Potential(Pot%r_infty)

          if(Pot%eq_dist.GT.zero)then
             Pot%diss_ener = - Potential(Pot%eq_dist)
          else
             Pot%radius = -one
             Pot%eq_dist = -one
             Pot%diss_ener = -one
          end if

       end IF

    CASE DEFAULT

       IF(Pot%is_dip_mom)THEN
          CALL Non_Consistent_Info("type of dipole moment","dipole moment definition")
       ELSE
          CALL Non_Consistent_Info("type of potential","potential definition")
       END IF

    END SELECT

    CALL Check_Potential ()

    ! Print potential information in .log file
    WRITE(uni_log,FMT='("---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/---")')
    IF(Pot%is_dip_mom)THEN
       WRITE(uni_log,FMT='("Dipole moment ",A,", of type ",A,", in atomic units:")') TRIM(Pot%name),TRIM(Pot%type)
       Pot_or_Mom_dip = "M"
    ELSE
       WRITE(uni_log,FMT='("Potential ",A,", of type ",A,", in atomic units:")') TRIM(Pot%name),TRIM(Pot%type)
       Pot_or_Mom_dip = "V"
    END IF

    SELECT CASE(Pot%type)

    CASE("HardSphere")

       WRITE(uni_log,FMT='(A,"(r) = infinity, if r < ",'//TRIM(formats%distance)//',".")') TRIM(Pot_or_Mom_dip),Pot%radius
       WRITE(uni_log,FMT='(A,"(r) = 0       , se r > ",'//TRIM(formats%distance)//',".")') TRIM(Pot_or_Mom_dip),Pot%radius

    CASE("SoftSphere")

       WRITE(uni_log,FMT='(A,"(r) = ",'//TRIM(formats%potential)//',"*((",'//TRIM(formats%distance)//&
            ',"/r)^",I2," -1), if r < ",'//TRIM(formats%distance)//',".")') TRIM(Pot_or_Mom_dip),Pot%epsilon,Pot%radius,&
            Pot%rep_exp,Pot%radius
       WRITE(uni_log,FMT='(A,"(r) = 0                                  , if r > ",'//TRIM(formats%distance)//',".")') &
            TRIM(Pot_or_Mom_dip),Pot%radius

    CASE("Square")

       WRITE(uni_log,FMT='(A,"(r) = ",'//TRIM(formats%potential)//',", if r < ",'//TRIM(formats%distance)//',".")') &
            &TRIM(Pot_or_Mom_dip),Pot%epsilon,Pot%radius
       WRITE(uni_log,FMT='(A,"(r) = 0          , if r > ",'//TRIM(formats%distance)//',".")') TRIM(Pot_or_Mom_dip),Pot%radius

    CASE("MurSorbie")

       WRITE(uni_log,FMT='(A,"(r) = -",'//TRIM(formats%potential)//',"*exp(-",'//TRIM(formats%distance)//',"*ro)*")')&
            TRIM(Pot_or_Mom_dip),Pot%epsilon,Pot%coef(1)
       WRITE(uni_log,FMT='("(")',Advance='NO')
       DO aux_int = 1,SIZE(Pot%coef)-1,1
          WRITE(uni_log,FMT='('//TRIM(formats%potential)//',"*ro^",I3,"+")',Advance='NO') Pot%coef(aux_int),aux_int

          IF(MOD(aux_int,3).EQ.0)THEN
             WRITE(uni_log,*)
             WRITE(uni_log,FMT='(1X)',Advance='NO')
          ENDIF

       END DO
       WRITE(uni_log,FMT='('//TRIM(formats%potential)//',"*ro^",I3,"+")') Pot%Coef(SIZE(Pot%coef)),SIZE(Pot%coef)
       WRITE(uni_log,FMT='(", ro = r - ",'//TRIM(formats%distance)//')') Pot%radius

    CASE("LenJon","LSLenJon")

       WRITE(uni_log,FMT='(A,"(r) = ",'//TRIM(formats%potential)//',"*((",'//TRIM(formats%distance)//&
            ',"/r)^12 - 2*(",'//TRIM(formats%distance)//',"/r)**6)")') TRIM(Pot_or_Mom_dip),Pot%epsilon,Pot%radius,Pot%radius

       IF(Pot%type.EQ."LSLenJon")THEN

          WRITE(uni_log,FMT='("Whose spectroscopic parameters were obtained by least square"'//&
               '" over the following points:")')
          
          WRITE(uni_log,FMT='("     R               ",A,"            ")') TRIM(Pot_or_Mom_dip)
          DO aux_int = 1,SIZE(Pot%r),1
             WRITE(uni_log,FMT='('//TRIM(formats%distance)//'," ",'//TRIM(formats%potential)//')') Pot%R(aux_int),Pot%V(aux_int)
          END DO
       END IF


    CASE("LenJonNM","LSLenJonNM")

       WRITE(uni_log,FMT='(A,"(r) = ",'//TRIM(formats%potential)//',"*((",'//TRIM(formats%distance)//&
            ',"/r)^",I2," - (",'//TRIM(formats%distance)//',"/r)^",I1,")")') &
            TRIM(Pot_or_Mom_dip),Pot%epsilon,Pot%radius,Pot%n_lenjon+Pot%m_lenjon,Pot%radius,Pot%n_lenjon

       IF(Pot%type.EQ."LSLenJonNM")THEN
          
          WRITE(uni_log,FMT='("Whose spectroscopic parameters were obtained by least square"'//&
               '" over the following points:")')
          
          WRITE(uni_log,FMT='("     R               ",A,"            ")') TRIM(Pot_or_Mom_dip)
          DO aux_int = 1,SIZE(Pot%r),1
             WRITE(uni_log,FMT='('//TRIM(formats%distance)//'," ",'//TRIM(formats%potential)//')') Pot%R(aux_int),Pot%V(aux_int)
          END DO
       END IF


    CASE("Morse","AppMorse")

       WRITE(uni_log,FMT='(A,"(r) = ",'//TRIM(formats%potential)//',"*((1 - exp(-",'//TRIM(formats%potential)//&
            ',"*(r-",'//TRIM(formats%distance)//',")))^2 - 1)")') TRIM(Pot_or_Mom_dip),Pot%epsilon,Pot%Beta,Pot%radius

       IF(Pot%type.EQ."AppMorse")THEN
          WRITE(uni_log,FMT='("Whose spectroscopic parameters were obtained by the Natural Cubic Splines"'//&
               '" over the following points:")')
          WRITE(uni_log,FMT='("     R               ",A,"            ")') TRIM(Pot_or_Mom_dip)
          DO aux_int = 1,SIZE(Pot%r),1
             WRITE(uni_log,FMT='('//TRIM(formats%distance)//'," ",'//TRIM(formats%potential)//')') Pot%R(aux_int),Pot%V(aux_int)
          END DO
       END IF

    CASE("LinearSpl","NatCubSpl","CubSplr6","CubSplr6r8","CubSplFitr6","CubSplFitr6r8","CubSplFitr6r8r10")

       SELECT CASE(Pot%type)

       CASE("LinearSpl")
          WRITE(uni_log,FMT='("Linear spline over the following points:")')
       CASE("NatCubSpl","CubSplr6","CubSplr6r8","CubSplFitr6","CubSplFitr6r8","CubSplFitr6r8r10")
          WRITE(uni_log,FMT='("Natural cubic spline over the following points:")')
       CASE DEFAULT
          IF(Pot%is_dip_mom)THEN
             CALL Non_Consistent_Info("type of potential","writing potential in .log file")
          ELSE
             CALL Non_Consistent_Info("type of dipole moment","writing potential in .log file")
          END IF
       END SELECT

       WRITE(uni_log,FMT='("     R               ",A,"            ")') TRIM(Pot_or_Mom_dip)
       DO aux_int = 1,SIZE(Pot%r),1
          WRITE(uni_log,FMT='('//TRIM(formats%distance)//'," ",'//TRIM(formats%potential)//')') Pot%R(aux_int),Pot%V(aux_int)
       END DO
       WRITE(uni_log,FMT='("---")')

       IF(Pot%type.EQ."CubSplr6".OR.Pot%type.EQ."CubSplFitr6")THEN
          WRITE(uni_log,FMT='("joined with ",'//TRIM(formats%potential)//',"/r⁶ tail")') Pot%coef_6
       END IF
       IF(Pot%type.EQ."CubSplr6r8".OR.Pot%type.EQ."CubSplFitr6r8")THEN
          WRITE(uni_log,FMT='("joined with ",'//TRIM(formats%potential)//',"/r⁶ + ",'//TRIM(formats%potential)//',"/r⁸ tail")')&
               Pot%coef_6, Pot%coef_8
       END IF
       IF(Pot%type.EQ."CubSplFitr6r8r10")THEN
          WRITE(uni_log,FMT='("joined with ",'//TRIM(formats%potential)//',"/r⁶ + ",'//TRIM(formats%potential)//',"/r⁸  + ",'&
               //TRIM(formats%potential)//',"/r¹⁰ tail")')&
               Pot%coef_6, Pot%coef_8, Pot%coef_10
       END IF

    CASE("Direct")

       IF(Pot%is_dip_mom)THEN
          WRITE(uni_log,FMT='("Direct potential. Using the following input:")')
       ELSE
          WRITE(uni_log,FMT='("Direct dipole moment. Using the following input:")')
       END IF

       current_string=>Pot%input
       DO WHILE(ASSOCIATED(current_string))
          WRITE(uni_log,FMT='(A)') TRIM(current_string%value)
          current_string=>current_string%next
       END DO

       WRITE(uni_log,FMT='("Where the keyword """,A,""" is replaced by the distance.")') TRIM(Pot%keyword)
       WRITE(uni_log,FMT='("The following command is used to run the electronic energy calculation:")')
       WRITE(uni_log,FMT='(A)') TRIM(Pot%ES_command)
       WRITE(uni_log,FMT='(", the following to get the energy:")')
       WRITE(uni_log,FMT='(A)') TRIM(Pot%get_ener_command_c)
       WRITE(uni_log,FMT='(" and cleaning for the next run with")')
       WRITE(uni_log,FMT='(A)') TRIM(Pot%rm_command)
       WRITE(uni_log,FMT='("Formats used for write internuclear distance: ",A)') TRIM(Pot%format_r)
       WRITE(uni_log,FMT='("Formats used for get energy: ",A)') TRIM(Pot%format_ener)

    CASE DEFAULT

       IF(Pot%is_dip_mom)THEN
          CALL Non_Consistent_Info("type of dipole moment","writing the potential in .log file")
       ELSE
          CALL Non_Consistent_Info("type of potential","writing the potential in .log file")
       END IF

    END SELECT

    IF(.NOT.Pot%is_dip_mom)THEN
       WRITE(uni_log,'("General properties:")') 
       WRITE(uni_log,'("Zero energy classical turning point: ",'//TRIM(formats%distance)//',1X,A)') &
            Pot%zero_energy_ret_p, at_uni_dist
       WRITE(uni_log,'("Absolute minimum energy: ",'//TRIM(formats%energy)//',1X,A)') Pot%abs_min_energy, at_uni_ener
       WRITE(uni_log,'("Lower bound: ",'//TRIM(formats%distance)//',1X,A)') Pot%lower_bound, at_uni_dist
       if(Pot%upper_bound.GT.zero)then
          WRITE(uni_log,'("Upper bound: ",'//TRIM(formats%distance)//',1X,A)') Pot%upper_bound, at_uni_dist
       else
          WRITE(uni_log,'("Potential has an undefined upperb bound")') 
       end if
       if(Pot%eq_dist.GT.zero)then
          WRITE(uni_log,'("Equilibrium distance: ",'//TRIM(formats%distance)//',1X,A)') Pot%eq_dist, at_uni_dist
          WRITE(uni_log,'("Dissociation energy: ",'//TRIM(formats%energy)//',1X,A)') Pot%diss_ener, at_uni_ener
          WRITE(uni_log,'("Force constant: ",'//TRIM(formats%energy)//',1X,A)') Pot%force_const, "u.a."
       else
          WRITE(uni_log,'("Potential without minimum")') 
       end if
    END IF

    WRITE(uni_log,*) 

    IF(print_pot) THEN
       IF(Pot%is_dip_mom)THEN
          WRITE(uni_log,FMT='("Printing dipole moment:")')
       ELSE
          WRITE(uni_log,FMT='("Printing potential:")')
       END IF

       IF(Pot%is_dip_mom)THEN
          WRITE(uni_log,FMT='("Distance/",A,"     Dipole Moment/",A)') TRIM(dist_uni_print), TRIM(ener_uni_print)
       ELSE
          WRITE(uni_log,FMT='("Distance/",A,"     Energy/",A)') TRIM(dist_uni_print), TRIM(ener_uni_print)
       END IF

       DO aux_int = 1,NINT((r_max_print-r_min_print)/delta_print)+1,1

          r_print = r_min_print + (aux_int-1)*delta_print

          IF(lin_log_print.EQ."Log") r_print = 10**r_print

          WRITE(uni_log,FMT='('//TRIM(formats%distance)//'," ")',ADVANCE='NO') r_print

          IF(Pot%is_dip_mom)THEN
             pot_print = Unit_Converter(at_uni_dip_mom,ener_uni_print,Potential(Unit_Converter(dist_uni_print,at_uni_dist,r_print)))
          ELSE
             pot_print = Unit_Converter(at_uni_ener,ener_uni_print,Potential(Unit_Converter(dist_uni_print,at_uni_dist,r_print)))
          END IF

          WRITE(uni_log,FMT='('//TRIM(formats%potential)//')') pot_print

       END DO

       WRITE(uni_log,FMT='("---")')

    END IF

  CONTAINS
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
    SUBROUTINE CheckIos()

      IMPLICIT NONE

      IF(ios_var.NE.0)THEN
         error_msg = "Error reading potential information."
         CALL Error(0)
      END IF

    END SUBROUTINE CheckIos
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
    SUBROUTINE CheckGreaterThanZero(information,int,float)

      IMPLICIT NONE

      CHARACTER(len = *), INTENT(in) :: information
      INTEGER, INTENT(in), OPTIONAL :: int
      REAL(kind = real_kind), INTENT(in), OPTIONAL :: float

      IF((PRESENT(int).AND.int.LE.0).OR.(PRESENT(float).AND.float.LE.0))THEN
         error_msg = "The value of "//TRIM(information)//" should be positive."
         CALL Error(0)
      END IF

    END SUBROUTINE CheckGreaterThanZero
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  END SUBROUTINE Define_Potential
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  SUBROUTINE Select_Potential(potential_name)

    IMPLICIT NONE

    CHARACTER(LEN = 40), INTENT(in) :: potential_name

    INTEGER :: no_pot

    IF(potential_name.EQ.free_particle_name)THEN
       Pot=>NULL()
       RETURN
    END IF

    no_pot = 0
    Pot=>Pot_list
    
    DO WHILE(ASSOCIATED(Pot))
       no_pot = no_pot+1
       IF(Pot%name.EQ.potential_name) RETURN
       Pot=>Pot%next
    END DO

    IF(LEN_TRIM(potential_name).GT.0)THEN    
       error_msg = TRIM(potential_name)//" is not a name of a given potential"
       CALL Error(0)
    END IF

    Pot=>Pot_list
    IF(no_pot.GE.2) THEN
       error_msg = "Pass the name of the potential for each method when more than one potential is defined."
       CALL Error(0)
    END IF

  END SUBROUTINE Select_Potential
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  SUBROUTINE Check_Potential()

    IMPLICIT NONE

    REAL(kind = real_kind) :: aux_real
    INTEGER :: aux_int
    REAL(kind = real_kind), parameter :: conv_thrsh = 0.0001_real_kind
    REAL(kind = real_kind), parameter :: back_step = 0.1_real_kind
    REAL(kind = real_kind), parameter :: cent_dif_step = 0.005_real_kind
    integer, parameter :: max_iter = 50

    REAL(kind = real_kind) :: Diff, pot_tmp_prev, pot_tmp
    INTEGER :: ind_min
    REAL(kind = real_kind) :: h_cspl, a_cspl, b_cspl, c_cspl, d_cspl
    REAL(kind = real_kind) :: raiz1, raiz2
    LOGICAL :: first_attempt
    REAL(kind = real_kind) :: tmp
    REAL(kind = real_kind) :: DerivFim
    REAL(kind = real_kind) :: last_energy
    INTEGER :: no_Pts

    ! Unit conversions
    Pot%radius = Unit_Converter(Pot%dist_uni,at_uni_dist,Pot%radius)
    Pot%eq_dist = Unit_Converter(Pot%dist_uni,at_uni_dist,Pot%eq_dist)
    Pot%r_infty = Unit_Converter(Pot%dist_uni,at_uni_dist,Pot%r_infty)
    Pot%zero_energy_ret_p = Unit_Converter(Pot%dist_uni,at_uni_dist,Pot%zero_energy_ret_p)
    Pot%lower_bound = Unit_Converter(Pot%dist_uni,at_uni_dist,Pot%lower_bound)
    Pot%upper_bound = Unit_Converter(Pot%dist_uni,at_uni_dist,Pot%upper_bound)
    IF(Pot%is_dip_mom)THEN
       Pot%force_const = Unit_Converter(Pot%energy_uni,at_uni_dip_mom,Pot%force_const)
    ELSE
       Pot%force_const = Unit_Converter(Pot%energy_uni,at_uni_ener,Pot%force_const)
    END IF

    ! distance conversion factor for inverse units
    aux_real = Unit_Converter(Pot%dist_uni,at_uni_dist,one)

    Pot%force_const = Pot%force_const/aux_real**2
    Pot%Beta = Pot%Beta/aux_real

    IF(Pot%is_dip_mom)THEN
       Pot%epsilon = Unit_Converter(Pot%energy_uni,at_uni_dip_mom,Pot%epsilon)
    ELSE
       Pot%epsilon = Unit_Converter(Pot%energy_uni,at_uni_ener,Pot%epsilon)
    END IF

    IF(Pot%is_dip_mom)THEN
       Pot%ener_infty = Unit_Converter(Pot%energy_uni,at_uni_dip_mom,Pot%ener_infty)
    ELSE
       Pot%ener_infty = Unit_Converter(Pot%energy_uni,at_uni_ener,Pot%ener_infty)
    END IF
    Pot%r_infty = Unit_Converter(Pot%dist_uni,at_uni_dist,Pot%r_infty)

    IF(Pot%is_dip_mom)THEN
       Pot%coef_6 = Unit_Converter(Pot%energy_uni,at_uni_dip_mom,Pot%coef_6)
    ELSE
       Pot%coef_6 = Unit_Converter(Pot%energy_uni,at_uni_ener,Pot%coef_6)
    END IF
    Pot%coef_6 = Pot%coef_6*aux_real**6

    IF(Pot%is_dip_mom)THEN
       Pot%coef_8 = Unit_Converter(Pot%energy_uni,at_uni_dip_mom,Pot%coef_8)
    ELSE
       Pot%coef_8 = Unit_Converter(Pot%energy_uni,at_uni_ener,Pot%coef_8)
    END IF
    Pot%coef_8 = Pot%coef_8*aux_real**8

    IF(ALLOCATED(Pot%coef))THEN
       DO aux_int = 1,SIZE(Pot%coef)
          Pot%coef(aux_int) = Pot%coef(aux_int)/aux_real**aux_int
       END DO
    END IF

    IF(ALLOCATED(Pot%r))THEN
       no_pts = SIZE(Pot%r)
    END IF

    IF(ALLOCATED(Pot%V))THEN
       ! last energy in user units
       IF(Pot%is_dip_mom)THEN
          last_energy = zero
       ELSE
          last_energy = Pot%V(no_pts)
       END IF
       
       IF(TRIM(Pot%type).EQ."SplCubFLJ68")THEN ! already a solution for equation for tail r-6 and r-8
          last_energy = Pot%ener_infty
       END IF
       
       DO aux_int = 1,no_pts,1
          Pot%V(aux_int) = Pot%V(aux_int) - last_energy
          Pot%r(aux_int) = Unit_Converter(Pot%dist_uni,at_uni_dist,Pot%r(aux_int))
          IF(Pot%is_dip_mom)THEN
             Pot%V(aux_int) = Unit_Converter(Pot%energy_uni,at_uni_dip_mom,Pot%V(aux_int))
          ELSE
             Pot%V(aux_int) = Unit_Converter(Pot%energy_uni,at_uni_ener,Pot%V(aux_int))
          END IF
       END DO
    END IF

    SELECT CASE(Pot%type)

       ! Potentials based on tabulated points
    CASE("LinearSpl","NatCubSpl","CubSplr6","CubSplr6r8","CubSplFitr6","CubSplFitr6r8","CubSplFitr6r8r10",&
         "LSLenJon","LSLenJonNM","AppMorse")
       
       CALL Get_Z_Nat_Cub_Spl (Pot%r,Pot%V,Pot%z)

       ind_min = 1

       SELECT CASE(Pot%type)

       CASE("LinearSpl","NatCubSpl","CubSplr6","CubSplr6r8","CubSplFitr6","CubSplFitr6r8","CubSplFitr6r8r10")

          Pot%lower_bound = Pot%R(1)

          Pot%upper_bound = Pot%R(SIZE(Pot%r))

          IF(Pot%type.EQ."CubSplFitr6") THEN

             aux_real = Pot%coef_6/Pot%R(no_Pts)**6
             DO aux_int = 1,no_Pts,1
                Pot%V(aux_int) = Pot%V(aux_int) + aux_real
                ! Since Pot%z depends only of the difference of ajacent points in V, it is not needed call Get_Z_Nat_Cub_Spl again
             END DO

          END IF

          IF(Pot%type.EQ."CubSplFitr6r8") THEN

             aux_real = Pot%coef_6/Pot%R(no_Pts)**6 + Pot%coef_8/Pot%R(no_Pts)**8
             DO aux_int = 1,no_Pts,1
                Pot%V(aux_int) = Pot%V(aux_int) + aux_real
                ! Since Pot%z depends only of the difference of ajacent points in V, it is not needed call Get_Z_Nat_Cub_Spl again
             END DO

          END IF

          IF(Pot%type.EQ."CubSplFitr6r8r10") THEN

             aux_real = Pot%R(no_Pts)**2
             aux_real = Pot%coef_6/aux_real**3 + Pot%coef_8/aux_real**4 + Pot%coef_10/aux_real**5
             DO aux_int = 1,no_Pts,1
                Pot%V(aux_int) = Pot%V(aux_int) + aux_real
                ! Since Pot%z depends only of the difference of ajacent points in V, it is not needed call Get_Z_Nat_Cub_Spl again
             END DO

          END IF

          IF(Pot%type.EQ."CubSplr6") THEN

             tmp = Natural_Cubic_Spline(Pot%R(no_Pts),Pot%R,Pot%V,Pot%z,DerivFim)

             Pot%coef_6=-DerivFim*Pot%R(no_Pts)**7/6

             DO aux_int = 1,no_Pts,1
                Pot%V(aux_int) = Pot%V(aux_int) - (DerivFim*Pot%R(no_Pts) + 6*Pot%V(no_Pts))/6
                ! Since Pot%z depends only of the difference of ajacent points in V, it is not needed call Get_Z_Nat_Cub_Spl again
             END DO

          END IF

          IF(Pot%type.EQ."CubSplr6r8") THEN

             tmp = Natural_Cubic_Spline(Pot%R(no_Pts),Pot%R,Pot%V,Pot%z,DerivFim)

             tmp = Pot%R(no_Pts)**6
             Pot%coef_6 = tmp*(8*Pot%V(no_Pts) + DerivFim*Pot%R(no_Pts))/2

             tmp = tmp*(Pot%R(no_Pts)**2)
             Pot%coef_8 = -tmp*(6*Pot%V(no_Pts) + DerivFim*Pot%R(no_Pts))/2

          END IF

       CASE ("LSLenJon","LSLenJonNM")

          if(Pot%type.eq."LSLenJon")then
             Pot%n_lenjon = 6
             Pot%m_lenjon = 6
          end if

          CALL Least_Square_Len_Jon()

       CASE("AppMorse")

       CASE DEFAULT

          IF(Pot%is_dip_mom)THEN
             CALL Non_Consistent_Info("type of dipole moment","internal definition of data")
          ELSE
             CALL Non_Consistent_Info("type of potential","internal definition of data")
          END IF

       END SELECT

       IF(.NOT.Pot%is_dip_mom)THEN ! Find the minimum and zero energy turnig point for data points potentials

          DO ind_min = 1,no_pts-1,1
             IF(Pot%V(ind_min).LT.Pot%V(ind_min+1)) exit
          END DO
          
          IF(ind_min.GE.no_Pts.OR.ind_min.LE.1)THEN ! No minimum
             Pot%radius = -one
             Pot%epsilon = -one

             Pot%eq_dist = -one
             Pot%diss_ener = zero
             Pot%force_const = -one
             
             Pot%abs_min_energy = Pot%V(ind_min) + last_energy
             
          ELSE
             
             IF((Pot%type).EQ."LinearSpl")THEN
                
                Pot%radius = Pot%R(ind_min)
                Pot%epsilon = -Pot%V(ind_min)

                Pot%eq_dist = Pot%radius
                Pot%diss_ener = Pot%epsilon

                ! Force constant = second derivative, difecences method
                Pot%force_const = (Pot%V(ind_min+1) - Pot%V(ind_min))/(Pot%r(ind_min+1) - Pot%r(ind_min))
                Pot%force_const = Pot%force_const - (Pot%V(ind_min) - Pot%V(ind_min-1))/(Pot%r(ind_min) - Pot%r(ind_min-1))
                Pot%force_const = 2*Pot%force_const/(Pot%r(ind_min+1) - Pot%r(ind_min-1))

                Pot%abs_min_energy = Pot%V(ind_min) + last_energy
                
             ELSE
                
                ind_min = ind_min-1
                first_attempt = .TRUE.

                do

                   h_cspl = Pot%r(ind_min+1) - Pot%r(ind_min)
                   a_cspl = (Pot%z(ind_min+1) - Pot%z(ind_min))/(6*h_cspl)
                   b_cspl = Pot%z(ind_min)/2
                   c_cspl = -(h_cspl/6)*(Pot%z(ind_min+1) + 2*Pot%z(ind_min)) + (Pot%V(ind_min+1) - Pot%V(ind_min))/h_cspl
                   d_cspl =  Pot%V(ind_min)
                   
                   raiz1 = (-b_cspl + SQRT(b_cspl**2 - 3*a_cspl*c_cspl))/(3*a_cspl) + Pot%r(ind_min)
                   raiz2 = (-b_cspl - SQRT(b_cspl**2 - 3*a_cspl*c_cspl))/(3*a_cspl) + Pot%r(ind_min)
                   
                   IF(raiz1.GT.Pot%r(ind_min).AND.raiz1.LT.Pot%r(ind_min+1))THEN

                      IF(raiz2.GT.Pot%r(ind_min).AND.raiz2.LT.Pot%r(ind_min+1))THEN
                         error_msg = "Two roots in the interval of Cubic Spline, for calculate eq_distance."
                         CALL Error(1)
                      END IF

                      Pot%radius = raiz1
                      exit

                   ELSE

                      IF(raiz2.GT.Pot%R(ind_min).AND.raiz2.LT.Pot%R(ind_min+1))THEN
                         Pot%radius = raiz2
                         exit
                      ELSE

                         IF(first_attempt)THEN
                            first_attempt = .FALSE.
                            ind_min = ind_min+1
                         else
                            error_msg = "No roots in the interval of Cubic Spline, for calculate eq_distance."
                            CALL Error(2)
                         END IF

                      END IF
                   END IF
                   
                end do
                
                ! Force constant = cubic spline second derivative
                Pot%force_const = 2*(3*a_cspl*(Pot%radius - Pot%r(ind_min)) + b_cspl)

                ! Cubic spline minimum energy
                Pot%abs_min_energy = b_cspl + (Pot%eq_dist - Pot%r(ind_min))*a_cspl
                Pot%abs_min_energy = c_cspl + (Pot%eq_dist - Pot%r(ind_min))*Pot%abs_min_energy
                Pot%abs_min_energy = d_cspl + (Pot%eq_dist - Pot%r(ind_min))*Pot%abs_min_energy + last_energy
                
                ! Variables for AppMorse
                Pot%epsilon = b_cspl + a_cspl*(Pot%radius - Pot%R(ind_min))
                Pot%epsilon = c_cspl + Pot%epsilon*(Pot%radius - Pot%R(ind_min))
                Pot%epsilon = -(d_cspl + Pot%epsilon*(Pot%radius - Pot%R(ind_min)))
                
                Pot%eq_dist = Pot%radius
                Pot%diss_ener = Pot%epsilon

                Pot%Beta = SQRT((6*a_cspl*(Pot%radius-Pot%R(ind_min))+2*b_cspl)/(2*Pot%epsilon))
                
             END IF
             
             
          END IF

          if(Pot%V(1).LE.Pot%V(no_pts))then ! first point below the last
             
             Pot%zero_energy_ret_p = Pot%r(1)
             
          ELSE
             
             DO aux_int = 1,no_pts-1,1
                IF(Pot%V(aux_int+1).LT.0) EXIT
             END DO

             IF((Pot%type).EQ."LinearSpl")THEN
                
                Pot%zero_energy_ret_p = (Pot%r(aux_int)*Pot%V(aux_int+1) - Pot%V(aux_int)*Pot%r(aux_int+1))/&
                     (Pot%V(aux_int+1) - Pot%V(aux_int))
       
             ELSE
                
                Pot%zero_energy_ret_p = Newton_Raphson(Potential,(Pot%r(aux_int)+Pot%r(aux_int+1))/2)
                IF(status.GT.0)THEN
                   WRITE(uni_log,FMT = '("Maximum number of iterations reached, calculating zero energy return point "'//&
                        '"for ",A," potential.")') trim(Pot%type)//"-"//trim(Pot%name)
                END IF


             END IF
             
          END IF

       END IF
       
    END SELECT

    SELECT CASE(Pot%type)

    CASE("HardSphere")

       Pot%zero_energy_ret_p = Pot%radius
       Pot%lower_bound = Pot%radius
       Pot%upper_bound = Pot%radius
       Pot%eq_dist = -one
       Pot%abs_min_energy = zero
       Pot%diss_ener = -one
       Pot%force_const = -one

    CASE("SoftSphere")
       
       Pot%zero_energy_ret_p = Pot%radius
       Pot%lower_bound = zero
       Pot%upper_bound = Pot%radius
       Pot%eq_dist = -one
       Pot%abs_min_energy = zero
       Pot%diss_ener = -one
       Pot%force_const = -one

    CASE("Square")

       IF(Pot%epsilon.LT.0) THEN
          Pot%zero_energy_ret_p = zero
          Pot%eq_dist = Pot%radius
          Pot%abs_min_energy = Pot%epsilon
          Pot%diss_ener = -Pot%epsilon
       ELSE
          Pot%zero_energy_ret_p = Pot%radius
          Pot%eq_dist = -one
          Pot%abs_min_energy = zero
          Pot%diss_ener = -one
       END IF

       Pot%lower_bound = zero
       Pot%upper_bound = Pot%radius

       Pot%force_const=zero

    CASE("LenJon","LSLenJon")

       Pot%zero_energy_ret_p = Pot%radius/(two**(one/6))
       Pot%lower_bound = zero
       Pot%upper_bound = -one
       Pot%eq_dist = Pot%radius
       Pot%abs_min_energy = -Pot%epsilon
       Pot%diss_ener = Pot%epsilon
       Pot%force_const = 72*Pot%epsilon/(Pot%radius**2)
       
    CASE("LenJonNM","LSLenJonNM")

       if(Pot%type.EQ."LenJonNM")then
          ! Passed values 
          Pot%eq_dist = Pot%radius
          Pot%diss_ener = Pot%epsilon
          
          ! Variables used in calculation
          Pot%radius = Pot%radius*((Pot%n_lenjon*one)/(Pot%n_lenjon+Pot%m_lenjon))**(one/Pot%m_lenjon)
          Pot%epsilon = (Pot%epsilon*(Pot%n_lenjon+Pot%m_lenjon))/(Pot%m_lenjon)*&
               ((Pot%n_lenjon+Pot%m_lenjon)*one/(Pot%n_lenjon))**(one*Pot%n_lenjon/Pot%m_lenjon)
       else

          ! The values from least square are the working radius and epsilon variables
          Pot%eq_dist = Pot%radius*((Pot%n_lenjon+Pot%m_lenjon)*one/(Pot%n_lenjon))**(one/Pot%m_lenjon)
          Pot%diss_ener = (Pot%epsilon*Pot%m_lenjon)/(Pot%n_lenjon+Pot%m_lenjon)*&
               (Pot%n_lenjon*one/(Pot%n_lenjon+Pot%m_lenjon))**(one*Pot%n_lenjon/Pot%m_lenjon)

       end if

       Pot%zero_energy_ret_p = Pot%radius
       Pot%lower_bound = zero
       Pot%upper_bound = -one
       Pot%abs_min_energy = -Pot%diss_ener
       Pot%force_const = Pot%epsilon*Pot%m_lenjon*(2*Pot%m_lenjon + Pot%m_lenjon + 1)/(Pot%radius**2)

    CASE("Morse","AppMorse")

       Pot%zero_energy_ret_p = Pot%radius - LOG(two)/Pot%beta
       Pot%lower_bound = zero
       Pot%upper_bound = -one

       Pot%eq_dist = Pot%radius
       Pot%abs_min_energy = -Pot%epsilon
       Pot%diss_ener = Pot%epsilon
       Pot%force_const = 2*Pot%epsilon*(Pot%beta**2)

    CASE("MurSorbie")

       Pot%eq_dist = Pot%radius
       Pot%diss_ener = Pot%epsilon
       Pot%abs_min_energy = -Pot%diss_ener

       ! Numerical search for zero point return point
       IF(.NOT.Pot%is_dip_mom)THEN

          Pot%zero_energy_ret_p=Pot%radius
          do
             Pot%zero_energy_ret_p = Pot%zero_energy_ret_p - back_step
             if(Potential(Pot%zero_energy_ret_p).GT.zero.OR.Pot%zero_energy_ret_p.LT.zero) exit
          end do

          Pot%zero_energy_ret_p = Newton_Raphson(Potential,Pot%zero_energy_ret_p)
          IF(status.GT.0)THEN
             WRITE(uni_log,FMT = '("Maximum number of iterations reached, calculating zero energy return point "'//&
                  '"for ",A," potential.")') trim(Pot%type)//"-"//trim(Pot%name)
          END IF

          ! Numerical search for point in repulsive region of potential where derivative becomes positive (incorrect 
          ! description of Murrel sorbie potential
          Diff = -one
          pot_tmp_prev = Potential(Pot%zero_energy_ret_p)
          Pot%lower_bound = Pot%zero_energy_ret_p
          DO WHILE(Diff.LT.0.AND.Pot%lower_bound.GT.0)
             
             Pot%lower_bound = Pot%lower_bound - back_step
             pot_tmp = Potential(Pot%lower_bound)
             
             Diff = pot_tmp_prev - pot_tmp
             
             pot_tmp_prev = pot_tmp
             
          END DO

          Pot%lower_bound = Pot%lower_bound + back_step
          Pot%upper_bound = -one

          Pot%force_const = (Potential(Pot%eq_dist+cent_dif_step) - 2*Potential(Pot%eq_dist) + &
               Potential(Pot%eq_dist-cent_dif_step))/(cent_dif_step**2)

       END IF

    CASE("Direct")
       
    case("LinearSpl","NatCubSpl","CubSplr6","CubSplr6r8","CubSplFitr6","CubSplFitr6r8","CubSplFitr6r8r10")
       
    CASE DEFAULT

       IF(Pot%is_dip_mom)THEN
          CALL Non_Consistent_Info("type of dipole moment","external definition of data")
       ELSE
          CALL Non_Consistent_Info("type of potential","external definition of data")
       END IF

    END SELECT

!!$    write(*,FMT='("      DATA (X_TAB(I),I=1,",I2,")/")') size(Pot%r)
!!$    do aux_int=1,size(Pot%r)-1
!!$       write(*,FMT='("     ",I1,5X,F20.13,",")/")') MOD(aux_int,10),Pot%r(aux_int)
!!$    end do
!!$    aux_int=size(Pot%r)
!!$    write(*,FMT='("     ",I1,5X,F20.13,"/")/")') MOD(aux_int,10),Pot%r(aux_int)
!!$
!!$    write(*,FMT='("      DATA (Y_TAB(I),I=1,",I2,")/")') size(Pot%V)
!!$    do aux_int=1,size(Pot%V)-1
!!$       write(*,FMT='("     ",I1,5X,F20.13,",")/")') MOD(aux_int,10),Pot%V(aux_int)
!!$    end do
!!$    aux_int=size(Pot%V)
!!$    write(*,FMT='("     ",I1,5X,F20.13,"/")/")') MOD(aux_int,10),Pot%V(aux_int)
!!$
!!$    write(*,FMT='("      DATA (D2YDX2_TAB(I),I=1,",I2,")/")') size(Pot%z)
!!$    do aux_int=1,size(Pot%z)-1
!!$       write(*,FMT='("     ",I1,5X,F20.13,",")/")') MOD(aux_int,10),Pot%z(aux_int)
!!$    end do
!!$    aux_int=size(Pot%z)
!!$    write(*,FMT='("     ",I1,5X,F20.13,"/")/")') MOD(aux_int,10),Pot%z(aux_int)

  CONTAINS
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
    SUBROUTINE Get_Z_Nat_Cub_Spl(x_tab,y_tab,d2ydx2_tab)

      IMPLICIT NONE

      REAL(kind = real_kind), DIMENSION(:), INTENT(in) :: x_tab, y_tab
      REAL(kind = real_kind), ALLOCATABLE, DIMENSION(:), INTENT(out) :: d2ydx2_tab

      REAL(kind = 8), ALLOCATABLE, DIMENSION(:) :: h, b, u, v
      INTEGER :: i
      INTEGER :: no_Pts

      ! Verifica tamanho
      IF(SIZE(x_tab).NE.SIZE(y_tab))THEN
         error_msg="Error in vector sizes in Cubic Spline."
         CALL Error(2)
      END IF
      no_Pts=SIZE(x_tab)

      IF(ALLOCATED(d2ydx2_tab)) THEN
         DEALLOCATE(d2ydx2_tab)
      END IF

      ALLOCATE(d2ydx2_tab(no_Pts))

      ALLOCATE(h(no_Pts-1),b(no_Pts-1),u(2:no_Pts-1),v(2:no_Pts-1))

      DO i = 1,no_Pts-1,1
         h(i) = x_tab(i+1) - x_tab(i)
         b(i) = (y_tab(i+1) - y_tab(i))/h(i)
      END DO

      u(2) = 2*(h(1) + h(2))
      v(2) = 6*(b(2) - b(1))

      DO i = 3,no_Pts-1,1
         u(i) = 2*(h(i)+h(i-1)) - h(i-1)**2/u(i-1)
         v(i) = 6*(b(i)-b(i-1)) - h(i-1)*v(i-1)/u(i-1)
      END DO

      d2ydx2_tab(no_Pts) = 0

      DO i = no_Pts-1,2,-1
         d2ydx2_tab(i) = (v(i)-h(i)*d2ydx2_tab(i+1))/u(i)
      END DO

      d2ydx2_tab(1) = 0

      DEALLOCATE(h,b,u,v)

    END SUBROUTINE Get_Z_Nat_Cub_Spl
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  END SUBROUTINE Check_Potential
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  SUBROUTINE Check_Unit(unit,info)

    IMPLICIT NONE

    CHARACTER (LEN = units_len), INTENT(in) :: unit
    CHARACTER(len = *), INTENT(in) :: info

    IF(status.EQ.1.OR.status.EQ.17)THEN
       error_msg = TRIM(unit)//" has an unknow prefix."
       CALL Error(0)
    END IF
    IF(status.EQ.20)THEN
       error_msg = TRIM(unit)//" is an unknow unit."
       CALL Error(0)
    END IF
    IF(status.EQ.21)THEN
       error_msg = TRIM(unit)//" is an unknow unit with an unknow prefix."
       CALL Error(0)
    END IF
    IF(status.EQ.16)THEN
       error_msg = TRIM(unit)//" is not unit of "//TRIM(info)//"."
       CALL Error(0)
    END IF
    IF(status.NE.0)THEN
       WRITE(unit=error_msg, FMT='("Inconsistency of unit: ",A," - ",A,". Unit converter status:",I2)')&
            TRIM(unit), TRIM(info), status
       CALL Error(2)
    END IF

  END SUBROUTINE Check_Unit
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  SUBROUTINE Least_Square_Len_Jon()

    IMPLICIT NONE

    INTEGER :: aux_int
    CHARACTER(LEN=1) :: TRANS_lapack
    INTEGER :: M_lapack, N_lapack, NRHS_lapack, LDA_lapack, LDB_lapack, LWORK_lapack, INFO_lapack
    REAL(kind = real_kind), ALLOCATABLE, DIMENSION(:) :: B_lapack, WORK_lapack
    REAL(kind = real_kind), ALLOCATABLE, DIMENSION(:,:) :: A_lapack

    m_lapack = size(Pot%V)
    n_lapack=2
    ALLOCATE(A_lapack(M_lapack,N_lapack),B_lapack(M_lapack))

    DO aux_int = 1,m_lapack,1
       B_lapack(aux_int) = Pot%V(aux_int)

       A_lapack(aux_int,2) = -1/Pot%R(aux_int)**Pot%n_lenjon
       A_lapack(aux_int,1) = -A_lapack(aux_int,2)/Pot%R(aux_int)**Pot%m_lenjon

    END DO

    TRANS_lapack = 'N'

    NRHS_lapack = 1
    LDA_lapack = MAX(1,M_lapack)
    LDB_lapack = MAX(1,M_lapack,N_lapack)
    LWORK_lapack = max( 1, M_lapack*N_lapack + max( M_lapack*N_lapack, NRHS_lapack ) )
    ALLOCATE(WORK_lapack(MAX(1,LWORK_lapack)))
    
    CALL DGELS( TRANS_lapack, M_lapack, N_lapack, NRHS_lapack, A_lapack, LDA_lapack, B_lapack, LDB_lapack, &
         WORK_lapack, LWORK_lapack, INFO_lapack )

    if(Pot%type.eq."LSLenJon")then
       Pot%radius = EXP(LOG(2*B_lapack(2)/B_lapack(1))/6)
       Pot%epsilon = (B_lapack(1)**2)/(4*B_lapack(2))
    else
       Pot%radius = EXP(LOG(B_lapack(2)/B_lapack(1))/6)
       Pot%epsilon = (B_lapack(1)**2)/B_lapack(2)
    end if

  END SUBROUTINE Least_Square_Len_Jon
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  REAL(kind = real_kind) FUNCTION Potential(r)

    IMPLICIT NONE

    REAL(kind = real_kind), INTENT(in) :: r       ! Local onde será calculado o potencial

    INTEGER :: aux_int                              ! Contador
    REAL(kind = real_kind) :: tmp, tmp2, ro          ! Auxiliar para os cálculos
    CHARACTER(LEN=30) :: saidaPad                 ! Saida Padrão do system_call
    TYPE(TypeStringLinkedList), POINTER :: current_string =>NULL()
    INTEGER :: ind, ind_ant

!!$real(kind = real_kind) :: a_inten,b_inten,c_inten,d_inten,r_inten

!    INTEGER :: i
!    REAL( kind = real_kind) :: V_sh, V_as

    IF(use_Pot_vector)THEN
       
       Potential=curr_pot_knot%value
       curr_pot_knot => curr_pot_knot%next
       status = 0
       
       RETURN
    END IF
    
!!Intensity potential
!!! Dont forget to define de diss energy
!!$    if(Pot%name.eq."estado_1")then
!!$       if(r.lt.0.26000000000D+01)then
!!$          error_msg = "Lower than lower bound - intensity"
!!$          CALL Error(2)
!!$       end if
!!$       ro = r
!!$       if(r.lt.0.31000000000D+01)then
!!$          r_inten = 0.26000000000D+01
!!$          a_inten = -0.16450669751D+01
!!$          b_inten = -0.56048531557D+00
!!$          c_inten = 0.77854248634D+00
!!$          d_inten = -0.38357310668D+00
!!$          
!!$       elseif(r.lt.0.36090021065D+01)then
!!$          r_inten = 0.31000000000D+01
!!$          a_inten = -0.17786206497D+01
!!$          b_inten = -0.69622659236D-01
!!$          c_inten = 0.20318282633D+00
!!$          d_inten = -0.10528820918D+00
!!$          
!!$       elseif(r.lt.0.45080000000D+01)then
!!$          r_inten = 0.36090021065D+01
!!$          a_inten = -0.17753022475D+01
!!$          b_inten = 0.55383113060D-01
!!$          c_inten = 0.42407065520D-01
!!$          d_inten = -0.26729071769D-01
!!$          
!!$       else
!!$          r_inten = 0.45080000000D+01
!!$          a_inten = -0.17106601507D+01
!!$          b_inten = 0.66823754314D-01
!!$          c_inten = -0.29681072124D-01
!!$          d_inten = 0.44781826775D-02
!!$          
!!$          if(r.gt.0.70000000000D+01) ro = 0.70000000000D+01
!!$       end if
!!$       Potential = a_inten + (ro-r_inten)*(b_inten + (ro-r_inten)*(c_inten + (ro-r_inten)*d_inten)) + 1.711196626183661_real_kind
!!$       return
!!$!(+ -1.7753022475 (* 0.8909979 0.055383113060 ) (* 0.793877258 0.042407065520 ) (* 0.70734297 -0.026729071769 ) )
!!$    end if




    IF(.NOT.ASSOCIATED(Pot))THEN
       Potential = 0.0_real_kind
    ELSE

       IF(r.LE.Pot%lower_bound-Pot%TolLimInf)THEN
          if(stop_error)then
             error_msg = "Attempt to calculate the potential below its lower bound."
             CALL Error(1)
          else
             status = 1
             return
          end if
       END IF

       SELECT CASE(Pot%type)

          ! HardSphere
          !
          ! V(r) = infty   if r < radius
          !        0       if r > radius
          !  
       CASE("HardSphere")

          Potential = 0.0_real_kind

          ! SoftSphere
          !
          ! V(r) = epsilon*((radius/r)^rep_exp -1)  if r < radius
          !        0                                if r > radius
          !
       CASE("SoftSphere")

          IF(r.GE.Pot%radius)THEN
             Potential = 0.0_real_kind
          ELSE
             Potential = Pot%epsilon*((Pot%radius/r)**Pot%rep_exp - 1)
          END IF

          ! Square
          ! 
          ! V(r) =  epsilon   if r < radius
          !         0         if r > radius
          !
       CASE("Square")

          IF(r.GE.Pot%radius)THEN
             Potential = 0.0_real_kind
          ELSE
             Potential = Pot%epsilon
          END IF

          ! Lennard Jones (12,6)
          !
          ! V(r) = epsilon*((radius/r)^12 - 2*(radius/r)**6)
          !
       CASE("LenJon","LSLenJon")

          tmp = (Pot%radius/r)**6
          tmp = tmp**2 - 2*tmp
          Potential = Pot%epsilon*tmp

          ! Lennard Jones (N+M,N)
          !
          ! V(r) = epsilon*((radius/r)^(N+M) - (radius/r)**N)
          !
       CASE("LenJonNM","LSLenJonNM")

          tmp = (Pot%radius/r)**Pot%n_lenjon
          tmp2 = (Pot%radius/r)**Pot%m_lenjon

          Potential = Pot%epsilon*(tmp*tmp2 - tmp)

          ! Morse
          !
          ! V(r) = epsilon*((1-exp(-beta*(r-radius)))^2 - 1)
          !
       CASE("Morse","AppMorse")

          Potential = Pot%epsilon*((1 - EXP(-Pot%Beta*(r - Pot%radius)))**2 - 1)

          ! Murrell Sorbie
          !
          ! V(r) = -epsilon(1+sum_{i=1}^{n} a_i(r-Re)^i)exp(-a_1(r-Re))
          !
       CASE("MurSorbie")

!!! Potential of Zhang et al
! alpha = epsilon
! beta = radius
! gamma = coef(1)
! lambda = coef(2)
! a_0 = coef(3)
! a_1 = coef(4)
! a_2 = coef(5)
! a_3 = coef(6)
! a_4 = coef(7)
! a_5 = coef(8)
! a_6 = coef(9)
! a_7 = coef(10)
! a_8 = coef(11)
! C_6 = coef(12)
! C_8 = coef(13)
! C_10 = coef(14)
! C_12 = coef(15)
! C_14 = coef(16)
!
!!$          V_sh = zero
!!$
!!$          ro = one
!!$          do i=0,8,1
!!$             V_sh = V_sh + Pot%coef(i+3)*ro*exp(-Pot%epsilon*(r-Pot%radius))
!!$             ro = ro*r
!!$          end do
!!$
!!$!          V_sh = V_sh
!!$
!!$          V_as = zero
!!$          
!!$          ro = r**6
!!$          do i=0,4,1
!!$             V_as = V_as + Pot%coef(i+12)/ro
!!$             ro = ro*(r**2)
!!$          end do
!!$
!!$          V_as = (one + tanh(Pot%coef(1) + Pot%coef(2)*R))*V_as/2
!!$
!!$          Potential = V_sh + V_as
!!$
          ro = r - Pot%radius

          tmp = 1+Pot%coef(1)*ro

          DO aux_int = 2,SIZE(Pot%coef),1
             tmp = tmp + Pot%coef(aux_int)*ro**aux_int
          END DO

          Potential = -Pot%epsilon*tmp*EXP(-Pot%coef(1)*ro)

          !
          ! Linear spline
          !
       CASE("LinearSpl")

          IF (r.GE.Pot%R(SIZE(Pot%r))) THEN
             Potential = Pot%V(SIZE(Pot%r))
          ELSE
             Potential = Linear_Spline(r,Pot%R,Pot%V)
          END IF

          !
          ! Natural Cubic Spline, Natural Cubic Spline joined with a c_6/r⁶ tail, and Natural Cubic Spline joined with a 
          ! c_6/r⁶ + c_8/r⁸ tail
          !
       CASE("NatCubSpl","CubSplr6","CubSplr6r8","CubSplFitr6","CubSplFitr6r8","CubSplFitr6r8r10")

          IF (r.GE.Pot%R(SIZE(Pot%r))) THEN

             SELECT CASE(Pot%type)

             CASE("NatCubSpl")

                Potential = Pot%V(SIZE(Pot%r))

             CASE("CubSplr6","CubSplFitr6")
                Potential = Pot%coef_6/r**6

             CASE("CubSplr6r8","CubSplFitr6r8")
                tmp = r**2
                Potential = Pot%coef_6/tmp**3 + Pot%coef_8/tmp**4

             CASE("CubSplFitr6r8r10")
                tmp = r**2
                Potential = Pot%coef_6/tmp**3 + Pot%coef_8/tmp**4 + Pot%coef_10/tmp**5

             CASE DEFAULT

                CALL Non_Consistent_Info("type of potential","SplCubcalculation")

             END SELECT

          ELSE
             Potential = Natural_Cubic_Spline(r,Pot%R,Pot%V,Pot%z)
          END IF

          !
          ! On the fly calculation of potential
          !
       CASE ("Direct")

          OPEN(Unit=Pot%uni_EE_input,File='INPUT',status='NEW')
          
          IF(r.GT.Pot%r_infty)THEN
             Potential=zero
             CLOSE(Pot%uni_EE_input)
          ELSE
             
             current_string=>Pot%input
             DO WHILE(ASSOCIATED(current_string))
                
                ind_ant = -LEN_TRIM(Pot%keyword)+1
                DO
                   ind = INDEX(current_string%value(ind_ant+LEN_TRIM(Pot%keyword):),TRIM(Pot%keyword)) &
                        &+ ind_ant + LEN_TRIM(Pot%keyword) - 1
                   IF(ind.GT.ind_ant+LEN_TRIM(Pot%keyword)-1) THEN
                      WRITE(Pot%uni_EE_input,FMT=TRIM(Pot%format_r),advance='no') &
                           TRIM(current_string%value(ind_ant+LEN_TRIM(Pot%keyword):ind-1)),r
                   ELSE
                      WRITE(Pot%uni_EE_input,FMT='(A)') TRIM(current_string%value(ind_ant+LEN_TRIM(Pot%keyword):))
                      EXIT
                   END IF
                   ind_ant = ind
                END DO
                
                current_string=>current_string%next
             END DO
             
             CLOSE(Pot%uni_EE_input)
             
             CALL system(TRIM(Pot%ES_command))
             
             CALL system_call_stdout_c2f(TRIM(Pot%get_ener_command_c),saidaPad)
             READ(saidaPad,FMT=TRIM(Pot%format_ener),iostat=ios_var) Potential
             IF(ios_var.NE.0.OR.LEN_TRIM(saidaPad).EQ.0)THEN
                error_msg = "Error in getting energy. Check commands, input and formats."
                CALL Error(0)
             END IF
             
             Potential = Potential - Pot%ener_infty

          END IF

          CALL system(TRIM(Pot%rm_command))
          
       CASE DEFAULT

          CALL Non_Consistent_Info("type of potential","calculation of potential")

       END SELECT

    END IF

    IF(save_Potential)then
       curr_pot_knot%value = Potential
       allocate(curr_pot_knot%next)
       curr_pot_knot => curr_pot_knot%next
    end IF

    status = 0

  END FUNCTION Potential
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  REAL(kind = real_kind) FUNCTION Linear_Spline(x,x_tab,y_tab,dydx)

    REAL(kind = real_kind), INTENT(in) :: x
    REAL(kind = real_kind), DIMENSION(:), INTENT(in) :: x_tab, y_tab
    REAL(kind = real_kind), INTENT(out), OPTIONAL :: dydx

    REAL(kind = real_kind) :: der
    INTEGER :: i
    INTEGER :: no_Pts

    no_Pts=SIZE(x_tab)

    i=Ordered_Search(x,x_tab)

    IF(i.EQ.SIZE(x_tab)) i=i-1
    IF(i.EQ.0) i=i+1

    ! Cálculo do spline
    der = (y_tab(i+1) - y_tab(i))/(x_tab(i+1) - x_tab(i))
    Linear_Spline = y_tab(i) + (x - x_tab(i))*der

    IF (PRESENT(dydx)) dydx = der

  END FUNCTION Linear_Spline
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  REAL(kind = real_kind) FUNCTION Natural_Cubic_Spline(x,x_tab,y_tab,d2ydx2_tab,dydx,d2ydx2,d3ydx3)

    IMPLICIT NONE

    REAL(kind = real_kind), INTENT(in) :: x
    REAL(kind = real_kind), DIMENSION(:), INTENT(in), TARGET :: x_tab, y_tab, d2ydx2_tab
    REAL(kind = real_kind), INTENT(out), OPTIONAL :: dydx, d2ydx2, d3ydx3

    REAL(kind = real_kind) :: h,D,B,TEMP
    INTEGER :: i
    INTEGER :: no_Pts

    no_Pts=SIZE(x_tab)

    i=Ordered_Search(x,x_tab)

    IF(i.EQ.SIZE(x_tab)) i=i-1
    IF(i.EQ.0) i=i+1

    ! Cálculo do spline
    h = x_tab(i+1) - x_tab(i)

    D = (d2ydx2_tab(i+1) - d2ydx2_tab(i))/h
    B = -(h/6)*(d2ydx2_tab(i+1) + 2*d2ydx2_tab(i)) + (y_tab(i+1) - y_tab(i))/h

    TEMP = (d2ydx2_tab(i)/2) + (x - x_tab(i))*D/6
    TEMP = B + (x - x_tab(i))*TEMP
    Natural_Cubic_Spline = y_tab(i) + (x - x_tab(i))*TEMP

    ! Calcula derivadas, se pedido
    IF(PRESENT(dydx)) THEN
       TEMP = d2ydx2_tab(i) + (x - x_tab(i))*D/2
       dydx = B + (x - x_tab(i))*TEMP
    END IF
    IF(PRESENT(d2ydx2)) THEN
       d2ydx2 = d2ydx2_tab(i) + (x - x_tab(i))*D
    END IF
    IF(PRESENT(d3ydx3)) THEN
       d3ydx3 = D
    END IF

  END FUNCTION Natural_Cubic_Spline
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  SUBROUTINE Classical_Turning_Points(Energy,in,out)

    IMPLICIT NONE

    REAL(kind = real_kind), INTENT(in) :: Energy
    REAL(kind = real_kind), INTENT(out) :: in, out
    
    REAL(kind = real_kind), PARAMETER :: conv_thrsh=1.0E-5_real_kind
    INTEGER, PARAMETER :: max_iter=200

    if(Pot%type.eq."Square")then
       if(Energy_TMP_module.lt.zero)then
          in = zero
          out = pot%radius
       else
          status = 1
       end if
       return
    end if

    Energy_TMP_module = Energy
    in = Newton_Raphson(Dif_Ener_Pot,(Pot%lower_bound+Pot%eq_dist)/2)

    IF(status.NE.0) THEN
       RETURN
    END IF
    if(energy.gt.zero)then
       out = -one
       status = 0
       return
    end if

    out = Newton_Raphson(Dif_Ener_Pot,2*Pot%eq_dist-in)

    if(status.eq.2)then
       out = Newton_Raphson(Dif_Ener_Pot,Pot%upper_bound-eps_NR_default)
    end if

    IF(status.NE.0) THEN
       RETURN
    END IF

    status = 0

  END SUBROUTINE Classical_Turning_Points
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  REAL(kind = real_kind) FUNCTION Dif_Ener_Pot(x)
    
    IMPLICIT NONE
    
    REAL(kind = real_kind), INTENT(in) :: x
    
    Dif_Ener_Pot = Potential(x)-Energy_TMP_module
    
  END FUNCTION Dif_Ener_Pot
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  SUBROUTINE Print_Potential_Data()

    IMPLICIT NONE

    INTEGER :: i
    TYPE(TypeStringLinkedList), POINTER :: current_string

    PRINT*,"Potential Data:"
    Pot=>Pot_list
    DO WHILE(ASSOCIATED(Pot))

       PRINT*,"----------------------------------------------------"
       PRINT*,"name = ",TRIM(Pot%name)
       PRINT*,"type = ",TRIM(Pot%type)
       PRINT*,"is_dip_mom = ",Pot%is_dip_mom

       PRINT*,"VARIABLES FOR EACH POTENTIAL:"
       PRINT*,"rep_exp = ",Pot%rep_exp
       PRINT*,"n_lenjon = ",Pot%n_lenjon
       PRINT*,"m_lenjon = ",Pot%m_lenjon
       PRINT*,"radius = ",Pot%radius
       PRINT*,"Beta = ",Pot%Beta
       PRINT*,"epsilon = ",Pot%epsilon
       PRINT*,"coef_6 = ",Pot%coef_6
       PRINT*,"coef_8 = ",Pot%coef_8
       PRINT*,"coef_10 = ",Pot%coef_10
       PRINT*,"r_infty = ",Pot%ener_infty
       PRINT*,"TolLimInf = ",Pot%TolLimInf

       IF(ALLOCATED(Pot%coef))THEN
          PRINT*,"coef:"
          DO i=1,SIZE(Pot%coef),1
             PRINT*,Pot%coef(i)
          END DO
       ELSE
          PRINT*,"coef: not allocated"
       END IF

       IF(ALLOCATED(Pot%r))THEN
          PRINT*,"r:"
          DO i=1,SIZE(Pot%r),1
             PRINT*,Pot%r(i)
          END DO
       ELSE
          PRINT*,"r: not allocated"
       END IF

       IF(ALLOCATED(Pot%V))THEN
          PRINT*,"V:"
          DO i=1,SIZE(Pot%V),1
             PRINT*,Pot%V(i)
          END DO
       ELSE
          PRINT*,"V: not allocated"
       END IF

       IF(ALLOCATED(Pot%z))THEN
          PRINT*,"z:"
          DO i=1,SIZE(Pot%z),1
             PRINT*,Pot%z(i)
          END DO
       ELSE
          PRINT*,"z: not allocated"
       END IF

       PRINT*,"lambda = ",Pot%lambda

       PRINT*,"VARIABLES FOR DIRECT CALCULATION:"
       PRINT*,"ES_command = ",TRIM(Pot%ES_command)
       PRINT*,"rm_command = ",TRIM(Pot%rm_command)
       PRINT*,"format_ener = ",TRIM(Pot%format_ener)
       PRINT*,"format_r = ",TRIM(Pot%format_r)
       PRINT*,"keyword = ",TRIM(Pot%keyword)
       PRINT*,"get_ener_commant_c = ",TRIM(Pot%get_ener_command_c)
       PRINT*,"uni_EE_input = ",Pot%uni_EE_input
       PRINT*,"input:"
       current_string=>Pot%input
       DO WHILE(ASSOCIATED(current_string))
          PRINT*,TRIM(current_string%value)
          current_string=>current_string%next
       END DO

       PRINT*,"VARIABLES FOR GENERAL PROPERTIES:"
       PRINT*,"zero_energy_ret_p = ",Pot%zero_energy_ret_p
       PRINT*,"lower_bound = ",Pot%lower_bound
       PRINT*,"upper_bound = ",Pot%upper_bound
       PRINT*,"eq_dist = ",Pot%eq_dist
       PRINT*,"diss_ener = ",Pot%eq_dist
       PRINT*,"force_const = ",Pot%force_const

       PRINT*,"UNITS USED IN POTENTIAL DEFINITION:"
       PRINT*,"enegy_uni = ",TRIM(Pot%energy_uni)
       PRINT*,"dist_uni = ",TRIM(Pot%dist_uni)

       PRINT*,"----------------------------------------------------"

       Pot=>Pot%next
    END DO

  END SUBROUTINE Print_Potential_Data
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  SUBROUTINE Non_Consistent_Info(info,locale)
    
    IMPLICIT NONE
    
    CHARACTER(len = *), INTENT(in) :: info, locale
    
    error_msg = "Inconsistency of "//TRIM(info)//" in "//TRIM(locale)//"."
    CALL Error(2)
    
  END SUBROUTINE Non_Consistent_Info
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  SUBROUTINE Free_Potential_Memory()
    
    IMPLICIT NONE

    TYPE(TypePot_DipMomCurve), POINTER :: next=>NULL()

    Pot => Pot_list
    do while(associated(Pot))
       next => Pot%next

       IF(ALLOCATED(Pot%coef)) DEALLOCATE(Pot%coef)
       IF(ALLOCATED(Pot%r)) DEALLOCATE(Pot%r)
       IF(ALLOCATED(Pot%V)) DEALLOCATE(Pot%V)
       IF(ALLOCATED(Pot%z)) DEALLOCATE(Pot%z)
       if(ASSOCIATED(Total_pot)) CALL Free_RLL(Total_pot)
       IF(ASSOCIATED(Pot%input)) CALL Free_SLL(Pot%input)

       DEALLOCATE(Pot)

       Pot=>next

    END do

  END SUBROUTINE Free_Potential_Memory
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
END MODULE ModPot
