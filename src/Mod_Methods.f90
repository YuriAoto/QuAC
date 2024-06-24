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
! Module for job calculation
!
! This module get the data from a given knot of the jobs linked list, print an header, by Header subroutine, print the parameters
! to be used in the calculation, by Print_Parameter subroutine, and finally do the calculation, in Calculte subroutine, calling the
! necessary subprograms of ModCalc module. There is a subroutine for each kind of job. All transformation from the user units to
! atomic units, used in calculation, is carried here, along with the eventuality unit conversion of the final result
!
MODULE ModMethods

  USE ModUtil, ONLY : real_kind,status,at_uni_dist,at_uni_ener,at_uni_dip_mom,at_uni_mass,at_uni_vel,error_msg,uni_out,uni_log,&
       uni_inp,ios_var,ios_ok,quad_met,Unit_Converter,Error,pi,zero, one, two, Double_Factorial, Sph_Bessel_J, Sph_Bessel_N
  USE ModPot, ONLY : stop_error, Pot, Potential, formats, TypeIntegerLinkedList, TypeRealLinkedList, ILL_Size, RLL_Size, &
       Free_RLL, Classical_Turning_Points, Select_Potential, Total_pot, curr_pot_knot, save_Potential, use_Pot_vector
  USE ModInput, ONLY : Job_Data, TypeJobSchedule, Jobs_Data_list, PHASE_SHIFTS_kw, SCATT_LENGTH_kw, BOUND_STATES_kw, &
       RES_ANALYSIS_kw, MASSES_kw, RED_MASS_kw, ATOMS_kw, FORMATS_kw
  USE ModCalc, ONLY : TypeIntParameters, Integration, Calc_Cross_Sec, Calc_Dif_Cross_Sec,&
       Calc_Phase_Shift, Calc_Scatt_Length, Analytical_LenJonN_Nm2_Scatt_Length, Analytical_Square_Phase_Shift,&
       Analytical_HardSphere_Phase_Shift
  USE ModSystem, ONLY : particles, Sec_Deriv, Sec_Deriv_tilde, r_to_y_meshkov, lambda

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: Header, Print_Parameters, Calculate

  TYPE TypeVibrationalResults
     INTEGER :: J
     REAL(kind = real_kind) :: Energy
     INTEGER :: level
     INTEGER :: no_iter
     REAL(kind = real_kind), ALLOCATABLE, DIMENSION(:,:) :: wave_function
     REAL(kind = real_kind) :: eq_dist

     TYPE(TypeVibrationalResults), POINTER :: next=>NULL()
  END TYPE TypeVibrationalResults

  TYPE TypeVibLevelsTransitions
     TYPE(TypeVibrationalResults), POINTER :: state_1=>NULL()
     TYPE(TypeVibrationalResults), POINTER :: state_2=>NULL()
     REAL(kind = real_kind) :: DipMom
     REAL(kind = real_kind) :: FrankCondon_Fac
     REAL(kind = real_kind) :: Einstein_A_Coeff
     
     TYPE(TypeVibLevelsTransitions), POINTER :: next=>NULL()
  END TYPE TypeVibLevelsTransitions
  
  TYPE TypeResultsList
     CHARACTER(LEN = 100) :: name

     TYPE(TypeVibrationalResults), POINTER :: Vib=>NULL()
     TYPE(TypeVibLevelsTransitions), POINTER :: VibTrans=>NULL()

     REAL(kind = real_kind), ALLOCATABLE, DIMENSION(:,:) :: phase_shifts
     REAL(kind = real_kind), ALLOCATABLE, DIMENSION(:) :: Dunhan_parameters
     REAL(kind = real_kind), ALLOCATABLE, DIMENSION(:) :: wave_vector ! tamp

     REAL(kind = real_kind), ALLOCATABLE, DIMENSION(:,:) :: trans_dip_mom_matrix
     REAL(kind = real_kind), ALLOCATABLE, DIMENSION(:,:) :: cross_section

     TYPE(TypeResultsList), POINTER :: next=>NULL()
  END TYPE TypeResultsList

  TYPE(TypeResultsList), PROTECTED, SAVE, POINTER :: Results_list=>NULL()
  TYPE(TypeResultsList), POINTER :: Results=>NULL()

CONTAINS
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  ! Header of this job
  SUBROUTINE Header()

    IMPLICIT NONE

    WRITE(uni_out,FMT='("---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/---")')
    WRITE(uni_out,*)

    SELECT CASE(Job_Data%job_method)

    CASE(MASSES_kw, RED_MASS_kw, ATOMS_kw)
       
       particles => Job_Data%job_particles
       WRITE(uni_out,FMT='("Particles Definition")') 
       WRITE(uni_out,*)

    CASE(FORMATS_kw)

       formats => Job_Data%job_formats

    CASE(PHASE_SHIFTS_kw,SCATT_LENGTH_kw,BOUND_STATES_kw,RES_ANALYSIS_kw)

       WRITE(uni_log,FMT='("---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/---")')
       WRITE(uni_log,*)

       IF(.NOT.ASSOCIATED(Results))THEN
          ALLOCATE(Results_list)
          Results => Results_list
       ELSE
          ALLOCATE(Results%next)
          Results => Results%next
       END IF
       Results%name = Job_Data%job_name

       WRITE(uni_out,FMT='(A)') TRIM(Job_Data%complete_job_name)
       WRITE(uni_out,*)
       WRITE(uni_out,FMT='("Additional data printed in .log file.")')
       WRITE(uni_out,*)
       
       WRITE(uni_log,FMT='(A)') TRIM(Job_Data%complete_job_name)
       WRITE(uni_log,*)
       
       SELECT CASE(Job_Data%job_method)
          
       CASE(PHASE_SHIFTS_kw,SCATT_LENGTH_kw,BOUND_STATES_kw)
          
          CALL Select_Potential(Job_Data%job%pot_name)

          if(Job_Data%job_method.eq.BOUND_STATES_kw)then
             lambda = Pot%lambda
          else
             lambda = 0
          end if
          
       CASE(RES_ANALYSIS_kw)
       
       CASE DEFAULT
          error_msg = "Inconsistency of method in Header, level 1: "//TRIM(Job_Data%job_method)
          CALL Error(2)
          
       END SELECT

    CASE DEFAULT
       error_msg = "Inconsistency of method in Header, level 2: "//TRIM(Job_Data%job_method)
       CALL Error(2)

    END SELECT

  END SUBROUTINE Header
!!$!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  SUBROUTINE Print_Parameters()

    IMPLICIT NONE
    
    SELECT CASE(Job_Data%job_method)

    CASE(MASSES_kw,RED_MASS_kw,ATOMS_kw)
       
       IF(Job_Data%job_method.EQ.ATOMS_kw)THEN
          WRITE(uni_out,FMT='("Atomic data for ",A," and ",A,":")') TRIM(particles%name1),TRIM(particles%name2)
          WRITE(uni_out,FMT='("Atomic numbers: ",I0," and ",I0)') particles%Z1,particles%Z2
          WRITE(uni_out,FMT='("Isotopic abundances: ",F0.4," and ",F0.4)') particles%abundance1,particles%abundance2
       END IF
       
       IF(Particles%mass1.GT.zero)THEN
          WRITE(uni_out,FMT='("Atomic masses: ",F0.8," and ",F0.8,1X,A)') &
               Particles%mass1,Particles%mass2,TRIM(Particles%mass_uni)
       END IF
       
       WRITE(uni_out,FMT='("Reduced mass: ",F0.8,1X,A," = ")',advance='NO') &
            Particles%reduced_mass,TRIM(Particles%mass_uni)
       
       ! Print particle information before unit conversion
       IF(Particles%mass1.GT.zero) Particles%mass1 = Unit_Converter(Particles%mass_uni,at_uni_mass,Particles%mass1)
       IF(Particles%mass2.GT.zero) Particles%mass2 = Unit_Converter(Particles%mass_uni,at_uni_mass,Particles%mass2)
       Particles%reduced_mass = Unit_Converter(Particles%mass_uni,at_uni_mass,Particles%reduced_mass)
       
       WRITE(uni_out,FMT='(F0.8,1X,A)') Particles%reduced_mass,TRIM(at_uni_mass)

       WRITE(uni_out,*)
       
    CASE(FORMATS_kw)
       
       if(formats%default_values) WRITE(uni_out,FMT='("Formats restored to default values.")')
       
       if(formats%energy_def)        WRITE(uni_out,FMT='("Format of energy set to ",A,".")')        trim(formats%energy)
       if(formats%cross_section_def) WRITE(uni_out,FMT='("Format of cross_section set to ",A,".")') trim(formats%cross_section)
       if(formats%angle_def)         WRITE(uni_out,FMT='("Format of angle set to ",A,".")')         trim(formats%angle)
       if(formats%ang_momentum_def)  WRITE(uni_out,FMT='("Format of ang_momentum set to ",A,".")')  trim(formats%ang_momentum)
       if(formats%distance_def)      WRITE(uni_out,FMT='("Format of distance set to ",A,".")')      trim(formats%distance)
       if(formats%wave_function_def) WRITE(uni_out,FMT='("Format of wave_function set to ",A,".")') trim(formats%wave_function)
       if(formats%potential_def)     WRITE(uni_out,FMT='("Format of potential set to ",A,".")')     trim(formats%potential)
       if(formats%energy_cm1_def)    WRITE(uni_out,FMT='("Format of energy_cm1 set to ",A,".")')    trim(formats%energy_cm1)
       if(formats%einstein_coef_def) WRITE(uni_out,FMT='("Format of Einstein_coef set to ",A,".")') trim(formats%einstein_coef)
       if(formats%dip_mom_def)       WRITE(uni_out,FMT='("Format of dip_mom set to ",A,".")')       trim(formats%dip_mom)
       if(formats%phase_shift_def)   WRITE(uni_out,FMT='("Format of phase_shift set to ",A,".")')   trim(formats%phase_shift)
       if(formats%dunhan_def)        WRITE(uni_out,FMT='("Format of Dunhan set to ",A,".")')        trim(formats%dunhan)
       if(formats%vib_level_def)     WRITE(uni_out,FMT='("Format of vib_level set to ",A,".")')     trim(formats%vib_level)
       
       WRITE(uni_out,*)
       
    CASE(PHASE_SHIFTS_kw,SCATT_LENGTH_kw,BOUND_STATES_kw)

       WRITE(uni_out,FMT='("-- Parameters used in this job --")')
       WRITE(uni_out,*)
       
       ! Integration method
       WRITE(uni_out,FMT='("Integration method: ",A)',Advance='NO') TRIM(Job_Data%job%int_method)
       CALL User_Default(Job_Data%job%int_method_def)
       
       ! Save steps
       IF(Job_Data%job%save_steps)THEN
          WRITE(uni_out,FMT='("Points saved from numerical integration writen in log file.")',Advance='NO')
          CALL User_Default(Job_Data%job%save_steps_def)
       END IF
       
       ! Integration step size
       WRITE(uni_out,FMT='("Integration step size = ",'//trim(formats%distance)//')',Advance='NO') Job_Data%job%int_step
       IF(.NOT.ASSOCIATED(Pot))THEN
          WRITE(uni_out,FMT='(" (dimension less unit)")',Advance='NO')
       ELSE
          WRITE(uni_out,FMT='(1X,A)',Advance='NO') TRIM(Pot%dist_uni)
          IF(.NOT.Job_Data%job%int_step_def.AND.Pot%dist_uni.NE.at_uni_dist) THEN
             WRITE(uni_out,FMT='("= ",'//trim(formats%distance)//',1X,A)',Advance='NO') &
                  Unit_Converter(Pot%dist_uni,at_uni_dist,Job_Data%job%int_step),TRIM(at_uni_dist)
          END IF
       END IF
       CALL User_Default(Job_Data%job%int_step_def)

       IF(Job_Data%job_method.EQ.PHASE_SHIFTS_kw.OR.Job_Data%job_method.EQ.SCATT_LENGTH_kw)THEN

          ! Final integration step size
          WRITE(uni_out,FMT='("Final integration step size = ",'//trim(formats%distance)//')',Advance='NO') &
               Job_Data%job%final_int_step
          IF(.NOT.ASSOCIATED(Pot))THEN
             WRITE(uni_out,FMT='(" (dimension less unit)")',Advance='NO')
          ELSE
             WRITE(uni_out,FMT='(1X,A)',Advance='NO') TRIM(Pot%dist_uni)
             IF(.NOT.Job_Data%job%final_int_step_def.AND.Pot%dist_uni.NE.at_uni_dist) THEN
                WRITE(uni_out,FMT='("= ",'//trim(formats%distance)//',1X,A)',Advance='NO') &
                     &Unit_Converter(Pot%dist_uni,at_uni_dist,Job_Data%job%final_int_step),&
                     &TRIM(at_uni_dist)
             END IF
          END IF
          CALL User_Default(Job_Data%job%final_int_step_def)

       END IF
       
       ! Print
       IF(Job_Data%job%print_at.GT.0)THEN
          WRITE(uni_out,FMT='("Integration printed in log file at each ",I3," steps.")',Advance='NO') Job_Data%job%print_at
          CALL User_Default(Job_Data%job%print_at_def)
       END IF

       IF(Job_Data%job%print_vib_iterations)THEN
          WRITE(uni_out,FMT='("Integration also printed in log file for each vibrational iteration.")',Advance='NO') 
          CALL User_Default(Job_Data%job%print_vib_iterations_def)
       END IF
       
       IF(Job_Data%job_method.EQ.PHASE_SHIFTS_kw.OR.Job_Data%job_method.EQ.SCATT_LENGTH_kw)THEN
          
          ! r_min
          WRITE(uni_out,FMT='("Integration starting at ",'//trim(formats%distance)//')',Advance='NO') Job_Data%job%r_min
          IF(ASSOCIATED(Pot))THEN
             WRITE(uni_out,FMT='(" (dimension less unit)")',Advance='NO')
          ELSE
             WRITE(uni_out,FMT='(1X,A)',Advance='NO') TRIM(Pot%dist_uni)
          END IF
          CALL User_Default(Job_Data%job%r_min_def,Job_Data%job%r_min_default_explanation)
          
          ! r_max
          WRITE(uni_out,FMT='("Integration finishing at ",'//trim(formats%distance)//')',Advance='NO') Job_Data%job%r_max
          IF(.NOT.ASSOCIATED(Pot))THEN
             WRITE(uni_out,FMT='("(dimension less unit)")',Advance='NO')
          ELSE
             WRITE(uni_out,FMT='(1X,A)',Advance='NO') TRIM(Pot%dist_uni)
          END IF
          CALL User_Default(Job_Data%job%r_max_def,Job_Data%job%r_max_default_explanation)

       ELSE

          IF(Job_Data%job%vib_fixed_mesh_points)then
             WRITE(uni_out,FMT='("Vibrational integration with fixed mesh points.")')
             WRITE(uni_out,FMT='("Mesh points between ",'//trim(formats%distance)//'," and ",'//trim(formats%distance)//',1X,A)') &
                  Job_Data%job%r_min,Job_Data%job%r_max,TRIM(Pot%dist_uni)
             CALL User_Default(Job_Data%job%r_max_def,Job_Data%job%r_max_default_explanation)
          else
          end IF
          
          ! Rmin and rmax information for BOUND_STATES
   
       END IF

       ! Wave function scale
       WRITE(uni_out,FMT='("Wave function scale: ",ES10.3)',Advance='NO') Job_Data%job%wf_scale
       CALL User_Default(Job_Data%job%wf_scale_def)
       
       ! Initial condition
       IF(Job_Data%job%ini_cond_free_part)THEN
          WRITE(uni_out,FMT='("Initial condition of free particle. ")',Advance='NO')
       ELSE
          WRITE(uni_out,FMT='("Initial null function and derivative of free particle. ")',Advance='NO')
       END IF
       CALL User_Default(Job_Data%job%ini_cond_def)

       IF(Job_Data%job_method.EQ.PHASE_SHIFTS_kw.OR.Job_Data%job_method.EQ.SCATT_LENGTH_kw)THEN

          ! Phase Shift method
          WRITE(uni_out,FMT='("Phase shift calculation method: ",A)',Advance='NO') TRIM(Job_Data%job%phase_shift_method)
          CALL User_Default(Job_Data%job%phase_shift_method_def)

          ! Additional steps
          WRITE(uni_out,FMT='("Integrating ",I5," extra steps to get wave function information.")',Advance='NO') &
               Job_Data%job%add_steps
          CALL User_Default(Job_Data%job%add_steps_def)

          ! Save at each
          WRITE(uni_out,FMT='("Integration points saved at each ",I3," steps.")',Advance='NO') Job_Data%job%save_at
          CALL User_Default(Job_Data%job%save_at_def)

       END IF

       WRITE(uni_out,*)

    CASE(RES_ANALYSIS_kw)

       WRITE(uni_out,FMT='("-- Parameters used in this analysis --")')
       WRITE(uni_out,*)

       WRITE(uni_out,*)

    CASE DEFAULT
       
       error_msg = "Inconsistency of method in Print_Parameters: "//TRIM(Job_Data%job_method)
       CALL Error(2)
       
    END SELECT

  CONTAINS
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
    SUBROUTINE User_Default(user_defined,explanation)
      
      IMPLICIT NONE
      
      LOGICAL, INTENT(in) :: user_defined
      CHARACTER(len = *), OPTIONAL :: explanation

      IF(user_defined)THEN
         WRITE(uni_out,FMT='(" (user defined).")')
      ELSE
         WRITE(uni_out,FMT='(" (default")',advance='no')
         IF(PRESENT(explanation).AND.LEN_TRIM(explanation).GT.0)THEN
            WRITE(uni_out,FMT='(" - ",A,")")') TRIM(explanation)
         ELSE
            WRITE(uni_out,FMT='(")")')
         END IF
      END IF
      
    END SUBROUTINE User_Default
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  END SUBROUTINE Print_Parameters
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  SUBROUTINE Calculate()
    
    IMPLICIT NONE
    
    SELECT CASE(Job_Data%job_method)

    CASE(PHASE_SHIFTS_kw)
       CALL PHASE_SHIFTS()

    CASE(SCATT_LENGTH_kw)
       CALL SCATT_LENGTH()

    CASE(BOUND_STATES_kw)
       CALL BOUND_STATES()

    CASE(RES_ANALYSIS_kw)
       CALL RES_ANALYSIS()

    CASE(MASSES_kw,RED_MASS_kw,ATOMS_kw,FORMATS_kw)

    CASE DEFAULT

       error_msg = "Inconsistency of method in Calculate: "//TRIM(Job_Data%job_method)
       CALL Error(2)

    END SELECT

    if(associated(Total_pot)) CALL Free_RLL(Total_pot)

    SELECT CASE(Job_Data%job_method)
    CASE(PHASE_SHIFTS_kw,SCATT_LENGTH_kw,BOUND_STATES_kw,RES_ANALYSIS_kw)
       
       WRITE(uni_out,FMT='("END OF ",A," CALCULATION")') TRIM(Job_Data%job_method)
       WRITE(uni_out,*)

    CASE(MASSES_kw,RED_MASS_kw,ATOMS_kw,FORMATS_kw)

    CASE DEFAULT
       error_msg = "Inconsistency of method in Calculate: "//TRIM(Job_Data%job_method)
       CALL Error(2)

    END SELECT


  END SUBROUTINE Calculate
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  SUBROUTINE PHASE_SHIFTS()
    
    IMPLICIT NONE

    TYPE(TypeIntParameters) :: arg
    REAL(kind = real_kind) :: energy, wave_vector
    REAL(kind = real_kind), ALLOCATABLE, DIMENSION(:,:) :: points
    LOGICAL :: is_first_int
    TYPE(TypeRealLinkedList), POINTER :: current_real => NULL()
    TYPE(TypeIntegerLinkedList), POINTER :: current_int => NULL()
    INTEGER :: i_l, i_E

    ALLOCATE(Results%phase_shifts(RLL_Size(Job_Data%job%energy),ILL_Size(Job_Data%job%l)))
    ALLOCATE(Results%wave_vector(RLL_Size(Job_Data%job%energy)))
    
    is_first_int = .TRUE.

    current_real=>Job_Data%job%energy
    DO i_E = 1,SIZE(Results%phase_shifts,1),1
       
       energy = Set_Energy(current_real%value,wave_vector)
       Results%wave_vector(i_E) = wave_vector
       
       CALL Print_Current_Calculation(Job_Data%job%save_steps.OR.Job_Data%job%print_at.GT.0,energy=energy)
       
       current_int=>Job_Data%job%l
       DO i_l = 1,SIZE(Results%phase_shifts,2),1
          
          CALL Set_Integration_Parameters(arg,current_int%value,is_first_int,wave_vector=wave_vector,&
               save_global_pot=is_first_int,use_global_pot=(.NOT.is_first_int))

          CALL Print_Current_Calculation(Job_Data%job%save_steps.OR.Job_Data%job%print_at.GT.0,l=current_int%value)

          CALL Integration(arg,current_int%value,energy,wave_vector,points)

          is_first_int = .FALSE.
          
          Results%phase_shifts(i_E,i_l) = calc_phase_shift(current_int%value,points,wave_vector) !energy,wave_vector)
          
          DEALLOCATE(points)
          
          current_int=>current_int%next
       END DO
       current_real=>current_real%next
    END DO
    
    ! Print the final result
    IF(Job_Data%job%print_ps) then
       
       WRITE(uni_log,FMT='("+----------------------------------------------+")')
       WRITE(uni_log,FMT='("|                  Phase Shifts                |")')
       WRITE(uni_log,*)
       
       current_real=>Job_Data%job%energy 
       DO i_E = 1,SIZE(Results%phase_shifts,1),1
          
          WRITE(uni_log,FMT='(14X,"Energy = ",'//TRIM(formats%energy)//',1X,A)') &
               current_real%value,TRIM(Job_Data%job%energy_uni)
          WRITE(uni_log,FMT='(21X,"l    phase shift")')
          
          current_int=>Job_Data%job%l
          DO i_l = 1,SIZE(Results%phase_shifts,2),1
             WRITE(uni_log,FMT='(21X,'//TRIM(formats%ang_momentum)//',"  ",'//trim(formats%phase_shift)//')') &
                  current_int%value,Results%phase_shifts(i_E,i_l)
             current_int=>current_int%next
          END DO
          
          WRITE(uni_log,*)
          current_real=>current_real%next
       END DO
       WRITE(uni_log,FMT='("+----------------------------------------------+")')
       WRITE(uni_log,*)
       
       WRITE(uni_out,FMT='("Phase shifts writen in log file.")')
       WRITE(uni_out,*)
       
    END IF

  END SUBROUTINE PHASE_SHIFTS
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  SUBROUTINE SCATT_LENGTH()
    
    IMPLICIT NONE

    TYPE(TypeIntParameters) :: arg
    REAL(kind = real_kind) :: energy, wave_vector
    INTEGER, PARAMETER :: l_zero = 0
    LOGICAL :: is_first_int
    REAL(kind = real_kind), ALLOCATABLE, DIMENSION(:,:) :: points
    REAL(kind = real_kind) :: scatt_length_result, scatt_length_result2, gamma, Rich_correction

    ALLOCATE(Results%phase_shifts(1,1))
    
    is_first_int = .TRUE.

    IF(ASSOCIATED(Job_Data%job%energy))THEN
       energy = Set_Energy(Job_Data%job%energy%value,wave_vector)
       CALL Print_Current_Calculation(Job_Data%job%save_steps.OR.Job_Data%job%print_at.GT.0,energy=energy)
       CALL Set_Integration_Parameters(arg,l_zero,is_first_int,wave_vector)
    ELSE
       CALL Set_Integration_Parameters(arg,l_zero,is_first_int)
    END IF
    
    IF(ASSOCIATED(Job_Data%job%energy))THEN
       CALL Integration(arg,l_zero,energy,wave_vector,points)
    ELSE
       CALL Integration(arg,l_zero,points=points)
    END IF
    
    is_first_int = .FALSE.

    IF(ASSOCIATED(Job_Data%job%energy))THEN
       scatt_length_result = Calc_Scatt_Length(points,Job_Data%job%scatt_len_method,wave_vector)
    ELSEIF(Job_Data%job%scatt_len_method.NE."MESHKOV")THEN
       scatt_length_result = Calc_Scatt_Length(points,Job_Data%job%scatt_len_method)
    ELSE
       scatt_length_result = Calc_Scatt_Length(points,Job_Data%job%scatt_len_method,r_bar=arg%r_bar_mesh,beta=arg%beta_mesh)
    END IF
    
    DEALLOCATE(points)
    
    Results%phase_shifts(1,1) = scatt_length_result
    
    IF(Job_Data%job%Rich_extr)THEN
       
       IF(Job_Data%job%scatt_len_method.NE."MESHKOV")THEN
          error_msg = "Attempt to calculate Richardson extrapolation by "//TRIM(Job_Data%job%scatt_len_method)
          CALL Error(2)
       END IF
       
       ! add_steps must be multiple of 4.
       Job_Data%job%add_steps = Job_Data%job%add_steps/2
       
       CALL Set_Integration_Parameters(arg,l_zero,is_first_int)
       
       IF(ASSOCIATED(Job_Data%job%energy))THEN
          CALL Integration(arg,0,energy,wave_vector,points)
       ELSE
          CALL Integration(arg,0,points=points)
       END IF
       
       IF(ASSOCIATED(Job_Data%job%energy))THEN
          error_msg = "Attempt to calculate Richardson extrapolation by a non-zero energy method."
          CALL Error(2)             
       ELSE
          scatt_length_result2 = Calc_Scatt_Length(points,Job_Data%job%scatt_len_method,r_bar=arg%r_bar_mesh,beta=arg%beta_mesh)
       END IF
       DEALLOCATE(points)
       
       gamma = two
       Rich_correction = (scatt_length_result - scatt_length_result2)/(gamma**4 - 1)
       
       Results%phase_shifts(1,1) = Results%phase_shifts(1,1) + Rich_correction
       
       Rich_correction = Unit_Converter(at_uni_dist,Job_Data%job%cross_sec_uni,Rich_correction)
       
    END IF
    
    scatt_length_result = Unit_Converter(at_uni_dist,Job_Data%job%cross_sec_uni,scatt_length_result)

    ! Impressão do resultado final no output
    WRITE(uni_out,FMT='("+----------------------------------------------+")')
    WRITE(uni_out,FMT='("|              Scattering Length               |")')
    WRITE(uni_out,*)
    
    WRITE(uni_out,FMT='("a_scatt/",A,"           = ",'//TRIM(formats%cross_section)//')')&
         &TRIM(Job_Data%job%cross_sec_uni),scatt_length_result
    IF(Job_Data%job%Rich_extr)THEN
       WRITE(uni_out,FMT='("a_scatt(Rich Ext)/",A," = ",'//TRIM(formats%cross_section)//')')&
            &TRIM(Job_Data%job%cross_sec_uni),scatt_length_result+Rich_correction
    END IF
    WRITE(uni_out,FMT='("Cross section(k=0)/",A,"² = 4 pi a_scatt²            = ",'//TRIM(formats%cross_section)//')')&
         TRIM(Job_Data%job%cross_sec_uni),4*pi*scatt_length_result**2
    IF(Job_Data%job%Rich_extr)THEN
       WRITE(uni_out,FMT='("Cross section(k=0; Rich Extr)/",A,"² = 4 pi a_scatt² = ",'//TRIM(formats%cross_section)//')')&
            TRIM(Job_Data%job%cross_sec_uni),4*pi*(scatt_length_result+Rich_correction)**2
    END IF
    
    WRITE(uni_out,*)
    
    WRITE(uni_out,FMT='("+----------------------------------------------+")')
    WRITE(uni_out,*)
    
    ! Analytical solution !!
    
  END SUBROUTINE SCATT_LENGTH
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  SUBROUTINE BOUND_STATES()

    IMPLICIT NONE

    TYPE(TypeVibrationalResults), POINTER :: current_calculation => NULL()
    TYPE(TypeIntegerLinkedList), POINTER :: current_l => NULL()
    TYPE(TypeIntegerLinkedList), POINTER :: current_level => NULL()
    TYPE(TypeRealLinkedList), POINTER :: current_real => NULL()
    REAL(kind = real_kind) :: q_new, q_sup, q_inf, harmonic_freq, dif_energy, energy_out, wave_vector
    REAL(kind = real_kind) :: vib_threshold_hartree
    INTEGER :: i_iter, level, no_knots
    LOGICAL :: levels_def, first_level, final_level
    LOGICAL :: very_first_int
    
    very_first_int=.true.

    if(ASSOCIATED(Job_Data%job%energy))then
       vib_threshold_hartree = Unit_Converter(Job_Data%job%energy_uni,at_uni_ener,Job_Data%job%vib_threshold)
    else
       vib_threshold_hartree = Unit_Converter("cm-1           ",at_uni_ener,Job_Data%job%vib_threshold)
    end if

    IF(ASSOCIATED(Job_Data%job%energy))THEN
       
       current_real=>Job_Data%job%energy

       DO WHILE(ASSOCIATED(current_real))

          current_l=>Job_Data%job%l
          DO WHILE(ASSOCIATED(current_l))

             i_iter = 0
             q_new = Set_Energy(current_real%value, wave_vector)

             IF(.NOT.ASSOCIATED(current_calculation))THEN
                ALLOCATE(Results%vib)
                current_calculation=>Results%vib
             ELSE
                ALLOCATE(current_calculation%next)
                current_calculation => current_calculation%next
             END IF
             
             current_calculation%J = current_l%value

             CALL Sec_Order_It(current_l%value,q_new,Job_Data%job%vib_max_iter,vib_threshold_hartree,energy_out,i_iter,&
                  current_calculation%wave_function)
             if(status.eq.0)then

                current_calculation%Energy = Pot%diss_ener+energy_out
                current_calculation%level = no_knots
                current_calculation%no_iter = i_iter
             
             end if

             current_l=>current_l%next

          END DO
          current_real=>current_real%next
       END DO
       
    ELSE
       if(Pot%force_const.eq.zero)then
          harmonic_freq = pi/Pot%radius
       else
          harmonic_freq = SQRT(Pot%force_const/Particles%reduced_mass)
       end if

       current_l=>Job_Data%job%l
       DO WHILE(ASSOCIATED(current_l))
          levels_def = .FALSE.
          level = -1
          first_level = .TRUE.
          IF(ASSOCIATED(Job_Data%job%vib_lev))THEN
             levels_def = .TRUE.
             current_level=>Job_Data%job%vib_lev
          END IF
          
          DO !Over vibrational levels

             i_iter = 0

             IF(levels_def)THEN
                IF(.NOT.ASSOCIATED(current_level)) EXIT
                level = current_level%value
             ELSE
                level = level+1
             END IF

             If(Job_Data%job%print_vib_iterations.OR.Job_Data%job%print_at.GT.0) WRITE(uni_log,FMT='("-------------")')
             If(Job_Data%job%print_vib_iterations.OR.Job_Data%job%print_at.GT.0) &
                  WRITE(uni_log,FMT='("Vibrational level: ",'//trim(formats%vib_level)//')') level
             IF(first_level)THEN
                if(Pot%force_const.eq.zero)then
                   q_inf = -Pot%diss_ener - harmonic_freq/5
                else
                   q_inf = -Pot%diss_ener + harmonic_freq/5
                end if
                q_sup = -Pot%diss_ener + harmonic_freq
             ELSE
                q_inf = energy_out
                q_sup = q_new + harmonic_freq
                IF(q_sup.GT.zero)THEN
                   q_sup = energy_out/2
                END IF
             END IF

             ! Determine q_inf
             i_iter = i_iter+1
             If(Job_Data%job%print_vib_iterations)then
                WRITE(uni_log,FMT='("Iteration number (q_inf): ",I3)') i_iter
                CALL Print_energy_cm_and_ua(q_inf,relative_energy=.true.)
             end If

             CALL VibrationalIntegration(current_l%value,q_inf,dif_energy)
             if(status.ne.0) exit

             IF(.NOT.is_lower(no_knots,level,dif_energy))THEN
                error_msg = "Q_inf is not lower than true energy."
                CALL Error(2)
             END IF
             If(Job_Data%job%print_vib_iterations) WRITE(uni_log,FMT='("---")')

             ! Determine q_sup
             i_iter = i_iter+1
             If(Job_Data%job%print_vib_iterations)then
                WRITE(uni_log,FMT='("Iteration number (q_sup): ",I3)') i_iter
                CALL Print_energy_cm_and_ua(q_sup,relative_energy=.true.)
             end If

             CALL VibrationalIntegration(current_l%value,q_sup,dif_energy)

             if(status.ne.0) exit
             If(Job_Data%job%print_vib_iterations) WRITE(uni_log,FMT='("---")')

             final_level = .FALSE.
             DO WHILE(is_lower(no_knots,level,dif_energy))
                IF(ABS(q_sup).LT.1.0E-10) THEN
                   final_level = .TRUE.
                   EXIT
                END IF
                q_sup = q_sup/2

                i_iter=i_iter+1

                If(Job_Data%job%print_vib_iterations)then
                   WRITE(uni_log,FMT='("Iteration number (q_sup): ",I3)') i_iter
                   CALL Print_energy_cm_and_ua(q_sup,relative_energy=.true.)
                end If

                CALL VibrationalIntegration(current_l%value,q_sup,dif_energy)
                if(status.ne.0) exit 
                If(Job_Data%job%print_vib_iterations) WRITE(uni_log,FMT='("---")')

             END DO
             if(status.ne.0) exit

             IF(final_level) EXIT

             dif_energy = one
             DO WHILE(i_iter.LE.Job_Data%job%vib_1st_order_max_iter) ! First order knot counter algorithm

                q_new = (q_sup+q_inf)/2

                i_iter = i_iter+1
                If(Job_Data%job%print_vib_iterations)then
                   WRITE(uni_log,FMT='("Iteration number: ",I3)') i_iter
                   WRITE(uni_log,FMT='("Energy inf: ")',ADVANCE='NO')
                   CALL Print_energy_cm_and_ua(q_inf,relative_energy=.true.)
                   WRITE(uni_log,FMT='("Energy med: ")',ADVANCE='NO')
                   CALL Print_energy_cm_and_ua(q_new,relative_energy=.true.)
                   WRITE(uni_log,FMT='("Energy sup: ")',ADVANCE='NO')
                   CALL Print_energy_cm_and_ua(q_sup,relative_energy=.true.)
                end If

                CALL VibrationalIntegration(current_l%value,q_new,dif_energy)
                if(status.ne.0) exit

                If(Job_Data%job%print_vib_iterations)then
                   WRITE(uni_log,FMT='("no_knots = ",'//trim(formats%vib_level)//')') no_knots
                   WRITE(uni_log,FMT='("Energy diff: ")',ADVANCE='NO')
                   CALL Print_energy_cm_and_ua(dif_energy)
                   WRITE(uni_log,FMT='("Energy med + diff: ")',ADVANCE='NO')
                   CALL Print_energy_cm_and_ua(q_new+dif_energy,relative_energy=.true.)
                   WRITE(uni_log,FMT='("---")')
                end If

                IF(no_knots.EQ.level.AND.(q_new + dif_energy.LT.q_sup).AND.(q_new + dif_energy.GT.q_inf))THEN
                   EXIT
                ELSE
                   IF(is_lower(no_knots,level,dif_energy))THEN
                      q_inf = q_new
                   ELSE
                      q_sup = q_new
                   END IF
                END IF
                
             END DO
             if(status.ne.0) exit

             IF(i_iter.GE.Job_Data%job%vib_1st_order_max_iter)THEN
                WRITE(uni_log,FMT='("==> Maximum number of iterations (",I3,") reached in fist order knot counter algorithm")') &
                     i_iter
             END IF
             
             If(Job_Data%job%print_vib_iterations)then
                WRITE(uni_log,*)
                WRITE(uni_log,FMT='("Starting second order algorithm")')
                WRITE(uni_log,FMT='("---")')
             end If

             q_new = q_new + dif_energy

             IF(.NOT.ASSOCIATED(current_calculation))THEN
                ALLOCATE(Results%vib)
                current_calculation=>Results%vib
             ELSE
                ALLOCATE(current_calculation%next)
                current_calculation => current_calculation%next
             END IF

             current_calculation%J = current_l%value

             CALL Sec_Order_It(current_l%value,q_new,Job_Data%job%vib_max_iter,Job_Data%job%vib_threshold,energy_out,i_iter,&
                  current_calculation%wave_function)
             if(status.ne.0) exit
             
             ! Energy relative to the bottom of the well
             current_calculation%Energy = Pot%diss_ener+energy_out
             current_calculation%level = no_knots
             current_calculation%no_iter = i_iter

             first_level = .FALSE.
             IF(levels_def) current_level=>current_level%next
          END DO
          if(status.ne.0) exit
          current_l=>current_l%next
       END DO
    END IF

    WRITE(uni_out,FMT='("+-------------------------------------------------------+")')
    WRITE(uni_out,FMT='("|                  Vibrational Levels                   |")')

    WRITE(uni_out,FMT='(2X,"J",1X,"level  Energy/",A,"   Rel_Energy/cm-1  no_iter" )') TRIM(at_uni_ener)

    current_calculation=>Results%vib
    DO WHILE(ASSOCIATED(current_calculation))

       WRITE(uni_out,FMT='(1X,'//trim(formats%ang_momentum)//',2X,'//trim(formats%vib_level)//&
            ',3X,'//TRIM(formats%energy)//',5X,'//TRIM(formats%energy_cm1)//',5X,I3)') &
            current_calculation%J,&
            current_calculation%level,&
            current_calculation%energy-Pot%diss_ener,&
            Unit_Converter(at_uni_ener,"cm-1           ",current_calculation%energy),&
            current_calculation%no_iter

       current_calculation=>current_calculation%next
    END DO
    WRITE(uni_out,FMT='("+-------------------------------------------------------+")')

  CONTAINS
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
    SUBROUTINE Sec_Order_It(l,energy_in,max_iter_SO,eps_conv,energy_out,i_iter,points_final)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(in) :: l
      REAL(kind = real_kind), INTENT(in) :: energy_in
      INTEGER, INTENT(in) :: max_iter_SO
      REAL(kind = real_kind), INTENT(in) :: eps_conv
      REAL(kind = real_kind), INTENT(out) :: energy_out
      INTEGER, INTENT(inout) :: i_iter
      REAL(kind = real_kind), ALLOCATABLE, DIMENSION(:,:), INTENT(out), OPTIONAL :: points_final

      REAL(kind = real_kind) :: dif_ener, pot_print
      Integer :: J, i
      
      energy_out = energy_in
      dif_ener = 2*eps_conv
      
      DO WHILE(i_iter.LE.max_iter_SO.AND.ABS(dif_ener).GT.eps_conv)
         i_iter = i_iter+1

         If(Job_Data%job%print_vib_iterations)then
            WRITE(uni_log,FMT='("Iteration number: ",I3)') i_iter         
            WRITE(uni_log,FMT='("Energy: ")',ADVANCE='NO')
            CALL Print_energy_cm_and_ua(energy_out,relative_energy=.true.)
         end If

         CALL VibrationalIntegration(l,energy_out,dif_ener)

         if(status.ne.0) return

         If(Job_Data%job%print_vib_iterations)then
            WRITE(uni_log,FMT='("Diff energy: ")',ADVANCE='NO')
            CALL Print_energy_cm_and_ua(dif_ener)
            WRITE(uni_log,FMT='("New energy: ")',ADVANCE='NO')
            CALL Print_energy_cm_and_ua(energy_out+dif_ener,relative_energy=.true.)
         end If
                
         energy_out = energy_out+dif_ener
      END DO

      CALL VibrationalIntegration(l,energy_out,dif_ener,points_final)

      if(status.ne.0) return

      if(Job_Data%job%print_at.GT.0)then

         Write(uni_log,*)
         Write(uni_log,'("  Vibrational wave function:")')
         Write(uni_log,'("  distance (a0)   radial function   radial potential  effective potential")')
         J = l*(l+1)
         do i=1,size(points_final,1),Job_Data%job%print_at
            
            pot_print = Potential(points_final(i,1))

            Write(uni_log,FMT='('//TRIM(formats%distance)//'," ",'//TRIM(formats%wave_function)//'," ",'&
                 //TRIM(formats%potential)//'," ",'//TRIM(formats%potential)//')')&
                 points_final(i,1),points_final(i,2),pot_print,pot_print+J/points_final(i,1)**2
         end do

      end if

    END SUBROUTINE Sec_Order_It
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
    SUBROUTINE VibrationalIntegration(l,energy,delta_energy,points_total)
      
      IMPLICIT NONE    
      
      INTEGER, INTENT(in) :: l
      REAL(kind = real_kind), INTENT(in) :: energy
      REAL(kind = real_kind), INTENT(out) :: delta_energy
      REAL(kind = real_kind), ALLOCATABLE, DIMENSION(:,:), INTENT(out), OPTIONAL :: points_total
      
      TYPE(TypeIntParameters) :: arg
      REAL(kind = real_kind), ALLOCATABLE, DIMENSION(:,:) :: points
      REAL(kind = real_kind) :: A_1, A_2, B_1, B_2
      REAL(kind = real_kind) :: wave_vector
      
      REAL(kind = real_kind) :: deriv_from_left, deriv_from_right
      REAL(kind = real_kind) :: beg_in, beg_out
      REAL(kind = real_kind) :: int_of_sq, int_of_sq_new
      
      LOGICAL :: save_pot, use_pot

      LOGICAL, SAVE :: step_size_def=.FALSE.
      REAL(real_kind), SAVE :: step_size_at_uni, pass_by_point

      INTEGER :: J, n_middle_from_min, n_steps, no_knots_new, i

      wave_vector = SIGN(SQRT(ABS(2*Particles%reduced_mass*energy)),energy)

      J = l*(l+1)

      stop_error = .FALSE.

      CALL Get_Initial_Points(energy,beg_in,beg_out)

      if(status.ne.0) return
      stop_error = .TRUE.

      IF(.NOT.step_size_def)THEN
         n_steps = NINT((beg_out - beg_in)/Unit_converter(pot%dist_uni,at_uni_dist,Job_Data%job%int_step))
         step_size_at_uni = (beg_out - beg_in)/n_steps
         pass_by_point = Pot%eq_dist
         step_size_def = .TRUE.
      END IF

      beg_in = pass_by_point + NINT((beg_in - pass_by_point)/step_size_at_uni)*step_size_at_uni
      if(beg_in.lt.Pot%lower_bound)then
         beg_in = beg_in+step_size_at_uni
      end if
      beg_out = pass_by_point + NINT((beg_out - pass_by_point)/step_size_at_uni)*step_size_at_uni

      n_steps = NINT((beg_out - beg_in)/step_size_at_uni)
!      if(MOD(n_steps,2).EQ.1)then ! Even number of intervals for Simpson's rule
!         beg_out = beg_out - step_size_at_uni
!         n_steps = NINT((beg_out - beg_in)/step_size_at_uni)
!      end if

      n_middle_from_min = NINT((Pot%eq_dist - beg_in)/step_size_at_uni)
!      if(MOD(n_middle_from_min,2).EQ.1) n_middle_from_min = n_middle_from_min+1 ! Even number of intervals for Simpson's rule

      IF(PRESENT(points_total)) ALLOCATE(points_total(n_steps+1,2))

      if(Job_Data%job%use_potential_vector)then
         if(very_first_int)then
            save_pot=.true.
            use_pot=.false.
         else
            save_pot=.false.
            use_pot=.true.
         end if
         
      else
         save_pot=.false.
         use_pot=.false.
      end if

     ! Outward integration
      CALL Set_Integration_Parameters(arg,l,.true.,wave_vector,vib_step_size=step_size_at_uni,beg_vib_int=beg_in&
           ,vib_n_steps=n_middle_from_min,save_vib_wf=PRESENT(points_total)&
           ,save_global_pot=save_pot,use_global_pot=use_pot)

      IF(arg%print.GT.0) WRITE(uni_log,FMT='("Outward integration:")') 

      CALL Integration(arg,l,energy,wave_vector,points,Int_of_sq_new,no_knots)

      IF(PRESENT(points_total))THEN
         DO i=1,n_middle_from_min+1
            points_total(i,:) = points(i,:)
            points_total(i,1) = points_total(i,1)/wave_vector
         END DO
         
      ELSE
         IF(arg%save_last_deriv)THEN
            deriv_from_left = wave_vector*points(2,2)/points(1,2)
         ELSE
            A_1 = (points(4,2) - points(2,2))/2
            A_2 = (points(5,2) - points(1,2))/2
            
            B_1 = (arg%step_size**2)*(Sec_Deriv(J,points(4,1),wave_vector,energy)*points(4,2) -&
                 Sec_Deriv(J,points(2,1),wave_vector,energy)*points(2,2))/12
            B_2 = (arg%step_size**2)*(Sec_Deriv(J,points(5,1),wave_vector,energy)*points(5,2) -&
                 Sec_Deriv(J,points(1,1),wave_vector,energy)*points(1,2))/12
            
            deriv_from_left = 16*(-A_1 + (37*A_2)/32 - (37*B_1)/5 - (17*B_2)/40)/(21*ABS(arg%step_size/wave_vector))
            deriv_from_left = deriv_from_left/points(3,2)
         END IF

      END IF
      
      IF(arg%save_last_deriv)THEN
         Int_of_sq = Int_of_sq_new/points(1,2)**2
      else
         Int_of_sq = Int_of_sq_new/points(SIZE(points,1)-2,2)**2
      end IF

      DEALLOCATE(points)

      ! Inward integration
      CALL Set_Integration_Parameters(arg,l,.false.,wave_vector,beg_vib_int=beg_out,vib_n_steps=(n_steps - n_middle_from_min),&
           save_vib_wf=PRESENT(points_total),save_global_pot=save_pot,use_global_pot=use_pot)

      IF(arg%print.GT.0) WRITE(uni_log,FMT='("Inward integration:")')
      CALL Integration(arg,l,energy,wave_vector,points,Int_of_sq_new,no_knots_new)

      IF(PRESENT(points_total))THEN

         IF(arg%save_last_deriv)THEN
            
            DO i=n_middle_from_min+2,n_steps+1
               points_total(i,:) = points(n_steps + 2 - i,:)
               points_total(i,2) = points_total(i,2)*points_total(n_middle_from_min+1,2)/points(n_steps + 1 - n_middle_from_min,2)
               points_total(i,1) = points_total(i,1)/wave_vector
            END DO
            
         ELSE
            DO i=n_middle_from_min+2,n_steps+1
               points_total(i,:) = points(n_steps + 2 - i,:)
               points_total(i,2) = points_total(i,2)*points_total(n_middle_from_min+1,2)/points(n_steps + 1 - n_middle_from_min,2)
               points_total(i,1) = points_total(i,1)/wave_vector
            END DO
            
         END IF

      ELSE
         IF(arg%save_last_deriv)THEN
            deriv_from_right = wave_vector*points(2,2)/points(1,2)
         ELSE
            A_1 = (points(2,2) - points(4,2))/2
            A_2 = (points(1,2) - points(5,2))/2
            B_1 = (arg%step_size**2)*(Sec_Deriv(J,points(2,1),wave_vector,energy)*points(2,2) -&
                 Sec_Deriv(J,points(4,1),wave_vector,energy)*points(4,2))/12
            B_2 = (arg%step_size**2)*(Sec_Deriv(J,points(1,1),wave_vector,energy)*points(1,2) -&
                 Sec_Deriv(J,points(5,1),wave_vector,energy)*points(5,2))/12
            
            deriv_from_right = 16*(-A_1 + (37*A_2)/32 - (37*B_1)/5 - (17*B_2)/40)/(21*ABS(arg%step_size/wave_vector))
            deriv_from_right = deriv_from_right/points(3,2)
         END IF
         
      END IF
      
      IF(arg%save_last_deriv)THEN
         Int_of_sq = Int_of_sq + Int_of_sq_new/points(1,2)**2
      else
         Int_of_sq = Int_of_sq + Int_of_sq_new/points(SIZE(points,1)-2,2)**2
      end IF

      no_knots = no_knots + no_knots_new
      
      IF(PRESENT(points_total))THEN
         delta_energy = zero

         Int_of_sq = Int_of_sq*(points_total(n_middle_from_min+1,2)**2)
         points_total(:,2) = points_total(:,2)/SQRT(Int_of_sq)

         IF(points_total(2,2).LT.zero) points_total(:,2) = -points_total(:,2)
         
      ELSE
         delta_energy = (deriv_from_left - deriv_from_right)/(2*Particles%reduced_mass*Int_of_sq)
      END IF

      very_first_int=.false.
      
      DEALLOCATE(points)

    END SUBROUTINE VibrationalIntegration
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
    SUBROUTINE Get_Initial_Points(vib_energy,ini_point_internal,ini_point_external)

      IMPLICIT NONE

      REAL(kind = real_kind), intent(IN) :: vib_energy
      REAL(kind = real_kind), intent(OUT) :: ini_point_internal,ini_point_external

      if(Job_Data%job%vib_fixed_mesh_points)then
         if(Job_Data%job%r_min.LE.zero)then
            ini_point_internal=Unit_converter(pot%dist_uni,at_uni_dist,Pot%zero_energy_ret_p) - &
                 Unit_Converter(pot%dist_uni,at_uni_dist,Job_Data%job%internal_ini_point_shift)
            if(ini_point_internal.LE.zero)then
               ini_point_internal=Unit_converter(pot%dist_uni,at_uni_dist,Pot%zero_energy_ret_p)
            end if
         else
            ini_point_internal=Unit_converter(pot%dist_uni,at_uni_dist,Job_Data%job%r_min)
         end if
         ini_point_external=Unit_converter(pot%dist_uni,at_uni_dist,Job_Data%job%r_max)

      else
         CALL Classical_Turning_Points(vib_energy,ini_point_internal,ini_point_external)
         if(status.ne.0) return

         ini_point_internal = ini_point_internal - Unit_Converter(pot%dist_uni,at_uni_dist,Job_Data%job%internal_ini_point_shift)
         IF(ini_point_internal.LT.Pot%lower_bound)THEN
            ini_point_internal = Pot%lower_bound
         END IF
         ini_point_external = ini_point_external + Unit_Converter(pot%dist_uni,at_uni_dist,Job_Data%job%external_ini_point_shift)
         
      end if

    END SUBROUTINE Get_Initial_Points
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
    ! True if the energy of the last calculation (represented by no_knots and dif_energy) is lower than the true eigenenergy
    LOGICAL FUNCTION is_lower(no_knots,level,dif_energy)
      
      IMPLICIT NONE
      
      INTEGER, INTENT(in) :: no_knots, level
      REAL(kind = real_kind), INTENT(in) :: dif_energy
      
      IF(no_knots.LT.level)THEN
         is_lower = .TRUE.
         RETURN
      END IF

      IF(no_knots.GT.level)THEN
         is_lower = .FALSE.
         RETURN
      END IF

      IF(dif_energy.GT.zero)THEN
         is_lower = .TRUE.
         RETURN
      END IF

      is_lower = .FALSE.

    END FUNCTION is_lower
    
  END SUBROUTINE BOUND_STATES
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  SUBROUTINE RES_ANALYSIS()
    
    IMPLICIT NONE

    TYPE(TypeResultsList), POINTER :: curr_results=>NULL(), curr_results_2=>NULL()
    TYPE(TypeJobSchedule), POINTER :: That_Data=>NULL(), That_Data_2=>NULL()

    TYPE(TypeIntegerLinkedList), POINTER :: current_ang_mom
    TYPE(TypeRealLinkedList), POINTER :: current_energy => NULL()
    TYPE(TypeRealLinkedList), POINTER :: current_angle => NULL()
    Integer :: i, i_l, i_E, i_A, i_min_1, i_min_2, i_level, i_param

    integer :: total_size, n_rot, n_vib, n_vib_rot

    TYPE(TypeVibrationalResults), POINTER :: current_level_1, current_level_2
    TYPE(TypeVibLevelsTransitions), POINTER :: current_trans=>NULL()

    Real(kind = real_kind) :: r_min, r_max, step_size, sum, energy_min, aux_real, energy_dif!, nu_cm1

    REAL(kind = real_kind), ALLOCATABLE, DIMENSION(:,:) :: points
    CHARACTER(LEN=12) :: param_name

!    REAL(kind = real_kind), parameter :: A_coef_factor = 3.313750913E-7 ! nu in cm-1 and mu in Deb, gives A in s-1
!    REAL(kind = real_kind), parameter :: A_coef_factor = 7.2356E-6 ! nu in cm-1 and mu in e*ANG, gives A in s-1
    REAL(kind = real_kind), parameter :: A_coef_factor = 2.1420067E+10 ! nu in hartree and mu in auDip, gives A in s-1
    INTEGER :: trans_deg_factor

    REAL(kind = real_kind) :: botton_energy_dif
    LOGICAL :: vibronic_trans

    CHARACTER(LEN=1) :: TRANS_lapack
    INTEGER :: M_lapack, N_lapack, NRHS_lapack, LDA_lapack, LDB_lapack, LWORK_lapack, INFO_lapack
    REAL(kind = real_kind), ALLOCATABLE, DIMENSION(:) :: B_lapack, WORK_lapack
    REAL(kind = real_kind), ALLOCATABLE, DIMENSION(:,:) :: A_lapack

    curr_results=>Select_Results(Job_Data%analysis%job_name_1)
    That_Data=>Select_Data(Job_Data%analysis%job_name_1)
    
    if(Job_Data%analysis%job_name_2_def)then
       curr_results_2=>Select_Results(Job_Data%analysis%job_name_2)
       That_Data_2=>Select_Data(Job_Data%analysis%job_name_2)
    end if

    SELECT CASE(Job_Data%analysis%analysis_type)

    CASE("PartialCrossSec")

       ALLOCATE(Results%cross_section(size(curr_results%phase_shifts,1),size(curr_results%phase_shifts,2)))

       WRITE(uni_out,FMT='("+-------------------------------------------------------+")')
       WRITE(uni_out,FMT='("|               Partial Cross Sections                  |")')

       write(uni_out,FMT='("Energy ")',ADVANCE='NO') 
       current_ang_mom=>That_Data%job%l
       do i_l=1,size(curr_results%phase_shifts,2),1
          
          write(uni_out,FMT='(4X,'//TRIM(formats%ang_momentum)//',4X)',ADVANCE='NO') current_ang_mom%value
          
          current_ang_mom=>current_ang_mom%next
       end do
       write(uni_out,FMT='()')

       current_energy=>That_Data%job%energy
       do i_E=1,size(curr_results%phase_shifts,1),1
          current_ang_mom=>That_Data%job%l

          write(uni_out,FMT='('//TRIM(formats%energy)//',1X)',ADVANCE='NO') current_energy%value

          do i_l=1,size(curr_results%phase_shifts,2),1
             
             Results%cross_section(i_E,i_l) = Calc_Cross_Sec((/curr_results%phase_shifts(i_E,i_l)/),curr_results%wave_vector(i_E),&
                  l_single=current_ang_mom%value)
             
             write(uni_out,FMT='('//TRIM(formats%cross_section)//')',ADVANCE='NO') Results%cross_section(i_E,i_l)
          
             current_ang_mom=>current_ang_mom%next
          end do
          write(uni_out,FMT='()')
          current_energy=>current_energy%next
       end do


       WRITE(uni_out,FMT='("+-------------------------------------------------------+")')       
       
    CASE("TotalCrossSec")

       ALLOCATE(Results%cross_section(size(curr_results%phase_shifts,1),1))

       WRITE(uni_out,FMT='("+-------------------------------------------------------+")')
       WRITE(uni_out,FMT='("|                Total Cross Sections                   |")')
       
       current_energy=>That_Data%job%energy
       do i_E=1,size(curr_results%phase_shifts,1),1
          
          Results%cross_section(i_E,1) = Calc_Cross_Sec(curr_results%phase_shifts(i_E,:),curr_results%wave_vector(i_E),&
               l_list=That_Data%job%l)
          
          write(uni_out,FMT='('//TRIM(formats%energy)//',1X,'//TRIM(formats%cross_section)//')')&
               current_energy%value,Results%cross_section(i_E,1)
          
          current_energy=>current_energy%next
          
       end do
       
       WRITE(uni_out,FMT='("+-------------------------------------------------------+")')

    CASE("DifferentialCrossSec")

       ALLOCATE(Results%cross_section(size(curr_results%phase_shifts,1),RLL_Size(Job_Data%analysis%angle)))

       WRITE(uni_out,FMT='("+-------------------------------------------------------+")')
       WRITE(uni_out,FMT='("|             Differential Cross Sections               |")')

       current_energy=>That_Data%job%energy
       do i_E=1,size(Results%cross_section,1),1

          write(uni_out,FMT='("Energy = ",'//TRIM(formats%energy)//',":")') current_energy%value

          current_angle=>Job_Data%analysis%angle
          do i_A=1,size(Results%cross_section,2),1

             Results%cross_section(i_E,i_A) = Calc_Dif_Cross_Sec(curr_results%phase_shifts(i_E,:),current_angle%value,&
                  curr_results%wave_vector(i_E),Job_Data%analysis%is_center_of_mass_ref)
             
             write(uni_out,FMT='('//TRIM(formats%angle)//',1X,'//TRIM(formats%cross_section)//')') &
                  current_angle%value,Results%cross_section(i_E,i_A)

             current_angle=>current_angle%next
          end DO
          
          current_energy=>current_energy%next
       end do

       WRITE(uni_out,FMT='("+-------------------------------------------------------+")')

    CASE("VibrationalIntensity")

       if(.NOT.Job_Data%analysis%Frank_Condon)then
          CALL Select_Potential(Job_Data%analysis%dipole_moment_name)
       end if

       trans_deg_factor = 1



       ! Same electronic state
       if(.not.associated(curr_results_2))then

          curr_results_2=>curr_results
          vibronic_trans = .FALSE.
          
          botton_energy_dif = zero

       else ! Different electronic state
          vibronic_trans = .TRUE.

          if(Job_Data%analysis%elec_trans_ener_def)then
             botton_energy_dif = Unit_Converter(Job_Data%analysis%elec_trans_ener_uni,at_uni_ener,Job_Data%analysis%elec_trans_ener)
          else

             CALL Select_Potential(That_Data_2%job%pot_name)
             botton_energy_dif = Pot%ener_infty - Pot%diss_ener

             if(Pot%lambda.eq.0) trans_deg_factor = 2
             
             CALL Select_Potential(That_Data%job%pot_name)
             botton_energy_dif = botton_energy_dif - Pot%ener_infty + Pot%diss_ener

             if(Pot%lambda.eq.0) trans_deg_factor = 1

          end if

       end if


       if(Job_Data%analysis%Frank_Condon)then
          WRITE(uni_out,FMT='("+----------------------------------------------------------------------------------+")')
          WRITE(uni_out,FMT='("|                            Transition intensities                                |")')
          WRITE(uni_out,FMT='(" Transition   Energy   Frank_Condon_fac    dipole_int          Einstein_Aij_coeff")')
       else
          WRITE(uni_out,FMT='("+------------------------------------------------------------------------+")')
          WRITE(uni_out,FMT='("|                       Transition intensities                           |")')
          WRITE(uni_out,FMT='(" Transition   Energy         dipole_int         Einstein_Aij_coeff")')
       end if

       current_level_1 => curr_results%vib
       do while(associated(current_level_1))
          current_level_2 => curr_results_2%vib
          do while(associated(current_level_2))
             
             if(vibronic_trans.OR.current_level_1%level.lt.current_level_2%level)then

                IF(.NOT.ASSOCIATED(current_trans))THEN
                   ALLOCATE(Results%VibTrans)
                   current_trans=>Results%VibTrans
                ELSE
                   ALLOCATE(current_trans%next)
                   current_trans => current_trans%next
                END IF

                current_trans%state_1 => current_level_1
                current_trans%state_2 => current_level_2

                energy_dif = current_level_2%energy - current_level_1%energy + botton_energy_dif
                
                write(uni_out,FMT='('//trim(formats%vib_level)//'," -> ",'//trim(formats%vib_level)//&
                     ',5X,'//TRIM(formats%energy_cm1)//',3X)',advance='no') &
                     current_level_2%level,current_level_1%level, &
                     Unit_Converter(at_uni_ener,"cm-1           ",energy_dif)
                
                step_size = current_level_1%wave_function(2,1) - current_level_1%wave_function(1,1)
                r_min = max(current_level_1%wave_function(1,1),&
                            current_level_2%wave_function(1,1))
                r_max = min(current_level_1%wave_function(size(current_level_1%wave_function,1),1),&
                            current_level_2%wave_function(size(current_level_2%wave_function,1),1))
                total_size = NINT((r_max - r_min)/step_size) + 1
                i_min_1 = NINT((current_level_1%wave_function(1,1) - r_min)/step_size) + 1
                i_min_2 = NINT((current_level_2%wave_function(1,1) - r_min)/step_size) + 1

                if(allocated(points)) deallocate(points)
                
                allocate(points(total_size,2))
                
                do i=1,total_size,1
                   if(i-i_min_1.ge.1.and.i-i_min_1.le.size(current_level_1%wave_function,1))then
                      points(i,1) = current_level_1%wave_function(i-i_min_1+1,1)

                   else
                      points(i,1) = current_level_2%wave_function(i-i_min_2+1,1)
                   end if
                   
                   if(i-i_min_1.ge.1.and.i-i_min_1.le.size(current_level_1%wave_function,1).and.&
                        i-i_min_2.ge.1.and.i-i_min_2.le.size(current_level_2%wave_function,1))then
                      points(i,2) = current_level_1%wave_function(i-i_min_1+1,2)*current_level_2%wave_function(i-i_min_2+1,2)
                      if(.NOT.Job_Data%analysis%Frank_Condon)then
                         points(i,2) = points(i,2)*Potential(points(i,1))
                      end if
                   else
                      points(i,2) = zero
                   end if

                end do

                sum = zero

                SELECT CASE(quad_met) ! Numerical integration
                   
                case("TRAPEZOIDAL")

                   do i=2,total_size,1
                      sum = sum + (points(i,2) + points(i-1,2))*(points(i,1)-points(i-1,1))/2
                   end do

                CASE("SIMPSON")
                   
                   do i=3,total_size,2
                      sum = sum + (points(i,2) + 4*points(i-1,2) + points(i-2,2))*(points(i,1)-points(i-2,1))/6
                   end do
                   if(i.eq.total_size+1)then
                      sum = sum + (points(i-1,2) + points(i-2,2))*(points(i-1,1) - points(i-2,1))/2 ! Trapezoidal for last interval
                   end if

                CASE DEFAULT
                   error_msg = trim(quad_met)//" is not an know quadrature name."
                   CALL Error(2)
                   
                END SELECT
                
                
                if(Job_Data%analysis%Frank_Condon)then
                   current_trans%FrankCondon_Fac = sum
                   sum = sum*Job_Data%analysis%ave_elec_dip_int
                   write(uni_out,FMT='(5X,'//trim(formats%dip_mom)//',5X)',advance='NO') current_trans%FrankCondon_Fac
                end if

                current_trans%DipMom = sum
                write(uni_out,FMT='(5X,'//trim(formats%dip_mom)//',5X)',advance='NO') current_trans%DipMom

                current_trans%Einstein_A_Coeff = trans_deg_factor*A_coef_factor*(sum**2)*(ABS(energy_dif)**3)
                
                write(uni_out,FMT='(3X,'//trim(formats%einstein_coef)//')') current_trans%Einstein_A_Coeff
                
                deallocate(points)
               
             end if
             
             current_level_2 => current_level_2%next
          end do
          current_level_1 => current_level_1%next
       end do

       WRITE(uni_out,FMT='("+------------------------------------------------------------------------+")')
       
    CASE("DunhanExpansion")

       CALL Select_Potential(That_Data%job%pot_name)

       n_vib = Job_Data%analysis%Dunhan_n_vib
       n_rot = Job_Data%analysis%Dunhan_n_rot
       n_vib_rot = Job_Data%analysis%Dunhan_n_vib_rot

       n_lapack = n_vib + n_rot + n_vib_rot
       m_lapack=0
       current_level_1 => curr_results%vib
       do while(associated(current_level_1))
          current_level_1 => current_level_1%next
          m_lapack = m_lapack+1
       end do

       ALLOCATE(Results%Dunhan_parameters(n_lapack))

       energy_min = Unit_Converter(at_uni_ener,"cm-1           ",curr_results%vib%energy)

       ALLOCATE(A_lapack(M_lapack,n_lapack),B_lapack(m_lapack))
       i_level=1
       current_level_1 => curr_results%vib
       do while(associated(current_level_1))

          ! The vibrational analysis is carried only for J=0
          B_lapack(i_level) = Unit_Converter(at_uni_ener,"cm-1           ",current_level_1%energy)
          
          if(current_level_1%J.eq.0.and.B_lapack(i_level).lt.energy_min) energy_min = B_lapack(i_level)
          
          ! Vibrational expansion
          do i_param=1,n_vib,1
             if(i_param.eq.1)then
                A_lapack(i_level,i_param) = (current_level_1%level + 0.5_real_kind)
             else
                A_lapack(i_level,i_param) = -ABS(A_lapack(i_level,i_param-1)*(current_level_1%level + 0.5_real_kind))
             end if
          end do
          
          ! Rotational expansion
          do i_param=n_vib+1,n_vib+n_rot,1
             if(i_param.EQ.n_vib+1)then
                A_lapack(i_level,i_param) = current_level_1%J*(current_level_1%J + one)
             else
                A_lapack(i_level,i_param) = -ABS(A_lapack(i_level,i_param-1)*current_level_1%J*(current_level_1%J + one))
             end if
          end do
          
          ! Vibro-rotational expansion: assuming n_vib and n_rot different from zero
          do i_param=n_vib+n_rot+1,n_vib+n_rot+n_vib_rot
             if(i_param.EQ.n_vib+n_rot+1)then
                A_lapack(i_level,i_param) = -A_lapack(i_level,1)*A_lapack(i_level,n_vib+1)
             else
                A_lapack(i_level,i_param) = A_lapack(i_level,i_param-1)*A_lapack(i_level,1)
             end if
          end do
          
          i_level=i_level+1
          
          current_level_1 => current_level_1%next
       end do
       
       if(n_vib.eq.0)then
          do i_level=1,size(B_lapack)
             B_lapack(i_level) = B_lapack(i_level)-energy_min
          end do
       end if

       TRANS_lapack = 'N'
       M_lapack = size(A_lapack,1)
       N_lapack = size(A_lapack,2)
       NRHS_lapack = 1
       LDA_lapack = MAX(1,M_lapack)
       LDB_lapack = MAX(1,M_lapack,N_lapack)
       LWORK_lapack = max( 1, M_lapack*N_lapack + max( M_lapack*N_lapack, NRHS_lapack ) )
       ALLOCATE(WORK_lapack(MAX(1,LWORK_lapack)))

       CALL DGELS( TRANS_lapack, M_lapack, N_lapack, NRHS_lapack, A_lapack, LDA_lapack, B_lapack, LDB_lapack, &
            WORK_lapack, LWORK_lapack, INFO_lapack )
       
       If(INFO_lapack.ne.0)then
          error_msg = "Error in LAPACK-DGELS - INFO_lapack not zero."
          CALL Error(2)
       end If

       Results%Dunhan_parameters = B_lapack(1:N_lapack)

       WRITE(uni_out,FMT='("+-------------------------------------------------------+")')
       WRITE(uni_out,FMT='("|                  Dunhan Parameters                    |")')
       
       do i_param=1,N_lapack
          CALL Set_Param_Name(i_param,param_name)
          write(uni_out,FMT='(A10," = ",'//trim(formats%dunhan)//')',advance='NO') TRIM(param_name),B_lapack(i_param)
          
          if(param_name.eq."B_e")then
             aux_real=(Unit_Converter("cm-1           ",at_uni_ener,B_lapack(i_param))*Particles%reduced_mass*2)
             if(aux_real.gt.zero)then
                write(uni_out,FMT='("  (r_0 = ",'//TRIM(formats%distance)//',1X,A,")")') &
                     sqrt(1/aux_real),&
                     TRIM(at_uni_dist)
             else
                write(uni_out,FMT='()')
             end if
          else
             write(uni_out,FMT='()')
          end if
          
       end do

       WRITE(uni_out,FMT='("+-------------------------------------------------------+")')
       
    CASE("EquilibriumDistance")

       WRITE(uni_out,FMT='("+-------------------------------------------------------+")')
       WRITE(uni_out,FMT='("|      Vibrationally averaged equilibrium distances     |")')

       current_level_1 => curr_results%vib
       do while(associated(current_level_1))
          
          write(uni_out,FMT='("Level ",'//trim(formats%vib_level)//',", J ",'//trim(formats%vib_level)//',":")',advance='no')&
               current_level_1%level,current_level_1%J
          
          if(allocated(points)) deallocate(points)
          
          allocate(points(size(current_level_1%wave_function,1),2))
          
          do i=1,size(points,1),1
             points(i,1) = current_level_1%wave_function(i,1)
             points(i,2) = points(i,1)*(current_level_1%wave_function(i,2)**2)
          end do
                
          current_level_1%eq_dist = zero
          do i=2,size(points,1),1
             current_level_1%eq_dist = current_level_1%eq_dist + (points(i,2)+points(i-1,2))*(points(i,1)-points(i-1,1))/2
          end do
          write(uni_out,FMT='(2X,'//TRIM(formats%distance)//')') current_level_1%eq_dist
          
          deallocate(points)
               
          current_level_1 => current_level_1%next
       end do
       
       WRITE(uni_out,FMT='("+-------------------------------------------------------+")')

    CASE DEFAULT

       error_msg = "Inconsistency analysis type: "//TRIM(Job_Data%analysis%analysis_type)
       CALL Error(2)
       
    END SELECT

    curr_results=>NULL()
    curr_results_2=>NULL()
    That_Data=>NULL()
    That_Data_2=>NULL()
    current_energy => NULL()
    current_angle => NULL()
    current_trans=>NULL()

  Contains
    FUNCTION Select_Results(job_name)

      IMPLICIT NONE

      TYPE(TypeResultsList), POINTER :: Select_Results
      CHARACTER(LEN = 100) :: job_name

      Select_Results=>Results_list
      DO WHILE(ASSOCIATED(Select_Results))
         if(Select_Results%name.eq.job_name) return
         Select_Results=>Select_Results%next
      END DO

      error_msg = "Unknown result name."
      CALL Error(0)

    END FUNCTION Select_Results

    FUNCTION Select_Data(job_name)

      IMPLICIT NONE

      TYPE(TypeJobSchedule), POINTER :: Select_Data
      CHARACTER(LEN = 100) :: job_name

      Select_Data=>Jobs_Data_list
      DO WHILE(ASSOCIATED(Select_Data))
         if(Select_Data%job_name.eq.job_name) return
         Select_Data=>Select_Data%next
      END DO

      error_msg = "Unknown result name."
      CALL Error(0)

    END FUNCTION Select_Data

    SUBROUTINE Set_Param_Name(index,name)

      IMPLICIT NONE

      INTEGER, intent(in) :: index
      character(len=*) :: name

      INTEGER :: index_shifted

      if(index.le.n_vib)then

         select case(index)
         case(1)
            name="W_e"
         case(2)
            name="W_e.X_e"
         case(3)
            name="W_e.Y_e"
         case(4)
            name="W_e.Z_e"
         case default
            write(name,fmt='("W_e.X_e(",I0,")")') index
         end select

         return
            
      end if

      if(index.le.n_vib+n_rot)then

         index_shifted=index-n_vib

         select case(index_shifted)
         case(1)
            name="B_e"
         case(2)
            name="D_e"
         case(3)
            name="H_e"
         case default
            write(name,fmt='("B_e(",I0,")")') index_shifted
         end select

         return
            
      end if

      if(index.le.n_vib+n_rot+n_vib_rot)then

         index_shifted=index-(n_vib+n_rot)

         select case(index_shifted)
         case(1)
            name="alpha_e"
         case(2)
            name="gamma_e"
         case default
            write(name,fmt='("alpha_e(",I0,")")') index_shifted
         end select

         return
            
      end if

      error_msg = "Trying to set parameter name for a index greater than the number of parameters asked."
      CALL Error(2)

    END SUBROUTINE Set_Param_Name

  END SUBROUTINE RES_ANALYSIS
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
    SUBROUTINE Print_energy_cm_and_ua(print_energy,advance,relative_energy)

      IMPLICIT NONE

      Real(kind = real_kind), intent(in) :: print_energy
      Logical, intent(in), optional :: advance
      Logical, intent(in), optional :: relative_energy

      Real(kind = real_kind) :: print_energy_cm

      if(present(relative_energy).and.relative_energy)then
         print_energy_cm = Unit_Converter(at_uni_ener,"cm-1           ",Pot%diss_ener+print_energy)
      else
         print_energy_cm = Unit_Converter(at_uni_ener,"cm-1           ",print_energy)
      end if

      if(present(advance).and..NOT.advance)then
         write(uni_log,FMT='('//TRIM(formats%energy_cm1)//'," cm-1 (",'//TRIM(formats%energy)//',1X,A,")")',ADVANCE='NO')&
              print_energy_cm,print_energy,TRIM(at_uni_ener)
      else
         write(uni_log,FMT='('//TRIM(formats%energy_cm1)//'," cm-1 (",'//TRIM(formats%energy)//',1X,A,")")')&
              print_energy_cm,print_energy,TRIM(at_uni_ener)
      end if

    END SUBROUTINE Print_energy_cm_and_ua
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
!!$  TYPE(TypeIntParameters) FUNCTION Set_General_Parameters()
!!$
!!$    IMPLICIT NONE
!!$
!!$    Set_General_Parameters%method = Job_Data%job%int_method
!!$    Set_General_Parameters%r_bar_mesh = Unit_Converter(Pot%dist_uni,at_uni_dist,Job_Data%job%r_bar_meshkov)
!!$    Set_General_Parameters%beta_mesh = Job_Data%job%beta_meshkov
!!$    Set_General_Parameters%save_at = Job_Data%job%save_at
!!$    Set_General_Parameters%save_points = Job_Data%job%save_steps
!!$    Set_General_Parameters%print = Job_Data%job%print_at
!!$
!!$    IF(Job_Data%job_method.EQ."SCATT_LENGTH".AND.Job_Data%job%scatt_len_method.EQ."MESHKOV")THEN
!!$       Set_General_Parameters%inc_step = -1
!!$       Set_General_Parameters%end_inc_step = -1
!!$       Set_General_Parameters%end_integration = Job_Data%job%add_steps
!!$       Set_General_Parameters%start_save = Set_General_Parameters%end_integration
!!$       Set_General_Parameters%save_at = 1
!!$    ELSE
!!$!!Set_General_Parameters%end_integration = NINT((data%r_max-data%r_min)/data%int_step)
!!$       Set_General_Parameters%inc_step = NINT((Job_Data%job%r_max-Job_Data%job%r_min)/Job_Data%job%int_step)
!!$       IF(Job_Data%job%int_step.GT.zero.AND.Job_Data%job_method.NE."BOUND_STATES")THEN
!!$          Set_General_Parameters%end_inc_step = Set_General_Parameters%inc_step + &
!!$               2*NINT(LOG(Job_Data%job%final_int_step/Job_Data%job%int_step)/LOG(two))
!!$       ELSE
!!$          Set_General_Parameters%end_inc_step=-1
!!$       END IF
!!$       Set_General_Parameters%start_save = Set_General_Parameters%end_inc_step + 1
!!$       Set_General_Parameters%end_integration = Set_General_Parameters%inc_step + Job_Data%job%add_steps 
!!$    END IF
!!$
!!$  END FUNCTION Set_General_Parameters
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
!!$  SUBROUTINE Set_Specific_Parameters(arg,l,wave_vector,save_global_pot,use_global_pot)
!!$
!!$    IMPLICIT NONE
!!$
!!$    TYPE(TypeIntParameters), INTENT(inout) :: arg
!!$    INTEGER, INTENT(in) :: l
!!$    REAL(kind = real_kind), INTENT(in), OPTIONAL :: wave_vector
!!$    LOGICAL, INTENT(in), OPTIONAL :: save_global_pot
!!$    LOGICAL, INTENT(in), OPTIONAL :: use_global_pot
!!$
!!$    REAL(kind = real_kind) :: q_tilde
!!$
!!$    IF(Job_Data%job%use_potential_vector)THEN
!!$       IF(PRESENT(save_global_pot))THEN
!!$          Arg%save_global_potential = save_global_pot
!!$       ELSE
!!$          Arg%save_global_potential = .FALSE.
!!$       END IF
!!$       IF(PRESENT(use_global_pot))THEN
!!$          Arg%use_global_potential = use_global_pot
!!$       ELSE
!!$          Arg%use_global_potential = .FALSE.
!!$       END IF
!!$    ELSE
!!$       Arg%use_global_potential = .FALSE.
!!$       Arg%save_global_potential = .FALSE.
!!$    END IF
!!$
!!$    IF(Job_Data%job_method.EQ."SCATT_LENGTH".AND.Job_Data%job%scatt_len_method.EQ."MESHKOV")THEN
!!$       arg%r_min = r_to_y_meshkov(Job_Data%job%r_min,Job_Data%job%r_bar_meshkov,Job_Data%job%beta_meshkov)
!!$
!!$       q_tilde = Sec_Deriv_tilde(l*(l+1),arg%r_min,arg%r_bar_mesh,arg%beta_mesh)
!!$       arg%initial_wf(1) = arg%step_size*(SQRT(-Q_tilde) - Q_tilde*arg%step_size/3)
!!$
!!$       arg%step_size = (one-arg%r_min)/arg%end_integration
!!$    ELSE
!!$       arg%r_min = Unit_converter(pot%dist_uni,at_uni_dist,Job_Data%job%r_min)
!!$       arg%step_size = Unit_Converter(Pot%dist_uni,at_uni_dist,Job_Data%job%int_step)
!!$
!!$       IF(PRESENT(wave_vector)) arg%r_min = arg%r_min*wave_vector
!!$       IF(PRESENT(wave_vector)) arg%step_size = arg%step_size*wave_vector
!!$       arg%initial_wf = Free_Part_Ini_Cond(arg%r_min,l)
!!$       IF(.NOT.Job_Data%job%ini_cond_free_part) arg%initial_wf(1) = zero
!!$    END IF
!!$    arg%initial_wf = arg%initial_wf*Job_Data%job%wf_scale
!!$
!!$  END SUBROUTINE Set_Specific_Parameters
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  SUBROUTINE Set_Integration_Parameters(arg,l,first_integration,wave_vector,vib_step_size,beg_vib_int,vib_n_steps&
       ,save_vib_wf,save_global_pot,use_global_pot)
    
    IMPLICIT NONE
    
    TYPE(TypeIntParameters), INTENT(inout) :: arg
    INTEGER, INTENT(in) :: l
    LOGICAL, INTENT(IN) :: first_integration
    REAL(kind = real_kind), INTENT(in), OPTIONAL :: wave_vector
    REAL(kind = real_kind), INTENT(in), OPTIONAL :: vib_step_size
    REAL(kind = real_kind), INTENT(in), OPTIONAL :: beg_vib_int
    INTEGER, INTENT(in), OPTIONAL :: vib_n_steps
    LOGICAL, INTENT(in), OPTIONAL :: save_vib_wf
    LOGICAL, INTENT(in), OPTIONAL :: save_global_pot
    LOGICAL, INTENT(in), OPTIONAL :: use_global_pot

    REAL(kind = real_kind) :: q_tilde

    if(first_integration)then
       
       arg%method = Job_Data%job%int_method
       arg%print = Job_Data%job%print_at
       arg%save_points = Job_Data%job%save_steps
       arg%save_last_deriv = .FALSE.
       
       IF((Job_Data%job_method.EQ.SCATT_LENGTH_kw.AND.Job_Data%job%scatt_len_method.NE."MESHKOV")&
            .OR.Job_Data%job_method.EQ.PHASE_SHIFTS_kw)THEN
          
          arg%save_at = Job_Data%job%save_at
          
          arg%inc_step = NINT((Job_Data%job%r_max-Job_Data%job%r_min)/Job_Data%job%int_step)
          
          arg%end_inc_step = arg%inc_step + &
               2*NINT(LOG(Job_Data%job%final_int_step/Job_Data%job%int_step)/LOG(two))
          
          arg%start_save = arg%end_inc_step + 1
          arg%end_integration = arg%inc_step + Job_Data%job%add_steps 
          
       ELSEIF(Job_Data%job_method.EQ.SCATT_LENGTH_kw.AND.Job_Data%job%scatt_len_method.EQ."MESHKOV")THEN
          
          arg%r_bar_mesh = Unit_Converter(Pot%dist_uni,at_uni_dist,Job_Data%job%r_bar_meshkov)
          arg%beta_mesh = Job_Data%job%beta_meshkov
          arg%inc_step = -1
          arg%end_inc_step = -1
          arg%save_at = 1
          
       ELSEIF(Job_Data%job_method.EQ.BOUND_STATES_kw)THEN
          
          arg%inc_step = -1
          arg%end_inc_step = -1
          arg%save_at = 1
          
          if(Job_Data%job%print_vib_wf_iterations) then
             arg%print=abs(Job_Data%job%print_at)
          else
             arg%print=-abs(Job_Data%job%print_at)
          end if
          
          if(Job_Data%job%print_vib_wf_iterations) then
             arg%print=abs(Job_Data%job%print_at)
          else
             arg%print=-abs(Job_Data%job%print_at)
          end if

          arg%initial_wf(2) = one
          arg%initial_wf(1) = zero
          
          arg%initial_wf = arg%initial_wf*Job_Data%job%wf_scale

          arg%step_size = vib_step_size*wave_vector

          IF(Job_Data%job%int_method.NE."NUMEROV".AND..NOT.(present(save_vib_wf).and.save_vib_wf)) arg%save_last_deriv = .TRUE.

       ELSE

          error_msg = "Inconsistency in Set Integration Parameters."
          CALL Error(2)

       END IF
          
    end if

    IF((Job_Data%job_method.EQ.SCATT_LENGTH_kw.AND.Job_Data%job%scatt_len_method.NE."MESHKOV")&
         .OR.Job_Data%job_method.EQ.PHASE_SHIFTS_kw)THEN
       
       arg%r_min = Unit_converter(pot%dist_uni,at_uni_dist,Job_Data%job%r_min)
       arg%step_size = Unit_Converter(Pot%dist_uni,at_uni_dist,Job_Data%job%int_step)
       
       IF(PRESENT(wave_vector)) arg%r_min = arg%r_min*wave_vector
       IF(PRESENT(wave_vector)) arg%step_size = arg%step_size*wave_vector
       
       arg%initial_wf = Free_Part_Ini_Cond(arg%r_min,l)
       IF(.NOT.Job_Data%job%ini_cond_free_part) arg%initial_wf(1) = zero
       
    ELSEIF(Job_Data%job_method.EQ."SCATT_LENGTH".AND.Job_Data%job%scatt_len_method.EQ."MESHKOV")THEN
       arg%r_min = r_to_y_meshkov(Job_Data%job%r_min,Job_Data%job%r_bar_meshkov,Job_Data%job%beta_meshkov)

       q_tilde = Sec_Deriv_tilde(l*(l+1),arg%r_min,arg%r_bar_mesh,arg%beta_mesh)
       arg%initial_wf(1) = arg%step_size*(SQRT(-Q_tilde) - Q_tilde*arg%step_size/3)
       
       arg%end_integration = Job_Data%job%add_steps
       arg%start_save = arg%end_integration

       arg%step_size = (one-arg%r_min)/arg%end_integration
       
    ELSEIF(Job_Data%job_method.EQ.BOUND_STATES_kw)THEN
       
       arg%r_min = beg_vib_int
       arg%r_min = arg%r_min*wave_vector
       
       IF(Job_Data%job%int_method.NE."NUMEROV".AND..NOT.(present(save_vib_wf).and.save_vib_wf)) THEN
          arg%start_save = vib_n_steps
          arg%end_integration = vib_n_steps
       else
          arg%start_save = vib_n_steps - 2
          arg%end_integration = arg%start_save + 4
       end IF

       if(.not.first_integration)then
          arg%step_size = -arg%step_size
       end if

       if(present(save_vib_wf).and.save_vib_wf) arg%start_save = 0

    END IF

    IF(Job_Data%job%use_potential_vector)THEN
       IF(PRESENT(save_global_pot).AND.save_global_pot)THEN
          
          save_Potential = .true.
          if(.not.associated(Total_pot))then
             allocate(Total_pot)
             curr_pot_knot => Total_pot
          end if

       ELSE

          save_Potential = .false.
          
       END IF
       
       IF(PRESENT(use_global_pot).AND.use_global_pot)THEN
          use_Pot_vector = .TRUE.
          if(first_integration) curr_pot_knot => Total_pot
       ELSE
          use_Pot_vector = .FALSE.
       END IF

    ELSE
       save_Potential = .FALSE.
       use_Pot_vector = .FALSE.
    END IF

    arg%initial_wf = arg%initial_wf*Job_Data%job%wf_scale

  END SUBROUTINE Set_Integration_Parameters
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  REAL(kind = real_kind) FUNCTION Set_Energy(energy, wave_vector)

    IMPLICIT NONE

    REAL(kind = real_kind), INTENT(in), OPTIONAL :: energy
    REAL(kind = real_kind), INTENT(out) :: wave_vector
    
    IF(PRESENT(energy))THEN
       Set_energy = Unit_Converter(Job_Data%job%energy_uni,at_uni_ener,energy)
       IF(status.NE.0)THEN
          Set_energy = Particles%reduced_mass*Unit_Converter(Job_Data%job%energy_uni,at_uni_vel,energy)**2/2
          IF(status.NE.0)THEN
             error_msg = "Inconsistent unit in converting velocity to energy."
             CALL Error(2)
          END IF
       END IF

       wave_vector = SIGN(SQRT(ABS(2*Particles%reduced_mass*Set_Energy)),Set_Energy)
!       wave_vector = SQRT(abs(2*Particles%reduced_mass*Set_Energy))

    ELSE
       Set_energy = one
       wave_vector = one
    END IF

  END FUNCTION Set_Energy
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  SUBROUTINE Print_Current_Calculation(is_printable,l,energy,angle)
    
    IMPLICIT NONE

    LOGICAL, INTENT(in) :: is_printable
    INTEGER, INTENT(in), OPTIONAL :: l
    REAL(kind = real_kind), INTENT(in), OPTIONAL :: energy
    REAL(kind = real_kind), INTENT(in), OPTIONAL :: angle

    IF(is_printable.AND.PRESENT(l))THEN
       WRITE(uni_log,FMT='("l = ",'//TRIM(formats%ang_momentum)//')') l
    END IF

    IF(is_printable.AND.PRESENT(energy))THEN
       WRITE(uni_log,FMT='("Energy = ",'//TRIM(formats%energy)//',1X,A)') energy,TRIM(at_uni_ener)
    END IF

    IF(is_printable.AND.PRESENT(angle))THEN
       WRITE(uni_log,FMT='("angle = ",'//TRIM(formats%angle)//')') angle
    END IF

  END SUBROUTINE Print_Current_Calculation
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  FUNCTION Free_Part_Ini_Cond(r,l)

    IMPLICIT NONE

    REAL(kind = real_kind), DIMENSION(2) :: Free_Part_Ini_Cond
    REAL(kind = real_kind), INTENT(in) :: r

    REAL(kind = real_kind) :: DF
    INTEGER, INTENT(in) :: l
    
    IF(l.LE.149)THEN
       DF = Double_Factorial(2*l+1)
    ELSE
       DF = 1.0E+200_real_kind
    END IF

    Free_Part_Ini_Cond(1) = (1/DF)*(r**(l+1))
    Free_Part_Ini_Cond(2) = (l+1)*(1/DF)*(r**l)

  END FUNCTION Free_Part_Ini_Cond
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
END MODULE ModMethods
