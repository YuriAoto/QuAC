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
! Module for input manipulation
!
! The input information of each job are stored in an linked list of TypeJobData defined type. All these variables are stored 
! in user units and are converted to atomic unit only locally in the ModMethods module. Angular momentum, energies and angles
! are also stored in linked lists. For (almost) all variables in TypeJobData there is a corresponding one with _def suffix
! that set .true. if the variable was passed by the user. The input is read by Get_Input function, where some consistency 
! between job and parameter and units checks are carried. Some additional checks and definition of default values are carried
! out by Check_Input_Data subroutine.
!
MODULE ModInput

  USE ModUtil, ONLY : real_kind, zero, at_uni_dist, at_uni_ener, at_uni_vel, at_uni_mass, error_msg, ios_var, ios_ok, ios_eof,&
       ios_eol, units_len, uni_inp_orig, uni_inp, uni_out, uni_log, Error, Unit_Converter, Atomic_Data
  USE ModSystem, ONLY : TypeParticles, particles
  USE ModPot, ONLY : Pot, TypeIntegerLinkedList, TypeRealLinkedList, Free_RLL, POTENTIAL_kw, DIP_MOM_kw, &
       Free_ILL, Non_Consistent_Info, Check_Unit, TypePrintFormats, formats, Select_Potential

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: end_of_input, command, TypeJobSchedule, particles, Job_Data, Jobs_Data_list, PHASE_SHIFTS_kw, SCATT_LENGTH_kw, &
       BOUND_STATES_kw, RES_ANALYSIS_kw, MASSES_kw, RED_MASS_kw, ATOMS_kw, FORMATS_kw, Print_Header, Get_Input, &
       Print_Input_Data, Check_Input_Data, Free_Input_Data_Memory, global_job_name

  ! Data for integration
  TYPE TypeJobData

     ! ANGULAR MOMENTUM VARIABLES:
     TYPE(TypeIntegerLinkedList), POINTER :: l=>NULL()

     ! ENERGY VARIABLES:
     TYPE(TypeRealLinkedList), POINTER :: energy=>NULL()
     CHARACTER (LEN = units_len) :: energy_uni = ""

     ! INTEGRATION VARIABLES:
     CHARACTER(LEN = 15) :: int_method = "NUMEROV"
     LOGICAL :: int_method_def = .FALSE.

     REAL(kind = real_kind) :: int_step = 0.005_real_kind
     LOGICAL :: int_step_def = .FALSE.

     INTEGER :: print_at = -1
     LOGICAL :: print_at_def = .FALSE.

     REAL(kind = real_kind) :: wf_scale = 1.0_real_kind
     LOGICAL :: wf_scale_def = .FALSE.

     CHARACTER(LEN = 40) :: pot_name = ""
     LOGICAL :: pot_name_def = .FALSE.

     LOGICAL :: ini_cond_free_part = .FALSE.
     LOGICAL :: ini_cond_def = .FALSE.

     LOGICAL :: save_steps = .FALSE.
     LOGICAL :: save_steps_def = .FALSE.

     ! PHASE SHIFT KEYWORDS:

     CHARACTER(LEN = 20) :: phase_shift_method = "REGRESSION"
     LOGICAL :: phase_shift_method_def = .FALSE.

     REAL(kind = real_kind) :: final_int_step = 0.005_real_kind
     LOGICAL :: final_int_step_def = .FALSE.

     REAL(kind = real_kind) :: r_min = -1
     LOGICAL :: r_min_def = .FALSE.
     CHARACTER(len = 100) :: r_min_default_explanation = ""
     
     REAL(kind = real_kind) :: r_max = -1
     LOGICAL :: r_max_def = .FALSE.
     CHARACTER(len = 100) :: r_max_default_explanation = ""

     INTEGER :: add_steps = 1000
     LOGICAL :: add_steps_def = .FALSE.

     INTEGER :: save_at = 1
     LOGICAL :: save_at_def = .FALSE.

     LOGICAL :: print_ps = .FALSE.
     LOGICAL :: print_ps_def = .FALSE.

     LOGICAL :: use_potential_vector = .FALSE.
     LOGICAL :: use_potential_vector_def = .FALSE.

     ! SCATTERING LENGHT VARIABLES:

     CHARACTER(LEN = 20) :: scatt_len_method = "DIRECT"
     LOGICAL :: scatt_len_method_def = .FALSE.

     REAL(kind = real_kind) :: r_bar_meshkov = -1
     REAL(kind = real_kind) :: beta_meshkov = -1
     LOGICAL :: meshkov_param_def = .FALSE.

     LOGICAL :: Rich_extr = .FALSE.
     LOGICAL :: Rich_extr_def = .FALSE.

     CHARACTER (LEN = units_len) :: cross_sec_uni = "a0"
     LOGICAL :: cross_sec_uni_def = .FALSE.

     ! VIBRATIONAL VARIABLES:

     TYPE(TypeIntegerLinkedList), POINTER :: vib_lev=>NULL()

     LOGICAL :: vib_fixed_mesh_points = .FALSE.
     LOGICAL :: vib_fixed_mesh_points_def = .FALSE.

     LOGICAL :: print_vib_iterations = .FALSE.
     LOGICAL :: print_vib_iterations_def = .FALSE.

     LOGICAL :: print_vib_wf_iterations = .FALSE.
     LOGICAL :: print_vib_wf_iterations_def = .FALSE.

     INTEGER :: vib_1st_order_max_iter = 20
     LOGICAL :: vib_1st_order_max_iter_def = .FALSE.

     INTEGER :: vib_max_iter = 50
     LOGICAL :: vib_max_iter_def = .FALSE.

     REAL(kind = real_kind) :: vib_threshold = 1.0E-5_real_kind
     LOGICAL :: vib_threshold_def = .FALSE.

     REAL(kind = real_kind) :: external_ini_point_shift = 1.0_real_kind
     REAL(kind = real_kind) :: internal_ini_point_shift = 0.7_real_kind
     LOGICAL :: ini_point_shift_def = .FALSE.
     
  END TYPE TypeJobData

  ! Data for analysis
  TYPE TypeAnalysisData
     CHARACTER(LEN = 30) :: analysis_type
     LOGICAL :: analysis_type_def = .FALSE.

     CHARACTER(LEN = 100) :: job_name_1
     LOGICAL :: job_name_1_def = .FALSE.

     CHARACTER(LEN = 100) :: job_name_2
     LOGICAL :: job_name_2_def = .FALSE.

     CHARACTER(LEN = 100) :: dipole_moment_name
     LOGICAL :: dipole_moment_name_def = .FALSE.

     LOGICAL :: Frank_Condon = .FALSE.
     LOGICAL :: Frank_Condon_def = .FALSE.
     REAL(kind = real_kind) :: ave_elec_dip_int
     CHARACTER (LEN = units_len) :: ave_elec_dip_int_uni = ""

     REAL(kind = real_kind) :: elec_trans_ener
     CHARACTER (LEN = units_len) :: elec_trans_ener_uni = ""
     LOGICAL :: elec_trans_ener_def = .FALSE.

     INTEGER :: Dunhan_n_vib = 0
     INTEGER :: Dunhan_n_rot = 0
     INTEGER :: Dunhan_n_vib_rot = 0
     LOGICAL :: Dunhan_def = .FALSE.

     TYPE(TypeRealLinkedList), POINTER :: angle=>NULL()

     LOGICAL :: is_center_of_mass_ref = .FALSE.
     LOGICAL :: is_center_of_mass_ref_def = .FALSE.

     CHARACTER (LEN = units_len) :: cross_sec_uni = "a0"
     LOGICAL :: cross_sec_uni_def = .FALSE.

  END TYPE TypeAnalysisData

  ! Schedule knot
  TYPE TypeJobSchedule
     CHARACTER(LEN = 30) :: job_method
     CHARACTER(LEN = 130) :: complete_job_name
     CHARACTER(LEN = 100) :: job_name

     TYPE(TypeParticles), POINTER :: job_particles=>NULL()
     TYPE(TypePrintFormats), POINTER :: job_formats=>NULL()

     TYPE(TypeJobData), POINTER :: job=>NULL()
     TYPE(TypeAnalysisData), POINTER :: analysis=>NULL()

     TYPE(TypeJobSchedule), POINTER :: next=>NULL()
  END TYPE TypeJobSchedule

  ! Name of the entire job
  CHARACTER (LEN = 100) :: global_job_name

  ! Schedule
  TYPE(TypeJobSchedule), PROTECTED, SAVE, POINTER :: Jobs_Data_list=>NULL()
  TYPE(TypeJobSchedule), POINTER :: Job_Data=>NULL()

  ! end of input
  LOGICAL :: end_of_input = .FALSE.
  
  CHARACTER(LEN=100) :: command

  CHARACTER(len = 30), PARAMETER :: PHASE_SHIFTS_kw="PHASE_SHIFTS"
  CHARACTER(len = 30), PARAMETER :: SCATT_LENGTH_kw="SCATT_LENGTH"
  CHARACTER(len = 30), PARAMETER :: BOUND_STATES_kw="BOUND_STATES"
  CHARACTER(len = 30), PARAMETER :: RES_ANALYSIS_kw="RES_ANALYSIS"

  CHARACTER(len = 30), PARAMETER :: MASSES_kw="MASSES"
  CHARACTER(len = 30), PARAMETER :: RED_MASS_kw="RED_MASS"
  CHARACTER(len = 30), PARAMETER :: ATOMS_kw="ATOMS"

  CHARACTER(len = 30), PARAMETER :: FORMATS_kw="FORMATS"

CONTAINS

!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  SUBROUTINE Print_Header()
    
    WRITE(uni_out,FMT='("    QuAC - Quantum Atomic Colision")')
    WRITE(uni_out,*)
    WRITE(uni_out,FMT='("    This software study  the continuum  spectra of two particle species,")')
    WRITE(uni_out,FMT='("    integrating the radial  Schrödinger equation, obtaining phase shifts")')
    WRITE(uni_out,FMT='("    and total and differential cross sections. It also does calculations")')
    WRITE(uni_out,FMT='("    of the bound levels.")')
    WRITE(uni_out,*)
    WRITE(uni_out,FMT='("    QuAC  Copyright (C) 2009, 2010, 2011, 2012  Yuri Alexandre Aoto")')
    WRITE(uni_out,FMT='("    This program comes with ABSOLUTELY NO WARRANTY; for details type `",A," -c''.")') TRIM(command)
    WRITE(uni_out,FMT='("    This is free software, and you are welcome to redistribute it")')
    WRITE(uni_out,FMT='("    under certain conditions; type `",A," -c'' for details.")') TRIM(command)
    WRITE(uni_out,*)
    WRITE(uni_out,FMT='("---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/---")')
    WRITE(uni_out,*)
    WRITE(uni_out,FMT='("Input file:")')
    WRITE(uni_out,FMT='("---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/---")')

    ! Print input file and take out the comments
    CALL Comments_Take_Out()

    WRITE(uni_out,FMT='("---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/---")')
    WRITE(uni_out,*)
    
    WRITE(uni_log,FMT='("    QuAC - Quantum Atomic Colision        ")')
    WRITE(uni_log,FMT='("    Copyright, 2009, 2010, Yuri Alexandre Aoto  ")')
    WRITE(uni_log,*)
  
  END SUBROUTINE Print_Header
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  SUBROUTINE Comments_Take_Out()
    
    INTEGER :: ios_var_prev=0
    CHARACTER(LEN=1) :: char, char_prev
    LOGICAL :: f90_comment = .FALSE., C_comment = .FALSE.
    
    READ(uni_inp_orig,'(A1)',advance='NO',iostat=ios_var) char_prev
    DO
       
       READ(uni_inp_orig,'(A1)',advance='NO',iostat=ios_var) char
       
       ! f90-style comment
       IF(.NOT.f90_comment.AND.char_prev.EQ."!".AND..NOT.C_comment)THEN
          f90_comment = .TRUE.
          WRITE(uni_out,'(25X,A1,"COMMENT:")',advance='NO') char_prev
          char_prev = char
          READ(uni_inp_orig,'(A1)',advance='NO',iostat=ios_var) char
       END IF

       ! C-style comment 
       IF(.NOT.f90_comment)THEN
          
          IF(.NOT.C_comment.AND.char_prev.EQ."/".AND.char.EQ."*") THEN
             C_comment=.TRUE.
             WRITE(uni_out,'(A1,A1,"BEGIN OF COMMENT:")',advance='NO') char_prev,char
             READ(uni_inp_orig,'(A1)',advance='NO',iostat=ios_var) char_prev
             READ(uni_inp_orig,'(A1)',advance='NO',iostat=ios_var) char
          END IF
          
          IF(C_comment.AND.char_prev.EQ."*".AND.char.EQ."/") THEN
             WRITE(uni_out,'(";END OF COMMENT",A1,A1)',advance='NO') char_prev,char
             C_comment=.FALSE.
             READ(uni_inp_orig,'(A1)',advance='NO',iostat=ios_var) char_prev
             READ(uni_inp_orig,'(A1)',advance='NO',iostat=ios_var) char
          END IF
         
       END IF

       IF(ios_var_prev.EQ.ios_ok)THEN ! Normal character

          WRITE(uni_out,'(A1)',advance='NO') TRIM(char_prev)

          IF(.NOT.C_comment.AND..NOT.f90_comment)THEN
             WRITE(uni_inp,'(A1)',advance='NO') TRIM(char_prev)
          END IF

       END IF

       IF(ios_var_prev.EQ.ios_eol)THEN ! End of line

          WRITE(uni_out,'("")')

          f90_comment = .FALSE.
          IF(.NOT.C_comment)THEN
             WRITE(uni_inp,'("")')
          END IF

       END IF

       IF(ios_var_prev.EQ.ios_eof)THEN ! End of file

          EXIT

       END IF

       ios_var_prev=ios_var
       char_prev=char

    END DO
    
  END SUBROUTINE Comments_Take_Out
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  SUBROUTINE Get_Input()
    
    IMPLICIT NONE

    LOGICAL, SAVE :: first_read = .TRUE.
    TYPE(TypeJobSchedule), POINTER :: current_sched_knot=>NULL()

    TYPE(TypeIntegerLinkedList), POINTER :: current_int => NULL()
    TYPE(TypeRealLinkedList), POINTER :: current_real =>NULL()

    CHARACTER(LEN = 30) :: info
    INTEGER :: aux_int, aux_int2
    REAL(kind = real_kind) :: aux_real
    CHARACTER (LEN = units_len) :: cross_sec_uni_tmp
    REAL(kind = real_kind) :: ini,fim,delta

    ! First read of input file
    IF(first_read)THEN
       
       ! Get job name
       READ(uni_inp,FMT='(A100)',Iostat=ios_var) global_job_name
       IF(ios_var.NE.0)THEN
          error_msg = "Empty input file."
          CALL Error(0)
       END IF
       
       ! To avoid confusion, global_job_name cannot begin with a keyword
       IF(  INDEX(ADJUSTL(global_job_name),TRIM(PHASE_SHIFTS_kw)).EQ.1.OR.&
            INDEX(ADJUSTL(global_job_name),TRIM(SCATT_LENGTH_kw)).EQ.1.OR.&
            INDEX(ADJUSTL(global_job_name),TRIM(BOUND_STATES_kw)).EQ.1.OR.&
            INDEX(ADJUSTL(global_job_name),TRIM(RES_ANALYSIS_kw)).EQ.1.OR.&
            INDEX(ADJUSTL(global_job_name),TRIM(POTENTIAL_kw))   .EQ.1.OR.&
            INDEX(ADJUSTL(global_job_name),TRIM(DIP_MOM_kw))     .EQ.1.OR.&
            INDEX(ADJUSTL(global_job_name),TRIM(MASSES_kw))      .EQ.1.OR.&
            INDEX(ADJUSTL(global_job_name),TRIM(RED_MASS_kw))    .EQ.1.OR.&
            INDEX(ADJUSTL(global_job_name),TRIM(ATOMS_kw))       .EQ.1.OR.&
            INDEX(ADJUSTL(global_job_name),TRIM(FORMATS_kw))     .EQ.1)THEN
          error_msg = TRIM(global_job_name)//" is not an allowed input name."
          CALL Error(0)
       END IF

       ! Allocate memory and point Job_Data to this new allocated area, where the datas will be stored
       ALLOCATE(Jobs_Data_list)
       Job_Data => Jobs_Data_list

       ! The first "job" is the print formats
       ALLOCATE(Job_Data%job_formats)
       Job_Data%job_method = FORMATS_kw
       Job_Data%complete_job_name = FORMATS_kw
       formats => Job_Data%job_formats
       CALL Set_default_formats()

    END IF

    READ(uni_inp,*,Iostat=ios_var) info

    IF(ios_var.NE.0)THEN
       end_of_input = .TRUE.
       first_read=.FALSE.
       RETURN
    END IF

    ! Allocate memory for next job, mass of formats
    SELECT CASE(info)

    CASE(POTENTIAL_kw,DIP_MOM_kw)
       end_of_input = .TRUE.
       !!?? RETURN

    CASE(PHASE_SHIFTS_kw,SCATT_LENGTH_kw,BOUND_STATES_kw,RES_ANALYSIS_kw,MASSES_kw,RED_MASS_kw,ATOMS_kw,FORMATS_kw)

       ALLOCATE(Job_Data%next)
       Job_Data => Job_Data%next
       Job_Data%job_method = info

       SELECT CASE(Job_Data%job_method)

       CASE(FORMATS_kw)
          Job_Data%complete_job_name = Job_Data%job_method
          ALLOCATE(Job_Data%job_formats)
          Job_Data%job_formats = formats
          formats => Job_Data%job_formats
       
       CASE(MASSES_kw,RED_MASS_kw,ATOMS_kw)
          Job_Data%complete_job_name = Job_Data%job_method
          ALLOCATE(Job_Data%job_particles)
          particles => Job_Data%job_particles

       CASE(PHASE_SHIFTS_kw,SCATT_LENGTH_kw,BOUND_STATES_kw)
          ALLOCATE(Job_Data%job)
          
       CASE(RES_ANALYSIS_kw)
          ALLOCATE(Job_Data%analysis)
         
       CASE DEFAULT
          error_msg = "Inconsistecy of method - "//TRIM(Job_Data%job_method)
          CALL Error(2)

       END SELECT
       
       ! Read jobs, masses or formats information
       SELECT CASE(info)
          
       CASE(FORMATS_kw)

          DO ! Loop over formats
             
             READ(uni_inp,*) info

             SELECT CASE(info)
             
             CASE("energy")
                BACKSPACE(uni_inp)
                READ(uni_inp,*,Iostat=ios_var) info, formats%energy
                formats%energy_def = .TRUE.
                
             CASE("energy_cm-1")
                BACKSPACE(uni_inp)
                READ(uni_inp,*,Iostat=ios_var) info, formats%energy_cm1
                formats%energy_cm1_def = .TRUE.
                
             CASE("cross_sec")
                BACKSPACE(uni_inp)
                READ(uni_inp,*,Iostat=ios_var) info, formats%cross_section
                formats%cross_section_def = .TRUE.
                
             CASE("angle")
                BACKSPACE(uni_inp)
                READ(uni_inp,*,Iostat=ios_var) info, formats%angle
                formats%angle_def = .TRUE.
                
             CASE("ang_momentum")
                BACKSPACE(uni_inp)
                READ(uni_inp,*,Iostat=ios_var) info, formats%ang_momentum
                formats%ang_momentum_def = .TRUE.
                
             CASE("distance")
                BACKSPACE(uni_inp)
                READ(uni_inp,*,Iostat=ios_var) info, formats%distance
                formats%distance_def = .TRUE.
                
             CASE("wave_function")
                BACKSPACE(uni_inp)
                READ(uni_inp,*,Iostat=ios_var) info, formats%wave_function
                formats%wave_function_def = .TRUE.
                
             CASE("potential")
                BACKSPACE(uni_inp)
                READ(uni_inp,*,Iostat=ios_var) info, formats%potential
                formats%potential_def = .TRUE.
                
             CASE("Einstein_coef")
                BACKSPACE(uni_inp)
                READ(uni_inp,*,Iostat=ios_var) info, formats%einstein_coef
                formats%einstein_coef_def = .TRUE.
                
             CASE("dipole_mom")
                BACKSPACE(uni_inp)
                READ(uni_inp,*,Iostat=ios_var) info, formats%dip_mom
                formats%dip_mom_def = .TRUE.
                
             CASE("phase_shift")
                BACKSPACE(uni_inp)
                READ(uni_inp,*,Iostat=ios_var) info, formats%phase_shift
                formats%phase_shift_def = .TRUE.
                
             CASE("Dunhan_param")
                BACKSPACE(uni_inp)
                READ(uni_inp,*,Iostat=ios_var) info, formats%dunhan
                formats%dunhan_def = .TRUE.
                
             CASE("vib_level")
                BACKSPACE(uni_inp)
                READ(uni_inp,*,Iostat=ios_var) info, formats%vib_level
                formats%vib_level_def = .TRUE.
                
             CASE("default")
                
                CALL Set_default_formats()
                formats%default_values = .TRUE.
                
             CASE default

                BACKSPACE(uni_inp)
                EXIT
                
             END SELECT

          END DO

       CASE(MASSES_kw,RED_MASS_kw,ATOMS_kw)

          BACKSPACE(uni_inp)

          SELECT CASE(info)
             
          CASE(MASSES_kw)
             
             READ(uni_inp,*,Iostat=ios_var) info, Particles%mass1, Particles%mass2, Particles%mass_uni
             IF(ios_var.NE.0)THEN
                error_msg = "Error in reading mass."
                CALL Error(0)
             END IF
             
             IF(Particles%mass1.LE.0)THEN
                error_msg = "Masses should be positive."
                CALL Error(0)
             END IF
             IF(Particles%mass2.LE.0)THEN
                error_msg = "Masses should be positive."
                CALL Error(0)
             END IF
             
             aux_real = Unit_Converter(Particles%mass_uni,at_uni_mass)
             CALL Check_Unit(Particles%mass_uni,"mass")
             
             Particles%reduced_mass = Particles%mass1*Particles%mass2/(Particles%mass1 + Particles%mass2)
             
          CASE(RED_MASS_kw)
             
             READ(uni_inp,*,Iostat=ios_var) info, Particles%reduced_mass, Particles%mass_uni
             IF(ios_var.NE.0)THEN
                error_msg = "Error in reading mass."
                CALL Error(0)
             END IF
             
             IF(Particles%reduced_mass.LE.0)THEN
                error_msg = "Reduced mass should be positive."
                CALL Error(0)
             END IF
             
             aux_real = Unit_Converter(Particles%mass_uni,at_uni_mass)
             CALL Check_Unit(Particles%mass_uni,"mass")
             
          CASE(ATOMS_kw)
             
             Particles%mass_uni = "amu"
             READ(uni_inp,*) info, Particles%atom1, Particles%atom2
             
             CALL Atomic_Data(Particles%atom1, mass=Particles%mass1, name=Particles%name1, Z=Particles%Z1, &
                  abundance=Particles%abundance1)
             IF(Particles%mass1.LT.zero)THEN
                error_msg = TRIM(Particles%atom1)//" is an unknow atomic symbol."
                CALL Error(0)
             END IF
             
             CALL Atomic_Data(Particles%atom2, mass=Particles%mass2, name=Particles%name2, Z=Particles%Z2, &
                  abundance=Particles%abundance2)
             IF(Particles%mass2.LT.zero)THEN
                error_msg = TRIM(Particles%atom2)//" is an unknow atomic symbol."
                CALL Error(0)
             END IF
             
             Particles%reduced_mass = Particles%mass1*Particles%mass2/(Particles%mass1 + Particles%mass2)
             
          CASE DEFAULT
             
             error_msg = "Inconsistecy of type of mass input - "//TRIM(Job_Data%job_method)
             CALL Error(2)
             
          END SELECT
          
       CASE(PHASE_SHIFTS_kw,SCATT_LENGTH_kw,BOUND_STATES_kw,RES_ANALYSIS_kw)

          BACKSPACE(uni_inp)
          READ(uni_inp,FMT='(A130)') Job_Data%complete_job_name
          READ(Job_Data%complete_job_name(LEN_TRIM(Job_Data%job_method)+1:),FMT='(A100)',iostat=ios_var) Job_Data%job_name
          Job_Data%job_name = ADJUSTL(Job_Data%job_name)
          
          current_sched_knot=>Jobs_Data_list
          DO WHILE(ASSOCIATED(current_sched_knot%next))
             IF(current_sched_knot%job_name.NE.''.AND.current_sched_knot%job_name.EQ.Job_Data%job_name)THEN
                error_msg = "The job name """//TRIM(Job_Data%job_name)//""" was already used in "//&
                     TRIM(current_sched_knot%complete_job_name)
                CALL Error(0)
             END IF
             current_sched_knot=>current_sched_knot%next
          END DO
          
          DO ! Loop over job informations
             READ(uni_inp,*,Iostat=ios_var) info
             IF(ios_var.NE.0)THEN
                EXIT
             END IF
             
             SELECT CASE(info)
!
!       CASE("KEYWORD")
!          
!          verify if KEYWORD was already passed in input for this job, with CheckAlreadyPassed subroutine
!          verify consistency of job_method, with CheckjobMethodConsistency subroutine
!          allocate memmory for new data, if needed
!          read parameters from this keyword
!          if number, verify iostat status, with CheckIostatVar subroutine
!          verify consistency of values (like negative values). Use CheckGreaterThanZero and CheckNegativeNumber subroutines
!          verify consistency of units. Use Unit_Converter and Check_Unit subprograms
!          set this keyword as already passed
!          let the uni_inp ready to read next keyword
!
                
                ! ANGULAR MOMENTUM KEYWORDS:

                ! List of values for l
             CASE("l")

                CALL CheckjobMethodConsistency((/BOUND_STATES_kw,PHASE_SHIFTS_kw/))
                CALL CheckAlreadyPassed(ASSOCIATED(Job_Data%job%l))

                current_int => NULL()

                DO

                   READ(uni_inp,*,Iostat=ios_var) aux_int
                   IF(ios_var.NE.0) EXIT
                   CALL CheckNegativeNumber(int=aux_int)

                   IF(.NOT.ASSOCIATED(current_int))THEN
                      ALLOCATE(Job_Data%job%l)
                      current_int=>Job_Data%job%l
                   ELSE
                      ALLOCATE(current_int%next)
                      current_int=>current_int%next
                   END IF

                   current_int%value = aux_int

                END DO

                CALL CheckPassedValue(ASSOCIATED(Job_Data%job%l))

                BACKSPACE(uni_inp)

                ! max value of l
             CASE("lmax")

                CALL CheckjobMethodConsistency((/BOUND_STATES_kw,PHASE_SHIFTS_kw/))
                CALL CheckAlreadyPassed(ASSOCIATED(Job_Data%job%l),"l")

                BACKSPACE(uni_inp)
                READ(uni_inp,*,Iostat=ios_var) info,aux_int
                CALL CheckIostatVar()
                CALL CheckNegativeNumber(int=aux_int)

                current_int => NULL()

                DO aux_int2=0,aux_int,1

                   IF(.NOT.ASSOCIATED(current_int))THEN
                      ALLOCATE(Job_Data%job%l)
                      current_int=>Job_Data%job%l
                   ELSE
                      ALLOCATE(current_int%next)
                      current_int=>current_int%next
                   END IF

                   current_int%value = aux_int2

                END DO

                ! ENERGY KEYWORDS:

                ! List of values for energy or velocity
             CASE("Energy","Velocity")

                select case(info)
                case("Energy")
                   CALL CheckjobMethodConsistency((/PHASE_SHIFTS_kw,SCATT_LENGTH_kw,BOUND_STATES_kw/))
                case("Velocity")
                   CALL CheckjobMethodConsistency((/PHASE_SHIFTS_kw,SCATT_LENGTH_kw/))
                end select

                CALL CheckAlreadyPassed(ASSOCIATED(Job_Data%job%energy),"Energy/Velocity")

                BACKSPACE(uni_inp)
                READ(uni_inp,*) info,Job_Data%job%energy_uni

                SELECT CASE(info)

                CASE("Energy")
                   aux_real = Unit_Converter(Job_Data%job%energy_uni,at_uni_ener)
                   CALL Check_Unit(Job_Data%job%energy_uni,"energy")

                CASE("Velocity")
                   aux_real = Unit_Converter(Job_Data%job%energy_uni,at_uni_vel)
                   CALL Check_Unit(Job_Data%job%energy_uni,"velocity")

                CASE DEFAULT

                   CALL Non_Consistent_Info("Energy/Velocity","unit check in Energy/Velocity in "//TRIM(Job_Data%complete_job_name))

                END SELECT

                current_real => NULL()
                aux_int=0

                DO

                   READ(uni_inp,*,Iostat=ios_var) aux_real
                   IF(ios_var.NE.0) EXIT
                   CALL CheckGreaterThanZero(float=aux_real)

                   IF(.NOT.ASSOCIATED(current_real))THEN
                      ALLOCATE(Job_Data%job%energy)
                      current_real=>Job_Data%job%energy
                   ELSE
                      ALLOCATE(current_real%next)
                      current_real=>current_real%next
                   END IF

                   aux_int = aux_int + 1
                   current_real%value = aux_real

                END DO

                CALL CheckPassedValue(ASSOCIATED(Job_Data%job%energy))

                BACKSPACE(uni_inp)

                ! Interval of expoents for energy, equally spaced in log scale   
             CASE("ExpInterEnergy","ExpInterVelocity")

                select case(info)
                case("ExpInterEnergy")
                   CALL CheckjobMethodConsistency((/PHASE_SHIFTS_kw,SCATT_LENGTH_kw,BOUND_STATES_kw/))
                case("ExpInterVelocity")
                   CALL CheckjobMethodConsistency((/PHASE_SHIFTS_kw,SCATT_LENGTH_kw/))
                end select

                CALL CheckAlreadyPassed(ASSOCIATED(Job_Data%job%energy),"Energy/Velocity")

                BACKSPACE(uni_inp)
                READ(uni_inp,*,Iostat=ios_var) info,ini,fim,delta,Job_Data%job%energy_uni
                CALL CheckIostatVar()

                SELECT CASE(info)

                CASE("ExpInterEnergy")
                   aux_real = Unit_Converter(Job_Data%job%energy_uni,at_uni_ener)
                   CALL Check_Unit(Job_Data%job%energy_uni,"energy")

                CASE("ExpInterVelocity")
                   aux_real = unit_Converter(Job_Data%job%energy_uni,at_uni_vel)
                   CALL Check_Unit(Job_Data%job%energy_uni,"velocity")

                CASE DEFAULT

                   CALL Non_Consistent_Info("ExpInterEnergy/ExpInterVelocity",&
                        "unit check in Energy/Velocity in "//TRIM(Job_Data%complete_job_name))

                END SELECT

                IF((fim-ini)*delta.LT.0)THEN
                   error_msg = "The variation of the Energy/Velocity shoud be consistent with the initial and final values, in "&
                        &//TRIM(Job_Data%complete_job_name)
                   CALL Error(0)
                END IF

                IF(ABS(delta).GT.ABS(fim-ini))THEN
                   error_msg = "Variation of Energy/Velocity greater than interval in "//TRIM(Job_Data%complete_job_name)
                   CALL Error(0)
                END IF

                current_real => NULL()

                DO aux_int = 1,NINT((fim-ini)/delta)+1,1

                   IF(.NOT.ASSOCIATED(current_real))THEN
                      ALLOCATE(Job_Data%job%energy)
                      current_real=>Job_Data%job%energy
                   ELSE
                      ALLOCATE(current_real%next)
                      current_real=>current_real%next
                   END IF

                   current_real%value = 10**(ini + (aux_int-1)*delta)

                END DO

                ! Interval for energy, equally spaced
             CASE("InterEnergy","InterVelocity")

                select case(info)
                case("InterEnergy")
                   CALL CheckjobMethodConsistency((/PHASE_SHIFTS_kw,SCATT_LENGTH_kw,BOUND_STATES_kw/))
                case("InterVelocity")
                   CALL CheckjobMethodConsistency((/PHASE_SHIFTS_kw,SCATT_LENGTH_kw/))
                end select

                CALL CheckAlreadyPassed(ASSOCIATED(Job_Data%job%energy),"Energy/Velocity")

                BACKSPACE(uni_inp)
                READ(uni_inp,*,Iostat=ios_var) info,ini,fim,delta,Job_Data%job%energy_uni
                CALL CheckIostatVar()
                CALL CheckGreaterThanZero(float=ini)
                CALL CheckGreaterThanZero(float=fim)

                SELECT CASE(info)

                CASE("InterEnergy")
                   aux_real = Unit_Converter(Job_Data%job%energy_uni,at_uni_ener)
                   CALL Check_Unit(Job_Data%job%energy_uni,"energy")

                CASE("InterVelocity")
                   aux_real = Unit_Converter(Job_Data%job%energy_uni,at_uni_vel)
                   CALL Check_Unit(Job_Data%job%energy_uni,"velocity")

                CASE DEFAULT

                   CALL Non_Consistent_Info("InterEnergy/InterVelocity",TRIM(Job_Data%complete_job_name))

                END SELECT

                IF((fim-ini)*delta.LT.0)THEN
                   error_msg = "The variation of the Energy/Velocity should be consistent with the initial and final values, in"&
                        &//TRIM(Job_Data%complete_job_name)
                   CALL Error(0)
                END IF

                IF(ABS(delta).GT.ABS(fim-ini))THEN
                   error_msg = "Variation of Energy/Velocity greater than interval in "//TRIM(Job_Data%complete_job_name)
                   CALL Error(0)
                END IF

                current_real => NULL()

                DO aux_int = 1,NINT((fim-ini)/delta)+1,1

                   IF(.NOT.ASSOCIATED(current_real))THEN
                      ALLOCATE(Job_Data%job%energy)
                      current_real=>Job_Data%job%energy
                   ELSE
                      ALLOCATE(current_real%next)
                      current_real=>current_real%next
                   END IF

                   current_real%value = ini + (aux_int-1)*delta

                END DO

                ! INTEGRATION KEYWORDS:

                ! Method for numerical integration
             CASE("IntMeth")

                CALL CheckjobMethodConsistency((/PHASE_SHIFTS_kw,SCATT_LENGTH_kw,BOUND_STATES_kw/))
                CALL CheckAlreadyPassed(Job_Data%job%int_method_def)

                BACKSPACE(uni_inp)
                READ(uni_inp,*) info,Job_Data%job%int_method

                IF(Job_Data%job%int_method.NE."NUMEROV".AND.Job_Data%job%int_method.NE."JOHNSON".AND.&
                     Job_Data%job%int_method.NE."ABM".AND.Job_Data%job%int_method.NE."RK4")THEN
                   error_msg = "Unknown method for numerical integration, "//TRIM(Job_Data%job%int_method)//", in "&
                        //TRIM(Job_Data%complete_job_name)
                   CALL Error(0)
                END IF

                Job_Data%job%int_method_def = .TRUE.

                ! Step size, in the same dimension of potential
             CASE("IntStep")

                CALL CheckjobMethodConsistency((/PHASE_SHIFTS_kw,SCATT_LENGTH_kw,BOUND_STATES_kw/))
                CALL CheckAlreadyPassed(Job_Data%job%int_step_def)

                BACKSPACE(uni_inp)
                READ(uni_inp,*,Iostat=ios_var) info,Job_Data%job%int_step
                CALL CheckIostatVar()
                !          CALL CheckGreaterThanZero(float=Job_Data%job%int_step)

                Job_Data%job%int_step_def = .TRUE.

                ! Print the steps in the numerical integration, at each Job_Data%job%print_at steps
             CASE("Print")

                CALL CheckjobMethodConsistency((/PHASE_SHIFTS_kw,SCATT_LENGTH_kw,BOUND_STATES_kw/))
                CALL CheckAlreadyPassed(Job_Data%job%print_at_def)

                BACKSPACE(uni_inp)
                READ(uni_inp,*,Iostat=ios_var) info,Job_Data%job%print_at
                CALL CheckIostatVar()

                Job_Data%job%print_at_def = .TRUE.

                ! scale for wave function
             CASE("WFScale")

                CALL CheckjobMethodConsistency((/PHASE_SHIFTS_kw,SCATT_LENGTH_kw,BOUND_STATES_kw/))
                CALL CheckAlreadyPassed(Job_Data%job%wf_scale_def)

                BACKSPACE(uni_inp)
                READ(uni_inp,*,Iostat=ios_var) info,Job_Data%job%wf_scale
                CALL CheckIostatVar()

                Job_Data%job%wf_scale_def = .TRUE.

                ! Potential name for this calculation
             CASE("PotName")

                CALL CheckAlreadyPassed(Job_Data%job%pot_name_def)

                BACKSPACE(uni_inp)
                READ(uni_inp,*) info,Job_Data%job%pot_name

                Job_Data%job%pot_name_def = .TRUE.

                ! Free particle initial condition
             CASE("FPIniCond")

                CALL CheckjobMethodConsistency((/PHASE_SHIFTS_kw,SCATT_LENGTH_kw,BOUND_STATES_kw/))
                CALL CheckAlreadyPassed(Job_Data%job%ini_cond_def)

                Job_Data%job%ini_cond_free_part = .TRUE.
                Job_Data%job%ini_cond_def = .TRUE.

                ! Print chosen points to be used in phase shift/vabrational calculation in .log file
             CASE("SaveSteps")

                CALL CheckjobMethodConsistency((/PHASE_SHIFTS_kw,SCATT_LENGTH_kw,BOUND_STATES_kw/))
                CALL CheckAlreadyPassed(Job_Data%job%save_steps_def)

                Job_Data%job%save_steps = .TRUE.
                Job_Data%job%save_steps_def = .TRUE.

                ! PHASE SHIFT KEYWORDS:

                ! Method for calculation of phase shift
             CASE("PhaseShiftMeth")

                CALL CheckjobMethodConsistency((/PHASE_SHIFTS_kw/))
                CALL CheckAlreadyPassed(Job_Data%job%phase_shift_method_def)

                BACKSPACE(uni_inp)
                READ(uni_inp,*,Iostat=ios_var) info,Job_Data%job%phase_shift_method

                IF(Job_Data%job%phase_shift_method.NE."REGRESSION")THEN
                   error_msg = "Unknown method for phase shift, "//TRIM(Job_Data%job%phase_shift_method)//", in "&
                        //TRIM(Job_Data%complete_job_name)
                   CALL Error(0)
                END IF

                Job_Data%job%phase_shift_method_def = .TRUE.

                ! Step size, in the same dimension of potential, for the final of integration
             CASE("FinalIntStep")

                CALL CheckjobMethodConsistency((/PHASE_SHIFTS_kw,SCATT_LENGTH_kw/))
                CALL CheckAlreadyPassed(Job_Data%job%final_int_step_def)

                BACKSPACE(uni_inp)
                READ(uni_inp,*,Iostat=ios_var) info,Job_Data%job%final_int_step
                CALL CheckIostatVar()
                CALL CheckGreaterThanZero(float=Job_Data%job%final_int_step)

                Job_Data%job%final_int_step_def = .TRUE.

                ! Initial distance for integration, in the same dimension of potential
             CASE("Rmin")

                CALL CheckjobMethodConsistency((/PHASE_SHIFTS_kw, SCATT_LENGTH_kw, BOUND_STATES_kw/))
                CALL CheckAlreadyPassed(Job_Data%job%r_min_def)

                BACKSPACE(uni_inp)
                READ(uni_inp,*,Iostat=ios_var) info,Job_Data%job%r_min
                CALL CheckIostatVar()
                CALL CheckGreaterThanZero(float=Job_Data%job%r_min)

                Job_Data%job%r_min_def = .TRUE.

                ! Final distance for integration, in the same dimension of potential, just before increase the step size
             CASE("Rmax")

                CALL CheckjobMethodConsistency((/PHASE_SHIFTS_kw, SCATT_LENGTH_kw, BOUND_STATES_kw/))
                CALL CheckAlreadyPassed(Job_Data%job%r_max_def)

                BACKSPACE(uni_inp)
                READ(uni_inp,*,Iostat=ios_var) info,Job_Data%job%r_max
                CALL CheckIostatVar()
                CALL CheckGreaterThanZero(float=Job_Data%job%r_max)

                Job_Data%job%r_max_def = .TRUE.

                ! Number of additional iterations to get integrations steps for subsequent calculations
             CASE("AddSteps")

                CALL CheckjobMethodConsistency((/PHASE_SHIFTS_kw,SCATT_LENGTH_kw/))
                CALL CheckAlreadyPassed(Job_Data%job%add_steps_def)

                BACKSPACE(uni_inp)
                READ(uni_inp,*,Iostat=ios_var) info,Job_Data%job%add_steps
                CALL CheckIostatVar()
                CALL CheckGreaterThanZero(int=Job_Data%job%add_steps)

                Job_Data%job%add_steps_def = .TRUE.

                ! Save points at each Job_Data%job%save_at steps
             CASE("SaveAtEach")

                CALL CheckjobMethodConsistency((/PHASE_SHIFTS_kw,SCATT_LENGTH_kw/))
                CALL CheckAlreadyPassed(Job_Data%job%save_at_def)

                BACKSPACE(uni_inp)
                READ(uni_inp,*,iostat=ios_var) info,Job_Data%job%save_at
                CALL CheckIostatVar()
                CALL CheckGreaterThanZero(int=Job_Data%job%save_at)

                Job_Data%job%save_at_def = .TRUE.

                ! Print the phase shifts in log file
             CASE("PrintPhaseShifts")

                CALL CheckjobMethodConsistency((/PHASE_SHIFTS_kw/))
                CALL CheckAlreadyPassed(Job_Data%job%print_ps_def)

                Job_Data%job%print_ps = .TRUE.
                Job_Data%job%print_ps_def = .TRUE.

                ! Use potential vector generated in first integration
             CASE("SaveGlobalPot")

                CALL CheckjobMethodConsistency((/PHASE_SHIFTS_kw,BOUND_STATES_kw/))
                CALL CheckAlreadyPassed(Job_Data%job%use_potential_vector_def)

                Job_Data%job%use_potential_vector = .TRUE.
                Job_Data%job%use_potential_vector_def = .TRUE.

                ! SCATTERING LENGHT KEYWORDS:

                ! Method for calculation of scattering lenght
             CASE("ScattLenMeth")

                CALL CheckjobMethodConsistency((/SCATT_LENGTH_kw/))
                CALL CheckAlreadyPassed(Job_Data%job%scatt_len_method_def)

                BACKSPACE(uni_inp)
                READ(uni_inp,*) info,Job_Data%job%scatt_len_method

                IF(Job_Data%job%scatt_len_method.NE."DIRECT".AND.Job_Data%job%scatt_len_method.NE."MESHKOV")THEN
                   error_msg = "Unknown method for scattering length, "//TRIM(Job_Data%job%scatt_len_method)//", in "&
                        //TRIM(Job_Data%complete_job_name)
                   CALL Error(0)
                END IF

                Job_Data%job%scatt_len_method_def = .TRUE.

                ! Meshkov parameters
             CASE ("MeshkovParam")

                CALL CheckjobMethodConsistency((/SCATT_LENGTH_kw/))
                CALL CheckAlreadyPassed(Job_Data%job%meshkov_param_def)

                BACKSPACE(uni_inp)
                READ(uni_inp,*,Iostat=ios_var) info, Job_Data%job%r_bar_meshkov, Job_Data%job%beta_meshkov
                CALL CheckIostatVar()
                CALL CheckGreaterThanZero(float=Job_Data%job%r_bar_meshkov)
                CALL CheckGreaterThanZero(float=Job_Data%job%beta_meshkov)

                Job_Data%job%meshkov_param_def = .TRUE.

                ! Richardson Extrapolation
             CASE ("RichardsonExtr")

                CALL CheckjobMethodConsistency((/SCATT_LENGTH_kw/))
                CALL CheckAlreadyPassed(Job_Data%job%Rich_extr_def)

                Job_Data%job%Rich_extr=.TRUE.
                Job_Data%job%Rich_extr_def = .TRUE.

                ! VIBRATIONAL KEYWORDS

                ! Vibrational levels
             CASE("VibLevel")

                CALL CheckjobMethodConsistency((/BOUND_STATES_kw/))
                CALL CheckAlreadyPassed(ASSOCIATED(Job_Data%job%vib_lev))

                current_int => NULL()

                aux_int2 = -1
                DO

                   READ(uni_inp,*,Iostat=ios_var) aux_int
                   IF(ios_var.NE.0) EXIT
                   CALL CheckNegativeNumber(int=aux_int)
                   IF(aux_int2.GT.aux_int)THEN
                      error_msg = "Pass the vibrational levels in crescent order."
                      CALL Error(0)
                   END IF
                   aux_int2 = aux_int

                   IF(.NOT.ASSOCIATED(current_int))THEN
                      ALLOCATE(Job_Data%job%vib_lev)
                      current_int=>Job_Data%job%vib_lev
                   ELSE
                      ALLOCATE(current_int%next)
                      current_int=>current_int%next
                   END IF

                   current_int%value = aux_int

                END DO

                CALL CheckPassedValue(ASSOCIATED(Job_Data%job%vib_lev))

                BACKSPACE(uni_inp)

                ! max value of vibrational levels
             CASE("VibLevelMax")

                CALL CheckjobMethodConsistency((/BOUND_STATES_kw/))
                CALL CheckAlreadyPassed(ASSOCIATED(Job_Data%job%vib_lev),"VibLevel")

                BACKSPACE(uni_inp)
                READ(uni_inp,*,Iostat=ios_var) info,aux_int
                CALL CheckIostatVar()
                CALL CheckNegativeNumber(int=aux_int)

                current_int => NULL()

                DO aux_int2=0,aux_int,1

                   IF(.NOT.ASSOCIATED(current_int))THEN
                      ALLOCATE(Job_Data%job%vib_lev)
                      current_int=>Job_Data%job%vib_lev
                   ELSE
                      ALLOCATE(current_int%next)
                      current_int=>current_int%next
                   END IF

                   current_int%value = aux_int2

                END DO

                ! Print wave function at each iteration
             CASE("VibFixedMesh")

                CALL CheckjobMethodConsistency((/BOUND_STATES_kw/))
                CALL CheckAlreadyPassed(Job_Data%job%vib_fixed_mesh_points_def)

                Job_Data%job%vib_fixed_mesh_points = .TRUE.
                Job_Data%job%vib_fixed_mesh_points_def = .TRUE.

                ! Print wave function at each iteration
             CASE("PrintVibIterations")

                CALL CheckjobMethodConsistency((/BOUND_STATES_kw/))
                CALL CheckAlreadyPassed(Job_Data%job%print_vib_iterations_def)

                Job_Data%job%print_vib_iterations = .TRUE.
                Job_Data%job%print_vib_iterations_def = .TRUE.

                ! Print wave function at each iteration
             CASE("PrintVibWFIterations")

                CALL CheckjobMethodConsistency((/BOUND_STATES_kw/))
                CALL CheckAlreadyPassed(Job_Data%job%print_vib_iterations_def)

                Job_Data%job%print_vib_wf_iterations = .TRUE.
                Job_Data%job%print_vib_iterations = .TRUE.
                Job_Data%job%print_vib_wf_iterations_def = .TRUE.

                ! Vibrational maximun number of iterations in the first order knots counter algorithm
             CASE("Vib1stOrderMaxIter")

                CALL CheckjobMethodConsistency((/BOUND_STATES_kw/))
                CALL CheckAlreadyPassed(Job_Data%job%vib_1st_order_max_iter_def)

                BACKSPACE(uni_inp)
                READ(uni_inp,*,Iostat=ios_var) info,Job_Data%job%vib_1st_order_max_iter
                CALL CheckIostatVar()
                CALL CheckGreaterThanZero(int=Job_Data%job%vib_1st_order_max_iter)

                Job_Data%job%vib_1st_order_max_iter_def = .TRUE.

                ! Vibrational maximun number of iterations
             CASE("VibMaxIter")

                CALL CheckjobMethodConsistency((/BOUND_STATES_kw/))
                CALL CheckAlreadyPassed(Job_Data%job%vib_max_iter_def)

                BACKSPACE(uni_inp)
                READ(uni_inp,*,Iostat=ios_var) info,Job_Data%job%vib_max_iter
                CALL CheckIostatVar()
                CALL CheckGreaterThanZero(int=Job_Data%job%vib_max_iter)

                Job_Data%job%vib_max_iter_def = .TRUE.

                ! Vibrational convergence threshold
             CASE("VibThreshold")

                CALL CheckjobMethodConsistency((/BOUND_STATES_kw/))
                CALL CheckAlreadyPassed(Job_Data%job%vib_threshold_def)

                BACKSPACE(uni_inp)
                READ(uni_inp,*,Iostat=ios_var) info,Job_Data%job%vib_threshold
                CALL CheckIostatVar()
                CALL CheckGreaterThanZero(float=Job_Data%job%vib_threshold)

                Job_Data%job%vib_threshold_def = .TRUE.

                ! Vibrational point shifts toward non classical region
             CASE("VibRiniShift")

                CALL CheckjobMethodConsistency((/BOUND_STATES_kw/))
                CALL CheckAlreadyPassed(Job_Data%job%ini_point_shift_def)

                BACKSPACE(uni_inp)
                READ(uni_inp,*,Iostat=ios_var) info,Job_Data%job%internal_ini_point_shift,Job_Data%job%external_ini_point_shift
                CALL CheckIostatVar()
                CALL CheckGreaterThanZero(float=Job_Data%job%external_ini_point_shift)
                CALL CheckGreaterThanZero(float=Job_Data%job%internal_ini_point_shift)

                Job_Data%job%ini_point_shift_def = .TRUE.

                ! ANALYSIS KEYWORDS

                ! Analysis type
             CASE("AnalysisType")

                CALL CheckjobMethodConsistency((/RES_ANALYSIS_kw/))
                CALL CheckAlreadyPassed(Job_Data%analysis%analysis_type_def)

                BACKSPACE(uni_inp)
                READ(uni_inp,*,Iostat=ios_var) info,Job_Data%analysis%analysis_type

                Job_Data%analysis%analysis_type_def = .TRUE.


                ! List of values for angle
             CASE("Angle")

                CALL CheckjobMethodConsistency((/RES_ANALYSIS_kw/))
                CALL CheckAlreadyPassed(ASSOCIATED(Job_Data%analysis%angle))

                current_real => Job_Data%analysis%angle

                current_real => NULL()

                DO

                   READ(uni_inp,*,Iostat=ios_var) aux_real
                   IF(ios_var.NE.0) EXIT
                   CALL CheckNegativeNumber(float=aux_real)

                   IF(.NOT.ASSOCIATED(current_real))THEN
                      ALLOCATE(Job_Data%analysis%angle)
                      current_real=>Job_Data%analysis%angle
                   ELSE
                      ALLOCATE(current_real%next)
                      current_real=>current_real%next
                   END IF

                   current_real%value = aux_real

                END DO

                CALL CheckPassedValue(ASSOCIATED(Job_Data%analysis%angle))

                BACKSPACE(uni_inp)

                ! Interval for angle, equally spaced
             CASE("InterAngle")

                CALL CheckjobMethodConsistency((/RES_ANALYSIS_kw/))
                CALL CheckAlreadyPassed(ASSOCIATED(Job_Data%analysis%angle),"Angle")

                BACKSPACE(uni_inp)
                READ(uni_inp,*,Iostat=ios_var) info,ini,fim,delta
                CALL CheckIostatVar()
                CALL CheckNegativeNumber(float=ini)
                CALL CheckNegativeNumber(float=fim)

                IF((fim-ini)*delta.LT.0)THEN
                   error_msg = "The variation of the angle should be consistent with the initial and final values, in "&
                        &//TRIM(Job_Data%complete_job_name)
                   CALL Error(0)
                END IF

                IF(ABS(delta).GT.ABS(fim-ini))THEN
                   error_msg = "Variation of angle greater than interval in "//TRIM(Job_Data%complete_job_name)
                   CALL Error(0)
                END IF

                current_real => NULL()

                DO aux_int = 1,NINT((fim-ini)/delta)+1,1

                   IF(.NOT.ASSOCIATED(current_real))THEN
                      ALLOCATE(Job_Data%analysis%angle)
                      current_real=>Job_Data%analysis%angle
                   ELSE
                      ALLOCATE(current_real%next)
                      current_real=>current_real%next
                   END IF

                   current_real%value = ini + (aux_int-1)*delta

                END DO

                ! Job used in analysis
             CASE("Job")

                CALL CheckjobMethodConsistency((/RES_ANALYSIS_kw/))
                CALL CheckAlreadyPassed(Job_Data%analysis%job_name_2_def)

                BACKSPACE(uni_inp)
                if(Job_Data%analysis%job_name_1_def)then
                   READ(uni_inp,FMT='(4X,A100)',Iostat=ios_var) Job_Data%analysis%job_name_2
                   Job_Data%analysis%job_name_2_def = .TRUE.
                else
                   READ(uni_inp,FMT='(4X,A100)',Iostat=ios_var) Job_Data%analysis%job_name_1
                   Job_Data%analysis%job_name_1_def = .TRUE.
                End if

                ! Jobs used in analysis - WARNING: names with space are not correct read
             CASE("Jobs")

                CALL CheckjobMethodConsistency((/RES_ANALYSIS_kw/))
                CALL CheckAlreadyPassed(Job_Data%analysis%job_name_1_def)
                CALL CheckAlreadyPassed(Job_Data%analysis%job_name_2_def)

                BACKSPACE(uni_inp)
                READ(uni_inp,*,Iostat=ios_var) info, Job_Data%analysis%job_name_1, Job_Data%analysis%job_name_2

                Job_Data%analysis%job_name_1_def = .TRUE.
                Job_Data%analysis%job_name_2_def = .TRUE.

                ! Electric Dipole, either expected or transition value
             CASE("ElecDipoleFunction")

                CALL CheckjobMethodConsistency((/RES_ANALYSIS_kw/))
                CALL CheckAlreadyPassed(Job_Data%analysis%dipole_moment_name_def)

                BACKSPACE(uni_inp)
                READ(uni_inp,*,Iostat=ios_var) info,Job_Data%analysis%dipole_moment_name

                Job_Data%analysis%dipole_moment_name_def = .TRUE.

                ! Use Frank Condon approximation
             CASE("FrankCondon")

                CALL CheckjobMethodConsistency((/RES_ANALYSIS_kw/))
                CALL CheckAlreadyPassed(Job_Data%analysis%Frank_Condon_def)

                BACKSPACE(uni_inp)
                READ(uni_inp,*,Iostat=ios_var) info,Job_Data%analysis%ave_elec_dip_int,Job_Data%analysis%ave_elec_dip_int_uni

                Job_Data%analysis%Frank_Condon = .TRUE.
                Job_Data%analysis%Frank_Condon_def = .TRUE.

                ! The T0 - difference of energy between the two minima
             CASE("ElecTransEnergy")

                CALL CheckjobMethodConsistency((/RES_ANALYSIS_kw/))
                CALL CheckAlreadyPassed(Job_Data%analysis%elec_trans_ener_def)

                BACKSPACE(uni_inp)
                READ(uni_inp,*,Iostat=ios_var) info,Job_Data%analysis%elec_trans_ener,Job_Data%analysis%elec_trans_ener_uni

                aux_real = Unit_Converter(Job_Data%analysis%elec_trans_ener_uni,at_uni_ener)
                CALL Check_Unit(Job_Data%analysis%elec_trans_ener_uni,"energy")

                Job_Data%analysis%elec_trans_ener_def = .TRUE.


                ! Differential cross sectio in the center of mass system of reference
             CASE("CenterOfMassSys")

                CALL CheckjobMethodConsistency((/RES_ANALYSIS_kw/))
                CALL CheckAlreadyPassed(Job_Data%analysis%is_center_of_mass_ref_def)

                Job_Data%analysis%is_center_of_mass_ref = .TRUE.
                Job_Data%analysis%is_center_of_mass_ref_def = .TRUE.

             CASE("DunhanNumberOfParam")

                CALL CheckjobMethodConsistency((/RES_ANALYSIS_kw/))
                CALL CheckAlreadyPassed(Job_Data%analysis%Dunhan_def)

                BACKSPACE(uni_inp)
                READ(uni_inp,*,Iostat=ios_var) info,Job_Data%analysis%Dunhan_n_vib,Job_Data%analysis%Dunhan_n_rot,&
                     Job_Data%analysis%Dunhan_n_vib_rot
                IF(ios_var.NE.0)THEN
                   error_msg = "Error in reading DunhanNumberOfParam."
                   CALL Error(0)
                END IF

                Job_Data%analysis%Dunhan_def = .TRUE.

             CASE("DunhanPureVibrational")

                CALL CheckjobMethodConsistency((/RES_ANALYSIS_kw/))
                CALL CheckAlreadyPassed(Job_Data%analysis%Dunhan_def)

                BACKSPACE(uni_inp)
                READ(uni_inp,*,Iostat=ios_var) info,Job_Data%analysis%Dunhan_n_vib
                IF(ios_var.NE.0)THEN
                   error_msg = "Error in reading DunhanPureVibrational."
                   CALL Error(0)
                END IF
                Job_Data%analysis%Dunhan_n_rot = 0
                Job_Data%analysis%Dunhan_n_vib_rot = 0

                Job_Data%analysis%Dunhan_def = .TRUE.

             CASE("DunhanPureRotational")

                CALL CheckjobMethodConsistency((/RES_ANALYSIS_kw/))
                CALL CheckAlreadyPassed(Job_Data%analysis%Dunhan_def)

                BACKSPACE(uni_inp)
                READ(uni_inp,*,Iostat=ios_var) info,Job_Data%analysis%Dunhan_n_rot
                IF(ios_var.NE.0)THEN
                   error_msg = "Error in reading DunhanPureRotational."
                   CALL Error(0)
                END IF
                Job_Data%analysis%Dunhan_n_vib = 0
                Job_Data%analysis%Dunhan_n_vib_rot = 0

                Job_Data%analysis%Dunhan_def = .TRUE.

                ! Dimension for print cross section
             CASE("CrossSecUni")

                CALL CheckjobMethodConsistency((/RES_ANALYSIS_kw,SCATT_LENGTH_kw/))

                BACKSPACE(uni_inp)
                READ(uni_inp,*,Iostat=ios_var) info,cross_sec_uni_tmp

                aux_real = Unit_Converter(cross_sec_uni_tmp,at_uni_dist)
                CALL Check_Unit(cross_sec_uni_tmp,"distance")

                IF(ASSOCIATED(Job_Data%analysis))THEN
                   CALL CheckAlreadyPassed(Job_Data%analysis%cross_sec_uni_def)
                   Job_Data%analysis%cross_sec_uni = cross_sec_uni_tmp
                   Job_Data%analysis%cross_sec_uni_def = .TRUE.
                ELSEIF(ASSOCIATED(Job_Data%job))THEN
                   CALL CheckAlreadyPassed(Job_Data%job%cross_sec_uni_def)
                   Job_Data%job%cross_sec_uni = cross_sec_uni_tmp
                   Job_Data%job%cross_sec_uni_def = .TRUE.
                ELSE
                   error_msg = "Neither analysis nor job associated in CrossSecUni reading."
                   CALL Error(2)
                END IF

             CASE DEFAULT
                BACKSPACE(uni_inp)
                EXIT

             END SELECT

          END DO

       CASE DEFAULT

          error_msg = "Inconsistecy of method - "//TRIM(Job_Data%job_method)
          CALL Error(2)
          
       END SELECT

    CASE DEFAULT
       
       IF(first_read)THEN
          error_msg = "Unknow expression: "//TRIM(info)
       ELSE
          error_msg = "Unknow expression: "//TRIM(info)//", in "//TRIM(Job_Data%complete_job_name)
       END IF
       CALL Error(0)
       
    END SELECT
    
    first_read=.FALSE.

  CONTAINS
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
    SUBROUTINE Set_default_formats()

      IMPLICIT NONE
      
      formats%energy = "ES15.8"
      formats%cross_section = "ES15.8"
      formats%angle = "F0.10"
      formats%ang_momentum = "I0"
      formats%distance = "F0.6"
      formats%wave_function = "ES15.8"
      formats%potential = "ES17.10"
      formats%energy_cm1 = "F0.5"
      formats%einstein_coef = "F8.5"
      formats%dip_mom = "F8.5"
      formats%phase_shift = "ES15.8"
      formats%dunhan = "F0.7"
      formats%vib_level = "I2"
      
    END SUBROUTINE Set_default_formats
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
    SUBROUTINE CheckAlreadyPassed(def,new_info)

      IMPLICIT NONE

      LOGICAL, INTENT(in) :: def
      CHARACTER(LEN=*), INTENT(in), OPTIONAL :: new_info
      
      IF(def)THEN
         IF(PRESENT(new_info)) THEN
            error_msg = TRIM(new_info)//" already passed in "//TRIM(Job_Data%complete_job_name)
         ELSE
            error_msg = TRIM(info)//" already passed in "//TRIM(Job_Data%complete_job_name)
         END IF
         CALL Error(0)
      END IF

    END SUBROUTINE CheckAlreadyPassed
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
    SUBROUTINE CheckPassedValue(is_assoc)

      IMPLICIT NONE

      LOGICAL, INTENT(in) :: is_assoc

      IF(.NOT.is_assoc)THEN
         error_msg = "No value for "//TRIM(info)//" passed, in "//TRIM(Job_Data%complete_job_name)
         CALL Error(0)
      END IF

    END SUBROUTINE CheckPassedValue
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
    SUBROUTINE CheckGreaterThanZero(int,float)

      IMPLICIT NONE

      INTEGER, INTENT(in), OPTIONAL :: int
      REAL(kind = real_kind), INTENT(in), OPTIONAL :: float
      
      IF((PRESENT(int).AND.Int.LE.0).OR.(PRESENT(float).AND.float.LE.0))THEN
         error_msg = "The value of "//TRIM(info)//" should be positive, in "//TRIM(Job_Data%complete_job_name)
         CALL Error(0)
      END IF

    END SUBROUTINE CheckGreaterThanZero
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
    SUBROUTINE CheckNegativeNumber(int,float)

      IMPLICIT NONE

      INTEGER, INTENT(in), OPTIONAL :: int
      REAL(kind = real_kind), INTENT(in), OPTIONAL :: float
      
      IF((PRESENT(int).AND.Int.LT.0).OR.(PRESENT(float).AND.float.LT.0))THEN
         error_msg = "Found negative value for "//TRIM(info)//", in "//TRIM(Job_Data%complete_job_name)
         CALL Error(0)
      END IF

    END SUBROUTINE CheckNegativeNumber
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
    SUBROUTINE CheckIostatVar()

      IMPLICIT NONE

      IF(ios_var.NE.0)THEN
         error_msg = "Error in reading "//TRIM(info)//" in "//TRIM(Job_Data%complete_job_name)
         CALL Error(0)
      END IF

    END SUBROUTINE CheckIostatVar
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
    SUBROUTINE CheckjobMethodConsistency(InconsitentMethods)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(in), DIMENSION(:) :: InconsitentMethods

      IF(.NOT.ANY(InconsitentMethods.EQ.Job_Data%job_method)) THEN
         error_msg = "Do not pass "//TRIM(info)//" for "//TRIM(Job_Data%job_method)
         CALL Error(0)
      END IF

    END SUBROUTINE CheckjobMethodConsistency

  END SUBROUTINE Get_Input
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  SUBROUTINE Print_Input_Data()

    IMPLICIT NONE

    TYPE(TypeIntegerLinkedList), POINTER :: current_int => NULL()
    TYPE(TypeRealLinkedList), POINTER :: current_real =>NULL()

    PRINT*,"Input data:"
    Job_Data=>Jobs_Data_list
    DO WHILE(ASSOCIATED(Job_Data))

       PRINT*,"----------------------------------------------------"
       PRINT*,"job_method = ",TRIM(Job_Data%job_method)
       PRINT*,"complete_job_name = ",TRIM(Job_Data%complete_job_name)
       PRINT*,"job_name = ",TRIM(Job_Data%job_name)

       IF(ASSOCIATED(Job_Data%job))THEN

          ! ANGULAR MOMENTUM VARIABLES:
          PRINT*,"ANGULAR MOMENTUM VARIABLES:"

          PRINT*,"l = "
          current_int=>Job_Data%job%l
          DO WHILE(ASSOCIATED(current_int))
             PRINT*,current_int%value
             current_int=>current_int%next
          END DO

          ! ENERGY VARIABLES:
          PRINT*,"ENERGY VARIABLES:"

          PRINT*,"energy = "
          current_real=>Job_Data%job%energy
          DO WHILE(ASSOCIATED(current_real))
             PRINT*,current_real%value
             current_real=>current_real%next
          END DO
          PRINT*,"energy_uni = ",Job_Data%job%energy_uni

          ! INTEGRATION VARIABLES:
          PRINT*,"INTEGRATION VARIABLES:"

          PRINT*,"int_method = ",TRIM(Job_Data%job%int_method)
          PRINT*,"int_method_def = ",Job_Data%job%int_method_def

          PRINT*,"int_step = ",Job_Data%job%int_step
          PRINT*,"int_step_def = ",Job_Data%job%int_step_def

          PRINT*,"print_at = ",Job_Data%job%print_at
          PRINT*,"print_at_def = ",Job_Data%job%print_at_def

          PRINT*,"wf_scale = ",Job_Data%job%wf_scale
          PRINT*,"wf_scale_def = ",Job_Data%job%wf_scale_def
          
          PRINT*,"pot_name = ",TRIM(Job_Data%job%pot_name)
          PRINT*,"pot_name_def = ",Job_Data%job%pot_name_def

          PRINT*,"ini_cond_free_part = ",Job_Data%job%ini_cond_free_part
          PRINT*,"ini_cond_def = ",Job_Data%job%ini_cond_def

          PRINT*,"save_steps = ",Job_Data%job%save_steps
          PRINT*,"save_steps_def = ",Job_Data%job%save_steps_def

          ! PHASE SHIFT KEYWORDS:
          PRINT*,"PHASE SHIFT KEYWORDS:"

          PRINT*,"phase_shift_method = ",TRIM(Job_Data%job%phase_shift_method)
          PRINT*,"phase_shift_method_def = ",Job_Data%job%phase_shift_method_def

          PRINT*,"final_int_step = ",Job_Data%job%final_int_step
          PRINT*,"final_int_step_def = ",Job_Data%job%final_int_step_def

          PRINT*,"r_min = ",Job_Data%job%r_min
          PRINT*,"r_min_def = ",Job_Data%job%r_min_def
          PRINT*,"r_min_default_explanation = ",Job_Data%job%r_min_default_explanation

          PRINT*,"r_max = ",Job_Data%job%r_max
          PRINT*,"r_max_def = ",Job_Data%job%r_max_def
          PRINT*,"r_max_default_explanation = ",Job_Data%job%r_min_default_explanation

          PRINT*,"add_steps = ",Job_Data%job%add_steps
          PRINT*,"add_steps_def = ",Job_Data%job%add_steps_def

          PRINT*,"save_at = ",Job_Data%job%save_at
          PRINT*,"save_at_def = ",Job_Data%job%save_at_def

          PRINT*,"print_ps = ",Job_Data%job%print_ps
          PRINT*,"print_ps_def = ",Job_Data%job%print_ps_def

          PRINT*,"use_potential_vector = ",Job_Data%job%use_potential_vector
          PRINT*,"use_potential_vector_def = ",Job_Data%job%use_potential_vector_def

          ! SCATTERING LENGHT VARIABLES:
          PRINT*,"SCATTERING LENGHT VARIABLES:"

          PRINT*,"scatt_len_method = ",TRIM(Job_Data%job%scatt_len_method)
          PRINT*,"scatt_len_method_def = ",Job_Data%job%scatt_len_method_def

          PRINT*,"r_bar_meshkov = ",Job_Data%job%r_bar_meshkov
          PRINT*,"beta_meshkov", Job_Data%job%beta_meshkov
          PRINT*,"meshkov_param_def = ",Job_Data%job%meshkov_param_def
          
          PRINT*,"Rich_extr = ",Job_Data%job%Rich_extr
          PRINT*,"Rich_extr_def = ",Job_Data%job%Rich_extr_def

          ! VIBRATIONAL VARIABLES:
          PRINT*,"VIBRATIONAL VARIABLES:"
          
          PRINT*,"VibLevel = "
          current_int=>Job_Data%job%vib_lev
          DO WHILE(ASSOCIATED(current_int))
             PRINT*,current_int%value
             current_int=>current_int%next
          END DO

          PRINT*,"vib_fixed_mesh_points = ",Job_Data%job%vib_fixed_mesh_points
          PRINT*,"vib_fixed_mesh_points_def = ",Job_Data%job%vib_fixed_mesh_points_def

          PRINT*,"print_vib_iterations = ",Job_Data%job%print_vib_iterations
          PRINT*,"print_vib_iterations_def = ",Job_Data%job%print_vib_iterations_def

          PRINT*,"vib_1st_order_max_iter = ",Job_Data%job%vib_1st_order_max_iter
          PRINT*,"vib_1st_order_max_iter_def = ",Job_Data%job%vib_1st_order_max_iter_def

          PRINT*,"vib_max_iter_iter = ",Job_Data%job%vib_max_iter
          PRINT*,"vib_max_iter_def = ",Job_Data%job%vib_max_iter_def

          PRINT*,"vib_threshold = ",Job_Data%job%vib_threshold
          PRINT*,"vib_threshold_def = ",Job_Data%job%vib_threshold_def

          PRINT*,"internal_ini_point_shift = ",Job_Data%job%internal_ini_point_shift
          PRINT*,"external_ini_point_shift = ",Job_Data%job%external_ini_point_shift
          PRINT*,"ini_point_shift_def = ",Job_Data%job%ini_point_shift_def
          
       END IF

       IF(ASSOCIATED(Job_Data%analysis))THEN

          PRINT*,"analysis_type = ",Job_Data%analysis%analysis_type
          PRINT*,"analysis_type_def = ",Job_Data%analysis%analysis_type_def

          PRINT*,"job_name_1 = ",Job_Data%analysis%job_name_1
          PRINT*,"job_name_1_def = ",Job_Data%analysis%job_name_1_def

          PRINT*,"job_name_2 = ",Job_Data%analysis%job_name_2
          PRINT*,"job_name_2_def = ",Job_Data%analysis%job_name_2_def

          PRINT*,"dipole_moment_name = ",Job_Data%analysis%dipole_moment_name
          PRINT*,"dipole_moment_name_def = ",Job_Data%analysis%dipole_moment_name_def

          PRINT*,"Frank_Condon = ",Job_Data%analysis%dipole_moment_name
          PRINT*,"Frank_Condon_def = ",Job_Data%analysis%dipole_moment_name_def

          PRINT*,"ave_elec_dip_int = ",Job_Data%analysis%ave_elec_dip_int
          PRINT*,"ave_elec_dip_int_uni = ",Job_Data%analysis%ave_elec_dip_int_uni

          PRINT*,"angle = "
          current_real=>Job_Data%analysis%angle
          DO WHILE(ASSOCIATED(current_real))
             PRINT*,current_real%value
             current_real=>current_real%next
          END DO

          PRINT*,"Dunhan_n_vib = ",Job_Data%analysis%Dunhan_n_vib
          PRINT*,"Dunhan_n_rot = ",Job_Data%analysis%Dunhan_n_rot
          PRINT*,"Dunhan_n_vib_rot = ",Job_Data%analysis%Dunhan_n_vib_rot
          PRINT*,"Dunhan_def = ",Job_Data%analysis%Dunhan_def

          PRINT*,"is_center_of_mass_ref = ",Job_Data%analysis%is_center_of_mass_ref
          PRINT*,"is_center_of_mass_ref_def = ",Job_Data%analysis%is_center_of_mass_ref_def

          PRINT*,"cross_sec_uni = ",Job_Data%analysis%cross_sec_uni
          PRINT*,"cross_sec_uni_def = ",Job_Data%analysis%cross_sec_uni_def

       END IF
       
       PRINT*,"----------------------------------------------------"

       Job_Data => Job_Data%next
    END DO

  END SUBROUTINE Print_Input_Data
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  SUBROUTINE Check_Input_Data()

    IMPLICIT NONE

    TYPE(TypeRealLinkedList), POINTER :: current_real
    TYPE(TypeJobSchedule), POINTER :: curr_job_1=>NULL()
    TYPE(TypeJobSchedule), POINTER :: curr_job_2=>NULL()

    REAL(kind = real_kind) :: aux_real

    Job_Data => Jobs_Data_list

    DO WHILE(ASSOCIATED(Job_Data))

       select case(Job_Data%job_method)

       CASE(FORMATS_kw)

          formats => Job_Data%job_formats

       CASE(MASSES_kw,RED_MASS_kw,ATOMS_kw)

          particles => Job_Data%job_particles

       CASE(PHASE_SHIFTS_kw,SCATT_LENGTH_kw,BOUND_STATES_kw)

          if(.not.associated(particles))then
             error_msg = "Pass the particle definition before jobs."
             CALL Error(0)
          end if
          
          ! Select the potential
          CALL Select_Potential(Job_Data%job%pot_name)

          if(Job_Data%job_method.EQ.BOUND_STATES_kw.AND.Pot%eq_dist.LE.zero)then
             error_msg = "Vibrational states can be calculated only for potentials with minimum, in "&
                  //TRIM(Job_Data%complete_job_name)
             CALL Error(0)
          end if

          ! Allocation of l
          
          IF(.NOT.ASSOCIATED(Job_Data%job%l))Then
             
             IF(Job_Data%job_method.EQ.BOUND_STATES_kw.OR.Job_Data%job_method.EQ.SCATT_LENGTH_kw)THEN
                ALLOCATE(Job_Data%job%l)
                Job_Data%job%l%value = 0
             ELSEIF(Job_Data%job_method.EQ.PHASE_SHIFTS_kw)THEN
                error_msg = "Value of l not passed in "//TRIM(Job_Data%complete_job_name)
                CALL Error(0)
             ELSE
                error_msg = "Inconsistency of job, "//TRIM(Job_Data%complete_job_name)//", in checking l."
                CALL Error(2)
             END IF
             
          END IF
          
          ! Allocation of energy
          IF(.NOT.ASSOCIATED(Job_Data%job%energy))THEN
             IF(Job_Data%job_method.EQ.PHASE_SHIFTS_kw)then
                error_msg = "Value of energy/velocity not passed in "//TRIM(Job_Data%complete_job_name)
                CALL Error(0) 
             END IF
          END IF
          IF(ASSOCIATED(Job_Data%job%energy).AND..NOT.associated(Particles))THEN
             error_msg = "Define particles information for a calculations with defined energy."
             CALL Error(0)
          END IF
          IF(Job_Data%job_method.EQ.BOUND_STATES_kw.AND.associated(Job_Data%job%energy))THEN
             current_real=>Job_Data%job%energy
             DO WHILE(ASSOCIATED(current_real))
                current_real%value = current_real%value - Unit_Converter(at_uni_ener,Job_Data%job%energy_uni,Pot%diss_ener)
                current_real=>current_real%next
             END DO
          END IF
          
          ! Integration method: related only with Meshkov calculation. Checked below.
          ! CHANGE POINTS SAVED IN BOUND IF INT METHOD IS OTHER THAN NUMEROV!!!
          
          ! Integration step: change to user units
          IF(.NOT.Job_Data%job%int_step_def)THEN
             Job_Data%job%int_step = Unit_Converter(at_uni_dist,Pot%dist_uni,Job_Data%job%int_step)
          END IF
          
          ! Print information: no check needed
          
          ! Scaling for wave function: no check needed
          
          ! Initial condition - more information with Meshkov calculation. Check below
          IF((.NOT.Job_Data%job%ini_cond_def).AND.(.NOT.ASSOCIATED(Pot).OR.Pot%type.EQ."Square"))THEN
             Job_Data%job%ini_cond_free_part = .TRUE.
          END IF
          
          ! Save points: some information related with Meshkov calculation. Checked below.
          
          ! Phase shift method: no check needed          
          
          ! Final integration step: change to user units and check if is greater than or equal Job_Data%job%int_step_def
          IF(.NOT.Job_Data%job%final_int_step_def)THEN
             Job_Data%job%final_int_step = Unit_Converter(at_uni_dist,Pot%dist_uni,Job_Data%job%final_int_step)
          END IF
          IF(Job_Data%job%final_int_step.LT.Job_Data%job%int_step)THEN
             IF(Job_Data%job%final_int_step_def.AND.Job_Data%job%int_step_def)THEN
                error_msg = "Final step size should be greater or equal the step size, in "//TRIM(Job_Data%job_method)
                CALL Error(0)
             ELSE
                Job_Data%job%final_int_step = 20*Job_Data%job%int_step
             END IF
          END IF
          
          ! R min.
          aux_real=Unit_Converter(at_uni_dist,pot%dist_uni,Pot%lower_bound)
          IF(Job_Data%job%r_min_def.AND.Job_Data%job_method.EQ.BOUND_STATES_kw.AND..NOT.job_Data%job%vib_fixed_mesh_points)then
             error_msg = "Pass Rmin for "//trim(BOUND_STATES_kw)//" only in a fixed mesh points calculation."
             CALL Error(2)
          end IF
          IF(.NOT.Job_Data%job%r_min_def)THEN
             IF(.NOT.ASSOCIATED(Pot).OR.Pot%type.EQ."Square")THEN
                Job_Data%job%r_min = Job_Data%job%int_step
             ELSE
                Job_Data%job%r_min = aux_real
             END IF
          END IF
          IF(Job_Data%job%r_min.LT.aux_real) THEN
             WRITE(unit=error_msg,fmt='("R min cannot be lower than ",'//TRIM(formats%distance)//'," ",A,'//&
                  '" the lower bound for this potential, in ",A)') aux_real, pot%dist_uni, TRIM(Job_Data%complete_job_name)
             CALL Error(0)
          END IF
          
          if(pot%upper_bound.gt.zero)then
             aux_real=Unit_Converter(at_uni_dist,pot%dist_uni,Pot%upper_bound)
          else
             aux_real=10*Unit_Converter(at_uni_dist,pot%dist_uni,Pot%eq_dist)
          end if
          ! R max. Some information related with Meshkov calculation. Checked below.
          IF(Job_Data%job%r_max_def.AND.Job_Data%job_method.EQ.BOUND_STATES_kw.AND..NOT.Job_Data%job%vib_fixed_mesh_points)then
             error_msg = "Pass Rmin for "//trim(BOUND_STATES_kw)//" only in a fixed mesh points calculation."
             CALL Error(2)
          end IF
          IF(.NOT.Job_Data%job%r_max_def)THEN
             IF(.NOT.ASSOCIATED(Pot))THEN
                error_msg = "Rmax should be passed for free particle calculation"
                CALL Error(0)
             ELSE
                Job_Data%job%r_max = aux_real
             END IF
          END IF
          
          ! Additional steps. Some information related with Meshkov calculation. Checked below.
          
          ! Save points at how much iterations. Some information related with Meshkov calculation. Checked below.

          ! Print phase shifts: no check needed

          ! Use potential vector: no check needed
          IF(Job_Data%job%use_potential_vector.AND.Job_Data%job_method.EQ.BOUND_STATES_kw&
               .AND..NOT.Job_Data%job%vib_fixed_mesh_points)THEN
             error_msg = "Pass SaveGlobalPot for "//trim(BOUND_STATES_kw)//" only for a fixed mesh points calculation, in "&
                  //TRIM(Job_Data%complete_job_name)
             CALL Error(0)
          END IF

          ! Scattering length method: no check needed
          
          ! Meshkov calculation related checks
          IF(Job_Data%job_method.EQ.SCATT_LENGTH_kw.AND.Job_Data%job%scatt_len_method.EQ."MESHKOV")THEN

             ! Integration step: change to user units
             IF(Job_Data%job%int_step_def)THEN
                error_msg = "Do not pass the integration step for Meshkov procedure. Use AddSteps instead."
                CALL Error(0)
             END IF
             
             ! Integration method.             
             IF(Job_Data%job%int_method_def.AND.Job_Data%job%int_method.NE."JOHNSON")THEN
                error_msg = "Calculation of scattering length by Meshkov method require JOHNSON method of integration, in "&
                     &//TRIM(Job_Data%complete_job_name)
                CALL Error(0)
             END IF
             Job_Data%job%int_method="JOHNSON"
             
             ! R max
             IF(Job_Data%job%r_max_def) THEN
                error_msg = "Do not pass the end of integration for scattering length calculation by MESHKOV, in "&
                     &//TRIM(Job_Data%complete_job_name)
                CALL Error(0)
             END IF
             
             ! Additional steps - multiple of 4, to allow painless Richardson extrapolation
             Job_Data%job%add_steps = 4*(Job_Data%job%add_steps/4)
             
             ! Save points at how much iterations.
             Job_Data%job%save_at = -1
             
             ! Meshkov parameters
             IF(.NOT.Job_Data%job%meshkov_param_def)THEN
                error_msg = "Pass the Meshkov parameters for calculation of scattering length by MESHKOV method, in "&
                     &//TRIM(Job_Data%complete_job_name)
                CALL Error(0)
             END IF
             
             ! Initial condition.
             IF(Job_Data%job%ini_cond_def) THEN
                error_msg = "Do not define initial condition for for scattering length calculation by MESHKOV, in "&
                     &//TRIM(Job_Data%complete_job_name)
                CALL Error(0)
             END IF

          ELSE

             ! R max and Rmin
             if(Job_Data%job%r_max.lt.Job_Data%job%r_min)then
                error_msg = "R max cannot be lower than R min, in "//TRIM(Job_Data%complete_job_name)
                CALL error(0)
             end If
             
             ! Meshkov parameters
             IF(Job_Data%job%meshkov_param_def) THEN
                error_msg = "Do not pass the Meshkov parameters in a calculation other than scattering length by MESHKOV, in "&
                     &//TRIM(Job_Data%complete_job_name)
                CALL Error(0)
             END IF
             
             ! Richardson extrapolation
             IF(Job_Data%job%Rich_extr_def) THEN
                error_msg = "Do not ask for Richardson extrapolation in a calculation other than scattering length by MESHKOV, in "&
                     &//TRIM(Job_Data%complete_job_name)
                CALL Error(0)
             END IF

          END IF

          ! Vibrational levels: no check needed

          ! Print vibrational iterations: no check needed

          ! Vibrational first order maximum iterations: no check needed

          ! Vibrational maximum iterations: no check needed

          ! Vibrational thereshold - convert to hartree
          if(ASSOCIATED(Job_Data%job%energy))then
             Job_Data%job%vib_threshold = Unit_Converter(Job_Data%job%energy_uni,at_uni_ener,Job_Data%job%vib_threshold)
          else
             Job_Data%job%vib_threshold = Unit_Converter("cm-1           ",at_uni_ener,Job_Data%job%vib_threshold)
          end if
          
          ! External initial point shift: change to user units
          IF(.NOT.Job_Data%job%ini_point_shift_def)THEN
             Job_Data%job%external_ini_point_shift = &
                  Unit_Converter(at_uni_dist,Pot%dist_uni,Job_Data%job%external_ini_point_shift)
             Job_Data%job%internal_ini_point_shift = &
                  Unit_Converter(at_uni_dist,Pot%dist_uni,Job_Data%job%internal_ini_point_shift)
          END IF

          ! Internal initial point shift: no check needed

          ! Unit for print cross sections.
          if(Job_Data%job%cross_sec_uni_def.AND.Job_Data%job_method.EQ.BOUND_STATES_kw)then
             error_msg = "Do not pass cross section unit for "//TRIM(BOUND_STATES_kw)//", in "//TRIM(Job_Data%complete_job_name)
             CALL Error(0)
          end if

       CASE(RES_ANALYSIS_kw)

          if(.NOT.Job_Data%analysis%analysis_type_def)then
             error_msg = "Pass the analysis type in "//TRIM(Job_Data%complete_job_name)
             CALL Error(0)
          end if

          ! Analysis type
          SELECT CASE(Job_Data%analysis%analysis_type)

          CASE("PartialCrossSec","TotalCrossSec","DifferentialCrossSec","VibrationalIntensity",&
               "DunhanExpansion","EquilibriumDistance")

          CASE DEFAULT
             error_msg = "Unknown analysis type: "//TRIM(Job_Data%analysis%analysis_type)
             CALL Error(0)
          END SELECT

          ! Job names
          if(Job_Data%analysis%job_name_1_def)then
             curr_job_1=>Jobs_Data_list
             do while(associated(curr_job_1))
                
                if(associated(curr_job_1, target=Job_Data))then
                   error_msg = "The job "//TRIM(Job_Data%analysis%job_name_1)//" must be done before its analysis."
                   CALL Error(0)
                End if
                
                if(curr_job_1%job_name.EQ.Job_Data%analysis%job_name_1) exit
                
                curr_job_1 => curr_job_1%next
             end do
          end if

          if(Job_Data%analysis%job_name_2_def)then
             curr_job_2=>Jobs_Data_list
             do while(associated(curr_job_2))
                
                if(associated(curr_job_2,target=Job_Data))then
                   error_msg = "The job "//TRIM(Job_Data%analysis%job_name_2)//" must be done before its analysis, in "&
                        //TRIM(Job_Data%complete_job_name)
                   CALL Error(0)
                End if
                
                if(curr_job_2%job_name.EQ.Job_Data%analysis%job_name_2) exit
                
                curr_job_2 => curr_job_2%next
             end do
          end if

          ! Job method and analysis type consistency
          if(.NOT.Job_Data%analysis%job_name_1_def)THEN
             error_msg = "Pass the name of the job to be analysed, in "//TRIM(Job_Data%complete_job_name)
             CALL Error(0)
          End if
          
          SELECT CASE(Job_Data%analysis%analysis_type)
             
          CASE("PartialCrossSec","TotalCrossSec","DifferentialCrossSec")

             if(Job_Data%analysis%job_name_2_def)then
                error_msg = "Pass just one job to be analysed for "//TRIM(Job_Data%analysis%analysis_type)//&
                     ", in "//TRIM(Job_Data%complete_job_name)
                CALL Error(0)
             end if
             
             if(curr_job_1%job_method.NE.PHASE_SHIFTS_kw)then
                error_msg = "The results analysis "//TRIM(Job_Data%analysis%analysis_type)//" can be done only on "&
                     //TRIM(PHASE_SHIFTS_kw)//" method, in "//TRIM(Job_Data%complete_job_name)
                CALL Error(0)
             end if

          CASE("VibrationalIntensity","DunhanExpansion","EquilibriumDistance")

             if(curr_job_1%job_method.NE.BOUND_STATES_kw)then
                error_msg = "The results analysis "//TRIM(Job_Data%analysis%analysis_type)//" can be done only on "&
                     //TRIM(BOUND_STATES_kw)//" method."
                CALL Error(0)
             end if

             if(Job_Data%analysis%cross_sec_uni_def)then
                error_msg = "Do not pass cross section unit for "//TRIM(Job_Data%analysis%analysis_type)//&
                     ", in "//TRIM(Job_Data%complete_job_name)
                CALL Error(0)
             end if

             SELECT CASE(Job_Data%analysis%analysis_type)

             CASE("DunhanExpansion","EquilibriumDistance")

                if(Job_Data%analysis%job_name_2_def)then
                   error_msg = "Pass just one job to be analysed for "//TRIM(Job_Data%analysis%analysis_type)//&
                        ", in "//TRIM(Job_Data%complete_job_name)
                   CALL Error(0)
                end if

             END SELECT
             
          CASE DEFAULT
             
             error_msg = "Unknown analysis type: "//TRIM(Job_Data%analysis%analysis_type)
             CALL Error(0)
             
          END SELECT
          
          ! Dipole Moment Name
          If(Job_Data%analysis%dipole_moment_name_def.AND.Job_Data%analysis%analysis_type.NE."VibrationalIntensity")THEN
             error_msg = "Pass the dipole moment only for VibrationalIntensity analysis type, in "//TRIM(Job_Data%complete_job_name)
             CALL Error(0)
          END If

          ! Frank Condon
          If(Job_Data%analysis%Frank_Condon_def)THEN

             IF(Job_Data%analysis%analysis_type.NE."VibrationalIntensity")then
                error_msg = "Ask for Frank Condon approximation only for VibrationalIntensity analysis type, in "&
                     //TRIM(Job_Data%complete_job_name)
                CALL Error(0)
             End IF

             If(Job_Data%analysis%dipole_moment_name_def)THEN
                error_msg = "Pass the dipole moment or ask for Frank Condon approximation."
                CALL Error(0)
             END If

             IF(.NOT.Job_Data%analysis%job_name_2_def)then
                error_msg = "Ask for Frank Condon only when two BOUND_STATES are being analysed, in "&
                     //TRIM(Job_Data%complete_job_name)
                CALL Error(0)
             End IF

          END If

          If(Job_Data%analysis%analysis_type.EQ."VibrationalIntensity".AND.&
               .NOT.Job_Data%analysis%dipole_moment_name_def.AND.&
               .NOT.Job_Data%analysis%Frank_Condon_def)THEN
             error_msg = "Pass the dipole moment or ask for Frank Condon approximation for"&
                  //" VibrationalIntensity analysis type, in "//TRIM(Job_Data%complete_job_name)
             CALL Error(0)
          END If


          ! Dunhan parameters
          If(Job_Data%analysis%Dunhan_def)then
             if(Job_Data%analysis%analysis_type.NE."DunhanExpansion")then
                error_msg = "Pass the number of Dunhan parameters only in a DunhanExpansion analysis type, in "&
                     //TRIM(Job_Data%complete_job_name)
                CALL Error (0)
             end if
          else
             if(Job_Data%analysis%analysis_type.EQ."DunhanExpansion")then
                error_msg = "Pass the number of Dunhan parameters in a DunhanExpansion analysis type, in "&
                     //TRIM(Job_Data%complete_job_name)
                CALL Error (0)
             end if
          end If

          ! Angle
          IF(ASSOCIATED(Job_Data%analysis%angle))Then
             if(Job_Data%analysis%analysis_type.NE."DifferentialCrossSec")then
                error_msg = "Pass angle only in a DifferentialCrossSec analysis type, in "//TRIM(Job_Data%complete_job_name)
                CALL Error (0)
             end if
          else
             if(Job_Data%analysis%analysis_type.EQ."DifferentialCrossSec")then
                error_msg = "Pass angle in a DifferentialCrossSec analysis type, in "//TRIM(Job_Data%complete_job_name)
                CALL Error (0)
             end if
          end IF

          ! Center of mass reference system
          if(Job_Data%analysis%is_center_of_mass_ref_def)then
             if(Job_Data%analysis%analysis_type.NE."DifferentialCrossSec")then
                error_msg = "Center of mass reference system only makes sense in a DifferentialCrossSec analysis type, in "&
                     //TRIM(Job_Data%complete_job_name)
                CALL Error(0)
             end if
          end if

          ! Center of mass calculation of the differenial cross section requires the mass of the individual particles
          IF(Job_Data%analysis%analysis_type.EQ."DifferentialCrossSec".AND.&
               .NOT.Job_Data%analysis%is_center_of_mass_ref.AND.Particles%mass1.LT.zero)THEN
             error_msg = "Individual masses should be known for diferential cross section in laboratory system of reference, in "//&
                  TRIM(Job_Data%job_method)
             CALL Error(0)
          END IF

       CASE DEFAULT

          error_msg = "Inconsistecy of method - "//TRIM(Job_Data%job_method)
          CALL Error(2)

       END select

       Job_Data => Job_Data%next

    END DO

  END SUBROUTINE Check_Input_Data
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  SUBROUTINE Free_Input_Data_Memory()

    IMPLICIT NONE

    TYPE(TypeJobSchedule), POINTER :: next => NULL()

    Job_Data => Jobs_Data_list
    DO WHILE(ASSOCIATED(Job_Data))

       next => Job_Data%next

       IF(ASSOCIATED(Job_Data%job))THEN
          CALL Free_ILL(Job_Data%job%l)
          CALL Free_ILL(Job_Data%job%vib_lev)
          CALL Free_RLL(Job_Data%job%energy)
       END IF

       IF(ASSOCIATED(Job_Data%analysis))THEN
          CALL Free_RLL(Job_Data%analysis%angle)
       END IF

       DEALLOCATE(Job_Data)

       Job_Data => next
    END DO

  END SUBROUTINE Free_Input_Data_Memory
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
END MODULE ModInput
