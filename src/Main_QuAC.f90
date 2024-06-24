
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
! Main program of QuAC
!
! Get the input name from command line, get input and potential, check the input data and loop over the linked list Jobs_Data_list
!
PROGRAM QuAC

  USE ModUtil, ONLY: verbose, FileDeleteStatus, error_msg, at_data_file, uni_file, ios_var, uni_inp_orig, uni_inp, uni_out, &
       uni_log, Timestamp, Error
  USE ModPot, ONLY: Pot, Set_uni_EE_input, Define_Potential, Print_Potential_Data, Free_Potential_Memory
  USE ModInput, ONLY: end_of_input, command, Jobs_Data_list, Job_Data, Print_Header, Get_Input, Print_Input_Data, &
       Check_Input_Data, Free_Input_Data_Memory, global_job_name
  USE ModMethods, only: Header, Print_Parameters, Calculate

  IMPLICIT NONE

  CHARACTER(LEN = 100) :: argument
  CHARACTER(LEN = 100) :: nomeArq
  CHARACTER(LEN = 35) :: replace
  CHARACTER(LEN = 10) :: FileStatus="NEW"

  CHARACTER(LEN = 50), parameter :: view_command="less"
  CHARACTER(LEN = 50), parameter :: info_command="info"
  CHARACTER(LEN = 50), parameter :: perl = "perl"
  CHARACTER(LEN = 100), parameter :: copying_file="/home/yuriaoto/Programacao/QuAC/COPYING"
  CHARACTER(LEN = 100), parameter :: info_file="/home/yuriaoto/Programacao/QuAC/doc/QuAC.info"
  CHARACTER(LEN = 100), parameter :: script_at_data="/home/yuriaoto/Programacao/QuAC/src/Show_at_data.pl"
  CHARACTER(LEN = 100), parameter :: script_uni="/home/yuriaoto/Programacao/QuAC/src/Show_units.pl"
  CHARACTER(LEN = 100), parameter :: no_lines_command="echo $((\`tput lines\` - 1))"

  CHARACTER(LEN = 10) :: tmp_suffix, out_suffix, log_suffix
  Integer :: suffix_index, i_arg
  Logical :: tmp_inp_file_ex

  CALL GET_COMMAND_ARGUMENT (0,command)

  CALL GET_COMMAND_ARGUMENT (1,nomeArq)
  CALL GET_COMMAND_ARGUMENT (2,replace)

  do i_arg=1,command_argument_count(),1

     CALL GET_COMMAND_ARGUMENT (i_arg,argument)

     select case(argument)

     case("-V","--version")
        if(verbose.gt.0)then
           print*,"QuAC 4.0"
        end if
        stop

     case("-h", "--help")
        if(verbose.gt.0)then
           print*,"Usage: ",trim(command)," [options] [input file]"
           print*,"Options:"
           print*,"-V, --version             print version information and exit"
           print*,"-h, --help                print this help and exit"
           print*,"-i, --info                info documentation (press q to exit)"
           print*,"-c, --copying             show copying information (by less, press q to exit) and exit"
           print*,"-a, --atoms               list atomic data information (by less, press q to exit) and exit"
           print*,"-u, --units               list units conversion information (by less, press q to exit) and exit"
           print*,"-q, --quiet               suppress any messages to standard output (set verbose level at 0)"
           print*,"-v, --verbose             increment verbose level (defaut=1)"
           print*,"-r, --replace             replace the output and log files"
           print*,"-k, --keepTempInput       do not delete the temporary input file"
        end if
        STOP

     case("-i","--info")
        CALL SYSTEM(TRIM(info_command)//" "//TRIM(info_file))
        STOP

     case("-c","--copying")
        CALL SYSTEM(TRIM(view_command)//" "//TRIM(copying_file))
        STOP

     case("-a","--atoms")
        CALL SYSTEM( TRIM(perl)//" "//TRIM(script_at_data)//" "//trim(at_data_file)//" `"//&
             trim(no_lines_command)//"` | "//trim(view_command))
        STOP

     case("-u","--units")
        CALL SYSTEM( TRIM(perl)//" "//TRIM(script_uni)//" "//trim(uni_file)//" `"//&
             trim(no_lines_command)//"` | "//trim(view_command))
        STOP

     case("-q", "--quiet")
        verbose=0

     case("-v","--verbose")
        verbose = verbose+1

     case("-r", "--replace")
        FileStatus="REPLACE"

     case("-k", "--keepTempInput")
        FileDeleteStatus="KEEP"

     case default
        nomeArq=argument

     end select

  end do
  
  out_suffix = ".out"
  log_suffix = ".log"

  uni_inp_orig = 20
  uni_inp = 25
  uni_out = 30
  uni_log = 35
  CALL set_uni_EE_input(50)

  ! Calculate tmp_suffix
  suffix_index=1
  write(tmp_suffix,fmt='(I3)') suffix_index
  tmp_suffix = ".tmp_"//ADJUSTL(tmp_suffix)
  inquire(file=TRIM(nomeArq)//TRIM(tmp_suffix), exist=tmp_inp_file_ex)
  do while(tmp_inp_file_ex)
     suffix_index=suffix_index+1
     if(suffix_index.ge.100)then
        error_msg = "More than 100 temporary input files! Please, clean your working directory!"
        CALL Error(0)
     end if
     write(tmp_suffix,fmt='(I3)') suffix_index
     tmp_suffix = ".tmp_"//ADJUSTL(tmp_suffix)
     inquire(file=TRIM(nomeArq)//TRIM(tmp_suffix), exist=tmp_inp_file_ex)
  end do

  ! Open files
  if(verbose.gt.1) write(6,fmt='(A)',advance='no') "Opening files... "
  ! Original Input
  OPEN(Unit = uni_inp_orig, File = TRIM(nomeArq),STATUS='OLD', Iostat=ios_var)
  IF(ios_var.NE.0)THEN
     error_msg="Error opening "//TRIM(nomeArq)//". It exists?"
     CALL Error(0)
  END IF
  ! Input without comments
  OPEN(Unit = uni_inp, File = TRIM(nomeArq)//TRIM(tmp_suffix),STATUS='NEW', Iostat=ios_var)
  IF(ios_var.NE.0)THEN
     error_msg="Error opening "//TRIM(nomeArq)//TRIM(tmp_suffix)//". File already exists."
     CALL Error(2)
  END IF
  ! Output file
  OPEN(Unit = uni_out, File = TRIM(nomeArq)//TRIM(out_suffix), STATUS = TRIM(FileStatus), Iostat=ios_var)
  IF(ios_var.NE.0)THEN
     error_msg="Error opening "//TRIM(nomeArq)//TRIM(out_suffix)//". File already exists."
     CALL Error(0)
  END IF
  ! Log file
  OPEN(Unit = uni_log, File = TRIM(nomeArq)//TRIM(log_suffix), STATUS = TRIM(FileStatus), Iostat=ios_var)
  IF(ios_var.NE.0)THEN
     error_msg="Error opening "//TRIM(nomeArq)//TRIM(log_suffix)//". File already exists."
     CALL Error(0)
  END IF
  
  if(verbose.gt.1) write(6,fmt='(A)') "done."
  if(verbose.gt.1) write(6,fmt='(A)',advance='no') "Printing header and creating temporary input without comments... "

  CALL Print_Header()

  if(verbose.gt.1) write(6,fmt='(A)') "done."
  
  Close(uni_inp)
  Close(uni_inp_orig)

  ! Read input without comments
  OPEN(Unit = uni_inp, File = TRIM(nomeArq)//TRIM(tmp_suffix),STATUS='OLD', Iostat=ios_var)
  IF(ios_var.NE.0)THEN
     error_msg="Error opening "//TRIM(nomeArq)//TRIM(tmp_suffix)//". Error in creating?"
     CALL Error(2)
  END IF

  ! Get data for jobs from input file
  if(verbose.gt.1) write(6,fmt='(A)',advance='no') "Getting data from input file... "
  DO
     CALL Get_Input()
     IF(end_of_input) EXIT
  END DO
  if(verbose.gt.1) write(6,fmt='(A)') "done."

  CALL Timestamp()
  WRITE(uni_out,*)

  ! Get potentials from input file
  if(verbose.gt.1) write(6,fmt='(A)',advance='no') "Getting potential information from input file... "
  DO
     CALL Define_Potential()
     IF(.NOT.associated(Pot)) EXIT
  END DO
  if(verbose.gt.1) write(6,fmt='(A)') "done."

  CLOSE(uni_inp,STATUS=TRIM(FileDeleteStatus))

  ! Check input
  if(verbose.gt.1) write(6,fmt='(A)',advance='no') "Checking input data... "
  CALL Check_Input_Data()
  if(verbose.gt.1) write(6,fmt='(A)') "done."

  if(verbose.gt.2)then
     CALL Print_Input_Data()
     CALL Print_Potential_Data()
  end if

  WRITE(uni_out,FMT='("Consistent input.")')
  WRITE(uni_out,*)
  WRITE(uni_out,FMT='("~~~~~~~~~~> ",A," <~~~~~~~~~~")') TRIM(global_job_name)
  WRITE(uni_out,FMT='("Output file:")')
  WRITE(uni_out,*)
  
  WRITE(uni_log,FMT='("~~~~~~~~~~> ",A," <~~~~~~~~~~")') TRIM(global_job_name)
  WRITE(uni_log,FMT='("log file:")')
  WRITE(uni_log,*)


  ! Run the jobs
  if(verbose.gt.1) write(6,fmt='(A)') "Jobs calculation:"
  Job_Data=>Jobs_Data_list
  Do While(associated(Job_Data))

     if(verbose.gt.1) write(6,fmt='(A)') "Job """//TRIM(Job_Data%complete_job_name)//""""
     CALL Header()

     if(verbose.gt.1) write(6,fmt='(A)',advance='no') "Printing parameters... "
     CALL Print_Parameters()
     if(verbose.gt.1) write(6,fmt='(A)') "done."

     if(verbose.gt.1) write(6,fmt='(A)',advance='no') "Calculation... "
     CALL Calculate()
     if(verbose.gt.1) write(6,fmt='(A)') "done."

     Job_Data => Job_Data%next
  END Do

  if(verbose.gt.1) write(6,fmt='(A)',advance='no') "Freeing memory... "
  CALL Free_Input_Data_Memory()
  CALL Free_Potential_Memory()
  CALL Timestamp()
  if(verbose.gt.1) write(6,fmt='(A)') "done."
  
  WRITE(uni_out,'("End of program.")')
  if(verbose.gt.0) PRINT*,"QuAC finshed"

END PROGRAM QuAC
