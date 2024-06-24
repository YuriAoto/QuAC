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
! Module for general utilities.
!
MODULE ModUtil

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: real_kind, pi, zero, one, two, three, at_uni_dist, at_uni_ener, at_uni_dip_mom, at_uni_vel, at_uni_mass, quad_met, &
       error_msg, ios_var, ios_ok, ios_eof, ios_eol, verbose, FileDeleteStatus, status, units_len, uni_inp_orig, uni_inp, uni_out, &
       uni_log, at_data_file, uni_file, conv_thrsh_NR_default, eps_NR_default, max_iter_NR_default,&
       Unit_Converter, Atomic_Data, Newton_Raphson, Ordered_Search, Timestamp, Error, &
       Sph_Bessel_J, Sph_Bessel_N, Double_Factorial

  INTEGER, PARAMETER :: real_kind = SELECTED_REAL_KIND(15,300)

  ! Some constants
  REAL(kind = real_kind), PARAMETER :: pi = 3.141592653589793238462643_real_kind
  REAL(kind = real_kind), PARAMETER :: zero = ABS(0.0_real_kind)
  REAL(kind = real_kind), PARAMETER :: one = 1.0_real_kind
  REAL(kind = real_kind), PARAMETER :: two = 2.0_real_kind
  REAL(kind = real_kind), PARAMETER :: three = 3.0_real_kind

  ! The atomic units
  INTEGER, parameter :: units_len = 15
  CHARACTER(LEN = units_len), PARAMETER :: at_uni_dist = "a0"
  CHARACTER(LEN = units_len), PARAMETER :: at_uni_ener = "hartree"
  CHARACTER(LEN = units_len), PARAMETER :: at_uni_vel = "auv"
  CHARACTER(LEN = units_len), PARAMETER :: at_uni_mass = "me"
  CHARACTER(LEN = units_len), PARAMETER :: at_uni_dip_mom = "auDip"

  CHARACTER(LEN = 15), parameter :: quad_met = "SIMPSON" ! "TRAPEZOIDAL"

  ! Error messages and variable for use in Iostat statement
  CHARACTER (LEN = 200) :: error_msg
  INTEGER :: ios_var
  INTEGER, parameter :: ios_ok=0, ios_eof=-1, ios_eol=-2 
  INTEGER :: status

  INTEGER :: verbose=1
  CHARACTER(LEN = 10) :: FileDeleteStatus="DELETE"

  ! Units for inpu, output and log file
  Integer :: uni_inp_orig
  INTEGER :: uni_inp
  INTEGER :: uni_out
  INTEGER :: uni_log

  ! Units for atomic data and units libraries
  INTEGER, PARAMETER :: uni_at_data = 130
  INTEGER, PARAMETER :: uni_uni = 120

  ! Files of atomic data and units libraries
  CHARACTER (LEN = 75), PARAMETER :: at_data_file = '/home/yuriaoto/Programacao/modulos/atData.lib'
  CHARACTER (LEN = 75), PARAMETER :: uni_file = '/home/yuriaoto/Programacao/modulos/units.lib'
  !  include "./localInstINCLUDE"

  REAL(kind = real_kind), parameter :: conv_thrsh_NR_default = 1.0E-10_real_kind
  REAL(kind = real_kind), parameter :: eps_NR_default = 0.01_real_kind
  INTEGER, parameter :: max_iter_NR_default = 50

CONTAINS
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  REAL(kind = real_kind) FUNCTION Unit_Converter(uni_in,uni_out,quanti_in)
    !
    ! Structure of unit for string: <prefix>.<unit>, where <prefix> stands for p=pico=10E-12,
    ! M=Mega=10E6, etc. and <unit> means J=Joule, m=meter, etc. Check the units.lib ASCII file for
    ! all units.
    !
    ! status:
    ! 1 -> prefix_in not found
    ! 2 -> prefix_out not found
    ! 4 -> uni_in not_found
    ! 8 -> uni_out not_found
    ! 16 -> uni_in and uni_out of diferent signal
    !
    ! So, status = 0, means OK and status = 7 = 4+2+1
    ! means prefix_in, prefix_out and uni_in not found, whereas 
    ! status = 18 = 16 + 2 means uni_in and uni_out of diferent signals and 
    ! prefix_out not found and so on.
    !
    IMPLICIT NONE

    CHARACTER(LEN=units_len), INTENT(in) :: uni_in, uni_out
    REAL(kind = real_kind), INTENT(in), OPTIONAL :: quanti_in

    INTEGER :: status_ini
    INTEGER :: ind_in, ind_out
    INTEGER :: prefix_factor_in, prefix_factor_out
    LOGICAL :: found_pre_in, found_pre_out
    INTEGER :: signal_uni_in, signal_uni_out
    CHARACTER(LEN=units_len) :: info = ""
    REAL(kind = real_kind) :: aux_real
    INTEGER :: aux_int
    INTEGER :: signal
    REAL(kind = real_kind) :: factor_in, factor_out

    status_ini = 0

    OPEN(Unit = uni_uni, File = TRIM(uni_file), Status = 'OLD')

    ind_in = INDEX(uni_in,".")
    ind_out = INDEX(uni_out,".")

    prefix_factor_in = 0
    prefix_factor_out = 0

    IF(ind_in.GT.0)THEN
       found_pre_in = .FALSE.
    ELSE
       found_pre_in = .TRUE.
    END IF

    IF(ind_out.GT.0)THEN
       found_pre_out = .FALSE.
    ELSE
       found_pre_out = .TRUE.
    END IF

    ! Get expoent of prefix
    DO WHILE(info.NE."---")

       READ(uni_uni,*,Iostat=ios_var) info,aux_int

       IF(info.EQ.uni_in(:ind_in-1))THEN
          prefix_factor_in = aux_int
          found_pre_in = .TRUE.
       END IF

       IF(info.EQ.uni_out(:ind_out-1))THEN
          prefix_factor_out = aux_int
          found_pre_out = .TRUE.
       END IF

    END DO
    BACKSPACE(uni_uni)

    IF(.NOT.found_pre_in) status_ini = status_ini + 1
    IF(.NOT.found_pre_out) status_ini = status_ini + 2

    signal_uni_in = -1
    signal_uni_out = -1
    ! Get factor
    DO WHILE(signal_uni_in.LT.0.OR.signal_uni_out.LT.0)

       READ(uni_uni,*,Iostat=ios_var) info,aux_real,signal
       IF(ios_var.NE.0) EXIT

       IF(info.EQ.uni_in(ind_in+1:))THEN
          factor_in = aux_real
          signal_uni_in = signal
       END IF

       IF(info.EQ.uni_out(ind_out+1:))THEN
          factor_out = aux_real
          signal_uni_out = signal
       END IF

    END DO

    CLOSE(uni_uni)

    IF(signal_uni_in.LT.0) status_ini = status_ini + 4
    IF(signal_uni_out.LT.0) status_ini = status_ini + 8
    IF(signal_uni_in.NE.signal_uni_out) status_ini = status_ini + 16

    status = status_ini

    IF(PRESENT(quanti_in).AND.status_ini.EQ.0)THEN
       Unit_Converter = 10**(one*(prefix_factor_in-prefix_factor_out))*quanti_in*factor_in/factor_out
    ELSE
       Unit_Converter = zero
    END IF

  END FUNCTION Unit_Converter
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!  
  SUBROUTINE Atomic_Data(symbol, name, Z, mass, abundance)
    !
    ! Get atomic data from atData.lib ASCII file
    !
    ! If the atomic symbol is not found in atData.lib file, name take the value "-empty-" and the numeric quantities -1, if passed
    !

    IMPLICIT NONE

    CHARACTER(LEN = 7), INTENT(in) :: symbol
    CHARACTER(LEN = 15), INTENT(out), optional :: name
    INTEGER, INTENT(out), optional :: Z
    REAL(kind = real_kind), INTENT(out), optional :: mass, abundance

    CHARACTER(LEN = 7) :: info
    CHARACTER(len = 15) :: aux_name
    REAL(kind = real_kind) :: aux_real
    INTEGER :: aux_Z

    ! Abre biblioteca de dados atômicos
    OPEN(Unit = uni_at_data, File = TRIM(at_data_file), Status = "OLD")

    ios_var = 0
    DO WHILE(ios_var.EQ.0)
       READ(uni_at_data,*,Iostat=ios_var) info
       IF(info.EQ.symbol)THEN
          
          READ(uni_at_data,*,Iostat=ios_var) aux_name
          IF(ios_var.NE.0) EXIT
          if(present(name)) name = aux_name

          READ(uni_at_data,*,Iostat=ios_var) aux_Z
          IF(ios_var.NE.0) EXIT
          if(present(Z)) Z = aux_Z
          
          READ(uni_at_data,*,Iostat=ios_var) aux_real
          IF(ios_var.NE.0) EXIT
          if(present(mass)) mass = aux_real
          
          READ(uni_at_data,*,Iostat=ios_var) aux_real
          IF(ios_var.NE.0) EXIT
          if(present(abundance)) abundance = aux_real
          
          EXIT
          
       ENDIF
    END DO

    IF(ios_var.NE.0.and.present(name)) name = "-empty-"
    IF(ios_var.NE.0.and.present(Z)) Z = -1
    IF(ios_var.NE.0.and.present(mass)) mass = -one
    IF(ios_var.NE.0.and.present(abundance)) abundance = -one

    CLOSE(uni_at_data)

  END SUBROUTINE Atomic_Data
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  INTEGER FUNCTION Ordered_Search (point,vector,initial_index)
    !
    ! Given the ordered (increasing or decreasing) vector, return the lower index of
    ! the two index that define the interval of teh vector in which te point belongs.
    ! The optional parameter, initial_index, shift the result by (initial_index_effective - 1)
    ! and should be used if a vector with a non fortran-standard range is passed and correspond
    ! to the origninal initial index of the vector.
    !

    REAL(kind = real_kind), INTENT(in) :: point
    REAL(kind = real_kind), DIMENSION(:), INTENT(in) :: vector
    INTEGER, INTENT(in), OPTIONAL :: initial_index

    INTEGER :: initial_index_effective
    INTEGER :: i_min,i_max,i_middle,no_pts
    LOGICAL :: is_crescent

    IF (PRESENT(initial_index)) THEN
       initial_index_effective=initial_index
    ELSE
       initial_index_effective=1
    END IF       

    no_pts=SIZE(vector)
    is_crescent=vector(no_pts).GT.vector(1)

    !
    ! Extreme cases. Necessary, since, for increasing ordered points, the intervals are
    ! [vector(1),vector(2)),[vector(2),vector(3)),...,[vector(n-1),vector(n))
    ! and the last point keeps out. For decreasing ordered points, the intervals are
    ! (vector(1),vector(2)],(vector(2),vector(3)],...,(vector(n-1),vector(n)]
    ! and the first point keeps out.
    !
    ! Study general way to do this
    !
    IF(point.EQ.vector(1)) THEN
       Ordered_Search = 1 + (initial_index_effective - 1)
       RETURN
    ELSE
       IF(point.EQ.vector(no_pts)) THEN
          Ordered_Search = no_pts-1 + (initial_index_effective - 1)
          RETURN
       END IF
    END IF
    IF(is_crescent.EQV.(point.LT.vector(1))) THEN
       Ordered_Search = 0 + (initial_index_effective - 1)
       RETURN
    ELSE
       IF(is_crescent.EQV.(point.GT.vector(no_pts))) THEN
       Ordered_Search = no_pts + (initial_index_effective - 1)
          RETURN
       END IF
    END IF

    i_min=0
    i_max=no_pts+1

    DO WHILE(i_max-i_min.GT.1)
       i_middle=(i_min+i_max)/2
       IF(is_crescent.EQV.(point.GE.vector(i_middle)))THEN
          i_min=i_middle
       ELSE
          i_max=i_middle
       END IF
    END DO

    Ordered_Search = i_min + (initial_index_effective - 1)

  END FUNCTION Ordered_Search
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  Real(kind = real_kind) function Newton_Raphson(Func,x_ini,conv_thrsh,eps,max_iter)

    implicit none

    Interface
       Function Func(x)
         Import :: real_kind
         real(kind = real_kind) :: Func
         real(kind = real_kind), intent(in) :: x
       end Function Func
    end Interface
    real(kind = real_kind), intent(in) :: x_ini
    real(kind = real_kind), optional, intent(in) :: conv_thrsh
    real(kind = real_kind), optional, intent(in) :: eps
    Integer, optional, intent(in) :: max_iter

    real(kind = real_kind) :: conv_thrsh_used = 1.0E-10_real_kind,eps_used = 0.01_real_kind
    Integer:: max_iter_used = 50

    real(kind = real_kind) :: x_old
    real(kind = real_kind) :: f, f_der
    Integer :: i_iter

! status = 0 -> ok
! status = 1 -> potential error
! status = 2 -> zero derivative
! status = 3 -> number of iterations reached

    conv_thrsh_used = conv_thrsh_NR_default
    eps_used = eps_NR_default
    max_iter_used = max_iter_NR_default

    if(present(conv_thrsh)) conv_thrsh_used = conv_thrsh
    if(present(eps)) eps_used = eps
    if(present(max_iter)) max_iter_used = max_iter

    Newton_Raphson = x_ini
    x_old = Newton_Raphson + 2*conv_thrsh_used

    f = Func(Newton_Raphson)
    i_iter = 1
    do while (abs(x_old-Newton_Raphson).gt.conv_thrsh_used.and.i_iter.lt.max_iter_used)

       f_der = (Func(Newton_Raphson+(eps_used/2)) - Func(Newton_Raphson-(eps_used/2)))/eps_used
       if(f_der.eq.zero)then
          status = 2
          return
       end if
       x_old = Newton_Raphson

       Newton_Raphson = Newton_Raphson - f/f_der
       f = Func(Newton_Raphson)

       if(status.ne.0) return
       
       i_iter = i_iter + 1
       
    end do

    if(i_iter.ge.max_iter_used) then
       status=3
    else
       status=0
    end if

  end function Newton_Raphson
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  SUBROUTINE Timestamp()
    !
    ! TIMESTAMP prints the current YMDHMS date as a time stamp, like 01 January 1970 00:00:00.000, in unit uni_out
    !
    !  Author:
    !
    !    John Burkardt (http://orion.math.iastate.edu/burkardt/f_src/f_src.html)
    !    Adapted by Yuri Alexandre Aoto (29 March 2009)
    !

    IMPLICIT NONE

    INTEGER :: year, month, day, hour, minute, second, ms
    INTEGER, DIMENSION(8) :: values
    CHARACTER(len = 8) :: date
    CHARACTER(len=5) :: zone
    CHARACTER(len=10) :: time
    CHARACTER(len = 9), PARAMETER, DIMENSION(12) :: month_list = (/ &
         'January  ', 'February ', 'March    ', 'April    ', &
         'May      ', 'June     ', 'July     ', 'August   ', &
         'September', 'October  ', 'November ', 'December ' /)

    CALL DATE_AND_TIME ( date, time, zone, values )

    year = values(1)
    month = values(2)
    day = values(3)
    hour = values(5)
    minute = values(6)
    second = values(7)
    ms = values(8)

    WRITE(uni_out, '(i2,1x,a,1x,i4,1x,a1,1x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
         day, TRIM(month_list(month)), year, "-", hour, ':', minute, ':', second, '.', ms

  END SUBROUTINE Timestamp
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  SUBROUTINE Error(error_type)
    !
    ! Print an error information, stored in the module variable error_msg, and stop the execution of the program.
    ! The error type stands for 0 = Input error, 1 = Calculation error and 2 = Coding error.
    !
    IMPLICIT NONE

    INTEGER, INTENT(in) :: error_type

    LOGICAL :: already_opened

    inquire(unit=uni_out,opened=already_opened)
    if(already_opened) WRITE(uni_out,*)

    SELECT CASE(error_type)

    CASE(0)
       
       if(already_opened) WRITE(uni_out,FMT='("Input error. ",A)') TRIM(error_msg)
       if(verbose.gt.0) PRINT*,"Input error. ",TRIM(error_msg)

    CASE(1)
       
       if(already_opened) WRITE(uni_out,FMT='("Calculation error. ",A)') TRIM(error_msg)
       if(verbose.gt.0) PRINT*,"Calculation error. ",TRIM(error_msg)

    CASE(2)
       
       if(already_opened) WRITE(uni_out,FMT='("Coding error. ",A)') TRIM(error_msg)
       if(already_opened) WRITE(uni_out,FMT='("Please, report this bug.")')
       
       if(verbose.gt.0) PRINT*,"Coding error. ",TRIM(error_msg)
       if(verbose.gt.0) PRINT*,"Please, report this bug."

    CASE DEFAULT
       
       if(already_opened) WRITE(uni_out,FMT='("Unknow error: ",I2," - ",A)') error_type,TRIM(error_msg)
       if(verbose.gt.0) PRINT*,"Unknow error: ",error_type,"-",TRIM(error_msg)

    END SELECT

    if(already_opened) WRITE(uni_out,FMT='("Program finished.")')
    if(verbose.gt.0) PRINT*,"Program finished."
    CLOSE(uni_inp,STATUS=TRIM(FileDeleteStatus))
    STOP

  END SUBROUTINE Error
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  REAL(kind = real_kind) FUNCTION Sph_Bessel_J(l,r)
    !
    ! Calculate the spherical Bessel function of first kind of order l
    !
    IMPLICIT NONE

    INTEGER, INTENT(in) :: l
    REAL(kind = real_kind), INTENT(in) :: r

    INTEGER :: i
    REAL(kind = real_kind) :: bessJ_0, bessJ_1, bessJ_tmp

    IF(ABS(r).LT.1.0E-150_real_kind)THEN
       SELECT CASE(l)
       CASE(0)
          Sph_Bessel_J = 1
          RETURN
       CASE DEFAULT
          Sph_Bessel_J = 0
          RETURN
       END SELECT
    ENDIF

    bessJ_0 = SIN(r)

    bessJ_1 = -COS(r) + SIN(r)/r

    IF(l.EQ.0)THEN
       Sph_Bessel_J = bessJ_0
       RETURN
    END IF

    IF(l.EQ.1)THEN
       Sph_Bessel_J = bessJ_1
       RETURN
    END IF

    ! Compute r*j_l(r) using the following recurrence formula
    !
    !     (2l+1)*j_l(r) = r[j_l+1(r)+j_l-1(r)]
    !
    ! i.e.:
    !
    !     (2l+1)*{r*j_l(r)} = r[ {r*j_l+1(r)} + {r*j_l-1(r)} ]
    !
    ! If (l(l+1)/r²) > 100, uses the following approximation
    !
    !     j_l(r) = r^l/(2l+1)!!
    !
    ! i.e:
    !
    !     {r*j_l(r)} = r^(l+1)/(2l+1)!!
    !
    IF((l*(l+1))/r**2.GT.100)THEN
       bessJ_tmp = (1/(Double_Factorial(2*l+1)))*(r**(l+1))
    ELSE
       DO i=2,l,1

          bessJ_tmp = ((2*i-1)/r)*BessJ_1 - BessJ_0

          bessJ_0 = bessJ_1
          bessJ_1 = bessJ_tmp
       END DO
    END IF

    Sph_Bessel_J = bessJ_tmp

  END FUNCTION Sph_Bessel_J
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  REAL(kind = real_kind) FUNCTION Sph_Bessel_N(l,r)
    !
    ! Calculate the spherical Bessel function of second kind of order l
    !
    IMPLICIT NONE

    INTEGER, INTENT(in) :: l
    REAL(kind = real_kind), INTENT(in) :: r

    INTEGER :: i
    REAL(kind = real_kind) :: bessN_0, bessN_1, bessN_tmp

    bessN_0 = -COS(r)

    bessN_1 = -SIN(r) - COS(r)/r

    IF(l.EQ.0)THEN
       Sph_Bessel_N = bessN_0
       RETURN
    END IF

    IF(l.EQ.1)THEN
       Sph_Bessel_N = bessN_1
       RETURN
    END IF

    ! Compute r*n_l(r) using the following recurrence formula
    !
    !     (2l+1)*n_l(r) = r[n_l+1(r)+n_l-1(r)]
    !
    ! i.e.:
    !
    !     (2l+1)*{r*n_l(r)} = r[ {r*n_l+1(r)} + {r*n_l-1(r)} ]
    !
    DO i=2,l,1

       bessN_tmp = ((2*i-1)/r)*bessN_1 - bessN_0

       bessN_0 = bessN_1
       bessN_1 = bessN_tmp
    END DO

    Sph_Bessel_N = bessN_tmp

  END FUNCTION Sph_Bessel_N
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
  REAL(kind = real_kind) FUNCTION Double_Factorial(k)
    !
    ! Calculate the double factorial function:
    ! k!! = k*(k-2)*(k-4)*...*1 for odd k
    ! k!! = k*(k-2)*(k-4)*...*2 for even k
    !
    IMPLICIT NONE

    INTEGER, INTENT(in) :: k
    integer :: i

    if(k.eq.0.or.k.eq.1)then
       Double_Factorial = one
       return
    end if

    Double_Factorial = k*one
    i = k-2
    do while(i.gt.0)
       Double_Factorial = Double_Factorial*i
       i = i-2
    end do

  END FUNCTION Double_Factorial
!!!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$!!
END MODULE ModUtil
