F90 = /usr/bin/gfortran -std=gnu -Wall -pedantic -llapack -std=gnu -ffpe-trap=invalid,zero,overflow #-g -O0

srcdir = /home/yuriaoto/Programacao/QuAC/src/

QuAC.x: function_call_system.o Mod_Util.o Mod_DiatPot.o Mod_Input.o Mod_System.o Mod_Calc.o Mod_Methods.o Main_QuAC.o
	$(F90) -o QuACexec function_call_system.o Mod_Util.o Mod_DiatPot.o Mod_Input.o Mod_System.o Mod_Calc.o Mod_Methods.o Main_QuAC.o

Main_QuAC.o: $(srcdir)Main_QuAC.f90 Mod_Util.o Mod_Input.o Mod_Methods.o Mod_DiatPot.o
	$(F90) -c $(srcdir)Main_QuAC.f90

Mod_Methods.o: $(srcdir)Mod_Methods.f90 Mod_Util.o Mod_System.o Mod_DiatPot.o Mod_Input.o Mod_Calc.o
	$(F90) -c $(srcdir)Mod_Methods.f90

Mod_Input.o: $(srcdir)Mod_Input.f90 Mod_Util.o Mod_DiatPot.o Mod_System.o
	$(F90) -c $(srcdir)Mod_Input.f90

Mod_Calc.o: $(srcdir)Mod_Calc.f90 Mod_Util.o Mod_DiatPot.o Mod_System.o
	$(F90) -c $(srcdir)Mod_Calc.f90

Mod_System.o: $(srcdir)Mod_System.f90 Mod_Util.o Mod_DiatPot.o
	$(F90) -c $(srcdir)Mod_System.f90

Mod_DiatPot.o: $(srcdir)Mod_DiatPot.f90 Mod_Util.o
	$(F90) -c $(srcdir)Mod_DiatPot.f90

function_call_system.o: $(srcdir)function_call_system.c
	$(CC) -c $(srcdir)function_call_system.c

Mod_Interpolacao.o: $(srcdir)Mod_Interpolacao.f90  Mod_Util.o
	$(F90) -c $(srcdir)Mod_Interpolacao.f90

Mod_Util.o: $(srcdir)Mod_Util.f90
	$(F90) -c $(srcdir)Mod_Util.f90

clean:
	rm Mod_Util.o modutil.mod Mod_System.o modsystem.mod Mod_Methods.o modmethods.mod Mod_DiatPot.o modpot.mod Mod_Input.o modinput.mod Mod_Calc.o modcalc.mod function_call_system.o Main_QuAC.o QuACexec
