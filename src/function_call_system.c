/* 
  QuAC - Quantum Atomic Collisions
  Copyright (C) 2009, 2010, 2011, 2012  Yuri Aoto

  This file is part of QuAC.

  QuAC is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

   Function to execute a system command and return its standard output.
   To be called from Fortran in the following way:

   CALL System_call_c2f(command,stdout)

   with command and stdout declared like:

   Character(LEN=m,kind=c_char) :: command
   Character(LEN=n) :: stdout

   with n and m positive integers.
*/
#include <stdio.h>

int system_call_stdout_c2f_(char *command, char *stdout, long int command_len, long int stdout_len){

  FILE *fp;
  int i;
  char c;

  fp = popen(command,"r");

  c=fgetc(fp);
  for (i=0;i<stdout_len;i++){
    if (c != EOF){
      stdout[i]=c;
      c=fgetc(fp);
    }else{
      stdout[i]=' ';
    }
  }

  return pclose(fp);
}
