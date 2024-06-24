!!$
          ! Valor máximo que a função pode atingir.
       CASE("MaxFuncao")

          If(Def(i_MET)%MaxFun)Then
             erroMSG = "Foi passado o valor máximo para função mais de uma vez em "//TRIM(Dados(i_MET)%Metodo)//" - "&
                  &//TRIM(Dados(i_MET)%nome)
             CALL Erro(0)
          End If

          Backspace(uniINP)
          Read(uniINP,*,iostat=iosVAR) info,RealAUX
          If(iosVAR.NE.0)Then
             erroMSG = "Falha ao ler MaxFuncao em "//TRIM(Dados(i_MET)%Metodo)//" - "//TRIM(Dados(i_MET)%nome)
             CALL Erro(0)
          End If

          If(RealAUX.LT.0)Then
             erroMSG = "Encontrado valor negativo para MaxFuncao em "//TRIM(Dados(i_MET)%Metodo)//" - "//TRIM(Dados(i_MET)%nome)
             CALL Erro(0)
          End If

          Do i_E = 1,Dados(i_MET)%Energia%n,1
             Do i_l = 0,Dados(i_MET)%l%n,1
                Dados(i_MET)%Arg(i_E,i_l)%Max_Psy = RealAUX
             End Do
          End Do

          Def(i_MET)%MaxFun = .true.

          ! Condição inicial para a integração numérica
       CASE("FuncaoInicial")

          If(Def(i_MET)%PsyIni)Then
             erroMSG = "Foi passada a condição inicial mais de uma vez em "//TRIM(Dados(i_MET)%Metodo)//" - "&
                  &//TRIM(Dados(i_MET)%nome)
             CALL Erro(0)
          End If

          Do i_E = 1,Dados(i_MET)%Energia%n,1
             Do i_l = 0,Dados(i_MET)%l%n,1

                Read(uniINP,*,Iostat=iosVAR) RealAUX,RealAUX2
                If(iosVAR.EQ.0)Then
                   Dados(i_MET)%Arg(i_E,i_l)%PsyIni(1) = RealAUX
                   Dados(i_MET)%Arg(i_E,i_l)%PsyIni(2) = RealAUX2
                Else
                   erroMSG = "Deve ser passado valor de FunçãoInicial para cada valor de l e E, em "//&
                        &TRIM(Dados(i_MET)%Metodo)//" - "//TRIM(Dados(i_MET)%nome)
                   CALL Erro(0)
                End If

             End Do
          End Do

          Def(i_MET)%PsyIni = .true.
          ! A integração numérica inicia em IniciaMultEner*Ener. Passado nas unidades do potencial e armazenado em unidades adimensionais ----> Isso muda da última versão
       CASE("InicialMultEner")

          If(Def(i_MET)%Rmin)Then
             erroMSG = "Foi passado o valor inicial de integração mais de uma vez em "//TRIM(Dados(i_MET)%Metodo)//" - "&
                  &//TRIM(Dados(i_MET)%nome)
             CALL Erro(0)
          End If

          Backspace(uniINP)
          Read(uniINP,*,iostat=iosVAR) info,Def(i_MET)%ValorRmin
          If(iosVAR.NE.0)Then
             erroMSG = "Falha ao ler InicialMultEner em "//TRIM(Dados(i_MET)%Metodo)//" - "//TRIM(Dados(i_MET)%nome)
             CALL Erro(0)
          End If

          If(RealAUX.LT.0)Then
             erroMSG = "Encontrado valor negativo para InicialMultEner em "//TRIM(Dados(i_MET)%Metodo)//" - "&
                  &//TRIM(Dados(i_MET)%nome)
             CALL Erro(0)
          End If

          If(TRIM(Dados(i_MET)%Metodo).EQ."PARTLIVRE")Then
             erroMSG = "Deve ser passado Rinicial e não InicialMultEner para a partícula livre."
             CALL Erro(0)
          End If

          SELECT CASE(TRIM(Pot%Tipo))

          CASE("EsfRigida","PotQuad")

             erroMSG = "Deve ser passado RInicial e não InicialMultEner para potenciais do tipo "//TRIM(Pot%Tipo)
             CALL Erro(0)

          CASE("LenJon","LenJon610","MurSorbie","EsfMole","SplLinear","SplCubNat","SplCubFLJ","SplCubFLJ68")

             CALL Unidades(EnerTMP,unidadeEner(i_MET),uniEner)
             If(UnidadesIguais())Then

                Do i_E = 1,Dados(i_MET)%Energia%n,1

                   CALL Unidades(Dados(i_MET)%Energia%Entry(i_E),unidadeEner(i_MET),uniEner,EnerTMP)
                   CALL VerErroUni(unidadeEner(i_MET))
                   
                   RealAUX = RaizPotEner(EnerTMP,Def(i_MET)%ValorRmin)
                   Do i_l = 0,Dados(i_MET)%l%n,1
                      if(Dados(i_MET)%vetor_onda%Entry(i_E).lt.0)then
                         Dados(i_MET)%Arg(i_E,i_l)%Rmin = RealAUX
                      else
                         Dados(i_MET)%Arg(i_E,i_l)%Rmin = RealAUX*Dados(i_MET)%vetor_onda%Entry(i_E)
                      end if
                   End Do
                      
                End Do
                   
             Else
                
                CALL Unidades(EnerTMP,unidadeEner(i_MET),uniVel)
                If(UnidadesIguais())Then
                   
                   Do i_E = 1,Dados(i_MET)%Energia%n,1
                      
                      CALL Unidades(Dados(i_MET)%Energia%Entry(i_E),unidadeEner(i_MET),uniVel,EnerTMP)
                      CALL VerErroUni(unidadeEner(i_MET))

                      EnerTMP = (MassaRed*EnerTMP**2)/2

                      RealAUX = RaizPotEner(EnerTMP,Def(i_MET)%ValorRmin)
                      Do i_l = 0,Dados(i_MET)%l%n,1
                         if(Dados(i_MET)%vetor_onda%Entry(i_E).lt.0)then
                            Dados(i_MET)%Arg(i_E,i_l)%Rmin = RealAUX
                         else
                            Dados(i_MET)%Arg(i_E,i_l)%Rmin = RealAUX*Dados(i_MET)%vetor_onda%Entry(i_E)
                         end if
                      End Do
                         
                   End Do
                      
                Else
                   erroMSG = "Inconsistência de unidades. Nem energia nem velocidade."
                   CALL Erro(0)
                End If

             End If

          CASE DEFAULT

             erroMSG = "Inconsistência de potencial na leitura de InicialMultEner." 
             CALL Erro(0)

          END SELECT

          Def(i_MET)%Rmin = .true.
          Def(i_MET)%IniMEner = .true.










    coef = pi*r_bar/(2*beta)
    tang = beta*(r/r_bar - 1)
    sec_sqr = 1 + tang**2
    
    g = coef*sec_sqr
    coef = coef*pi

    g_fst_deriv = coef*sec_sqr*tang
    coef = coef*pi/2

    g_snd_deriv = coef*sec_sqr*(3*sec_sqr - 2)

    SecDiffer_tilde = g**2*SecDiffer(J,r,wave_vector,Ener) + g_snd_deriv/(2*g) - 3*(g_fst_deriv/g)**2/4




    CASE("NUMEROV")

       i = 1 ! Passo onde a função está sendo calculada
       r = Arg%Rmin + i*Arg%passo ! r onde a função está sendo calculada

       psi2 = Arg%PsyIni(1)
       ! Primeiro passo para Numerov
       psi1 = PrimeiroPassoNum(Arg%Rmin,Arg%PsyIni(1),Arg%PsyIni(2),Arg%passo)

       ! psi1 -> função que acabou de ser calculada
       ! psi2 -> função no passo imediatamente anterior

       funcCalc = psi1
       CALL testesSobreFuncao ()
       If(FunAlta)Then
          RETURN
       End If

       multPasso=MIN(MAX(NINT(5.0_real_kind/wave_vector),1),comecoSalva/10)

       passoAtual=Arg%passo

       r = Arg%Rmin + 2*passoAtual ! Passo onde a função está sendo calculada

       ! NUMEROV
       do i=2,FimDaInt,1

          If(i.eq.ComecoSalva+1)then
             passoAtual=multPasso*passoAtual
             psi2=psiANTpnovo
          End If

          r = r+passoAtual ! Passo onde a função está sendo calculada

          psiTemp = psi1
          psi1 = PassoNumerov(r,psi1,psi2,passoAtual)

          psi2 = psiTemp

          if (i.eq.ComecoSalva-multPasso) then
             psiANTpnovo=psi1
             print*,"qwddassda"
          End if

          funcCalc = psi1
          CALL testesSobreFuncao ()
          If(FunAlta)Then
             RETURN
          End If

       End do

!!$       passoNOVO=multPasso*Arg%passo
!!$       i = FimDaInt + 1 ! Passo onde a função está sendo calculada
!!$       r = Arg%Rmin + i*Arg%passo ! r onde a função está sendo calculada
!!$
!!$       psi2 = psiANTpnovo
!!$
!!$       do i=ComecoSalva+1,FimDaInt,1
!!$
!!$          r = Arg%Rmin + i*Arg%passo ! Passo onde a função está sendo calculada
!!$
!!$          psiTemp = psi1
!!$          psi1 = PassoNumerov(r,psi1,psi2,passoNOVO)
!!$          psi2 = psiTemp
!!$
!!$          funcCalc = psi1
!!$          CALL testesSobreFuncao ()
!!$          If(FunAlta)Then
!!$             RETURN
!!$          End If
!!$
!!$       End do


       ! Previsor-Corretor de Adams-Bashforth-Moulton. retirado de:
       !


       ! Primeiro passo para Numerov
       V0 = SecDiffer(J,Arg%Rmin,wave_vector,Ener)
       V1 = SecDiffer(J,Arg%Rmin + Arg%passo,wave_vector,Ener)
       V2 = SecDiffer(J,Arg%Rmin + 2*Arg%passo,wave_vector,Ener)

       tmp = Arg%PsyIni(1)*(1 - V2*Arg%passo**2/24) + Arg%passo*Arg%PsyIni(2)*(1 - V2*Arg%passo**2/12)
       tmp = tmp + (Arg%passo**2/24)*7*V0*Arg%psyIni(1) - (V2*Arg%passo**4/36)*V0*Arg%PsyIni(1)
       psi1 = tmp/(1 - V1*Arg%passo**2/4 + V1*V2*Arg%passo**4/18)

       ! psi1 -> função que acabou de ser calculada
       ! psi2 -> função no passo imediatamente anterior





  Real(kind = real_kind) function Root_Finder(Func,x_min,x_max,conv_thrsh,max_iter,status)

    implicit none

    Interface
       Function Func(x)
         Import :: real_kind
         real(kind = real_kind) :: Func
         real(kind = real_kind), intent(in) :: x
       end Function Func
    end Interface
    real(kind = real_kind), intent(in) :: x_min, x_max
    real(kind = real_kind), intent(in) :: conv_thrsh
    Integer, intent(in) :: max_iter
    Integer, intent(out) :: status

    real(kind = real_kind) :: x_inf, x_sup, x_middle
    real(kind = real_kind) :: func_inf, func_sup, func_middle
    Integer :: i_iter

    x_sup = x_max
    x_inf = x_min

    func_sup=Func(x_max)
    func_inf=Func(x_min)

    if(func_inf*func_sup.gt.0)then
       status = -1
       return
    end if

    i_iter = 1
    do while (abs(x_sup-x_inf).gt.conv_thrsh.and.i_iter.lt.max_iter)
       x_middle = (x_sup+x_inf)/2
       func_middle = Func(x_middle)

       if(func_middle*func_inf.gt.zero)then
          x_inf = x_middle
          func_inf = func_middle
       else
          x_sup = x_middle
          func_sup = func_middle
       end if

       i_iter = i_iter + 1

    end do

    if(i_iter.ge.max_iter) then
       status=1
    else
       status=0
    end if

    Root_Finder = x_middle

  end function Root_Finder




  SUBROUTINE Least_Square_Len_Jon()

    IMPLICIT NONE

    REAL(kind = real_kind) :: g1g1=0, g1g2=0, g2g2=0, fg1=0, fg2=0
    REAL(kind = real_kind) :: invDe_r6, Det, Det1, Det2
    INTEGER :: aux_int

    DO aux_int = 1,SIZE(Pot%r),1
       invDe_r6 = 1/Pot%R(aux_int)**6

       g1g1 = g1g1 + invDe_r6**2
       g1g2 = g1g2 - invDe_r6**3
       g2g2 = g2g2 + invDe_r6**4

       fg1 = fg1 - Pot%V(aux_int)*invDe_r6
       fg2 = fg2 + Pot%V(aux_int)*invDe_r6**2
    END DO

    Det = g1g1*g2g2 - g1g2*g1g2
    Det1 = fg1*g2g2 - fg2*g1g2
    Det2 = g1g1*fg2 - g1g2*fg1

    Pot%radius = EXP(LOG(2*Det2/Det1)/6)
    Pot%epsilon = (Det1**2)/(4*Det*Det2)

  END SUBROUTINE Least_Square_Len_Jon
