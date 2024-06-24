c=====================================================================
      PROGRAM  JLD_RE_test
c Driver program calling subroutine JDL for calculating scattering 
c lengths a_s for three model 15-level Lennard-Jones(2n,n) (n=4, 5 & 6)
c potentials using, in turn, both the
c  (a) one-parameter (Rbar) mapping function  y(r)= (r - Rbar)/(r + Rbar)
c                          and 
c  (b) two-parameters (Rbar,beta) mapping function
c                  y(r) = 2*arctan[beta*(r/Rbar-1)]/Pi 
c The main points illustrated here are the convergence tests for 
c   a_s-values calculated using Johnson-s log-derivative method alone 
c   (denoted `SL-JLD') and the results (labeled `SL-JLD-RE') obtained by 
c   also applying a (Mes,Mes/2)-Richardson extrapolation procedure to 
c   those results with respect to number of grid points (Mes)
c
c This code may be used as a template for creating an interface to JLD. 
c The main thing that must be done in order to use this routine 
c on your problem is to modify subroutine-function which calculates 
c the potential U(r) and apply the Richardson Extrapolation procedure. 
c after performing two calls to subroutine JDL with different numbers 
c of mesh points (as illustratrd below)
c-----------------------------------------------------------------------
c         COPYRIGHT 2011  by  V.V. Meshkov and A.V. Stolyarov
c-----------------------------------------------------------------------
      INTEGER*4  i, Mes, Mes1
      DOUBLE PRECISION  Scale, Rmin
      DOUBLE PRECISION  Rbar, beta, gamma
      DOUBLE PRECISION  Yl, Yl2, SL_err, SL_JLD, SL_JLD2, SL_JLD_RE
      INTEGER*4  n, Imap
      DOUBLE PRECISION  De, Re, C4
      COMMON/LJ/De,Re,n
      DOUBLE PRECISION  u_LJ_2n_n 
      EXTERNAL          u_LJ_2n_n 
c************************** Initialization *****************************
c Set the scaling factor of the energy
      Scale = 1.d0 
c Set the left boundary point
      Rmin = 0.3d0 
c Set the equilibrium distance of LJ(2n,n) potentials 
      Re = 1.d0 
c----------------------------------------------------------------------
      open(10,file='LJ(2n,n)_testing_REAL8.out') ! Name of output file
c======================================================================
      write(10,*)' n map    Mes       SL-JLD           SL-JLD-RE'//
     + '       RE-correction'
c----------------------------------------------------------------------
      DO  n= 6,6,1
c** Loop over cases of three 15-level LJ(2n,n) potentials for n= 4, 5 
c  and 6 with well depths = 1000, 2165 and 3761 cm-1
c* Set non-zero long range C_4 coefficient for the case n=4 
          C4 = 0.d0
c          IF(n.EQ.4) C4 = 2.d3
c----------------------------------------------------------------------
c          DO  Imap = 1,1
c  Setup of the mapping parameter Rbar
              Rbar = Re 
              mes = 16000
c----------------------------------------------------------------------
c              DO i = 5,5
c Setup of number of grid points 
c                  Mes = 500*(2**i) ! Mes must be an even integer
c                  if(Mes.NE.2*(Mes/2)) then
c                      write(*,*) ' No. mesh points Mes must be even !!!'
c                      Mes = Mes + 1
c                      endif  
c                  Mes1 = Mes/2
c                  gamma = DFLOAT(Mes)/DFLOAT(Mes1)
c======================================================================
c                  IF(Imap.EQ.0) THEN 
c ... Integrate using one-parameter mapping function 
c         CALL JLD(u_LJ_2n_n,Scale,C4,Rmin,Mes,Rbar,0.d0,Yl,SL_JLD)
c         CALL JLD(u_LJ_2n_n,Scale,C4,Rmin,Mes1,Rbar,0.d0,Yl2,SL_JLD2)
c                    ELSE
c ... Integrate using two-parameter mapping function 
c     Set the mapping parameter beta
                      beta = dfloat(n)/2.d0-1.d0 
         CALL JLD(u_LJ_2n_n,Scale,C4,Rmin,Mes,Rbar,beta,Yl,SL_JLD) 
c         CALL JLD(u_LJ_2n_n,Scale,C4,Rmin,Mes1,Rbar,beta,Yl2,SL_JLD2)
c                    ENDIF
c======================================================================
c        Apply Richardson's extrapolation to zero step size  h->0 
c======================================================================
c                  SL_err = (SL_JLD - SL_JLD2)/(gamma**4-1.d0)
c                  SL_JLD_RE = SL_JLD + SL_err
c======================================================================         
c           Output of the resulting a_s-values 
c======================================================================
c                  write(10,'(2I3,I8,3F18.12)') n,Imap,Mes,SL_JLD,
c     1                                        SL_JLD_RE,SL_err
c                  ENDDO
c              ENDDO
          ENDDO
      close(10)
      stop
      END
c======================================================================

c**********************************************************************
           DOUBLE PRECISION FUNCTION U_LJ_2n_n(r)
c U_LJ_2n_n is a sample user-supplied potential energy function 
c subroutine which calculates the Lennard-Jones(2n,n) potentials for 
c cases   n = 4, 5 and 6
c
c Re= 1 Angst  is the equilibrium distance 
c De is the dissociation energy    
      INTEGER*4 n
      DOUBLE PRECISION r
      DOUBLE PRECISION De,Re
      DOUBLE PRECISION t
      COMMON/LJ/De,Re,n
c=====================================================================
c     Setup of dissociation energy of LJ(2n,n) potentials 
      if(n.eq.4) De = 1000.d0 
      if(n.eq.5) De = 2165.d0 
      if(n.eq.6) De = 3761.d0   
c======================================================================
      t = (Re/r)**n
      U_LJ_2n_n = De*t*(t-2.d0)
      RETURN
      END
c======================================================================

c**********************************************************************
      SUBROUTINE JLD(U, Scale, C4, Rmin, Mes, Rbar, beta, Yl, SL)
c Subroutine JLD uses either the one- or two-parameter mapping function
c to propagate the logarithmic derivative of the modified wave function
c to the right (outer) boundary point y=1 corresponding to r->infinity
c    (a) with the one-parameter (Rbar) mapping function: 
c            y(r) = (r - Rbar)/(r + Rbar) defined on the domain [-1,1] 
c    (b) with the two-parameter (Rbar&beta) mapping function: 
c            y(r) = 2*arctan[beta*(r/Rbar-1)]/Pi   
c        defined on the finite domain [(2/Pi)arctan(-beta),1] 
c      
c If nput value of beta=0 the one-parameter mapping function (a) is used 
c Otherwise (beta.ne.0) the two-parameter mapping function (b) is used 
c         Rbar should be close to equilibrium distance Re 
c                            Rbar=Re 
c         When  beta > 0, beta should be close to the magic number n/2-1
c                          beta=-n/2-1
c======================================================================
c The user must provide a subroutine-function which calculates the 
c interaction potential U(r) defined on the interval
c                      r\in [Rmin, infinity) where
c     r    is the conventional radial coordinate and
c     Rmin is the left boundary point 
c     It is assumed that the origin of the potential is defined as 
c                         U(r->infinity)->0  
c
c The outward finite-difference propagation of the modified Ricatti 
c equation is performed by Johnson's log-derivative method: see
c        [1] B.R.Johnson, J.Comp.Phys.v.13,pp.445-449 (1973)
c        [2] B.R.Johnson, J.Chem.Phys.v.67,pp.4086-4093 (1977) 
c======================================================================
c                            Input variables 
c======================================================================    
c  U is the name of the user-supplied subroutine function that
c  calculates the potential. 
c  U must be declared in an external statement in the user calling
c  program, and it should be written as follows.
c
c  double precision function U(r)
c  double precision r
c  ----------
c  It calculates the potential at r and return the value in U
c  ----------
c  where 
c  U is the potential defined on the interval [Rmin, infinity] 
c    it is assumed that U(r->infinity)->0
c  r is the radial coordinate   
c  ----------
c  return
c  end
c----------------------------------------------------------------------
c  Scale is the scaled factor of the energy:    Scale = hbar^2/2mu
c  where mu is the reduced mass and hbar is Planck constant 
c----------------------------------------------------------------------
c  C4 is the dispersion coefficient for the case of C4/r^4 potential
c  Over wise (when n>4) C4 must be defined to be zero
c----------------------------------------------------------------------
c  Rmin is the left boundary point
c----------------------------------------------------------------------
c  Mes is the number of grid points; it must be an EVEN integer 
c----------------------------------------------------------------------
c  Rbar & beta are the parameters of the i2-parameter mapping function
c----------------------------------------------------------------------
c                      Output variables  
c  Yl is the logarithmic derivative of the modified wave function at 
c     the right (outer) boundary point y=1 (r->infinity)
c  SL is the resulting value of the scattering length
c----------------------------------------------------------------------
      INTEGER*4 Mes 
      DOUBLE PRECISION U, Scale, C4, Rmin, Rbar, beta, Yl, SL
c----------------------------------------------------------------------
      INTEGER*4 j
      DOUBLE PRECISION Pi, Pi2, Pi4, Ymin, y, r, w, h, h3, h6, Q, s,t,z
      PARAMETER ( Pi=3.1415926535897932384626433832795D0 )  
c  ************************** Initialization **************************
      IF(Rbar.LE.Rmin)THEN 
          WRITE (*,*) ' STOP Rbar<Rmin'
          STOP
          ENDIF
c----------------------------------------------------------------------
      IF(beta.EQ.0.d0)THEN
c----------------------------------------------------------------------         
c         r(y) = Rbar*(1+y)/(1-y)
c----------------------------------------------------------------------
          t = 4.d0*Rbar*Rbar
          s = -Scale*t
          Ymin = (Rmin - Rbar)/(Rmin + Rbar)
          Q = s*U(Rmin)/(1.d0-Ymin)**4
        ELSE
c---------------------------------------------------------------------- 
c         r(y) = Rbar*(1+tan(Pi*y/2)/beta)
c----------------------------------------------------------------------  
          Pi2 = Pi/2.d0
          Pi4 = Pi2**2
          s = -Scale*Pi4*(Rbar/beta)**2
          Ymin = DATAN(beta*(Rmin/Rbar-1.d0))/Pi2
          Q = Pi4 + s*U(Rmin)/(DCOS(Pi2*Ymin)**4)
        ENDIF
c======================================================================    
      h = (1.d0 - Ymin)/DFLOAT(Mes)
      h3 = h*h/3.D0
      h6 = h3/2.D0
c  *********** Estimation of initial z0-value by WKB method ************
      z = h*(DSQRT(-Q) - h*Q/3.D0)
c  ***** Outward propagation by Johnson's log-derivative method *****
      DO j = 1,Mes-1
          y = Ymin + h*DFLOAT(j)
          IF(beta.EQ.0.d0 )THEN
              r = Rbar * (1.d0 + y)/(1.d0 - y)
              Q = s * U(r) /(1.d0 - y)**4
            ELSE
              r = Rbar * (1.d0 + DTAN(Pi2*y)/beta)
              Q = Pi4 + s * U(r)/(DCOS(Pi2*y)**4)         
            ENDIF
          w = 2.d0*Q
          if(j.NE.2*(j/2)) w = 2.d0*w/(1.d0+h6*Q) 
          z = z/(1.d0 + z) - h3*w
          print*,r,z
          ENDDO
      z = z/(1.d0 + z)
c  *********** Final step at the right bound point *********************
      IF(beta.EQ.0.d0 )THEN
          Yl = (z - h3*Scale*C4/t)/h 
          SL = Rbar*(2.d0*Yl - 1.d0)
        ELSE
          Yl = (z - h3*(Pi4-s*C4))/h
          write(10,*),Yl
          SL = Rbar*(Yl/Pi2/beta + 1.d0)
        ENDIF
c  *********************************************************************
      return
c last line of subroutine JLD
      END
c=======================================================================
