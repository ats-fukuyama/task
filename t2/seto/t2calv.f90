!--------------------------------------------------------------------
!
! MODULE T2CALV
! 
!
!      CALCULATION OF PLASMA PARAMETERS  
!       FOR MULTI-FLUID TRANSPORT MODEL 
!      (ANOMALOUS TRANSPORT MODEL IS BASED ON TWO-FLUID MODEL)
!
!                       2014-05-31 H.Seto
!
!-------------------------------------------------------------------
MODULE T2CALV
  
  USE T2CNST, ONLY:ikind,rkind
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC T2CALV_EXECUTE
  
CONTAINS
  
  
  !---------------------------------------------------------
  !
  !       CALCULATION OF FUNDAMENTAL PHYSICAL QUANTITIES 
  ! 
  !
  !                     LAST UPDATE 2014-05-31 H.Seto
  !
  !
  !---------------------------------------------------------
  
  SUBROUTINE T2CALV_EXECUTE(i_m)
    
    INTEGER,INTENT(IN)::i_m
    
    CALL T2CALV_SETUP(i_m)
    
    CALL T2CALV_FRICTION_COEFFICIENT  ! L11, L12, L21, L22

    CALL T2CALV_VISCOUS_COEFFICIENT   ! Mu1, Mu2, Mu3
    
    CALL T2CALV_ANOMALOUS_COEFFICIENT !
    
    CALL T2CALV_MODIFY_COEFFICIENT    ! Lx1 ~ Lx4, Mux1 ~ Mux4
    RETURN

  END SUBROUTINE T2CALV_EXECUTE
  
  SUBROUTINE T2CALV_SETUP(i_m)
    
    USE T2CNST,ONLY: Amp,Aee,Pi,Eps0
    
    USE T2COMM,ONLY:&
         & NSMAX,BpNF,BtNF,&
         & NnNF,FrNF,FbNF,FtNF,FpNF,&
         & PpNF,QrNF,QbNF,QtNF,QpNF,&
         & i2crt,d2rzm,d2mfc1,d2jm1,&
         !
         & XvecIn,&
         & Pa, Pz, R_rz,  R_mc,&
         & GRt  , G11xCo, G12Co, G22Co,  G33Co,&
         &        G11Ct,  G12Ct, G22xCt, G33Ct,&  
         & UgrCt, UgpCt,&
         & BtCo, BtCt, BtSq, BpCo, BpCt, BpSq, Bb, BbSq,&
         & Mm, Ee,Vv,BaseNu,X,Y,Z,&
         & Nn, FrCt, Fb, FtCo, FpCt, UrCt, Ub, UtCo, UpCt, UuSq,&
         & Pp, QrCt, Qb, QtCo, QpCt, WrCt, Wb, WtCo, WpCt, Tt
    
    
    INTEGER(ikind),INTENT(IN)::i_m
    INTEGER(ikind)::&
         i_m1d,i_m2d,i_s,j_s,vOffset
    REAL(   rkind),DIMENSION(1:NSMAX,1:NSMAX)::&
         lnLamb
    REAL(   rkind)::&
         & nnA,frCtA,fbA,ftCoA,fpCtA,urCtA,ubA,utCoA,upCtA,&
         & ppA,qrCtA,qbA,qtCoA,qpCtA,wrCtA,wbA,wtCoA,wpCtA,&
         !
         & ttA,ttB,mmA,mmB,vvA,vvB,eeA,eeB,nnB,&
         & nnE,ttE_eV,dPsiDr,lnLambAB
    
    ! R_rz   : Local major radisu: R
    ! R_mc   : radial flux label : r, \rho, \sigma etc...
    ! GRt    : \sqrt{g} 
    ! G11xCo : g_{\sig\sig}*R_mc
    ! G12Co  : g_{\sig\chi}
    ! G22Co  : g_{\chi\chi}
    ! G33Co  : g_{\zeta\zeta}
    !
    ! G11Ct  : g^{\sig\sig}
    ! G12Ct  : g^{\sig\chi}
    ! G22xCt : g^{\chi\chi}*R_mc
    ! G33Ct  : g^{\zeta\zeta} 
    !
    ! UgrCt  : CT RADIAL   FLAME MOVING VELOCITY: default 0
    ! UgrCt  : CT POLOIDAL FLAME MOVING VELOCITY: default 0
    
    ! BpCt   : CT POLOIDAL MAGNETIC FIELD          : B^{\chi}
    ! BpCo   : Co POLOIDAL MAGNETIC FIELD          : B_{\chi}
    ! BpSq   : SQUARED INTENSITY 
    !                  of POLOIDAL MAGENTIC FIELD  : Bp^{2}
    !
    ! BtCt   : CT Toroidal magnetic field          : B^{\zeta}
    ! BtCo   : CO Toroidal magnetic field          : B_{\zeta}
    ! BtSq   : SQUARED INTENSITY 
    !                  of TOROIDAL MAGENTIC FIELD  : Bt^{2}
    ! Bb     : INTENSITY OF MAGNETIC FIELD         : B 
    ! BbSq   : SQUARED INTENSITY of MAGNETIC FIELD : B^{2} 
    
    ! initialization
    i_m1d = i2crt( 3,i_m)
    i_m2d = i2crt( 2,i_m)
    R_rz  = d2rzm( 1,i_m)
    R_mc  = d2mfc1(1,i_m)
        
    GRt      = d2jm1(1,i_m)
    G11xCo   = d2jm1(2,i_m)
    G12Co    = d2jm1(3,i_m)
    G22Co    = d2jm1(4,i_m)
    G33Co    = d2jm1(5,i_m)
    G11Ct    = d2jm1(6,i_m)
    G12Ct    = d2jm1(7,i_m)
    G22xCt   = d2jm1(8,i_m)
    G33Ct    = d2jm1(9,i_m)
    
    UgrCt  = 0.D0
    UgpCt  = 0.D0

    DO i_s = 1,NSMAX
       Mm(i_s) = Pa(i_s)*Amp ! mass 
       Ee(i_s) = Pz(i_s)*Aee ! electric charge
    ENDDO
    
    dPsiDr = BpNF*XvecIn(1,i_m1d)
    
    BtCo = BtNF*XvecIn(2,i_m1d)
    BtCt = BtCo*G33Ct
    
    BpCt = dPsiDr/GRt
    BpCo = BpCt*G22Co
    
    UgrCt = 0.D0
    UgpCt = 0.D0
    
    BpSq = bpCo*bpCt
    BtSq = btCo*btCt
    
    BbSq = BpSq + BtSq
    Bb   = SQRT(BbSq)
    
    DO i_s = 1, NSMAX
       
       vOffset =  10*(i_s - 1) + 5
       
       nnA    = NnNF*XvecIn(vOffset+ 1,i_m2d)
       ppA    = PpNF*XvecIn(vOffset+ 6,i_m2d)
       IF((nnA.LT.0.D0).OR.(ppA.LT.0.D0))THEN
          WRITE(6,*)'NEGATIVE  DENSITY or PRESSURE'
          WRITE(6,*)'SPECIS=',i_s,'NODE=',i_m2d,&
               'N=',nnA*1.D-20,'1.D20/m3',&
               'P=',ppA*1.D-23/Aee,'keV*1.D20/m3'
          STOP       
       ENDIF
       
       frCtA = FrNF*XvecIn(vOffset+ 2,i_m2d)*R_mc
       fbA   = FbNF*XvecIn(vOffset+ 3,i_m2d)
       ftCoA = FtNF*XvecIn(vOffset+ 4,i_m2d)
       fpCtA = FpNF*XvecIn(vOffset+ 5,i_m2d)
       
       urCtA = frCtA/nnA
       ubA   = fbA  /nnA
       utCoA = ftCoA/nnA
       upCtA = fpCtA/nnA

       qrCtA = QrNF*XvecIn(vOffset+ 7,i_m2d)*R_mc
       qbA   = QbNF*XvecIn(vOffset+ 8,i_m2d)
       qtCoA = QtNF*XvecIn(vOffset+ 9,i_m2d)
       qpCtA = QpNF*XvecIn(vOffset+10,i_m2d)
       
       wrCtA = qrCtA/ppA
       wbA   = qbA  /ppA
       wtCoA = qtCoA/ppA
       wpCtA = qpCtA/ppA
       
       Nn(  i_s) = nnA
       FrCt(i_s) = frCtA
       Fb(  i_s) = fbA
       FtCo(i_s) = ftCoA
       FpCt(i_s) = fpCtA

       UrCt(i_s) = urCtA
       Ub(  i_s) = ubA
       UtCo(i_s) = utCoA
       UpCt(i_s) = upCtA
       UuSq(i_s) = upCtA*upCtA*G22Co + utCoA*utCoA*G33Ct
       
       Pp(  i_s) = ppA
       QrCt(i_s) = qrCtA
       Qb(  i_s) = qbA
       QtCo(i_s) = qtCoA
       QpCt(i_s) = qpCtA

       WrCt(i_s) = wrCtA
       Wb(  i_s) = wbA
       WtCo(i_s) = wtCoA
       WpCt(i_s) = wpCtA

       Tt(i_s) = ppA/nnA

    ENDDO
    
    ! THERMAL VELOCITY [m/s]
    ! v_{ta} = \sqrt{2T_{a}/M_{a}}
    DO i_s = 1, NSMAX
       ttA = Tt(i_s)
       mmA = Mm(i_s)
       Vv(i_s) = SQRT(2.0D0*ttA/mmA)
    ENDDO
    
    ! COULOMB LOGARITHM lnLamb 
    ! ref: COLLISIONAL TRANSPORT IN MAGNETIZED PLASMAS
    !      P. HELANDER AND D.J. SIGMAR (2002)  P.4
    nnE    = Nn(1)
    ttE_eV = Tt(1)/Aee
    DO j_s = 1, NSMAX
    DO i_s = 1, NSMAX
       lnLamb(i_s,j_s) = 18.4D0 -1.15D0*LOG10(nnE)+2.3D0*LOG10(ttE_eV)
    ENDDO
    ENDDO
    
    !  basic collision frequency [Hz]
    !  ref: COLLISIONAL TRANSPORT IN MAGNETIZED PLASMAS
    !      P. HELANDER AND D.J. SIGMAR (2002)  P.277
    DO j_s = 1, NSMAX
       eeB   = Ee(j_s)
       nnB   = Nn(j_s)
       DO i_s = 1, NSMAX
          eeA = Ee(i_s)
          mmA = Mm(i_s)
          vvA = Vv(i_s)
          lnLambAB = lnLamb(i_s,j_s)
          BaseNu(i_s,j_s) = (nnB*(eeA**2)*(eeB**2)*lnLambAB) &
               &          / (4.D0*Pi*(Eps0**2)*(vvA**3)*(mmA**2))
       ENDDO
    ENDDO
    
    ! DEFINITION OF DIMENSIONLESS PARAMETERS 
    !
    ! x(i_s,j_s) = x_ab = v_{tb}/v_{ta}
    ! y(i_s,j_s) = y_ab = m_{ a}/m_{ b}
    ! z(i_s,j_s) = z_ab = t_{ a}/t_{ b}

    DO j_s = 1, NSMAX
       vvB = Vv(j_s)
       mmB = Mm(j_s)
       ttB = Tt(j_s)
       DO i_s = 1, NSMAX 
          vvA = Vv(i_s)
          mmA = Mm(i_s)
          ttA = Tt(i_s)
          X(i_s,j_s) = vvB/vvA
          Y(i_s,j_s) = mmA/mmB
          Z(i_s,j_s) = ttA/ttB
       ENDDO
    ENDDO
   
!    ! set known variables for T2EXEC: KnownVar
!    !
!    ! KnownVar(1   ,:): MAGNETIC FIELD INTENSITY      : B   
!    ! KnownVar(2   ,:): MAJOR RADIUS                  : R   
!    ! KnownVar(2N+1,:): NORMALISED PARALLEL FLOW      : u_{a\para}
!    ! KnownVar(2N+2,:): NORMALIZED PARALLEL HEAT FLOW : w_{a\para}
!    
!    KnownVar(1,i_m) = bb
!    KnownVar(2,i_m) = r_rz
!    
!    DO i_s = 0, NSMAX-1
!       i_w = 2*i_s + 2
!       KnownVar(i_w+1,i_m) = ub(i_s)
!       KnownVar(i_w+2,i_m) = wb(i_s)
!    ENDDO

    RETURN

  END SUBROUTINE T2CALV_SETUP

  ! BRAGINSKII'S MATRIX ELEMENT OF COLLISION OPERATOR
  ! REF: COLLISIONAL TRANSPORT IN MAGNETIZED PLASMAS
  !      P. HELANDER AND D.J. SIGMAR (2002)
  !      P.239, P.278
  !
  ! temp = 1/SQRT(1 + x_{ab}^{2})
  !
  ! DEFINITION OF VARIABLRS FOR BRAGINSKII MATRIX ELEMENTS
  ! 
  ! m00[i_s,j_s] : M_{ab}^{00} [ND]
  ! m01[i_s,j_s] : M_{ab}^{01} [ND]
  ! m10[i_s,j_s] : M_{ab}^{10} [ND]
  ! m11[i_s,j_s] : M_{ab}^{11} [ND]
  ! n00[i_s,j_s] : N_{ab}^{00} [ND]
  ! n01[i_s,j_s] : N_{ab}^{01} [ND]
  ! n10[i_s,j_s] : N_{ab}^{10} [ND]
  ! n11[i_s,j_s] : N_{ab}^{11} [ND]
  ! l11[i_s,j_s] : l^{ab}_{11} [kg/m2*s2]
  ! l12[i_s,j_s] : l^{ab}_{12} [kg/m2*s2]
  ! l21[i_s,j_s] : l^{ab}_{21} [kg/m2*s2]
  ! l22[i_s,j_s] : l^{ab}_{22} [kg/m2*s2] 
  SUBROUTINE T2CALV_FRICTION_COEFFICIENT
    
    USE T2CNST,ONLY: Eps0,Pi
    
    USE T2COMM,ONLY:&
         & NSMAX, R_rz, R_mc, BaseNu, X, Y, Z, Mm, Nn,&
         & L11,L12,L21,L22,Hex
    
    INTEGER(ikind)::i_s,j_s,k_s
    REAL(   rkind)::&
         & temp,xAB,yAB,zAB,zBA,&
         & nuAB,m00AB,m01AB,m10AB,m11AB, &
         & nuAC,m00AC,m01AC,m10AC,m11AC, &
         &      n00AB,n01AB,n10AB,n11AB, &
         &      l11AB,l12AB,l21AB,l22AB
    
    REAL(   rkind),DIMENSION(1:NSMAX,1:NSMAX)::&
         nu,m00,m01,m10,m11,n00,n01,n10,n11
    
    ! COLLISION FREQUENCY [1/s]
    ! REF: COLLISIONAL TRANSPORT IN MAGNETIZED PLASMAS
    !      P. HELANDER AND D.J. SIGMAR (2002)  P.277
    DO j_s = 1, NSMAX
    DO i_s = 1, NSMAX
       nu(i_s,j_s) = BaseNu(i_s,j_s)/(0.75D0*SQRT(Pi))
    ENDDO
    ENDDO
    
    ! Braginsikii matrix elements: M^{ab]_{ij}, N^{ab}_{ij}
    DO j_s = 1, NSMAX
    DO i_s = 1, NSMAX
       
       xAB  = X(i_s,j_s) 
       yAB  = Y(i_s,j_s)
       zAB  = Z(i_s,j_s)
       
       temp = 1.D0/SQRT(1.D0+xAB**2)
       m00AB = -(temp**3)*(1.D0+yAB)
       m01AB = -(temp**5)*(1.D0+yAB)*1.5D0
       m11AB = -(temp**5)*(3.25D0+4.D0*(xAB**2)+7.50D0*(xAB**4))
       
       m00(i_s,j_s) =  m00AB
       m01(i_s,j_s) =  m01AB
       m10(i_s,j_s) =  m01AB
       m11(i_s,j_s) =  m11AB
       
       n00(i_s,j_s) = -m00AB
       n01(j_s,i_s) = -m01AB*xAB/zAB
       n10(i_s,j_s) = -m01AB
       n11(i_s,j_s) = (temp**5)*6.75D0*zAB*(xAB**2)
       
    ENDDO
    ENDDO
    
    ! PARALLEL FRICTION COEFFICIENTS: l^{ab]_{ij}
    l11AB = 0.D0
    l12AB = 0.D0
    l21AB = 0.D0
    l22AB = 0.D0
    
    DO j_s = 1, NSMAX
    DO i_s = 1, NSMAX
       nuAB  = nu( i_s,j_s)

       n00AB = n00(i_s,j_s)
       n01AB = n01(i_s,j_s)
       n10AB = n10(i_s,j_s)
       n11AB = n11(i_s,j_s)
       
       l11AB = n00AB*nuAB
       l12AB = n01AB*nuAB
       l21AB = n10AB*nuAB
       l22AB = n11AB*nuAB
       
       IF(i_s.EQ.j_s)THEN
          DO k_s = 1, NSMAX
             nuAC  = nu( i_s,k_s)
             
             m00AC = m00(i_s,k_s)
             m01AC = m01(i_s,k_s)
             m10AC = m10(i_s,k_s)
             m11AC = m11(i_s,k_s)
             
             l11AB = l11AB + m00AC*nuAC
             l12AB = l12AB + m01AC*nuAC
             l21AB = l21AB + m10AC*nuAC
             l22AB = l22AB + m11AC*nuAC
          ENDDO
       ENDIF

       temp = Mm(i_s)*Nn(i_s)
       L11(i_s,j_s) = l11AB*temp ! l_{11}^{ab}
       L12(i_s,j_s) = l12AB*temp ! l_{12}^{ab}
       L21(i_s,j_s) = l21AB*temp ! l_{21}^{ab}
       L22(i_s,j_s) = l22AB*temp ! l_{22}^{ab}
    ENDDO
    ENDDO       
    
    ! HEAT EXCHANGE COEFFICIENT
    DO i_s = 1, NSMAX
       Hex(i_s) = 0.D0
       DO j_s = 1, NSMAX
          zBA  = Z( j_s,i_s)
          nuAB = nu(i_s,j_s)
          Hex(i_s) = Hex(i_s) + 1.5D0*(1.D0 - zBA)*nuAB
       ENDDO
    ENDDO
        
    RETURN
    
  END SUBROUTINE T2CALV_FRICTION_COEFFICIENT

  ! NEOCLASSICAL PARALLEL COEFFICIENTS
  ! Ref: COLLISIONAL TRANSPORT IN MAGNETIZED PLASMAS
  !      P. HELANDER AND D.J. SIGMAR (2002)  CHAP.12 AND P.279
  !
  !      trapped particle fraction rate is defined by
  !
  !      f_t = 1.46*(r/R0)^{1/2} - 0.46*(r/R0)^{3/2}
  !      
  ! Ref: Y.B. Kim et al., Phys. Fluids B 3 (1991) 2050
  !
  ! In TASK/T2 Kij is defined as follows 
  !
  ! On magnetic axis and outside LCFS
  !
  !            Kij = Kij_PS  (1)
  !
  ! inside LCFS except for magnetic axis
  !
  !            Kij = Kij_PS*Kij_BN/(Kij_PS+Kij_BN) (2)
  ! 
  ! (2) become (1) in the limit where r -> 0 or  r -> 1, 
  ! since Kij_BN -> +infty and K_ij = Kij_PS
  !
  ! iar  : INVERSE ASPECT RATIO: r/R0
  ! iarRt: Square root of inverse aspect ratio: SQRT(r/R0)  
  ! f_t  : trapped particle fraction rate
  ! q    : safety factor (in TASK/T2 q is defined 
  !        by local magnetic slope
  !                                   q =  B^{\zeta}/B^{\chi}
  !
  ! k11ps: MODIFIED VISCOSITY COEFFICIENT IN PS REGIME: K^{11}_{aPS}
  ! k12ps: MODIFIED VISCOSITY COEFFICIENT IN PS REGIME: K^{12}_{aPS}
  ! k22ps: MODIFIED VISCOSITY COEFFICIENT IN PS REGIME: K^{22}_{aPS}
  !
  ! k11bn: MODIFIED VISCOSITY COEFFICIENT IN BN REGIME: K^{11}_{aBN}
  ! k12bn: MODIFIED VISCOSITY COEFFICIENT IN BN REGIME: K^{12}_{aBN}
  ! k22bn: MODIFIED VISCOSITY COEFFICIENT IN BN REGIME: K^{22}_{aBN}
  !
  ! k11  : MODIFIED VISCOSITY COEFFICIENT             : K^{11}_{a}
  ! k12  : MODIFIED VISCOSITY COEFFICIENT             : K^{12}_{a}
  ! k22  : MODIFIED VISCOSITY COEFFICIENT             : K^{22}_{a}
  !
  ! mu1 : NEOCLASSICAL VISCOSITY COEFFICIENT          : mu_{1a}
  ! mu2 : NEOCLASSICAL VISCOSITY COEFFICIENT          : mu_{2a} 
  ! mu3 : NEOCLASSICAL VISCOSITY COEFFICIENT          : mu_{3a}
  SUBROUTINE T2CALV_VISCOUS_COEFFICIENT

    USE T2CNST,ONLY: DeCoef
    USE T2COMM,ONLY: NSMAX,RA,RR,R_mc,R_rz,I_xa,&
         & Mm,Nn,Pp,BtCt,BpCt,&
         & Mu1,Mu2,Mu3

    USE LIBT2, ONLY:integ_f
    

    REAL(   rkind)::&
         & iar,iarRt,f_t,f_c,q,temp,&
         & mmA,nnA,ppA,&
         & psCoef,k11ps,k12ps,k22ps,&
         & bnCoef,k11bn,k12bn,k22bn,&
         &        k11,  k12,  k22,derr
    
    IF((R_mc.GT.0.D0).OR.(R_mc.LT.1.D0))THEN
       
       iar   = RA*SQRT(R_mc)/RR
       iarRt = SQRT(iar)
       f_t   = 1.46D0*iarRt - 0.46D0*(iarRt**3)
       f_c   = 1.D0 - f_t
       q     = btCt/bpCt 
       temp  = (2.D0*(q**2)*(r_rz**2)*f_t*DeCoef) /(3.D0*(iar**2)*f_c)
       
       DO I_xa = 1, NSMAX
          mmA = mm(I_xa)
          nnA = nn(I_xa)
          ppA = pp(I_xa)
          
          !C MODIFIED VISCOSITY COEFFICIENT IN PS REGIME
          CALL INTEG_F(FUNC_k11ps,k11ps,derr,EPS=1.D-8,ILST=0)
          CALL INTEG_F(FUNC_k12ps,k12ps,derr,EPS=1.D-8,ILST=0)
          CALL INTEG_F(FUNC_k22ps,k22ps,derr,EPS=1.D-8,ILST=0)
          
          psCoef = 0.4D0*ppA*DeCoef
          
          k11ps = psCoef*k11ps
          k12ps = psCoef*k12ps
          k22ps = psCoef*k22ps
          
          ! MODIFIED VISCOSITY COEFFICIENT IN BN REGIME
          CALL INTEG_F(FUNC_k11bn,k11bn,derr,EPS=1.D-8,ILST=0)
          CALL INTEG_F(FUNC_k12bn,k12bn,derr,EPS=1.D-8,ILST=0)
          CALL INTEG_F(FUNC_k22bn,k22bn,derr,EPS=1.D-8,ILST=0)
          
          bnCoef = mmA*nnA*temp
          
          k11bn = bnCoef*k11bn
          k12bn = bnCoef*k12bn
          k22bn = bnCoef*k22bn
          
          ! MODIFIED VISCOSITY COEFFICIENT IN INTER-REGIMES 
          k11 = k11bn*k11ps/(k11bn + k11ps)
          k12 = k12bn*k12ps/(k12bn + k12ps)
          k22 = k22bn*k22ps/(k22bn + k22ps)
          
          ! NEOCLASSICAL VISCOSITY COEFFICIENTS: \mu_{ai}
          Mu1(I_xa) = k11
          Mu2(I_xa) = k12 - 2.5D0*k11
          Mu3(I_xa) = k22 - 5.0D0*k12 + 6.25D0*k11
       ENDDO
    ELSE       
       DO I_xa = 1, NSMAX
          
          ppA = pp(I_xa)
          
          !C MODIFIED VISCOSITY COEFFICIENT IN PS REGIME
          CALL INTEG_F(FUNC_k11ps,k11ps,derr,EPS=1.D-8,ILST=0)
          CALL INTEG_F(FUNC_k12ps,k12ps,derr,EPS=1.D-8,ILST=0)
          CALL INTEG_F(FUNC_k22ps,k22ps,derr,EPS=1.D-8,ILST=0)
          
          psCoef = 0.4D0*ppA*DeCoef
          
          k11ps = psCoef*k11ps
          k12ps = psCoef*k12ps
          k22ps = psCoef*k22ps
          
          ! MODIFIED VISCOSITY COEFFICIENT IN INTER-REGIMES 
          k11 = k11ps
          k12 = k12ps
          k22 = k22ps
          
          ! NEOCLASSICAL VISCOSITY COEFFICIENTS: \mu_{ai}
          Mu1(I_xa) = k11
          Mu2(I_xa) = k12 - 2.5D0*k11
          Mu3(I_xa) = k22 - 5.0D0*k12 + 6.25D0*k11

       ENDDO
    ENDIF
    
    RETURN
    
  END SUBROUTINE T2CALV_VISCOUS_COEFFICIENT

! CALCULATE K11_PS 
  !        BY DOUBLE EXPONENTIAL INTEGRATION FORMULA
  FUNCTION FUNC_k11ps(xaSq)
 
    USE LIBT2, ONLY:func_phi,func_G
    USE T2COMM,ONLY:NSMAX,I_xa,X,Y,Z,BaseNu
    INTEGER(ikind)::i_xb
    REAL(   rkind),INTENT(IN)::xaSq
    REAL(   rkind)::FUNC_k11ps
    REAL(   rkind)::&
         xA,xB,baseNuAB,xBA,yBA,zAB,&
         nuS,nuB,nuD,nuT
    
    xA = SQRT(xASq) ! V/Vt_a
    
    nuT = 0.D0
    
    DO i_xb  = 1, NSMAX
       xBA      = X(  i_xb,I_xa)    ! Vt_a/Vt_b
       yBA      = Y(  i_xb,I_xa)    ! m_b/m_a
       zAB      = Z(  I_xa,i_xb)    ! T_a/T_b
       baseNuAB = BaseNu(I_xa,i_xb) 
       xB  = xA*xBA            ! (V/Vt_a) * (Vt_a/Vt_b)
       nuD =      baseNuAB*(func_phi(xB)-func_G(xB))/(xA**3)
       nuS = 2.D0*baseNuAB*zAB*(1.D0+yBA)*func_G(xB)/xA 
       nuB = 2.D0*baseNuAB*func_G(xB)/(Xa**3)
       
       nuT = nuT + 2.D0*nuS - nuB + nuD
    ENDDO
    
    IF(    nuT.GT.0.D0)THEN 
       nuT = 1.D0/nuT
    !ELSEIF(nuT.EQ.0.D0)THEN
    !   d0nut = 0.D0
    ELSE
       WRITE(6,*)'ERROR IN FUNC_K11PS'
       STOP
    END IF
    
    FUNC_K11PS = EXP(-xA**2)*(xA**5)*nuT
    
    RETURN
    
  END FUNCTION FUNC_k11ps
  
  ! CALCULATE K12_PS 
  !        BY DOUBLE EXPONENTIAL INTEGRATION FORMULA
  FUNCTION FUNC_k12ps(xaSq)
 
    USE LIBT2, ONLY:func_phi,func_G
    USE T2COMM,ONLY:NSMAX,I_xa,X,Y,Z,BaseNu
    INTEGER(ikind)::i_xb
    REAL(   rkind),INTENT(IN)::xaSq
    REAL(   rkind)::FUNC_k12ps
    REAL(   rkind)::&
         xA,xB,baseNuAB,xBA,yBA,zAB,&
         nuS,nuB,nuD,nuT
    
    xA = SQRT(xASq) ! V/Vt_a
    
    nuT = 0.D0
    
    DO i_xb  = 1, NSMAX
       xBA      = X(  i_xb,I_xa)    ! Vt_a/Vt_b
       yBA      = Y(  i_xb,I_xa)    ! m_b/m_a
       zAB      = Z(  I_xa,i_xb)    ! T_a/T_b
       BaseNuAB = baseNu(I_xa,i_xb) 
       xB  = xA*xBA            ! (V/Vt_a) * (Vt_a/Vt_b)
       nuD =      baseNuAB*(func_phi(xB)-func_G(xB))/(xA**3)
       nuS = 2.D0*baseNuAB*zAB*(1.D0+yBA)*func_G(xB)/xA 
       nuB = 2.D0*baseNuAB*func_G(xB)/(Xa**3)
       nuT = nuT + 2.D0*nuS - nuB + nuD
    ENDDO
    
    IF(    nuT.GT.0.D0)THEN 
       nuT = 1.D0/nuT
    !ELSEIF(nuT.EQ.0.D0)THEN
    !   d0nut = 0.D0
    ELSE
       WRITE(6,*)'ERROR IN FUNC_K11PS'
       STOP
    END IF
    
    FUNC_K12PS = EXP(-xA**2)*(xA**7)*nuT
    
    RETURN
    
  END FUNCTION FUNC_k12ps
  
  ! CALCULATE K22_PS 
  !        BY DOUBLE EXPONENTIAL INTEGRATION FORMULA
  FUNCTION FUNC_k22ps(xaSq)
    
    USE LIBT2, ONLY:func_phi,func_G
    USE T2COMM,ONLY:NSMAX,I_xa,X,Y,Z,BaseNu
    INTEGER(ikind)::i_xb
    REAL(   rkind),INTENT(IN)::xaSq
    REAL(   rkind)::FUNC_k22ps
    REAL(   rkind)::&
         xA,xB,baseNuAB,xBA,yBA,zAB,&
         nuS,nuB,nuD,nuT
    
    xA = SQRT(xASq) ! V/Vt_a
    
    nuT = 0.D0
    
    DO i_xb  = 1, NSMAX
       xBA      = X(  i_xb,I_xa)    ! Vt_a/Vt_b
       yBA      = Y(  i_xb,I_xa)    ! m_b/m_a
       zAB      = Z(  I_xa,i_xb)    ! T_a/T_b
       baseNuAB = BaseNu(I_xa,i_xb) 
       xB  = xA*xBA            ! (V/Vt_a) * (Vt_a/Vt_b)
       nuD =      baseNuAB*(func_phi(xB)-func_G(xB))/(xA**3)
       nuS = 2.D0*baseNuAB*zAB*(1.D0+yBA)*func_G(xB)/xA 
       nuB = 2.D0*baseNuAB*func_G(xB)/(Xa**3)
       nuT = nuT + 2.D0*nuS - nuB + nuD
    ENDDO
    
    IF(    nuT.GT.0.D0)THEN 
       nuT = 1.D0/nuT
       !ELSEIF(nuT.EQ.0.D0)THEN
    !   d0nut = 0.D0
    ELSE
       WRITE(6,*)'ERROR IN FUNC_K11PS'
       STOP
    END IF
    
    FUNC_K22PS = EXP(-xA**2)*(xA**9)*nuT
    
    RETURN
    
  END FUNCTION FUNC_k22ps

  ! CALCULATE K11_BN 
  !        BY DOUBLE EXPONENTIAL INTEGRATION FORMULA
  FUNCTION FUNC_k11bn(xaSq)
        
    USE LIBT2, ONLY:func_phi,func_G
    USE T2COMM,ONLY:NSMAX,I_xa,X,BaseNu
    INTEGER(ikind)::i_xb
    REAL(   rkind),INTENT(IN)::xaSq
    REAL(   rkind)::FUNC_k11bn
    REAL(   rkind)::&
         xA,xB,baseNuAB,xBA,nuD,nuT
    
    xA = SQRT(XaSq) !C V/Vt_a

    nuT = 0.D0

    DO i_xb  = 1, NSMAX
       
       xBA      = X(  i_xb,I_xa)  ! Vt_a/Vt_b
       baseNuAB = BaseNu(I_xa,i_xb)
       xB       = xA*xBA          ! V/Vt_b = (V/Vt_a)*(Vt_a/Vt_b)
       nuD      = baseNuAB*(func_phi(xB)-func_G(xB))/(xA**3)
       nuT      = nuT + nuD
    ENDDO
    
    FUNC_k11bn = EXP(-xA**2)*(xA**3)*nuT
    
    RETURN
    
  END FUNCTION FUNC_k11bn
  
  ! CALCULATE K12_BN 
  !        BY DOUBLE EXPONENTIAL INTEGRATION FORMULA
  FUNCTION FUNC_k12bn(xaSq)
        
    USE LIBT2, ONLY:func_phi,func_G
    USE T2COMM,ONLY:NSMAX,I_xa,X,BaseNu
    INTEGER(ikind)::i_xb
    REAL(   rkind),INTENT(IN)::xaSq
    REAL(   rkind)::FUNC_k12bn
    REAL(   rkind)::&
         xA,xB,baseNuAB,xBA,nuD,nuT
    
    xA = SQRT(XaSq) !C V/Vt_a

    nuT = 0.D0

    DO i_xb  = 1, NSMAX
       
       xBA      = X(  i_xb,I_xa)  ! Vt_a/Vt_b
       baseNuAB = BaseNu(I_xa,i_xb)
       xB       = xA*xBA          ! V/Vt_b = (V/Vt_a)*(Vt_a/Vt_b)
       nuD      = baseNuAB*(func_phi(xB)-func_G(xB))/(xA**3)
       nuT      = nuT + nuD
    ENDDO
    
    FUNC_k12bn = EXP(-xA**2)*(xA**5)*nuT
    
    RETURN
    
  END FUNCTION FUNC_k12bn


  ! CALCULATE K22_BN 
  !        BY DOUBLE EXPONENTIAL INTEGRATION FORMULA
  FUNCTION FUNC_k22bn(xaSq)
    
    USE LIBT2, ONLY:func_phi,func_G
    USE T2COMM,ONLY:NSMAX,I_xa,X,BaseNu
    INTEGER(ikind)::i_xb
    REAL(   rkind),INTENT(IN)::xaSq
    REAL(   rkind)::FUNC_k22bn
    REAL(   rkind)::&
         xA,xB,baseNuAB,xBA,nuD,nuT
    
    xA = SQRT(XaSq) !C V/Vt_a

    nuT = 0.D0

    DO i_xb  = 1, NSMAX
       
       xBA      = X(     i_xb,I_xa)  ! Vt_a/Vt_b
       baseNuAB = BaseNu(I_xa,i_xb)
       xB       = xA*xBA          ! V/Vt_b = (V/Vt_a)*(Vt_a/Vt_b)
       nuD      = baseNuAB*(func_phi(xB)-func_G(xB))/(xA**3)
       nuT      = nuT + nuD
    ENDDO
    
    FUNC_k22bn = EXP(-xA**2)*(xA**7)*nuT
    
    RETURN
    
  END FUNCTION FUNC_k22bn

  SUBROUTINE T2CALV_ANOMALOUS_COEFFICIENT
    
    USE T2CNST,ONLY:Aee
    
    USE T2COMM,ONLY:&
         & NSMAX,R_mc,Mm,Nn,Tt,BpCt,BtCo,BbSq,Bb,GRt,G11Ct, &
         & FtAnom1,FtAnom2,FtAnom3,FtAnom4,FtAnom5,  &
         & GtAnom1,GtAnom2,GtAnom3,GtAnom4,GtAnom5

    INTEGER(ikind)::&
         i_s
    REAL(   rkind)::&
         ttE,&
         mmA,eeA,eeB,vvA,&
         d_anom,m_anom,x_anom,temp,&
         ftAnom1E,ftAnom2E,ftAnom3E,ftAnom4E
    
    !
    ! COEFFICIENTS OF ANORMALOUS TRANSPORT BY QUASI-LINEAR THEORY
    ! REF: S.I. Itoh, Phys. Fluids B, 4 796         (1992)
    !      K.C. Shaing, Phys. Plasmas, 31 2249      (1988)
    !      M. Honda et. al, Nucl. Fusion, 50 095012 (2010)
    !
    IF(R_mc.LE.1.D0)THEN
       d_anom = 0.45D0*R_mc+0.05D0
       m_anom = 0.45D1*R_mc+0.05D1
       x_anom = 0.45D1*R_mc+0.05D1
    ELSE
       d_anom = 0.5D0
       m_anom = 0.5D1
       x_anom = 0.5D1 
    ENDIF
    
    ! Toroidal Force by Anomalous Transport (1)
    ttE = Tt(1)
    temp = (Aee**2)*d_anom/ttE
    
    ftAnom1E = -temp*(1.5D0-m_anom/d_anom)*GRt*BpCt*G11Ct*ttE/Aee !ne
    ftAnom2E =  temp*(1.5D0-m_anom/d_anom)*GRt*BpCt*G11Ct    /Aee !pe
    ftAnom3E =  temp*BtCo*Bb !Gamma_{e\para}
    ftAnom4E = -temp*BbSq    !Gamma_{e\zeta}
    
    FtAnom1(1) = ftAnom1E
    FtAnom2(1) = ftAnom2E
    FtAnom3(1) = ftAnom3E
    FtAnom4(1) = ftAnom4E

    DO i_s =2, NSMAX
       FtAnom1(i_s) = -ftAnom1E
       FtAnom2(i_s) = -ftAnom2E
       FtAnom3(i_s) = -ftAnom3E
       FtAnom4(i_s) = -ftAnom4E
    ENDDO

    ! Toroidal Force by Anomalous Transport (2)
    DO i_s = 1, NSMAX
       FtAnom5(i_s) = Mm(i_s)*Nn(i_s)*m_anom
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2CALV_ANOMALOUS_COEFFICIENT


  SUBROUTINE T2CALV_MODIFY_COEFFICIENT

    ! PARALLEL FRICTION COEFFICIENTS 
    !          WITH RESPECT TO MOMENTUM AND TOTAL HEAT FLUX
    !
    ! DEFINITION OF VARIABLRS 
    !            FOR NEOCLASSICAL FRICTION COEFFICIENTS
    !
    ! l01: \bar{l}^{ab}_{01} [kg/m2*s2/m3]
    ! l02: \bar{l}^{ab}_{02} [kg/m2*s2/m3]
    ! l03: \bar{l}^{ab}_{03} [kg/m2*s2/J ]
    ! l04: \bar{l}^{ab}_{04} [kg/m2*s2/J ]

    USE T2COMM,ONLY:&
         & NSMAX,L11,L12,L21,L22,Mu1, Mu2, Mu3,&
         &       Lx1,Lx2,Lx3,Lx4,Mux1,Mux2,Mux3,Mux4

    INTEGER(ikind)::i_s,j_s
    REAL(   rkind)::&
         & l11AB,l12AB,l21AB,l22AB,lx1AB,lx2AB,lx3AB,lx4AB,&
         & mu1A, mu2A, mu3A,       mux1A,mux2A,mux3A,mux4A

    DO j_s = 1, NSMAX
    DO i_s = 1, NSMAX
       
       l11AB = L11(i_s,j_s)
       l12AB = L12(i_s,j_s)
       l21AB = L21(i_s,j_s)
       l22AB = L22(i_s,j_s)
       
       lx1AB =  l11AB + l12AB
       lx2AB = -0.4D0 * l12AB
       lx3AB =  l21AB + l22AB
       lx4AB = -0.4D0 * l22AB

       lx3AB = 2.5D0*lx1AB - lx3AB
       lx4AB = 2.5D0*lx2AB - lx4AB
       
       Lx1(i_s,j_s) = lx1AB
       Lx2(i_s,j_s) = lx2AB 
       Lx3(i_s,j_s) = lx3AB
       Lx4(i_s,j_s) = lx4AB
       
    ENDDO
    ENDDO

    ! mux1 : NEOCLASSICAL VISCOSITY COEFFICIENT : \bar{\mu}_{1a}
    ! mux2 : NEOCLASSICAL VISCOSITY COEFFICIENT : \bar{\mu}_{2a} 
    ! mux3 : NEOCLASSICAL VISCOSITY COEFFICIENT : \bar{\mu}_{3a} 
    ! mux3 : NEOCLASSICAL VISCOSITY COEFFICIENT : \bar{\mu}_{3a} 
    
    DO i_s = 1, NSMAX
       
       mu1A = Mu1(i_s)
       mu2A = Mu2(i_s)
       mu3A = Mu3(i_s)
       
       mux1A = mu1A - mu2A
       mux2A = 0.4D0* mu2A 
       mux3A = mu2A - mu3A
       mux4A = 0.4D0* mu3A
       
       mux3A = 2.5D0*mux1A + mux3A
       mux4A = 2.5D0*mux2A + mux4A
       
       Mux1(i_s) = mux1A
       Mux2(i_s) = mux2A
       Mux3(i_s) = mux3A
       Mux4(i_s) = mux4A

    ENDDO
    
    RETURN
  
  END SUBROUTINE T2CALV_MODIFY_COEFFICIENT
  
END MODULE T2CALV
