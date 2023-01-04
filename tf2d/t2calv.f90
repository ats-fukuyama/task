!--------------------------------------------------------------------
!
! MODULE T2CALV
! 
!
!      CALCULATION OF PLASMA PARAMETERS  
!       FOR MULTI-FLUID TRANSPORT MODEL 
!      (ANOMALOUS TRANSPORT MODEL IS BASED ON TWO-FLUID MODEL)
!
!      VALIDITY OF VISCOUS AND FRICTION COEFFICIENTS 
!      ARE PARTLY CONFIRMED BY THE COMPARISON WITH 
!      TASK/T2 RESULTS AND ANALITICAL VALUES IN CTMP.
!
!                   LAST UPDATE    2014-07-01 H.Seto
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
  !                     LAST UPDATE 2014-07-01 H.Seto
  !
  !
  !---------------------------------------------------------
  
  SUBROUTINE T2CALV_EXECUTE(i_m)
    
    INTEGER,INTENT(IN)::i_m

    CALL T2CALV_SETUP(i_m)
    
    CALL T2CALV_FRICTION_COEFFICIENT  ! L11, L12, L21, L22
    
    CALL T2CALV_VISCOUS_COEFFICIENT   ! Mu1, Mu2, Mu3

    CALL T2CALV_ANOMALOUS_COEFFICIENT ! F^{QL}_{a\zeta}, G^{QL}_{a\zeta}
    
    CALL T2CALV_ADDITIONAL_COEFFICIENT
    
    RETURN

  END SUBROUTINE T2CALV_EXECUTE
  
  SUBROUTINE T2CALV_SETUP(i_m)
    
    USE T2CNST,ONLY: Amp,Aee,Pi,Eps0
    USE T2VOUT,ONLY: T2VOUT_EXECUTE    
    USE T2GOUT,ONLY: T2_GOUT
    USE T2COMM,ONLY:&
         & NSMAX,EtNF,BpNF,BtNF,NFMAX,&
         & NnNF,FrNF,FbNF,FtNF,FpNF,&
         & PpNF,QrNF,QbNF,QtNF,QpNF,&
         & i2crt,d2rzm,GlobalCrd,Metric,d2xout,&
         !
         & XvecIn,Xvec,&
         & Pa, Pz, R_rz,  R_mc,&
         & GRt  , G11xCo, G12Co, G22Co,  G33Co,&
         &        G11Ct,  G12Ct, G22xCt, G33Ct,&  
         & UgrCt, UgpCt,&
         & BtCo, BtCt, BtSq, BpCo, BpCt, BpSq, Bb, BbSq,EtCo,&
         & Mm, Ee,Vv,BaseNu,Nu,Nuh,X,Y,Z,&
         & Nn, FrCt, Fb, FtCo, FtCt, FpCo, FpCt, UrCt, Ub, UtCo, UpCt, UuSq,&
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
         & nnE,ttE_eV,dPsiDr,lnLambAB,temp
    
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
    !
    ! initialization
    !
    i_m1d = i2crt(    3,i_m)
    i_m2d = i2crt(    2,i_m)
    R_rz  = d2rzm(    1,i_m)
    R_mc  = GlobalCrd(1,i_m)
        
    GRt      = Metric(1,i_m)
    G11xCo   = Metric(2,i_m)
    G12Co    = Metric(3,i_m)
    G22Co    = Metric(4,i_m)
    G33Co    = Metric(5,i_m)
    G11Ct    = Metric(6,i_m)
    G12Ct    = Metric(7,i_m)
    G22xCt   = Metric(8,i_m)
    G33Ct    = Metric(9,i_m)
    
    UgrCt  = 0.D0
    UgpCt  = 0.D0

    DO i_s = 1,NSMAX
       Mm(i_s) = Pa(i_s)*Amp ! mass 
       Ee(i_s) = Pz(i_s)*Aee ! electric charge
    ENDDO
    
    dpsidr = XvecIn(1,i_m1d)*BpNF
    BtCo   = XvecIn(2,i_m1d)*BtNF
    EtCo   = XvecIn(3,i_m1d)*EtNF

    BtCt   = BtCo*G33Ct
    
    BpCt   = dpsidr/GRt
    BpCo   = BpCt*G22Co
    
    BpSq   = BpCo*BpCt
    BtSq   = BtCo*BtCt
    
    BbSq   = BpSq + BtSq
    Bb     = SQRT(BbSq)

    DO i_s = 1, NSMAX
       
       vOffset =  10*(i_s - 1) + NFMAX
       nnA   = XvecIn(vOffset+ 1,i_m2d)*NnNF
       frCtA = XvecIn(vOffset+ 2,i_m2d)*FrNF
       fbA   = XvecIn(vOffset+ 3,i_m2d)*FbNF
       ftCoA = XvecIn(vOffset+ 4,i_m2d)*FtNF
       fpCtA = XvecIn(vOffset+ 5,i_m2d)*FpNF
       ppA   = XvecIn(vOffset+ 6,i_m2d)*PpNF
       qrCtA = XvecIn(vOffset+ 7,i_m2d)*QrNF
       qbA   = XvecIn(vOffset+ 8,i_m2d)*QbNF
       qtCoA = XvecIn(vOffset+ 9,i_m2d)*QtNF
       qpCtA = XvecIn(vOffset+10,i_m2d)*QpNF
       
       IF((nnA.LT.0.D0).OR.(ppA.LT.0.D0))THEN
          WRITE(6,*)'NEGATIVE  DENSITY or PRESSURE'
          WRITE(6,*)'SPECIS=',i_s,'NODE=',i_m2d,&
               'N=',nnA*1.D-20,'1.D20/m3',&
               'P=',ppA*1.D-23/Aee,'keV*1.D20/m3'
          d2xout = XvecIn
          !CALL T2VOUT_EXECUTE
          CALL T2_GOUT
          STOP       
       ENDIF
       
       urCtA = frCtA/nnA
       ubA   = fbA  /nnA
       utCoA = ftCoA/nnA
       upCtA = fpCtA/nnA
              
       wrCtA = qrCtA/ppA
       wbA   = qbA  /ppA
       wtCoA = qtCoA/ppA
       wpCtA = qpCtA/ppA
       
       Nn(  i_s) = nnA
       FrCt(i_s) = frCtA
       Fb(  i_s) = fbA
       FtCo(i_s) = ftCoA
       FtCt(i_s) = ftCoA*G33Ct
       FpCt(i_s) = fpCtA
       FpCo(i_s) = frCtA*G12Co + fpCtA*G22Co

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
    ! ref: Tokamaks 3rd Ed. J.Wesson et al., P.740
    !nnE_E20 = Nn(1)*1.D-20
    !ttE_keV = Tt(1)*1.D-3/Aee
    ! electron-electron collision
    !lnLamb(1,1)        = 14.9D0 - 0.5D0*LOG(nnE_E20) + LOG(ttE_keV)
    !! electron-ion collision
    !DO i_s = 1, NSMAX!
    !   lnLamb(i_s,  1) = 15.2D0 - 0.5D0*LOG(nnE_E20) + LOG(ttE_keV)
    !   lnLamb(  1,i_s) = lnLamb(  1,i_s)
    !ENDDO

    DO j_s = 1, NSMAX
    DO i_s = 1, NSMAX
       lnLamb(i_s,j_s) = 17.D0! for debug
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
               &          / (4.D0*Pi*(Eps0**2)*(mmA**2)*(vvA**3))
       ENDDO
    ENDDO

    ! COLLISION FREQUENCY [1/s]
    ! REF: COLLISIONAL TRANSPORT IN MAGNETIZED PLASMAS
    !      P. HELANDER AND D.J. SIGMAR (2002)  P.277

    DO j_s = 1, NSMAX
    DO i_s = 1, NSMAX
       Nu(i_s,j_s) = BaseNu(i_s,j_s)/(0.75D0*SQRT(Pi))
    ENDDO
    ENDDO

    ! HEAT EXCHANGE FREQUENCY [Hz]

    temp = 3.D0*SQRT(2.D0)*(SQRT(PI)**3)*(EPS0**2)
    DO j_s = 1,NSMAX
    DO i_s = 1,NSMAX
       Nuh(i_s,j_s) = temp*Mm(i_s)*Mm(j_s)&
            &  * (SQRT(Tt(i_s)/Mm(i_s) + Tt(j_s)/Mm(j_s))**3) &
            &  /(Nn(j_s)*(Ee(i_s)**2)*(Ee(j_s)**2)*lnLamb(i_s,j_s))
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

  !
  ! RESULT OF TEST CALCULATION
  !
  !  TASK/T2 with lnLamb=17.D0   H and S .P242 with Z=1.D0
  !
  ! l^{ee}_{11}/Coef = 0.9986       1.00
  ! l^{ee}_{12}/Coef = 1.4994       1.50
  ! l^{ee}_{21}/Coef = 1.4994       1.50
  ! l^{ee}_{22}/Coef = 4.6631       4.66
  !
  !  TASK/T2 with lnLamb=17.D0   H and S .P242 with Z=1.D0
  !
  ! l^{ei}_{11}/Coef = 0.9986       1.00
  ! l^{ei}_{12}/Coef = 4.0830D-4    0
  ! l^{ei}_{21}/Coef = 1.4994       1.50
  ! l^{ei}_{22}/Coef = 1.8368D-3    0
  !
  SUBROUTINE T2CALV_FRICTION_COEFFICIENT
    
    USE T2CNST,ONLY: Eps0,Pi,Aee
    
    USE T2COMM,ONLY:&
         & NSMAX, R_rz, R_mc, BaseNu,Nu,X, Y, Z, Mm, Nn,&
         & L11,L12,L21,L22,Hex
    
    INTEGER(ikind)::i_s,j_s,k_s
    REAL(   rkind)::&
         & temp,xAB,yAB,zAB,zBA,&
         & nuAB,m00AB,m01AB,m10AB,m11AB, &
         & nuAC,m00AC,m01AC,m10AC,m11AC, &
         &      n00AB,n01AB,n10AB,n11AB, &
         &      l11AB,l12AB,l21AB,l22AB
    
    REAL(   rkind),DIMENSION(1:NSMAX,1:NSMAX)::&
         & m00,m01,m10,m11,n00,n01,n10,n11,nuh
    
    !REAL(   rkind):: mmE,nnE,nuEI,temp2! for debug
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
       nuAB  = Nu( i_s,j_s)

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
             nuAC  = Nu( i_s,k_s)
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
    
    !>>>>> ********** for debug ************* >>>>>
    !mmE   = Mm(1)
    !nnE   = Nn(1)
    !nuEI  = Nu(1,2)
    !temp2 = -mmE*nnE*nuEI
    !print*,'r=',SQRT(R_mc)
    !print*,'L^{ee}_{11}=',L11(1,1)/temp2,'L^{ee}_{12}/C =',L12(1,1)/temp2
    !print*,'L^{ee}_{21}=',L21(1,1)/temp2,'L^{ee}_{12}/C =',L22(1,1)/temp2
    !temp2 = mmE*nnE*nuEI
    !print*,'r=',SQRT(R_mc)
    !print*,'L^{ei}_{11}=',L11(1,2)/temp2,'L^{ei}_{12}/C =',L12(1,2)/temp2
    !print*,'L^{ei}_{21}=',L21(1,2)/temp2,'L^{ei}_{12}/C =',L22(1,2)/temp2
    !stop
    ! <<<<< ********** for debug ************* <<<<<
   
    ! HEAT EXCHANGE TERM
    DO i_s = 1, NSMAX
       Hex(i_s) = 0.D0
       DO j_s = 1, NSMAX
          Hex(i_s) = Hex(i_s) + 1.5D0*(1.D0 - Z(j_s,i_s))*nuh(i_s,j_s)
       ENDDO
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2CALV_FRICTION_COEFFICIENT

  ! NEOCLASSICAL PARALLEL COEFFICIENTS
  ! Ref: COLLISIONAL TRANSPORT IN MAGNETIZED PLASMAS
  !      P. HELANDER AND D.J. SIGMAR (2002)  CHAP.12 AND P.279
  !
  !                 LAST UPDATE 2014-06-09 H.Seto 
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
  !
  !  *****  RESULT OF TEST CALCULATION *****
  ! in Banana Regime 
  !   TASK/T2 with lnLamb= 17.D0    H & S P241 (12.48)
  !  \hat{mu}_i1 = Coef*( 0.554)    Coef*( 0.533) 
  !  \hat{mu}_i2 = Coef*(-0.634)    Coef*(-0.625)
  !  \hat{mu}_i3 = Coef*( 1.413)    Coef*( 1.387)
  !
  ! in P-S Regime 
  !
  !   TASK/T2 with lnLamb= 17.D0    H & S P241 (12.48)
  !   K_i11 = Coef*( 1.19)          Coef*( 1.26) 
  !   K_i12 = Coef*( 5.56)          Coef*( 5.99)
  !   K_i22 = Coef*( 31.8)          Coef*( 39.4)
  !
  SUBROUTINE T2CALV_VISCOUS_COEFFICIENT

    USE T2CNST,ONLY: DeCoef
    USE T2COMM,ONLY: NSMAX,RA,RR,R_mc,R_rz,I_xa,&
         &           Mm,Nn,Pp,BtCt,BpCt,Mu1,Mu2,Mu3,Nu
    USE t2lib, ONLY:integ_f
    
    
    REAL(   rkind)::&
         & iar,iarRt,f_t,f_c,q,temp,&
         & mmA,nnA,ppA,&
         & psCoef,k11ps,k12ps,k22ps,&
         & bnCoef,k11bn,k12bn,k22bn,&
         &        k11,  k12,  k22,derr
    REAL(   rkind)::&
         temp2,nuII,mu1I,mu2I,mu3I
    
    IF((R_mc.GT.0.D0).AND.(R_mc.LT.1.D0))THEN
       
       iar   = RA*SQRT(R_mc)/RR
       iarRt = SQRT(iar)
       f_t   = 1.46D0*iarRt - 0.46D0*(iarRt**3)
       f_c   = 1.D0 - f_t
       q     = BtCt/BpCt 
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
          
          !>>>>> ********** for debug ************* >>>>>
          !IF(i_xa.eq.2)THEN
          !   temp2 = Pp(2)/Nu(2,2)
          !   print*,'r=',SQRT(R_mc),'K11PS=',k11ps/temp2,&
          !        'K12PS=',k12ps/temp2,'K22PS=',k22ps/temp2
          !ENDIF
          ! <<<<< ********** for debug ************* <<<<<

          ! MODIFIED VISCOSITY COEFFICIENT IN BN REGIME
          CALL INTEG_F(FUNC_k11bn,k11bn,derr,EPS=1.D-8,ILST=0)
          CALL INTEG_F(FUNC_k12bn,k12bn,derr,EPS=1.D-8,ILST=0)
          CALL INTEG_F(FUNC_k22bn,k22bn,derr,EPS=1.D-8,ILST=0)
          !>>>>> ********** for debug ************* >>>>>
          !IF(i_xa.eq.2)THEN
          !   nuII = Nu(2,2)
          !   mu1i = DeCoef/nuII* k11bn
          !   mu2i = DeCoef/nuII*(k12bn-2.50D0*k11bn)
          !   mu3i = DeCoef/nuII*(k22bn-5.00D0*k12bn + 6.25D0*k11bn)
          !   print*,'r=',SQRT(R_mc),'mu1i=',mu1i,'mu2i=',mu2i,'mu3i=',mu3i
          !ENDIF
          ! <<<<< ********** for debug ************* <<<<<
          
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
 
    USE t2lib, ONLY:func_phi,func_G
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
       xB  = xA*xBA                 ! (V/Vt_a) * (Vt_a/Vt_b)
       nuD =      baseNuAB*(func_phi(xB)-func_G(xB))/(xA**3)
       nuS = 2.D0*baseNuAB*zAB*(1.D0+yBA)*func_G(xB)/xA 
       nuB = 2.D0*baseNuAB*func_G(xB)/(Xa**3)
       
       nuT = nuT + 2.D0*nuS - nuB + nuD
    ENDDO
    
    IF(    nuT.GT.0.D0)THEN 
       nuT = 1.D0/nuT
    ELSE
       WRITE(6,*)'ERROR IN FUNC_K11PS'
       WRITE(6,*)nut
       STOP
    END IF
    
    FUNC_k11ps = EXP(-xA**2)*(xA**5)*nuT
    
    RETURN
    
  END FUNCTION FUNC_k11ps
  
  ! CALCULATE K12_PS 
  !        BY DOUBLE EXPONENTIAL INTEGRATION FORMULA
  FUNCTION FUNC_k12ps(xaSq)
 
    USE t2lib, ONLY:func_phi,func_G
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
       xB  = xA*xBA                 ! (V/Vt_a) * (Vt_a/Vt_b)
       nuD =      baseNuAB*(func_phi(xB)-func_G(xB))/(xA**3)
       nuS = 2.D0*baseNuAB*zAB*(1.D0+yBA)*func_G(xB)/xA 
       nuB = 2.D0*baseNuAB*func_G(xB)/(Xa**3)
       nuT = nuT + 2.D0*nuS - nuB + nuD
    ENDDO
    
    IF(    nuT.GT.0.D0)THEN 
       nuT = 1.D0/nuT
    ELSE
       WRITE(6,*)'ERROR IN FUNC_K11PS'
       STOP
    END IF
    
    FUNC_k12ps = EXP(-xA**2)*(xA**7)*nuT
    
    RETURN
    
  END FUNCTION FUNC_k12ps
  
  ! CALCULATE K22_PS 
  !        BY DOUBLE EXPONENTIAL INTEGRATION FORMULA
  FUNCTION FUNC_k22ps(xaSq)
    
    USE t2lib, ONLY:func_phi,func_G
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
    ELSE
       WRITE(6,*)'ERROR IN FUNC_K11PS'
       STOP
    END IF
    
    FUNC_k22ps = EXP(-xA**2)*(xA**9)*nuT
    
    RETURN
    
  END FUNCTION FUNC_k22ps

  ! CALCULATE K11_BN 
  !        BY DOUBLE EXPONENTIAL INTEGRATION FORMULA
  FUNCTION FUNC_k11bn(xaSq)
    
    USE t2lib, ONLY:func_phi,func_G
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
        
    USE t2lib, ONLY:func_phi,func_G
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
    
    USE t2lib, ONLY:func_phi,func_G
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
         & NSMAX,R_mc,Mm,Ee,Nn,Tt,BpCt,BtCo,BbSq,Bb,GRt,G11Ct, &
         & FtAnom1,FtAnom2,FtAnom3,FtAnom4,FtAnom5,  &
         & GtAnom1,GtAnom2,GtAnom3,GtAnom4,GtAnom5

    INTEGER(ikind)::&
         i_s
    REAL(   rkind)::&
         eeE,ttE,eeA,ttA,&
         d_anom,m_anom,x_anom,temp1,temp2
    
    !
    ! COEFFICIENTS OF ANORMALOUS TRANSPORT BY QUASI-LINEAR THEORY
    ! REF: S.I. Itoh, Phys. Fluids B, 4 796         (1992)
    !      K.C. Shaing, Phys. Plasmas, 31 2249      (1988)
    !      M. Honda et. al, Nucl. Fusion, 50 095012 (2010)
    !
    IF(R_mc.LE.1.D0)THEN
       !d_anom = 0.45D0*R_mc+0.05D0
       d_anom = 0.03D0
       m_anom = 0.45D1*R_mc+0.05D1
       x_anom = 0.45D1*R_mc+0.05D1
    ELSE
       d_anom = 0.5D0
       m_anom = 0.5D1
       x_anom = 0.5D1 
    ENDIF
    
    ! Toroidal Force by Anomalous Transport (1)
    ttE = Tt(1)
    eeE = Ee(1)
    DO i_s = 1,NSMAX
       eeA = Ee(i_s)
       temp1 = (eeA/ABS(eeA))*((eeE**2)*d_anom/ttE)
       temp2 = 0.5D0*G11Ct*GRt*BpCt/eeE
       !temp2 = (1.5D0 - (m_anom/d_anom))*G11Ct*GRt*BpCt/eeE
       FtAnom1(i_s) = -temp1*temp2*ttE
       FtAnom2(i_s) =  temp1*temp2
       !FtAnom1(i_s) = 0.D0
       !FtAnom2(i_s) = 0.D0
       FtAnom3(i_s) =  temp1*BtCo*Bb
       FtAnom4(i_s) = -temp1*BbSq
       
    ENDDO
    ! Toroidal Force by Anomalous Transport (2)
    DO i_s = 1, NSMAX
       FtAnom5(i_s) = Mm(i_s)*Nn(i_s)*m_anom*G11Ct
    ENDDO
    
    ! EW Toroidal Force by Anomalous Transport (1)
    DO i_s =1, NSMAX
       ttA   = Tt(i_s)
       eeA   = Ee(i_s)
       temp1 =  (eeA**2)*m_anom/ttA
       temp2 = -(2.5D0-(x_anom/m_anom))*GRt*BpCt*G11Ct*ttA/eeA
       GtAnom1(i_s) =  temp1*temp2*ttA
       GtAnom2(i_s) = -temp1*temp2
       GtAnom3(i_s) = -temp1*BtCo*Bb*0.4D0
       GtAnom4(i_s) =  temp1*BbSq   *0.4D0
    ENDDO
    
    ! Toroidal Force by Anomalous Transport (2)
    DO i_s = 1, NSMAX
       GtAnom5(i_s) = 0.D0
    ENDDO
    
    RETURN
    
  END SUBROUTINE T2CALV_ANOMALOUS_COEFFICIENT
  
  SUBROUTINE T2CALV_ADDITIONAL_COEFFICIENT
    
    USE T2COMM,ONLY:&
         & NSMAX,&
         & GRt,BpCt,BtCo,BtCt,BtSq,Bb,BbSq,&
         & Nn,Pp,Tt,Ee,UrCt,Ub,UtCo,UpCt,WtCo,WpCt,&
         & Mu1,Mu2,Mu3,L11,L12,L21,L22,&
         & BNCXb1,BNCXb2,BNCXb3,BNCXb4,BNCXb5,BNCXb6,&
         & BNCXt1,BNCXt2,BNCXt3,&
         & BNCPp1,BNCPp2,BNCPp3,BNCPp4,BNCPp5,BNCPp6,&
         & BNCPp7,BNCPp8,BNCPp9,&
         & BNCQb1,BNCQb2,BNCQb3,BNCQb4,&
         & BNCQt1,BNCQt2,BNCQt3,BNCQt4,&
         & CNCV01,CNCV02,CNCV03,CNCV04,CNCV05,CNCV06,CNCV07,&
         & CNCV08,CNCV09,CNCV10,CNCV11,CNCV12,CNCV13,&
         & CNCF01,CNCF02,CNCF03,CNCF04

    INTEGER(ikind)::i_s,j_s    
    REAL(   rkind)::b01,b02,b03,b04,b05,b06,b07,b08,b09,b10,b11,b12
    REAL(   rkind)::&
         & usrCt(1:NSMAX),uspCt(1:NSMAX),udpCt(1:NSMAX),wdpCt(1:NSMAX),&
         & mux1( 1:NSMAX),mux2( 1:NSMAX),mux3( 1:NSMAX),mux4( 1:NSMAX),&
         & lx1(  1:NSMAX,1:NSMAX),lx2(  1:NSMAX,1:NSMAX),&
         & lx3(  1:NSMAX,1:NSMAX),lx4(  1:NSMAX,1:NSMAX)

    b01 = 3.D0*GRt
    b02 = BtSq/(Bb**3)
    b03 = 2.D0*BpCt/3.D0
    b04 = BpCt/BtCo
    b05 = BpCt/Bb 
    b06 = b05*b05
    b07 = BpCt*BtCo/BbSq
    b08 = b02*b02
    b09 = 2.D0*BtCt/3.D0
    b10 = 2.D0*BpCt/3.D0
    b11 = BtSq/BbSq-1.D0/3.D0
    b12 = BpCt*BtCo/BbSq

    !        b01*b02*b03*b04*b05*b06*b07*b08*b09*b10*b11*b12
    BNCXb1 = b01*b02*b03*b04
    BNCXb2 = b01*b02*b03
    BNCXb3 = b01    *b03    *b05
    BNCXb4 = b01                *b06
    BNCXb5 = b01*b02    *b04*b05
    BNCXb6 = b01*b02        *b05
    !        b01*b02*b03*b04*b05*b06*b07*b08*b09*b10*b11*b12
    BNCXt1 = b01*b02    *b04        *b07
    BNCXt2 = b01*b02                *b07
    BNCXt3 = b01            *b05    *b07
    !        b01*b02*b03*b04*b05*b06*b07*b08*b09*b10*b11*b12
    BNCPp1 = b01*b02    *b04
    BNCPp2 = b01*b02
    BNCPp3 = b01            *b05
    BNCPp4 = b01*b02        *b05
    BNCPp5 = b01        *b04            *b08
    BNCPp6 = b01                        *b08
    BNCPp7 = b01                *b06        
    BNCPp8 = b01*b02    *b04*b05                    
    BNCPp9 = b01*b02        *b05              
    !        b01*b02*b03*b04*b05*b06*b07*b08*b09*b10*b11*b12
    BNCQb1 = b01*b02                        *b09
    BNCQb2 = b01            *b05            *b09
    BNCQb3 = b01*b02                            *b10
    BNCQb4 = b01            *b05                *b10
    !        b01*b02*b03*b04*b05*b06*b07*b08*b09*b10*b11*b12
    BNCQt1 = b01*b02                                *b11
    BNCQt2 = b01            *b05                    *b11
    BNCQt3 = b01*b02                                    *b12
    BNCQt4 = b01            *b05                        *b12
    
    DO i_s = 1, NSMAX       
       usrCt(i_s) = UrCt(i_s)/3.D0 
       uspCt(i_s) = UpCt(i_s)/3.D0 - Ub(  i_s)*b05
       udpCt(i_s) = UtCo(i_s)*b04  - UpCt(i_s)
       wdpCt(i_s) = WtCo(i_s)*b04  - WpCt(i_s) 
    ENDDO
    
    ! mux1 : NEOCLASSICAL VISCOSITY COEFFICIENT : \bar{\mu}_{1a}
    ! mux2 : NEOCLASSICAL VISCOSITY COEFFICIENT : \bar{\mu}_{2a} 
    ! mux3 : NEOCLASSICAL VISCOSITY COEFFICIENT : \bar{\mu}_{3a} 
    ! mux3 : NEOCLASSICAL VISCOSITY COEFFICIENT : \bar{\mu}_{3a} 
    DO i_s = 1, NSMAX       
       mux1(i_s) = Mu1(i_s) - Mu2(i_s)
       mux2(i_s) =    0.4D0 * Mu2(i_s)
       mux3(i_s) = Mu2(i_s) - Mu3(i_s) + 2.5D0*(Mu1(i_s) - Mu2(i_s))
       mux4(i_s) =    0.4D0 * Mu3(i_s) + Mu2(i_s)
    ENDDO
    
    DO i_s = 1,NSMAX
       CNCV01(i_s) = mux1(i_s)        /Nn(i_s) 
       CNCV02(i_s) = mux2(i_s)        /Pp(i_s)
       CNCV03(i_s) = mux3(i_s)*Tt(i_s)/Nn(i_s)
       CNCV04(i_s) = mux4(i_s)*Tt(i_s)/Pp(i_s)
       ! for viscous heating
       CNCV05(i_s) = mux1(i_s)*usrCt(i_s)/Nn(i_s) 
       CNCV06(i_s) = mux1(i_s)*uspCt(i_s)/Nn(i_s) 
       CNCV07(i_s) = mux2(i_s)*usrCt(i_s)/Pp(i_s)
       CNCV08(i_s) = mux2(i_s)*uspCt(i_s)/Pp(i_s)
       CNCV09(i_s) = mux1(i_s)*udpCt(i_s)/Nn(i_s) 
       CNCV10(i_s) = mux2(i_s)*udpCt(i_s)/Pp(i_s) 
       ! for lorentz
       CNCV11(i_s) =  Ee(i_s)*mux1(i_s)*udpCt(i_s) &
            &        +Ee(i_s)*mux2(i_s)*wdpCt(i_s)
       CNCV12(i_s) =  Ee(i_s)*mux1(i_s)
       CNCV13(i_s) =  Ee(i_s)*mux2(i_s)
    ENDDO
    
    DO j_s = 1, NSMAX
    DO i_s = 1, NSMAX
       lx1(i_s,j_s) =         L11(i_s,j_s) + L12(i_s,j_s) 
       lx2(i_s,j_s) =               -0.4D0 * L12(i_s,j_s) 
       lx3(i_s,j_s) =  2.5D0*(L11(i_s,j_s) + L12(i_s,j_s))&
            &               -(L21(i_s,j_s) + L22(i_s,j_s))
       lx4(i_s,j_s) = -L12(i_s,j_s) +0.4D0 * L22(i_s,j_s)
    ENDDO
    ENDDO
    
    DO j_s=1,NSMAX
    DO i_s=1,NSMAX
       CNCF01(i_s,j_s) = GRt*lx1(i_s,j_s)        /Nn(j_s)
       CNCF02(i_s,j_s) = GRt*lx2(i_s,j_s)        /Pp(j_s)
       CNCF03(i_s,j_s) = GRt*lx3(i_s,j_s)*Tt(i_s)/Nn(j_s)
       CNCF04(i_s,j_s) = GRt*lx4(i_s,j_s)*Tt(i_s)/Pp(j_s)
    ENDDO
    ENDDO

    RETURN
    
  END SUBROUTINE T2CALV_ADDITIONAL_COEFFICIENT

  !
  !
  !
  !
  SUBROUTINE T2CALV_COEF_CHECK

    USE T2COMM,ONLY:&
         & NSMAX,&
         & GRt,G11xCo,G12Co,G22Co, G33Co,&
         &     G11Ct, G12Ct,G22xCt,G33Ct,&  
         &     UgrCt, UgpCt,R_rz,R_mc,   &
         !
         & BtCo,BtCt,BtSq,BpCo,BpCt,BpSq,Bb,BbSq,&
         !
         & Ee,Mm,Nn,Vv,FrCt,Fb,FtCo,FpCt,UrCt,Ub,&
         & UtCo,UpCt,UuSq,Pp,QrCt,Qb,QtCo,QpCt,WrCt,&
         & Wb,WtCo,WpCt,Tt,&
         !
         & X,Y,Z,BaseNu,Nu,Nuh,&
         !
         & L11,L12,L21,L22,Mu1,Mu2,Mu3,Hex,&   
         & FtAnom1,FtAnom2,FtAnom3,FtAnom4,FtAnom5,&
         & GtAnom1,GtAnom2,GtAnom3,GtAnom4,GtAnom5,&
         & BNCXb1,BNCXb2,BNCXb3,BNCXb4,BNCXb5,BNCXb6,&
         & BNCXt1,BNCXt2,BNCXt3,&
         & BNCPp1,BNCPp2,BNCPp3,BNCPp4,BNCPp5,BNCPp6,&
         & BNCPp7,BNCPp8,BNCPp9,&
         & BNCQb1,BNCQb2,BNCQb3,BNCQb4,&
         & BNCQt1,BNCQt2,BNCQt3,BNCQt4,& 
         !
         & CNCV01,CNCV02,CNCV03,CNCV04,CNCV05,&
         & CNCV06,CNCV07,CNCV08,CNCV09,CNCV10,&
         & CNCV11,CNCV12,CNCV13,&
         & CNCF01,CNCF02,CNCF03,CNCF04


    INTEGER(ikind)::i_s,j_s

    OPEN(30,FILE="TEST_T2CALV.txt")
100 FORMAT(A10,1X,'val=',1X,D15.8)
110 FORMAT(A10,1X,'i=',I1,1X,'val=',D15.8)
120 FORMAT(A10,1X,'i=',I1,1X,'j=',I1,1X,'val=',D15.8)
    
    WRITE(30,100)"GRt    ",GRt
    WRITE(30,100)"G11xCo ",G11xCo
    WRITE(30,100)"G12Co  ",G12Co
    WRITE(30,100)"G22Co  ",G22Co
    WRITE(30,100)"G33Co  ",G33Co
    WRITE(30,100)"G11Ct  ",G11Ct
    WRITE(30,100)"G12Ct  ",G12Ct
    WRITE(30,100)"G22xCt ",G22xCt
    WRITE(30,100)"G33Ct  ",G33Ct
    WRITE(30,100)"UgrCt  ",UgrCt
    WRITE(30,100)"UgpCt  ",UgpCt
    WRITE(30,100)"R_rz   ",R_rz
    WRITE(30,100)"R_mc   ",R_mc
    
    WRITE(30,100)"BtCo   ",BtCo
    WRITE(30,100)"BtCt   ",BtCt
    WRITE(30,100)"BtSq   ",BtSq
    WRITE(30,100)"BpCo   ",BpCo
    WRITE(30,100)"BpCt   ",BpCt
    WRITE(30,100)"BpSq   ",BpSq
    WRITE(30,100)"Bb     ",Bb
    WRITE(30,100)"BbSq   ",BbSq
    
    DO i_s = 1, NSMAX
       WRITE(30,110)"Ee     ",i_s,Ee(  i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"Mm     ",i_s,Mm(  i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"Tt     ",i_s,Tt(  i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"Vv     ",i_s,Vv(  i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"Nn     ",i_s,Nn(  i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"FrCt   ",i_s,FrCt(i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"Fb     ",i_s,Fb(  i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"FtCo   ",i_s,FtCo(i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"FpCt   ",i_s,FpCt(i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"UrCt   ",i_s,UrCt(i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"Ub     ",i_s,Ub(  i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"UtCo   ",i_s,UtCo(i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"UpCt   ",i_s,UpCt(i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"UuSq   ",i_s,UuSq(i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"Pp     ",i_s,Pp(  i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"QrCt   ",i_s,QrCt(i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"Qb     ",i_s,Qb(  i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"QtCo   ",i_s,QtCo(i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"QtCo   ",i_s,QtCo(i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"QpCt   ",i_s,QpCt(i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"WrCt   ",i_s,WrCt(i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"Wb     ",i_s,Wb(  i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"WtCo   ",i_s,WtCo(i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"WpCt   ",i_s,WpCt(i_s)
    ENDDO

    DO i_s = 1, NSMAX
       WRITE(30,110)"Mu1    ",i_s,Mu1(  i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"Mu2    ",i_s,Mu2(  i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"Mu3    ",i_s,Mu3(  i_s)
    ENDDO

    DO i_s = 1, NSMAX
    DO j_s = 1, NSMAX
       WRITE(30,120)"X      ",i_s,j_s,X(     i_s,j_s)
    ENDDO
    ENDDO
    DO i_s = 1, NSMAX
    DO j_s = 1, NSMAX
       WRITE(30,120)"Y      ",i_s,j_s,Y(     i_s,j_s)
    ENDDO
    ENDDO
    DO i_s = 1, NSMAX
    DO j_s = 1, NSMAX
        WRITE(30,120)"Z     ",i_s,j_s,Z(     i_s,j_s)
    ENDDO
    ENDDO
    DO i_s = 1, NSMAX
    DO j_s = 1, NSMAX
       WRITE(30,120)"BaseNu ",i_s,j_s,BaseNu(i_s,j_s)
    ENDDO
    ENDDO
    DO i_s = 1, NSMAX
    DO j_s = 1, NSMAX
       WRITE(30,120)"Nu     ",i_s,j_s,Nu(    i_s,j_s)
    ENDDO
    ENDDO
    DO i_s = 1, NSMAX
    DO j_s = 1, NSMAX
       WRITE(30,120)"Nuh    ",i_s,j_s,Nuh(   i_s,j_s)
    ENDDO
    ENDDO

    DO i_s = 1, NSMAX
    DO j_s = 1, NSMAX
       WRITE(30,120)"L11    ",i_s,j_s,L11(i_s,j_s)
    ENDDO
    ENDDO    
    DO i_s = 1, NSMAX
    DO j_s = 1, NSMAX
       WRITE(30,120)"L12    ",i_s,j_s,L12(i_s,j_s)    
    ENDDO
    ENDDO    
    DO i_s = 1, NSMAX
    DO j_s = 1, NSMAX
       WRITE(30,120)"L21    ",i_s,j_s,L21(i_s,j_s)
    ENDDO
    ENDDO    
    DO i_s = 1, NSMAX
    DO j_s = 1, NSMAX
       WRITE(30,120)"L22    ",i_s,j_s,L22(i_s,j_s)
    ENDDO
    ENDDO    
    DO i_s = 1, NSMAX
       WRITE(30,110)"Hex    ",i_s,Hex(i_s)
    ENDDO    

    DO i_s = 1, NSMAX
       WRITE(30,110)"FtAnom1",i_s,FtAnom1(i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"FtAnom2",i_s,FtAnom2(i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"FtAnom3",i_s,FtAnom3(i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"FtAnom4",i_s,FtAnom4(i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"FtAnom5",i_s,FtAnom5(i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"GtAnom1",i_s,GtAnom1(i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"GtAnom2",i_s,GtAnom2(i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"GtAnom3",i_s,GtAnom3(i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"GtAnom4",i_s,GtAnom4(i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"GtAnom5",i_s,GtAnom5(i_s)
    ENDDO

    WRITE(30,100)"BNCXb1 ",BNCXb1
    WRITE(30,100)"BNCXb2 ",BNCXb2
    WRITE(30,100)"BNCXb3 ",BNCXb3
    WRITE(30,100)"BNCXb4 ",BNCXb4
    WRITE(30,100)"BNCXb5 ",BNCXb5
    WRITE(30,100)"BNCXb6 ",BNCXb6

    WRITE(30,100)"BNCXt1 ",BNCXt1
    WRITE(30,100)"BNCXt2 ",BNCXt2
    WRITE(30,100)"BNCXt3 ",BNCXt3

    WRITE(30,100)"BNCPp1 ",BNCPp1
    WRITE(30,100)"BNCPp2 ",BNCPp2
    WRITE(30,100)"BNCPp3 ",BNCPp3
    WRITE(30,100)"BNCPp4 ",BNCPp4
    WRITE(30,100)"BNCPp5 ",BNCPp5
    WRITE(30,100)"BNCPp6 ",BNCPp6
    WRITE(30,100)"BNCPp7 ",BNCPp7
    WRITE(30,100)"BNCPp8 ",BNCPp8
    WRITE(30,100)"BNCPp9 ",BNCPp9

    WRITE(30,100)"BNCQb1 ",BNCQb1
    WRITE(30,100)"BNCQb2 ",BNCQb2
    WRITE(30,100)"BNCQb3 ",BNCQb3
    WRITE(30,100)"BNCQb4 ",BNCQb4

    WRITE(30,100)"BNCQt1 ",BNCQt1
    WRITE(30,100)"BNCQt2 ",BNCQt2
    WRITE(30,100)"BNCQt3 ",BNCQt3
    WRITE(30,100)"BNCQt4 ",BNCQt4

    DO i_s = 1, NSMAX
       WRITE(30,110)"CNCV01 ",i_s,CNCV01(i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"CNCV02 ",i_s,CNCV02(i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"CNCV03 ",i_s,CNCV03(i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"CNCV04 ",i_s,CNCV04(i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"CNCV05 ",i_s,CNCV05(i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"CNCV06 ",i_s,CNCV06(i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"CNCV07 ",i_s,CNCV07(i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"CNCV08 ",i_s,CNCV08(i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"CNCV09 ",i_s,CNCV09(i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"CNCV10 ",i_s,CNCV10(i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"CNCV11 ",i_s,CNCV11(i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"CNCV12 ",i_s,CNCV12(i_s)
    ENDDO
    DO i_s = 1, NSMAX
       WRITE(30,110)"CNCV13 ",i_s,CNCV13(i_s)
    ENDDO

    DO i_s = 1, NSMAX
    DO j_s = 1, NSMAX
       WRITE(30,120)"CNCF01 ",i_s,j_s,CNCF01(i_s,j_s)
    ENDDO
    ENDDO
    DO i_s = 1, NSMAX
    DO j_s = 1, NSMAX
       WRITE(30,120)"CNCF02 ",i_s,j_s,CNCF02(i_s,j_s)
    ENDDO
    ENDDO
    DO i_s = 1, NSMAX
    DO j_s = 1, NSMAX
       WRITE(30,120)"CNCF03 ",i_s,j_s,CNCF03(i_s,j_s)
    ENDDO
    ENDDO
    DO i_s = 1, NSMAX
    DO j_s = 1, NSMAX
       WRITE(30,120)"CNCF04 ",i_s,j_s,CNCF04(i_s,j_s)
    ENDDO
    ENDDO    

    CLOSE(30)
    
    RETURN

  END SUBROUTINE T2CALV_COEF_CHECK
END MODULE T2CALV
