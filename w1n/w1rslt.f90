MODULE w1rslt

CONTAINS

!     ****** POWER ABSORPTION AS A FUNCTION OF KZ ******

  SUBROUTINE W1CLPW(NZ)
    USE w1comm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NZ
    INTEGER:: NX,NS,NXANT1,NXANT2,I
    REAL(rkind):: FACT

    IF(NSYM.NE.0.AND.NZ.NE.1.AND.NZ.NE.(NZPMAX/2+1)) THEN
       FACT = 2.0*RZ
    ELSEIF(NSYM.EQ.-1.AND.NZ.EQ.(NZPMAX/2+1)) THEN
       FACT = 0.0
    ELSE
       FACT =     RZ
    ENDIF

    DO NS=1,NSMAX
       DO NX=1,NXPMAX
          PABSX(NX,NS)=PABSX(NX,NS)+PABS(NX,NS)*FACT
       END DO
    END DO

    DO NX=1,NXPMAX
       FLUXX(NX)=FLUXX(NX)+FLUX(NX)*FACT
    END DO

    DO NS=1,NSMAX
       PAK(NZ,NS)=0.D0
       DO NX=1,NXPMAX
          PAK(NZ,NS)=PAK(NZ,NS)+PABS(NX,NS)*RZ
       END DO
    END DO

    NXANT1=NXPMAX+NXVMAX+NXVMAX/2+1
    NXANT2=NXPMAX       +NXVMAX/2
    PANTK(NZ)=-RZ*DCONJG(CE2DA(NZ,NXANT1,2))*CFJY1 &
              -RZ*DCONJG(CE2DA(NZ,NXANT2,2))*CFJY2 &
              -RZ*DCONJG(CE2DA(NZ,NXANT1,3))*CFJZ1 &
              -RZ*DCONJG(CE2DA(NZ,NXANT2,3))*CFJZ2

    PAKT(NZ,3)=0.D0
    DO NS=1,NSMAX
       PAKT(NZ,3)=PAKT(NZ,3)+PAK(NZ,NS)
    END DO

    PAKT(NZ,2)=0.
    DO I=1,3
       PAKT(NZ,2)=PAKT(NZ,2)+( FLUXX(NXPMAX+NXVMAX+NXVMAX) &
                              -FLUXX(1          ) &
                              +FLUXX(NXPMAX        ) &
                              -FLUXX(NXPMAX+1      ))
    END DO
    PAKT(NZ,1)=PAKT(NZ,2)+PAKT(NZ,3)

    RETURN
  END SUBROUTINE W1CLPW

!     ****** INITIALIZE POWER ABSORPTION ******

  SUBROUTINE W1PWRI
    USE w1comm
    IMPLICIT NONE
    INTEGER:: NX,NS

    DO NX=1,NXTMAX
       FLUXX(NX)=0.D0
    END DO

    DO NS=1,NSMAX
       DO NX=1,NXPMAX
          PABSX(NX,NS)=0.D0
       END DO
    END DO

    DO NX=1,NXPMAX
       AJCDX(NX)=0.D0
    END DO

    RETURN
  END SUBROUTINE W1PWRI

!     ******* POWER ABSORPTION AND ENERGY FLUXX *******

  SUBROUTINE W1PWRS
    USE w1comm
    IMPLICIT NONE
    INTEGER:: NS,NX,NXANT1,NXANT2,NZ,NA
    REAL(rkind):: PHASE
    COMPLEX(rkind):: CPHASE,CJYH,CJZH,CJYL,CJZL
    
    AJCDT=0.D0
    DO NS=1,NSMAX
       PABSXZ(NS)=0.D0
    END DO

    DO NS=1,NSMAX
       DO NX=1,NXPMAX
          PABSXZ(NS)  =PABSXZ(NS)  +PABSX(NX,NS)
       END DO
    END DO

    DO NS=1,NSMAX
       PABSX(     1,NS)=2.D0*PABSX(     1,NS)/(XA(2)-XA(1))
       PABSX(NXPMAX,NS)=2.D0*PABSX(NXPMAX,NS)/(XA(NXPMAX)-XA(NXPMAX-1))
       DO NX=2,NXPMAX-1
          PABSX( NX,NS)=2.D0*PABSX( NX,NS)/(XA(NX+1)-XA(NX-1))
       END DO
    END DO

    PABSTT=0.D0
    DO NS=1,NSMAX
       PABSTT=PABSTT+PABSXZ(NS)
    END DO

    AJCDT=0.D0
    DO NX=1,NXPMAX
       AJCDT =AJCDT +AJCDX(NX)
    END DO
    AJCDX(  1   )=2.D0*AJCDX(  1   )/(XA(2)-XA(1))
    AJCDX(NXPMAX)=2.D0*AJCDX(NXPMAX)/(XA(NXPMAX)-XA(NXPMAX-1))
    DO NX=2,NXPMAX-1
       AJCDX( NX   )=2.D0*AJCDX( NX   )/(XA(NX+1)-XA(NX-1))
    END DO

    NXANT1=NXPMAX+NXVMAX+NXVMAX/2+1
    NXANT2=NXPMAX       +NXVMAX/2
    PANT1 = 0.D0
    PANT2 = 0.D0
    DO NZ = 1 , NZPMAX
       PANT1 = PANT1 + DCONJG(CE2DA(NZ,NXANT1,2))*CJ1(NZ) &
                     + DCONJG(CE2DA(NZ,NXANT1,3))*CJ3(NZ)
       PANT2 = PANT2 + DCONJG(CE2DA(NZ,NXANT2,2))*CJ2(NZ) &
                     + DCONJG(CE2DA(NZ,NXANT2,3))*CJ4(NZ)
    END DO
    PANT1 = - PANT1 * DZ
    PANT2 = - PANT2 * DZ
    PANT  = PANT1 + PANT2

    PWALL1 = -FLUXX(NXPMAX+NXVMAX+1)
    PWALL2 =  FLUXX(NXPMAX+NXVMAX  )
    PWALL  =  PWALL1 + PWALL2

!    PIBW1 = FLUXX(NXPMAX+NXVMAX+NXVMAX) - FLUXX(       1)
    PIBW1 = 0.D0
!    PIBW2 = FLUXX(NXPMAX              ) - FLUXX(NXPMAX+1)
    PIBW2 = 0.D0
    PIBW  = PIBW1 + PIBW2

    PERR  = PANT - ( PABSTT + PWALL + PIBW )

    DO NA = 1 , NAMAX
       RANT1(NA)=0.D0
       XANT1(NA)=0.D0
       CJYH=CJ1(NZANTYH(NA))
       IF( ABS( CJYH ) .GT. 1.E-6 ) THEN
          RANT1(NA) = RANT1(NA)-DBLE( CE2DA(NZANTYH(NA),NXANT1,2)/ CJYH)
          XANT1(NA) = XANT1(NA)+IMAG( CE2DA(NZANTYH(NA),NXANT1,2)/ CJYH)
       END IF
       DO NZ=NZANTZH(NA),NZANTLH(NA)
          CJZH=CJ3(NZ)
          IF( ABS( CJZH ) .GT. 1.E-6 ) THEN
             RANT1(NA) = RANT1(NA)-DBLE( CE2DA(NZ,NXANT1,3)/ CJZH)
             XANT1(NA) = XANT1(NA)+IMAG( CE2DA(NZ,NXANT1,3)/ CJZH)
          END IF
       END DO
       RANT2(NA)=0.D0
       XANT2(NA)=0.D0
       CJYL=CJ2(NZANTYL(NA))
       IF( ABS( CJYL ) .GT. 1.E-6 ) THEN
          RANT2(NA) = RANT2(NA)-DBLE( CE2DA(NZANTYL(NA),NXANT2,2)/ CJYL)
          XANT2(NA) = XANT2(NA)+IMAG( CE2DA(NZANTYL(NA),NXANT2,2)/ CJYL)
       END IF
       DO NZ=NZANTZL(NA),NZANTLL(NA)
          CJZL=CJ4(NZ)
          IF( ABS( CJZL ) .GT. 1.E-6 ) THEN
             RANT2(NA) = RANT2(NA)-DBLE( CE2DA(NZ,NXANT2,3)/ CJZL)
             XANT2(NA) = XANT2(NA)+IMAG( CE2DA(NZ,NXANT2,3)/ CJZL)
          END IF
       END DO
    END DO
    RETURN
  END SUBROUTINE W1PWRS

!     ****** CALCULATE DRIVEN CURRENT ******

  SUBROUTINE W1CLCD(NZ)
    USE w1comm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NZ
    INTEGER:: NX,NS
    REAL(rkind):: XR,VTE,VPH,W,EFCD
    REAL(rkind):: Y0,Y1,Y2,Y3,E0,E1,E2,E3,RLNLMD,AJCD

    AJCDK(NZ)=0.D0

    DO NX=1,NXPMAX
       XR=XAM(NX)/RR
       DO NS=1,NSMAX
          IF(IELEC(NS).EQ.1) THEN
             VTE=SQRT(PROFTR(NX,NS)*AEE*1.D3/AME)
             VPH=2.D0*PI*RF*1.D6/RKZ
             W=ABS(VPH/VTE)
             IF(ABS(VPH).GT.VC) THEN
                EFCD=0.D0
             ELSE
                IF(WVYSIZ.LE.0.D0) THEN
                   EFCD=W1CDEF(W,ZEFF,XR,0.D0,NCDTYP)
                ELSE
                   Y0=0.00D0*WVYSIZ/RR
                   Y1=0.25D0*WVYSIZ/RR
                   Y2=0.50D0*WVYSIZ/RR
                   Y3=0.75D0*WVYSIZ/RR
                   E0=W1CDEF(W,ZEFF,XR,Y0,NCDTYP)
                   E1=W1CDEF(W,ZEFF,XR,Y1,NCDTYP)
                   E2=W1CDEF(W,ZEFF,XR,Y2,NCDTYP)
                   E3=W1CDEF(W,ZEFF,XR,Y3,NCDTYP)
!     *** WEIGHTING WITH PARABOLIC ELECTRIC FIELD PROFILE ***
                   EFCD=(256.D0*E0+450.D0*E1+288.D0*E2+98.D0*E3)/1092.D0
!     *** WEIGHTING WITH PARABOLIC POWER DENSITY PROFILE ***
!                   EFCD=(16.D0*E0+30.D0*E1+24.D0*E2+14.D0*E3)/84.D0
                ENDIF
             ENDIF
             IF(PROFPN(NX,1).LE.0.D0) THEN
                AJCD=0.D0
             ELSE
                RLNLMD=16.1D0 - 1.15D0*LOG10(PROFPN(NX,1)) &
                              + 2.30D0*LOG10(PROFTR(NX,1))
                AJCD=0.384D0*PROFTR(NX,NS)*EFCD/(PROFPN(NX,1)*RLNLMD) &
                     *PABS(NX,NS)
             END IF
             IF(VPH.LT.0.D0) AJCD=-AJCD
             AJCDX(NX)=AJCDX(NX)+AJCD
             AJCDK(NZ)=AJCDK(NZ)+AJCD
          ENDIF
       END DO
    END DO
    RETURN
  END SUBROUTINE W1CLCD

!     ****** CURRENT DRIVE EFFICIENCY ******

!      WT = V / VT : PHASE VELOCITY NORMALIZED BY THERMAL VELOCITY
!      Z  = ZEFF   : EFFECTIVE Z
!      XR = X / RR : NORMALIZED X
!      YR = Y / RR : NORMALIZED Y
!      ID : 0 : LANDAU DAMPING
!           1 : TTMP

  FUNCTION W1CDEF(WT,Z,XR,YR,ID)
    USE w1comm,ONLY: rkind
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: WT,Z,XR,YR
    INTEGER,INTENT(IN):: ID
    REAL(rkind):: W1CDEF,R,D,C,A,RM,RC,W,EFF0,EFF1,Y2,Y1,EFF2,YT,ARG,EFF3

    R=SQRT(XR*XR+YR*YR)
    IF(ID.EQ.0) THEN
       D=3.D0/Z
       C=3.83D0
       A=0.D0
       RM=1.38D0
       RC=0.389D0
    ELSE
       D=11.91D0/(0.678D0+Z)
       C=4.13D0
       A=12.3D0
       RM=2.48D0
       RC=0.0987D0
    ENDIF
    IF(WT.LE.1.D-20) THEN
       W=1.D-20
    ELSE
       W=WT
    ENDIF
    EFF0=D/W+C/Z**0.707D0+4.D0*W*W/(5.D0+Z)
    EFF1=1.D0-R**0.77D0*SQRT(3.5D0**2+W*W)/(3.5D0*R**0.77D0+W)

    Y2=(R+XR)/(1.D0+R)
    IF(Y2.LT.0.D0) Y2=0.D0
    Y1=SQRT(Y2)
    EFF2=1.D0+A*(Y1/W)**3

    IF(Y2.LE.1.D-20) THEN
       YT=(1.D0-Y2)*WT*WT/1.D-60
    ELSE
       YT=(1.D0-Y2)*WT*WT/Y2
    ENDIF
    IF(YT.GE.0.D0.AND.RC*YT.LT.40.D0) THEN
       ARG=(RC*YT)**RM
       IF(ARG.LE.100.D0) THEN
          EFF3=1.D0-MIN(EXP(-ARG),1.D0)
       ELSE
          EFF3=1.D0
       ENDIF
    ELSE
       EFF3=1.D0
    ENDIF

    W1CDEF=EFF0*EFF1*EFF2*EFF3
    RETURN
  END FUNCTION W1CDEF

!     ****** CALCULATION OF FORCE AND HELICITY ******

  SUBROUTINE W1HELD(MODE)
    USE w1comm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: MODE
    INTEGER:: NX,K,KK,I,II,NS
    REAL(rkind):: RKV,RNPR,RW,DX,FWP,WP2I,TAU,VPH,FACTC,AWEC
    REAL(rkind):: XR,YR,EP,RLNLMD,VTE,TAU0,FRAC,W,TAU1,TAU2,TAU3,TAU4
    REAL(rkind):: TAUE,TAUI,A6,A7
    COMPLEX(rkind):: CAL,CSKXL,CSPXL,CSPZL,CPH,CEX,CEY,CEZ,CPT
    COMPLEX(rkind):: CJX,CDX,CJY,CDY,CJZ,CDZ,CVX,CVZ
    
    RKV=2.D6*PI*RF/VC
    RNPR=RKZ/RKV
    RW=2.D6*PI*RF

    DO NX = 1 , NXPMAX
       DX = RKV * ( XA( NX+1 ) - XA( NX ) )
       K  = (NX-1)*MODE + 2
       KK = (NX-1)*MODE/2

       DO I=1,3
          CEF(NX,I)=0.D0
          CBF(NX,I)=0.D0
          CAF(NX,I)=0.D0
          CED(NX,I)=0.D0
       END DO
       CAF(NX,4)=0.D0

       DO NS=1,NSMAX
          DO I=1,7
             CJF(NX,NS,I)=0.D0
          END DO
          DO I=1,3
             CJD(NX,NS,I)=0.D0
          END DO
       END DO

       DO I=1,MODE
          CAL=CA(K+I)
          II=(I-1)/2+1
          CSKXL=CSKX(KK+II)
          CSPXL=CSPX(KK+II)
          CSPZL=CSPZ(KK+II)
          IF(MOD(I,2).EQ.0) THEN
             CSKXL=-CSKXL
             CSPZL=-CSPZL
          ENDIF

          CPH  = CDEXP(0.5D0*CI*DX*CSKXL)
          CEX  = CAL*CPH*CSPXL
          CEY  = CAL*CPH
          CEZ  = CAL*CPH*CSPZL
          CPT  = CI*(CSKXL*CEX+RNPR*CEZ)/(CSKXL**2+RNPR**2)

          CEF(NX,1) = CEF(NX,1)+CEX
          CEF(NX,2) = CEF(NX,2)+CEY
          CEF(NX,3) = CEF(NX,3)+CEZ

          CED(NX,1) = CED(NX,1)+CI*CSKXL*CEX
          CED(NX,2) = CED(NX,2)+CI*CSKXL*CEY
          CED(NX,3) = CED(NX,3)+CI*CSKXL*CEZ

          CBF(NX,1) = CBF(NX,1)-RNPR*CEY
          CBF(NX,2) = CBF(NX,2)+RNPR*CEX-CSKXL*CEZ
          CBF(NX,3) = CBF(NX,3)+CSKXL*CEY

          CAF(NX,1) = CAF(NX,1)-CI*CEX+CSKXL*CPT
          CAF(NX,2) = CAF(NX,2)-CI*CEY
          CAF(NX,3) = CAF(NX,3)-CI*CEZ+RNPR*CPT
          CAF(NX,4) = CAF(NX,4)+CPT

          DO NS=1,NSMAX 
             CJX=-CI*( CM0(1,NX,NS)         *CEX &
                      +CM0(2,NX,NS)         *CEY &
                      +CM1(1,NX,NS)   *CSKXL*CEZ)
             CDX=-CI*(+CM2(1,NX,NS)*CI*CSKXL*CEX &
                      +CM2(2,NX,NS)*CI*CSKXL*CEY)
             CJY=-CI*(-CM0(2,NX,NS)         *CEX &
                      +CM0(3,NX,NS)         *CEY)
             CDY=-CI*(+CM1(2,NX,NS)*CI      *CEZ &
                      -CM2(2,NX,NS)*CI*CSKXL*CEX &
                      +CM2(3,NX,NS)*CI*CSKXL*CEY)
             CJZ=-CI*( CM0(4,NX,NS)         *CEZ &
                      -CM1(2,NX,NS)   *CSKXL*CEY)
             CDZ=-CI*(+CM1(1,NX,NS)*CI      *CEX &
                      +CM2(4,NX,NS)*CI*CSKXL*CEZ)
             CVX=-CI*( CM0(1,NX,NS)            *CEX &
                      +CM0(2,NX,NS)            *CEY &
                      +CM1(1,NX,NS)*CSKXL      *CEZ &
                      +CM2(1,NX,NS)*CSKXL*CSKXL*CEX &
                      +CM2(2,NX,NS)*CSKXL*CSKXL*CEY)
             CVZ=-CI*( CM0(4,NX,NS)            *CEZ &
                      +CM1(1,NX,NS)*CSKXL      *CEX &
                      -CM1(2,NX,NS)*CSKXL      *CEY &
                      +CM2(4,NX,NS)*CSKXL*CSKXL*CEZ)

             CJF(NX,NS,1)=CJF(NX,NS,1)+CJX
             CJF(NX,NS,2)=CJF(NX,NS,2)+CJY
             CJF(NX,NS,3)=CJF(NX,NS,3)+CJZ
             CJF(NX,NS,4)=CJF(NX,NS,4)+CVX
             CJF(NX,NS,5)=CJF(NX,NS,5)+CI*CSKXL*CVX
             CJF(NX,NS,6)=CJF(NX,NS,6)+CVZ
             CJF(NX,NS,7)=CJF(NX,NS,7)+CI*CSKXL*CVZ
             CJD(NX,NS,1)=CJD(NX,NS,1)+CDX
             CJD(NX,NS,2)=CJD(NX,NS,2)+CDY
             CJD(NX,NS,3)=CJD(NX,NS,3)+CDZ
          END DO
       END DO
    END DO

    DO NX=1,NXPMAX
       RHL(NX,1)=DCONJG(CAF(NX,1))*CBF(NX,1) &
                +DCONJG(CAF(NX,2))*CBF(NX,2) &
                +DCONJG(CAF(NX,3))*CBF(NX,3)
       RHL(NX,2)=DCONJG(CAF(NX,2))*CAF(NX,3)*(-CI) &
                -DCONJG(CAF(NX,3))*CAF(NX,2)*(-CI) &
                +2.D0*DCONJG(CAF(NX,4))*CBF(NX,1)
       RHL(NX,3)=DCONJG(CEF(NX,1))*CBF(NX,1) &
                +DCONJG(CEF(NX,2))*CBF(NX,2) &
                +DCONJG(CEF(NX,3))*CBF(NX,3)

       RHL(NX,4)=0.D0
       RHL(NX,5)=0.D0

       DO NS=1,NSMAX
          FWP= 1.D20*AEE*AEE*PZ(NS)*PZ(NS)/(AMP*PA(NS)*EPS0*RW*RW)
          IF(PROFPN(NX,NS).LE.0.D0) THEN
             WP2I=1.D0
          ELSE
             WP2I=1.D0/(FWP*PROFPN(NX,NS))
          ENDIF
          FHL(NX,NS,1)=RNPR*(DCONJG(CEF(NX,1))*CJF(NX,NS,1) &
                            +DCONJG(CEF(NX,2))*CJF(NX,NS,2) &
                            +DCONJG(CEF(NX,3))*CJF(NX,NS,3) &
                            +DCONJG(CED(NX,1))*CJD(NX,NS,1) &
                            +DCONJG(CED(NX,2))*CJD(NX,NS,2) &
                            +DCONJG(CED(NX,3))*CJD(NX,NS,3))
          FHL(NX,NS,2)=CI*(DCONJG(CJF(NX,NS,4))*CED(NX,3) &
                          +DCONJG(CJF(NX,NS,5))*CEF(NX,3))
!XX                   +CI*(WP2I*DCONJG(CJF(NX,NS,4))*CJF(NX,NS,7) &
!XX                       +WP2I*DCONJG(CJF(NX,NS,5))*CJF(NX,NS,6))

          RHL(NX,4)=RHL(NX,4)+       FHL(NX,NS,2)
          RHL(NX,5)=RHL(NX,5)-PZ(NS)*FHL(NX,NS,2)
       END DO
    END DO

    TAU=(4.D0*PI*(EPS0/AEE)**2*(AME/AEE)**2)/(1.D20)
    VPH=2.D0*PI*RF*1.D6/RKZ
    FACTC=3.D0*SQRT(0.5D0*PI)/(ZEFF*(0.29+0.46/(1.08+ZEFF)))
    AWEC=2.D0*PI*RF*1.D6*EPS0*AEE/(VC*AMP)

    DO NX=1,NXPMAX
       XR=XAM(NX)/RR
       YR=0.1D0/RR
       EP=SQRT(XR*XR+YR*YR)
       RLNLMD=16.1D0 - 1.15D0*LOG10(PROFPN(NX,1)) &
                     + 2.30D0*LOG10(PROFTR(NX,1))
       VTE=SQRT(PROFTR(NX,1)*AEE*1.D3/AME)
       TAU0=TAU*VTE**3/(PROFPN(NX,1)*RLNLMD)
       FRAC=1.D0-1.9D0*SQRT(EP)+ABS(EP)
 
       NS=1
       W=VPH/VTE
       TAU1=TAU0*W*W1CDEF(W,ZEFF,0.D0,0.D0,NCDTYP)
       TAU2=TAU0*FACTC
       TAU3=TAU0*W*W1CDEF(W,ZEFF,XR,YR,NCDTYP)
       TAU4=TAU0*FACTC*FRAC
       AHL(NX,NS,1)=-AWEC*(PZ(NS)/PA(NS))*FHL(NX,NS,1)*TAU1
       AHL(NX,NS,2)=-AWEC*(PZ(NS)/PA(NS))*FHL(NX,NS,2)*TAU2
       AHL(NX,NS,3)=-AWEC*(PZ(NS)/PA(NS))*FHL(NX,NS,1)*TAU3
       AHL(NX,NS,4)=-AWEC*(PZ(NS)/PA(NS))*FHL(NX,NS,2)*TAU4

       DO NS=2,NSMAX
          TAUE=TAU0*FACTC*SQRT(PROFTR(NX,NS)/PROFTR(NX,1))**3 &
               *SQRT(2.D0*PA(NS)*AMP/AME)/(ZEFF*PZ(NS)**2) 
          TAUI=TAU0*FACTC*PA(NS)*AMP/(AME*PZ(NS)**2)
          TAU1=(1.D0-PZ(NS)/ZEFF)/(1.D0/TAUE+1.D0/TAUI)
          TAU2=TAU1
          TAU3=(1.D0-PZ(NS)/ZEFF)/(1.D0/TAUE+1.D0/TAUI)
          TAU4=TAU1
          AHL(NX,NS,1)=-AWEC*(PZ(NS)/PA(NS))*FHL(NX,NS,1)*TAU1
          AHL(NX,NS,2)=-AWEC*(PZ(NS)/PA(NS))*FHL(NX,NS,2)*TAU2
          AHL(NX,NS,3)=-AWEC*(PZ(NS)/PA(NS))*FHL(NX,NS,1)*TAU3
          AHL(NX,NS,4)=-AWEC*(PZ(NS)/PA(NS))*FHL(NX,NS,2)*TAU4
       END DO
    END DO

    DO NS=1,NSMAX
       AHLT(NS,1)=0.D0
       AHLT(NS,2)=0.D0
       AHLT(NS,3)=0.D0
       AHLT(NS,4)=0.D0
       DO NX=1,NXPMAX
          AHLT(NS,1)=AHLT(NS,1)+AHL(NX,NS,1)*(XA(NX+1)-XA(NX))
          AHLT(NS,2)=AHLT(NS,2)+AHL(NX,NS,2)*(XA(NX+1)-XA(NX))
          AHLT(NS,3)=AHLT(NS,3)+AHL(NX,NS,3)*(XA(NX+1)-XA(NX))
          AHLT(NS,4)=AHLT(NS,4)+AHL(NX,NS,4)*(XA(NX+1)-XA(NX))
       END DO
    END DO

!     RHL(,6)=EPARA=HDOT*C/(B0*C)    V/M FOR MW
!     RHL(,7)=EPARA=FE2*NE*ME/QE     V/M FOR MW

    A6=1.D6/(BB*VC)
    A7=1.D6*2.D0*PI*RF*1.D6*EPS0/(VC*1.D20*AEE)
    DO NX=1,NXPMAX
       RHL(NX,6)=A6*RHL(NX,3)/PROFB(NX)
       RHL(NX,7)=A7*FHL(NX,1,2)/PROFPN(NX,1)
    END DO
    RETURN
  END SUBROUTINE W1HELD
END MODULE w1rslt
