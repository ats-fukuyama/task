MODULE w1disp

CONTAINS

!     ******* LOCAL PLASMA PARAMETERS *******

  SUBROUTINE W1DSPA
    USE w1comm
    USE w1alfa
    USE libdsp
    USE w1lib
    IMPLICIT NONE

!      CM0 , CD0         CM1 , CD1          CM2 , CD2
!      1   2   0         0   0   1          1   2   0
!    -(2)  3   0         0   0   2        -(2)  3   0
!      0   0   4         1 -(2)  0          0   0   4

    COMPLEX(rkind),DIMENSION(:,:),ALLOCATABLE:: CGZ,CZ,CDZ,CDDZ,CDDDZ
    REAL(rkind),DIMENSION(:),ALLOCATABLE:: RKPP,ALAM
    INTEGER:: NX,NC,NN,IHMIN,NCMAXS,NS
    REAL(rkind):: RT2,RW,VALPHA,RKPR,RNPR,FWP,FWC,FVT,WP,WC,RNPP2
    REAL(rkind):: VTE,VCRIT3,VCRIT,VRES,VLIM,Y1,Y2,Y3,COEF,UD,AKPR,ARG
    REAL(rkind):: RT,RL,ALAMc,ALAMM,ALAMP,ALAMD,ALAMN,TPR,RNT,TP
    REAL(rkind):: CC1,CC2,CC3,CC4,CC5,CC6,CC7,CC8,CC9,CC0
    COMPLEX(rkind):: CW1,CW2,CW3,CW4,CW5,CF2,CF3,CF4,CG2,CG3,CG4
    COMPLEX(rkind):: CDT,CDX,CA1,CA2,CA3

    ALLOCATE(RKPP(NXTMAX))
    ALLOCATE(CGZ(NXPMAX,NCMAX),CZ(NXPMAX,NCMAX),CDZ(NXPMAX,NCMAX))
    ALLOCATE(CDDZ(NXPMAX,NCMAX),CDDDZ(NXPMAX,NCMAX))
    ALLOCATE(ALAM(0:NHMAX+1))

    RT2= SQRT(2.D0)
    RW  = 2.D6*PI*RF
    VALPHA=SQRT(2.D0*3.5D6*AEE/(4.D0*AMP))

    IF(NMODEL.LE.1) THEN
       RKPR=RKZ
       RNPR=VC*RKPR/RW
       DO NX=1,NXPMAX
          RKPP(NX)=0.D0
       END DO
    ELSE
       RKPR=RKZ
       RNPR=VC*RKPR/RW
       DO NX=1,NXPMAX
          CDT =1.D0
          CDX =0.D0
          DO NS=1,NSMAX
             FWP= 1.D20*AEE*AEE*PZ(NS)*PZ(NS)/(AMP*PA(NS)*EPS0)
             FWC= AEE*PZ(NS)*BB/(AMP*PA(NS))
             WP = FWP*PROFPN(NX,NS)
             WC = FWC*PROFB(NX)
             CDT=CDT-WP/(RW*RW-WC*WC)
             CDX=CDX-WP*WC/((RW*RW-WC*WC)*RW)
          END DO
          RNPP2=((CDT-RNPR*RNPR)**2-CDX**2)/(CDT-RNPR*RNPR)
          RKPP(NX)= RW*SQRT(MAX(RNPP2,0.D0))/VC
       END DO
    END IF

    DO NS=1,NSMAX
       FWP= 1.D20*AEE*AEE*PZ(NS)*PZ(NS)/(AMP*PA(NS)*EPS0*RW*RW)
       FWC= AEE*PZ(NS)*BB/(AMP*PA(NS))
       FVT= AEE*1.D3/(AMP*PA(NS))
       DO NX=1,NXPMAX
          CM0(1,NX,NS) = 0.D0
          CM0(2,NX,NS) = 0.D0
          CM0(3,NX,NS) = 0.D0
          CM0(4,NX,NS) = 0.D0
          CM1(1,NX,NS) = 0.D0
          CM1(2,NX,NS) = 0.D0
          CM2(1,NX,NS) = 0.D0
          CM2(2,NX,NS) = 0.D0
          CM2(3,NX,NS) = 0.D0
          CM2(4,NX,NS) = 0.D0
       END DO

       IF(MOD(NALPHA,2).EQ.1.AND.NS.EQ.4) THEN
          DO NX=1,NXPMAX
             RKPR= RKZ
             WP  = FWP*PROFPN(NX,NS)
             WC  = FWC*PROFB(NX)
             VTE   =SQRT(2.D0*PROFTR(NX,1)*1.D3*AEE/AME)
             VCRIT3=.75D0*SQRT(PI)*AME*(PROFPN(NX,2)/(AMP*PA(2)) &
                                       +PROFPN(NX,3)/(AMP*PA(3))) &
                                      /PROFPN(NX,1)
             VCRIT =VTE*VCRIT3**(1.D0/3.D0)
             DO NC=1,2*ABS(IHARM(4))+1
                NN  = NC-ABS(IHARM(4))-1
                VRES=(RW-NN*WC)/RKPR
                IF(VALPHA**2-VRES**2.GT.0.D0.AND. &
                   RKPP(NX).GT.0.D0) THEN
                   VLIM=SQRT(VALPHA**2-VRES**2)
                   Y1=RKPP(NX)*VCRIT/WC
                   Y2=VRES/VCRIT
                   Y3=VLIM/VCRIT
                   COEF=4.5D0*PI*WP*RW &
                        /(ABS(RKPR)*VCRIT*LOG(1.D0+(VALPHA/VCRIT)**3))
                   CM0(1,NX,NS)=CM0(1,NX,NS) &
                               +CI*COEF*W1FALF(Y1,Y2,Y3,NN,1)
                   CM0(2,NX,NS)=CM0(2,NX,NS) &
                               -   COEF*W1FALF(Y1,Y2,Y3,NN,2)
                   CM1(1,NX,NS)=CM1(1,NX,NS) &
                               +CI*COEF*W1FALF(Y1,Y2,Y3,NN,3)*RW/(RKPP(NX)*VC)
                   CM0(3,NX,NS)=CM0(3,NX,NS) &
                               +CI*COEF*W1FALF(Y1,Y2,Y3,NN,4)
                   CM1(2,NX,NS)=CM1(2,NX,NS) &
                               +   COEF*W1FALF(Y1,Y2,Y3,NN,5)*RW/(RKPP(NX)*VC)
                   CM0(4,NX,NS)=CM0(4,NX,NS) &
                               +CI*COEF*W1FALF(Y1,Y2,Y3,NN,6)
                END IF
             END DO
          END DO
       ELSE
          IF(NMODEL.LE.3.OR.IHARM(NS).LT.0.OR.IHARM(NS).GE.3) THEN
             IF(NMODEL.LE.3.OR.IHARM(NS).LT.0) THEN
                IHMIN=0
             ELSE
                IHMIN=3
             END IF
             NCMAXS=2*ABS(IHARM(NS))+1
             DO NX=1,NXPMAX
                RKPR= RKZ
                WC  = FWC*PROFB(NX)
                UD  = SQRT(FVT*PROFPU(NX,NS))
                AKPR= RT2*ABS(RKPR)*SQRT(FVT*PROFTR(NX,NS))
                DO NC=1,NCMAXS
                   NN= NC-ABS(IHARM(NS))-1
                   ARG = (RW-NN*WC)/AKPR
                   CGZ(NX,NC)= ARG
                   CZ(NX,NC)= ARG
                END DO
             END DO
             DO NC=1,NCMAXS
                CALL DSPFNV(NXPMAX,CGZ(1:NXPMAX,NC), &
                                   CZ(1:NXPMAX,NC),CDZ(1:NXPMAX,NC), &
                                   CDDZ(1:NXPMAX,NC),CDDDZ(1:NXPMAX,NC))
             END DO
             DO NX=1,NXPMAX
                RKPR= RKZ
                WP  = FWP*PROFPN(NX,NS)
                WC  = FWC*PROFB(NX)
                UD  = SQRT(FVT*PROFPU(NX,NS))
                AKPR= RT2*ABS(RKPR)*SQRT(FVT*PROFTR(NX,NS))
                RT  = PROFTP(NX,NS)/PROFTR(NX,NS)
                RL  = RKPP(NX)*RKPP(NX)*FVT*PROFTP(NX,NS)/(WC*WC)
                CALL LAMBDA(ABS(IHARM(NS))+1,RL,ALAM)
                DO NC=1,NCMAXS
                   NN=NC-ABS(IHARM(NS))-1
                   IF(ABS(NN).GE.IHMIN) THEN
                      CA1 = RW*CZ(NX,NC)/AKPR+0.5D0*(1.D0-RT)*CDZ(NX,NC)
                      CA2 = 0.5D0*(RT*RW/WC-NN*(RT-1.D0))*CDZ(NX,NC)
                      CA3 =-(RW-NN*WC*(1.D0-1.D0/RT))*CGZ(NX,NC)*CDZ(NX,NC) &
                            /AKPR
                      ALAMC=ALAM(ABS(NN))
                      ALAMM=ALAM(ABS(NN-1))
                      ALAMP=ALAM(ABS(NN+1))
                      ALAMD=0.5D0*(ALAMM+ALAMP)-ALAMC
                      ALAMN=0.5D0*(ALAMM-ALAMP)
                      CM0(1,NX,NS)=CM0(1,NX,NS) &
                                  +WP*NN*ALAMN*CA1
                      CM0(2,NX,NS)=CM0(2,NX,NS) &
                                  +WP*CI*NN*ALAMD*CA1
                      CM1(1,NX,NS)=CM1(1,NX,NS) &
                                  -WP*ALAMN*CA2*RW/(RKPR*VC)
                      CM0(3,NX,NS)=CM0(3,NX,NS) &
                                  +WP*(NN*ALAMN-2.0D0*RL*ALAMD)*CA1
                      CM1(2,NX,NS)=CM1(2,NX,NS) &
                                  +WP*CI*ALAMD*CA2*RW/(RKPR*VC)
                      CM0(4,NX,NS)=CM0(4,NX,NS) &
                                  +WP*ALAMC*CA3
                   END IF
                END DO
             END DO
          END IF
          IF(NMODEL.GT.3.AND.IHARM(NS).GE.0) THEN
             DO NX = 1 , NXPMAX
                RKPR = RKZ
                WC   = FWC*PROFB(NX)
                UD   = SQRT(FVT*PROFPU(NX,NS))
                AKPR = RT2*ABS(RKPR)*SQRT(FVT*PROFTR(NX,NS))
                CGZ(NX,1) = (RW-(-2.D0)*WC-RKPR*UD)/AKPR
                CGZ(NX,2) = (RW-(-1.D0)*WC-RKPR*UD)/AKPR
                CGZ(NX,3) = (RW           -RKPR*UD)/AKPR
                CGZ(NX,4) = (RW-( 1.D0)*WC-RKPR*UD)/AKPR
                CGZ(NX,5) = (RW-( 2.D0)*WC-RKPR*UD)/AKPR
             END DO

             DO NC=1,5
                CALL DSPFNV(NXPMAX,CGZ(1:NXPMAX,NC), &
                                   CZ(1:NXPMAX,NC),CDZ(1:NXPMAX,NC), &
                                   CDDZ(1:NXPMAX,NC),CDDDZ(1:NXPMAX,NC))
             END DO
             DO NX=1,NXPMAX
                RKPR= RKZ
                RNPR= VC*RKPR/RW
                WP  = FWP*PROFPN(NX,NS)
                WC  = FWC*PROFB(NX)
                UD  = SQRT(FVT*PROFPU(NX,NS))
                AKPR= RT2*ABS(RKPR)*SQRT(FVT*PROFTR(NX,NS))
                TPR = FVT*PROFTR(NX,NS)
                RNT = PROFTP(NX,NS)/PROFTR(NX,NS)
                TP  = WP*FVT*PROFTP(NX,NS)*(RW/(VC*WC))**2

                CC1 = (RW-RKPR*UD)/AKPR
                CC2 = 0.5D0*(1.D0-RNT)
                CC3 =-RKPR*UD/AKPR
                CC4 = 0.5D0*RNT*RW/WC
                CC5 = RW/AKPR
                CC6 =-WC/AKPR*(1.D0-(1.D0+UD*UD/TPR)/RNT)
                CC7 =-UD*RW/(2.D0*RKPR*TPR)
                CC8 = UD*WC*(1.D0-2.D0/RNT)/(2.D0*RKPR*TPR)
                CC9 = 0.5D0*RW/AKPR
                CC0 =-0.5D0*WC*(1.D0-1.D0/RNT)/AKPR

                CW1 = CC1*CZ(NX,1) + CC2*CDZ(NX,1)
                CW2 = CC1*CZ(NX,2) + CC2*CDZ(NX,2)
                CW3 = CC1*CZ(NX,3) + CC2*CDZ(NX,3)
                CW4 = CC1*CZ(NX,4) + CC2*CDZ(NX,4)
                CW5 = CC1*CZ(NX,5) + CC2*CDZ(NX,5)

                CF2 =-1.D0*CC3*CZ(NX,2) + (CC4-CC2)*CDZ(NX,2)
                CF3 =                     (CC4    )*CDZ(NX,3)
                CF4 = 1.D0*CC3*CZ(NX,4) + (CC4+CC2)*CDZ(NX,4)

                CG2 = (CC5-CC6)*CZ(NX,2) + (CC7-CC8)*CDZ(NX,2) &
                    + (CC9-CC0)*CDDZ(NX,2)
                CG3 = (CC5    )*CZ(NX,3) + (CC7    )*CDZ(NX,3) &
                    + (CC9    )*CDDZ(NX,3)
                CG4 = (CC5+CC6)*CZ(NX,4) + (CC7+CC8)*CDZ(NX,4) &
                    + (CC9+CC0)*CDDZ(NX,4)

                CM0(1,NX,NS)=CM0(1,NX,NS) &
                            +    WP*0.5D0*(CW4+CW2)
                CM0(2,NX,NS)=CM0(2,NX,NS) &
                            + CI*WP*0.5D0*(CW4-CW2)
                CM0(3,NX,NS)=CM0(3,NX,NS) &
                            +    WP*0.5D0*(CW4+CW2)
                CM0(4,NX,NS)=CM0(4,NX,NS) &
                            +    WP      * CG3
                CM1(1,NX,NS)=CM1(1,NX,NS) &
                            -    WP*0.5D0*(CF4-CF2)/RNPR
                CM1(2,NX,NS)=CM1(2,NX,NS) &
                            + CI*WP*0.5D0*(CF4-2.D0*CF3+CF2)/RNPR
                CM2(1,NX,NS)=CM2(1,NX,NS) &
                            +    TP*0.5D0*(CW5-CW4-CW2+CW1)
                CM2(2,NX,NS)=CM2(2,NX,NS) &
                            + CI*TP*0.5D0*(CW5-2.D0*CW4+2.D0*CW2-CW1)
                CM2(3,NX,NS)=CM2(3,NX,NS) &
                            +    TP*0.5D0*(CW5-3.D0*CW4+4.D0*CW3-3.D0*CW2+CW1)
                CM2(4,NX,NS)=CM2(4,NX,NS) &
                            +    TP*0.5D0*(CG4-2.D0*CG3+CG2)
             END DO
          END IF
       END IF
    END DO
    
    DEALLOCATE(RKPP)
    DEALLOCATE(CGZ,CZ,CDZ,CDDZ,CDDDZ)
    DEALLOCATE(ALAM)

    DO NX = 1 , NXPMAX
       CD0( 1 , NX ) =  1.D0 - RNPR*RNPR
       CD0( 2 , NX ) =  0.D0
       CD0( 3 , NX ) =  1.D0 - RNPR*RNPR
       CD0( 4 , NX ) =  1.D0
       CD1( 1 , NX ) =  RNPR
       CD1( 2 , NX ) =  0.D0
       CD2( 1 , NX ) =  0.D0
       CD2( 2 , NX ) =  0.D0
       CD2( 3 , NX ) = -1.D0
       CD2( 4 , NX ) = -1.D0
    END DO

    DO NS = 1 , NSMAX
       DO NX = 1 , NXPMAX
          CD0( 1 , NX ) = CD0( 1 , NX ) + CM0( 1 , NX , NS )
          CD0( 2 , NX ) = CD0( 2 , NX ) + CM0( 2 , NX , NS )
          CD0( 3 , NX ) = CD0( 3 , NX ) + CM0( 3 , NX , NS )
          CD0( 4 , NX ) = CD0( 4 , NX ) + CM0( 4 , NX , NS )
          CD1( 1 , NX ) = CD1( 1 , NX ) + CM1( 1 , NX , NS )
          CD1( 2 , NX ) = CD1( 2 , NX ) + CM1( 2 , NX , NS )
          CD2( 1 , NX ) = CD2( 1 , NX ) + CM2( 1 , NX , NS )
          CD2( 2 , NX ) = CD2( 2 , NX ) + CM2( 2 , NX , NS )
          CD2( 3 , NX ) = CD2( 3 , NX ) + CM2( 3 , NX , NS )
          CD2( 4 , NX ) = CD2( 4 , NX ) + CM2( 4 , NX , NS )
       END DO
    END DO
    RETURN
  END SUBROUTINE W1DSPA
END MODULE w1disp
