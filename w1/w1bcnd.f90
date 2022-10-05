MODULE w1bcnd

CONTAINS

!     ****** SET BOUNDARY CONDITIONS AT R=RA ******

  SUBROUTINE W1_BCND
    USE w1comm
    IMPLICIT NONE
    INTEGER:: I,J
    REAL(rkind):: RW,RKV,RKPR,RNPR,RCE
    COMPLEX(rkind):: CKKV2,CKKV,CFCTD,CFCTB,CKKW2,CKKW,CKKU2,CRE,CRM
    COMPLEX(rkind):: CRJ,CFCJ,CJY1,CJY2,CJZ1,CJZ2

    WRITE(6,*) '@@@ point 2400'
    RW=2.D6*PI*RF
    RKV=RW/VC
    RKPR=RKZ
    RNPR=VC*RKPR/RW
    RCE=VC*EPS0
    CKKV2=RNPR*RNPR-1.D0
    CKKV =SQRT(CKKV2)
    CFCTD=EXP(-2.D0*CKKV*RKV*(RD-RA))
    CFCTB=EXP(-2.D0*CKKV*RKV*(RB-RA))
    IF(ABS(WALLR).GT.1.D-12) THEN
       CKKW2=RNPR*RNPR-1.D0-DCMPLX(0.D0,1.D0/(RCE*RKV*WALLR))
       CKKW =SQRT(CKKW2)
       CKKU2=RNPR*RNPR-CKKW2
       CRE=CFCTB*(CKKV-CKKW)/(CKKV+CKKW)
       CRM=CFCTB*(CKKW-CKKV*CKKU2)/(CKKW+CKKV*CKKU2)
    ELSE 
       CRE=-CFCTB
       CRM=-CFCTB
    ENDIF
    CRJ=-CFCTD
    CFCJ=(0.D0,-0.5D0)*EXP(CKKV*RKV*(RD-RA))/(RCE*CKKV)
    CJY1=CFCJ*CFJY1
    CJY2=CFCJ*CFJY2
    CJZ1=CFCJ*CFJZ1
    CJZ2=CFCJ*CFJZ2

    WRITE(6,*) '@@@ point 2401'

    DO I=1,3
       DO J=1,5
          CGIN(I,J)=(0.D0,0.D0)
          CGOT(I,J)=(0.D0,0.D0)
       END DO
    END DO

    WRITE(6,*) '@@@ point 2402'
    IF(MOD(NMODEL,2).EQ.1) THEN
       CGIN(1,1)=     1.D0+CRE
       CGIN(1,2)=-CI*(1.D0-CRE)*CKKV
       CGIN(2,3)=     1.D0+CRM
       CGIN(2,4)=-CI*(1.D0-CRM)/CKKV
       CGIN(3,1)=    (1.D0+CRJ)     *CJY1
       CGIN(3,2)=-CI*(1.D0-CRJ)*CKKV*CJY1
       CGIN(3,3)=    (1.D0+CRJ)     *CJZ1
       CGIN(3,4)=-CI*(1.D0-CRJ)/CKKV*CJZ1
       CGOT(1,1)=     1.D0+CRE
       CGOT(1,2)= CI*(1.D0-CRE)*CKKV
       CGOT(2,3)=     1.D0+CRM
       CGOT(2,4)= CI*(1.D0-CRM)/CKKV
       CGOT(3,1)=    (1.D0+CRJ)     *CJY2
       CGOT(3,2)= CI*(1.D0-CRJ)*CKKV*CJY2
       CGOT(3,3)=    (1.D0+CRJ)     *CJZ2
       CGOT(3,4)= CI*(1.D0-CRJ)*CKKV*CJZ2
    ELSE
       CGIN(1,1)=  1.D0+CRE
       CGIN(1,2)=-(1.D0-CRE)*CKKV
       CGIN(2,3)=  1.D0+CRM
       CGIN(2,4)=-(1.D0-CRM)/CKKV
       CGIN(3,1)= (1.D0+CRJ)     *CJY1
       CGIN(3,2)=-(1.D0-CRJ)*CKKV*CJY1
       CGIN(3,3)= (1.D0+CRJ)     *CJZ1
       CGIN(3,4)=-(1.D0-CRJ)/CKKV*CJZ1
       CGOT(1,1)=  1.D0+CRE
       CGOT(1,2)=-(1.D0-CRE)*CKKV
       CGOT(2,3)=  1.D0+CRM
       CGOT(2,4)=-(1.D0-CRM)/CKKV
       CGOT(3,1)= (1.D0+CRJ)     *CJY2
       CGOT(3,2)=-(1.D0-CRJ)*CKKV*CJY2
       CGOT(3,3)= (1.D0+CRJ)     *CJZ2
       CGOT(3,4)=-(1.D0-CRJ)*CKKV*CJZ2
    END IF
    RETURN
  END SUBROUTINE W1_BCND

!     ******* ELECTROMAGNETIC FIELD IN VACUUM *******

  SUBROUTINE W1EVAC(NZ)
    USE w1comm
    IMPLICIT NONE
    INTEGER,INTENT(IN):: NZ
    INTEGER:: NSF,NXVH,I,NX1,NX2
    REAL(rkind):: RW,RKV,RCE,FACT,RKPR,RNPR,X,RJ,DX
    COMPLEX(rkind):: CKKV2,CKKV,CFCTD,CFCTB,CKKW2,CKKW,CKKU2,CRE,CRM
    COMPLEX(rkind):: CRJ,CFCJ,CJY1,CJY2,CJZ1,CJZ2,CEAM,CEAP
    COMPLEX(rkind):: CAY1,CAY2,CAZ1,CAZ2
    COMPLEX(rkind):: CEX1,CEY1,CEZ1,CEX1DX,CEY1DX,CEZ1DX
    COMPLEX(rkind):: CEX2,CEY2,CEZ2,CEX2DX,CEY2DX,CEZ2DX
    COMPLEX(rkind):: CJ2DX,CJ2DY,CJ2DZ

    RW=2.D6*PI*RF
    RKV=RW/VC
    RCE=VC*EPS0

    IF(NSYM.NE.0.AND.NZ.NE.1.AND.NZ.NE.(NZMAX/2+1)) THEN
       FACT = 2.0*RZ
    ELSEIF(NSYM.EQ.-1.AND.NZ.EQ.(NZMAX/2+1)) THEN
       FACT = 0.0
    ELSE
       FACT =     RZ
    ENDIF
    RKPR=RKZ
    RNPR=VC*RKPR/RW
    CKKV2=RNPR*RNPR-1.D0
    CKKV =SQRT(CKKV2)
    CFCTD=EXP(-2.D0*CKKV*RKV*(RD-RA))
    CFCTB=EXP(-2.D0*CKKV*RKV*(RB-RA))
    IF(ABS(WALLR).GT.1.D-12) THEN
       CKKW2=RNPR*RNPR-1.D0-DCMPLX(0.D0,1.D0/(RCE*RKV*WALLR))
       CKKW =SQRT(CKKW2)
       CKKU2=RNPR*RNPR-CKKW2
       CRE=CFCTB*(CKKV-CKKW)/(CKKV+CKKW)
       CRM=CFCTB*(CKKW-CKKV*CKKU2)/(CKKW+CKKV*CKKU2)
    ELSE
       CRE=-CFCTB
       CRM=-CFCTB
    ENDIF
    CRJ=-CFCTD
    CFCJ=(0.D0,-0.5D0)*EXP(CKKV*RKV*(RD-RA))/(RCE*CKKV)
    CJY1=CFCJ*CFJY1
    CJY2=CFCJ*CFJY2
    CJZ1=CFCJ*CFJZ1
    CJZ2=CFCJ*CFJZ2

    IF(NMODEL.EQ.1.OR.NMODEL.EQ.3) THEN
       NSF=4*NXPMAX+4
    ELSEIF(NMODEL.EQ.5) THEN
       NSF=6*NXPMAX+4
    ELSE
       NSF=3*NXPMAX+4
    ENDIF
    NXVH=NXVMAX/2

    DO I=1,NXVMAX
       NX1=NXPMAX+NXVMAX+I
       NX2=NXPMAX+NXVMAX+1-I
       X=XAM(NX1)
       CEAM=EXP( CKKV*RKV*(X+RA))
       CEAP=EXP(-CKKV*RKV*(X+RA))

       IF(I.LE.NXVH) THEN
          DX=RKV*(RB-RD)/DBLE(NXVH-1)
          CAY1=(0.D0,0.D0)
          CAY2=(0.D0,0.D0)
          CAZ1=(0.D0,0.D0)
          CAZ2=(0.D0,0.D0)
       ELSE
          DX=RKV*(RD-RA)/DBLE(NXVH-1)
          CAY1=CJY1
          CAY2=CJY2
          CAZ1=CJZ1
          CAZ2=CJZ2
       ENDIF

       CEX1  =CA(2    )*( 0.D0,-1.D0)*RNPR/CKKV*(CEAM-CRM*CEAP) &
             +CAZ1     *( 0.D0,-1.D0)*RNPR/CKKV*(CEAM-CRJ*CEAP)
       CEY1  =CA(1    )                        *(CEAM+CRE*CEAP) &
             +CAY1                             *(CEAM+CRJ*CEAP)
       CEZ1  =CA(2    )                        *(CEAM+CRM*CEAP) &
             +CAZ1                             *(CEAM+CRJ*CEAP)
       CEX1DX=CA(2    )*(-1.D0, 0.D0)*RNPR     *(CEAM+CRM*CEAP) &
             +CAZ1     *(-1.D0, 0.D0)*RNPR     *(CEAM+CRJ*CEAP)
       CEY1DX=CA(1    )*( 0.D0,-1.D0)*CKKV     *(CEAM-CRE*CEAP) &
             +CAY1     *( 0.D0,-1.D0)*CKKV     *(CEAM-CRJ*CEAP)
       CEZ1DX=CA(2    )*( 0.D0,-1.D0)*CKKV     *(CEAM-CRM*CEAP) &
             +CAZ1     *( 0.D0,-1.D0)*CKKV     *(CEAM-CRJ*CEAP)
       CEX2  =CA(NSF  )*( 0.D0, 1.D0)*RNPR/CKKV*(CEAM-CRM*CEAP) &
             +CAZ2     *( 0.D0, 1.D0)*RNPR/CKKV*(CEAM-CRJ*CEAP)
       CEY2  =CA(NSF-1)                        *(CEAM+CRE*CEAP) &
             +CAY2                             *(CEAM+CRJ*CEAP)
       CEZ2  =CA(NSF  )                        *(CEAM+CRM*CEAP) &
             +CAZ2                             *(CEAM+CRJ*CEAP)
       CEX2DX=CA(NSF  )*(-1.D0, 0.D0)*RNPR     *(CEAM+CRM*CEAP) &
             +CAZ2     *(-1.D0, 0.D0)*RNPR     *(CEAM+CRJ*CEAP)
       CEY2DX=CA(NSF-1)*( 0.D0, 1.D0)*CKKV     *(CEAM-CRE*CEAP) &
             +CAY2     *( 0.D0, 1.D0)*CKKV     *(CEAM-CRJ*CEAP)
       CEZ2DX=CA(NSF  )*( 0.D0, 1.D0)*CKKV     *(CEAM-CRM*CEAP) &
             +CAZ2     *( 0.D0, 1.D0)*CKKV     *(CEAM-CRJ*CEAP)

       CE2DA(NZ,NX1,1)=CEX1
       CE2DA(NZ,NX1,2)=CEY1
       CE2DA(NZ,NX1,3)=CEZ1
       CJ2DX          = 0.5D0*RNPR*CEZ1
       CJ2DY          =                -CEY1DX
       CJ2DZ          = 0.5D0*RNPR*CEX1-CEZ1DX
       FLUXX(NX1)=FLUXX(NX1) - FACT*( DCONJG(CEX1)*CJ2DX &
                                    + DCONJG(CEY1)*CJ2DY &
                                    + DCONJG(CEZ1)*CJ2DZ ) * RCE

       CE2DA(NZ,NX2,1)=CEX2
       CE2DA(NZ,NX2,2)=CEY2
       CE2DA(NZ,NX2,3)=CEZ2
       CJ2DX          = 0.5D0*RNPR*CEZ2
       CJ2DY          =                -CEY2DX
       CJ2DZ          = 0.5D0*RNPR*CEX2-CEZ2DX
       FLUXX(NX2)=FLUXX(NX2) - FACT*( DCONJG(CEX2)*CJ2DX &
                                    + DCONJG(CEY2)*CJ2DY &
                                    + DCONJG(CEZ2)*CJ2DZ ) * RCE
    END DO

    RETURN
  END SUBROUTINE W1EVAC
END MODULE w1bcnd
