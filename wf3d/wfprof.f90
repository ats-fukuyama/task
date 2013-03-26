!     $Id$

!     ****** PSI ******

SUBROUTINE WFSPSI(X,Y,Z,PSI)

  use wfcomm
  implicit none
  real(8) :: PSI,X,Y,Z

  IF(MODELB.EQ.0.OR.&
 &   MODELB.EQ.1.OR.&
 &   MODELB.EQ.2.OR.&
 &   MODELB.EQ.3) THEN
     PSI=(X*X+Y*Y)/(RA*RA)
  ENDIF
  RETURN
END SUBROUTINE WFSPSI

!     ****** MAGNETIC FIELD PROFILE ******

SUBROUTINE WFSMAG(IN,BABS,AL)

  use wfcomm
  implicit none
  integer :: IN,I
  real(8) :: BLO(3),AL(3),XC,YC,ZC,BABS
  real(8) :: XT,YT,ZT
  real(8) :: R,Q,PHI,AC(3),BLOT(3)

  XC=XND(IN)
  YC=YND(IN)
  ZC=ZND(IN)
  call RCTORT(XC,YC,ZC,XT,YT,ZT)

  R=sqrt(XT*XT+YT*YT)
  Q0=1.d0
  QA=3.d0
  Q = Q0+(QA-Q0)*(R/RA)**2

  IF(MODELB.EQ.0) THEN
     BLO(1)=BB
     BLO(2)=0.D0
     BLO(3)=0.D0
  ELSEIF(MODELB.EQ.1) THEN
     BLO(1)=0.D0
     BLO(2)=BB
     BLO(3)=0.D0
  ELSEIF(MODELB.EQ.2) THEN
     BLO(1)=0.D0
     BLO(2)=0.D0
     BLO(3)=BB
  ELSEIF(MODELB.eq.3) THEN
     if(MODELG.eq.1) then
        BLOT(1)=-BB*YT/(Q*RR)
        BLOT(2)= BB*XT/(Q*RR)
        BLOT(3)= BB*RR/(RR+XT)
        CALL ATTOAC(BLOT(1),BLOT(2),BLOT(3),ZT,BLO(1),BLO(2),BLO(3))
     else
        BLO(1)=-BB*YT/(Q*RR)
        BLO(2)= BB*XT/(Q*RR)
        BLO(3)= BB*RR/(RR+XT)
     end if
  ENDIF
  
  BABS=0.D0
  DO I=1,3
     BABS=BABS+BLO(I)*BLO(I)
  ENDDO
  BABS=SQRT(BABS)

  DO I=1,3
     if(BABS.eq.0.0) then
        AL(I)=0.0
     else
        AL(I)=BLO(I)/BABS
     end if
  ENDDO
  
  RETURN
END SUBROUTINE WFSMAG

!     ****** DENSITY & TEMPERATURE PROFILE ******

SUBROUTINE WFSDEN(IN,RN,RTPR,RTPP,RZCL)

  use wfcomm
  implicit none
  integer :: IN,NS,NSI
  real(8) :: RN(NSM),RTPR(NSM),RTPP(NSM),RZCL(NSM),FACT,PSI,Z
  real(8) :: TE,TI,RLAMEE,RLAMEI,RLAMII,SN,PNN0,VTE,RNUEE,RNUEI,RNUEN
  real(8) :: RNUE,VTI,RNUIE,RNUII,RNUIN,RNUI
  real(8) :: XT,YT,ZT

  IF(MODELP.EQ.0) THEN
     FACT=1.D0
  ELSE
     if(MODELG.eq.1) then
        CALL RCTORT(XND(IN),YND(IN),ZND(IN),XT,YT,ZT)
        CALL WFSPSI(XT,YT,ZT,PSI)
     else
        CALL WFSPSI(XND(IN),YND(IN),ZND(IN),PSI)
     end if
     IF(PSI.LT.1.D0) THEN
        IF(MODELP.EQ.1) THEN
           FACT=1.D0
        ELSEIF(MODELP.EQ.2) THEN
           FACT=1.D0-PSI
        ELSEIF(MODELP.EQ.3) THEN
           FACT=EXP(-(ZPMAX-ZND(IN))/0.10D0)
        ELSE
           WRITE(6,*) 'XX WDSDEN: UNKNOWN MODELP = ',MODELP
        ENDIF
     ELSE
        FACT=0.D0
     ENDIF
     
     Z=ZND(IN)
     IF(Z.LT.ZPMIN.OR.Z.GT.ZPMAX) THEN
        FACT=0.D0
     ENDIF
  ENDIF
  
  DO NS=1,NSMAX
     RN(NS)  =(PN(NS)  -PNS(NS))*FACT+PNS(NS)
     RTPR(NS)=PTPR(NS)
     RTPP(NS)=PTPP(NS)
  ENDDO
  
  IF(RN(1).GT.0.D0) THEN
     
     TE=(RTPR(1)+2.D0*RTPP(1))*1.D3/3.D0
     TI=(RTPR(2)+2.D0*RTPP(2))*1.D3/3.D0
     RLAMEE= 8.0D0+2.3D0*(LOG10(TE)-0.5D0*LOG10(RN(1)))
     RLAMEI= RLAMEE+0.3D0
     RLAMII=12.1D0+2.3D0*(LOG10(TI)-0.5D0*LOG10(RN(1)))
     SN=1.D-20
     PNN0=PPN0/(PTN0*AEE)
     
     DO NS=1,NSMAX
        IF(PZCL(NS).EQ.0) THEN
           IF(NS.EQ.1) THEN
              TE=(RTPR(1)+2.D0*RTPP(1))*1.D3/3.D0
              VTE=SQRT(2.D0*TE*AEE/AME)
              RNUEE=RN(1)*RLAMEE &
                   &                 /(1.24D-4*SQRT(TE*1.D-3)**3)
              RNUEI=0.D0
              DO NSI=2,NSMAX
                 RNUEI=RNUEI+PZ(NSI)**2*RN(NSI)
              ENDDO
              RNUEI=RNUEI*RLAMEI/(1.51D-4*SQRT(TE*1.D-3)**3)
              RNUEN=PNN0*SN*0.88D0*VTE
              RNUE=RNUEE+RNUEI+RNUEN
              RZCL(NS)=RNUE/(2.D6*PI*RF)
           ELSE
              TI=(RTPR(NS)+2.D0*RTPP(NS))*1.D3/3.D0
              VTI=SQRT(2.D0*TI*AEE/(PA(NS)*AMP))
              RNUIE=PZ(NS)**2*RN(1)*RLAMEI &
                   &                 /(2.00D-1*SQRT(TE*1.D-3)**3*PA(NS))
              RNUII=0.D0
              DO NSI=2,NSMAX
                 RNUII=RNUII+PZ(NSI)**2*RN(NSI)
              ENDDO
              RNUII=RNUII*RLAMII &
                   &                 /(5.31D-3*SQRT(TI*1.D-3)**3*SQRT(PA(NS)))
              RNUIN=PNN0*SN*0.88D0*VTI
              RNUI=RNUIE+RNUII+RNUIN
              RZCL(NS)=RNUI/(2.D6*PI*RF)
           ENDIF
        ELSE
           RZCL(NS)=PZCL(NS)
        ENDIF
     ENDDO
     
  ELSE
     DO NS=1,NSMAX
        RZCL(NS)=0.D0
     ENDDO
  ENDIF
  
!  WRITE(6,*) 'ZND= ',ZND(IN)
!  WRITE(6,*) 'RN = ',RN(1),RN(2)
!  WRITE(6,*) 'RT = ',RTPR(1),RTPR(2)
!  WRITE(6,*) 'RZ = ',RZCL(1),RZCL(2)
!  WRITE(6,*) 'E  = ',RNUE,RNUEE,RNUEI,RNUEN
!  WRITE(6,*) 'I  = ',RNUI,RNUIE,RNUII,RNUIN
!  STOP

  RETURN
END SUBROUTINE WFSDEN
