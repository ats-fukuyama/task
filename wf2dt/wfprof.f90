!     $Id: wfprof.f90,v 1.6 2011/12/15 19:06:55 maruyama Exp $

!     ****** set psi ******

SUBROUTINE WFSPSI(LR,LZ,PSI)

  use wfcomm
  implicit none
  real(8),intent(in) :: LR,LZ
  real(8),intent(out):: PSI

!  MODELB=0
  PSI=(LR*LR+LZ*LZ)/(RA*RA)

  RETURN
END SUBROUTINE WFSPSI

!     ****** set magnetic field ******

SUBROUTINE WFSMAG(R,Z,BABS,AL)

  use wfcomm
  implicit none
  integer :: I
  real(8),intent(in) :: R,Z
  real(8),intent(out):: BABS,AL(3)
  real(8) :: BLO(3),LR,LZ
  real(8) :: L,Q

! L : distance from the center of plasma
! Q : safety factor

! --- initialize ---
  LR = R-RR
  LZ = Z
  L  = sqrt(LR*LR+LZ*LZ)

  Q0 = 1.d0
  QA = 3.d0
  Q  = Q0+(QA-Q0)*(L/RA)**2

! --- set B field at NODE NN ---
  if (MODELB.eq.0) then
     BLO(1) =-BB*LZ/(Q*RR)  ! r   direction
     BLO(2) = BB*RR/(RR+LR) ! phi direction
     BLO(3) = BB*LR/(Q*RR)  ! z   direction
  else
     BLO=0.d0
  end if
     
  BABS=0.d0
  DO I=1,3
     BABS=BABS+BLO(I)*BLO(I)
  ENDDO
  BABS=SQRT(BABS)

  DO I=1,3
     if(BABS.eq.0.d0) then
        AL(I)=0.0
     else
        AL(I)=BLO(I)/BABS
     end if
  ENDDO
  
  RETURN
END SUBROUTINE WFSMAG

!     ****** set density & collision frequency ******

SUBROUTINE WFSDEN(R,Z,RN,RTPR,RTPP,RZCL)

  use wfcomm
  implicit none
  integer :: NS,NSI
  real(8),intent(in) :: R,Z
  real(8),intent(out):: RN(NSM),RTPR(NSM),RTPP(NSM),RZCL(NSM)
  real(8) :: LR,LZ
  real(8) :: TE,TI,RLAMEE,RLAMEI,RLAMII,SN,PNN0,VTE,RNUEE,RNUEI,RNUEN
  real(8) :: RNUE,VTI,RNUIE,RNUII,RNUIN,RNUI,FACT,PSI

  ! --- set FACT ---

  if(MODELP.eq.0) then
     FACT=1.D0
  else
     LR=R-RR
     LZ=Z
     call WFSPSI(LR,LZ,PSI)
     if(PSI.lt.1.D0) then
        if(MODELP.eq.1) then
           FACT=1.D0
        elseif(MODELP.eq.2) then
           FACT=1.D0-PSI
        else
           write(6,*) 'XX WDSDEN: UNKNOWN MODELP = ',MODELP
        endif
     else
        FACT=0.D0
     endif
  end if

  ! --- set density at NODE NN ---

  do NS=1,NSMAX
     RN(NS)  =(PN(NS)-PNS(NS))*FACT+PNS(NS)
     RTPR(NS)=PTPR(NS)
     RTPP(NS)=PTPP(NS)
  enddo

  ! --- set collision frequency ---
  
  if(RN(1).gt.0.D0) then
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
                   &    /(1.24D-4*SQRT(TE*1.D-3)**3)
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
                   &    /(2.00D-1*SQRT(TE*1.D-3)**3*PA(NS))
              RNUII=0.D0
              DO NSI=2,NSMAX
                 RNUII=RNUII+PZ(NSI)**2*RN(NSI)
              ENDDO
              RNUII=RNUII*RLAMII &
                   &    /(5.31D-3*SQRT(TI*1.D-3)**3*SQRT(PA(NS)))
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
