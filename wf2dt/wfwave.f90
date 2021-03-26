!     $Id: wfwave.f90,v 1.36 2012/03/05 06:29:02 maruyama Exp $

!     ********** WF WAVE SOLVER **********

subroutine WFWAVE

  use wfcomm
  implicit none
  integer :: IERR
  real :: GTMAIN,GTSOLV,GCPUT0,GCPUT1,GCPUT2,GCPUT3

  GTMAIN=0.0
  GTSOLV=0.0
  
  call GUTIME(GCPUT0)
  
  if (nrank.eq.0) write(6,*) '--- WFWPRE start ---'
  call WFWPRE(IERR)
  if(IERR.ne.0) goto 9000

  if (nrank.eq.0) write(6,*) '--- CVCALC start ---'
  call CVCALC
  
  call GUTIME(GCPUT1)
  
  if (nrank.eq.0) write(6,*) '--- CVSOLV start ---'
  call CVSOLV
  
  call GUTIME(GCPUT2)

  if (nrank.eq.0) write(6,*) '--- CALFLD start ---'
  call CALFLD
  call PWRABS
  call PWRRAD
!  CALL TERMEP
!  CALL WFCALB
  CALL LPEFLD

  call GUTIME(GCPUT3)
  GTSOLV=GTSOLV+GCPUT2-GCPUT1
  GTMAIN=GTMAIN+GCPUT3-GCPUT2+GCPUT1-GCPUT0

!  if (nrank.eq.0) write (6,'(A/5F12.3)') &
!       "GCPUT0,GCPUT1,GCPUT2,GCPUT3,GTSOLV=", &
!        GCPUT0,GCPUT1,GCPUT2,GCPUT3,GTSOLV
!  if (nrank.eq.0) CALL LPEFLD

  if (nrank.eq.0) write(6,100) GTMAIN,GTSOLV
100 format(' ','****** CPU TIME : MAIN = ',F10.3,' SEC',5X,&
         &                     ': SOLV = ',F10.3,' SEC ******')
  
9000 continue

  return
end subroutine WFWAVE

!     ********** WF WAVE PREPARATION **********

subroutine WFWPRE(IERR)

  USE wfcomm
  USE wfparm
  USE plload,ONLY: pl_load
  USE femmeshprep
  USE feminterpolate
  USE libbes
  IMPLICIT NONE
  REAL(rkind):: RGAMMA
  INTEGER,INTENT(OUT) :: IERR

  IERR=0

  CALL fem_meshprep

  CALL fem_setup_zone

  CALL pl_load(ierr)
  if(IERR.ne.0) return

  IF(MODELG.EQ.11) THEN
     RGAMMA=ABS(Hpitch1*RRCH)
     IF(RGAMMA.LT.1.D-5) THEN
        HA1=0.D0
     ELSE
        HA1=RGAMMA*BESKNX(1,2.D0*RGAMMA)+BESKNX(2,2.D0*RGAMMA)
     ENDIF
     RKAP=SQRT((1.D0+2.D0*HA1)/(1.D0-2.D0*HA1))
     WRITE(6,'(A,1P2E12.4)') 'HA1,RKAP=',HA1,RKAP
  ENDIF
  

  SELECT CASE(MODELG)
  CASE(0,12)
     CALL WFBPSI(RA,0.D0,PSIA)
  CASE(1:10,13)
     CALL WFBPSI(RR+RA,0.D0,PSIA)
  END SELECT

  call LPELMT

  if (nrank.eq.0) write(6,*) '----- SETBDY start ---'
  call SETBDY(IERR)
  if(IERR.ne.0) return

  if (nrank.eq.0) write(6,*) '----- SETLSD start ---'
  call SETLSD
  
  if (nrank.eq.0) write(6,*) '----- MODANT start ---'
  call MODANT(IERR)
  if(IERR.ne.0) return
  
  if (nrank.eq.0) write(6,*) '----- SETEWG start ---'
  call SETEWG
  if(IERR.ne.0) return
  
  if (nrank.eq.0) write(6,*) '----- DEFMLEN start ---'
  call DEFMLEN
  
  call wffld_allocate

  if (nrank.eq.0) call WF_VIEW

  return
end subroutine WFWPRE

!     ****** DIELECTRIC TENSOR ******

subroutine DTENSR(NE,DTENS)

  use wfcomm
  implicit none
  integer,intent(in) :: NE
  integer :: IN,I,J,ID
  complex(rkind),intent(out):: DTENS(NSM,3,3,3)
  integer    :: NS,NN
  real(rkind)    :: R,Z,WW,WP(NSM),WC(NSM),BABS,AL(3),RN(NSM),RTPR(NSM)
  real(rkind)    :: RTPP(NSM),RZCL(NSM),FP,FR,FZ,DR,DZ,F
  complex(rkind) :: CWP,CWC,CDT0,CDX0,CDP0,CDT,CDP,CDX,CDAMP
  complex(rkind) :: CRR,CRP,CRZ,CPR,CPP,CPZ,CZR,CZP,CZZ

  ! ----- initialize -----  

  do J=1,3
     do I=1,3
        do IN=1,3
           do NS=1,NSMAX
              DTENS(NS,IN,I,J)=(0.d0,0.d0)
           end do
        end do
     end do
  end do

  WW=2.D0*PI*RF*1.D6

  DO NS=1,NSMAX
!     WRITE(6,'(I8,1P2E12.4)') NS,PA(NS),PZ(NS)
     WP(NS)=PZ(NS)*PZ(NS)*AEE*AEE*1.D20/(PA(NS)*AMP*EPS0*WW*WW)
     WC(NS)=PZ(NS)*AEE/(PA(NS)*AMP*WW)
  ENDDO

  ! ----- collisional cold plasma model -----

  do IN=1,3
     
     NN=NDELM(IN,NE)
     R=RNODE(NN)
     Z=ZNODE(NN)
     
     CALL WFSMAG(R,Z,BABS,AL)
     FR=AL(1)
     FP=AL(2)
     FZ=AL(3)
     
     CALL WFSDEN(R,Z,RN,RTPR,RTPP,RZCL)

     do NS=1,NSMAX
        
        CWP = WP(NS)*RN(NS)/(1.D0+CII*RZCL(NS))
        CWC = WC(NS)*BABS  /(1.D0+CII*RZCL(NS))
        CDT0= CWP/(1.D0-CWC**2)
        CDX0= CII*CWP*CWC/(1.D0-CWC**2)
        CDP0= CWP
        
        CDT=CDT0       
        CDP=CDP0-CDT0
        CDX=CDX0

        CRR= CDT   +CDP*FR*FR
        CRP= CDX*FZ+CDP*FR*FP
        CRZ=-CDX*FP+CDP*FR*FZ
        CPR=-CDX*FZ+CDP*FP*FR
        CPP= CDT   +CDP*FP*FP
        CPZ= CDX*FR+CDP*FP*FZ
        CZR= CDX*FP+CDP*FZ*FR
        CZP=-CDX*FR+CDP*FZ*FP
        CZZ= CDT   +CDP*FZ*FZ
        
        DTENS(NS,IN,1,1)=DTENS(NS,IN,1,1)-CRR
        DTENS(NS,IN,1,2)=DTENS(NS,IN,1,2)-CRP
        DTENS(NS,IN,1,3)=DTENS(NS,IN,1,3)-CRZ
        DTENS(NS,IN,2,1)=DTENS(NS,IN,2,1)-CPR
        DTENS(NS,IN,2,2)=DTENS(NS,IN,2,2)-CPP
        DTENS(NS,IN,2,3)=DTENS(NS,IN,2,3)-CPZ
        DTENS(NS,IN,3,1)=DTENS(NS,IN,3,1)-CZR
        DTENS(NS,IN,3,2)=DTENS(NS,IN,3,2)-CZP
        DTENS(NS,IN,3,3)=DTENS(NS,IN,3,3)-CZZ

     END do


     IF(WDAMP.GT.0.D0) THEN
        CDAMP=CII*PZCL(NSMAX)
        F=FDAMP
        IF(R-BDRMIN.LT.WDAMP) THEN
           ID=1
           IF(MDAMP.EQ.1.AND. &
              Z.GT.ZDAMP_MIN.AND.Z.LT.ZDAMP_MAX) ID=0
           IF(ID.EQ.1) THEN
              DR=R-BDRMIN
!              DTENS(NSMAX,IN,1,1)=DTENS(NSMAX,IN,1,1)+F*(WDAMP-DR)/(DR-CDAMP)
              DTENS(NSMAX,IN,2,2)=DTENS(NSMAX,IN,2,2)+F*(WDAMP-DR)/(DR-CDAMP)
              DTENS(NSMAX,IN,3,3)=DTENS(NSMAX,IN,3,3)+F*(WDAMP-DR)/(DR-CDAMP)
           END IF
        END IF
        IF(BDRMAX-R.LT.WDAMP) THEN
           ID=1
           IF(MDAMP.EQ.2.AND. &
              Z.GT.ZDAMP_MIN.AND.Z.LT.ZDAMP_MAX) ID=0
           IF(ID.EQ.1) THEN
              DR=BDRMAX-R
!              DTENS(NSMAX,IN,1,1)=DTENS(NSMAX,IN,1,1)+F*(WDAMP-DR)/(DR-CDAMP)
              DTENS(NSMAX,IN,2,2)=DTENS(NSMAX,IN,2,2)+F*(WDAMP-DR)/(DR-CDAMP)
              DTENS(NSMAX,IN,3,3)=DTENS(NSMAX,IN,3,3)+F*(WDAMP-DR)/(DR-CDAMP)
           END IF
        END IF
        IF(Z-BDZMIN.LT.WDAMP) THEN
           ID=1
           IF(MDAMP.EQ.3.AND. &
              R.GT.RDAMP_MIN.AND.R.LT.RDAMP_MAX) ID=0
           IF(ID.EQ.1) THEN
              DZ=Z-BDZMIN
              DTENS(NSMAX,IN,1,1)=DTENS(NSMAX,IN,1,1)+F*(WDAMP-DZ)/(DZ-CDAMP)
              DTENS(NSMAX,IN,2,2)=DTENS(NSMAX,IN,2,2)+F*(WDAMP-DZ)/(DZ-CDAMP)
!              DTENS(NSMAX,IN,3,3)=DTENS(NSMAX,IN,3,3)+F*(WDAMP-DZ)/(DZ-CDAMP)
           END IF
        END IF
        IF(BDZMAX-Z.LT.WDAMP) THEN
           ID=1
           IF(MDAMP.EQ.4.AND. &
              R.GT.RDAMP_MIN.AND.R.LT.RDAMP_MAX) ID=0
           IF(ID.EQ.1) THEN
              DZ=BDZMAX-Z
              DTENS(NSMAX,IN,1,1)=DTENS(NSMAX,IN,1,1)+F*(WDAMP-DZ)/(DZ-CDAMP)
              DTENS(NSMAX,IN,2,2)=DTENS(NSMAX,IN,2,2)+F*(WDAMP-DZ)/(DZ-CDAMP)
!              DTENS(NSMAX,IN,3,3)=DTENS(NSMAX,IN,3,3)+F*(WDAMP-DZ)/(DZ-CDAMP)
           END IF
        END IF
     END IF
  end do

  return
end subroutine DTENSR

!     ***** INTERPOLATION TENSOR *****

subroutine MUTENSR(NE,MU)

  use wfcomm
  implicit none
  integer,intent(in) :: NE
  integer :: ISD,NSD,I,J,K
  real(rkind),intent(out)::MU(3,3,6)
  real(rkind) :: A(3),B(3),C(3),L(3)

  do ISD=1,3
     NSD=ABS(NSDELM(ISD,NE))
     IF(MODELWF.EQ.0) THEN
        L(ISD)=LSID(NSD)
     ELSE
        IF(NSDELM(ISD,NE).GT.0.D0) THEN
           L(ISD)=LSID(NSD)
        ELSE
           L(ISD)=-LSID(NSD)
        END IF
     END IF
  end do

  call WFABC(NE,A,B,C)

  do K=1,6
     do J=1,3
        do I=1,3
           MU(I,J,K)=0.d0
        end do
     end do
  end do

  MU(1,1,1)= L(1)*B(2)
  MU(1,3,1)= L(1)*C(2)
  MU(1,2,4)= 1.D0
  MU(1,1,3)=-L(3)*B(3)
  MU(1,3,3)=-L(3)*C(3)

  MU(2,1,1)=-L(1)*B(1)
  MU(2,3,1)=-L(1)*C(1)
  MU(2,2,5)= 1.D0
  MU(2,1,2)= L(2)*B(3)
  MU(2,3,2)= L(2)*C(3)

  MU(3,1,2)=-L(2)*B(2)
  MU(3,3,2)=-L(2)*C(2)
  MU(3,2,6)= 1.D0
  MU(3,1,3)= L(3)*B(1)
  MU(3,3,3)= L(3)*C(1)

  return
end subroutine MUTENSR

!     ****** CURRENT COEFFICIENT VECTOR CALCULATION ******
!     LIF: Line Integral of interpolation Function
SUBROUTINE CVCALC

  use wfcomm
  implicit none
  integer    :: NE,NA,IJ,IV,I
  real(rkind)    :: RW,PHASE,MU(3,3,6),A(3),B(3),C(3)
  real(rkind)    :: R1,Z1,R2,Z2,LIF(3),R21,Z21
  complex(rkind) :: CJ(3),CVJ

  RW=2.D0*PI*RF*1.D6

  DO NE=1,NEMAX
     DO IV=1,6
        CVTOT(IV,NE)=(0.d0,0.d0)
     ENDDO
  ENDDO
  
  DO NA=1,NAMAX
     PHASE =APH(NA)*PI/180.D0
     CVJ=CII*RW*RMU0*AJ(NA)*EXP(CII*(PHASE))
     IF(JNUM(NA).EQ.1) THEN
        NE=JELMT(1,NA)
        CALL WFABC(NE,A,B,C)
        R1=RJ(1,NA)
        Z1=ZJ(1,NA)
        CJ(1)=0.D0
        CJ(2)=CVJ*RR
        CJ(3)=0.D0
        SELECT CASE(MODELG)
        CASE(0,12)
           do I=1,3
              LIF(I)=A(I)*RR &
                    +B(I)*R1*RR &
                    +C(I)*Z1*RR
           end do
        CASE(1:10,13)
           do I=1,3
              LIF(I)=A(I)*R1 &
                    +B(I)*R1*R1 &
                    +C(I)*Z1*R1
           end do
        END SELECT
        call MUTENSR(NE,MU)
        do I=1,3
           do IV=1,6
              CVTOT(IV,NE)= CVTOT(IV,NE)&
                           +LIF(I)*( MU(I,1,IV)*CJ(1)&
                                    +MU(I,2,IV)*CJ(2)&
                                    +MU(I,3,IV)*CJ(3))
           end do
        end do
     ELSE
     DO IJ=2,JNUM(NA)
        NE=JELMT(IJ,NA)
        CALL WFABC(NE,A,B,C)
        R1=RJ(IJ-1,NA)
        Z1=ZJ(IJ-1,NA)
        R2=RJ(IJ,NA)
        Z2=ZJ(IJ,NA)
        R21=R2-R1
        Z21=Z2-Z1

        CJ(1)=CVJ*R21  
        CJ(2)=(0.d0,0.d0) 
        CJ(3)=CVJ*Z21  

        SELECT CASE(MODELG)
        CASE(0,12)
           do I=1,3
              LIF(I)=A(I)*RR &
                    +B(I)*(R1+R2)*RR/2.D0 &
                    +C(I)*(Z1+Z2)*RR/2.D0
           end do
        CASE(1:10,13)
           do I=1,3
              LIF(I)=A(I)*(R1+R2)/2.d0 &
                    +B(I)*(R2**2+R1*R2+R1**2)/3.d0 &
                    +C(I)*(R2*Z1+Z2*R1+2.d0*R2*Z2+2.d0*R1*Z1)/6.d0
           end do
        END SELECT
        call MUTENSR(NE,MU)

!        WRITE(16,*) NE
!        DO I=1,3
!           DO J=1,3
!              WRITE(16,'(2I3,1P6E12.4)') I,J,MU(I,J,1:6)
!           END DO
!        END DO

        do I=1,3
           do IV=1,6
              CVTOT(IV,NE)= CVTOT(IV,NE)&
                           +LIF(I)*( MU(I,1,IV)*CJ(1)&
                                    +MU(I,2,IV)*CJ(2)&
                                    +MU(I,3,IV)*CJ(3))
!                         TEMP=LIF(I)*( MU(I,1,IV)*CJ(1)&
!                                    +MU(I,2,IV)*CJ(2)&
!                                    +MU(I,3,IV)*CJ(3))
!                         IF(ABS(TEMP).ne.0.D0) THEN
!                           WRITE(16,'(2I6,1P3E12.4)') I,IV,LIF(I),TEMP
!                           WRITE(16,'(1P3E12.4)') MU(I,1,IV),CJ(1)
!                           WRITE(16,'(1P3E12.4)') MU(I,2,IV),CJ(2)
!                           WRITE(16,'(1P3E12.4)') MU(I,3,IV),CJ(3)
!                        END IF
           end do
        end do
     end DO
     END IF
  end DO

!  do NE=1,NEMAX
!     do IV=1,6
!        if(nrank.eq.0.and.CVTOT(IV,NE).ne.(0.d0,0.d0)) &
!                                   & write(16,*) NE,IV,CVTOT(IV,NE)
!     end do
!  end do

  RETURN
END SUBROUTINE CVCALC

!     ****** LOCAL ELEMENT MATRIX  ******

SUBROUTINE CMCALC(NE)

  USE wfcomm
  IMPLICIT NONE
  INTEGER,INTENT(IN):: NE
  INTEGER:: I,J

  ! --- initialize CM ---

  DO J=1,6
     DO I=1,6
        CM(I,J)=(0.D0,0.D0)
     END DO
  END DO

  CALL CMCALCV(NE)

  CALL CMCALCS(NE)

  CALL CMCALCP(NE)

  RETURN
END SUBROUTINE CMCALC


!     ****** LOCAL ELEMENT MATRIX : volume integral ******

SUBROUTINE CMCALCV(NE)

  use wfcomm
  implicit none
  integer,intent(in) :: NE 
  integer :: I,J,K,M,N,ISD,NSD
  complex(rkind) :: CM1(3,3)
  real(rkind) :: S,L(3)
  real(rkind) :: A(3),B(3),C(3),AW(3),BW(3),CW(3)
  real(rkind) :: R(3),Z(3)

  ! --- initialize ---

  S=SELM(NE)
  call WFABC(NE,A,B,C)
  call WFNODE(NE,R,Z)

  do ISD=1,3
     NSD=ABS(NSDELM(ISD,NE))
     IF(MODELWF.EQ.0) THEN
        L(ISD)=LSID(NSD)
     ELSE
        IF(NSDELM(ISD,NE).GT.0.D0) THEN
           L(ISD)=LSID(NSD)
        ELSE
           L(ISD)=-LSID(NSD)
        END IF
     END IF
  end do

  do ISD=1,3
     M=ISD
     N=ISD+1
     if(N.gt.3) N=N-3
     AW(ISD)=L(ISD)*(A(M)*B(N)-A(N)*B(M))
     BW(ISD)=L(ISD)*(B(M)*C(N)-B(N)*C(M))
     CW(ISD)=L(ISD)*(C(M)*A(N)-C(N)*A(M))
  end do

  ! ----- rotErotF term -----

  ! --- E1F1 ---

  do J=1,3
     do I=1,3
        CM1(I,J)=(0.d0,0.d0)
     end do
  end do

  SELECT CASE(MODELG)
  CASE(0,12)
     do K=1,3 
        do J=1,3
           do I=1,3
              CM1(I,J)=CM1(I,J) &
                      +(RKZ**2)*RR &
                       *( AW(I)*AW(J)-(AW(I)*BW(J)+AW(J)*BW(I))*Z(K) &
                         +CW(I)*CW(J)-(BW(I)*CW(J)+BW(J)*CW(I))*R(K) &
                         +BW(I)*BW(J)*(R(K)**2+Z(K)**2)) &
                       *S*AIF1(K) &
                      +4.d0*BW(I)*BW(J)*RR*S*AIF1(K)
           end do
        end do
     end do
  CASE(1:10,13)
     do K=1,3 
        do J=1,3
           do I=1,3
              CM1(I,J)=CM1(I,J) &
                      +(real(NPH)**2)/R(K) &
                       *( AW(I)*AW(J)-(AW(I)*BW(J)+AW(J)*BW(I))*Z(K) &
                         +CW(I)*CW(J)-(BW(I)*CW(J)+BW(J)*CW(I))*R(K) &
                         +BW(I)*BW(J)*(R(K)**2+Z(K)**2)) &
                       *S*AIF1(K) &
                      +4.d0*BW(I)*BW(J)*R(K)*S*AIF1(K)
           end do
        end do
     end do
  END SELECT
  do J=1,3
     do I=1,3
        CM(I,J)=CM(I,J)+CM1(I,J)
     end do
  end do

  ! --- E1F2 --- 

  do J=1,3
     do I=1,3
        CM1(I,J)=(0.d0,0.d0)
     end do
  end do

  SELECT CASE(MODELG)
  CASE(0,12)
     do K=1,3
        do J=1,3
           do I=1,3
              CM1(I,J)=CM1(I,J) &
                      +(CII*RKZ*RR) &
                       *(-B(I) &
                          *(AW(J)-BW(J)*Z(K)) &
                         +C(I) &
                          *(CW(J)-BW(J)*R(K))) &
                       *S*AIF1(K)
           end do
        end do
     end do
  CASE(1:10,13)
     do K=1,3
        do J=1,3
           do I=1,3
              CM1(I,J)=CM1(I,J) &
!
                       +(CII*real(NPH)) &
                       *(-B(I) &
                          *(AW(J)-BW(J)*Z(K)) &
                         +C(I) &
                          *(CW(J)-BW(J)*R(K))) &
                       *S*AIF1(K)&
!
                      +(CII*real(NPH))&
                       *(-(AW(J)-BW(J)*Z(K))/R(K))&
                       *S*AIF2(I,K)
           end do
        end do
     end do
  END SELECT

  do J=1,3
     do I=1,3
        CM(I+3,J)=CM(I+3,J)+CM1(I,J)
     end do
  end do

  ! --- E2F1 ---

  do J=1,3
     do I=1,3
        CM1(I,J)=(0.d0,0.d0)
     end do
  end do

  SELECT CASE(MODELG)
  CASE(0,12)
     do K=1,3
        do J=1,3
           do I=1,3
              CM1(I,J)=CM1(I,J) &
                      -(CII*RKZ*RR) &
                       *(-B(J) &
                          *(AW(I)-BW(I)*Z(K)) &
                         +C(J)&
                          *(CW(I)-BW(I)*R(K))) &
                        *S*AIF1(K)
           end do
        end do
     end do
  CASE(1:10,13)
     do K=1,3
        do J=1,3
           do I=1,3
              CM1(I,J)=CM1(I,J) &
!
                      -(CII*real(NPH)) &
                       *(-B(J) &
                          *(AW(I)-BW(I)*Z(K)) &
                         +C(J)&
                          *(CW(I)-BW(I)*R(K))) &
                        *S*AIF1(K) &
!
                      -(CII*real(NPH)) &
                       *(-(AW(I)-BW(I)*Z(K))/R(K)) &
                        *S*AIF2(J,K)
           end do
        end do
     end do
  END SELECT
  do J=1,3
     do I=1,3
        CM(I,J+3)=CM(I,J+3)+CM1(I,J)
     end do
  end do

  ! --- E2F2 ---

  do J=1,3
     do I=1,3
        CM1(I,J)=(0.d0,0.d0)
     end do
  end do

  SELECT CASE(MODELG)
  CASE(0,12)
     do K=1,3
        do J=1,3
           do I=1,3
              CM1(I,J)=CM1(I,J) &
                      +(B(I)*B(J)+C(I)*C(J))*RR*S*AIF1(K)
!                      +B(J)*S*AIF2(I,K) &
!                      +B(I)*S*AIF2(J,K) &
!                      +1.D0/RR*S*AIF3(I,J,K)
           end do
        end do
     end do
  CASE(1:10,13)
     do K=1,3
        do J=1,3
           do I=1,3
              CM1(I,J)=CM1(I,J) &
                      +(B(I)*B(J)+C(I)*C(J))*R(K)*S*AIF1(K) &
                      +B(J)*S*AIF2(I,K) &
                      +B(I)*S*AIF2(J,K) &
                      +1.D0/(R(K))*S*AIF3(I,J,K)
           end do
        end do
     end do
  END SELECT
  do J=1,3
     do I=1,3
        CM(I+3,J+3)=CM(I+3,J+3)+CM1(I,J)
     end do
  end do

  RETURN
END SUBROUTINE CMCALCV

!     ****** LOCAL ELEMENT MATRIX : Surface integral ******

SUBROUTINE CMCALCS(NE)

  use wfcomm
  implicit none
  integer,intent(in) :: NE 
  integer :: I,J,K,I1,J1,K1,ID,ISD,NSD,IND(2),ND1,ND2
  real(rkind) :: S,L(3),RN,ZN
  real(rkind) :: A(3),B(3),C(3),AW(3),BW(3),CW(3)
  real(rkind) :: R(3),Z(3)
  COMPLEX(rkind):: CTEMP

  ! --- if no boundary side, return ---

  ID=0
  DO ISD=1,3
     NSD=ABS(NSDELM(ISD,NE))
     IF(KASID(NSD).EQ.1) ID=1
  END DO
  IF(ID.EQ.0) RETURN

  ! --- initialize ---

  S=SELM(NE)
  call WFABC(NE,A,B,C)
  call WFNODE(NE,R,Z)

  DO ISD=1,3
     NSD=ABS(NSDELM(ISD,NE))
     IF(KASID(NSD).EQ.1) THEN
        IF(MODELWF.EQ.0) THEN
           L(ISD)=LSID(NSD)
        ELSE
           IF(NSDELM(ISD,NE).GT.0.D0) THEN
              L(ISD)=LSID(NSD)
           ELSE
              L(ISD)=-LSID(NSD)
           END IF
        END IF

        ND1=ISD
        ND2=ISD+1
        IF(ND2.GT.3) ND2=ND2-3
        AW(ISD)=L(ISD)*(A(ND1)*B(ND2)-A(ND2)*B(ND1))
        BW(ISD)=L(ISD)*(B(ND1)*C(ND2)-B(ND2)*C(ND1))
        CW(ISD)=L(ISD)*(C(ND1)*A(ND2)-C(ND2)*A(ND1))
     END IF
  END DO

  DO ISD=1,3
     NSD=ABS(NSDELM(ISD,NE))
     IF(KASID(NSD).EQ.1) THEN

        IND(1)=ISD
        IND(2)=ISD+1
        IF(IND(2).GT.3) IND(2)=IND(2)-3
        RN= (Z(IND(2))-Z(IND(1)))/L(ISD)
        ZN=-(R(IND(2))-R(IND(1)))/L(ISD)

  ! ----- rotErotF term -----

  ! --- E1F1 ---

        I=ISD
        J=ISD
        SELECT CASE(MODELG)
        CASE(0,12)
           DO K1=1,2
              K=IND(K1)
              CTEMP=CM(I,J)
              CM(I,J)=CM(I,J) &
                      -2.D0*BW(J)*(ZN*(AW(I)-BW(I)*Z(K)) &
                                  +RN*(CW(I)-BW(I)*R(K)))*RR*L(ISD)*AIE1(K)
!                    WRITE(21,'(A,I10,3I5,1P4E12.4)') &
!                         'CM:',NE,I,J,K,CM(I,J),CM(I,J)-CTEMP
!                    WRITE(21,'(5X,1P5E12.4)') &
!                         BW(J),ZN,AW(I),BW(I),Z(K)
!                    WRITE(21,'(5X,1P5E12.4)') &
!                         RN,CW(I),R(K),RR,L(ISD),AIE1(K)
           END DO
        CASE(1:10,13)
           DO K1=1,2
              K=IND(K1)
              CM(I,J)=CM(I,J) &
                      -2.D0*BW(J)*(ZN*(AW(I)-BW(I)*Z(K)) &
                                  +RN*(CW(I)-BW(I)*R(K)))*R(K)*L(ISD)*AIE1(K)
           END DO
        END SELECT

  ! --- E2F1 --- 

        I=ISD
        DO J1=1,2
           J=IND(J1)
           SELECT CASE(MODELG)
           CASE(0,12)
              DO K1=1,2
                 K=IND(K1)
                 CTEMP=CM(I,J+3)
                 CM(I,J+3)=CM(I,J+3) &
                          +CII*RKZ*RR &
                           *(-ZN*(-CW(J)+BW(J)*R(K)) &
                             -RN*( AW(J)-BW(J)*Z(K)))*L(ISD)*AIE2(J,K)
                    WRITE(21,'(A,I10,3I5,1P4E12.4)') &
                         'CM:',NE,I,J+3,K,CM(I,J+3),CM(I,J+3)-CTEMP
              END DO
           CASE(1:10,13)
              DO K1=1,2
                 K=IND(K1)
                 CM(I,J+3)=CM(I,J+3) &
                          +CII*NPH &
                           *(-ZN*(-CW(J)+BW(J)*R(K)) &
                             -RN*( AW(J)-BW(J)*Z(K)))*L(ISD)*AIE2(J,K)
              END DO
           END SELECT
        END DO

  ! --- E2F2 ---

        DO J1=1,2
           J=IND(J1)
           DO I1=1,2
              I=IND(I1)
              SELECT CASE(MODELG)
              CASE(0,12)
                 DO K1=1,2
                    K=IND(K1)
                    CTEMP=CM(I+3,J+3)
                    CM(I+3,J+3)=CM(I+3,J+3) &
                               +(RN*B(J)+ZN*C(J))*RR &
                                *L(ISD)*AIE3(I,J,K)
                    WRITE(21,'(A,I10,3I5,1P4E12.4)') &
                         'CM:',NE,I+3,J+3,K,CM(I+3,J+3),CM(I+3,J+3)-CTEMP
                 END DO
              CASE(1:10,13)
                 DO K1=1,2
                    K=IND(K1)
                    CM(I+3,J+3)=CM(I+3,J+3) &
                               +((RN*B(J)+ZN*C(J))*R(K) &
                                 +RN*(A(J)+B(J)*R(K)+C(J)*Z(K))) &
                                *L(ISD)*AIE3(I,J,K)
                 END DO
              END SELECT
           END DO
        END DO
     END IF
  END DO

  RETURN
END SUBROUTINE CMCALCS

!     ****** LOCAL ELEMENT MATRIX : Surface integral ******

SUBROUTINE CMCALCP(NE)

  use wfcomm
  implicit none
  integer,intent(in) :: NE 
  integer :: I,J,K,NS,IN,JJ,II
  real(rkind) :: RW,WC,WC2
  real(rkind) :: S
  real(rkind) :: A(3),B(3),C(3)
  real(rkind) :: R(3),Z(3),MU(3,3,6)
  COMPLEX(rkind):: CM2(6,6)
  complex(rkind) :: DTENS(NSM,3,3,3)
  complex(rkind) :: DTENST(3,3,3)

  ! --- initialize ---

  RW=2.D0*PI*RF*1.D6
  WC=RW/VC
  WC2=WC**2

  S=SELM(NE)  
  call WFABC(NE,A,B,C)
  call WFNODE(NE,R,Z)

  ! ----- dielectric tensor term -----

  call DTENSR(NE,DTENS)
  call MUTENSR(NE,MU)

  DO J=1,6
     DO I=1,6
        CM2(I,J)=(0.D0,0.D0)
     END DO
  END DO

  do J=1,3
     do I=1,3
        do IN=1,3
           if(I.eq.J) then
              DTENST(IN,I,J)=(1.d0,0.d0)
           else
              DTENST(IN,I,J)=(0.d0,0.d0)
           end if
        end do
     end do
  end do

  ! --- assemble dielectric tensor ---

  do J=1,3
     do I=1,3
        do IN=1,3
           do NS=1,NSMAX
              DTENST(IN,I,J)=DTENST(IN,I,J)+DTENS(NS,IN,I,J)
           end do
        end do
     end do
  end do

  SELECT CASE(MODELG)
  CASE(0,12)
     do JJ=1,6
        do II=1,6
           do K=1,3
              do J=1,3
                 do I=1,3
                    CM2(II,JJ)= CM2(II,JJ)&
                                +((MU(I,1,II)*DTENST(J,1,1)&
                                  +MU(I,2,II)*DTENST(J,2,1)&
                                  +MU(I,3,II)*DTENST(J,3,1))&
                                  *MU(K,1,JJ)&
                                 +(MU(I,1,II)*DTENST(J,1,2)&
                                  +MU(I,2,II)*DTENST(J,2,2)&
                                  +MU(I,3,II)*DTENST(J,3,2))&
                                  *MU(K,2,JJ)&
                                 +(MU(I,1,II)*DTENST(J,1,3)&
                                  +MU(I,2,II)*DTENST(J,2,3)&
                                  +MU(I,3,II)*DTENST(J,3,3))&
                                  *MU(K,3,JJ))&
                                 *RR*S*AIF3(I,J,K)
                 end do
              end do
           end do
        end do
     end do
  CASE(1:10,13)
     do JJ=1,6
        do II=1,6
           do K=1,3
              do J=1,3
                 do I=1,3
                    CM2(II,JJ)= CM2(II,JJ)&
                                +((MU(I,1,II)*DTENST(J,1,1)&
                                  +MU(I,2,II)*DTENST(J,2,1)&
                                  +MU(I,3,II)*DTENST(J,3,1))&
                                  *MU(K,1,JJ)&
                                 +(MU(I,1,II)*DTENST(J,1,2)&
                                  +MU(I,2,II)*DTENST(J,2,2)&
                                  +MU(I,3,II)*DTENST(J,3,2))&
                                  *MU(K,2,JJ)&
                                 +(MU(I,1,II)*DTENST(J,1,3)&
                                  +MU(I,2,II)*DTENST(J,2,3)&
                                  +MU(I,3,II)*DTENST(J,3,3))&
                                  *MU(K,3,JJ))&
                                 *R(J)*S*AIF3(I,J,K)
                 end do
              end do
           end do
        end do
     end do
  END SELECT

  do J=1,6
     do I=1,6
        CM(I,J)=CM(I,J)-WC2*CM2(I,J)
     end do
  end do

!  write(*,*) 'NE=',NE
!  do I=1,6
!     write(6,'(6(A1,E9.3,A1,E9.3,A1))')&
!          "(",real(CM(I,1)),',',aimag(CM(I,1)),")",&
!          "(",real(CM(I,2)),',',aimag(CM(I,2)),")",&
!          "(",real(CM(I,3)),',',aimag(CM(I,3)),")",&
!          "(",real(CM(I,4)),',',aimag(CM(I,4)),")",&
!          "(",real(CM(I,5)),',',aimag(CM(I,5)),")",&
!          "(",real(CM(I,6)),',',aimag(CM(I,6)),")"
!  end do

  return
END SUBROUTINE CMCALCP

!     ******* ELECTRIC FIELD CALCULATION *******

SUBROUTINE CALFLD

  use wfcomm
  implicit none
  integer :: NN,NSD,NV

  DO NSD=1,NSDMAX
     CESD(NSD)=(0.d0,0.d0)
  ENDDO
  DO NN=1,NNMAX
     CEND(NN) =(0.d0,0.d0)
  END DO

  DO NSD=1,NSDMAX
     NV=NVNSD(NSD)
     if (NV.eq.0) then
        IF(KBSID(NSD).NE.0) THEN
           CESD(NSD)=CEBSD(KBSID(NSD))
        ELSE
           CESD(NSD)=(0.d0,0.d0)
        END IF
     else
        CESD(NSD)=CSV(NV)
     end if
!     if(nrank.eq.0) write(6,*) NSD,CESD(NSD),KASID(NSD)
  END DO
  DO NN=1,NNMAX
     NV=NVNN(NN)
     if (NV.eq.0) then
        IF(KBNOD(NN).NE.0) THEN
           CEND(NN)=CEBND(KBNOD(NN))
        ELSE
           CEND(NN)=(0.d0,0.d0)
        END IF
     else
        CEND(NN)=CSV(NV)
     end if
!     if(nrank.eq.0) write(6,*) NN,CEND(NN),KANOD(NN)
  END DO
 
  RETURN
END SUBROUTINE CALFLD

!     ******* POWER ABSORPTION *******

SUBROUTINE PWRABS

  use wfcomm
  USE libmpi
  implicit none

  integer    :: NE,IN,NN,NSD,NS
  integer    :: I,J,K,II,JJ
  real(rkind),dimension(:,:),ALLOCATABLE:: PABS
  real(rkind)    :: RW,S,MU(3,3,6),R(3),Z(3)
  complex(rkind) :: DTENS(NSM,3,3,3),CTENS(NSM,3,3,3)
  complex(rkind) :: CIWE,CINT(NSM,6,6),CE(6)
  INTEGER,ALLOCATABLE:: nelm_len_nrank(:),nelm_pos_nrank(:)
  REAL(rkind),ALLOCATABLE:: rdata(:),rdata_tot(:)
  INTEGER:: ipos,n,nblk,ndata,nelm1,nelm2,nsize_high,nsize_low

  ALLOCATE(nelm_len_nrank(0:nsize-1),nelm_pos_nrank(0:nsize-1))
  nblk=nemax/nsize
  nsize_high=nemax-nblk*nsize
  nsize_low=nsize-nsize_high
  ipos=0
  DO n=0,nsize_high-1
     nelm_len_nrank(n)=nblk+1
     nelm_pos_nrank(n)=ipos
     ipos=ipos+nelm_len_nrank(n)
  END DO
  DO n=nsize_high,nsize-1
     nelm_len_nrank(n)=nblk
     nelm_pos_nrank(n)=ipos
     ipos=ipos+nelm_len_nrank(n)
  END DO
  IF(ipos.NE.nemax) THEN
     WRITE(6,'(A,2I8)') 'XX ne parallel error: nemax,ipos=',nemax,ipos
     STOP
  END IF

  ! --- initialize ---
  
  allocate(PABS(NSMAX,NEMAX))
  
  RW=2.D0*PI*RF*1.D6
  CIWE=CII*RW*EPS0

  do NE=nelm_pos_nrank(nrank)+1,nelm_pos_nrank(nrank)+nelm_len_nrank(nrank)
     S=SELM(NE)

     ! --- calculate conductivity tensor ---

     call DTENSR(NE,DTENS)
     do NS=1,NSMAX
        do IN=1,3
           do J=1,3
              do I=1,3
                 CTENS(NS,IN,I,J)=-CIWE*DTENS(NS,IN,I,J)
              end do
           end do
        end do
     end do

     call WFNODE(NE,R,Z)
     call MUTENSR(NE,MU)
     
     CINT=0.d0

     do NS=1,NSMAX
        SELECT CASE(MODELG)
        CASE(0,12)
           do JJ=1,6
              do II=1,6
                 do K=1,3
                    do J=1,3
                       do I=1,3
                          CINT(NS,II,JJ)= CINT(NS,II,JJ) &
                                         +((MU(I,1,II)*CTENS(NS,J,1,1) &
                                           +MU(I,2,II)*CTENS(NS,J,2,1) &
                                           +MU(I,3,II)*CTENS(NS,J,3,1)) &
                                           *MU(K,1,JJ)&
                                          +(MU(I,1,II)*CTENS(NS,J,1,2)&
                                           +MU(I,2,II)*CTENS(NS,J,2,2)&
                                           +MU(I,3,II)*CTENS(NS,J,3,2))&
                                           *MU(K,2,JJ)&
                                          +(MU(I,1,II)*CTENS(NS,J,1,3)&
                                           +MU(I,2,II)*CTENS(NS,J,2,3)&
                                           +MU(I,3,II)*CTENS(NS,J,3,3))&
                                           *MU(K,3,JJ))&
                                          *RR*S*AIF3(I,J,K)
                       end do
                    end do
                 end do
              end do
           end do
        CASE(1:10,13)
           do JJ=1,6
              do II=1,6
                 do K=1,3
                    do J=1,3
                       do I=1,3
                          CINT(NS,II,JJ)= CINT(NS,II,JJ) &
                                         +((MU(I,1,II)*CTENS(NS,J,1,1) &
                                           +MU(I,2,II)*CTENS(NS,J,2,1) &
                                           +MU(I,3,II)*CTENS(NS,J,3,1)) &
                                           *MU(K,1,JJ)&
                                          +(MU(I,1,II)*CTENS(NS,J,1,2)&
                                           +MU(I,2,II)*CTENS(NS,J,2,2)&
                                           +MU(I,3,II)*CTENS(NS,J,3,2))&
                                           *MU(K,2,JJ)&
                                          +(MU(I,1,II)*CTENS(NS,J,1,3)&
                                           +MU(I,2,II)*CTENS(NS,J,2,3)&
                                           +MU(I,3,II)*CTENS(NS,J,3,3))&
                                           *MU(K,3,JJ))&
                                          *R(J)*S*AIF3(I,J,K)
                       end do
                    end do
                 end do
              end do
           end do
        END SELECT
     end do

     do I=1,3
        NSD=NSDELM(I,NE)
        if(NSD.lt.0) then
           NSD=-NSD
           CE(I)=-CESD(NSD)
        else
           CE(I)=CESD(NSD)
        end if
     end do
     do I=1,3
        NN=NDELM(I,NE)
        CE(I+3)=CEND(NN)
     end do

     do NS=1,NSMAX
        PABS(NS,NE)=0.d0
        do JJ=1,6
           do II=1,6
              PABS(NS,NE)=PABS(NS,NE)&
                            +0.5d0*real(CONJG(CE(II))*CINT(NS,II,JJ)*CE(JJ))
           end do
        end do
     end do

  end do

  nelm1=nelm_pos_nrank(nrank)+1
  nelm2=nelm_pos_nrank(nrank)+nelm_len_nrank(nrank)
  ndata=nelm_len_nrank(nrank)
  ALLOCATE(rdata(ndata),rdata_tot(nemax))
  DO ns=1,nsmax
     rdata(1:ndata)=pabs(ns,nelm1:nelm2)
     CALL mtx_allgatherv_real8(rdata,ndata,rdata_tot,nemax, &
          nelm_len_nrank,nelm_pos_nrank)
     pabs(ns,1:nemax)=rdata_tot(1:nemax)
  END DO

  do NS=1,NSMAX
     PABST(NS)=0.d0
     do NE=1,NEMAX
        PABST(NS)=PABST(NS)+PABS(NS,NE)
     end do
  end do

  PABSTT=0.D0
  DO NS=1,NSMAX
     PABSTT=PABSTT+PABST(NS)
  END DO

  deallocate(PABS)

  RETURN
END SUBROUTINE PWRABS

!     ******* POWER RADIATION *******

SUBROUTINE PWRRAD

  use wfcomm
  implicit none

  integer    :: NE,NA,I,NN,IV
  integer    :: IJ,NSD
  real(rkind)    :: PHASE,RW,LIF(3),A(3),B(3),C(3)
  real(rkind)    :: R1,R2,Z1,Z2,R21,Z21,MU(3,3,6)
  complex(rkind) :: CE(6),CJ(3),CVJ

  ! --- initialize ---

  RW=2.D0*PI*RF*1.D6

  do NA=1,NAMAX
     PHASE =APH(NA)*PI/180.D0
     CVJ=AJ(NA)*EXP(CII*(PHASE))
     CIMP(NA)=(0.d0,0.d0)
     do IJ=2,JNUM(NA)
        NE=JELMT(IJ,NA)
        R1=RJ(IJ-1,NA)
        Z1=ZJ(IJ-1,NA)
        R2=RJ(IJ,NA)
        Z2=ZJ(IJ,NA)
        R21=R2-R1
        Z21=Z2-Z1
        CJ(1)=CVJ*R21
        CJ(2)=(0.d0,0.d0)
        CJ(3)=CVJ*Z21

        call WFABC(NE,A,B,C)

        SELECT CASE(MODELG)
        CASE(0,12)
           do I=1,3
              LIF(I)= A(I)*RR &
                     +B(I)*RR*(R2+R1)/2.d0 &
                     +C(I)*RR*(Z1+Z2)/2.d0           
           end do
        CASE(1:10,13)
           do I=1,3
              LIF(I)= A(I)*(R1+R2)/2.d0 &
                     +B(I)*(R2**2+R1*R2+R1**2)/3.d0 &
                     +C(I)*(R2*Z1+Z2*R1+2.d0*R2*Z2+2.d0*R1*Z1)/6.d0           
           end do
        END SELECT
        do I=1,3
        NSD=NSDELM(I,NE)
        if(NSD.lt.0) then
           NSD=-NSD
           CE(I)=-CESD(NSD)
           else
           CE(I)=CESD(NSD)
           end if
        end do
        do I=1,3
           NN=NDELM(I,NE)
           CE(I+3)=CEND(NN)
        end do

        call MUTENSR(NE,MU)

        do I=1,3
           do IV=1,6
              CIMP(NA)=CIMP(NA)&
                         -0.5d0*LIF(I)*CONJG(CE(IV))&
                                      *( MU(I,1,IV)*CJ(1)&
                                        +MU(I,2,IV)*CJ(2)&
                                        +MU(I,3,IV)*CJ(3))
           end do
        end do

     end do
  end do

  CTIMP=(0.d0,0.d0)

  do NA=1,NAMAX
     CTIMP=CTIMP+CIMP(NA)
  end do

  RETURN
END SUBROUTINE PWRRAD

!     ******* OUTPUT FIELD DATA *******

  SUBROUTINE LPEFLD

    USE wfcomm
    IMPLICIT NONE
    INTEGER:: I,NR1,NR2,NR3,NR4,NS,NR,NA,NZ,J
    REAL(rkind):: P
    REAL(rkind):: P1(NSM),P2(NSM),P3(NSM),P4(NSM)

    IF(NPRINT.LT.1) RETURN
    IF(nrank.NE.0) RETURN
    
!    WRITE(6,110) (EMAX(I),I=1,3),ETMAX,PNMAX
!110 FORMAT(1H ,'EXMAX  =',1PE12.4 &
!         ,3X ,'EYMAX  =',1PE12.4 &
!         ,3X ,'EZMAX  =',1PE12.4/ &
!         1H ,'EMAX   =',1PE12.4 &
!         ,3X ,'PNMAX  =',1PE12.4)

    WRITE(6,120) DBLE(CTIMP),PABSTT
120 FORMAT(1H ,'RADIATED POWER =',1PE12.4/ &
         1H ,'ABSORBED POWER =',1PE12.4)

    DO NS=1,NSMAX
       WRITE(6,126) NS,PABST(NS)
126    FORMAT(1H ,'      PABS(',I2,') =',1PE12.4)
    END DO

    WRITE(6,130)
130 FORMAT(1H ,' I JNUM', '  AJ(I)','  APH(I)','  AWD(I)', &
            ' APOS(I)',' XJ(I)','  YJ(I)', &
            8X,'LOADING IMP.[ohm]')
    DO NA=1,NAMAX
       WRITE(6,140) NA,JNUM(NA),AJ(NA),APH(NA),AWD(NA),APOS(NA), &
                               RJ(1,NA),ZJ(1,NA),CIMP(NA)
140    FORMAT(1H ,I2,I3,0PF8.2,F8.2,1X,4F7.4,2X,'(',1P2E12.4,')')
    END DO

    IF(NPRINT.LT.2) RETURN

    ! field output

    RETURN
  END SUBROUTINE LPEFLD
  
!     ******* OUTPUT ELEMENT DATA *******

SUBROUTINE LPELMT

  use wfcomm
  implicit none

  integer :: I,J,NA

  IF(NPRINT.LT.3) RETURN
  IF(nrank.NE.0) RETURN
     
  WRITE(6,110) NNMAX
110 FORMAT(/' ','NODE DATA     : #### NNMAX =',I5,' ####'/&
         &       ' ',2('  NNMAX',' KANOD',&
         &       9X,'R',14X,'Z',9X))
  WRITE(6,115) (I,KANOD(I),RNODE(I),ZNODE(I),&
       &              I=1,NNMAX)
115 FORMAT((' ',2(2I6,2X,1P2E15.7,2X)))
  
  WRITE(6,120) NEMAX,(I,(NDELM(J,I),J=1,3),I=1,NEMAX)
120 FORMAT(/' ','ELEMENT DATA  : #### NEMAX =',I5,' ####'/&
         &      (' ',4(I6,'(',3I5,')',2X)))
  
  WRITE(6,125) NEMAX,(I,(NSDELM(J,I),J=1,3),I=1,NEMAX)
125 FORMAT(/' ','SIDE    DATA  : #### NEMAX =',I5,' ####'/&
         &      (' ',2(I8,'(',3I8,')',2X)))
  
  DO NA=1,NAMAX
     WRITE(6,140) NA,JNUM0(NA)
140  FORMAT(/' ','ORIGINAL ANTENNA DATA : NA =',I5,' JNUM0 =',I5/&
          &          ' ',2('  NO.',13X,' RJ0',11X,' ZJ0',6X))
     WRITE(6,150) (I,RJ0(I,NA),ZJ0(I,NA),I=1,JNUM0(NA))
150  FORMAT((' ',2(I5,8X,1P2E15.7)))
     
     WRITE(6,154) NA,JNUM(NA)
154  FORMAT(/' ','MODIFIED ANTENNA DATA : NA =',I5,' JNUM  =',I5/&
          &          ' ',2('  NO.',' JELM',8X,' JR ',11X,' JZ ',6X))
     WRITE(6,156) (I,JELMT(I,NA),RJ(I,NA),ZJ(I,NA),I=1,JNUM(NA))
156  FORMAT((' ',2(2I5,3X,1P2E15.7)))
  ENDDO
  
  RETURN
END SUBROUTINE LPELMT
