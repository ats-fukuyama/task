! wfwave.f90

MODULE wfwave

  PRIVATE
  PUBLIC wf_wave
  PUBLIC wf_wpre
  PUBLIC wf_calfld
  PUBLIC wf_pwrabs
  PUBLIC wf_pwrrad
  PUBLIC wf_lpefld
  PUBLIC wf_lpelmt

CONTAINS

!     ********** WF WAVE SOLVER **********

  SUBROUTINE wf_wave

    use wfcomm
    USE wfsolv
  implicit none
  integer :: IERR
  real :: GTMAIN,GTSOLV,GCPUT0,GCPUT1,GCPUT2,GCPUT3

  GTMAIN=0.0
  GTSOLV=0.0
  
  call GUTIME(GCPUT0)

  if (nrank.eq.0) write(6,*) '--- wf_wpre start ---'
  call wf_wpre(IERR)
  if(IERR.ne.0) goto 9000

  CALL wf_field_allocate
  
  if (nrank.eq.0) write(6,*) '--- wf_cvcalc start ---'
  call wf_cvcalc
  
  call GUTIME(GCPUT1)
  
  if (nrank.eq.0) write(6,*) '--- wf_cvsolv start ---'
  call wf_cvsolv
  
  call GUTIME(GCPUT2)

  if (nrank.eq.0) write(6,*) '--- wf_calfld start ---'
  call wf_calfld
  call wf_pwrabs
  call wf_pwrrad
!  CALL TERMEP
!  CALL WFCALB
  CALL wf_lpefld

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
  end subroutine wf_wave

!     ********** WF WAVE PREPARATION **********

  subroutine wf_wpre(IERR)

  USE wfcomm
  USE wfparm
  USE wfprof,ONLY: wf_bpsi
  USE plload,ONLY: pl_load
  USE wfant
  USE wfsolv
  USE wfindex
  USE wfsub
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
        Hpitch2=0.D0
     ELSE
        Hpitch2=RGAMMA*BESKNX(1,2.D0*RGAMMA)+BESKNX(2,2.D0*RGAMMA)
     ENDIF
     RKAP=SQRT((1.D0+2.D0*Hpitch2)/(1.D0-2.D0*Hpitch2))
     WRITE(6,'(A,1P2E12.4)') 'Hpitch2,RKAP=',Hpitch2,RKAP
  ENDIF
  

  SELECT CASE(MODELG)
  CASE(0,12)
     CALL wf_bpsi(RA,0.D0,PSIA)
  CASE(1:10,13)
     CALL wf_bpsi(RR+RA,0.D0,PSIA)
  END SELECT

  call wf_lpelmt

  if (nrank.eq.0) write(6,*) '----- wf_set_bdy start ---'
  call wf_set_bdy(IERR)
  if(IERR.ne.0) return

  if (nrank.eq.0) write(6,*) '----- wf_set_lside start ---'
  call wf_set_lside
  
  if (nrank.eq.0) write(6,*) '----- wf_modant start ---'
  call wf_modant(IERR)
  if(IERR.ne.0) return
  
  if (nrank.eq.0) write(6,*) '----- wf_set_ewg start ---'
  call wf_set_ewg
  if(IERR.ne.0) return
  
  if (nrank.eq.0) write(6,*) '----- wf_defmlen start ---'
  call wf_defmlen
  
  if (nrank.eq.0) call WF_VIEW

  return
  end subroutine wf_wpre

!     ****** CURRENT COEFFICIENT VECTOR CALCULATION ******
!     LIF: Line Integral of interpolation Function

  SUBROUTINE wf_cvcalc

    use wfcomm
    USE wfsolv
    USE wfsub
  implicit none
  integer    :: NE,NA,IJ,IV,I
  real(rkind)    :: RW,PHASE,MU(3,3,6),A(3),B(3),C(3)
  real(rkind)    :: R1,Z1,R2,Z2,LIF(3),R21,Z21
  complex(rkind) :: CJ(3),CVJ

  RW=2.D0*PI*RF*1.D6

  DO NE=1,nelm_max
     DO IV=1,6
        CVTOT(IV,NE)=(0.d0,0.d0)
     ENDDO
  ENDDO
  
  DO NA=1,nant_max
     PHASE =APH(NA)*PI/180.D0
     CVJ=CII*RW*RMU0*AJ(NA)*EXP(CII*(PHASE))
     IF(JNUM(NA).EQ.1) THEN
        NE=JELMT(1,NA)
        CALL wf_set_abc(NE,A,B,C)
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
        call wf_mutensr(NE,MU)
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
        CALL wf_set_abc(NE,A,B,C)
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
        call wf_mutensr(NE,MU)

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

!  do NE=1,nelm_max
!     do IV=1,6
!        if(nrank.eq.0.and.CVTOT(IV,NE).ne.(0.d0,0.d0)) &
!                                   & write(16,*) NE,IV,CVTOT(IV,NE)
!     end do
!  end do

  RETURN
  END SUBROUTINE wf_cvcalc



!     ******* ELECTRIC FIELD CALCULATION *******

  SUBROUTINE wf_calfld

  use wfcomm
  implicit none
  integer :: NN,NSD,NV

  DO NSD=1,nseg_max
     CESD(NSD)=(0.d0,0.d0)
  ENDDO
  DO NN=1,node_max
     CEND(NN) =(0.d0,0.d0)
  END DO

  DO NSD=1,nseg_max
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
  DO NN=1,node_max
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
END SUBROUTINE wf_calfld

!     ******* POWER ABSORPTION *******

  SUBROUTINE wf_pwrabs

    use wfcomm
    USE wfsolv
    USE wfsub
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
  nblk=nelm_max/nsize
  nsize_high=nelm_max-nblk*nsize
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
  IF(ipos.NE.nelm_max) THEN
     WRITE(6,'(A,2I8)') 'XX ne parallel error: nelm_max,ipos=',nelm_max,ipos
     STOP
  END IF

  ! --- initialize ---
  
  allocate(PABS(NSMAX,nelm_max))
  
  RW=2.D0*PI*RF*1.D6
  CIWE=CII*RW*EPS0

  do NE=nelm_pos_nrank(nrank)+1,nelm_pos_nrank(nrank)+nelm_len_nrank(nrank)
     S=SELM(NE)

     ! --- calculate conductivity tensor ---

     call wf_dtensr(NE,DTENS)
     do NS=1,NSMAX
        do IN=1,3
           do J=1,3
              do I=1,3
                 CTENS(NS,IN,I,J)=-CIWE*DTENS(NS,IN,I,J)
              end do
           end do
        end do
     end do

     call wf_set_node(NE,R,Z)
     call wf_mutensr(NE,MU)
     
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
        NSD=nseg_nside_nelm(I,NE)
        if(NSD.lt.0) then
           NSD=-NSD
           CE(I)=-CESD(NSD)
        else
           CE(I)=CESD(NSD)
        end if
     end do
     do I=1,3
        NN=node_nside_nelm(I,NE)
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
  ALLOCATE(rdata(ndata),rdata_tot(nelm_max))
  DO ns=1,nsmax
     rdata(1:ndata)=pabs(ns,nelm1:nelm2)
     CALL mtx_allgatherv_real8(rdata,ndata,rdata_tot,nelm_max, &
          nelm_len_nrank,nelm_pos_nrank)
     pabs(ns,1:nelm_max)=rdata_tot(1:nelm_max)
  END DO

  do NS=1,NSMAX
     PABST(NS)=0.d0
     do NE=1,nelm_max
        PABST(NS)=PABST(NS)+PABS(NS,NE)
     end do
  end do

  PABSTT=0.D0
  DO NS=1,NSMAX
     PABSTT=PABSTT+PABST(NS)
  END DO

  deallocate(PABS)

  RETURN
END SUBROUTINE wf_pwrabs

!     ******* POWER RADIATION *******

  SUBROUTINE wf_pwrrad

    use wfcomm
    USE wfsolv
    USE wfsub
  implicit none

  integer    :: NE,NA,I,NN,IV
  integer    :: IJ,NSD
  real(rkind)    :: PHASE,RW,LIF(3),A(3),B(3),C(3)
  real(rkind)    :: R1,R2,Z1,Z2,R21,Z21,MU(3,3,6)
  complex(rkind) :: CE(6),CJ(3),CVJ

  ! --- initialize ---

  RW=2.D0*PI*RF*1.D6

  do NA=1,nant_max
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

        call wf_set_abc(NE,A,B,C)

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
        NSD=nseg_nside_nelm(I,NE)
        if(NSD.lt.0) then
           NSD=-NSD
           CE(I)=-CESD(NSD)
           else
           CE(I)=CESD(NSD)
           end if
        end do
        do I=1,3
           NN=node_nside_nelm(I,NE)
           CE(I+3)=CEND(NN)
        end do

        call wf_mutensr(NE,MU)

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

  do NA=1,nant_max
     CTIMP=CTIMP+CIMP(NA)
  end do

  RETURN
END SUBROUTINE wf_pwrrad

!     ******* OUTPUT FIELD DATA *******

  SUBROUTINE wf_lpefld

    USE wfcomm
    IMPLICIT NONE
    INTEGER:: NS,NA

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
    DO NA=1,nant_max
       WRITE(6,140) NA,JNUM(NA),AJ(NA),APH(NA),AWD(NA),APOS(NA), &
                               RJ(1,NA),ZJ(1,NA),CIMP(NA)
140    FORMAT(1H ,I2,I3,0PF8.2,F8.2,1X,4F7.4,2X,'(',1P2E12.4,')')
    END DO

    IF(NPRINT.LT.2) RETURN

    ! field output

    RETURN
  END SUBROUTINE wf_lpefld
  
!     ******* OUTPUT ELEMENT DATA *******

  SUBROUTINE wf_lpelmt

  use wfcomm
  implicit none

  integer :: I,J,NA

  IF(NPRINT.LT.3) RETURN
  IF(nrank.NE.0) RETURN
     
  WRITE(6,110) node_max
110 FORMAT(/' ','NODE DATA     : #### node_max =',I5,' ####'/&
         &       ' ',2('  node_max',' KANOD',&
         &       9X,'R',14X,'Z',9X))
  WRITE(6,115) (I,KANOD(I),xnode(I),ynode(I),&
       &              I=1,node_max)
115 FORMAT((' ',2(2I6,2X,1P2E15.7,2X)))
  
  WRITE(6,120) nelm_max,(I,(node_nside_nelm(J,I),J=1,3),I=1,nelm_max)
120 FORMAT(/' ','ELEMENT DATA  : #### nelm_max =',I5,' ####'/&
         &      (' ',4(I6,'(',3I5,')',2X)))
  
  WRITE(6,125) nelm_max,(I,(nseg_nside_nelm(J,I),J=1,3),I=1,nelm_max)
125 FORMAT(/' ','SIDE    DATA  : #### nelm_max =',I5,' ####'/&
         &      (' ',2(I8,'(',3I8,')',2X)))
  
  DO NA=1,nant_max
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
  END SUBROUTINE wf_lpelmt

END MODULE wfwave
