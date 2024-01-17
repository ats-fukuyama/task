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
  integer    :: nelm,nant,IV,I,np
  real(rkind)    :: RW,PHASE,MU(3,3,6),A(3),B(3),C(3)
  real(rkind)    :: x1,y1,x2,y2,LIF(3),x21,y21
  complex(rkind) :: CJ(3),CVJ

  RW=2.D0*PI*RF*1.D6

  DO nelm=1,nelm_max
     DO IV=1,6
        CVTOT(IV,nelm)=(0.d0,0.d0)
     ENDDO
  ENDDO
  
  DO nant=1,nant_max
     PHASE =APH(nant)*PI/180.D0
     CVJ=CII*RW*RMU0*AJ(nant)*EXP(CII*(PHASE))
     IF(np_max_nant(nant).EQ.1) THEN
        nelm=nelm_np_nant(1,nant)
        CALL wf_set_abc(nelm,A,B,C)
        x1=x_np_nant(1,nant)
        y1=y_np_nant(1,nant)
        CJ(1)=0.D0
        CJ(2)=CVJ*RR
        CJ(3)=0.D0
        SELECT CASE(MODELG)
        CASE(0,12)
           do I=1,3
              LIF(I)=A(I)*RR &
                    +B(I)*x1*RR &
                    +C(I)*y1*RR
           end do
        CASE(1:10,13)
           do I=1,3
              LIF(I)=A(I)*x1 &
                    +B(I)*x1*x1 &
                    +C(I)*y1*x1
           end do
        END SELECT
        call wf_mutensr(nelm,MU)
        do I=1,3
           do IV=1,6
              CVTOT(IV,nelm)= CVTOT(IV,nelm)&
                           +LIF(I)*( MU(I,1,IV)*CJ(1)&
                                    +MU(I,2,IV)*CJ(2)&
                                    +MU(I,3,IV)*CJ(3))
           end do
        end do
     ELSE
     DO np=2,np_max_nant(nant)
        nelm=nelm_np_nant(np,nant)
        CALL wf_set_abc(nelm,A,B,C)
        x1=x_np_nant(np-1,nant)
        y1=y_np_nant(np-1,nant)
        x2=x_np_nant(np,nant)
        y2=y_np_nant(np,nant)
        x21=x2-x1
        y21=y2-y1

        CJ(1)=CVJ*x21  
        CJ(2)=(0.d0,0.d0) 
        CJ(3)=CVJ*y21  

        SELECT CASE(MODELG)
        CASE(0,12)
           do I=1,3
              LIF(I)=A(I)*RR &
                    +B(I)*(x1+x2)*RR/2.D0 &
                    +C(I)*(y1+y2)*RR/2.D0
           end do
        CASE(1:10,13)
           do I=1,3
              LIF(I)=A(I)*(x1+x2)/2.d0 &
                    +B(I)*(x2**2+x1*x2+x1**2)/3.d0 &
                    +C(I)*(x2*y1+y2*x1+2.d0*x2*y2+2.d0*x1*y1)/6.d0
           end do
        END SELECT
        call wf_mutensr(nelm,MU)

!        WRITE(16,*) NE
!        DO I=1,3
!           DO J=1,3
!              WRITE(16,'(2I3,1P6E12.4)') I,J,MU(I,J,1:6)
!           END DO
!        END DO

        do I=1,3
           do IV=1,6
              CVTOT(IV,nelm)= CVTOT(IV,nelm)&
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
  integer :: node,nseg,nvar

  DO nseg=1,nseg_max
     CESD_nseg(nseg)=(0.d0,0.d0)
  ENDDO
  DO node=1,node_max
     CEND_node(node) =(0.d0,0.d0)
  END DO

  DO nseg=1,nseg_max
     nvar=nvar_nseg(nseg)
     if (nvar.eq.0) then
        IF(nbdy_nseg(nseg).NE.0) THEN
           CESD_nseg(nseg)=CESD_nbdy(nbdy_nseg(nseg))
        ELSE
           CESD_nseg(nseg)=(0.d0,0.d0)
        END IF
     else
        CESD_nseg(nseg)=CSV(nvar)
     end if
!     if(nrank.eq.0) write(6,*) nseg,CESD_nseg(nseg),nbdy_nseg(nseg)
  END DO
  DO node=1,node_max
     nvar=nvar_node(node)
     if (nvar.eq.0) then
        IF(nbdy_node(node).NE.0) THEN
           CEND_node(node)=CEND_nbdy(nbdy_node(node))
        ELSE
           CEND_node(node)=(0.d0,0.d0)
        END IF
     else
        CEND_node(node)=CSV(nvar)
     end if
!     if(nrank.eq.0) write(6,*) node,CEND_node(node),mode_node(node)
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

  integer    :: nelm,IN,node,nseg,ns,nside
  integer    :: I,J,K,II,JJ
  real(rkind)    :: RW,S,MU(3,3,6),x(3),y(3)
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
     WRITE(6,'(A,2I8)') 'XX wf_pwrabs error: nelm_max,ipos=',nelm_max,ipos
     STOP
  END IF

  ! --- initialize ---
  
  RW=2.D0*PI*RF*1.D6
  CIWE=CII*RW*EPS0

  do nelm=nelm_pos_nrank(nrank)+1,nelm_pos_nrank(nrank)+nelm_len_nrank(nrank)
     S=area_nelm(nelm)

     ! --- calculate conductivity tensor ---

     call wf_dtensr(nelm,DTENS)
     do NS=1,NSMAX
        do IN=1,3
           do J=1,3
              do I=1,3
                 CTENS(NS,IN,I,J)=-CIWE*DTENS(NS,IN,I,J)
              end do
           end do
        end do
     end do

     call wf_set_node(nelm,x,y)
     call wf_mutensr(nelm,MU)
     
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
                                          *x(J)*S*AIF3(I,J,K)
                       end do
                    end do
                 end do
              end do
           end do
        END SELECT
     end do

     do nside=1,3
        nseg=nseg_nside_nelm(nside,nelm)
        if(nseg.lt.0) then
           nseg=-nseg
           CE(I)=-CESD_nseg(nseg)
        else
           CE(I)=CESD_nseg(nseg)
        end if
     end do
     do nside=1,3
        node=node_nside_nelm(nside,nelm)
        CE(I+3)=CEND_node(node)
     end do

     do ns=1,NSMAX
        pabs_ns_nelm(ns,nelm)=0.d0
        do JJ=1,6
           do II=1,6
              pabs_ns_nelm(ns,nelm)=pabs_ns_nelm(ns,nelm) &
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
     rdata(1:ndata)=pabs_ns_nelm(ns,nelm1:nelm2)
     CALL mtx_allgatherv_real8(rdata,ndata,rdata_tot,nelm_max, &
          nelm_len_nrank,nelm_pos_nrank)
     pabs_ns_nelm(ns,1:nelm_max)=rdata_tot(1:nelm_max)
  END DO

  do ns=1,NSMAX
     pabs_ns(ns)=0.d0
     do nelm=1,nelm_max
        pabs_ns(ns)=pabs_ns(ns)+pabs_ns_nelm(ns,nelm)
     end do
  end do

  pabs_tot=0.D0
  DO ns=1,NSMAX
     pabs_tot=pabs_tot+pabs_ns(ns)
  END DO

  RETURN
END SUBROUTINE wf_pwrabs

!     ******* POWER RADIATION *******

  SUBROUTINE wf_pwrrad

    use wfcomm
    USE wfsolv
    USE wfsub
  implicit none

  integer    :: nelm,nant,node,nseg
  integer    :: nside,I,IV,np
  real(rkind)    :: PHASE,RW,LIF(3),A(3),B(3),C(3)
  real(rkind)    :: x1,x2,y1,y2,x21,y21,MU(3,3,6)
  complex(rkind) :: CE(6),CJ(3),CVJ

  ! --- initialize ---

  RW=2.D0*PI*RF*1.D6

  do nant=1,nant_max
     PHASE =APH(nant)/180.D0
     CVJ=AJ(nant)*EXP(CII*(PHASE))
     cimp_nant(nant)=(0.d0,0.d0)
     do np=2,np_max_nant(nant)
        nelm=nelm_np_nant(np,nant)
        x1=x_np_nant(np-1,nant)
        y1=y_np_nant(np-1,nant)
        x2=x_np_nant(np,nant)
        y2=y_np_nant(np,nant)
        x21=x2-x1
        y21=y2-y1
        CJ(1)=CVJ*x21
        CJ(2)=(0.d0,0.d0)
        CJ(3)=CVJ*y21

        call wf_set_abc(nelm,A,B,C)

        SELECT CASE(MODELG)
        CASE(0,12)
           do I=1,3
              LIF(I)= A(I)*RR &
                     +B(I)*RR*(x2+x1)/2.d0 &
                     +C(I)*RR*(y1+y2)/2.d0           
           end do
        CASE(1:10,13)
           do I=1,3
              LIF(I)= A(I)*(x1+x2)/2.d0 &
                     +B(I)*(x2**2+x1*x2+x1**2)/3.d0 &
                     +C(I)*(x2*y1+y2*x1+2.d0*x2*y2+2.d0*x1*y1)/6.d0           
           end do
        END SELECT
        do nside=1,3
           nseg=nseg_nside_nelm(nside,nelm)
           if(nseg.lt.0) then
              nseg=-nseg
              CE(I)=-CESD_nseg(nseg)
           else
              CE(I)= CESD_nseg(nseg)
           end if
        end do
        do nside=1,3
           node=node_nside_nelm(nside,nelm)
           CE(I+3)=CEND_node(node)
        end do

        call wf_mutensr(nelm,MU)

        do I=1,3
           do IV=1,6
              cimp_nant(nant)=cimp_nant(nant)&
                         -0.5d0*LIF(I)*CONJG(CE(IV))&
                                      *( MU(I,1,IV)*CJ(1)&
                                        +MU(I,2,IV)*CJ(2)&
                                        +MU(I,3,IV)*CJ(3))
           end do
        end do

     end do
  end do

  cimp_tot=(0.d0,0.d0)
  do nant=1,nant_max
     cimp_tot=cimp_tot+cimp_nant(nant)
  end do

  RETURN
END SUBROUTINE wf_pwrrad

!     ******* OUTPUT FIELD DATA *******

  SUBROUTINE wf_lpefld

    USE wfcomm
    IMPLICIT NONE
    INTEGER:: NS,nant

    IF(NPRINT.LT.1) RETURN
    IF(nrank.NE.0) RETURN
    
!    WRITE(6,110) (EMAX(I),I=1,3),ETMAX,PNMAX
!110 FORMAT(1H ,'EXMAX  =',1PE12.4 &
!         ,3X ,'EYMAX  =',1PE12.4 &
!         ,3X ,'EZMAX  =',1PE12.4/ &
!         1H ,'EMAX   =',1PE12.4 &
!         ,3X ,'PNMAX  =',1PE12.4)

    WRITE(6,120) DBLE(cimp_tot),pabs_tot
120 FORMAT(1H ,'RADIATED POWER =',1PE12.4/ &
         1H ,'ABSORBED POWER =',1PE12.4)

    DO ns=1,NSMAX
       WRITE(6,126) NS,pabs_ns(ns)
126    FORMAT(1H ,'      PABS(',I2,') =',1PE12.4)
    END DO

    WRITE(6,130)
130 FORMAT(1H ,' I JNUM', '  AJ(I)','  APH(I)','  AWD(I)', &
            ' APOS(I)',' XJ(I)','  YJ(I)', &
            8X,'LOADING IMP.[ohm]')
    DO nant=1,nant_max
       WRITE(6,140) nant,np_max_nant(nant), &
            AJ(nant),APH(nant),AWD(nant),APOS(nant), &
            x_np_nant(1,nant),y_np_nant(1,nant),cimp_nant(nant)
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

  integer :: nelm,np,np0,nside,nant,node

  IF(NPRINT.LT.3) RETURN
  IF(nrank.NE.0) RETURN
     
  WRITE(6,110) node_max
110 FORMAT(/' ','NODE DATA     : #### node_max =',I5,' ####'/&
            ' ',2('  node_max',' mode_node',&
                9X,'X',14X,'Y',9X))
  WRITE(6,115) (node,mode_node(node),xnode(node),ynode(node),&
       &              node=1,node_max)
115 FORMAT((' ',2(2I6,2X,1P2E15.7,2X)))
  
  WRITE(6,120) nelm_max,(nelm,(node_nside_nelm(nside,nelm),nside=1,3), &
       nelm=1,nelm_max)
120 FORMAT(/' ','ELEMENT DATA  : #### nelm_max =',I5,' ####'/&
         &      (' ',4(I6,'(',3I5,')',2X)))
  
  WRITE(6,125) nelm_max,(nelm,(nseg_nside_nelm(nside,nelm),nside=1,3), &
       nelm=1,nelm_max)
125 FORMAT(/' ','SIDE    DATA  : #### nelm_max =',I5,' ####'/&
         &      (' ',2(I8,'(',3I8,')',2X)))
  
  DO nant=1,nant_max
     WRITE(6,140) nant,np0_max_nant(nant)
140  FORMAT(/' ','ORIGINAL ANTENNA DATA : nant =',I5,' JNUM0 =',I5/&
          &          ' ',2('  NO.',13X,' RJ0',11X,' ZJ0',6X))
     WRITE(6,150) (np0,x_np0_nant(np0,nant),y_np0_nant(np0,nant), &
          np0=1,np0_max_nant(nant))
150  FORMAT((' ',2(I5,8X,1P2E15.7)))
     
     WRITE(6,154) nant,np_max_nant(nant)
154  FORMAT(/' ','MODIFIED ANTENNA DATA : nant =',I5,' JNUM  =',I5/&
          &          ' ',2('  NO.',' JELM',8X,' JR ',11X,' JZ ',6X))
     WRITE(6,156) (np,np_max_nant(nant), &
          x_np_nant(np,nant),y_np_nant(np,nant),np=1,np_max_nant(nant))
156  FORMAT((' ',2(2I5,3X,1P2E15.7)))
  ENDDO
  
  RETURN
  END SUBROUTINE wf_lpelmt

END MODULE wfwave
