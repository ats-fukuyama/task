!     $Id$

! *****************************
!     PREPARATION OF FPLOOP
! *****************************
      MODULE fpprep

      USE fpcomm
      USE fpinit
      USE fpsave
      USE fpcoef
      USE equnit_mod

      contains

!     ***** create mesh quantities *****

      SUBROUTINE fp_mesh(IERR)

      USE libmtx
      USE plprof
      USE fpbroadcast
      USE fpwrin
      USE fpwmin

      Implicit none
      integer :: ierr,NSA,NSB,NS,NR,NP,NTH,id,NP2
      character(LEN=80)::line 
      real(8)::rhon,rhol,rhol1,rhol2,A1,epsl,fact,rint0,es,rint2,ql,BT
      real(8),DIMENSION(:),POINTER:: work,workg
      real(8)::suml, delh, etal, x, psib, pcos
      integer:: NG

!     ----- exec EQ -----

      IF(MODELG.EQ.3) THEN
         IF(nrank.EQ.0) THEN
            CALL eq_load(MODELG,KNAMEQ,IERR)
            IF(IERR.NE.0) THEN
               write(6,*) 'XX FPMESH:EQLOAD:IERR=',IERR
               RETURN
            ENDIF
         ENDIF
         CALL fp_eq_broadcast
!         write(LINE,'(A,I5)') 'nrmax=',51
!         call eq_parm(2,line,ierr)
!         write(LINE,'(A,I5)') 'nthmax=',64
!         call eq_parm(2,line,ierr)
!         write(LINE,'(A,I5)') 'nsumax=',64
!         call eq_parm(2,line,ierr)
         CALL eqcalq(IERR)
         CALL eqgetb(BB,RR,RIP,RA,RKAP,RDLT,RB)
      ENDIF

!     ----- set radial mesh -----

      IF(NRMAX.EQ.1) THEN
         DELR=DELR1
      ELSE
         DELR=(RMAX-RMIN)/NRMAX
      ENDIF

      IF(NRMAX.EQ.1) THEN
         RM(1)=R1
         RG(1)=R1-0.5D0*DELR
         RG(2)=R1+0.5D0*DELR
      ELSE
         DO NR=1,NRMAX
            RM(NR)=RMIN+DELR*(NR-1)+0.5D0*DELR
            RG(NR)=RMIN+DELR*(NR-1)
         ENDDO
         RG(NRMAX+1)=RMAX
      ENDIF

!     ----- load WR resluts -----

      ID=0
      DO NSA=1,NSAMAX
         NS=NS_NSA(NSA)
         IF(MODELW(NS).EQ.1.OR.MODELW(NS).EQ.2) ID=1
      ENDDO
      IF(ID.EQ.1) THEN
         CALL fp_wr_read(IERR)
         IF(IERR.NE.0) RETURN
      ENDIF

!     ----- load WM resluts -----

      ID=0
      DO NSA=1,NSAMAX
         NS=NS_NSA(NSA)
         IF(MODELW(NS).EQ.4) ID=1
      ENDDO
      IF(ID.EQ.1) THEN
         CALL fp_wm_read(IERR)
         IF(IERR.NE.0) RETURN
      ENDIF


!     ----- set approximate poloidal magneticl field -----

      DO NR=1,NRMAX
         RHON=RM(NR)
         CALL pl_qprf(RHON,QL)
         BT=BB
         BP= RSRHON(RHON)*BT/(RR*QL)
         EPSRM(NR)=RSRHON(RHON)/RR
      ENDDO
      RHON=RG(NRMAX+1)
      CALL pl_qprf(RHON,QL)
      BT=BB
      BP= RSRHON(RHON)*BT/(RR*QL)
      EPSRM(NRMAX+1)=RSRHON(RHON)/RR

      DO NR=1,NRMAX+1
         RHON=RG(NR)
         CALL pl_qprf(RHON,QL)
         BT=BB
         BP(NR)= RSRHON(RHON)*BT/(RR*QL)
         EPSRG(NR)=RSRHON(RHON)/RR
      ENDDO

!     ----- set parallel current density -----

      DO NR=1,NRMAX
         RJ1(NR)=(RG(NR+1)*BP(NR+1)-RM(NR)*BP(NR)) &
               /(RMU0*RM(NR)*DELR)
      ENDDO

!     ----- set parallel electric field -----

      DO NR=1,NRMAX
         E1(NR)=E0
         E2(NR)=E0
         RJ2(NR)=RJ1(NR)
      ENDDO

!     ----- set momentum space mesh -----

      DELTH=PI/NTHMAX
      DO NSB=1,NSBMAX
         DELP(NSB)=PMAX(NSB)/NPMAX
         DO NP=1,NPMAX
            PG(NP,NSB)=DELP(NSB)*(NP-1)
            PM(NP,NSB)=DELP(NSB)*(NP-0.5D0)
         ENDDO
!         DO NP2=1,NP2MAX
!            PG2(NP2,NSB)=PG(2,NSB)/DBLE(NP2MAX)*(NP2-1)
!            PM2(NP2,NSB)=PG(2,NSB)/DBLE(NP2MAX)*(NP2-0.5D0)
!         END DO
!         PG2(NP2MAX+1,NSB)=PG(2,NSB)
         PG(NPMAX+1,NSB)=PMAX(NSB)
      ENDDO

      DO NTH=1,NTHMAX
         THG(NTH)=DELTH*(NTH-1)
         THM(NTH)=DELTH*(NTH-0.5D0)

         SINM(NTH)=SIN(THM(NTH))
         COSM(NTH)=COS(THM(NTH))
         SING(NTH)=SIN(THG(NTH))
         COSG(NTH)=COS(THG(NTH))
      ENDDO
      THG(NTHMAX+1)=PI
      SING(NTHMAX+1)=0.D0
      COSG(NTHMAX+1)=-1.D0

      DO NSB=1,NSBMAX
      DO NP=1,NPMAX
      DO NTH=1,NTHMAX
         VOLP(NTH,NP,NSB)=2.D0*PI*SINM(NTH)*PM(NP,NSB)**2*DELP(NSB)*DELTH
      ENDDO
      ENDDO
      ENDDO

      DO NR=1,NRMAX
         RHOL=RM(NR)
         RHOL1=RG(NR)
         RHOL2=RG(NR+1)
         VOLR(NR)=2.D0*PI*RSRHON(RHOL)*(RSRHON(RHOL2)-RSRHON(RHOL1)) &
                *2.D0*PI*RR
      ENDDO
      TVOLR=0.D0
      DO NR=1,NRMAX
         TVOLR=TVOLR+VOLR(NR)
      ENDDO

!     ----- set bounce-average parameters -----

      IF (MODELA.EQ.0) THEN
         DO NR=1,NRMAX+1
            ITL(NR)=0
            ITU(NR)=0
!         ENDDO
!         DO NR=1,NRMAX+1
            ITLG(NR)=0
            ITUG(NR)=0
         ENDDO
      ELSE
         DO NR=1,NRMAX+1
            A1=ACOS(SQRT(2.D0*EPSRM(NR)/(1.D0+EPSRM(NR))))
            DO NTH=1,NTHMAX/2
               IF (THG(NTH).LE.A1.AND.THG(NTH+1).GE.A1) GOTO 201
            ENDDO
            NTH=NTHMAX/2-1
  201       CONTINUE

            ITL(NR)=NTH
            ITU(NR)=NTHMAX-NTH+1

            EPSL=COSM(ITL(NR))**2/(2.D0-COSM(ITL(NR))**2)
            IF(nprocs.gt.1.and.NRANK.eq.1) &
                 WRITE(6,'(A,2I5,1P2E12.4)') 'NR,NTHC,EPSRM=',NR,NTH,EPSRM(NR),EPSL
!            WRITE(6,'(A,1P2E12.4)') 'EPSRM(NR)=',EPSRM(NR),EPSL
!            WRITE(6,'(A,2I8)') 'ITL,ITU=',ITL(NR),ITU(NR)
            EPSRM(NR)=EPSL
         ENDDO

         IF(NRANK.eq.1) WRITE(6,*) " "

         DO NR=1,NRMAX+1
            A1=ACOS(SQRT(2.D0*EPSRG(NR)/(1.D0+EPSRG(NR))))
            DO NTH=1,NTHMAX/2
               IF (THG(NTH).LE.A1.AND.THG(NTH+1).GE.A1) GOTO 202
            ENDDO
            NTH=NTHMAX/2-1
  202       CONTINUE

            ITLG(NR)=NTH
            ITUG(NR)=NTHMAX-NTH+1

            EPSL=COSM(ITLG(NR))**2/(2.D0-COSM(ITLG(NR))**2)
            IF(nprocs.gt.1.and.NRANK.eq.1) &
                 WRITE(6,'(A,2I5,1P2E12.4)') 'NR,NTHC,EPSRG=',NR,NTH,EPSRG(NR),EPSL

            EPSRG(NR)=EPSL
         ENDDO

      ENDIF

      IF (MODELA.EQ.0) THEN
         DO NR=NRSTART,NREND
            DO NTH=1,NTHMAX
               ETAM(NTH,NR)=PI/2.D0
               RLAMDA(NTH,NR)=1.D0
               RLAMDC(NTH,NR)=1.D0
            ENDDO
            DO NTH=1,NTHMAX+1
               ETAG(NTH,NR)=PI/2.D0
            ENDDO
         ENDDO

         DO NR=NRSTART,NREND
            DO NTH=1,NTHMAX
               ETAM_G(NTH,NR)=PI/2.D0
               RLAMDA_G(NTH,NR)=1.D0
               RLAMDC_G(NTH,NR)=1.D0
            ENDDO
            DO NTH=1,NTHMAX+1
               ETAG_G(NTH,NR)=PI/2.D0
            ENDDO
         ENDDO
      ELSE
         DO NR=NRSTART,NREND
            CALL SET_RLAMDA(NR)
         ENDDO
!         CALL SET_RLAMDA(NRMAX+1)
         DO NR=NRSTART,NREND
            CALL SET_RLAMDA_G(NR)
         ENDDO
      END IF

      allocate(work(nrstart:nrendx),workg(NRMAX))

      DO NTH=1,NTHMAX
         DO NR=NRSTART,NRENDX
            work(NR)=RLAMDA(NTH,NR)
         ENDDO
         CALL mtx_gatherv_real8(work(NRSTART:NRENDX),MTXLEN(NRANK+1), &
                                workg,NRMAX,MTXLEN,MTXPOS)
         CALL mtx_broadcast_real8(workg,NRMAX)
         DO NR=1,NRMAX
            RLAMDAG(NTH,NR)=workg(NR)
         ENDDO
         RLAMDAG(NTH,NRMAX+1)=RLAMDAG(NTH,NRMAX)
      ENDDO

      DO NTH=1,NTHMAX
         DO NR=NRSTART,NRENDX
            work(NR)=ETAM(NTH,NR)
         ENDDO
         CALL mtx_gatherv_real8(work(NRSTART:NRENDX),MTXLEN(NRANK+1), &
                                workg,NRMAX,MTXLEN,MTXPOS)
         CALL mtx_broadcast_real8(workg,NRMAX)
         DO NR=1,NRMAX
            ETAMG(NTH,NR)=workg(NR)
         ENDDO
      ENDDO

      DO NTH=1,NTHMAX
         DO NR=NRSTART,NRENDX
            work(NR)=RLAMDA_G(NTH,NR)
         ENDDO
         CALL mtx_gatherv_real8(work(NRSTART:NRENDX),MTXLEN(NRANK+1), &
                                workg,NRMAX,MTXLEN,MTXPOS)
         CALL mtx_broadcast_real8(workg,NRMAX)
         DO NR=1,NRMAX
            RLAMDA_GG(NTH,NR)=workg(NR)
         ENDDO
         RLAMDA_GG(NTH,NRMAX+1)=RLAMDA_GG(NTH,NRMAX)
      ENDDO

      DO NTH=1,NTHMAX
         DO NR=NRSTART,NRENDX
            work(NR)=ETAM_G(NTH,NR)
         ENDDO
         CALL mtx_gatherv_real8(work(NRSTART:NRENDX),MTXLEN(NRANK+1), &
                                workg,NRMAX,MTXLEN,MTXPOS)
         CALL mtx_broadcast_real8(workg,NRMAX)
         DO NR=1,NRMAX
            ETAM_GG(NTH,NR)=workg(NR)
         ENDDO
      ENDDO

      deallocate(work,workg)

      IERR=0
      RETURN
      END subroutine fp_mesh

!
!-----------------------------
!
      SUBROUTINE SET_RLAMDA(NR)

      USE libmtx
      USE plprof
!      USE fpbroadcast
!      USE fpwrin
!      USE fpwmin

      IMPLICIT NONE
      integer:: NR, NTH, NG
      real(8):: EPSL, FACT, A1, RINT0, ES, DELH
      real(8):: SUML, ETAL, X, PSIB, PCOS, RINT2

      EPSL=EPSRM(NR)
      FACT=(1.D0+EPSL)/(2.D0*EPSL)
      
      DO NTH=1,ITL(NR)
         ETAM(NTH,NR)=PI/2.D0
      ENDDO
      
      DO NTH=ITL(NR)+1,ITU(NR)-1
         A1=FACT*COSM(NTH)**2
         ETAM(NTH,NR)=0.5D0*DACOS(1.D0-2.D0*A1)
      ENDDO
      
      DO NTH=ITU(NR),NTHMAX
         ETAM(NTH,NR)=PI/2.D0
      ENDDO
      
      DO NTH=1,ITL(NR)
         ETAG(NTH,NR)=PI/2.D0
      ENDDO
      
      DO NTH=ITL(NR)+1,ITU(NR)
         A1=FACT*COSG(NTH)**2
         ETAG(NTH,NR)=0.5D0*DACOS(1.D0-2.D0*A1)
      ENDDO
      
      DO NTH=ITU(NR)+1,NTHMAX+1
         ETAG(NTH,NR)=PI/2.D0
      ENDDO
      
      NRX=NR
      DO NTH=1,NTHMAX/2
         NTHX=NTH
         IF(NTH.LT.ITL(NR)) THEN
            CALL DEFT(RINT0,ES,H0DE,EPSDE,0,FPFN0U)
         ELSEIF(NTH.ge.ITL(NR)) THEN
!                  CALL DEFT(RINT0,ES,H0DE,EPSDE,0,FPFN0T)
            DELH=2.D0*ETAM(NTH,NR)/NAVMAX
            suml=0.D0
            DO NG=1,NAVMAX
               ETAL=DELH*(NG-0.5D0)
               X=EPSRM(NR)*COS(ETAL)*RR
               PSIB=(1.D0+EPSRM(NR))/(1.D0+X/RR)
               IF(COSM(NTH).ge.0.D0)THEN
                  PCOS=SQRT(1.D0-PSIB*SINM(NTH)**2)
               ELSE
                  PCOS=-SQRT(1.D0-PSIB*SINM(NTH)**2)
               END IF
               suml=suml+COSM(NTH)/PCOS
            END DO
            RINT0=suml*DELH/ABS(COSM(NTH))
         ELSE
            RINT0=0.D0
         ENDIF
         CALL DEFT(RINT2,ES,H0DE,EPSDE,0,FPFN2A)
         RLAMDA(NTH,NR)=RINT0*ABS(COSM(NTH))/PI
         RLAMDC(NTH,NR)=RINT2/(PI*(1.D0+EPSRM(NR))*(COSG(NTH)))
      ENDDO
      IF(ITL(NR).ne.NTHMAX/2)THEN
         RLAMDA(ITL(NR),NR)=0.5D0*(RLAMDA(ITL(NR)-1,NR) &
              +RLAMDA(ITL(NR)+1,NR))
      ELSE!IF(ITL(NR).eq.NTHMAX/2)THEN
         RLAMDA(ITL(NR),NR)=RINT0*ABS(COSM(ITL(NR)))/PI
      END IF
      DO NTH=1,NTHMAX/2
         RLAMDA(NTHMAX-NTH+1,NR)=RLAMDA(NTH,NR)
         RLAMDC(NTHMAX-NTH+2,NR)=RLAMDC(NTH,NR)
      ENDDO
      RLAMDC(NTHMAX/2+1,NR)=0.D0

      END SUBROUTINE SET_RLAMDA

!
!-----------------------------
!
      SUBROUTINE SET_RLAMDA_G(NR)

      USE libmtx
      USE plprof
!      USE fpbroadcast
!      USE fpwrin
!      USE fpwmin

      IMPLICIT NONE
      integer:: NR, NTH, NG
      real(8):: EPSL, FACT, A1, RINT0, ES, DELH
      real(8):: SUML, ETAL, X, PSIB, PCOS, RINT2

      EPSL=EPSRG(NR)
      FACT=(1.D0+EPSL)/(2.D0*EPSL)
      
      DO NTH=1,ITLG(NR)
         ETAM_G(NTH,NR)=PI/2.D0
      ENDDO
      
      DO NTH=ITLG(NR)+1,ITUG(NR)-1
         A1=FACT*COSM(NTH)**2
         ETAM_G(NTH,NR)=0.5D0*DACOS(1.D0-2.D0*A1)
      ENDDO
      
      DO NTH=ITUG(NR),NTHMAX
         ETAM_G(NTH,NR)=PI/2.D0
      ENDDO
      
      DO NTH=1,ITLG(NR)
         ETAG_G(NTH,NR)=PI/2.D0
      ENDDO
      
      DO NTH=ITLG(NR)+1,ITUG(NR)
         A1=FACT*COSG(NTH)**2
         ETAG_G(NTH,NR)=0.5D0*DACOS(1.D0-2.D0*A1)
      ENDDO
      
      DO NTH=ITUG(NR)+1,NTHMAX+1
         ETAG_G(NTH,NR)=PI/2.D0
      ENDDO
      
      NRX=NR
      DO NTH=1,NTHMAX/2
         NTHX=NTH
         IF(NTH.LT.ITLG(NR)) THEN
            CALL DEFT(RINT0,ES,H0DE,EPSDE,0,FPFN0U_G)
         ELSEIF(NTH.ge.ITLG(NR)) THEN
            DELH=2.D0*ETAM_G(NTH,NR)/NAVMAX
            suml=0.D0
            DO NG=1,NAVMAX
               ETAL=DELH*(NG-0.5D0)
               X=EPSRG(NR)*COS(ETAL)*RR
               PSIB=(1.D0+EPSRG(NR))/(1.D0+X/RR)
               IF(COSM(NTH).ge.0.D0)THEN
                  PCOS=SQRT(1.D0-PSIB*SINM(NTH)**2)
               ELSE
                  PCOS=-SQRT(1.D0-PSIB*SINM(NTH)**2)
               END IF
               suml=suml+COSM(NTH)/PCOS
            END DO
            RINT0=suml*DELH/ABS(COSM(NTH))
         ELSE
            RINT0=0.D0
         ENDIF
         CALL DEFT(RINT2,ES,H0DE,EPSDE,0,FPFN2A_G)
         RLAMDA_G(NTH,NR)=RINT0*ABS(COSM(NTH))/PI
         RLAMDC_G(NTH,NR)=RINT2/(PI*(1.D0+EPSRG(NR))*(COSG(NTH)))
      ENDDO
      IF(ITLG(NR).ne.NTHMAX/2)THEN
         RLAMDA_G(ITLG(NR),NR)=0.5D0*(RLAMDA_G(ITLG(NR)-1,NR) &
              +RLAMDA_G(ITLG(NR)+1,NR))
      ELSE!IF(ITL(NR).eq.NTHMAX/2)THEN
         RLAMDA_G(ITLG(NR),NR)=RINT0*ABS(COSM(ITLG(NR)))/PI
      END IF
      DO NTH=1,NTHMAX/2
         RLAMDA_G(NTHMAX-NTH+1,NR)=RLAMDA_G(NTH,NR)
         RLAMDC_G(NTHMAX-NTH+2,NR)=RLAMDC_G(NTH,NR)
      ENDDO
      RLAMDC_G(NTHMAX/2+1,NR)=0.D0

      END SUBROUTINE SET_RLAMDA_G

! *************************
!     INITIAL DATA SAVE
! *************************

      SUBROUTINE FPCINI

      IMPLICIT NONE
      integer:: NSA,NR, NTH, NP

      ISAVE=0

      DO NSA=1,NSAMAX
      DO NR=NRSTART,NREND
         DO NP=1,NPMAX+1
         DO NTH=1,NTHMAX
            DPP(NTH,NP,NR,NSA)=0.D0
            DPT(NTH,NP,NR,NSA)=0.D0
            FPP(NTH,NP,NR,NSA)=0.D0
         ENDDO
         ENDDO
         DO NP=1,NPMAX
         DO NTH=1,NTHMAX+1
            DTP(NTH,NP,NR,NSA)=0.D0
            DTT(NTH,NP,NR,NSA)=0.D0
            FTH(NTH,NP,NR,NSA)=0.D0
         ENDDO
         ENDDO
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE FPCINI

! ***************************************************************

!                       SET OF INTEGRAND

! ***************************************************************

! ============================================================

      REAL*8 FUNCTION  FPFN0U(X,XM,XP)

      IMPLICIT NONE
      real(8):: X, XM, XP, XX, A0, A1, A2

      XX=X
      XX=XM
      A0=ETAM(NTHX,NRX)
      A1=1.D0+EPSRM(NRX)*COS(A0*XP)
      A2=(1.D0+EPSRM(NRX))*SINM(NTHX)**2
      FPFN0U=A0*SQRT(A1/(A1-A2))

      RETURN
      END FUNCTION  FPFN0U

! ============================================================

      REAL*8 FUNCTION  FPFN0U_G(X,XM,XP)

      IMPLICIT NONE
      real(8):: X, XM, XP, XX, A0, A1, A2

      XX=X
      XX=XM
      A0=ETAM_G(NTHX,NRX)
      A1=1.D0+EPSRG(NRX)*COS(A0*XP)
      A2=(1.D0+EPSRG(NRX))*SINM(NTHX)**2
      FPFN0U_G=A0*SQRT(A1/(A1-A2))

      RETURN
      END FUNCTION  FPFN0U_G

! ============================================================

      REAL*8 FUNCTION FPFN0T(X,XM,XP)

      IMPLICIT NONE
      real(8):: X, XM, XP, XX, A0, A1, A2

      XX=X
      A0=ETAM(NTHX,NRX)
      A1=1.D0+EPSRM(NRX)*COS(A0*XP)
      A2=2.D0*EPSRM(NRX)*SIN(0.5D0*A0*XM)*SIN(0.5D0*A0*(XP+2.D0))
      FPFN0T=A0*SQRT(A1/A2)

      RETURN
      END FUNCTION FPFN0T

! ============================================================

      REAL*8 FUNCTION FPFN1A(X,XM,XP)

      IMPLICIT NONE
      real(8):: X, XM, XP, XX, A0, A1, A2

      XX=X
      XX=XM
      A0=ETAM(NTHX,NRX)
      A1=1.D0+EPSRM(NRX)*COS(A0*XP)
      FPFN1A=A0*SQRT(A1)

      RETURN
      END FUNCTION FPFN1A

! ============================================================

      REAL*8 FUNCTION FPFN2A(X,XM,XP)

      IMPLICIT NONE
      real(8):: X, XM, XP, XX, A0, A1, A2

      XX=X
      XX=XM
      A0=ETAG(NTHX,NRX)
      A1=1.D0+EPSRM(NRX)*COS(A0*XP)
      A2=A1-(1.D0+EPSRM(NRX))*SING(NTHX)**2
      FPFN2A=A0*SQRT(A1*A2)

      RETURN
      END FUNCTION FPFN2A
! ============================================================

      REAL*8 FUNCTION FPFN2A_G(X,XM,XP)

      IMPLICIT NONE
      real(8):: X, XM, XP, XX, A0, A1, A2

      XX=X
      XX=XM
      A0=ETAG_G(NTHX,NRX)
      A1=1.D0+EPSRG(NRX)*COS(A0*XP)
      A2=A1-(1.D0+EPSRG(NRX))*SING(NTHX)**2
      FPFN2A_G=A0*SQRT(A1*A2)

      RETURN
      END FUNCTION FPFN2A_G
!--------------------------------
      REAL*8 FUNCTION  FPFN2U(X,XM,XP)
                           
      IMPLICIT NONE
      real(8):: X, XM, XP, XX, A0, A1, A2

      XX=X
      XX=XM
      A0=ETAG(NTHX,NRX)
      A1=1.D0+EPSRM(NRX)*COS(A0*XP)
      A2=(1.D0+EPSRM(NRX))*SING(NTHX)**2
      FPFN2U=A0*SQRT(A1/(A1-A2))

      RETURN
      END FUNCTION  FPFN2U
!------------------------------
      REAL*8 FUNCTION FPFN2T(X,XM,XP)

      IMPLICIT NONE
      real(8):: X, XM, XP, XX, A0, A1, A2

      XX=X
      A0=ETAG(NTHX,NRX)
      A1=1.D0+EPSRM(NRX)*COS(A0*XP)
      A2=2.D0*EPSRM(NRX)*SIN(0.5D0*A0*XM)*SIN(0.5D0*A0*(XP+2.D0))
      FPFN2T=A0*SQRT(A1/A2)

      RETURN
      END FUNCTION FPFN2T

!-------------------------------------------------

      SUBROUTINE fp_prep(ierr)

      USE plprof
      USE fpnfrr
      USE libmtx
      USE libmtx
      Implicit none

      integer :: ierr,NSA,NSB,NS,NR,NP,NTH,NSFP,NSFD,NSBA,N,NREND1
      real(8) :: FL, RSUM1, RSUM2, RTFD0L, RHON, RNE, RTE
      real(8) :: RLNRL, FACT, RSUM, RSUM11, rsum3, rsum4
      TYPE(pl_plf_type),DIMENSION(NSMAX):: PLF
      INTEGER,DIMENSION(nprocs):: ima1,ima2,nra1,nra2,nma1,nma2
      real(8),DIMENSION(:),POINTER:: work,workg

!     ----- Initialize time counter -----

      TIMEFP=0.D0
      NTG1=0
      NTG2=0
!     ----- Check nprocs -----
      IF(nprocs.GT.nrmax) THEN
         IF(nrank.EQ.0) THEN
            WRITE(6,*) 'XX fp_prep: nrmax must be greater than nprocs.'
            WRITE(6,*) 'XX          nrmax,nprocs=',nrmax,nprocs
         ENDIF
         ierr=1
         RETURN
      ENDIF

!     ----- Set matrix size -----

      imtxsize=nthmax*npmax*nrmax
      IF(modeld.EQ.0) THEN
         imtxwidth=4*nthmax-1
      ELSE
         imtxwidth=4*nthmax*npmax-1
      ENDIF
      CALL mtx_setup(imtxsize,imtxstart,imtxend,imtxwidth)
      nrstart=(imtxstart-1)/(nthmax*npmax)+1
      nrend=  (imtxend  -1)/(nthmax*npmax)+1
      nrend1= (imtxend    )/(nthmax*npmax)+1
      IF(nrend1.EQ.nrend) THEN
         NRENDX=NREND-1
      ELSE
         NRENDX=NREND
      ENDIF
      nmstart=nthmax*npmax*(nrstart-1)+1
      nmend  =nthmax*npmax* nrend
      CALL mtx_gather_integer(imtxstart,ima1)
      CALL mtx_gather_integer(imtxend,  ima2)
      CALL mtx_gather_integer(nrstart,  nra1)
      CALL mtx_gather_integer(nrend,    nra2)
      CALL mtx_gather_integer(nmstart,  nma1)
      CALL mtx_gather_integer(nmend,    nma2)
      IF(nrank.EQ.0) THEN
         write(6,'(A,2I10)') '  imtxsize,imtxwidth=',imtxsize,imtxwidth
         write(6,'(A,A)') '     nrank   imtxstart   imtxend   nrstart',&
                          '     nrend   nmstart     nmend'
         DO N=1,nprocs
            write(6,'(7I10)') N,ima1(N),ima2(N),nra1(N),nra2(N),nma1(N),nma2(N)
         ENDDO
      ENDIF
      CALL mtx_cleanup

!     ----- Allocate variables -----

      CALL fp_allocate
      call fp_allocate_ntg1
      call fp_allocate_ntg2

!     ----- Get mtxlen and mtxpos -----

      CALL mtx_allgather_integer(nrend-nrstart+1,mtxlen)
      CALL mtx_allgather_integer(nrstart-1,mtxpos)
      if(nrank.eq.0) then
         DO N=1,NPROCS
            WRITE(6,'(A,3I10)') '  nrank,mtxpos,mtxlen = ', &
                                   n-1,mtxpos(n),mtxlen(n)
         ENDDO
!     ----- Check NS_NSA and NS_NSB -----

         DO NSA=1,NSAMAX
            IF(NS_NSA(NSA).EQ.0) THEN
               WRITE(6,*) 'XX NS_NSA(NSA)=0 for NSA=', NSA
               WRITE(6,*) '   NS_NSA(NSA)=NSA substituted'
               NS_NSA(NSA)=NSA
            ENDIF
         ENDDO

         DO NSB=1,NSBMAX
            IF(NS_NSB(NSB).EQ.0) THEN
               WRITE(6,*) 'XX NS_NSB(NSB)=0 for NSB=', NSB
               WRITE(6,*) '   NS_NSB(NSB)=NSB substituted'
               NS_NSB(NSB)=NSB
            ENDIF
         ENDDO
      ENDIF

!     ----- NSB_NSA and NSA_NSB array -----

      DO NSA=1,NSAMAX
         NSB_NSA(NSA)=0
      ENDDO
      DO NSB=1,NSBMAX
         NSA_NSB(NSB)=0
      ENDDO
      DO NSA=1,NSAMAX
         DO NSB=1,NSBMAX
            IF(NS_NSA(NSA).EQ.NS_NSB(NSB)) THEN
               NSB_NSA(NSA)=NSB
               NSA_NSB(NSB)=NSA
            ENDIF
         ENDDO
      ENDDO

      IF(nrank.EQ.0) THEN
         DO NSA=1,NSAMAX
            WRITE(6,'(A,3I5)') 'NSA,NS_NSA,NSB_NSA=', &
                                NSA,NS_NSA(NSA),NSB_NSA(NSA)
         ENDDO
         DO NSB=1,NSBMAX
            WRITE(6,'(A,3I5)') 'NSB,NS_NSB,NSA_NSB=', &
                                NSB,NS_NSB(NSB),NSA_NSB(NSB)
         ENDDO
         ierr=0

         DO NSA=1,NSAMAX
            IF(NSB_NSA(NSA).EQ.0) THEN
               WRITE(6,*) 'XX NS_NSA has no correponding NS_NSB for NSA=',NSA
               IERR=1
            ENDIF
         ENDDO
         if(ierr.ne.0) return

         WRITE(6,*) '--------------------'
         DO NSA=1,NSAMAX
            WRITE(6,'(A,2I3)') 'NSA,NS=',NSA,NS_NSA(NSA)
         ENDDO
         WRITE(6,*) '--------------------'
         DO NSB=1,NSBMAX
            WRITE(6,'(A,2I3)') 'NSB,NS=',NSB,NS_NSB(NSB)
         ENDDO
         WRITE(6,*) '--------------------'
      ENDIF

!     ----- create meches -----

      CALL fp_mesh(ierr)
!     ----- Initialize velocity distribution function of all species -----

      DO NSB=1,NSBMAX
         NS=NS_NSB(NSB)
         DO NR=1,NRMAX
            DO NP=1,NPMAX
               FL=FPMXWL(PM(NP,NSB),NR,NS)
               DO NTH=1,NTHMAX
                  FNS(NTH,NP,NR,NSB)=FL
               END DO
            ENDDO
         END DO
      END DO

!--------- normalize bounce average parameter ---------

      IF(MODELA.eq.1)THEN
         NSB=1 ! arbitrary
         DO NR=NRSTART,NREND
            RSUM1=0.D0
            RSUM2=0.D0
            RSUM3=0.D0
            RSUM4=0.D0
            DO NP=1,NPMAX
               DO NTH=1,NTHMAX
                  RSUM1 = RSUM1+VOLP(NTH,NP,NSB)&!*FNS(NTH,NP,NR,NSB) &
                       *RLAMDAG(NTH,NR)
                  RSUM2 = RSUM2+VOLP(NTH,NP,NSB)!*FNS(NTH,NP,NR,NSB)
                  RSUM3 = rsum3+VOLP(NTH,NP,NSB)*RLAMDA_GG(NTH,NR)
                  RSUM4 = rsum4+VOLP(NTH,NP,NSB)
               END DO
            END DO
            IF(RSUM1.EQ.0.D0) &
                 WRITE(6,'(1P3E12.4)') VOLP(1,1,NSB),FNS(1,1,1,NSB),RLAMDA(1,1)
            RCOEF(NR)=RSUM2/RSUM1
            RCOEF_G(NR)=RSUM4/RSUM3
!            write(*,*) rsum2/rsum1, rsum4/rsum3
         END DO
      ELSE
         DO NSB=1,NSBMAX
            DO NR=NRSTART,NREND
               RCOEF(NR)=1.D0
               RCOEF_G(NR)=1.D0
           ENDDO
         ENDDO
      END IF


!     ----- set boundary distribution functions -----

      IF(NRSTART.EQ.1) THEN
         DO NSA=1,NSAMAX
            NS=NS_NSA(NSA)
            NSBA=NSB_NSA(NSA)
            DO NP=1,NPMAX
               FL=FPMXWL(PM(NP,NSBA),0,NS)
               DO NTH=1,NTHMAX
                  FS1(NTH,NP,NSA)=FL ! at r=0
               ENDDO
            ENDDO
         END DO
      ELSE
         DO NSA=1,NSAMAX
            DO NP=1,NPMAX
               DO NTH=1,NTHMAX
                  FS1(NTH,NP,NSA)=0.D0
               ENDDO
            ENDDO
         ENDDO
      ENDIF

!      IF(NREND.EQ.NRMAX) THEN
         DO NSA=1,NSAMAX
            NS=NS_NSA(NSA)
            NSBA=NSB_NSA(NSA)
            DO NP=1,NPMAX
               CALL FPMXWL_EDGE(NP,NSA,FL)
!               FL=FPMXWL(PM(NP,NSBA),NRMAX+1,NS)
               DO NTH=1,NTHMAX
                  FS2(NTH,NP,NSA)=FL ! at R=1.0+DELR/2
               ENDDO
            ENDDO
            IF(MODELA.eq.1)THEN
               RSUM1=0.D0
               RSUM2=0.D0
               RSUM3=0.D0
               RSUM4=0.D0
               DO NP=1,NPMAX
                  DO NTH=1,NTHMAX
                     RSUM1 = RSUM1+VOLP(NTH,NP,NSBA)&!*FS2(NTH,NP,NSA) &
                                  *RLAMDAG(NTH,NRMAX+1)
                     RSUM2 = RSUM2+VOLP(NTH,NP,NSBA)!*FS2(NTH,NP,NSA)
                     RSUM3 = RSUM3+VOLP(NTH,NP,NSBA)*RLAMDA_GG(NTH,NRMAX+1)
                     RSUM4 = RSUM4+VOLP(NTH,NP,NSBA)
                  END DO
               END DO
               RCOEF2(1)=RSUM2/RSUM1
               RCOEF2_G(1)=RSUM4/RSUM3
            ELSE
               RCOEF2(1)=1.D0
               RCOEF2_G(1)=1.D0
            END IF

         ENDDO
!      END IF

      allocate(work(nrstart:nrendx),workg(NRMAX))

      DO NR=NRSTART,NRENDX
         work(NR)=RCOEF(NR)
      ENDDO
      CALL mtx_gatherv_real8(work(NRSTART:NRENDX),MTXLEN(NRANK+1), &
           workg,NRMAX,MTXLEN,MTXPOS)
      CALL mtx_broadcast_real8(workg,NRMAX)
      DO NR=1,NRMAX
         RCOEFG(NR)=workg(NR)
      ENDDO
      RCOEFG(NRMAX+1)=RCOEF2(1)
!      RCOEFG(NRMAX+1)=RCOEFG(NRMAX)

      DO NR=NRSTART,NRENDX
         work(NR)=RCOEF_G(NR)
      ENDDO
      CALL mtx_gatherv_real8(work(NRSTART:NRENDX),MTXLEN(NRANK+1), &
           workg,NRMAX,MTXLEN,MTXPOS)
      CALL mtx_broadcast_real8(workg,NRMAX)
      DO NR=1,NRMAX
         RCOEF_GG(NR)=workg(NR)
      ENDDO
      RCOEF_GG(NRMAX+1)=RCOEF2_G(1)

      deallocate(work,workg)

      DO NTH=1,NTHMAX
         DO NR=1,NRMAX+1
            RLAMDAG(NTH,NR)=RLAMDAG(NTH,NR)*RCOEFG(NR)
         END DO
         DO NR=NRSTART, NREND
            RLAMDA(NTH,NR)=RLAMDA(NTH,NR)*RCOEFG(NR)
         END DO

         DO NR=1,NRMAX+1
            RLAMDA_GG(NTH,NR)=RLAMDA_GG(NTH,NR)*RCOEF_GG(NR)
         END DO
         DO NR=NRSTART, NREND
            RLAMDA_G(NTH,NR)=RLAMDA_G(NTH,NR)*RCOEF_GG(NR)
         END DO
      ENDDO

!      IF(NRANK.eq.1)THEN
!      DO NR=1,NRMAX+1
!         RSUM4=0.D0
!         DO NP=1,NPMAX
!            DO NTH=1,NTHMAX
!               RSUM4 = rsum4+VOLP(NTH,NP,1)*RLAMDAG(NTH,NR)
!            END DO
!         END DO
!         WRITE(*,'(I3,1P4E14.6)') NR,RM(NR),RLAMDAG(1,NR),RCOEFG(NR),RSUM4
!      END DO
!      write(*,*) " "
!      write(*,*) " "
!         RSUM4=0.D0
!         DO NP=1,NPMAX
!            DO NTH=1,NTHMAX
!               RSUM4 = rsum4+VOLP(NTH,NP,1)*RLAMDA_GG(NTH,NR)
!            END DO
!         END DO
!      DO NR=1,NRMAX+1
!         WRITE(*,'(I3,1P4E14.6)') NR,RG(NR),RLAMDA_GG(1,NR),RCOEF_GG(NR),rsum4
!      END DO
!      END IF

!     ----- set parameters for target species -----

      DO NSA=1,NSAMAX
         NS=NS_NSA(NSA)
         AEFP(NSA)=PZ(NS)*AEE
         AMFP(NSA)=PA(NS)*AMP
         RNFP0(NSA)=PN(NS)
         RNFPS(NSA)=PNS(NS)
         RTFP0(NSA)=(PTPR(NS)+2.D0*PTPP(NS))/3.D0
         RTFPS(NSA)=PTS(NS)

         PTFP0(NSA)=SQRT(RTFP0(NSA)*1.D3*AEE*AMFP(NSA))
         VTFP0(NSA)=SQRT(RTFP0(NSA)*1.D3*AEE/AMFP(NSA))
      ENDDO

      DO NSB=1,NSBMAX
         NS=NS_NSB(NSB)
         AEFD(NSB)=PZ(NS)*AEE
         AMFD(NSB)=PA(NS)*AMP
         RTFD0L=(PTPR(NS)+2.D0*PTPP(NS))/3.D0
         RNFD0(NSB)=PN(NS)

         PTFD0(NSB)=SQRT(RTFD0L*1.D3*AEE*AMFD(NSB))
         VTFD0(NSB)=SQRT(RTFD0L*1.D3*AEE/AMFD(NSB))
      ENDDO

!     ----- set profile data -----

      DO NR=NRSTART,NREND+1

         RHON=RM(NR)
         CALL PL_PROF(RHON,PLF)

         DO NSA=1,NSAMAX
            NS=NS_NSA(NSA)
            RNFP(NR,NSA)=PLF(NS)%RN
            RTFP(NR,NSA)=(PLF(NS)%RTPR+2.D0*PLF(NS)%RTPP)/3.D0
            PTFP(NR,NSA)=SQRT(RTFP(NR,NSA)*1.D3*AEE*AMFP(NSA))
            VTFP(NR,NSA)=SQRT(RTFP(NR,NSA)*1.D3*AEE/AMFP(NSA))
         ENDDO

         DO NSB=1,NSBMAX
            NS=NS_NSB(NSB)
            AEFD(NSB)=PZ(NS)*AEE
            AMFD(NSB)=PA(NS)*AMP
            RNFD(NR,NSB)=PLF(NS)%RN
            RTFD(NR,NSB)=(PLF(NS)%RTPR+2.D0*PLF(NS)%RTPP)/3.D0
            PTFD(NR,NSB)=SQRT(RTFD(NR,NSB)*1.D3*AEE*AMFD(NSB))
            VTFD(NR,NSB)=SQRT(RTFD(NR,NSB)*1.D3*AEE/AMFD(NSB))
         ENDDO

         RNE=PLF(1)%RN
         RTE=(PLF(1)%RTPR+2.D0*PLF(1)%RTPP)/3.D0

         DO NSA=1,NSAMAX
            NSFP=NS_NSB(NSA)
            DO NSB=1,NSBMAX
               NSFD=NS_NSB(NSB)

               IF(NSFP.EQ.1.AND.NSFD.EQ.1) THEN
                  RLNRL=14.9D0-0.5D0*LOG(RNE)+LOG(RTE)
               ELSEIF(NSFP.EQ.1.OR.NSFD.EQ.1) THEN
                  RLNRL=15.2D0-0.5D0*LOG(RNE)+LOG(RTE)
               ELSE
                  RLNRL=17.3D0-0.5D0*LOG(RNE)+1.5D0*LOG(RTFD(NR,NSB))
               ENDIF
               FACT=AEFP(NSA)**2*AEFD(NSB)**2*RLNRL/(4.D0*PI*EPS0**2)
               RNUD(NR,NSB,NSA)=FACT*RNFP0(NSA)*1.D20 &
                      /(SQRT(2.D0)*VTFD(NR,NSB)*PTFP0(NSA)**2)
               RNUF(NR,NSB,NSA)=FACT*RNFP0(NSA)*1.D20 &
                      /(2*AMFD(NSB)*VTFD(NR,NSB)**2*PTFP0(NSA))
            ENDDO
         ENDDO
      ENDDO

!     ----- set relativistic parameters -----

      IF (MODELR.EQ.0) THEN
         DO NSA=1,NSAMAX
            THETA0(NSA)=0.D0
            DO NR=NRSTART,NREND
               THETA(NR,NSA)=0.D0
            ENDDO
         ENDDO
      ELSE
         DO NSA=1,NSAMAX
            THETA0(NSA)=RTFP0(NSA)*1.D3*AEE/(AMFP(NSA)*VC*VC)
            DO NR=NRSTART,NREND
               THETA(NR,NSA)=THETA0(NSA)*RTFP(NR,NSA)/RTFP0(NSA)
            ENDDO
         ENDDO
      ENDIF

!$$$      DO NSA=1,NSAMAX
!$$$         NS=NS_NSA(NSA)
!$$$      DO NR=NRSTART,NREND
!$$$         NTH=1
!$$$         DO NP=2,NPMAX-1
!$$$            IF(FNS(NTH,NP+1,NR,NSA).LE.1.D-70) THEN
!$$$               DFDP=0.D0
!$$$            ELSE
!$$$               DFDP=(FNS(NTH,NP+1,NR,NS)-FNS(NTH,NP-1,NR,NS))
!$$$     &                    /(2.D0*DELP*FNS(NTH,NP,NR,NS))
!$$$            ENDIF
!$$$            IF(FNS(NTH,NP,NR,NSA).LE.1.D-70) THEN
!$$$               DFDP1=0.D0
!$$$            ELSE
!$$$               DFDP1=(LOG(FNS(NTH,NP+1,NR,NS))-LOG(FNS(NTH,NP-1,NR,NS)))
!$$$     &                    /(2.D0*DELP)
!$$$            ENDIF
!$$$            DFDP2=-PM(NP)*RTFP0(NSA)/RTFP(NR,NSA)
!$$$            write(6,'(I5,1P5E12.4)')
!$$$     &           NP,PM(NP),DFDP,DFDP1,DFDP2,FNS(NTH,NP,NR,NS)
!$$$         ENDDO
!$$$      ENDDO
!$$$      ENDDO
      N_IMPL=0
      CALL NF_REACTION_COEF
!      IF(nprocs.gt.1.and.NRANK.eq.1)THEN 
!      open(8,file='rlamdag_new.dat')
!         DO NR=1,NRMAX
!         DO NTH=1,NTHMAX
!            WRITE(8,'(2I5,1PE12.4)') NR,NTH,RLAMDAG(NTH,NR)
!         END DO
!         WRITE(8,*)" "
!         WRITE(8,*)" "
!         END DO
!      CLOSE(8)
!      END IF
!      NCALCNR=0
      DO NSA=1,NSAMAX
         CALL FP_COEF(NSA)
         NSBA=NSB_NSA(NSA)
         DO NTH=1,NTHMAX
            DO NP=1,NPMAX
!               DO NR=1,NRMAX
               DO NR=NRSTART,NREND
                  F(NTH,NP,NR)=FNS(NTH,NP,NR,NSBA)
               END DO
            END DO
         END DO
         CALL FPWEIGHT(NSA,IERR)
!         IF(MODELR.eq.1.and.MODELC.eq.4.and.NCALCNR.eq.2)THEN
!            NCALCNR=1
!            CALL FP_COEF(NSA)
!         END IF
      END DO
      ISAVE=0
      IF(NTG1.eq.0) CALL FPWAVE_CONST ! all nrank must have RPWT  
      CALL FPSSUB
      IF(nrank.EQ.0) THEN
         CALL FPSGLB
         CALL FPWRTGLB
         CALL FPSPRF
         CALL FPWRTPRF
      ENDIF
      IERR=0
      RETURN
      END subroutine fp_prep

      END MODULE fpprep
