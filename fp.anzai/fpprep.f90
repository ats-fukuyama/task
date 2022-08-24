!     $Id: fpprep.f90,v 1.41 2013/02/08 07:36:24 nuga Exp $

! *****************************
!     PREPARATION OF FPLOOP
! *****************************
      MODULE fpprep

      USE fpcomm
      USE fpinit
      USE fpsave
      USE fpcoef
      USE fpcalw
      USE fpbounce
      USE equnit
      USE fpmpi
      USE libmpi
      USE fpcaleind
      USE fpdisrupt

      contains

!     ***** create mesh quantities *****

      SUBROUTINE fp_mesh(IERR)

      USE libmtx
      USE plprof
      USE fpbroadcast
      USE fpwrin
      USE fpwmin

      Implicit none
      integer :: ierr,NSA,NS,NR,NP,NTH,id
!      character(LEN=80)::line 
      real(kind8)::rhon,rhol,rhol1,rhol2,A1,epsl,ql,BT
      real(kind8),DIMENSION(:),POINTER:: work,workg
      real(kind8):: Rmass, RRTFP, RPTFP,RVTFP, sumEmax

!     ----- define upper boundary of p from Emax-----
      sumEmax=0.D0
      DO NS=1,NSMAX
         sumEmax=sumEmax+EMAX(ns)
      END DO
      IF(sumEmax.ne.0.D0)THEN ! IF Emax is given by a input file --
         IF(NRANK.eq.0) WRITE(6,'(A,1P6E12.4)') "old pmax= ", (pmax(ns),ns=1,nsmax)
         DO NS=1, NSMAX
            IF(EMAX(NS).ne.0.D0)THEN
               Rmass=PA(NS)*AMP
               RRTFP=(PTPR(NS)+2.D0*PTPP(NS))/3.D0
               RPTFP=SQRT(RRTFP*1.D3*AEE*Rmass)
               RVTFP=SQRT(RRTFP*1.D3*AEE/Rmass)
               pmax(ns)=SQRT(2.D0*AEE*Rmass*1.D3*EMAX(NS))/RPTFP
            END IF
         END DO
         IF(NRANK.eq.0) WRITE(6,'(A,1P6E12.4)') "new pmax= ", (pmax(ns),ns=1,nsmax)
      END IF
!      IF(NRANK.eq.0) WRITE(6,'(A,1P6E12.4)') "new vmax[m/s] ", (pmax(nsb)*RVTFP,nsb=1,nsbmax)

!     ----- exec EQ -----

      IF(MODELG.EQ.3.OR.MODELG.EQ.5.OR.MODELG.EQ.8) THEN
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
      WRITE(6,*) 'BB  =',BB
      WRITE(6,*) 'RR  =',RR
      WRITE(6,*) 'RIP =',RIP
      WRITE(6,*) 'RA  =',RA
      RKAP=1.D0

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
         QLM(NR)=QL
         BT=BB
         BP= RSRHON(RHON)*BT/(RR*QL)
         EPSRM(NR)=RSRHON(RHON)/RR
         BPM(NR)= RSRHON(RHON)*BT/(RR*QL)
!         IF(NRANK.eq.0) &
!         write(6,'(A,I5,1P5E12.4)') 'nr,rm,rsrhon,epsrm,bpm,ql=', &
!              NR,RM(NR),RSRHON(RHON),EPSRM(NR),BPM(NR),QL
      ENDDO
!      RHON=RG(NRMAX+1)
      RHON=RM(NRMAX)+DELR
      CALL pl_qprf(RHON,QL)
      QLM(NRMAX+1)=QL
      BT=BB
      BP= RSRHON(RHON)*BT/(RR*QL)
      EPSRM(NRMAX+1)=RSRHON(RHON)/RR
      BPM(NRMAX+1)= RSRHON(RHON)*BT/(RR*QL)
      IF(NRANK.eq.0) &
           write(6,'(A,I5,1P5E12.4)') 'nr,rm,rsrhon,epsrm,bpm,ql=', &
           NRMAX+1,RM(NRMAX)+DELR,RSRHON(RHON),EPSRM(NRMAX+1),BPM(NRMAX+1),QLM(NRMAX+1)

!      IF(NRANK.eq.0) WRITE(*,*) "BP=", BP
      DO NR=1,NRMAX+1
         RHON=RG(NR)
         CALL pl_qprf(RHON,QL)
         QLG(NR)=QL
         BT=BB
         BP(NR)= RSRHON(RHON)*BT/(RR*QL)
         EPSRG(NR)=RSRHON(RHON)/RR
         BPG(NR)= RSRHON(RHON)*BT/(RR*QL)
      ENDDO
!     ----- set parallel current density -----

      DO NR=1,NRMAX
         RJ1(NR)=(RG(NR+1)*BP(NR+1)-RM(NR)*BP(NR)) &
               /(RMU0*RM(NR)*DELR)
      ENDDO

!     ----- set momentum space mesh -----

      DELTH=PI/NTHMAX
      DO NS=1,NSMAX
         DELP(NS)=PMAX(NS)/NPMAX
         DO NP=1,NPMAX
            PG(NP,NS)=DELP(NS)*(NP-1)
            PM(NP,NS)=DELP(NS)*(NP-0.5D0)
         ENDDO
         PG(NPMAX+1,NS)=PMAX(NS)
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

!     set momentum volume element: VOLP
      DO NS=1,NSMAX
      DO NP=NPSTART,NPEND
      DO NTH=1,NTHMAX
         VOLP(NTH,NP,NS)=2.D0*PI*SINM(NTH)*PM(NP,NS)**2*DELP(NS)*DELTH
      ENDDO
      ENDDO
      ENDDO

!     set volume element: VOLR= int dchi*dzeta (Jacob*dr)
      SELECT CASE(MODELG)
      CASE(0:2) ! cylinder
         DO NR=1,NRMAX
            RHOL=RM(NR)
            RHOL1=RG(NR)
            RHOL2=RG(NR+1)
            VOLR(NR)=2.D0*PI*RSRHON(RHOL)*(RSRHON(RHOL2)-RSRHON(RHOL1)) &
                 *2.D0*PI*RR
         ENDDO
      CASE(3,5,8) ! toroidal
         DO NR=1,NRMAX
            RHOL=RM(NR)
            RHOL1=RG(NR)
            RHOL2=RG(NR+1)
            VOLR(NR)=2.D0*PI*RSRHON(RHOL)*(RSRHON(RHOL2)-RSRHON(RHOL1)) &
                 *2.D0*PI*RR
         ENDDO
      END SELECT

      TVOLR=0.D0
      DO NR=1,NRMAX
         TVOLR=TVOLR+VOLR(NR)
      ENDDO

      IF(NRANK.eq.0) THEN
         WRITE(6,'(A,3E14.6)') "DEVICE, RR, RA, BB", RR, RA, BB
      END IF
!     ----- set bounce-average parameters -----

      IF (MODELA.EQ.0) THEN
         DO NR=1,NRMAX+1
            ITL(NR)=0
            ITU(NR)=0
            ITLG(NR)=0
            ITUG(NR)=0
            EPSRM2(NR) = EPSRM(NR)
            EPSRG2(NR) = EPSRG(NR)
         ENDDO
      ELSE
! THM(ITL(NR)) is always trapped region 
! Exact boundary A1 lays on between THM(ITL-1) and THM(ITL)
         IF(nsize.gt.1.and.NRANK.eq.1) & 
              WRITE(6,'(A)') '# NR,ITL(NR),EPSRM2,EPSRM,THM(ITL-1)<THC<THM(ITL),THG(ITL) on RM(NR)'
         DO NR=1,NRMAX+1
            A1=ACOS(SQRT(2.D0*EPSRM(NR)/(1.D0+EPSRM(NR))))
            NTH=0
!            DO WHILE (THG(NTH+1).le.A1)
            DO WHILE (THM(NTH+1).le.A1)
               NTH = NTH+1
            END DO
            ITL(NR)=NTH+1
!            ITL(NR)=NTH
            ITU(NR)=NTHMAX-ITL(NR)+1

            EPSL=COSM(ITL(NR))**2/(2.D0-COSM(ITL(NR))**2)

            EPSRM2(NR) = EPSRM(NR)
            EPSRM(NR)=EPSL
            IF(nsize.gt.1.and.NRANK.eq.1) &
                 WRITE(6,'(2I5,1P8E12.4)') NR,ITL(NR),EPSRM2(NR),EPSRM(NR),THM(ITL(NR)-1), A1, THM(ITL(NR)),THG(ITL(NR))

            IF(A1.ge.THG(ITL(NR)))THEN ! The mesh point including A1. RLAMDA(NTH) on such a point can not calculate numerically.
               ITL_judge(NR)=ITL(NR)
!               ITL_judge(NR)=ITL(NR)-1
            ELSE
               ITL_judge(NR)=ITL(NR)-1
            END IF

         ENDDO

         IF(NRANK.eq.1) WRITE(6,*) " "
         IF(nsize.gt.1.and.NRANK.eq.1) &
              WRITE(6,'(A)') '# NR,ITLG(NR),EPSRG2,EPSRG,THG(ITL),THC on RG(NR)'

         DO NR=1,NRMAX+1
            A1=ACOS(SQRT(2.D0*EPSRG(NR)/(1.D0+EPSRG(NR))))
            NTH=0
!            DO WHILE (THG(NTH+1).le.A1)
            DO WHILE (THM(NTH+1).le.A1)
               NTH = NTH+1
            END DO
            ITLG(NR)=NTH+1
            ITUG(NR)=NTHMAX-ITLG(NR)+1

            EPSL=COSM(ITLG(NR))**2/(2.D0-COSM(ITLG(NR))**2)

            EPSRG2(NR) = EPSRG(NR)
            EPSRG(NR)=EPSL
            IF(nsize.gt.1.and.NRANK.eq.1) &
                 WRITE(6,'(2I5,1P4E12.4)') NR,NTH,EPSRG2(NR),EPSRG(NR),THG(NTH), A1

            IF(A1.ge.THG(ITLG(NR)))THEN ! The mesh point including A1. RLAMDA(NTH) on such a point can not calculate numerically.
               ITLG_judge(NR)=ITLG(NR)
            ELSE
               ITLG_judge(NR)=ITLG(NR)-1
            END IF
         ENDDO

      ENDIF

      IF (MODELA.EQ.0) THEN
         DO NR=NRSTART,NREND
            DO NTH=1,NTHMAX
               ETAM(NTH,NR)=PI*0.5D0
               RLAMDA(NTH,NR)=1.D0
               RLAMDC(NTH,NR)=1.D0
               RLAMDC_G(NTH,NR)=1.D0
            ENDDO
            DO NTH=1,NTHMAX+1
               ETAG_RG(NTH,NR)=PI/2.D0
            ENDDO
         ENDDO
!
         DO NR=NRSTART,NRENDWG
            DO NTH=1,NTHMAX
               ETAM_RG(NTH,NR)=PI*0.5D0
               RLAMDA_RG(NTH,NR)=1.D0
            END DO
            DO NTH=1,NTHMAX+1
               ETAG_RG(NTH,NR)=PI/2.D0
            ENDDO
         END DO

         DO NTH=1, NTHMAX
            RLAMDAG_RG(NTH,NRMAX+1)=1.D0
            RLAMDA_NRMAXP1(NTH)=1.D0
         END DO
         DO NR=1,NRMAX+1
            RFSADG(NR)=1.D0
            RFSADG_RG(NR)=1.D0
         END DO
      ELSE
         RLAMDA_NRMAXP1(:)=0.D0
         DO NR=1,NRMAX+1
            CALL SET_RFSAD(NR)
         END DO
         CALL SET_BOUNCE_PARAM

      END IF ! MODELA

      allocate(work(nrstart:nrend),workg(NRMAX))

      CALL mtx_set_communicator(comm_nr)
      DO NTH=1,NTHMAX
         DO NR=NRSTART,NREND
            work(NR)=ETAM(NTH,NR)
         ENDDO
         CALL mtx_allgatherv_real8(work(NRSTART:NREND),MTXLEN(NRANK+1), &
                                   workg,NRMAX,MTXLEN,MTXPOS)
         DO NR=1,NRMAX
            ETAMG(NTH,NR)=workg(NR)
         ENDDO
      ENDDO

      DO NTH=1,NTHMAX
         DO NR=NRSTART,NREND
            work(NR)=RLAMDA(NTH,NR)
         ENDDO
         CALL mtx_allgatherv_real8(work(NRSTART:NREND),MTXLEN(NRANK+1), &
                                   workg,NRMAX,MTXLEN,MTXPOS)
         DO NR=1,NRMAX
            RLAMDAG(NTH,NR)=workg(NR)
         ENDDO
         RLAMDAG(NTH,NRMAX+1)=RLAMDA_NRMAXP1(NTH)
      ENDDO
      
      DO NTH=1,NTHMAX
         DO NR=NRSTART,NREND
            work(NR)=ETAM_RG(NTH,NR)
         ENDDO
         CALL mtx_allgatherv_real8(work(NRSTART:NREND),MTXLEN(NRANK+1), &
                                   workg,NRMAX,MTXLEN,MTXPOS)
         DO NR=1,NRMAX
            ETAMG_RG(NTH,NR)=workg(NR)
         ENDDO
      ENDDO

      DO NTH=1,NTHMAX
         DO NR=NRSTART,NREND
            work(NR)=RLAMDA_RG(NTH,NR)
         ENDDO
         CALL mtx_allgatherv_real8(work(NRSTART:NREND),MTXLEN(NRANK+1), &
                                   workg,NRMAX,MTXLEN,MTXPOS)
         DO NR=1,NRMAX
            RLAMDAG_RG(NTH,NR)=workg(NR)
         ENDDO
      ENDDO

      CALL mtx_reset_communicator

!      IF(MODELA.eq.1)THEN
!      IF(NRANK.eq.0)THEN
!      open(8,file='RLAMDAG.dat')
!      DO NR =1, NRMAX+1
!         DO NTH=1,NTHMAX
!            WRITE(8,'(2I4, 4E14.6)') NR, NTH, COSM(NTH), RLAMDAG(NTH,NR), RFSADG(NR)
!         END DO
!         WRITE(8,*) " "
!         WRITE(8,*) " "
!      END DO
!      close(8)
!      open(8,file='RLAMDAG_RG.dat')
!      DO NR =1, NRMAX+1
!         DO NTH=1,NTHMAX
!            WRITE(8,'(2I4, 4E14.6)') NR, NTH, COSM(NTH), RLAMDAG_RG(NTH,NR), RFSADG_RG(NR)
!         END DO
!         WRITE(8,*) " "
!         WRITE(8,*) " "
!      END DO
!      close(8)
!      END IF
!      END IF

!      IF(NREND.eq.NRMAX)THEN
!         open(8,file='RLAMDA_RG_MAX.dat')
!         DO NTH=1, NTHMAX
!            WRITE(8,'(2I4, 4E14.6)') NRMAX+1, NTH, COSM(NTH), RLAMDA_RG(NTH,NRMAX+1), RFSADG_RG(NRMAX+1)
!         END DO
!         CLOSE(8)
!      END IF

!      IF(NRANK.eq.0)THEN
!      DO NR=1,NRMAX
!         WRITE(*,*) NR, "NRG", ITLG(NR)
!         WRITE(*,*) NR, "NRM ", ITL(NR)
!      END DO
!      END IF

      deallocate(work,workg)
      IERR=0
      RETURN
      END subroutine fp_mesh

! *************************
!     INITIAL DATA SAVE
! *************************

      SUBROUTINE FPCINI

      IMPLICIT NONE
      integer:: NSA, NR, NTH, NP

      ISAVE=0

      DO NSA=1,NSAMAX
      DO NR=NRSTART,NREND
         DO NP=NPSTART,NPENDWG
         DO NTH=1,NTHMAX
            DPP(NTH,NP,NR,NSA)=0.D0
            DPT(NTH,NP,NR,NSA)=0.D0
            FPP(NTH,NP,NR,NSA)=0.D0

            DCPP(NTH,NP,NR,NSA)=0.D0
            DCPT(NTH,NP,NR,NSA)=0.D0
            FCPP(NTH,NP,NR,NSA)=0.D0

            FEPP(NTH,NP,NR,NSA)=0.D0
            FSPP(NTH,NP,NR,NSA)=0.D0

            DWPP(NTH,NP,NR,NSA)=0.D0
            DWPT(NTH,NP,NR,NSA)=0.D0

            DLPP(NTH,NP,NR,NSA)=0.D0
            FLPP(NTH,NP,NR,NSA)=0.D0
         ENDDO
         ENDDO

         DO NP=NPSTARTW,NPENDWM
         DO NTH=1,NTHMAX+1
            DTP(NTH,NP,NR,NSA)=0.D0
            DTT(NTH,NP,NR,NSA)=0.D0
            FTH(NTH,NP,NR,NSA)=0.D0

            DCTP(NTH,NP,NR,NSA)=0.D0
            DCTT(NTH,NP,NR,NSA)=0.D0
            FCTH(NTH,NP,NR,NSA)=0.D0

            FETH(NTH,NP,NR,NSA)=0.D0
            FSTH(NTH,NP,NR,NSA)=0.D0

            DWTP(NTH,NP,NR,NSA)=0.D0
            DWTT(NTH,NP,NR,NSA)=0.D0
         ENDDO
         ENDDO
      ENDDO
      ENDDO

      RETURN
      END SUBROUTINE FPCINI

! ============================================================
      SUBROUTINE fp_comm_setup

      USE libmtx
      IMPLICIT NONE
      INTEGER:: ierr, NREND1, keys,N
      INTEGER,DIMENSION(nsize):: ima1,ima2,npa1,npa2,nra1,nra2,nma1,nma2,insa1,insa2

!     ----- Check nsize -----
      IF(N_partition_s*N_partition_r*N_partition_p.ne.nsize)THEN
         IF(NRANK.eq.0) THEN
            WRITE(6,*) 'XX fp_prep: N_partition_s*N_partition_r != nsize.'
            WRITE(6,'(A,3I4)') 'XX N_partition_s, N_partition_r, nsize=', &
                                   N_partition_s, N_partition_r, nsize
         END IF
         ierr=1
         RETURN
      END IF

      IF(MOD(NPMAX,N_partition_p).ne.0)THEN
         IF(NRANK.eq.0) THEN
            WRITE(6,*) 'XX fp_prep: MOD(NPMAX,N_partition_p)!=0'
            WRITE(6,'(A,2I4)') 'XX NPMAX, N_partition_p', &
                                   NPMAX, N_partition_p
         END IF
         ierr=1
         RETURN
      END IF

      IF(MOD(NRMAX,N_partition_r).ne.0)THEN
         IF(NRANK.eq.0) THEN
            WRITE(6,*) 'XX fp_prep: MOD(NRMAX,N_partition_r)!=0'
            WRITE(6,'(A,2I4)') 'XX NRMAX, N_partition_r', &
                                   NRMAX, N_partition_r
         END IF
         ierr=1
         RETURN
      END IF

!     ----- Set commpunicator -----
!      CALL mtx_comm_split2D(N_partition_s,N_partition_r,comm_nr,comm_nsa) !2D
      CALL mtx_comm_split3D(N_partition_s,N_partition_r,N_partition_p, &
           comm_nsa,comm_nr,comm_np,comm_nrnp,comm_nsanp,comm_nsanr)
      IF(N_partition_p.ne.1)THEN
         IF(NRMAX.ne.N_partition_r)THEN
            WRITE(*,*) "IF N_partition_p does not equal 1, N_partition_r must be NRMAX"
            STOP
         END IF
      END IF

!      nranks=comm_nr%rank
!      nsizes=comm_nr%size
!      colors=comm_nr%rankg
!      keys=  comm_nr%rankl
!      NSASTART = (NSAMAX/N_partition_s)*colors+1 !2D
!      NSAEND =   (NSAMAX/N_partition_s)*(colors+1)


      IF(N_partition_s.GT.NSAMAX) THEN
         IF(NRANK.EQ.0) &
              write(6,'(A,2I5)') 'XX N_partition_s.GT.NSAMAX:', &
              N_partition_s,NSAMAX
         STOP
      END IF
      keys=comm_nsa%rankl
      NSASTART = (NSAMAX/N_partition_s)*keys+1 !3D
      NSAEND =   (NSAMAX/N_partition_s)*(keys+1)

      imtxsize=nthmax*npmax*nrmax
      IF(modeld.EQ.0) THEN
         imtxwidth=4*nthmax-1
      ELSE
         imtxwidth=4*nthmax*npmax-1
      ENDIF
      imtxwidth=npmax*nthmax*nrmax*2-1
!      CALL mtx_set_communicator(comm_nr)! 2D
      CALL mtx_set_communicator(comm_nrnp)! 3D
      CALL mtx_setup(imtxsize,imtxstart,imtxend,imtxwidth)

      nrstart=(imtxstart-1)/(nthmax*npmax)+1
      nrend=  (imtxend  -1)/(nthmax*npmax)+1
      nrend1= (imtxend    )/(nthmax*npmax)+1
!      IF(nrend1.EQ.nrend) THEN
!         NRENDX=NREND-1
!      ELSE
!         NRENDX=NREND
!      ENDIF

      NPSTART=( imtxstart-1- (nrstart-1)*nthmax*npmax )/nthmax +1
      NPEND  =( imtxend - (nrend-1)*nthmax*npmax )/nthmax

!---- SET OF SHADOWS
      IF(NRSTART.eq.1)THEN
         NRSTARTW=1
      ELSE
         NRSTARTW=NRSTART-1
      END IF
      IF(NREND.eq.NRMAX)THEN
         NRENDWM=NRMAX
         NRENDWG=NRMAX+1
      ELSE
         NRENDWM=NREND+1
         NRENDWG=NREND+1
      END IF

      IF(NPSTART.eq.1)THEN ! SET SHADOW OF NP
         NPSTARTW=1
      ELSE
         NPSTARTW=NPSTART-1
      END IF
      IF(NPEND.eq.NPMAX)THEN
         NPENDWM=NPMAX   ! for PM(NP)
         NPENDWG=NPMAX+1 ! for PG(NP)
      ELSE
         NPENDWM=NPEND+1
         NPENDWG=NPEND+1
      END IF 
!---- OND OF SHADOW

!      nmstart=nthmax*npmax*(nrstart-1)+1 !2D
!      nmend  =nthmax*npmax* nrend
      nmstart= imtxstart !3D
      nmend  = imtxend

      CALL mtx_reset_communicator

      CALL mtx_gather1_integer(imtxstart,ima1)
      CALL mtx_gather1_integer(imtxend,  ima2)
      CALL mtx_gather1_integer(npstart,  npa1)
      CALL mtx_gather1_integer(npend,    npa2)
      CALL mtx_gather1_integer(nrstart,  nra1)
      CALL mtx_gather1_integer(nrend,    nra2)
      CALL mtx_gather1_integer(nmstart,  nma1)
      CALL mtx_gather1_integer(nmend,    nma2)
      CALL mtx_gather1_integer(nsastart,insa1)
      CALL mtx_gather1_integer(nsaend,  insa2)

      IF(nrank.EQ.0) THEN
         write(6,'(A,2I10)') '  imtxsize,imtxwidth=',imtxsize,imtxwidth
         write(6,'(A,A,A)') '  nrank imtxstart end  npstart  end', &
                                   ' nrstart   end  nmstart  end', &
                                   ' nsastart  end'
         DO N=1,nsize
            write(6,'(11I7)') N,ima1(N),ima2(N),npa1(N),npa2(N), &
                                nra1(N),nra2(N),nma1(N),nma2(N), &
                                insa1(N),insa2(N)
         ENDDO
      ENDIF

!      IF(NRANK.eq.0) WRITE(6,*) "      NRANK, nsa_rank, nr_rank, np_rank, nrnp_rank, nsanp_rank, nsanr_rank"
!      write(6,'(7I10)') NRANK, comm_nsa%rankl, comm_nr%rankl, &
!           comm_np%rankl, comm_nrnp%rankl, comm_nsanp%rankl, comm_nsanr%rankl

!      IF(NRANK.eq.0) WRITE(6,*) "RANK, NSAS, NSAE,  NRS,  NRE,  NPS,  NPE"
!      WRITE(6,'(7I6)') NRANK, NSASTART, NSAEND, NRSTART, NREND, NPSTART, NPEND

      CALL mtx_cleanup

      CALL mtx_reset_communicator

      END SUBROUTINE fp_comm_setup
!-------------------------------------------------------------
      SUBROUTINE fp_set_nsa_nsb

      IMPLICIT NONE
      INTEGER:: NSA, NSB, IERR

      if(nrank.eq.0) then
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
      ELSE
         DO NSA=1,NSAMAX
            IF(NS_NSA(NSA).EQ.0) THEN
               NS_NSA(NSA)=NSA
            ENDIF
         ENDDO
         DO NSB=1,NSBMAX
            IF(NS_NSB(NSB).EQ.0) THEN
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


      END SUBROUTINE fp_set_nsa_nsb
!-------------------------------------------------------------
      SUBROUTINE FNSP_INIT

        USE fpsub
      IMPLICIT NONE
      INTEGER:: NTH,NP,NR,NSA,NS,NSB
      REAL(8):: FL

      DO NSA=NSASTART,NSAEND
         NS=NS_NSA(NSA)
         DO NR=NRSTARTW,NRENDWM
            IF(NR.ge.1.and.NR.le.NRMAX)THEN
               DO NP=NPSTARTW,NPENDWM
                  FL=FPMXWL(PM(NP,NS),NR,NS)
!                  IF(NRANK.EQ.0) WRITE(6,'(A,4I5,1PE12.4)') &
!                       'NSA,NS,NR,NP,FL=',NSA,NS,NR,NP,FL
                  DO NTH=1,NTHMAX
                     FNSP(NTH,NP,NR,NSA)=FL
                     FNS0(NTH,NP,NR,NSA)=FL
                  END DO
               ENDDO
            END IF
         END DO
      END DO

      DO NSA=NSASTART,NSAEND
         NS=NS_NSA(NSA)
         IF(MODEL_DELTA_F(NS).eq.1)THEN
            DO NR=NRSTARTW,NRENDWM
               IF(NR.ge.1.and.NR.le.NRMAX)THEN
                  DO NP=NPSTARTW,NPENDWM
                     FL=FPMXWL(PM(NP,NS),NR,NS)
                     DO NTH=1,NTHMAX
                        FNSP_MXWL(NTH,NP,NR,NSA)=FL
                        FNSP_DEL(NTH,NP,NR,NSA)=FL*1.D-30
                     END DO
                  ENDDO
               END IF
            END DO
         END IF
      END DO

!      DO NS=1, NSMAX
!         DO NR=NRSTART,NREND
!            DO NP=NPSTART,NPEND
!               FL=FPMXWL(PM(NP,NS),NR,NS)
!               DO NTH=1,NTHMAX
!                  FNSB(NTH,NP,NR,NS)=FL
!               END DO
!            ENDDO
!         END DO
!      END DO

      CALL update_fnsb_maxwell

      END SUBROUTINE FNSP_INIT
!-------------------------------------------------------------
      SUBROUTINE FNSP_INIT_EDGE

      USE fpsub
      IMPLICIT NONE
      INTEGER:: NTH,NP,NR,NSA,NS
      REAL(8):: FL

      DO NSA=NSASTART,NSAEND
         NS=NS_NSA(NSA)
         NR=NRMAX+1
         DO NP=NPSTARTW,NPENDWM
            FL=FPMXWL(PM(NP,NS),NR,NS)
            DO NTH=1,NTHMAX
               FS1(NTH,NP,NSA)=FL ! rho=1.0 fixed value
            END DO
         ENDDO
      END DO
!     ----- set boundary distribution functions -----
      
      DO NSA=NSASTART,NSAEND
         NS=NS_NSA(NSA)
         DO NP=NPSTARTW,NPENDWM
            FL=FPMXWL(PM(NP,NS),0,NS)
            DO NTH=1,NTHMAX
               FS0(NTH,NP,NSA)=FL ! at r=0 fixed value
            ENDDO
         ENDDO
      END DO

      IF(MODELD_boundary.eq.0)THEN
         DO NSA=NSASTART,NSAEND
            NS=NS_NSA(NSA)
            DO NP=NPSTARTW,NPENDWM
               CALL FPMXWL_EDGE(NP,NSA,FL)
               DO NTH=1,NTHMAX
!                  FS2(NTH,NP,NS)=FL ! at R=1.0+DELR/2
                  FS2(NTH,NP,NSA)=FL ! at R=1.0+DELR/2
               ENDDO
            ENDDO
         ENDDO
      ELSEIF(MODELD_boundary.eq.1)THEN
         IF(NREND.eq.NRMAX)THEN
            DO NSA=NSASTART,NSAEND
               NS=NS_NSA(NSA)
               DO NP=NPSTARTW,NPENDWM
                  DO NTH=1,NTHMAX
                     FS2(NTH,NP,NSA) = 2.D0*FS1(NTH,NP,NSA)-FNSP(NTH,NP,NRMAX,NSA) ! linear
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            FS2(:,:,:)=0.D0
         END IF
      END IF

      END SUBROUTINE FNSP_INIT_EDGE
!-------------------------------------------------------------
      SUBROUTINE fp_set_normalize_param

      USE plprof
      USE EG_READ
      IMPLICIT NONE
      INTEGER:: NSA, NSB, NS, NSFP, NSFD, NR, ISW_CLOG, i, j
      TYPE(pl_prf_type),DIMENSION(NSMAX):: PLF
      real(kind8):: RTFD0L, RHON, RNE, RTE, RLNRL, FACT, RNA, RTA, RNB, RTB, SUM, AMFDL
      real(kind8):: A_D, tau_se_E0, k_energy, log_energy, sigma_cx0, sigma_cx, tau_cx_E1
      real(kind8):: tau_se_E0E1, k_energy1, log10_neu0, log10_neus, alpha, beta, N_NEUT, E_CR

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
      RT_E=RTFPS(1)*1.D-2
      RN_E=RNFPS(1)*1.D-1

      DO NSB=1,NSBMAX
         NS=NS_NSB(NSB)
         AEFD(NSB)=PZ(NS)*AEE
         AMFD(NSB)=PA(NS)*AMP
         RTFD0L=(PTPR(NS)+2.D0*PTPP(NS))/3.D0
         RNFD0(NSB)=PN(NS)
         RTFD0(NSB)=(PTPR(NS)+2.D0*PTPP(NS))/3.D0
         RTFDS(NSB)=PTS(NS)

         PTFD0(NSB)=SQRT(RTFD0L*1.D3*AEE*AMFD(NSB))
         VTFD0(NSB)=SQRT(RTFD0L*1.D3*AEE/AMFD(NSB))
      ENDDO

!     ----- set initial profile data -----
      DO NR=NRSTART, NRENDWG ! fixed initial vaule
         RHON=RG(NR)
         CALL PL_PROF(RHON,PLF)
         DO NSA=1, NSAMAX
            NS=NS_NSA(NSA)
            RNFP_G(NR,NSA)=PLF(NS)%RN
            RTFP_G(NR,NSA)=(PLF(NS)%RTPR+2.D0*PLF(NS)%RTPP)/3.D0
         END DO
      END DO

      DO NR=1,NRMAX ! change in time
         RHON=RM(NR)
         CALL PL_PROF(RHON,PLF)
         DO NS=1,NSMAX
            RT_TEMP(NR,NS)=(PLF(NS)%RTPR+2.D0*PLF(NS)%RTPP)/3.D0
            RN_TEMP(NR,NS)=PLF(NS)%RN
         END DO
      END DO

      IF(MODEL_EX_READ_Tn.ne.0)THEN
         CALL READ_EXP_DATA
         CALL MAKE_EXP_PROF(timefp)
         IF(NRANK.eq.0) WRITE(*,'(A,E14.6)') "time_exp_offset= ", time_exp_offset
      END IF
!      WRITE(*,*) NR, NS, RT_READ(NR,NS), RN_READ(NR,NS), RT_TEMP(NR,NS), RN_TEMP(NR,NS)

      DO NR=NRSTART,NRENDWM
         RHON=RM(NR)
         CALL PL_PROF(RHON,PLF)

         DO NSA=1,NSAMAX ! fixed initial value
            NS=NS_NSA(NSA)
            RNFP(NR,NSA)=PLF(NS)%RN
            RTFP(NR,NSA)=(PLF(NS)%RTPR+2.D0*PLF(NS)%RTPP)/3.D0
            PTFP(NR,NSA)=SQRT(RTFP(NR,NSA)*1.D3*AEE*AMFP(NSA))
            VTFP(NR,NSA)=SQRT(RTFP(NR,NSA)*1.D3*AEE/AMFP(NSA))
         ENDDO

         DO NSB=1,NSBMAX
            NS=NS_NSB(NSB)
            RNFD(NR,NSB)=PLF(NS)%RN
            RTFD(NR,NSB)=(PLF(NS)%RTPR+2.D0*PLF(NS)%RTPP)/3.D0
            PTFD(NR,NSB)=SQRT(RTFD(NR,NSB)*1.D3*AEE*AMFD(NSB))
            VTFD(NR,NSB)=SQRT(RTFD(NR,NSB)*1.D3*AEE/AMFD(NSB))
         ENDDO
!         WRITE(*,'(3I4, 4E14.6)') NRANK, NR, NPSTART, RTFD(NR,1), RTFD(NR,2)

         IF(MODEL_EX_READ_Tn.ne.0)THEN
            DO NSA=1,NSAMAX ! temporal
               NS=NS_NSA(NSA)
               RNFP(NR,NSA)=RN_TEMP(NR,NS)
               RTFP(NR,NSA)=RT_TEMP(NR,NS)
               PTFP(NR,NSA)=SQRT(RTFP(NR,NSA)*1.D3*AEE*AMFP(NSA))
               VTFP(NR,NSA)=SQRT(RTFP(NR,NSA)*1.D3*AEE/AMFP(NSA))
            END DO
         END IF

         IF(MODEL_IMPURITY.eq.1)THEN
            zeff_imp= (target_zeff*SPITOT - RNFP0(1)) &
                 /(SPITOT-RNFP0(1))
            PZ(N_impu)=zeff_imp
            NS=NS_NSA(N_impu)
            AEFP(N_impu)=PZ(NS)*AEE 
            AEFD(N_impu)=PZ(NS)*AEE 
         END IF

         RNE=PLF(1)%RN
         RTE=(PLF(1)%RTPR+2.D0*PLF(1)%RTPP)/3.D0
         E_EDGEM=0.D0

!-----   Coulomb log

         ISW_CLOG=0 ! =0 Wesson, =1 NRL
         DO NSA=1,NSAMAX
            NSFP=NS_NSB(NSA)
            RNA=RN_TEMP(NR,NSFP)
            RTA=RT_TEMP(NR,NSFP)
            DO NSB=1,NSBMAX
               NSFD=NS_NSB(NSB)
               RNB=RN_TEMP(NR,NSFD)
               RTB=RT_TEMP(NR,NSFD)
               IF(ISW_CLOG.eq.0)THEN
                  IF(PZ(NSFP).eq.-1.and.PZ(NSFD).eq.-1) THEN !e-e
                     RLNRL=14.9D0-0.5D0*LOG(RNE)+LOG(RTE) 
                  ELSEIF(PZ(NSFP).eq.-1.OR.PZ(NSFD).eq.-1) THEN
                     RLNRL=15.2D0-0.5D0*LOG(RNE)+LOG(RTE) ! e-i T>10eV
                  ELSE
                     RLNRL=17.3D0-0.5D0*LOG(RNE)+1.5D0*LOG(RTB) ! i-i T < m_i/m_p*10 keV, single charge
                  ENDIF
               ELSEIF(ISW_CLOG.eq.1)THEN
                  IF(PZ(NSFP).eq.-1.and.PZ(NSFD).eq.-1) THEN !e-e
                     RLNRL=23.5D0-LOG(SQRT(RNA*1.D14)*(RTA*1.D3)**(-1.25D0) ) - &
                          SQRT(1.D-5+(LOG(RTA*1.D3)-2.D0)**2/16.D0 )
                  ELSEIF(PZ(NSFP).eq.-1.OR.PZ(NSFD).eq.-1) THEN
                     RLNRL=24.D0-LOG(SQRT(RNA*1.D14)*(RTA*1.D3)**(-1) ) ! Ti*(me/mi) < 10eV < Te
                  ELSE
                     RLNRL=23.D0-LOG(PZ(NSFP)*PZ(NSFD)*(PA(NSFP)+PA(NSFD)) &
                          /(PA(NSFP)*(RTB*1.D3)+PA(NSFD)*(RTA*1.D3))* &
                          SQRT((RNA*1.D14)*PZ(NSFP)**2/(RTA*1.D3) + (RNB*1.D14)*PZ(NSFD)**2/(RTB*1.D3) ) )
                  ENDIF
               END IF

!               IF(NRSTART.eq.NRMAX) WRITE(*,*) NR,"Coulomb log",NSA,NSB,RLNRL
               LNLAM(NR,NSB,NSA)=RLNRL
               FACT=AEFP(NSA)**2*AEFD(NSB)**2*RLNRL/(4.D0*PI*EPS0**2)
               RNUD(NR,NSB,NSA)=FACT*RNFP0(NSA)*1.D20 &
                      /(SQRT(2.D0)*VTFD(NR,NSB)*PTFP0(NSA)**2)
               RNUF(NR,NSB,NSA)=FACT*RNFP0(NSA)*1.D20 &
                      /(2*AMFD(NSB)*VTFD(NR,NSB)**2*PTFP0(NSA))

            ENDDO
            IF(NRANK.eq.0)THEN
               tau_ta0(NSA)=(4.D0*PI*EPS0**2)*PTFP0(NSA)**3 &
                    /( AMFP(NSA)*AEFP(NSA)**4*LNLAM(1,NSA,NSA)*RNFP0(NSA)*1.D20 )
            END IF
         ENDDO
      ENDDO

      IF(NRANK.eq.0) WRITE(6,'(A,1P5E14.6)') "tau_ta0 = ",(tau_ta0(NSA),NSA=1,NSAMAX)
      IF(NRANK.eq.0.and.NSBMAX.ge.2)THEN
         SUM=0.D0
         DO NSB=2,NSBMAX
            SUM = SUM + RNFD0(NSB)*AEFD(NSB)**2
         END DO
         ZEFF = SUM/(RNFD0(1)*AEFD(1)**2)
         WRITE(6,'(A,1PE12.4)') " ZEFF = ", ZEFF
         WRITE(6,'(A,1P2E12.4)') "LN_LAMBDA = ", LNLAM(1,1,1), LNLAM(1,2,1)
      END IF

      CALL mtx_broadcast1_real8(ZEFF)
      call mtx_broadcast_real8(tau_ta0,nsamax)
!     ----- set relativistic parameters -----

      IF (MODELR.EQ.0) THEN
         DO NS=1,NSMAX
            THETA0(NS)=0.D0
            DO NR=NRSTART,NREND
               THETA(NR,NS)=0.D0
            ENDDO
         ENDDO
      ELSE
         DO NS=1,NSMAX
            RTFD0L=(PTPR(NS)+2.D0*PTPP(NS))/3.D0
            AMFDL=PA(NS)*AMP
            THETA0(NS)=RTFD0L*1.D3*AEE/(AMFDL*VC*VC)
            DO NR=NRSTART,NREND
               THETA(NR,NS)=THETA0(NS)*RT_TEMP(NR,NS)/RTFD0L
            ENDDO
         ENDDO
      ENDIF

! ---- NBI energy loss time
!      IF(MODEL_NBI.ne.0.and.NRANK.eq.0.and.MODEL_CX_LOSS.eq.1)THEN ! assuming AMFD(NSA_F1)=D, D beam, ZEFF=1

      IF(MODEL_NBI.ne.0.and.MODEL_CX_LOSS.eq.1)THEN ! assuming AMFD(NSA_F1)=D, D beam, ZEFF=1
         NS=NS_NSA(NSA_F1)
         A_D=RN_TEMP(NRSTART,1)*1.D20*AEFD(1)**4*LNLAM(NRSTART,1,NSA_F1)/(2.D0*PI*EPS0**2*AMFP(NSA_F1)**2)
         tau_se_E0 = (3.D0*sqrt(2.D0*PI)*(RT_TEMP(NRSTART,1)*AEE*1.D3)**1.5)/(sqrt(AMFD(1))*AMFP(NSA_F1)*A_D) 

         Ebeam0= 140.D3
         log_energy = dlog10(Ebeam0)
         sigma_cx0 = 0.6937D-14*(1.D0-0.155D0*log_energy)**2/ &
              (1.D0+0.1112D-14*Ebeam0**3.3D0)*1.D-4 * SQRT(Ebeam0*AEE/AMFD(NSA_F1))

         j=0
         DO i=1, 100
            Ebeam1= (140.D0-1.D0*i)*1.D3
            log_energy = dlog10(Ebeam1)
            sigma_cx = 0.6937D-14*(1.D0-0.155D0*log_energy)**2/ &
                 (1.D0+0.1112D-14*Ebeam1**3.3D0)*1.D-4 *SQRT(Ebeam1*AEE/AMFD(NSA_F1))
            IF(sigma_cx.le.sigma_cx0*exp(1.D0))THEN
               j=i
            END IF
         END DO
         Ebeam1= (140.D0-1.D0*j)*1.D3
         log_energy = dlog10(Ebeam1)
         sigma_cx = 0.6937D-14*(1.D0-0.155D0*log_energy)**2/ &
              (1.D0+0.1112D-14*Ebeam1**3.3D0)*1.D-4 *SQRT(Ebeam1*AEE/AMFD(NSA_F1))

         E_CR=14.8D0*PA(NS)/PA(NS)**(2.D0/3.D0)*RT_TEMP(NRSTART,1)*1.D3
!         tau_se_E0E1=tau_se_E0*log(Ebeam0/Ebeam1)*0.5D0
         tau_se_E0E1=tau_se_E0*&
              log( (Ebeam0**(1.5D0)+E_CR**(1.5D0) )/( Ebeam1**(1.5D0)+E_CR**(1.5D0) ) )&
              /3.D0

         log10_neu0=LOG10(RN_NEU0)
         log10_neus=LOG10(RN_NEUS)
!     neutral gas profile
         alpha=2.5D0
         beta=0.8D0
         N_NEUT=10**( (log10_neu0-log10_neus)*(1.D0-RM(NRSTART)**alpha)**beta+log10_neus )

         tau_cx_E1 = 1.D0/(N_NEUT*1.D20*EXP(1.D0)*sigma_cx0)

         IF(NRANK.eq.0)THEN
            WRITE(*,'(A,E14.6)') "E_CR on axis=   ", E_CR
            WRITE(*,'(A,1PE14.6,2E14.6)') "E0, tau_se_E0 on axis=   ", Ebeam0, tau_se_E0
            WRITE(*,'(A,1PE14.6,2E14.6)') "E1, tau_se_E0E1, wo E_CR on axis= ",&
                 Ebeam1, tau_se_E0E1, tau_se_E0*log(Ebeam0/Ebeam1)*0.5D0
            WRITE(*,'(A,E14.6)') "tau_cx_E1 on axis=   ", tau_cx_E1
         END IF
         IF(NRSTART.eq.NRMAX.and.NPSTART.eq.1.and.NSASTART.eq.1)THEN
            WRITE(*,'(A,1PE14.6,2E14.6)') "E0, tau_se_E0, on edge=   ", Ebeam0, tau_se_E0
            WRITE(*,'(A,1PE14.6,2E14.6)') "E1, tau_se_E0E1, on edge= ", Ebeam1, tau_se_E0E1, tau_se_E0*log(Ebeam0/Ebeam1)*0.5D0
            WRITE(*,'(A,E14.6)') "tau_cx_E1 on edge=   ", tau_cx_E1
         END IF
      END IF
! ----set runaway
      IF(MODEL_DISRUPT.ne.0.and.nt_init.eq.0)THEN 
         IF(MODEL_IMPURITY.eq.1)THEN
            DO NSB=1,NSBMAX
               DO NR=NRSTART,NREND
                  RN_MGI(NR,NSB)=RNFD(NR,NSB)
               END DO
               RN0_MGI(NSB)=RNFD0(NSB)
            END DO
            IF(NRANK.eq.0) WRITE(*,'(A,1P2E14.6)') "zeff_imp, pz = ",zeff_imp, PZ(n_impu)
         END IF

         CALL set_initial_disrupt_param
         CALL set_post_disrupt_Clog_f
         CALL set_post_disrupt_Clog
         call mtx_broadcast_real8(POST_tau_ta0_f,nsamax)
      END IF
      RWS_PARA(:,:)=0.D0
      RWS_PERP(:,:)=0.D0

      END SUBROUTINE fp_set_normalize_param
!==============================================================
      SUBROUTINE Coulomb_log

      IMPLICIT NONE
      INTEGER:: NR, NSA, NSFP, NSFD, NSB, ISW_CLOG
      DOUBLE PRECISION:: RTA,RTB,RNA,RNB, RLNRL, FACT,RNE
      double precision,dimension(NSAMAX,NSBMAX):: CLOG
      double precision:: VTFDL, PTFDL

      DO NR=NRSTART,NRENDWM
         ISW_CLOG=0 ! =0 Wesson, =1 NRL
         DO NSA=1,NSAMAX
            NSFP=NS_NSA(NSA)
            RNE=RN_TEMP(NR,1)
            RNA=RN_TEMP(NR,NSFP)
            RTA=RT_TEMP(NR,NSFP)
            DO NSB=1,NSBMAX
               NSFD=NS_NSB(NSB)
               IF(MODEL_disrupt.eq.0)THEN
!                  IF(NSFD.NE.0)THEN
                     RNB=RN_TEMP(NR,NSFD)
                     RTB=RT_TEMP(NR,NSFD)
!                  ELSE
!                     RNB=RNFD(NR,NSB)
!                     RTB=RTFD(NR,NSB)
!                  END IF
               ELSE
                  RNB=RNFD(NR,NSB)
                  RTB=RT_quench(NR)
               END IF

               vtfdl=SQRT(RTB*1.D3*AEE/AMFD(NSB))
               ptfdl=SQRT(RTB*1.D3*AEE*AMFD(NSB))

               IF(ISW_CLOG.eq.0)THEN
                  IF(PZ(NSFP).eq.-1.and.PZ(NSFD).eq.-1) THEN !e-e
                     RLNRL=14.9D0-0.5D0*LOG(RNA)+LOG(RTA) 
                  ELSEIF(PZ(NSFP).eq.-1.and.PZ(NSFD).ne.-1) THEN
                     RLNRL=15.2D0-0.5D0*LOG(RNA)+LOG(RTA) ! e-i T>10eV
                  ELSEIF(PZ(NSFP).ne.-1.and.PZ(NSFD).eq.-1) THEN
                     RLNRL=15.2D0-0.5D0*LOG(RNB)+LOG(RTB) ! e-i T>10eV
                  ELSE
                     RLNRL=17.3D0-0.5D0*LOG(RNE)+1.5D0*LOG(RTB) ! i-i T < m_i/m_p*10 keV, single charge
                  ENDIF
               ELSEIF(ISW_CLOG.eq.1)THEN
                  IF(PZ(NSFP).eq.-1.and.PZ(NSFD).eq.-1) THEN !e-e
                     RLNRL=23.5D0-LOG(SQRT(RNA*1.D14)*(RTA*1.D3)**(-1.25D0) ) - &
                          SQRT(1.D-5+(LOG(RTA*1.D3)-2.D0)**2/16.D0 )
                  ELSEIF(PZ(NSFP).eq.-1.and.PZ(NSFD).ne.-1) THEN
                     RLNRL=24.D0-LOG(SQRT(RNA*1.D14)*(RTA*1.D3)**(-1) ) ! Ti*(me/mi) < 10eV < Te
                  ELSEIF(PZ(NSFP).ne.-1.and.PZ(NSFD).eq.-1) THEN
                     RLNRL=24.D0-LOG(SQRT(RNB*1.D14)*(RTB*1.D3)**(-1) ) ! Ti*(me/mi) < 10eV < Te
                  ELSE
                     RLNRL=23.D0-LOG(PZ(NSFP)*PZ(NSFD)*(PA(NSFP)+PA(NSFD)) &
                          /(PA(NSFP)*(RTB*1.D3)+PA(NSFD)*(RTA*1.D3))* &
                          SQRT((RNA*1.D14)*PZ(NSFP)**2/(RTA*1.D3) + (RNB*1.D14)*PZ(NSFD)**2/(RTB*1.D3) ) )
                  ENDIF
               END IF
               LNLAM(NR,NSB,NSA)=RLNRL
               CLOG(NSA,NSB)=RLNRL
               FACT=AEFP(NSA)**2*AEFD(NSB)**2*RLNRL/(4.D0*PI*EPS0**2)
               RNUD(NR,NSB,NSA)=FACT*RNFP0(NSA)*1.D20 &
                      /(SQRT(2.D0)*VTFD(NR,NSB)*PTFP0(NSA)**2)
               RNUF(NR,NSB,NSA)=FACT*RNFP0(NSA)*1.D20 &
                      /(2*AMFD(NSB)*VTFD(NR,NSB)**2*PTFP0(NSA))
            ENDDO
         ENDDO
      ENDDO

      END SUBROUTINE Coulomb_log
!==============================================================
      SUBROUTINE fp_continue(ierr)

      USE fpnfrr
      USE libmtx
      USE fpnflg
      IMPLICIT NONE
      integer :: ierr,NSA,NR,NP,NTH,NSB

      IF(NRANK.eq.0) &
      WRITE(6,*) "----- SET COEFFICIENTS AND DISTRIBUTION FUNCTIONS -----"

      N_IMPL=0
      CALL fusion_source_init

      CALL Define_Bulk_NP
      CALL FP_COEF(0)

      DO NSA=NSASTART,NSAEND
         DO NR=NRSTART,NREND
            DO NP=NPSTARTW,NPENDWM
               DO NTH=1,NTHMAX
                  F(NTH,NP,NR)=FNSP(NTH,NP,NR,NSA)
               END DO
            END DO
         END DO
         CALL FPWEIGHT(NSA,IERR)
      END DO
      IF(MODELS.ne.0)THEN
         CALL mtx_set_communicator(comm_nsa)
         CALL source_allreduce(SPPF)
         CALL mtx_reset_communicator
      END IF

      ISAVE=0

      IERR=0
 
      RETURN
      END SUBROUTINE fp_continue
!-------------------------------------------------------------
      SUBROUTINE fp_set_initial_value_from_f

      USE libmtx
      IMPLICIT NONE

      RNS(:,:)=0.D0
      IF(nrank.EQ.0) WRITE(*,*) "SET INITIAL VALUE FROM f"
      IF(NRANK.eq.0.and.MODEL_disrupt.ne.0)THEN
         CALL display_disrupt_initials
      END IF

      CALL FPSSUB

      IF(nrank.EQ.0) THEN
         CALL FPSGLB
         CALL FPSPRF
         CALL FPWRTGLB
         CALL FPWRTPRF
      ENDIF
!      CALL mtx_broadcast_real8(RT_T,NRMAX*NSBMAX)
      CALL mtx_broadcast_real8(RNS,NRMAX*NSAMAX)
      CALL mtx_broadcast1_integer(NTG1)
      CALL mtx_broadcast1_integer(NTG2)

      END SUBROUTINE fp_set_initial_value_from_f
!-------------------------------------------------------------
      SUBROUTINE fp_prep(ierr)

      USE plprof
      USE fpnfrr
      USE libmtx
      USE fpnflg
      USE FP_READ_FIT
      USE FPOUTDATA

      Implicit none

      integer :: ierr,NSA,NS,NR,N,NSW,i,NSFP,NSB
      real:: gut1, gut2, gut_prep
      real(8):: SIGMA
      real(8),dimension(:),allocatable:: conduct_temp, E1_temp
      integer,dimension(6):: idata
      integer,dimension(6*nsize):: idata2

      CALL GUTIME(gut1)

!     ----- Initialize time counter -----

      TIMEFP=0.D0
      NTG1=0
      NTG2=0
      NT_init=0
      gut_comm(:)=0.0

      CALL fp_comm_setup

!     ----- Allocate variables -----

      CALL fp_allocate
      call fp_allocate_ntg1
      call fp_allocate_ntg2

      CALL mtx_set_communicator(comm_nr)
      allocate(MTXLEN(nsize),MTXPOS(nsize))

      CALL mtx_set_communicator(comm_nsanr)
      allocate(SAVLEN(nsize)) 
      allocate(SAVPOS(nsize,NSAEND-NSASTART+1)) 
      CALL mtx_reset_communicator
      allocate(Rank_Partition_Data(6,0:nsize-1))

      idata(1)=NPSTARTW
      idata(2)=NPENDWM
      idata(3)=NRSTARTW
      idata(4)=NRENDWM
      idata(5)=NSASTART
      idata(6)=NSAEND

      CALL mtx_gather_integer(idata,6,idata2)

      IF(NRANK.eq.0)THEN
         DO N=0, nsize-1
            DO I=1,6
               Rank_Partition_Data(I,N)=idata2(I+N*6)
            END DO
         END DO
      END IF

!     ----- Get mtxlen and mtxpos -----
!     MTXLEN(NRANK+1): the number of NR grid points for each RANK
!     MTXPOS(NRANK): 
!
      CALL mtx_set_communicator(comm_nr)
      CALL mtx_allgather1_integer(nrend-nrstart+1,mtxlen)
      CALL mtx_allgather1_integer(nrstart-1,mtxpos)

      CALL mtx_set_communicator(comm_nsanr)
      CALL mtx_allgather1_integer(nrend-nrstart+1,savlen)
      NSW=NSAEND-NSASTART+1
      DO N=1,NSW
         NSA=N+NSASTART-1
         CALL mtx_allgather1_integer((nsa-1)*NRMAX+NRSTART,savpos(1:nsize,N))
      END DO
      CALL mtx_reset_communicator

      CALL fp_set_nsa_nsb

!     ----- create meches -----
      CALL fp_mesh(ierr)

!     ----- Initialize diffusion coef. -----
      call FPCINI
      RNS_DELF(:,:)=0.D0

!     ----- set parameters for target species -----
      CALL fp_set_normalize_param

!     ----- Initialize velocity distribution function of all species -----

      CALL FNSP_INIT     
      CALL FNSP_INIT_EDGE
      IF(NRANK.EQ.0) WRITE(6,*) 'END INIT'

!     ----- set background f

      CALL mtx_set_communicator(comm_nsa)
      CALL update_fnsb
      CALL mtx_reset_communicator

!     ----- READ FIT3D result for NBI -----
      IF(MODEL_NBI.eq.2)THEN
         DO NS=1, NBEAMMAX
            NSFP=NSSPB(NS)
            IF(PA(NSFP).eq.1)THEN
               CALL READ_FIT3D_H
            ELSEIF(PA(NSFP).eq.2)THEN
               CALL READ_FIT3D_D
            END IF
            CALL SV_WEIGHT_R
        END DO
      END IF

!     ----- set parallel electric field -----

      IF(MODEL_DISRUPT.eq.0)THEN
         DO NR=1,NRMAX
            E1(NR)=E0!*E_drei0(1)
         END DO
         DO NR=NRSTART,NREND
            EP(NR)=E1(NR) ! plus
            EM(NR)=0.D0 ! minus
         END DO
         allocate(conduct_temp(NRSTART:NREND))
         DO NR=NRSTART,NREND
            CALL SPITZER_SIGMA(NR,SIGMA)
            conduct_temp(NR)=sigma
         END DO
         CALL mtx_set_communicator(comm_nr)
         call mtx_allgather_real8(conduct_temp,NREND-NRSTART+1,conduct_sp)
         CALL mtx_reset_communicator
      ELSE
         allocate(conduct_temp(NRSTART:NREND))
         allocate(E1_temp(NRSTART:NREND))
         DO NR=NRSTART,NREND
            CALL SPITZER_SIGMA(NR,SIGMA)
            conduct_temp(NR)=sigma
            IF(MODELE.eq.0)THEN
               E1_temp(NR)=RJ_ohm(NR)/SIGMA*1.D6 ! fit to initial current prof
            ELSEIF(MODELE.eq.1)THEN
               E1_temp(NR)=E0*E_drei0(1) ! uniform
            ELSEIF(MODELE.eq.2)THEN
               E1_temp(NR)=E0*(1.D0-RM(NR)**1.5)**1 ! arbitrary profile
            END IF
            EP(NR)=0.D0
         END DO
         CALL mtx_set_communicator(comm_nr)
         call mtx_allgather_real8(conduct_temp,NREND-NRSTART+1,conduct_sp)
         call mtx_allgather_real8(E1_temp,NREND-NRSTART+1,E1)
         CALL mtx_reset_communicator
      END IF

!continue start
      IF(MODELS.eq.3) CALL NF_LG_FUNCTION
      IF(MODELS.ne.0) CALL NF_REACTION_COEF

      CALL fp_continue(ierr)
      IF(MODELS.ge.2)THEN
         CALL ALLREDUCE_NF_RATE
         CALL PROF_OF_NF_REACTION_RATE(1)
      END IF
      CALL fp_set_initial_value_from_f

      CALL GUTIME(gut2)
      gut_prep=gut2-gut1
      IF(NRANK.eq.0) WRITE(6,'(A,E14.6)') "---------------PREP_TIME=", gut_prep

      IF(OUTPUT_TXT_DELTA_F.eq.1.and.NRANK.eq.0) CALL OUT_TXT_FNS_DEL
 
      RETURN
      END subroutine fp_prep
!-----
      END MODULE fpprep
