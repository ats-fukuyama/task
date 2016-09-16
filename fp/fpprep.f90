!     $Id: fpprep.f90,v 1.41 2013/02/08 07:36:24 nuga Exp $

! *****************************
!     PREPARATION OF FPLOOP
! *****************************
      MODULE fpprep

      USE fpcomm
      USE fpinit
      USE fpsave
      USE fpcoef
      USE fpbounce
      USE equnit_mod
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
      integer :: ierr,NSA,NSB,NS,NR,NP,NTH,id,NP2
      character(LEN=80)::line 
      real(kind8)::rhon,rhol,rhol1,rhol2,A1,epsl,fact,rint0,es,rint2,ql,BT
      real(kind8),DIMENSION(:),POINTER:: work,workg
      real(kind8)::suml, delh, etal, x, psib, pcos
      integer:: NG

!     ----- exec EQ -----

      IF(MODELG.EQ.3.OR.MODELG.EQ.8) THEN
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
!      WRITE(6,*) 'RKAP=',RKAP,' set to 1.0'
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
         IF(NRANK.eq.0) &
         write(6,'(A,I5,1P5E12.4)') 'nr,rm,rsrhon,epsrm,bpm,ql=', &
              NR,RM(NR),RSRHON(RHON),EPSRM(NR),BPM(NR),QL
      ENDDO
!      RHON=RG(NRMAX+1)
      RHON=RM(NRMAX)+DELR
      CALL pl_qprf(RHON,QL)
      QLM(NRMAX+1)=QL
      BT=BB
      BP= RSRHON(RHON)*BT/(RR*QL)
      EPSRM(NRMAX+1)=RSRHON(RHON)/RR
      BPM(NRMAX+1)= RSRHON(RHON)*BT/(RR*QL)
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
      DO NSB=1,NSBMAX
         DELP(NSB)=PMAX(NSB)/NPMAX
         DO NP=1,NPMAX
            PG(NP,NSB)=DELP(NSB)*(NP-1)
            PM(NP,NSB)=DELP(NSB)*(NP-0.5D0)
         ENDDO
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

!     set momentum volume element: VOLP
      DO NSB=1,NSBMAX
      DO NP=NPSTART,NPEND
      DO NTH=1,NTHMAX
         VOLP(NTH,NP,NSB)=2.D0*PI*SINM(NTH)*PM(NP,NSB)**2*DELP(NSB)*DELTH
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
      CASE(3:4) ! toroidal
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
               ETAM_G(NTH,NR)=PI*0.5D0
               RLAMDA_G(NTH,NR)=1.D0
               RLAMDC_G(NTH,NR)=1.D0
            ENDDO
            DO NTH=1,NTHMAX+1
               ETAG(NTH,NR)=PI/2.D0
               ETAG_G(NTH,NR)=PI/2.D0
            ENDDO
         ENDDO

         DO NTH=1, NTHMAX
            RLAMDA_GG(NTH,NRMAX+1)=1.D0
         END DO
         DO NR=1,NRMAX+1
            RFSADG(NR)=1.D0
            RFSAD_GG(NR)=1.D0
         END DO
      ELSE
         DO NR=1,NRMAX+1
            CALL SET_RFSAD(NR)
         END DO
         CALL SET_BOUNCE_PARAM

      END IF ! MODELA

      allocate(work(nrstart:nrend),workg(NRMAX))

      CALL mtx_set_communicator(comm_nr)
      DO NTH=1,NTHMAX
         DO NR=NRSTART,NREND
            work(NR)=RLAMDA(NTH,NR)
         ENDDO
         CALL mtx_allgatherv_real8(work(NRSTART:NREND),MTXLEN(NRANK+1), &
                                   workg,NRMAX,MTXLEN,MTXPOS)
         DO NR=1,NRMAX
            RLAMDAG(NTH,NR)=workg(NR)
         ENDDO
         RLAMDAG(NTH,NRMAX+1)=RLAMDAG(NTH,NRMAX)
      ENDDO
      
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
            work(NR)=RLAMDA_G(NTH,NR)
         ENDDO
         CALL mtx_allgatherv_real8(work(NRSTART:NREND),MTXLEN(NRANK+1), &
                                   workg,NRMAX,MTXLEN,MTXPOS)
         DO NR=1,NRMAX
            RLAMDA_GG(NTH,NR)=workg(NR)
         ENDDO
      ENDDO

      DO NTH=1,NTHMAX
         DO NR=NRSTART,NREND
            work(NR)=ETAM_G(NTH,NR)
         ENDDO
         CALL mtx_allgatherv_real8(work(NRSTART:NREND),MTXLEN(NRANK+1), &
                                   workg,NRMAX,MTXLEN,MTXPOS)
         DO NR=1,NRMAX
            ETAM_GG(NTH,NR)=workg(NR)
         ENDDO
      ENDDO
      CALL mtx_reset_communicator

      IF(MODELA.eq.1)THEN
      IF(NRANK.eq.0)THEN
      open(8,file='RLAMDAG.dat')
      DO NR =1, NRMAX
         DO NTH=1,NTHMAX
            WRITE(8,'(2I4, 4E14.6)') NR, NTH, COSM(NTH), RLAMDAG(NTH,NR), RFSADG(NR)
         END DO
         WRITE(8,*) " "
         WRITE(8,*) " "
      END DO
      close(8)
      END IF
      END IF

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
      integer:: NSA,NR, NTH, NP

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
      INTEGER:: ierr, colors, keys, N, NREND1, NPEND1, I
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


      keys=comm_nsa%rankl
      NSASTART = (NSAMAX/N_partition_s)*keys+1 !3D
      NSAEND =   (NSAMAX/N_partition_s)*(keys+1)

      imtxsize=nthmax*npmax*nrmax
      IF(modeld.EQ.0) THEN
         imtxwidth=4*nthmax-1
      ELSE
         imtxwidth=4*nthmax*npmax-1
      ENDIF
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

!      IF(nrank.EQ.0) THEN
!         write(6,'(A,2I10)') '  imtxsize,imtxwidth=',imtxsize,imtxwidth
!         write(6,'(A,A,A)') '     nrank   imtxstart   imtxend   npstart    npend', &
!                          '    nrstart     nrend   nmstart     nmend', &
!                          '      nsastart      nsaend'
!         DO N=1,nsize
!            write(6,'(11I10)') N,ima1(N),ima2(N),npa1(N),npa2(N),nra1(N),nra2(N), &
!                              nma1(N),nma2(N),insa1(N),insa2(N)
!         ENDDO
!      ENDIF

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

      IMPLICIT NONE
      INTEGER:: NTH,NP,NR,NSA,NSB,NS,NSBA
      REAL(8):: FL

      DO NSA=NSASTART,NSAEND
         NS=NS_NSB(NSA)
         DO NR=NRSTARTW,NRENDWM
            IF(NR.ge.1.and.NR.le.NRMAX)THEN
               DO NP=NPSTARTW,NPENDWM
                  FL=FPMXWL(PM(NP,NSA),NR,NS)
                  DO NTH=1,NTHMAX
                     FNSP(NTH,NP,NR,NSA)=FL
                     FNS0(NTH,NP,NR,NSA)=FL
                  END DO
               ENDDO
            END IF
         END DO
      END DO

      END SUBROUTINE FNSP_INIT
!-------------------------------------------------------------
      SUBROUTINE FNSP_INIT_EDGE

      IMPLICIT NONE
      INTEGER:: NTH,NP,NR,NSA,NSB,NS,NSBA
      integer:: isw
      REAL(8):: FL

      DO NSA=NSASTART,NSAEND
         NS=NS_NSB(NSA)
         NR=NRMAX+1
         DO NP=NPSTARTW,NPENDWM
            FL=FPMXWL(PM(NP,NSA),NR,NS)
            DO NTH=1,NTHMAX
               FS1(NTH,NP,NSA)=FL ! rho=1.0 fixed value
            END DO
         ENDDO
      END DO
!     ----- set boundary distribution functions -----
      
      DO NSA=NSASTART,NSAEND
         NS=NS_NSA(NSA)
         NSBA=NSB_NSA(NSA)
         DO NP=NPSTARTW,NPENDWM
            FL=FPMXWL(PM(NP,NSBA),0,NS)
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
                  FS2(NTH,NP,NS)=FL ! at R=1.0+DELR/2
               ENDDO
            ENDDO
         ENDDO
      ELSEIF(MODELD_boundary.eq.1)THEN
         IF(NREND.eq.NRMAX)THEN
            DO NSA=NSASTART,NSAEND
               NS=NS_NSA(NSA)
               DO NP=NPSTARTW,NPENDWM
                  DO NTH=1,NTHMAX
                     FS2(NTH,NP,NS) = 2.D0*FS1(NTH,NP,NS)-FNSP(NTH,NP,NRMAX,NSA) ! linear
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            FS2(:,:,:)=0.D0
         END IF
      END IF

      END SUBROUTINE FNSP_INIT_EDGE
!-------------------------------------------------------------
      SUBROUTINE fp_set_bounceaverage_param

      USE libmtx
      IMPLICIT NONE
      INTEGER:: NTH, NP, NR, NSA, NSB, NSBA
      real(kind8):: RSUM1,RSUM2,RSUM3,RSUM4,RSUM5,RSUM6
      real(kind8),DIMENSION(:),POINTER:: work,workg

      CALL mtx_set_communicator(comm_np) 
      IF(MODELA.eq.1)THEN
         NSB=NSASTART ! arbitrary
         DO NR=NRSTART,NREND
            RSUM1=0.D0
            RSUM2=0.D0
            RSUM3=0.D0
            RSUM4=0.D0
            DO NTH=1,NTHMAX
               RSUM1 = RSUM1+RLAMDAG(NTH,NR)/RFSADG(NR)*SINM(NTH)*DELTH
               RSUM2 = RSUM2+SINM(NTH)*DELTH
               RSUM3 = rsum3+RLAMDA_GG(NTH,NR)/RFSAD_GG(NR)*SINM(NTH)*DELTH
               RSUM4 = rsum4+SINM(NTH)*DELTH
            END DO
            IF(RSUM1.EQ.0.D0) &
                 WRITE(6,'(1P3E12.4)') VOLP(1,1,NSB),FNSP(1,1,1,NSB),RLAMDA(1,1)

            RCOEFN(NR)=1.D0
            RCOEFN_G(NR)=1.D0
         END DO
      ELSE
         DO NR=NRSTART,NREND
            RCOEFN(NR)=1.D0
            RCOEFN_G(NR)=1.D0
         ENDDO
      END IF

!
      IF(MODELA.eq.1)THEN
         NSBA=NSASTART
         RSUM1=0.D0
         RSUM2=0.D0
         RSUM3=0.D0
         RSUM4=0.D0
         RSUM5=0.D0
         RSUM6=0.D0
         DO NTH=1,NTHMAX
            RSUM1 = RSUM1+RLAMDA_GG(NTH,1)/RFSAD_GG(1)*SINM(NTH)*DELTH
            RSUM2 = RSUM2+SINM(NTH)*DELTH
            RSUM3 = RSUM3+RLAMDA_GG(NTH,NRMAX+1)/RFSAD_GG(NRMAX+1)*SINM(NTH)*DELTH
!            RSUM4 = RSUM4+SINM(NTH)*DELTH
            RSUM5 = RSUM5+RLAMDAG(NTH,NRMAX+1)/RFSADG(NRMAX+1)*SINM(NTH)*DELTH
!            RSUM6 = RSUM6+SINM(NTH)*DELTH
         END DO

!         RCOEFN_GG(1) = RSUM2/RSUM1
!         RCOEFN_GG(NRMAX+1) = RSUM2/RSUM3
!         RCOEFNG(NRMAX+1) = RSUM2/RSUM5

         RCOEFN_GG(1) = 1.D0
         RCOEFN_GG(NRMAX+1) = 1.D0
         RCOEFNG(NRMAX+1) = 1.D0
      ELSE
         RCOEFN_GG(1) = 1.D0
         RCOEFN_GG(NRMAX+1) = 1.D0
         RCOEFNG(NRMAX+1) = 1.D0
      END IF
!      CALL mtx_reset_communicator
!!!!
      IF(MODELA.eq.1)THEN
         DO NR=NRSTART,NREND
            RSUM1=0.D0
            RSUM2=0.D0
            DO NTH=1,NTHMAX/2
               RSUM1 = RSUM1 + RLAMDA(NTH,NR)*COSM(NTH)*SINM(NTH)*DELTH
               RSUM2 = RSUM2 + COSM(NTH)*SINM(NTH)*DELTH
            END DO
            RCOEFJ(NR)=RSUM2/RSUM1*RFSADG(NR)
         END DO
      ELSE
         DO NR=NRSTART,NREND
            RCOEFJ(NR)=1.D0
         END DO
      END IF
!!!!!!!!!!!!!!!!!!!!!

      allocate(work(nrstart:nrend),workg(NRMAX))

      CALL mtx_set_communicator(comm_nr)

      DO NR=NRSTART,NREND
         work(NR)=RCOEFN(NR)
      ENDDO
      CALL mtx_allgatherv_real8(work(NRSTART:NREND),MTXLEN(NRANK+1), &
                                workg,NRMAX,MTXLEN,MTXPOS)
      DO NR=1,NRMAX
         RCOEFNG(NR)=workg(NR)
      ENDDO

      DO NR=NRSTART,NREND
         work(NR)=RCOEFN_G(NR)
      ENDDO
      CALL mtx_allgatherv_real8(work(NRSTART:NREND),MTXLEN(NRANK+1), &
                                workg,NRMAX,MTXLEN,MTXPOS)
      DO NR=2,NRMAX
         RCOEFN_GG(NR)=workg(NR)
      ENDDO

      DO NR=NRSTART,NREND
         work(NR)=RCOEFJ(NR)
      ENDDO
      CALL mtx_allgatherv_real8(work(NRSTART:NREND),MTXLEN(NRANK+1), &
                                workg,NRMAX,MTXLEN,MTXPOS)
      DO NR=1,NRMAX
         RCOEFJG(NR)=workg(NR)
      ENDDO
      CALL mtx_reset_communicator

      deallocate(work,workg)

      DO NR=1,NRMAX
         DO NTH=1,NTHMAX
            RLAMDAG(NTH,NR)=RLAMDAG(NTH,NR)*RCOEFNG(NR)
         END DO
      END DO
      DO NR=NRSTART,NREND
         DO NTH=1,NTHMAX
            RLAMDA(NTH,NR)=RLAMDA(NTH,NR)*RCOEFNG(NR)
         END DO
      END DO

!      IF(NRANK.eq.0)THEN
!         DO NR=1,NRMAX
!            WRITE(*,'(I3,3E16.8)') NR, RM(NR), RCOEFNG(NR), RCOEFJG(NR)
!         END DO
!      END IF


      END SUBROUTINE fp_set_bounceaverage_param
!-------------------------------------------------------------
      SUBROUTINE fp_set_normalize_param

      USE plprof
      IMPLICIT NONE
      INTEGER:: NSA, NSB, NS, NSBA, NSFP, NSFD, NR, ISW_CLOG, NP
      TYPE(pl_plf_type),DIMENSION(NSMAX):: PLF
      real(kind8):: RTFD0L, RHON, RNE, RTE, RLNRL, FACT, RNA, RTA, RNB, RTB, SUM

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

!     ----- set profile data -----
      DO NR=NRSTART, NRENDWG
         RHON=RG(NR)
         CALL PL_PROF(RHON,PLF)
         DO NSA=1, NSAMAX
            RNFP_G(NR,NSA)=PLF(NS)%RN
            RTFP_G(NR,NSA)=(PLF(NS)%RTPR+2.D0*PLF(NS)%RTPP)/3.D0
         END DO
      END DO
      DO NR=NRSTART,NRENDWM

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
!            AEFD(NSB)=PZ(NS)*AEE
!            AMFD(NSB)=PA(NS)*AMP
            RNFD(NR,NSB)=PLF(NS)%RN
            RTFD(NR,NSB)=(PLF(NS)%RTPR+2.D0*PLF(NS)%RTPP)/3.D0
            PTFD(NR,NSB)=SQRT(RTFD(NR,NSB)*1.D3*AEE*AMFD(NSB))
            VTFD(NR,NSB)=SQRT(RTFD(NR,NSB)*1.D3*AEE/AMFD(NSB))
         ENDDO

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
            RNA=RNFP(NR,NSA)
            RTA=RTFP(NR,NSA)
            DO NSB=1,NSBMAX
               NSFD=NS_NSB(NSB)
               RNB=RNFD(NR,NSB)
               RTB=RTFD(NR,NSB)
               IF(ISW_CLOG.eq.0)THEN
                  IF(PZ(NSFP).eq.-1.and.PZ(NSFD).eq.-1) THEN !e-e
                     RLNRL=14.9D0-0.5D0*LOG(RNE)+LOG(RTE) 
                  ELSEIF(PZ(NSFP).eq.-1.OR.PZ(NSFD).eq.-1) THEN
                     RLNRL=15.2D0-0.5D0*LOG(RNE)+LOG(RTE) ! e-i T>10eV
                  ELSE
                     RLNRL=17.3D0-0.5D0*LOG(RNE)+1.5D0*LOG(RTFD(NR,NSB)) ! i-i T < m_i/m_p*10 keV, single charge
                  ENDIF
               ELSEIF(ISW_CLOG.eq.1)THEN
                  IF(PZ(NSFP).eq.-1.and.PZ(NSFD).eq.-1) THEN !e-e
                     RLNRL=23.5D0-LOG(SQRT(RNA*1.D14)*(RTA*1.D3)**(-1.25D0) ) - &
                          SQRT(1.D-5+(LOG(RTA*1.D3)-2.D0)**2/16.D0 )
                  ELSEIF(PZ(NSFP).eq.-1.OR.PZ(NSFD).eq.-1) THEN
                     RLNRL=24.D0-LOG(SQRT(RNA*1.D14)*(RTA*1.D3)**(-1) ) ! Ti*(me/mi) < 10eV < Te
                  ELSE
                     RLNRL=23.D0-LOG(PZ(NSA)*PZ(NSB)*(PA(NSA)+PA(NSB)) &
                          /(PA(NSA)*(RTB*1.D3)+PA(NSB)*(RTA*1.D3))* &
                          SQRT((RNA*1.D14)*PZ(NSA)**2/(RTA*1.D3) + (RNB*1.D14)*PZ(NSB)**2/(RTB*1.D3) ) )
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
         DO NSB=1,NSBMAX
            THETA0(NSB)=0.D0
            DO NR=NRSTART,NREND
               THETA(NR,NSB)=0.D0
            ENDDO
         ENDDO
      ELSE
         DO NSB=1,NSBMAX
            THETA0(NSB)=RTFD0(NSB)*1.D3*AEE/(AMFD(NSB)*VC*VC)
            DO NR=NRSTART,NREND
               THETA(NR,NSB)=THETA0(NSB)*RTFD(NR,NSB)/RTFD0(NSB)
            ENDDO
         ENDDO
      ENDIF

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

      END SUBROUTINE fp_set_normalize_param
!==============================================================
      SUBROUTINE Coulomb_log

      IMPLICIT NONE
      INTEGER:: NTH, NP, NR, NSA, NSFP, NSFD, NSB, ISW_CLOG
      DOUBLE PRECISION:: RTE,NTE,RTA,RTB,RNA,RNB, RLNRL, FACT,RNE
      double precision,dimension(NSAMAX,NSBMAX):: CLOG
      double precision:: VTFDL, PTFDL

!      CALL FPSSUB2
!      CALL FPSPRF2

      DO NR=NRSTART,NRENDWM
         ISW_CLOG=0 ! =0 Wesson, =1 NRL
         DO NSA=1,NSAMAX
            NSFP=NS_NSB(NSA)
            IF(MODEL_disrupt.eq.0)THEN
               RNA=RN_IMPL(NR,NSA)
               RTA=RT_T(NR,NSA)
               RNE=RN_IMPL(NR,1)
            ELSE
               RNE=RNFP(NR,1)
               RNA=RNS(NR,NSA)
               RTA=RT_quench(NR)
            END IF

!            RNA=RN_IMPL(NR,NSA)
!            RTA=RT_IMPL(NR,NSA)
            DO NSB=1,NSBMAX
               NSFD=NS_NSB(NSB)
               IF(MODEL_disrupt.eq.0)THEN
                  IF(NSB.le.NSAMAX)THEN
                     RNB=RN_IMPL(NR,NSB)
                     RTB=RT_IMPL(NR,NSB)
                  ELSE
                     RNB=RNFD(NR,NSB)
                     RTB=RTFD(NR,NSB)
                  END IF
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
                     RLNRL=23.D0-LOG(PZ(NSA)*PZ(NSB)*(PA(NSA)+PA(NSB)) &
                          /(PA(NSA)*(RTB*1.D3)+PA(NSB)*(RTA*1.D3))* &
                          SQRT((RNA*1.D14)*PZ(NSA)**2/(RTA*1.D3) + (RNB*1.D14)*PZ(NSB)**2/(RTB*1.D3) ) )
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
!      WRITE(*,'(I3,A,10E16.8)') NR, " Coulomb log=", LNLAM(NR,1,1), LNLAM(NR,1,2), LNLAM(NR,2,1), LNLAM(NR,2,2), RN_IMPL(NR,2), RT_T(NR,2), RT_IMPL(NR,2) 
      ENDDO
!      IF(NR.eq.NRMAX) WRITE(*,'(A,1P4E16.8)') "Coulomb log",CLOG(1,1), CLOG(1,2), CLOG(2,1), CLOG(2,2)

      END SUBROUTINE Coulomb_log
!==============================================================
      SUBROUTINE fp_continue(ierr)

      USE fpnfrr
      USE libmtx
      USE fpnflg
      IMPLICIT NONE
      integer :: ierr,NSA,NSB,NS,NR,NP,NTH,NSBA,N,NSW,j
      INTEGER:: NSEND, NSWI
      real:: gut1, gut2, gut_prep

      IF(NRANK.eq.0) &
      WRITE(6,*) "----- SET COEFFICIENTS AND DISTRIBUTION FUNCTIONS -----"

      N_IMPL=0
      IF(MODELS.eq.3) CALL NF_LG_FUNCTION
      IF(MODELS.ne.0) CALL NF_REACTION_COEF
      CALL fusion_source_init

      DO NSA=NSASTART,NSAEND
         CALL FP_COEF(NSA)
         NSBA=NSB_NSA(NSA)
         DO NR=NRSTART,NREND
            DO NP=NPSTARTW,NPENDWM
               DO NTH=1,NTHMAX
                  F(NTH,NP,NR)=FNSP(NTH,NP,NR,NSBA)
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
      IF(NTG1.eq.0.and.MODEL_WAVE.ne.0) CALL FPWAVE_CONST ! all nrank must have RPWT  

      ISAVE=0

      IERR=0
 
      RETURN
      END SUBROUTINE fp_continue
!-------------------------------------------------------------
      SUBROUTINE fp_set_initial_value_from_f

      USE libmtx
      IMPLICIT NONE
      INTEGER:: NR, NSB

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
      CALL mtx_broadcast_real8(RT_T,NRMAX*NSAMAX)
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

      Implicit none

      integer :: ierr,NSA,NSB,NS,NR,NP,NTH,NSBA,N,NSW,j,i
      INTEGER:: NSEND, NSWI
      real:: gut1, gut2, gut_prep
      real(8):: alp, z_i, h_alpha_z, lambda_alpha, gamma_alpha_z, G_conner
      real(8):: G_conner_nr, G_conner_lm, SIGMA
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
!      WRITE(6,*) "START MESH"
      CALL fp_mesh(ierr)
!      WRITE(6,*) "END MESH"
!     ----- Initialize diffusion coef. -----
      call FPCINI
!     ----- set parameters for target species -----
      CALL fp_set_normalize_param
!     ----- Initialize velocity distribution function of all species -----
      CALL FNSP_INIT
      CALL FNSP_INIT_EDGE
!      WRITE(6,*) "END INIT"
!     ----- normalize bounce average parameter ---------
!      CALL fp_set_bounceaverage_param ! RCOEF
!     ----- set background f
      CALL mtx_set_communicator(comm_nsa)
      CALL update_fnsb
      CALL mtx_reset_communicator
!     ----- set parallel electric field -----
      IF(MODEL_DISRUPT.eq.0)THEN
         DO NR=1,NRMAX
            E1(NR)=E0!*E_drei0(1)
         END DO
         DO NR=NRSTART,NREND
            EP(NR)=E1(NR) ! plus
            EM(NR)=0.D0 ! minus
         END DO
         DO NR=NRSTART,NREND
            allocate(conduct_temp(NRSTART:NREND))
            CALL SPITZER_SIGMA(NR,SIGMA)
            conduct_temp(NR)=sigma
         END DO
         CALL mtx_set_communicator(comm_nr)
         call mtx_allgather_real8(conduct_temp,NREND-NRSTART+1,conduct_sp)
         CALL mtx_reset_communicator
      ELSE
         DO NR=NRSTART,NREND
            allocate(conduct_temp(NRSTART:NREND))
            allocate(E1_temp(NRSTART:NREND))
            CALL SPITZER_SIGMA(NR,SIGMA)
            conduct_temp(NR)=sigma
            IF(MODELE.eq.0)THEN
               E1_temp(NR)=RJ_ohm(NR)/SIGMA*1.D6 ! fit to initial current prof
            ELSEIF(MODELE.eq.1)THEN
               E1_temp(NR)=E0*E_drei0(1) ! uniform
            ELSEIF(MODELE.eq.2)THEN
               E1_temp(NR)=E0*(1.D0-RM(NR)**1.5)**1 ! arbitrary profile
            END IF
         END DO
         CALL mtx_set_communicator(comm_nr)
         call mtx_allgather_real8(conduct_temp,NREND-NRSTART+1,conduct_sp)
         call mtx_allgather_real8(E1_temp,NREND-NRSTART+1,E1)
         CALL mtx_reset_communicator
      END IF

!continue start
      CALL fp_continue(ierr)
      CALL fp_set_initial_value_from_f

      CALL GUTIME(gut2)
      gut_prep=gut2-gut1
      IF(NRANK.eq.0) WRITE(6,'(A,E14.6)') "---------------PREP_TIME=", gut_prep
 
      RETURN
      END subroutine fp_prep
!-----
      END MODULE fpprep
