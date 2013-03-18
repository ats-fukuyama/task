!     $Id$

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
!         write(6,'(A,I5,1P4E12.4)') 'nr,rm,epsrm,bpm=', &
!              NR,RM(NR),RSRHON(RHON),EPSRM(NR),BPM(NR)
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
      DO NR=NRSTART,NREND
         DO NTH=1,NTHMAX
            PSIPM_P(NTH,NR)=0.D0
         END DO
      END DO
      DO NR=NRSTART,NREND
         DO NTH=1,NTHMAX+1
            PSIPG_P(NTH,NR)=0.D0
         END DO
      END DO
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

      DO NSB=1,NSBMAX
!      DO NP=1,NPMAX
      DO NP=NPSTART,NPEND
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
            ITLG(NR)=0
            ITUG(NR)=0
            EPSRM2(NR) = EPSRM(NR)
            EPSRG2(NR) = EPSRG(NR)
         ENDDO
      ELSE
         DO NR=1,NRMAX+1
            A1=ACOS(SQRT(2.D0*EPSRM(NR)/(1.D0+EPSRM(NR))))
            NTH=0
            DO WHILE (THG(NTH+1).le.A1)
               NTH = NTH+1
            END DO
            ITL(NR)=NTH
            ITU(NR)=NTHMAX-NTH+1

            EPSL=COSM(ITL(NR))**2/(2.D0-COSM(ITL(NR))**2)

!            WRITE(6,'(A,3I5,1P2E12.4)') 'NR,ITL,ITU,EPSRM=', &
!                          NR,ITL(NR),ITU(NR),EPSRM(NR),EPSL
            EPSRM2(NR) = EPSRM(NR)
            EPSRM(NR)=EPSL
         ENDDO

         IF(NRANK.eq.1) WRITE(6,*) " "

         DO NR=1,NRMAX+1
            A1=ACOS(SQRT(2.D0*EPSRG(NR)/(1.D0+EPSRG(NR))))
            NTH=0
            DO WHILE (THG(NTH+1).le.A1)
               NTH = NTH+1
            END DO
            ITLG(NR)=NTH
            ITUG(NR)=NTHMAX-NTH+1

            EPSL=COSM(ITLG(NR))**2/(2.D0-COSM(ITLG(NR))**2)
            IF(nsize.gt.1.and.NRANK.eq.1) &
                 WRITE(6,'(A,2I5,1P2E12.4)') 'NR,NTHC,EPSRG=',NR,NTH,EPSRG(NR),EPSL
            EPSRG2(NR) = EPSRG(NR)
            EPSRG(NR)=EPSL
!            EPSRG2(NR) = EPSRG(NR)
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
         DO NR=NRSTART,NREND
            CALL SET_BOUNCE_PARAM(NR)
         END DO
         CALL SET_BOUNCE_PARAM(NRMAX+1)
         DO NR=1,NRMAX+1
            CALL SET_RFSAD(NR)
         END DO
      END IF ! MODELA

      allocate(work(nrstart:nrendx),workg(NRMAX))

      CALL mtx_set_communicator(comm_nr)
      DO NTH=1,NTHMAX
         DO NR=NRSTART,NRENDX
            work(NR)=RLAMDA(NTH,NR)
         ENDDO
         CALL mtx_allgatherv_real8(work(NRSTART:NRENDX),MTXLEN(NRANK+1), &
                                   workg,NRMAX,MTXLEN,MTXPOS)
         DO NR=1,NRMAX
            RLAMDAG(NTH,NR)=workg(NR)
         ENDDO
         RLAMDAG(NTH,NRMAX+1)=RLAMDAG(NTH,NRMAX)
      ENDDO

      DO NTH=1,NTHMAX
         DO NR=NRSTART,NRENDX
            work(NR)=ETAM(NTH,NR)
         ENDDO
         CALL mtx_allgatherv_real8(work(NRSTART:NRENDX),MTXLEN(NRANK+1), &
                                   workg,NRMAX,MTXLEN,MTXPOS)
         DO NR=1,NRMAX
            ETAMG(NTH,NR)=workg(NR)
         ENDDO
      ENDDO

      DO NTH=1,NTHMAX
         DO NR=NRSTART,NRENDX
            work(NR)=RLAMDA_G(NTH,NR)
         ENDDO
         CALL mtx_allgatherv_real8(work(NRSTART:NRENDX),MTXLEN(NRANK+1), &
                                   workg,NRMAX,MTXLEN,MTXPOS)
         DO NR=1,NRMAX
            RLAMDA_GG(NTH,NR)=workg(NR)
         ENDDO
      ENDDO

      DO NTH=1,NTHMAX
         DO NR=NRSTART,NRENDX
            work(NR)=ETAM_G(NTH,NR)
         ENDDO
         CALL mtx_allgatherv_real8(work(NRSTART:NRENDX),MTXLEN(NRANK+1), &
                                   workg,NRMAX,MTXLEN,MTXPOS)
         DO NR=1,NRMAX
            ETAM_GG(NTH,NR)=workg(NR)
         ENDDO
      ENDDO
      CALL mtx_reset_communicator

!      IF(NRANK.eq.0)THEN
!      open(8,file='RLAMDAG100_tpb_ex_killeen_fine.dat')
!      DO NR =1, NRMAX
!      DO NTH=1,NTHMAX
!         WRITE(8,'(2I4, 4E14.6)') NR, NTH, NTH-0.5D0, COSM(NTH), RLAMDAG(NTH,NR), RLAMDA_GG(NTH,NR)
!      END DO
!      WRITE(8,*) " "
!      WRITE(8,*) " "
!      END DO
!      close(8)
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
      integer:: NSA,NR, NTH, NP

      ISAVE=0

      DO NSA=1,NSAMAX
      DO NR=NRSTART,NREND
!         DO NP=1,NPMAX+1
         DO NP=NPSTART,NPENDWG
         DO NTH=1,NTHMAX
            DPP(NTH,NP,NR,NSA)=0.D0
            DPT(NTH,NP,NR,NSA)=0.D0
            FPP(NTH,NP,NR,NSA)=0.D0
         ENDDO
         ENDDO
!         DO NP=1,NPMAX
         DO NP=NPSTARTW,NPENDWM
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

! ============================================================
      SUBROUTINE fp_comm_setup

      USE libmtx
      IMPLICIT NONE
      INTEGER:: ierr, colors, keys, N, NREND1, NPEND1
      INTEGER,DIMENSION(nsize):: ima1,ima2,npa1,npa2,nra1,nra2,nma1,nma2,insa1,insa2

!     ----- Check nsize -----
      IF(N_partition_s*N_partition_r*N_partition_p.ne.nsize)THEN
         IF(NRANK.eq.0) THEN
            WRITE(6,*) 'XX fp_prep: N_partition_s*N_partition_r != nsize.'
            WRITE(6,'(A,3I4)') 'XX N_partition_s, N_partiton_r, nsize=', &
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
      CALL mtx_reset_communicator

      nrstart=(imtxstart-1)/(nthmax*npmax)+1
      nrend=  (imtxend  -1)/(nthmax*npmax)+1
      nrend1= (imtxend    )/(nthmax*npmax)+1
!      IF(nrend1.EQ.nrend) THEN
!         NRENDX=NREND-1
!      ELSE
         NRENDX=NREND
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


      IF(nrank.EQ.0) THEN
         write(6,'(A,2I10)') '  imtxsize,imtxwidth=',imtxsize,imtxwidth
         write(6,'(A,A,A)') '     nrank   imtxstart   imtxend   npstart    npend', &
                          '    nrstart     nrend   nmstart     nmend', &
                          '      nsastart      nsaend'
         DO N=1,nsize
            write(6,'(11I10)') N,ima1(N),ima2(N),npa1(N),npa2(N),nra1(N),nra2(N), &
                              nma1(N),nma2(N),insa1(N),insa2(N)
         ENDDO
      ENDIF

!      IF(NRANK.eq.0) WRITE(6,*) "      NRANK, nsa_rank, nr_rank, np_rank, nrnp_rank, nsanp_rank, nsanr_rank"
!      write(6,'(7I10)') NRANK, comm_nsa%rankl, comm_nr%rankl, &
!           comm_np%rankl, comm_nrnp%rankl, comm_nsanp%rankl, comm_nsanr%rankl

      CALL mtx_cleanup


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



      DO NSB=NSASTART,NSAEND
         NS=NS_NSB(NSB)
         DO NR=NRSTARTW,NRENDWM
            IF(NR.ge.1.and.NR.le.NRMAX)THEN
!               DO NP=1,NPMAX
               DO NP=NPSTARTW,NPENDWM
                  FL=FPMXWL(PM(NP,NSB),NR,NS)
                  DO NTH=1,NTHMAX
                     FNSP(NTH,NP,NR,NSB)=FL
                     FNS0(NTH,NP,NR,NSB)=FL
                  END DO
               ENDDO
            END IF
         END DO
         NR=NRMAX+1
!         DO NP=1,NPMAX
         DO NP=NPSTARTW,NPENDWM
            FL=FPMXWL(PM(NP,NSB),NR,NS)
            DO NTH=1,NTHMAX
               FS3(NTH,NP,NSB)=FL ! rho=1.0
            END DO
         ENDDO
      END DO
!     ----- set boundary distribution functions -----
      
!      DO NSA=1,NSAMAX
      DO NSA=NSASTART,NSAEND
         NS=NS_NSA(NSA)
         NSBA=NSB_NSA(NSA)
!         DO NP=1,NPMAX
         DO NP=NPSTARTW,NPENDWM
            FL=FPMXWL(PM(NP,NSBA),0,NS)
            DO NTH=1,NTHMAX
               FS1(NTH,NP,NSA)=FL ! at r=0
            ENDDO
         ENDDO
      END DO
!
!      DO NSA=1,NSAMAX
      DO NSA=NSASTART,NSAEND
         NS=NS_NSA(NSA)
         DO NP=NPSTARTW,NPENDWM
            CALL FPMXWL_EDGE(NP,NSA,FL)
            DO NTH=1,NTHMAX
               FS2(NTH,NP,NS)=FL ! at R=1.0+DELR/2
            ENDDO
         ENDDO
      ENDDO

      END SUBROUTINE FNSP_INIT
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
!            DO NP=NPSTART,NPEND
!               DO NTH=1,NTHMAX
!                  RSUM1 = RSUM1+VOLP(NTH,NP,NSB)*RLAMDAG(NTH,NR)/RFSADG(NR)*FNSP(NTH,NP,NR,NSB)
!                  RSUM2 = RSUM2+VOLP(NTH,NP,NSB)*FNSP(NTH,NP,NR,NSB)
!                  RSUM3 = rsum3+VOLP(NTH,NP,NSB)*RLAMDA_GG(NTH,NR)/RFSAD_GG(NR)
!                  RSUM4 = rsum4+VOLP(NTH,NP,NSB)
!               END DO
!            END DO
            DO NTH=1,NTHMAX
               RSUM1 = RSUM1+RLAMDAG(NTH,NR)/RFSADG(NR)*SINM(NTH)*DELTH
               RSUM2 = RSUM2+SINM(NTH)*DELTH
               RSUM3 = rsum3+RLAMDA_GG(NTH,NR)/RFSAD_GG(NR)*SINM(NTH)*DELTH
               RSUM4 = rsum4+SINM(NTH)*DELTH
            END DO
            IF(RSUM1.EQ.0.D0) &
                 WRITE(6,'(1P3E12.4)') VOLP(1,1,NSB),FNSP(1,1,1,NSB),RLAMDA(1,1)
!            CALL p_theta_integration(RSUM1)
!            CALL p_theta_integration(RSUM2)
!            CALL p_theta_integration(RSUM3)
!            CALL p_theta_integration(RSUM4)

            RCOEFN(NR)=RSUM2/RSUM1
            RCOEFN_G(NR)=RSUM4/RSUM3
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
!         DO NP=NPSTART,NPEND
!            DO NTH=1,NTHMAX
!               RSUM1 = RSUM1+VOLP(NTH,NP,NSBA)*RLAMDA_GG(NTH,1) &
!                    /RFSAD_GG(1)*FS1(NTH,NP,NSBA)
!               RSUM2 = RSUM2+VOLP(NTH,NP,NSBA)*FS1(NTH,NP,NSBA)
!               RSUM3 = RSUM3+VOLP(NTH,NP,NSBA)*RLAMDA_GG(NTH,NRMAX+1) &
!                    /RFSAD_GG(NRMAX+1)*FS3(NTH,NP,NSBA)
!               RSUM4 = RSUM4+VOLP(NTH,NP,NSBA)*FS3(NTH,NP,NSBA)
!               RSUM5 = RSUM5+VOLP(NTH,NP,NSBA)*RLAMDAG(NTH,NRMAX+1) &
!                    /RFSADG(NRMAX+1)*FS2(NTH,NP,NSBA)
!               RSUM6 = RSUM6+VOLP(NTH,NP,NSBA)*FS2(NTH,NP,NSBA)
!            END DO
!         END DO
!         CALL p_theta_integration(RSUM1)
!         CALL p_theta_integration(RSUM2)
!         CALL p_theta_integration(RSUM3)
!         CALL p_theta_integration(RSUM4)
!         CALL p_theta_integration(RSUM5)
!         CALL p_theta_integration(RSUM6)
         DO NTH=1,NTHMAX
            RSUM1 = RSUM1+RLAMDA_GG(NTH,1)/RFSAD_GG(1)*SINM(NTH)*DELTH
            RSUM2 = RSUM2+SINM(NTH)*DELTH
            RSUM3 = RSUM3+RLAMDA_GG(NTH,NRMAX+1)/RFSAD_GG(NRMAX+1)*SINM(NTH)*DELTH
!            RSUM4 = RSUM4+SINM(NTH)*DELTH
            RSUM5 = RSUM5+RLAMDAG(NTH,NRMAX+1)/RFSADG(NRMAX+1)*SINM(NTH)*DELTH
!            RSUM6 = RSUM6+SINM(NTH)*DELTH
         END DO

         RCOEFN_GG(1) = RSUM2/RSUM1
         RCOEFN_GG(NRMAX+1) = RSUM2/RSUM3
         RCOEFNG(NRMAX+1) = RSUM2/RSUM5
!         RCOEFN_GG(NRMAX+1) = RSUM4/RSUM3
!         RCOEFNG(NRMAX+1) = RSUM6/RSUM5
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

      allocate(work(nrstart:nrendx),workg(NRMAX))

      CALL mtx_set_communicator(comm_nr)

      DO NR=NRSTART,NRENDX
         work(NR)=RCOEFN(NR)
      ENDDO
      CALL mtx_allgatherv_real8(work(NRSTART:NRENDX),MTXLEN(NRANK+1), &
                                workg,NRMAX,MTXLEN,MTXPOS)
      DO NR=1,NRMAX
         RCOEFNG(NR)=workg(NR)
      ENDDO

      DO NR=NRSTART,NRENDX
         work(NR)=RCOEFN_G(NR)
      ENDDO
      CALL mtx_allgatherv_real8(work(NRSTART:NRENDX),MTXLEN(NRANK+1), &
                                workg,NRMAX,MTXLEN,MTXPOS)
      DO NR=2,NRMAX
         RCOEFN_GG(NR)=workg(NR)
      ENDDO

      DO NR=NRSTART,NRENDX
         work(NR)=RCOEFJ(NR)
      ENDDO
      CALL mtx_allgatherv_real8(work(NRSTART:NRENDX),MTXLEN(NRANK+1), &
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
      INTEGER:: NSA, NSB, NS, NSBA, NSFP, NSFD, NR
      TYPE(pl_plf_type),DIMENSION(NSMAX):: PLF
      real(kind8):: RTFD0L, RHON, RNE, RTE, RLNRL, FACT

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


      END SUBROUTINE fp_set_normalize_param
!==============================================================
      SUBROUTINE fp_continue(ierr)

      USE fpnfrr
      USE libmtx
      IMPLICIT NONE
      integer :: ierr,NSA,NSB,NS,NR,NP,NTH,NSBA,N,NSW,j
      INTEGER:: NSEND, NSWI
      real:: gut1, gut2, gut_prep

      CALL GUTIME(gut1)
      IF(NRANK.eq.0) &
      WRITE(6,*) "----- RESET COEFFICIENTS FOR NEW PARAMETERS -----"

!     ----- set parallel electric field -----
      DO NR=1,NRMAX
         E1(NR)=E0/(1.D0+EPSRM(NR))
         IF(MODELE.eq.0)THEN
            EP(NR)=0.D0 ! plus
            EM(NR)=0.D0 ! minus
         ELSEIF(MODELE.eq.1)THEN
            EP(NR)=E1(NR) ! plus
            EM(NR)=E1(NR) ! minus
         END IF
      ENDDO
      N_IMPL=0
      CALL NF_REACTION_COEF
      NCALCNR=0
      CALL fusion_source_init
      IF(NRANK.eq.0) WRITE(*,*) "source_init"
      DO NSA=NSASTART,NSAEND
         CALL FP_COEF(NSA)
         IF(NRANK.eq.0) WRITE(*,*) "coef",nsa
         CALL FPWEIGHT(NSA,IERR)
         IF(NRANK.eq.0) WRITE(*,*) "weight",nsa
      END DO
!      CALL mtx_set_communicator(comm_nr)
!      CALL source_allreduce(SPPF)
!      CALL mtx_reset_communicator
      ISAVE=0
!      IF(NTG1.eq.0) CALL FPWAVE_CONST ! all nrank must have RPWT  
      IERR=0
      CALL GUTIME(gut2)
      gut_prep=gut2-gut1
      IF(NRANK.eq.0) WRITE(6,'(A,E14.6)') "-------- CONTINUE_TIME=", gut_prep
 
      RETURN
      END SUBROUTINE fp_continue
!-------------------------------------------------------------
      SUBROUTINE fp_prep(ierr)

      USE plprof
      USE fpnfrr
      USE libmtx

      Implicit none

      integer :: ierr,NSA,NSB,NS,NR,NP,NTH,NSBA,N,NSW,j
      INTEGER:: NSEND, NSWI
      real:: gut1, gut2, gut_prep

      CALL GUTIME(gut1)
!     ----- Initialize time counter -----

      TIMEFP=0.D0
      NTG1=0
      NTG2=0
      gut_comm(:)=0.0

      CALL fp_comm_setup

!     ----- Allocate variables -----

      CALL fp_allocate
      call fp_allocate_ntg1
      call fp_allocate_ntg2

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
!     ----- Initialize velocity distribution function of all species -----
      CALL FNSP_INIT
!     ----- normalize bounce average parameter ---------
      CALL fp_set_bounceaverage_param
!     ----- set parameters for target species -----
      CALL fp_set_normalize_param
!     ----- set background f
      CALL mtx_set_communicator(comm_nsa)
      CALL update_fnsb
      CALL mtx_reset_communicator
!     ----- set parallel electric field -----
      DO NR=1,NRMAX
!         E1(NR)=E0/(1.D0+EPSRM(NR))
         E1(NR)=E0
         IF(MODELE.eq.0)THEN
            EP(NR)=0.D0 ! plus
            EM(NR)=0.D0 ! minus
         ELSEIF(MODELE.eq.1)THEN
            EP(NR)=E1(NR) ! plus
            EM(NR)=E1(NR) ! minus
         END IF
         RJ_M(NR)=0.D0
      ENDDO
!      IF(NRANK.eq.0)THEN
!         DO NR=1,NRMAX
!            WRITE(*,*) NR, RM(NR), E1(NR)
!         END DO
!      END IF

      EP_PHIM(:,:,:)=0.D0
      EP_PHIG(:,:,:)=0.D0

      N_IMPL=0
      CALL NF_REACTION_COEF
      NCALCNR=0
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
      CALL mtx_set_communicator(comm_nr)
!      CALL source_allreduce(SPPF)
      CALL mtx_reset_communicator
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
      CALL GUTIME(gut2)
      gut_prep=gut2-gut1
      IF(NRANK.eq.0) WRITE(6,'(A,E14.6)') "---------------PREP_TIME=", gut_prep
 
      RETURN
      END subroutine fp_prep
!-----
      END MODULE fpprep
