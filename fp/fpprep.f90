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
         QLM(NR)=QL
         BT=BB
         BP= RSRHON(RHON)*BT/(RR*QL)
         EPSRM(NR)=RSRHON(RHON)/RR
         BPM(NR)= RSRHON(RHON)*BT/(RR*QL)
      ENDDO
!      RHON=RG(NRMAX+1)
      RHON=RM(NRMAX)+DELR
      CALL pl_qprf(RHON,QL)
      QLM(NRMAX+1)=QL
      BT=BB
      BP= RSRHON(RHON)*BT/(RR*QL)
      EPSRM(NRMAX+1)=RSRHON(RHON)/RR
      BPM(NRMAX+1)= RSRHON(RHON)*BT/(RR*QL)

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
            ITLG(NR)=0
            ITUG(NR)=0
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
            IF(nprocs.gt.1.and.NRANK.eq.1) &
                 WRITE(6,'(A,3I5,1P2E12.4)') 'NR,ITL,ITU,EPSRM=',NR,ITL(NR),ITU(NR),EPSRM(NR),EPSL
            EPSRM2(NR) = EPSRM(NR)
            EPSRM(NR)=EPSL
!            EPSRM2(NR) = EPSRM(NR)
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
            IF(nprocs.gt.1.and.NRANK.eq.1) &
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

      DO NTH=1,NTHMAX
         DO NR=NRSTART,NRENDX
            work(NR)=RLAMDA(NTH,NR)
         ENDDO
         CALL mtx_gatherv_real8_local(work(NRSTART:NRENDX),MTXLEN(NRANK+1), &
                                workg,NRMAX,MTXLEN,MTXPOS,ncoms)
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
         CALL mtx_gatherv_real8_local(work(NRSTART:NRENDX),MTXLEN(NRANK+1), &
                                workg,NRMAX,MTXLEN,MTXPOS,ncoms)
         CALL mtx_broadcast_real8(workg,NRMAX)
         DO NR=1,NRMAX
            ETAMG(NTH,NR)=workg(NR)
         ENDDO
      ENDDO

      DO NTH=1,NTHMAX
         DO NR=NRSTART,NRENDX
            work(NR)=RLAMDA_G(NTH,NR)
         ENDDO
         CALL mtx_gatherv_real8_local(work(NRSTART:NRENDX),MTXLEN(NRANK+1), &
                                workg,NRMAX,MTXLEN,MTXPOS,ncoms)
         CALL mtx_broadcast_real8(workg,NRMAX)
         DO NR=1,NRMAX
            RLAMDA_GG(NTH,NR)=workg(NR)
         ENDDO
      ENDDO

      DO NTH=1,NTHMAX
         DO NR=NRSTART,NRENDX
            work(NR)=ETAM_G(NTH,NR)
         ENDDO
         CALL mtx_gatherv_real8_local(work(NRSTART:NRENDX),MTXLEN(NRANK+1), &
                                workg,NRMAX,MTXLEN,MTXPOS,ncoms)
         CALL mtx_broadcast_real8(workg,NRMAX)
         DO NR=1,NRMAX
            ETAM_GG(NTH,NR)=workg(NR)
         ENDDO
      ENDDO

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

      IF(NRANK.eq.0)THEN
      DO NR=1,NRMAX
         WRITE(*,*) NR, "NRG", ITLG(NR)
         WRITE(*,*) NR, "NRM ", ITL(NR)
      END DO
      END IF

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

! ============================================================

      SUBROUTINE fp_prep(ierr)

      USE plprof
      USE fpnfrr
      USE libmtx

      Implicit none

      integer :: ierr,NSA,NSB,NS,NR,NP,NTH,NSFP,NSFD,NSBA,N,NREND1,NSW,j
      real(kind8) :: FL, RSUM1, RSUM2, RTFD0L, RHON, RNE, RTE
      real(kind8) :: RLNRL, FACT, RSUM, RSUM11, rsum3, rsum4, rsum5, rsum6
      TYPE(pl_plf_type),DIMENSION(NSMAX):: PLF
      INTEGER,DIMENSION(nprocs):: ima1,ima2,nra1,nra2,nma1,nma2,insa1,insa2
      real(kind8),DIMENSION(:),POINTER:: work,workg

!     ----- Initialize time counter -----

      TIMEFP=0.D0
      NTG1=0
      NTG2=0
!     ----- Check nprocs -----
!      IF(nprocs.GT.nrmax) THEN
!         IF(nrank.EQ.0) THEN
!            WRITE(6,*) 'XX fp_prep: nrmax must be greater than nprocs.'
!            WRITE(6,*) 'XX          nrmax,nprocs=',nrmax,nprocs
!         ENDIF
!         ierr=1
!         RETURN
!      ENDIF
!      N_partition_s=2
!      N_partition_r=25
      IF(N_partition_s*N_partition_r.ne.NPROCS)THEN
         IF(NRANK.eq.0) THEN
            WRITE(6,*) 'XX fp_prep: The partition number must be N_partition_s*N_partition_r=NPROCS.'
            WRITE(6,'(A,3I4)') 'XX N_partition_s, N_partiton_r, NPROCS=', N_partition_s, N_partition_r, NPROCS
         END IF
         ierr=1
         RETURN
      END IF

!     ----- Set matrix size -----
      call mtx_comm_split_s(N_partition_s,colors,keys,ncoms)
      call mtx_comm_split_r(N_partition_r,colorr,keyr,ncomr)
      NSASTART = (NSAMAX/N_partition_s)*colors+1
      NSAEND =   (NSAMAX/N_partition_s)*(colors+1)

!      WRITE(*,'(A,5I4)') "NRANK, colors, NRANKS, colorr, NRANKR", &
!           nrank, colors, nranks, colorr, nrankr

      imtxsize=nthmax*npmax*nrmax
      IF(modeld.EQ.0) THEN
         imtxwidth=4*nthmax-1
      ELSE
         imtxwidth=4*nthmax*npmax-1
      ENDIF
      CALL mtx_setup(imtxsize,imtxstart,imtxend,imtxwidth,ncoms)
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
      CALL mtx_gather_integer(nsastart,insa1)
      CALL mtx_gather_integer(nsaend,  insa2)
      IF(nrank.EQ.0) THEN
         write(6,'(A,2I10)') '  imtxsize,imtxwidth=',imtxsize,imtxwidth
         write(6,'(A,A)') '     nrank   imtxstart   imtxend   nrstart',&
                          '     nrend   nmstart     nmend      nsastart      nsaend'
         DO N=1,nprocs
            write(6,'(9I10)') N,ima1(N),ima2(N),nra1(N),nra2(N),nma1(N),nma2(N),insa1(N),insa2(N)
         ENDDO
      ENDIF
      CALL mtx_cleanup

!     ----- Allocate variables -----

      CALL fp_allocate
      call fp_allocate_ntg1
      call fp_allocate_ntg2

!     ----- Get mtxlen and mtxpos -----
!     MTXLEN(NRANK+1): the number of NR grid points for each RANK
!     MTXPOS(NRANK):
      CALL mtx_allgather_integer(nrend-nrstart+1,mtxlen)
      CALL mtx_allgather_integer(nrstart-1,mtxpos)
      CALL mtx_allgather_integer(nrend-nrstart+1,savlen)
      NSW=NSAEND-NSASTART+1
      DO N=1,NSW
         NSA=N+NSASTART-1
         CALL mtx_allgather_integer_sav( (nsa-1)*NRMAX+NRSTART,N,NSW,ncomw)
      END DO

      if(nrank.eq.0) then
!         DO N=1,NPROCS
!            WRITE(6,'(A,5I8)') '  nrank,mtxpos,mtxlen = ', &
!                                   n-1,mtxpos(n),mtxlen(n)
!         ENDDO
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

!      DO NSB=1,NSBMAX
      DO NSB=NSASTART,NSAEND
         NS=NS_NSB(NSB)
         DO NR=NRSTART-1,NREND+1
            IF(NR.ge.1.and.NR.le.NRMAX)THEN
               DO NP=1,NPMAX
                  FL=FPMXWL(PM(NP,NSB),NR,NS)
                  DO NTH=1,NTHMAX
                     FNSP(NTH,NP,NR,NSB)=FL
                  END DO
               ENDDO
            END IF
         END DO
         NR=NRMAX+1
         DO NP=1,NPMAX
            FL=FPMXWL(PM(NP,NSB),NR,NS)
            DO NTH=1,NTHMAX
               FS3(NTH,NP,NSB)=FL ! rho=1.0
            END DO
         ENDDO
      END DO

!--------- normalize bounce average parameter ---------

      IF(MODELA.eq.1)THEN
         NSB=NSASTART ! arbitrary
         DO NR=NRSTART,NREND
            RSUM1=0.D0
            RSUM2=0.D0
            RSUM3=0.D0
            RSUM4=0.D0
            DO NP=1,NPMAX
               DO NTH=1,NTHMAX
                  RSUM1 = RSUM1+VOLP(NTH,NP,NSB)*RLAMDAG(NTH,NR)/RFSADG(NR)*FNSP(NTH,NP,NR,NSB)
                  RSUM2 = RSUM2+VOLP(NTH,NP,NSB)*FNSP(NTH,NP,NR,NSB)
                  RSUM3 = rsum3+VOLP(NTH,NP,NSB)*RLAMDA_GG(NTH,NR)/RFSAD_GG(NR)
                  RSUM4 = rsum4+VOLP(NTH,NP,NSB)
               END DO
            END DO
            IF(RSUM1.EQ.0.D0) &
                 WRITE(6,'(1P3E12.4)') VOLP(1,1,NSB),FNSP(1,1,1,NSB),RLAMDA(1,1)

            RCOEFN(NR)=RSUM2/RSUM1
            RCOEFN_G(NR)=RSUM4/RSUM3
!            RCOEF(NR) = ( QLM(NR)*RR )
!            RCOEF_G(NR) = ( QLG(NR)*RR )
         END DO
      ELSE
         DO NR=NRSTART,NREND
!            RCOEF(NR)=1.D0
!            RCOEF_G(NR)=1.D0
            RCOEFN(NR)=1.D0
            RCOEFN_G(NR)=1.D0
         ENDDO
      END IF

!     ----- set boundary distribution functions -----
      
!      DO NSA=1,NSAMAX
      DO NSA=NSASTART,NSAEND
         NS=NS_NSA(NSA)
         NSBA=NSB_NSA(NSA)
         DO NP=1,NPMAX
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
         NSBA=NSB_NSA(NSA)
         DO NP=1,NPMAX
            CALL FPMXWL_EDGE(NP,NS,FL)
            DO NTH=1,NTHMAX
               FS2(NTH,NP,NS)=FL ! at R=1.0+DELR/2
            ENDDO
         ENDDO
      ENDDO
!
      IF(MODELA.eq.1)THEN
         NSBA=NSASTART
         RSUM1=0.D0
         RSUM2=0.D0
         RSUM3=0.D0
         RSUM4=0.D0
         RSUM5=0.D0
         RSUM6=0.D0
         DO NP=1,NPMAX
            DO NTH=1,NTHMAX
               RSUM1 = RSUM1+VOLP(NTH,NP,NSBA)*RLAMDA_GG(NTH,1) &
                    /RFSAD_GG(1)*FS1(NTH,NP,NSBA)
               RSUM2 = RSUM2+VOLP(NTH,NP,NSBA)*FS1(NTH,NP,NSBA)
               RSUM3 = RSUM3+VOLP(NTH,NP,NSBA)*RLAMDA_GG(NTH,NRMAX+1) &
                    /RFSAD_GG(NRMAX+1)*FS3(NTH,NP,NSBA)
               RSUM4 = RSUM4+VOLP(NTH,NP,NSBA)*FS3(NTH,NP,NSBA)
               RSUM5 = RSUM5+VOLP(NTH,NP,NSBA)*RLAMDAG(NTH,NRMAX+1) &
                    /RFSADG(NRMAX+1)*FS2(NTH,NP,NSBA)
               RSUM6 = RSUM6+VOLP(NTH,NP,NSBA)*FS2(NTH,NP,NSBA)
            END DO
         END DO
         RCOEFN_GG(1) = RSUM2/RSUM1
         RCOEFN_GG(NRMAX+1) = RSUM4/RSUM3
         RCOEFNG(NRMAX+1) = RSUM6/RSUM5
      ELSE
         RCOEFN_GG(1) = 1.D0
         RCOEFN_GG(NRMAX+1) = 1.D0
         RCOEFNG(NRMAX+1) = 1.D0
      END IF
!!!!
      allocate(work(nrstart:nrendx),workg(NRMAX))

!      DO NR=NRSTART,NRENDX
!         work(NR)=RCOEF(NR)
!      ENDDO
!      CALL mtx_gatherv_real8_local(work(NRSTART:NRENDX),MTXLEN(NRANK+1), &
!           workg,NRMAX,MTXLEN,MTXPOS,ncoms)
!      CALL mtx_broadcast_real8(workg,NRMAX)
!      DO NR=1,NRMAX
!         RCOEFG(NR)=workg(NR)
!      ENDDO
!      RCOEFG(NRMAX+1)=(QLM(NRMAX)*RR)

!      DO NR=NRSTART,NRENDX
!         work(NR)=RCOEF_G(NR)
!      ENDDO
!      CALL mtx_gatherv_real8_local(work(NRSTART:NRENDX),MTXLEN(NRANK+1), &
!           workg,NRMAX,MTXLEN,MTXPOS,ncoms)
!      CALL mtx_broadcast_real8(workg,NRMAX)
!      DO NR=1,NRMAX
!         RCOEF_GG(NR)=workg(NR)
!      ENDDO
!      RCOEF_GG(NRMAX+1)=1.D0/(QLG(NRMAX+1)*RR)

      DO NR=NRSTART,NRENDX
         work(NR)=RCOEFN(NR)
      ENDDO
      CALL mtx_gatherv_real8_local(work(NRSTART:NRENDX),MTXLEN(NRANK+1), &
           workg,NRMAX,MTXLEN,MTXPOS,ncoms)
      CALL mtx_broadcast_real8(workg,NRMAX)
      DO NR=1,NRMAX
         RCOEFNG(NR)=workg(NR)
      ENDDO

      DO NR=NRSTART,NRENDX
         work(NR)=RCOEFN_G(NR)
      ENDDO
      CALL mtx_gatherv_real8_local(work(NRSTART:NRENDX),MTXLEN(NRANK+1), &
           workg,NRMAX,MTXLEN,MTXPOS,ncoms)
      CALL mtx_broadcast_real8(workg,NRMAX)
      DO NR=2,NRMAX
         RCOEFN_GG(NR)=workg(NR)
      ENDDO

!      IF(NRANK.eq.0)THEN
!         open(8,file='rcoefng_norfsad.dat')
!         DO NR=1,NRMAX
!            WRITE(8,'(7E14.6)') RM(NR), RCOEFNG(NR), RFSADG(NR) &
!                 , RG(NR), RCOEFN_GG(NR), RFSAD_GG(NR), QLM(NR)
!         END DO
!         close(8)
!         open(8,file='volp_r_tpb_ex_killeen_fine.dat')
!         DO NTH=1,NTHMAX/2
!            DO NR=1,NRMAX
!               WRITE(8,'(2I4,E14.6)') NTH, NR, RLAMDAG(NTH,NR)/RFSADG(NR)
!            END DO
!            WRITE(8,*) " "
!            WRITE(8,*) " "
!         END DO
!         close(8)
!         open(8,file='r_ram_ram_ex_killeen_fine.dat')
!         DO NTH = 1, NTHMAX
!            DO NR=1,NRMAX
!               WRITE(8,'(2I4,2E14.6)') NTH, NR, RLAMDAG(NTH,NR), RFSADG(NR)
!            END DO
!            WRITE(8,*) " "
!            WRITE(8,*) " "
!         END DO
!         close(8)
!      END IF

      deallocate(work,workg)

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

      N_IMPL=0
      CALL NF_REACTION_COEF
      NCALCNR=0
!      DO NSA=1,NSAMAX
      CALL fusion_source_init
      DO NSA=NSASTART,NSAEND
         CALL FP_COEF(NSA)
         NSBA=NSB_NSA(NSA)
         DO NR=NRSTART,NREND
            DO NP=1,NPMAX
               DO NTH=1,NTHMAX
                  F(NTH,NP,NR)=FNSP(NTH,NP,NR,NSBA)
               END DO
            END DO
         END DO
         CALL FPWEIGHT(NSA,IERR)
      END DO
      CALL source_allreduce(SPPF,ncomr)
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
!-----
      END MODULE fpprep
