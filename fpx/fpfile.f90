!     $Id: fpfile.f90,v 1.2 2010/03/19 09:27:53 nuga Exp $
!
!     ***********************************************************
!
!           SAVE VELOCITY DISTRIBUTION DATA
!
!     ***********************************************************
!
      MODULE fpfile

      USE fpcomm

      contains
!----------------------------------------------------------------
      SUBROUTINE FP_SAVE

      USE libfio

      CALL FWOPEN(21,KNAMFP,0,MODEFW,'FP',IERR)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX FPSAVE: FWOPEN: IERR=',IERR
         RETURN
      ENDIF

      REWIND(21)
      WRITE(21) NRMAX,NPMAX,NTHMAX,NSAMAX
      WRITE(21) DELR,DELTH,RMIN,RMAX
      DO NSA=1,NSAMAX
         NS=NS_NSA(NSA)
         WRITE(21) NS_NSA(NSA)
         WRITE(21) DELP(NS)
         WRITE(21) AEFP(NSA),AMFP(NSA),RNFP0(NSA),RTFP0(NSA)
         WRITE(21) (((FNS(NTH,NP,NR,NSA),NTH=1,NTHMAX), &
                      NP=1,NPMAX),NR=1,NRMAX)
      ENDDO
      CLOSE(21)

      WRITE(6,*) '# DATA WAS SUCCESSFULLY SAVED TO THE FILE.'

  900 RETURN
      END SUBROUTINE FP_SAVE
!
!     ***********************************************************
!
!           LOAD VELOCITY DISTRIBUTION DATA
!
!     ***********************************************************
!
      SUBROUTINE FP_LOAD

      USE libfio

      CALL FROPEN(21,KNAMFP,0,MODEFR,'FP',IERR)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX FP_LOAD: FROPEN: IERR=',IERR
         RETURN
      ENDIF

      REWIND(21)
      READ(21) NRMAX,NPMAX,NTHMAX,NSAMAX
      READ(21) DELR,DELTH,RMIN,RMAX
      DO NSA=1,NSAMAX
         NS=NS_NSA(NSA)
         READ(21) NS_NSA(NSA)
         READ(21) DELP(NS)
         READ(21) AEFP(NSA),AMFP(NSA),RNFP0(NSA),RTFP0(NSA)
         READ(21) (((FNS(NTH,NP,NR,NSA),NTH=1,NTHMAX), &
                      NP=1,NPMAX),NR=1,NRMAX)
      ENDDO
      CLOSE(21)

      WRITE(6,*) '# DATA WAS SUCCESSFULLY LOADED FROM THE FILE.'

  900 RETURN
      END SUBROUTINE FP_LOAD
!------------------------------------------      
      SUBROUTINE FP_SAVE2

        USE libfio
      IMPLICIT NONE
      INTEGER:: NSA, NR, NP, NTH, NSB, NS, IERR

      CALL FWOPEN(21,KNAMFP,0,MODEFW,'FP',IERR)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX FPSAVE: FWOPEN: IERR=',IERR
         RETURN
      ENDIF

      REWIND(21)

      WRITE(21) NRMAX,NPMAX,NTHMAX,NSAMAX,NSBMAX,NSMAX
      WRITE(21) DELR,DELTH,RMIN,RMAX
      DO NSA=1,NSAMAX
         WRITE(21) NS_NSA(NSA)
         WRITE(21) AEFP(NSA),AMFP(NSA),RNFP0(NSA),RTFP0(NSA)
      ENDDO
      DO NS=1,NSMAX
         WRITE(21) DELP(NS),PMAX(NS),EMAX(NS)
      ENDDO

      DO NSA=1, NSAMAX
         DO NR=1, NRMAX
            DO NP=1, NPMAX
               WRITE(21) (FNS(NTH,NP,NR,NSA), NTH=1, NTHMAX) ! -> FNSP, FNSB
            END DO
         END DO
      END DO

! required for FP
! the values that update in each step

      WRITE(21) TIMEFP, NT_init

      WRITE(21) ( (RN_TEMP(NR,NS), NR=1,NRMAX), NS=1,NSMAX)
      WRITE(21) ( (RT_TEMP(NR,NS), NR=1,NRMAX), NS=1,NSMAX)
      WRITE(21) ( (RT_BULK(NR,NS), NR=1,NRMAX), NS=1,NSAMAX)
      WRITE(21) ( (RNS_DELF(NR,NS), NR=1,NRMAX), NS=1,NSMAX)
      WRITE(21) (E1(NR),NR=1,NRMAX) ! -> EP

      IF(MODELD.ne.0)THEN
         WRITE(21) ((((WEIGHR_G(NTH,NP,NR,NSA), NTH=1,NTHMAX),NP=1,NPMAX),NR=1, NRMAX+1),NSA=1,NSAMAX)
      END IF

      CLOSE(21)

      WRITE(6,*) '# DATA WAS SUCCESSFULLY SAVED TO THE FILE.'

!      DO NR=1, NRMAX
!         WRITE(*,'(A,I4,3E14.6)') "TEST 1 ", NR, RT_BULK(NR,2), RN_TEMP(NR,3), FNS(1,10,NR,2)
!      END DO

  900 RETURN
      END SUBROUTINE FP_SAVE2
!------------------------------------------      
      SUBROUTINE FP_LOAD2

        USE libfio
      IMPLICIT NONE
      INTEGER:: NSA, NR, NP, NTH, NSB, IERR, NS

      CALL FROPEN(21,KNAMFP,0,MODEFR,'FP',IERR)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX FP_LOAD2: FROPEN: IERR=',IERR
         RETURN
      ENDIF

! required for FP
! the values that update in each step
      REWIND(21)

      READ(21) NRMAX,NPMAX,NTHMAX,NSAMAX,NSBMAX,NSMAX
      READ(21) DELR,DELTH,RMIN,RMAX
      DO NSA=1,NSAMAX
         READ(21) NS_NSA(NSA)
         READ(21) AEFP(NSA),AMFP(NSA),RNFP0(NSA),RTFP0(NSA)
      ENDDO
      DO NS=1,NSMAX
         READ(21) DELP(NS),PMAX(NS),EMAX(NS)
      ENDDO

      DO NSA=1, NSAMAX
         DO NR=1, NRMAX
            DO NP=1, NPMAX
               READ(21) (FNS(NTH,NP,NR,NSA), NTH=1, NTHMAX) ! -> FNSP, FNSB
            END DO
         END DO
      END DO

! ------

      READ(21) TIMEFP, NT_init

      READ(21) ( (RN_TEMP(NR,NS), NR=1,NRMAX), NS=1,NSMAX)
      READ(21) ( (RT_TEMP(NR,NS), NR=1,NRMAX), NS=1,NSMAX)
      READ(21) ( (RT_BULK(NR,NS), NR=1,NRMAX), NS=1,NSAMAX)
      READ(21) ( (RNS_DELF(NR,NS), NR=1,NRMAX), NS=1,NSMAX)
      READ(21) (E1(NR),NR=1,NRMAX) ! -> EP

      IF(MODELD.ne.0)THEN
         READ(21) ((((WEIGHR_G(NTH,NP,NR,NSA), NTH=1,NTHMAX),NP=1,NPMAX),NR=1, NRMAX+1),NSA=1,NSAMAX)
      END IF

      CLOSE(21)

      WRITE(6,*) '# DATA WAS SUCCESSFULLY LOADED FROM THE FILE.'

  900 RETURN
      END SUBROUTINE FP_LOAD2
!------------------------------------------      
      SUBROUTINE FP_PRE_LOAD

      USE fpprep
      USE libmpi
      USE fpmpi
      IMPLICIT NONE
      integer:: NSW, N, NSA, ierr, i, NR
!      real:: gut1, gut2
      integer,dimension(1:6):: idata
      integer,dimension(6*nsize):: idata2

!      CALL GUTIME(gut1) 
      CALL fp_comm_setup
      CALL fp_allocate
      call fp_allocate_ntg1
      call fp_allocate_ntg2

      CALL mtx_set_communicator(comm_nr)
      allocate(MTXLEN(nsize),MTXPOS(nsize))
      CALL mtx_set_communicator(comm_nsanr)
      allocate(SAVLEN(nsize)) 
      allocate(SAVPOS(nsize,NSAEND-NSASTART+1)) 
      CALL mtx_reset_communicator
      allocate(Rank_Partition_Data(1:6,0:nsize-1))

      CALL mtx_reset_communicator

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

      END SUBROUTINE FP_PRE_LOAD
!------------------------------------------      
      SUBROUTINE FP_POST_LOAD

      USE fpprep
      USE libmpi
      USE fpmpi
      USE plprof
      USE fpnfrr
      USE fpnflg
      USE fpsub
      IMPLICIT NONE
      integer:: NSA, ierr, NR, NS, NP, NTH, NBEAM
      double precision:: FL, RHON
      real(8),dimension(NRMAX,NSMAX):: tempt, tempn 
      TYPE(pl_plf_type),DIMENSION(NSMAX):: PLF

!     ----- set parameters for target species -----
      CALL fp_set_normalize_param
!      CALL GUTIME(gut2) 
!      WRITE(*,'(A,E14.6)') "PRE LOAD GUT= ", gut2-gut1

      CALL mtx_reset_communicator
      CALL Bcast_loaded_data
      CALL scatter_fns_to_fns0

      CALL mtx_set_communicator(comm_nsa)
      CALL update_fnsb
      CALL mtx_reset_communicator

      DO NR=NRSTART,NREND
         EP(NR)=E1(NR)
         EM(NR)=0.D0
      END DO

      CALL FNSP_INIT_EDGE
      IF(MODELS.eq.3) CALL NF_LG_FUNCTION
      IF(MODELS.ne.0) CALL NF_REACTION_COEF
      CALL fp_continue(ierr)
      IF(MODELS.ge.2)THEN
         CALL ALLREDUCE_NF_RATE
         CALL PROF_OF_NF_REACTION_RATE(1)
      END IF
      CALL fp_set_initial_value_from_f

      END SUBROUTINE FP_POST_LOAD
!------------------------------------------      
      SUBROUTINE Bcast_loaded_data

      USE libmpi
      IMPLICIT NONE
      INTEGER,DIMENSION(99):: idata
      real(8),DIMENSION(99):: rdata
      integer:: NR, NSB, NSA, NTH, NP

      idata(1)=NT_init
      CALL mtx_broadcast_integer(idata,1)
      NT_init = idata(1)

      DO NR=1,NRMAX
         rdata(NR)=E1(NR)
      END DO
      rdata(NRMAX+1)=TIMEFP
      CALL mtx_broadcast_real8(rdata,NRMAX+1)
      DO NR=1,NRMAX
         E1(NR)=rdata(NR)
      END DO
      TIMEFP=rdata(NRMAX+1)

      CALL mtx_broadcast_real8(RN_TEMP,NRMAX*NSMAX)
      CALL mtx_broadcast_real8(RT_TEMP,NRMAX*NSMAX)
      CALL mtx_broadcast_real8(RT_BULK,NRMAX*NSAMAX)
      CALL mtx_broadcast_real8(RNS_DELF,NRMAX*NSMAX)

      IF(MODELD.ne.0)THEN
         CALL mtx_broadcast_real8(WEIGHR_G,NTHMAX*NPMAX*(NRMAX+1)*NSAMAX)
         DO NSA=NSASTART, NSAEND
            DO NR=NRSTART, NREND+1
               DO NP=NPSTART, NPEND
                  DO NTH=1, NTHMAX
                     WEIGHR(NTH,NP,NR,NSA)=WEIGHR_G(NTH,NP,NR,NSA)
                  END DO
               END DO
            END DO
         END DO
      END IF

      END SUBROUTINE Bcast_loaded_data
!------------------------------------------      
! GATHER DATA to RANK 0
      SUBROUTINE FP_PRE_SAVE

      USE libmpi
      IMPLICIT NONE
      double precision,dimension(NRMAX)::temp
      double precision,dimension(NRSTART:NREND)::temp_l
      double precision,dimension(NTHMAX,NPSTART:NPEND,NRSTART:NREND):: dsend
      double precision,dimension(nthmax,npmax,nrmax):: drecv
      double precision,dimension(NTHMAX,NPMAX):: dsend_max, drecv_max
      integer,dimension(NTHMAX,NPMAX)::vloc_max

      double precision,dimension(nthmax,npmax,nrmax+1,nsastart:nsaend):: temp_l2
      integer:: nsend
      INTEGER:: NR, NSB, NSA, NTH, NP!, dest, source, tag, nn, Isum
!!!!!!!!!!!!

      IF(MODELD.ne.0)THEN
         DO NSA=NSASTART, NSAEND
            CALL mtx_set_communicator(comm_nrnp)
            nsend=NTHMAX*(NPEND-NPSTART+1)*(NREND-NRSTART+1)
            DO NR=NRSTART, NREND
               DO NP=NPSTART, NPEND
                  DO NTH=1, NTHMAX
                     dsend(nth,np,nr)=WEIGHR(NTH,NP,NR,NSA)
                  END DO
               END DO
            END DO
            CALL mtx_allgather_real8(dsend,nsend,drecv)
            DO NR=1, NRMAX
               DO NP=1, NPMAX
                  DO NTH=1, NTHMAX
                     temp_l2(nth,np,nr,nsa)=drecv(NTH,NP,NR)
                  END DO
               END DO
            END DO
            IF(NREND.eq.NRMAX)THEN
               CALL mtx_set_communicator(comm_np)
               DO NP=NPSTART, NPEND
                  DO NTH=1, NTHMAX
                     dsend_max(nth,np)=WEIGHR(NTH,NP,NRMAX+1,NSA)
                  END DO
               END DO
               CALL mtx_allreduce_real8(dsend_max,NTHMAX*NPMAX,3, &
                                        drecv_max,vloc_max)
               DO NP=1, NPMAX
                  DO NTH=1, NTHMAX
                     temp_l2(nth,np,nrmax+1,nsa)=drecv_max(NTH,NP)
                  END DO
               END DO
            END IF
         END DO
         CALL mtx_set_communicator(comm_nsa)
         nsend=NTHMAX*NPMAX*(NRMAX+1)*(NSAEND-NSASTART+1) 
         CALL mtx_gather_real8(temp_l2,nsend,WEIGHR_G)
      END IF ! END MODELD

      CALL mtx_reset_communicator

      END SUBROUTINE FP_PRE_SAVE
!------------------------------------------      
      SUBROUTINE OPEN_EVOLVE_DATA_OUTPUT

      USE libmpi      
      IMPLICIT NONE

      IF(OUTPUT_TXT_F1.eq.1) OPEN(9,file="f1_1.dat") 
      IF(OUTPUT_TXT_BEAM_WIDTH.eq.1)THEN
         OPEN(24,file="time_evol_f1.dat") 
      END IF
      IF(OUTPUT_TXT_BEAM_DENS.eq.1) OPEN(31,file="t-beam_count.dat")
      IF(OUTPUT_TXT_HEAT_PROF.eq.1)THEN
         OPEN(29,file="RPCS2_NSB_H.dat")
         OPEN(30,file="RPCS2_NSB_H_DEL.dat")
         OPEN(32,file="RSPL_NSB.dat")
      END IF
      IF(MODELS.eq.2)THEN
         open(25,file='fusion_reaction_rate_prof.dat')
         open(26,file='fusion_reaction_rate_tot.dat')
      END IF
      IF(OUTPUT_TXT_DELTA_F.eq.1) OPEN(33,file='output_fns_del.dat')

!      open(19,file='p-T_bulk.dat')
!      open(20,file='DCPP.dat')
!      open(16,file='err_message_for_RT_BULK.dat')

      END SUBROUTINE OPEN_EVOLVE_DATA_OUTPUT
!------------------------------------------      
      SUBROUTINE CLOSE_EVOLVE_DATA_OUTPUT

      USE libmpi      
      IMPLICIT NONE

      IF(OUTPUT_TXT_F1.eq.1) close(9)
      IF(OUTPUT_TXT_BEAM_WIDTH.eq.1)THEN
         close(24)
!         close(34)
      END IF
      IF(OUTPUT_TXT_BEAM_DENS.eq.1) close(31)
      IF(OUTPUT_TXT_HEAT_PROF.eq.1)THEN
         close(29)
         close(30)
         close(32)
      END IF

      IF(MODELS.eq.2)THEN
         close(25)
         close(26)
      END IF
      IF(OUTPUT_TXT_DELTA_F.eq.1) close(33)
!      close(19)
!      close(20)  
!      close(16)  

      END SUBROUTINE CLOSE_EVOLVE_DATA_OUTPUT
!------------------------------------------      

      END MODULE fpfile
