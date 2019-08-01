!     ***********************************************************

!           SET INITIAL PROFILE

!     ***********************************************************

      SUBROUTINE TRPROF(ierr)

      USE TRCOMM, ONLY : &
           ABRHOG, AJ, AJNB, AJNBU, AJOH, AJTOR, AJU, ALP, &
           ANC, ANFE, ANNU, AR1RHOG, ARRHOG, BB, BP, BPRHO, DR, DVRHO, &
           DVRHOG, ETA, GRG, GRM, MDLEQ0, MDLJQ, MDNCLS, MDNI, &
           MODELG, MODEP, NFM, NGR, NGST, NGT, NRAMAX, NRMAX, NROMAX, NSM, &
           NT, NTMAX, PBM, PBMU, PEX, PI, PICU, PN, PNBU, PNC, PNFE, PNS, &
           PNSA, PNSS, PNSSA, PRF, PROFJ1, PROFJ2, PROFN1, PROFN2, PROFT1, &
           PROFT2, PROFU1, PROFU2, PT, PTS, PZ, PZC, PZFE, Q0, QP, QPU, &
           RDP, RG, RHOA, RIP, RIPA, RIPS, RM, RMJRHO, RMJRHOU, RMU0, RN, &
           RNF, RNFU, RNU, RPSI, RR, RT, RTF, RTU, RU, RW, SEX, SNBU, &
           SUMPBM, SWLU, T, TPRE, TST, TTRHO, TTRHOG, VPAR, VPOL, &
           VPRP, VSEC, VTOR, WROT, WROTU, RDPS, KUFDIR, KUFDCG, KUFDEV, &
           NTMAX_SAVE,ALLOCATE_TRCOMM, ABVRHOG, RDPVRHOG, &
           ABVRHO,abrho
      USE TRCOM0, ONLY : NSTM
      USE TRCOM1, ONLY : NTAMAX,KDIRX
      IMPLICIT NONE
      INTEGER(4):: IERR, NR, NS, NF
      REAL(8)   :: ANEAVE, ANI, ANZ, DILUTE
      REAL(8)   :: FACT, FACTOR0, FACTORM, FACTORP, FCTR, PROF
      REAL(8)   :: SUML
      REAL(8), DIMENSION(NRMAX) :: DSRHO

      CALL ALLOCATE_TRCOMM(IERR)
      IF(IERR.NE.0) RETURN

      CALL TR_EQS_SELECT(0)

      NRAMAX=INT(RHOA*NRMAX)
      DR = 1.D0/DBLE(NRMAX)

      NT    = 0
      T     = 0.D0
      TPRE  = 0.D0
      TST   = 0.D0
      VSEC  = 0.D0
      NGR   = 0
      NGT   = 0
      NGST  = 0
      RIP   = RIPS

      CALL TR_EDGE_DETERMINER(0)
      CALL TR_EDGE_SELECTOR(0)
      IF(RHOA.NE.1.D0) NRMAX=NROMAX
      DO NR=1,NRMAX
         RG(NR) = DBLE(NR)*DR
         RM(NR) =(DBLE(NR)-0.5D0)*DR
         VTOR(NR)=0.D0
         VPAR(NR)=0.D0
         VPRP(NR)=0.D0
         VPOL(NR)=0.D0
         WROT(NR)=0.D0
         DO NS=1,NSTM
            RN(NR,NS)=0.D0
            RT(NR,NS)=0.D0
            RU(NR,NS)=0.D0
         END DO
         DO NF=1,NFM
            RW(NR,NF)=0.D0
            RNF(NR,NF)=0.D0
            RTF(NR,NF)=0.D0
         END DO

         PROF   = (1.D0-(ALP(1)*RM(NR))**PROFN1)**PROFN2
         RN(NR,1:NSM) = (PN(1:NSM)-PNS(1:NSM))*PROF+PNS(1:NSM)

         PROF   = (1.D0-(ALP(1)*RM(NR))**PROFT1)**PROFT2
         RT(NR,1:NSM) = (PT(1:NSM)-PTS(1:NSM))*PROF+PTS(1:NSM)

         PEX(NR,1:NSM) = 0.D0
         SEX(NR,1:NSM) = 0.D0
         PRF(NR,1:NSM) = 0.D0
         RNF(NR,1:NFM) = 0.D0
         PBM(NR)=0.D0
         WROT(NR)=0.D0
         VTOR(NR)=0.D0

         IF(MDLEQ0.EQ.1) THEN
            PROF   = (1.D0-(ALP(1)*RM(NR))**PROFU1)**PROFU2
            RN(NR,7) = (PN(7)-PNS(7))*PROF+PNS(7)
            RN(NR,8) = (PN(8)-PNS(8))*PROF+PNS(8)
            ANNU(NR) = RN(NR,7)+RN(NR,8)
         ENDIF

         RW(NR,1:NFM) = 0.D0

         SUMPBM=SUMPBM+PBM(NR)
      ENDDO
      CALL TR_EDGE_DETERMINER(1)
      CALL TR_EDGE_SELECTOR(1)

!     *** CALCULATE GEOMETRIC FACTOR ***

      CALL TRSTGF
      CALL TRGFRG

!     *** CALCULATE PZC,PZFE ***

      CALL TRZEFF

!     *** CALCULATE ANEAVE ***

      ANEAVE=SUM(RN(1:NRMAX,1)*RM(1:NRMAX))*2.D0*DR

!     *** CALCULATE IMPURITY DENSITY
!                ACCORDING TO ITER PHYSICS DESIGN GUIDELINE ***

      DO NR=1,NRMAX
         ANC (NR)= (0.9D0+0.60D0*(0.7D0/ANEAVE)**2.6D0)*PNC *1.D-2*RN(NR,1)
         ANFE(NR)= (0.0D0+0.05D0*(0.7D0/ANEAVE)**2.3D0)*PNFE*1.D-2*RN(NR,1)
         ANI = SUM(PZ(2:NSM)*RN(NR,2:NSM))
         ANZ = PZFE(NR)*ANFE(NR)+PZC(NR)*ANC(NR)
         DILUTE = 1.D0-ANZ/ANI
         RN(NR,2:NSM) = RN(NR,2:NSM)*DILUTE
      ENDDO
      PNSS(1)=PNS(1)
      PNSS(2:NSM)=PNS(2:NSM)*DILUTE
      PNSS(7)=PNS(7)
      PNSS(8)=PNS(8)
      IF(RHOA.NE.1.D0) THEN
         PNSSA(1)=PNSA(1)
         PNSSA(2:NSM)=PNSA(2:NSM)*DILUTE
         PNSSA(7)=PNSA(7)
         PNSSA(8)=PNSA(8)
      ENDIF

!     *** CALCULATE PROFILE OF AJ(R) ***

!     *** THIS MODEL ASSUMES GIVEN JZ PROFILE ***

      DO NR=1,NRMAX
         IF((1.D0-RM(NR)**ABS(PROFJ1)).LE.0.D0) THEN
            PROF=0.D0
         ELSE
            PROF= (1.D0-RM(NR)**ABS(PROFJ1))**ABS(PROFJ2)
         ENDIF
         AJOH(NR)= PROF
         AJ(NR)  = PROF
      ENDDO

      NR=1
         FACTOR0=RMU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
         FACTORP=ABVRHOG(NR  )/TTRHOG(NR  )
         RDPVRHOG(NR)=FACTOR0*DR/FACTORP
         RDP(NR)=RDPVRHOG(NR)*DVRHOG(NR)
         BP(NR)=AR1RHOG(NR)*RDP(NR)/RR
      DO NR=2,NRMAX
         FACTOR0=RMU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
         FACTORM=ABVRHOG(NR-1)/TTRHOG(NR-1)
         FACTORP=ABVRHOG(NR  )/TTRHOG(NR  )
         RDPVRHOG(NR)=(FACTORM*RDPVRHOG(NR-1)+FACTOR0*DR)/FACTORP
         RDP(NR)=RDPVRHOG(NR)*DVRHOG(NR)
         BP(NR)=AR1RHOG(NR)*RDP(NR)/RR
      ENDDO
      NR=1
         FACTOR0=RR/(RMU0*DVRHO(NR))
         FACTORP=ABVRHOG(NR  )
         AJTOR(NR) =FACTOR0*FACTORP*RDPVRHOG(NR)/DR
      DO NR=2,NRMAX
         FACTOR0=RR/(RMU0*DVRHO(NR))
         FACTORM=ABVRHOG(NR-1)
         FACTORP=ABVRHOG(NR  )
         AJTOR(NR) =FACTOR0*(FACTORP*RDPVRHOG(NR)-FACTORM*RDPVRHOG(NR-1))/DR
      ENDDO

      RDPS=2.D0*PI*RMU0*RIP*1.D6*DVRHOG(NRMAX)/ABVRHOG(NRMAX)
      FACT=RDPS/RDP(NRMAX)
      RDP(1:NRMAX)=FACT*RDP(1:NRMAX)
      RDPVRHOG(1:NRMAX)=FACT*RDPVRHOG(1:NRMAX)
      AJOH(1:NRMAX)=FACT*AJOH(1:NRMAX)
      AJ(1:NRMAX)  =AJOH(1:NRMAX)
      BP(1:NRMAX)  =FACT*BP(1:NRMAX)
      QP(1:NRMAX)  =TTRHOG(1:NRMAX)*ARRHOG(1:NRMAX) &
                    /(4.D0*PI**2*RDPVRHOG(1:NRMAX))

!     *** calculate q_axis ***
      Q0=FCTR(RG(1),RG(2),QP(1),QP(2))

!     *** THIS MODEL ASSUMES CONSTANT EZ ***

      IF(PROFJ1.LE.0.D0.OR.MDNCLS.EQ.1) THEN
         CALL TRZEFF
         CALL TRCFET
         IF(PROFJ1.GT.0.D0.AND.MDNCLS.EQ.1) GOTO 2000

         AJOH(1:NRMAX)=1.D0/ETA(1:NRMAX)
         AJ(1:NRMAX)  =1.D0/ETA(1:NRMAX)

         NR=1
            FACTOR0=RMU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORP=ABVRHOG(NR  )/TTRHOG(NR  )
            RDPVRHOG(NR)=FACTOR0*DR/FACTORP
            RDP(NR)=RDPVRHOG(NR)*DVRHOG(NR)
         DO NR=2,NRMAX
            FACTOR0=RMU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORM=ABVRHOG(NR-1)/TTRHOG(NR-1)
            FACTORP=ABVRHOG(NR  )/TTRHOG(NR  )
            RDPVRHOG(NR)=(FACTORM*RDPVRHOG(NR-1)+FACTOR0*DR)/FACTORP
            RDP(NR)=RDPVRHOG(NR)*DVRHOG(NR)
         ENDDO
         NR=1
            FACTOR0=RR/(RMU0*DVRHO(NR))
            FACTORP=ABVRHOG(NR  )
            AJTOR(NR)=FACTOR0*FACTORP*RDPVRHOG(NR)/DR
         DO NR=2,NRMAX
            FACTOR0=RR/(RMU0*DVRHO(NR))
            FACTORM=ABVRHOG(NR-1)
            FACTORP=ABVRHOG(NR  )
            AJTOR(NR)=FACTOR0*(FACTORP*RDPVRHOG(NR)-FACTORM*RDPVRHOG(NR-1))/DR
         ENDDO
         BP(1:NRMAX)=AR1RHOG(1:NRMAX)*RDP(1:NRMAX)/RR

         RDPS=2.D0*PI*RMU0*RIP*1.D6*DVRHOG(NRMAX)/ABVRHOG(NRMAX)
         FACT=RDPS/RDP(NRMAX)
         RDP(1:NRMAX)  =FACT*RDP(1:NRMAX)
         RDPVRHOG(1:NRMAX)=FACT*RDPVRHOG(1:NRMAX)
         AJOH(1:NRMAX) =FACT*AJOH(1:NRMAX)
         AJ(1:NRMAX)   =AJOH(1:NRMAX)
         AJTOR(1:NRMAX)=FACT*AJTOR(1:NRMAX)
         BP(1:NRMAX)   =FACT*BP(1:NRMAX)
         QP(1:NRMAX)   =TTRHOG(1:NRMAX)*ARRHOG(1:NRMAX) &
                        /(4.D0*PI**2*RDPVRHOG(1:NRMAX))
      ENDIF
 2000 CONTINUE
      SUML=0.D0
      DO NR=1,NRMAX
         SUML=SUML+RDP(NR)*DR
         RPSI(NR)=SUML
         BPRHO(NR)=BP(NR)
      ENDDO

      GRG(1)=0.0
      GRM(1:NRMAX)  =SNGL(RM(1:NRMAX))
      GRG(2:NRMAX+1)=SNGL(RG(1:NRMAX))

      call trsetg(ierr)

      RETURN
      END SUBROUTINE TRPROF

!     ***********************************************************

!           SET GEOMETRICAL FACTOR

!     ***********************************************************

      SUBROUTINE trsetg(ierr)

      USE trcomm, ONLY : modelg, nrmax, knameq, knameq2, RR, RA
      USE tr_bpsd, ONLY: tr_bpsd_init,tr_bpsd_set,tr_bpsd_get
      USE equnit_mod, ONLY: eq_parm,eq_prof,eq_calc,eq_load
      USE pl_vmec_mod, ONLY: pl_vmec
      IMPLICIT NONE
      INTEGER,INTENT(OUT):: ierr
      CHARACTER(len=80):: line

      CALL tr_bpsd_init
      CALL tr_bpsd_set(ierr)

      IF(modelg.EQ.3.OR.modelg.EQ.5.OR.modelg.EQ.8) THEN
         WRITE(line,'(A,I5)') 'nrmax=',nrmax+1
         CALL eq_parm(2,line,ierr)
         WRITE(line,'(A,I5)') 'nthmax=',64
         CALL eq_parm(2,line,ierr)
         WRITE(line,'(A,I5)') 'nsumax=',0
         CALL eq_parm(2,line,ierr)
         IF(modelg.EQ.8) THEN
            WRITE(line,'(A,A,A,A)') 'knameq2=','"',TRIM(knameq2),'"'
            WRITE(6,'(A,A)') 'line=',line
            CALL eq_parm(2,line,ierr)
         END IF
         CALL eq_load(modelg,knameq,ierr) ! load eq data and calculate eq
         IF(ierr.NE.0) THEN
            WRITE(6,*) 'XX eq_load: ierr=',ierr
            RETURN
         ENDIF
         CALL tr_bpsd_get(ierr)  ! 
         IF(ierr.NE.0) WRITE(6,*) 'XX tr_bpsd_get: ierr=',ierr
!         call trgout
      ELSEIF(modelg.EQ.7) THEN
         CALL pl_vmec(knameq,ierr) ! load vmec data
         CALL tr_bpsd_get(ierr)  ! 
!         call trgout
      ELSEIF(modelg.EQ.9) THEN
         CALL eq_prof ! initial calculation of eq
         CALL eq_calc         ! recalculate eq
         CALL tr_bpsd_get(ierr)  ! 
!         call trgout
      ENDIF

      RETURN
      END SUBROUTINE trsetg
      
!     ***********************************************************

!           SET GEOMETRIC FACTOR AT HALF MESH

!     ***********************************************************

      SUBROUTINE TRSTGF

      USE TRCOMM, ONLY : &
           ABRHO, ABRHOU, AR1RHO, AR1RHOU, AR2RHO, AR2RHOU, ARRHO, &
           ARRHOU, BB, BP, BPRHO, DVRHO, DVRHOU, EPSRHO, &
           MDPHIA, MODELG, NRMAX, PHIA, PI, QP, QRHO, RA, RG, &
           RHOG, RHOM, RJCB, RKAP, RKPRHO, RKPRHOU, RM, RMJRHO, &
           RMJRHOU, RMNRHO, RMNRHOU, RR, TTRHO, TTRHOU, VOLAU, &
           ABVRHO, ABVRHOG, PVOLRHOG, PSURRHOG, ABB1RHO
      IMPLICIT NONE
      INTEGER(4) :: NR
      REAL(8)    :: RKAPS, RHO_A

      RKAPS=SQRT(RKAP)

      DO NR=1,NRMAX
         BPRHO(NR)=BP(NR)
         QRHO(NR)=QP(NR)

         TTRHO(NR)=BB*RR
         DVRHO(NR)=2.D0*PI*RKAP*RA*RA*2.D0*PI*RR*RM(NR)
         ABRHO(NR)=1.D0/(RKAPS*RA*RR)**2
         ABVRHO(NR)=DVRHO(NR)**2*ABRHO(NR)
         ARRHO(NR)=1.D0/RR**2
         AR1RHO(NR)=1.D0/(RKAPS*RA)
         AR2RHO(NR)=1.D0/(RKAPS*RA)**2
         RMJRHO(NR)=RR
         RMNRHO(NR)=RA*RG(NR)
         RKPRHO(NR)=RKAP
         RJCB(NR)=1.D0/(RKAPS*RA)
         RHOM(NR)=RM(NR)/RJCB(NR)
         RHOG(NR)=RG(NR)/RJCB(NR)
         EPSRHO(NR)=RMNRHO(NR)/RMJRHO(NR)
         ABB1RHO(NR)=BB*(1.D0+0.25D0*EPSRHO(NR)**2)
         PVOLRHOG(NR)=PI*RKAP*(RA*RG(NR))**2*2.D0*PI*RR
         PSURRHOG(NR)=PI*(RKAP+1.D0)*RA*RG(NR)*2.D0*PI*RR
      ENDDO

      RETURN
      END SUBROUTINE TRSTGF


!     ***********************************************************

!           CALCULATING FLUX FROM SOURCE TERM FILES

!     ***********************************************************

      SUBROUTINE FLUX

      USE TRCOMM, ONLY : DR, DVRHO, MDLFLX, NRMAX, NSM, PZ, RGFLX, SNBU, SWLU
      IMPLICIT NONE
      INTEGER(4):: NR
      REAL(8),DIMENSION(NRMAX)::  SALEL, SALIL


      IF(MDLFLX.EQ.0) THEN
         RGFLX(1:NRMAX,1:NSM)=0.D0
      ELSE
         DO NR=1,NRMAX
            SALEL(NR)=SNBU(1,NR,1)+SWLU(1,NR)/PZ(2)
            RGFLX(NR,1)=SUM(SALEL(1:NR)*DVRHO(1:NR))*DR
            SALIL(NR)=SNBU(1,NR,2)+SWLU(1,NR)
            RGFLX(NR,2)=SUM(SALIL(1:NR)*DVRHO(1:NR))*DR
            RGFLX(NR,3:NSM)=0.D0
         ENDDO
      ENDIF

      RETURN
      END SUBROUTINE FLUX

!     ***********************************************************

!           GEOMETRIC QUANTITIES AT GRID MESH

!     ***********************************************************

      SUBROUTINE TRGFRG

      USE TRCOMM, ONLY : &
           ABB2RHOG, ABRHO, ABRHOG, AIB2RHOG, AR1RHO, AR1RHOG, AR2RHO, &
           AR2RHOG, ARHBRHOG, ARRHO, ARRHOG, BB, DVRHO, DVRHOG, EPSRHO, &
           NRMAX, RG, RKPRHO, RKPRHOG, RM, TTRHO, TTRHOG, ABVRHO, ABVRHOG
      IMPLICIT NONE
      INTEGER(4) :: NR
      REAL(8)    :: RGL

      DO NR=1,NRMAX-1
         AR1RHOG(NR)=0.5D0*(AR1RHO(NR)+AR1RHO(NR+1))
         AR2RHOG(NR)=0.5D0*(AR2RHO(NR)+AR2RHO(NR+1))
         RKPRHOG(NR)=0.5D0*(RKPRHO(NR)+RKPRHO(NR+1))
         TTRHOG (NR)=0.5D0*(TTRHO (NR)+TTRHO (NR+1))
         DVRHOG (NR)=0.5D0*(DVRHO (NR)+DVRHO (NR+1))
         ARRHOG (NR)=0.5D0*(ARRHO (NR)+ARRHO (NR+1))
         ABRHOG (NR)=0.5D0*(ABRHO (NR)+ABRHO (NR+1))

         ABVRHOG(NR)=DVRHOG(NR)**2*ABRHOG(NR)
         ABB2RHOG(NR)=BB**2*(1.D0+0.5D0*EPSRHO(NR)**2)
         AIB2RHOG(NR)=(1.D0+1.5D0*EPSRHO(NR)**2)/BB**2
      ENDDO
      NR=NRMAX
         RGL=RG(NR)

         CALL AITKEN(RGL,AR1RHOG(NR),RM,AR1RHO,2,NRMAX)
         CALL AITKEN(RGL,AR2RHOG(NR),RM,AR2RHO,2,NRMAX)
         CALL AITKEN(RGL,RKPRHOG(NR),RM,RKPRHO,2,NRMAX)
         CALL AITKEN(RGL,TTRHOG (NR),RM,TTRHO ,2,NRMAX)
         CALL AITKEN(RGL,DVRHOG (NR),RM,DVRHO ,2,NRMAX)
         CALL AITKEN(RGL,ARRHOG (NR),RM,ARRHO ,2,NRMAX)
         CALL AITKEN(RGL,ABRHOG (NR),RM,ABRHO ,2,NRMAX)

         ABVRHOG(NR)=DVRHOG(NR)**2*ABRHOG(NR)
         ABB2RHOG(NR)=BB**2*(1.D0+0.5D0*EPSRHO(NR)**2)
         AIB2RHOG(NR)=(1.D0+1.5D0*EPSRHO(NR)**2)/BB**2
         ARHBRHOG(NR)=AR2RHOG(NR)*AIB2RHOG(NR)

      RETURN
      END SUBROUTINE TRGFRG

!     ***********************************************************

!           EDGE VALUE SELECTOR

!     ***********************************************************

      SUBROUTINE TR_EDGE_SELECTOR(NSW)

!        NSW = 0: store edge value; substitute rhoa value
!              1: restore original edge value

      USE TRCOMM, ONLY : &
           NRAMAX, NSM, NSTM, PNSS, PNSSA, PTS, PTSA, RHOA, RN, RT, &
           PNSSO,PTSO,PNSSAO,PTSAO
      IMPLICIT NONE
      INTEGER(4), INTENT(IN) :: NSW
      INTEGER(4) :: NS

      IF(RHOA.EQ.1.D0) RETURN

      IF(NSW.EQ.0) THEN
         DO NS=1,NSM
            PNSSO(NS)=PNSS(NS)
            PTSO (NS)=PTS (NS)

            PNSS (NS)=PNSSAO(NS)
            PTS  (NS)=PTSAO (NS)
         ENDDO
      ELSE
         DO NS=1,NSM
            PNSS (NS)=PNSSO(NS)
            PTS  (NS)=PTSO (NS)
         ENDDO
      ENDIF
      RETURN
    END SUBROUTINE TR_EDGE_SELECTOR
!
    SUBROUTINE TR_EDGE_DETERMINER(NSW)

!        NSW = 0: store edge value; substitute rhoa value
!              1: restore original edge value

      USE TRCOMM, ONLY : &
           NRAMAX, NSM, NSTM, PNSS, PNSSA, PTS, PTSA, RHOA, RN, RT, &
           PNSSO,PTSO,PNSSAO,PTSAO
      IMPLICIT NONE
      INTEGER(4), INTENT(IN) :: NSW
      INTEGER(4) :: NS

      IF(RHOA.NE.1.D0) THEN
         IF(NSW.EQ.0) THEN
            DO NS=1,NSM
               PNSSAO(NS)=PNSS(NS)
               PTSAO (NS)=PTS (NS)
            ENDDO
         ELSE
            DO NS=1,NSM
               PNSSAO(NS)=RN(NRAMAX,NS)
               PTSAO (NS)=RT(NRAMAX,NS)
               PNSSA (NS)=PNSSAO(NS)
               PTSA  (NS)=PTSAO (NS)
            ENDDO
         ENDIF
      ENDIF

      RETURN
    END SUBROUTINE TR_EDGE_DETERMINER
