!     ***********************************************************

!           SET INITIAL PROFILE

!     ***********************************************************

      SUBROUTINE TRPROF

      USE TRCOMM, ONLY : ABRHOG, AJ, AJNB, AJNBU, AJOH, AJTOR, AJU, ALP, ANC, ANFE, ANNU, AR1RHOG, ARRHOG, BB, BP, BPRHO, &
     &                   DR, DVRHO, DVRHOG, ETA, GRG, GRM, MDLEQ0, MDLJQ, MDLUF, MDNCLS, MDNI, MODELG, MODEP, NFM, NGR,   &
     &                   NGST, NGT, NRAMAX, NRM, NRMAX, NROMAX, NSM, NT, NTMAX, PBM, PBMU, PEX, PI, PICU, PN, PNBU, PNC,  &
     &                   PNFE, PNS, PNSA, PNSS, PNSSA, PRF, PROFJ1, PROFJ2, PROFN1, PROFN2, PROFT1, PROFT2, PROFU1, PROFU2,&
     &                   PT, PTS, PZ, PZC, PZFE, Q0, QP, QPU, RDP, RG, RHOA, RIP, RIPA, RIPS, RM, RMJRHO, RMJRHOU, RMU0,  &
     &                   RN, RNF, RNFU, RNU, RPSI, RR, RT, RTU, RW, SEX, SNBU, SUMPBM, SWLU, T, TPRE, TST, TTRHO, TTRHOG, &
     &                   VPAR, VPOL, VPRP, VSEC, VTOR, WROT, WROTU, RDPS, &
           & KUFDIR, KUFDCG, KUFDEV, MDLXP, NTMAX_SAVE,ALLOCATE_TRCOMM, PECU
      USE TRCOM1, ONLY : NTAMAX,KDIRX
      IMPLICIT NONE
      INTEGER(4):: IERR, NF, NR, NS
      REAL(8)   :: ANEAVE, ANI, ANZ, DILUTE, FACT, FACTOR0, FACTORM, FACTORP, FCTR, PROF
      REAL(8)   :: SUML
      REAL(8), DIMENSION(NRMAX) :: DSRHO
      
         if(nrmax>NRM)then
            print*,'CAUTION!!'
            print*,'NRMAX=',NRMAX,'<','NRM=',NRM
            print*,'Code is Stopped and you must check parameters.'
            stop
         endif 
        CALL ALLOCATE_TRCOMM(IERR)
        IF(IERR.NE.0) RETURN
        IF(MDLUF.NE.0.AND.MDLXP.NE.0) CALL IPDB_OPEN(KUFDEV, KUFDCG)
        IF(MDLUF.NE.0) CALL UFILE_INTERFACE(KDIRX,KUFDIR,KUFDEV,KUFDCG,0)

         CALL TR_EQS_SELECT(0)

         IF(MDLUF.EQ.1) THEN
!            IF(INIT.EQ.2.AND.NT.NE.0) THEN
            IF(NT.NE.0) THEN
               NT=0
               NTMAX=NTMAX_SAVE
            ENDIF
            CALL TR_UFILE_CONTROL(1)
         ELSEIF(MDLUF.EQ.2) THEN
            CALL TR_UFILE_CONTROL(2)
         ELSEIF(MDLUF.EQ.3) THEN
!            IF(INIT.EQ.2.AND.NT.NE.0) THEN
            IF(NT.NE.0) THEN
               NT=0
               NTMAX=NTMAX_SAVE
            ENDIF
            CALL TR_UFILE_CONTROL(3)
         ELSE
            CALL TR_UFILE_CONTROL(0)
         ENDIF
         

      NT    = 0
      T     = 0.D0
      TPRE  = 0.D0
      TST   = 0.D0
      VSEC  = 0.D0
      NGR   = 0
      NGT   = 0
      NGST  = 0
      RIP   = RIPS

      IF(MDLUF.EQ.1.OR.MDLUF.EQ.3) THEN
         IF(NTMAX.GT.NTAMAX) NTMAX=NTAMAX
      ENDIF
      
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

         select case(MDLUF)
         case(1)
            IF(MDNI.EQ.0) THEN
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFN1)**PROFN2
               RN(NR,1) = RNU(1,NR,1)
               RN(NR,2) = RNU(1,NR,2)
               RN(NR,3) = (PN(3)-PNS(3))*PROF+PNS(3)
               RN(NR,4) = (PN(4)-PNS(4))*PROF+PNS(4)

               IF(RHOA.EQ.1.D0.OR.(RHOA.NE.1.D0.AND.NR.LE.NRAMAX)) THEN
                  PROF   = (1.D0-(ALP(1)*RM(NR)/RHOA)**PROFT1)**PROFT2
                  RT(NR,1:2) = (PT(1:2)-PTS(1:2))*PROF+PTS(1:2)
               ELSEIF(RHOA.NE.1.D0.AND.NR.GT.NRAMAX) THEN
                  RT(NR,1) = RTU(1,NR,1)
                  RT(NR,2) = RTU(1,NR,2)
               ENDIF
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFT1)**PROFT2
               RT(NR,3) = (RTU(1,NR,2)-RTU(1,NRMAX,2))*PROF +RTU(1,NRMAX,2)
               RT(NR,4) = (RTU(1,NR,2)-RTU(1,NRMAX,2))*PROF +RTU(1,NRMAX,2)
            ELSE
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFN1)**PROFN2
               RN(NR,1) = RNU(1,NR,1)
               RN(NR,2) = RNU(1,NR,2)
               RN(NR,3) = RNU(1,NR,3)
               RN(NR,4) = (PN(4)-PNS(4))*PROF+PNS(4)

!               IF(RHOA.EQ.1.D0.OR.(RHOA.NE.1.D0.AND.NR.LE.NRAMAX)) THEN
!                  PROF   = (1.D0-(ALP(1)*RM(NR)/RHOA)**PROFT1)**PROFT2
!                  DO NS=1,3
!                     RT(NR,NS) = (PT(NS)-PTS(NS))*PROF+PTS(NS)
!                  ENDDO
!               ELSEIF(RHOA.NE.1.D0.AND.NR.GT.NRAMAX) THEN
!                  RT(NR,1) = RTU(1,NR,1)
!                  RT(NR,2) = RTU(1,NR,2)
!                  RT(NR,3) = RTU(1,NR,3)
!               ENDIF
               RT(NR,1) = RTU(1,NR,1)
               RT(NR,2) = RTU(1,NR,2)
               RT(NR,3) = RTU(1,NR,3)
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFT1)**PROFT2
               RT(NR,4) = (RTU(1,NR,2)-RTU(1,NRMAX,2))*PROF +RTU(1,NRMAX,2)
            ENDIF

            PEX(NR,1) = PNBU(1,NR,1)
            PEX(NR,2) = PNBU(1,NR,2)
            PEX(NR,3) = 0.D0
            PEX(NR,4) = 0.D0
            SEX(NR,1:NSM) = 0.D0
            PRF(NR,1) = PICU(1,NR,1)+PECU(1,NR)
            PRF(NR,2) = PICU(1,NR,2)
            PRF(NR,3) = 0.D0
            PRF(NR,4) = 0.D0

            RNF(NR,1) = RNFU(1,NR)
            RNF(NR,2:NFM) = RNFU(1,NR)
            PBM(NR)   = PBMU(1,NR)
            WROT(NR)  = WROTU(1,NR)
            VTOR(NR)  = WROTU(1,NR)*RMJRHOU(1,NR)
         case(2)
            IF(MDNI.EQ.0) THEN ! ** MDNI **
            IF(MODEP.EQ.1) THEN
               IF(RHOA.EQ.1.D0.OR.(RHOA.NE.1.D0.AND.NR.LE.NRAMAX)) THEN
                  PROF   = (1.D0-(ALP(1)*RM(NR)/RHOA)**PROFN1)**PROFN2
                  RN(NR,1:2) = (PN(1:2)-PNS(1:2))*PROF+PNS(1:2)
               ELSEIF(RHOA.NE.1.D0.AND.NR.GT.NRAMAX) THEN
                  RN(NR,1) = RNU(1,NR,1)
                  RN(NR,2) = RNU(1,NR,2)
               ENDIF
               PROF   = (1.D0-(ALP(1)*RM(NR)/RHOA)**PROFN1)**PROFN2
               RN(NR,3) = (RNU(1,NR,2)-RNU(1,NRMAX,2))*PROF +RNU(1,NRMAX,2)
               RN(NR,4) = (RNU(1,NR,2)-RNU(1,NRMAX,2))*PROF +RNU(1,NRMAX,2)

               IF(RHOA.EQ.1.D0.OR.(RHOA.NE.1.D0.AND.NR.LE.NRAMAX)) THEN
                  PROF   = (1.D0-(ALP(1)*RM(NR)/RHOA)**PROFT1)**PROFT2
                  RT(NR,1:2) = (PT(1:2)-PTS(1:2))*PROF+PTS(1:2)
               ELSEIF(RHOA.NE.1.D0.AND.NR.GT.NRAMAX) THEN
                  RT(NR,1) = RTU(1,NR,1)
                  RT(NR,2) = RTU(1,NR,2)
               ENDIF
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFT1)**PROFT2
               RT(NR,3) = (RTU(1,NR,2)-RTU(1,NRMAX,2))*PROF +RTU(1,NRMAX,2)
               RT(NR,4) = (RTU(1,NR,2)-RTU(1,NRMAX,2))*PROF +RTU(1,NRMAX,2)
            ELSEIF(MODEP.EQ.2) THEN
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFN1)**PROFN2
               RN(NR,1) = RNU(1,NR,1)
               RN(NR,2) = RNU(1,NR,2)
               RN(NR,3) = (PN(3)-PNS(3))*PROF+PNS(3)
               RN(NR,4) = (PN(4)-PNS(4))*PROF+PNS(4)

               RT(NR,3) = RTU(1,NR,2)
               RT(NR,4) = RTU(1,NR,2)
            ELSEIF(MODEP.EQ.3) THEN
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFN1)**PROFN2
               RN(NR,1) = RNU(1,NR,1)
               RN(NR,2) = RNU(1,NR,2)
               RN(NR,3) = (PN(3)-PNS(3))*PROF+PNS(3)
               RN(NR,4) = (PN(4)-PNS(4))*PROF+PNS(4)

               IF(RHOA.EQ.1.D0.OR.(RHOA.NE.1.D0.AND.NR.LE.NRAMAX)) THEN
                  PROF   = (1.D0-(ALP(1)*RM(NR)/RHOA)**PROFT1)**PROFT2
                  RT(NR,1:2) = (PT(1:2)-PTS(1:2))*PROF+PTS(1:2)
               ELSEIF(RHOA.NE.1.D0.AND.NR.GT.NRAMAX) THEN
                  RT(NR,1) = RTU(1,NR,1)
                  RT(NR,2) = RTU(1,NR,2)
               ENDIF
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFT1)**PROFT2
               RT(NR,3) = (RTU(1,NR,2)-RTU(1,NRMAX,2))*PROF +RTU(1,NRMAX,2)
               RT(NR,4) = (RTU(1,NR,2)-RTU(1,NRMAX,2))*PROF +RTU(1,NRMAX,2)
            ENDIF
            ELSE ! ** MDNI **
            IF(MODEP.EQ.1) THEN
               IF(RHOA.EQ.1.D0.OR.(RHOA.NE.1.D0.AND.NR.LE.NRAMAX)) THEN
                  PROF   = (1.D0-(ALP(1)*RM(NR)/RHOA)**PROFN1)**PROFN2
                  RN(NR,1:3) = (PN(1:3)-PNS(1:3))*PROF+PNS(1:3)
               ELSEIF(RHOA.NE.1.D0.AND.NR.GT.NRAMAX) THEN
                  RN(NR,1:3) = RNU(1,NR,1:3)
               ENDIF
               PROF   = (1.D0-(ALP(1)*RM(NR)/RHOA)**PROFN1)**PROFN2
               RN(NR,4) = (RNU(1,NR,2)-RNU(1,NRMAX,2))*PROF +RNU(1,NRMAX,2)

               IF(RHOA.EQ.1.D0.OR.(RHOA.NE.1.D0.AND.NR.LE.NRAMAX)) THEN
                  PROF   = (1.D0-(ALP(1)*RM(NR)/RHOA)**PROFT1)**PROFT2
                  RT(NR,1:3) = (PT(1:3)-PTS(1:3))*PROF+PTS(1:3)
               ELSEIF(RHOA.NE.1.D0.AND.NR.GT.NRAMAX) THEN
                  RT(NR,1:3) = RTU(1,NR,1:3)
               ENDIF
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFT1)**PROFT2
               RT(NR,4) = (RTU(1,NR,2)-RTU(1,NRMAX,2))*PROF +RTU(1,NRMAX,2)
            ELSEIF(MODEP.EQ.2) THEN
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFN1)**PROFN2
               RN(NR,1:3) = RNU(1,NR,1:3)
               RN(NR,4) = (PN(4)-PNS(4))*PROF+PNS(4)

               RT(NR,4) = RTU(1,NR,2)
            ELSEIF(MODEP.EQ.3) THEN
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFN1)**PROFN2
               RN(NR,1:3) = RNU(1,NR,1:3)
               RN(NR,4) = (PN(4)-PNS(4))*PROF+PNS(4)

               IF(RHOA.EQ.1.D0.OR.(RHOA.NE.1.D0.AND.NR.LE.NRAMAX)) THEN
                  PROF   = (1.D0-(ALP(1)*RM(NR)/RHOA)**PROFT1)**PROFT2
                  RT(NR,1:3) = (PT(1:3)-PTS(1:3))*PROF+PTS(1:3)
               ELSEIF(RHOA.NE.1.D0.AND.NR.GT.NRAMAX) THEN
                  RT(NR,1:3) = RTU(1,NR,1:3)
               ENDIF
               PROF   = (1.D0-(ALP(1)*RM(NR))**PROFT1)**PROFT2
               RT(NR,4) = (RTU(1,NR,2)-RTU(1,NRMAX,2))*PROF +RTU(1,NRMAX,2)
            ENDIF
            ENDIF ! ** MDNI **

            PEX(NR,1)=PNBU(1,NR,1)
            PEX(NR,2)=PNBU(1,NR,2)
            PEX(NR,3)=0.D0
            PEX(NR,4)=0.D0

            SEX(NR,1)=SNBU(1,NR,1)+SWLU(1,NR)/PZ(2)
            SEX(NR,2)=SNBU(1,NR,2)+SWLU(1,NR)
            SEX(NR,3)=0.D0
            SEX(NR,4)=0.D0
!            PRF(NR,1:NSM)=0.D0
            PRF(NR,1) = PICU(1,NR,1)+PECU(1,NR)
!            print*,"NR:PRF:PECU",NR,PRF(NR,1),PECU(1,NR)
!            PRF(NR,1)=PRF(NR,1)*0.0d0
!            print*,"PRF????????",NR,PRF(NR,1),PECU(1,NR)
            PRF(NR,2) = PICU(1,NR,2)
            RNF(NR,1)=RNFU(1,NR)
            RNF(NR,2:NFM)=0.D0
            PBM(NR)=0.D0
            WROT(NR) =WROTU(1,NR)
            VTOR(NR) =WROTU(1,NR)*RMJRHOU(1,NR)
         case(3)
            PROF   = (1.D0-(ALP(1)*RM(NR))**PROFN1)**PROFN2
            RN(NR,1) = RNU(1,NR,1)
            RN(NR,2) = RNU(1,NR,2)
            RN(NR,3) = RNU(1,NR,3)
            RN(NR,4) = (PN(4)-PNS(4))*PROF+PNS(4)
            RT(NR,1) = RTU(1,NR,1)
            RT(NR,2) = RTU(1,NR,2)
            RT(NR,3) = RTU(1,NR,3)
            PROF   = (1.D0-(ALP(1)*RM(NR))**PROFT1)**PROFT2
            RT(NR,4) = (RTU(1,NR,2)-RTU(1,NRMAX,2))*PROF +RTU(1,NRMAX,2)

!$$$            PROF   = (1.D0-(ALP(1)*RM(NR))**PROFN1)**PROFN2
!$$$            DO NS=1,NSM
!$$$               RN(NR,NS) = (PN(NS)-PNS(NS))*PROF+PNS(NS)
!$$$            ENDDO
!$$$C
!$$$            PROF   = (1.D0-(ALP(1)*RM(NR))**PROFT1)**PROFT2
!$$$            DO NS=1,NSM
!$$$               RT(NR,NS) = (PT(NS)-PTS(NS))*PROF+PTS(NS)
!$$$            ENDDO

            PEX(NR,1) = PNBU(1,NR,1)
            PEX(NR,2) = PNBU(1,NR,2)
            PEX(NR,3) = 0.D0
            PEX(NR,4) = 0.D0
            SEX(NR,1:NSM) = 0.D0
            PRF(NR,1) = PICU(1,NR,1)
            PRF(NR,2) = PICU(1,NR,2)
            PRF(NR,3) = 0.D0
            PRF(NR,4) = 0.D0
            RNF(NR,1:NFM) = 0.D0
            PBM(NR)=0.D0
            PRF(NR,4) = 0.D0
            WROT(NR)  = WROTU(1,NR)
            VTOR(NR)  = WROTU(1,NR)*RMJRHOU(1,NR)
         case default
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
         end select

!!-------------------------------------------------
!! 2009/07/09
!         if (modelg==20 .or. modelg==22) then
!            call t3d_input_prof('nt',1)
!         elseif(modelg==21 .or. modelg==23) then
!            call t3d_input_prof('nt',2)
!         end if
!!--------------------------------------------------

         IF(MDLEQ0.EQ.1) THEN
            PROF   = (1.D0-(ALP(1)*RM(NR))**PROFU1)**PROFU2
            RN(NR,7) = (PN(7)-PNS(7))*PROF+PNS(7)
            RN(NR,8) = (PN(8)-PNS(8))*PROF+PNS(8)
            ANNU(NR) = RN(NR,7)+RN(NR,8)
         ENDIF

         RW(NR,1:NFM) = 0.D0

         SUMPBM=SUMPBM+PBM(NR)
      ENDDO
!      CALL PLDATA_SETR(RG,RM)
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

      IF(MDLUF.NE.3) THEN
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
      ELSE
         PNSS(1:NSM)=PNS(1:NSM)
         PNSS(7)=PNS(7)
         PNSS(8)=PNS(8)
         IF(RHOA.NE.1.D0) THEN
            PNSSA(1:NSM)=PNSA(1:NSM)
            PNSSA(7)=PNSA(7)
            PNSSA(8)=PNSA(8)
         ENDIF
      ENDIF

!     *** CALCULATE PROFILE OF AJ(R) ***

!     *** THIS MODEL ASSUMES GIVEN JZ PROFILE ***


      IF(MDLUF.EQ.1) THEN
         IF(MDLJQ.EQ.0) THEN ! *** MDLJQ ***
         NR=1
            AJ(NR)=AJU(1,NR)
            AJNB(NR)=AJNBU(1,NR)
            AJOH(NR)=AJ(NR)-AJNB(NR)
            FACTOR0=RMU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
            RDP(NR)=FACTOR0*DR/FACTORP
            BP(NR) =AR1RHOG(NR)*RDP(NR)/RR
         DO NR=2,NRMAX
            AJ(NR)=AJU(1,NR)
            AJNB(NR)=AJNBU(1,NR)
            AJOH(NR)=AJ(NR)-AJNB(NR)
            FACTOR0=RMU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)/TTRHOG(NR-1)
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
            RDP(NR)=(FACTOR0*DR+FACTORM*RDP(NR-1))/FACTORP
            BP(NR) =AR1RHOG(NR)*RDP(NR)/RR
         ENDDO
         NR=1
            FACTOR0=RR/(RMU0*DVRHO(NR))
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
            AJTOR(NR)=FACTOR0*FACTORP*RDP(NR)/DR
         DO NR=2,NRMAX
            FACTOR0=RR/(RMU0*DVRHO(NR))
            FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
            AJTOR(NR)=FACTOR0*(FACTORP*RDP(NR)-FACTORM*RDP(NR-1))/DR
         ENDDO
         QP(1:NRMAX)=TTRHOG(1:NRMAX)*ARRHOG(1:NRMAX)*DVRHOG(1:NRMAX) &
     &              /(4.D0*PI**2*RDP(1:NRMAX))
         ELSE IF(MDLJQ.EQ.1) THEN ! *** MDLJQ ***
            QP(1:NRMAX) =QPU(1,1:NRMAX)
            RDP(1:NRMAX)=TTRHOG(1:NRMAX)*ARRHOG(1:NRMAX)*DVRHOG(1:NRMAX)/(4.D0*PI**2*QP(1:NRMAX))
            BP(1:NRMAX) =AR1RHOG(1:NRMAX)*RDP(1:NRMAX)/RR

         NR=1
            FACTOR0=TTRHO(NR)**2/(RMU0*BB*DVRHO(NR))
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
            AJ(NR) =FACTOR0*FACTORP*RDP(NR)/DR
            AJOH(NR)=AJ(NR)
         DO NR=2,NRMAX
            FACTOR0=TTRHO(NR)**2/(RMU0*BB*DVRHO(NR))
            FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)/TTRHOG(NR-1)
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
            AJ(NR) =FACTOR0*(FACTORP*RDP(NR)-FACTORM*RDP(NR-1))/DR
            AJOH(NR)=AJ(NR)
         ENDDO
         NR=1
            FACTOR0=RR/(RMU0*DVRHO(NR))
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
            AJTOR(NR) =FACTOR0*FACTORP*RDP(NR)/DR
         DO NR=2,NRMAX
            FACTOR0=RR/(RMU0*DVRHO(NR))
            FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
            AJTOR(NR) =FACTOR0*(FACTORP*RDP(NR)-FACTORM*RDP(NR-1))/DR
         ENDDO
         ELSE IF(MDLJQ.EQ.2) THEN ! *** MDLJQ ***
            QP(1:NRMAX) =QPU(1,1:NRMAX)
            RDP(1:NRMAX)=TTRHOG(1:NRMAX)*ARRHOG(1:NRMAX)*DVRHOG(1:NRMAX)/(4.D0*PI**2*QP(1:NRMAX))
            BP(1:NRMAX) =AR1RHOG(1:NRMAX)*RDP(1:NRMAX)/RR

         NR=1
            AJ(NR) =1.d-6
            AJOH(NR)=AJ(NR)
         DO NR=2,NRMAX
            AJ(NR) =1.D-6
            AJOH(NR)=AJ(NR)
         ENDDO
         NR=1
            AJTOR(NR) =0.d0
         DO NR=2,NRMAX
            AJTOR(NR) =0.d0
         ENDDO
         ENDIF ! *** MDLJQ ***
      ELSEIF(MDLUF.EQ.2) THEN
         IF(MDLJQ.EQ.0) THEN  ! *** MDLJQ ***
            NR=1
            AJ(NR)=AJU(1,NR)
            AJNB(NR)=AJNBU(1,NR)
            AJOH(NR)=AJ(NR)-AJNB(NR)
            FACTOR0=RMU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
            RDP(NR)=FACTOR0*DR/FACTORP
            BP(NR) =AR1RHOG(NR)*RDP(NR)/RR
         DO NR=2,NRMAX
            AJ(NR)=AJU(1,NR)
            AJNB(NR)=AJNBU(1,NR)
            AJOH(NR)=AJ(NR)-AJNB(NR)
            FACTOR0=RMU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)/TTRHOG(NR-1)
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
            RDP(NR)=(FACTOR0*DR+FACTORM*RDP(NR-1))/FACTORP
            BP(NR) =AR1RHOG(NR)*RDP(NR)/RR
         ENDDO

!$$$         DO NR=1,NRMAX
!$$$            AJ(NR)=AJU(1,NR)
!$$$            AJNB(NR)=AJNBU(1,NR)
!$$$            AJOH(NR)=AJ(NR)-AJNB(NR)
!$$$            RDP(NR)=RMU0/(RR*DVRHOG(NR)*ABRHOG(NR))*SUM(AJ(1:NRMAX)*DVRHO(1:NRMAX)*DR
!$$$            BP(NR) =AR1RHOG(NR)*RDP(NR)/RR
!$$$         ENDDO


         NR=1
            FACTOR0=RR/(RMU0*DVRHO(NR))
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
            AJTOR(NR)=FACTOR0*FACTORP*RDP(NR)/DR
         DO NR=2,NRMAX
            FACTOR0=RR/(RMU0*DVRHO(NR))
            FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
            AJTOR(NR)=FACTOR0*(FACTORP*RDP(NR)-FACTORM*RDP(NR-1))/DR
         ENDDO
         QP(1:NRMAX)=TTRHOG(1:NRMAX)*ARRHOG(1:NRMAX)*DVRHOG(1:NRMAX)/(4.D0*PI**2*RDP(1:NRMAX))

         ELSEIF(MDLJQ.EQ.1) THEN ! *** MDLJQ ***

            QP(1:NRMAX) =QPU(1,1:NRMAX)
            RDP(1:NRMAX)=TTRHOG(1:NRMAX)*ARRHOG(1:NRMAX)*DVRHOG(1:NRMAX)/(4.D0*PI**2*QP(1:NRMAX))
            BP(1:NRMAX) =AR1RHOG(1:NRMAX)*RDP(1:NRMAX)/RR

            NR=1
               FACTOR0=TTRHO(NR)**2/(RMU0*BB*DVRHO(NR))
               FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
               AJ(NR) =FACTOR0*FACTORP*RDP(NR)/DR
               AJOH(NR)=AJ(NR)
            DO NR=2,NRMAX
               FACTOR0=TTRHO(NR)**2/(RMU0*BB*DVRHO(NR))
               FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)/TTRHOG(NR-1)
               FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
               AJ(NR) =FACTOR0*(FACTORP*RDP(NR)-FACTORM*RDP(NR-1))/DR
               AJOH(NR)=AJ(NR)
            ENDDO
            NR=1
               FACTOR0=RR/(RMU0*DVRHO(NR))
               FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
               AJTOR(NR) =FACTOR0*FACTORP*RDP(NR)/DR
            DO NR=2,NRMAX
               FACTOR0=RR/(RMU0*DVRHO(NR))
               FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)
               FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
               AJTOR(NR) =FACTOR0*(FACTORP*RDP(NR)-FACTORM*RDP(NR-1))/DR
            ENDDO
         ELSEIF(MDLJQ.EQ.2) THEN ! *** MDLJQ ***

            QP(1:NRMAX) =QPU(1,1:NRMAX)
            RDP(1:NRMAX)=TTRHOG(1:NRMAX)*ARRHOG(1:NRMAX)*DVRHOG(1:NRMAX)/(4.D0*PI**2*QP(1:NRMAX))
            BP(1:NRMAX) =AR1RHOG(1:NRMAX)*RDP(1:NRMAX)/RR

            NR=1
               AJ(NR) =1.D-6
               AJOH(NR)=AJ(NR)
            DO NR=2,NRMAX
               AJ(NR) =1.D-6
               AJOH(NR)=AJ(NR)
            ENDDO
            NR=1
               AJTOR(NR) =0.D0
            DO NR=2,NRMAX
               AJTOR(NR) =0.D0
            ENDDO
         ENDIF ! *** MDLJQ ***
         RIPA=DVRHOG(NRAMAX)*ABRHOG(NRAMAX)*RDP(NRAMAX)*1.D-6 /(2.D0*PI*RMU0)
      ELSEIF(MDLUF.EQ.3) THEN
         DO NR=1,NRMAX
            AJOH(NR)=AJU(1,NR)
            AJ(NR)  =AJU(1,NR)
         ENDDO

         NR=1
            FACTOR0=RMU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
            RDP(NR)=FACTOR0*DR/FACTORP
         DO NR=2,NRMAX
            FACTOR0=RMU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)/TTRHOG(NR-1)
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
            RDP(NR)=(FACTORM*RDP(NR-1)+FACTOR0*DR)/FACTORP
         ENDDO
         BP(1:NRMAX)=AR1RHOG(1:NRMAX)*RDP(1:NRMAX)/RR

         RDPS=2.D0*PI*RMU0*RIP*1.D6/(DVRHOG(NRMAX)*ABRHOG(NRMAX))
         FACT=RDPS/RDP(NRMAX)
         AJOH(1:NRMAX)=FACT*AJOH(1:NRMAX)
         AJ(1:NRMAX)  =AJOH(1:NRMAX)
         BP(1:NRMAX)  =FACT*BP(1:NRMAX)
         QP(1:NRMAX)  =TTRHOG(1:NRMAX)*ARRHOG(1:NRMAX)*DVRHOG(1:NRMAX) /(4.D0*PI**2*RDP(1:NRMAX))
      ELSE
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
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
            RDP(NR)=FACTOR0*DR/FACTORP
            BP(NR)=AR1RHOG(NR)*RDP(NR)/RR
         DO NR=2,NRMAX
            FACTOR0=RMU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)/TTRHOG(NR-1)
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
            RDP(NR)=(FACTORM*RDP(NR-1)+FACTOR0*DR)/FACTORP
            BP(NR)=AR1RHOG(NR)*RDP(NR)/RR
         ENDDO
         NR=1
            FACTOR0=RR/(RMU0*DVRHO(NR))
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
            AJTOR(NR) =FACTOR0*FACTORP*RDP(NR)/DR
         DO NR=2,NRMAX
            FACTOR0=RR/(RMU0*DVRHO(NR))
            FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
            AJTOR(NR) =FACTOR0*(FACTORP*RDP(NR)-FACTORM*RDP(NR-1))/DR
         ENDDO

         RDPS=2.D0*PI*RMU0*RIP*1.D6/(DVRHOG(NRMAX)*ABRHOG(NRMAX))
         FACT=RDPS/RDP(NRMAX)
         RDP(1:NRMAX)=FACT*RDP(1:NRMAX)
         AJOH(1:NRMAX)=FACT*AJOH(1:NRMAX)
         AJ(1:NRMAX)  =AJOH(1:NRMAX)
         BP(1:NRMAX)  =FACT*BP(1:NRMAX)
         QP(1:NRMAX)  =TTRHOG(1:NRMAX)*ARRHOG(1:NRMAX)*DVRHOG(1:NRMAX)/(4.D0*PI**2*RDP(1:NRMAX))
      ENDIF

!!-------------------------------------------------
!! 2009/07/09
!      if (modelg==20 .or. modelg==22) then
!         call t3d_input_prof('iota',1)
!      elseif(modelg==21 .or. modelg==23) then
!         call t3d_input_prof('iota',2)
!      end if
!!--------------------------------------------------      

!      write(6,*) 'in trprof'
!      write(6,'(1P5E12.4)') (qp(nr),nr=1,nrmax)
!      pause
!     *** calculate q_axis ***
      Q0=FCTR(RG(1),RG(2),QP(1),QP(2))

!     calculate plasma current inside the calucated region (rho <= rhoa)
!     necessary for MDLEQB = 1 and MDLUF /= 0
      IF(MDLUF.EQ.1.OR.MDLUF.EQ.3) THEN
         DSRHO(1:NRAMAX)=DVRHO(1:NRAMAX)/(2.D0*PI*RMJRHO(1:NRAMAX))
         RIPA=SUM(AJ(1:NRAMAX)*DSRHO(1:NRAMAX))*DR/1.D6
      ENDIF

!     *** THIS MODEL ASSUMES CONSTANT EZ ***

      IF(PROFJ1.LE.0.D0.OR.MDNCLS.EQ.1) THEN
         CALL TRZEFF
         CALL TRCFET
         IF(PROFJ1.GT.0.D0.AND.MDNCLS.EQ.1) GOTO 2000

         AJOH(1:NRMAX)=1.D0/ETA(1:NRMAX)
         AJ(1:NRMAX)  =1.D0/ETA(1:NRMAX)

         NR=1
            FACTOR0=RMU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
            RDP(NR)=FACTOR0*DR/FACTORP
         DO NR=2,NRMAX
            FACTOR0=RMU0*BB*DVRHO(NR)*AJ(NR)/TTRHO(NR)**2
            FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)/TTRHOG(NR-1)
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )/TTRHOG(NR  )
            RDP(NR)=(FACTORM*RDP(NR-1)+FACTOR0*DR)/FACTORP
         ENDDO
         NR=1
            FACTOR0=RR/(RMU0*DVRHO(NR))
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
            AJTOR(NR)=FACTOR0*FACTORP*RDP(NR)/DR
         DO NR=2,NRMAX
            FACTOR0=RR/(RMU0*DVRHO(NR))
            FACTORM=DVRHOG(NR-1)*ABRHOG(NR-1)
            FACTORP=DVRHOG(NR  )*ABRHOG(NR  )
            AJTOR(NR)=FACTOR0*(FACTORP*RDP(NR)-FACTORM*RDP(NR-1))/DR
         ENDDO
         BP(1:NRMAX)=AR1RHOG(1:NRMAX)*RDP(1:NRMAX)/RR

         RDPS=2.D0*PI*RMU0*RIP*1.D6/(DVRHOG(NRMAX)*ABRHOG(NRMAX))
         FACT=RDPS/RDP(NRMAX)
         RDP(1:NRMAX)  =FACT*RDP(1:NRMAX)
         AJOH(1:NRMAX) =FACT*AJOH(1:NRMAX)
         AJ(1:NRMAX)   =AJOH(1:NRMAX)
         AJTOR(1:NRMAX)=FACT*AJTOR(1:NRMAX)
         BP(1:NRMAX)   =FACT*BP(1:NRMAX)
         QP(1:NRMAX)   =TTRHOG(1:NRMAX)*ARRHOG(1:NRMAX)*DVRHOG(1:NRMAX)/(4.D0*PI**2*RDP(1:NRMAX))
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

!      IF(RHOA.NE.1.D0) NRMAX=NROMAX

      call trsetg(ierr)

      RETURN
      END SUBROUTINE TRPROF

!     ***********************************************************

!           SET GEOMETRICAL FACTOR

!     ***********************************************************

      subroutine trsetg(ierr)

      use trcomm, only : &
           & modelg, nrmax, knameq
      use tr_bpsd, only: tr_bpsd_get
      use equnit_mod, only: eq_parm,eq_prof,eq_calc,eq_load
      use equunit_mod, only: equ_prof,equ_calc
      use pl_vmec_mod, only: pl_vmec
      integer, intent(out):: ierr
      character(len=80):: line

         if(modelg.eq.3) then
            write(line,'(A,I5)') 'nrmax=',nrmax+1
            call eq_parm(2,line,ierr)
            write(line,'(A,I5)') 'nthmax=',64
            call eq_parm(2,line,ierr)
            write(line,'(A,I5)') 'nsumax=',0
            call eq_parm(2,line,ierr)
            call eq_load(modelg,knameq,ierr) ! load eq data and calculate eq
            call tr_bpsd_get(ierr)  ! 
            if(ierr.ne.0) write(6,*) 'XX2 ierr=',ierr
         elseif(modelg.eq.7) then
            call pl_vmec(knameq,ierr) ! load vmec data
            call tr_bpsd_get(ierr)  ! 
            call trgout
         elseif(modelg.eq.8) then
            call equ_prof ! initial calculation of eq
            call equ_calc         ! recalculate eq
            call tr_bpsd_get(ierr)  ! 
            call trgout
         elseif(modelg.eq.9) then
            call eq_prof ! initial calculation of eq
            call eq_calc         ! recalculate eq
            call tr_bpsd_get(ierr)  ! 
            call trgout
!!---------------------
!! 2009/07/10
         elseif(modelg.eq.22 .or. modelg.eq.23 .or. modelg.eq.24) then
            call t3d_vmec2tr(knameq)
!! --------------------
         end if
      return
      end subroutine trsetg

!     ***********************************************************

!           SET GEOMETRIC FACTOR AT HALF MESH

!     ***********************************************************

      SUBROUTINE TRSTGF

      USE TRCOMM, ONLY : ABRHO, ABRHOU, AR1RHO, AR1RHOU, AR2RHO, AR2RHOU, ARRHO, ARRHOU, BB, BP, BPRHO, DVRHO, DVRHOU, &
     &                   EPSRHO, MDLUF, MDPHIA, MODELG, NRMAX, PHIA, PI, QP, QRHO, RA, RG, RHOG, RHOM, RJCB, RKAP, RKPRHO, &
     &                   RKPRHOU, RM, RMJRHO, RMJRHOU, RMNRHO, RMNRHOU, RR, TTRHO, TTRHOU, VOLAU
      IMPLICIT NONE
      INTEGER(4) :: NR
      REAL(8)    :: RKAPS, RHO_A

      RKAPS=SQRT(RKAP)
      IF(MDLUF.NE.0) THEN
         DO NR=1,NRMAX
            TTRHO(NR)=TTRHOU(1,NR)
            DVRHO(NR)=DVRHOU(1,NR)
            ABRHO(NR)=ABRHOU(1,NR)
            ARRHO(NR)=ARRHOU(1,NR)
            AR1RHO(NR)=AR1RHOU(1,NR)
            AR2RHO(NR)=AR2RHOU(1,NR)
            RMJRHO(NR)=RMJRHOU(1,NR)
            RMNRHO(NR)=RMNRHOU(1,NR)
            RKPRHO(NR)=RKPRHOU(1,NR)
            IF(MDPHIA.EQ.0) THEN
!     define rho_a from phi_a data
               RHO_A=SQRT(PHIA/(PI*BB))
               RJCB(NR)=1.D0/RHO_A
               RHOM(NR)=RM(NR)*RHO_A
               RHOG(NR)=RG(NR)*RHO_A
            ELSE
               RHO_A=SQRT(VOLAU(1)/(2.D0*PI**2*RMJRHOU(1,NRMAX)))
               RJCB(NR)=1.D0/RHO_A
               RHOM(NR)=RM(NR)/RJCB(NR)
               RHOG(NR)=RG(NR)/RJCB(NR)
            ENDIF
            EPSRHO(NR)=RMNRHO(NR)/RMJRHO(NR)
         ENDDO
         CALL FLUX
      ELSE
         DO NR=1,NRMAX
            BPRHO(NR)=BP(NR)
            QRHO(NR)=QP(NR)

            TTRHO(NR)=BB*RR
            DVRHO(NR)=2.D0*PI*RKAP*RA*RA*2.D0*PI*RR*RM(NR)
            ABRHO(NR)=1.D0/(RKAPS*RA*RR)**2
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
         ENDDO
      ENDIF

      RETURN
      END SUBROUTINE TRSTGF

!     ***********************************************************

!           GEOMETRIC QUANTITIES AT GRID MESH

!     ***********************************************************

      SUBROUTINE TRGFRG

      USE TRCOMM, ONLY : ABB2RHOG, ABRHO, ABRHOG, AIB2RHOG, AR1RHO, AR1RHOG, AR2RHO, AR2RHOG, ARHBRHOG, ARRHO, ARRHOG, &
     &                   BB, DVRHO, DVRHOG, EPSRHO, NRMAX, RG, RKPRHO, RKPRHOG, RM, TTRHO, TTRHOG
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

         ABB2RHOG(NR)=BB**2*(1.D0+0.5D0*EPSRHO(NR)**2)
         AIB2RHOG(NR)=(1.D0+1.5D0*EPSRHO(NR)**2)/BB**2
         ARHBRHOG(NR)=AR2RHOG(NR)*AIB2RHOG(NR)
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

      USE TRCOMM, ONLY : MDLUF, NRAMAX, NSM, NSTM, PNSS, PNSSA, PTS, PTSA, RHOA, RN, RT, PNSSO,PTSO,PNSSAO,PTSAO
      IMPLICIT NONE
      INTEGER(4), INTENT(IN) :: NSW
      INTEGER(4) :: NS

      IF(RHOA.EQ.1.D0) RETURN

      IF(MDLUF.EQ.0) THEN
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
      ELSE
         IF(NSW.EQ.0) THEN
            DO NS=1,NSM
               PNSSO(NS)=PNSS (NS)
               PTSO (NS)=PTS  (NS)

               PNSS (NS)=PNSSA(NS)
               PTS  (NS)=PTSA (NS)
            ENDDO
         ELSE
            DO NS=1,NSM
               PNSS (NS)=PNSSO(NS)
               PTS  (NS)=PTSO (NS)
            ENDDO
         ENDIF
      ENDIF
      RETURN
!
      ENTRY TR_EDGE_DETERMINER(NSW)

      IF(MDLUF.EQ.0.AND.RHOA.NE.1.D0) THEN
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
      END SUBROUTINE TR_EDGE_SELECTOR
