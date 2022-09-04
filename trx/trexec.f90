! trexec.f90

MODULE trexec

  PRIVATE
  PUBLIC tr_exec,tr_eval,tr_coef_decide

CONTAINS

!     ***********************************************************

!           MAIN ROUTINE FOR TRANSPORT CALCULATION

!     ***********************************************************

      SUBROUTINE tr_exec(DT,IERR)

      USE TRCOMM, ONLY : &
     &   AJ, AJU, AMM, ANC, ANFE, ANNU, AX, AY, AZ, BB, &
     &   DR, DVRHO, DVRHOG, EPSLTR, LDAB, LMAXTR, MDLEQB, MDLEQE, &
     &   MDLEQN, MDLPCK, MLM, NEQMAX, NFM, &
     &   NRAMAX, NRMAX, NROMAX, NSM, NSMAX, NSS, NST, NSV, &
     &   NTUM, &
     &   PA, PZ, PZC, PZFE, RDP, RHOA, RIP, RIPU, &
     &   RMU0, RN, RT, RU, RW, T, TTRHO, TTRHOG, &
     &   VLOOP, VSEC, X, XV, Y, YV, Z, ZV ,NEQMAXM, DIPDT, &
         RDPVRHOG, abvrhog, rkind
      USE TRCOM1, ONLY : TMU, TMU1, NTXMAX, NTXMAX1
      USE libbnd
      USE libitp
      IMPLICIT NONE
      REAL(rkind),INTENT(IN) :: DT
      INTEGER,INTENT(OUT) :: IERR
      INTEGER:: I, ICHCK, INFO, J, L, LDB, M, MWRMAX, &
           N, NEQ, NEQ1, NEQRMAX, NR, NRHS, NSSN, NSSN1, &
           NSTN, NSTN1, NSVN, NSVN1, KL, KU
      REAL(rkind)   :: AJL, FACTOR0, FACTORM, FACTORP, TSL
      INTEGER,DIMENSION(NEQMAXM*NRMAX) :: IPIV
      REAL(rkind),DIMENSION(NEQMAXM*NRMAX)    :: XX
      REAL(rkind),DIMENSION(2,NRMAX)  :: YY
      REAL(rkind),DIMENSION(NSMAX,NRMAX):: ZZ

      IERR=0
      ICHCK=0

      L=0
!     /* Setting New Variables */
      CALL TRATOX

!     /* Store Variables for Convergence Check */
      forall(J=1:NEQMAX,NR=1:NRMAX) XX(NEQMAX*(NR-1)+J) = XV(J,NR)
      YY(1:NFM,1:NRMAX) = YV(1:NFM,1:NRMAX)

 2000 CONTINUE

!      CALL TR_EDGE_SELECTOR(0)

!     /* Calcualte matrix coefficients */

      CALL TRMTRX(NEQRMAX)

!     /* Solve matrix equation */

      IF(MDLPCK.EQ.0) THEN
         MWRMAX=4*NEQRMAX-1
         CALL BANDRD(AX,X,NEQRMAX*NRMAX,MWRMAX,LDAB,IERR)
         IF(IERR.EQ.30000) THEN
            WRITE(6,*) 'XX ERROR IN TRLOOP : MATRIX AA IS SINGULAR  AT T = ',T
            IERR=9000
            RETURN
         ENDIF
      ELSEIF(MDLPCK.EQ.1) THEN
         M=NEQRMAX*NRMAX
         N=NEQRMAX*NRMAX
         KL=2*NEQRMAX-1
         KU=2*NEQRMAX-1
         NRHS=1
         LDB=MLM
         CALL LAPACK_DGBTRF(M,N,KL,KU,AX,LDAB,IPIV,INFO)
         IF(INFO.NE.0) THEN
            WRITE(6,*) 'XX ERROR IN TRLOOP : DGBTRF, INFO = ',INFO
            IERR=9001
            RETURN
         ENDIF
         CALL LAPACK_DGBTRS('N',N,KL,KU,NRHS,AX,LDAB,IPIV,X,LDB,INFO)
         IF(INFO.NE.0) THEN
            WRITE(6,*) 'XX ERROR IN TRLOOP : DGBTRS, INFO = ',INFO
         ENDIF
      ELSE
         N=NEQRMAX*NRMAX
         KL=2*NEQRMAX-1
         KU=2*NEQRMAX-1
         NRHS=1
         LDB=MLM
         CALL LAPACK_DGBSV(N,KL,KU,NRHS,AX,LDAB,IPIV,X,LDB,INFO)
         IF(INFO.NE.0) THEN
            WRITE(6,*) 'XX ERROR IN TRLOOP : DGBSV, INFO = ',INFO
         ENDIF
      ENDIF

!    /* Solve equation for fast particl

      Y(1:NFM,1:NRMAX) = Y(1:NFM,1:NRMAX)/AY(1:NFM,1:NRMAX)

!     /* Convergence check */

      DO I=1,NEQRMAX*NRMAX
         IF (ABS(X(I)-XX(I)).GT.EPSLTR*ABS(X(I))) GOTO 3000
      ENDDO
      DO J=1,NFM
      DO NR=1,NRMAX
         IF (ABS(Y(J,NR)-YY(J,NR)).GT.EPSLTR*ABS(Y(J,NR))) GOTO 3000
      ENDDO
      ENDDO

      GOTO 4000

 3000 L=L+1
      IF(L.GE.LMAXTR) GOTO 4000

!     /* Stored Variables for Convergence Check */
      DO I=1,NEQRMAX*NRMAX
         XX(I) = X(I)
      ENDDO
      DO J=1,NFM
      DO NR=1,NRMAX
         YY(J,NR) = Y(J,NR)
      ENDDO
      ENDDO

      DO NR=1,NRMAX
         DO NEQ=1,NEQMAX
            NSSN=NSS(NEQ)
            NSVN=NSV(NEQ)
            NSTN=NST(NEQ)
            IF(NSVN.EQ.0) THEN
               IF(NSTN.EQ.0) THEN
                  RDP(NR) = XV(NEQ,NR)
               ELSE
                  RDP(NR) = 0.5D0*(XV(NEQ,NR)+X(NEQRMAX*(NR-1)+NSTN))
!                  if(nr==nrmax) write(6,*) &
!                         NEQMAX,NEQRMAX,NEQRMAX*(NR-1)+NSTN, &
!                         NEQMAX*(NR-1)+NSTN,NSTN
!                  if(nr==nrmax) write(6,*) &
!                         "NSTN/=0",XV(NEQ,NR),X(NEQRMAX*(NR-1)+NSTN), &
!                         X(NEQMAX*(NR-1)+NSTN)
               ENDIF
               RDPVRHOG(NR) = RDP(NR) / DVRHOG(NR)
            ELSEIF(NSVN.EQ.1) THEN
               IF(MDLEQN.NE.0) THEN
                  IF(NSSN.EQ.1.AND.MDLEQE.EQ.0) THEN
                     RN(NR,NSSN) = 0.D0
                     DO NEQ1=1,NEQMAX
                        NSSN1=NSS(NEQ1)
                        NSVN1=NSV(NEQ1)
                        NSTN1=NST(NEQ1)
                        IF(NSVN1.EQ.1.AND.NSSN1.NE.1) THEN
                           IF(NSTN1.EQ.0) THEN
                              RN(NR,NSSN) = RN(NR,NSSN) &
                                   + PZ(NSSN1)*XV(NEQ1,NR)
                           ELSE
                              RN(NR,NSSN) = RN(NR,NSSN) &
                                   + PZ(NSSN1)*0.5D0*(XV(NEQ1,NR) &
                                                   +X(NEQRMAX*(NR-1)+NSTN1))
                           ENDIF
                        ENDIF
                     ENDDO
                     RN(NR,NSSN) = RN(NR,NSSN) &
                          +PZFE(NR)*ANFE(NR)+PZC(NR)*ANC(NR)
                  ELSEIF(NSSN.EQ.1.AND.MDLEQE.EQ.1) THEN
                     RN(NR,NSSN) = 0.5D0*(XV(NEQ,NR)+X(NEQRMAX*(NR-1)+NEQ))
                  ELSEIF(NSSN.EQ.2.AND.MDLEQE.EQ.2) THEN
                     RN(NR,NSSN) = 0.D0
                     DO NEQ1=1,NEQMAX
                        NSSN1=NSS(NEQ1)
                        NSVN1=NSV(NEQ1)
                        NSTN1=NST(NEQ1)
                        IF(NSVN1.EQ.1.AND.NSSN1.NE.2) THEN
                           IF(NSTN1.EQ.0) THEN
                              RN(NR,NSSN) = RN(NR,NSSN) &
                                   +ABS(PZ(NSSN1))*XV(NEQ1,NR)
                           ELSE
                              RN(NR,NSSN) = RN(NR,NSSN) &
                                   +ABS(PZ(NSSN1))*0.5D0*(XV(NEQ1,NR) &
                                   +X(NEQRMAX*(NR-1)+NSTN1))
                           ENDIF
                        ENDIF
                     ENDDO
                     RN(NR,NSSN) = RN(NR,NSSN) &
                          +PZFE(NR)*ANFE(NR)+PZC(NR)*ANC(NR)
                  ELSE
                     IF(NSTN.EQ.0) THEN
                        RN(NR,NSSN) = XV(NEQ,NR)
                     ELSE
                        RN(NR,NSSN) = 0.5D0*(XV(NEQ,NR) &
                                            + X(NEQRMAX*(NR-1)+NSTN))
                     ENDIF
                  ENDIF
               ELSE
                  IF(NSTN.EQ.0) THEN
                     RN(NR,NSSN) = XV(NEQ,NR)
                  ELSE
                     RN(NR,NSSN) = 0.5D0*(XV(NEQ,NR) +  X(NEQRMAX*(NR-1)+NSTN))
                  ENDIF
               ENDIF
            ELSEIF(NSVN.EQ.2) THEN
               IF(NSTN.EQ.0) THEN
                  IF(NSSN.NE.NSM) THEN
                     RT(NR,NSSN) = XV(NEQ,NR)/RN(NR,NSSN)
!                     write(6,'(A,2I5,1P3E12.4)') '-1- ',NR,NSSN, &
!                     & XV(NEQ,NR),RN(NR,NSSN),RT(NR,NSSN)
                  ELSE
                     IF(RN(NR,NSM).LT.1.D-70) THEN
                        RT(NR,NSM) = 0.D0
                     ELSE
                        RT(NR,NSM) = XV(NEQ,NR)/RN(NR,NSM)
                     ENDIF
                  ENDIF
               ELSE
                  IF(NSSN.NE.NSM) THEN
                     RT(NR,NSSN) = 0.5D0*(XV(NEQ,NR) &
                                         +X(NEQRMAX*(NR-1)+NSTN))/RN(NR,NSSN)
!                     write(6,'(A,3I5,1P4E12.4)') &
!                             '-2- ',NR,NSSN,NSTN, &
!                              XV(NEQ,NR),RN(NR,NSSN),RT(NR,NSSN),&
!                              X(NEQRMAX*(NR-1)+NSTN)
                  ELSE
                     IF(RN(NR,NSM).LT.1.D-70) THEN
                        RT(NR,NSM) = 0.D0
                     ELSE
                        RT(NR,NSM) = 0.5D0*(XV(NEQ,NR) &
                             +X(NEQRMAX*(NR-1)+NSTN))/RN(NR,NSM)
!                     write(6,'(A,2I5,1P3E12.4)') '-3- ',NR,NSM, &
!                     & XV(NEQ,NR),RN(NR,NSM),RT(NR,NSM)
                     ENDIF
                  ENDIF
               ENDIF
            ELSEIF(NSVN.EQ.3) THEN
               IF(NSTN.EQ.0) THEN
                  RU(NR,NSSN) = XV(NEQ,NR)/(PA(NSSN)*AMM*RN(NR,NSSN))
               ELSE
                  RU(NR,NSSN) = 0.5D0*(XV(NEQ,NR) &
                       + X(NEQRMAX*(NR-1)+NSTN))/(PA(NSSN)*AMM*RN(NR,NSSN))
               ENDIF
            ENDIF
         END DO
         IF(NSMAX.GE.8) THEN
            ANNU(NR)=RN(NR,7)+RN(NR,8)
         END IF
!!!         BP(NR)=AR1RHOG(NR)*RDP(NR)/RR
         RW(NR,1)  = 0.5D0*(YV(1,NR)+Y(1,NR))
         RW(NR,2)  = 0.5D0*(YV(2,NR)+Y(2,NR))
      ENDDO

      CALL TRCHCK(ICHCK)
      IF(ICHCK.EQ.1) THEN
         CALL TRGLOB
         CALL TRATOT
         CALL TRATOTN
         CALL TRATOG
         IERR=1
         RETURN
      END IF

!      CALL TR_EDGE_SELECTOR(1)

!     /* update matrix coefficients */

      CALL TRCALC(IERR)
      IF(IERR.NE.0) RETURN
      GOTO 2000

 4000 T=T+DT
      VSEC=VSEC+VLOOP*DT
      RIP=RIP+DIPDT*DT

!      write(6,'(A,1P4E12.4)') "T,RIP,RIPE,DIP=",T,RIP,RIPE,DIPDT*DT

!     /* Making new XV(NEQ,NR) and YV(NF,NR) */
      DO NEQ=1,NEQMAX
         NSTN=NST(NEQ)
         IF(NSTN.NE.0) THEN
            DO NR=1,NRMAX
               XV(NEQ,NR) = X(NEQRMAX*(NR-1)+NSTN)
            ENDDO
         ENDIF
      ENDDO

      YV(1:NFM,1:NRMAX) = Y(1:NFM,1:NRMAX)

!     /* Making New Physical Variables */

      CALL TRXTOA

!      CALL TR_EDGE_SELECTOR(1)

!     /* Check negative temperature */

      IF(ICHCK.EQ.0) CALL TRCHCK(ICHCK)

      IF(ICHCK.EQ.1) THEN
         CALL TRGLOB
         CALL TRATOT
         CALL TRATOTN
         CALL TRATOG
         IERR=1
         RETURN
      ENDIF

!     /* calculate physical quauntities */

      CALL TRCALC(IERR)
      IF(IERR.NE.0) RETURN

! !     *** DATA ACQUISITION FOR SHOWING GRAPH AND STATUS ***

!       IDGLOB=0
!       IF(MOD(NT,NTSTEP).EQ.0) THEN
!          IF(IDGLOB.EQ.0) CALL TRGLOB
!          IDGLOB=1
!          CALL TRSNAP
!       ENDIF
!       IF(MOD(NT,NGTSTP).EQ.0) THEN
!          IF(IDGLOB.EQ.0) CALL TRGLOB
!          IDGLOB=1
!          CALL TRATOT
!          CALL TRATOTN
!       ENDIF
!       IF(MOD(NT,NGRSTP).EQ.0) THEN
!          IF(IDGLOB.EQ.0) CALL TRGLOB
!          IDGLOB=1
!          CALL TRATOG
!       ENDIF

      RETURN
    END SUBROUTINE tr_exec

      SUBROUTINE tr_eval(NT,IERR)

      USE TRCOMM, ONLY : NGRSTP, NGTSTP, NTSTEP
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: NT
      INTEGER,INTENT(OUT) :: IERR
      integer:: IDGLOB

      CALL TRCALC(IERR)
      IF(IERR.ne.0) RETURN

      IDGLOB=0
      IF(MOD(NT,NTSTEP).EQ.0) THEN
         IF(IDGLOB.EQ.0) CALL TRGLOB
         IDGLOB=1
         CALL TRSNAP
      ENDIF
      IF(MOD(NT,NGTSTP).EQ.0) THEN
         IF(IDGLOB.EQ.0) CALL TRGLOB
         IDGLOB=1
         CALL TRATOT
         CALL TRATOTN
      ENDIF
      IF(MOD(NT,NGRSTP).EQ.0) THEN
         IF(IDGLOB.EQ.0) CALL TRGLOB
         IDGLOB=1
         CALL TRATOG
      ENDIF
      IERR=0
      RETURN
    END SUBROUTINE tr_eval

!     ***********************************************************

!           COMPUTE MATRIX ELEMENTS

!     ***********************************************************

      SUBROUTINE TRMTRX(NEQRMAX)

      USE TRCOMM, ONLY : AEE, AKTB, AMZ, AX, AY, AZ, DD, DI, DT, &
     &                   DVRHOG, EPS0, LDAB, MDLCD, MDLEQB, MDLPCK, &
     &                   MLM, NEA, NEQM, NEQMAX, NRMAX, NSMAX, &
     &                   NSS, NST, NTUM, NVM, PI, PNB, PNF, RA, RIP, RIPA, &
     &                   RKEV, RMU0, RN, RR, RT, RTM, TAUB, TAUF, TAUK, VI, &
     &                   VV, X, XV, Y, YV, Z, ZV, RDPS, &
     &                   ABVRHOG, rkind
      USE TRCOM1, ONLY : A, B, C, D, PPA, PPB, PPC, RD
      IMPLICIT NONE
      INTEGER, INTENT(INOUT):: NEQRMAX
      INTEGER:: KL, MV, MVV, MW, MWMAX, NEQ, NEQ1, NR, NS, NS1, NSTN, NSW, &
     &             NV, NW
      REAL(rkind)   :: ADV, C1, COEF, COULOG, DV53, FADV, PRV, RDPA, RLP

! Boundary condition for magnetic diffusion equation

      IF(MDLEQB.NE.0) THEN
         IF(MDLCD.EQ.0) THEN
            RDPS=2.D0*PI*RMU0*RIP*1.D6*DVRHOG(NRMAX)/ABVRHOG(NRMAX)
         ELSE
            NEQ=1
            IF(NSS(NEQ).EQ.0) THEN
               RDPA=XV(NEQ,NRMAX)
            ELSE
               RDPA=0.D0
            ENDIF
            RLP=RA*(LOG(8.D0*RR/RA)-2.D0)
            RDPS= RDPA -4.D0*PI*PI*RMU0*RA*DVRHOG(NRMAX)/(RLP*ABVRHOG(NRMAX))
         ENDIF
      ENDIF

      COEF = AEE**4*1.D20/(3.D0*SQRT(2.D0*PI)*PI*EPS0**2)

      A(1:NEQMAX,1:NEQMAX,1:NRMAX)=0.D0
      B(1:NEQMAX,1:NEQMAX,1:NRMAX)=0.D0
      C(1:NEQMAX,1:NEQMAX,1:NRMAX)=0.D0

      D(1:NEQMAX,1:NRMAX)=0.D0
      PPA(1:NEQMAX,1:NRMAX)=0.D0
      PPB(1:NEQMAX,1:NRMAX)=0.D0
      PPC(1:NEQMAX,1:NRMAX)=0.D0

      VV(1:NEQMAX,1:NEQMAX,1:4,1:3)=0.D0
      DD(1:NEQMAX,1:NEQMAX,1:4,1:3)=0.D0
      VI(1:NEQMAX,1:NEQMAX,1:2,1:3)=0.D0
      DI(1:NEQMAX,1:NEQMAX,1:2,1:3)=0.D0

      MWMAX=4*NEQMAX-1
      AX(1:MWMAX,1:NEQMAX*NRMAX) = 0.D0

!     FADV = 0.D0  : Explicit scheme
!            0.5D0 : Crank-Nicolson scheme
!            1.D0  : Full implicit scheme

      FADV=1.D0

      PRV=(1.D0-FADV)*DT
      ADV=FADV*DT

!          +----------+
!    ***   |   NR=1   |   ***
!          +----------+

      NR=1
      NSW=1
      CALL TR_COEF_DECIDE(NR,NSW,DV53)

      DO NV=1,NEQMAX
      DO NW=1,NEQMAX
         A(NV,NW,NR) = 0.5D0*VI(NV,NW,1,NSW)+DI(NV,NW,1,NSW)
         B(NV,NW,NR) = 0.5D0*VI(NV,NW,1,NSW)-DI(NV,NW,1,NSW) &
     &                -0.5D0*VI(NV,NW,2,NSW)-DI(NV,NW,2,NSW)
         C(NV,NW,NR) =-0.5D0*VI(NV,NW,2,NSW)+DI(NV,NW,2,NSW)
      ENDDO
      ENDDO

      DO NS=1,NSMAX
         NEQ=NEA(NS,2)
         DO NS1=1,NSMAX
            IF(NS1.NE.NS) THEN
               NEQ1=NEA(NS1,2)
               C1=COEF/((RTM(NS)+RTM(NS1))**1.5D0*AMZ(NS)*AMZ(NS1)) &
     &                  *DV53 *COULOG(NS,NS1,RN(NR,1),RT(NR,NS))
               B(NEQ,NEQ, NR)=B(NEQ,NEQ, NR)-RN(NR,NS1)*C1
               B(NEQ,NEQ1,NR)=B(NEQ,NEQ1,NR)+RN(NR,NS )*C1
            ENDIF
         ENDDO
      ENDDO

      CALL TR_IONIZATION(NR)
      CALL TR_CHARGE_EXCHANGE(NR)

!     ***** RHS Vector *****

      DO NEQ=1,NEQMAX
         X(NEQMAX*(NR-1)+NEQ) = RD(NEQ,NR)*XV(NEQ,NR)+DT*D(NEQ,NR)
      ENDDO

      DO NW=1,NEQMAX
      DO NV=1,NEQMAX
         X(NEQMAX*(NR-1)+NV) = X(NEQMAX*(NR-1)+NV) &
     &                        +PRV*(+B(NV,NW,NR)*XV(NW,NR  ) &
     &                              +C(NV,NW,NR)*XV(NW,NR+1))
      ENDDO
      ENDDO

!     ***** Evolution of fast ion components *****

      Y(1,NR)=(1.D0-PRV/TAUB(NR))*YV(1,NR)+PNB(NR)*DT/(RKEV*1.D20)
      Y(2,NR)=(1.D0-PRV/TAUF(NR))*YV(2,NR)+PNF(NR)*DT/(RKEV*1.D20)
      AY(1,NR)=1.D0+ADV/TAUB(NR)
      AY(2,NR)=1.D0+ADV/TAUF(NR)

!          +---------------------+
!    ***   |   NR=2 to NRMAX-1   |   ***
!          +---------------------+

      NSW=2
      DO NR=2,NRMAX-1
         CALL TR_COEF_DECIDE(NR,NSW,DV53)

         DO NV=1,NEQMAX
         DO NW=1,NEQMAX
            A(NV,NW,NR) = 0.5D0*VI(NV,NW,1,NSW)+DI(NV,NW,1,NSW)
            B(NV,NW,NR) = 0.5D0*VI(NV,NW,1,NSW)-DI(NV,NW,1,NSW) &
     &                   -0.5D0*VI(NV,NW,2,NSW)-DI(NV,NW,2,NSW)
            C(NV,NW,NR) =-0.5D0*VI(NV,NW,2,NSW)+DI(NV,NW,2,NSW)
         ENDDO
         ENDDO

         DO NS=1,NSMAX
            NEQ=NEA(NS,2)
            DO NS1=1,NSMAX
               IF(NS1.NE.NS) THEN
                  NEQ1=NEA(NS1,2)
                  C1=COEF/((RTM(NS)+RTM(NS1))**1.5D0*AMZ(NS)*AMZ(NS1)) &
     &                     *DV53*COULOG(NS,NS1,RN(NR,1),RT(NR,NS))
                  B(NEQ,NEQ, NR)=B(NEQ,NEQ, NR)-RN(NR,NS1)*C1
                  B(NEQ,NEQ1,NR)=B(NEQ,NEQ1,NR)+RN(NR,NS )*C1
               ENDIF
            ENDDO
         ENDDO

         CALL TR_IONIZATION(NR)
         CALL TR_CHARGE_EXCHANGE(NR)

!     ***** RHS Vector *****

         DO NEQ=1,NEQMAX
            X(NEQMAX*(NR-1)+NEQ) = RD(NEQ,NR)*XV(NEQ,NR)+DT*D(NEQ,NR)
         ENDDO

         DO NW=1,NEQMAX
         DO NV=1,NEQMAX
            X(NEQMAX*(NR-1)+NV) = X(NEQMAX*(NR-1)+NV) &
     &                          +PRV*(A(NV,NW,NR)*XV(NW,NR-1) &
     &                               +B(NV,NW,NR)*XV(NW,NR  ) &
     &                               +C(NV,NW,NR)*XV(NW,NR+1))
         ENDDO
         ENDDO

!     ***** Evolution of fast ion components *****

         Y(1,NR)=(1.D0-PRV/TAUB(NR))*YV(1,NR)+PNB(NR)*DT/(RKEV*1.D20)
         Y(2,NR)=(1.D0-PRV/TAUF(NR))*YV(2,NR)+PNF(NR)*DT/(RKEV*1.D20)
         AY(1,NR)=1.D0+ADV/TAUB(NR)
         AY(2,NR)=1.D0+ADV/TAUF(NR)
!
      ENDDO

!          +--------------+
!    ***   |   NR=NRMAX   |   ***
!          +--------------+

      NR=NRMAX
      NSW=3
      CALL TR_COEF_DECIDE(NR,NSW,DV53)

      DO NV=1,NEQMAX
      DO NW=1,NEQMAX
         A(NV,NW,NR) = 0.5D0*VI(NV,NW,1,NSW)+DI(NV,NW,1,NSW) &
     &                +DI(NV,NW,2,NSW)/3.D0
         B(NV,NW,NR) = 0.5D0*VI(NV,NW,1,NSW)-DI(NV,NW,1,NSW) &
     &                -DI(NV,NW,2,NSW)*3.D0
         C(NV,NW,NR) =-0.0D0
      ENDDO
      ENDDO

      DO NS=1,NSMAX
         NEQ=NEA(NS,2)
         DO NS1=1,NSMAX
            IF(NS1.NE.NS) THEN
               NEQ1=NEA(NS1,2)
               C1=COEF/((RTM(NS)+RTM(NS1))**1.5D0*AMZ(NS)*AMZ(NS1)) &
     &                 *DV53*COULOG(NS,NS1,RN(NR,1),RT(NR,NS))
               B(NEQ,NEQ, NR)=B(NEQ,NEQ, NR)-RN(NR,NS1)*C1
               B(NEQ,NEQ1,NR)=B(NEQ,NEQ1,NR)+RN(NR,NS )*C1
            ENDIF
         ENDDO
      ENDDO

      CALL TR_IONIZATION(NR)
      CALL TR_CHARGE_EXCHANGE(NR)

!     ***** RHS Vector *****

      DO NEQ=1,NEQMAX
         X(NEQMAX*(NR-1)+NEQ) = RD(NEQ,NR)*XV(NEQ,NR)+DT*D(NEQ,NR)
      ENDDO

      DO NW=1,NEQMAX
      DO NV=1,NEQMAX
         X(NEQMAX*(NR-1)+NV) = X(NEQMAX*(NR-1)+NV) &
     &                        +PRV*(A(NV,NW,NR)*XV(NW,NR-1) &
     &                             +B(NV,NW,NR)*XV(NW,NR  ) )
      ENDDO
      ENDDO

!     ***** Evolution of fast ion components *****

      Y(1,NR)=(1.D0-PRV/TAUB(NR))*YV(1,NR)+PNB(NR)*DT/(RKEV*1.D20)
      Y(2,NR)=(1.D0-PRV/TAUF(NR))*YV(2,NR)+PNF(NR)*DT/(RKEV*1.D20)
      AY(1,NR)=1.D0+ADV/TAUB(NR)
      AY(2,NR)=1.D0+ADV/TAUF(NR)

!     ***** Band Matrix *****

      CALL TR_BAND_GEN(NEQRMAX,ADV)
      IF(MDLPCK.NE.0) CALL TR_BAND_LAPACK(A,B,C,AX,NRMAX,NEQRMAX,NVM,LDAB,MLM)

!     ***** RHS Vector Reduce *****

      DO NR=1,NRMAX
      DO NEQ=1,NEQMAX
         NSTN=NST(NEQ)
         IF(NSTN.NE.0) THEN
            X(NEQRMAX*(NR-1)+NSTN) = X(NEQMAX*(NR-1)+NEQ) &
                 +PPA(NEQ,NR)+PPB(NEQ,NR)+PPC(NEQ,NR)
         ENDIF
      ENDDO
      ENDDO

!     ***** Surface Boundary Condition for Bp *****

      NEQ=1
      IF(MDLEQB.NE.0) THEN
         IF(NEA(0,0).EQ.1) THEN
            MVV=NEQRMAX*(NRMAX-1)+NEQ
            IF(MDLPCK.EQ.0) THEN
               DO MW=1,MWMAX
                  AX(MW,MVV)=0.D0
               ENDDO
               AX(2*NEQRMAX,MVV)=RD(NEQ,NRMAX)
            ELSE
               KL=2*NEQRMAX-1
               DO MW=3*NEQRMAX,NEQRMAX+1,-1
                  DO MV=NEQRMAX*NRMAX-2*NEQRMAX+KL,NEQRMAX*NRMAX
                     AX(MW,MV)=0.D0
                  ENDDO
               ENDDO
               AX(2*NEQRMAX+KL,MVV)=RD(NEQ,NRMAX)
            ENDIF
            X(MVV)=RDPS
         ENDIF
      ENDIF

!      nr=nrmax
!      WRITE(6,'(A,I6)') '** NRMAX:',nrmax
!      DO NV=1,NEQMAX
!         DO NW=1,NEQMAX
!            IF(NW.EQ.1) THEN
!               WRITE(6,'(A,2I6,4ES12.4)') 'NV,NW,A,B,C,X=',NV,NW, &
!                    A(NV,NW,NR),B(NV,NW,NR),C(NV,NW,NR),X(NEQMAX*(NR-1)+NV)
!            ELSE
!               WRITE(6,'(A,2I6,3ES12.4)') 'NV,NW,A,B,C  =',NV,NW, &
!                    A(NV,NW,NR),B(NV,NW,NR),C(NV,NW,NR)
!            END IF
!         END DO
!      END DO
               
      RETURN
      END SUBROUTINE TRMTRX

!     ***********************************************************

!           BAND MATRIX GENERATOR

!     ***********************************************************

      SUBROUTINE TR_BAND_GEN(NEQRMAX,ADV)

      USE TRCOMM, ONLY : AX, MDLPCK, NEQM, NEQMAX, NRMAX, XV, rkind
      USE TRCOM1, ONLY : A, B, C, PPA, PPB, PPC, RD
      IMPLICIT NONE
      INTEGER,INTENT(INOUT):: NEQRMAX
      REAL(rkind)   ,INTENT(IN)   :: ADV
      INTEGER:: NEQ, NR, NV, NW

      DO NR=1,NRMAX

         A(1:NEQMAX,1:NEQMAX,NR) = -ADV*A(1:NEQMAX,1:NEQMAX,NR)
         B(1:NEQMAX,1:NEQMAX,NR) = -ADV*B(1:NEQMAX,1:NEQMAX,NR)
         C(1:NEQMAX,1:NEQMAX,NR) = -ADV*C(1:NEQMAX,1:NEQMAX,NR)

         DO NEQ=1,NEQMAX
            B(NEQ,NEQ,NR)=B(NEQ,NEQ,NR)+RD(NEQ,NR)
         ENDDO

         IF(NR.EQ.1) THEN
            CALL TR_BANDREDUCE(B,XV,PPB,NR,NR,  NEQRMAX)
            CALL TR_BANDREDUCE(C,XV,PPC,NR,NR+1,NEQRMAX)
         ELSEIF(NR.EQ.NRMAX) THEN
            CALL TR_BANDREDUCE(A,XV,PPA,NR,NR-1,NEQRMAX)
            CALL TR_BANDREDUCE(B,XV,PPB,NR,NR,  NEQRMAX)
         ELSE
            CALL TR_BANDREDUCE(A,XV,PPA,NR,NR-1,NEQRMAX)
            CALL TR_BANDREDUCE(B,XV,PPB,NR,NR,  NEQRMAX)
            CALL TR_BANDREDUCE(C,XV,PPC,NR,NR+1,NEQRMAX)
         ENDIF

         IF(MDLPCK.EQ.0) THEN
!            forall(NV=1:NEQRMAX,NW=1:NEQRMAX)
!               AX(  NEQRMAX+NW-NV,NEQRMAX*(NR-1)+NV) = A(NV,NW,NR)
!               AX(2*NEQRMAX+NW-NV,NEQRMAX*(NR-1)+NV) = B(NV,NW,NR)
!               AX(3*NEQRMAX+NW-NV,NEQRMAX*(NR-1)+NV) = C(NV,NW,NR)
!            end forall
            DO NV=1,NEQRMAX
               DO NW=1,NEQRMAX
                  AX(  NEQRMAX+NW-NV,NEQRMAX*(NR-1)+NV) = A(NV,NW,NR)
                  AX(2*NEQRMAX+NW-NV,NEQRMAX*(NR-1)+NV) = B(NV,NW,NR)
                  AX(3*NEQRMAX+NW-NV,NEQRMAX*(NR-1)+NV) = C(NV,NW,NR)
               ENDDO
            ENDDO
         ENDIF

      ENDDO

      RETURN
      END SUBROUTINE TR_BAND_GEN

!     ***********************************************************

!           ADDITIONAL BAND MATRIX GENERATOR FOR LAPACK

!     ***********************************************************

      SUBROUTINE TR_BAND_LAPACK(A,B,C,AX,NRMAX,NEQRMAX,NVM,LDAB,MLM)

      USE TRCOMM,ONLY: rkind
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: NRMAX, NEQRMAX, NVM, LDAB, MLM
      REAL(rkind),DIMENSION(NVM,NVM,NRMAX),INTENT(IN) :: A, B, C
      REAL(rkind),DIMENSION(LDAB,MLM)  ,INTENT(OUT) :: AX
      INTEGER:: KL, NR, NV, NW

      KL=2*NEQRMAX-1
      NR=1
      forall(NV=1:NEQRMAX,NW=1:NEQRMAX)
         AX(3*NEQRMAX+NW-NV+KL,NEQRMAX*(NR-1)+NV) = A(NW,NV,NR+1)
         AX(2*NEQRMAX+NW-NV+KL,NEQRMAX*(NR-1)+NV) = B(NW,NV,NR  )
      end forall
      forall(NR=2:NRMAX-1,NV=1:NEQRMAX,NW=1:NEQRMAX)
         AX(3*NEQRMAX+NW-NV+KL,NEQRMAX*(NR-1)+NV) = A(NW,NV,NR+1)
         AX(2*NEQRMAX+NW-NV+KL,NEQRMAX*(NR-1)+NV) = B(NW,NV,NR  )
         AX(  NEQRMAX+NW-NV+KL,NEQRMAX*(NR-1)+NV) = C(NW,NV,NR-1)
      end forall
      NR=NRMAX
      forall(NV=1:NEQRMAX,NW=1:NEQRMAX)
         AX(2*NEQRMAX+NW-NV+KL,NEQRMAX*(NR-1)+NV) = B(NW,NV,NR  )
         AX(  NEQRMAX+NW-NV+KL,NEQRMAX*(NR-1)+NV) = C(NW,NV,NR-1)
      end forall

      RETURN
      END SUBROUTINE TR_BAND_LAPACK

!     ***********************************************************

!           BAND REDUCER

!     ***********************************************************

      SUBROUTINE TR_BANDREDUCE(A,XR,P,NR,NR1,NEQRMAX)


      USE TRCOMM, ONLY : NEQM, NEQMAX, NEQMAXM, NNS, NRMAX, NST, rkind
      IMPLICIT NONE
      REAL(rkind),DIMENSION(NEQMAXM,NEQMAXM,NRMAX),INTENT(INOUT):: A
      REAL(rkind),DIMENSION(NEQMAXM,NRMAX)     ,INTENT(IN)   :: XR
      REAL(rkind),DIMENSION(NEQMAXM,NRMAX)     ,INTENT(INOUT):: P
      INTEGER                     ,INTENT(IN)  :: NR,NR1
      INTEGER                     ,INTENT(OUT)  :: NEQRMAX
      INTEGER :: L, LOOP, NBSIZE, NEQ, NNSN, NNSN1, NNSOLD, NV, NW
      REAL(rkind),DIMENSION(NEQMAXM,NEQMAXM) :: AA, AL
      REAL(rkind),DIMENSION(NEQMAXM)         :: AM

      DO NV=1,NEQMAX
         DO NW=1,NEQMAX
            AA(NV,NW)=A(NV,NW,NR)
         ENDDO
      ENDDO
      AM(1:NEQMAX)=0.D0

      NBSIZE=NEQMAX
      LOOP=0
      DO NEQ=1,NEQM
         NNSN=NNS(NEQ)
         NNSN1=NNS(NEQ)
         IF(NNSN.EQ.0) THEN
            NEQRMAX=NEQMAX-LOOP
            GOTO 1000
         ENDIF
         IF(LOOP.EQ.0) THEN
            NNSOLD=NNSN
         ELSE
            IF(NNSN.GT.NNSOLD) THEN
               NNSN=NNSN-LOOP
               NNSOLD=NNSN
            ELSE
               STOP 'XX TR_BANDREDUCE: ERROR'
            ENDIF
         ENDIF
         LOOP=LOOP+1

!     /* Obtaining correction terms */

         DO NV=1,NEQMAX
            AL(LOOP,NV)=AA(NV,NNSN1)*XR(NNSN1,NR1)
         ENDDO

!     /* Band reducer*/

         IF(NNSN.NE.NBSIZE) THEN
            NBSIZE=NBSIZE-1
            DO NV=1,NBSIZE+1
               DO NW=NNSN,NBSIZE
                  A(NV,NW,NR)=A(NV,NW+1,NR)
               ENDDO
            ENDDO
            A(1:NBSIZE+1,NBSIZE+1,NR)=0.D0

            DO NW=1,NBSIZE+1
               DO NV=NNSN,NBSIZE
                  A(NV,NW,NR)=A(NV+1,NW,NR)
               ENDDO
            ENDDO
            A(NBSIZE+1,1:NBSIZE+1,NR)=0.D0
         ELSE
            NBSIZE=NBSIZE-1
            A(1:NBSIZE+1,NBSIZE+1,NR)=0.D0
            A(NBSIZE+1,1:NBSIZE+1,NR)=0.D0
         ENDIF
      ENDDO

 1000 CONTINUE

!     /* Consummation */

      DO L=1,LOOP
         AM(1:NEQMAX)=AM(1:NEQMAX)+AL(L,1:NEQMAX)
      ENDDO
      DO NEQ=1,NEQMAX
         IF(NST(NEQ).NE.0) P(NEQ,NR)=P(NEQ,NR)-AM(NEQ)
      ENDDO

      RETURN
      END SUBROUTINE TR_BANDREDUCE

!     ***********************************************************

!           CONVERSION FROM ANT TO XV

!     ***********************************************************

      SUBROUTINE TRATOX

      USE TRCOMM, ONLY : AKTB, AMM, NEQMAX, NRMAX, NSMAX, NSS, NSV, PA, RDP, RN, RT, RU, RW, XV, YV, ZV
      IMPLICIT NONE
      INTEGER:: NEQ, NR, NS, NSSN, NSVN


      DO NR=1,NRMAX
         DO NEQ=1,NEQMAX
            NSVN=NSV(NEQ)
            NSSN=NSS(NEQ)
            IF(NSVN.EQ.0) THEN
               XV(NEQ,NR) = RDP(NR)
            ELSEIF(NSVN.EQ.1) THEN
               XV(NEQ,NR) = RN(NR,NSSN)
            ELSEIF(NSVN.EQ.2) THEN
               XV(NEQ,NR) = RN(NR,NSSN)*RT(NR,NSSN)
            ELSEIF(NSVN.EQ.3) THEN
               XV(NEQ,NR) = PA(NSSN)*AMM*RN(NR,NSSN)*RU(NR,NSSN)
            ENDIF
         ENDDO
         YV(  1,NR) = RW(NR,1)
         YV(  2,NR) = RW(NR,2)
      ENDDO

      RETURN
      END SUBROUTINE TRATOX

!     ***********************************************************

!           CONVERT FROM XV TO ANT

!     ***********************************************************

      SUBROUTINE TRXTOA

      USE TRCOMM, ONLY : &
           AKTB, AMM, ANC, ANFE, ANNU, DR, MDLEQE, MDLEQN, &
           NEQMAX, NRMAX, NROMAX, NSM, NSMAX, NSS, NSV, &
           PA, PZ, PZC, PZFE, RDP, RN, RPSI, RT, RU, RW, &
           XV, YV, ZV, RDPVRHOG, DVRHOG, rkind
      IMPLICIT NONE
      INTEGER:: N, NEQ, NEQ1, NR, NS, NSSN, NSSN1, NSVN, NSVN1
      REAL(rkind)   :: SUM


      DO NR=1,NRMAX
         DO NEQ=1,NEQMAX
            NSSN=NSS(NEQ)
            NSVN=NSV(NEQ)
            IF(NSVN.EQ.0) THEN
               RDP(NR) = XV(NEQ,NR)
               RDPVRHOG(NR)=RDP(NR) / DVRHOG(NR)
            ELSEIF(NSVN.EQ.1) THEN
               IF(MDLEQN.NE.0) THEN
                  IF(NSSN.EQ.1.AND.MDLEQE.EQ.0) THEN
                     RN(NR,NSSN) = 0.D0
                     DO NEQ1=1,NEQMAX
                        NSSN1=NSS(NEQ1)
                        NSVN1=NSV(NEQ1)
                        IF(NSVN1.EQ.1.AND.NSSN1.NE.1) THEN
                           RN(NR,NSSN) = RN(NR,NSSN)+ABS(PZ(NSSN1))*XV(NEQ1,NR)
                        ENDIF
                     ENDDO
                     RN(NR,NSSN) = RN(NR,NSSN)+PZFE(NR)*ANFE(NR)+PZC(NR)*ANC(NR)
!     I think the above expression is obsolete.
                  ELSEIF(NSSN.EQ.2.AND.MDLEQE.EQ.2) THEN
                     RN(NR,NSSN) = 0.D0
                     DO NEQ1=1,NEQMAX
                        NSSN1=NSS(NEQ1)
                        NSVN1=NSV(NEQ1)
                        IF(NSVN1.EQ.1.AND.NSSN1.NE.2) THEN
                           RN(NR,NSSN) = RN(NR,NSSN)+ABS(PZ(NSSN1))*XV(NEQ1,NR)
                        ENDIF
                     ENDDO
                     RN(NR,NSSN) = RN(NR,NSSN)+PZFE(NR)*ANFE(NR)+PZC(NR)*ANC(NR)
                  ELSE
                     RN(NR,NSSN) = XV(NEQ,NR)
                  ENDIF
               ELSE
                  RN(NR,NSSN) = XV(NEQ,NR)
               ENDIF
            ELSEIF(NSVN.EQ.2) THEN
               IF(NSSN.NE.NSM) THEN
                  RT(NR,NSSN) = XV(NEQ,NR)/RN(NR,NSSN)
               ELSE
                  IF(RN(NR,NSM).LT.1.D-70) THEN
                     RT(NR,NSM) = 0.D0
                  ELSE
                     RT(NR,NSM) = XV(NEQ,NR)/RN(NR,NSM)
                  ENDIF
               ENDIF
            ELSEIF(NSVN.EQ.3) THEN
               RU(NR,NSSN) = XV(NEQ,NR)/(PA(NSSN)*AMM*RN(NR,NSSN))
            ENDIF
         ENDDO
         RW(NR,1) = YV(1,NR)
         RW(NR,2) = YV(2,NR)
         IF(NSMAX.GE.8) THEN
            ANNU(NR) = RN(NR,7)+RN(NR,8)
         END IF
!!!         BP(NR)   = AR1RHOG(NR)*RDP(NR)/RR
      ENDDO

      SUM=0.D0
      DO NR=1,NROMAX
         SUM=SUM+RDP(NR)*DR
         RPSI(NR)=SUM
      ENDDO

      N=NRMAX
!      CALL PLDATA_SETP(N,RN,RT,RU)

      RETURN
      END SUBROUTINE TRXTOA

!     ***********************************************************

!           SAVE TRANSIENT DATA FOR GRAPHICS

!     ***********************************************************

      SUBROUTINE TRATOTN

      USE TRCOMM, ONLY : GTS, GVT, IZERO, MDLST, NGPST, NGST, NGT, NRMAX, NTM, RT, T, TSST
      IMPLICIT NONE
      INTEGER:: K, L, IX
      REAL   :: GUCLIP


      IF(NGT.GE.NTM) RETURN

      IF(MDLST.EQ.0) RETURN

      IF(TSST.GT.T) RETURN

      NGST=NGST+1
      K=0
      L=IZERO
      IX=INT((NRMAX+1-IZERO)/(NGPST-1))

      GTS(NGST) = GUCLIP(T)

  100 GVT(NGST,K+89) = GUCLIP(RT(L,1))
      GVT(NGST,K+90) = GUCLIP(RT(L,2))
      GVT(NGST,K+91) = GUCLIP(RT(L,3))
      GVT(NGST,K+92) = GUCLIP(RT(L,4))

      K=K+4
      L=L+IX
      IF(L.LE.IZERO+(NGPST-1)*IX) GOTO 100

      RETURN
      END SUBROUTINE TRATOTN

!     ***********************************************************

!           CHECK NEGATIVE DENSITY OR TEMPERATURE

!     ***********************************************************

      SUBROUTINE TRCHCK(ICHCK)

      USE TRCOMM, ONLY : NEQMAX, NRMAX, NSS, NSV, NT, RT
      IMPLICIT NONE
      INTEGER,INTENT(OUT):: ICHCK
      INTEGER :: IND, NEQ, NR, NSSN

      ICHCK = 0
      DO NEQ=1,NEQMAX
         IF(NSV(NEQ).EQ.2) THEN
            DO NR=1,NRMAX
               IF(RT(NR,NSS(NEQ)).LT.0.D0) &
                    & write(6,*) NT,NR,NEQ,NSS(NEQ),RT(NR,NSS(NEQ))
               IF(RT(NR,NSS(NEQ)).LT.0.D0) GOTO 100
            ENDDO
         ENDIF
      ENDDO
      GOTO 9000

  100 CONTINUE
      IND=0
      WRITE(6,*) 'XX ERROR : NEGATIVE TEMPERATURE AT STEP ',NT
      DO NEQ=1,NEQMAX
         NSSN=NSS(NEQ)
         IF(NSSN.EQ.1.AND.IND.NE.1) THEN
            WRITE(6,*) '     TE (',NR,')=',RT(NR,NSSN)
            IND=1
         ELSEIF(NSSN.EQ.2.AND.IND.NE.2) THEN
            WRITE(6,*) '     TD (',NR,')=',RT(NR,NSSN)
            IND=2
         ELSEIF(NSSN.EQ.3.AND.IND.NE.3) THEN
            WRITE(6,*) '     TT (',NR,')=',RT(NR,NSSN)
            IND=3
         ELSEIF(NSSN.EQ.4.AND.IND.NE.4) THEN
            WRITE(6,*) '     TA (',NR,')=',RT(NR,NSSN)
            IND=4
         ELSEIF(NSSN.EQ.5.AND.IND.NE.5) THEN
            WRITE(6,*) '     TI1(',NR,')=',RT(NR,NSSN)
            IND=5
         ELSEIF(NSSN.EQ.6.AND.IND.NE.6) THEN
            WRITE(6,*) '     TI2(',NR,')=',RT(NR,NSSN)
            IND=6
!         ELSEIF(NSSN.EQ.7.AND.IND.NE.7) THEN
!            WRITE(6,*) '     TNC(',NR,')=',RT(NR,NSSN)
!            IND=7
!         ELSEIF(NSSN.EQ.8.AND.IND.NE.8) THEN
!            WRITE(6,*) '     TNH(',NR,')=',RT(NR,NSSN)
!            IND=8
         ENDIF
      ENDDO
      ICHCK=1

 9000 RETURN
      END SUBROUTINE TRCHCK

!     ***********************************************************

!           DECIDE COEEFICIENTS FOR EQUATIONS

!     ***********************************************************

      SUBROUTINE TR_COEF_DECIDE(NR,NSW,DV53)

      USE trcomm
      USE TRCOM1, ONLY : D, RD
      USE libitp
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: NR, NSW
      REAL(rkind)   , INTENT(OUT)::DV53
      INTEGER:: &
           IND, NEQ, NEQ1, NEQLMAX, NF, NI, NJ, NMK, NMKL, NO, NRF, NRJ, &
           NSSN, NSSN1, NSVN, NSVN1, NV, NW
      REAL(rkind)   :: &
           C83, CC, DISUMN, DISUMT1, DISUMT2, DV11, DV23, VISUMN, &
           VISUMT1, VISUMT2
      REAL(rkind),DIMENSION(2) :: F2C
      REAL(rkind),DIMENSION(4) :: SIG, FCB


      DV11=DVRHO(NR)
      DV23=DVRHO(NR)**(2.D0/3.D0)
      DV53=DVRHO(NR)**(5.D0/3.D0)
!!      CC=2.5D0
      CC=1.5D0
      C83=8.D0/3.D0

      DO NEQ=1,NEQMAX  ! *** NEQ ***

         NSSN=NSS(NEQ)
         NSVN=NSV(NEQ)
         IF(NSSN.NE.0) RTM(NSSN)=RT(NR,NSSN)*RKEV/(PA(NSSN)*AMM)

!     /* Coefficients of left hand side variables */

         IF(NSVN.EQ.0) THEN
            RD(NEQ,NR)=1.D0
         ELSEIF(NSVN.EQ.1) THEN
            RD(NEQ,NR)=DV11
         ELSEIF(NSVN.EQ.2) THEN
            RD(NEQ,NR)=1.5D0*DV53
         ELSEIF(NSVN.EQ.3) THEN
            RD(NEQ,NR)=DV11
         ENDIF

!     /* Coefficients of main variables */

         IF(NSVN.EQ.0) THEN ! for mag. diff. eq.
            IF(NSW.NE.3) THEN ! NR/=NRMAX
               DO NF=1,4
                  NRF=NR+(NF-2)
                  IF(NRF.EQ.0.OR.NRF.GT.NRMAX) THEN
                     FCB(NF)=0.D0
                  ELSE
                     FCB(NF)=DVRHO(NRF)*ABRHO(NRF)/TTRHO(NRF)
                  ENDIF
               ENDDO
            ELSE ! NR=NRMAX
               DO NF=1,4
                  FCB(NF)=0.D0
               ENDDO
            ENDIF
            IF(NR.EQ.NRMAX-1) FCB(4)=AITKEN2P(RM(NR+1)+DR, &
     &                        DVRHO(NR+1)*ABRHO(NR+1)/TTRHO(NR+1), &
     &                        DVRHO(NR  )*ABRHO(NR  )/TTRHO(NR  ), &
     &                        DVRHO(NR-1)*ABRHO(NR-1)/TTRHO(NR-1), &
     &                        RM(NR+1),RM(NR),RM(NR-1))
         ENDIF
!
         IF(NSVN.NE.0) THEN ! for no mag. diff. eqs.
            DO NMK=1,3
               IF (NSW.EQ.2.OR.NSW.NE.NMK) THEN
                  NMKL=NMK-2
                  FA(NMK,NSW)=DVRHO(NR+NMKL)*AR1RHO(NR+NMKL)/DR
                  FB(NMK,NSW)=DVRHO(NR+NMKL)*AR2RHO(NR+NMKL)/(DR*DR)
               ELSE
                  FA(NMK,NSW)=0.D0
                  FB(NMK,NSW)=0.D0
               ENDIF
            ENDDO
         ELSE ! for mag. diff. eq.
            DO NMK=1,3
               FC(NMK,NSW)=0.5D0*(FCB(NMK)+FCB(NMK+1))/(DR*DR)
            ENDDO
         ENDIF
!
         SIG(1)= 2.D0
         SIG(2)= 2.D0
         SIG(3)=-2.D0
         SIG(4)=-2.D0
         DO NMK=1,4  ! *** NMK ***
            IF(NMK.EQ.1) THEN
               NI=1
               NJ=1
            ELSEIF(NMK.EQ.2) THEN
               NI=2
               NJ=2
            ELSEIF(NMK.EQ.3) THEN
               NI=2
               NJ=1
            ELSE
               NI=3
               NJ=2
            ENDIF
            NRJ=NR+(NJ-2)

            IF(MDDIAG.EQ.0) THEN  ! *** MDDIAG ***

            IF(NSVN.EQ.0) THEN
               IF(NSW.NE.3) THEN
                  F2C(NJ)=ETA(NRJ+1)*TTRHO(NRJ+1)/(RMU0*ARRHO(NRJ+1)*DVRHO(NRJ+1))
                  VV(NEQ,NEQ  ,NMK,NSW)= SIG(NMK)*F2C(NJ)*FC(NI,NSW)
                  DD(NEQ,NEQ  ,NMK,NSW)=          F2C(NJ)*FC(NI,NSW)
               ENDIF
            ELSEIF(NSVN.EQ.1) THEN
               IF(NSW.EQ.2.OR.NSW.NE.NI.AND.NRJ.NE.0) THEN
                  VV(NEQ,NEQ  ,NMK,NSW)= FA(NI,NSW)*AVP(NRJ,NSSN)
                  DD(NEQ,NEQ  ,NMK,NSW)= FB(NI,NSW)*AD(NRJ,NSSN)
               ENDIF
            ELSEIF(NSVN.EQ.2) THEN
               IF(NSW.EQ.2.OR.NSW.NE.NI.AND.NRJ.NE.0) THEN
                  IF(MDLFLX.EQ.0) THEN
                  VV(NEQ,NEQ  ,NMK,NSW)= FA(NI,NSW)*(AVK(NRJ,NSSN)+AVP(NRJ,NSSN)*CC)*DV23
                  ELSE
                  VV(NEQ,NEQ  ,NMK,NSW)=(FA(NI,NSW)*AVK(NRJ,NSSN)+CC*RGFLX(NRJ,NSSN)/(RNV(NRJ,NSSN)*DR))*DV23
                  ENDIF
                  DD(NEQ,NEQ  ,NMK,NSW)= FB(NI,NSW)*AK(NRJ,NSSN)*DV23
                  VV(NEQ,NEQ-1,NMK,NSW)= 0.D0*DV23
                  DD(NEQ,NEQ-1,NMK,NSW)= FB(NI,NSW)*(AD(NRJ,NSSN)*CC-AK(NRJ,NSSN))*RTV(NRJ,NSSN)*DV23
               ENDIF
            ELSEIF(NSVN.EQ.3) THEN
               IF(NSW.EQ.2.OR.NSW.NE.NI.AND.NRJ.NE.0) THEN
                  VV(NEQ,NEQ  ,NMK,NSW)= FA(NI,NSW)*VOID
                  DD(NEQ,NEQ  ,NMK,NSW)= FB(NI,NSW)*VOID
               ENDIF
            ENDIF

            ELSE  ! *** MDDIAG ***

!     *** OFF-DIAGONAL TRANSPORT ***
!     *
!     This expression below is the interim way to solve the problem.
!     That is because this simple expression cannot be used if we
!       consider impurity transports.
            IF(MDLEQ0.EQ.1) THEN
               NEQLMAX=NEQMAX-2 ! This means "Local NEQMAX".
            ELSE
               NEQLMAX=NEQMAX
            ENDIF

            IF(NSVN.EQ.0) THEN
               IF(NSW.NE.3) THEN
                  F2C(NJ)=ETA(NRJ+1)*TTRHO(NRJ+1)/(RMU0*ARRHO(NRJ+1)*DVRHO(NRJ+1))
                  VV(NEQ,NEQ  ,NMK,NSW)= SIG(NMK)*F2C(NJ)*FC(NI,NSW)
                  DD(NEQ,NEQ  ,NMK,NSW)=          F2C(NJ)*FC(NI,NSW)
               ENDIF
            ELSEIF(NSVN.EQ.1) THEN
               IF(NSSN.EQ.7.OR.NSSN.EQ.8) THEN
                  IF(NSW.EQ.2.OR.NSW.NE.NI.AND.NRJ.NE.0) THEN
                     VV(NEQ,NEQ  ,NMK,NSW)= FA(NI,NSW)*AVP(NRJ,NSSN)
                     DD(NEQ,NEQ  ,NMK,NSW)= FB(NI,NSW)*AD(NRJ,NSSN)
                  ENDIF
               ELSE
               IF(NSW.EQ.2.OR.NSW.NE.NI.AND.NRJ.NE.0) THEN
                  VV(NEQ,NEQ  ,NMK,NSW)= FA(NI,NSW)*(AVP(NRJ,NSSN)+RGFLS(NRJ,5,NSSN)/RNV(NRJ,NSSN))
                  DO NEQ1=1,NEQLMAX
                     NSSN1=NSS(NEQ1)
                     NSVN1=NSV(NEQ1)
                     IF(NSVN1.EQ.1) THEN
                        DD(NEQ,NEQ1,NMK,NSW)= FB(NI,NSW)*RNV(NRJ,NSSN)*ADLD(NRJ,NSSN1,NSSN)/RNV(NRJ,NSSN1)
                     ELSEIF(NSVN1.EQ.2) THEN
                        DD(NEQ,NEQ1,NMK,NSW)= FB(NI,NSW)*RNV(NRJ,NSSN)*ADLP(NRJ,NSSN1,NSSN)/RPV(NRJ,NSSN1)
                     ENDIF
                  ENDDO
               ENDIF
               ENDIF
            ELSEIF(NSVN.EQ.2) THEN
               IF(NSW.EQ.2.OR.NSW.NE.NI.AND.NRJ.NE.0) THEN
                  IF(MDLFLX.EQ.0) THEN
                     VV(NEQ,NEQ  ,NMK,NSW)= FA(NI,NSW)*DV23 &
     &                    *( AVK(NRJ,NSSN)+AVP(NRJ,NSSN)*CC &
     &                      +RQFLS(NRJ,5,NSSN)   /RPV(NRJ,NSSN) &
     &                      +RGFLS(NRJ,5,NSSN)*CC/RNV(NRJ,NSSN))
                  ELSE
                     VV(NEQ,NEQ  ,NMK,NSW)= DV23*(FA(NI,NSW) &
     &                    *( AVK(NRJ,NSSN) &
     &                    + RQFLS(NRJ,5,NSSN)   / RPV(NRJ,NSSN) &
     &                    + RGFLS(NRJ,5,NSSN)*CC/ RNV(NRJ,NSSN)) &
     &                    +(RGFLX(NRJ,NSSN)  *CC/(RNV(NRJ,NSSN)*DR)))
                  ENDIF
                  DO NEQ1=1,NEQLMAX
                     NSSN1=NSS(NEQ1)
                     NSVN1=NSV(NEQ1)
                     IF(NSVN1.EQ.1) THEN
                        DD(NEQ,NEQ1,NMK,NSW)= FB(NI,NSW)*DV23 &
     &                       *  RPV(NRJ,NSSN) &
     &                       *( AKLD(NRJ,NSSN1,NSSN) &
     &                         +ADLD(NRJ,NSSN1,NSSN)*CC) &
     &                       /  RNV(NRJ,NSSN1)
                     ELSEIF(NSVN1.EQ.2) THEN
                        DD(NEQ,NEQ1,NMK,NSW)= FB(NI,NSW)*DV23 &
     &                       *  RPV(NRJ,NSSN) &
     &                       *( AKLP(NRJ,NSSN1,NSSN) &
     &                         +ADLP(NRJ,NSSN1,NSSN)*CC) &
     &                       /  RPV(NRJ,NSSN1)
                     ENDIF
                  ENDDO
               ENDIF
            ELSEIF(NSVN.EQ.3) THEN
               IF(NSW.EQ.2.OR.NSW.NE.NI.AND.NRJ.NE.0) THEN
                  VV(NEQ,NEQ  ,NMK,NSW)= FA(NI,NSW)*VOID
                  DD(NEQ,NEQ  ,NMK,NSW)= FB(NI,NSW)*VOID
               ENDIF
            ENDIF
!     *
!     **************

            ENDIF  ! *** MDDIAG ***
         ENDDO  ! *** NMK ***

      ENDDO  ! *** NEQ ***

      IF(MDLEQE.EQ.1) THEN
         IND=0
         DO NEQ=1,NEQMAX
            NSSN=NSS(NEQ)
            NSVN=NSV(NEQ)
            IF(NSSN.EQ.1.AND.NSVN.EQ.1) THEN
               DO NMK=1,4
                  VV(NEQ,NEQ,NMK,NSW)=0.D0
                  DD(NEQ,NEQ,NMK,NSW)=0.D0
               ENDDO
               DO NEQ1=1,NEQMAX
                  NSSN1=NSS(NEQ1)
                  NSVN1=NSV(NEQ1)
                  IF(NSSN1.NE.1.AND.NSVN1.EQ.1) THEN
                     DO NMK=1,4
                        VV(NEQ,NEQ1,NMK,NSW) = VV(NEQ1,NEQ1,NMK,NSW)*PZ(NSSN1)
                        DD(NEQ,NEQ1,NMK,NSW) = DD(NEQ1,NEQ1,NMK,NSW)*PZ(NSSN1)
                     ENDDO
                  ENDIF
               ENDDO
               IND=IND+1
            ELSEIF(NSSN.EQ.1.AND.NSVN.EQ.2) THEN
               DO NMK=1,4
                  VV(NEQ,NEQ-1,NMK,NSW)=0.D0
                  DD(NEQ,NEQ-1,NMK,NSW)=0.D0
               ENDDO
               DO NEQ1=1,NEQMAX
                  NSSN1=NSS(NEQ1)
                  NSVN1=NSV(NEQ1)
                  IF(NSSN1.NE.1.AND.NSVN1.EQ.2) THEN
                     DO NMK=1,4
                        VV(NEQ,NEQ1-1,NMK,NSW) = VV(NEQ1,NEQ1-1,NMK,NSW) *PZ(NSSN1)
                        DD(NEQ,NEQ1-1,NMK,NSW) = DD(NEQ1,NEQ1-1,NMK,NSW) *PZ(NSSN1)
                     ENDDO
                  ENDIF
               ENDDO
               IND=IND+1
            ENDIF
            IF(IND.EQ.2) GOTO 1000
         ENDDO
      ENDIF
 1000 CONTINUE

!     /* Interim Parameter */

      IF(NR.EQ.NRMAX) THEN
      DO NV=1,NEQMAX
      DO NW=1,NEQMAX
         NO=1
         VI(NV,NW,NO,NSW)=0.5D0*(VV(NV,NW,NO,NSW)+VV(NV,NW,NO+2,NSW))
         DI(NV,NW,NO,NSW)=0.5D0*(DD(NV,NW,NO,NSW)+DD(NV,NW,NO+2,NSW))
         NO=2
         VI(NV,NW,NO,NSW)=VV(NV,NW,NO,NSW)
         DI(NV,NW,NO,NSW)=DD(NV,NW,NO,NSW)
      ENDDO
      ENDDO
      ELSE
      DO NV=1,NEQMAX
      DO NW=1,NEQMAX
      DO NO=1,2
         VI(NV,NW,NO,NSW)=0.5D0*(VV(NV,NW,NO,NSW)+VV(NV,NW,NO+2,NSW))
         DI(NV,NW,NO,NSW)=0.5D0*(DD(NV,NW,NO,NSW)+DD(NV,NW,NO+2,NSW))
      ENDDO
      ENDDO
      ENDDO
      ENDIF

!     /* Coefficients of source term */

      DO NEQ=1,NEQMAX
         NSSN=NSS(NEQ)
         NSVN=NSV(NEQ)
         IF(NR.EQ.NRMAX.AND.MDDIAG.EQ.0) THEN
            IF(NSSN.EQ.0) THEN
               D(NEQ,NR) = 0.D0
            ELSEIF(NSSN.EQ.5) THEN
               D(NEQ,NR) = VOID
            ELSEIF(NSSN.EQ.6) THEN
               D(NEQ,NR) = VOID
            ELSEIF(NSSN.EQ.7.OR.NSSN.EQ.8) THEN
               D(NEQ,NR) = SSIN(NR,NSSN)*DV11 &
                    +(-VI(NEQ,NEQ,2,NSW)+C83*DI(NEQ,NEQ,2,NSW))*RNV(NR,NSSN)
            ELSE
               IF(NSVN.EQ.1) THEN
                  D(NEQ,NR) = (SSIN(NR,NSSN)+SPE(NR,NSSN)/DT)*DV11 &
                       +(-VI(NEQ,NEQ,2,NSW)+C83*DI(NEQ,NEQ,2,NSW))*RNV(NR,NSSN)
               ELSEIF(NSVN.EQ.2) THEN
                  D(NEQ,NR) = (PIN(NR,NSSN)/(RKEV*1.D20)  )*DV53 &
                       +(-VI(NEQ,NEQ-1,2,NSW)+C83*DI(NEQ,NEQ-1,2,NSW)) &
                       *RNV(NR,NSSN) &
                       +(-VI(NEQ,NEQ  ,2,NSW)+C83*DI(NEQ,NEQ  ,2,NSW)) &
                       *RPV(NR,NSSN)
               ELSEIF(NSVN.EQ.3) THEN
                  D(NEQ,NR) = VOID
               ELSE
                  STOP 'XX TR_COEF_DECIDE: must be NSSV=0 if NSSN=0'
               ENDIF
            ENDIF

!     *** OFF-DIAGONAL TRANSPORT ***
!     *
         ELSEIF(NR.EQ.NRMAX.AND.MDDIAG.NE.0) THEN
            IF(NSSN.EQ.0) THEN
               D(NEQ,NR) = 0.D0
            ELSEIF(NSSN.EQ.5) THEN
               D(NEQ,NR) = VOID
            ELSEIF(NSSN.EQ.6) THEN
               D(NEQ,NR) = VOID
            ELSEIF(NSSN.EQ.7.OR.NSSN.EQ.8) THEN
               D(NEQ,NR) = SSIN(NR,NSSN)*DV11
               DO NEQ1=1,NEQMAX
                  NSSN1=NSS(NEQ1)
                  NSVN1=NSV(NEQ1)
!                  IF((NSSN1.EQ.7.OR.NSSN1.EQ.8).AND.NSVN1.EQ.1) THEN
                  IF(NSVN1.EQ.1) THEN
                     D(NEQ,NR) = D(NEQ,NR)+(-VI(NEQ,NEQ1,2,NSW)+C83*DI(NEQ,NEQ1,2,NSW))*RNV(NR,NSSN1)
                  ENDIF
               ENDDO
            ELSE
               IF(NSVN.EQ.1) THEN
                  D(NEQ,NR) = (SSIN(NR,NSSN)+SPE(NR,NSSN)/DT)*DV11
                  DO NEQ1=1,NEQMAX
                     NSSN1=NSS(NEQ1)
                     NSVN1=NSV(NEQ1)
                     IF(NSVN1.EQ.1) THEN
                        D(NEQ,NR) = D(NEQ,NR)+(-VI(NEQ,NEQ1,2,NSW)+C83*DI(NEQ,NEQ1,2,NSW))*RNV(NR,NSSN1)
                     ELSEIF(NSVN1.EQ.2) THEN
                        D(NEQ,NR) = D(NEQ,NR)+(-VI(NEQ,NEQ1,2,NSW)+C83*DI(NEQ,NEQ1,2,NSW))*RPV(NR,NSSN1)
                     ENDIF
                  ENDDO
               ELSEIF(NSVN.EQ.2) THEN
                  D(NEQ,NR) = (PIN(NR,NSSN)/(RKEV*1.D20)  )*DV53
                  DO NEQ1=1,NEQMAX
                     NSSN1=NSS(NEQ1)
                     NSVN1=NSV(NEQ1)
                     IF(NSVN1.EQ.1) THEN
                        D(NEQ,NR) = D(NEQ,NR)+(-VI(NEQ,NEQ1,2,NSW)+C83*DI(NEQ,NEQ1,2,NSW))*RNV(NR,NSSN1)
                     ELSEIF(NSVN1.EQ.2) THEN
                        D(NEQ,NR) = D(NEQ,NR)+(-VI(NEQ,NEQ1,2,NSW)+C83*DI(NEQ,NEQ1,2,NSW))*RPV(NR,NSSN1)
                     ENDIF
                  ENDDO
               ELSEIF(NSVN.EQ.3) THEN
                  D(NEQ,NR) = VOID
               ELSE
                  STOP 'XX TR_COEF_DECIDE: must be NSSV=0 if NSSN=0'
               ENDIF
            ENDIF
!     *
!     *************

         ELSE
            IF(NSSN.EQ.0) THEN
               D(NEQ,NR) = ETA(NR  )*BB/(TTRHO(NR  )*ARRHO(NR  )*DR)*(AJ(NR  )-AJOH(NR  )) &
     &                    -ETA(NR+1)*BB/(TTRHO(NR+1)*ARRHO(NR+1)*DR)*(AJ(NR+1)-AJOH(NR+1))
            ELSEIF(NSSN.EQ.5) THEN
               D(NEQ,NR) = VOID
            ELSEIF(NSSN.EQ.6) THEN
               D(NEQ,NR) = VOID
            ELSEIF(NSSN.EQ.7.OR.NSSN.EQ.8) THEN
               D(NEQ,NR) = SSIN(NR,NSSN)*DV11
            ELSE
               IF(NSVN.EQ.1) THEN
                  D(NEQ,NR) = (SSIN(NR,NSSN)+SPE(NR,NSSN)/DT)*DV11
               ELSEIF(NSVN.EQ.2) THEN
                  D(NEQ,NR) = (PIN(NR,NSSN)/(RKEV*1.D20)    )*DV53
               ELSEIF(NSVN.EQ.3) THEN
                  D(NEQ,NR) = VOID
               ELSE
                  STOP 'XX TR_COEF_DECIDE: must be NSSV=0 if NSSN=0'
               ENDIF
            ENDIF
         ENDIF
      ENDDO

      IF(MDLEQE.EQ.1) THEN

      IF(NR.EQ.NRMAX) THEN
         IND=0
         VISUMN=0.D0
         DISUMN=0.D0
         VISUMT1=0.D0
         DISUMT1=0.D0
         VISUMT2=0.D0
         DISUMT2=0.D0
         DO NEQ=1,NEQMAX
            NSSN=NSS(NEQ)
            NSVN=NSV(NEQ)
            IF(NSSN.EQ.1.AND.NSVN.EQ.1) THEN
               D(NEQ,NR) = 0.D0
               DO NEQ1=1,NEQMAX
                  NSSN1=NSS(NEQ1)
                  NSVN1=NSV(NEQ1)
                  IF(NSSN1.NE.1.AND.NSVN1.EQ.1) THEN
                     D(NEQ,NR) = D(NEQ,NR)+PZ(NSSN1)*(SSIN(NR,NSSN1)+SPE(NR,NSSN1)/DT)*DV11
                     VISUMN = VISUMN+PZ(NSSN1)*VI(NEQ1,NEQ1,2,NSW)*RNV(NR,NSSN1)
                     DISUMN = DISUMN+PZ(NSSN1)*DI(NEQ1,NEQ1,2,NSW)*RNV(NR,NSSN1)
                  ENDIF
               ENDDO
               D(NEQ,NR) = D(NEQ,NR)+(-VISUMN+C83*DISUMN)
               IND=IND+1
            ELSEIF(NSSN.EQ.1.AND.NSVN.EQ.2) THEN
               D(NEQ,NR) = 0.D0
               DO NEQ1=1,NEQMAX
                  NSSN1=NSS(NEQ1)
                  NSVN1=NSV(NEQ1)
                  IF(NSSN1.NE.1.AND.NSVN1.EQ.2) THEN
                     D(NEQ,NR) = D(NEQ,NR)+PZ(NSSN1)*(PIN(NR,NSSN1)/(RKEV*1.D20)  )*DV53
                     VISUMT1 = VISUMT1+PZ(NSSN1)*VI(NEQ1,NEQ1-1,2,NSW)*RNV(NR,NSSN1)
                     DISUMT1 = DISUMT1+PZ(NSSN1)*DI(NEQ1,NEQ1-1,2,NSW)*RNV(NR,NSSN1)
                     VISUMT2 = VISUMT2+PZ(NSSN1)*VI(NEQ1,NEQ1  ,2,NSW)*RPV(NR,NSSN1)
                     DISUMT2 = DISUMT2+PZ(NSSN1)*DI(NEQ1,NEQ1  ,2,NSW)*RPV(NR,NSSN1)
                  ENDIF
               ENDDO
               D(NEQ,NR) = D(NEQ,NR)+(-VISUMT1+C83*DISUMT1)+(-VISUMT2+C83*DISUMT2)
               IND=IND+1
            ENDIF
           IF(IND.EQ.2) GOTO 2000
         ENDDO
      ELSE
         DO NEQ=1,NEQMAX
            NSSN=NSS(NEQ)
            NSVN=NSV(NEQ)
            IF(NSSN.EQ.1.AND.NSVN.EQ.1) THEN
               D(NEQ,NR) = 0.D0
               DO NEQ1=1,NEQMAX
                  NSSN1=NSS(NEQ1)
                  NSVN1=NSV(NEQ1)
                  IF(NSSN1.NE.1.AND.NSVN1.EQ.1) THEN
                     D(NEQ,NR)=D(NEQ,NR)+PZ(NSSN1)*D(NEQ1,NR)
                  ENDIF
               ENDDO
               GOTO 2000
            ENDIF
         ENDDO
      ENDIF

      ENDIF
 2000 CONTINUE

      RETURN
      END SUBROUTINE TR_COEF_DECIDE

!     ***********************************************************

!           IONIZATION

!     ***********************************************************

      SUBROUTINE TR_IONIZATION(NR)

      USE TRCOMM, ONLY : DVRHO, MDLEQ0, NEA, NSMAX, PN, TSIE
      USE TRCOM1, ONLY : B
      IMPLICIT NONE
!      INCLUDE 'trcomm.inc'
      INTEGER,INTENT(IN):: NR
      INTEGER:: NV, NS, NW, NEQ

      IF(MDLEQ0.EQ.1) THEN
!     *** electron density ***
         NV=NEA(1,1)
         DO NS=7,8
            NW=NEA(NS,1)
            B(NV,NW,NR)=B(NV,NW,NR)+TSIE(NR)*DVRHO(NR)
         ENDDO
         IF(NSMAX.EQ.2) THEN
!     *** deuterium density ***
            NV=NEA(2,1)
            DO NS=7,8
               NW=NEA(NS,1)
               B(NV,NW,NR)=B(NV,NW,NR)+TSIE(NR)*DVRHO(NR)
            ENDDO
         ELSE
!     *** deuterium density ***
            NV=NEA(2,1)
            DO NS=7,8
               NW=NEA(NS,1)
               B(NV,NW,NR)=B(NV,NW,NR)+(PN(2)/(PN(2)+PN(3)))*TSIE(NR)*DVRHO(NR)
            ENDDO
!     *** tritium density ***
            NV=NEA(3,1)
            DO NS=7,8
               NW=NEA(NS,1)
               B(NV,NW,NR)=B(NV,NW,NR)+(PN(3)/(PN(2)+PN(3)))*TSIE(NR)*DVRHO(NR)
            ENDDO
         ENDIF
!     *** neutrals ***
         DO NS=7,8
            NEQ=NEA(NS,1)
            B(NEQ,NEQ,NR)=B(NEQ,NEQ,NR)-TSIE(NR)*DVRHO(NR)
         ENDDO
      ENDIF

      RETURN
      END SUBROUTINE TR_IONIZATION

!     ***********************************************************

!           CHARGE EXCHANGE

!     ***********************************************************

      SUBROUTINE TR_CHARGE_EXCHANGE(NR)

      USE TRCOMM, ONLY : DVRHO, MDLEQ0, NEA, TSCX
      USE TRCOM1, ONLY : B
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: NR
      INTEGER:: NEQ, NV, NW


      IF(MDLEQ0.EQ.1) THEN
         NEQ=NEA(7,1)
         B(NEQ,NEQ,NR)=B(NEQ,NEQ,NR)-TSCX(NR)*DVRHO(NR)
         NV=NEA(8,1)
         NW=NEA(7,1)
         B(NV,NW,NR)=B(NV,NW,NR)+TSCX(NR)*DVRHO(NR)
      ENDIF

      RETURN
      END SUBROUTINE TR_CHARGE_EXCHANGE

!     ***************************************************************

!           LOCAL VARIABLES FOR N, T AND P ON GRID

!     **************************************************************

      FUNCTION RNV(NR,NS)

      USE TRCOMM, ONLY : NRMAX, PNSS, PNSSA, RHOA, RN, rkind
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NR,NS
      REAL(rkind):: RNV

      IF(NR.EQ.NRMAX) THEN
         IF(RHOA.EQ.1.D0) THEN
            RNV=PNSS(NS)
         ELSE
            RNV=PNSSA(NS)
         ENDIF
      ELSEIF(NR.NE.0.AND.NR.NE.NRMAX) THEN
         RNV=0.5D0*(RN(NR,NS)+RN(NR+1,NS))
      ENDIF

      RETURN
      END FUNCTION RNV

      FUNCTION RTV(NR,NS)

      USE TRCOMM, ONLY : NRMAX, PTS, PTSA, RHOA, RT, rkind
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NR,NS
      REAL(rkind):: RTV

      IF(NR.EQ.NRMAX) THEN
         IF(RHOA.EQ.1.D0) THEN
            RTV=PTS(NS)
         ELSE
            RTV=PTSA(NS)
         ENDIF
      ELSEIF(NR.NE.0.AND.NR.NE.NRMAX) THEN
         RTV=0.5D0*(RT(NR,NS)+RT(NR+1,NS))
      ENDIF

      RETURN
      END FUNCTION RTV

      FUNCTION RPV(NR,NS)

      USE TRCOMM, ONLY : NRMAX, PNSS, PNSSA, PTS, PTSA, RHOA, RN, RT, rkind
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NR, NS
      REAL(rkind):: RPV


      IF(NR.EQ.NRMAX) THEN
         IF(RHOA.EQ.1.D0) THEN
            RPV=PNSS(NS)*PTS(NS)
         ELSE
            RPV=PNSSA(NS)*PTSA(NS)
         ENDIF
      ELSEIF(NR.NE.0.AND.NR.NE.NRMAX) THEN
         RPV=0.5D0*(RN(NR,NS)*RT(NR,NS)+RN(NR+1,NS)*RT(NR+1,NS))
      ENDIF

!      IF(NR.EQ.NRMAX.AND.NS.EQ.1) WRITE(6,'(A,1PE12.4)') 'RPV=',RPV

      RETURN
      END FUNCTION RPV

END MODULE trexec
