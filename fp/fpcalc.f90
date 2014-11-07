!     $Id$
!
! ******************************************************************
!
!     CALCULATION OF COLLISINAL DIFUSION COEFFICIENTS: DC AND FC
!
! ******************************************************************
!

      MODULE fpcalc

      USE fpcomm
      USE fpcalcn
      USE fpcalcnr
      USE fpmpi
      USE libmpi
      USE libspf, ONLY: ERF0,ERF1
      real(8):: PMAXC
      integer:: NSB_ISO
      real(8):: RTFD0L, RTFDL
      integer,parameter:: ISW_NOTAIL=0

      contains

!----------------------------------

      SUBROUTINE FP_CALC(NSA)

      USE libmtx
      IMPLICIT NONE
      integer:: NSA,NSB,NR, NP, NTH
      integer:: nsrc,nsend
      real(8):: RGAMH, RGAMH2, RZI, RTE, PFPL, VFPL, U, DCTTL, RGAMA, DFDP, DFDTH

      DO NR=NRSTART,NREND
         DO NP=NPSTART,NPENDWG
         DO NTH=1,NTHMAX
            DCPP(NTH,NP,NR,NSA)=0.D0
            DCPT(NTH,NP,NR,NSA)=0.D0
            FCPP(NTH,NP,NR,NSA)=0.D0
         ENDDO
         ENDDO
         DO NP=NPSTARTW,NPENDWM
         DO NTH=1,NTHMAX+1
            DCTP(NTH,NP,NR,NSA)=0.D0
            DCTT(NTH,NP,NR,NSA)=0.D0
            FCTH(NTH,NP,NR,NSA)=0.D0
         ENDDO
         ENDDO
      ENDDO

      DO NSB=1,NSBMAX
      DO NR=NRSTART,NREND
         DO NP=NPSTART,NPENDWG
         DO NTH=1,NTHMAX
            DCPP2(NTH,NP,NR,NSB,NSA)=0.D0
            DCPT2(NTH,NP,NR,NSB,NSA)=0.D0
            FCPP2(NTH,NP,NR,NSB,NSA)=0.D0
         ENDDO
         ENDDO
         DO NP=NPSTARTW,NPENDWM
         DO NTH=1,NTHMAX+1
            DCTP2(NTH,NP,NR,NSB,NSA)=0.D0
            DCTT2(NTH,NP,NR,NSB,NSA)=0.D0
            FCTH2(NTH,NP,NR,NSB,NSA)=0.D0
         ENDDO
         ENDDO
      ENDDO
      ENDDO

!      DO NR=NRSTART,NREND
!         DO NP=1, NPMAX+1
!         DO NTH=1,4
!            DCPPB(NTH,NP,NR,NSA)=0.D0
!            DCPTB(NTH,NP,NR,NSA)=0.D0
!            FCPPB(NTH,NP,NR,NSA)=0.D0
!         END DO
!         END DO
!      END DO

!      DO NSB=1,NSBMAX
!      DO NR=NRSTART,NREND
!         DO NP=1, NPMAX+1
!         DO NTH=1,4
!            DCPP2B(NTH,NP,NR,NSB,NSA)=0.D0
!            DCPT2B(NTH,NP,NR,NSB,NSA)=0.D0
!            FCPP2B(NTH,NP,NR,NSB,NSA)=0.D0
!         END DO
!         END DO
!      END DO
!      END DO

      DO NR=NRSTART,NREND
         DO NSB=1,NSBMAX
!            if(nr.eq.1) write(6,'(A,I8,1P2E12.4)') 
!     &           ' NSB,RN,RNFD=',NSB,RN(NSB),RNFD(NR,NSB)
!
            IF(MODELC.EQ.0.or.MODELC.eq.1) THEN
               CALL FPCALC_L(NR,NSB,NSA)
            ELSEIF(MODELC.EQ.2.or.MODELC.eq.3) THEN
               IF(NS_NSB(NSB).EQ.NS_NSA(NSA)) THEN
                  CALL FPCALC_NL(NR,NSB,NSA)
               ELSE
                  CALL FPCALC_L(NR,NSB,NSA)
               ENDIF
            ELSEIF(MODELC.EQ.4) THEN
               IF(MODELR.eq.0)THEN
                  CALL FPCALC_NL(NR,NSB,NSA)
               ELSE IF(MODELR.eq.1)THEN
                  CALL FPCALC_NLR(NR,NSB,NSA)
               END IF
            ELSEIF(MODELC.EQ.5) THEN
               IF(NS_NSB(NSB).NE.NS_NSA(NSA)) THEN
                  CALL FPCALC_L(NR,NSB,NSA)
               ENDIF
            ELSEIF(MODELC.EQ.6) THEN
               IF(NS_NSB(NSB).NE.NS_NSA(NSA)) THEN
                  CALL FPCALC_NL(NR,NSB,NSA)
               ENDIF
            ELSEIF(MODELC.EQ.7) THEN
               IF(MODELR.eq.0)THEN
                  IF(NS_NSB(NSB).EQ.NS_NSA(NSA)) THEN
                     CALL FPCALC_NL(NR,NSB,NSA)
                  ELSE
                     CALL FPCALC_L(NR,NSB,NSA)
                  ENDIF
               ELSE IF(MODELR.eq.1)THEN
                  IF(NS_NSB(NSB).EQ.NS_NSA(NSA)) THEN
                     CALL FPCALC_NLR(NR,NSB,NSA)
                  ELSE
                     CALL FPCALC_L(NR,NSB,NSA)
                  ENDIF
               END IF
! For conductivity check. Karney adopt simple calculation for e-i
            ELSEIF(MODELC.EQ.-1) THEN
               IF(NS_NSB(NSB).EQ.NS_NSA(NSA)) THEN
                  MODELC=0
                  CALL FPCALC_L(NR,NSB,NSA)
                  MODELC=-1
               ENDIF
            ELSEIF(MODELC.EQ.-2) THEN 
               IF(NS_NSB(NSB).EQ.NS_NSA(NSA)) THEN
                  MODELC=1
                  CALL FPCALC_L(NR,NSB,NSA)
                  MODELC=-2
               ENDIF
            ELSEIF(MODELC.EQ.-4) THEN
               IF(NS_NSB(NSB).EQ.NS_NSA(NSA)) THEN
                  IF(MODELR.eq.0)THEN
                     CALL FPCALC_NL(NR,NSB,NSA)
                  ELSE IF(MODELR.eq.1)THEN
                     CALL FPCALC_NLR(NR,NSB,NSA)
                  END IF
               ELSEIF(NS_NSB(NSB).eq.1.and.NS_NSA(NSA).eq.2)THEN
                  IF(MODELR.eq.0)THEN
                     CALL FPCALC_NL(NR,NSB,NSA)
                  ELSE IF(MODELR.eq.1)THEN
                     CALL FPCALC_NLR(NR,NSB,NSA)
                  END IF
               ENDIF
            ENDIF
         ENDDO
!        ----- Simple electron-ion collision term using ZEFF ----- neglect i-e
         IF(MODELC.LT.0) THEN
            IF(NS_NSA(NSA).EQ.1) THEN
               DO NSB=1,NSBMAX
                  IF(NS_NSB(NSB).EQ.2) THEN
                     RGAMH=RNUD(NR,NSA,NSA)*SQRT(2.D0)*VTFD(NR,NSA)*AMFP(NSA) &
                          /(RNFP0(NSA)*PTFP0(NSA)*1.D20)*RNFD0(NSA)
                     RGAMH2=RGAMH*RNFP(NR,NSA)*1.D20*PTFP0(NSA)/AMFP(NSA)/RNFD0(NSA)
                     rZI = -PZ(2)/PZ(1) &
!                     rZI = ( PZ(2)/PZ(1) )**2 * RNFD(NR,NSB)/RNFP(NR,NSA) &
                          /(14.9D0-0.5D0*LOG(RNFP(NR,NSA))+LOG(RTFP(NR,NSA)))* &
                          (15.2D0-0.5D0*LOG(RNFP(NR,NSA))+LOG(RTFP(NR,NSA)))
                     DO NP=NPSTARTW,NPENDWM
                        RGAMA =SQRT(1.D0+PM(NP,NSA)**2*THETA0(NSA))
                        PFPL=PM(NP,NSA)*PTFP0(NSA)
                        VFPL=PFPL/AMFP(NSA)/RGAMA
                        DCTTL=0.5D0*RZI*RGAMH2/VFPL
                        DO NTH=1,NTHMAX+1
                           DCTT2(NTH,NP,NR,NSB,NSA) &
                                =DCTT2(NTH,NP,NR,NSB,NSA)+DCTTL
                        ENDDO
                     ENDDO
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

!     ----- bounce average -----
         IF(MODELA.EQ.1) THEN
            IF(MODELC.EQ.0.or.MODELC.eq.1) THEN
               CALL FPCALC_LAV(NR,NSA)
            ELSEIF(MODELC.EQ.2.or.MODELC.eq.3) THEN
               CALL FPCALC_NLAV(NR,NSA)
            ELSEIF(MODELC.EQ.4) THEN
               CALL FPCALC_NLAV(NR,NSA)
            ELSEIF(MODELC.EQ.5) THEN
               CALL FPCALC_LAV(NR,NSA)
            ELSEIF(MODELC.EQ.6) THEN
               CALL FPCALC_NLAV(NR,NSA)
            ELSEIF(MODELC.EQ.-1) THEN
               CALL FPCALC_LAV(NR,NSA)
            ELSEIF(MODELC.EQ.-2) THEN
               CALL FPCALC_NLAV(NR,NSA)
            ENDIF
         ENDIF
!     sum up coefficients by species
         DO NSB=1,NSBMAX
            DO NP=NPSTART,NPENDWG
               DO NTH=1,NTHMAX
                  DCPP(NTH,NP,NR,NSA)=DCPP(NTH,NP,NR,NSA) &
                                     +DCPP2(NTH,NP,NR,NSB,NSA)
                  DCPT(NTH,NP,NR,NSA)=DCPT(NTH,NP,NR,NSA) &
                                     +DCPT2(NTH,NP,NR,NSB,NSA)
                  FCPP(NTH,NP,NR,NSA)=FCPP(NTH,NP,NR,NSA) &
                                     +FCPP2(NTH,NP,NR,NSB,NSA)
               END DO
!               DO NTH=1,4
!                  DCPPB(NTH,NP,NR,NSA)=DCPPB(NTH,NP,NR,NSA) &
!                                     +DCPP2B(NTH,NP,NR,NSB,NSA)
!                  DCPTB(NTH,NP,NR,NSA)=DCPTB(NTH,NP,NR,NSA) &
!                                     +DCPT2B(NTH,NP,NR,NSB,NSA)
!                  FCPPB(NTH,NP,NR,NSA)=FCPPB(NTH,NP,NR,NSA) &
!                                     +FCPP2B(NTH,NP,NR,NSB,NSA)
!               END DO
            END DO
!            DO NP=1,NPMAX
            DO NP=NPSTARTW,NPENDWM
               DO NTH=1,NTHMAX+1
                  DCTP(NTH,NP,NR,NSA)=DCTP(NTH,NP,NR,NSA) &
                                     +DCTP2(NTH,NP,NR,NSB,NSA)
                  DCTT(NTH,NP,NR,NSA)=DCTT(NTH,NP,NR,NSA) &
                                     +DCTT2(NTH,NP,NR,NSB,NSA)
                  FCTH(NTH,NP,NR,NSA)=FCTH(NTH,NP,NR,NSA) &
                                     +FCTH2(NTH,NP,NR,NSB,NSA)
               END DO
            END DO
         END DO
      ENDDO

!      IF(NSA.eq.1)THEN
!         open(8,file='T_c4.dat')
!         DO NP=2,NPMAX
!            WRITE(8,'(I3,2E18.10)') NP, PG(NP,1), PG(NP,1)*DCPP2(1,NP,1,1,1)/FCPP2(1,NP,1,1,1)
!         END DO
!         close(8)
!      ENDIF

!      IF(NRANK.eq.0)THEN
!      open(8,file='DCe-i.dat')
!      DO NP=1,NPMAX+1
!            WRITE(8,'(4E16.8)') &
!                 PG(NP,NSA), DCTT2(1,NP,1,2,1), DCPP2(1,NP,1,2,1), FCPP2(1,NP,1,2,1)
!      END DO
!      close(8)
!      END IF

!      open(9,file='flux_th.dat')
!      DO NP=1,NPMAX
!!         DO NTH=1,NTHMAX+1
!            WRITE(9,'(4E16.8)') &
!!PM(NP,NSA)*COSG(NTH),PM(NP,NSA)*SING(NTH),
!                 PM(NP,NSA), DCTP2(NTH,NP,NR,1,1),DCTT2(NTH,NP,NR,1,1),FCTH2(NTH,NP,NR,1,1)
!!         END DO
!!         WRITE(9,*) " "
!!         WRITE(9,*) " "
!      END DO
!      close(9)
!      END IF

!      IF(MODELR.eq.1)THEN
!         DO NP=2,NPMAX+1
!            RGAMA =SQRT(1.D0+PG(NP,1)**2*THETA0(1))
!            WRITE(*,'(I3,3E18.10)') NP, PG(NP,1),FCPP2(NTH,NP,NR,1,1)/DCPP2(NTH,NP,NR,1,1)*RGAMA/PG(NP,1), -THETA0(1)/THETA(NR,NSA)
!         END DO
!      ELSE
!         DO NP=2,NPMAX+1
!            WRITE(*,'(I3,3E18.10)') NP, PG(NP,1),FCPP2(NTH,NP,NR,1,1)/DCPP2(NTH,NP,NR,1,1)/PG(NP,1), -RTFP0(1)/RTFD(NR,1)
!         END DO
!      END IF
!      END IF

      IF(nrank.eq.0) THEN
      IF(MOD(IDBGFP/8,2).EQ.1) THEN
!
!     +++ plot of D_coll +++
!
         CALL PAGES
         CALL FPGRFA(1,DCPP,PG,1,'@DCPP(P)@',NTHMAX+1,NPMAX+1,NRMAX+1, &
                                             NTHMAX,NPMAX,NRMAX,NSA)
         CALL FPGRFA(2,DCTT,PM,2,'@DCTT(P)@',NTHMAX+1,NPMAX+1,NRMAX+1, &
                                             NTHMAX,NPMAX,NRMAX,NSA)
         CALL FPGRFA(3,FCPP,PG,1,'@FCPP(P)@',NTHMAX+1,NPMAX+1,NRMAX+1, &
                                             NTHMAX,NPMAX,NRMAX,NSA)
         CALL PAGEE
      ENDIF

      IF(MOD(IDBGFP/16,2).EQ.1) THEN
!
!     +++ plot of D_coll +++
!
         CALL PAGES
         CALL FPGRFB(1,DCPP,THM,1,'@DCPP(TH)@',NTHMAX+1,NPMAX+1,NRMAX+1, &
                                               NTHMAX,NPMAX,NRMAX,NSA)
         CALL FPGRFB(2,DCTT,THG,2,'@DCTT(TH)@',NTHMAX+1,NPMAX+1,NRMAX+1, &
                                               NTHMAX,NPMAX,NRMAX,NSA)
         CALL FPGRFB(3,FCPP,THM,1,'@FCPP(TH)@',NTHMAX+1,NPMAX+1,NRMAX+1, &
                                               NTHMAX,NPMAX,NRMAX,NSA)
         CALL PAGEE
      ENDIF
      ENDIF

      RETURN
      END SUBROUTINE FP_CALC

!--------------------------------------------------------------
      SUBROUTINE FPGRFA(ID,DATA,P,IND,TITLE,NTHM,NPM,NRM, &
                                            NTHMAX,NPMAX,NRMAX,NSA)

      USE libgrf,ONLY: grd1d
      implicit none
      integer, intent(IN)::  ID,IND,NTHM,NPM,NRM,NTHMAX,NPMAX,NRMAX,NSA
      real(8), dimension(NTHM,NPM,NRM,NSA), intent(IN):: DATA
      real(8), dimension(NPM), intent(IN):: P
      character(len=*):: title
      real(8), dimension(NPM,9):: WORK
      integer, dimension(9):: NTHG
      integer:: nth,np,nr,ng,npmaxg

      if(NRMAX.GT.1) then
    1    WRITE(6,*) '## NR ?'
         READ(5,*,err=1,end=9999)
         if(NR.LT.1.OR.NR.GT.NRMAX) THEN
            write(6,*) 'XX NR must be between 1 and NRMAX:',NRMAX
            goto 1
         endif
      else
         nr=1
      endif

      nthg(1)=1
      nthg(2)=NTHMAX/8
      nthg(3)=NTHMAX/4
      if(ind.eq.0) then
         nthg(4)=NTHMAX/2-NTHMAX/8
         nthg(5)=NTHMAX/2
         nthg(6)=NTHMAX/2+NTHMAX/8+1
         nthg(7)=NTHMAX-NTHMAX/4
         nthg(8)=NTHMAX-NTHMAX/8
         nthg(9)=NTHMAX
         npmaxg=npmax
      elseif(ind.eq.1) then
         nthg(4)=NTHMAX/2-NTHMAX/8
         nthg(5)=NTHMAX/2
         nthg(6)=NTHMAX/2+NTHMAX/8+1
         nthg(7)=NTHMAX-NTHMAX/4
         nthg(8)=NTHMAX-NTHMAX/8
         nthg(9)=NTHMAX
         npmaxg=npmax+1
      elseif(ind.eq.2) then
         nthg(4)=NTHMAX/2-NTHMAX/8+1
         nthg(5)=NTHMAX/2+1
         nthg(6)=NTHMAX/2+NTHMAX/8+1
         nthg(7)=NTHMAX+1-NTHMAX/4
         nthg(8)=NTHMAX+1-NTHMAX/8
         nthg(9)=NTHMAX+1
         npmaxg=npmax
      endif

      do ng=1,9
         nth=NTHG(ng)
         do np=1,npmaxg
            work(np,ng)=data(nth,np,nr,NSA)
         enddo
      enddo
            
      CALL GRD1D(ID,p,work,npm,npmaxg,9,TITLE,0)

 9999 return
      end SUBROUTINE FPGRFA
!--------------------------------------
      SUBROUTINE FPGRFB(ID,DATA,TH,IND,TITLE,NTHM,NPM,NRM, &
                                            NTHMAX,NPMAX,NRMAX,NSA)
!
      USE libgrf,ONLY: grd1d
      implicit none
      integer, intent(IN)::  ID,IND,NTHM,NPM,NRM,NTHMAX,NPMAX,NRMAX,NSA
      real(8), dimension(NTHM,NPM,NRM,NSA), intent(IN):: DATA
      real(8), dimension(NTHM), intent(IN):: TH
      character(len=*):: title
      real(8), dimension(NTHM,10):: WORK
      integer, dimension(10):: NPG
      integer:: nth,np,nr,ng,nthmaxg

      if(NRMAX.GT.1) then
    1    WRITE(6,*) '## NR ?'
         READ(5,*,err=1,end=9999)
         if(NR.LT.1.OR.NR.GT.NRMAX) THEN
            write(6,*) 'XX NR must be between 1 and NRMAX:',NRMAX
            goto 1
         endif
      else
         nr=1
      endif

      do ng=1,10
         npg(ng)=NINT(0.1*NPMAX*ng)
      enddo
      if(ind.eq.0) then
         nthmaxg=nthmax
      elseif(ind.eq.1) then
         nthmaxg=nthmax
      elseif(ind.eq.2) then
         nthmaxg=nthmax+1
      endif

      do ng=1,10
         np=npg(ng)
         do nth=1,nthmaxg
            work(nth,ng)=data(nth,np,nr,NSA)
         enddo
      enddo
            
      CALL GRD1D(ID,th,work,nthm,nthmaxg,10,TITLE,0)

 9999 return
      end SUBROUTINE FPGRFB

!
! ***************************************************************
!
!                       SET OF INTEGRAND
!
! ***************************************************************
!
      FUNCTION FPFN0R(X)
!
      IMPLICIT NONE
      real(8)::FPFN0R
      real(8),INTENT(IN):: X
      real(8)::PN
!
      PN=X
      FPFN0R=PN**2*FPRMXW(PN)
!
      RETURN
      END FUNCTION FPFN0R
!---------------------------------------------------------------
      FUNCTION FPFN1R(X,XM,XP)

      IMPLICIT NONE
      real(8)::FPFN1R
      real(8),INTENT(IN)::X, XM, XP
      real(8)::PN, A, B

      A=0.5D0*PNFP
      PN=A*XP
      B=PN**4/(1.D0+PN**2*TMC2FD0)
      FPFN1R=A*B*FPRMXW(PN)

      RETURN
      END FUNCTION FPFN1R
!
! ===============================================================
!
      FUNCTION FPFN2R(X)

      real(8):: FPFN2R
      real(8),INTENT(IN):: X
      real(8):: A, PN, B

      A=1.D0
      PN=A*(X+PNFP)
      B=PN*SQRT(1.D0+PN**2*TMC2FD0)
      FPFN2R=A*B*FPRMXW(PN)

      RETURN
      END FUNCTION FPFN2R
!
! ==============================================================
!
      FUNCTION FPFN3R(X,XM,XP)

      real(8):: FPFN3R
      real(8),INTENT(IN):: X, XM, XP
      real(8):: A, PN

      A=0.5D0*PNFP
      PN=A*XP
      FPFN3R=A*PN**2*FPRMXW(PN)

      RETURN
      END FUNCTION FPFN3R
!
! ===============================================================
!
      FUNCTION FPFN4R(X,XM,XP)

      real(8):: FPFN4R
      real(8),INTENT(IN):: X, XM, XP
      real(8):: A, PN, B

      A=0.5D0*PNFP
      PN=A*XP
      B=PN**2/SQRT(1.D0+PN**2*TMC2FD0)
      FPFN4R=A*B*FPRMXW(PN)

      RETURN
      END FUNCTION FPFN4R
!
! ===============================================================
!
      FUNCTION FPFN5R(X,XM,XP)

      real(8):: FPFN5R
      real(8),INTENT(IN):: X, XM, XP
      real(8):: A, PN, B

      A=0.5D0*PNFP
      PN=A*XP
      B=PN**4/(SQRT(1.D0+PN**2*TMC2FD0))**3
      FPFN5R=A*B*FPRMXW(PN)

      RETURN
      END FUNCTION FPFN5R
!
! ===============================================================
!
      FUNCTION FPFN6R(X)

      real(8):: FPFN6R
      real(8),INTENT(IN):: X
      real(8):: A, PN

      A=1.D0
      PN=A*(X+PNFP)
      FPFN6R=A*PN*FPRMXW(PN)

      RETURN
      END FUNCTION FPFN6R
!
! =============================================================== 
      FUNCTION FPFN7R(X,XM,XP)
!                            
      real(8):: FPFN7R
      real(8),INTENT(IN):: X, XM, XP
      real(8):: A, PN, B, PMAX2

      PMAX2=PMAXC
      A=0.5D0*(PNFP-PMAX2)
      PN=A*XP+0.5D0*(PNFP+PMAX2)
      B=PN**4/(1.D0+PN**2*TMC2FD0)
      FPFN7R=A*B*FPRMXW(PN)

      RETURN
      END FUNCTION FPFN7R
!
! =============================================================== 
!
      FUNCTION FPFN8R(X,XM,XP)

      real(8):: FPFN8R
      real(8),INTENT(IN):: X, XM, XP
      real(8):: A, PN, B, PMAX2

      PMAX2=PMAXC
      A=0.5D0*(PNFP-PMAX2)
      PN=A*XP+0.5D0*(PNFP+PMAX2)
      FPFN8R=A*PN**2*FPRMXW(PN)

      RETURN
      END FUNCTION FPFN8R
!
! =============================================================== 
!
      FUNCTION FPFN9R(X,XM,XP)

      real(8):: FPFN9R
      real(8),INTENT(IN):: X, XM, XP
      real(8):: A, PN, B, PMAX2

      PMAX2=PMAXC
      A=0.5D0*(PNFP-PMAX2)
      PN=A*XP+0.5D0*(PNFP+PMAX2)
      B=PN**2/SQRT(1.D0+PN**2*TMC2FD0)
      FPFN9R=A*B*FPRMXW(PN)

      RETURN
      END FUNCTION FPFN9R
!
! =============================================================== 
!
      FUNCTION FPFN10R(X,XM,XP)

      real(8):: FPFN10R
      real(8),INTENT(IN):: X, XM, XP
      real(8):: A, PN, B, PMAX2

      PMAX2=PMAXC
      A=0.5D0*(PNFP-PMAX2)
      PN=A*XP+0.5D0*(PNFP+PMAX2)
      B=PN**4/(SQRT(1.D0+PN**2*TMC2FD0))**3
      FPFN10R=A*B*FPRMXW(PN)

      RETURN
      END FUNCTION FPFN10R
!
! ===============================================================
!
      FUNCTION FPRMXW(PN)

      USE plprof
      real(8):: FPRMXW
      real(8),INTENT(IN):: PN
      real(8):: EX

      IF(MODELR.eq.1)THEN
         EX=(1.D0-SQRT(1.D0+PN**2*TMC2FD0))/TMC2FD
!         IF (EX.LT.-100.D0)THEN
!            FPRMXW=0.D0
!         ELSE
            FPRMXW=EXP(EX)
!         ENDIF
      ELSEIF(MODELR.eq.0)THEN
         EX=-PN**2/(2.D0*RTFDL/RTFD0L)
!         IF (EX.LT.-100.D0)THEN
!            FPRMXW=0.D0
!         ELSE
            FPRMXW=EXP(EX)
!         ENDIF
      ENDIF

      RETURN
      END FUNCTION FPRMXW
!
! ************************************************************
!
!       CALCULATION OF LINEAR COLLISIONAL OPERATOR
!
! ************************************************************
!
      SUBROUTINE FPCALC_L(NR,NSB,NSA)

      USE libde
      IMPLICIT NONE

      integer:: NSA, NSB, NR, NP, NTH, NSBA
      real(8):: RGAMH, RGAMH2, RZI, RTE, PFPL, VFPL, U, DCTTL
      real(8):: RNNL, RNUFL, RNUDL, DCPPL, FCPPL, V
      real(8):: PNFPL, RGAMA, vtatb, ptatb, PCRIT
      real(8):: RINT0, ES0, RINT1, ES1, RINT2, ES2, RINT4, ES4, RINT5, ES5
      real(8):: RINT6, ES6, RINT7, ES7, RINT8, ES8, RINT9, ES9
      real(8):: RINT3, ES3, p_thermal, v_thermal, pe_thermal, EX

!     ------ define --------
      RNNL=RNFD(NR,NSB)/RNFP0(NSA)
      RNUFL=RNUF(NR,NSB,NSA)*RNNL
      RNUDL=RNUD(NR,NSB,NSA)*RNNL

!      RGAMH=RNUD(NR,NSB,NSA)*SQRT(2.D0)*VTFD(NR,NSB)*AMFP(NSA) &
!              /(RNFP0(NSA)*PTFP0(NSA)*1.D20)
      RGAMH=AEFP(NSA)**2*AEFD(NSB)**2*LNLAM(NR,NSB,NSA)/(4.D0*PI*EPS0**2) &
           *AMFP(NSA)/PTFP0(NSA)**3 
!     RGAMH = \hat{\Gamma}/n_b
      NSBA=NSB_NSA(NSA)

!      WRITE(*,*) " "
!      WRITE(*,'(2I4,3E16.8)') NSA,NSB,RGAMH,RNFD(1,NSB),VTFD(1,NSB)
!
!     ----- Non-Relativistic -----
!
      IF(MODELR.EQ.0) THEN 
         IF(MODELC.eq.0)THEN ! maxwellian
            IF(MODEL_DISRUPT.eq.0)THEN
               RTFDL=RTFD(NR,NSB)
               RTFD0L=RTFD0(NSB)
               v_thermal=VTFD(NR,NSB)
            ELSEIF(MODEL_DISRUPT.ge.1)THEN
               RTFDL=RT_quench(NR) ! [keV]
               v_thermal=SQRT( RTFDL*1.D3*AEE/AMFD(NSB))
            END IF
         DO NP=NPSTART,NPENDWG
            IF(NP.EQ.1) THEN
               DCPPL=RGAMH*RNFD(NR,NSB)*1.D20*(2.D0/(3.D0*SQRT(PI))) &
                          *(VTFP0(NSA)/(SQRT(2.D0)*VTFD(NR,NSB))) 
               FCPPL=0.D0
            ELSE
               PFPL=PG(NP,NSBA)*PTFP0(NSA)
               VFPL=PFPL/AMFP(NSA)
               V=VFPL/VTFP0(NSA)
!               U=VFPL/(SQRT(2.D0)*VTFD(NR,NSB))
               U=VFPL/(SQRT(2.D0)*v_thermal)
               DCPPL= 0.5D0*RNUDL/U   *(ERF0(U)/U**2-ERF1(U)/U)
               FCPPL=-      RNUFL/U**2*(ERF0(U)-U*ERF1(U))
            ENDIF
            DO NTH=1,NTHMAX
               DCPP2(NTH,NP,NR,NSB,NSA)=DCPP2(NTH,NP,NR,NSB,NSA)+DCPPL
               FCPP2(NTH,NP,NR,NSB,NSA)=FCPP2(NTH,NP,NR,NSB,NSA)+FCPPL
            ENDDO
!            WRITE(*,'(3I4,3E16.8)') NP, NSA,NSB,PG(NP,NSBA), DCPPL, FCPPL
         ENDDO

!         DO NP=1,NPMAX
         DO NP=NPSTARTW,NPENDWM
            PFPL=PM(NP,NSA)*PTFP0(NSA)
            VFPL=PFPL/AMFP(NSA)
            V=VFPL/VTFP0(NSA)
!            U=VFPL/(SQRT(2.D0)*VTFD(NR,NSB))
            U=VFPL/(SQRT(2.D0)*v_thermal)
            DCTTL= 0.25D0*RNUDL/U &
                         *((2.D0-1.D0/U**2)*ERF0(U)+ERF1(U)/U)
!
            DO NTH=1,NTHMAX+1
               DCTT2(NTH,NP,NR,NSB,NSA)=DCTT2(NTH,NP,NR,NSB,NSA)+DCTTL
            ENDDO
         ENDDO

         ELSEIF(MODELC.eq.1)THEN ! isotoropic
            PMAXC=PMAX(NSBA)
            TMC2FD =0.D0
            TMC2FD0=0.D0
            RGAMA=1.D0
            NSB_ISO=NSB
            IF(MODEL_DISRUPT.eq.0)THEN
               RTFDL=RTFD(NR,NSB)
               RTFD0L=RTFD0(NSB)
            ELSEIF(MODEL_DISRUPT.ge.1)THEN
!               RTFDL=RTFD(NR,NSB)*1.D-1
               RTFDL=RT_quench(NR) ! [keV]
               RTFD0L=RTFD0(NSB)
            END IF
            DO NP=NPSTART,NPENDWG
               IF(NP.EQ.1) THEN
                  PNFPL=PG(NP,NSBA)
                  PNFP=PNFPL
                  PCRIT=0.D0
                  CALL DEHIFT(RINT0,ES0,H0DE,EPSDE,0,FPFN0R,"N0R_NP1")
                  PNFP=PCRIT
                  CALL DEHIFT(RINT2,ES2,H0DE,EPSDE,0,FPFN2R,"N2R_NP1")
                  DCPPL=RGAMH/(3.D0*RINT0)*( &
                       (AMFD(NSB)*PTFP0(NSA))                     &
                       /(AMFP(NSA)*PTFD0(NSB))*RINT2 )             &
                       *RNFD(NR,NSB)*1.D20
                  FCPPL=0.D0
               ELSE
                  PFPL=PG(NP,NSBA)*PTFP0(NSA)
                  VFPL=PFPL/SQRT(AMFP(NSA)**2+PFPL**2/VC**2)
                  PNFPL=PG(NP,NSBA)
                  PNFP=PNFPL
                  vtatb=(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))
                  ptatb=PG(NP,NSBA)/RGAMA
                  PCRIT=SQRT(vtatb**2/(1.D0-TMC2FD0*vtatb**2*ptatb**2)) &
                       *ptatb
                  IF(PCRIT.le.PMAX(NSB))THEN
                     CALL DEHIFT(RINT0,ES0,H0DE,EPSDE,0,FPFN0R,"N0R_DCPP")
                     PNFP=PCRIT
                     CALL DEFT  (RINT1,ES1,H0DE,EPSDE,0,FPFN1R,"N1R_DCPP")
                     CALL DEHIFT(RINT2,ES2,H0DE,EPSDE,0,FPFN2R,"N2R_DCPP")
                     PNFP=PNFPL
                     DCPPL=RGAMH/(3.D0*RINT0)*(                       &
                          (AMFP(NSA)**2*PTFD0(NSB)**2*RGAMA**3)       &
                          /(AMFD(NSB)**2*PTFP0(NSA)**2*PNFP**3)*RINT1 &
                          +(AMFD(NSB)*PTFP0(NSA))                     &
                          /(AMFP(NSA)*PTFD0(NSB))*RINT2 )             &
                          *RNFD(NR,NSB)*1.D20
                     PNFP=PCRIT
                     CALL DEFT  (RINT4,ES4,H0DE,EPSDE,0,FPFN4R,"N4R_FCPP")
                     CALL DEFT  (RINT5,ES5,H0DE,EPSDE,0,FPFN5R,"N5R_FCPP")
                     CALL DEHIFT(RINT6,ES6,H0DE,EPSDE,0,FPFN6R,"N6R_FCPP")
                     PNFP=PNFPL
                     FCPPL=-RGAMH/(3.D0*RINT0)*(                &
                          (AMFP(NSA)*RGAMA**2)/(AMFD(NSB)*PNFP**2) &
                          *(3.D0*RINT4-TMC2FD0*RINT5)           &
                          +2.D0*(PTFP0(NSA)*PNFP)                  &
                          /(PTFD0(NSB)*RGAMA)*THETA0(NSA)*RINT6 )   &
                          *RNFD(NR,NSB)*1.D20
                  ELSEIF(PCRIT.gt.PMAX(NSB))THEN
                     CALL DEHIFT(RINT0,ES0,H0DE,EPSDE,0,FPFN0R,"N0R_PMAX_DCPP")
                     PNFP=PMAX(NSBA)
                     CALL DEFT  (RINT1,ES1,H0DE,EPSDE,0,FPFN1R,"N1R_PMAX_DCPP")
                     PNFP=PCRIT
                     CALL DEFT  (RINT7,ES7,H0DE,EPSDE,0,FPFN7R,"N7R_PMAX_DCPP")
                     CALL DEHIFT(RINT2,ES2,H0DE,EPSDE,0,FPFN2R,"N2R_PMAX_DCPP")
                     PNFP=PNFPL
                     DCPPL=RGAMH/(3.D0*RINT0)*(                 &
                          (AMFP(NSA)**2*PTFD0(NSB)**2*RGAMA**3) &
                          /(AMFD(NSB)**2*PTFP0(NSA)**2*PNFP**3) &
                          *(RINT1+RINT7)                        &
                          +(AMFD(NSB)*PTFP0(NSA))                &
                          /(AMFP(NSA)*PTFD0(NSB))*RINT2 )       &
                          *RNFD(NR,NSB)*1.D20
                     PNFP=PMAX(NSBA)
                     CALL DEFT  (RINT4,ES4,H0DE,EPSDE,0,FPFN4R,"N4R_PMAX_FCPP")
                     CALL DEFT  (RINT5,ES5,H0DE,EPSDE,0,FPFN5R,"N5R_PMAX_FCPP")
                     PNFP=PCRIT
                     CALL DEFT  (RINT8,ES8,H0DE,EPSDE,0,FPFN9R,"N9R_PMAX_FCPP")
                     CALL DEFT  (RINT9,ES9,H0DE,EPSDE,0,FPFN10R,"N10R_PMAX_FCPP")
                     CALL DEHIFT(RINT6,ES6,H0DE,EPSDE,0,FPFN6R,"N6R_PMAX_FCPP")
                     PNFP=PNFPL
                     FCPPL=-RGAMH/(3.D0*RINT0)*(                &
                          (AMFP(NSA)*RGAMA**2)/(AMFD(NSB)*PNFP**2) &
                          *(3.D0*(RINT4+RINT8)-TMC2FD0*(RINT5+RINT9)) &
                          +2.D0*(PTFP0(NSA)*PNFP)/(PTFD0(NSB)*RGAMA) &
                          *THETA0(NSA)*RINT6 ) &
                          *RNFD(NR,NSB)*1.D20
                  ENDIF
               ENDIF
               DO NTH=1,NTHMAX
                  DCPP2(NTH,NP,NR,NSB,NSA)=DCPP2(NTH,NP,NR,NSB,NSA)+DCPPL
                  FCPP2(NTH,NP,NR,NSB,NSA)=FCPP2(NTH,NP,NR,NSB,NSA)+FCPPL
               ENDDO
            ENDDO

            DO NP=NPSTARTW,NPENDWM
               PFPL=PM(NP,NSA)*PTFP0(NSA)
               VFPL=PFPL/SQRT(AMFP(NSA)**2+PTFP(NR,NSA)**2/VC**2)
               PNFPL=PM(NP,NSBA)
               PNFP=PNFPL
               vtatb=(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))
               ptatb=PM(NP,NSA)/RGAMA
               PCRIT=SQRT(vtatb**2/(1.D0-TMC2FD0*vtatb**2*ptatb**2)) &
                    *ptatb
               IF(PCRIT.le.PMAX(NSBA))THEN
                  CALL DEHIFT(RINT0,ES0,H0DE,EPSDE,0,FPFN0R,"N0R_DCTT")
                  PNFP=PCRIT
                  CALL DEFT  (RINT1,ES1,H0DE,EPSDE,0,FPFN1R,"N1R_DCTT")
                  CALL DEHIFT(RINT2,ES2,H0DE,EPSDE,0,FPFN2R,"N2R_DCTT")
                  CALL DEFT  (RINT3,ES3,H0DE,EPSDE,0,FPFN3R,"N3R_DCTT")
                  PNFP=PNFPL
                  DCTTL=RGAMH/(3.D0*RINT0)*( &
                       1.5D0*RGAMA/PNFP*RINT3 &
                       -0.5D0*(AMFP(NSA)**2*PTFD0(NSB)**2*RGAMA**3) &
                       /(AMFD(NSB)**2*PTFP0(NSA)**2*PNFP**3)*RINT1 &
                       +(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))*RINT2 ) &
                       *RNFD(NR,NSB)*1.D20
               ELSEIF(PCRIT.gt.PMAX(NSB))THEN
                  CALL DEHIFT(RINT0,ES0,H0DE,EPSDE,0,FPFN0R,"N0R_PMAX_DCTT")
                  PNFP=PMAX(NSBA)
                  CALL DEFT  (RINT1,ES1,H0DE,EPSDE,0,FPFN1R,"N1R_PMAX_DCTT")
                  CALL DEFT  (RINT3,ES3,H0DE,EPSDE,0,FPFN3R,"N3R_PMAX_DCTT")
                  PNFP=PNFPL*PTFP0(NSA)*AMFD(NSB)/(PTFD0(NSB)*AMFP(NSA))
                  CALL DEFT  (RINT4,ES1,H0DE,EPSDE,0,FPFN7R,"N7R_PMAX_DCTT")
                  CALL DEFT  (RINT5,ES3,H0DE,EPSDE,0,FPFN8R,"N8R_PMAX_DCTT")
                  CALL DEHIFT(RINT2,ES2,H0DE,EPSDE,0,FPFN2R,"N2R_PMAX_DCTT")
                  PNFP=PNFPL
                  DCTTL=RGAMH/(3.D0*RINT0)*( &
                       1.5D0*RGAMA/PNFP*(RINT3+RINT5) &
                       -0.5D0*(AMFP(NSA)**2*PTFD0(NSB)**2*RGAMA**3) &
                       /(AMFD(NSB)**2*PTFP0(NSA)**2*PNFP**3)*(RINT1+RINT4) &
                       +(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))*RINT2 ) &
                       *RNFD(NR,NSB)*1.D20
              ENDIF
               
               DO NTH=1,NTHMAX+1
                  DCTT2(NTH,NP,NR,NSB,NSA)=DCTT2(NTH,NP,NR,NSB,NSA)+DCTTL
               ENDDO
            ENDDO
         ENDIF
!     ----- Relativistic -----
!
      ELSE
         PMAXC=PMAX(NSBA)
         IF(MODEL_DISRUPT.eq.0)THEN
            RTFDL=RTFD(NR,NSB)
            RTFD0L=RTFD0(NSB)
            P_thermal=SQRT(RTFDL*1.D3*AEE*AMFD(NSB))
         ELSEIF(MODEL_DISRUPT.ge.1)THEN
            RGAMH=AEFP(NSA)**2*AEFD(NSB)**2*POST_LNLAM(NR,NSB,NSA)/(4.D0*PI*EPS0**2) &
                 *AMFP(NSA)/PTFP0(NSA)**3 
            IF(ISW_NOTAIL.eq.0)THEN
               RTFDL=RT_quench(NR) ! [keV] ! thermal quench
            ELSE
               RTFDL=RTFD(NR,NSB) ! no tail effect
            END IF
            RTFD0L=RTFD0(NSB)
            P_thermal=SQRT(RTFDL*1.D3*AEE*AMFD(NSB))
         END IF
         pe_thermal=SQRT(RTFDL*1.D3*AEE*AMFD(NSA))
         TMC2FD =(P_thermal/(AMFD(NSB)*VC))**2
         TMC2FD0=(PTFD0(NSB)/(AMFD(NSB)*VC))**2
         v_thermal=SQRT( RTFDL*1.D3*AEE/AMFD(NSB)*(1.D0-2.5D0*TMC2FD+55.D0/8.D0*TMC2FD**2) )
         
         DO NP=NPSTART,NPENDWG
            IF(NP.EQ.1) THEN
               PNFPL=PG(NP,NSBA)
               PNFP=PNFPL
               RGAMA=SQRT(1.D0+PNFP**2*THETA0(NSA))
               PCRIT=0.D0
               CALL DEHIFT(RINT0,ES0,H0DE,EPSDE,0,FPFN0R,"DEHIFT_0R")
               PNFP=PCRIT
               CALL DEHIFT(RINT2,ES2,H0DE,EPSDE,0,FPFN2R,"DEHIFT_2R")
               DCPPL=RGAMH/(3.D0*RINT0)*( &
                    (AMFD(NSB)*PTFP0(NSA))                     &
                    /(AMFP(NSA)*PTFD0(NSB))*RINT2 )             &
                    *RNFD(NR,NSB)*1.D20
               FCPPL=0.D0
            ELSE
               PFPL=PG(NP,NSBA)*PTFP0(NSA)
               VFPL=PFPL/SQRT(AMFP(NSA)**2+PFPL**2/VC**2)
               PNFPL=PG(NP,NSBA)
               PNFP=PNFPL
               RGAMA=SQRT(1.D0+PNFP**2*THETA0(NSA))
               vtatb=(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))
               ptatb=PG(NP,NSBA)/RGAMA
               PCRIT=SQRT(vtatb**2/(1.D0-TMC2FD0*vtatb**2*ptatb**2)) &
                    *ptatb

               EX=(1.D0-SQRT(1.D0+PG(NP,NSA)**2*TMC2FD0))/TMC2FD

               IF(EX.ge.-500.D0.and.NSA.eq.NSB)THEN
!               IF(PTFP0(NSA)/Pe_thermal*PG(NP,NSBA).le.PMAX(NSA).and.NSA.eq.NSB) THEN ! normal velocity
!               IF(PTFP0(NSA)/Pe_thermal*PG(NP,NSBA).le.PMAX(NSA)) THEN ! normal velocity
!               IF(PTFP0(NSA)*PG(NP,NSBA).le.PTFP0(NSA)*PMAX(NSA)*2) THEN ! normal velocity
!                  IF(PCRIT.le.PMAX(NSB))THEN
                  IF(PCRIT.le.PMAX(NSBA))THEN
                     CALL DEHIFT(RINT0,ES0,H0DE,EPSDE,0,FPFN0R,"DEHIFT_0R_2")
                     PNFP=PCRIT
                     CALL DEFT  (RINT1,ES1,H0DE,EPSDE,0,FPFN1R,"DEFT1_le")
                     CALL DEHIFT(RINT2,ES2,H0DE,EPSDE,0,FPFN2R,"DEHIFT_2R_2")
                     PNFP=PNFPL
                     DCPPL=RGAMH/(3.D0*RINT0)*(                       &
                          (AMFP(NSA)**2*PTFD0(NSB)**2*RGAMA**3)       &
                          /(AMFD(NSB)**2*PTFP0(NSA)**2*PNFP**3)*RINT1 &
                          +(AMFD(NSB)*PTFP0(NSA))                     &
                          /(AMFP(NSA)*PTFD0(NSB))*RINT2 )             &
                          *RNFD(NR,NSB)*1.D20
                     PNFP=PCRIT
                     CALL DEFT  (RINT4,ES4,H0DE,EPSDE,0,FPFN4R,"DEFT4_le")
                     CALL DEFT  (RINT5,ES5,H0DE,EPSDE,0,FPFN5R,"DEFT5_le")
                     CALL DEHIFT(RINT6,ES6,H0DE,EPSDE,0,FPFN6R,"DEHIFT_6R")
                     PNFP=PNFPL
                     FCPPL=-RGAMH/(3.D0*RINT0)*(                &
                          (AMFP(NSA)*RGAMA**2)/(AMFD(NSB)*PNFP**2) &
                          *(3.D0*RINT4-TMC2FD0*RINT5)           &
                          +2.D0*(PTFP0(NSA)*PNFP)                  &
                          /(PTFD0(NSB)*RGAMA)*THETA0(NSA)*RINT6 )   &
                          *RNFD(NR,NSB)*1.D20
!                     WRITE(6,'(A,2I5,1P10E14.6)') "low v e-e ", NP, NSB, TMC2FD, RINT0, RINT1, PCRIT, DCPPL, FCPPL
!                  ELSEIF(PCRIT.gt.PMAX(NSB))THEN
                  ELSEIF(PCRIT.gt.PMAX(NSBA))THEN
                     CALL DEHIFT(RINT0,ES0,H0DE,EPSDE,0,FPFN0R,"DEHIFT_0R_3")
                     PNFP=PMAX(NSBA)
                     CALL DEFT  (RINT1,ES1,H0DE,EPSDE,0,FPFN1R,"DEFT1_gt")
                     PNFP=PCRIT
!                     CALL DEFT  (RINT7,ES7,H0DE,EPSDE,0,FPFN7R,"DEFT7_gt")
                     RINT7=0.D0
!                     CALL DEHIFT(RINT2,ES2,H0DE,EPSDE,0,FPFN2R,"DEHIFT_2R_3")
                     RINT2=0.D0
                     PNFP=PNFPL
                     DCPPL=RGAMH/(3.D0*RINT0)*(                 &
                          (AMFP(NSA)**2*PTFD0(NSB)**2*RGAMA**3) &
                          /(AMFD(NSB)**2*PTFP0(NSA)**2*PNFP**3) &
                          *(RINT1+RINT7)                        &
                          +(AMFD(NSB)*PTFP0(NSA))                &
                          /(AMFP(NSA)*PTFD0(NSB))*RINT2 )       &
                          *RNFD(NR,NSB)*1.D20
                     PNFP=PMAX(NSBA)
                     CALL DEFT  (RINT4,ES4,H0DE,EPSDE,0,FPFN4R,"DEFT4_gt")
                     CALL DEFT  (RINT5,ES5,H0DE,EPSDE,0,FPFN5R,"DEFT5_gt")
                     PNFP=PCRIT
                     RINT8=0.D0
                     RINT9=0.D0
                     RINT6=0.D0
!                     CALL DEFT  (RINT8,ES8,H0DE,EPSDE,0,FPFN9R,"DEFT9_gt")
!                     CALL DEFT  (RINT9,ES9,H0DE,EPSDE,0,FPFN10R,"DEFT10_gt")
!                     CALL DEHIFT(RINT6,ES6,H0DE,EPSDE,0,FPFN6R,"DEHIFT_6R_2")
                     PNFP=PNFPL
                     FCPPL=-RGAMH/(3.D0*RINT0)*(                &
                          (AMFP(NSA)*RGAMA**2)/(AMFD(NSB)*PNFP**2) &
                          *(3.D0*(RINT4+RINT8)-TMC2FD0*(RINT5+RINT9)) &
                          +2.D0*(PTFP0(NSA)*PNFP)/(PTFD0(NSB)*RGAMA) &
                          *THETA0(NSA)*RINT6 ) &
                          *RNFD(NR,NSB)*1.D20
!                     WRITE(6,'(A,2I5,1P10E14.6)') "low v e-i ", NP, NSB, TMC2FD, RINT0, RINT1, PCRIT, DCPPL, FCPPL
                  ENDIF
               ELSE! high velocity limit 
                  DCPPL = RGAMH*(AMFP(NSA)/PTFP0(NSA))**2*(RGAMA/PG(NP,NSBA))**3 &
                       *v_thermal**2 * RNFD(NR,NSB)*1.D20
                  FCPPL =-RGAMH*(RGAMA/PG(NP,NSBA))**2 &
                       *(1.D0-2.5D0*TMC2FD+55.D0/8.D0*TMC2FD**2) &
                       * RNFD(NR,NSB)*1.D20*AMFP(NSA)/AMFD(NSB)

!                  WRITE(6,'(A,2I5,1P4E14.6)') "high v ", NP, NSB, TMC2FD, v_thermal, DCPPL, FCPPL
               END IF
            END IF
            DO NTH=1,NTHMAX
               DCPP2(NTH,NP,NR,NSB,NSA)=DCPP2(NTH,NP,NR,NSB,NSA)+DCPPL
               FCPP2(NTH,NP,NR,NSB,NSA)=FCPP2(NTH,NP,NR,NSB,NSA)+FCPPL
            ENDDO
         ENDDO

         DO NP=NPSTARTW,NPENDWM
            PFPL=PM(NP,NSA)*PTFP0(NSA)
            VFPL=PFPL/SQRT(AMFP(NSA)**2+PTFP(NR,NSA)**2/VC**2)
            PNFPL=PM(NP,NSBA)
            PNFP=PNFPL
            RGAMA=SQRT(1.D0+PNFP**2*THETA0(NSA))
            vtatb=(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))
            ptatb=PM(NP,NSA)/RGAMA
            PCRIT=SQRT(vtatb**2/(1.D0-TMC2FD0*vtatb**2*ptatb**2)) &
                 *ptatb

            EX=(1.D0-SQRT(1.D0+PM(NP,NSA)**2*TMC2FD0))/TMC2FD
            IF(EX.ge.-500.D0.and.NSA.eq.NSB)THEN
!            IF(PTFP0(NSA)/Pe_thermal*PG(NP,NSBA).le.PMAX(NSA).and.NSA.eq.NSB) THEN ! normal velocity
!            IF(PTFP0(NSA)/P_thermal*PG(NP,NSBA).le.PMAX(NSB))THEN
!            IF(PTFP0(NSA)*PG(NP,NSBA).le.PE_thermal*PMAX(NSA)*2) THEN ! normal velocity
               IF(PCRIT.le.PMAX(NSBA))THEN
                  CALL DEHIFT(RINT0,ES0,H0DE,EPSDE,0,FPFN0R,"DEHIFT_0R_4")
                  PNFP=PCRIT
                  CALL DEFT  (RINT1,ES1,H0DE,EPSDE,0,FPFN1R,"DEFT1_le_tt")
                  CALL DEHIFT(RINT2,ES2,H0DE,EPSDE,0,FPFN2R,"DEHIFT_2R_4")
                  CALL DEFT  (RINT3,ES3,H0DE,EPSDE,0,FPFN3R,"DEFT3_le_tt")
                  PNFP=PNFPL
                  DCTTL=RGAMH/(3.D0*RINT0)*( &
                       1.5D0*RGAMA/PNFP*RINT3 &
                       -0.5D0*(AMFP(NSA)**2*PTFD0(NSB)**2*RGAMA**3) &
                       /(AMFD(NSB)**2*PTFP0(NSA)**2*PNFP**3)*RINT1 &
                       +(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))*RINT2 ) &
                       *RNFD(NR,NSB)*1.D20
!               ELSEIF(PCRIT.gt.PMAX(NSB))THEN
               ELSEIF(PCRIT.gt.PMAX(NSBA))THEN
                  CALL DEHIFT(RINT0,ES0,H0DE,EPSDE,0,FPFN0R,"DEHIFT_0R_5")
                  PNFP=PMAX(NSBA)
                  CALL DEFT  (RINT1,ES1,H0DE,EPSDE,0,FPFN1R,"DEFT1_gt_tt")
                  CALL DEFT  (RINT3,ES3,H0DE,EPSDE,0,FPFN3R,"DEFT3_gt_tt")
                  PNFP=PNFPL*PTFP0(NSA)*AMFD(NSB)/(PTFD0(NSB)*AMFP(NSA))
!                  CALL DEFT  (RINT4,ES1,H0DE,EPSDE,0,FPFN7R,"DEFT7_gt_tt")
!                  CALL DEFT  (RINT5,ES3,H0DE,EPSDE,0,FPFN8R,"DEFT8_gt_tt")
!                  CALL DEHIFT(RINT2,ES2,H0DE,EPSDE,0,FPFN2R,"DEHIFT_2R_5")
                  RINT4=0.D0
                  RINT5=0.D0
                  RINT2=0.D0
                  PNFP=PNFPL
                  DCTTL=RGAMH/(3.D0*RINT0)*( &
                       1.5D0*RGAMA/PNFP*(RINT3+RINT5) &
                       -0.5D0*(AMFP(NSA)**2*PTFD0(NSB)**2*RGAMA**3) &
                       /(AMFD(NSB)**2*PTFP0(NSA)**2*PNFP**3)*(RINT1+RINT4) &
                       +(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))*RINT2 ) &
                       *RNFD(NR,NSB)*1.D20
               ENDIF
            ELSE
!               DCTTL=RGAMH/PG(NP,NSBA)*(1.D0-P_thermal**2/PFPL**2)
               DCTTL=0.5D0*RGAMH*RGAMA/PM(NP,NSBA)* &
                    (1.D0-(v_thermal*AMFP(NSA)/PTFP0(NSA))**2*(RGAMA/PM(NP,NSBA))**2) &
                    * RNFD(NR,NSB)*1.D20
            END IF
            DO NTH=1,NTHMAX+1
               DCTT2(NTH,NP,NR,NSB,NSA)=DCTT2(NTH,NP,NR,NSB,NSA)+DCTTL
            ENDDO
         ENDDO
      ENDIF

      RETURN
      END SUBROUTINE FPCALC_L
!
! ************************************************************
!
!       BOUNCE AVERAGE OF LINEAR COLLISIONAL OPERATOR
!
! ************************************************************
!
      SUBROUTINE FPCALC_LAV(NR,NSA)

      IMPLICIT NONE

      integer:: NSA, NSB,NR, NP, NTH
      real(8):: RGAMH, RGAMH2, RZI, RTE, PFPL, VFPL, U, DCTTL
      real(8):: FACT
      real(8):: PNFP,TMC2FD,TMC2FD0 
      DOUBLE PRECISION:: DELH, sum, etal, psib, pcos, arg, x
      INTEGER:: NG

      DO NSB=1,NSBMAX
!         DO NP=1,NPMAX+1
         DO NP=NPSTART,NPENDWG
            DO NTH=1,NTHMAX
               FACT=RLAMDA(NTH,NR)
               DCPP2(NTH,NP,NR,NSB,NSA) &
                    =FACT*DCPP2(NTH,NP,NR,NSB,NSA)
               FCPP2(NTH,NP,NR,NSB,NSA) &
                    =FACT*FCPP2(NTH,NP,NR,NSB,NSA)
            ENDDO
            DCPP2(ITL(NR),NP,NR,NSB,NSA) &
                     =RLAMDA(ITL(NR),NR)/4.D0 &
                       *( DCPP2(ITL(NR)-1,NP,NR,NSB,NSA) &
                                   /RLAMDA(ITL(NR)-1,NR) &
                         +DCPP2(ITL(NR)+1,NP,NR,NSB,NSA) &
                                   /RLAMDA(ITL(NR)+1,NR) &
                         +DCPP2(ITU(NR)-1,NP,NR,NSB,NSA) &
                                   /RLAMDA(ITU(NR)-1,NR) &
                         +DCPP2(ITU(NR)+1,NP,NR,NSB,NSA) &
                                   /RLAMDA(ITU(NR)+1,NR))

            FCPP2(ITL(NR),NP,NR,NSB,NSA) &
                    =RLAMDA(ITL(NR),NR)/4.D0  &
                      *( FCPP2(ITL(NR)-1,NP,NR,NSB,NSA) &
                                    /RLAMDA(ITL(NR)-1,NR) &
                        +FCPP2(ITL(NR)+1,NP,NR,NSB,NSA) &
                                    /RLAMDA(ITL(NR)+1,NR) &
                        +FCPP2(ITU(NR)-1,NP,NR,NSB,NSA) &
                                    /RLAMDA(ITU(NR)-1,NR) &
                        +FCPP2(ITU(NR)+1,NP,NR,NSB,NSA) &
                                    /RLAMDA(ITU(NR)+1,NR)) 
            DCPP2(ITU(NR),NP,NR,NSB,NSA)=DCPP2(ITL(NR),NP,NR,NSB,NSA)
            FCPP2(ITU(NR),NP,NR,NSB,NSA)=FCPP2(ITL(NR),NP,NR,NSB,NSA)
         END DO
!         FACT=RLAMDC(NTH,NR)
!         DO NP=1,NPMAX
!            DO NTH=1,NTHMAX+1
!               DCTT2(NTH,NP,NR,NSB,NSA)=FACT*DCTT2(NTH,NP,NR,NSB,NSA)
!            ENDDO
!         END DO
!         DO NP=1,NPMAX
         DO NP=NPSTARTW,NPENDWM
            DO NTH=1,NTHMAX+1
               IF(NTH.NE.NTHMAX/2+1) THEN 
                  DELH = 2.D0*ETAG(NTH,NR)/NAVMAX
                  SUM=0.D0
                  DO NG=1,NAVMAX
                     ETAL = DELH*(NG-0.5D0)
                     X=EPSRM2(NR)*COS(ETAL)*RR 
                     PSIB=(1.D0+EPSRM2(NR))/(1.D0+X/RR)
                     ARG=1.D0-PSIB*SING(NTH)**2
                     PCOS = SQRT(ARG)
                     sum=sum + DCTT2(NTH,NP,NR,NSB,NSA)*PCOS/(PSIB*COSG(NTH)) 
                  END DO
                  DCTT2(NTH,NP,NR,NSB,NSA)=sum*DELH
               ELSE
                  DCTT2(NTH,NP,NR,NSB,NSA)=0.D0
               END IF
            ENDDO
            DO NTH=ITL(NR)+1,NTHMAX/2
               DCTT2(NTH,NP,NR,NSB,NSA)        &
                    =(DCTT2(NTH,NP,NR,NSB,NSA) &
                    +DCTT2(NTHMAX-NTH+2,NP,NR,NSB,NSA))*0.5D0
               DCTT2(NTHMAX-NTH+2,NP,NR,NSB,NSA) &
                    =DCTT2(NTH,NP,NR,NSB,NSA) 
            END DO
         END DO ! NP
      END DO ! NSB

      RETURN
      END SUBROUTINE FPCALC_LAV
      END MODULE fpcalc
