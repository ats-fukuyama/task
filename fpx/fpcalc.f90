! fpcalc.f90
!
! ******************************************************************
!
!     CALCULATION OF COLLISINAL DIFUSION COEFFICIENTS: DC AND FC
!
! ******************************************************************
!

      MODULE fpcalc

      USE fpcomm
      USE fpmpi
      USE libmpi
      USE libspf, ONLY: ERF0,ERF1
      REAL(rkind):: PMAXC
      integer:: NSB_ISO
      REAL(rkind):: THETA0L_C, THETAL_C
      REAL(rkind):: RTFD0L_C, RTFDL_C, PNFP_C
      REAL(rkind):: RNFD0L_C, RNFDL_C

      integer,parameter:: ISW_NOTAIL=0
      integer,parameter:: MODEL_DE=1

      contains

!----------------------------------

      SUBROUTINE FP_CALC

      USE libmtx
      USE fpcalcn, ONLY: FPCALC_NL
      USE fpcalcnr, ONLY: FPCALC_NLR 
      IMPLICIT NONE
      integer:: NSA,NSB,NSSB,NR, NP, NTH,NS
      REAL(rkind):: RGAMH, RGAMH2, RZI, PFPL, VFPL, DCTTL, RGAMA

      DO NSA=NSASTART,NSAEND
         NS=NS_NSA(NSA)
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
         
         DO NR=NRSTART,NREND
            DO NSB=1,NSBMAX
               NSSB=NS_NSB(NSB)
!            if(nr.eq.1) write(6,'(A,I8,1P2E12.4)') 
!     &           ' NSB,RN,RNFD=',NSB,RN(NSB),RNFD(NR,NSB)
!
               IF(MODELC(NSSB).EQ.0.or. &
                  MODELC(NSSB).eq.1.or.MODELC(NS).eq.2) THEN
                  CALL FPCALC_L(NR,NSB,NSA)
               ELSEIF(MODELC(NSSB).EQ.4) THEN
                  IF(MODELR.eq.0)THEN
                     CALL FPCALC_NL(NR,NSB,NSA)
                  ELSE IF(MODELR.eq.1)THEN
                     CALL FPCALC_NLR(NR,NSB,NSA)
                  END IF
               ELSEIF(MODELC(NSSB).EQ.5) THEN
                  IF(NS_NSB(NSB).NE.NS_NSA(NSA)) THEN
                     CALL FPCALC_L(NR,NSB,NSA)
                  ENDIF
               ELSEIF(MODELC(NSSB).EQ.6) THEN
                  IF(NS_NSB(NSB).NE.NS_NSA(NSA)) THEN
                     CALL FPCALC_NL(NR,NSB,NSA)
                  ENDIF
               ELSEIF(MODELC(NSSB).EQ.7) THEN
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
               ELSEIF(MODELC(NSSB).EQ.-1) THEN
                  IF(NS_NSB(NSB).EQ.NS_NSA(NSA)) THEN
                     MODELC(NSSB)=0
                     CALL FPCALC_L(NR,NSB,NSA)
                     MODELC(NSSB)=-1
                  ENDIF
               ELSEIF(MODELC(NSSB).EQ.-2) THEN 
                  IF(NS_NSB(NSB).EQ.NS_NSA(NSA)) THEN
                     MODELC(NSSB)=1
                     CALL FPCALC_L(NR,NSB,NSA)
                     MODELC(NSSB)=-2
                  ENDIF
               ELSEIF(MODELC(NSSB).EQ.-4) THEN
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
            IF(MODELC(NS).LT.0) THEN
               IF(NS_NSA(NSA).EQ.1) THEN
                  DO NSB=1,NSBMAX
                     IF(NS_NSB(NSB).EQ.2) THEN
                        RGAMH=RNUD(NR,NSA,NSA)*SQRT(2.D0)*VTFD(NR,NSA) &
                             *AMFP(NSA) &
                             /(RNFP0(NSA)*PTFP0(NSA)*1.D20)*RNFD0(NSA)
                        RGAMH2=RGAMH*RNFP(NR,NSA)*1.D20*PTFP0(NSA) &
                             /AMFP(NSA)/RNFD0(NSA)
                        rZI = -PZ(2)/PZ(1) &
!                     rZI = ( PZ(2)/PZ(1) )**2 * RNFD(NR,NSB)/RNFP(NR,NSA) &
                          /(14.9D0-0.5D0*LOG(RNFP(NR,NSA))+LOG(RTFP(NR,NSA))) &
                          *(15.2D0-0.5D0*LOG(RNFP(NR,NSA))+LOG(RTFP(NR,NSA)))
                        DO NP=NPSTARTW,NPENDWM
                           RGAMA =SQRT(1.D0+PM(NP,NS)**2*THETA0(NS))
                           PFPL=PM(NP,NS)*PTFP0(NSA)
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
               CALL FPCALC_LAV(NR,NSA)
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
               END DO
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
      END DO ! NSA

      IF(nrank.eq.0) THEN
      IF(MOD(IDBGFP/8,2).EQ.1) THEN
!
!     +++ plot of D_coll +++
!
         DO NSA=NSASTART,NSAEND
         CALL PAGES
         CALL FPGRFA(1,DCPP,PG,1,'@DCPP(P)@',NTHMAX+1,NPMAX+1,NRMAX+1, &
                                             NTHMAX,NPMAX,NRMAX,NSA)
         CALL FPGRFA(2,DCTT,PM,2,'@DCTT(P)@',NTHMAX+1,NPMAX+1,NRMAX+1, &
                                             NTHMAX,NPMAX,NRMAX,NSA)
         CALL FPGRFA(3,FCPP,PG,1,'@FCPP(P)@',NTHMAX+1,NPMAX+1,NRMAX+1, &
                                             NTHMAX,NPMAX,NRMAX,NSA)
         CALL PAGEE
         END DO
      ENDIF


      IF(MOD(IDBGFP/16,2).EQ.1) THEN
!
!     +++ plot of D_coll +++
!
         DO NSA=NSASTART,NSAEND
         CALL PAGES
         CALL FPGRFB(1,DCPP,THM,1,'@DCPP(TH)@',NTHMAX+1,NPMAX+1,NRMAX+1, &
                                               NTHMAX,NPMAX,NRMAX,NSA)
         CALL FPGRFB(2,DCTT,THG,2,'@DCTT(TH)@',NTHMAX+1,NPMAX+1,NRMAX+1, &
                                               NTHMAX,NPMAX,NRMAX,NSA)
         CALL FPGRFB(3,FCPP,THM,1,'@FCPP(TH)@',NTHMAX+1,NPMAX+1,NRMAX+1, &
                                               NTHMAX,NPMAX,NRMAX,NSA)
         CALL PAGEE
         END DO
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
      REAL(rkind), dimension(NTHM,NPM,NRM,NSA), intent(IN):: DATA
      REAL(rkind), dimension(NPM), intent(IN):: P
      character(len=*):: title
      REAL(rkind), dimension(NPM,9):: WORK
      integer, dimension(9):: NTHG
      integer:: nth,np,nr,ng,npmaxg

      if(NRMAX.GT.1) then
    1    WRITE(6,*) '## NR ?'
         READ(5,*,err=1,end=9999) NR
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
      REAL(rkind), dimension(NTHM,NPM,NRM,NSA), intent(IN):: DATA
      REAL(rkind), dimension(NTHM), intent(IN):: TH
      character(len=*):: title
      REAL(rkind), dimension(NTHM,10):: WORK
      integer, dimension(10):: NPG
      integer:: nth,np,nr,ng,nthmaxg

      if(NRMAX.GT.1) then
    1    WRITE(6,*) '## NR ?'
         READ(5,*,err=1,end=9999) NR
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
      REAL(rkind)::FPFN0R
      REAL(rkind),INTENT(IN):: X
      REAL(rkind)::PN, FACT
!
      IF(MODEL_DE.eq.0)THEN
         FACT=1.D0
      ELSE
         FACT=SQRT(RTFDL_C/RTFD0L_C)
      END IF

!      PN=X
      PN=X*FACT
      FPFN0R=PN**2*FPRMXW(PN)*FACT
!
      RETURN
      END FUNCTION FPFN0R
!---------------------------------------------------------------
      FUNCTION FPFN1R(X,XM,XP)

      IMPLICIT NONE
      REAL(rkind)::FPFN1R
      REAL(rkind),INTENT(IN)::X, XM, XP
      REAL(rkind)::PN, A, B
      REAL(RKIND):: DUMMY

      DUMMY=X
      DUMMY=XM
      A=0.5D0*PNFP_C
      PN=A*XP
      B=PN**4/(1.D0+PN**2*THETA0L_C)
      FPFN1R=A*B*FPRMXW(PN)

      RETURN
      END FUNCTION FPFN1R
!
! ===============================================================
!
      FUNCTION FPFN2R(X)

      REAL(rkind):: FPFN2R
      REAL(rkind),INTENT(IN):: X
      REAL(rkind):: A, PN, B, FACT

      IF(MODEL_DE.eq.0)THEN
         FACT=1.D0
      ELSE
         FACT=SQRT(RTFDL_C/RTFD0L_C)
      END IF

      A=1.D0
      PN=A*(X*FACT+PNFP_C)
      IF(MODEL_DE.eq.0)THEN
         B=PN*SQRT(1.D0+PN**2*THETA0L_C)
      ELSE
         B=PN*SQRT(1.D0+PN**2*THETA0L_C)
      END IF
      FPFN2R=A*B*FPRMXW(PN)*FACT

      RETURN
      END FUNCTION FPFN2R
!
! ==============================================================
!
      FUNCTION FPFN3R(X,XM,XP)

      REAL(rkind):: FPFN3R
      REAL(rkind),INTENT(IN):: X, XM, XP
      REAL(rkind):: A, PN
      REAL(RKIND):: DUMMY

      DUMMY=X
      DUMMY=XM
      A=0.5D0*PNFP_C
      PN=A*XP
      FPFN3R=A*PN**2*FPRMXW(PN)

      RETURN
      END FUNCTION FPFN3R
!
! ===============================================================
!
      FUNCTION FPFN4R(X,XM,XP)

      REAL(rkind):: FPFN4R
      REAL(rkind),INTENT(IN):: X, XM, XP
      REAL(rkind):: A, PN, B
      REAL(RKIND):: DUMMY

      DUMMY=X
      DUMMY=XM
      A=0.5D0*PNFP_C
      PN=A*XP
      B=PN**2/SQRT(1.D0+PN**2*THETA0L_C)
      FPFN4R=A*B*FPRMXW(PN)

      RETURN
      END FUNCTION FPFN4R
!
! ===============================================================
!
      FUNCTION FPFN5R(X,XM,XP)

      REAL(rkind):: FPFN5R
      REAL(rkind),INTENT(IN):: X, XM, XP
      REAL(rkind):: A, PN, B
      REAL(RKIND):: DUMMY

      DUMMY=X
      DUMMY=XM
      A=0.5D0*PNFP_C
      PN=A*XP
      B=PN**4/(SQRT(1.D0+PN**2*THETA0L_C))**3
      FPFN5R=A*B*FPRMXW(PN)

      RETURN
      END FUNCTION FPFN5R
!
! ===============================================================
!
      FUNCTION FPFN6R(X)

      REAL(rkind):: FPFN6R
      REAL(rkind),INTENT(IN):: X
      REAL(rkind):: A, PN, FACT

      IF(MODEL_DE.eq.0)THEN
         FACT=1.D0
      ELSE
         FACT=SQRT(RTFDL_C/RTFD0L_C)
      END IF

      A=1.D0
      PN=A*(X*FACT+PNFP_C)
      FPFN6R=A*PN*FPRMXW(PN)*FACT

      RETURN
      END FUNCTION FPFN6R
!
! =============================================================== 
      FUNCTION FPFN7R(X,XM,XP)!
!                            
      REAL(rkind):: FPFN7R
      REAL(rkind),INTENT(IN):: X, XM, XP
      REAL(rkind):: A, PN, B, PMAX2
      REAL(RKIND):: DUMMY


      DUMMY=X
      DUMMY=XM
      PMAX2=PMAXC
      A=0.5D0*(PNFP_C-PMAX2)
      PN=( A*XP+0.5D0*(PNFP_C+PMAX2) )
      B=PN**4/(1.D0+PN**2*THETA0L_C)
      FPFN7R=A*B*FPRMXW(PN)

      RETURN
      END FUNCTION FPFN7R
!
! =============================================================== 
!
      FUNCTION FPFN8R(X,XM,XP)!

      REAL(rkind):: FPFN8R
      REAL(rkind),INTENT(IN):: X, XM, XP
      REAL(rkind):: A, PN, PMAX2
      REAL(RKIND):: DUMMY

      DUMMY=X
      DUMMY=XM
      PMAX2=PMAXC
      A=0.5D0*(PNFP_C-PMAX2)
      PN=( A*XP+0.5D0*(PNFP_C+PMAX2) )
      FPFN8R=A*PN**2*FPRMXW(PN)

      RETURN
      END FUNCTION FPFN8R
!
! =============================================================== 
!
      FUNCTION FPFN9R(X,XM,XP)!

      REAL(rkind):: FPFN9R
      REAL(rkind),INTENT(IN):: X, XM, XP
      REAL(rkind):: A, PN, B, PMAX2
      REAL(RKIND):: DUMMY

      DUMMY=X
      DUMMY=XM
      PMAX2=PMAXC
      A=0.5D0*(PNFP_C-PMAX2)
      PN=( A*XP+0.5D0*(PNFP_C+PMAX2) )
      B=PN**2/SQRT(1.D0+PN**2*THETA0L_C)
      FPFN9R=A*B*FPRMXW(PN)

      RETURN
      END FUNCTION FPFN9R
!
! =============================================================== 
!
      FUNCTION FPFN10R(X,XM,XP)!

      REAL(rkind):: FPFN10R
      REAL(rkind),INTENT(IN):: X, XM, XP
      REAL(rkind):: A, PN, B, PMAX2
      REAL(rkind):: DUMMY

      DUMMY=X
      DUMMY=XM
      PMAX2=PMAXC
      A=0.5D0*(PNFP_C-PMAX2)
      PN=( A*XP+0.5D0*(PNFP_C+PMAX2) )
      B=PN**4/(SQRT(1.D0+PN**2*THETA0L_C))**3
      FPFN10R=A*B*FPRMXW(PN)

      RETURN
      END FUNCTION FPFN10R
!
! ===============================================================
!
      FUNCTION FPRMXW(PN)

      USE libbes,ONLY: beseknx
      USE plprof
      REAL(rkind):: FPRMXW
      REAL(rkind),INTENT(IN):: PN
      REAL(rkind):: EX
      REAL(rkind):: FACT, DKBSL, Z

      IF(MODELR.eq.1)THEN
         Z=1.D0/THETAL_C
         DKBSL=BESEKNX(2,Z)
         FACT=RNFDL_C*SQRT(THETA0L_C)/(4.D0*PI*RTFDL_C*DKBSL) &
              *RTFD0L_C

         EX=(1.D0-SQRT(1.D0+PN**2*THETA0L_C)) /THETAL_C

!         IF (EX.LT.-100.D0)THEN
!            FPRMXW=0.D0
!         ELSE
            FPRMXW=EXP(EX)*FACT
!         ENDIF
      ELSEIF(MODELR.eq.0)THEN
         EX=-PN**2/(2.D0*RTFDL_C/RTFD0L_C)
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

      integer:: NSA, NSB, NR, NP, NTH, NSSA, NSSB
      REAL(rkind):: RGAMH, PFPL, VFPL, U, DCTTL
      REAL(rkind):: RNNL, RNUFL, RNUDL, DCPPL, FCPPL, V
      REAL(rkind):: PNFPL, RGAMA, vtatb, ptatb, PCRIT
      REAL(rkind):: RINT0, ES0, RINT1, ES1, RINT2, ES2, RINT4, ES4, RINT5, ES5
      REAL(rkind):: RINT6, ES6, RINT7, ES7, RINT8, ES8, RINT9, ES9
      REAL(rkind):: RINT3, ES3, p_thermal, v_thermal, RNFDL

!     ------ define --------
      RNNL=RNFD(NR,NSB)/RNFP0(NSA)
      RNUFL=RNUF(NR,NSB,NSA)*RNNL
      RNUDL=RNUD(NR,NSB,NSA)*RNNL

!      RGAMH=RNUD(NR,NSB,NSA)*SQRT(2.D0)*VTFD(NR,NSB)*AMFP(NSA) &
!              /(RNFP0(NSA)*PTFP0(NSA)*1.D20)
      RGAMH=AEFP(NSA)**2*AEFD(NSB)**2*LNLAM(NR,NSB,NSA)/(4.D0*PI*EPS0**2) &
           *AMFP(NSA)/PTFP0(NSA)**3 
!     RGAMH = \hat{\Gamma}/n_b
      NSSA=NS_NSA(NSA)
      NSSB=NS_NSB(NSB)

!      WRITE(*,*) " "
!      WRITE(*,'(2I4,3E16.8)') NSA,NSB,RGAMH,RNFD(1,NSB),VTFD(1,NSB)
!
!     ----- Non-Relativistic -----
!
      IF(MODELR.EQ.0) THEN 
         IF(MODELC(NSSB).eq.0)THEN ! maxwellian
            IF(MODEL_DISRUPT.eq.0)THEN
               RTFDL_C=RTFD(NR,NSB)
               RTFD0L_C=RTFD0(NSB)
               v_thermal=VTFD(NR,NSB)
            ELSEIF(MODEL_DISRUPT.ge.1)THEN
               RTFDL_C=RT_quench(NR) ! [keV]
               v_thermal=SQRT( RTFDL_C*1.D3*AEE/AMFD(NSB))
            END IF
            DO NP=NPSTART,NPENDWG
               IF(NP.EQ.1) THEN
                  DCPPL=RGAMH*RNFD(NR,NSB)*1.D20*(2.D0/(3.D0*SQRT(PI))) &
                       *(VTFP0(NSA)/(SQRT(2.D0)*VTFD(NR,NSB))) 
                  FCPPL=0.D0
               ELSE
                  PFPL=PG(NP,NSSA)*PTFP0(NSA)
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

            DO NP=NPSTARTW,NPENDWM
               PFPL=PM(NP,NSSA)*PTFP0(NSA)
               VFPL=PFPL/AMFP(NSA)
               V=VFPL/VTFP0(NSA)
               U=VFPL/(SQRT(2.D0)*v_thermal)
               DCTTL= 0.25D0*RNUDL/U &
                    *((2.D0-1.D0/U**2)*ERF0(U)+ERF1(U)/U)
!
               DO NTH=1,NTHMAX+1
                  DCTT2(NTH,NP,NR,NSB,NSA)=DCTT2(NTH,NP,NR,NSB,NSA)+DCTTL
               ENDDO
            ENDDO
         ELSEIF(MODELC(NSSB).ge.2)THEN ! isotoropic
            PMAXC=PMAX(NSSB)
            THETAL_C =0.D0
            THETA0L_C=0.D0
            RGAMA=1.D0
            NSB_ISO=NSB
            IF(MODELC(NSSB).eq.1)THEN ! Collision term for initial temperature
               RTFDL_C=RTFD(NR,NSB)
               RTFD0L_C=RTFD0(NSB)
               RNFDL=RNFD(NR,NSB)
            ELSEIF(MODELC(NSSB).eq.2)THEN ! for updating temperature
               RTFDL_C=RT_TEMP(NR,NSB)
               RTFD0L_C=RTFD0(NSB)
               RNFDL=RN_TEMP(NR,NSB)
               RGAMH=AEFP(NSA)**2*AEFD(NSB)**2*LNLAM(NR,NSB,NSA) &
                    /(4.D0*PI*EPS0**2) &
                    *AMFP(NSA)/PTFP0(NSA)**3 
            END IF
            IF(MODEL_DISRUPT.eq.1)THEN
               RTFDL_C=RT_quench(NR) ! [keV]
               RTFD0L_C=RTFD0(NSB)
               IF(MODEL_IMPURITY.eq.0)THEN 
                  RNFDL_C=RNFD(NR,NSB)
               ELSE
                  RNFDL=RN_MGI(NR,NSB)
               END IF
               RGAMH=AEFP(NSA)**2*AEFD(NSB)**2*POST_LNLAM(NR,NSB,NSA) &
                    /(4.D0*PI*EPS0**2) &
                    *AMFP(NSA)/PTFP0(NSA)**3 
            END IF

            DO NP=NPSTART,NPENDWG
               IF(NP.EQ.1) THEN
                  PNFPL=PG(NP,NSSA)
                  PNFP_C=PNFPL
                  PCRIT=0.D0
                  CALL DEHIFT(RINT0,ES0,H0DE,EPSDE,0,FPFN0R,"N0R_NP1")
                  PNFP_C=PCRIT
                  CALL DEHIFT(RINT2,ES2,H0DE,EPSDE,0,FPFN2R,"N2R_NP1")
                  DCPPL=RGAMH/(3.D0*RINT0)*( &
                       (AMFD(NSB)*PTFP0(NSA))                     &
                       /(AMFP(NSA)*PTFD0(NSB))*RINT2 )             &
                       *RNFDL*1.D20
                  FCPPL=0.D0
               ELSE
!                  PFPL=PG(NP,NSBA)*PTFP0(NSA)
!                  VFPL=PFPL/SQRT(AMFP(NSA)**2+PFPL**2/VC**2)
                  PNFPL=PG(NP,NSSA)
                  PNFP_C=PNFPL
                  vtatb=(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))
                  ptatb=PG(NP,NSSA)/RGAMA
!                  PCRIT=SQRT(vtatb**2/(1.D0-THETA0L_C*vtatb**2*ptatb**2)) &
!                       *ptatb
                  PCRIT=(AMFD(NSB)*PTFP0(NSA)) &
                       /(AMFP(NSA)*PTFD0(NSB))*PG(NP,NSSA)
                  IF(PCRIT.le.PMAX(NSSB))THEN
                     CALL DEHIFT(RINT0,ES0,H0DE,EPSDE,0,FPFN0R,"N0R_DCPP")
                     PNFP_C=PCRIT
                     CALL DEFT  (RINT1,ES1,H0DE,EPSDE,0,FPFN1R,"N1R_DCPP")
                     CALL DEHIFT(RINT2,ES2,H0DE,EPSDE,0,FPFN2R,"N2R_DCPP")
                     PNFP_C=PNFPL
                     DCPPL=RGAMH/(3.D0*RINT0)*(                       &
                          (AMFP(NSA)**2*PTFD0(NSB)**2*RGAMA**3)       &
                          /(AMFD(NSB)**2*PTFP0(NSA)**2*PNFP_C**3)*RINT1 &
                          +(AMFD(NSB)*PTFP0(NSA))                     &
                          /(AMFP(NSA)*PTFD0(NSB))*RINT2 )             &
                          *RNFDL*1.D20
                     PNFP_C=PCRIT
                     CALL DEFT  (RINT4,ES4,H0DE,EPSDE,0,FPFN4R,"N4R_FCPP")
                     CALL DEFT  (RINT5,ES5,H0DE,EPSDE,0,FPFN5R,"N5R_FCPP")
                     CALL DEHIFT(RINT6,ES6,H0DE,EPSDE,0,FPFN6R,"N6R_FCPP")
                     PNFP_C=PNFPL
                     FCPPL=-RGAMH/(3.D0*RINT0)*(                &
                          (AMFP(NSA)*RGAMA**2)/(AMFD(NSB)*PNFP_C**2) &
                          *(3.D0*RINT4-THETA0L_C*RINT5)           &
                          +2.D0*(PTFP0(NSA)*PNFP_C)                  &
                          /(PTFD0(NSB)*RGAMA)*THETA0(NSSA)*RINT6 )   &
                          *RNFDL*1.D20
                  ELSEIF(PCRIT.gt.PMAX(NSSB))THEN
                     CALL DEHIFT(RINT0,ES0,H0DE,EPSDE,0,FPFN0R,"N0R_PMAX_DCPP")
!                     PNFP_C=PMAX(NSBA)
                     PNFP_C=PMAX(NSSB)
                     CALL DEFT  (RINT1,ES1,H0DE,EPSDE,0,FPFN1R,"N1R_PMAX_DCPP")
                     PNFP_C=PCRIT
                     CALL DEFT  (RINT7,ES7,H0DE,EPSDE,0,FPFN7R,"N7R_PMAX_DCPP")
                     CALL DEHIFT(RINT2,ES2,H0DE,EPSDE,0,FPFN2R,"N2R_PMAX_DCPP")
                     PNFP_C=PNFPL
                     DCPPL=RGAMH/(3.D0*RINT0)*(                 &
                          (AMFP(NSA)**2*PTFD0(NSB)**2*RGAMA**3) &
                          /(AMFD(NSB)**2*PTFP0(NSA)**2*PNFP_C**3) &
                          *(RINT1+RINT7)                        &
                          +(AMFD(NSB)*PTFP0(NSA))                &
                          /(AMFP(NSA)*PTFD0(NSB))*RINT2 )       &
                          *RNFDL*1.D20
!                     PNFP_C=PMAX(NSBA)
                     PNFP_C=PMAX(NSSB)
                     CALL DEFT(RINT4,ES4,H0DE,EPSDE,0,FPFN4R,"N4R_PMAX_FCPP")
                     CALL DEFT(RINT5,ES5,H0DE,EPSDE,0,FPFN5R,"N5R_PMAX_FCPP")
                     PNFP_C=PCRIT
                     CALL DEFT(RINT8,ES8,H0DE,EPSDE,0,FPFN9R,"N9R_PMAX_FCPP")
                     CALL DEFT(RINT9,ES9,H0DE,EPSDE,0,FPFN10R,"N10R_PMAX_FCPP")
                     CALL DEHIFT(RINT6,ES6,H0DE,EPSDE,0,FPFN6R,"N6R_PMAX_FCPP")
                     PNFP_C=PNFPL
                     FCPPL=-RGAMH/(3.D0*RINT0)*(                &
                          (AMFP(NSA)*RGAMA**2)/(AMFD(NSB)*PNFP_C**2) &
                          *(3.D0*(RINT4+RINT8)-THETA0L_C*(RINT5+RINT9)) &
                          +2.D0*(PTFP0(NSA)*PNFP_C)/(PTFD0(NSB)*RGAMA) &
                          *THETA0(NSSA)*RINT6 ) &
                          *RNFDL*1.D20
                  ENDIF
               ENDIF
               DO NTH=1,NTHMAX
                  DCPP2(NTH,NP,NR,NSB,NSA)=DCPP2(NTH,NP,NR,NSB,NSA)+DCPPL
                  FCPP2(NTH,NP,NR,NSB,NSA)=FCPP2(NTH,NP,NR,NSB,NSA)+FCPPL
               ENDDO
            ENDDO

            DO NP=NPSTARTW,NPENDWM
!               PNFPL=PM(NP,NSBA)
               PNFPL=PM(NP,NSSA)
               PNFP_C=PNFPL
               vtatb=(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))
               ptatb=PM(NP,NSSA)/RGAMA
!               PCRIT=SQRT(vtatb**2/(1.D0-THETA0L_C*vtatb**2*ptatb**2)) &
!                    *ptatb
               PCRIT=(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))*PM(NP,NSSA)
               IF(PCRIT.le.PMAX(NSSB))THEN
                  CALL DEHIFT(RINT0,ES0,H0DE,EPSDE,0,FPFN0R,"N0R_DCTT")
                  PNFP_C=PCRIT
                  CALL DEFT  (RINT1,ES1,H0DE,EPSDE,0,FPFN1R,"N1R_DCTT")
                  CALL DEHIFT(RINT2,ES2,H0DE,EPSDE,0,FPFN2R,"N2R_DCTT")
                  CALL DEFT  (RINT3,ES3,H0DE,EPSDE,0,FPFN3R,"N3R_DCTT")
                  PNFP_C=PNFPL
                  DCTTL=RGAMH/(3.D0*RINT0)*( &
                       1.5D0*RGAMA/PNFP_C*RINT3 &
                       -0.5D0*(AMFP(NSA)**2*PTFD0(NSB)**2*RGAMA**3) &
                       /(AMFD(NSB)**2*PTFP0(NSA)**2*PNFP_C**3)*RINT1 &
                       +(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))*RINT2 ) &
                       *RNFDL*1.D20
               ELSEIF(PCRIT.gt.PMAX(NSSB))THEN
                  CALL DEHIFT(RINT0,ES0,H0DE,EPSDE,0,FPFN0R,"N0R_PMAX_DCTT")
!                  PNFP_C=PMAX(NSBA)
                  PNFP_C=PMAX(NSSB)
                  CALL DEFT  (RINT1,ES1,H0DE,EPSDE,0,FPFN1R,"N1R_PMAX_DCTT")
                  CALL DEFT  (RINT3,ES3,H0DE,EPSDE,0,FPFN3R,"N3R_PMAX_DCTT")
                  PNFP_C=PNFPL*PTFP0(NSA)*AMFD(NSB)/(PTFD0(NSB)*AMFP(NSA))
                  CALL DEFT  (RINT4,ES1,H0DE,EPSDE,0,FPFN7R,"N7R_PMAX_DCTT")
                  CALL DEFT  (RINT5,ES3,H0DE,EPSDE,0,FPFN8R,"N8R_PMAX_DCTT")
                  CALL DEHIFT(RINT2,ES2,H0DE,EPSDE,0,FPFN2R,"N2R_PMAX_DCTT")
                  PNFP_C=PNFPL
                  DCTTL=RGAMH/(3.D0*RINT0)*( &
                       1.5D0*RGAMA/PNFP_C*(RINT3+RINT5) &
                       -0.5D0*(AMFP(NSA)**2*PTFD0(NSB)**2*RGAMA**3) &
                       /(AMFD(NSB)**2*PTFP0(NSA)**2*PNFP_C**3)*(RINT1+RINT4) &
                       +(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))*RINT2 ) &
                       *RNFDL*1.D20
              ENDIF
               
               DO NTH=1,NTHMAX+1
                  DCTT2(NTH,NP,NR,NSB,NSA)=DCTT2(NTH,NP,NR,NSB,NSA)+DCTTL
               ENDDO
            ENDDO
         ENDIF
!     ----- Relativistic -----
!
      ELSE
!         PMAXC=PMAX(NSBA)
         PMAXC=PMAX(NSSB)
         RTFD0L_C=RTFD0(NSB)
         RNFD0L_C=RNFD0(NSB)
         THETA0L_C=THETA0(NSSB)
         IF(MODEL_DISRUPT.eq.0)THEN
            IF(MODELC(NSSB).eq.1.or.MODELC(NSSB).eq.0)THEN ! constant T
               RNFDL_C=RNFD(NR,NSB) 
               RTFDL_C=RTFD(NR,NSB)
            ELSEIF(MODELC(NSSB).eq.2)THEN ! variable n, T
               RNFDL_C=RN_TEMP(NR,NS_NSB(NSB))
               RTFDL_C=RT_TEMP(NR,NS_NSB(NSB))
            END IF
         ELSEIF(MODEL_DISRUPT.ge.1)THEN
            IF(MODEL_IMPURITY.eq.0)THEN 
               RNFDL_C=RNFD(NR,NSB)
            ELSE
               RNFDL_C=RN_MGI(NR,NSB)
            END IF
            RGAMH=AEFP(NSA)**2*AEFD(NSB)**2*POST_LNLAM(NR,NSB,NSA)/(4.D0*PI*EPS0**2) &
                 *AMFP(NSA)/PTFP0(NSA)**3 
            RTFDL_C=RT_quench(NR)
         END IF
         THETAL_C =THETA0L_C*RTFDL_C/RTFD0L_C
         P_thermal=SQRT(RTFDL_C*1.D3*AEE*AMFD(NSB))
         v_thermal=SQRT(RTFDL_C*1.D3*AEE/AMFD(NSB) &
                        *(1.D0-2.5D0*THETAL_C+55.D0/8.D0*THETAL_C**2) )
         CALL DEHIFT(RINT0,ES0,H0DE,EPSDE,0,FPFN0R,"DEHIFT_0R")

         DO NP=NPSTART,NPENDWG
            IF(NP.EQ.1) THEN
!               PNFPL=PG(NP,NSBA)
               PNFPL=PG(NP,NSSA)
               PNFP_C=PNFPL
               RGAMA=SQRT(1.D0+PNFP_C**2*THETA0(NSSA))
               PCRIT=0.D0
               PNFP_C=PCRIT
               CALL DEHIFT(RINT2,ES2,H0DE,EPSDE,0,FPFN2R,"DEHIFT_2R")
               DCPPL=RGAMH/(3.D0*RINT0)*( &
                    (AMFD(NSB)*PTFP0(NSA))                     &
                    /(AMFP(NSA)*PTFD0(NSB))*RINT2 )             &
!                    *RNFD(NR,NSB)*1.D20
                    *RNFDL_C*1.D20 
               FCPPL=0.D0
            ELSE
!               PNFPL=PG(NP,NSBA)
               PNFPL=PG(NP,NSSA)
               PNFP_C=PNFPL
               RGAMA=SQRT(1.D0+PNFP_C**2*THETA0(NSSA))
               VFPL=PG(NP,NSSA)*PTFP0(NSA)/(AMFP(NSA)*RGAMA)
               vtatb=(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))
               ptatb=PG(NP,NSSA)/RGAMA
               PCRIT=(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))*PG(NP,NSSA)

               IF(VFPL.le.v_thermal*10)THEN
                  IF(PCRIT.le.PMAX(NSSB))THEN
!                     CALL DEHIFT(RINT0,ES0,H0DE,EPSDE,0,FPFN0R,"DEHIFT_0R_2")
                     PNFP_C=PCRIT
                     CALL DEFT  (RINT1,ES1,H0DE,EPSDE,0,FPFN1R,"DEFT1_le")
!                     IF(TIMEFP.ge.0.D-3) WRITE(*,'(3I5,A,E14.6)') NSA, NSB, NP, " RINT1=", RINT1
                     CALL DEHIFT(RINT2,ES2,H0DE,EPSDE,0,FPFN2R,"DEHIFT_2R_2")
!                     IF(TIMEFP.eq.5.D-3) WRITE(*,'(3I5,A,E14.6)') NSA, NSB, NP, " RINT2=", RINT2
                     PNFP_C=PNFPL
                     DCPPL=RGAMH/(3.D0*RINT0)*(                       &
                          (AMFP(NSA)**2*PTFD0(NSB)**2*RGAMA**3)       &
                          /(AMFD(NSB)**2*PTFP0(NSA)**2*PNFP_C**3)*RINT1 &
                          +(AMFD(NSB)*PTFP0(NSA))                     &
                          /(AMFP(NSA)*PTFD0(NSB))*RINT2 )             &
!                          *RNFD(NR,NSB)*1.D20
                          *RNFDL_C*1.D20
                     PNFP_C=PCRIT
                     CALL DEFT  (RINT4,ES4,H0DE,EPSDE,0,FPFN4R,"DEFT4_le")
!                     IF(TIMEFP.eq.5.D-3) WRITE(*,'(3I5,A,E14.6)') NSA, NSB, NP, "RINT4=", RINT4
                     CALL DEFT  (RINT5,ES5,H0DE,EPSDE,0,FPFN5R,"DEFT5_le")
!                     IF(TIMEFP.eq.5.D-3) WRITE(*,'(3I5,A,E14.6)') NSA, NSB, NP, "RINT5=", RINT5
                     CALL DEHIFT(RINT6,ES6,H0DE,EPSDE,0,FPFN6R,"DEHIFT_6R")
!                     IF(TIMEFP.eq.5.D-3) WRITE(*,'(3I5,A,E14.6)') NSA, NSB, NP, "RINT6=", RINT6
                     PNFP_C=PNFPL
                     FCPPL=-RGAMH/(3.D0*RINT0)*(                &
                          (AMFP(NSA)*RGAMA**2)/(AMFD(NSB)*PNFP_C**2) &
                          *(3.D0*RINT4-THETA0L_C*RINT5)           &
                          +2.D0*(PTFP0(NSA)*PNFP_C)                  &
                          /(PTFD0(NSB)*RGAMA)*THETA0(NSSA)*RINT6 )   &
!                          *RNFD(NR,NSB)*1.D20
                          *RNFDL_C*1.D20
!                     WRITE(6,'(A,2I5,1P10E14.6)') "low v e-e ", NP, NSB, THETAL_C, RINT0, RINT1, PCRIT, DCPPL, FCPPL
                  ELSEIF(PCRIT.gt.PMAX(NSSB))THEN
!                     CALL DEHIFT(RINT0,ES0,H0DE,EPSDE,0,FPFN0R,"DEHIFT_0R_3")
                     PNFP_C=PMAX(NSSB)
                     CALL DEFT  (RINT1,ES1,H0DE,EPSDE,0,FPFN1R,"DEFT1_gt")
                     PNFP_C=PCRIT
                     RINT7=0.D0
                     RINT2=0.D0
                     PNFP_C=PNFPL
                     DCPPL=RGAMH/(3.D0*RINT0)*(                 &
                          (AMFP(NSA)**2*PTFD0(NSB)**2*RGAMA**3) &
                          /(AMFD(NSB)**2*PTFP0(NSA)**2*PNFP_C**3) &
                          *(RINT1+RINT7)                        &
                          +(AMFD(NSB)*PTFP0(NSA))                &
                          /(AMFP(NSA)*PTFD0(NSB))*RINT2 )       &
!                          *RNFD(NR,NSB)*1.D20
                          *RNFDL_C*1.D20
                     PNFP_C=PMAX(NSSB)
                     CALL DEFT  (RINT4,ES4,H0DE,EPSDE,0,FPFN4R,"DEFT4_gt")
                     CALL DEFT  (RINT5,ES5,H0DE,EPSDE,0,FPFN5R,"DEFT5_gt")
                     PNFP_C=PCRIT
                     RINT8=0.D0
                     RINT9=0.D0
                     RINT6=0.D0
                     PNFP_C=PNFPL
                     FCPPL=-RGAMH/(3.D0*RINT0)*(                &
                          (AMFP(NSA)*RGAMA**2)/(AMFD(NSB)*PNFP_C**2) &
                          *(3.D0*(RINT4+RINT8)-THETA0L_C*(RINT5+RINT9)) &
                          +2.D0*(PTFP0(NSA)*PNFP_C)/(PTFD0(NSB)*RGAMA) &
                          *THETA0(NSSA)*RINT6 ) &
!                          *RNFD(NR,NSB)*1.D20
                          *RNFDL_C*1.D20
!                     WRITE(6,'(A,2I5,1P10E14.6)') "low v e-i ", NP, NSB, THETAL_C, RINT0, RINT1, PCRIT, DCPPL, FCPPL
                  ENDIF
               ELSE! high velocity limit 
                  DCPPL = RGAMH*(AMFP(NSA)/PTFP0(NSA))**2*(RGAMA/PG(NP,NSSA))**3 &
                       *v_thermal**2 &
!                       * RNFD(NR,NSB)*1.D20
                       * RNFDL_C*1.D20
                  FCPPL =-RGAMH*(RGAMA/PG(NP,NSSA))**2 &
                       *(1.D0-2.5D0*THETAL_C+55.D0/8.D0*THETAL_C**2) &
                       *AMFP(NSA)/AMFD(NSB) &
!                       * RNFD(NR,NSB)*1.D20 
                       * RNFDL_C*1.D20 

!                  WRITE(6,'(A,2I5,1P4E14.6)') "high v ", NP, NSB, THETAL_C, v_thermal, DCPPL, FCPPL
               END IF
            END IF
            DO NTH=1,NTHMAX
               DCPP2(NTH,NP,NR,NSB,NSA)=DCPP2(NTH,NP,NR,NSB,NSA)+DCPPL
               FCPP2(NTH,NP,NR,NSB,NSA)=FCPP2(NTH,NP,NR,NSB,NSA)+FCPPL
            ENDDO
         ENDDO

         DO NP=NPSTARTW,NPENDWM
            PNFP_C=PM(NP,NSSA)
            PNFP_C=PNFPL
            RGAMA=SQRT(1.D0+PNFP_C**2*THETA0(NSSA))
            VFPL=PM(NP,NSSA)*PTFP0(NSA)/(AMFP(NSA)*RGAMA)
            vtatb=(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))
            ptatb=PM(NP,NSSA)/RGAMA
            PCRIT=(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))*PM(NP,NSSA)

            IF(VFPL.le.v_thermal*10)THEN
               IF(PCRIT.le.PMAX(NSSB))THEN
                  CALL DEHIFT(RINT0,ES0,H0DE,EPSDE,0,FPFN0R,"DEHIFT_0R_4")
                  PNFP_C=PCRIT
                  CALL DEFT  (RINT1,ES1,H0DE,EPSDE,0,FPFN1R,"DEFT1_le_tt")
                  CALL DEHIFT(RINT2,ES2,H0DE,EPSDE,0,FPFN2R,"DEHIFT_2R_4")
                  CALL DEFT  (RINT3,ES3,H0DE,EPSDE,0,FPFN3R,"DEFT3_le_tt")
                  PNFP_C=PNFPL
                  DCTTL=RGAMH/(3.D0*RINT0)*( &
                       1.5D0*RGAMA/PNFP_C*RINT3 &
                       -0.5D0*(AMFP(NSA)**2*PTFD0(NSB)**2*RGAMA**3) &
                       /(AMFD(NSB)**2*PTFP0(NSA)**2*PNFP_C**3)*RINT1 &
                       +(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))*RINT2 ) &
!                       *RNFD(NR,NSB)*1.D20
                       *RNFDL_C*1.D20
               ELSEIF(PCRIT.gt.PMAX(NSSB))THEN
                  CALL DEHIFT(RINT0,ES0,H0DE,EPSDE,0,FPFN0R,"DEHIFT_0R_5")
                  PNFP_C=PMAX(NSSB)
                  CALL DEFT  (RINT1,ES1,H0DE,EPSDE,0,FPFN1R,"DEFT1_gt_tt")
                  CALL DEFT  (RINT3,ES3,H0DE,EPSDE,0,FPFN3R,"DEFT3_gt_tt")
                  PNFP_C=PNFPL*PTFP0(NSA)*AMFD(NSB)/(PTFD0(NSB)*AMFP(NSA))
                  RINT4=0.D0
                  RINT5=0.D0
                  RINT2=0.D0
                  PNFP_C=PNFPL
                  DCTTL=RGAMH/(3.D0*RINT0)*( &
                       1.5D0*RGAMA/PNFP_C*(RINT3+RINT5) &
                       -0.5D0*(AMFP(NSA)**2*PTFD0(NSB)**2*RGAMA**3) &
                       /(AMFD(NSB)**2*PTFP0(NSA)**2*PNFP_C**3)*(RINT1+RINT4) &
                       +(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))*RINT2 ) &
!                       *RNFD(NR,NSB)*1.D20
                       *RNFDL_C*1.D20
               ENDIF
            ELSE
               DCTTL=0.5D0*RGAMH*RGAMA/PM(NP,NSSA)* &
                    (1.D0-(v_thermal*AMFP(NSA)/PTFP0(NSA))**2*(RGAMA/PM(NP,NSSA))**2) &
!                    * RNFD(NR,NSB)*1.D20
                    * RNFDL_C*1.D20
            END IF
            DO NTH=1,NTHMAX+1
               DCTT2(NTH,NP,NR,NSB,NSA)=DCTT2(NTH,NP,NR,NSB,NSA)+DCTTL
            ENDDO
         ENDDO
      ENDIF

!      IF(TIMEFP.eq.5.D-3.and.NPSTART.eq.1)THEN
!         DO NP=NPSTART,NPEND
!            WRITE(*,'(A,3I4,4E14.6)') "TEST ", NSA, NSB, NP, PM(NP,NSSA), DCPP2(1,NP,1,NSB,NSA), DCTT2(1,NP,1,NSB,NSA), FCPP2(1,NP,1,NSB,NSA)
!         END DO
!      END IF

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
      REAL(rkind):: FACT
      REAL(RKIND):: DELH, sum, etal, psib, pcos, arg, x
      INTEGER:: NG, ITLB, ITUB, NSSB

! DCPP, FCPP
      DO NSB=1,NSBMAX
         DO NP=NPSTART,NPENDWG
            DO NTH=1,NTHMAX
               FACT=RLAMDA(NTH,NR)
               DCPP2(NTH,NP,NR,NSB,NSA) &
                    =FACT*DCPP2(NTH,NP,NR,NSB,NSA)!*2.D0
               FCPP2(NTH,NP,NR,NSB,NSA) &
                    =FACT*FCPP2(NTH,NP,NR,NSB,NSA)!*2.D0
            ENDDO
!            ITLB=ITL(NR)-1
!            ITUB=NTHMAX-ITLB+1
!            DCPP2(ITLB,NP,NR,NSB,NSA) &
!                     =RLAMDA(ITLB,NR)/4.D0 &
!                       *( DCPP2(ITLB-1,NP,NR,NSB,NSA) &
!                                   /RLAMDA(ITLB-1,NR) &
!                         +DCPP2(ITLB+1,NP,NR,NSB,NSA) &
!                                   /RLAMDA(ITLB+1,NR) &
!                         +DCPP2(ITUB-1,NP,NR,NSB,NSA) &
!                                   /RLAMDA(ITUB-1,NR) &
!                         +DCPP2(ITUB+1,NP,NR,NSB,NSA) &
!                                   /RLAMDA(ITUB+1,NR) )

!            FCPP2(ITLB,NP,NR,NSB,NSA) &
!                    =RLAMDA(ITLB,NR)/4.D0  &
!                      *( FCPP2(ITLB-1,NP,NR,NSB,NSA) &
!                                    /RLAMDA(ITLB-1,NR) &
!                        +FCPP2(ITLB+1,NP,NR,NSB,NSA) &
!                                    /RLAMDA(ITLB+1,NR) &
!                        +FCPP2(ITUB-1,NP,NR,NSB,NSA) &
!                                    /RLAMDA(ITUB-1,NR) &
!                        +FCPP2(ITUB+1,NP,NR,NSB,NSA) &
!                                    /RLAMDA(ITUB+1,NR) ) 
!            DCPP2(ITUB,NP,NR,NSB,NSA)=DCPP2(ITLB,NP,NR,NSB,NSA)
!            FCPP2(ITUB,NP,NR,NSB,NSA)=FCPP2(ITLB,NP,NR,NSB,NSA)
         END DO

! DCTT
         DO NP=NPSTARTW,NPENDWM
            DO NTH=1,NTHMAX+1
               IF(NTH.NE.NTHMAX/2+1) THEN 
                  DELH = 2.D0*ETAG(NTH,NR)/NAVMAX
                  SUM=0.D0
                  DO NG=1,NAVMAX
                     ETAL = DELH*(NG-0.5D0)
                     X=EPSRM2(NR)*COS(ETAL)
                     PSIB=(1.D0+EPSRM2(NR))/(1.D0+X)
                     ARG=1.D0-PSIB*SING(NTH)**2
                     PCOS = SQRT(ARG)
                     sum=sum + DCTT2(NTH,NP,NR,NSB,NSA)*PCOS/(PSIB*ABS(COSG(NTH))) 
                  END DO
                  DCTT2(NTH,NP,NR,NSB,NSA)=sum*DELH*Line_Element(NR)*A_chi0(NR)*2.D0
               ELSE
                  DCTT2(NTH,NP,NR,NSB,NSA)=0.D0
               END IF
            ENDDO

!            ITLB=ITL_judge(NR)+1
!            DO NTH=ITLB,NTHMAX/2
!               DCTT2(NTH,NP,NR,NSB,NSA)        &
!                    =(DCTT2(NTH,NP,NR,NSB,NSA) &
!                    +DCTT2(NTHMAX-NTH+2,NP,NR,NSB,NSA))*0.5D0
!               DCTT2(NTHMAX-NTH+2,NP,NR,NSB,NSA) &
!                    =DCTT2(NTH,NP,NR,NSB,NSA) 
!            END DO

         END DO ! NP
      END DO ! NSB

!     DPT, DTP, FTH
      DO NSB=1,NSBMAX
         NSSB=NS_NSB(NSB)
         IF(MODELC(NSSB).ge.3)THEN
         DO NP=NPSTART,NPENDWG
            DO NTH=1,NTHMAX
               DELH = 2.D0*ETAM(NTH,NR)/NAVMAX
               SUM=0.D0
               DO NG=1,NAVMAX
                  ETAL = DELH*(NG-0.5D0)
                  X=EPSRM2(NR)*COS(ETAL)
                  PSIB=(1.D0+EPSRM2(NR))/(1.D0+X)
                  sum=sum + DCPT2(NTH,NP,NR,NSB,NSA)/SQRT(PSIB)
               END DO
               DCPT2(NTH,NP,NR,NSB,NSA)=sum*DELH*Line_Element(NR)*A_chi0(NR)*2.D0
            ENDDO

            ITLB=ITL(NR)-1
            ITUB=NTHMAX-ITLB+1

            DCPT2(ITLB,NP,NR,NSB,NSA) &
                     =RLAMDA(ITLB,NR)/4.D0 &
                       *( DCPT2(ITLB-1,NP,NR,NSB,NSA) &
                                   /RLAMDA(ITLB-1,NR) &
                         +DCPT2(ITLB+1,NP,NR,NSB,NSA) &
                                   /RLAMDA(ITLB+1,NR) &
                         +DCPT2(ITUB-1,NP,NR,NSB,NSA) &
                                   /RLAMDA(ITUB-1,NR) &
                         +DCPT2(ITUB+1,NP,NR,NSB,NSA) &
                                   /RLAMDA(ITUB+1,NR) )

            DCPT2(ITUB,NP,NR,NSB,NSA)=DCPT2(ITLB,NP,NR,NSB,NSA)
         END DO

! DTP, FTH
         DO NP=NPSTARTW,NPENDWM
            DO NTH=1,NTHMAX+1
               DELH = 2.D0*ETAG(NTH,NR)/NAVMAX
               SUM=0.D0
               DO NG=1,NAVMAX
                  ETAL = DELH*(NG-0.5D0)
                  X=EPSRM2(NR)*COS(ETAL)
                  PSIB=(1.D0+EPSRM2(NR))/(1.D0+X)
                  sum=sum + 1.D0/SQRT(PSIB)
               END DO
               DCTP2(NTH,NP,NR,NSB,NSA)=sum*DCTP2(NTH,NP,NR,NSB,NSA)*DELH*Line_Element(NR)*A_chi0(NR)*2.D0
               FCTH2(NTH,NP,NR,NSB,NSA)=sum*FCTH2(NTH,NP,NR,NSB,NSA)*DELH*Line_Element(NR)*A_chi0(NR)*2.D0
            ENDDO

            ITLB=ITL_judge(NR)+1
            DO NTH=ITLB,NTHMAX/2
               DCTP2(NTH,NP,NR,NSB,NSA)        &
                    =(DCTP2(NTH,NP,NR,NSB,NSA) &
                    +DCTP2(NTHMAX-NTH+2,NP,NR,NSB,NSA))*0.5D0
               DCTP2(NTHMAX-NTH+2,NP,NR,NSB,NSA) &
                    =DCTP2(NTH,NP,NR,NSB,NSA) 

               FCTH2(NTH,NP,NR,NSB,NSA)        &
                    =(FCTH2(NTH,NP,NR,NSB,NSA) &
                    +FCTH2(NTHMAX-NTH+2,NP,NR,NSB,NSA))*0.5D0
               FCTH2(NTHMAX-NTH+2,NP,NR,NSB,NSA) &
                    =FCTH2(NTH,NP,NR,NSB,NSA) 
            END DO
         END DO ! NP
         END IF
      END DO ! NSB

      RETURN
      END SUBROUTINE FPCALC_LAV
!-----------------------------------------
      SUBROUTINE FPCALC_NLAV(NR,NSA)
!
      IMPLICIT NONE
      integer,intent(in):: NR, NSA
      integer:: NSB, NTH, NP, NG
      REAL(rkind):: DELH, ETAL, X, PSIB, PCOS, ARG
      REAL(rkind):: sum1, sum2, sum3, sum4, sum5, sum6
      INTEGER:: ISW_LAV
     
      ISW_LAV=0
! INTEGRATION OF BOUNCE AVERAGING
      DO NSB = 1,NSBMAX
         DO NP=NPSTART,NPENDWG
            DO NTH=1,NTHMAX
               IF(NTH.ne.ITL(NR).and.NTH.ne.ITU(NR))THEN
                  DELH=2.D0*ETAM(NTH,NR)/NAVMAX
                  sum1=0.D0
                  sum2=0.D0
                  sum3=0.D0
                  DO NG=1,NAVMAX
                     ETAL=DELH*(NG-0.5D0)
                     X=EPSRM2(NR)*COS(ETAL)
                     PSIB=(1.D0+EPSRM2(NR))/(1.D0+X)
                     PCOS=SQRT(1.D0-PSIB*SINM(NTH)**2)

                     sum1=sum1 + DCPP2(NTH,NP,NR,NSB,NSA)*ABS(COSM(NTH))/PCOS
                     sum2=sum2 + DCPT2(NTH,NP,NR,NSB,NSA)/SQRT(PSIB)
                     sum3=sum3 + FCPP2(NTH,NP,NR,NSB,NSA)*ABS(COSM(NTH))/PCOS
                  END DO ! END NAVMAX 
                  
                  DCPP2(NTH,NP,NR,NSB,NSA)=SUM1*DELH
                  DCPT2(NTH,NP,NR,NSB,NSA)=SUM2*DELH
                  FCPP2(NTH,NP,NR,NSB,NSA)=SUM3*DELH
               END IF
            END DO ! END NTH
!            IF(ISW_LAV.eq.1)THEN
!               INTH=0
!               DO NTH=ITL(NR),ITL(NR)+1
!                  INTH=INTH+1
!                  DELH=2.D0*ETAM(NTH,NR)/NAVMAX
!                  sum7=0.D0
!                  sum8=0.D0
!                  sum9=0.D0
!                  
!                  temp7 = DCPP2B(INTH,NP,NR,NSB,NSA)
!                  temp8 = FCPP2B(INTH,NP,NR,NSB,NSA)
!                  temp9 = DCPT2B(INTH,NP,NR,NSB,NSA)
!                  
!                  IF (COSM(NTH).GE.0.D0) THEN
!                     DO NG=1,NAVMAX
!                        ETAL=DELH*(NG-0.5D0)
!                        X=EPSRM(NR)*COS(ETAL)*RR
!                        PSIB=(1.D0+EPSRM(NR))/(1.D0+X/RR)
!                        PCOS=SQRT(1.D0-PSIB*SINM(NTH)**2)
!                        
!                        sum7=sum7 + temp7*COSM(NTH)/PCOS
!                        sum8=sum8 + temp8*COSM(NTH)/PCOS
!                        sum9=sum9 + temp9/SQRT(PSIB)
!                     END DO ! END NAVMAX 
!                  ELSE ! SIGN OF PCOS
!                     DO NG=1,NAVMAX
!                        ETAL=DELH*(NG-0.5D0)
!                        X=EPSRM(NR)*COS(ETAL)*RR
!                        PSIB=(1.D0+EPSRM(NR))/(1.D0+X/RR)
!                        PCOS=-SQRT(1.D0-PSIB*SINM(NTH)**2)
!                        
!                        sum7=sum7 + temp7*COSM(NTH)/PCOS
!                        sum8=sum8 + temp8*COSM(NTH)/PCOS
!                        sum9=sum9 + temp9/SQRT(PSIB)
!                     END DO ! END NAVMAX
!                  END IF
!                  
!                  DCPP2B(INTH,NP,NR,NSB,NSA)=SUM7*DELH
!                  FCPP2B(INTH,NP,NR,NSB,NSA)=SUM8*DELH
!                  DCPT2B(INTH,NP,NR,NSB,NSA)=SUM9*DELH
!               END DO ! END NTH
!               DO NTH=ITU(NR),ITU(NR)+1
!                  INTH=INTH+1
!                  DELH=2.D0*ETAM(NTH,NR)/NAVMAX
!                  sum7=0.D0
!                  sum8=0.D0
!                  sum9=0.D0
!                  
!                  temp7 = DCPP2B(INTH,NP,NR,NSB,NSA)
!                  temp8 = FCPP2B(INTH,NP,NR,NSB,NSA)
!                  temp9 = DCPT2B(INTH,NP,NR,NSB,NSA)
!                  
!                  IF (COSM(NTH).GE.0.D0) THEN
!                     DO NG=1,NAVMAX
!                        ETAL=DELH*(NG-0.5D0)
!                        X=EPSRM(NR)*COS(ETAL)*RR
!                        PSIB=(1.D0+EPSRM(NR))/(1.D0+X/RR)
!                        PCOS=SQRT(1.D0-PSIB*SINM(NTH)**2)
!                        
!                        sum7=sum7 + temp7*COSM(NTH)/PCOS
!                        sum8=sum8 + temp8*COSM(NTH)/PCOS
!                        sum9=sum9 + temp9/SQRT(PSIB)
!                     END DO ! END NAVMAX 
!                  ELSE ! SIGN OF PCOS
!                     DO NG=1,NAVMAX
!                        ETAL=DELH*(NG-0.5D0)
!                        X=EPSRM(NR)*COS(ETAL)*RR
!                        PSIB=(1.D0+EPSRM(NR))/(1.D0+X/RR)
!                        PCOS=-SQRT(1.D0-PSIB*SINM(NTH)**2)
!                        
!                        sum7=sum7 + temp7*COSM(NTH)/PCOS
!                        sum8=sum8 + temp8*COSM(NTH)/PCOS
!                        sum9=sum9 + temp9/SQRT(PSIB)
!                     END DO ! END NAVMAX
!                  END IF
!
!                  DCPP2B(INTH,NP,NR,NSB,NSA)=SUM7*DELH
!                  FCPP2B(INTH,NP,NR,NSB,NSA)=SUM8*DELH
!                  DCPT2B(INTH,NP,NR,NSB,NSA)=SUM9*DELH
!               END DO ! END NTH
!            END IF ! NEW GRID ISW
         END DO ! END NP
         DO NP=NPSTARTW,NPENDWM
            DO NTH=1,NTHMAX+1
               IF(NTH.NE.NTHMAX/2+1) THEN
                  DELH=2.D0*ETAG(NTH,NR)/NAVMAX
                  sum4=0.D0
                  sum5=0.D0
                  sum6=0.D0
                  DO NG=1,NAVMAX
                     ETAL=DELH*(NG-0.5D0)
                     X=EPSRM2(NR)*COS(ETAL)
                     PSIB=(1.D0+EPSRM2(NR))/(1.D0+X)
                     ARG=1.D0-PSIB*SING(NTH)**2
                     IF(ARG.GT.0.D0) THEN
                        IF (COSG(NTH).GE.0.D0) THEN
                           PCOS= SQRT(ARG)
                        ELSE
                           PCOS=-SQRT(ARG)
                        ENDIF
                     ELSE
                        PCOS=0.D0
                     ENDIF
                     sum4=sum4 + DCTP2(NTH,NP,NR,NSB,NSA)/SQRT(PSIB)
                     sum5=sum5 + DCTT2(NTH,NP,NR,NSB,NSA)*PCOS/(PSIB*COSG(NTH))
                     sum6=sum6 + FCTH2(NTH,NP,NR,NSB,NSA)/SQRT(PSIB)
                  END DO ! END NAVMAX
                  DCTP2(NTH,NP,NR,NSB,NSA)=sum4*DELH
                  DCTT2(NTH,NP,NR,NSB,NSA)=sum5*DELH
                  FCTH2(NTH,NP,NR,NSB,NSA)=sum6*DELH
               ELSE
                  DCTP2(NTH,NP,NR,NSB,NSA)=0.D0
                  DCTT2(NTH,NP,NR,NSB,NSA)=0.D0!?
                  FCTH2(NTH,NP,NR,NSB,NSA)=0.D0
               ENDIF ! END NTH!=pi/2
            END DO ! END NTH
         END DO ! END NP
      END DO ! END NSB
! END OF INTEGRATION

! BALANCE TRAPPED REGION for P direction
      DO NSB=1,NSBMAX
         DO NP=NPSTART,NPENDWG
            DO NTH=ITL(NR)+1,NTHMAX/2
               DCPP2(NTH,NP,NR,NSB,NSA) &
                    =(DCPP2(NTH,NP,NR,NSB,NSA) &
                    +DCPP2(NTHMAX-NTH+1,NP,NR,NSB,NSA))*0.5D0
               DCPP2(NTHMAX-NTH+1,NP,NR,NSB,NSA) &
                    =DCPP2(NTH,NP,NR,NSB,NSA)
            END DO ! END NTH
            IF(ISW_LAV.ne.1)THEN
               DCPP2(ITL(NR),NP,NR,NSB,NSA)=RLAMDA(ITL(NR),NR)*0.25D0     &
                    *( DCPP2(ITL(NR)-1,NP,NR,NSB,NSA)/RLAMDA(ITL(NR)-1,NR) &
                    +DCPP2(ITL(NR)+1,NP,NR,NSB,NSA)/RLAMDA(ITL(NR)+1,NR) &
                    +DCPP2(ITU(NR)-1,NP,NR,NSB,NSA)/RLAMDA(ITU(NR)-1,NR) &
                    +DCPP2(ITU(NR)+1,NP,NR,NSB,NSA)/RLAMDA(ITU(NR)+1,NR))
               DCPP2(ITU(NR),NP,NR,NSB,NSA)=DCPP2(ITL(NR),NP,NR,NSB,NSA)
            ELSE
!               DCPP2(ITL(NR),NP,NR,NSB,NSA)=RLAMDA(ITL(NR),NR)*0.25D0 &
!                    *( DCPP2B(1,NP,NR,NSB,NSA)/RLAMDA(ITL(NR)-1,NR) &
!                    +DCPP2B(2,NP,NR,NSB,NSA)/RLAMDA(ITL(NR)+1,NR)   &
!                    +DCPP2B(3,NP,NR,NSB,NSA)/RLAMDA(ITU(NR)-1,NR)   &
!                    +DCPP2B(4,NP,NR,NSB,NSA)/RLAMDA(ITU(NR)+1,NR))
!               DCPP2(ITU(NR),NP,NR,NSB,NSA)=DCPP2(ITL(NR),NP,NR,NSB,NSA)
            END IF
         END DO ! END NP
      END DO ! END NSB
      DO NSB=1,NSBMAX
         DO NP=NPSTART,NPENDWG
            DO NTH=ITL(NR)+1,NTHMAX/2
               FCPP2(NTH,NP,NR,NSB,NSA) &
                    =(FCPP2(NTH,NP,NR,NSB,NSA) &
                    +FCPP2(NTHMAX-NTH+1,NP,NR,NSB,NSA))*0.5D0
               FCPP2(NTHMAX-NTH+1,NP,NR,NSB,NSA) &
                    =FCPP2(NTH,NP,NR,NSB,NSA)
            END DO ! END NTH
            IF(ISW_LAV.ne.1)THEN
               FCPP2(ITL(NR),NP,NR,NSB,NSA)=RLAMDA(ITL(NR),NR)*0.25D0     &
                    *( FCPP2(ITL(NR)-1,NP,NR,NSB,NSA)/RLAMDA(ITL(NR)-1,NR) &
                    +FCPP2(ITL(NR)+1,NP,NR,NSB,NSA)/RLAMDA(ITL(NR)+1,NR) &
                    +FCPP2(ITU(NR)-1,NP,NR,NSB,NSA)/RLAMDA(ITU(NR)-1,NR) &
                    +FCPP2(ITU(NR)+1,NP,NR,NSB,NSA)/RLAMDA(ITU(NR)+1,NR))
               FCPP2(ITU(NR),NP,NR,NSB,NSA)=FCPP2(ITL(NR),NP,NR,NSB,NSA)
            ELSE
!               FCPP2(ITL(NR),NP,NR,NSB,NSA)=RLAMDA(ITL(NR),NR)*0.25D0     &
!                    *( FCPP2B(1,NP,NR,NSB,NSA)/RLAMDA(ITL(NR)-1,NR) &
!                    +FCPP2B(2,NP,NR,NSB,NSA)/RLAMDA(ITL(NR)+1,NR) &
!                    +FCPP2B(3,NP,NR,NSB,NSA)/RLAMDA(ITU(NR)-1,NR) &
!                    +FCPP2B(4,NP,NR,NSB,NSA)/RLAMDA(ITU(NR)+1,NR))
!               FCPP2(ITU(NR),NP,NR,NSB,NSA)=FCPP2(ITL(NR),NP,NR,NSB,NSA)
            END IF
         END DO ! END NP
      END DO ! END NSB
      DO NSB=1,NSBMAX
         DO NP=NPSTART,NPENDWG
            DO NTH=ITL(NR)+1,NTHMAX/2
               DCPT2(NTH,NP,NR,NSB,NSA) &
                    =(DCPT2(NTH,NP,NR,NSB,NSA) &
!                     +DCPT2(NTHMAX-NTH+1,NP,NR,NSB,NSA))*0.5D0 ! symmetry
                     -DCPT2(NTHMAX-NTH+1,NP,NR,NSB,NSA))*0.5D0 ! assymmetry
               DCPT2(NTHMAX-NTH+1,NP,NR,NSB,NSA) &
!                    =DCPT2(NTH,NP,NR,NSB,NSA) ! symmetry
                    =-DCPT2(NTH,NP,NR,NSB,NSA) ! assymmetry
            END DO ! END NTH
            IF(ISW_LAV.ne.1)THEN
               DCPT2(ITL(NR),NP,NR,NSB,NSA)=RLAMDA(ITL(NR),NR)*0.25D0     &
                    *( DCPT2(ITL(NR)-1,NP,NR,NSB,NSA)/RLAMDA(ITL(NR)-1,NR) &
                    +DCPT2(ITL(NR)+1,NP,NR,NSB,NSA)/RLAMDA(ITL(NR)+1,NR) &
!                    +DCPT2(ITU(NR)-1,NP,NR,NSB,NSA)/RLAMDA(ITU(NR)-1,NR) & ! symmetry
!                    +DCPT2(ITU(NR)+1,NP,NR,NSB,NSA)/RLAMDA(ITU(NR)+1,NR))
                    -DCPT2(ITU(NR)-1,NP,NR,NSB,NSA)/RLAMDA(ITU(NR)-1,NR) & ! assymetry
                    -DCPT2(ITU(NR)+1,NP,NR,NSB,NSA)/RLAMDA(ITU(NR)+1,NR))
!               DCPT2(ITU(NR),NP,NR,NSB,NSA)=DCPT2(ITL(NR),NP,NR,NSB,NSA) ! symmetry
               DCPT2(ITU(NR),NP,NR,NSB,NSA)=-DCPT2(ITL(NR),NP,NR,NSB,NSA) ! assymmetry
            ELSE
!               DCPT2(ITL(NR),NP,NR,NSB,NSA)=RLAMDA(ITL(NR),NR)*0.25D0     &
!                    *( DCPT2B(1,NP,NR,NSB,NSA)/RLAMDA(ITL(NR)-1,NR) &
!                    +DCPT2B(2,NP,NR,NSB,NSA)/RLAMDA(ITL(NR)+1,NR) &
!                    +DCPT2B(3,NP,NR,NSB,NSA)/RLAMDA(ITU(NR)-1,NR) &
!                    +DCPT2B(4,NP,NR,NSB,NSA)/RLAMDA(ITU(NR)+1,NR))
!               DCPT2(ITU(NR),NP,NR,NSB,NSA)=-DCPT2(ITL(NR),NP,NR,NSB,NSA)
            END IF
         END DO ! END NP
      END DO ! END NSB
! END OF BALANCE TRAPPED REGION for P direction

! BALANCE TRAPPED REGION for THETA direction
      DO NSB=1,NSBMAX
         DO NP=NPSTARTW,NPENDWM
            DO NTH=ITL(NR)+1,NTHMAX/2
               DCTT2(NTH,NP,NR,NSB,NSA)        &
                    =(DCTT2(NTH,NP,NR,NSB,NSA) &
                    +DCTT2(NTHMAX-NTH+2,NP,NR,NSB,NSA))*0.5D0
               DCTT2(NTHMAX-NTH+2,NP,NR,NSB,NSA) &
                    =DCTT2(NTH,NP,NR,NSB,NSA)
            END DO
         END DO
      END DO
      DO NSB=1,NSBMAX
         DO NP=NPSTARTW,NPENDWM
            DO NTH=ITL(NR)+1,NTHMAX/2
               FCTH2(NTH,NP,NR,NSB,NSA)        &
                    =(FCTH2(NTH,NP,NR,NSB,NSA) &
!                     +FCTH2(NTHMAX-NTH+2,NP,NR,NSB,NSA))*0.5D0 ! symmetry
                     -FCTH2(NTHMAX-NTH+2,NP,NR,NSB,NSA))*0.5D0 ! a
               FCTH2(NTHMAX-NTH+2,NP,NR,NSB,NSA) &
!                    =FCTH2(NTH,NP,NR,NSB,NSA) ! symmetry
                    =-FCTH2(NTH,NP,NR,NSB,NSA)
            END DO
         END DO
      END DO
      DO NSB=1,NSBMAX
         DO NP=NPSTARTW,NPENDWM
            DO NTH=ITL(NR)+1,NTHMAX/2
               DCTP2(NTH,NP,NR,NSB,NSA)        &
                    =(DCTP2(NTH,NP,NR,NSB,NSA) &
!                     +DCTP2(NTHMAX-NTH+2,NP,NR,NSB,NSA))*0.5D0 ! symmetry
                     -DCTP2(NTHMAX-NTH+2,NP,NR,NSB,NSA))*0.5D0
               DCTP2(NTHMAX-NTH+2,NP,NR,NSB,NSA) &
!                    =DCTP2(NTH,NP,NR,NSB,NSA) ! symmetry
                    =-DCTP2(NTH,NP,NR,NSB,NSA)
            END DO
         END DO
      END DO
! END OF BALANCE TRAPPED REGION for THETA direction
      
      RETURN
      END SUBROUTINE FPCALC_NLAV
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      END MODULE fpcalc
