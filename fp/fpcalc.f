C     $Id$
C
C ******************************************************************
C
C     CALCULATION OF COLLISINAL DIFUSION COEFFICIENTS: DC AND FC
C
C ******************************************************************
C
      SUBROUTINE FPCALC(NSA)
C
      INCLUDE 'fpcomm.inc'

      DO NR=1,NRMAX
         DO NP=1,NPMAX+1
         DO NTH=1,NTHMAX
            DCPP(NTH,NP,NR,NSA)=0.D0
            DCPT(NTH,NP,NR,NSA)=0.D0
            FCPP(NTH,NP,NR,NSA)=0.D0
         ENDDO
         ENDDO
         DO NP=1,NPMAX
         DO NTH=1,NTHMAX+1
            DCTP(NTH,NP,NR,NSA)=0.D0
            DCTT(NTH,NP,NR,NSA)=0.D0
            FCTH(NTH,NP,NR,NSA)=0.D0
         ENDDO
         ENDDO
      ENDDO

      DO NSB=1,NSBMAX
      DO NR=1,NRMAX
         DO NP=1,NPMAX+1
         DO NTH=1,NTHMAX
            DCPP2(NTH,NP,NR,NSB,NSA)=0.D0
            DCPT2(NTH,NP,NR,NSB,NSA)=0.D0
            FCPP2(NTH,NP,NR,NSB,NSA)=0.D0
         ENDDO
         ENDDO
         DO NP=1,NPMAX
         DO NTH=1,NTHMAX+1
            DCTP2(NTH,NP,NR,NSB,NSA)=0.D0
            DCTT2(NTH,NP,NR,NSB,NSA)=0.D0
            FCTH2(NTH,NP,NR,NSB,NSA)=0.D0
         ENDDO
         ENDDO
      ENDDO
      ENDDO
C
      DO NR=1,NRMAX
         DO NSB=1,NSBMAX
C
c            if(nr.eq.1) write(6,'(A,I8,1P2E12.4)') 
c     &           ' NSB,RN,RNFD=',NSB,RN(NSB),RNFD(NR,NSB)
C
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
            ELSEIF(MODELC.EQ.-1) THEN
               IF(NS_NSB(NSB).EQ.NS_NSA(NSA)) THEN
                  CALL FPCALC_L(NR,NSB,NSA)
               ENDIF
            ELSEIF(MODELC.EQ.-2) THEN
               IF(NS_NSB(NSB).EQ.NS_NSA(NSA)) THEN
                  CALL FPCALC_NL(NR,NSB,NSA)
               ENDIF
            ENDIF
         ENDDO

C        ----- Simple electron-ion collision term using ZEFF -----

         IF(MODELC.LT.0) THEN
            IF(NS_NSA(NSA).EQ.1) THEN
               DO NSB=1,NSBMAX
                  IF(NS_NSB(NSB).EQ.2) THEN
c                     RNNL=RNFP(NR,NSA)/RNFP0(NSA)
c                     RNUDL=RNUD(NR,NSB,NSA)*RNNL
      RGAMH=RNUD(NR,NSA,NSA)*SQRT(2.D0)*VTFD(NR,NSA)*AMFP(NSA)
     &        /(RNFP0(NSA)*PTFP0(NSA)*1.D20)
                   RGAMH2=RGAMH*RNFD(NR,NSA)*1.D20*PTFP0(NSA)/AMFP(NSA)
               RTE=(RTPR(1)+RTPP(1)*2.D0)/3.D0
              rZI = -PZ(2)/PZ(1)/(14.9D0-0.5D0*LOG(RN(1))+LOG(RTE))*
     &              (15.2D0-0.5D0*LOG(RN(1))+LOG(RTE))
c                     RZI=-PZ(2)/PZ(1)*RNUD(NR,NSB,NSA)/RNUD(NR,NSA,NSA)
                     DO NP=1,NPMAX
                        PFPL=PM(NP)*PTFP0(NSA)
                        VFPL=PFPL/AMFP(NSA)
                        U=VFPL/VTFP0(NSA)
                        DCTTL=0.5D0*RZI*RGAMH2/VFPL
c                        write(*,*)DCTTL
                        DO NTH=1,NTHMAX+1
                           DCTT2(NTH,NP,NR,NSB,NSA)
     &                          =DCTT2(NTH,NP,NR,NSB,NSA)+DCTTL
                        ENDDO
                     ENDDO
c         write(*,*) "test",DCTT2(2,2,1,NSB,NSA),NSA,NSB
                  ENDIF
               ENDDO
            ENDIF
         ENDIF

C     ----- bounce average -----

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

c     sum up coefficients by species

         DO NSB=1,NSBMAX
            DO NP=1,NPMAX+1
               DO NTH=1,NTHMAX
                  DCPP(NTH,NP,NR,NSA)=DCPP(NTH,NP,NR,NSA)
     &                               +DCPP2(NTH,NP,NR,NSB,NSA)
                  DCPT(NTH,NP,NR,NSA)=DCPT(NTH,NP,NR,NSA)
     &                               +DCPT2(NTH,NP,NR,NSB,NSA)
                  FCPP(NTH,NP,NR,NSA)=FCPP(NTH,NP,NR,NSA)
     &                               +FCPP2(NTH,NP,NR,NSB,NSA)
               END DO
            END DO
            DO NP=1,NPMAX
               DO NTH=1,NTHMAX+1
                  DCTP(NTH,NP,NR,NSA)=DCTP(NTH,NP,NR,NSA)
     &                               +DCTP2(NTH,NP,NR,NSB,NSA)
                  DCTT(NTH,NP,NR,NSA)=DCTT(NTH,NP,NR,NSA)
     &                               +DCTT2(NTH,NP,NR,NSB,NSA)
                  FCTH(NTH,NP,NR,NSA)=FCTH(NTH,NP,NR,NSA)
     &                               +FCTH2(NTH,NP,NR,NSB,NSA)
               END DO
            END DO
         END DO
      ENDDO

C
      IF(MOD(IDBGFP/8,2).EQ.1) THEN
C
C     +++ plot of D_coll +++
C
         CALL PAGES
         CALL FPGRFA(1,DCPP,PG,1,'@DCPP(P)@',NTHM,NPM,NRM,
     &                                       NTHMAX,NPMAX,NRMAX,NSA)
         CALL FPGRFA(2,DCTT,PM,2,'@DCTT(P)@',NTHM,NPM,NRM,
     &                                       NTHMAX,NPMAX,NRMAX,NSA)
         CALL FPGRFA(3,FCPP,PG,1,'@FCPP(P)@',NTHM,NPM,NRM,
     &                                       NTHMAX,NPMAX,NRMAX,NSA)
         CALL PAGEE
      ENDIF
C
      IF(MOD(IDBGFP/16,2).EQ.1) THEN
C
C     +++ plot of D_coll +++
C
         CALL PAGES
         CALL FPGRFB(1,DCPP,THM,1,'@DCPP(TH)@',NTHM,NPM,NRM,
     &                                         NTHMAX,NPMAX,NRMAX,NSA)
         CALL FPGRFB(2,DCTT,THG,2,'@DCTT(TH)@',NTHM,NPM,NRM,
     &                                         NTHMAX,NPMAX,NRMAX,NSA)
         CALL FPGRFB(3,FCPP,THM,1,'@FCPP(TH)@',NTHM,NPM,NRM,
     &                                         NTHMAX,NPMAX,NRMAX,NSA)
         CALL PAGEE
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE FPGRFA(ID,DATA,P,IND,TITLE,NTHM,NPM,NRM,
     &                                      NTHMAX,NPMAX,NRMAX,NSA)
C
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
      end
C
      SUBROUTINE FPGRFB(ID,DATA,TH,IND,TITLE,NTHM,NPM,NRM,
     &                                      NTHMAX,NPMAX,NRMAX,NSA)
C
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
      end
C
C ************************************************************
C
C       CALCULATION OF LINEAR COLLISIONAL OPERATOR
C
C ************************************************************
C
      SUBROUTINE FPCALC_L(NR,NSB,NSA)
C
      INCLUDE 'fpcomm.inc'
      COMMON /FPFNV1/ PNFP,TMC2FD,TMC2FD0,PMAXC
      EXTERNAL FPFN0R,FPFN1R,FPFN2R,FPFN3R,FPFN4R,FPFN5R,FPFN6R
     &     ,FPFN7R,FPFN8R,FPFN9R,FPFN10R

C     ------ define --------
      RNNL=RNFD(NR,NSB)/RNFP0(NSA)
      RNUFL=RNUF(NR,NSB,NSA)*RNNL
      RNUDL=RNUD(NR,NSB,NSA)*RNNL

      RGAMH=RNUD(NR,NSB,NSA)*SQRT(2.D0)*VTFD(NR,NSB)*AMFP(NSA)
     &        /(RNFP0(NSA)*PTFP0(NSA)*1.D20)
C
C     ----- Non-Relativistic -----
C
      IF(MODELR.EQ.0) THEN 
         DO NP=1,NPMAX+1
            IF(NP.EQ.1) THEN
               DCPPL=RGAMH*RNFD(NR,NSB)*1.D20*(2.D0/(3.D0*SQRT(PI))) 
     &                    *(VTFP0(NSA)/(SQRT(2.D0)*VTFD(NR,NSB))) 
               FCPPL=0.D0
            ELSE
               PFPL=PG(NP)*PTFP0(NSA)
               VFPL=PFPL/AMFP(NSA)
               V=VFPL/VTFP0(NSA)
               U=VFPL/(SQRT(2.D0)*VTFD(NR,NSB))
               DCPPL= 0.5D0*RNUDL/U   *(ERF0(U)/U**2-ERF1(U)/U)
               FCPPL=-      RNUFL/U**2*(ERF0(U)-U*ERF1(U))
            ENDIF
            DO NTH=1,NTHMAX
               DCPP2(NTH,NP,NR,NSB,NSA)=DCPP2(NTH,NP,NR,NSB,NSA)+DCPPL
               FCPP2(NTH,NP,NR,NSB,NSA)=FCPP2(NTH,NP,NR,NSB,NSA)+FCPPL
            ENDDO
         ENDDO

         DO NP=1,NPMAX
            PFPL=PM(NP)*PTFP0(NSA)
            VFPL=PFPL/AMFP(NSA)
            V=VFPL/VTFP0(NSA)
            U=VFPL/(SQRT(2.D0)*VTFD(NR,NSB))
            DCTTL= 0.25D0*RNUDL/U
     &                   *((2.D0-1.D0/U**2)*ERF0(U)+ERF1(U)/U)
C
            DO NTH=1,NTHMAX+1
               DCTT2(NTH,NP,NR,NSB,NSA)=DCTT2(NTH,NP,NR,NSB,NSA)+DCTTL
            ENDDO
         ENDDO
C     ----- Relativistic -----
C
      ELSE
         PMAXC=PMAX
         DO NP=1,NPMAX+1
            IF(NP.EQ.1) THEN
              DCPPL=RGAMH*RNFD(NR,NSB)*1.D20*(2.D0/(3.D0*SQRT(PI)))
     &                    *(VTFP0(NSA)/(SQRT(2.D0)*VTFD(NR,NSB)))
              FCPPL=0.D0
            ELSE
               PFPL=PG(NP)*PTFP0(NSA)
               VFPL=PFPL/SQRT(AMFP(NSA)**2+PFPL**2/VC**2)
               PNFPL=PG(NP)
               PNFP=PNFPL
               TMC2FD =(PTFD(NR,NSB)/(AMFD(NSB)*VC))**2
               TMC2FD0=(PTFD0(NSB)/(AMFD(NSB)*VC))**2
               TMC2FP0=(PTFP0(NSA)/(AMFP(NSA)*VC))**2
               RGAMA=SQRT(1.D0+PNFP**2*TMC2FP0)
               vtatb=(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))
               ptatb=PG(NP)/RGAMA
               PCRIT=SQRT(vtatb**2/(1.D0-TMC2FD0*vtatb**2*ptatb**2))
     &              *ptatb
              IF(PCRIT.le.PMAX)THEN
                  CALL DEHIFT(RINT0,ES0,H0DE,EPSDE,0,FPFN0R)
                  PNFP=PCRIT
                  CALL DEFT  (RINT1,ES1,H0DE,EPSDE,0,FPFN1R)
                  CALL DEHIFT(RINT2,ES2,H0DE,EPSDE,0,FPFN2R)
                  PNFP=PNFPL
                  DCPPL=RGAMH/(3.D0*RINT0)*(
     &                 (AMFP(NSA)**2*PTFD0(NSB)**2*RGAMA**3)
     &                 /(AMFD(NSB)**2*PTFP0(NSA)**2*PNFP**3)*RINT1
     &                +(AMFD(NSB)*PTFP0(NSA))
     &                 /(AMFP(NSA)*PTFD0(NSB))*RINT2 )
     &             *RNFD(NR,NSB)*1.D20
                  PNFP=PCRIT
                  CALL DEFT  (RINT4,ES4,H0DE,EPSDE,0,FPFN4R)
                  CALL DEFT  (RINT5,ES5,H0DE,EPSDE,0,FPFN5R)
                  CALL DEHIFT(RINT6,ES6,H0DE,EPSDE,0,FPFN6R)
                  PNFP=PNFPL
                  FCPPL=-RGAMH/(3.D0*RINT0)*(
     &              (AMFP(NSA)*RGAMA**2)/(AMFD(NSB)*PNFP**2)
     &                 *(3.D0*RINT4-TMC2FD0*RINT5)
     &              +2.D0*(PTFP0(NSA)*PNFP)
     &                 /(PTFD0(NSB)*RGAMA)*TMC2FP0*RINT6 )
     &             *RNFD(NR,NSB)*1.D20
               ELSEIF(PCRIT.gt.PMAX)THEN
                  CALL DEHIFT(RINT0,ES0,H0DE,EPSDE,0,FPFN0R)
                  PNFP=PMAX
                  CALL DEFT  (RINT1,ES1,H0DE,EPSDE,0,FPFN1R)
                  PNFP=PCRIT
                  CALL DEFT  (RINT7,ES7,H0DE,EPSDE,0,FPFN7R)
                  CALL DEHIFT(RINT2,ES2,H0DE,EPSDE,0,FPFN2R)
                  PNFP=PNFPL
                  DCPPL=RGAMH/(3.D0*RINT0)*( 
     &                 (AMFP(NSA)**2*PTFD0(NSB)**2*RGAMA**3)
     &                 /(AMFD(NSB)**2*PTFP0(NSA)**2*PNFP**3)
     &                 *(RINT1+RINT7)
     &                +(AMFD(NSB)*PTFP0(NSA))
     &                 /(AMFP(NSA)*PTFD0(NSB))*RINT2 )
     &                 *RNFD(NR,NSB)*1.D20
                  PNFP=PMAX
                  CALL DEFT  (RINT4,ES4,H0DE,EPSDE,0,FPFN4R)
                  CALL DEFT  (RINT5,ES5,H0DE,EPSDE,0,FPFN5R)
                  PNFP=PCRIT
                  CALL DEFT  (RINT8,ES8,H0DE,EPSDE,0,FPFN9R)
                  CALL DEFT  (RINT9,ES9,H0DE,EPSDE,0,FPFN10R)
                  CALL DEHIFT(RINT6,ES6,H0DE,EPSDE,0,FPFN6R)
                  PNFP=PNFPL
                  FCPPL=-RGAMH/(3.D0*RINT0)*(
     &                 (AMFP(NSA)*RGAMA**2)/(AMFD(NSB)*PNFP**2)
     &                 *(3.D0*(RINT4+RINT8)-TMC2FD0*(RINT5+RINT9))
     &                 +2.D0*(PTFP0(NSA)*PNFP)/(PTFD0(NSB)*RGAMA)
     &                 *TMC2FP0*RINT6 )
     &             *RNFD(NR,NSB)*1.D20
               ENDIF
            ENDIF
            DO NTH=1,NTHMAX
               DCPP2(NTH,NP,NR,NSB,NSA)=DCPP2(NTH,NP,NR,NSB,NSA)+DCPPL
               FCPP2(NTH,NP,NR,NSB,NSA)=FCPP2(NTH,NP,NR,NSB,NSA)+FCPPL
            ENDDO
         ENDDO

         DO NP=1,NPMAX
            PFPL=PM(NP)*PTFP0(NSA)
            VFPL=PFPL/SQRT(AMFP(NSA)**2+PTFP(NR,NSA)**2/VC**2)
            PNFPL=PM(NP)
            PNFP=PNFPL
            TMC2FD =PTFD(NR,NSB)**2/(AMFD(NSB)*VC)**2
            TMC2FD0=PTFD0(NSB)**2/(AMFD(NSB)*VC)**2
            TMC2FP0=PTFP0(NSA)**2/(AMFP(NSA)*VC)**2
            RGAMA=SQRT(1.D0+PNFP**2*TMC2FP0)
            vtatb=(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))
            ptatb=PM(NP)/RGAMA
            PCRIT=SQRT(vtatb**2/(1.D0-TMC2FD0*vtatb**2*ptatb**2))
     &           *ptatb
            IF(PCRIT.le.PMAX)THEN
               CALL DEHIFT(RINT0,ES0,H0DE,EPSDE,0,FPFN0R)
               PNFP=PNFPL*PTFP0(NSA)*AMFD(NSB)/(PTFD0(NSB)*AMFP(NSA))
               CALL DEFT  (RINT1,ES1,H0DE,EPSDE,0,FPFN1R)
               CALL DEHIFT(RINT2,ES2,H0DE,EPSDE,0,FPFN2R)
               CALL DEFT  (RINT3,ES3,H0DE,EPSDE,0,FPFN3R)
               PNFP=PNFPL
               DCTTL=RGAMH/(3.D0*RINT0)*(
     &            1.5D0*RGAMA/PNFP*RINT3
     &           -0.5D0*(AMFP(NSA)**2*PTFD0(NSB)**2*RGAMA**3)
     &           /(AMFD(NSB)**2*PTFP0(NSA)**2*PNFP**3)*RINT1
     &           +(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))*RINT2 )
     &          *RNFD(NR,NSB)*1.D20
            ELSEIF(PCRIT.gt.PMAX)THEN
               CALL DEHIFT(RINT0,ES0,H0DE,EPSDE,0,FPFN0R)
               PNFP=PMAX
               CALL DEFT  (RINT1,ES1,H0DE,EPSDE,0,FPFN1R)
               CALL DEFT  (RINT3,ES3,H0DE,EPSDE,0,FPFN3R)
               PNFP=PNFPL*PTFP0(NSA)*AMFD(NSB)/(PTFD0(NSB)*AMFP(NSA))
               CALL DEFT  (RINT4,ES1,H0DE,EPSDE,0,FPFN7R)
               CALL DEFT  (RINT5,ES3,H0DE,EPSDE,0,FPFN8R)
               CALL DEHIFT(RINT2,ES2,H0DE,EPSDE,0,FPFN2R)
               PNFP=PNFPL
               DCTTL=RGAMH/(3.D0*RINT0)*(
     &            1.5D0*RGAMA/PNFP*(RINT3+RINT5)
     &           -0.5D0*(AMFP(NSA)**2*PTFD0(NSB)**2*RGAMA**3)
     &           /(AMFD(NSB)**2*PTFP0(NSA)**2*PNFP**3)*(RINT1+RINT4)
     &           +(AMFD(NSB)*PTFP0(NSA))/(AMFP(NSA)*PTFD0(NSB))*RINT2 )
     &          *RNFD(NR,NSB)*1.D20
            ENDIF

            DO NTH=1,NTHMAX+1
               DCTT2(NTH,NP,NR,NSB,NSA)=DCTT2(NTH,NP,NR,NSB,NSA)+DCTTL
            ENDDO
         ENDDO
      ENDIF

      RETURN
      END
C
C ************************************************************
C
C       BOUNCE AVERAGE OF LINEAR COLLISIONAL OPERATOR
C
C ************************************************************
C
      SUBROUTINE FPCALC_LAV(NR,NSA)
C
      INCLUDE 'fpcomm.inc'
C
      DO NTH=1,NTHMAX
         FACT=RLAMDA(NTH,NR)
         DO NP=1,NPMAX+1
            DO NSB=1,NSBMAX
               DCPP2(NTH,NP,NR,NSB,NSA)
     &              =FACT*DCPP2(NTH,NP,NR,NSB,NSA)
               FCPP2(NTH,NP,NR,NSB,NSA)
     &              =FACT*FCPP2(NTH,NP,NR,NSB,NSA)
            END DO
         ENDDO
C
         FACT=RLAMDC(NTH,NR)
         DO NP=1,NPMAX
            DO NSB=1,NSBMAX
               DCTT2(NTH,NP,NR,NSB,NSA)=FACT*DCTT2(NTH,NP,NR,NSB,NSA)
            END DO
         ENDDO

         DO NP=1,NPMAX+1
            DO NSB=1,NSBMAX
               DCPP2(ITL(NR),NP,NR,NSB,NSA)
     &              =RLAMDA(ITL(NR),NR)/4.D0
     &                *( DCPP2(ITL(NR)-1,NP,NR,NSB,NSA)
     &                                          /RLAMDA(ITL(NR)-1,NR)
     &                  +DCPP2(ITL(NR)+1,NP,NR,NSB,NSA)
     &                                          /RLAMDA(ITL(NR)+1,NR)
     &                  +DCPP2(ITU(NR)-1,NP,NR,NSB,NSA)
     &                                          /RLAMDA(ITU(NR)-1,NR)
     &                  +DCPP2(ITU(NR)+1,NP,NR,NSB,NSA)
     &                                          /RLAMDA(ITU(NR)+1,NR))

               FCPP2(ITL(NR),NP,NR,NSB,NSA)
     &              =RLAMDA(ITL(NR),NR)/4.D0
     &                *( FCPP2(ITL(NR)-1,NP,NR,NSB,NSA)
     &                                          /RLAMDA(ITL(NR)-1,NR)
     &                  +FCPP2(ITL(NR)+1,NP,NR,NSB,NSA)
     &                                          /RLAMDA(ITL(NR)+1,NR)
     &                  +FCPP2(ITU(NR)-1,NP,NR,NSB,NSA)
     &                                          /RLAMDA(ITU(NR)-1,NR)
     &                  +FCPP2(ITU(NR)+1,NP,NR,NSB,NSA)
     &                                          /RLAMDA(ITU(NR)+1,NR))
               DCPP2(ITU(NR),NP,NR,NSB,NSA)=DCPP2(ITL(NR),NP,NR,NSB,NSA)
               FCPP2(ITU(NR),NP,NR,NSB,NSA)=FCPP2(ITL(NR),NP,NR,NSB,NSA)
            END DO
         END DO
      END DO

      RETURN
      END
C
C ***************************************************************
C
C                       SET OF INTEGRAND
C
C ***************************************************************
C
      REAL*8 FUNCTION FPFN0R(X)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-F,H,O-Z)
      COMMON /FPFNV1/ PNFP,TMC2FD,TMC2FD0,PMAXC
C
      PN=X
      FPFN0R=PN**2*FPRMXW(PN)
C
      RETURN
      END
C---------------------------------------------------------------
      REAL*8 FUNCTION FPFN1R(X,XM,XP)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-F,H,O-Z)
      COMMON /FPFNV1/ PNFP,TMC2FD,TMC2FD0,PMAXC
C
      A=0.5D0*PNFP
      PN=A*XP
      B=PN**4/(1.D0+PN**2*TMC2FD0)
      FPFN1R=A*B*FPRMXW(PN)
C
      RETURN
      END
C
C ===============================================================
C
      REAL*8 FUNCTION FPFN2R(X)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-F,H,O-Z)
      COMMON /FPFNV1/ PNFP,TMC2FD,TMC2FD0,PMAXC
C
      A=1.D0
      PN=A*(X+PNFP)
      B=PN*SQRT(1.D0+PN**2*TMC2FD0)
      FPFN2R=A*B*FPRMXW(PN)
C
      RETURN
      END
C
C ==============================================================
C
      REAL*8 FUNCTION FPFN3R(X,XM,XP)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-F,H,O-Z)
      COMMON /FPFNV1/ PNFP,TMC2FD,TMC2FD0,PMAXC
C
      A=0.5D0*PNFP
      PN=A*XP
      FPFN3R=A*PN**2*FPRMXW(PN)
C
      RETURN
      END
C
C ===============================================================
C
      REAL*8 FUNCTION FPFN4R(X,XM,XP)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-F,H,O-Z)
      COMMON /FPFNV1/ PNFP,TMC2FD,TMC2FD0,PMAXC
C
      A=0.5D0*PNFP
      PN=A*XP
      B=PN**2/SQRT(1.D0+PN**2*TMC2FD0)
      FPFN4R=A*B*FPRMXW(PN)
C
      RETURN
      END
C
C ===============================================================
C
      REAL*8 FUNCTION FPFN5R(X,XM,XP)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-F,H,O-Z)
      COMMON /FPFNV1/ PNFP,TMC2FD,TMC2FD0,PMAXC
C
      A=0.5D0*PNFP
      PN=A*XP
      B=PN**4/(SQRT(1.D0+PN**2*TMC2FD0))**3
      FPFN5R=A*B*FPRMXW(PN)
C
      RETURN
      END
C
C ===============================================================
C
      REAL*8 FUNCTION FPFN6R(X)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-F,H,O-Z)
      COMMON /FPFNV1/ PNFP,TMC2FD,TMC2FD0,PMAXC
C
      A=1.D0
      PN=A*(X+PNFP)
      FPFN6R=A*PN*FPRMXW(PN)
C
      RETURN
      END
C
C =============================================================== 
      REAL*8 FUNCTION FPFN7R(X,XM,XP)
C                            
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-F,H,O-Z)
      COMMON /FPFNV1/ PNFP,TMC2FD,TMC2FD0,PMAXC
C
      PMAX=PMAXC
c      write(*,*)PMAX
      A=0.5D0*(PNFP-PMAX)
      PN=A*XP+0.5D0*(PNFP+PMAX)
      B=PN**4/(1.D0+PN**2*TMC2FD0)
      FPFN7R=A*B*FPRMXW(PN)
C
      RETURN
      END
C
C =============================================================== 
C
      REAL*8 FUNCTION FPFN8R(X,XM,XP)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-F,H,O-Z)
      COMMON /FPFNV1/ PNFP,TMC2FD,TMC2FD0,PMAXC
C
      PMAX=PMAXC
      A=0.5D0*(PNFP-PMAX)
      PN=A*XP+0.5D0*(PNFP+PMAX)
      FPFN8R=A*PN**2*FPRMXW(PN)
C
      RETURN
      END
C
C =============================================================== 
C
      REAL*8 FUNCTION FPFN9R(X,XM,XP)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-F,H,O-Z)
      COMMON /FPFNV1/ PNFP,TMC2FD,TMC2FD0,PMAXC
C
      PMAX=PMAXC
      A=0.5D0*(PNFP-PMAX)
      PN=A*XP+0.5D0*(PNFP+PMAX)
      B=PN**2/SQRT(1.D0+PN**2*TMC2FD0)
      FPFN9R=A*B*FPRMXW(PN)
C
      RETURN
      END
C
C =============================================================== 
C
      REAL*8 FUNCTION FPFN10R(X,XM,XP)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-F,H,O-Z)
      COMMON /FPFNV1/ PNFP,TMC2FD,TMC2FD0,PMAXC
C
      PMAX=PMAXC
      A=0.5D0*(PNFP-PMAX)
      PN=A*XP+0.5D0*(PNFP+PMAX)
      B=PN**4/(SQRT(1.D0+PN**2*TMC2FD0))**3
      FPFN10R=A*B*FPRMXW(PN)
C
      RETURN
      END
C
C ===============================================================
C
      REAL*8 FUNCTION FPRMXW(PN)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-F,H,O-Z)
      COMMON /FPFNV1/ PNFP,TMC2FD,TMC2FD0,PMAXC
C
      EX=(1.D0-SQRT(1.D0+PN**2*TMC2FD0))/TMC2FD
      IF (EX.LT.-100.D0)THEN
         FPRMXW=0.D0
      ELSE
         FPRMXW=EXP(EX)
      ENDIF
C
      RETURN
      END
