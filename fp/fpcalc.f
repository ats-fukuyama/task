C     $Id$
C
C ******************************************************************
C
C     CALCULATION OF COLLISINAL DIFUSION COEFFICIENTS: DC AND FC
C
C ******************************************************************
C
      SUBROUTINE FPCALC
C
      INCLUDE 'fpcomm.inc'
C
      DO NR=1,NRMAX
         DO NP=1,NPMAX+1
         DO NTH=1,NTHMAX
            DCPP(NTH,NP,NR)=0.D0
            DCPT(NTH,NP,NR)=0.D0
            FCPP(NTH,NP,NR)=0.D0
            DO NS=1,NSMAX
               DCPP2(NTH,NP,NR,NS)=0.D0
               DCPT2(NTH,NP,NR,NS)=0.D0
               FCPP2(NTH,NP,NR,NS)=0.D0
            END DO
         ENDDO
         ENDDO
         DO NP=1,NPMAX
         DO NTH=1,NTHMAX+1
            DCTP(NTH,NP,NR)=0.D0
            DCTT(NTH,NP,NR)=0.D0
            FCTH(NTH,NP,NR)=0.D0
            Do NS=1,NSMAX
               DCTP2(NTH,NP,NR,NS)=0.D0
               DCTT2(NTH,NP,NR,NS)=0.D0
               FCTH2(NTH,NP,NR,NS)=0.D0
            END DO
         ENDDO
         ENDDO
      ENDDO
C
      DO NR=1,NRMAX
         DO NS=1,NSMAX
C
c            if(nr.eq.1) write(6,'(A,I8,1P2E12.4)') 
c     &           ' NS,RN,RNFD=',NS,RN(NS),RNFD(NR,NS)
C
            IF(MODELC.EQ.0.or.MODELC.eq.1) THEN
               CALL FPCALC_L(NR,NS)
            ELSEIF(MODELC.EQ.2.or.MODELC.eq.3) THEN
               IF(NS.EQ.NSFP) THEN
                  CALL FPCALC_NL(NR,NS)
               ELSE
                  CALL FPCALC_L(NR,NS)
               ENDIF
            ELSEIF(MODELC.EQ.4) THEN
               CALL FPCALC_NL(NR,NS)
            ELSEIF(MODELC.EQ.5) THEN
               IF(NS.NE.NSFP) THEN
                  CALL FPCALC_L(NR,NS)
               ENDIF
            ELSEIF(MODELC.EQ.6) THEN
               IF(NS.NE.NSFP) THEN
                  CALL FPCALC_NL(NR,NS)
               ENDIF
            ELSEIF(MODELC.EQ.-1) THEN
               IF(NS.EQ.NSFP) THEN
                  CALL FPCALC_L(NR,NS)
               ENDIF
            ELSEIF(MODELC.EQ.-2) THEN
               IF(NS.EQ.NSFP) THEN
                  CALL FPCALC_NL(NR,NS)
               ENDIF
            ENDIF
c            write(6,'(2I8,1P3E12.4)')
c     &           NR,NS,DCPP(5,5,NR),DCTT(5,5,NR),FCPP(5,5,NR)
c     divide coefficients by species
            DO NP=1,NPMAX
               DO NTH=1,NTHMAX
                  if(NS.eq.1)then
                     DCPP2(NTH,NP,NR,NS)=DCPP(NTH,NP,NR)
                     DCPT2(NTH,NP,NR,NS)=DCPT(NTH,NP,NR)
                     DCTP2(NTH,NP,NR,NS)=DCTP(NTH,NP,NR)
                     DCTT2(NTH,NP,NR,NS)=DCTT(NTH,NP,NR)
                     FCPP2(NTH,NP,NR,NS)=FCPP(NTH,NP,NR)
                     FCTH2(NTH,NP,NR,NS)=FCTH(NTH,NP,NR)
                  else
                     differ1=0.D0
                     differ2=0.D0
                     differ3=0.D0
                     differ4=0.D0
                     differ5=0.D0
                     differ6=0.D0
                     Do NS2=1,NS-1
                        differ1=differ1+DCPP2(NTH,NP,NR,NS2)
                        differ2=differ2+DCPT2(NTH,NP,NR,NS2)
                        differ3=differ3+DCTP2(NTH,NP,NR,NS2)
                        differ4=differ4+DCTT2(NTH,NP,NR,NS2)
                        differ5=differ5+FCPP2(NTH,NP,NR,NS2)
                        differ6=differ6+FCTH2(NTH,NP,NR,NS2)
                     END DO
                     DCPP2(NTH,NP,NR,NS)=DCPP(NTH,NP,NR)-differ1
                     DCPT2(NTH,NP,NR,NS)=DCPT(NTH,NP,NR)-differ2
                     DCTP2(NTH,NP,NR,NS)=DCTP(NTH,NP,NR)-differ3
                     DCTT2(NTH,NP,NR,NS)=DCTT(NTH,NP,NR)-differ4
                     FCPP2(NTH,NP,NR,NS)=FCPP(NTH,NP,NR)-differ5
                     FCTH2(NTH,NP,NR,NS)=FCTH(NTH,NP,NR)-differ6
                  end if
               END DO
c     check for collision coef.
               if(NS.eq.NSMAX+1.and.NP.ne.1)then
                  IF(MODELR.eq.0)THEN
                     TMC2FP0=0.D0
                  ELSE IF(MODELR.eq.1)THEN
                     TMC2FP0=(PTFP0/(AMFP*VC))**2
                  END IF
                  PNFP=PG(NP)
                  RGAMA=SQRT(1.D0+PNFP**2*TMC2FP0)
                  pbar=PNFP/RGAMA
c                  write(*,1274)FCPP2(1,NP,NR,1)/DCPP2(1,NP,NR,1)/pbar,
c     &                 FCPP2(1,NP,NR,2)/DCPP2(1,NP,NR,2)/pbar,
c     &                 FCPP(1,NP,NR)/DCPP(1,NP,NR)/pbar
                  write(*,1281)NP,FCPP2(1,NP,1,1)/DCPP2(1,NP,1,1)/pbar
     &                 ,FCPP2(1,NP,1,2)/DCPP2(1,NP,1,2)/pbar
     &                 ,FCPP2(1,NP,1,1),FCPP2(1,NP,1,2),FCPP(1,NP,1)
               end if
c     end of check
            END DO
         ENDDO
 1274    FORMAT(3e14.5)
 1281    FORMAT(I2,5e14.5)
C
C        ----- Simple ion collision term using ZEFF -----
C
         IF(MODELC.LT.0) THEN
            NS=1
            RNNL=RNFD(NR,NS)/RNFP0
            RNUDL=RNUD(NR,NS)*RNNL
            DO NP=1,NPMAX
               PFPL=PM(NP)*PTFP0
               VFPL=PFPL/AMFP
               U=VFPL/VTFP0
               DCTTL=0.5D0*ZEFF*RNUDL/U
               DO NTH=1,NTHMAX+1
                  DCTT(NTH,NP,NR)=DCTT(NTH,NP,NR)+DCTTL
               ENDDO
            ENDDO
         ENDIF

C
         IF(MODELA.EQ.1) THEN
            IF(MODELC.EQ.0.or.MODELC.eq.1) THEN
               CALL FPCALC_LAV(NR)
            ELSEIF(MODELC.EQ.2.or.MODELC.eq.3) THEN
               CALL FPCALC_NLAV(NR)
            ELSEIF(MODELC.EQ.4) THEN
               CALL FPCALC_NLAV(NR)
            ELSEIF(MODELC.EQ.5) THEN
               CALL FPCALC_LAV(NR)
            ELSEIF(MODELC.EQ.6) THEN
               CALL FPCALC_NLAV(NR)
            ELSEIF(MODELC.EQ.-1) THEN
               CALL FPCALC_LAV(NR)
            ELSEIF(MODELC.EQ.-2) THEN
               CALL FPCALC_NLAV(NR)
            ENDIF
         ENDIF
      ENDDO

C
      IF(MOD(IDBGFP/8,2).EQ.1) THEN
C
C     +++ plot of D_coll +++
C
         CALL PAGES
         CALL FPGRFA(1,DCPP,PG,1,'@DCPP(P)@',NTHM,NPM,NRM,
     &                                       NTHMAX,NPMAX,NRMAX)
         CALL FPGRFA(2,DCTT,PM,2,'@DCTT(P)@',NTHM,NPM,NRM,
     &                                       NTHMAX,NPMAX,NRMAX)
         CALL FPGRFA(3,FCPP,PG,1,'@FCPP(P)@',NTHM,NPM,NRM,
     &                                       NTHMAX,NPMAX,NRMAX)
         CALL PAGEE
      ENDIF
C
      IF(MOD(IDBGFP/16,2).EQ.1) THEN
C
C     +++ plot of D_coll +++
C
         CALL PAGES
         CALL FPGRFB(1,DCPP,THM,1,'@DCPP(TH)@',NTHM,NPM,NRM,
     &                                         NTHMAX,NPMAX,NRMAX)
         CALL FPGRFB(2,DCTT,THG,2,'@DCTT(TH)@',NTHM,NPM,NRM,
     &                                         NTHMAX,NPMAX,NRMAX)
         CALL FPGRFB(3,FCPP,THM,1,'@FCPP(TH)@',NTHM,NPM,NRM,
     &                                         NTHMAX,NPMAX,NRMAX)
         CALL PAGEE
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE FPGRFA(ID,DATA,P,IND,TITLE,NTHM,NPM,NRM,
     &                                      NTHMAX,NPMAX,NRMAX)
C
      implicit none
      integer, intent(IN)::  ID,IND,NTHM,NPM,NRM,NTHMAX,NPMAX,NRMAX
      real(8), dimension(NTHM,NPM,NRM), intent(IN):: DATA
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
            work(np,ng)=data(nth,np,nr)
         enddo
      enddo
            
      CALL GRD1D(ID,p,work,npm,npmaxg,9,TITLE,0)

 9999 return
      end
C
      SUBROUTINE FPGRFB(ID,DATA,TH,IND,TITLE,NTHM,NPM,NRM,
     &                                      NTHMAX,NPMAX,NRMAX)
C
      implicit none
      integer, intent(IN)::  ID,IND,NTHM,NPM,NRM,NTHMAX,NPMAX,NRMAX
      real(8), dimension(NTHM,NPM,NRM), intent(IN):: DATA
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
            work(nth,ng)=data(nth,np,nr)
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
      SUBROUTINE FPCALC_L(NR,NS)
C
      INCLUDE 'fpcomm.inc'
      COMMON /FPFNV1/ PNFP,TMC2FD,TMC2FD0
      EXTERNAL FPFN0R,FPFN1R,FPFN2R,FPFN3R,FPFN4R,FPFN5R,FPFN6R
     &     ,FPFN7R,FPFN8R,FPFN9R,FPFN10R

c            if(nr.eq.1) write(6,'(A,I8,1P2E12.4)') 
c     &           ':NS,RN,RNFD=',NS,RN(NS),RNFD(NR,NS)
C     ------ define --------
      AMFD=PA(NS)*AMP
      AEFD=PZ(NS)*AEE 
      PTFPL=PTFP(NR)
      VTFPL=VTFP(NR)
      RNNL=RNFD(NR,NS)/RNFP0
      RNUFL=RNUF(NR,NS)*RNNL
      RNUDL=RNUD(NR,NS)*RNNL
      PTFDL=PTFD(NR,NS)
      VTFDL=VTFD(NR,NS)

c      AEFP=PZ(NSFP)*AEE
      RGAMH=RNUD(NR,NS)*SQRT(2.D0)*VTFD(NR,NS)*AMFP
     &        /(RNFP0*PTFP0*1.D20)
      RNFD0=PN(NS)
      RNFDS=PNS(NS)
      RTFD0=(PTPR(NS)+2.D0*PTPP(NS))/3.D0
      RTFDS=PTS(NS)
C
C      IF(MODELR.EQ.0) THEN
         PTFD0=SQRT(RTFD0*1.D3*AEE*AMFD)
         VTFD0=SQRT(RTFD0*1.D3*AEE/AMFD)
C      ELSE
C         RKE=RTFD0*1.D3*AEE
C         PTFD0=SQRT(RKE*RKE+2.D0*RKE*AMFD*VC*VC)/VC
C         VTFD0=PTFD0/SQRT(AMFD**2+PTFD0**2/VC**2)
C      ENDIF

C
C     ----- Non-Relativistic -----
C
      IF(MODELR.EQ.0) THEN 
         DO NP=1,NPMAX+1
            IF(NP.EQ.1) THEN
               DCPPL=RGAMH*RNFD(NR,NS)*1.D20*(2.D0/(3.D0*SQRT(PI))) 
     &                    *(VTFP0/(SQRT(2.D0)*VTFD(NR,NS))) 
               FCPPL=0.D0
            ELSE
               PFPL=PG(NP)*PTFP0
               VFPL=PFPL/AMFP
               V=VFPL/VTFP0
               U=VFPL/(SQRT(2.D0)*VTFD(NR,NS))
               DCPPL= 0.5D0*RNUDL/U   *(ERF0(U)/U**2-ERF1(U)/U)
               FCPPL=-      RNUFL/U**2*(ERF0(U)-U*ERF1(U))
c               if(NSFP.eq.1.and.NS.eq.2) DCPPL=0.D0
c               if(NSFP.eq.1.and.NS.eq.2) FCPPL=0.D0
C               WRITE(6,'(I5,1P4E12.4)')NP,DCPPL,FCPPL,FCPPL/DCPPL/PG(NP)
            ENDIF
            DO NTH=1,NTHMAX
               DCPP(NTH,NP,NR)=DCPP(NTH,NP,NR)+DCPPL
               FCPP(NTH,NP,NR)=FCPP(NTH,NP,NR)+FCPPL
            ENDDO
         ENDDO

         DO NP=1,NPMAX
            PFPL=PM(NP)*PTFP0
            VFPL=PFPL/AMFP
            V=VFPL/VTFP0
            U=VFPL/(SQRT(2.D0)*VTFD(NR,NS))
            DCTTL= 0.25D0*RNUDL/U
     &                   *((2.D0-1.D0/U**2)*ERF0(U)+ERF1(U)/U)
C
c     for check conductivity
c            if(NSFP.eq.1.and.NS.eq.2)then
c               RTE=(RTPR(1)+RTPP(1)*2.D0)/3.D0
c              rZI = -PZ(2)/PZ(1)/(14.9D0-0.5D0*LOG(RN(1))+LOG(RTE))*
c     &              (15.2D0-0.5D0*LOG(RN(1))+LOG(RTE))
c              DCTTL = RNUD(NR,1)*SQRT(2.D0)*VTFD(NR,1)*AMFP
c     &        /(RNFP0*PTFP0*1.D20)*RNFD(1,NS)*1.D20*PTFP0/AMFP
c     &             *rZI/VFPL*0.5D0
c            endif
            DO NTH=1,NTHMAX+1
               DCTT(NTH,NP,NR)=DCTT(NTH,NP,NR)+DCTTL
            ENDDO
         ENDDO

C     ----- Relativistic -----
C
      ELSE
         DO NP=1,NPMAX+1
            IF(NP.EQ.1) THEN
C              DCPPL=RNUDL*(2.D0/(3.D0*SQRT(PI)))
C     &                    *(VTFP0/(SQRT(2.D0)*VTFD(NR,NS)))  
              DCPPL=RGAMH*RNFD(NR,NS)*1.D20*(2.D0/(3.D0*SQRT(PI)))
     &                    *(VTFP0/(SQRT(2.D0)*VTFD(NR,NS)))
              FCPPL=0.D0
C               WRITE(6,'(I5,1P4E12.4)') NP,DCPPL,FCPPL
            ELSE
               PFPL=PG(NP)*PTFP0
               VFPL=PFPL/SQRT(AMFP**2+PFPL**2/VC**2)
               PNFPL=PG(NP)
               PNFP=PNFPL
               TMC2FD =(PTFDL/(AMFD*VC))**2
               TMC2FD0=(PTFD0/(AMFD*VC))**2
               TMC2FP0=(PTFP0/(AMFP*VC))**2
               RGAMA=SQRT(1.D0+PNFP**2*TMC2FP0)
               vtatb=(AMFD*PTFP0)/(AMFP*PTFD0)
               ptatb=PG(NP)/RGAMA
               PCRIT=SQRT(vtatb**2/(1.D0-TMC2FD0*vtatb**2*ptatb**2))
     &              *ptatb
              IF(PCRIT.le.PMAX)THEN
                  CALL DEHIFT(RINT0,ES0,H0DE,EPSDE,0,FPFN0R)
C                  PNFP=PNFPL*PTFP0*AMFD/(PTFD0*AMFP)
                  PNFP=PCRIT
                  CALL DEFT  (RINT1,ES1,H0DE,EPSDE,0,FPFN1R)
                  CALL DEHIFT(RINT2,ES2,H0DE,EPSDE,0,FPFN2R)
                  PNFP=PNFPL
                  DCPPL=RGAMH/(3.D0*RINT0)*( (AMFP**2*PTFD0**2*RGAMA**3)
     &              /(AMFD**2*PTFP0**2*PNFP**3)*RINT1+
     &              (AMFD*PTFP0)/(AMFP*PTFD0)*RINT2 )
     &             *RNFD(NR,NS)*1.D20
C                  PNFP=PNFPL*PTFP0*AMFD/(PTFD0*AMFP)
                  PNFP=PCRIT
                  CALL DEFT  (RINT4,ES4,H0DE,EPSDE,0,FPFN4R)
                  CALL DEFT  (RINT5,ES5,H0DE,EPSDE,0,FPFN5R)
                  CALL DEHIFT(RINT6,ES6,H0DE,EPSDE,0,FPFN6R)
                  PNFP=PNFPL
                  FCPPL=-RGAMH/(3.D0*RINT0)*(
     &              (AMFP*RGAMA**2)/(AMFD*PNFP**2)*(3.D0*RINT4
     &              -TMC2FD0*RINT5)+2.D0*(PTFP0*PNFP)/(PTFD0*RGAMA)
     &              *TMC2FP0*RINT6 )
     &             *RNFD(NR,NS)*1.D20
               ELSEIF(PCRIT.gt.PMAX)THEN
                  CALL DEHIFT(RINT0,ES0,H0DE,EPSDE,0,FPFN0R)
                  PNFP=PMAX
                  CALL DEFT  (RINT1,ES1,H0DE,EPSDE,0,FPFN1R)
C                  PNFP=PNFPL*PTFP0*AMFD/(PTFD0*AMFP)
                  PNFP=PCRIT
                  CALL DEFT  (RINT7,ES7,H0DE,EPSDE,0,FPFN7R)
                  CALL DEHIFT(RINT2,ES2,H0DE,EPSDE,0,FPFN2R)
                  PNFP=PNFPL
                  DCPPL=RGAMH/(3.D0*RINT0)*( (AMFP**2*PTFD0**2*RGAMA**3)
     &              /(AMFD**2*PTFP0**2*PNFP**3)*(RINT1+RINT7)+
     &              (AMFD*PTFP0)/(AMFP*PTFD0)*RINT2 )
     &             *RNFD(NR,NS)*1.D20
                  PNFP=PMAX
                  CALL DEFT  (RINT4,ES4,H0DE,EPSDE,0,FPFN4R)
                  CALL DEFT  (RINT5,ES5,H0DE,EPSDE,0,FPFN5R)
C                  PNFP=PNFPL*PTFP0*AMFD/(PTFD0*AMFP)
                  PNFP=PCRIT
                  CALL DEFT  (RINT8,ES8,H0DE,EPSDE,0,FPFN9R)
                  CALL DEFT  (RINT9,ES9,H0DE,EPSDE,0,FPFN10R)
                  CALL DEHIFT(RINT6,ES6,H0DE,EPSDE,0,FPFN6R)
                  PNFP=PNFPL
                  FCPPL=-RGAMH/(3.D0*RINT0)*(
     &              (AMFP*RGAMA**2)/(AMFD*PNFP**2)*(3.D0*(RINT4+RINT8)
     &           -TMC2FD0*(RINT5+RINT9))+2.D0*(PTFP0*PNFP)/(PTFD0*RGAMA)
     &              *TMC2FP0*RINT6 )
     &             *RNFD(NR,NS)*1.D20
               ENDIF
C               WRITE(6,'(I5,1P4E12.4)') NP,DCPPL,FCPPL,FCPPL/DCPPL/PNFP
C     &              *RGAMA,PNFP/RGAMA
            ENDIF
            DO NTH=1,NTHMAX
               DCPP(NTH,NP,NR)=DCPP(NTH,NP,NR)+DCPPL
               FCPP(NTH,NP,NR)=FCPP(NTH,NP,NR)+FCPPL
            ENDDO
         ENDDO

         DO NP=1,NPMAX
            PFPL=PM(NP)*PTFP0
            VFPL=PFPL/SQRT(AMFP**2+PTFPL**2/VC**2)
            PNFPL=PM(NP)
            PNFP=PNFPL
            TMC2FD =PTFDL**2/(AMFD*VC)**2
            TMC2FD0=PTFD0**2/(AMFD*VC)**2
            TMC2FP0=PTFP0**2/(AMFP*VC)**2
            RGAMA=SQRT(1.D0+PNFP**2*TMC2FP0)
            vtatb=(AMFD*PTFP0)/(AMFP*PTFD0)
            ptatb=PM(NP)/RGAMA
            PCRIT=SQRT(vtatb**2/(1.D0-TMC2FD0*vtatb**2*ptatb**2))
     &           *ptatb
            IF(PCRIT.le.PMAX)THEN
               CALL DEHIFT(RINT0,ES0,H0DE,EPSDE,0,FPFN0R)
               PNFP=PNFPL*PTFP0*AMFD/(PTFD0*AMFP)
               CALL DEFT  (RINT1,ES1,H0DE,EPSDE,0,FPFN1R)
               CALL DEHIFT(RINT2,ES2,H0DE,EPSDE,0,FPFN2R)
               CALL DEFT  (RINT3,ES3,H0DE,EPSDE,0,FPFN3R)
               PNFP=PNFPL
               DCTTL=RGAMH/(3.D0*RINT0)*(
     &            1.5D0*RGAMA/PNFP*RINT3
     &           -0.5D0*(AMFP**2*PTFD0**2*RGAMA**3)
     &           /(AMFD**2*PTFP0**2*PNFP**3)*RINT1
     &           +(AMFD*PTFP0)/(AMFP*PTFD0)*RINT2 )
     &          *RNFD(NR,NS)*1.D20
            ELSEIF(PCRIT.gt.PMAX)THEN
               CALL DEHIFT(RINT0,ES0,H0DE,EPSDE,0,FPFN0R)
               PNFP=PMAX
               CALL DEFT  (RINT1,ES1,H0DE,EPSDE,0,FPFN1R)
               CALL DEFT  (RINT3,ES3,H0DE,EPSDE,0,FPFN3R)
               PNFP=PNFPL*PTFP0*AMFD/(PTFD0*AMFP)
               CALL DEFT  (RINT4,ES1,H0DE,EPSDE,0,FPFN7R)
               CALL DEFT  (RINT5,ES3,H0DE,EPSDE,0,FPFN8R)
               CALL DEHIFT(RINT2,ES2,H0DE,EPSDE,0,FPFN2R)
               PNFP=PNFPL
               DCTTL=RGAMH/(3.D0*RINT0)*(
     &            1.5D0*RGAMA/PNFP*(RINT3+RINT5)
     &           -0.5D0*(AMFP**2*PTFD0**2*RGAMA**3)
     &           /(AMFD**2*PTFP0**2*PNFP**3)*(RINT1+RINT4)
     &           +(AMFD*PTFP0)/(AMFP*PTFD0)*RINT2 )
     &          *RNFD(NR,NS)*1.D20
            ENDIF

            DO NTH=1,NTHMAX+1
               DCTT(NTH,NP,NR)=DCTT(NTH,NP,NR)+DCTTL
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
      SUBROUTINE FPCALC_LAV(NR)
C
      INCLUDE 'fpcomm.inc'
C
C         DO NR=1,NRMAX
            DO NTH=1,NTHMAX
               FACT=RLAMDA(NTH,NR)
               DO NP=1,NPMAX+1
                  DCPP(NTH,NP,NR)=FACT*DCPP(NTH,NP,NR)
                  FCPP(NTH,NP,NR)=FACT*FCPP(NTH,NP,NR)
c-----------------divide
                  DO NS=1,NSMAX
                     DCPP2(NTH,NP,NR,NS)=FACT*DCPP2(NTH,NP,NR,NS)
                     FCPP2(NTH,NP,NR,NS)=FACT*FCPP2(NTH,NP,NR,NS)
                  END DO
c-----------------end of divide
               ENDDO
            ENDDO
C
            DO NTH=1,NTHMAX+1
               FACT=RLAMDC(NTH,NR)
               DO NP=1,NPMAX
                  DCTT(NTH,NP,NR)=FACT*DCTT(NTH,NP,NR)
c-----------------divide
                  DO NS=1,NSMAX
                     DCTT2(NTH,NP,NR,NS)=FACT*DCTT2(NTH,NP,NR,NS)                     
                  END DO
c-----------------end of divide
               ENDDO
            ENDDO
C         END DO
C         DO NR=1,NRMAX
            DO NP=1,NPMAX+1
               DCPP(ITL(NR),NP,NR)
     &              =RLAMDA(ITL(NR),NR)/4.D0
     &                    *( DCPP(ITL(NR)-1,NP,NR)/RLAMDA(ITL(NR)-1,NR)
     &                      +DCPP(ITL(NR)+1,NP,NR)/RLAMDA(ITL(NR)+1,NR)
     &                      +DCPP(ITU(NR)-1,NP,NR)/RLAMDA(ITU(NR)-1,NR)
     &                      +DCPP(ITU(NR)+1,NP,NR)/RLAMDA(ITU(NR)+1,NR))

C
               FCPP(ITL(NR),NP,NR)
     &              =RLAMDA(ITL(NR),NR)/4.D0
     &                    *( FCPP(ITL(NR)-1,NP,NR)/RLAMDA(ITL(NR)-1,NR)
     &                      +FCPP(ITL(NR)+1,NP,NR)/RLAMDA(ITL(NR)+1,NR)
     &                      +FCPP(ITU(NR)-1,NP,NR)/RLAMDA(ITU(NR)-1,NR)
     &                      +FCPP(ITU(NR)+1,NP,NR)/RLAMDA(ITU(NR)+1,NR))
               DCPP(ITU(NR),NP,NR)=DCPP(ITL(NR),NP,NR)
               FCPP(ITU(NR),NP,NR)=FCPP(ITL(NR),NP,NR)
c--------------divide coeff. by species
               Do NS=1,NSMAX
               DCPP2(ITL(NR),NP,NR,NS)
     &              =RLAMDA(ITL(NR),NR)/4.D0
     &                *( DCPP2(ITL(NR)-1,NP,NR,NS)/RLAMDA(ITL(NR)-1,NR)
     &                  +DCPP2(ITL(NR)+1,NP,NR,NS)/RLAMDA(ITL(NR)+1,NR)
     &                  +DCPP2(ITU(NR)-1,NP,NR,NS)/RLAMDA(ITU(NR)-1,NR)
     &                  +DCPP2(ITU(NR)+1,NP,NR,NS)/RLAMDA(ITU(NR)+1,NR))

               FCPP2(ITL(NR),NP,NR,NS)
     &              =RLAMDA(ITL(NR),NR)/4.D0
     &                *( FCPP2(ITL(NR)-1,NP,NR,NS)/RLAMDA(ITL(NR)-1,NR)
     &                  +FCPP2(ITL(NR)+1,NP,NR,NS)/RLAMDA(ITL(NR)+1,NR)
     &                  +FCPP2(ITU(NR)-1,NP,NR,NS)/RLAMDA(ITU(NR)-1,NR)
     &                  +FCPP2(ITU(NR)+1,NP,NR,NS)/RLAMDA(ITU(NR)+1,NR))
               DCPP2(ITU(NR),NP,NR,NS)=DCPP2(ITL(NR),NP,NR,NS)
               FCPP2(ITU(NR),NP,NR,NS)=FCPP2(ITL(NR),NP,NR,NS)
               END DO
c--------------end of divide
            ENDDO
C         ENDDO

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
      COMMON /FPFNV1/ PNFP,TMC2FD,TMC2FD0
C
C      A=PNFD
      PN=X
      FPFN0R=PN**2*FPRMXW(PN)
C
      RETURN
      END
C---------------------------------------------------------------
      REAL*8 FUNCTION FPFN1R(X,XM,XP)
C
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-F,H,O-Z)
      COMMON /FPFNV1/ PNFP,TMC2FD,TMC2FD0
C
      XX=XM
      XX=X
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
      COMMON /FPFNV1/ PNFP,TMC2FD,TMC2FD0
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
      COMMON /FPFNV1/ PNFP,TMC2FD,TMC2FD0
C
      XX=X
      XX=XM
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
      COMMON /FPFNV1/ PNFP,TMC2FD,TMC2FD0
C
      XX=X
      XX=XM
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
      COMMON /FPFNV1/ PNFP,TMC2FD,TMC2FD0
C
      XX=X
      XX=XM
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
      COMMON /FPFNV1/ PNFP,TMC2FD,TMC2FD0
C
      A=1.D0
      PN=A*(X+PNFP)
      FPFN6R=A*PN*FPRMXW(PN)
C      FPFN6R=A*PN**2*FPRMXW(PN)
C
      RETURN
      END
C
C =============================================================== 
      REAL*8 FUNCTION FPFN7R(X,XM,XP)
C                            
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-F,H,O-Z)
      COMMON /FPFNV1/ PNFP,TMC2FD,TMC2FD0
C
      XX=XM
      XX=X
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
      COMMON /FPFNV1/ PNFP,TMC2FD,TMC2FD0
C
      XX=X
      XX=XM
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
      COMMON /FPFNV1/ PNFP,TMC2FD,TMC2FD0
C
      XX=X
      XX=XM
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
      COMMON /FPFNV1/ PNFP,TMC2FD,TMC2FD0
C
      XX=X
      XX=XM
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
      COMMON /FPFNV1/ PNFP,TMC2FD,TMC2FD0
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
