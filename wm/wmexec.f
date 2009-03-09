C     $Id$
C
C     ****** CALCULATE ANTENNA EXCITATION ******
C
      SUBROUTINE WMEXEC(IERR)
C
      INCLUDE 'wmcomm.inc'
C
      IERR=0
      MODEEG=0
C
      CALL DPCHEK(IERR)
      IF(IERR.NE.0) RETURN
C
      if(mdlwmf.eq.0) then
C
         CALL WMSETG(IERR)
         IF(IERR.NE.0) RETURN
         CALL WMSETJ(IERR)
         IF(IERR.NE.0) RETURN
         CALL WMSOLV
         CALL WMEFLD
         CALL WMBFLD
         CALL WMPABS
         if(myrank.eq.0) then
            CALL WMPFLX
            CALL WMPANT
         endif
      else
         nrmax=nrmax+1
         call wmfem_main
         call wmfem_post
         nrmax=nrmax-1
      endif
C
      IF(MYRANK.EQ.0) THEN
         CALL WMPOUT
         IF(MODELW.EQ.1) CALL WMDOUT(IERR)
      ENDIF
      RETURN
      END
C
C     ****** DEBUG: NO CALCULATION ******
C
      SUBROUTINE WMDEBUG(IERR)
C
      INCLUDE 'wmcomm.inc'
C
      IERR=0
      MODEEG=0
      CALL WMSETG(IERR)
      IF(IERR.NE.0) RETURN
      CALL WMSETJ(IERR)
      IF(IERR.NE.0) RETURN
C
C      CALL WMSOLV
C      CALL WMEFLD
C      CALL WMBFLD
C      CALL WMPABS
C
      IF(MYRANK.EQ.0) THEN
C         CALL WMPFLX
C         CALL WMPANT
C         CALL WMPOUT
C         IF(MODELW.EQ.1) CALL WMDOUT(IERR)
      ENDIF
      RETURN
      END
C
C     ****** CALCULATE ANTENNA CURRENT ******
C     
      SUBROUTINE WMSETJ(IERR)
C
      INCLUDE 'wmcomm.inc'
C
      DIMENSION RHO(0:3,1:3,2:3)
C
      DATA RHO/ 3.832, 1.841, 3.054, 4.201,
     &          7.016, 5.332, 6.706, 8.015,
     &         10.173, 8.536, 9.969,11.346,
     &          2.405, 3.832, 5.136, 6.380,
     &          5.520, 7.016, 8.417, 9.761,
     &          8.654,10.173,11.620,13.015/
C
         IF(RD.LE.RA.OR.RD.GE.RB) THEN
            IF(MYRANK.EQ.0) WRITE(6,*) '!! WMSETJ: RD = (RA+RB)/2'
            RD=0.5D0*(RA+RB)
            IF(MYRANK.EQ.0) 
     &           WRITE(6,'(A,1P3E12.4)') 'RA,RB,RD=',RA,RB,RD
         ENDIF
C
      DO NDX=1,NDSIZ
      DO MDX=1,MDSIZ
         CJANT(1,MDX,NDX)=(0.D0,0.D0)
         CJANT(2,MDX,NDX)=(0.D0,0.D0)
         CJANT(3,MDX,NDX)=(0.D0,0.D0)
      ENDDO
      ENDDO
C
      MODELWG=0
      IF(MODELJ.EQ.0) THEN
         CALL WMCANT
      ELSEIF(MODELJ.EQ.1) THEN
         MODELWG=1
      ELSEIF(MODELJ.EQ.2) THEN
         CJANT(2,1,1)= 1.D0
      ELSEIF(MODELJ.EQ.3) THEN
         CJANT(3,1,1)= 1.D0
      ELSE
         NMODE=MOD(MODELJ,10)
         IMODE=MODELJ/10
         IF(NTH0.LT.0.OR.NTH0.GT.3) GOTO 9000
         IF(NMODE.LT.1.OR.NMODE.GT.3) GOTO 9000
         IF(IMODE.LT.2.OR.IMODE.GT.3) GOTO 9000
         RF =0.5D0*VC/PI
     &            *SQRT((NPH0/RR)**2+(RHO(NTH0,NMODE,IMODE)/RB)**2)
     &            *1.D-6
         RFI=0.D0
         CRF=DCMPLX(RF,RFI)
         CJANT(IMODE,1,1)=(1.D0,0.D0)
      ENDIF
C
      IERR=0
      RETURN
C
 9000 WRITE(6,*) 'XX WMCALJ ERROR'
      IERR=1
      RETURN
      END
C
C     ****** CALCULATE ANTENNA CURRENT ******
C
      SUBROUTINE WMCANT
C
      INCLUDE 'wmcomm.inc'
C
      DIMENSION CJT(MDM,NDM,NAM)
      DIMENSION CJZ(MDM,NDM,NAM)
C
      DO NA=1,NAMAX
         TH1=THJ1(NA)*PI/180.D0
         TH2=THJ2(NA)*PI/180.D0
         PH1=PHJ1(NA)*PI/180.D0
         PH2=PHJ2(NA)*PI/180.D0
C
         CAJ=AJ(NA)*EXP(DCMPLX(0.D0,APH(NA)*PI/180.D0))
C   
      DO ND=NDMIN,NDMAX
         NDX=ND-NDMIN+1
         NN=NPH0+NHC*ND
         IF(NN.EQ.0.OR.ABS(PH2-PH1).LE.1.D-15) THEN
            CJN=-CI
         ELSE
            CJN=(EXP(-CI*NN*PH2)-EXP(-CI*NN*PH1))/(NN*(PH2-PH1))
         ENDIF
      DO MD=MDMIN,MDMAX
         MDX=MD-MDMIN+1
         MM=NTH0+MD
         IF(ABS(MM+BETAJ).LE.0.D0) THEN
            CJMP=-CI*(TH2-TH1)
         ELSE
            CJMP=(EXP(-CI*(MM+BETAJ)*TH2)-EXP(-CI*(MM+BETAJ)*TH1))
     &           /(MM+BETAJ)
         ENDIF
         IF(ABS(MM-BETAJ).LE.0.D0) THEN
            CJMM=-CI*(TH2-TH1)
         ELSE
            CJMM=(EXP(-CI*(MM-BETAJ)*TH2)-EXP(-CI*(MM-BETAJ)*TH1))
     &           /(MM-BETAJ)
         ENDIF
         CJTEMP=CAJ/(8*PI**2)*CJN*(CJMP+CJMM)
         CJT(MDX,NDX,NA)=CJTEMP*COS(2*PI*ANTANG/360.D0)
         CJZ(MDX,NDX,NA)=CJTEMP*SIN(2*PI*ANTANG/360.D0)
      ENDDO
      ENDDO
      ENDDO
C
      DO ND=NDMIN,NDMAX
         NDX=ND-NDMIN+1
      DO MD=MDMIN,MDMAX
         MDX=MD-MDMIN+1
      DO NA=1,NAMAX
         CJANT(2,MDX,NDX)=CJANT(2,MDX,NDX)+CJT(MDX,NDX,NA)
         CJANT(3,MDX,NDX)=CJANT(3,MDX,NDX)+CJZ(MDX,NDX,NA)
         IF(MYRANK.EQ.0) THEN
            IF(NPRINT.GE.3) WRITE(6,'(A,2I4,1P2E15.7)') 
     &                   'NN,MM,CJANT=',
     &                   NPH0+NHC*ND,NTH0+MD,CJANT(2,MDX,NDX)
         ENDIF
      ENDDO
      ENDDO
      ENDDO
C
      RETURN
      END
C
C     ****** CALCULATE WAVEGUIDE FIELD ******
C
      SUBROUTINE WMCWGF
C
      INCLUDE 'wmcomm.inc'
C
      DIMENSION CJT(MDM,NDM,NAM)
      DIMENSION CJZ(MDM,NDM,NAM)
C
      DO NA=1,NAMAX
         TH1=THJ1(NA)*PI/180.D0
         TH2=THJ2(NA)*PI/180.D0
         WGH=RB*ABS(TH2-TH1)
         TH0=0.5D0*(TH1+TH2)
         RKWG=PI/WGH
         DO ND=NDMIN,NDMAX
            NDX=ND-NDMIN+1
            NN=NPH0+NHC*ND
            RKPHI=NN/(RR+RB*COS(TH0))
         ENDDO
C
         CAJ=AJ(NA)*EXP(DCMPLX(0.D0,APH(NA)*PI/180.D0))
C   
      DO ND=NDMIN,NDMAX
         NDX=ND-NDMIN+1
         NN=NPH0+NHC*ND
         IF(NN.EQ.0.OR.ABS(PH2-PH1).LE.1.D-15) THEN
            CJN=-CI
         ELSE
            CJN=(EXP(-CI*NN*PH2)-EXP(-CI*NN*PH1))/(NN*(PH2-PH1))
         ENDIF
      DO MD=MDMIN,MDMAX
         MDX=MD-MDMIN+1
         MM=NTH0+MD
         IF(ABS(MM+BETAJ).LE.0.D0) THEN
            CJMP=-CI*(TH2-TH1)
         ELSE
            CJMP=(EXP(-CI*(MM+BETAJ)*TH2)-EXP(-CI*(MM+BETAJ)*TH1))
     &           /(MM+BETAJ)
         ENDIF
         IF(ABS(MM-BETAJ).LE.0.D0) THEN
            CJMM=-CI*(TH2-TH1)
         ELSE
            CJMM=(EXP(-CI*(MM-BETAJ)*TH2)-EXP(-CI*(MM-BETAJ)*TH1))
     &           /(MM-BETAJ)
         ENDIF
         CJTEMP=CAJ/(8*PI**2)*CJN*(CJMP+CJMM)
         CJT(MDX,NDX,NA)=CJTEMP*COS(2*PI*ANTANG/360.D0)
         CJZ(MDX,NDX,NA)=CJTEMP*SIN(2*PI*ANTANG/360.D0)
      ENDDO
      ENDDO
      ENDDO
C
      DO ND=NDMIN,NDMAX
         NDX=ND-NDMIN+1
      DO MD=MDMIN,MDMAX
         MDX=MD-MDMIN+1
      DO NA=1,NAMAX
         CJANT(2,MDX,NDX)=CJANT(2,MDX,NDX)+CJT(MDX,NDX,NA)
         CJANT(3,MDX,NDX)=CJANT(3,MDX,NDX)+CJZ(MDX,NDX,NA)
         IF(MYRANK.EQ.0) THEN
            IF(NPRINT.GE.3) WRITE(6,'(A,2I4,1P2E15.7)') 
     &                   'NN,MM,CJANT=',
     &                   NPH0+NHC*ND,NTH0+MD,CJANT(2,MDX,NDX)
         ENDIF
      ENDDO
      ENDDO
      ENDDO
C
      RETURN
      END
