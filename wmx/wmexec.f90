! wmexec.f90

MODULE wmexec

  PRIVATE
  PUBLIC wm_exec

CONTAINS

  !     ****** CALCULATE ANTENNA EXCITATION ******
  
  SUBROUTINE wm_exec(IERR)

    USE wmcomm
    USE wmsetg
    USE wmsolv
    USE wmemfp
    USE wmpout
    USE dpparm
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr
    INTEGER:: nnr

    IERR=0
    
    MODEEG=0
    MODELK=1

    CALL wm_setg(IERR)
    IF(IERR.NE.0) RETURN
    
    CALL dpprep(NTHMAX,NRMAX+1,XRHO(1),XRHO(NRMAX+1),IERR)
    IF(IERR.NE.0) RETURN

    CALL wm_setj(IERR)
    IF(IERR.NE.0) RETURN

    CALL wm_setew
     
    CALL wm_solv
    
    CALL wm_efield
    CALL wm_bfield
    CALL wm_pabs
    IF(nrank.EQ.0) THEN
       CALL wm_pwrflux
       CALL wm_pwrant
    ENDIF
    IF(nrank.EQ.0) THEN
       CALL wm_pout
!       IF(nprint.GE.5) CALL wm_dout(ierr)
    ENDIF
    RETURN
  END SUBROUTINE wm_exec
  
!     ****** CALCULATE ANTENNA CURRENT ******

  SUBROUTINE wm_setj(ierr)
    USE wmcomm
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: ierr
    INTEGER:: NDX,MDX,NMODE,IMODE,NA
    REAL(rkind),DIMENSION(0:3,1:3,2:3):: RHO
    DATA RHO/   3.832, 1.841, 3.054, 4.201, &
                7.016, 5.332, 6.706, 8.015, &
               10.173, 8.536, 9.969,11.346, &
                2.405, 3.832, 5.136, 6.380, &
                5.520, 7.016, 8.417, 9.761, &
                8.654,10.173,11.620,13.015/
    
    IF(RD.LE.RA.OR.RD.GE.RB) THEN
       IF(NRANK.EQ.0) WRITE(6,*) '!! WMSETJ: RD = (RA+RB)/2'
       RD=0.5D0*(RA+RB)
       IF(NRANK.EQ.0) &
            WRITE(6,'(A,1P3E12.4)') 'RA,RB,RD=',RA,RB,RD
    ENDIF

    modewg=0
    DO NA=1,NAMAX
       IF(AEWGT(NA)**2+AEWGZ(NA)**2.GT.0.D0) modewg=modewg+1
    END DO
          
    DO NDX=1,NDSIZ
       DO MDX=1,MDSIZ
          CJANT(1,MDX,NDX)=(0.D0,0.D0)
          CJANT(2,MDX,NDX)=(0.D0,0.D0)
          CJANT(3,MDX,NDX)=(0.D0,0.D0)
       ENDDO
    ENDDO

    IF(MODELJ.EQ.0) THEN
       CALL WMCANT
    ELSEIF(MODELJ.EQ.1) THEN
       modewg=0
       DO NA=1,NAMAX
          IF(AEWGT(NA)**2+AEWGZ(NA)**2.GT.0.D0) modewg=modewg+1
       END DO
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
       RF =0.5D0*VC/PI &
            *SQRT((NPH0/RR)**2+(RHO(NTH0,NMODE,IMODE)/RB)**2) &
            *1.D-6
       RFI=0.D0
       CJANT(IMODE,1,1)=(1.D0,0.D0)
    ENDIF
    IERR=0
    RETURN
   
9000 WRITE(6,*) 'XX WMCALJ ERROR'
    IERR=1
    RETURN
  END SUBROUTINE wm_setj

  !     ****** CALCULATE ANTENNA CURRENT ******

  SUBROUTINE WMCANT
    USE wmcomm
    IMPLICIT NONE
    COMPLEX(rkind):: CJT(nthmax,nhhmax,namax),CJZ(nthmax,nhhmax,namax)
    INTEGER:: NA,ND,NDX,NN,MD,MDX,MM
    REAL(rkind):: TH1,TH2,PH1,PH2
    COMPLEX(rkind):: CAJ,CJN,CJMP,CJMM,CJTEMP
    
    DO NA=1,NAMAX
       TH1=THJ1(NA)*PI/180.D0
       TH2=THJ2(NA)*PI/180.D0
       PH1=PHJ1(NA)*PI/180.D0
       PH2=PHJ2(NA)*PI/180.D0
       CAJ=AJ(NA)*EXP(DCMPLX(0.D0,APH(NA)*PI/180.D0))
!   
       DO ND=NDMIN,NDMAX
          NDX=ND-NDMIN+1
          NN=NPH0+NHC*ND
          IF(NN.EQ.0.OR.ABS(PH2-PH1).LE.1.D-15) THEN
             CJN=-CI*EXP(-CI*NN*PH1)
          ELSE
             CJN=(EXP(-CI*NN*PH2)-EXP(-CI*NN*PH1))/(NN*(PH2-PH1))
          ENDIF
!
          DO MD=MDMIN,MDMAX
             MDX=MD-MDMIN+1
             MM=NTH0+MD
             IF(ABS(MM+BETAJ(NA)).LE.0.D0) THEN
                CJMP=-CI*(TH2-TH1)
             ELSE
                CJMP=(EXP(-CI*(MM+BETAJ(NA))*TH2)-EXP(-CI*(MM+BETAJ(NA))*TH1)) &
                     /(MM+BETAJ(NA))
             ENDIF
             IF(ABS(MM-BETAJ(NA)).LE.0.D0) THEN
                CJMM=-CI*(TH2-TH1)
             ELSE
                CJMM=(EXP(-CI*(MM-BETAJ(NA))*TH2)-EXP(-CI*(MM-BETAJ(NA))*TH1)) &
                     /(MM-BETAJ(NA))
             ENDIF
             CJTEMP=CAJ/(8*PI**2)*CJN*(CJMP+CJMM)
             CJT(MDX,NDX,NA)=CJTEMP*COS(2*PI*ANTANG(NA)/360.D0)
             CJZ(MDX,NDX,NA)=CJTEMP*SIN(2*PI*ANTANG(NA)/360.D0)
          ENDDO
       ENDDO
    ENDDO
    DO ND=NDMIN,NDMAX
       NDX=ND-NDMIN+1
       DO MD=MDMIN,MDMAX
          MDX=MD-MDMIN+1
          DO NA=1,NAMAX
             CJANT(2,MDX,NDX)=CJANT(2,MDX,NDX)+CJT(MDX,NDX,NA)
             CJANT(3,MDX,NDX)=CJANT(3,MDX,NDX)+CJZ(MDX,NDX,NA)
             IF(NRANK.EQ.0) THEN
                IF(NPRINT.GE.3) WRITE(6,'(A,2I4,1P2E15.7)')  &
                                        'NN,MM,CJANT=', &
                                        NPH0+NHC*ND,NTH0+MD,CJANT(2,MDX,NDX)
             ENDIF
          ENDDO
       ENDDO
    ENDDO
    RETURN
  END SUBROUTINE WMCANT
  
!     ****** CALCULATE WAVEGUIDE FIELD (simple) ******

  SUBROUTINE wm_setew
    USE wmcomm
    USE wmsub
    IMPLICIT NONE
    COMPLEX(rkind):: cetemp2(nthmax),cetemp3(nthmax)
    INTEGER:: NHH,NTH,NA,ND,NN
    REAL(rkind):: TH1,TH2,PH1,PH2,ANG,RWPH,WTH,TH0,DTH,TH
    COMPLEX(rkind):: CW,CAJ,COEF,CJA,CJB,CETH,CEPH

    DO NHH=1,NHHMAX
       DO NTH=1,NTHMAX
          CEWALL(1,NTH,NHH)=(0.D0,0.D0)
          CEWALL(2,NTH,NHH)=(0.D0,0.D0)
          CEWALL(3,NTH,NHH)=(0.D0,0.D0)
       END DO
    END DO

    CW=2.D0*PI*CMPLX(RF,RFI)*1.D6
    DO NA=1,NAMAX
       TH1=THJ1(NA)*PI/180.D0
       TH2=THJ2(NA)*PI/180.D0
       PH1=PHJ1(NA)*PI/180.D0
       PH2=PHJ2(NA)*PI/180.D0
       ANG=ANTANG(NA)*PI/180.D0
       CAJ=EXP(CMPLX(0.D0,APH(NA)*PI/180.D0))
       RWPH=RR+RB
       WTH=ABS(TH2-TH1)
       TH0=0.5D0*(TH1+TH2)
       DTH=2.D0*PI/NTHMAX
       DO ND=NDMIN,NDMAX
          NHH=ND-NDMIN+1
          NN=NPH0+NHC*ND
          CJA=EXP(-CI*NN*PH1)/(2.D0*PI)
          COEF=-RWPH*CW*SIN(ANG)/VC-NN
          IF(ABS(COEF).EQ.0.D0) THEN
             CJB=PH2-PH1
          ELSE
             CJB=-CI*(EXP(CI*COEF*(PH2-PH1))-1.D0)/COEF
          ENDIF
          DO NTH=1,NTHMAX
             TH=DTH*(NTH-1)
             IF(TH.GE.TH1.AND.TH.LE.TH2)THEN
                CETH=CJA*CJB*AEWGT(NA)
                CEPH=CJA*CJB*AEWGZ(NA)*COS((TH-TH0)/WTH*PI)
                DO NHH=1,NHHMAX
                   CEWALL(2,NTH,NHH)=CEWALL(2,NTH,NHH)+CETH
                   CEWALL(3,NTH,NHH)=CEWALL(3,NTH,NHH)+CEPH
                END DO
             ELSE IF(TH >= TH1+2.D0*PI .AND. TH <= TH2+2.D0*PI) THEN
                CETH=CJA*CJB*CAJ*AEWGT(NA)
                CEPH=CJA*CJB*CAJ*AEWGZ(NA)*COS((TH-TH0-2.D0*PI)/WTH*PI)
                DO NHH=1,NHHMAX
                   CEWALL(2,NTH,NHH)=CEWALL(2,NTH,NHH)+CETH
                   CEWALL(3,NTH,NHH)=CEWALL(3,NTH,NHH)+CEPH
                END DO
             END IF
          END DO
       END DO
    END DO

    DO nhh=1,nhhmax
       DO nth=1,nthmax
          cetemp2(nth)=cewall(2,nth,nhh)
          cetemp3(nth)=cewall(3,nth,nhh)
       END DO
       CALL WMSUBC(cetemp2)
       CALL WMSUBC(cetemp3)
       DO nth=1,nthmax
          cewall(2,nth,nhh)=cetemp2(nth)
          cewall(3,nth,nhh)=cetemp3(nth)
       END DO
    END DO
    RETURN
  END SUBROUTINE wm_setew
END MODULE wmexec
