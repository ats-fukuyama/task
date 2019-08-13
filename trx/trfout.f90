!     ***********************************************************

!          FILE Output PROFILE DATA

!     ***********************************************************

      SUBROUTINE trfout

      USE TRCOMM,ONLY : GT,GRM,GVT,GVRT,GVR,NGT,NGR,NRMAX,NCTM,NCRTM,NCGM, &
          RM,RN,RT,ANC,ANFE,AD,AV,AK,PIN,POH,PNB,PNF,PEX,PRF,PFCL, &
          SSIN,AJBS,AJ,QP,ZEFF,PRB,PRC,PRL,AR1RHOG,AR2RHOG,PVOLRHOG, &
          BB,BP,EZOH,VTOR,VPOL,AJOH,ER
      USE libfio
      IMPLICIT NONE
      CHARACTER(LEN=80):: LINE
      CHARACTER(LEN=80),SAVE:: FILENAME='stdout'
      INTEGER,SAVE:: NFL=6
      CHARACTER(LEN=2):: KID
      CHARACTER(LEN=8):: KNCTM,KNCGM,KNCRTM
      INTEGER:: I,ID,IERR,N,NR
      REAL(8):: VC(32)

      WRITE(6,'(A,I3,A,A)') 'NFL=',NFL, &
                            '  FILENAME=',FILENAME(1:LEN_TRIM(FILENAME))
1     CALL KKINT(NCTM,KNCTM)
      CALL KKINT(NCGM,KNCGM)
      CALL KKINT(NCRTM,KNCRTM)
      WRITE(6,'(13A)') &
           '# FOUT: GT[0-',KNCTM(1:LEN_TRIM(KNCTM)),'], ', &
                   'RT[0-',KNCRTM(1:LEN_TRIM(KNCRTM)),'], ', &
                   'RN[0-',KNCRTM(1:LEN_TRIM(KNCRTM)),'], ', &
                   'RG[0-',KNCGM(1:LEN_TRIM(KNCGM)),'], ', &
                   'A:all,CP,CN:csv,F:filename,?:help,X:exit'
      READ(5,'(A80)',END=9000,ERR=1) LINE
      KID=LINE(1:2)
      CALL GUCPTL(KID(1:1))
      CALL GUCPTL(KID(2:2))
      IF(KID(1:1).EQ.'X') GOTO 9000
      IF(KID(1:1).EQ.'F') THEN
3        WRITE(6,'(A)') '# output file name?'
         READ(5,*,END=1,ERR=3) FILENAME
         IF(FILENAME(1:6).EQ.'stdout') THEN
            NFL=6
         ELSE
            NFL=21
            CALL FWOPEN(NFL,FILENAME,1,0,'FOUT',IERR)
            IF(IERR.NE.0) GOTO 1
         ENDIF
         GOTO 1
      ENDIF
      IF(KID(1:1).EQ.'?') THEN
         CALL VIEWGTLIST(NCTM)
         CALL VIEWRTLIST(NCRTM)
         GOTO 1
      ENDIF
      IF(KID(1:1).EQ.'A') THEN
         DO I=1,NCTM
            CALL TRF1DGT(NFL,GT,GVT,NGT,I)
         ENDDO
         DO I=1,NCRTM
            CALL TRF2DRT(NFL,GRM,GT,GVRT,NRMAX,NGT,I)
         ENDDO
         DO I=1,NCRTM
            CALL TRF2DRN(NFL,GRM,GVRT,NRMAX,NGT,I)
         ENDDO
         DO I=1,NCGM
            CALL TRF2DRG(NFL,GRM,GT,GVR,NRMAX,NGR,I)
         ENDDO
         GOTO 1
      ENDIF

!     ----- csv file for ITPA Particle transport benchmark -----

      IF(KID.EQ.'CP') THEN
         WRITE(NFL,'(4A)') &
              ',10^19/m^3,10^19/m^3,10^19/m^3,10^19/m^3,10^19/m^3,m^2/s,', &
              'm/s,m^2/s,keV,keV,MW/m^3,MW/m^3,MW/m^3,MW/m^3,MW/m^3,MW/m^3,', &
              'MW/m^3, 10^19 /m^3/s, 10^19 /m^3/s,MA/m^2,MA/m^2,,,MW/m^3,',&
              'MW/m^3,MW/m^3,10^19/m^3,m^3, , ,m/s'
         WRITE(NFL,'(3A)') &
              'rho,ne,ni,nBe,nAr,nHe,D,V,ChiI=ChiE,Te,Ti,Pi,Pe,Paux,', &
              'HICe,HICi,HALe,HALi,Sedge,Spel,Jbs,Jtot,q,Zeff,Hbre,Hcyc,', &
              'Hlin,nD+nT,pvol,<d rho>,<((d rho)^2>,V_'
         DO NR=1,NRMAX
            VC( 1)=RM(NR) !rho
            VC( 2)=RN(NR,1)*10.D0 !ne
            VC( 3)=(RN(NR,2)+RN(NR,3)+RN(NR,4)+ANC(NR)+ANFE(NR))*10.D0 !ni
            VC( 4)=ANC(NR)*10.D0 !nBe
            VC( 5)=ANFE(NR)*10.D0 !nAr
            VC( 6)=RN(NR,4)*10.D0 !nHe
!            IF(MDL_ITPA.EQ.2) THEN
               VC( 7)=AD(NR,1) !D
               VC( 8)=AV(NR,1) !V
!            ELSE
!               VC( 7)=AD(NR,2) !D
!               VC( 8)=AV(NR,2) !V
!            END IF
            VC( 9)=AK(NR,2) !Chi
            VC(10)=RT(NR,1) !Te
            VC(11)=RT(NR,2) !Ti
            VC(12)=(PIN(NR,2)+PIN(NR,3)+PIN(NR,4))*1.D-6 !Pi
            VC(13)=PIN(NR,1)*1.D-6 !Pe
            VC(14)=(POH(NR)+PNB(NR)+PNF(NR) &
                   +PEX(NR,1)+PEX(NR,2)+PEX(NR,3)+PEX(NR,4) &
                   +PRF(NR,1)+PRF(NR,2)+PRF(NR,3)+PRF(NR,4))*1.D-6 ! Paux
            VC(15)= PRF(NR,1)*1.D-6                         ! HICe
            VC(16)=(PRF(NR,2)+PRF(NR,3)+PRF(NR,4))*1.D-6    ! HICi
            VC(17)= PFCL(NR,1)*1.D-6                        ! HALe
            VC(18)=(PFCL(NR,2)+PFCL(NR,3)+PFCL(NR,4))*1.D-6 ! HALi
!            VC(19)= SPT(NR)*1.D1                            ! Sedge
!            VC(20)= SPL(NR)*1.D1                            ! Spel
            VC(19)= 0.D0
            VC(20)= 0.D0
            VC(21)= AJBS(NR)*1.D-6                          ! Jbs
            VC(22)= AJ(NR)*1.D-6                            ! Jtot
            VC(23)= QP(NR)                                  ! q
            VC(24)= ZEFF(NR)                                ! Zeff
            VC(25)= PRB(NR)*1.D-6                           ! Hbre
            VC(26)= PRC(NR)*1.D-6                           ! Hcyc
            VC(27)= PRL(NR)*1.D-6                           ! Hlin
            VC(28)= (RN(NR,2)+RN(NR,3))*10.D0               ! nD
            VC(29)= PVOLRHOG(NR)                            ! pvol
!            VC(30)= AR1RHOG(NR)*RHOTA                       ! <(nabla rho)>
!            VC(31)= AR2RHOG(NR)*RHOTA**2                    ! <(nabla rho)^2>
            VC(30)=0.D0
            VC(31)=0.D0
            VC(32)=0.D0
 !           IF(MDL_ITPA.EQ.2) THEN
 !              VC(32)= AV(NR,1)*AR2RHOG(NR)*RHOTA/AR1RHOG(NR)  ! V var
 !           ELSE
 !              VC(32)= AV(NR,2)*AR2RHOG(NR)*RHOTA/AR1RHOG(NR)  ! V var
 !           END IF
            WRITE(NFL,'(31(1PE13.6,","),1PE13.6)') (VC(N),N=1,32)
         END DO
         GO TO 1
      END IF
            
!     ----- csv file for Nopparit work on HT6M -----

      IF(KID.EQ.'CN') THEN
         WRITE(NFL,'(2A)') &
              ',10^20/m^3,10^20/m^3,keV,keV,T,T,V/m,V/m,m/s,m/s,A/m^2,,',&
              'W/m^3,W/m^3'
         WRITE(NFL,'(A)') &
              'rho,ne,ni,Te,Ti,BT,BP,EZOH,ER,VTOR,VPOL,AJOH,QP,POH,PNB'
         DO NR=1,NRMAX
            VC( 1)=RM(NR) !rho
            VC( 2)=RN(NR,1) !ne
            VC( 3)=RN(NR,2) !ni
            VC( 4)=RT(NR,1) !Te
            VC( 5)=RT(NR,2) !Ti
            VC( 6)=BB
            VC( 7)=BP(NR)
            VC( 8)=EZOH(NR)
            VC( 9)=ER(NR)
            VC(10)=VTOR(NR)
            VC(11)=VPOL(NR)
            VC(12)=AJOH(NR)
            VC(13)=QP(NR)
            VC(14)=POH(NR)
            VC(15)=PNB(NR)
            WRITE(NFL,'(14(1PE13.6,","),1PE13.6)') (VC(N),N=1,15)
         END DO
         GO TO 1
      END IF
            
!     ----- Default for gt,rt,rn,rg -----

      READ(LINE(3:),*,END=1,ERR=1) ID
      SELECT CASE(KID)
!     ----- history of global data -----
      CASE('GT')
         IF(ID.EQ.0) then
            DO I=1,NCTM
               CALL TRF1DGT(NFL,GT,GVT,NGT,I)
            ENDDO
         ELSE
            CALL TRF1DGT(NFL,GT,GVT,NGT,ID)
         ENDIF
!     ----- history of profile data -----
      CASE('RT')
         IF(ID.EQ.0) then
            DO I=1,NCRTM
               CALL TRF2DRT(NFL,GRM,GT,GVRT,NRMAX,NGT,I)
            ENDDO
         ELSE
            CALL TRF2DRT(NFL,GRM,GT,GVRT,NRMAX,NGT,ID)
         ENDIF
!     ----- snapshot of profile data -----
      CASE('RN')
         IF(ID.EQ.0) then
            DO I=1,NCRTM
               CALL TRF2DRN(NFL,GRM,GVRT,NRMAX,NGT,I)
            ENDDO
         ELSE
            CALL TRF2DRN(NFL,GRM,GVRT,NRMAX,NGT,ID)
         ENDIF
!     ----- history of profile data -----
      CASE('RG')
         IF(ID.EQ.0) then
            DO I=1,NCGM
               CALL TRF2DRG(NFL,GRM,GT,GVR,NRMAX,NGR,I)
            ENDDO
         ELSE
            CALL TRF2DRG(NFL,GRM,GT,GVR,NRMAX,NGR,ID)
         ENDIF
      END SELECT
      GOTO 1

9000  RETURN
      END SUBROUTINE TRFOUT

!     ===== 1D file output ====

      SUBROUTINE TRF1DGT(NFL,GT,GF,NTMAX,ID)
      USE TRCOMM,ONLY: NTM,NCTM
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NFL,NTMAX,ID
      REAL(4),DIMENSION(NTMAX),INTENT(IN):: GT
      REAL(4),DIMENSION(NTM,NCTM),INTENT(IN):: GF
      CHARACTER(LEN=4):: KSTR
      INTEGER:: NT

      CALL GETKGT(ID,KSTR)
      WRITE(NFL,'(I4,A,A,A)') ID,':',KSTR,'(t)'
      WRITE(NFL,'(A,I8)') 'DIM=',1
      WRITE(NFL,'(A,I8)') 'NUM=',NTMAX
      WRITE(NFL,'(1P2E15.7)') (GT(NT),GF(NT,ID),NT=1,NTMAX)
      IF(NFL.NE.6) WRITE(6,'(I4,A,A,A,I6,A)') ID,':',KSTR,'(',NTMAX,'): fout'
      RETURN
      END SUBROUTINE TRF1DGT

!     ===== 2D file output ====

      SUBROUTINE TRF2DRT(NFL,GR,GT,GF,NRMAX,NTMAX,ID)
      USE TRCOMM,ONLY: NCRTM,NTM
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NFL,NRMAX,NTMAX,ID
      REAL(4),DIMENSION(NRMAX),INTENT(IN):: GR
      REAL(4),DIMENSION(NTMAX),INTENT(IN):: GT
      REAL(4),DIMENSION(NRMAX,NTM,NCRTM),INTENT(IN):: GF
      CHARACTER(LEN=4):: KSTR
      INTEGER:: NR,NT
      
      CALL GETKRT(ID,KSTR)
      WRITE(NFL,'(I4,A,A,A)') ID,':',KSTR,'(r,t)'
      WRITE(NFL,'(A,I8)') 'DIM=',2
      WRITE(NFL,'(A,2I8)') 'NUM=',NRMAX,NTMAX
      WRITE(NFL,'(1P5E15.7)') (GR(NR),NR=1,NRMAX)
      WRITE(NFL,'(1P5E15.7)') (GT(NT),NT=1,NTMAX)
      DO NT=1,NTMAX
         WRITE(NFL,'(1P5E15.7)') (GF(NR,NT,ID),NR=1,NRMAX)
      ENDDO
      IF(NFL.NE.6) WRITE(6,'(I4,A,A,A,I4,A,I6,A)') &
           ID,':',KSTR,'(',NRMAX,',',NTMAX,'): fout'
      RETURN
      END SUBROUTINE TRF2DRT

!     ===== 2D file last data output ====

      SUBROUTINE TRF2DRN(NFL,GR,GF,NRMAX,NTMAX,ID)
      USE TRCOMM,ONLY: NCRTM,NTM
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NFL,NRMAX,NTMAX,ID
      REAL(4),DIMENSION(NRMAX),INTENT(IN):: GR
      REAL(4),DIMENSION(NRMAX,NTM,NCRTM),INTENT(IN):: GF
      CHARACTER(LEN=4):: KSTR
      INTEGER:: NR

      CALL GETKRT(ID,KSTR)
      WRITE(NFL,'(I4,A,A,A)') ID,':',KSTR,'(r,tmax)'
      WRITE(NFL,'(A,I8)') 'DIM=',1
      WRITE(NFL,'(A,I8)') 'NUM=',NRMAX
      DO NR=1,NRMAX
         WRITE(NFL,'(1P2E15.7)') GR(NR),GF(NR,NTMAX,ID)
      ENDDO
      IF(NFL.NE.6) WRITE(6,'(I4,A,A,A,I4,A)') &
           ID,':',KSTR,'(',NRMAX,'): fout'
      RETURN
      END SUBROUTINE TRF2DRN

!     ===== 2D file output ====

      SUBROUTINE TRF2DRG(NFL,GR,GT,GF,NRMAX,NTMAX,ID)
      USE TRCOMM,ONLY: NCGM,NGM
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NFL,NRMAX,NTMAX,ID
      REAL(4),DIMENSION(NRMAX),INTENT(IN):: GR
      REAL(4),DIMENSION(NTMAX),INTENT(IN):: GT
      REAL(4),DIMENSION(NRMAX+1,NGM,NCGM),INTENT(IN):: GF
      CHARACTER(LEN=4):: KSTR
      INTEGER:: NR,NT
      
      CALL GETKRT(ID,KSTR)
      WRITE(NFL,'(I4,A,A,A)') ID,':',KSTR,'(r,t)'
      WRITE(NFL,'(A,I8)') 'DIM=',2
      WRITE(NFL,'(A,2I8)') 'NUM=',NRMAX,NTMAX
      WRITE(NFL,'(1P5E15.7)') (GR(NR),NR=1,NRMAX)
      WRITE(NFL,'(1P5E15.7)') (GT(NT),NT=1,NTMAX)
      DO NT=1,NTMAX
         WRITE(NFL,'(1P5E15.7)') (GF(NR,NT,ID),NR=1,NRMAX)
      ENDDO
      IF(NFL.NE.6) WRITE(6,'(I4,A,A,A,I4,A,I6,A)') &
           ID,':',KSTR,'(',NRMAX,',',NTMAX,'): fout'
      RETURN
      END SUBROUTINE TRF2DRG

