! trfout.f90

MODULE trfout

  PRIVATE
  PUBLIC tr_fout

CONTAINS

!     ***********************************************************

!          FILE Output PROFILE DATA

!     ***********************************************************

      SUBROUTINE tr_fout

      USE TRCOMM
      USE libchar
      USE libfio
      IMPLICIT NONE
      CHARACTER(LEN=80):: LINE
      CHARACTER(LEN=80),SAVE:: FILENAME='stdout'
      INTEGER,SAVE:: NFL=6
      CHARACTER(LEN=2):: KID
      CHARACTER(LEN=8):: KNCTM,KNCGM,KNCRTM
      CHARACTER(LEN=30):: KFORM
      INTEGER:: I,ID,IERR,N,NR,NTL
      REAL(rkind):: VCL(32)
      REAL,DIMENSION(:,:),ALLOCATABLE:: GF

      WRITE(6,'(A,I3,A,A)') 'NFL=',NFL, &
                            '  FILENAME=',FILENAME(1:LEN_TRIM(FILENAME))
1     CALL KKINT(NCTM,KNCTM)
      CALL KKINT(NCGM,KNCGM)
      CALL KKINT(NCRTM,KNCRTM)
      WRITE(6,'(13A)') &
           '# FOUT: GT[0-',KNCTM(1:LEN_TRIM(KNCTM)),'], ', &
                   'RT[0-',KNCRTM(1:LEN_TRIM(KNCRTM)),'], ', &
                   'RN[0-',KNCRTM(1:LEN_TRIM(KNCRTM)),'], ', &
                   'RG[0-',KNCGM(1:LEN_TRIM(KNCGM)),']:txt, ', &
                   'AR:all,'
      WRITE(6,'(4A)') &
           '        CR[1-',KNCRTM(1:LEN_TRIM(KNCRTM)),'], ', &
                   'AC for all CR,CP,CN:csv,F:filename,?:help,X:exit'
      READ(5,'(A80)',END=9000,ERR=1) LINE
      KID=LINE(1:2)
      CALL toupper(KID(1:1))
      CALL toupper(KID(2:2))
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
      IF(KID.EQ.'AR') THEN
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
            VCL( 1)=RM(NR) !rho
            VCL( 2)=RN(NR,1)*10.D0 !ne
            VCL( 3)=(RN(NR,2)+RN(NR,3)+RN(NR,4)+ANC(NR)+ANFE(NR))*10.D0 !ni
            VCL( 4)=ANC(NR)*10.D0 !nBe
            VCL( 5)=ANFE(NR)*10.D0 !nAr
            VCL( 6)=RN(NR,4)*10.D0 !nHe
!            IF(MDL_ITPA.EQ.2) THEN
               VCL( 7)=AD(NR,1) !D
               VCL( 8)=AV(NR,1) !V
!            ELSE
!               VCL( 7)=AD(NR,2) !D
!               VCL( 8)=AV(NR,2) !V
!            END IF
            VCL( 9)=AK(NR,2) !Chi
            VCL(10)=RT(NR,1) !Te
            VCL(11)=RT(NR,2) !Ti
            VCL(12)=(PIN(NR,2)+PIN(NR,3)+PIN(NR,4))*1.D-6 !Pi
            VCL(13)=PIN(NR,1)*1.D-6 !Pe
            VCL(14)=(POH(NR)+PNB_NR(NR)+PNF_NR(NR) &
                   +PEX(NR,1)+PEX(NR,2)+PEX(NR,3)+PEX(NR,4) &
                   +PRF(NR,1)+PRF(NR,2)+PRF(NR,3)+PRF(NR,4))*1.D-6 ! Paux
            VCL(15)= PRF(NR,1)*1.D-6                         ! HICe
            VCL(16)=(PRF(NR,2)+PRF(NR,3)+PRF(NR,4))*1.D-6    ! HICi
            VCL(17)= PFCL(NR,1)*1.D-6                        ! HALe
            VCL(18)=(PFCL(NR,2)+PFCL(NR,3)+PFCL(NR,4))*1.D-6 ! HALi
!            VCL(19)= SPT(NR)*1.D1                            ! Sedge
!            VCL(20)= SPL(NR)*1.D1                            ! Spel
            VCL(19)= 0.D0
            VCL(20)= 0.D0
            VCL(21)= AJBS(NR)*1.D-6                          ! Jbs
            VCL(22)= AJ(NR)*1.D-6                            ! Jtot
            VCL(23)= QP(NR)                                  ! q
            VCL(24)= ZEFF(NR)                                ! Zeff
            VCL(25)= PRB(NR)*1.D-6                           ! Hbre
            VCL(26)= PRC(NR)*1.D-6                           ! Hcyc
            VCL(27)= PRL(NR)*1.D-6                           ! Hlin
            VCL(28)= (RN(NR,2)+RN(NR,3))*10.D0               ! nD
            VCL(29)= PVOLRHOG(NR)                            ! pvol
!            VCL(30)= AR1RHOG(NR)*RHOTA                       ! <(nabla rho)>
!            VCL(31)= AR2RHOG(NR)*RHOTA**2                    ! <(nabla rho)^2>
            VCL(30)=0.D0
            VCL(31)=0.D0
            VCL(32)=0.D0
 !           IF(MDL_ITPA.EQ.2) THEN
 !              VCL(32)= AV(NR,1)*AR2RHOG(NR)*RHOTA/AR1RHOG(NR)  ! V var
 !           ELSE
 !              VCL(32)= AV(NR,2)*AR2RHOG(NR)*RHOTA/AR1RHOG(NR)  ! V var
 !           END IF
            WRITE(NFL,'(31(1PE13.6,","),1PE13.6)') (VCL(N),N=1,32)
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
            VCL( 1)=RM(NR) !rho
            VCL( 2)=RN(NR,1) !ne
            VCL( 3)=RN(NR,2) !ni
            VCL( 4)=RT(NR,1) !Te
            VCL( 5)=RT(NR,2) !Ti
            VCL( 6)=BB
            VCL( 7)=BP(NR)
            VCL( 8)=EZOH(NR)
            VCL( 9)=ER(NR)
            VCL(10)=VTOR(NR)
            VCL(11)=VPOL(NR)
            VCL(12)=AJOH(NR)
            VCL(13)=QP(NR)
            VCL(14)=POH(NR)
            VCL(15)=PNB_NR(NR)
            WRITE(NFL,'(14(1PE13.6,","),1PE13.6)') (VCL(N),N=1,15)
         END DO
         GO TO 1
      END IF
            
!     ----- csv file for profile data -----

      IF(KID.EQ.'CR') THEN
         ALLOCATE(GF(NRMAX,NGT))
         READ(LINE(3:),*,END=1,ERR=1) ID
         DO NTL=1,NGT
            DO NR=1,NRMAX
               GF(NR,NTL)=GVRT(NR,NTL,ID)
            END DO
         END DO
         CALL TRF2DCRT(NFL,GRM,GT,GF,NRMAX,NGT,KVRT(ID))
         WRITE(6,'(A,A)') '# data saved in ',TRIM(FILENAME)
         DEALLOCATE(GF)
         GO TO 1
      END IF
            
!     ----- csv file for all profile data -----

      IF(KID.EQ.'AC') THEN
         ALLOCATE(GF(NRMAX,NGT))
         DO ID=1,NCRTM
            DO NTL=1,NGT
               DO NR=1,NRMAX
                  GF(NR,NTL)=GVRT(NR,NTL,ID)
               END DO
            END DO
            CALL TRF2DCRT(NFL,GRM,GT,GF,NRMAX,NGT,KVRT(ID))
         END DO
         WRITE(6,'(A,A)') '# data saved in ',TRIM(FILENAME)
         DEALLOCATE(GF)
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
      END SUBROUTINE tr_fout

!     ===== 1D file output ====

      SUBROUTINE TRF1DGT(NFL,GT,GF,NTLMAX,ID)
      USE TRCOMM,ONLY: NTM,NCTM
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NFL,NTLMAX,ID
      REAL,DIMENSION(NTLMAX),INTENT(IN):: GT
      REAL,DIMENSION(NTM,NCTM),INTENT(IN):: GF
      CHARACTER(LEN=4):: KSTR
      INTEGER:: NTL

      CALL GETKGT(ID,KSTR)
      WRITE(NFL,'(I4,A,A,A)') ID,':',KSTR,'(t)'
      WRITE(NFL,'(A,I8)') 'DIM=',1
      WRITE(NFL,'(A,I8)') 'NUM=',NTLMAX
      WRITE(NFL,'(ES15.7,",",ES15.7)') (GT(NTL),GF(NTL,ID),NTL=1,NTLMAX)
      IF(NFL.NE.6) WRITE(6,'(I4,A,A,A,I6,A)') ID,':',KSTR,'(',NTLMAX,'): fout'
      RETURN
      END SUBROUTINE TRF1DGT

!     ===== 2D file output ====

      SUBROUTINE TRF2DRT(NFL,GR,GT,GF,NRMAX,NTLMAX,ID)
      USE TRCOMM,ONLY: NCRTM,NTM
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NFL,NRMAX,NTLMAX,ID
      REAL,DIMENSION(NRMAX),INTENT(IN):: GR
      REAL,DIMENSION(NTLMAX),INTENT(IN):: GT
      REAL,DIMENSION(NRMAX,NTM,NCRTM),INTENT(IN):: GF
      CHARACTER(LEN=4):: KSTR
      INTEGER:: NR,NTL
      CHARACTER(LEN=80):: FMT
      
      CALL GETKRT(ID,KSTR)
      WRITE(NFL,'(I4,A,A,A)') ID,':',KSTR,'(r,t)'
      WRITE(NFL,'(A,I8)') 'DIM=',2
      WRITE(NFL,'(A,2I8)') 'NUM=',NRMAX,NTLMAX
      WRITE(FMT,'(A,I8,A)') '(",",ES15.7,',NRMAX,'(",",ES15.7))' 
      WRITE(NFL,FMT) (GR(NR),NR=1,NRMAX)
      WRITE(FMT,'(A,I8,A)') '(ES15.7,",",ES15.7,',NRMAX,'(",",ES15.7))'
      DO NTL=1,NTLMAX
         WRITE(NFL,FMT) GT(NTL),(GF(NR,NTL,ID),NR=1,NRMAX)
      ENDDO
      IF(NFL.NE.6) WRITE(6,'(I4,A,A,A,I4,A,I6,A)') &
           ID,':',KSTR,'(',NRMAX,',',NTLMAX,'): fout'
      RETURN
      END SUBROUTINE TRF2DRT

!     ===== 2D file last data output ====

      SUBROUTINE TRF2DRN(NFL,GR,GF,NRMAX,NTLMAX,ID)
      USE TRCOMM,ONLY: NCRTM,NTM
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NFL,NRMAX,NTLMAX,ID
      REAL,DIMENSION(NRMAX),INTENT(IN):: GR
      REAL,DIMENSION(NRMAX,NTM,NCRTM),INTENT(IN):: GF
      CHARACTER(LEN=4):: KSTR
      INTEGER:: NR

      CALL GETKRT(ID,KSTR)
      WRITE(NFL,'(I4,A,A,A)') ID,':',KSTR,'(r,tmax)'
      WRITE(NFL,'(A,I8)') 'DIM=',1
      WRITE(NFL,'(A,I8)') 'NUM=',NRMAX
      DO NR=1,NRMAX
         WRITE(NFL,'(1P2E15.7)') GR(NR),GF(NR,NTLMAX,ID)
      ENDDO
      IF(NFL.NE.6) WRITE(6,'(I4,A,A,A,I4,A)') &
           ID,':',KSTR,'(',NRMAX,'): fout'
      RETURN
      END SUBROUTINE TRF2DRN

!     ===== 2D file output ====

      SUBROUTINE TRF2DRG(NFL,GR,GT,GF,NRMAX,NTLMAX,ID)
      USE TRCOMM,ONLY: NCGM,NGM
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NFL,NRMAX,NTLMAX,ID
      REAL,DIMENSION(NRMAX),INTENT(IN):: GR
      REAL,DIMENSION(NTLMAX),INTENT(IN):: GT
      REAL,DIMENSION(NRMAX+1,NGM,NCGM),INTENT(IN):: GF
      CHARACTER(LEN=4):: KSTR
      INTEGER:: NR,NT
      
      CALL GETKRT(ID,KSTR)
      WRITE(NFL,'(I4,A,A,A)') ID,':',KSTR,'(r,t)'
      WRITE(NFL,'(A,I8)') 'DIM=',2
      WRITE(NFL,'(A,2I8)') 'NUM=',NRMAX,NTLMAX
      WRITE(NFL,'(1P5E15.7)') (GR(NR),NR=1,NRMAX)
      WRITE(NFL,'(1P5E15.7)') (GT(NT),NT=1,NTLMAX)
      DO NT=1,NTLMAX
         WRITE(NFL,'(1P5E15.7)') (GF(NR,NT,ID),NR=1,NRMAX)
      ENDDO
      IF(NFL.NE.6) WRITE(6,'(I4,A,A,A,I4,A,I6,A)') &
           ID,':',KSTR,'(',NRMAX,',',NTLMAX,'): fout'
      RETURN
      END SUBROUTINE TRF2DRG

!     ===== 2D csv output ====

      SUBROUTINE TRF2DCRT(NFL,GR,GT,GF,NRMAX,NTLMAX,TITLE)
        IMPLICIT NONE
        INTEGER,INTENT(IN):: NFL,NRMAX,NTLMAX
        REAL,INTENT(IN):: GR(NRMAX)
        REAL,INTENT(IN):: GT(NTLMAX)
        REAL,INTENT(IN):: GF(NRMAX,NTLMAX)
        CHARACTER(LEN=*),INTENT(IN):: TITLE
        CHARACTER(LEN=30):: KFORM
        INTEGER:: NR,NTL
        
        WRITE(NFL,'(A,A,I6,A,I6)') TITLE,', NRMAX=',NRMAX,', NTLMAX= ',NTLMAX
        WRITE(KFORM,'(A,I6,A)') '(A,',NTLMAX,'(",",ES14.6))'
        WRITE(6,'(A)') KFORM
        WRITE(NFL,KFORM) TITLE,(GT(NTL),NTL=1,NTLMAX)
        WRITE(KFORM,'(A,I6,A)') '(ES14.6,',NTLMAX,'(",",ES14.6))'
        WRITE(6,'(A)') KFORM
        DO NR=1,NRMAX
           WRITE(NFL,KFORM) GR(NR),(GF(NR,NTL),NTL=1,NTLMAX)
        END DO
      END SUBROUTINE TRF2DCRT
      
    END MODULE trfout
