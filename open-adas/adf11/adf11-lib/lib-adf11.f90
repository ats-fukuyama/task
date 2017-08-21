MODULE ADF11

  INTEGER, PARAMETER :: dp=SELECTED_REAL_KIND(12,100)

!     Data size required by XXDATA_11
  INTEGER,PARAMETER:: ISDIMD = 200   ! maximum number of blocks
  INTEGER,PARAMETER:: ITDIMD = 50    ! maximum number of temp values
  INTEGER,PARAMETER:: IDDIMD = 40    ! maximum number of dens values
  INTEGER,PARAMETER:: IZDIMD = 120   ! maximum number of IZ0s
  INTEGER,ALLOCATABLE,DIMENSION(:):: IZ0A,ICLASSA
  CHARACTER(LEN=256),ALLOCATABLE,DIMENSION(:):: KFNAMA
  INTEGER,ALLOCATABLE,DIMENSION(:):: IDMAXA,ITMAXA,ISMAXA,IS1MINA,IS1MAXA
  REAL(dp),ALLOCATABLE,DIMENSION(:,:):: DDENSA,DTEMPA
  REAL(dp),ALLOCATABLE,DIMENSION(:,:,:,:):: DRCOFA
  REAL(dp),ALLOCATABLE,DIMENSION(:,:,:,:,:,:):: UDRCOFA
  INTEGER:: ND_TABLE(IZDIMD,12)
  INTEGER:: NDMAX,ISMAX,IDMAX,ITMAX
  
CONTAINS

  SUBROUTINE  READ_ADF11(IERR)

    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR
    INTEGER:: LUN1,LUN2,ND,IZ0,IC,IZ,ID,IT,IS
    REAL(dp):: DDENS(IDDIMD),DTEMP(ITDIMD)
    REAL(dp):: DRCOF(ISDIMD,ITDIMD,IDDIMD)
    REAL(dp),ALLOCATABLE,DIMENSION(:):: DDENSL,DTEMPL
    REAL(dp),ALLOCATABLE,DIMENSION(:,:):: DRCOFL
    REAL(dp),ALLOCATABLE,DIMENSION(:,:,:,:):: UDRCOFL
    REAL(dp):: FX(IDDIMD,ITDIMD)
    REAL(dp):: FY(IDDIMD,ITDIMD)
    REAL(dp):: FXY(IDDIMD,ITDIMD)

    LUN1=10
    LUN2=11
    DO IC=1,12
       DO IZ=1,IZDIMD
          ND_TABLE(IZ,IC)=0
       END DO
    END DO
    IF(ALLOCATED(IZ0A))    DEALLOCATE(IZ0A)
    IF(ALLOCATED(ICLASSA)) DEALLOCATE(ICLASSA)
    IF(ALLOCATED(KFNAMA))  DEALLOCATE(KFNAMA)
    IF(ALLOCATED(IDMAXA))  DEALLOCATE(IDMAXA)
    IF(ALLOCATED(ITMAXA))  DEALLOCATE(ITMAXA)
    IF(ALLOCATED(ISMAXA))  DEALLOCATE(ISMAXA)
    IF(ALLOCATED(IS1MINA)) DEALLOCATE(IS1MINA)
    IF(ALLOCATED(IS1MAXA)) DEALLOCATE(IS1MAXA)
    IF(ALLOCATED(DDENSA))  DEALLOCATE(DDENSA)
    IF(ALLOCATED(DTEMPA))  DEALLOCATE(DTEMPA)
    IF(ALLOCATED(DRCOFA))  DEALLOCATE(DRCOFA)
    IF(ALLOCATED(UDRCOFA)) DEALLOCATE(UDRCOFA)

    CALL FROPEN(LUN1,'ADF11-FILE-LIST',1,0,'ADF11-LIST',IERR)
    IF(IERR.NE.0) THEN
       WRITE(6,*) 'XX READ_ADF11: FROPEN(LIST): IERR=',IERR
       IERR=1
       RETURN
    END IF
    READ(LUN1,*,ERR=9002,END=9003) NDMAX
    ALLOCATE(IZ0A(NDMAX),ICLASSA(NDMAX),KFNAMA(NDMAX))
    ALLOCATE(IDMAXA(NDMAX),ITMAXA(NDMAX),ISMAXA(NDMAX))
    ALLOCATE(IS1MINA(NDMAX),IS1MAXA(NDMAX))

    DO ND=1,NDMAX
       READ(LUN1,*,ERR=9004,END=9005) IZ0A(ND),ICLASSA(ND),KFNAMA(ND)
       CALL FROPEN(LUN2,KFNAMA(ND),1,0,'ADF11-DATA',IERR)
       IF(IERR.NE.0) THEN
          WRITE(6,*) 'XX READ_ADF11: FROPEN(ADF11): IERR,ND=',IERR,ND
          IERR=6
          RETURN
       END IF
       CALL FUNC_ADF11(LUN2,ICLASSA(ND),IZ0,IS1MINA(ND),IS1MAXA(ND), &
                  ISMAXA(ND),IDMAXA(ND),ITMAXA(ND),DDENS,DTEMP,DRCOF,0)
       IF(IZ0.NE.IZ0A(ND)) THEN
          WRITE(6,*) 'XX READ_ADF11: INCONSISTENT IZ0: ND,IZ0,IZ0A(ND)=', &
                      ND,IZ0,IZ0A(ND)
          IERR=7
          RETURN
       END IF
       REWIND LUN2
       CLOSE(LUN2)
       REWIND LUN2
       ND_TABLE(IZ0,ICLASSA(ND))=ND
    END DO
    CLOSE(LUN1)

    ISMAX=ISMAXA(1)
    IDMAX=IDMAXA(1)
    ITMAX=ITMAXA(1)
    DO ND=2,NDMAX
       ISMAX=MAX(ISMAX,ISMAXA(ND))
       IDMAX=MAX(IDMAX,IDMAXA(ND))
       ITMAX=MAX(ITMAX,ITMAXA(ND))
    END DO
    ALLOCATE(DDENSL(IDMAX))
    ALLOCATE(DTEMPL(ITMAX))
    ALLOCATE(DRCOFL(IDMAX,ITMAX))
    ALLOCATE(UDRCOFL(4,4,IDMAX,ITMAX))
    ALLOCATE(DDENSA(IDMAX,NDMAX))
    ALLOCATE(DTEMPA(ITMAX,NDMAX))
    ALLOCATE(DRCOFA(IDMAX,ITMAX,ISMAX,NDMAX))
    ALLOCATE(UDRCOFA(4,4,IDMAX,ITMAX,ISMAX,NDMAX))

    DO ND=1,NDMAX
       CALL FROPEN(LUN2,KFNAMA(ND),1,0,'ADF11-DATA',IERR)
       IF(IERR.NE.0) THEN
          WRITE(6,*) 'XX READ_ADF11: FROPEN(ADF11+): IERR,ND=',IERR,ND
          IERR=8
          RETURN
       END IF
       CALL FUNC_ADF11(LUN2,ICLASSA(ND),IZ0,IS1MINA(ND),IS1MAXA(ND), &
                  ISMAXA(ND),IDMAXA(ND),ITMAXA(ND),DDENS,DTEMP,DRCOF,0)
       REWIND LUN2
       CLOSE(LUN2)

       WRITE(6,'(A,A68)')   'AF11-READ: ',TRIM(KFNAMA(ND))
       WRITE(6,'(A,7I6)') 'ND,PZ0,PZMIN,PZMAX,ICLASS,IDMAX,ITMAX=', &
                        ND,IZ0,IS1MINA(ND)-1,IS1MAXA(ND)-1,ICLASSA(ND), &
                        IDMAXA(ND),ITMAXA(ND)
       DO ID=1,IDMAXA(ND)
          DDENSA(ID,ND)=DDENS(ID)
       END DO
       DO IT=1,ITMAXA(ND)
          DTEMPA(IT,ND)=DTEMP(IT)
       END DO
       DO IS=1,ISMAXA(ND)
          DO IT=1,ITMAXA(ND)
             DO ID=1,IDMAXA(ND)
                DRCOFA(ID,IT,IS,ND)=DRCOF(IS,IT,ID)
             END DO
          END DO
       END DO
    END DO

    DO ND=1,NDMAX
       DDENSL(1:IDMAXA(ND))=DDENSA(1:IDMAXA(ND),ND)
       DTEMPL(1:ITMAXA(ND))=DTEMPA(1:ITMAXA(ND),ND)
       DO IS=1,ISMAXA(ND)
          DO IT=1,ITMAXA(ND)
             FX(1         ,IT)=0.D0
             FX(IDMAXA(ND),IT)=0.D0
          END DO
          DO ID=1,IDMAXA(ND)
             FX(ID,1         )=0.D0
             FX(ID,ITMAXA(ND))=0.D0
          END DO
          FXY(1,1)=0.D0
          FXY(IDMAXA(ND),1)=0.D0
          FXY(1,ITMAXA(ND))=0.D0
          FXY(IDMAXA(ND),ITMAXA(ND))=0.D0
          DRCOFL(1:IDMAXA(ND),1:ITMAXA(ND)) &
               =DRCOFA(1:IDMAXA(ND),1:ITMAXA(ND),IS,ND)
          CALL SPL2D(DDENSL,DTEMPL,DRCOFL, &
                     FX,FY,FXY, &
                     UDRCOFL, &
                     IDMAX,IDMAXA(ND),ITMAXA(ND),3,3,IERR)
          IF(IERR.NE.0) THEN
             WRITE(6,*) 'XX READ_ADF11: SPLD2d ERROR:'
             WRITE(6,*) '   IERR,IS,ND=',IERR,IS,ND
             IERR=8
             RETURN
          END IF
          UDRCOFA(1:4,1:4,1:IDMAXA(ND),1:ITMAXA(ND),IS,ND) &
               = UDRCOFL(1:4,1:4,1:IDMAXA(ND),1:ITMAXA(ND))
       END DO
    END DO
    IERR=0
    RETURN

9002 WRITE(6,*) 'XX READ_ADF11: READ ERROR in reading NDMAX'
    IERR=2
    RETURN
9003 WRITE(6,*) 'XX READ_ADF11: END OF FILE in reading NDMAX'
    IERR=3
    RETURN
9004 WRITE(6,*) 'XX READ_ADF11: READ ERROR in reading file list: ND=',ND
    IERR=4
    RETURN
9005 WRITE(6,*) 'XX READ_ADF11: END OF FILE in reading file list: ND=',ND
    IERR=5
    RETURN
  END SUBROUTINE READ_ADF11

! --- Calculate ADF11 data for density PN and temperatrure PT ---
!        ND: Data number
!        

  SUBROUTINE CALC_ADF11(ND,IZ,PN,PT,DR,IERR)

    IMPLICIT NONE
    INTEGER,INTENT(IN):: ND,IZ
    REAL(dp),INTENT(IN):: PN,PT
    REAL(dp),INTENT(OUT):: DR
    INTEGER,INTENT(OUT):: IERR
    INTEGER:: IS
    REAL(dp):: DENS,TEMP

    IS=IZ-IS1MINA(ND)           ! IS=IZ+1 
    DENS=LOG10(PN)+20.D0-6.D0   ! 10^{20} m^{-3} -> cm^{-3}
    WRITE(6,'(A,1P3E12.4)') 'PN,LOG10(PN),DENS=',PN,LOG10(PN),DENS
    WRITE(6,'(A,1P3E12.4)') 'DDENSA(1,ND),DDENSA(IDMAXA(ND),ND)=', &
                             DDENSA(1,ND),DDENSA(IDMAXA(ND),ND)
    IF(DENS.LT.DDENSA(1,ND))          DENS=DDENSA(1,ND)
    IF(DENS.GT.DDENSA(IDMAXA(ND),ND)) DENS=DDENSA(IDMAXA(ND),ND)
    TEMP=LOG10(PT)+3.D0   ! kev -> ev
    IF(TEMP.LT.DTEMPA(1,ND))          TEMP=DTEMPA(1,ND)
    IF(TEMP.GT.DTEMPA(ITMAXA(ND),ND)) TEMP=DTEMPA(ITMAXA(ND),ND)
    WRITE(6,'(A,1P4E12.4)') 'PN,PT,DENS,TEMP   =',PN,PT,DENS,TEMP

    CALL SPL2DF(DENS,TEMP,DR, &
                DDENSA(1:IDMAXA(ND),ND), &
                DTEMPA(1:ITMAXA(ND),ND), &
                UDRCOFA(1:4,1:4,1:IDMAX,1:ITMAXA(ND),IS,ND), &
                IDMAX,IDMAXA(ND),ITMAXA(ND),IERR)
    WRITE(6,'(A,1P3E12.4)') 'DENS,TEMP,DR=',DENS,TEMP,DR
    RETURN
  END SUBROUTINE CALC_ADF11


!--- ADF11 Read routine for unresolved data ---
        
  SUBROUTINE FUNC_ADF11(IUNIT,ICLASS,IZ0,IS1MIN,IS1MAX, &
                        ISMAX,IDMAX,ITMAX,DDENS,DTEV,DRCOF,ILIST)

    IMPLICIT NONE
!-----------------------------------------------------------------------
!     Data size required by XXDATA_11
!     --------------------------------        
      INTEGER,PARAMETER:: NDPTN  = 128 ! maximum no. of partitions in one level
      INTEGER,PARAMETER:: NDPTNL = 4   ! maximum level of partitions
      INTEGER,PARAMETER:: NDPTNC = 256 ! maximum no. of components in a parti,
      INTEGER,PARAMETER:: NDCNCT = 100 ! maximum no. of elem. in connect. vec.
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     Input variables (i.e. non-parameter) required by XXDATA_11
!     ----------------------------------------------------------
      INTEGER:: IUNIT          ! unit to which input file is allocated
      INTEGER:: ICLASS         ! class of data (1 - 12 ):
!                                  1-acd, 2-scd, 3-ccd, 4-prb, 5-prc
!                                  6-qcd, 7-xcd, 8-plt, 9-pls,10-zcd
!                                 11-ycd,12-ecd
      INTEGER:: ILIST ! 0: no print out
!                       1: print data summary
!                       2: print resolved data summary
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     Output variables required by XXDATA_11
!     --------------------------------------
      INTEGER::   IZ0       , IS1MIN    , IS1MAX
      INTEGER::   IBLMX     , ISMAX     , ITMAX     , IDMAX
      INTEGER::   NPTNL     , NCNCT
      REAL(dp)::  DNR_AMS
      INTEGER::   NPTN(NDPTNL)          , NPTNC(NDPTNL,NDPTN)
      INTEGER::   IPTNLA(NDPTNL)        , IPTNA(NDPTNL,NDPTN) 
      INTEGER::   IPTNCA(NDPTNL,NDPTN,NDPTNC)
      INTEGER::   ICNCTV(NDCNCT)
      INTEGER::   ISPPR(ISDIMD)   , ISPBR(ISDIMD)   , ISSTGR(ISDIMD)
      REAL(dp)::  DDENS(IDDIMD)         , DTEV(ITDIMD)
      REAL(dp)::  DRCOF(ISDIMD,ITDIMD,IDDIMD)
      LOGICAL::   LRES    , LSTAN     , LPTN 
      CHARACTER(LEN=12):: DNR_ELE
!-----------------------------------------------------------------------

!     Variables only used by this test program
!     ----------------------------------------
      INTEGER:: I, J, K 
!-----------------------------------------------------------------------
 
!-----------------------------------------------------------------------
!     Pass unit number and dimension parameters to XXDATA_11 and
!     get back contents of file.
      CALL XXDATA_11( IUNIT  , ICLASS , &
                      ISDIMD , IDDIMD , ITDIMD , &
                      NDPTNL , NDPTN  , NDPTNC , NDCNCT , &
                      IZ0    , IS1MIN , IS1MAX , &
                      NPTNL  , NPTN   , NPTNC  , & 
                      IPTNLA , IPTNA  , IPTNCA , &
                      NCNCT  , ICNCTV , &
                      IBLMX  , ISMAX  , DNR_ELE, DNR_AMS, &
                      ISPPR  , ISPBR  , ISSTGR , &
                      IDMAX  , ITMAX  , &
                      DDENS  , DTEV   , DRCOF  , &
                      LRES   , LSTAN  , LPTN &
                     )
!-----------------------------------------------------------------------
   
      IF(ILIST.EQ.0) RETURN
      PRINT *,'XXDATA_11 Test Program'      
      PRINT *,'----------------------'
      PRINT *,' '
      PRINT *,'The following is some information about the'
      PRINT *,' '
      PRINT 1001,IZ0,IS1MIN,IS1MAX, &
                             LRES,LSTAN,LPTN, &
                             NPTNL
      IF(ILIST.GE.2) THEN
         PRINT 1002
         DO I=1,NPTNL
            PRINT 1003,I,IPTNLA(I),NPTN(I)
            DO J=1,NPTN(I)
               PRINT '(30X,I3,4X,I2,8X,I2,8X,10I3)',J,IPTNA(I,J),NPTNC(I,J), &
                              (IPTNCA(I,J,K),K=1,NPTNC(I,J))
            ENDDO
         ENDDO
      
         PRINT 1006,NCNCT
         PRINT 1007,(ICNCTV(I),I=1,NCNCT)

         PRINT 1008,IBLMX,ISMAX
         PRINT 1009,(ISPPR(I),I=1,IBLMX)
         PRINT 1010,(ISPBR(I),I=1,IBLMX)
         PRINT 1011,(ISSTGR(I),I=1,IBLMX)
      END IF

      PRINT 1012,IDMAX,ITMAX
      PRINT 1013,DDENS(1),DDENS(IDMAX)
      PRINT 1014,DTEV(1),DTEV(ITMAX)
      IF(ICLASS.LE.10) THEN
          PRINT 1015,DRCOF(1,1,1), &
                                DRCOF(IBLMX,1,1), &
                                DRCOF(1,ITMAX,1), &
                                DRCOF(1,1,IDMAX), &
                                DRCOF(IBLMX,ITMAX,IDMAX)
      ELSE 
          PRINT 1016,DRCOF(1,1,1), &
                                DRCOF(IBLMX,1,1), &
                                DRCOF(1,ITMAX,1), &
                                DRCOF(1,1,IDMAX), &
                                DRCOF(IBLMX,ITMAX,IDMAX)
      ENDIF 
      IF(ILIST.EQ.1) RETURN
         
 1001 FORMAT(7X,'IZ0    = ',I3,7X,'IS1MIN = ',I3,7X,'IS1MAX = ',I3/ &
             7X,'LRES   = ',L3,7X,'LSTAN  = ',L3,7X,'LPTN   = ',L3/ &
             7X,'NPTNL  = ',I3)
 1002 FORMAT('  Partition information:'/ &
             8X,' I IPTNLA(I) NPTN(I)', &
             2X,' J IPTNA(I,J) NPTNC(I,J)', &
             2X,'  IPTNCA(I,J,K)')
 1003 FORMAT(8X,I2,I5,I10)
 1004 FORMAT(30X,I2,I5,I11,8X,10I3/56X,10I3/56X,10I3/56X,10I3/56X,10I3) 
! 1004 FORMAT(30X,I2,I5,I11)
 1005 FORMAT(56X,I2,I5)
 1006 FORMAT('  Connection information:'/ &
             7X,'NCNCT  = ',I3,3X,'(ICNCTV(I),I=1,NCNCT)')
 1007 FORMAT(22X,20I3)
 1008 FORMAT(7X,'IBLMX  = ',I3,7X,'ISMAX  = ',I3)
 1009 FORMAT(22X,'(ISPPR(I),I=1,IBLMX)'/(22X,20I3))
 1010 FORMAT(22X,'(ISPBR(I),I=1,IBLMX)'/(22X,20I3))
 1011 FORMAT(22X,'(ISSTGR(I),I=1,IBLMX)'/(22X,20I3))
 1012 FORMAT(7X,'IDMAX  = ',I3,7X,'ITMAX  = ',I3)
 1013 FORMAT(7X,'DDENS(1) = ',F10.5,3X,'DDENS(IDMAX) = ',F10.5)
 1014 FORMAT(7X,'DTEV(1)  = ',F10.5,3X,'DTEV(ITMAX)  = ',F10.5)
 1015 FORMAT(7X,'DRCOF(1,1,1)              = ',F10.5/ &
             7X,'DRCOF(IBLMX,1,1)          = ',F10.5/ &
             7X,'DRCOF(1,ITMAX,1)          = ',F10.5/ &
             7X,'DRCOF(1,1,IDMAX)          = ',F10.5/ &
             7X,'DRCOF(IBLMX,ITMAX,IDMAX)  = ',F10.5/)
 1016 FORMAT(7X,'DRCOF(1,1,1)              = ',1PD10.3/ &
             7X,'DRCOF(IBLMX,1,1)          = ',1PD10.3/ &
             7X,'DRCOF(1,ITMAX,1)          = ',1PD10.3/ &
             7X,'DRCOF(1,1,IDMAX)          = ',1PD10.3/ &
             7X,'DRCOF(IBLMX,ITMAX,IDMAX)  = ',1PD10.3/)
!-------------------------------------------------------------------
  END SUBROUTINE FUNC_ADF11

! --- save ADF11 data in a TASK-specific binary data file ---
!        ADF11 data files are specified in ADF11-FILE-LIST
!        binary data is written in ADF11-bin.data

  SUBROUTINE SAVE_ADF11_bin(IERR)

    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR
    INTEGER:: LUN,ND,ID,IT,IS,I,J

    LUN=20
    CALL FWOPEN(LUN,'ADF11-bin.data',0,0,'adf11-bin',IERR)
    IF(IERR.NE.0) THEN
       WRITE(6,*) 'XX SAVE_ADF11_bin: FWOPEN: IERR=',IERR
       IERR=1
       RETURN
    END IF
       
    WRITE(LUN) NDMAX,ISMAX,IDMAX,ITMAX
    DO ND=1,NDMAX
       WRITE(LUN) IZ0A(ND),ICLASSA(ND)
       WRITE(LUN) KFNAMA(ND)
       WRITE(LUN) IDMAXA(ND),ITMAXA(ND),ISMAXA(ND),IS1MINA(ND),IS1MAXA(ND)
       WRITE(LUN) (DDENSA(ID,ND),ID=1,IDMAXA(ND))
       WRITE(LUN) (DTEMPA(IT,ND),IT=1,ITMAXA(ND))
       WRITE(LUN) (((DRCOFA(ID,IT,IS,ND), &
                     ID=1,IDMAXA(ND)),IT=1,ITMAXA(ND)),IS=1,ISMAXA(ND))
       WRITE(LUN) (((((UDRCOFA(I,J,ID,IT,IS,ND),I=1,4),J=1,4), &
                     ID=1,IDMAXA(ND)),IT=1,ITMAXA(ND)),IS=1,ISMAXA(ND))
    END DO
    CLOSE(LUN)
    IERR=0
    RETURN
  END SUBROUTINE SAVE_ADF11_bin

! --- Load ADF11 data from a TASK-specific binary data file ---
!        binary data is loaded from ADF11-bin.data

  SUBROUTINE LOAD_ADF11_bin(IERR)

    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR
    INTEGER:: LUN,ND,ID,IT,IS,I,J,IST,IC,IZ

    IF(ALLOCATED(IZ0A))    DEALLOCATE(IZ0A)
    IF(ALLOCATED(ICLASSA)) DEALLOCATE(ICLASSA)
    IF(ALLOCATED(KFNAMA))  DEALLOCATE(KFNAMA)
    IF(ALLOCATED(IDMAXA))  DEALLOCATE(IDMAXA)
    IF(ALLOCATED(ITMAXA))  DEALLOCATE(ITMAXA)
    IF(ALLOCATED(ISMAXA))  DEALLOCATE(ISMAXA)
    IF(ALLOCATED(IS1MINA)) DEALLOCATE(IS1MINA)
    IF(ALLOCATED(IS1MAXA)) DEALLOCATE(IS1MAXA)
    IF(ALLOCATED(DDENSA))  DEALLOCATE(DDENSA)
    IF(ALLOCATED(DTEMPA))  DEALLOCATE(DTEMPA)
    IF(ALLOCATED(DRCOFA))  DEALLOCATE(DRCOFA)
    IF(ALLOCATED(UDRCOFA)) DEALLOCATE(UDRCOFA)

    LUN=20
    CALL FROPEN(LUN,'ADF11-bin.data',0,0,'adf11-bin',IERR)
    IF(IERR.NE.0) THEN
       WRITE(6,*) 'XX LOAD_ADF11_bin: FROPEN: IERR=',IERR
       IERR=1
       RETURN
    END IF
       
    READ(LUN,IOSTAT=IST,ERR=9002,END=9012) &
         NDMAX,ISMAX,IDMAX,ITMAX
    WRITE(6,'(A,4I10)') &
         'NDMAX,ISMAX,IDMAX,ITMAX=', &
          NDMAX,ISMAX,IDMAX,ITMAX
    WRITE(6,'(A)') &
         ' ND    IZ0 ICLASS  IDMAX  ITMAX  ISMAX IS1MIN IS1MAX filename'
    ALLOCATE(IZ0A(NDMAX),ICLASSA(NDMAX),KFNAMA(NDMAX))
    ALLOCATE(IDMAXA(NDMAX),ITMAXA(NDMAX),ISMAXA(NDMAX))
    ALLOCATE(IS1MINA(NDMAX),IS1MAXA(NDMAX))
    ALLOCATE(DDENSA(IDMAX,NDMAX))
    ALLOCATE(DTEMPA(ITMAX,NDMAX))
    ALLOCATE(DRCOFA(IDMAX,ITMAX,ISMAX,NDMAX))
    ALLOCATE(UDRCOFA(4,4,IDMAX,ITMAX,ISMAX,NDMAX))
    DO IC=1,12
       DO IZ=1,IZDIMD
          ND_TABLE(IZ,IC)=0
       END DO
    END DO

    DO ND=1,NDMAX
       READ(LUN,IOSTAT=IST,ERR=9003,END=9013) &
            IZ0A(ND),ICLASSA(ND)
       ND_TABLE(IZ0A(ND),ICLASSA(ND))=ND
       READ(LUN,IOSTAT=IST,ERR=9004,END=9014) &
            KFNAMA(ND)
       READ(LUN,IOSTAT=IST,ERR=9005,END=9015) &
            IDMAXA(ND),ITMAXA(ND),ISMAXA(ND),IS1MINA(ND),IS1MAXA(ND)
       WRITE(6,'(I3,7I7,A)') ND,IZ0A(ND),ICLASSA(ND),IDMAXA(ND),ITMAXA(ND), &
                        ISMAXA(ND),IS1MINA(ND),IS1MAXA(ND),TRIM(KFNAMA(ND))

       READ(LUN,IOSTAT=IST,ERR=9006,END=9016) &
            (DDENSA(ID,ND),ID=1,IDMAXA(ND))
       READ(LUN,IOSTAT=IST,ERR=9007,END=9017) &
            (DTEMPA(IT,ND),IT=1,ITMAXA(ND))
!       WRITE(6,'(A,1P2E12.4)') 'DDENSA:',DDENSA(1,ND),DDENSA(IDMAXA(ND),ND)
!       WRITE(6,'(A,1P2E12.4)') 'DTEMPA:',DTEMPA(1,ND),DTEMPA(ITMAXA(ND),ND)

       READ(LUN,IOSTAT=IST,ERR=9008,END=9018) &
            (((DRCOFA(ID,IT,IS,ND), &
                    ID=1,IDMAXA(ND)),IT=1,ITMAXA(ND)),IS=1,ISMAXA(ND))

       READ(LUN,IOSTAT=IST,ERR=9009,END=9019) &
                 (((((UDRCOFA(I,J,ID,IT,IS,ND),I=1,4),J=1,4), &
                    ID=1,IDMAXA(ND)),IT=1,ITMAXA(ND)),IS=1,ISMAXA(ND))

    END DO
    CLOSE(LUN)
    IERR=0
    RETURN

9002 WRITE(6,*) 'XX LOAD_ADF11_bin: READ FILE ERROR 2: IOSTAT=',IST
    IERR=2
    RETURN
9003 WRITE(6,*) 'XX LOAD_ADF11_bin: READ FILE ERROR 3: IOSTAT=',IST
    IERR=3
    RETURN
9004 WRITE(6,*) 'XX LOAD_ADF11_bin: READ FILE ERROR 4: IOSTAT=',IST
    IERR=4
    RETURN
9005 WRITE(6,*) 'XX LOAD_ADF11_bin: READ FILE ERROR 5: IOSTAT=',IST
    IERR=5
    RETURN
9006 WRITE(6,*) 'XX LOAD_ADF11_bin: READ FILE ERROR 6: IOSTAT=',IST
    IERR=6
    RETURN
9007 WRITE(6,*) 'XX LOAD_ADF11_bin: READ FILE ERROR 7: IOSTAT=',IST
    IERR=7
    RETURN
9008 WRITE(6,*) 'XX LOAD_ADF11_bin: READ FILE ERROR 8: IOSTAT=',IST
    IERR=8
    RETURN
9009 WRITE(6,*) 'XX LOAD_ADF11_bin: READ FILE ERROR 9: IOSTAT=',IST
    IERR=9
    RETURN
9012 WRITE(6,*) 'XX LOAD_ADF11_bin: READ END OF FILE 12: IOSTAT=',IST
    IERR=12
    RETURN
9013 WRITE(6,*) 'XX LOAD_ADF11_bin: READ END OF FILE 13: IOSTAT=',IST
    IERR=13
    RETURN
9014 WRITE(6,*) 'XX LOAD_ADF11_bin: READ END OF FILE 14: IOSTAT=',IST
    IERR=14
    RETURN
9015 WRITE(6,*) 'XX LOAD_ADF11_bin: READ END OF FILE 15: IOSTAT=',IST
    IERR=15
    RETURN
9016 WRITE(6,*) 'XX LOAD_ADF11_bin: READ END OF FILE 16: IOSTAT=',IST
    IERR=16
    RETURN
9017 WRITE(6,*) 'XX LOAD_ADF11_bin: READ END OF FILE 17: IOSTAT=',IST
    IERR=17
    RETURN
9018 WRITE(6,*) 'XX LOAD_ADF11_bin: READ END OF FILE 18: IOSTAT=',IST
    IERR=18
    RETURN
9019 WRITE(6,*) 'XX LOAD_ADF11_bin: READ END OF FILE 19: IOSTAT=',IST
    IERR=19
    RETURN
  END SUBROUTINE LOAD_ADF11_bin

END MODULE ADF11
