MODULE ADF11

  INTEGER, PARAMETER :: dp=SELECTED_REAL_KIND(12,100)

!     Data size required by XXDATA_11
  INTEGER,PARAMETER:: ISDIMD = 200   ! maximum number of blocks
  INTEGER,PARAMETER:: ITDIMD = 60    ! maximum number of temp values
  INTEGER,PARAMETER:: IDDIMD = 40    ! maximum number of dens values
  INTEGER,PARAMETER:: IZDIMD = 160   ! maximum number of IZ0s
  INTEGER,ALLOCATABLE,DIMENSION(:):: IZ0A,ICLASSA,IM0A
  CHARACTER(LEN=256),ALLOCATABLE,DIMENSION(:):: KFNAMA
  INTEGER,ALLOCATABLE,DIMENSION(:):: IDMAXA,ITMAXA,IZTOTA,IZMINA,IZMAXA
  REAL(dp),ALLOCATABLE,DIMENSION(:,:):: DDENSA,DTEMPA
  REAL(dp),ALLOCATABLE,DIMENSION(:,:,:,:):: DRCOFA
  REAL(dp),ALLOCATABLE,DIMENSION(:,:,:,:,:,:):: UDRCOFA
  REAL(dp),ALLOCATABLE,DIMENSION(:,:,:,:):: UDRCOFL
  INTEGER:: ND_TABLE(IZDIMD,12)
  INTEGER:: NDMAX,ISMAX,IDMAX,ITMAX

! NX=1:NXMAX  grid points
! IZ=IZMINA:IZMAXA=IZMINA+IZTOTA-1
! IS=1,IZTOTA
! IZ=IZMINA+IS-1
! IS=IZ-IZMINA+1

CONTAINS

  SUBROUTINE  READ_ADF11(adas_adf11_filename,IERR)

    USE libspl2d
    USE libfio
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN):: adas_adf11_filename
    INTEGER,INTENT(OUT):: IERR
    INTEGER:: LUN1,LUN2,ND,IZ0,IC,ID,IT,IS,IZ
    REAL(dp),ALLOCATABLE:: DDENS(:),DTEMP(:),DRCOF(:,:,:)
    REAL(dp),ALLOCATABLE,DIMENSION(:):: DDENSL,DTEMPL
    REAL(dp),ALLOCATABLE,DIMENSION(:,:):: DRCOFL
    REAL(dp),ALLOCATABLE:: FX(:,:),FY(:,:),FXY(:,:)
    LOGICAL:: FLAG

    LUN1=10
    LUN2=11
    DO IC=1,12
       DO IZ=1,IZDIMD
          ND_TABLE(IZ,IC)=0
       END DO
    END DO
    IF(ALLOCATED(IZ0A))    DEALLOCATE(IZ0A)
    IF(ALLOCATED(IM0A))    DEALLOCATE(IM0A)
    IF(ALLOCATED(ICLASSA)) DEALLOCATE(ICLASSA)
    IF(ALLOCATED(KFNAMA))  DEALLOCATE(KFNAMA)
    IF(ALLOCATED(IDMAXA))  DEALLOCATE(IDMAXA)
    IF(ALLOCATED(ITMAXA))  DEALLOCATE(ITMAXA)
    IF(ALLOCATED(IZTOTA))  DEALLOCATE(IZTOTA)
    IF(ALLOCATED(IZMINA)) DEALLOCATE(IZMINA)
    IF(ALLOCATED(IZMAXA)) DEALLOCATE(IZMAXA)
    IF(ALLOCATED(DDENSA))  DEALLOCATE(DDENSA)
    IF(ALLOCATED(DTEMPA))  DEALLOCATE(DTEMPA)
    IF(ALLOCATED(DRCOFA))  DEALLOCATE(DRCOFA)
    IF(ALLOCATED(UDRCOFA)) DEALLOCATE(UDRCOFA)
    IF(ALLOCATED(UDRCOFL)) DEALLOCATE(UDRCOFL)
    IF(ALLOCATED(FX)) DEALLOCATE(FX)
    IF(ALLOCATED(FY)) DEALLOCATE(FY)
    IF(ALLOCATED(FXY)) DEALLOCATE(FXY)

    CALL FROPEN(LUN1,adas_adf11_filename,1,0,'ADF11-LIST',IERR)
    IF(IERR.NE.0) THEN
       WRITE(6,*) 'XX READ_ADF11: FROPEN(LIST): IERR=',IERR
       IERR=1
       RETURN
    END IF
    READ(LUN1,*,ERR=9002,END=9003) NDMAX
    ALLOCATE(IZ0A(NDMAX),IM0A(NDMAX),ICLASSA(NDMAX),KFNAMA(NDMAX))
    ALLOCATE(IDMAXA(NDMAX),ITMAXA(NDMAX),IZTOTA(NDMAX))
    ALLOCATE(IZMINA(NDMAX),IZMAXA(NDMAX))
    ALLOCATE(DDENS(IDDIMD),DTEMP(ITDIMD))
    ALLOCATE(DRCOF(ISDIMD,ITDIMD,IDDIMD))

    DO ND=1,NDMAX
       READ(LUN1,*,ERR=9004,END=9005) IZ0A(ND),ICLASSA(ND),KFNAMA(ND)
       CALL FROPEN(LUN2,KFNAMA(ND),1,0,'ADF11-DATA',IERR)
       IF(IERR.NE.0) THEN
          WRITE(6,*) 'XX READ_ADF11: FROPEN(ADF11): IERR,ND=',IERR,ND
          IERR=6
          RETURN
       END IF
       CALL FUNC_ADF11(LUN2,ICLASSA(ND),IZ0,IZMINA(ND),IZMAXA(ND), &
                       IZTOTA(ND),IDMAXA(ND),ITMAXA(ND),DDENS,DTEMP,DRCOF,0)
       IF(IZ0A(ND).EQ.151) THEN
          IM0A(ND)=2
          IZ0A(ND)=1
       ELSE IF(IZ0A(ND).EQ.152) THEN
          IM0A(ND)=3
          IZ0A(ND)=1
       ELSE
          IM0A(ND)=0
       END IF

       IF(IZ0.NE.IZ0A(ND)) THEN
          WRITE(6,*) 'XX READ_ADF11: INCONSISTENT IZ0: ND,IZ0,IZ0A(ND)=', &
                      ND,IZ0,IZ0A(ND)
          IERR=7
          RETURN
       END IF
       REWIND LUN2
       CLOSE(LUN2)
       REWIND LUN2
       IF(IM0A(ND).EQ.2) THEN
          ND_TABLE(151,ICLASSA(ND))=ND
       ELSE IF(IM0A(ND).EQ.3) THEN
          ND_TABLE(152,ICLASSA(ND))=ND
       ELSE
          ND_TABLE(IZ0A(ND),ICLASSA(ND))=ND
       END IF
       IF(IM0A(ND).EQ.0) ND_TABLE(IZ0,ICLASSA(ND))=ND
    END DO
    CLOSE(LUN1)

    ISMAX=IZTOTA(1)
    IDMAX=IDMAXA(1)
    ITMAX=ITMAXA(1)
    DO ND=2,NDMAX
       ISMAX=MAX(ISMAX,IZTOTA(ND))
       IDMAX=MAX(IDMAX,IDMAXA(ND))
       ITMAX=MAX(ITMAX,ITMAXA(ND))
    END DO

    WRITE(6,'(A,4I10)') &
         'NDMAX,ISMAX,IDMAX,ITMAX=', &
          NDMAX,ISMAX,IDMAX,ITMAX

    ALLOCATE(DDENSL(IDMAX))
    ALLOCATE(DTEMPL(ITMAX))
    ALLOCATE(DRCOFL(IDMAX,ITMAX))
    ALLOCATE(UDRCOFL(4,4,IDMAX,ITMAX))
    ALLOCATE(DDENSA(IDMAX,NDMAX))
    ALLOCATE(DTEMPA(ITMAX,NDMAX))
    ALLOCATE(DRCOFA(IDMAX,ITMAX,ISMAX,NDMAX))
    ALLOCATE(UDRCOFA(4,4,IDMAX,ITMAX,ISMAX,NDMAX))
    ALLOCATE(FX(IDMAX,ITMAX))
    ALLOCATE(FY(IDMAX,ITMAX))
    ALLOCATE(FXY(IDMAX,ITMAX))

    DO ND=1,NDMAX
       CALL FROPEN(LUN2,KFNAMA(ND),1,0,'ADF11-DATA',IERR)
       IF(IERR.NE.0) THEN
          WRITE(6,*) 'XX READ_ADF11: FROPEN(ADF11+): IERR,ND=',IERR,ND
          IERR=8
          RETURN
       END IF
       CALL FUNC_ADF11(LUN2,ICLASSA(ND),IZ0,IZMINA(ND),IZMAXA(ND), &
                  IZTOTA(ND),IDMAXA(ND),ITMAXA(ND),DDENS,DTEMP,DRCOF,0)
       REWIND LUN2
       CLOSE(LUN2)

       DO ID=1,IDMAXA(ND)
          DDENSA(ID,ND)=DDENS(ID)
       END DO
       DO IT=1,ITMAXA(ND)
          DTEMPA(IT,ND)=DTEMP(IT)
       END DO
       DO IS=1,IZTOTA(ND)
          DO IT=1,ITMAXA(ND)
             DO ID=1,IDMAXA(ND)
                DRCOFA(ID,IT,IS,ND)=DRCOF(IS,IT,ID)
             END DO
          END DO
       END DO
    END DO

    WRITE(6,'(A)') &
       '  ND IZ0 IM0 ICL ID  IT  IS  IS1 IS1 filename'
    WRITE(6,'(A)') &
       '             ASS MAX MAX MAX MIN MAX'
    DO ND=1,NDMAX
       WRITE(6,'(9I4,1X,A)') ND,IZ0A(ND),IM0A(ND),ICLASSA(ND),IDMAXA(ND), &
             ITMAXA(ND),IZTOTA(ND),IZMINA(ND),IZMAXA(ND),TRIM(KFNAMA(ND))
    END DO

    DO ND=1,NDMAX
       DDENSL(1:IDMAXA(ND))=DDENSA(1:IDMAXA(ND),ND)
       DTEMPL(1:ITMAXA(ND))=DTEMPA(1:ITMAXA(ND),ND)
       DO IS=1,IZTOTA(ND)
          DO IT=1,ITMAXA(ND)
             FX(1         ,IT)=0.D0
             FX(IDMAXA(ND),IT)=0.D0
          END DO
          DO ID=1,IDMAXA(ND)
             FY(ID,1         )=0.D0
             FY(ID,ITMAXA(ND))=0.D0
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

    WRITE(6,'(A)') &
       ' IZ0   acd scd ccd prb prc qcd xcd plt pls zcd ycd ecd'
    DO IZ=1,IZDIMD
       FLAG=.FALSE.
       DO IC=1,12
          IF(ND_TABLE(IZ,IC).NE.0) FLAG=.TRUE.
       END DO
       IF(FLAG) WRITE(6,'(I4,2X,12I4)') IZ,(ND_TABLE(IZ,IC),IC=1,12)
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

  SUBROUTINE CALC_ADF11(ND,NZ,PN,PT,DR,IERR)

    USE libspl2d
    USE commpi
    IMPLICIT NONE
    INTEGER,INTENT(IN):: ND,NZ   ! 1 <= NZ <=NZMAX
    REAL(dp),INTENT(IN):: PN,PT  ! [n20], [keV]
    REAL(dp),INTENT(OUT):: DR    ! [n20/s]
    INTEGER,INTENT(OUT):: IERR
    INTEGER:: IS
    REAL(dp):: DENS,TEMP

    IS=NZ-IZMINA(ND)+1
!    WRITE(6,*) 'IS,IZ,IZMINA(ND)=',IS,IZ,IZMINA(ND)
    IF(IS.LT.1) IS=1
    IF(IS.GT.IZTOTA(ND)) IS=IZTOTA(ND)
    DENS=PN+20.D0-6.D0   ! 10^{20} m^{-3} -> cm^{-3}
!    WRITE(6,'(A,1P2E12.4)') 'LOG_10 PN,DENS=',PN,DENS
!    WRITE(6,'(A,1P2E12.4)') 'DDENSA(1,ND),DDENSA(IDMAXA(ND),ND)=', &
!                             DDENSA(1,ND),DDENSA(IDMAXA(ND),ND)
    IF(DENS.LT.DDENSA(1,ND))          DENS=DDENSA(1,ND)
    IF(DENS.GT.DDENSA(IDMAXA(ND),ND)) DENS=DDENSA(IDMAXA(ND),ND)
    TEMP=PT+3.D0   ! kev -> ev
    IF(TEMP.LT.DTEMPA(1,ND))          TEMP=DTEMPA(1,ND)
    IF(TEMP.GT.DTEMPA(ITMAXA(ND),ND)) TEMP=DTEMPA(ITMAXA(ND),ND)
!    WRITE(6,'(A,1P4E12.4)') 'PN,PT,DENS,TEMP   =',PN,PT,DENS,TEMP

    UDRCOFL(1:4,1:4,1:IDMAXA(ND),1:ITMAXA(ND)) &
       =UDRCOFA(1:4,1:4,1:IDMAXA(ND),1:ITMAXA(ND),IS,ND)
    CALL SPL2DF(DENS,TEMP,DR, &
                DDENSA(1:IDMAXA(ND),ND),DTEMPA(1:ITMAXA(ND),ND), &
                UDRCOFL,IDMAX,IDMAXA(ND),ITMAXA(ND),IERR)
    DR=DR+14.D0
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
      INTEGER,ALLOCATABLE::   NPTN(:)          , NPTNC(:,:)
      INTEGER,ALLOCATABLE::   IPTNLA(:)        , IPTNA(:,:) 
      INTEGER,ALLOCATABLE::   IPTNCA(:,:,:)
      INTEGER,ALLOCATABLE::   ICNCTV(:)
      INTEGER,ALLOCATABLE::   ISPPR(:)   , ISPBR(:)   , ISSTGR(:)
      REAL(dp)::  DDENS(IDDIMD)         , DTEV(ITDIMD)
      REAL(dp)::  DRCOF(ISDIMD,ITDIMD,IDDIMD)
      LOGICAL::   LRES    , LSTAN     , LPTN 
      CHARACTER(LEN=12):: DNR_ELE
!-----------------------------------------------------------------------

!     Variables only used by this test program
!     ----------------------------------------
      INTEGER:: I, J, K 
!-----------------------------------------------------------------------
 
      ALLOCATE(NPTN(NDPTNL)          , NPTNC(NDPTNL,NDPTN))
      ALLOCATE(IPTNLA(NDPTNL)        , IPTNA(NDPTNL,NDPTN))
      ALLOCATE(IPTNCA(NDPTNL,NDPTN,NDPTNC))
      ALLOCATE(ICNCTV(NDCNCT))
      ALLOCATE(ISPPR(ISDIMD)   , ISPBR(ISDIMD)   , ISSTGR(ISDIMD))

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
! 1004 FORMAT(30X,I2,I5,I11,8X,10I3/56X,10I3/56X,10I3/56X,10I3/56X,10I3) 
! 1004 FORMAT(30X,I2,I5,I11)
! 1005 FORMAT(56X,I2,I5)
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

  SUBROUTINE SAVE_ADF11_bin(adas_bin_filename,IERR)

    USE libfio
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN):: adas_bin_filename
    INTEGER,INTENT(OUT):: IERR
    INTEGER:: LUN,ND,ID,IT,IS,I,J

    LUN=20
    CALL FWOPEN(LUN,adas_bin_filename,0,0,'adf11-bin',IERR)
    IF(IERR.NE.0) THEN
       WRITE(6,*) 'XX SAVE_ADF11_bin: FWOPEN: IERR=',IERR
       IERR=1
       RETURN
    END IF
       
    WRITE(LUN) NDMAX,ISMAX,IDMAX,ITMAX
    DO ND=1,NDMAX
       WRITE(LUN) IZ0A(ND),IM0A(ND),ICLASSA(ND)
       WRITE(LUN) KFNAMA(ND)
       WRITE(LUN) IDMAXA(ND),ITMAXA(ND),IZTOTA(ND),IZMINA(ND),IZMAXA(ND)
       WRITE(LUN) (DDENSA(ID,ND),ID=1,IDMAXA(ND))
       WRITE(LUN) (DTEMPA(IT,ND),IT=1,ITMAXA(ND))
       WRITE(LUN) (((DRCOFA(ID,IT,IS,ND), &
                     ID=1,IDMAXA(ND)),IT=1,ITMAXA(ND)),IS=1,IZTOTA(ND))
!       WRITE(6,'(1P5E12.4)') (((DRCOFA(ID,IT,IS,ND), &
!                     ID=1,IDMAXA(ND)),IT=1,ITMAXA(ND)),IS=1,IZTOTA(ND))
       WRITE(LUN) (((((UDRCOFA(I,J,ID,IT,IS,ND),I=1,4),J=1,4), &
                     ID=1,IDMAXA(ND)),IT=1,ITMAXA(ND)),IS=1,IZTOTA(ND))
    END DO
    CLOSE(LUN)
    IERR=0
    RETURN
  END SUBROUTINE SAVE_ADF11_bin

! --- Load ADF11 data from a TASK-specific binary data file ---
!        binary data is loaded from ADF11-bin.data

  SUBROUTINE LOAD_ADF11_bin(adas_bin_filename,IERR)

    USE libfio
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN):: adas_bin_filename
    INTEGER,INTENT(OUT):: IERR
    INTEGER:: LUN,ND,ID,IT,IS,I,J,IST,IC,IZ
    LOGICAL:: FLAG

    IF(ALLOCATED(IZ0A))    DEALLOCATE(IZ0A)
    IF(ALLOCATED(IM0A))    DEALLOCATE(IM0A)
    IF(ALLOCATED(ICLASSA)) DEALLOCATE(ICLASSA)
    IF(ALLOCATED(KFNAMA))  DEALLOCATE(KFNAMA)
    IF(ALLOCATED(IDMAXA))  DEALLOCATE(IDMAXA)
    IF(ALLOCATED(ITMAXA))  DEALLOCATE(ITMAXA)
    IF(ALLOCATED(IZTOTA))  DEALLOCATE(IZTOTA)
    IF(ALLOCATED(IZMINA)) DEALLOCATE(IZMINA)
    IF(ALLOCATED(IZMAXA)) DEALLOCATE(IZMAXA)
    IF(ALLOCATED(DDENSA))  DEALLOCATE(DDENSA)
    IF(ALLOCATED(DTEMPA))  DEALLOCATE(DTEMPA)
    IF(ALLOCATED(DRCOFA))  DEALLOCATE(DRCOFA)
    IF(ALLOCATED(UDRCOFA)) DEALLOCATE(UDRCOFA)
    IF(ALLOCATED(UDRCOFL)) DEALLOCATE(UDRCOFL)

    LUN=20
    CALL FROPEN(LUN,adas_bin_filename,0,0,'adf11-bin',IERR)
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
       '  ND IZ0 IM0 ICL ID  IT  IS  IS1 IS1 filename'
    WRITE(6,'(A)') &
       '             ASS MAX MAX MAX MIN MAX'

    ALLOCATE(IZ0A(NDMAX),IM0A(NDMAX),ICLASSA(NDMAX),KFNAMA(NDMAX))
    ALLOCATE(IDMAXA(NDMAX),ITMAXA(NDMAX),IZTOTA(NDMAX))
    ALLOCATE(IZMINA(NDMAX),IZMAXA(NDMAX))
    ALLOCATE(DDENSA(IDMAX,NDMAX))
    ALLOCATE(DTEMPA(ITMAX,NDMAX))
    ALLOCATE(DRCOFA(IDMAX,ITMAX,ISMAX,NDMAX))
    ALLOCATE(UDRCOFA(4,4,IDMAX,ITMAX,ISMAX,NDMAX))
    ALLOCATE(UDRCOFL(4,4,IDMAX,ITMAX))

    DO IC=1,12
       DO IZ=1,IZDIMD
          ND_TABLE(IZ,IC)=0
       END DO
    END DO

    DO ND=1,NDMAX
       READ(LUN,IOSTAT=IST,ERR=9003,END=9013) &
            IZ0A(ND),IM0A(ND),ICLASSA(ND)
       IF(IM0A(ND).EQ.2) THEN
          ND_TABLE(151,ICLASSA(ND))=ND
       ELSE IF(IM0A(ND).EQ.3) THEN
          ND_TABLE(152,ICLASSA(ND))=ND
       ELSE
          ND_TABLE(IZ0A(ND),ICLASSA(ND))=ND
       END IF

       READ(LUN,IOSTAT=IST,ERR=9004,END=9014) &
            KFNAMA(ND)
       READ(LUN,IOSTAT=IST,ERR=9005,END=9015) &
            IDMAXA(ND),ITMAXA(ND),IZTOTA(ND),IZMINA(ND),IZMAXA(ND)
       WRITE(6,'(9I4,1X,A)') ND,IZ0A(ND),IM0A(ND),ICLASSA(ND),IDMAXA(ND), &
             ITMAXA(ND),IZTOTA(ND),IZMINA(ND),IZMAXA(ND),TRIM(KFNAMA(ND))

       READ(LUN,IOSTAT=IST,ERR=9006,END=9016) &
            (DDENSA(ID,ND),ID=1,IDMAXA(ND))
       READ(LUN,IOSTAT=IST,ERR=9007,END=9017) &
            (DTEMPA(IT,ND),IT=1,ITMAXA(ND))
!       WRITE(6,'(A,1P2E12.4)') 'DDENSA:',DDENSA(1,ND),DDENSA(IDMAXA(ND),ND)
!       WRITE(6,'(A,1P2E12.4)') 'DTEMPA:',DTEMPA(1,ND),DTEMPA(ITMAXA(ND),ND)

       READ(LUN,IOSTAT=IST,ERR=9008,END=9018) &
            (((DRCOFA(ID,IT,IS,ND), &
                    ID=1,IDMAXA(ND)),IT=1,ITMAXA(ND)),IS=1,IZTOTA(ND))
!       WRITE(6,'(1P5E12.4)') (((DRCOFA(ID,IT,IS,ND), &
!                     ID=1,IDMAXA(ND)),IT=1,ITMAXA(ND)),IS=1,IZTOTA(ND))
!
       READ(LUN,IOSTAT=IST,ERR=9009,END=9019) &
                 (((((UDRCOFA(I,J,ID,IT,IS,ND),I=1,4),J=1,4), &
                    ID=1,IDMAXA(ND)),IT=1,ITMAXA(ND)),IS=1,IZTOTA(ND))

    END DO
    CLOSE(LUN)

    WRITE(6,'(A)') &
       ' ICL   (1) (2) (3) (4) (5) (6) (7) (8) (9)(10)(11)(12)'
    WRITE(6,'(A)') &
       ' IZ0   acd scd ccd prb prc qcd xcd plt pls zcd ycd ecd'
    DO IZ=1,IZDIMD
       FLAG=.FALSE.
       DO IC=1,12
          IF(ND_TABLE(IZ,IC).NE.0) FLAG=.TRUE.
       END DO
       IF(FLAG) WRITE(6,'(I4,2X,12I4)') IZ,(ND_TABLE(IZ,IC),IC=1,12)
    END DO

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

! --- broadcast loaded ADF11-bin-data ---

  SUBROUTINE broadcast_ADF11_bin(IERR)

    USE libmpi
    USE commpi
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: IERR
    INTEGER:: ND,ID,IT,IS,nmax,n,IC,I,J
    REAL(dp),ALLOCATABLE:: work(:)

    CALL mtx_broadcast1_integer(NDMAX)
    CALL mtx_broadcast1_integer(ISMAX)
    CALL mtx_broadcast1_integer(IDMAX)
    CALL mtx_broadcast1_integer(ITMAX)

    IF(nrank.NE.0) THEN
       IF(ALLOCATED(IZ0A))    DEALLOCATE(IZ0A)
       IF(ALLOCATED(IM0A))    DEALLOCATE(IM0A)
       IF(ALLOCATED(ICLASSA)) DEALLOCATE(ICLASSA)
       IF(ALLOCATED(KFNAMA))  DEALLOCATE(KFNAMA)
       IF(ALLOCATED(IDMAXA))  DEALLOCATE(IDMAXA)
       IF(ALLOCATED(ITMAXA))  DEALLOCATE(ITMAXA)
       IF(ALLOCATED(IZTOTA))  DEALLOCATE(IZTOTA)
       IF(ALLOCATED(IZMINA)) DEALLOCATE(IZMINA)
       IF(ALLOCATED(IZMAXA)) DEALLOCATE(IZMAXA)
       IF(ALLOCATED(DDENSA))  DEALLOCATE(DDENSA)
       IF(ALLOCATED(DTEMPA))  DEALLOCATE(DTEMPA)
       IF(ALLOCATED(DRCOFA))  DEALLOCATE(DRCOFA)
       IF(ALLOCATED(UDRCOFA)) DEALLOCATE(UDRCOFA)
       IF(ALLOCATED(UDRCOFL)) DEALLOCATE(UDRCOFL)

       ALLOCATE(IZ0A(NDMAX),IM0A(NDMAX),ICLASSA(NDMAX),KFNAMA(NDMAX))
       ALLOCATE(IDMAXA(NDMAX),ITMAXA(NDMAX),IZTOTA(NDMAX))
       ALLOCATE(IZMINA(NDMAX),IZMAXA(NDMAX))
       ALLOCATE(DDENSA(IDMAX,NDMAX))
       ALLOCATE(DTEMPA(ITMAX,NDMAX))
       ALLOCATE(DRCOFA(IDMAX,ITMAX,ISMAX,NDMAX))
       ALLOCATE(UDRCOFA(4,4,IDMAX,ITMAX,ISMAX,NDMAX))
       ALLOCATE(UDRCOFL(4,4,IDMAX,ITMAX))
    END IF

    CALL mtx_broadcast_integer(IZ0A,NDMAX)
    CALL mtx_broadcast_integer(IM0A,NDMAX)
    CALL mtx_broadcast_integer(ICLASSA,NDMAX)
    DO ND=1,NDMAX
       CALL mtx_broadcast_character(KFNAMA(ND),256)
    END DO
    CALL mtx_broadcast_integer(IDMAXA,NDMAX)
    CALL mtx_broadcast_integer(ITMAXA,NDMAX)
    CALL mtx_broadcast_integer(IZTOTA,NDMAX)
    CALL mtx_broadcast_integer(IZMINA,NDMAX)
    CALL mtx_broadcast_integer(IZMAXA,NDMAX)

    DO ND=1,NDMAX
       ALLOCATE(work(IDMAXA(ND)))
       IF(nrank.EQ.0) THEN
          DO ID=1,IDMAXA(ND)
             work(ID)=DDENSA(ID,ND)
          END DO
       END IF
       CALL mtx_broadcast_real8(work,IDMAXA(ND))
       IF(nrank.NE.0) THEN
          DO ID=1,IDMAXA(ND)
             DDENSA(ID,ND)=work(ID)
          END DO
          DEALLOCATE(work)
          ALLOCATE(work(ITMAXA(ND)))
       ELSE
          DEALLOCATE(work)
          ALLOCATE(work(ITMAXA(ND)))
          DO IT=1,ITMAXA(ND)
             work(IT)=DTEMPA(IT,ND)
          END DO
       END IF
       CALL mtx_broadcast_real8(work,ITMAXA(ND))
       IF(nrank.NE.0) THEN
          DO IT=1,ITMAXA(ND)
             DTEMPA(IT,ND)=work(IT)
          END DO
       END IF
       DEALLOCATE(work)
    END DO

    DO ND=1,NDMAX
       nmax=IZTOTA(ND)*ITMAXA(ND)*IDMAXA(ND)
       ALLOCATE(work(nmax))
       IF(nrank.EQ.0) THEN
          n=0
          DO IS=1,IZTOTA(ND)
             DO IT=1,ITMAXA(ND)
                DO ID=1,IDMAXA(ND)
                   n=n+1
                   work(n)=DRCOFA(ID,IT,IS,ND)
                END DO
             END DO
          END DO
       END IF
       CALL mtx_broadcast_real8(work,nmax)
       IF(nrank.NE.0) THEN
          n=0
          DO IS=1,IZTOTA(ND)
             DO IT=1,ITMAXA(ND)
                DO ID=1,IDMAXA(ND)
                   n=n+1
                   DRCOFA(ID,IT,IS,ND)=work(n)
                END DO
             END DO
          END DO
       END IF
       DEALLOCATE(work)
    END DO
    DO ND=1,NDMAX
       nmax=IZTOTA(ND)*ITMAXA(ND)*IDMAXA(ND)*4*4
       ALLOCATE(work(nmax))
       IF(nrank.EQ.0) THEN
          n=0
          DO IS=1,IZTOTA(ND)
             DO IT=1,ITMAXA(ND)
                DO ID=1,IDMAXA(ND)
                   DO J=1,4
                      DO I=1,4
                         n=n+1
                         work(n)=UDRCOFA(I,J,ID,IT,IS,ND)
                      END DO
                   END DO
                END DO
             END DO
          END DO
       END IF
       CALL mtx_broadcast_real8(work,nmax)
       IF(nrank.NE.0) THEN
          n=0
          DO IS=1,IZTOTA(ND)
             DO IT=1,ITMAXA(ND)
                DO ID=1,IDMAXA(ND)
                   DO J=1,4
                      DO I=1,4
                         n=n+1
                         UDRCOFA(I,J,ID,IT,IS,ND)=work(n)
                      END DO
                   END DO
                END DO
             END DO
          END DO
       END IF
       DEALLOCATE(work)
    END DO

    DO IC=1,12
       CALL mtx_broadcast_integer(ND_TABLE(1:IZDIMD,IC),IZDIMD)
    END DO

    IERR=0
    RETURN
  END SUBROUTINE broadcast_ADF11_bin

  SUBROUTINE IONIZE_EQ2(IZ0,PN,PT,PZAV,PZ2AV,PWBAV,PWLAV,IND,IERR)

!   IND=0 No plot
!       1 Plot PNZ(PZ) 
!       2 Plot DRI(PZ),DRR(PZ),PNZ(PZ) 

    USE libgrf
    IMPLICIT NONE
    INTEGER,INTENT(IN):: IZ0,IND
    INTEGER,INTENT(OUT):: IERR
    REAL(dp),INTENT(IN):: PN,PT
    REAL(dp),INTENT(OUT):: PZAV,PZ2AV,PWBAV,PWLAV
    REAL(dp),DIMENSION(:),ALLOCATABLE:: PNZ
    INTEGER:: NXMAX,NX,IZ,NDI,NDR,NDB,NDL
    REAL(dp):: PWB,PWL

    NDI=ND_TABLE(IZ0,2)  ! ionization rate
    NDR=ND_TABLE(IZ0,1)  ! recombination rate
    NDB=ND_TABLE(IZ0,4)  ! brems power rate
    NDL=ND_TABLE(IZ0,8)  ! line power rate

    NXMAX=MAX(IZTOTA(NDI),IZTOTA(NDR))
    IF(ALLOCATED(PNZ)) DEALLOCATE(PNZ)
    ALLOCATE(PNZ(0:NXMAX))

    CALL IONIZE_EQ1(IZ0,PN,PT,PNZ,IND,IERR)

    PZAV=0.D0
    PZ2AV=0.D0
    PWBAV=0.D0
    PWLAV=0.D0
    DO NX=1,NXMAX
       IZ=NX-1+IZMINA(NDI)
       PZAV =PZAV +IZ*PNZ(NX)
       PZ2AV=PZ2AV+IZ*IZ*PNZ(NX)
       IF(NDB.EQ.0) THEN
          PWB=0.D0
       ELSE
          CALL CALC_ADF11(NDB,IZ,PN,PT,PWB,IERR)
       END IF
       PWBAV=PWBAV+10.D0**PWB*PNZ(NX)
       IF(NDL.EQ.0) THEN
          PWL=0.D0
       ELSE
          CALL CALC_ADF11(NDL,IZ,PN,PT,PWL,IERR)
       END IF
       PWLAV=PWLAV+10.D0**PWL*PNZ(NX)
    END DO

    RETURN
  END SUBROUTINE IONIZE_EQ2

  SUBROUTINE IONIZE_EQ1(IZ0,PN,PT,PNZ,IND,IERR)

!   IND=0 No plot
!       1 Plot PNZ(PZ) 
!       2 Plot DRI(PZ),DRR(PZ),PNZ(PZ) 

    USE libbnd
    USE libgrf
    use, intrinsic :: ieee_arithmetic, only : IEEE_SELECTED_REAL_KIND
    use, intrinsic :: ieee_exceptions
    IMPLICIT NONE
    
    INTEGER,INTENT(IN):: IZ0,IND
    INTEGER,INTENT(OUT):: IERR
    REAL(dp),INTENT(IN):: PN,PT
    REAL(dp),ALLOCATABLE,DIMENSION(:), INTENT(OUT):: PNZ
    INTEGER:: NDI,NDR,NXMAX,NX,IZ
    REAL(dp):: DR,PNZTOT,RATIO
    REAL(dp),ALLOCATABLE,DIMENSION(:):: DRI,DRR,XDATA,CVEC
    REAL(dp),ALLOCATABLE,DIMENSION(:,:):: FDATA,CMAT
!    REAL(dp),ALLOCATABLE,DIMENSION(:,:):: CMAT_SAVE
    LOGICAL:: should_halt, was_flagged
    
    NDI=ND_TABLE(IZ0,2)  ! ionization rate
    NDR=ND_TABLE(IZ0,1)  ! recombination rate

    NXMAX=IZTOTA(NDI)
!    WRITE(6,'(A,3I5)') &
!         'IZTOTA,IZMINA,IZMAXA=', &
!          IZTOTA(NDI),IZMINA(NDI),IZMAXA(NDI)
    IF(ALLOCATED(DRI)) DEALLOCATE(DRI)
    IF(ALLOCATED(DRR)) DEALLOCATE(DRR)
    IF(ALLOCATED(PNZ)) DEALLOCATE(PNZ)
    IF(ALLOCATED(XDATA)) DEALLOCATE(XDATA)
    IF(ALLOCATED(FDATA)) DEALLOCATE(FDATA)
    ALLOCATE(DRI(NXMAX),DRR(NXMAX),PNZ(0:NXMAX))
    ALLOCATE(XDATA(NXMAX+1),FDATA(NXMAX+1,2))

    DO NX=1,NXMAX
       IZ=NX-1+IZMINA(NDI)
       CALL CALC_ADF11(NDI,IZ,PN,PT,DR,IERR)
       DRI(NX)=10.D0**DR
       IF(IERR.NE.0) THEN
          WRITE(6,*) 'XX test-adf11: ionize_eq1: IERR=',IERR
          GOTO 9000
       END IF
       IZ=NX-1+IZMINA(NDR)
       CALL CALC_ADF11(NDR,IZ,PN,PT,DR,IERR)
       DRR(NX)=10.D0**DR
       IF(IERR.NE.0) THEN
          WRITE(6,*) 'XX test-adf11: ionize_eq2: IERR=',IERR
          GOTO 9000
       END IF
    END DO

    IF(NXMAX.EQ.1) THEN
       PNZ(0)=DRR(1)/(DRR(1)+DRI(1))
       PNZ(1)=DRI(1)/(DRR(1)+DRI(1))
       IERR=0
       RETURN       
    END IF

    IF(IND.EQ.2) THEN
       DO NX=1,NXMAX
          XDATA(NX)=DBLE(NX)
          FDATA(NX,1)=LOG10(DRI(NX))
          FDATA(NX,2)=LOG10(DRR(NX))
          write(6,'(A,I5,1P4E12.4)') 'DR:',NX,DRI(NX),DRR(NX)
       END DO

       CALL PAGES
       CALL GRD1D(0,XDATA,FDATA,NXMAX+1,NXMAX,2,'@DRI,DRR vs Z@',2)
       CALL PAGEE
    END IF

    ALLOCATE(CMAT(3,NXMAX+1),CVEC(NXMAX+1))
    CMAT(1,1)= 0.D0
    CMAT(2,1)=-DRI(1)
    CMAT(3,1)= DRR(1)
    DO NX=2,NXMAX
       CMAT(1,NX)=  DRI(NX-1)
       CMAT(2,NX)=-(DRI(NX)+DRR(NX-1))
       CMAT(3,NX)=  DRR(NX)
    END DO
    CMAT(1,NXMAX+1)=  DRI(NXMAX)
    CMAT(2,NXMAX+1)= -DRR(NXMAX)
    CMAT(3,NXMAX+1)=  0.D0

    DO NX=1,NXMAX+1
       CVEC(NX)=0.D0
    END DO

    NX=1
    DO IZ=1,NXMAX
       RATIO=DRR(IZ)/DRI(IZ)
       IF(RATIO.LT.1.D0) NX=IZ
    END DO

    CMAT(1,NX)=0.D0
    CMAT(2,NX)=1.D0
    CMAT(3,NX)=0.D0
    CVEC(NX)=1.D0

!    CMAT_SAVE=CMAT

    CALL BANDRD(CMAT,CVEC,NXMAX+1,3,3,IERR)

!    CMAT=CMAT_SAVE

    DO IZ=0,NXMAX
       PNZ(IZ)=CVEC(IZ+1)
    END DO
!       IF(IZ.LT.NXMAX) THEN
!          RATIO=DRR(IZ)/DRI(IZ)
!       ELSE
!          RATIO=0.D0
!       END IF
!       write(6,'(A,I5,1P2E12.4)') 'IZ,PNZ,R/I=',IZ,PNZ(IZ),RATIO

    DO IZ=0,NXMAX
       IF(IZ.EQ.0) THEN
          CVEC(IZ+1)=CMAT(2,IZ+1)*PNZ(IZ) &
                    +CMAT(3,IZ+1)*PNZ(IZ+1)
       ELSE IF(IZ.EQ.NXMAX) THEN
          CVEC(IZ+1)=CMAT(1,IZ+1)*PNZ(IZ-1) &
                    +CMAT(2,IZ+1)*PNZ(IZ)
       ELSE
          CVEC(IZ)=CMAT(1,IZ)*PNZ(IZ-1) &
                  +CMAT(2,IZ)*PNZ(IZ) &
                  +CMAT(3,IZ)*PNZ(IZ+1)
       END IF
!        write(6,'(A,I5,1P7E10.2)') 'EQ: ',IZ, &
!              CMAT(1,IZ+1),CMAT(2,IZ+1),CMAT(3,IZ+1),PNZ(IZ),CVEC(IZ+1), &
!              DRI(MIN(IZ,NXMAX-1)),DRR(MIN(IZ,NXMAX-1))
    END DO
    DEALLOCATE(CMAT,CVEC)
!    DEALLOCATE(CMAT_SAVE)

    PNZTOT=0.D0
    DO IZ=0,NXMAX
       PNZTOT=PNZTOT+PNZ(IZ)
    END DO

    ! Get the original halting mode and signal state
    call ieee_get_halting_mode(IEEE_UNDERFLOW, should_halt)
    call ieee_get_flag(IEEE_UNDERFLOW, was_flagged)

    ! Ensure we aren't going to halt on underflow
    call ieee_set_halting_mode(IEEE_UNDERFLOW, .FALSE.)
    
     DO IZ=0,NXMAX
!        WRITE(6,'(A,1P2E12.4)') 'PNZ:',PNZ(IZ),PNZTOT
       IF(PNZ(IZ).LT.1.D-8) THEN
          PNZ(IZ)=1.D-8/PNZTOT
       ELSE
          PNZ(IZ)=PNZ(IZ)/PNZTOT
       END IF
    END DO

    ! And restore our old state
    call ieee_set_halting_mode(IEEE_UNDERFLOW, should_halt)
    call ieee_set_flag(IEEE_UNDERFLOW, was_flagged)

    IF(IND.GE.1) THEN    
       DO NX=0,NXMAX
          XDATA(NX+1)=DBLE(NX)
          FDATA(NX+1,1)=MAX(PNZ(NX),1.D-8)
       END DO

       CALL PAGES
!       CALL GRD1D(0,XDATA,FDATA,NXMAX+1,NXMAX+1,1,'@PNZ vs Z@',2)
       CALL GRD1D(0,XDATA,FDATA,NXMAX+1,NXMAX+1,1,'@PNZ vs Z@',0)
       CALL PAGEE
    END IF

9000 CONTINUE
    RETURN
  END SUBROUTINE IONIZE_EQ1

END MODULE ADF11
