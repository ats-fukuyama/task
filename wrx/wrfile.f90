! wrfile.f90

MODULE wrfile

  PRIVATE
  PUBLIC wr_save
  PUBLIC wr_load
  PUBLIC wr_write
  PUBLIC eccd_write

CONTAINS

!***********************************************************************
!     save ray data
!***********************************************************************

  SUBROUTINE wr_save

    USE wrcomm
    USE libfio
    IMPLICIT NONE
    INTEGER:: NFL,IERR,NRAY,I,NSTP

      NFL=21
      CALL FWOPEN(NFL,KNAMWR,0,MODEFW,'WR',IERR)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX WRSAVE: FWOPEN ERROR: IERR=',IERR
         RETURN
      END IF

      WRITE(NFL,ERR=9) NRAYMAX
      DO NRAY=1,NRAYMAX
         WRITE(NFL,ERR=9) NSTPMAX_NRAY(NRAY)
      END DO
      WRITE(6,*,ERR=9) 'NRAYMAX=',NRAYMAX
      DO NRAY=1,NRAYMAX
         WRITE(6,*,ERR=9) 'NSTPMAX=',NSTPMAX_NRAY(NRAY)
      END DO
      DO NRAY=1,NRAYMAX
         WRITE(NFL,ERR=9) (RAYIN(I,NRAY),I=1,8)
         WRITE(NFL,ERR=9) (CEXS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         WRITE(NFL,ERR=9) (CEYS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         WRITE(NFL,ERR=9) (CEZS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         WRITE(NFL,ERR=9) (RKXS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         WRITE(NFL,ERR=9) (RKYS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         WRITE(NFL,ERR=9) (RKZS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         WRITE(NFL,ERR=9) (RXS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         WRITE(NFL,ERR=9) (RYS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         WRITE(NFL,ERR=9) (RZS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         WRITE(NFL,ERR=9) (RAYRB1(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         WRITE(NFL,ERR=9) (RAYRB2(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         DO I=0,8
            WRITE(NFL,ERR=9) (RAYS(I,NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         ENDDO
         WRITE(NFL,ERR=9) (BNXS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         WRITE(NFL,ERR=9) (BNYS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         WRITE(NFL,ERR=9) (BNZS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         WRITE(NFL,ERR=9) (BABSS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
!         DO NSTP=NSTPMAX(NRAY)-10,NSTPMAX_NRAY(NRAY)
!            WRITE(6,'(A,I5,1PE12.4)') 'NSTP,BABSS=',NSTP,BABSS(NSTP,NRAY)
!         END DO
      ENDDO
      CLOSE(NFL)

      WRITE(6,*) '# DATA WAS SUCCESSFULLY SAVED TO THE FILE: ',TRIM(KNAMWR)
      RETURN

    9 WRITE(6,*) 'XX WRLOAD: File IO error detected: KNAMFR= ',TRIM(KNAMWR)
    RETURN
  END SUBROUTINE wr_save

!***********************************************************************
!     load ray data
!***********************************************************************

  SUBROUTINE wr_load(NSTAT)

    USE wrcomm
    USE wrcalpwr,ONLY: wr_calc_pwr
    USE libfio
    IMPLICIT NONE
    INTEGER,INTENT(OUT):: NSTAT
    INTEGER:: NFL,IERR,I,NRAY,NSTP,NSTPMAX_temp
    REAL(rkind):: RF,RP,ZP,PHI,RNK,ANGP,ANGT,UU
    INTEGER,ALLOCATABLE:: NTEMP(:)

      NSTAT=0

      NFL=21
      CALL FROPEN(NFL,KNAMWR,0,MODEFR,'WR',IERR)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX WRLOAX: FROPEN ERROR: IERR=',IERR
         RETURN
      END IF

      READ(NFL,END=8,ERR=9) NRAYMAX
      ALLOCATE(NTEMP(NRAYMAX))
      NSTPMAX_temp=0
      DO NRAY=1,NRAYMAX
         READ(NFL,END=8,ERR=9) NTEMP(NRAY)
         IF(NTEMP(NRAY).GT.NSTPMAX_temp) NSTPMAX_temp=NTEMP(NRAY)
      END DO
      WRITE(6,*) '## NRAYMAX,NSTPMAX=',NRAYMAX,NSTPMAX_temp
      IF(NSTPMAX.LT.NSTPMAX_temp) NSTPMAX=NSTPMAX_temp
      CALL wr_allocate
      DO NRAY=1,NRAYMAX
         NSTPMAX_NRAY(NRAY)=NTEMP(NRAY)
      END DO
      DEALLOCATE(NTEMP)
      DO NRAY=1,NRAYMAX
         READ(NFL,END=8,ERR=9) (RAYIN(I,NRAY),I=1,8)
         READ(NFL,END=8,ERR=9) (CEXS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         READ(NFL,END=8,ERR=9) (CEYS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         READ(NFL,END=8,ERR=9) (CEZS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         READ(NFL,END=8,ERR=9) (RKXS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         READ(NFL,END=8,ERR=9) (RKYS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         READ(NFL,END=8,ERR=9) (RKZS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         READ(NFL,END=8,ERR=9) (RXS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         READ(NFL,END=8,ERR=9) (RYS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         READ(NFL,END=8,ERR=9) (RZS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         READ(NFL,END=8,ERR=9) (RAYRB1(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         READ(NFL,END=8,ERR=9) (RAYRB2(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         DO I=0,8
            READ(NFL,END=8,ERR=9) (RAYS(I,NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         ENDDO
         READ(NFL,END=8,ERR=9) (BNXS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         READ(NFL,END=8,ERR=9) (BNYS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         READ(NFL,END=8,ERR=9) (BNZS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
         READ(NFL,END=8,ERR=9) (BABSS(NSTP,NRAY),NSTP=0,NSTPMAX_NRAY(NRAY))
      ENDDO
      CLOSE(NFL)

      WRITE(6,*) '# DATA WAS SUCCESSFULLY SAVED TO THE FILE: ',TRIM(KNAMWR)
      IF(RAYRB1(1,1).EQ.0.D0.AND.RAYRB2(1,1).EQ.0.D0) THEN
         NSTAT=1
      ELSE
         NSTAT=2
      ENDIF
      RF=RAYIN(1,1)
      RP=RAYIN(2,1)
      ZP=RAYIN(3,1)
      PHI=RAYIN(4,1)
      RNK=RAYIN(5,1)
      ANGP=RAYIN(6,NRAY)
      ANGT=RAYIN(7,NRAY)
      UU=RAYIN(8,NRAY)

      CALL wr_calc_pwr

      RETURN

    8 WRITE(6,*) 'XX WRLOAD: End of file detected: KNAMFR= ',TRIM(KNAMWR)
      RETURN
    9 WRITE(6,*) 'XX WRLOAD: File IO error detected: KNAMFR= ',TRIM(KNAMWR)
    RETURN
  END SUBROUTINE wr_load

!***********************************************************************
!     write ray data as ascii format
!***********************************************************************

  SUBROUTINE wr_write

    USE wrcomm
    USE libfio
    IMPLICIT NONE
    INTEGER:: NFL,IERR,NRAY,I,NSTP

    NFL=22
    CALL FWOPEN(NFL,KNAMWRW,1,MODEFW,'WR',IERR)
    IF(IERR.NE.0) THEN
       WRITE(6,*) 'XX wr_write: FWOPEN ERROR: IERR=',IERR
       RETURN
    END IF

    WRITE(NFL,'(I8)') NRAYMAX
    DO NRAY=1,NRAYMAX
       WRITE(NFL,'(8ES16.8)') (RAYIN(I,NRAY),I=1,8)
       WRITE(NFL,'(I8)') NSTPMAX_NRAY(NRAY)
       DO NSTP=1,NSTPMAX_NRAY(NRAY)
          WRITE(NFL,'(27ES16.8)') &
               (RAYS(I,NSTP,NRAY),I=1,8), &
               CEXS(NSTP,NRAY),CEYS(NSTP,NRAY),CEZS(NSTP,NRAY), &
               RKXS(NSTP,NRAY),RKYS(NSTP,NRAY),RKZS(NSTP,NRAY), &
               RXS(NSTP,NRAY),RYS(NSTP,NRAY),RZS(NSTP,NRAY), &
               RAYRB1(NSTP,NRAY),RAYRB2(NSTP,NRAY), &
               BNXS(NSTP,NRAY),BNYS(NSTP,NRAY),BNZS(NSTP,NRAY), &
               BABSS(NSTP,NRAY),RAYS(0,NSTP,NRAY)
       END DO
    END DO
    CLOSE(NFL)

    WRITE(6,*) '# DATA WAS SUCCESSFULLY WRITTEN TO THE FILE: ',TRIM(KNAMWR)
    RETURN
  END SUBROUTINE wr_write
  
!駆動電流計算  
  SUBROUTINE eccd_write

    USE wrcomm
    USE libfio
    USE plcomm !
    !USE wrexer !
    USE wrcalpwr !
    USE plprof
    IMPLICIT NONE
    REAL(rkind) :: x, y, z !
    REAL :: time1, time2 !
    REAL(rkind) :: RK, PABSN !
    TYPE(pl_prf_type),DIMENSION(NSMAX) :: plf !次元をNSMAX→NRAYMAXに変更
    INTEGER:: NFL,IERR,NRAY,I,NSTP
    INTEGER :: NS !  
    INTEGER :: nsa, nray_exec !  

    NFL=26
    CALL FWOPEN(NFL,KNAMWRW,1,MODEFW,'WR',IERR)
    IF(IERR.NE.0) THEN
       WRITE(6,*) 'XX wr_write: FWOPEN ERROR: IERR=',IERR
       RETURN
    END IF
          
    WRITE(NFL,'(I8)') NRAYMAX
    
    !吸収パワー
    !nsamax_dp = nsamax_wr !
    !CALL GUTIME(TIME1) !
    DO NRAY=1,NRAYMAX
       !nray_exec = nray !
       !omega = 2.D6*PI*RFIN(nray)
       !rkv = omega / VC
       !rnv = VC/omega
       !CALL wr_setup_start_point(NRAY,RAYS(0,0,NRAY),nstp,IERR)
       !nstpmax_nray(nray)=nstp
       !IF (IERR.NE.0) CYCLE
       !CALL wr_exec_single_ray(NRAY,RAYS(0,0,NRAY),nstp,IERR)
       !nstpmax_nray(nray)=nstp
       !IF (IERR.NE.0) CYCLE

       !DO nsa=1,nsmax_wr
          !DO nstp=0,nstpmax_nray(nray)
             !pwr_nsa_nstp_nray(nsa,nstp,nray)=pwr_nsa_nstp(nsa,nstp,nray)
             !WRITE(NFL,'(2I8,F12.4)') nstp, nray, pwr_nsa_nstp_nray(nsa,nstp,nray)
          !END DO
       !END DO 
       !=========
       
       !WRITE(NFL,'(8ES16.8)') (RAYIN(I,NRAY),I=1,8)
       WRITE(NFL,'(I8)') NSTPMAX_NRAY(NRAY)
       DO NSTP=1,NSTPMAX_NRAY(NRAY)
       
          !温度と密度用
          x = RXS(NSTP, NRAY) !
          y = RYS(NSTP, NRAY) !
          z = RZS(NSTP, NRAY) !
          CALL pl_prof3d(x,y,z,plf) !      
       
          WRITE(NFL,'(11ES16.8)') &
               !(RAYS(I,NSTP,NRAY),I=1,8), &
               RKZS(NSTP,NRAY), &
               RXS(NSTP,NRAY),RYS(NSTP,NRAY),RZS(NSTP,NRAY), &
               BNXS(NSTP,NRAY),BNYS(NSTP,NRAY),BNZS(NSTP,NRAY), &
               BABSS(NSTP,NRAY), &
               plf(1)%RN, plf(1)%RTPR, plf(1)%RTPP
               
       END DO
    END DO
    CLOSE(NFL)

    WRITE(6,*) '# DATA WAS SUCCESSFULLY WRITTEN TO THE FILE: ',TRIM(KNAMWR)
    RETURN
  END SUBROUTINE eccd_write
  

END MODULE wrfile
