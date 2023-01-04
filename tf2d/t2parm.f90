!C
!C Input parapmeter interface
!C
MODULE T2PARM
  
  USE T2CNST, ONLY: ikind,rkind
  
  PRIVATE
  PUBLIC T2_PARM,T2_VIEW
  
CONTAINS

  !C
  !C 
  !C

  SUBROUTINE T2_PARM(MODE,KIN,IERR)
!     MODE=0 : standard namelinst input
!     MODE=1 : namelist file input
!     MODE=2 : namelist line input

!     IERR=0 : normal end
!     IERR=1 : namelist standard input error
!     IERR=2 : namelist file does not exist
!     IERR=3 : namelist file open error
!     IERR=4 : namelist file read error
!     IERR=5 : namelist file abormal end of file
!     IERR=6 : namelist line input error
!     IERR=7 : unknown MODE
!     IERR=10X : input parameter out of range

    USE libkio
    IMPLICIT NONE
    INTEGER,INTENT(IN):: MODE
    CHARACTER(LEN=*),INTENT(IN):: KIN
    INTEGER,INTENT(OUT):: IERR

1   CALL TASK_PARM(MODE,'T2',KIN,T2_NLIN,T2_PLST,IERR)
    IF(IERR.NE.0 .AND. IERR.NE.2) RETURN

    CALl T2_CHECK(IERR)
    IF(MODE.EQ.0.AND.IERR.NE.0) GOTO 1
    IF(IERR.NE.0) IERR=IERR+100

    RETURN
  END SUBROUTINE T2_PARM

!     ****** INPUT NAMELIST ******

  SUBROUTINE T2_NLIN(NID,IST,IERR)


    USE T2COMM,ONLY:&
         c10rname,i0fnum,&
         NNMAX,NQMAX,NDMAX,NSMAX,NPMIN,NLMAX,&
         !
         CoordinateSwitch,TestCase,EqSet,&
         !
         i1mlvl,i1rdn2,d1rec,d0rw,RR,RA,&
         i0nm,i0nn,i0tm,i0tn,d0qc,d0qs,d0bc,&
         d1nc,d1ns,d1nw,d1tc,d1ts,d1tw,Pa,Pz,&
         ! AF
         dt,time_init,eps_conv, &
         ntmax,ntstep,nt0dmax,nt0dstep,nt2dmax,nt2dstep,nconvmax, &
         idfile,idprint,idplot,idmode,idebug,&
         !
         UsePotentialDescription,UseNormalization,&
         UseSUPGFEM,             UseCoefficientCheck,&
         UseAnomalousTransportFT,UseAnomalousTransportGT,&
         !
         SolveElectron,SolveIons,&
         SolveBp,SolveBt,SolveEt,SolveEp,SolveEr,&
         SolveNn,SolveFr,SolveFb,SolveFt,SolveFp,&
         SolvePp,SolveQr,SolveQb,SolveQt,SolveQp,&
         !
         LockBpOnAxis,LockBtOnAxis,LockEtOnAxis,LockEpOnAxis,&
         LockErOnAxis,LockNnOnAxis,LockFrOnAxis,LockFbOnAxis,&
         LockFtOnAxis,LockFpOnAxis,LockPpOnAxis,LockQrOnAxis,&
         LockQbOnAxis,LockQtOnAxis,LockQpOnAxis,&
         !
         LockBpOnWall,LockBtOnWall,LockEtOnWall,LockEpOnWall,&
         LockErOnWall,LockNnOnWall,LockFrOnWall,LockFbOnWall,&
         LockFtOnWall,LockFpOnWall,LockPpOnWall,LockQrOnWall,&
         LockQbOnWall,LockQtOnWall,LockQpOnWall,&
         !
         TestMS,TestAV,TestAT,TestDT,TestGV,TestGT,&
         TestES,TestEV,TestET,TestSS,TestLEQ,TestLAX,TestLWL
         
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: NID
    INTEGER,INTENT(OUT) :: IST,IERR

    NAMELIST /T2/ &
         c10rname,i0fnum,&
         NNMAX,NQMAX,NDMAX,NSMAX,NPMIN,NLMAX,&
         !
         CoordinateSwitch,TestCase,EqSet,&
         !
         i1mlvl,i1rdn2,d1rec,d0rw,RR,RA,&
         i0nm,i0nn,i0tm,i0tn,d0qc,d0qs,d0bc,&
         d1nc,d1ns,d1nw,d1tc,d1ts,d1tw,Pa,Pz,&
         ! AF
         dt,time_init,eps_conv, &
         ntmax,ntstep,nt0dmax,nt0dstep,nt2dmax,nt2dstep,nconvmax, &
         idfile,idprint,idplot,idmode,idebug,&
         !
         UsePotentialDescription,UseNormalization,&
         UseSUPGFEM,             UseCoefficientCheck,&
         UseAnomalousTransportFT,UseAnomalousTransportGT,&
         !
         SolveElectron,SolveIons,&
         SolveBp,SolveBt,SolveEt,SolveEp,SolveEr,&
         SolveNn,SolveFr,SolveFb,SolveFt,SolveFp,&
         SolvePp,SolveQr,SolveQb,SolveQt,SolveQp,&
         !
         LockBpOnAxis,LockBtOnAxis,LockEtOnAxis,LockEpOnAxis,&
         LockErOnAxis,LockNnOnAxis,LockFrOnAxis,LockFbOnAxis,&
         LockFtOnAxis,LockFpOnAxis,LockPpOnAxis,LockQrOnAxis,&
         LockQbOnAxis,LockQtOnAxis,LockQpOnAxis,&
         !
         LockBpOnWall,LockBtOnWall,LockEtOnWall,LockEpOnWall,&
         LockErOnWall,LockNnOnWall,LockFrOnWall,LockFbOnWall,&
         LockFtOnWall,LockFpOnWall,LockPpOnWall,LockQrOnWall,&
         LockQbOnWall,LockQtOnWall,LockQpOnWall,&
         !
         TestMS,TestAV,TestAT,TestDT,TestGV,TestGT,&
         TestES,TestEV,TestET,TestSS,TestLEQ,TestLAX,TestLWL

    READ(NID,T2,IOSTAT=IST,ERR=9800,END=9900)

    IERR=0
    RETURN

9800 IERR=8
    RETURN
9900 IERR=9
    RETURN
  END SUBROUTINE T2_NLIN

!     ***** INPUT PARAMETER LIST *****

  SUBROUTINE T2_PLST

    IMPLICIT NONE
    WRITE(6,'(A)') '# NAMELIST T2 contains followings'
    WRITE(6,'(A)') 'c10rname,i0fnum'
    WRITE(6,'(A)') 'NNMAX,NQMAX,NDMAX,NSMAX,NPMIN,NLMAX'
    WRITE(6,'(A)') 'CoordinateSwitch'
         !
             !
    WRITE(6,'(A)') 'i1mlvl,i1rdn2,d1rec,d0rw,RR,RA'
    WRITE(6,'(A)') 'i0nm,i0nn,i0tm,i0tn,d0qc,d0qs,d0bc'
    WRITE(6,'(A)') 'd1nc,d1ns,d1nw,d1tc,d1ts,d1tw,Pa,Pz'
    ! AF
    WRITE(6,'(A)') 'dt,time_init,eps_conv'
    WRITE(6,'(A)') 'ntmax,ntstep,nt0dmax,nt0dstep,nt2dmax,'
    WRITE(6,'(A)') 'nt2dstep,nconvmaxidfile,idprint,idplot'
    WRITE(6,'(A)') 'idmode,idebug'
         !
    WRITE(6,'(A)') 'UsePotentialDescription,UseNormalization'
    WRITE(6,'(A)') 'UseSUPGFEM,UseCoefficientCheck'
    WRITE(6,'(A)') 'UseAnomalousTransportFT,UseAnomalousTransportGT'
    !
    WRITE(6,'(A)') 'SolveField,  SolveElectron,SolveIons'
    WRITE(6,'(A)') 'SolveDensity,SolveFlux,SolvePressure,SolveHeatFlux'
         !
    WRITE(6,'(A)') 'LockPoloidalMageticFieldOnAxis'
    WRITE(6,'(A)') 'LockToroidalMageticFieldOnAxis'
    WRITE(6,'(A)') 'LockRadialElectricFieldOnAxis'
    WRITE(6,'(A)') 'LockPoloidalElectricFieldOnAxis'
    WRITE(6,'(A)') 'LockToroidalElectricFieldOnAxis'
    WRITE(6,'(A)') 'LockDensityOnAxis'
    WRITE(6,'(A)') 'LockRadialFluxOnAxis'
    WRITE(6,'(A)') 'LockParallelFluxOnAxis'
    WRITE(6,'(A)') 'LockToroidalFluxOnAxis'
    WRITE(6,'(A)') 'LockPoroidalFluxOnAxis'
    WRITE(6,'(A)') 'LockPressureOnAxis'
    WRITE(6,'(A)') 'LockRadialHeatFluxOnAxis'
    WRITE(6,'(A)') 'LockParallelHeatFluxOnAxis'
    WRITE(6,'(A)') 'LockToroidalHeatFluxOnAxis'
    WRITE(6,'(A)') 'LockPoroidalHeatFluxOnAxis'
         !
    WRITE(6,'(A)') 'LockPoloidalMageticFieldOnWall'
    WRITE(6,'(A)') 'LockToroidalMageticFieldOnWall'
    WRITE(6,'(A)') 'LockRadialElectricFieldOnWall'
    WRITE(6,'(A)') 'LockPoloidalElectricFieldOnWall'
    WRITE(6,'(A)') 'LockToroidalElectricFieldOnWall'
    WRITE(6,'(A)') 'LockDensityOnWall'
    WRITE(6,'(A)') 'LockRadialFluxOnWall'
    WRITE(6,'(A)') 'LockParallelFluxOnWall'
    WRITE(6,'(A)') 'LockToroidalFluxOnWall'
    WRITE(6,'(A)') 'LockPoroidalFluxOnWall'
    WRITE(6,'(A)') 'LockPressureOnWall'
    WRITE(6,'(A)') 'LockRadialHeatFluxOnWall'
    WRITE(6,'(A)') 'LockParallelHeatFluxOnWall'
    WRITE(6,'(A)') 'LockToroidalHeatFluxOnWall'
    WRITE(6,'(A)') 'LockPoroidalHeatFluxOnWall'
    
    RETURN

  END SUBROUTINE T2_PLST

  !     ****** CHECK INPUT PARAMETER ******
  
  SUBROUTINE T2_CHECK(IERR)

    USE T2COMM,ONLY: NPMIN
    IMPLICIT NONE
    INTEGER:: IERR

    IERR=0

    IF(NPMIN <= 0) THEN
       WRITE(6,'(A,I12)') &
            'XX T2_check: INVALID NPMIN: ',NPMIN
       IERR=1
    ENDIF
    RETURN
  END SUBROUTINE T2_CHECK

!     ****** SHOW PARAMETERS ******
  SUBROUTINE T2_VIEW
    RETURN
  END SUBROUTINE T2_VIEW
END MODULE T2PARM
