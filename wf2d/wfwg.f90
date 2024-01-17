! wfwg.f90

MODULE wfwg

  PRIVATE
  PUBLIC wf_set_ewg

CONTAINS

!     ***** SET BOUNDARY ELECTRIC FIELD *****

  SUBROUTINE wf_set_ewg

    USE wfcomm
    USE wfload,ONLY: wf_read_wg 
    implicit none
    INTEGER:: NSD,NN1,NN2,NN,NBSD,NBND,IERR
    REAL(rkind):: ANGLE,R,Z,PHASE,PROD,FACTOR,SN,th_wg,thn1,thn2
    COMPLEX(rkind):: CEX,CEY,CEZ,CEBSD0,CEBND0
    REAL(rkind),PARAMETER:: EPSWG=1.D-12

    ANGLE=ANGWG*PI/180.D0

    ! --- WG Electric field on boundary side ---  

    NBSID=0
    DO NSD=1,NSDMAX
       IF(KASID(NSD).EQ.1) NBSID=NBSID+1
    END DO
    IF(ALLOCATED(NSDBS)) DEALLOCATE(NSDBS)
    IF(ALLOCATED(CEBSD)) DEALLOCATE(CEBSD)
    ALLOCATE(NSDBS(NBSID))
    ALLOCATE(CEBSD(NBSID))
    NBSID=0
    DO NSD=1,NSDMAX
       IF(KASID(NSD).EQ.1) THEN
          NBSID=NBSID+1
          NSDBS(NBSID)=NSD
          KBSID(NSD)=NBSID
       ELSE
          KBSID(NSD)=0
       END IF
    END DO
    DO NBSD=1,NBSID
       NSD=NSDBS(NBSD)
       NN1=NDSID(1,NSD)
       NN2=NDSID(2,NSD)
       R=0.5D0*(RNODE(NN1)+RNODE(NN2))
       Z=0.5D0*(ZNODE(NN1)+ZNODE(NN2))

!       WRITE(6,'(A,I6,2ES12.4)') 'wg: modewg,R,Z:',modelwg,R,Z

       SELECT CASE(MODELWG)
       CASE(0,1)
          IF((R.GE.R1WG).AND.(R.LE.R2WG).AND. &
               (Z.GE.Z1WG).AND.(Z.LE.Z2WG)) THEN
             WRITE(6,'(A,2ES12.4)') 'R,Z in WG:',R,Z
             PROD=(R2WG-R1WG)*(RNODE(NN2)-RNODE(NN1)) &
                  +(Z2WG-Z1WG)*(ZNODE(NN2)-ZNODE(NN1))
             IF(ABS(R1WG-R2WG).LT.1.D-8) THEN
                FACTOR=(Z-0.5D0*(Z1WG+Z2WG))**2/(Z1WG-Z2WG)**2
             ELSE IF(ABS(Z1WG-Z2WG).LT.1.D-8) THEN
                FACTOR=(R-0.5D0*(R1WG+R2WG))**2/(R1WG-R2WG)**2
             ELSE
                FACTOR=(R-0.5D0*(R1WG+R2WG))**2/(R1WG-R2WG)**2 &
                      +(Z-0.5D0*(Z1WG+Z2WG))**2/(Z1WG-Z2WG)**2
             END IF
             SN=SQRT((R   -R1WG)**2+(Z   -Z1WG)**2) &
                  /SQRT((R2WG-R1WG)**2+(Z2WG-Z1WG)**2) ! SN=0 at 1, 1 at 2
             PHASE=(PH1WG+(PH2WG-PH1WG)*SN+DPHWG*4.D0*SN*(1.D0-SN))*PI/180.D0
             CEBSD0= AMPWG*EXP(CII*PHASE) &
                  *(COS(ANGLE)+CII*ELPWG*SIN(ANGLE))
             IF(MODELWG.EQ.1) CEBSD0=CEBSD0*EXP(-10.D0*FACTOR)
             ! IF(PROD.GT.0.D0) CEBSD0=-CEBSD0
             CEBSD(NBSD)=CEBSD0
             IF(nrank.EQ.0.AND.idebuga(41).EQ.1) &
                  WRITE(6,'(A,2I8,1P5E12.4)') &
                  'SD:',NSD,NBSD,CEBSD(NBSD),factor,PHASE,ANGLE
          ELSE
             CEBSD(NBSD)=(0.D0,0.D0)
          END IF
       CASE(5,6)
          th_wg=ATAN2(Z,R-RR)*180.D0/PI
          IF(th_wg.GE.th_wg_min.AND.th_wg.LE.th_wg_max) THEN
             thn1=ATAN2(znode(NN1),(rnode(NN1)-RR))
             thn2=ATAN2(znode(NN2),(rnode(NN2)-RR))
             WRITE(6,'(A,I6,3ES12.4)') &
                  'thn1,thn2,thn1-thn2:',NSD,thn1,thn2,thn1-thn2
             factor=(th_wg-th_wg_min)/(th_wg_max-th_wg_min)
             phase=phase_wg_min+(phase_wg_max-phase_wg_min)*factor &
                  +(phase_wg_cen-0.5D0*(phase_wg_min+phase_wg_max)) &
                  *4.D0*factor*(1.D0-factor)
             CEBSD0= AMPWG*EXP(CII*PHASE*PI/180.D0) &
                  *(COS(ANGLE)+CII*ELPWG*SIN(ANGLE))
             IF(MODELWG.EQ.6) THEN
                factor=4.D0*factor*(1.D0-factor)
                CEBSD0=CEBSD0*EXP(-gauss_wg*4.D0*factor*(1.D0-factor))
             END IF
             CEBSD(NBSD)=CEBSD0
             IF(nrank.EQ.0.AND.idebuga(41).EQ.1) THEN
!                WRITE(6,'(A,3ES12.4)') '*th_wg:',th_wg,th_wg_min,th_wg_max
                WRITE(6,'(A,2I8,1P5E12.4)') &
                     'SD:',NSD,NBSD,CEBSD(NBSD),factor,PHASE, &
                     EXP(-gauss_wg*factor)
             END IF
          ELSE
             CEBSD(NBSD)=(0.D0,0.D0)
          END IF
       CASE(12)
          IF((R.GE.R1WG-EPSWG).AND.(R.LE.R2WG+EPSWG).AND. &
               (Z.GE.Z1WG-EPSWG).AND.(Z.LE.Z2WG+EPSWG)) THEN
             PROD=(R2WG-R1WG)*(RNODE(NN2)-RNODE(NN1)) &
                  +(Z2WG-Z1WG)*(ZNODE(NN2)-ZNODE(NN1))
             CALL wf_read_wg(Z,CEX,CEY,CEZ,IERR)
             IF(nrank.EQ.0.AND.idebuga(41).EQ.1) &
                  write(6,'(A,1P6E12.4)') 'R,Z,CEY=', &
                  R,Z,CEY,PROD,ZNODE(NN2)-ZNODE(NN1)
!!!           IF(PROD.GT.0.D0) CEY=-CEY
             CEBSD(NBSD)=AMPWG*CEY
          ELSE
             CEBSD(NBSD)=(0.D0,0.D0)
          END IF
       END SELECT
    END DO

    ! --- WG Electric field on boundary node ---  

    NBNOD=0
    DO NN=1,NNMAX
       IF(KANOD(NN).EQ.1) NBNOD=NBNOD+1
    END DO
    IF(ALLOCATED(NNDBS)) DEALLOCATE(NNDBS)
    IF(ALLOCATED(CEBND)) DEALLOCATE(CEBND)
    ALLOCATE(NNDBS(NBNOD))
    ALLOCATE(CEBND(NBNOD))
    NBNOD=0
    DO NN=1,NNMAX
       IF(KANOD(NN).EQ.1) THEN
          NBNOD=NBNOD+1
          NNDBS(NBNOD)=NN
          KBNOD(NN)=NBNOD
       ELSE
          KBNOD(NN)=0
       END IF
    END DO
    DO NBND=1,NBNOD
       NN=NNDBS(NBND)
       R=RNODE(NN)
       Z=ZNODE(NN)
       SELECT CASE(MODELWG)
       CASE(0,1)
          IF((R.GE.R1WG).AND.(R.LE.R2WG).AND. &
               (Z.GE.Z1WG).AND.(Z.LE.Z2WG)) THEN
             IF(ABS(R1WG-R2WG).LT.1.D-8) THEN
                FACTOR=(Z-0.5D0*(Z1WG+Z2WG))**2/(Z1WG-Z2WG)**2
             ELSE IF(ABS(Z1WG-Z2WG).LT.1.D-8) THEN
                FACTOR=(R-0.5D0*(R1WG+R2WG))**2/(R1WG-R2WG)**2
             ELSE
                FACTOR=(R-0.5D0*(R1WG+R2WG))**2/(R1WG-R2WG)**2 &
                      +(Z-0.5D0*(Z1WG+Z2WG))**2/(Z1WG-Z2WG)**2
             END IF
             SN=SQRT((R   -R1WG)**2+(Z   -Z1WG)**2) &
                  /SQRT((R2WG-R1WG)**2+(Z2WG-Z1WG)**2) ! SN=0 at 1, 1 at 2
             PHASE=(PH1WG+(PH2WG-PH1WG)*SN+DPHWG*4.D0*SN*(1.D0-SN))*PI/180.D0
             CEBND0= AMPWG*EXP(CII*PHASE) &
                  *(SIN(ANGLE)-CII*ELPWG*COS(ANGLE))
             IF(MODELWG.EQ.1) CEBND0=CEBND0*EXP(-10.D0*FACTOR)
             CEBND(NBND)=CEBND0
             IF(nrank.EQ.0.AND.idebuga(41).EQ.1) &
                  WRITE(6,'(A,2I8,1P5E12.4)') &
                  'ND:',NN,NBND,CEBND(NBND),factor,PHASE,ANGLE
          ELSE
             CEBND(NBND)=(0.D0,0.D0)
          END IF
       CASE(5,6)
          th_wg=ATAN2(Z,R-RR)*180.D0/PI
          IF(th_wg.GE.th_wg_min.AND.th_wg.LE.th_wg_max) THEN
             factor=(th_wg-th_wg_min)/(th_wg_max-th_wg_min)
             phase=phase_wg_min+(phase_wg_max-phase_wg_min)*factor &
                  +(phase_wg_cen-0.5D0*(phase_wg_min+phase_wg_max)) &
                  *4.D0*factor*(1.D0-factor)
             CEBND0= AMPWG*EXP(CII*PHASE*PI/180.D0) &
                  *(SIN(ANGLE)-CII*ELPWG*COS(ANGLE))
             IF(MODELWG.EQ.6) THEN
                factor=4.D0*factor*(1.D0-factor)
                CEBND0=CEBND0*EXP(-gauss_wg*4.D0*factor*(1.D0-factor))
             END IF
             CEBND(NBND)=CEBND0
             IF(nrank.EQ.0.AND.idebuga(41).EQ.1) THEN
!                WRITE(6,'(A,3ES12.4)') '*th_wg:',th_wg,th_wg_min,th_wg_max
                WRITE(6,'(A,2I8,1P5E12.4)') &
                     'ND:',NN,NBND,CEBND(NBND),factor,PHASE, &
                     EXP(-gauss_wg*factor)
             END IF
          ELSE
             CEBND(NBND)=(0.D0,0.D0)
          END IF
       CASE(12)
          IF((R.GE.R1WG-EPSWG).AND.(R.LE.R2WG+EPSWG).AND. &
               (Z.GE.Z1WG-EPSWG).AND.(Z.LE.Z2WG+EPSWG)) THEN
             CALL wf_read_wg(Z,CEX,CEY,CEZ,IERR)
             IF(nrank.EQ.0.AND.idebuga(41).EQ.1) &
                  write(6,'(A,1P4E12.4)') 'R,Z,CEZ=',R,Z,CEZ
             CEBND(NBND)=AMPWG*CEZ
          ELSE
             CEBND(NBND)=(0.D0,0.D0)
          END IF
       END SELECT
    END DO
    RETURN
  END SUBROUTINE wf_set_ewg

END MODULE wfwg
