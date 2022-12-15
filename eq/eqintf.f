C     $Id$
C
C     ***** GET PSIN AND MAGNETIC FIELD *****
C
C           PSIN=0 ON MAGNETIC AXIS
C           PSIN=1 ON PLASMA BOUNDARY
C
C           PSI=PSI0 ON MAGNETIC AXIS
C           PSI=0    ON PLASMA BOUNDARY
C
      SUBROUTINE GETRZ(RP,ZP,PHIP,BR,BZ,BT,RHON)
C
C     *** Input ***
C       RP, ZP : (R,Z) at which one would like to know magnetic fields
C       PHIP   : null
C     *** Output ***
C       BR     : major radius component of the magnetic field
C       BZ     : vertical component of the magnetic field
C       BT     : toroidal magnetic field
C       RHON   : normalized radial coordinate corresponding to (R,Z)
C
      USE libspl2d
      INCLUDE '../eq/eqcomq.inc'
C
      CALL SPL2DD(RP,ZP,PSI,DPSIR,DPSIZ,
     &            RG,ZG,UPSIRZ,NRGM,NRGMAX,NZGMAX,IERR)
C
      PSIN=1.D0-PSI/PSI0
      IF(PSIN.LE.0.D0) THEN
         PSIN=0.D0
         BT=FNTTS(0.D0)/(2.D0*PI*RR)
         BR=0.D0
         BZ=0.D0
      ELSE
         BT=FNTTS(SQRT(PSIN))/(2.D0*PI*RP)
         BR=-DPSIZ/(2.D0*PI*RP)
         BZ= DPSIR/(2.D0*PI*RP)
C         write(6,'(A,1P3E12.4)') 'PSI,DPSIZ,DPSIR      =',PSI,DPSIR,DPSIZ
      ENDIF
      RHON=FNRHON(PSIN)
C
      RETURN
      END
C
C     ***** GET PARAMETERS *****
C
      SUBROUTINE EQGETB(BB1,RR1,RIP1,RA1,RKAP1,RDEL1,RB1)
C
      INCLUDE '../eq/eqcomq.inc'
C
      BB1  =BB
      RR1  =RR
      RIP1 =RIP
      RA1  =RA
      RKAP1=RKAP
      RDEL1=RDLT
      RB1  =RB
      RETURN
      END
C
C     ***** GET PRESSURE AS A FUNCTION OF PSIN *****
C
      SUBROUTINE GETPP(RHON,PP)
C
      INCLUDE '../eq/eqcomq.inc'
C
      PP=FNPPS(RHON)
      RETURN
      END
C
C     ***** GET SAFETY FACTOR AS A FUNCTION OF PSIN *****
C
      SUBROUTINE GETQP(RHON,QP)
C
      INCLUDE '../eq/eqcomq.inc'
C
      QP=FNQPS(RHON)
      RETURN
      END
C
C     ***** GET MINIMUM R AS A FUNCTION OF PSIN *****
C
      SUBROUTINE GETRMN(RHON,RRMINL)
C
      INCLUDE '../eq/eqcomq.inc'
C
      RRMINL=FNRRMN(RHON)
      RETURN
      END
C
C     ***** GET MAXIMUM R AS A FUNCTION OF PSIN *****
C
      SUBROUTINE GETRMX(RHON,RRMAXL)
C
      INCLUDE '../eq/eqcomq.inc'
C
      RRMAXL=FNRRMX(RHON)
      RETURN
      END
C
C     ***** GET MAGNETIC AXIS *****
C
      SUBROUTINE GETAXS(RAXIS1,ZAXIS1)
C
      INCLUDE '../eq/eqcomq.inc'
C
      RAXIS1=RAXIS
      ZAXIS1=ZAXIS
      RETURN
      END
C
C     ***** GET PLASMA BOUNDARY POSITION *****
C
      SUBROUTINE eqget_rzsu(rsu1,zsu1,nsumax1)
C
      INCLUDE '../eq/eqcomq.inc'
      REAL(rkind),ALLOCATABLE,INTENT(OUT):: rsu1(:),zsu1(:)
      INTEGER,INTENT(OUT):: nsumax1
      INTEGER:: nsu
C
      WRITE(6,'(A)') '@@@ point 31'
      WRITE(6,*) ALLOCATED(rsu1)
      WRITE(6,'(A)') '@@@ point 310'
      IF(ALLOCATED(rsu1)) DEALLOCATE(rsu1)
      WRITE(6,'(A)') '@@@ point 311'
      IF(ALLOCATED(zsu1)) DEALLOCATE(zsu1)
      WRITE(6,'(A)') '@@@ point 312'
      ALLOCATE(rsu1(nsumax),zsu1(nsumax))
      WRITE(6,'(A)') '@@@ point 32'
      DO nsu=1,nsumax
         rsu1(nsu)=rsu(nsu)
         zsu1(nsu)=zsu(nsu)
      enddo
      WRITE(6,'(A)') '@@@ point 33'
      nsumax1=nsumax
      return
      end
C
C     ***** GET PLASMA BOUNDARY POSITION *****
C
      SUBROUTINE GETRSU(RSU1,ZSU1,N,NSUMAX1)
C
      INCLUDE '../eq/eqcomq.inc'
      DIMENSION RSU1(N),ZSU1(N)
C
      NSUMAX1=NSUMAX
      DO NSU=1,MIN(N,NSUMAX)
         RSU1(NSU)=RSU(NSU)
         ZSU1(NSU)=ZSU(NSU)
      ENDDO
      RETURN
      END
C
C     ***** GET R and Z for rhot and th *****
C
      SUBROUTINE GET_RZ(rhon_,chip_,R_,Z_)
C
      USE libspl2d
      INCLUDE '../eq/eqcomq.inc'
      REAL(rkind):: chip_
C
      CALL SPL2DF(chip_,rhon_,R_,
     &                  CHIP,RHOT,URPS,NTHMP,NTHMAX+1,NRMAX,IERR)
      CALL SPL2DF(chip_,rhon_,Z_,
     &                  CHIP,RHOT,UZPS,NTHMP,NTHMAX+1,NRMAX,IERR)
      RETURN
      END
C
C     ***** GET magnetic field for rhot and th *****
C
      SUBROUTINE GET_RZB(rhon_,chip_,R_,Z_,BR_,BZ_,BT_,BB_)
C
      USE bpsd,ONLY: rkind
      IMPLICIT NONE
      REAL(rkind),INTENT(IN):: rhon_,chip_
      REAL(rkind),INTENT(OUT):: R_,Z_,BR_,BZ_,BT_,BB_
      REAL(rkind):: rhon_dummy
C
      CALL GET_RZ(rhon_,chip_,R_,Z_)
      CALL GETRZ(R_,Z_,0.D0,BR_,BZ_,BT_,rhon_dummy)
      BB_=SQRT(BR_**2+BZ_**2+BT_**2)
      RETURN
      END
C
C     ***** GET total magnetic field for rhot and th *****
C
      SUBROUTINE GET_B(rhon_,chip_,BB_)
C
      USE bpsd,ONLY: rkind
      IMPLICIT NONE
      REAL(rkind),INTENT(IN):: rhon_,chip_
      REAL(rkind),INTENT(OUT):: BB_
      REAL(rkind):: R_,Z_,BR_,BZ_,BT_,RHON
C
      CALL GET_RZ(rhon_,chip_,R_,Z_)
      CALL GETRZ(R_,Z_,0.D0,BR_,BZ_,BT_,RHON)
      BB_=SQRT(BR_**2+BZ_**2+BT_**2)
      RETURN
      END
C
C     ***** GET BBMIN and BBMAX for rhot *****
C
      SUBROUTINE GET_BMINMAX(rhon_,BBMIN_,BBMAX_)
C
      USE bpsd,ONLY: rkind
      USE libspl1d
      INCLUDE '../eq/eqcomq.inc'
      REAL(rkind),INTENT(IN):: rhon_
      REAL(rkind),INTENT(OUT):: BBMIN_,BBMAX_
C
      CALL SPL1DF(rhon_,BBMIN_,RHOT,UBBMIN,NRMAX,IERR)
      CALL SPL1DF(rhon_,BBMAX_,RHOT,UBBMAX,NRMAX,IERR)
      RETURN
      END
C
C     ***** GET DVDRHO for rhot *****
C
      SUBROUTINE GET_DVDRHO(rhon_,DVDRHO_)
C
      USE bpsd,ONLY: rkind
      USE libspl1d
      INCLUDE '../eq/eqcomq.inc'
      REAL(rkind),INTENT(IN):: rhon_
      REAL(rkind),INTENT(OUT):: DVDRHO_

      CALL SPL1DF(FNPSIP(rhon_),DVDRHO_,PSIP,UDVDRHO,NRMAX,IERR)

      RETURN
      END
      
