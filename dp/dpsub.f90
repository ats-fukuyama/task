! dpsub.f90

MODULE dpsub

  PRIVATE
  PUBLIC dkaijou,CFQ,CFQZ,CFQZ_Z,CFQZ_EXP

CONTAINS

!     ****** function kaijou ******

  FUNCTION dkaijou(n)

      USE bpsd_kinds,ONLY: rkind
      USE bpsd_constants,ONLY: CI,PI
      IMPLICIT NONE
      INTEGER,INTENT(IN):: n
      REAL(rkind):: dkaijou
      INTEGER,PARAMETER::  &
          ka(0:10)=(/1,1,2,6,24,120,720,5040,40320,362880,3628800/)
      REAL(rkind):: d
      INTEGER:: i

      IF(n.LT.0) THEN
         WRITE(6,*) 'XX DKAIJO: WRONG ARGUMENT: ',n
      ELSEIF(n.LE.10) THEN
         dkaijou=DBLE(ka(n))
      ELSE
         d=DBLE(ka(10))
         DO i=11,n
            d=d*DBLE(i)
         ENDDO
         dkaijou=d
      ENDIF
      RETURN
    END FUNCTION dkaijou

!     ****** CALCULATE F ******

  FUNCTION CFQ(Q,CZ,CNPR,RMU)

      USE bpsd_kinds,ONLY: rkind
      IMPLICIT NONE
      REAL(rkind),INTENT(IN):: Q,RMU
      COMPLEX(rkind),INTENT(IN):: CZ,CNPR
      COMPLEX(rkind):: CFQ
      COMPLEX(rkind):: CFQ0

      CFQ0=CFQZ(Q,RMU-CZ)
      CFQ=CFQ0 &
         +RMU*CNPR**2/2.D0*(   CFQZ(Q-1,RMU-CZ) &
                            -2*CFQ0 &
                            +  CFQZ(Q+1,RMU-CZ))
      RETURN
  END FUNCTION CFQ
!     
!     ****** CALCULATE SHKAROFSKY ******
!           *** WITH Z-FUNCTION ***

      FUNCTION CFQZ(Q,CZ)

      USE bpsd_kinds,ONLY: rkind
      USE bpsd_constants,ONLY: CI,PI
      USE libdsp,ONLY: DSPFN
      IMPLICIT NONE
      REAL(rkind),INTENT(IN):: Q
      COMPLEX(rkind),INTENT(IN):: CZ
      COMPLEX(rkind):: CFQZ,CZ2,CDZ2
      COMPLEX(rkind):: CSUM,CTERM,CGZ
      REAL(rkind):: TEMP
      INTEGER:: NUMAX,NU

      IF(ABS(CZ).GT.15.D0) THEN
         NUMAX=20
         CSUM=0.D0
         DO NU=0,NUMAX
            CTERM=-dgamma(Q+NU)/(-CZ)**(NU+1)
            CSUM=CSUM+CTERM
            IF(ABS(CTERM).LE.1.D-12) GOTO 100
         ENDDO
  100    CONTINUE
         TEMP=DBLE(SQRT(CZ))
         IF(ABS(TEMP).LT.1.D-12) THEN
            CSUM=CSUM-  CI*PI*(-CZ)**(Q-1)*EXP(CZ)
         ELSEIF(TEMP.LT.0.D0) THEN
            CSUM=CSUM-2*CI*PI*(-CZ)**(Q-1)*EXP(CZ)
         ENDIF
         CFQZ=CSUM/dgamma(Q)
      ELSE
         NUMAX=NINT(Q-3.D0/2.D0)
         CSUM=(0.D0,0.D0)
         DO NU=0,NUMAX 
            CTERM=(-CZ)**NU*dgamma(Q-1-NU)
            CSUM=CSUM+CTERM
         ENDDO
         CGZ=CI*SQRT(CZ)
         CALL DSPFN(CGZ,CZ2,CDZ2)
         CSUM=CSUM+SQRT(PI)*(-CZ)**NUMAX*(CI*SQRT(CZ)*CZ2)
         CFQZ=CSUM/dgamma(Q)
      ENDIF
      RETURN
  END FUNCTION CFQZ
!     
!     ****** CALCULATE SHKAROFSKY ******
!           *** WITH Z-FUNCTION ***

  FUNCTION CFQZ_Z(Q,CZ)

      USE bpsd_kinds,ONLY: rkind
      USE bpsd_constants,ONLY: CI,PI
      USE libdsp
      IMPLICIT NONE
      REAL(rkind),INTENT(IN):: Q
      COMPLEX(rkind),INTENT(IN):: CZ
      COMPLEX(rkind):: CFQZ_Z
      COMPLEX(rkind):: CSUM,CTERM,CGZ,CZ2,CDZ2
      INTEGER:: NUMAX,NU

      NUMAX=NINT(Q-3.D0/2.D0)
      CSUM=(0.D0,0.D0)
      DO NU=0,NUMAX 
         CTERM=(-CZ)**NU*dgamma(Q-1-NU)
         CSUM=CSUM+CTERM
      ENDDO
      CGZ=CI*SQRT(CZ)
      CALL DSPFN(CGZ,CZ2,CDZ2)
      CSUM=CSUM+SQRT(PI)*(-CZ)**NUMAX*(CI*SQRT(CZ)*CZ2)
      CFQZ_Z=CSUM/dgamma(Q)
      RETURN
  END FUNCTION CFQZ_Z
!     
!     ****** CALCULATE SHKAROFSKY ******
!     *** WITH ASYMPTOTIC EXPANSION ***

  FUNCTION CFQZ_EXP(Q,CZ)

      USE bpsd_kinds,ONLY: rkind
      USE bpsd_constants,ONLY: PI
      IMPLICIT NONE
      REAL(rkind),INTENT(IN):: Q
      COMPLEX(rkind),INTENT(IN):: CZ
      COMPLEX(rkind):: CFQZ_EXP
      COMPLEX(rkind):: CFQZ,CSUM
      REAL(rkind):: RREZ,RIMZ,SIG
      INTEGER:: NU

      CFQZ=(0.D0,0.D0)

         RREZ=DBLE(CZ)
         IF(RREZ.GE.-(11.D0+(Q-5.D0/2.D0)*1.35D0).AND.RREZ.LE.11.D0)THEN
            DO NU=0,50
               CSUM=(-CZ)**NU*dgamma(Q-1-NU)
               CFQZ=CFQZ+CSUM
            ENDDO
            CFQZ=(CFQZ-PI*(-CZ)**(Q-3.D0/2.D0)*SQRT(CZ)*EXP(CZ)) &
                 /dgamma(Q)
         ELSE
            RIMZ=DIMAG(CZ)
            IF(RIMZ.EQ.0.D0.AND.RREZ.LT.0.D0)THEN
               SIG=1.D0
            ELSE
               SIG=0.D0
            ENDIF
            DO NU=0,10
               CSUM=dgamma(Q+NU)*((-CZ)**(-1-NU))
               CFQZ=CFQZ+CSUM
            ENDDO
            CFQZ=(CFQZ+(CFQZ+dgamma(Q+(NU+1))*(-CZ)**(-1-NU-1)))/2.D0
            CFQZ_EXP &
                 =-(CFQZ-PI*SIG*(-CZ)**(Q-3.D0/2.D0)*SQRT(CZ)*EXP(CZ)) &
                 /dgamma(Q)
         ENDIF
      RETURN
  END FUNCTION CFQZ_EXP

!****************************************************
!                   gamma function                  *
!                in double precision                *
!      COPYRIGHT : M.Mori  JUNE 30 1989  V.1        *
!****************************************************

  FUNCTION dgamma(X)

      USE bpsd_kinds,ONLY: rkind
      IMPLICIT NONE
      REAL(rkind),INTENT(IN):: X
      REAL(rkind):: dgamma
      INTEGER,PARAMETER:: IN=19
      REAL(rkind),PARAMETER:: &
           C(0:19)=(/  1.0D0, &
                      -0.4227843350984671D0, &
                      -0.2330937364217867D0, &
                       0.1910911013876915D0, &
                      -0.2455249000540002D-1, &
                      -0.1764524455014432D-1, &
                       0.8023273022267347D-2, &
                      -0.8043297756042470D-3, &
                      -0.3608378162548D-3, &
                       0.1455961421399D-3, &
                      -0.175458597517D-4, &
                      -0.25889950224D-5, &
                       0.13385015466D-5, &
                      -0.2054743152D-6, &
                      -0.1595268D-9, &
                       0.62756218D-8, &
                      -0.12736143D-8, &
                       0.923397D-10, &
                       0.120028D-10, &
                      -0.42202D-11 /)
      REAL(rkind):: XX,A,FCTR,Z,Y
      INTEGER:: M,MG,I

      IF (X .GT. 57.0D0) GO TO 901

      XX = X
      IF (XX .LE. 1.5D0) THEN
        IF (XX .GE. 0.5D0) THEN
          A = XX - 1.0D0
          FCTR = 1.0D0
        ELSE
          M = INT(XX)
          A = XX - M
          IF (A .EQ. 0.0D0) THEN
            GO TO 902
          ELSE IF (A .GE. -0.5D0) THEN
            MG = IABS(M) + 1
          ELSE
            MG = IABS(M) + 2
            A = A + 1.0D0
          END IF
          Z = 1.0D0
          DO I = 1, MG
            Z = Z * XX
            XX = XX + 1.0D0
          ENDDO
          FCTR = 1.0D0 / Z
        END IF

      ELSE
        M = INT (XX)
        A = XX - M
        IF (A .LE. 0.5D0) THEN
          MG = M - 1
        ELSE
          MG = M
          A = A - 1.0D0
        END IF
        Z = 1.0D0
        DO I = 1, MG
          Z = Z * (XX - 1.0D0)
          XX = XX - 1.0D0
       ENDDO
        FCTR = Z
      END IF

      Y = C(IN)
      DO I = IN - 1, 0, -1
        Y = C(I) + A * Y
      ENDDO

      dgamma = FCTR / ((1.0D0 + A) * Y)
      RETURN

  901 CONTINUE
      WRITE (6,2001) X
 2001 FORMAT (' (FUNC.dgamma) X(=',D23.16,')', &
              ' must be smaller than 57.0')
      dgamma = 1.0D75
      RETURN

  902 CONTINUE
      WRITE (6,2002) X
 2002 FORMAT (' (FUNC.dgamma) invalid argument', &
              ' X =',D23.16)
      dgamma = 1.0D75
      RETURN

  END FUNCTION dgamma
END MODULE dpsub
