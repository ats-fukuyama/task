! fpfout.f90
!
! ***********************************************************
!
!                FILE OUTPUT of GRAPH DATA
!
! ***********************************************************

      MODULE FPFOUT

      USE fpcomm

      contains

! ***********************************************************
!
!                        TIME EVOLUTION
!
! ***********************************************************

      SUBROUTINE FPFOTT(STRING,FT)

      REAL(rkind),DIMENSION(NTG1M):: FT
      CHARACTER(LEN=*),INTENT(IN):: STRING

      CALL FPFOUX1(STRING,NTG1,1,FT,NTG1M)
      RETURN
      END SUBROUTINE FPFOTT

! ***********************************************************
!
!                        RADIAL PROFILE
!
! ***********************************************************

      SUBROUTINE FPFOTR(STRING,FR)

      REAL(rkind),DIMENSION(NRMAX+1,NTG2M):: FR
      CHARACTER(LEN=*),INTENT(IN):: STRING

      CALL FPFOUX2(STRING,NRMAX,NTG2,FR,NRMAX+1)

      RETURN
      END SUBROUTINE FPFOTR

! ***********************************************************
!
!                        Momentum Dependence
!
! ***********************************************************

      SUBROUTINE FPFOTP(STRING,FG)

      REAL(rkind),DIMENSION(NTHMAX+1,NPMAX+1,NRMAX+1):: FG
      REAL(rkind),dimension(NPMAX+1,NRMAX+1):: FL
      CHARACTER(LEN=*),INTENT(IN):: STRING

   1  WRITE(6,'(A,I4,A)') '# INPUT NTH (1..',NTHMAX,' or 0) :'
      READ(5,*,ERR=1,END=9000) NTH
      IF(NTH.LT.1) GOTO 9000
      IF(NTH.GT.NTHMAX) GOTO 1

      DO NR=1,NRMAX
         DO NP=1,NPMAX
            FL(NP,NR)=FG(NTH,NP,NR)
         ENDDO
      ENDDO

      CALL FPFOUX3(STRING,NPMAX,NRMAX,FL,NPMAX+1)
      GOTO 1

 9000 RETURN
      END SUBROUTINE FPFOTP

! ***********************************************************
!
!                Two-dimensional data (np,nth)
!
! ***********************************************************

      SUBROUTINE FPFOTC(STRING,FG,MODE)
!
      REAL(rkind),DIMENSION(NTHMAX+1,NPMAX+1,NRMAX+1):: FG
      REAL(rkind),DIMENSION(NPMAX+1,NTHMAX+1):: FL
      CHARACTER(LEN=*),INTENT(IN):: STRING
!
    1 CONTINUE
      IF(NRMAX.GT.1) THEN
         IF(MODE.EQ.0) THEN
            NRG=NRMAX+1
         ELSE
            NRG=NRMAX
         ENDIF
         WRITE(6,'(A,I4,A)') '# INPUT NR (1..',NRG,' or 0) :'
         READ(5,*,ERR=1,END=9000) NR
         IF(NR.LT.1) GOTO 9000
         IF(NR.GT.NRG) GOTO 1
      ELSE
         NR=1
      END IF

      IF(MOD(MODE,2).EQ.0) THEN
         NPG=NPMAX
      ELSE
         NPG=NPMAX+1
      ENDIF

      IF(MOD(MODE/2,2).EQ.0) THEN
         NTHG=NTHMAX
      ELSE
         NTHG=NTHMAX+1
      ENDIF

      DO NTH=1,NTHG
         DO NP=1,NPG
            FL(NP,NTH)=FG(NTH,NP,NR)
         ENDDO
      ENDDO

      CALL FPFOUX4(STRING,NPG,NTHG,FL,NPMAX+1)
      IF(NRMAX.GT.1) GOTO 1

 9000 RETURN
      END SUBROUTINE FPFOTC
!-------------------------------------
      SUBROUTINE FPFOTC_2(STRING,FG,MODE)
!
      REAL(rkind),DIMENSION(NTHMAX,NPMAX,NRMAX):: FG
      REAL(rkind),DIMENSION(NPMAX,NTHMAX):: FL
      CHARACTER(LEN=*),INTENT(IN):: STRING
!
    1 CONTINUE
      IF(NRMAX.GT.1) THEN
         IF(MODE.EQ.0) THEN
            NRG=NRMAX+1
         ELSE
            NRG=NRMAX
         ENDIF
         WRITE(6,'(A,I4,A)') '# INPUT NR (1..',NRG,' or 0) :'
         READ(5,*,ERR=1,END=9000) NR
         IF(NR.LT.1) GOTO 9000
         IF(NR.GT.NRG) GOTO 1
      ELSE
         NR=1
      END IF

      IF(MOD(MODE,2).EQ.0) THEN
         NPG=NPMAX
      ELSE
         NPG=NPMAX+1
      ENDIF

      IF(MOD(MODE/2,2).EQ.0) THEN
         NTHG=NTHMAX
      ELSE
         NTHG=NTHMAX+1
      ENDIF

      DO NTH=1,NTHG
         DO NP=1,NPG
            FL(NP,NTH)=FG(NTH,NP,NR)
         ENDDO
      ENDDO

      CALL FPFOUX4(STRING,NPG,NTHG,FL,NPMAX+1)
      IF(NRMAX.GT.1) GOTO 1

 9000 RETURN
      END SUBROUTINE FPFOTC_2

!-------------------------------------
      SUBROUTINE FPFOUX1(STRING,N1MAX,N2MAX,FL,N1M)

      REAL(rkind),DIMENSION(N1M):: FL
      CHARACTER(LEN=*),INTENT(IN):: STRING

      WRITE(22,'(A)') 'PTG'
      WRITE(22,'(2I10)') N1MAX,1
      WRITE(22,'(1P5E15.7)') (PTG(N1),N1=1,N1MAX)
      WRITE(22,'(A)') STRING
      WRITE(22,'(2I10)') N1MAX,1
      DO N2=1,N2MAX
         WRITE(22,'(1P5E15.7)') (FL(N1),N1=1,N1MAX)
      ENDDO
      RETURN
      END SUBROUTINE FPFOUX1
!----------------------------------------
      SUBROUTINE FPFOUX2(STRING,N1MAX,N2MAX,FL,N1M)

      REAL(rkind),DIMENSION(N1M,N2MAX):: FL
      CHARACTER(LEN=*),INTENT(IN):: STRING

      WRITE(22,'(A)') 'RTG'
      WRITE(22,'(2I10)') 1,N2MAX
      WRITE(22,'(1P5E15.7)') (RTG(N2),N2=1,N2MAX)
      WRITE(22,'(A)') STRING
      WRITE(22,'(2I10)') N1MAX,N2MAX
      DO N2=1,N2MAX
         WRITE(22,'(1P5E15.7)') (FL(N1,N2),N1=1,N1MAX)
      ENDDO
      RETURN
      END SUBROUTINE FPFOUX2
!----------------------------------------
      SUBROUTINE FPFOUX3(STRING,N1MAX,N2MAX,FL,N1M)

      REAL(rkind),DIMENSION(N1M,N2MAX):: FL
      CHARACTER(LEN=*),INTENT(IN):: STRING

      WRITE(22,'(A)') STRING
      WRITE(22,'(2I10)') N1MAX,N2MAX
      DO N2=1,N2MAX
         WRITE(22,'(1P5E15.7)') (FL(N1,N2),N1=1,N1MAX)
      ENDDO
      RETURN
      END SUBROUTINE FPFOUX3
!----------------------------------------
      SUBROUTINE FPFOUX4(STRING,N1MAX,N2MAX,FL,N1M)

      REAL(rkind),DIMENSION(N1M,N2MAX):: FL
      CHARACTER(LEN=*),INTENT(IN):: STRING

      WRITE(22,'(A)') STRING
      WRITE(22,'(2I10)') N1MAX,N2MAX
      DO N2=1,N2MAX
         WRITE(22,'(1P5E15.7)') (FL(N1,N2),N1=1,N1MAX)
      ENDDO
      RETURN
      END SUBROUTINE FPFOUX4

!-----------------------------------------------
      END MODULE FPFOUT
