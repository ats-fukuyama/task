C     $Id$

      INCLUDE 'fpcomm.inc'
      CHARACTER*80 KNAMFP0

      MODEFR=0
      MODEFW=0

    1 WRITE(6,*) '## How many input fpdatas?'
      READ(5,*,ERR=1,END=9000) NFMAX

    2 WRITE(6,*) '## First input fpdata file name?'
      READ(5,'(A)',ERR=2,END=9000) KNAMFP0
      CALL DPLDFP(KNAMFP0,NR0,RMIN,RMAX1)

      DO NFI=2,NFMAX
    3    WRITE(6,'(A,I2,A)') '## #',NFI,' input fpdata file name?'
         READ(5,'(A)',ERR=3,END=9000) KNAMFP0
         CALL DPLDFP(KNAMFP0,NR0,RMIN1,RMAX1)
      ENDDO
      RMAX=RMAX1
      NRMAX=NR0

      MODEFW=0
    4 WRITE(6,*) '## Output fpdata file name?'
      READ(5,'(A)',ERR=4,END=9000) KNAMFP
      CALL FPSAVE
 9000 CONTINUE
      STOP
      END
C
C****** LOAD VELOCITY DISTRIBUTION DATA ******
C
      SUBROUTINE DPLDFP(KNAMFP0,NR0,RMIN1,RMAX1)
C
      INCLUDE 'fpcomm.inc'
      CHARACTER*80 KNAMFP0
C      
      CALL FROPEN(21,KNAMFP0,0,MODEFR,'FP',IERR)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX DPLDFP: FROPEN: IERR=',IERR
         RETURN
      ENDIF
      REWIND(21)

      READ(21) NRMAX,NPMAX,NTHMAX,NSAMAX
      READ(21) DELR,DELP,DELTH,RMIN1,RMAX1
      DO NSA=1,NSAMAX
         READ(21) NS_NSA(NSA)
         READ(21) AEFP(NSA),AMFP(NSA),RNFP0(NSA),RTFP0(NSA)
         READ(21) (((FNS(NTH,NP,NR+NR0,NSA),NTH=1,NTHMAX),
     &                NP=1,NPMAX),NR=1,NRMAX)
      ENDDO
      CLOSE(21)
C
      WRITE(6,*) '# DATA WAS SUCCESSFULLY LOADED FROM THE FILE.'
      write(6,'(A,2I5,1P2E12.4,2I5)') 
     &     'NRMIN/MAX,RMIN,RMAX,NPMAX,NTHMAX=',
     &     NR0+1,NR0+NRMAX,RMIN1,RMAX1,NPMAX,NTHMAX

      NR0=NR0+NRMAX
      RETURN
      END
C
C     ***********************************************************
C
C           SAVE VELOCITY DISTRIBUTION DATA
C
C     ***********************************************************
C
      SUBROUTINE FPSAVE
C
      INCLUDE 'fpcomm.inc'
C
      CALL FWOPEN(21,KNAMFP,0,MODEFW,'FP',IERR)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX FPSAVE: FWOPEN: IERR=',IERR
         RETURN
      ENDIF
C
      REWIND(21)
      WRITE(21) NRMAX,NPMAX,NTHMAX,NSAMAX
      WRITE(21) DELR,DELP,DELTH,RMIN,RMAX
      DO NSA=1,NSAMAX
         WRITE(21) NS_NSA(NSA)
         WRITE(21) AEFP(NSA),AMFP(NSA),RNFP0(NSA),RTFP0(NSA)
         WRITE(21) (((FNS(NTH,NP,NR,NSA),NTH=1,NTHMAX),
     &                NP=1,NPMAX),NR=1,NRMAX)
      ENDDO
      CLOSE(21)
C
      WRITE(6,*) '# DATA WAS SUCCESSFULLY SAVED TO THE FILE.'
C
  900 RETURN
      END
