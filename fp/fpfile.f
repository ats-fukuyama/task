C     $Id$
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
C
C     ***********************************************************
C
C           LOAD VELOCITY DISTRIBUTION DATA
C
C     ***********************************************************
C
      SUBROUTINE FPLOAD
C
      INCLUDE 'fpcomm.inc'
C
      CALL FROPEN(21,KNAMFP,0,MODEFR,'FP',IERR)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX DPLDFP: FROPEN: IERR=',IERR
         RETURN
      ENDIF
C
      REWIND(21)
      READ(21) NRMAX,NPMAX,NTHMAX,NSAMAX
      READ(21) DELR,DELP,DELTH,RMIN,RMAX
      DO NSA=1,NSAMAX
         READ(21) NS_NSA(NSA)
         READ(21) AEFP(NSA),AMFP(NSA),RNFP0(NSA),RTFP0(NSA)
         READ(21) (((FNS(NTH,NP,NR,NSA),NTH=1,NTHMAX),
     &                NP=1,NPMAX),NR=1,NRMAX)
      ENDDO
      CLOSE(21)
C
      WRITE(6,*) '# DATA WAS SUCCESSFULLY LOADED FROM THE FILE.'
C
  900 RETURN
      END
