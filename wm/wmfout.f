C     $Id$
C
C***********************************************************************
C
C     Save field data
C
C***********************************************************************
C
      SUBROUTINE WMSAVE
C
      INCLUDE '../wm/wmcomm.inc'
C
      CALL FWOPEN(21,KNAMWM,0,5,'WM',IERR)
      IF(IERR.NE.0) RETURN
C
      REWIND(21)
      WRITE(21) MDMAX,NDMAX,NRMAX,NSMAX
      WRITE(21) ((((CEFLD(I,MD,ND,NR),I=1,3),MD=1,MDMAX),ND=1,NDMAX),
     &           NR=1,NRMAX)
      WRITE(21) ((((CBFLD(I,MD,ND,NR),I=1,3),MD=1,MDMAX),ND=1,NDMAX),
     &           NR=1,NRMAX)
      WRITE(21) ((((PABS(MD,ND,NR,NS),       MD=1,MDMAX),ND=1,NDMAX),
     &           NR=1,NRMAX),NS=1,NSMAX)
      CLOSE(21)
C
      WRITE(6,*) '# DATA WAS SUCCESSFULLY SAVED TO THE FILE.'
C
      RETURN
      END
