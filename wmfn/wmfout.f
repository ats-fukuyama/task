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
      USE libfio
      INCLUDE 'wmcomm.inc'
C
      CALL FWOPEN(21,KNAMWM,0,MODEFW,'WM',IERR)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX WMSAVE: FWOPEN: IERR=',IERR
         RETURN
      ENDIF
C
C      WRITE(6,*) MDSIZ,NDSIZ,NRMAX
C      WRITE(6,'(4I5,1P2E12.4)') 
C     &           ((((NR,ND,MD,I,CEFLD(I,MD,ND,NR),I=1,3),MD=1,MDSIZ),
C     &                  ND=1,NDSIZ),NR=1,NRMAX)
C
      REWIND(21)
      WRITE(21) MDSIZ,NDSIZ,NRMAX,NSMAX
      WRITE(21) RA,RR,BB,CRF,NPH0,NTH0
      WRITE(21) ((((CEFLD(I,MD,ND,NR),I=1,3),MD=1,MDSIZ),ND=1,NDSIZ),
     &           NR=1,NRMAX)
      WRITE(21) ((((CBFLD(I,MD,ND,NR),I=1,3),MD=1,MDSIZ),ND=1,NDSIZ),
     &           NR=1,NRMAX)
      WRITE(21) ((((PABS(MD,ND,NR,NS),       MD=1,MDSIZ),ND=1,NDSIZ),
     &           NR=1,NRMAX),NS=1,NSMAX)
      CLOSE(21)
C
      WRITE(6,*) '# DATA WAS SUCCESSFULLY SAVED TO THE FILE.'
C
      RETURN
      END
