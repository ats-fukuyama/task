C
      SUBROUTINE WMLOOP(IERR)
C
      INCLUDE 'wmcomm.inc'
      REAL*8 NPHS0
C
      NPHSMAX_SAVE=NPHSMAX
      DO NPHS=1,NPHSM
         DO 

      SELECT CASE(MDLWM_NFS) THEN
      CASE(0)
         NPHSMAX=1
         
         CALL WMEXEC(IERR)
      CASE(1)
         DO NPHS=1,NPHSMAX-1
            NPH0S(NPHS)=-NPHSMAX/2+(NPHS-1)
         
