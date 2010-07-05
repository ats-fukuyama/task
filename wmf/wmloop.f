C
      SUBROUTINE WMLOOP(IERR)
C
      INCLUDE 'wmcomm.inc'
      REAL*8 NPHS0
C
      NPHSMAX_SAVE=NPHSMAX
      DO NPHS=1,NPHSM
         NPH0S_SAVE(NPHS)=NPH0S(NPHS)
         PFRACS_SAVE(NPHS)=PFRAC(NPHS)
      ENDDO

      SELECT CASE(MDLWM_NFS) THEN
      CASE(0)
         NPHSMAX=1
         NPH0S(1)=NPH0
         PFRACS(1)=1.D0
      CASE(2:3)
         DO NPHS=1,NPHSMAX
            NPH0S(NPHS)=-NPHSMAX/2+(NPHS-1)
            PFRACS(NPHS)=1.D0
         ENDDO
      END SELECT

      DO 
         
