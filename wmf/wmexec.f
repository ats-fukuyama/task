C     $Id$
C
C     ****** CALCULATE ANTENNA EXCITATION ******
C
      SUBROUTINE WMEXEC(IERR)
C
      INCLUDE 'wmcomm.inc'
C
      IERR=0
      MODEEG=0
C
C
C      if(mdlwmf.eq.0) then
C
C         CALL WMSETG(IERR)
C         IF(IERR.NE.0) RETURN
C         CALL DPCHEK(NTHMAX,NRMAX+1,XRHO(1),XRHO(NRMAX+1),RR,IERR)
C         IF(IERR.NE.0) RETURN
C         CALL WMSETJ(IERR)
C         IF(IERR.NE.0) RETURN

C         CALL WMSOLV
C         CALL WMEFLD
C         CALL WMBFLD
C         CALL WMPABS
C         if(NRANK.eq.0) then
C            CALL WMPFLX
C            CALL WMPANT
C         endif
C      else
         CALL DPCHEK(NTHMAX,NRMAX+1,XRHO(1),XRHO(NRMAX+1),RR,IERR)
         IF(IERR.NE.0) RETURN
         nrmax=nrmax+1
         call wmfem_pre
         call wmfem_allocate
         call wmfem_main
         call wmfem_post
         nrmax=nrmax-1
C      endif
C
      RETURN
      END
C
C     ****** DISPLAY OUTPUT DATA ******
C
      SUBROUTINE WMPOUT
C
      INCLUDE 'wmcomm.inc'
C
      IF(NPRINT.LT.1) RETURN
C
      IF(NRANK.EQ.0) THEN
         WRITE(6,601) NPH0,CRADTT,PABSTT,PCURT
         WRITE(6,602) (PABST(NS),NS=1,NSMAX)
      ENDIF
C
      IF(NPRINT.LT.2) RETURN
C
      IF(NRANK.EQ.0) THEN
C         WRITE(6,*) '   MD   IMP                        PABSKT'
C         DO ND=NDMIN,NDMAX
C            NDX=ND-NDMIN+1
C            NN=NPH0+ND
C         DO MD=MDMIN,MDMAX
C            MDX=MD-MDMIN+1
C            MM=NTH0+MD
C            WRITE(6,603) NN,MM,CRADKT(MDX,NDX),
C     &                   (PABSKT(MDX,NDX,NS),NS=1,NSMAX)
C         ENDDO
C         ENDDO
      ENDIF
C
      IF(NPRINT.LT.4) RETURN
C
      DO NR=1,NRMAX+1
         IF(NRANK.EQ.0) 
     &        WRITE(6,'(A,I3,1P6E10.2)') 'NR,E1,E2,E3=',NR,
     &               CEFLD(1,1,1,NR),CEFLD(2,1,1,NR),CEFLD(3,1,1,NR)
      ENDDO
C
      RETURN
C
  601 FORMAT(' ',5X,3X,'NPH0=',I4/
     &       ' ',5X,3X,'RANT=',1PE12.4,10X,'LANT=',1PE12.4/
     &       ' ',5X,3X,'PABS=',1PE12.4,10X,'IDRV=',1PE12.4)
  602 FORMAT(' ',5X,3X,27X,'PABS=',1P3E12.4)
C  603 FORMAT(' ',I5,I5,3X,1P2E12.4,3X,1P3E12.4)
C  604 FORMAT(' ',8X,3F8.4,24X,2F8.4)
C  605 FORMAT(' ',9F8.4)
C
      END
