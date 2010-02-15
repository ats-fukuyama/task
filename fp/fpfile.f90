!     $Id$
!
!     ***********************************************************
!
!           SAVE VELOCITY DISTRIBUTION DATA
!
!     ***********************************************************
!
      MODULE fpfile

      USE fpcomm

      contains
!----------------------------------------------------------------
      SUBROUTINE FP_SAVE

      CALL FWOPEN(21,KNAMFP,0,MODEFW,'FP',IERR)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX FPSAVE: FWOPEN: IERR=',IERR
         RETURN
      ENDIF

      REWIND(21)
      WRITE(21) NRMAX,NPMAX,NTHMAX,NSAMAX
      WRITE(21) DELR,DELTH,RMIN,RMAX
      DO NSA=1,NSAMAX
         NSBA=NSB_NSA(NSA)
         WRITE(21) NS_NSA(NSA)
         WRITE(21) DELP(NSBA)
         WRITE(21) AEFP(NSA),AMFP(NSA),RNFP0(NSA),RTFP0(NSA)
         WRITE(21) (((FNS(NTH,NP,NR,NSBA),NTH=1,NTHMAX), &
                      NP=1,NPMAX),NR=1,NRMAX)
      ENDDO
      CLOSE(21)

      WRITE(6,*) '# DATA WAS SUCCESSFULLY SAVED TO THE FILE.'

  900 RETURN
      END SUBROUTINE FP_SAVE
!
!     ***********************************************************
!
!           LOAD VELOCITY DISTRIBUTION DATA
!
!     ***********************************************************
!
      SUBROUTINE FP_LOAD

      CALL FROPEN(21,KNAMFP,0,MODEFR,'FP',IERR)
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX DPLDFP: FROPEN: IERR=',IERR
         RETURN
      ENDIF

      REWIND(21)
      READ(21) NRMAX,NPMAX,NTHMAX,NSAMAX
      READ(21) DELR,DELTH,RMIN,RMAX
      DO NSA=1,NSAMAX
         READ(21) NS_NSA(NSA)
         READ(21) DELP(NSBA)
         READ(21) AEFP(NSA),AMFP(NSA),RNFP0(NSA),RTFP0(NSA)
         READ(21) (((FNS(NTH,NP,NR,NSA),NTH=1,NTHMAX), &
                      NP=1,NPMAX),NR=1,NRMAX)
      ENDDO
      CLOSE(21)

      WRITE(6,*) '# DATA WAS SUCCESSFULLY LOADED FROM THE FILE.'

  900 RETURN
      END SUBROUTINE FP_LOAD
!------------------------------------------      
      END MODULE fpfile
