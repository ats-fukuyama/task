C     $Id$
      SUBROUTINE GSOPEN
      RETURN
      END
C
      SUBROUTINE GSCLOS
      RETURN
      END
C
      SUBROUTINE PAGES
      RETURN
      END
C
      SUBROUTINE PAGEE
      RETURN
      END
C
      SUBROUTINE WFGPPC(NW,NWMAX,KWD)
      RETURN
      END
C
      SUBROUTINE WFGPFC(NW,NWMAX,KWD)
      RETURN
      END
C
      SUBROUTINE WFGPBC(NW,NWMAX,KWD)
      RETURN
      END
C
      SUBROUTINE WFGPFR(NW,NWMAX,KWD)
      RETURN
      END
C
      SUBROUTINE WFGPRM
      RETURN
      END
C
      SUBROUTINE WFGNAS(ID)
      RETURN
      END
C
      SUBROUTINE WFGDIV
      RETURN
      END
C
      SUBROUTINE WFPLTA
      RETURN
      END
C
      SUBROUTINE GUFLSH
      RETURN
      END
C
C     ****** GET CPUTIME ******
C
      SUBROUTINE GUTIME(T)
C
      CALL DVTIME(NT,NTICK)
      T=REAL(DBLE(NT)/DBLE(NTICK))
      RETURN
      END
C
C     ++++++ CONVERT FROM CHARACTER TO INTEGER (ASCII) ++++++
C
      SUBROUTINE CHRASC(KTEXT,IASC,NBCHAR)
C
      CHARACTER KTEXT*256
      DIMENSION IASC(256)
      IF(NBCHAR.LE.0) RETURN
      DO 1000 I=1,NBCHAR
         IASC(I)=ICHAR(KTEXT(I:I))
 1000 CONTINUE
      RETURN
      END
C
C     ++++++ CONVERT FROM INTEGER TO CHARACTER (ASCII) ++++++
C
      SUBROUTINE ASCCHR(IASC,KTEXT,NBCHAR)
C
      CHARACTER KTEXT*256
      DIMENSION IASC(256)
      IF(NBCHAR.LE.0) RETURN
      DO 1000 I=1,NBCHAR
         KTEXT(I:I)=CHAR(IASC(I))
 1000 CONTINUE
      RETURN
      END
