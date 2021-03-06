!     $Id$
SUBROUTINE GSOPEN
  RETURN
END SUBROUTINE GSOPEN

SUBROUTINE GSCLOS
  RETURN
END SUBROUTINE GSCLOS

SUBROUTINE PAGES
  RETURN
END SUBROUTINE PAGES

SUBROUTINE PAGEE
  RETURN
END SUBROUTINE PAGEE

SUBROUTINE WFGPPC(NW,NWMAX,KWD)
  RETURN
END SUBROUTINE WFGPPC

SUBROUTINE WFGPFC(NW,NWMAX,KWD)
  RETURN
END SUBROUTINE WFGPFC

SUBROUTINE WFGPBC(NW,NWMAX,KWD)
  RETURN
END SUBROUTINE WFGPBC

SUBROUTINE WFGPFR(NW,NWMAX,KWD)
  RETURN
END SUBROUTINE WFGPFR

SUBROUTINE WFGPRM
  RETURN
END SUBROUTINE WFGPRM

SUBROUTINE WFGNAS(ID)
  RETURN
END SUBROUTINE WFGNAS

SUBROUTINE WFGDIV
  RETURN
END SUBROUTINE WFGDIV

SUBROUTINE WFPLTA
  RETURN
END SUBROUTINE WFPLTA

SUBROUTINE GUFLSH
  RETURN
END SUBROUTINE GUFLSH

!     ****** GET CPUTIME ******

SUBROUTINE GUTIME(T)
  
  CALL DVTIME(NT,NTICK)
  T=REAL(DBLE(NT)/DBLE(NTICK))
  RETURN
END SUBROUTINE GUTIME

!     ++++++ CONVERT FROM CHARACTER TO INTEGER (ASCII) ++++++

SUBROUTINE CHRASC(KTEXT,IASC,NBCHAR)
  
  CHARACTER KTEXT*256
  DIMENSION IASC(256)
  IF(NBCHAR.LE.0) RETURN
  DO I=1,NBCHAR
     IASC(I)=ICHAR(KTEXT(I:I))
  END DO
  RETURN
END DO
  
!     ++++++ CONVERT FROM INTEGER TO CHARACTER (ASCII) ++++++

SUBROUTINE ASCCHR(IASC,KTEXT,NBCHAR)
  
  CHARACTER KTEXT*256
  DIMENSION IASC(256)
  IF(NBCHAR.LE.0) RETURN
  DO I=1,NBCHAR
     KTEXT(I:I)=CHAR(IASC(I))
  END DO
  RETURN
END SUBROUTINE ASCCHR
