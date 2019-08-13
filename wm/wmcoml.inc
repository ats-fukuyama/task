
C     $Id$
C
      USE plcomm
      USE pllocal
      USE commpi
      USE plprof,ONLY: pl_mag_old,pl_prof_old
      USE dpcomm
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-F,H,O-Z)
C     
      include '../wm/wmcom0.inc'
      include '../wm/wmcom1.inc'
