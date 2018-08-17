
C     $Id: wmcoml.inc,v 1.2 2013/01/22 16:21:46 fukuyama Exp $
C
      USE plcomm
      USE pllocal
      USE commpi
      USE plprof,ONLY: pl_mag_old,pl_prof_old
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-F,H,O-Z)
C     
      include '../dpseki/dpcom1.inc'
      include '../wmseki/wmcom0.inc'
      include '../wmseki/wmcom1.inc'
