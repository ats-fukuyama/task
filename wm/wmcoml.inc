
C     $Id$
C
      USE plcomm
      USE pllocal
      USE plprof,ONLY: pl_mag_old,pl_prof_old
      IMPLICIT COMPLEX*16(C),REAL*8(A-B,D-F,H,O-Z)
C     
      include '../mpi/mpicom.inc'
      include '../dp/dpcom1.inc'
      include '../wm/wmcom0.inc'
      include '../wm/wmcom1.inc'
