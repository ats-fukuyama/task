C     $Id$
C     Fortran MDSplus access to International Profile Database.
C     Making connection to server
      subroutine IPDB_OPEN(kdev, kdcg)
      implicit none
      character*80 kdev, kdcg
      return
      end
C
C
C     Closing connection to server
      subroutine IPDB_CLOSE
      implicit none
      return
      end
C
C     for 1D profiles
C
      subroutine IPDB_MDS1(kdev, kdcg, kfid, ntm,
     &                     t, val, ntmax, ierr)
      implicit none
C
      character*80 kdev, kdcg, kfid*10
      integer ntm, ntmax, ierr
      real*8 t(ntm), val(ntm)
C
      ierr=1
      return
      end
C
C     for 2D profiles
C
      subroutine IPDB_MDS2(kdev, kdcg, kfid, nrm, ntm,
     &                     r, t, val, nrmax, ntmax, ierr)
      implicit none
      character*80 kdev, kdcg, kfid*10
      integer nrm, ntm, ntmax, nrmax, ierr
      real*8 t(ntm), r(nrm), val(ntm,nrm)
C
      ierr=1
      return
      end
C
C     for 2D profiles
C
      subroutine IPDB_RAW2(kdev, kdcg, kfid, nrm, ntm,
     &                     r, t, val, nrmax, ntmax, ierr)
      implicit none
      character*80 kdev, kdcg, kfid*10
      integer nrm, ntm, ntmax, nrmax, ierr
      real*8 t(ntm), r(nrm), val(ntm,nrm)
C
      ierr=1
      return
      end
C
      integer function mds_errstr(ierr, str)
      implicit none
      integer ierr
      character*(*) str
      integer dsc, size, status
      character*19 cmd
      integer mdsfdescr, mdsfvalue, mds__cstring
      write(cmd, '(''GETMSG('',I11,'')'')') ierr
c     dsc = mdsfdescr(mds__cstring(), str, 0, len(str))
c     status = mdsfvalue(cmd //char(0), dsc, 0, size)
c     call mds_addblanks(str)
      mds_errstr = 0
      return
      end
