C     $Id$
C     Fortran MDSplus access to International Profile Database.
C     Making connection to server
      subroutine IPDB_OPEN(kdev, kdcg)
      implicit none
      include 'mdslib.inc'
      integer ishot, status
      character*40 server
      character*40 tree
      character*80 kdev, kdcg
      integer ikdev, ikdcg
      common /IPDB1/ server
      common /IPDB2/ tree, ishot
C
      server = 'tokamak-profiledb.ukaea.org.uk'//CHAR(0)
      call ktrim(kdev,ikdev)
      tree   = kdev(1:ikdev)//char(0)
      call ktrim(kdcg,ikdcg)
      read(kdcg,*) ishot
C     Make connection
      if (server.ne.'localhost') then
         status = MdsConnect(server)
         if (status/2*2 .eq. status) then 
            write(6,*) 'Mdsconnected to server ',server
            write(6,*) 'Mdsconnect returns status= ',status
         else
            write(6,*) 'Mds SUCCESSFULLY connected to server ',server
         endif
      endif
      end
C
C     Closing connection to server
      subroutine IPDB_CLOSE
      implicit none
      include 'mdslib.inc'
      integer status
      character*40 server
      common /IPDB1/ server
C
      if (server.ne.'localhost') status = MdsDisconnect()
      if (mod(status,2) .eq. 0) then
         write(6,*) 'Mdsconnect returns status= ',status
      else
         write(6,*) 'Mds SUCCESSFULLY closed connection to server'
      endif
      end
C
C     for 1D profiles
C
      subroutine IPDB_MDS1(kdev, kdcg, kfid, ntm,
     &                     t, val, ntmax, ierr)
      implicit none
      include 'mdslib.inc'
      external descr
      integer status
      character cresult*25
C     
      integer nptsx,nptsy,ierror,ittyrd,ittywr
      integer lenstr
      parameter (lenstr=30)
      character*40 server
      character*40 tree
      character*40 signal
      character*80 errmsg
      integer ishot, irank
      integer ianswer, iretlen, L
      character*(lenstr) xlab, ylab
      integer MAXPTS,  idate, iudate, ish
      parameter (MAXPTS=81920)
      real x(MAXPTS), y(MAXPTS), pgasa,ip
C
      character*80 kdev, kdcg, kfid*10
      integer ikdev, ikdcg, ikfid, ntm, ntmax, ierr
      real*8 t(ntm), val(ntm)
      common /IPDB2/ tree, ishot
C
      ierror=0
      ittyrd=5
      ittywr=6
      iretlen=0
c$$$      server = 'tokamak-profiledb.ukaea.org.uk'//CHAR(0)
c$$$      call ktrim(kdev,ikdev)
c$$$      tree   = kdev(1:ikdev)//char(0)
c$$$      call ktrim(kdcg,ikdcg)
c$$$      read(kdcg,*) ishot
c$$$C     Make connection
c$$$      if (server.ne.'localhost') then
c$$$         status = MdsConnect(server)
c$$$         if (status/2*2 .eq. status) then 
c$$$            write(6,*) 'Mdsconnected to server ',server
c$$$            write(6,*) 'Mdsconnect returns status= ',status
c$$$         else
c$$$            write(6,*) 'Mds SUCCESSFULLY connected to server ',server
c$$$         endif
c$$$      endif
C
C     Open tree
      status = MdsOpen(tree, ishot)
      if (status/2*2 .eq. status) then 
         write(6,*) 'MdsOpen FAILS for (tree, ishot) :',tree,ishot
         write(6,*) 'MdsOpen returns status= ',status
C      else
C         write (6,*) 'MdsOpen SUCCESS for (tree, ishot): ',tree,ishot
      endif
C     
C     Initialization
      ierr=0
      do l=1,ntm
         val(l) = 0.d0
      enddo
C
C     Get integer, string, and real 0-D values
c$$$      ianswer = descr(IDTYPE_LONG,ish,0)
c$$$      status = MdsValue('.zerod:shot'//CHAR(0),ianswer, 0, iretlen)
c$$$      write (6,*) ' SHOT: ',ish, status
c$$$      ianswer = descr(IDTYPE_LONG,idate,0)
c$$$      status = MdsValue('.zerod:date'//CHAR(0),ianswer, 0, iretlen)
c$$$      write (6,*) ' DATE: ',idate, status
c$$$      ianswer = descr(IDTYPE_LONG,iudate,0)
c$$$      status = MdsValue('.zerod:update'//CHAR(0),ianswer, 0, iretlen)
c$$$      write (6,*) ' UPDATE: ',iudate, status
c$$$      ianswer = descr(IDTYPE_FLOAT,pgasa,0)
c$$$      status = MdsValue('.zerod:pgasa'//CHAR(0),ianswer, 0, iretlen)
c$$$      write (6,*) ' PGASA: ',pgasa, status
c$$$C     
c$$$      ianswer = descr(IDTYPE_CSTRING,cresult,0,LEN(cresult))
c$$$      status = MdsValue('.zerod:phase'//CHAR(0),ianswer, 0, iretlen)
c$$$      write (6,*) ' PHASE: ',cresult, status
C     
c$$$      ianswer = descr(IDTYPE_FLOAT,ip,0)
c$$$      status = MdsValue('.zerod:ip'//CHAR(0),ianswer, 0, iretlen)
c$$$      write (6,*) ' Plasma current: ',ip, status
C     
      call ktrim(kfid,ikfid)
      SIGNAL='.oned:'//kfid(1:ikfid)
C
C     Get the rank of a signal:
      ianswer = descr(IDTYPE_LONG,irank,0)
      status = MdsValue('rank('//signal//')'//CHAR(0),ianswer, 0, 
     &     iretlen)
      if (mod(status,2) .eq. 0) then
         ierr = 1
         goto 9000
      endif
C      write (6,*)
C      write (6,*) ' Signal=',signal
C      write (6,*) ' Rank of signal: ', irank
C     Get x units
      ianswer = descr(IDTYPE_CSTRING,xlab,0,LEN(xlab))
      status = MdsValue('units_of(dim_of('//signal//',0))'//CHAR(0),
     &     ianswer, 0, iretlen)
C      write (6,*) ' X units: ',xlab
C     Get number of x points
      ianswer = descr(IDTYPE_LONG,nptsx,0)
      status = MdsValue('size(dim_of('//signal//',0))'//CHAR(0),
     &     ianswer, 0, iretlen)
C      write (6,*) ' Number of X points: ', nptsx
C     Data access for X
      ianswer = descr(IDTYPE_FLOAT,X,MAXPTS,0)
      status = MdsValue('dim_of('//signal//',0)'//CHAR(0),ianswer, 
     &     0, iretlen)
C      write (6,*) ' X array length: ',iretlen
C      write (6,*) ' X array: ',(x(L),L=1,nptsx)
      write(6,'(2A10,A7,I2,A2,I4,A1)') 
     &     ' Signal = ',KFID, 'Rank = ',irank,' (',nptsx,')'
      ntmax=nptsx
      do L=1,ntmax
         t(L) = dble(x(L))
      enddo

      if(irank.eq.1) then

C     Get Y units (if a 1-D signal)
         ianswer = descr(IDTYPE_CSTRING,ylab,0,LEN(ylab))
         status = MdsValue('units_of('//signal//')'//CHAR(0),
     &        ianswer, 0, iretlen)
C         write (6,*) ' Y units: ',ylab
C     Get number of y points (if a 1-D signal)
         ianswer = descr(IDTYPE_LONG,nptsy,0)
         status = MdsValue('size('//signal//')'//CHAR(0),ianswer, 
     &        0, iretlen)
C         write (6,*) ' Number of Y points: ', nptsy
C     Data access for Y (if a 1-D signal)
         ianswer = descr(IDTYPE_FLOAT,Y,MAXPTS,0)
         status = MdsValue(signal//CHAR(0),ianswer, 0, iretlen)
C     write (6,*) ' Y array length: ',iretlen
C     write (6,*) ' Y array: ',(y(L),L=1,nptsy)
         if(ntmax.eq.nptsy) then
            do L=1,ntmax
               val(L) = dble(y(L))
            enddo
         else
            stop 'error occurs'
         endif

      else

         write(ittywr,'('' Not prepared for rank ='',I2)') irank
         stop

      endif

C 9000 if (server.ne.'localhost') status = MdsDisconnect()
 9000 end

C     ******

C     for 2D profiles
C
      subroutine IPDB_MDS2(kdev, kdcg, kfid, nrm, ntm,
     &                     r, t, val, nrmax, ntmax, ierr)
      implicit none
      include 'mdslib.inc'
      external descr
      integer status
      character cresult*25
C     
      integer nptsx,nptsy,nptsz,ierror,ittyrd,ittywr
      integer lenstr
      parameter (lenstr=30)
      character*40 server
      character*40 tree
      character*40 signal
      character*80 errmsg
      integer ishot, irank
      integer ianswer, iretlen, L,LL
      character*(lenstr) xlab, ylab, zlab
      integer MAXPTS,  idate, iudate, ish
      parameter (MAXPTS=81920)
      real x(MAXPTS), y(MAXPTS), z(MAXPTS), pgasa, ip
C
      character*80 kdev, kdcg, kfid*10
      integer ikdev, ikdcg, ikfid, nrm, ntm, ntmax, nrmax, ierr, md
      real*8 t(ntm), r(nrm), val(ntm,nrm)
      real*8 valctr(ntm), valedg(ntm)
      real*8 fctr, aitken2p
      common /IPDB2/ tree, ishot
C
      ierror=0
      ittyrd=5
      ittywr=6
      iretlen=0
c$$$      server = 'tokamak-profiledb.ukaea.org.uk'//CHAR(0)
c$$$      call ktrim(kdev,ikdev)
c$$$      tree   = kdev(1:ikdev)//char(0)
c$$$      call ktrim(kdcg,ikdcg)
c$$$      read(kdcg,*) ishot
c$$$C     Make connection
c$$$      if (server.ne.'localhost') then
c$$$         status = MdsConnect(server)
c$$$         if (status/2*2 .eq. status) then 
c$$$            write(6,*) 'Mdsconnected to server ',server
c$$$            write(6,*) 'Mdsconnect returns status= ',status
c$$$         else
c$$$            write(6,*) 'Mds SUCCESSFULLY connected to server ',server
c$$$         endif
c$$$      endif
C
C     Open tree
      status = MdsOpen(tree, ishot)
      if (status/2*2 .eq. status) then 
         write(6,*) 'MdsOpen FAILS for (tree, ishot) :',tree,ishot
         write(6,*) 'MdsOpen returns status= ',status
      else
C         write (6,*) 'MdsOpen SUCCESS for (tree, ishot): ',tree,ishot
      endif
C
C     Initialization
      ierr=0
      do L=1,nrm
         do LL=1,ntm
            val(LL,L) = 0.d0
         enddo
      enddo
C
C     Get integer, string, and real 0-D values
c$$$      ianswer = descr(IDTYPE_LONG,ish,0)
c$$$      status = MdsValue('.zerod:shot'//CHAR(0),ianswer, 0, iretlen)
c$$$      write (6,*) ' SHOT: ',ish, status
c$$$      ianswer = descr(IDTYPE_LONG,idate,0)
c$$$      status = MdsValue('.zerod:date'//CHAR(0),ianswer, 0, iretlen)
c$$$      write (6,*) ' DATE: ',idate, status
c$$$      ianswer = descr(IDTYPE_LONG,iudate,0)
c$$$      status = MdsValue('.zerod:update'//CHAR(0),ianswer, 0, iretlen)
c$$$      write (6,*) ' UPDATE: ',iudate, status
c$$$      ianswer = descr(IDTYPE_FLOAT,pgasa,0)
c$$$      status = MdsValue('.zerod:pgasa'//CHAR(0),ianswer, 0, iretlen)
c$$$      write (6,*) ' PGASA: ',pgasa, status
c$$$C     
c$$$      ianswer = descr(IDTYPE_CSTRING,cresult,0,LEN(cresult))
c$$$      status = MdsValue('.zerod:phase'//CHAR(0),ianswer, 0, iretlen)
c$$$      write (6,*) ' PHASE: ',cresult, status
C     
c$$$      ianswer = descr(IDTYPE_FLOAT,ip,0)
c$$$      status = MdsValue('.zerod:ip'//CHAR(0),ianswer, 0, iretlen)
c$$$      write (6,*) ' Plasma current: ',ip, status
C     
C     Choose 1 and 2-D signals
      call ktrim(kfid,ikfid)
      SIGNAL='.twod:'//kfid(1:ikfid)
C     Get the rank of a signal:
      ianswer = descr(IDTYPE_LONG,irank,0)
      status = MdsValue('rank('//signal//')'//CHAR(0),ianswer, 0, 
     &     iretlen)
      if (mod(status,2) .eq. 0) then
         ierr = 1
         goto 9000
      endif
C      write (6,*)
C      write (6,*) ' Signal=',signal
C      write (6,*) ' Rank of signal: ', irank
C     Get x units
      ianswer = descr(IDTYPE_CSTRING,xlab,0,LEN(xlab))
      status = MdsValue('units_of(dim_of('//signal//',0))'//CHAR(0),
     &     ianswer, 0, iretlen)
C      write (6,*) ' X units: ',xlab
C     Get number of x points
      ianswer = descr(IDTYPE_LONG,nptsx,0)
      status = MdsValue('size(dim_of('//signal//',0))'//CHAR(0),
     &     ianswer, 0, iretlen)
C      write (6,*) ' Number of X points: ', nptsx
C     Data access for X
      ianswer = descr(IDTYPE_FLOAT,X,MAXPTS,0)
      status = MdsValue('dim_of('//signal//',0)'//CHAR(0),ianswer, 
     &     0, iretlen)
C      write (6,*) ' X array length: ',iretlen
C      write (6,*) ' X array: ',(x(L),L=1,nptsx)
      nrmax=nptsx
      do L=1,nrmax
         r(L) = dble(x(L))
      enddo

      if(irank.eq.2) then

C     Get Y units (if a 2-D signal)
         ianswer = descr(IDTYPE_CSTRING,ylab,0,LEN(zlab))
         status = MdsValue('units_of(dim_of('//signal//',1))'
     &        //CHAR(0),ianswer, 0, iretlen)
C         write (6,*) ' Y units: ',ylab
C     Get number of Y points
         ianswer = descr(IDTYPE_LONG,nptsy,0)
         status = MdsValue('size(dim_of('//signal//',1))'//CHAR(0),
     &        ianswer, 0, iretlen)
C         write (6,*) ' Number of Y points: ', nptsy
C     Data access for Y (if a 2-D signal)
         ianswer = descr(IDTYPE_FLOAT,Y,MAXPTS,0)
         status = MdsValue('dim_of('//signal//',1)'//CHAR(0),
     &        ianswer, 0, iretlen)
C         write (6,*) ' Y array length: ',iretlen
C         write (6,*) ' Y array: ',(y(L),L=1,nptsy)
         ntmax=nptsy
         do L=1,ntmax
            t(L) = dble(y(L))
         enddo

C     Get Z units (if a 2-D signal)
         ianswer = descr(IDTYPE_CSTRING,zlab,0,LEN(xlab))
         status = MdsValue('units_of('//signal//')'//CHAR(0),
     &        ianswer, 0, iretlen)
C         write (6,*) ' Z units: ',zlab
C     Get number of z points (if a 2-D signal)
         ianswer = descr(IDTYPE_LONG,nptsz,0)
         status = MdsValue('size('//signal//')'//CHAR(0),
     &        ianswer, 0, iretlen)
C         write (6,*) ' Number of Z points: ', nptsz
C     Data access for Z (if a 2-D signal)
         ianswer = descr(IDTYPE_FLOAT,Z,MAXPTS,0)
         status = MdsValue(signal//CHAR(0),
     &        ianswer, 0, iretlen)
C         write (6,*) ' Z array length: ',iretlen
C         write (6,*) ' Z array: ',(z(L),L=1,nptsz)
         write(6,'(2A10,A7,I2,A2,I4,A1,I4,A1)') 
     &       ' Signal = ',KFID, 'Rank = ',irank,' (',nptsx,',',nptsy,')'
         if (nptsz.eq.nptsx*nptsy) then
            do L=1,ntmax
               do LL=1,nrmax
                  val(L,LL) = dble(z(nrmax*(L-1)+LL))
               enddo
            enddo
         else
            stop 'error occurs'
         endif
C     
C     *****
C     For interpolation, center and/or edge value is added
C     if necessary. The variables with center value having to
C     be zero such as BPOL, RMINOR and SURF were beforehand 
C     manupulated and thus in this section MD=4 will be selected.
C     *****
C     
         MD=0
         IF(r(1).NE.0.D0.AND.r(nrmax).NE.1.D0) THEN
            DO L=1,ntmax
               valctr(L)=FCTR(r(1),r(2),val(L,1),val(L,2))
               valedg(L)=AITKEN2P(1.D0,val(L,nrmax),val(L,nrmax-1),
     &                   val(L,nrmax-2),r(nrmax),r(nrmax-1),r(nrmax-2))
            ENDDO
            MD=1
         ELSEIF(r(1).NE.0.D0.AND.r(nrmax).EQ.1.D0) THEN
            DO L=1,ntmax
               valctr(L)=FCTR(r(1),r(2),val(L,1),val(L,2))
            ENDDO
            MD=2
         ELSEIF(r(1).EQ.0.D0.AND.r(nrmax).NE.1.D0) THEN
            DO L=1,ntmax
               valedg(L)=AITKEN2P(1.D0,val(L,nrmax),val(L,nrmax-1),
     &                   val(L,nrmax-2),r(nrmax),r(nrmax-1),r(nrmax-2))
            ENDDO
            MD=3
         ELSEIF(r(1).EQ.0.D0.AND.r(nrmax).EQ.1.D0) THEN
            MD=4
         ELSE
            STOP 'error occurs.'
         ENDIF
C     
         CALL DATA_ERROR_CORRECT(kdev,kdcg,kfid,r,val,
     &        ntmax,nrm,ntm)
C     
         IF(MD.EQ.1) THEN
            nrmax=nrmax+2
            DO LL=nrmax-1,2,-1
               r(LL)=r(LL-1)
            ENDDO
            r(1)=0.D0
            r(nrmax)=1.D0
            DO L=1,ntmax
               DO LL=nrmax-1,2,-1
                  val(L,LL)=val(L,LL-1)
               ENDDO
               val(L,1)=valctr(L)
               val(L,nrmax)=valedg(L)
            ENDDO
         ELSEIF(MD.EQ.2) THEN
            nrmax=nrmax+1
            DO LL=nrmax,2,-1
               r(LL)=r(LL-1)
            ENDDO
            r(1)=0.D0
            DO L=1,ntmax
               DO LL=nrmax,2,-1
                  val(L,LL)=val(L,LL-1)
               ENDDO
               val(L,1)=valctr(L)
            ENDDO
         ELSEIF(MD.EQ.3) THEN
            nrmax=nrmax+1
            DO L=1,ntmax
               val(L,nrmax)=valedg(L)
            ENDDO
         ENDIF

      else

         write(ittywr,'('' Not prepared for rank ='',I2)') irank
         stop

      endif

C 9000 if (server.ne.'localhost') status = MdsDisconnect()
 9000 end
C     ******

C     for 2D profiles
C
      subroutine IPDB_RAW2(kdev, kdcg, kfid, nrm, ntm,
     &                     r, t, val, nrmax, ntmax, ierr)
      implicit none
      include 'mdslib.inc'
      external descr
      integer status
      character cresult*25
C     
      integer nptsx,nptsy,nptsz,ierror,ittyrd,ittywr
      integer lenstr
      parameter (lenstr=30)
      character*40 server
      character*40 tree
      character*40 signal
      character*80 errmsg
      integer ishot, irank
      integer ianswer, iretlen, L,LL
      character*(lenstr) xlab, ylab, zlab
      integer MAXPTS,  idate, iudate, ish
      parameter (MAXPTS=81920)
      real x(MAXPTS), y(MAXPTS), z(MAXPTS), pgasa, ip
C
      character*80 kdev, kdcg, kfid*10
      integer ikdev, ikdcg, ikfid, nrm, ntm, ntmax, nrmax, ierr
      real*8 t(ntm), r(nrm), val(ntm,nrm)
      common /IPDB2/ tree, ishot
C
      ierror=0
      ittyrd=5
      ittywr=6
      iretlen=0
c$$$      server = 'tokamak-profiledb.ukaea.org.uk'//CHAR(0)
c$$$      call ktrim(kdev,ikdev)
c$$$      tree   = kdev(1:ikdev)//char(0)
c$$$      call ktrim(kdcg,ikdcg)
c$$$      read(kdcg,*) ishot
c$$$C     Make connection
c$$$      if (server.ne.'localhost') then
c$$$         status = MdsConnect(server)
c$$$         if (status/2*2 .eq. status) then 
c$$$            write(6,*) 'Mdsconnected to server ',server
c$$$            write(6,*) 'Mdsconnect returns status= ',status
c$$$         else
c$$$            write(6,*) 'Mds SUCCESSFULLY connected to server ',server
c$$$         endif
c$$$      endif
C
C     Open tree
      status = MdsOpen(tree, ishot)
      if (status/2*2 .eq. status) then 
         write(6,*) 'MdsOpen FAILS for (tree, ishot) :',tree,ishot
         write(6,*) 'MdsOpen returns status= ',status
      else
C         write (6,*) 'MdsOpen SUCCESS for (tree, ishot): ',tree,ishot
      endif
C
C     Initialization
      ierr=0
      do L=1,nrm
         do LL=1,ntm
            val(LL,L) = 0.d0
         enddo
      enddo
C
C     Get integer, string, and real 0-D values
c$$$      ianswer = descr(IDTYPE_LONG,ish,0)
c$$$      status = MdsValue('.zerod:shot'//CHAR(0),ianswer, 0, iretlen)
c$$$      write (6,*) ' SHOT: ',ish, status
c$$$      ianswer = descr(IDTYPE_LONG,idate,0)
c$$$      status = MdsValue('.zerod:date'//CHAR(0),ianswer, 0, iretlen)
c$$$      write (6,*) ' DATE: ',idate, status
c$$$      ianswer = descr(IDTYPE_LONG,iudate,0)
c$$$      status = MdsValue('.zerod:update'//CHAR(0),ianswer, 0, iretlen)
c$$$      write (6,*) ' UPDATE: ',iudate, status
c$$$      ianswer = descr(IDTYPE_FLOAT,pgasa,0)
c$$$      status = MdsValue('.zerod:pgasa'//CHAR(0),ianswer, 0, iretlen)
c$$$      write (6,*) ' PGASA: ',pgasa, status
c$$$C     
c$$$      ianswer = descr(IDTYPE_CSTRING,cresult,0,LEN(cresult))
c$$$      status = MdsValue('.zerod:phase'//CHAR(0),ianswer, 0, iretlen)
c$$$      write (6,*) ' PHASE: ',cresult, status
C     
c$$$      ianswer = descr(IDTYPE_FLOAT,ip,0)
c$$$      status = MdsValue('.zerod:ip'//CHAR(0),ianswer, 0, iretlen)
c$$$      write (6,*) ' Plasma current: ',ip, status
C     
C     Choose 1 and 2-D signals
      call ktrim(kfid,ikfid)
      SIGNAL='.twod:'//kfid(1:ikfid)
C     Get the rank of a signal:
      ianswer = descr(IDTYPE_LONG,irank,0)
      status = MdsValue('rank('//signal//')'//CHAR(0),ianswer, 0, 
     &     iretlen)
      if (mod(status,2) .eq. 0) then
         ierr = 1
         goto 9000
      endif
C      write (6,*)
C      write (6,*) ' Signal=',signal
C      write (6,*) ' Rank of signal: ', irank
C     Get x units
      ianswer = descr(IDTYPE_CSTRING,xlab,0,LEN(xlab))
      status = MdsValue('units_of(dim_of('//signal//',0))'//CHAR(0),
     &     ianswer, 0, iretlen)
C      write (6,*) ' X units: ',xlab
C     Get number of x points
      ianswer = descr(IDTYPE_LONG,nptsx,0)
      status = MdsValue('size(dim_of('//signal//',0))'//CHAR(0),
     &     ianswer, 0, iretlen)
C      write (6,*) ' Number of X points: ', nptsx
C     Data access for X
      ianswer = descr(IDTYPE_FLOAT,X,MAXPTS,0)
      status = MdsValue('dim_of('//signal//',0)'//CHAR(0),ianswer, 
     &     0, iretlen)
C      write (6,*) ' X array length: ',iretlen
C      write (6,*) ' X array: ',(x(L),L=1,nptsx)
      nrmax=nptsx
      do L=1,nrmax
         r(L) = dble(x(L))
      enddo

      if(irank.eq.2) then

C     Get Y units (if a 2-D signal)
         ianswer = descr(IDTYPE_CSTRING,ylab,0,LEN(zlab))
         status = MdsValue('units_of(dim_of('//signal//',1))'
     &        //CHAR(0),ianswer, 0, iretlen)
C         write (6,*) ' Y units: ',ylab
C     Get number of Y points
         ianswer = descr(IDTYPE_LONG,nptsy,0)
         status = MdsValue('size(dim_of('//signal//',1))'//CHAR(0),
     &        ianswer, 0, iretlen)
C         write (6,*) ' Number of Y points: ', nptsy
C     Data access for Y (if a 2-D signal)
         ianswer = descr(IDTYPE_FLOAT,Y,MAXPTS,0)
         status = MdsValue('dim_of('//signal//',1)'//CHAR(0),
     &        ianswer, 0, iretlen)
C         write (6,*) ' Y array length: ',iretlen
C         write (6,*) ' Y array: ',(y(L),L=1,nptsy)
         ntmax=nptsy
         do L=1,ntmax
            t(L) = dble(y(L))
         enddo

C     Get Z units (if a 2-D signal)
         ianswer = descr(IDTYPE_CSTRING,zlab,0,LEN(xlab))
         status = MdsValue('units_of('//signal//')'//CHAR(0),
     &        ianswer, 0, iretlen)
C         write (6,*) ' Z units: ',zlab
C     Get number of z points (if a 2-D signal)
         ianswer = descr(IDTYPE_LONG,nptsz,0)
         status = MdsValue('size('//signal//')'//CHAR(0),
     &        ianswer, 0, iretlen)
C         write (6,*) ' Number of Z points: ', nptsz
C     Data access for Z (if a 2-D signal)
         ianswer = descr(IDTYPE_FLOAT,Z,MAXPTS,0)
         status = MdsValue(signal//CHAR(0),
     &        ianswer, 0, iretlen)
C         write (6,*) ' Z array length: ',iretlen
C         write (6,*) ' Z array: ',(z(L),L=1,nptsz)
         write(6,'(2A10,A7,I2,A2,I4,A1,I4,A1)') 
     &       ' Signal = ',KFID, 'Rank = ',irank,' (',nptsx,',',nptsy,')'
         if (nptsz.eq.nptsx*nptsy) then
            do L=1,ntmax
               do LL=1,nrmax
                  val(L,LL) = dble(z(nrmax*(L-1)+LL))
               enddo
            enddo
         else
            stop 'error occurs'
         endif
C     
         IF(r(nrmax).GT.1.D0) THEN
            DO LL=1,nrmax
               IF(r(LL).GT.1.D0) THEN
                  nrmax=LL-1
                  GOTO 1000
               ENDIF
            ENDDO
         ENDIF
 1000    CONTINUE
         IF(KFID.EQ.'NEXP'.OR.KFID.EQ.'NEXPEB') THEN
            DO L=1,ntmax
               DO LL=1,nrm
                  val(L,LL)=val(L,LL)*1.D-20
               ENDDO
            ENDDO
         ELSE
            DO L=1,ntmax
               DO LL=1,nrm
                  val(L,LL)=val(L,LL)*1.D-3
               ENDDO
            ENDDO
         ENDIF

      else

         write(ittywr,'('' Not prepared for rank ='',I2)') irank
         stop

      endif

C 9000 if (server.ne.'localhost') status = MdsDisconnect()
 9000 end

C     ******

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
      mds_errstr = status
      return
      end
