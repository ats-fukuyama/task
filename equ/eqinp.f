      module eqinp_mod
      use eqinit_mod
      public
      contains
c
c=======================================================================
      subroutine eqinp
c=======================================================================
c     read namelist &equ
c=======================================================================
      use aaa_mod
      use par_mod
      use geo_mod
      use equ_mod
      use vac_mod
      use eqv_mod
      use cnt_mod
      implicit none
! local variables
      character aname*10
      integer   kstepk, istop, ist, ierr, n, nn, il
      data kstepk / 0 /
c=======================================================================
c:restart
      if(iread.gt.0)then
      call eqrest
      else
c=======================================================================
c:initial
      call eqinit
c-----
      if(kstepk.eq.0)then
      rewind ft05
      kstepk=1
      endif
c-----
      istop=0
c      read(ft05,equ,end=888,err=999)
      call eqnlin(ft05,ist,ierr)
         write(6,*) 'eqnlin done',ierr
      if(ierr.eq.9) goto 888
      if(ierr.eq.8) goto 999
      call eqflin(ft05,ierr)
         write(6,*) 'eqflin done',ierr
      istop=1
 888  if(istop.eq.0)then
      write(ft06,*)'>>&equ : end<< :IOSTAT=', ist
      rewind ft05
C      write(ft06,equ)
      call eqnlin(-ft06,ist,ierr)
      stop
      endif
 999  if(istop.eq.0)then
      write(ft06,*)'>>&equ : err<< :IOSTAT=', ist
      rewind ft05
C      write(ft06,equ)
      call eqnlin(-ft06,ist,ierr)
      stop
      endif
c=======================================================================
c     input data check
c=======================================================================
      call eqchek(ierr)
      if(ierr.ne.0) stop
      endif
c=======================================================================
c     print out
c=======================================================================
      call eqview
c-----
      return
      end subroutine eqinp
c
c=======================================================================
      subroutine eqrest
c=======================================================================
c     restart namelist &equrst
c=======================================================================
      use aaa_mod
      use par_mod
      use geo_mod
      use equ_mod
      use vac_mod
      use eqv_mod
      use cnt_mod
      implicit none
! local variables
      integer   i, k, n, nc, nfix, nnc, istop
c=======================================================================
c=====NAMELIST FOR RESTART
      namelist/equrst/title,
     >rvac,zvac,msfx,
     >ivac,cvac,ncoil,ivtab,ivgrp,isep,dsep,
     >ieqmax,eeqmax,iodmax,eodmax,iadmax,eadmax,bavmax,bavmin,nsumax,
     >ieqout
c=======================================================================
      istop=0
      read(ft05,equrst,end=888,err=999)
      istop=1
 888  if(istop.eq.0)then
      write(6,*)'-----  STANDARD DATA FOR &EQURST --'
      rewind ft05
      istop=1
      endif 
 999  if(istop.eq.0)then
      write(6,*)'-----  ERROR AT EQINP : &EQURST ------'
      rewind ft05
      endif
c-----
      if(nr.gt.irdm)nr=irdm
      if(nz.gt.izdm)nz=izdm
      if(nv.gt.ivdm)nv=ivdm
c-----
      if(nr*nz*nv.eq.0)then
         write(6,*) '========== DATA ERROR : NR*NZ*NV=0 ====='
         stop
      endif
      return
      end subroutine eqrest
c
      end module eqinp_mod
