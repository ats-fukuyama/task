!@fortran double
MODULE mcnbi_all
  CONTAINS

!======================================================================c
!                                                                      c
!                                                                      c
!                                                                      c
!      PROGRAM    prof_test
      SUBROUTINE mcnbi
!                                                                      c
!----------------------------------------------------------------------c
!                                                                      c
!                                           by  s.murakami
! 04.05.27
!    The B<0 case with a B>0 file is added.
! 04.05.20
!    The pich up routine of out9 file is added.
! 07.12
!    To Fortran90, SUBROUTINE
!                                                                      c
!======================================================================c

  USE mcnmod, ONLY :a, aveng, avengs, ckb, cpi, depepa, depera, depeta, depppa,    &
    depths, depehs, depihs, deppra, deppta, deppxa, depvha, depvra, depvva, dtime, &
    e0int, eng, etlos, fl, fmx0, ich, iend, iit, iout, isw, it, it0, itmax, mm2,   &
    mpmax3, mrmax, mrmax3, mtmax3, mvmax, mvmax2, mxmax, mymax3, nloss, npt, ntest,&
    prf, r, rtptcl, tmass, wmxt, yc
  USE mcnmod, ONLY : allocate_b, deallocate_b, allocate_restrt2, deallocate_restrt2, &
    deallocate_restrt22, allocate_restrt1, allocate_through, deallocate_through
  IMPLICIT NONE
  INTEGER(4)::i, iir, iirm, iirp, im, iph, ir, ith, ixx, iyy, jl, n      ! , ivh, ivs
  REAL(8):: csr, cvol, delr, dpeirm, dpeirp, dpeirt, dprirm, dprirp, dprirt, dvhirm,&
    dvhirp, dvhirt, dvrirm, dvrirp, dvrirt, dvvmax, dw, etmax, feng, rra, rrr, srcrt,&
    vvmax          !  ,vpp
  REAL(8), ALLOCATABLE::  depep(:), deper(:), depet(:), deppp(:), depprf(:), deppt(:),&
    depvh(:), depvr(:), depex(:,:), deppx(:,:)

!----------------------------------------------------------------------c
!----------------------------------------------------------------------c
!     open(10,
!    & file='/ktmp1/mura/data2/VRHFXTOR.AXM150.V66D43P0.BM100N12.FMT'
!    & ,status='old',form='formatted')
!  out20 for radial distribution
!     open(20,file='test.out20',STATUS='UNKNOWN')
!  out21 for the Green function
!     open(21,file='test.out21',STATUS='UNKNOWN')
!  out22 for loss distribution
!     open(22,file='test.out22',STATUS='UNKNOWN')
!     open(10,file='ptcl.dat',status='unknown')


      CALL allocate_through
      ALLOCATE(depep(mpmax3), deper(mrmax3), depet(mtmax3), deppp(mpmax3))
      ALLOCATE(depprf(mrmax3),deppt(mtmax3), depvh(mrmax3), depvr(mrmax3))
      ALLOCATE(depex(mxmax,mymax3), deppx(mxmax,mymax3))

!---------------------------------------c
!            set parameters
!---------------------------------------c

      CALL setprm

!---------------------------------------c
!            set field data
!---------------------------------------c

      CALL setfld

!---------------------------------------c
!       set particle distribution
!---------------------------------------c

      CALL setptl

      CALL allocate_b

!---------------------------------------c
!       set Maxwell distribution
!---------------------------------------c

      CALL fmxwll

!---------------------------------------c
!      set charge exchange data
!---------------------------------------c
      CALL setcx

!=======================================c
!=======================================c
!=====      start main routine     =====c
!=======================================c
!=======================================c

!...for the case without orbit calculation
!      call mvptl

      DO 10 it = it0, itmax
!
!---------------------------------------c
!           move particle
!---------------------------------------c

         CALL mvptl

!...for the case without orbit calculation
!           time = time+dt
!
!$$$         ipsum = 0
!$$$         do 11 ip=1, ntest
!$$$            if (isw(ip).eq.1) ipsum = ipsum +1
!$$$   11    continue
!$$$c
!$$$         if (ipsum.eq.0) then
!$$$            write(6,*) '   ***** ipsum = 0 '
!$$$            go to 555
!$$$         end if

!---------------------------------------c
!         particle collsion
!---------------------------------------c

         CALL colptl

!---------------------------------------c
!      store the memories
!               for plotter routine
!---------------------------------------c

         IF (ich.eq.1) GO TO 10


!---------------------------------------c
!             diagnosis
!---------------------------------------c

         IF (it/iout*iout.eq.it) THEN
            CALL diagi
         ENDIF
!
!---------------------------------------c

 10   CONTINUE
!
  555 CONTINUE

!=======================================c
!=======================================c
!======      end main routine      =====c
!=======================================c
!=======================================c

      CALL deallocate_b

      CALL restrt(0)

      CALL deallocate_restrt2

      IF(ich.eq.1) GO TO 850

!---------------------------------------c
!---------------------------------------c

      DO n=1,ntest
        IF (isw(n).eq.0) WRITE(6,1002) n,eng(n),yc(n,1)
        IF (isw(n).lt.0) WRITE(6,1003) n,eng(n),yc(n,1),isw(n)
      END DO

      DO i=1,iend+1
        fl(i)=1.d0
      END DO

      dw=e0int*1.5d0/dfloat(iend)
!      dw=tiev*0.1d0

      DO n=1,ntest
        IF (isw(n).eq.0) THEN
          IF (jl.gt.iend) THEN
            jl=iend+1
          ELSE
            jl=eng(n)/dw+1
          ENDIF
          fl(jl)=fl(jl)+1
        ENDIF
      END DO

      PRINT *, ' '
      PRINT *, ' '
      PRINT *, ' '
      WRITE(6,665)
      WRITE(6,667)
      DO i=1,iend+1
        fl(i)=dlog(fl(i))
        feng  = (dfloat(i-1)+0.5d0)*dw
        WRITE(6,676) feng, fl(i)
      END DO

 665  FORMAT(/,5x,' lost particle distribution')
 667  FORMAT(7x,'eng[keV]',12x,'alog(fl)' )
 676  FORMAT(2e12.4)
      PRINT *, ' '
      PRINT *, ' '
      PRINT *, ' '

!      write(6,*) 'mplot,maxnpl',mplot,maxnpl

850   CONTINUE

!      write(12) rgc,zgc,pgc,tang,pang

!ppp
      PRINT *, ' ckb,tmass=', ckb,tmass

!     etmax : maxmum value of energy in grid

      etmax= e0int*2.d0
      vvmax= sqrt(2.0d0*etmax*ckb/tmass)
      dvvmax=vvmax/dfloat(mvmax-1)
!ppp
      PRINT *, ' Output of distribution function'
      PRINT *, ' vvmax=', vvmax

      DO im=1,mm2
!CDIR NODEP
        DO ir=1,mrmax3
!ppp
!!            irm = mod(ir,mrmax)                      ! 200710
!!            ssr=(dfloat(irm)-0.5d0)/dfloat(mrmax-1)  ! 200710

!!            DO ivs=1,mvmax                           ! 200712
!!               vpp=(dfloat(ivs)-0.5d0)*dvvmax
!!               DO ivh=1,mvmax2
!!                  depvv(ivh,ivs,ir)=depvv(ivh,ivs,ir) &
!!     &                 +depvva(ivh,ivs,ir,im)/vpp &
!!     &                 +wmxt(ir,im)*fmx0(ivh,ivs,ir)/vpp
!! !c     &                 +depvva(ivh,ivs,ir,im)/vpp/ssr
!! !c     &                 +wmxt(ir,im)*fmx0(ivh,ivs,ir)/vpp/ssr
!!               END DO
!!            END DO

          depprf(ir) = depprf(ir) +deppra(ir,im)
          deper(ir)  = deper(ir)  +depera(ir,im)
          depvh(ir) =  depvh(ir) +depvha(ir,im)
          depvr(ir) =  depvr(ir) +depvra(ir,im)
        END DO

!CDIR NODEP
        DO ith = 1, mtmax3
          deppt(ith)=deppt(ith)+deppta(ith,im)
          depet(ith)=depet(ith)+depeta(ith,im)
        END DO

!CDIR NODEP
        DO iph = 1, mpmax3
          deppp(iph)=deppp(iph)+depppa(iph,im)
          depep(iph)=depep(iph)+depepa(iph,im)
        END DO

        DO ixx = 1, mxmax
          DO iyy=1,mymax3
            deppx(ixx,iyy)=deppx(ixx,iyy)+deppxa(ixx,iyy,im)
            depex(ixx,iyy)=depex(ixx,iyy)+deppxa(ixx,iyy,im)
          END DO
        END DO

      END DO


!!      vmcm = 1.d2               ! 200710
      delr = 1.d0/dfloat(mrmax-1)
      DO iir=1,mrmax-1

        rrr = (float(iir)-0.5d0)*delr
        rra = rrr*a
        csr = cpi*(1.d2*a*delr)**2*(2.d0*dfloat(iir)-1.d0)
        cvol= 2.d0*cpi*R*1.d2*csr
        srcrt = rtptcl/cvol/itmax

        dprirt= depprf(iir)*srcrt
        dpeirt= deper(iir)*srcrt
        dvhirt= depvh(iir)*srcrt
        dvrirt= depvr(iir)*srcrt

        iirp=iir+mrmax
        dprirp= depprf(iirp)*srcrt
        dpeirp= deper(iirp)*srcrt
        dvhirp= depvh(iirp)*srcrt
        dvrirp= depvr(iirp)*srcrt

        iirm=iirp+mrmax
        dprirm= depprf(iirm)*srcrt
        dpeirm= deper(iirm)*srcrt
        dvhirm= depvh(iirm)*srcrt
        dvrirm= depvr(iirm)*srcrt

        WRITE(20,11000) iir, rrr, dpeirt, dvhirt, rra, cvol
        WRITE(30,11000) iir, rrr, dpeirp, dvhirp, rra, cvol
        WRITE(40,11000) iir, rrr, dpeirm, dvhirm, rra, cvol

!        write(30,11000) iir, rrr, dprirp, dpeirp, dvhirp, dvrirp
!        write(40,11000) iir, rrr, dprirm, dpeirm, dvhirm, dvrirm

      END DO

!$$$      do 853 iir =1,mrmax
!$$$      do 852 ivpl=1,mvmax2
!$$$c
!$$$         write(21,11001) (depvv(ivpl,j,iir),j=1,32)
!$$$c
!$$$         iirp=iir+mrmax
!$$$         write(31,11001) (depvv(ivpl,j,iirp),j=1,32)
!$$$c
!$$$         iirm=iirp+mrmax
!$$$         write(41,11001) (depvv(ivpl,j,iirm),j=1,32)
!$$$c
!$$$ 852  continue
!$$$ 853  continue
!$$$c
!$$$      do 854 i=1,iit
!$$$         write(22,11002) tmhst(i),enghst(i),vhhst(i)
!$$$     &        ,enghsp(i),enghsm(i),vhhstp(i),vhhstm(i)
!$$$ 854  continue
!$$$c
!$$$*VDIR NODEP
!$$$         do 200 i=1,50
!$$$            do 202 j=1,3
!$$$               seir(i,j)=0.d0
!$$$               seer(i,j)=0.d0
!$$$ 202        continue
!$$$ 200     continue
!$$$c
!$$$      do 210 m=1,mm2
!$$$*VDIR NODEP
!$$$         do 212 i=1,50
!$$$            do 214 j=1,3
!$$$               seir(i,j)=seir(i,j)+eir(i,j,m)
!$$$               seer(i,j)=seer(i,j)+eer(i,j,m)
!$$$ 214        continue
!$$$ 212     continue
!$$$c
!$$$ 210  continue
!$$$c
!$$$      do 220 i=1,50
!$$$c
!$$$         rrr = (float(i)-0.5d0)/50.d0
!$$$         csr = cpi*(1.d2*a/50.d0)**2*(2.d0*dfloat(i)-1.d0)
!$$$         cvol= 2.d0*cpi*R*1.d2*csr
!$$$         engrt = rtptcl/cvol
!$$$c
!$$$         eirss1 = seir(i,1)*engrt
!$$$         eerss1 = seer(i,1)*engrt
!$$$         eirss2 = seir(i,2)*engrt
!$$$         eerss2 = seer(i,2)*engrt
!$$$         eirss3 = seir(i,3)*engrt
!$$$         eerss3 = seer(i,3)*engrt
!$$$         write(23,11003) rrr, eerss1, eirss1, eerss2,
!$$$     &                        eirss2, eerss3, eirss3
!$$$ 220  continue
!$$$c
!$$$      do 230 ith=1,mtmax
!$$$         tht = dfloat(ith-1)/dfloat(mtmax-1)
!$$$         write(24,11004) tht, deppt(ith), depet(ith)
!$$$         ithp = ith+mtmax
!$$$         write(34,11004) tht, deppt(ithp), depet(ithp)
!$$$         ithm = ithm+mtmax
!$$$         write(44,11004) tht, deppt(ithm), depet(ithm)
!$$$ 230  continue
!$$$c
!$$$      do 240 iph=1,mpmax
!$$$         pht = dfloat(iph-1)/dfloat(mpmax-1)
!$$$         write(25,11004) pht, deppp(iph), depep(iph)
!$$$         iphp=iph+mpmax
!$$$         write(35,11004) pht, deppp(iphp), depep(iphp)
!$$$         iphm=iphp+mpmax
!$$$         depppm=-deppp(iphm)
!$$$         depepm=-depep(iphm)
!$$$         write(45,11004) pht, depppm, depepm
!$$$ 240  continue
!$$$c
!$$$      do 250 iyy=1,mymax
!$$$         write(26,11005) (deppx(j,iyy),j=1,100)
!$$$         write(27,11005) (depex(j,iyy),j=1,100)
!$$$         iyyp=iyy+mymax
!$$$         write(36,11005) (deppx(j,iyyp),j=1,100)
!$$$         write(37,11005) (depex(j,iyyp),j=1,100)
!$$$         iyym=iyyp+mymax
!$$$         write(46,11005) (deppx(j,iyym),j=1,100)
!$$$         write(47,11005) (depex(j,iyym),j=1,100)
!$$$ 250  continue
!$$$c---------------------------------------c


      WRITE(6,999)

      WRITE(6,1001) (npt(i),depths(i),depehs(i),depihs(i),aveng(i), &
          avengs(i),dtime(i),nloss(i),etlos(i),prf(i),i=1,iit)


  666 FORMAT(10e12.4)
  668 FORMAT(/,5x,'**distribution of loss particles**',/)
  999 FORMAT(/,7x,'npt',5x,'egy dep',6x,'egy dep(e)',5x,'egy dep(i)', &
     &4x,'ave egy(test)',2x,'ave egys(test)',4x,'dtime',7x,'nloss', &
     &5x,'etloss',9x,'prf',/)

 1001 FORMAT(5x,i5,6e15.5,i5,2e15.5)
 1002 FORMAT(10x,'n=',i5,5x,'eng at loss =',e12.3, &
     &5x,'flux surf. at stop=',e12.3)
 1003 FORMAT(10x,'n=',i5,5x,'eng at singular=',e12.3, &
     &5x,'flux surf. at singular=',e12.3,5x,'type=',i3)
 1000 FORMAT(1h ,i7,1p5e13.4)
11000 FORMAT(1h ,i4,5e13.4)
11001 FORMAT(1h ,32e13.4)
11002 FORMAT(1h ,7e13.4)
11003 FORMAT(1h ,7e13.4)
11004 FORMAT(1h ,3e13.4)
11005 FORMAT(1h ,100e13.4)

      DEALLOCATE(depep,deper,depet,deppp,depprf,deppt,depvh,depvr,depex,deppx)
      CALL deallocate_restrt22
      CALL deallocate_through
!      STOP
      RETURN
      END SUBROUTINE


!======================================================================c
      SUBROUTINE setprm
!======================================================================c

  USE mcnmod, ONLY : enorm, ioptn, istrt, iswrd, it0, ix, iy, jplot, maxnpl
  IMPLICIT NONE

      NAMELIST/nmoptn/ioptn,jplot,istrt,enorm

!---------------------------------------c
!     read  nmoptn
!         ioptn,  jplot, istrt, enorm
!---------------------------------------c

      READ (5,nmoptn)
      WRITE(6,nmoptn)

!---------------------------------------c
!            set parameters
!---------------------------------------c

      CALL consts
      CALL input(ioptn)


      maxnpl = 0
      iswrd  = 0
      ix=959321
      iy=19825

      it0=1

      RETURN
      END SUBROUTINE


!======================================================================c
      SUBROUTINE setfld
!======================================================================c
!     set field data
!----------------------------------------------------------------------c

  USE mcnmod, ONLY : a, b0, bb0, bco, c1bf, c1et, c1g, c1i, c2bf, c2et, c2g, c2i,  &
    c3i, c3bf, c3et, c3g, cug, cui, depepa, depera, depeta, depexa, depppa, deppra,&
    deppta, deppxa, depvha, depvra, depvva, eot, kmsh, mm2, mmx, mpmax3, mrmax3,   &
    mtmax3, mvmax, mvmax2, mxmax, mymax3, psi, psia, wmxt
  IMPLICIT NONE
  INTEGER(4):: i, icon, im, iph, ir, ith, iv1, iv2, ixx, iyy, j
  REAL(8):: ds1, ds2, ds3, ds4, dds12, dds34,  &
    s(kmsh+1), wkbf(kmsh+1), wkc1(kmsh+1),wkc2(kmsh+1),wkc3(kmsh+1), dy(2)

!-------------------------------------c
!          read vmec-nwboz field
!-------------------------------------c

      CALL iodisk

!-------------------------------------c

      b0   = bb0
      psia = psi(kmsh+1)
!!      sa   = sqrt(psi(kmsh+1))  ! 200710
      IF (psia.gt.0.d0) THEN
        a = sqrt(psia*2.d0/b0)
      ELSE
        a = sqrt(-psia*2.d0/b0)
      END IF

!-------------------------------------c

  555 CONTINUE

      DO j=1,kmsh+1
!       s(j) = sqrt( psi(j)/psia )
        s(j) = sqrt( psi(j) )
      END DO

!     ds1 = sqrt(psi(2)/psia)  -sqrt(psi(1)/psia)
!     ds2 = sqrt(psi(3)/psia) -sqrt(psi(2)/psia)
      ds1 = sqrt(psi(2)) -sqrt(psi(1))
      ds2 = sqrt(psi(3)) -sqrt(psi(2))
      dds12 = (ds1 +ds2)/2.d0

!     ds3 = sqrt(psi(kmsh)/psia) -sqrt(psi(kmsh -1)/psia)
!     ds4 = sqrt(psi(kmsh +1)/psia) -sqrt(psi(kmsh)/psia)

      ds3 = sqrt(psi(kmsh)) -sqrt(psi(kmsh -1))
      ds4 = sqrt(psi(kmsh +1)) -sqrt(psi(kmsh))
      dds34 = (ds3 +ds4)/2.d0

      DO i=1,mmx+1

        DO j=1,kmsh+1
          wkbf(j)=bco(j,i)
        END DO

        dy(1) = ( (wkbf(3) -wkbf(2))/ds2 -(wkbf(2) -wkbf(1))/ds1 )/dds12
        dy(2) = ( (wkbf(kmsh+1) -wkbf(kmsh))/ds4   -(wkbf(kmsh) -wkbf(kmsh-1))/ds3 )/dds34

        CALL dinspl( s,wkbf, dy, kmsh+1, wkc1, wkc2, wkc3, icon)

        DO j=1, kmsh+1
          c1bf(j,i) =  wkc1(j)
          c2bf(j,i) =  wkc2(j)
          c3bf(j,i) =  wkc3(j)
        END DO

      END DO

!--------------------------------------------c

!     make the spline table

!--------------------------------------------c

      dy(1) = ((eot(3) -eot(2))/ds2 -(eot(2) -eot(1))/ds1 )/dds12
      dy(2) = ((eot(kmsh+1) -eot(kmsh))/ds4  -(eot(kmsh) -eot(kmsh-1))/ds3 )/dds34

      CALL dinspl( s, eot, dy, kmsh+1, c1et, c2et, c3et, icon)

      dy(1) = ((cui(3) -cui(2))/ds2 -(cui(2) -cui(1))/ds1 )/dds12
      dy(2) = ((cui(kmsh+1) -cui(kmsh))/ds4  -(cui(kmsh) -cui(kmsh-1))/ds3 )/dds34

      CALL dinspl( s, cui, dy, kmsh+1, c1i, c2i, c3i, icon)

      dy(1) = ((cug(3) -cug(2))/ds2 -(cug(2) -cug(1))/ds1 )/dds12
      dy(2) = ((cug(kmsh+1) -cug(kmsh))/ds4  -(cug(kmsh) -cug(kmsh-1))/ds3 )/dds34

      CALL dinspl( s, cug, dy, kmsh+1, c1g, c2g, c3g, icon)

!--------------------------------------------c
!      clear   deppra(mrmax,mm2),
!              depvha(mrmax,mm2),
!              depvva(mvmax2,mvmax,mrmax,mm2)
!--------------------------------------------c


      DO im=1,mm2
        DO ir=1,mrmax3
          deppra(ir,im)=0.d0
          depera(ir,im)=0.d0
          depvha(ir,im)=0.d0
          depvra(ir,im)=0.d0
          wmxt  (ir,im)=0.d0

          DO iv1=1,mvmax
            DO iv2=1,mvmax2
              depvva(iv2,iv1,ir,im)=0.d0
            END DO
          END DO
        END DO

        DO ith=1,mtmax3
          deppta(ith,im)=0.d0
          depeta(ith,im)=0.d0
        END DO

        DO iph=1,mpmax3
          depppa(iph,im)=0.d0
          depepa(iph,im)=0.d0
        END DO

        DO ixx=1,mxmax
          DO iyy=1,mymax3
            deppxa(ixx,iyy,im)=0.d0
            depexa(ixx,iyy,im)=0.d0
          END DO
        END DO
!
      END DO

      RETURN
      END SUBROUTINE


!======================================================================c
      SUBROUTINE setptl
!======================================================================c

!     set particle distribution
!     --- particle position and velocity  ----

!----------------------------------------------------------------------c

  USE mcnmod, ONLY : istrt
  IMPLICIT NONE
!
!-------------------------------------c
!   set particle position and velocity
!-------------------------------------c

      IF (istrt.eq.0) THEN
        CALL vinit
      ELSE
        CALL restrt(1)
      END IF

!-------------------------------------c
      RETURN
      END SUBROUTINE

!======================================================================c
      SUBROUTINE setcx
!======================================================================c

!     set charge exchange data

!----------------------------------------------------------------------c

  USE mcnmod, ONLY : denn0, engn0, mxcx
  IMPLICIT NONE
  INTEGER(4):: i, j
  REAL(8):: zr, zrms, denein(mxcx), denion(mxcx), tein(mxcx), tiin(mxcx)


!-------------------------------------c
!   read cx data
!-------------------------------------c

      PRINT 1775
      PRINT 1776

      DO i=1,mxcx
        READ(16,1777) j,zr,denein(i),denion(i),tein(i),tiin(i), denn0(i),engn0(i),zrms
        PRINT 1777, j,zr,denein(i),denion(i),tein(i),tiin(i), denn0(i),engn0(i),zrms
      END DO

 1775 FORMAT(2x,5hpoint,5x,6hradius,8x,2hne,12x,2hni,12x,2hte, 12x,2hti, &
     &  12x,2hn0,12x,2he0,12x,8hrms dev.)
 1776 FORMAT(2x,117(1h-))
 1777 FORMAT(1x,i5,9(3x,1pe11.4))

!-------------------------------------c
      RETURN
      END SUBROUTINE


!======================================================================c
      SUBROUTINE mvptl
!======================================================================c

!     move particles

!----------------------------------------------------------------------c

  USE mcnmod, ONLY : a, cpi2, dt, iit, isw, it, jit, ntest, psia, time, yc
  IMPLICIT NONE
  INTEGER(4):: n
  REAL(8):: sr

      IF (it.eq.1) THEN
        iit=0
        jit=0
!         dw=tiev*4.d0
!!         iadsw = 0  ! 200710
      END IF
!
!------------------------------c
!     call Runge-Kutta-Huta
!------------------------------c

      CALL rkhn(yc,time,dt)


!--------------------------------------------c
!        check the escaping particles
!--------------------------------------------c
!CDIR NODEP
      DO n=1,ntest
        IF(yc(n,1).lt.0.d0) THEN
          yc(n,1)= -yc(n,1)
          yc(n,2)=  yc(n,2)+cpi2
        ENDIF

        sr=sqrt(yc(n,1)/psia)*a
!
        IF ((isw(n).ne.0).and.(sr.ge.a)) THEN
          isw(n) = 0
        ENDIF
      END DO

      RETURN
      END SUBROUTINE


!======================================================================c
      SUBROUTINE colptl
!======================================================================c
!     move particles
!----------------------------------------------------------------------c

  USE mcnmod, ONLY : adinv, ddt, dee, dei, der, dir, dt, eer, eir, eng, es, icoll, &
    imsh, isw, mm2, tdep, tdepe, tdepi, yc
  IMPLICIT NONE
  INTEGER(4):: i, ic, j, m

!-----------------------------------------------------------------c
!     collide n particles with background ions and electrons
!                by monte carlo method
!-----------------------------------------------------------------c
!
      IF(icoll.eq.0) RETURN
!
!CDIR NODEP
      DO m=1,mm2
        DO j=1,3
          DO i=1,imsh+1
            dir(i,j,m)=0.d0
            der(i,j,m)=0.d0
          END DO
        END DO
!
       dee(m)=0.d0
       dei(m)=0.d0
      END DO
!
      ddt = dt/ dfloat(icoll)
      DO ic=1,icoll
        CALL monten(yc,adinv,eng,es,isw,ddt,dee,dei,der,dir)
!        call montcx(yc,adinv,eng,es,isw,ddt,dee,dei,der,dir)
      END DO
!
!CDIR NODEP
      DO m=1,mm2
        DO j=1,3
          DO i=1,imsh+1
            eir(i,j,m) = eir(i,j,m) +dir(i,j,m)
            eer(i,j,m) = eer(i,j,m) +der(i,j,m)
          END DO
        END DO
      END DO


!CDIR NODEP
      DO m=1,mm2
        tdep(m)  = tdep(m)  +dee(m) +dei(m)
        tdepe(m) = tdepe(m) +dee(m)
        tdepi(m) = tdepi(m) +dei(m)
      END DO
!
      RETURN
      END SUBROUTINE


!======================================================================c
      SUBROUTINE consts
!======================================================================c
!        set physical constants
!----------------------------------------------------------------------c

  USE mcnmod, ONLY : cecgs, ceps0, charge, ckb, clight, cme, cmp, cmyu0, cpi, cpi2
  IMPLICIT NONE

      NAMELIST/nmcnst/cpi,charge,cecgs,ckb,clight,cme,cmp,cmyu0,ceps0

!---------------------------------------------c
!     pi
      cpi  = atan2(0.d0, -1.d0)
      cpi2 = 2.d0*cpi

!     cpi=3.141592653589793238462643383279

!     charge
      charge=1.6021892d-19

!     boltzman constant
      ckb=charge

!     elementaly charge in cgs unit
      cecgs=4.803242d-10

!     speed of light
      clight=2.99792458d+8

!     electron mass
      cme=9.109534d-31

!     proton mass
      cmp=1.6726485d-27

!     permeability
      cmyu0=4.*cpi*1.d-7

!     permittivity
      ceps0=1./(cmyu0*clight**2)

!---------------------------------------------c

      WRITE(6,nmcnst)

      RETURN
      END SUBROUTINE

!======================================================================c
      SUBROUTINE input(ioptn)
!======================================================================c

  USE mcnmod, ONLY: a, ai, ane, anecgs, anem1, anem2, anew, ani, anicgs, at, b0,   &
    bres, cecgs, cfb, cfe0, cfi0, charge, ckb, cme, cmi, cmp, cpi, crpe, ddt, dt,  &
    el, em, emin, enorm, eot0, eota, epsa, epsh, erf, esp0, g, ich, icoll, iout,   &
    itmax, jout, kout, nout, nsw, pna0, pna1, pne0, pne1, psia, r, sgm, tchg, te,  &
    teev,teevw, tem1, tem2, tew, ti, tiev, tievw, tim1, tim2, tiw, tmass,ve,vi,vnorm
  IMPLICIT NONE
  INTEGER(4),INTENT(IN):: ioptn
  REAL(8):: cfbe, cfbi, cgcme, cgcmi, cgczi, cgtchg, cgte, cgti, cgtmss, czi, rchg, &
    zeff, zi, zt

      NAMELIST/nmtest/at,zt,ai,zi
      NAMELIST/nminp1/b0,r,a,epsh,sgm,el,em,eot0,eota,esp0,erf,bres
      NAMELIST /nmcomp/ dt,emin,itmax,icoll,ddt, iout,jout,nout,kout,nsw,ich
!     namelist/nmsngl/e0ev,pitch0,sr0,thta0,fai0
      NAMELIST/inputd/teev,teevw,tem1,tem2,tiev,tievw,tim1,tim2,ane,anew,anem1,anem2,zeff
!      namelist/nmadam/ depsa, depsr

!----------------------------------------c

!     at : test particle mass ratio
!     zt : test particle charge number
!     ai : balk ion      mass ratio
!     zi : balk ion      charge number

!----------------------------------------c
!     set default parameters

      at=1.
      zt=1.

      ai=2.
      zi=1.

!----------------------------------------c
!     read  parameters

      READ(5,nmtest)
      WRITE(6,nmtest)

!----------------------------------------c

!     tmass : test particle mass
!     tchg  : test particle charge
!     cmi   : balk ion      mass
!     czi   : balk ion      charge

!     dt    : time step size
!     emin  : deth particle energy(ev)
!     itmax : maximum time step

!----------------------------------------c
!     set default parameters

!$$$c for electron test particle
!$$$c
!$$$      tmass= cme*at
!$$$      tchg= -charge*zt
!$$$      cgtmss =  cme*at*1.d3
!$$$      cgtchg = -cecgs*zt

! for ion test particle
      tmass  = cmp*at
      tchg   = charge*zt
      cgtmss = cmp*at*1.d3
      cgtchg = cecgs*zt

      cmi= cmp*ai
      czi= charge*zi
      cgcmi = cmp*ai*1.d3
      cgczi = cecgs*zi
      cgcme = cme*1.d3
      crpe=cmp/cme

      dt=1.d-8
      emin=50.d0
      itmax=1000
      icoll=1
      iout=100
      jout=100
      nout=100
      kout=100
      nsw=2

!----------------------------------------c
!     read paramters

      READ(5,nmcomp)
      WRITE(6,nmcomp)

!----------------------------------------c

!     teev : ele. temperature(ev)
!     tiev : ion  temperature(ev)
!     ane  : ele. density
!     zeff : effective charge

!----------------------------------------c
!     set default parameters

      teev=1.e3
      tiev=1.e3
      ane=5.e19
      zeff=1.

!----------------------------------------c
!     read paramters

      READ(5,inputd)
      WRITE(6,inputd)

!----------------------------------------c
!     error limit of adams methods

!     depsa : absolute error limit
!     depsr : relative error limit
!----------------------------------------c
!     read paramters

!      read(5,nmadam)
!      write(6,nmadam)

!----------------------------------------c
!     te     : ele. temperature(j)
!     ti     : ion  temperature(j)
!     anecgs : ele. density in cgs unit
!     ani    : ion  density
!     anicgs : ion  density in cgs unit

!----------------------------------------c

      te=charge*teev
      ti=charge*tiev
      tew = charge*teevw
      tiw = charge*tievw

      cgte=te*1.d7
      cgti=ti*1.d7

      anecgs=ane*1.d-6
      ani=ane/zi
      anicgs=ani*1.d-6

!----------------------------------------c
!     coulomb logarithm
!----------------------------------------c
!     if(teev.lt.50.)
!    &   clog=23.4-1.15*alog10(anecgs)+3.45*alog10(teev)
!     if(teev.ge.50.)
!    &   clog=25.3-1.15*alog10(anecgs)+2.3*alog10(teev)

!----------------------------------------c
!     thermal velocities
!----------------------------------------c

      ve=sqrt(2.d0*te/cme)
      vi=sqrt(2.d0*ti/cmi)
      vnorm=sqrt(2.*ckb*enorm/tmass)


!----------------------------------------c
!     collision frequencies
!----------------------------------------c
!ecgs
      rchg = cecgs/charge

      cfb=1.d0/3.d7*sqrt(2.d0/ai)*zi**4*anicgs/tiev**1.5d0

      cfbi=4.d0/3.d0*sqrt(cpi/cgcmi)*(cgcmi/cgtmss)**2 &
     &      *cgczi**2*cgtchg**2*anicgs/cgti**1.5d0
!     &                *cgtchg**2*anicgs/cgti**1.5d0 &
      cfi0=1.5d0*sqrt(0.5d0*cpi)*cfbi*zeff

      cfbe=4.d0/3.d0*sqrt(cpi/cgcme)*(cgcme/cgtmss)**2 &
     &      *cecgs**2*cgtchg**2*anecgs/cgte**1.5d0
      cfe0=1.5d0*sqrt(0.5d0*cpi)*cfbe

      PRINT 666, cfb,cfbi,cfi0,cfbe,cfe0

 666  FORMAT(' cfb,cfbi,cfi0,cfbe,cfe0=',5e12.3)
!----------------------------------------c
!     critical velocity and
!                 slowing down time
!----------------------------------------c

!     z1=zeff*at/ai
!     z2=zeff/z1
!     vc=(3.*sqrt(cpi)/4.*z1/(at*crpe))**(1./3.)*ve
!     ec=0.5*tmass*vc**2
!     ecte=ec/te
!     gte=2.*cpi*zt**2*cecgs**4*anecgs*clog/(tmass*1.e3)**2
!     taus=(vc*100.)**3/(2.*gte*z1)

!----------------------------------------c
!     parameters for a helical system
!----------------------------------------c
!     set default parameters

!     magnetic field strength
      b0=4.

!     major radii
      r=5.

!     minor radii
      a=0.5

      epsh=0.3
      sgm=0.

!     coil parameter
      el=2.
      em=10.

!     rotational transform number
      eot0=0.5
      eota=1.5

      esp0=0.

      pne0=0.
      pne1=0.
      pna0=0.
      pna1=0.

!     rf strength ======> erf [v/cm]
      erf=200.

!----------------------------------------c
!     read parameters

      READ (5,nminp1)
      WRITE(6,nminp1)

!----------------------------------------c

      epsa=a/r
      psia=0.5*b0*a**2
      g=b0*r
      IF (b0.lt.0.d0) THEN
         b0 = -b0
         PRINT *, '*** b0 is replaced by -b0.'
      END IF

!----------------------------------------c

!ccs  if(ioptn.ne.0) return

!         parameters for a single test particle

!     e0ev=1.0d+4
!     pitch0=90.
!     sr0=0.0
!     thta0=0.
!     fai0=0.

!     read (5,nmsngl)
!     write(6,nmsngl)

!      print *, '**** in sub. input ***'
!      print *, 'Ti,Te,Vi,Ve= ', ti,te,vi,ve

      RETURN
      END SUBROUTINE

!======================================================================c
      SUBROUTINE bfld2(y,b)
!======================================================================c
!     give the magnetic strength
!----------------------------------------------------------------------c

  USE mcnmod, ONLY : bco, c1bf, c2bf, c3bf, cm, cn, kmsh, mmx, psi, psia
  IMPLICIT NONE
  REAL(8),INTENT(INOUT):: y(4)
  REAL(8),INTENT(OUT)  :: b
  INTEGER(4):: i, jp
  REAL(8):: bci, cx, ds, sj, sp, th


      IF (y(1).le.0.d0) THEN
        y(1) = 0.d0
        jp   = 1
      ELSE IF(y(1).ge.psia) THEN
        y(1) = psia
        jp   = kmsh +1
      ELSE
!
        jp = int(y(1)*kmsh/psia +1.d0)
!
      END IF
!
      sp = sqrt(y(1))
      sj = sqrt( psi(jp))
      ds = sp -sj

      b  = 0.d0
!
      DO i=1, mmx+1
        bci = bco(jp,i) +(  c1bf(jp,i) +(c2bf(jp,i) +c3bf(jp,i)*ds)*ds)*ds
!
        th = cm(i)*y(2) -cn(i)*y(3)
        cx = cos(th)
        b  = bci*cx+b

      END DO

      RETURN
      END SUBROUTINE


!======================================================================c
      SUBROUTINE bfld2n(yc,b)
!======================================================================c


!     give the magnetic strength

!----------------------------------------------------------------------c

  USE mcnmod, ONLY : bco, c1bf, c2bf, c3bf, cm, cn, kmsh, mmx, ntest, psi, psia
  IMPLICIT NONE
  REAL(8),INTENT(IN) :: yc(ntest,4)
  REAL(8),INTENT(OUT):: b(ntest)
  INTEGER(4):: i, ip, jp(ntest)
  REAL(8):: bci(ntest), cx(ntest), ds(ntest), sj(ntest), sp(ntest),    &
            th(ntest), y1(ntest)

!CDIR NODEP
      DO ip=1,ntest

        IF(yc(ip,1).ge.0.d0) THEN
          y1(ip)=yc(ip,1)
        ELSE
          y1(ip)=0.d0
        END IF
!
        jp(ip) = int(y1(ip)*kmsh/psia) +1
!
        IF(jp(ip).gt.kmsh) THEN
          jp(ip) = kmsh
          y1(ip) = psia
        END IF
!
        sp(ip) = sqrt(y1(ip))
        sj(ip) = sqrt( psi(jp(ip)) )
        ds(ip) = sp(ip) -sj(ip)
!
        b(ip)  = 0.d0
!
      END DO

      DO i=1, mmx+1
!CDIR NODEP
        DO ip=1,ntest
          bci(ip)= bco(jp(ip),i) +(c1bf(jp(ip),i)  +( c2bf(jp(ip),i)  &
                  +c3bf(jp(ip),i)*ds(ip))*ds(ip))*ds(ip)
          th(ip) = cm(i)*yc(ip,2) -cn(i)*yc(ip,3)
          cx(ip) = cos(th(ip))
          b(ip)  = bci(ip)*cx(ip)+b(ip)
        END DO
!
      END DO

!
      RETURN
      END SUBROUTINE


!======================================================================c
      SUBROUTINE vinit
!======================================================================c

  USE mcnmod , ONLY : a, adinv, ane, ane, bres, ckb, clog, cpi, e0ev, e0int, ekin0, &
    eng, esp0, esp0, fai0, ips, ispsi, ispsi, isswt, isw, ntest, pecrh, pitch0, psia,&
    ram, rswt, rtptcl, sr0, swt, taus, tchg, teev, thta0, time, tmass, vel, y, yc
  IMPLICIT NONE
  INTEGER(4)::ichk, ipt9, ipt9a, ipt9b, ismpsi, ix, j, k1, n, nt1, nt2, nt3
  REAL(8)::  b, crate, e01, e02, engsr, espi, fai, psi1, rtp1, rtp2, rtp3, rtpt9,   &
    smamu, smeng, smlmd, smphi, smpsi, smroh, smtht, thta, v0, bhflrt(3), bhfprt(3)

      NAMELIST/nmix/  ix
      NAMELIST/ne0/   e0int, pecrh
      NAMELIST/ncrt/  crate


      time = 0.d0

!---------------------------------------c
!     read  ix
!---------------------------------------c

      READ (5,nmix)
      WRITE(6,nmix)

!---------------------------------------c
!     read  e0int
!---------------------------------------c

      READ (5,ne0)
      WRITE(6,ne0)

!---------------------------------------c
!     read  ncrt
!---------------------------------------c

      READ (5,ncrt)
      WRITE(6,ncrt)

!------------------------------------------c
!               ioptn = 0
!------------------------------------------c

!        call ranu2n(ix,rndn1,ntest,icon)
!        call ranu2n(ix,rndn2,ntest,icon)
!        call ranu2n(ix,rndn3,ntest,icon)

      engsr =0.d0
      ismpsi = 0

!!!      dntest = 180.0d0/ntest  ! 200710

      nt1 = 0
      nt2 = 0
      nt3 = 0

      e01 = 0.9d0*e0int
      e02 = 0.9d0*0.5d0*e0int

!-----------------------------------c
!     count the particle number in out9 file
!-----------------------------------c

      ipt9=0
 40   CONTINUE
      READ(9,7001,END=50) k1,smpsi,smtht,smphi,smroh,smamu,smeng,smlmd
      IF (smpsi.le.1.d-4) THEN
        GO TO 40
      END IF
      ipt9=ipt9+1
      GO TO 40
 50   CONTINUE

      rtpt9=dfloat(ipt9)/dfloat(ntest)
      PRINT *, '*** partile number in out9 = ', ipt9
      PRINT *, '*** reading rate = ', rtpt9

      REWIND(9)

      ipt9a=0
      DO n=1,ntest
!cccccccc

 12     CONTINUE
        READ(9,7001) k1,smpsi,smtht,smphi,smroh,smamu,smeng,smlmd
 7001   FORMAT(i4,7d24.6)

        IF (smpsi.le.1.d-4) THEN
          ismpsi = ismpsi +1
          GO TO 12
        END IF
        ipt9a=ipt9a+1
        ipt9b=int(rtpt9*dfloat(n))
        IF (ipt9a.ne.ipt9b) THEN
          GO TO 12
        END IF

!-------------------------------c
!   psi_a<0 case
!-------------------------------c
        IF (ispsi.eq.1) THEN
          smtht=-smtht
          smphi=-smphi
          smlmd=-smlmd
        ELSE IF (ispsi.eq.2) THEN
          smtht=-smtht
          smphi=-smphi
          smlmd=-smlmd
        END IF
!-------------------------------c
!-------------------------------c

        WRITE(6,7002) n,smpsi,smtht,smphi,smroh,smamu,smeng,smlmd
 7002   FORMAT(i4,7d12.3)

        e0ev(n)   = smeng
        vel(n)    = sqrt(2.d0*e0ev(n)*ckb/tmass)
!...co
!           ram(n)    = -smlmd
!ccc
!...cnt
        ram(n)    = smlmd
!cc
        sr0(n)    = sqrt(smpsi)*a
        pitch0(n) = acos(smlmd)*180.d0/cpi
        thta0(n)  = smtht
        fai0(n)   = smphi

        psi1 = smpsi*psia
        thta = smtht
        fai  = smphi
        v0   = vel(n)

!cccccccc
!           print *, 'n,vppn,vpln,swn=',n,vppn,vpln,swn
!cc
!
        swt(n) = 1.d0
        engsr   = engsr +swt(n)*e0ev(n)
!
        IF (e0ev(n).gt.e01) THEN
          isswt(n) = 1
          nt1 = nt1+1
        ELSE IF (e0ev(n).gt.e02) THEN
          isswt(n) = 2
          nt2 = nt2+1
        ELSE
          isswt(n) = 3
           nt3 = nt3+1
        END IF
!
!cc
!...for ion test particle
        clog=24.0d0-dlog(sqrt(ane*1.d-6)/(teev))
        taus = teev**(1.5d0)/(1.6d-15*ane*clog)
        rswt(n) = crate*(taus/2.d0)

        ichk=mod(n,200)
        IF (ichk.eq.0) THEN
          PRINT *, 'n, clog, taus,rswt=',n, clog,taus,rswt(n)
        END IF

!cccc
        yc(n,1) = psi1
        yc(n,2) = thta
        yc(n,3) = fai
!
        y(1) = psi1
        y(2) = thta
        y(3) = fai
!
        CALL bfld2(y,b)
!
        yc(n,4) = tmass*v0*ram(n)/tchg/b
        espi    = esp0*(1.d0 -yc(n,1)/psia) *(1.d0 -yc(n,1)/psia)
        adinv(n)= 0.5d0*tmass*v0*v0*(1.d0 -ram(n)**2)/b
        eng(n)  = e0ev(n)
        ekin0(n)= e0ev(n)
!
!     print *, 'eng(n)=', n, eng(n)
!     print *, yc(n,1),yc(n,2),yc(n,3),yc(n,4)
!
        isw(n)  = 1
!
        ips(n)  = -1
        IF (b.ge.bres) ips(n)= 1
!
      END DO
!
      READ(15,7018 ) (bhflrt(j), j=1,3)
      READ(15,7018 ) (bhfprt(j), j=1,3)
 7018 FORMAT(3f15.4)

!      rtptcl = pecrh/(engsr*ckb)
!     rtptcl =  (1.d0-bhflrt(3))*pecrh/ntest

      rtp1   =  (1.d0-bhflrt(3))*bhfprt(3)*pecrh/nt1
      if (nt2.ne.0) then
        rtp2   =  (1.d0-bhflrt(2))*bhfprt(2)*pecrh/nt2
      else
        rtp2 =0
      end if
      if (nt3.ne.0) then
        rtp3   =  (1.d0-bhflrt(1))*bhfprt(1)*pecrh/nt3
      else
        rtp3 =0
      end if
!      rtp2   =  (1.d0-bhflrt(2))*bhfprt(2)*pecrh/nt2
!      rtp3   =  (1.d0-bhflrt(1))*bhfprt(1)*pecrh/nt3

      rtptcl =  1.d0

      DO n=1,ntest
        IF (isswt(n).eq.1) THEN
          swt(n)=rtp1
        ELSE IF (isswt(n).eq.2) THEN
          swt(n)=rtp2
        ELSE
          swt(n)=rtp3
        END IF
      END DO
!
      PRINT *, '  nt1,  nt2,  nt3 = ' ,  nt1,  nt2,  nt3
      PRINT *, ' rtp1, rtp2, rtp3 = ' , rtp1, rtp2, rtp3
      PRINT *, ' pecrh, engsr, rtptcl = ',pecrh,engsr,rtptcl
      PRINT *, ' ismpsi = ', ismpsi
!
 7010 FORMAT(3e18.5)
!
      RETURN
      END SUBROUTINE


!=============================================c
      SUBROUTINE restrt(irs)
!=============================================c

  USE mcnmod , ONLY : adinv, adiv, aveng, aveng1, avengs, avens1, b, bf, db, dbfi,   &
    dbth, ddt, dee, dei, dep, dep1, depe, depe1, depehs, depepa, depera, depeta,     &
    depexa, depi, depihs, depppa, deppra, deppta, deppxa, depths, depvha, depvra,    &
    depvva, der, dir, dt, dtime, dtime1, e0ev, e0int, echeck, eeer, eer, eiir, einit,&
    eir, ekin0, eler, elir, eloss, emin, eng, enghsm, enghsp, enghst, es, etlos, f,  &
    fai0, fderf, fdrf, ffm, ffp, fl, fld, flrf, fmx0, fs, ich, icoll, ifout, iit,    &
    iout, ips, ips1, iran, isswt, isw, iswrd, it, it0, ix, iy, jit, jobno, jout,     &
    kout, loss, maxnpl, nloss, nout, npl, npt, nsw, ntot, pang, pecrh, pgc, pitch0,  &
    prf, prf1, ram, rgc, rswt, rtptcl, sr0, sr1, swt, tang, tchg, tdep, tdepe, tdepi,&
    thta0, time, tmass, tmhst, vel, vhhst, vhhstm, vhhstp, vx, vy, wmxt, y, yc, yeer,&
    yiir, zgc
  USE mcnmod , ONLY : allocate_restrt1, deallocate_restrt1
  IMPLICIT NONE
  INTEGER(4),INTENT(IN):: irs
  REAL(8):: vel0

      NAMELIST/ne0/   e0int, pecrh
!-----------------------------------------------c

      CALL allocate_restrt1

      IF (irs.eq.0) THEN

        WRITE(12) time, it

        WRITE(12) b,db,dbth,dbfi
        WRITE(12) einit
        WRITE(12) adinv
        WRITE(12) dt,emin,icoll,ddt,iout,jout,nout,kout,nsw,ich,ifout
        WRITE(12) ekin0, echeck
        WRITE(12) eir,eer, dir,der, dei,dee
        WRITE(12) eeer, eiir, yeer, yiir
        WRITE(12) eler ,elir
        WRITE(12) f,fs,fl, ffp,ffm
        WRITE(12) isw,ips,ips1
        WRITE(12) eloss, loss, ntot
        WRITE(12) maxnpl, npl, jit, iit
!
        WRITE(12)  etlos,   dep,depe  &
             , depi, dtime, prf       &
             ,aveng,avengs, npt       &
             ,nloss
        WRITE(12)  dep1, depe1,  prf1 &
             ,aveng1,avens1,dtime1    &
             ,  flrf,  fdrf, fderf    &
             ,   fld
!
        WRITE(12) adiv,tmass,tchg,rtptcl
        WRITE(12) vel,ram,bf,eng, es
        WRITE(12) tdep, tdepe, tdepi
        WRITE(12) e0ev,pitch0,sr0,thta0,fai0,jobno
!
!
        WRITE(12) rgc, zgc   &
               ,  pgc, tang  &
               , pang
        WRITE(12) vel0,ram
        WRITE(12) vx ,vy
        WRITE(12) yc,y
!
        WRITE(12) swt,rswt
        WRITE(12) isswt
        WRITE(12) deppra, depera, deppta, depppa, &
                  depeta, depepa, depvra, depvha, &
                  depvva, deppxa, depexa
        WRITE(12) fmx0, wmxt
        WRITE(12) sr1
        WRITE(12) ix,iy,iswrd,iran
        WRITE(12) tmhst,vhhst,enghst ,enghsp,enghsm &
                 ,vhhstp,vhhstm ,depths,depehs,depihs

      ELSE

        READ(14) time, it0

        READ(14) b,db,dbth,dbfi
        READ(14) einit
        READ(14) adinv
        READ(14) dt,emin,icoll,ddt ,iout,jout,nout,kout,nsw,ich,ifout
        READ(14) ekin0, echeck
        READ(14) eir,eer, dir,der, dei,dee
        READ(14) eeer, eiir, yeer,yiir
        READ(14) eler ,elir
        READ(14) f,fs,fl, ffp,ffm
        READ(14) isw,ips,ips1
        READ(14) eloss, loss, ntot
        READ(14) maxnpl, npl, jit, iit
!
        READ(14) etlos,   dep,depe   &
               , depi, dtime, prf    &
                ,aveng,avengs, npt   &
                ,nloss
        READ(14)  dep1, depe1,  prf1  &
               ,aveng1,avens1,dtime1  &
               ,  flrf,  fdrf, fderf  &
               ,   fld
!
        READ(14) adiv,tmass,tchg,rtptcl
        READ(14) vel,ram,bf ,eng, es
        READ(14) tdep, tdepe, tdepi
        READ(14) e0ev,pitch0,sr0,thta0,fai0,jobno
!
        READ(14) rgc, zgc   &
              ,  pgc, tang  &
              , pang
        READ(14) vel0,ram
        READ(14) vx ,vy
        READ(14) yc,y
!
        READ(14) swt,rswt
        READ(14) isswt
        READ(14) deppra, depera, deppta, depppa, &
                 depeta, depepa, depvra, depvha, &
                 depvva, deppxa, depexa
        READ(14) fmx0, wmxt
        READ(14) sr1
        READ(14) ix,iy,iswrd,iran
        READ(14) tmhst,vhhst,enghst,enghsp,enghsm     &
                ,vhhstp,vhhstm ,depths,depehs,depihs

!$$$         read(14) tmhs2,vhhs2,enghs2
!$$$     &        ,engh2p,engh2m
!$$$     &        ,vhhs2p,vhhs2m
!$$$     &        ,depth2,depeh2,depih2
!$$$
!$$$c
!$$$         do 300 iii=1,10000
!$$$          tmhst(iii)=tmhs2(iii)
!$$$          vhhst(iii)=vhhs2(iii)
!$$$          enghst(iii)=enghs2(iii)
!$$$          enghsp(iii)=engh2p(iii)
!$$$          enghsm(iii)=engh2m(iii)
!$$$          vhhstp(iii)=vhhs2p(iii)
!$$$          vhhstm(iii)=vhhs2m(iii)
!$$$          depths(iii)=depth2(iii)
!$$$          depehs(iii)=depeh2(iii)
!$$$          depihs(iii)=depih2(iii)
!$$$ 300   continue

!---------------------------------------c
!     read  ne0
!---------------------------------------c

        READ (5,ne0)
        WRITE(6,ne0)



        PRINT *, ' time = ',time, '   it0= ',it0
        it0=it0+1

      END IF

      CALL deallocate_restrt1

      RETURN
      END SUBROUTINE


!======================================================================c
      SUBROUTINE monten(yc,adinv,eng,es,isw,dtau,dee,dei,ee,ei)
!======================================================================c

  USE mcnmod, ONLY : a, ane, anecgs, anem1, anem2, anew, cfe0, cfi0, ckb, clog, cme,&
    cmi, cpi2, depera, depexa, deppra, deppxa, depvha, depvra, e0int, ekin0, imsh,   &
    isswt, it, mm, mm2, mpmax, mrmax, mtmax, mvmax, mxmax, mymax, ntest, psia, rswt, &
    sr1, swt, tchg, te, teev, tem1, tem2, tew, ti, tim1, tim2, time, tiw, tmass, ve, &
    vi, wmxt
  IMPLICIT NONE
  INTEGER(4), INTENT(IN) :: isw(ntest)
  REAL(8),INTENT(IN)   :: dtau
  REAL(8),INTENT(INOUT):: adinv(ntest), dee(mm2), dei(mm2), yc(ntest,4),             &
    ee(imsh+1,3,mm2),ei(imsh+1,3,mm2)
  REAL(8),INTENT(OUT)  :: eng(ntest),es(ntest)
  INTEGER(4):: idsgn, iip, iir, iir2, iit, iix, iiy, ivvh, ivvs, m, n, n2,           &
    isf0(ntest), isf1(ntest), isf2(ntest), nt(ntest)
  REAL(8):: cctm, cfd, cfde, cfdi, cfee, cfei, dcfee, dcfei, dder, ddpr, ddswt, ddvh,&
    ddvr, ddw, dee1, dei1, delph, delr, delr2, delth, delxx, delyy, dvvmax, e1, e1d, &
    epp, etmax, ph, ran, rt0e, rt0i, rt2e, rt2i, th, v0, v1, vvh, vvmax, wrate, xx,yy,&
    bb(ntest), cfe1(ntest), cfi1(ntest), e0(ntest), fxe(ntest), fxe1(ntest),         &
    fxe2(ntest), fxi(ntest), fxi1(ntest), fxi2(ntest), gxe(ntest), gxi(ntest),       &
    r0(ntest), r1(ntest), rt1e(ntest), rt1i(ntest), sr(ntest), velh(ntest),          &
    velr(ntest), vels(ntest), xe(ntest), xi(ntest)
  REAL(8):: espr1=1.0d-11, espr2=1.0d-10
      save

!      data espr1,espr2/1.0d-11,1.0d-10/
!!!      data espr1,espr2/1.0d-6,1.0d-5/


!     etmax : maxmum value of energy in grid
      etmax= e0int*2.d0
!     vvmax : maxmum value of velocity in grid
      vvmax= sqrt(2.0*etmax*ckb/tmass)
!     dvvmax: size of grid
      dvvmax=vvmax/dfloat(mvmax-1)

      delr  = 1.d0/dfloat(mrmax-1)
      delr2 = 1.d0/dfloat(imsh)
      delth = cpi2/dfloat(mtmax-1)
      delph = cpi2/dfloat(mpmax-1)
      delxx = 2.d0/dfloat(mxmax-1)
      delyy = 2.d0/dfloat(mymax-1)

      CALL iseifu(isf0)
      CALL iseifu(isf1)
      CALL iseifu(isf2)

      CALL bfld2n(yc,bb)

!CDIR NODEP
      DO n=1,ntest
        sr(n)=sqrt(yc(n,1)/psia)
        IF (sr(n).le.0.d0) THEN
          sr(n)=0.1d-16
        ELSE IF (sr(n).ge.1.d0) THEN
          sr(n) = 1.d0
        END IF
!
        rt0e=(1.d0-tew/te)*(1.d0-sr(n)**tem1)**tem2+tew/te
        rt0i=(1.d0-tiw/ti)*(1.d0-sr(n)**tim1)**tim2+tiw/ti
        ran =(1.d0-anew/ane)*(1.d0-sr(n)**anem1)**anem2+anew/ane
!
!     coulomb logarithm
!
        epp=50.d0*ckb/te
!
!     if(rt0e.lt.epp) then
!     clog=23.4-1.15*dlog10(anecgs*ran)+3.45*dlog10(teev*rt0e)
!     else
!     clog=25.3-1.15*dlog10(anecgs*ran)+2.3*dlog10(teev*rt0e)
!     end if
!
        clog=24.0d0-dlog(sqrt(anecgs*ran)/(teev*rt0e))
!
        rt1e(n)=sqrt(1.d0/rt0e)
        rt1i(n)=sqrt(1.d0/rt0i)
        rt2e=ran*clog/rt0e**1.5d0
        rt2i=ran*clog/rt0i**1.5d0
        cfe1(n)=cfe0*rt2e
        cfi1(n)=cfi0*rt2i
!
        velh(n)=yc(n,4)*tchg/tmass*bb(n)
        vels(n)=sqrt(adinv(n)*bb(n)*2.d0/tmass)
        v0=sqrt(velh(n)**2+vels(n)**2)
        e0(n)=0.5d0*tmass*v0**2
        r0(n)=velh(n)/v0
        xe(n)=v0*rt1e(n)/ve
        xi(n)=v0*rt1i(n)/vi

!...for electrons
!         if (xi(n).lt.20.d0) then
!            xi(n)=20.d0
!         end if
!...for ions
        IF (xi(n).lt.1.d0) THEN
          xi(n)=1.d0
        END IF
!         ee00=e0(n)/ckb
!         if (ee00.lt.0.5d3) then
!           swt(n)=0.d0
!        end if

      END DO
!
!     print *,'Start gosakn'
      CALL gosakn(xe,fxe,fxe1,fxe2,gxe)
      CALL gosakn(xi,fxi,fxi1,fxi2,gxi)
!     print *,'End   gosakn'
!
!     print *,'Start roop 200'

      DO n2=1,mm
!CDIR NODEP
        DO m=1,mm2
          n = n2+mm*(m-1)

!----------------------------c
!     evaluate the weight
!----------------------------c
          cctm   = -(time/rswt(n))**8

          IF (cctm.lt.-20.d0) THEN
            wrate = 0.d0
          ELSE
            wrate = exp(cctm)
          END IF

          cfdi=  cfi1(n)*(fxi(n)-gxi(n))/xi(n)**3
          cfde=  cfe1(n)*(fxe(n)-gxe(n))/xe(n)**3
          cfd = cfdi+cfde
!
          cfee=2.d0*cfe1(n)*(tmass/cme)*gxe(n)/xe(n)
          cfei=2.d0*cfi1(n)*(tmass/cmi)*gxi(n)/xi(n)
          dcfee=-2.d0*cfe1(n)*(tmass/cme)/sqrt(te*e0(n))       &
                *(fxe2(n)+6.d0*gxe(n)) /(4.d0*xe(n)**2)*rt1e(n)
          dcfei=-2.d0*cfi1(n)*(tmass/cmi)/sqrt(ti*e0(n))       &
                *(fxi2(n)+6.d0*gxi(n)) /(4.d0*xi(n)**2)*rt1i(n)
!
!     print 665, n,cfdi,cfde,cfei,cfee
 665     FORMAT(i4,' cfd,cfe = ',4e12.3)
!
          IF(abs(abs(r0(n))-1.d0).lt.espr1) THEN
            r0(n)=r0(n)/abs(r0(n))*(1.d0-espr2)
          END IF
!
         r1(n)=r0(n)*(1.d0-cfd*dtau) &
     &        +isf0(n)*sqrt((1.d0-r0(n)**2)*cfd*dtau)
!
          IF(r1(n).ge.1.d0) THEN
            r1(n)=  1.d0-(r1(n)-int(r1(n)))
          ELSE IF (r1(n).le.-1.d0) THEN
            r1(n)= -1.d0-(r1(n)+int(abs(r1(n))))
          END IF
!
!     print *,cfee,dtau,dcfee
          dee1=-2.d0*cfee*dtau*(e0(n)-(1.5d0+e0(n)*dcfee/cfee)*te*rt0e)  &
               +2.d0*isf1(n)*sqrt(te*rt0e*e0(n)*cfee*dtau)
!     print *,'dee is evaluated.'
!     print *,cfei,dtau,dcfei
          dei1=-2.d0*cfei*dtau*(e0(n)-(1.5d0+e0(n)*dcfei/cfei)*ti*rt0i) &
               +2.d0*isf2(n)*sqrt(ti*rt0i*e0(n)*cfei*dtau)
!     print *,'dei is evaluated.'
!
!
          dei1=dei1*wrate
          dee1=dee1*wrate

          e1d = dee1+dei1
          e1  = e0(n)+e1d

!     e1=e0(n)+dee1
!     e1=e0(n)
!     print 666, n,e0(n),dee1,dei1,r0(n),dtau,cfd
 666     FORMAT(i4,' e0,dee1,dei1,r0,dtau,cfd= ',6e12.3)
!
          eng(n)  = e1/ckb
!
!-----------------------------------------------------c
!        Energy deposition profile (radail)

!                       ei(n,1) <---ion total f
!                       ei(n,2) <---ion f+
!                       ei(n,3) <---ion f-
!                       ee(n,1) <---ele. total f
!                       ee(n,2) <---ele. f+
!                       ee(n,3) <---ele. f-
!-----------------------------------------------------c

          nt(n)=int(sr(n)/delr2)+1

          ei(nt(n),  1,m)    = ei(nt(n),  1,m)-swt(n)*dei1
          ee(nt(n),  1,m)    = ee(nt(n),  1,m)-swt(n)*dee1
!
          idsgn = 3-isswt(n)
          ei(nt(n),idsgn ,m) = ei(nt(n),idsgn  ,m) -dei1
          ee(nt(n),idsgn ,m) = ee(nt(n),idsgn  ,m) -dee1

!-----------------------------------------------------c
!
          ekin0(n) = ekin0(n) +(e1 -e0(n))/ckb

          dee(m)=dee(m)-swt(n)*dee1
          dei(m)=dei(m)-swt(n)*dei1

          v1=sqrt(2.d0*e1/tmass)
          velh(n)=r1(n)*v1
          vels(n)=v1*sqrt(1.d0-r1(n)**2)
          es(n)=0.5d0*tmass*vels(n)**2
          yc(n,4)=velh(n)*tmass/bb(n)/tchg

          adinv(n)=es(n)/bb(n)
          es(n)=es(n)/ckb
!
        END DO
      END DO
!
      IF (it.eq.1) THEN

!CDIR NODEP
        DO n=1,ntest
          velr(n)=0.d0
        END DO
      ELSE

!CDIR NODEP
        DO n=1,ntest
          velr(n)=a*(sr(n)-sr1(n))/dtau
        END DO
      END IF
!
      DO n2=1,mm
!CDIR NODEP
        DO m=1,mm2
          n = n2+mm*(m-1)
!
!------------ evaluate ivvh, ivvs and iir  ----->>>>>>>>>>>>
          vvh  = vvmax+velh(n)
          ivvh = int(vvh/dvvmax)+1
!
          ivvs = int(vels(n)/dvvmax)+1
          iir  = int(sr(n)  /delr  )+1

!------------ evaluate ith and iph ----->>>>>>>>>>>>

          th  = dmod(yc(n,2),cpi2)
          IF (th.lt.0.0d0) THEN
            th  = cpi2+th
          END IF

          iit = int(th/delth)+1

          ph  = dmod(yc(n,3),cpi2)
          IF (ph.lt.0.0d0) THEN
            ph  = cpi2+ph
          END IF

          iip = int(ph/delph)+1
!

!------------ evaluate ixx and iyy ----->>>>>>>>>>>>

          xx = sr(n)*cos(th)+1.d0
          yy = sr(n)*sin(th)+1.d0

          iix = int(xx/delxx)+1
          iiy = int(yy/delyy)+1

!------------------------------------------------------------c
!     evaluate the weight

          cctm   = -(time/rswt(n))**8
!
          IF (cctm.lt.-20.d0) THEN
            wrate  =0.d0
          ELSE
            wrate  =exp(cctm)
          END IF
!
!
          ddpr= swt(n)
          dder= swt(n)*r1(n)**2
          ddvh= swt(n)*(1-r1(n)**2)
          ddvr= wrate*swt(n)*velr(n)*dtau
!
          ddw =         wrate *swt(n)*1.d8*dtau/ntest
          ddswt = (1.d0-wrate)*swt(n)*1.d8*dtau/ntest
!

!$$$  ddw = wrate*swt(n)*1.d8*dtau/ntest/vels(n)/sr(n)
!$$$  ddw = swt(n)*1.d8*dtau/ntest/vels(n)/sr(n)
!$$$  ddvh= swt(n)*dtau/ntest/sr(n)*velh(n)
!$$$  ddswt = (1.d0-wrate)*swt(n)*1.d8*dtau/ntest/sr(n)
!
!-------------------------------------------------------------c
!           for total sum
!-------------------------------------------------------------c
!
!            depvva(ivvh  ,ivvs  ,iir,m) = depvva(ivvh  ,ivvs  ,iir,m)
!     &           + ddw

          idsgn = isswt(n) -1
          iir2 = idsgn*mrmax+iir

          deppra(iir2,  m) = deppra(iir2,  m)+ddpr
          depera(iir2,  m) = depera(iir2,  m)+dder
          depvha(iir2  ,m) = depvha(iir2  ,m)+ddvh
          depvra(iir2  ,m) = depvra(iir2  ,m)+ddvr

!            deppta(iit,  m) = deppta(iit,  m)+ddpr
!            depeta(iit,  m) = depeta(iit,  m)+dder
!            depppa(iip,  m) = depppa(iip,  m)+ddpr
!            depepa(iip,  m) = depepa(iip,  m)+dder

          deppxa(iix,iiy,m) = deppxa(iix,iiy,m)+ddpr
          depexa(iix,iiy,m) = depexa(iix,iiy,m)+dder

          wmxt(iir  ,m)   = wmxt(iir  ,m)+ddswt

!$$$c------------------------------------------------------------c
!$$$c     for positive or negative weight particle sum
!$$$c-------------------------------------------------------------c
!$$$c
!$$$            idsgn = 2-isswt(n)
!$$$            iir2  = idsgn*mrmax+iir
!$$$            iit2  = idsgn*mtmax+iit
!$$$            iip2  = idsgn*mpmax+iip
!$$$            iiy2  = idsgn*mymax+iiy
!$$$c
!$$$            depvva(ivvh  ,ivvs  ,iir2,m) = depvva(ivvh  ,ivvs  ,iir2,m)
!$$$     &           + ddw
!$$$c
!$$$            depvha(iir2  ,m) = depvha(iir2  ,m)+ddvh
!$$$            depvra(iir2  ,m) = depvra(iir2  ,m)+ddvr
!$$$c
!$$$            deppra(iir2,  m) = deppra(iir2,  m)+ddpr
!$$$            depera(iir2,  m) = depera(iir2,  m)+dder
!$$$            deppta(iit2,  m) = deppta(iit2,  m)+ddpr
!$$$            depeta(iit2,  m) = depeta(iit2,  m)+dder
!$$$            depppa(iip2,  m) = depppa(iip2,  m)+ddpr
!$$$            depepa(iip2,  m) = depepa(iip2,  m)+dder
!$$$c
!$$$            deppxa(iix,iiy2,m) = deppxa(iix,iiy2,m)+ddpr
!$$$            depexa(iix,iiy2,m) = depexa(iix,iiy2,m)+dder
!$$$c
!$$$            wmxt(iir2 ,m)   = wmxt(iir2 ,m)+ddswt
!$$$c

        END DO
      END DO

!CDIR NODEP
      DO  n=1,ntest
        sr1(n)=sr(n)
      END DO
!
      RETURN
      END SUBROUTINE


!======================================================================c
      SUBROUTINE montcx(yc,adinv,eng,es,isw,dtau,dee,dei,ee,ei)
!======================================================================c

  USE mcnmod, ONLY : a, ane, anecgs, anem1, anem2, anew, cfe0, cfi0, ckb, clog, cme, &
    cmi, cpi2, denn0, depepa, depera, depeta, depexa, depppa, deppra, deppta, deppxa,&
    depvha, depvra, depvva, e0int, ekin0, engn0, imsh, iran, isswt, iswrd, it, ix,   &
    iy, mm, mm2, mpmax, mrmax, mtmax, mvmax, mxcx, mxmax, mymax, ntest, psia, rswt,  &
    sr1, swt, tchg, te, teev, tem1, tem2, tew, ti, tim1, tim2, time, tiw, tmass,     &
    ve, vi, wmxt
  IMPLICIT NONE
  INTEGER(4), INTENT(IN) :: isw(ntest)
  REAL(8),INTENT(IN)   :: dtau
  REAL(8),INTENT(INOUT):: adinv(ntest),dee(mm2), dei(mm2), yc(ntest,4),              &
    ee(imsh+1,3,mm2),ei(imsh+1,3,mm2)
  REAL(8),INTENT(OUT)  :: eng(ntest),es(ntest)
  INTEGER(4):: icon, idsgn, iip, iip2, iir, iir2, iit2, iix, iiy, iiy2, ircx, ivvh,  &
    ivvs, iit, m, n, n2, isf0(ntest),isf1(ntest),isf2(ntest), nt(ntest)
  REAL(8)::aloge, cctm, cfd, cfde, cfdi, cfee, cfei, dcfee, dcfei, dder, ddpr, ddswt,&
    ddvh, ddvr, ddw, dee1, dei1, delph, delr, delr2, delth, delxx, delyy, dvvmax, e1,&
    e1d, epp, etmax, ph, r1, ran, rt0e, rt0i, rt2e, rt2i, rtcx, rteng, sigcx, th, v1,&
    vvh, vvmax, wrate, xx, yy, bb(ntest), cfe1(ntest), cfi1(ntest), e0(ntest),       &
    fxe(ntest), fxe1(ntest), fxe2(ntest), fxi(ntest), fxi1(ntest), fxi2(ntest),      &
    gxe(ntest), gxi(ntest), r0(ntest),rt1e(ntest), rt1i(ntest), sr(ntest), velh(ntest),&
    velr(ntest), vels(ntest), velt(ntest), xe(ntest), xi(ntest)
  REAL(4):: rndcx(ntest)
  REAL(8):: espr1=1.0d-11, espr2=1.0d-10
      save
!      data espr1,espr2/1.0d-11,1.0d-10/
!!!      data espr1,espr2/1.0d-6,1.0d-5/


!     etmax : maxmum value of energy in grid
      etmax= e0int*2.d0
!     vvmax : maxmum value of velocity in grid
      vvmax= sqrt(2.0*etmax*ckb/tmass)
!     dvvmax: size of grid
      dvvmax=vvmax/dfloat(mvmax-1)

      delr  = 1.d0/dfloat(mrmax-1)
      delr2 = 1.d0/dfloat(imsh)
      delth = cpi2/dfloat(mtmax-1)
      delph = cpi2/dfloat(mpmax-1)
      delxx = 2.d0/dfloat(mxmax-1)
      delyy = 2.d0/dfloat(mymax-1)

      CALL iseifu(isf0)
      CALL iseifu(isf1)
      CALL iseifu(isf2)

      CALL bfld2n(yc,bb)

!CDIR NODEP
      DO n=1,ntest
        sr(n)=sqrt(yc(n,1)/psia)
        IF (sr(n).le.0.d0) THEN
          sr(n)=0.1d-16
        ELSE IF (sr(n).ge.1.d0) THEN
          sr(n) = 1.d0
        END IF
!
        rt0e=(1.d0-tew/te)*(1.d0-sr(n)**tem1)**tem2+tew/te
        rt0i=(1.d0-tiw/ti)*(1.d0-sr(n)**tim1)**tim2+tiw/ti
        ran =(1.d0-anew/ane)*(1.d0-sr(n)**anem1)**anem2+anew/ane
!
!     coulomb logarithm
!
        epp=50.d0*ckb/te
!
!     if(rt0e.lt.epp) then
!     clog=23.4-1.15*dlog10(anecgs*ran)+3.45*dlog10(teev*rt0e)
!     else
!     clog=25.3-1.15*dlog10(anecgs*ran)+2.3*dlog10(teev*rt0e)
!     end if
!
        clog=24.0d0-dlog(sqrt(anecgs*ran)/(teev*rt0e))
!
        rt1e(n)=sqrt(1.d0/rt0e)
        rt1i(n)=sqrt(1.d0/rt0i)
        rt2e=ran*clog/rt0e**1.5d0
        rt2i=ran*clog/rt0i**1.5d0
        cfe1(n)=cfe0*rt2e
        cfi1(n)=cfi0*rt2i
!
        velh(n)=yc(n,4)*tchg/tmass*bb(n)
        vels(n)=sqrt(adinv(n)*bb(n)*2.d0/tmass)
        velt(n)=sqrt(velh(n)**2+vels(n)**2)
        e0(n)=0.5d0*tmass*velt(n)**2
        r0(n)=velh(n)/velt(n)
        xe(n)=velt(n)*rt1e(n)/ve
        xi(n)=velt(n)*rt1i(n)/vi

!...for electrons
!         if (xi(n).lt.20.d0) then
!            xi(n)=20.d0
!         end if

      END DO
!
!     print *,'Start gosakn'
      CALL gosakn(xe,fxe,fxe1,fxe2,gxe)
      CALL gosakn(xi,fxi,fxi1,fxi2,gxi)
!     print *,'End   gosakn'
!
!     print *,'Start roop 200'

      DO  n2=1,mm
!CDIR NODEP
        DO m=1,mm2
          n = n2+mm*(m-1)

!----------------------------c
!     evaluate the weight
!----------------------------c

          cctm   = -(time/rswt(n))**8

          IF (cctm.lt.-20.d0) THEN
            wrate = 0.d0
          ELSE
            wrate = exp(cctm)
          END IF

          cfdi=  cfi1(n)*(fxi(n)-gxi(n))/xi(n)**3
          cfde=  cfe1(n)*(fxe(n)-gxe(n))/xe(n)**3
          cfd = cfdi+cfde
!
          cfee=2.d0*cfe1(n)*(tmass/cme)*gxe(n)/xe(n)
          cfei=2.d0*cfi1(n)*(tmass/cmi)*gxi(n)/xi(n)
          dcfee=-2.d0*cfe1(n)*(tmass/cme)/sqrt(te*e0(n))      &
                *(fxe2(n)+6.d0*gxe(n))/(4.d0*xe(n)**2)*rt1e(n)
          dcfei=-2.d0*cfi1(n)*(tmass/cmi)/sqrt(ti*e0(n))      &
                *(fxi2(n)+6.d0*gxi(n))/(4.d0*xi(n)**2)*rt1i(n)
!
!     print 665, n,cfdi,cfde,cfei,cfee
 665      FORMAT(i4,' cfd,cfe = ',4e12.3)
!
          IF(abs(abs(r0(n))-1.d0).lt.espr1) THEN
            r0(n)=r0(n)/abs(r0(n))*(1.d0-espr2)
          END IF
!
          r1=r0(n)*(1.d0-cfd*dtau) +isf0(n)*sqrt((1.d0-r0(n)**2)*cfd*dtau)
!
          IF(r1.ge.1.d0) THEN
            r1=  1.d0-(r1-int(r1))
          ELSE IF (r1.le.-1.d0) THEN
            r1= -1.d0-(r1+int(abs(r1)))
          END IF
!
!     print *,cfee,dtau,dcfee
          dee1=-2.d0*cfee*dtau*(e0(n)-(1.5d0+e0(n)*dcfee/cfee)*te*rt0e) &
               +2.d0*isf1(n)*sqrt(te*rt0e*e0(n)*cfee*dtau)
!     print *,'dee is evaluated.'
!     print *,cfei,dtau,dcfei
          dei1=-2.d0*cfei*dtau*(e0(n)-(1.5d0+e0(n)*dcfei/cfei)*ti*rt0i) &
               +2.d0*isf2(n)*sqrt(ti*rt0i*e0(n)*cfei*dtau)
!     print *,'dei is evaluated.'
!
!
          dei1=dei1*wrate
          dee1=dee1*wrate

          e1d = dee1+dei1
          e1  = e0(n)+e1d

!     e1=e0(n)+dee1
!     e1=e0(n)
!     print 666, n,e0(n),dee1,dei1,r0(n),dtau,cfd
 666      FORMAT(i4,' e0,dee1,dei1,r0,dtau,cfd= ',6e12.3)
!
          eng(n)  = e1/ckb
!
!$$$  if(eng(n).lt.emin) then
!$$$  dee1=0.d0
!$$$  dei1=0.d0
!$$$  e1=e0(n)
!$$$  eng(n)=e1/ckb
!$$$  end if
!
!-----------------------------------------------------c
!        Energy deposition profile (radail)

!                       ei(n,1) <---ion total f
!                       ei(n,2) <---ion f+
!                       ei(n,3) <---ion f-
!                       ee(n,1) <---ele. total f
!                       ee(n,2) <---ele. f+
!                       ee(n,3) <---ele. f-
!-----------------------------------------------------c

          nt(n)=int(sr(n)/delr2)+1

          ei(nt(n),  1,m)    = ei(nt(n),  1,m)-swt(n)*dei1
          ee(nt(n),  1,m)    = ee(nt(n),  1,m)-swt(n)*dee1
!
          idsgn = 3-isswt(n)
          ei(nt(n),idsgn ,m) = ei(nt(n),idsgn  ,m) -dei1
          ee(nt(n),idsgn ,m) = ee(nt(n),idsgn  ,m) -dee1

!-----------------------------------------------------c
!
          ekin0(n) = ekin0(n) +(e1 -e0(n))/ckb

          dee(m)=dee(m)-swt(n)*dee1
          dei(m)=dei(m)-swt(n)*dei1

          v1=sqrt(2.d0*e1/tmass)
          velh(n)=r1*v1
          vels(n)=v1*sqrt(1.d0-r1**2)
          es(n)=0.5d0*tmass*vels(n)**2
          yc(n,4)=velh(n)*tmass/bb(n)/tchg

          adinv(n)=es(n)/bb(n)
          es(n)=es(n)/ckb
!
        END DO
      END DO
!
      IF (it.eq.1) THEN

!CDIR NODEP
        DO n=1,ntest
          velr(n)=0.d0
        END DO
      ELSE

!CDIR NODEP
        DO n=1,ntest
          velr(n)=a*(sr(n)-sr1(n))/dtau
        END DO
      END IF
!
      DO n2=1,mm
!CDIR NODEP
        DO m=1,mm2
          n = n2+mm*(m-1)
!
!------------ evaluate ivvh, ivvs and iir  ----->>>>>>>>>>>>
          vvh  = vvmax+velh(n)
          ivvh = int(vvh/dvvmax)+1
!
          ivvs = int(vels(n)/dvvmax)+1
          iir  = int(sr(n)  /delr  )+1

!------------ evaluate ith and iph ----->>>>>>>>>>>>

          th  = dmod(yc(n,2),cpi2)
          IF (th.lt.0.0d0) THEN
            th  = cpi2+th
          END IF

          iit = int(th/delth)+1

          ph  = dmod(yc(n,3),cpi2)
          IF (ph.lt.0.0d0) THEN
            ph  = cpi2+ph
          END IF

          iip = int(ph/delph)+1
!

!------------ evaluate ixx and iyy ----->>>>>>>>>>>>

          xx = sr(n)*cos(th)+1.d0
          yy = sr(n)*sin(th)+1.d0

          iix = int(xx/delxx)+1
          iiy = int(yy/delyy)+1

!------------------------------------------------------------c
!     evaluate the weight

          cctm   = -(time/rswt(n))**8
!
          IF (cctm.lt.-20.d0) THEN
            wrate=0.d0
          ELSE
            wrate  =exp(cctm)
          END IF
!
          ddvh= wrate*swt(n)*velh(n)*dtau
          ddvr= wrate*swt(n)*velr(n)*dtau
!
          ddpr= swt(n)*dtau
          dder= swt(n)*eng(n)*dtau
!
          ddw =         wrate *swt(n)*1.d8*dtau/ntest
          ddswt = (1.d0-wrate)*swt(n)*1.d8*dtau/ntest
!

!$$$  ddw = wrate*swt(n)*1.d8*dtau/ntest/vels(n)/sr(n)
!$$$  ddw = swt(n)*1.d8*dtau/ntest/vels(n)/sr(n)
!$$$  ddvh= swt(n)*dtau/ntest/sr(n)*velh(n)
!$$$  ddswt = (1.d0-wrate)*swt(n)*1.d8*dtau/ntest/sr(n)
!
!-------------------------------------------------------------c
!           for total sum
!-------------------------------------------------------------c
!
          depvva(ivvh  ,ivvs  ,iir,m) = depvva(ivvh  ,ivvs  ,iir,m) + ddw

          depvha(iir  ,m) = depvha(iir  ,m)+ddvh
          depvra(iir  ,m) = depvra(iir  ,m)+ddvr

          deppra(iir,  m) = deppra(iir,  m)+ddpr
          depera(iir,  m) = depera(iir,  m)+dder
          deppta(iit,  m) = deppta(iit,  m)+ddpr
          depeta(iit,  m) = depeta(iit,  m)+dder
          depppa(iip,  m) = depppa(iip,  m)+ddpr
          depepa(iip,  m) = depepa(iip,  m)+dder

          deppxa(iix,iiy,m) = deppxa(iix,iiy,m)+ddpr
          depexa(iix,iiy,m) = depexa(iix,iiy,m)+dder

          wmxt(iir  ,m)   = wmxt(iir  ,m)+ddswt

!------------------------------------------------------------c
!     for positive or negative weight particle sum
!-------------------------------------------------------------c

          idsgn = 2-isswt(n)
          iir2  = idsgn*mrmax+iir
          iit2  = idsgn*mtmax+iit
          iip2  = idsgn*mpmax+iip
          iiy2  = idsgn*mymax+iiy
!
          depvva(ivvh  ,ivvs  ,iir2,m) = depvva(ivvh  ,ivvs  ,iir2,m) + ddw

          depvha(iir2  ,m) = depvha(iir2  ,m)+ddvh
          depvra(iir2  ,m) = depvra(iir2  ,m)+ddvr

          deppra(iir2,  m) = deppra(iir2,  m)+ddpr
          depera(iir2,  m) = depera(iir2,  m)+dder
          deppta(iit2,  m) = deppta(iit2,  m)+ddpr
          depeta(iit2,  m) = depeta(iit2,  m)+dder
          depppa(iip2,  m) = depppa(iip2,  m)+ddpr
          depepa(iip2,  m) = depepa(iip2,  m)+dder
!
          deppxa(iix,iiy2,m) = deppxa(iix,iiy2,m)+ddpr
          depexa(iix,iiy2,m) = depexa(iix,iiy2,m)+dder

          wmxt(iir2 ,m)   = wmxt(iir2 ,m)+ddswt


        END DO
      END DO

!CDIR NODEP
      DO  n=1,ntest
        sr1(n)=sr(n)
      END DO

!-------------------------------------c
!     charge exchange with neutral
!-------------------------------------c
!
      CALL rjuflp(ntest,ix,iy,iran,rndcx,iswrd,icon)

!CDIR NODEP
      DO n=1,ntest

        IF (eng(n).ge.1.0d3) THEN
          aloge=dlog10(eng(n))
          sigcx=0.6937d-14*(1.d0-0.155d0*aloge)**2/(1.d0+0.1112d-14*eng(n)**3.3d0)
!
          ircx = int(sr(n)*dfloat(mxcx))+1
          rtcx = dtau*velt(n)*1.d2*sigcx*denn0(ircx)
!
          IF (rndcx(n).lt.rtcx) THEN
            rteng   = engn0(ircx)/eng(n)
            yc(n,4) =  yc(n,4)*sqrt(rteng)
            eng(n)  =   eng(n)*rteng
            es(n)   =    es(n)*rteng
            adinv(n)= adinv(n)*rteng
            ekin0(n)= ekin0(n)*rteng
          END IF
        END IF
!
      END DO
!
      RETURN
      END SUBROUTINE


!======================================================================c
      SUBROUTINE fmxwll
!======================================================================c

  USE mcnmod, ONLY : ckb, e0int, fmx0, mrmax, mvmax, mvmax2, te, tem1, tem2, tew, &
    tmass
  IMPLICIT NONE
  integer(4):: ir, ivpl, ivpp
  REAL(8):: dvvmax, etmax, ffsum, rte, sr, vpl, vpp, vthe, vvmax, xe, &
    fmxa(mvmax2,mvmax)
      save

!     etmax : maxmum value of energy in grid
      etmax= e0int*2.d0
!     vvmax : maxmum value of velocity in grid
      vvmax= sqrt(2.0d0*etmax*ckb/tmass)
!     dvvmax: size of grid
      dvvmax=vvmax/dfloat(mvmax-1)

      DO ir = 1, mrmax
!ppp
        PRINT *, 'ir =',ir

        IF (ir.eq.1) THEN
          sr = 1.d-4
        ELSE
          sr = dfloat(ir-1)/dfloat(mrmax-1)
        END IF

        rte=(te -tew )*(1.d0-sr**tem1 )**tem2 +tew
        vthe=sqrt(2.d0*rte/tmass)

!ppp
        PRINT *, 'sr,rte,vthe=',sr,rte,vthe

        ffsum = 0.d0
        DO ivpl=1,mvmax2
          vpl=dfloat(ivpl-32)*dvvmax
          DO ivpp=1,mvmax
            vpp  = dfloat(ivpp-1)*dvvmax
            xe  = (vpl**2+vpp**2)/vthe**2
            fmxa(ivpl,ivpp) = vpp*exp(-xe)
            ffsum= ffsum +fmxa(ivpl,ivpp)
          END DO
        END DO

!ppp
        PRINT *, 'ffsum=',ffsum
!     ffmxa is normalized by ffsum.
        IF (ffsum.gt.0.d0) THEN
!
          DO ivpl=1,mvmax2
            DO ivpp=1,mvmax
              fmx0(ivpl,ivpp,ir)= fmxa(ivpl,ivpp)/ffsum
            END DO
          END DO
        END IF
!
      END DO

      RETURN
      END SUBROUTINE

!======================================================================c
      SUBROUTINE trans(y,rg,zg,pg)
!======================================================================c

  USE mcnmod, ONLY : cpi
  IMPLICIT NONE
  REAL(8), INTENT(IN) :: y(4)
  REAL(8), INTENT(OUT):: rg, zg, pg
  REAL(8):: pol, tol

      pol = mod(y(2),2.*cpi)
      tol = mod(y(2),2.*cpi)
      rg  = y(1)*cos(pol)
      zg  = y(1)*sin(pol)
      pg  = tol

      RETURN
      END SUBROUTINE

!======================================================================c
      SUBROUTINE iseifu(isf)
!======================================================================c

  USE mcnmod, ONLY : iran, iswrd, ix, iy, ntest
  IMPLICIT NONE
  INTEGER(4),INTENT(OUT):: isf(ntest)
  INTEGER(4):: i, icon
  REAL(4):: ran(ntest)

!      call ranu2n(ran,icon)
!      call djuflp(ntest,ix,iy,iran,ran,iswrd,icon)
      CALL rjuflp(ntest,ix,iy,iran,ran,iswrd,icon)
!      call rjufsp(ntest,ix,iy,ran,icon)

!     if(icon.ne.0) write(6,100) icon

!  100 format(1h ,'********** error in sub. iseifu ***********')

!CDIR NODEP
      DO i=1,ntest
        isf(i)=int(2.d0*ran(i))*2-1
      END DO

      RETURN
      END SUBROUTINE

!======================================================================c
      SUBROUTINE gosakn(x,fx,fx1,fx2,gx)
!======================================================================c

  USE mcnmod, ONLY : ntest
  IMPLICIT NONE
  EXTERNAL derf
  REAL(8),INTENT(IN) :: x(ntest)
  REAL(8),INTENT(OUT):: fx(ntest), fx1(ntest), fx2(ntest), gx(ntest)
  INTEGER(4):: ip
  REAL(8):: f(ntest),derf
  REAL(8):: cpi=3.141592653589793d0
!      DATA RH / 0.70710678118654752440D+00/
!      data cpi/3.141592653589793d0/

!CDIR NODEP
      DO ip=1,ntest

        f(ip)= DERF(X(ip))*0.5D0

        fx(ip)=2.d0*f(ip)
        fx1(ip)=2.d0/sqrt(cpi)*exp(-x(ip)**2)
        fx2(ip)=-2.d0*x(ip)*fx1(ip)
        gx(ip)=0.d0

        IF(abs(x(ip)).gt.1.e-10) THEN
          gx(ip)=(fx(ip)-x(ip)*fx1(ip))/(2.d0*x(ip)**2)
        END IF

      END DO

      RETURN
      END SUBROUTINE

!======================================================================c
      SUBROUTINE diagi
!======================================================================c
!     computation of out put data,
!     such as energy deposition
!     and velocity distribution of minority ions.
!----------------------------------------------------------------------c

  USE mcnmod , ONLY:  adinv, aveng, avengs,  ckb, dep, depehs, depihs, depths, dt,  &
    dtime, e0int, ekin0, eloss, eng, enghsm, enghsp, enghst, es, etlos, f, ffm, ffp,&
    iend, iit, isw, it, loss, mm2, nloss, npt, ntest, ntot, prf, rswt, swt, tchg,   &
    tdep, tdepe, tdepi, teev, time, tmass, tmhst, vhhst, vhhstm, vhhstp, vx, vy, y, yc
  IMPLICIT NONE
  INTEGER(4):: i, iernum, ip, kf, n, nermax, nerptl, nprtcl, ntotm, ntotp, ierdis(20)
  REAL(8):: cctm, dlter, dw, eeror, engerr, ereng, ermax, errtot, estot, etot, etotm,&
    etotp, feng, pg, rg, stdep, stdepe, stdepi, teng, vh, vhtot, vhtotm, vhtotp, vs, &
    wrate, zg, b


      iit=iit+1

      nloss(iit) = loss
      etlos(iit) = eloss

      ntot   = 0
      ntotp  = 0
      ntotm  = 0
      nerptl = 0
      nprtcl = 0

      etot   = 0.d0
      estot  = 0.d0
      etotp  = 0.d0
      etotm  = 0.d0
      vhtot  = 0.d0
      vhtotp = 0.d0
      vhtotm = 0.d0
      ermax  = 0.d0
      teng   = 0.d0
      stdep  = 0.d0
      stdepe = 0.d0
      stdepi = 0.d0

      DO n=1,ntest

        IF (isw(n).eq.1) THEN
          y(1) = yc(n,1)
          y(2) = yc(n,2)
          y(3) = yc(n,3)
          y(4) = yc(n,4)
          CALL bfld2(y,b)
          vh = yc(n,4)*tchg*b/tmass
          vs = sqrt(adinv(n)*2.d0*b/tmass)
!
          nprtcl = nprtcl+1
          vx(n)  = vh
          vy(n)  = vs
          es(n)  = 0.5d0*tmass*vs**2/ckb
          eng(n) = es(n)+0.5d0*tmass*vh**2/ckb
!
          ntot   = ntot +1
!
          nerptl = nerptl+1
          engerr = abs((eng(n) -ekin0(n))/ekin0(n))
!
          IF (engerr.gt.ermax) THEN
            nermax = n
            ermax  = engerr
          END IF
!
          errtot = errtot+engerr

!
          cctm   = -(time/rswt(n))**8
!
          IF (cctm.lt.-20.d0) THEN
            wrate=0.d0
          ELSE
            wrate  =exp(cctm)
          END IF
!     wrate=1.d0
!
          IF (swt(n).gt.0.d0) THEN
            ntotp   = ntotp +1
            etotp   = etotp  +eng(n)*wrate +1.5d0*teev*(1.d0-wrate)
            vhtotp  = vhtotp +vh*wrate
          ELSE
            ntotm   = ntotm +1
            etotm   = etotm  +eng(n)*wrate +1.5d0*teev*(1.d0-wrate)
            vhtotm = vhtotm +vh*wrate
          END IF
!
!     etot   = etot  +swt(n)*eng(n)
          etot   = etot  +eng(n)*wrate +1.5d0*teev*(1.d0-wrate)
          estot  = estot +es(n)
          vhtot  = vhtot +vh*wrate*swt(n)
!
        ENDIF
!
        IF (isw(n).ge.0) THEN
          teng=teng+eng(n)
        ENDIF
!
      END DO

      DO i=1,mm2
        stdep  = stdep  +tdep(i)
        stdepe = stdepe +tdepe(i)
        stdepi = stdepi +tdepi(i)
      END DO

      prf(iit)    = teng +dep(iit)
      aveng(iit)  = etot/float(ntot)
      avengs(iit) = estot/float(ntot)
      npt(iit)    = ntot
      dtime(iit)  = dt*float(it)
      tmhst(iit)  = dtime(iit)
      vhhst(iit)  = vhtot/float(ntot)
      vhhstp(iit)  = vhtotp/float(ntot)
      vhhstm(iit)  = vhtotm/float(ntot)
      enghst(iit)  = etot/float(ntot)
      enghsp(iit)  = etotp/float(ntotp)

      IF (ntotm.gt.0) THEN
        enghsm(iit)  = etotm/float(ntotm)
      ELSE
        enghsm(iit)  = 0.d0
      END IF

      depths(iit)  = stdep/ckb
      depehs(iit)  = stdepe/ckb
      depihs(iit)  = stdepi/ckb

      DO i=1,iend+1
        f(i)  = 1.d0
        ffp(i) = 1.d0
        ffm(i)  = 1.d0
      END DO

!       dw=tiev/10.0
      dw=e0int*1.5d0/dfloat(iend)

      DO n=1,ntest
        IF (isw(n).eq.1) THEN

!            jf=int(es(n)/dw)+1
!            if (jf.gt.iend) jf=iend+1
!            fs(jf)=fs(jf)+1.d0

          kf=int(eng(n)/dw)+1
          IF(kf.gt.iend) kf=iend+1
          f(kf)=f(kf)+1.d0
!
          IF (swt(n).gt.0.d0) THEN
            ffp(kf)=ffp(kf)+1.d0
          ELSE
            ffm(kf)=ffm(kf)+1.d0
          END IF

        ENDIF
      END DO
!
!
      WRITE(6,665) time
      WRITE(6,667)
      DO i=1,iend+1
        f(i)  =  dlog(f(i))
        ffp(i)=  dlog(ffp(i))
        ffm(i)=  dlog(ffm(i))
        feng  = (dfloat(i-1)+0.5d0)*dw
        WRITE(6,676) feng, f(i), ffp(i),ffm(i)
      END DO

 665  FORMAT(/,5x,'time=',e10.3,/,15x,'alog(f)',/)
 667  FORMAT(7x,'eng[keV]',12x,'alog(ff)',12x, 'alog(ffp)',12x,'alog(ffm)' )
 676  FORMAT(4e12.4)
!
      eeror = errtot/ntot

      DO i=1,4
        y(i)=yc(nermax,i)
      END DO

      CALL trans(y,rg,zg,pg)

      WRITE(6,701) time
      WRITE(6,702) eeror
      WRITE(6,703) nerptl
      WRITE(6,704) ermax
      WRITE(6,705) eng(nermax)
      WRITE(6,706) rg
      WRITE(6,707) zg
      WRITE(6,708) pg
      WRITE(6,709) nermax

  701  FORMAT(/,/,' ****** particle energy error check ****** ' &
     &                ,/,'  time = ',e10.3)
  702 FORMAT(2x,' error averaged   = ',e10.3)
  703 FORMAT(2x,' count particle   = ',i6)
  704 FORMAT(2x,' maximum error    = ',e10.3)
  705 FORMAT(2x,' energy(maxerro)  = ',e10.3)
  706 FORMAT(2x,'     rg(maxerro)  = ',e10.3)
  707 FORMAT(2x,'     zg(maxerro)  = ',e10.3)
  708 FORMAT(2x,'     pg(maxerro)  = ',e10.3)
  709 FORMAT(2x,'        maxerro   = ',i5)

      dlter = ermax /20.d0

      DO i=1, 20
        ierdis(i) = 0
      END DO


      DO ip=1, ntest
        IF (isw(ip).eq.1) THEN
          ereng  = abs(eng(ip) -ekin0(ip))/ekin0(ip)
          iernum = int(ereng/dlter)+1
          ierdis(iernum) = ierdis(iernum) +1
        END IF
      END DO

!      add max error particle
      ierdis(20) = ierdis(20) +1

      WRITE(6,711)
      WRITE(6,712)
      WRITE(6,713) ierdis
      WRITE(6,714)

  711 FORMAT(/,1x,'error distribution')
  712 FORMAT(5x,'0----+----+----+----+----+----+----+----+----', &
     & '+----+----+----+----+----+----+----+----+----+----+----ermax')
  713 FORMAT(5x,20i5)
  714 FORMAT(5x,'0----+----+----+----+----+----+----+----+----', &
     & '+----+----+----+----+----+----+----+----+----+----+----ermax',/)

      RETURN
      END SUBROUTINE

END MODULE
