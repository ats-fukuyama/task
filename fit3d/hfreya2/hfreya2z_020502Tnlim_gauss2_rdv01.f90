MODULE hfreya_all



  CONTAINS
!***********************************************************************
!       program      h f r e y a  2
!                   last modified by S.Murakami (2007.8.22)
!***********************************************************************
!                   2007.12 To Fortran90 subroutine

      SUBROUTINE hfreya

  USE hfrmod, ONLY : icoef
  USE hfrmod, ONLY : allocate_a, deallocate_a, allocate_b, deallocate_b, &
                     allocate_through, deallocate_through
  IMPLICIT NONE

      CALL allocate_through
      CALL allocate_a

!-----------------------------------c
!     set the beam parameters
!-----------------------------------c

      CALL nbiinp

!-----------------------------------c
!     read the field data
!-----------------------------------c

      CALL dkread

!-----------------------------------c
!     normalize the data
!-----------------------------------c

      CALL dimles
      PRINT *, 'dimles is ended.'
      IF(icoef.eq.0) CALL coeff
      IF(icoef.eq.1) CALL coeff1

!-----------------------------------c
!     beam initialization
!-----------------------------------c

      CALL binit
!
      CALL deallocate_a
      CALL allocate_b

!-----------------------------------c
!     calculate the NBI deposition
!-----------------------------------c

      CALL nfreya
!
!-----------------------------------c
!     write the resuls
!-----------------------------------c

      CALL fryout

!-----------------------------------c

      CALL deallocate_b
      CALL deallocate_through

! 100  STOP
      RETURN
      END SUBROUTINE


!#####################################################################
      SUBROUTINE rzf(ps,theta,phi,r,z)
!#####################################################################
!
!     iterface for rzpfun
!
!--------------------------------------------------------------------c

  IMPLICIT NONE
  REAL(8),INTENT(IN) :: phi, ps, theta
  REAL(8),INTENT(OUT):: r, z
  REAL(8):: pp, rr, zz, x(3)

      x(1)=ps
      x(2)=theta
      x(3)=phi

      CALL rzpfun(x,rr,zz,pp)
      r=rr
      z=zz
      RETURN
      END SUBROUTINE


!#####################################################################
      SUBROUTINE dkread
!#####################################################################

!     Read the field data by calling iodisk.

!--------------------------------------------------------------------c

  USE hfrmod, ONLY : a0, ai, aiota, b2inm, bnm, gg, loop, mm, mn0, mnumbr, modmax,&
    nn, nnumbr, pnm, psi, r0, rnm, vol, vprime, znm
  IMPLICIT NONE
  INTEGER(4):: i, k, mx, nhel, mnum(mm), nnum(mm)
  REAL(8):: b00, bbnm(mm,0:nn), ppnm(mm,0:nn), psip(0:nn), rrnm(mm,0:nn), zznm(mm,0:nn),  &
     zbm(mm)

!----------------------------c
!     Read the field data
!----------------------------c

      CALL iodisk_hfr(modmax,mn0,nhel,loop,nnum,mnum,psi,gg,ai,   &
                      aiota,vol,bbnm,b2inm,rrnm,zznm,ppnm)

!-------------------------------------------c
!     Calculate the values at the axis
!-------------------------------------------c

!     b00=fcent(bbnm,psi,mn0,m)
!     r0 =fcent(rrnm,psi,mn0,m)

      b00 = bbnm(1,0)
      r0 = rrnm(1,0)
      a0 = sqrt(2.d0*psi(loop)/b00)


      DO k=1,modmax
        zbm(k)=0.0
        DO i=0,loop
          bbnm(k,i)=bbnm(k,i)/b00
          zbm(k)   =dmax1(zbm(k),dabs(bbnm(k,i)))
        END DO
      END DO

      DO i=0,loop
        ai(i)  = ai(i)/b00
        gg(i)  = gg(i)/b00
        psi(i) = psi(i)/b00
        psip(i)= psip(i)/b00
      END DO

      psi(0)  =0.0d0
      psip(0) =0.0d0
      ai(0)   =0.0d0

!     aiota(0)=fcent0(aiota,psi)
!     gg(0)   =fcent0(gg,psi)


!--------------------------------------------c
!    Print out the field data
!--------------------------------------------c

      PRINT 600,  modmax,loop
      PRINT 601, (nnum(k),mnum(k),k=1,modmax)
!
      PRINT 602, (i,psi(i),psip(i),gg(i),ai(i),aiota(i),vol(i),vprime(i),i=0,loop)
!
      PRINT 610, (nnum(k),mnum(k),k=1,11)
!
      DO i=0,loop
        PRINT 611, i,(bbnm(k,i),k=1,11)
      END DO
!
      PRINT 610, (nnum(k),mnum(k),k=12,22)

      DO i=0,loop
        PRINT 611, i,(bbnm(k,i),k=12,22)
      END DO


 600  FORMAT(1h ,'  *** modmax loop ***',2i5)
 601  FORMAT(1h //,'  *** n ,m ****',//(2i10))
 602  FORMAT(1h //,4x,'i',9x,'psi',8x,'psip',10x,'gg',10x,'ai', &
     &     7x,'aiota',9x,'vol',7x,'vprime'//(i5,1p7e12.3))
 610  FORMAT(1h ,//,' ****** fourier amplitude ** bnm **'//, &
     &     3x,11(2x,'(',i3,i5,')'),//)
 611  FORMAT(i3,1p11e12.2)
!
!--------------------------------------------c

      DO k=1,modmax
        mnumbr(k)=mnum(k)
        nnumbr(k)=nnum(k)
        DO i=0,loop
          bnm(k,i)=bbnm(k,i)
          rnm(k,i)=rrnm(k,i)
          znm(k,i)=zznm(k,i)
          pnm(k,i)=ppnm(k,i)
        END DO
      END DO
      mn0=1


!--------------------------------------------c
!    Print out the field data
!--------------------------------------------c

      mx=min(11,modmax)
      PRINT 620, (nnumbr(k),mnumbr(k),k=1,mx)

      DO i=0,loop
        PRINT 611, i,(bnm(k,i),k=1,mx)
      END DO

 620  FORMAT(1h ,//,' ****** fourier amplitude ** bnm **'//,3x,11(2x,'(',i3,i5,')'),//)

!--------------------------------------------c

      RETURN
      END SUBROUTINE


!#####################################################################
      SUBROUTINE findps(r,phi,z,ps,thet,sphi)
!#####################################################################
!
!     Find the position in the Boozer coordinates
!     from the position in the cylinder coordinates.

!--------------------------------------------------------------------c

  USE hfrmod, ONLY : r0, xmin

  IMPLICIT NONE
  REAL(8),INTENT(IN) :: r, phi, z
  REAL(8),INTENT(OUT):: ps,thet,sphi
  INTEGER(4):: i, j, imx, jmx, ix, iz
  INTEGER(4):: ifst=0
  REAL(8):: rr, sr, sz, zl, zmin, zpsi, zt, zz, x0(3), xxx(21,64),zzz(21,64)
!
!  Set the initial guess of phi
      sphi=phi

!  Normalization of  R and z
      sr=r/r0
      sz=z/r0


      imx=20
      jmx=64
      xmin=0.01/dfloat(imx)

!  Initialize the table

      IF(ifst.eq.0) THEN
        DO i=1,imx+1
          DO j=1,jmx
            zt=dfloat(j-1)/dfloat(jmx)*3.14159265*2.
!x             zt=atan2(z,r-r0)
            zpsi=dfloat(i)/dfloat(imx)
            CALL rzf(zpsi,zt,phi,rr,zz)
            xxx(i,j)=rr
            zzz(i,j)=zz
          END DO
        END DO
        ifst=1
      ENDIF

!  Obtain the nearest grid point
!
      zmin=1.d70
      DO 200 i=1,imx+1
        DO 200 j=1,jmx
          zl=(xxx(i,j)-sr)**2+(zzz(i,j)-sz)**2
          IF(zl.gt.zmin) GO TO 200
          zmin=zl
          ix=i
          iz=j
 200      CONTINUE
      thet=dfloat(iz-1)/dfloat(jmx)*3.14159265*2.
!x    thet=zt
      IF(ix.le.imx) THEN
        ps=dfloat(ix)/dfloat(imx)
      ELSE
         ps=1.5
      ENDIF
!x    write(6,600) ix,iz,ps,zmin
 600  FORMAT(' *** ix,iz,ps,zmin***',2i5,1p2e14.4)

!
!x    if(ps.gt.1.0) return
!
!x    if(ic05pc.eq.0) return
!
      x0(1)=ps
      x0(2)=thet
      x0(3)=phi

!     Find ps,thet,sphi
!     with the initial guess x0(n).
!
      CALL findpa(x0,sr,phi,sz,ps,thet,sphi)
!
      RETURN
      END SUBROUTINE


!#####################################################################
      SUBROUTINE findpa(x0,r,phi,zf,ps,thet,sphi)
!#####################################################################

  USE hfrmod, ONLY: pphi, rr, zz

  IMPLICIT NONE
  REAL(8),INTENT(IN) :: x0(*), r, phi, zf
  REAL(8),INTENT(OUT):: ps, thet, sphi
  INTEGER(4),PARAMETER:: n=3, ldfjac=n, lr=(n*(n+1))/2
  INTEGER(4):: ifail, iwk(n), nev(2)
  INTEGER(4):: maxfev=100, mode=0
  REAL(8):: er, x(n), fvec(n), wk(n*(n+3)), dwk(n*(2*n+2))
  REAL(8):: fjac(ldfjac,n),diag(n),qtf(n),rq(lr),w(n,4)
  REAL(8):: xtol,factor
  INTEGER(4):: nprint,nfev,njev

!!  n is the number of nonlinear algebrac equations to be solved by the
!!!     data x0/0.6,0.9,1.2,1.5,1.8,0.,0.,0.,0./

!C---NEC 01/02/19 start
!      real wk(n*(n+3)),dwk(n*(2*n+2))
!      integer iwk(n),nev(2)
!*     real wk(n*(3*n+7))
!*     integer iwk(2*n)
!       external fcn
!      external fcnv,fcnj
!C---NEC END
!
!x    write(6,601) rr,zz,pphi,r0
 601  FORMAT(' *** rr,zz,pphi,r0**',1p4e12.3)
      x(1)=x0(1)
      x(2)=x0(2)
      x(3)=x0(3)
      rr=r
      zz=zf
      pphi=phi

      nprint=1
      factor=1.d0
      xtol=0.d0

!C---NEC 01/02/19 start
      call c05pcf(fcn,n,x,fvec,fjac,ldfjac,xtol,maxfev,diag,mode, &
     &     factor,nprint,nfev,njev,rq,lr,qtf,w,ifail)
!      nev(1) = maxfev
!      nev(2) = 100*n
!      er     = 0.1d-10
!      CALL DLSRDS(fcnv,fcnj,x,n,er,nev,mode,iwk,wk,dwk,ifail)
!C---NEC END
!x    if(x(1).le.xmin) then
!x    x(1)=xmin
!x    x(2)=x0(2)
!x    x(3)=x0(3)
!x    endif
      IF(ifail.ne.0) THEN
        WRITE(7,6000) ifail,fvec(1),fvec(2),fvec(3),x(1)
      ENDIF
 6000 FORMAT(' ***ifail**fvec**psi*',i5,1p4e13.4)
      ps=x(1)
      thet=x(2)
      sphi=x(3)

      RETURN
      END SUBROUTINE


!#####################################################################
      SUBROUTINE fcn(n,x,fvec,fjac,ldfjac,iflag)
!#####################################################################

  USE hfrmod, ONLY: abmnum, cp, cr, cz, dp, dr, dz, ep, er, ez, loop, mnumbr,    &
    modmax, nnumbr, pnmn, pphi, psino, rnmn, rr, znmn, zz
  IMPLICIT NONE

  INTEGER(4),INTENT(IN) :: n, ldfjac, iflag
  REAL(8),INTENT(INOUT) :: x(n)
  REAL(8),INTENT(OUT):: fvec(n), fjac(ldfjac,n)
  INTEGER(4):: i, ii, j, k
  REAL(8):: cosf, dpharm, drharm, dxabm, dzharm, p, pharm, r, rharm, sinf, x1, x2,  &
    x3, xabm, xx, z, zharm

!x    if(x(1).le.xmin) then
!x    iflag=-1
!x    return
!x    endif

      IF(x(1).lt.0.0) THEN
        x(1)=-x(1)
        x(2)=x(2)+3.1415926535897
      ENDIF
!
      x1=x(1)
      x2=x(2)
      x3=x(3)

      IF(iflag.eq.2) GO TO 200
!
      CALL rzpfun(x,r,z,p)

      fvec(1)=r-rr
      fvec(2)=z-zz
      fvec(3)=p-pphi

      RETURN

 200  CONTINUE
      IF(x1.gt.1.0) GO TO 500
      IF(x1.eq.0.0) THEN
        ii=0
        xx=0.
        GO TO 50
      ENDIF

      DO i=1,loop
        IF(x1.gt.psino(i-1).and.x1.le.psino(i)) THEN
          ii=i-1
          xx=x1-psino(i-1)
          GO TO 50
        ENDIF
      END DO

 50   CONTINUE

      DO i=1,ldfjac
        DO j=1,n
          fjac(i,j)=0.
        END DO
      END DO

      DO k=1,modmax
        xabm=x1**abmnum(k)
        dxabm=x1**(abmnum(k)-1)*abmnum(k)
        cosf=dcos(mnumbr(k)*x2-nnumbr(k)*x3)
        sinf=dsin(mnumbr(k)*x2-nnumbr(k)*x3)
        rharm=rnmn(k,ii)+(cr(k,ii)+(dr(k,ii)+er(k,ii)*xx)*xx)*xx
        zharm=znmn(k,ii)+(cz(k,ii)+(dz(k,ii)+ez(k,ii)*xx)*xx)*xx
        pharm=pnmn(k,ii)+(cp(k,ii)+(dp(k,ii)+ep(k,ii)*xx)*xx)*xx
        drharm=cr(k,ii)+(2.d0*dr(k,ii)+3.d0*er(k,ii)*xx)*xx
        dzharm=cz(k,ii)+(2.d0*dz(k,ii)+3.d0*ez(k,ii)*xx)*xx
        dpharm=cp(k,ii)+(2.d0*dp(k,ii)+3.d0*ep(k,ii)*xx)*xx
        fjac(1,1)=fjac(1,1)+(drharm*xabm+rharm*dxabm)*cosf
        fjac(2,1)=fjac(2,1)+(dzharm*xabm+zharm*dxabm)*sinf
        fjac(3,1)=fjac(3,1)+(dpharm*xabm+pharm*dxabm)*sinf
        fjac(1,2)=fjac(1,2)-mnumbr(k)*rharm*xabm*sinf
        fjac(1,3)=fjac(1,3)+nnumbr(k)*rharm*xabm*sinf
        fjac(2,2)=fjac(2,2)+mnumbr(k)*zharm*xabm*cosf
        fjac(2,3)=fjac(2,3)-nnumbr(k)*zharm*xabm*cosf
        fjac(3,2)=fjac(3,2)+mnumbr(k)*pharm*xabm*cosf
        fjac(3,3)=fjac(3,3)-nnumbr(k)*pharm*xabm*cosf
      END DO
      fjac(3,3)=1.+fjac(3,3)
      RETURN

 500  CONTINUE
!x    write(6,600) x(1),x(2),x(3)
 600  FORMAT(1h //,' **** x **** out of range****',1p3e14.4)
 800  CONTINUE
      DO i=1,ldfjac
        DO j=1,n
          fjac(i,j)=0.
        END DO
      END DO
!
      DO k=1,modmax
        cosf=dcos(mnumbr(k)*x2-nnumbr(k)*x3)
        sinf=dsin(mnumbr(k)*x2-nnumbr(k)*x3)
        IF(k.eq.1) GO TO 910
        fjac(1,1)=fjac(1,1)+rnmn(k,loop)*cosf
        fjac(1,2)=fjac(1,2)-mnumbr(k)*rnmn(k,loop)*sinf*x1
        fjac(1,3)=fjac(1,3)+nnumbr(k)*rnmn(k,loop)*sinf*x1
 910    fjac(2,1)=fjac(2,1)+znmn(k,loop)*sinf
        fjac(2,2)=fjac(2,2)+mnumbr(k)*znmn(k,loop)*cosf*x1
        fjac(2,3)=fjac(2,3)-nnumbr(k)*znmn(k,loop)*cosf*x1
        fjac(3,2)=fjac(3,2)+mnumbr(k)*pnmn(k,loop)*cosf
        fjac(3,3)=fjac(3,3)-nnumbr(k)*pnmn(k,loop)*cosf
      END DO

      fjac(3,3)=fjac(3,3)+1.0
      RETURN
      END SUBROUTINE


!#####################################################################
      SUBROUTINE rzpfun(x,rr,zz,pp)

!     Find the position in the cylinder coordinates
!     form the position in the Boozer coordinates.

!--------------------------------------------------------------------c

  USE hfrmod, ONLY : abmnum, cp, cr, cz, dp, dr, dz, ep, er, ez, loop, mnumbr,    &
    modmax, nnumbr, pnmn, psino, rnmn, znmn
  IMPLICIT NONE
  REAL(8),INTENT(IN) :: x(3)
  REAL(8),INTENT(OUT):: rr, zz, pp
  INTEGER(4):: i, ii, k
  REAL(8):: cosf, sinf, xabm, xx, xx2, xx3

      IF(x(1).lt.0.d0) THEN
!!         iflag=0                             200710
      ELSE IF(x(1).eq.0.0d0) THEN
        pp=x(3)
        rr=rnmn(1,0)
        zz=0.0d0
      ELSE IF(x(1).gt.1.0d0) THEN
        rr=rnmn(1,loop)
        zz=0.0d0
        pp=x(3)
        DO k=1,modmax
          cosf=dcos(mnumbr(k)*x(2)-nnumbr(k)*x(3))
          sinf=dsin(mnumbr(k)*x(2)-nnumbr(k)*x(3))
          zz=zz+znmn(k,loop)*sinf*x(1)
          pp=pp+pnmn(k,loop)*sinf
          rr=rr+rnmn(k,loop)*cosf*x(1)
        END DO
        rr=rr-rnmn(1,loop)*dcos(mnumbr(1)*x(2)-nnumbr(1)*x(3))*x(1)
      ELSE
        DO i=1,loop
          IF(x(1).gt.psino(i-1).and.x(1).le.psino(i)) THEN
            ii=i-1
            xx=x(1)-psino(ii)
            GO TO 50
          ENDIF
        END DO
 50     CONTINUE

        rr=0.0d0
        zz=0.0d0
        pp=x(3)

        xx2=xx*xx
        xx3=xx2*xx
        DO k=1,modmax
          xabm=x(1)**abmnum(k)
          cosf=dcos(mnumbr(k)*x(2)-nnumbr(k)*x(3))
          sinf=dsin(mnumbr(k)*x(2)-nnumbr(k)*x(3))
          rr=rr+(rnmn(k,ii)+cr(k,ii)*xx+dr(k,ii)*xx2+er(k,ii)*xx3)*xabm*cosf
          zz=zz+(znmn(k,ii)+cz(k,ii)*xx+dz(k,ii)*xx2+ez(k,ii)*xx3)*xabm*sinf
          pp=pp+(pnmn(k,ii)+cp(k,ii)*xx+dp(k,ii)*xx2+ep(k,ii)*xx3)*xabm*sinf
        END DO
      ENDIF

      RETURN
      END SUBROUTINE
!#####################################################################c
      SUBROUTINE nbiinp
!#####################################################################

!     Input the beam parameters

!--------------------------------------------------------------------c

  USE hfrmod, ONLY: ambip, amf, ane00, anexm, angle1, angle2, angle3, angle4,     &
    anwall, aphght, apwdth, ate00, atewal, atexm, ati00, atiwal, atixm, awb, b0,  &
    bamp, bdvghz, bdvgvt, bhght, bhzfoc, blengt, bvtfoc, bwdth, cnel1, cnel2, cnel3,&
    crtc, crtfe, crthe, crto, cte1, cte2, dt0, ebkev, efract, eps, er1, er1mod, er2,&
    ibouce, ic05pc, icmax, icoef, icut, icut0, ideb, ipartp, irunge, jobno, lmonte, &
    lorbp, lpi2, lxy, mb, modmxn, nhght, nhis, npart, nplotp, nprntp, npshap, npskp,&
    nptf, nshape, nwdth, phi0, phipp, psibm, qmasb, r00, r1chab, r2chab, rbore,     &
    rmj0, rmjpvt, rpivot, rpp, tlim, zb, zchab, zf
  IMPLICIT NONE

!< &inpnbi >------------------------------------------------------------

  INTEGER(4):: i, iab, nbdim
  REAL(8):: ab, an0
  REAL(8):: tolnbh=0.2d0
  LOGICAL :: lfauc

      NAMELIST /inporb/ ipartp,tlim,dt0,eps,modmxn,nprntp,nplotp, &
     &                  irunge,lmonte,nhis,psibm,b0,lpi2, &
     &                  icmax,lorbp
      NAMELIST /bin/ icut,icut0,icoef,ideb
!      dimension inj(10)
!!!     data ifpout /0,0,0,0,0,0,0,0,0,0/
!      data lneutl/0/
!      data clbcx/0.5d0/
!      data tolnbi/0.2d0/,tolnbh/0.2d0/


!< &freyin >------------------------------------------------------------
!
      NAMELIST /inpnbi/ jobno,iab,amf,zf,zb,ambip

      NAMELIST /freyin/npart,npskp,mb,bamp,efract,ebkev,bhght, &
     &         bwdth,blengt,nhght,nwdth,nshape,aphght,apwdth,npshap, &
     &         rmjpvt,rpivot, &
     &         bhzfoc,bvtfoc,bdvghz,bdvgvt,angle1,angle2,angle3, &
     &         angle4, &
     &         an0,r00,rbore,r1chab,r2chab,zchab,ibouce, &
     &         nptf,tolnbh, &
     &         ate00,ane00,anwall,atewal,cnel1,cnel2,cnel3,cte1,cte2, &
     &         ati00,atiwal,phi0,lfauc
      NAMELIST /protp1/ rpp,phipp
      NAMELIST /wbound/ er1, er2, er1mod, rmj0
      NAMELIST /zqinp/ crthe, crtc, crto, crtfe

      NAMELIST /prof/  anexm, atixm, atexm


!< &freyin >------------------------------------------------------------

      lxy=.true.
      irunge=0
      modmxn=0
      lpi2=0
      lorbp=.true.
      icmax=0
      lmonte=.false.
      lfauc=.false.
      ic05pc=1
      nbdim=1
      b0=1.0
      npart=10000
      npskp=10
      mb=1
      r1chab=350.d0
      r2chab=630.d0
      zchab=250.d0
      DO i=1,2
        bamp(i)=333.d0
        efract(i,1)=0.65d0
        efract(i,2)=0.30d0
        efract(i,3)=0.05d0
        ebkev(i)=150.d0
        bhght(i)=50.d0
        bwdth(i)=50.d0
        blengt(i)=900.d0
        nhght(i)=1
        nwdth(i)=1
        nshape(i)=1
        aphght(i)=1.d+05
        apwdth(i)=1.d+05
        npshap(i)=1
        bhzfoc(i)=900.d0
        bvtfoc(i)=900.d0
        bdvghz(i)=0.9d0
        bdvgvt(i)=0.9d0
        angle1(i)=0.d0
        angle2(i)=0.d0
        angle3(i)=15.d0
        phi0(i)=.0
        angle4(i)=0.d0
!x    rpivot(i)=rlimc
!x    rmjpvt(i)=rmajc
      END DO
!< &inpnbi >------------------------------------------------------------

      icut=2
      icoef=0
      ideb=1
      READ(5,bin)
      WRITE(6,bin)

      iab=2
      READ(5,inpnbi)
      WRITE(6,inpnbi)


      awb=iab
      qmasb=awb*1.6726485d-24
      ab=iab
!x    if(iab.eq.0) kb=0
!x    if(iab.eq.1) kb=mprot
!x    if(iab.eq.2) kb=mdeut
!x    if(iab.eq.3) kb=mtrit
!x    mbeam=kb
      READ(5,freyin)

!     <<  cal bhzfoc,bvtfoc  >>
      IF(lfauc) THEN
        CALL calfoc( angle3,rmjpvt(1),rpivot(1),blengt,bhzfoc,bvtfoc )
      ENDIF

      WRITE(6,freyin)

!........  default value
      rpp   = 325.8d0
      phipp =  72.0d0

      READ (5,protp1)
      WRITE(6,protp1)

      phipp=phipp*0.017453293d0


      READ(5,inporb)
      WRITE(6,inporb)

      READ(5,wbound)
      WRITE(6,wbound)

      READ(5,zqinp)
      WRITE(6,zqinp)

      READ(5,prof)
      WRITE(6,prof)

      RETURN
      END SUBROUTINE

!#####################################################################
      SUBROUTINE calfoc( angle3,rmaj,rlim,blengt,bhzfoc,bvtfoc )
!#####################################################################

  IMPLICIT NONE
  REAL(8),INTENT(IN) :: angle3(2), rmaj, rlim, blengt(2)
  REAL(8),INTENT(OUT):: bhzfoc(2),bvtfoc(2)
  REAL(8):: a, b, c, foc, t, th0, th

      a = rmaj
      b = rlim
      c = blengt(1)
      t = 3.141592654d0*angle3(1)/180.d0

      th0 =  (a+b)*sin(t)/a
      IF( abs(th0).ge.1.0d0 ) THEN
         bhzfoc(1) = blengt(1)
         bvtfoc(1) = blengt(1)
         RETURN
      ENDIF

      th  = acos( th0 )
      foc = c+(a+b)*cos(t)-a*sin(th)

      bhzfoc(1) = foc
      bvtfoc(1) = foc

      RETURN
      END SUBROUTINE

!#####################################################################
      SUBROUTINE nfreya
!#####################################################################

!     Calculate the NBI depositions

!--------------------------------------------------------------------c
  USE hfrmod, ONLY : amiss, angle3, aplost, bmdep, bmtot, efract, hbeam, hdep, hener,&
    hrbeam, hrlose, hrr1, hrr2, mb, mfm1, mmf, mnengy, orlost, passed, pplost, psif, &
    psivol, sgxn, sgxnmi, sivol, volume, xnel, xte
  IMPLICIT NONE

  INTEGER(4):: i, j, m
  REAL(8):: zsum, znorm, sw1(101), sw2(101), ww1(4,101), zr2(101)

!---------------------c
      CALL hflux
!---------------------c

!---------------------c
      CALL hrcal
!---------------------c

      DO 10 i=1,mfm1
        DO 10 j=1,mnengy
          IF(bmtot(j).eq.0.d0) GO TO 10
          hrbeam(j,i)=hrbeam(j,i)*volume/(psivol(i)*bmtot(j))
          hrr2(j,i)=hrr2(j,i)*volume/(psivol(i)*bmtot(j))
          hrr1(j,i)=hrr1(j,i)*volume/(psivol(i)*bmtot(j))
          DO m=1,mb
            hbeam(i,j,m)=hbeam(i,j,m)*volume/(psivol(i)*bmtot(j))
         END DO
 10   CONTINUE

      DO i=1,mmf
        zr2(i)=dfloat(i-1)/dfloat(mmf-1)
      END DO

      DO j=1,mnengy
        hrbeam(j,mmf)=0.0d0
        DO m=1,mb
          hbeam(mmf,j,m)=0.0d0
        END DO
      END DO

      DO 20 j=1,mnengy
        hrlose(j)=0.d0
        IF(bmtot(j).eq.0.d0) GO TO 20
        hrlose(j)=1.d0-bmdep(j)/bmtot(j)
 20   CONTINUE

      DO j=1,mnengy
        DO m=1,mb
          DO i=1,mfm1
            sw1(i)=hbeam(i,j,m)
          END DO
          CALL dspln(sw2,zr2,mfm1,sw1,psif,mfm1,ww1,0)
          DO i=1,mfm1
            hdep(i,4-j,m)=dmax1(0.0d0,sw2(i))
          END DO
        END DO
      END DO

      DO m=1,mb
        DO j=1,mnengy
          zsum=0.0
          DO i=1,mfm1
            zsum=zsum+hdep(i,j,m)*sivol(i)
          END DO
          znorm=1.0
          IF(bmtot(j).ne.0.0d0.and.zsum.ne.0.0d0) THEN
            znorm=volume/zsum*bmdep(j)/bmtot(j)
          ENDIF
          WRITE(6,600) znorm
          DO i=1,mfm1
            hdep(i,j,m)=hdep(i,j,m)*znorm
          END DO
        END DO
      END DO
!
600   FORMAT(1h /,' ***** hfreya znorm=*****',1pe14.4)

      RETURN


!--------------------------c
      entry fryout
!--------------------------c
!
!x    if(time.lt.tnbon) return
!x    if(time.gt.tnboff) return
!     write(6,8005)
!     write(6,7000) volume,(psivol(i), i=1,mfm1)
 7000 FORMAT(9h volumes:,1p5e12.3,/(9x,1p5e12.3))
!x    call upage(0)
      WRITE(6,8005)
 8005 FORMAT(1h,///,'  **********  output nfreya *******',///)
!     write(6,8006)
 8006 FORMAT(1h //,5x,'te',15x,'ne')
      WRITE(6,8007) (i,xte(i),xnel(i),i=1,mfm1)
 8007 FORMAT(1h ,i3,0pf10.3,1pe15.4)
      WRITE(6,8008) (i,sgxn(1,i),sgxn(2,i),sgxn(3,i),i=1,mfm1)
 8008 FORMAT(1h //,4x,'i',8x,'sgxn(1)',8x,'sgxn(2)',8x,'sgxn(3)'/(i5,1p3e15.4))
      WRITE(6,8009) (sgxnmi(j),j=1,3)
 8009 FORMAT(1h //,'  sgxnm '/(1p3e15.4))
      WRITE(6,7002) (hener(i),i=1,mnengy)
 7002 FORMAT(1h1 // ,'energy components',2x,(3f15.2))
      WRITE(6,7003) (bmtot(j), j=1,mnengy)
 7003 FORMAT(1h ,'total particles',4x,(3f15.2))
      WRITE(6,7014) (orlost(j),j=1,mnengy)
 7014 FORMAT(1h ,'particle orbit loss',(3f15.2))
      WRITE(6,7015) (aplost(j),j=1,mnengy)
 7015 FORMAT(1h ,'aperature loss',5x,(3f15.2))
      WRITE(6,7016) (amiss(j),j=1,mnengy)
 7016 FORMAT(1h ,'missed chamber',5x,(3f15.2))
      WRITE(6,7017) (passed(j),j=1,mnengy)
 7017 FORMAT(1h ,'passed through',5x,(3f15.2))
!
      WRITE(6,6000) (pplost(j),j=1,mnengy)
 6000 FORMAT(1h ,'lost at protection plate',5x,(3f15.2))

      WRITE(6,7004) (bmdep(j), j=1,mnengy)
 7004 FORMAT(1h ,'particles deposited',(3f15.2))
      WRITE(6,7005) (hrlose(j), j=1,mnengy)
 7005 FORMAT(1h ,'fraction lost',6x,(3f15.4))
      WRITE(15,7018) (hrlose(j), j=1,mnengy)
 7018 FORMAT(3f15.4)
      WRITE(15,7019) (efract(1,j), j=mnengy,1,-1)
 7019 FORMAT(3f15.4)
      WRITE(6,7006)
 7006 FORMAT(1h ///,'  beam deposition ')
      WRITE(6,7001) ((hrbeam(j,i), j=1,mnengy), i=1,mmf)
 7001 FORMAT((20x,3f15.5))
      WRITE(6,7007)
 7007 FORMAT(///)
      DO m=1,mb
         WRITE(6,7020) angle3(m)
         WRITE(6,7001)((hbeam(i,j,m),j=1,mnengy),i=1,mmf)
      END DO
!
!x    write(6,7021)
!x    do 1010 m=1,mb
!x    write(6,7001) ((hdep(i,j,m),j=1,mnengy),i=1,mmf)
 1010 CONTINUE
 7021 FORMAT(1h //,   '  hdep')
 7020 FORMAT(1h //,'      angle3  =',f10.3,//)
!cc



!******************************************************************
      CALL ppout
!******************************************************************

      RETURN
      END SUBROUTINE


!#####################################################################
      SUBROUTINE binit
!#####################################################################
!     initializes freya type neutral beam injection
!--------------------------------------------------------------------c

  USE hfrmod, ONLY : angle1, angle2, angle3, angle4, awb, bamp, beff, bpart, bpartt,&
    cang1, cang2, cang3, cang4, ebeam, ebkev, efract, hbeff, hener, hrfac, iecomp,  &
    inumpt, lbeam, mnengy, npart, sang1, sang2, sang3, sang4, vbeam
  IMPLICIT NONE
  INTEGER(4),PARAMETER:: mb=2, mb3m=3*mb-1
  INTEGER(4):: i, j, ichang, mb2
  REAL(8):: awbi, bparti, eb, hrfaci, temp, xj, xnpart
  REAL(8):: pio180=0.017453293
!!cc
!      data pi,twopi,forpi,pio180/3.141592654,6.283185308, &
!     &     12.56637062,0.017453293/,nbdim/1/

!cc
!cc
!name                =====< nbi input data >======
!name                          call nbiiut
!name                =============================
!cc
      DO i=1,mb
        ebeam(i)=1.d+03*ebkev(i)
        DO j=1,3
          xj=j
          eb=ebkev(i)/(awb*xj)
          CALL logint(eb,beff(i,j))
          hbeff(4-j,i)=beff(i,j)
        END DO
      END DO
!cc
      DO i=1,mb
        cang1(i) = cos(angle1(i)*pio180)
        cang2(i) = cos(angle2(i)*pio180)
        cang3(i) = cos(angle3(i)*pio180)
        cang4(i) = cos(angle4(i)*pio180)
        sang1(i) = sin(angle1(i)*pio180)
        sang2(i) = sin(angle2(i)*pio180)
        sang3(i) = sin(angle3(i)*pio180)
        sang4(i) = sin(angle4(i)*pio180)
      END DO
!cc
      mb2=2*mb
      DO i=1,mb
        hener(i)    = ebeam(i)/3.d0
        hener(i+mb) = ebeam(i)/2.d0
        hener(i+mb2)= ebeam(i)
      END DO
!cc
!     mb3m=3*mb-1
!...see paramter

 21   CONTINUE
      ichang=0
      DO 22 j=1,mb3m
        IF(hener(j).le.hener(j+1)) GO TO 22
        temp=hener(j)
        hener(j)=hener(j+1)
        hener(j+1)=temp
        ichang=1
 22   CONTINUE
      IF(ichang.ge.1) GO TO 21
!cc
      mnengy=3*mb

      i=1
 23   CONTINUE
      IF(hener(i).le.0.999d0*hener(i+1).or.hener(i).ge.1.001d0*hener(i+1)) GO TO 24
      mnengy=mnengy-1
      IF(i.ge.mnengy) GO TO 25
      DO j=i,mnengy
        hener(j)=hener(j+1)
      END DO
      GO TO 23
 24   CONTINUE
      i=i+1
      IF(i.ge.mnengy) GO TO 25
      GO TO 23
 25   CONTINUE
!cc
      DO i=1,mnengy
        DO j=1,mb
          IF(ebeam(j).gt.0.999*hener(i).and.ebeam(j).lt.1.001*hener(i)) lbeam(j,1)=i
          temp=ebeam(j)/2.d0
          IF(temp.gt.0.999*hener(i).and.temp.lt.1.001*hener(i)) lbeam(j,2)=i
          temp=ebeam(j)/3.d0
          IF(temp.gt.0.999*hener(i).and.temp.lt.1.001*hener(i)) lbeam(j,3)=i
        END DO
      END DO
   30 CONTINUE
!cc
      awbi=1.d0/awb
      DO i=1,mnengy
        vbeam(i)=1.384d+06*sqrt(awbi*hener(i))
      END DO
!cc
      bpartt=0.
      DO i=1,mb
        bparti=6.2418d+18*bamp(i)
        DO j=1,3
!     bpart(i,j)=bparti*efract(i,j)*beff(i,j)*j
          bpart(i,j)=bparti*efract(i,j)*j
          bpartt=bpartt+bpart(i,j)
        END DO
      END DO

      xnpart = npart
      hrfac  = bpartt/xnpart
      hrfaci = 1.d0/hrfac

      DO i=1,mb
        inumpt(i)=0.d0
        DO j=1,3
          iecomp(i,j)=hrfaci*bpart(i,j)
          inumpt(i)=inumpt(i)+iecomp(i,j)
        END DO
      END DO
!cc
      RETURN
      END SUBROUTINE


!#####################################################################
      SUBROUTINE hflux
!#####################################################################

  USE hfrmod , ONLY : anexm, atexm, cnel1, cte1, dps, dpsii, loop, loop, mfm1, mmf,&
    nn, psi, psicut, psif, psimin, psivol, r1chab, r2chab, rmax, rmin, sivol, vol, &
    volume, xnel, xnion, xte, zchab, zmax, zmin
  IMPLICIT NONE
  INTEGER(4):: i,jj
  REAL(8):: dpsi,zps, ww(4,nn), zpsi(0:nn), zr2(61), zvol1(61), zvol(61)

      rmax=r2chab
      rmin=r1chab
      zmax=zchab
      zmin=-zmax
!....
!....
!ccccccccc
!      mmf=61
!      mfm1=mmf-1
!

      psimin=0.0d0
      psicut=1.0d0
      dpsi=(psicut-psimin)/float(mfm1)
      dps=dpsi
      dpsii=1.d0/dps

      DO i=1,mmf
        psif(i)=dpsi*float(i-1)+psimin
        zr2(i)=(dfloat(i-1)/dfloat(mmf-1))**2
      END DO
!
      IF(cnel1.lt.0.1d0) cnel1=0.1d0
      IF(cte1.lt.0.1d0) cte1=0.1d0

      DO i=1,mmf
        zps=sqrt(psif(i))

        IF (i.eq.1) THEN
          xnel(1)= anexm(1)
          xte(1) = atexm(1)
        ELSE
          DO jj=1,10
            xnel(i)=xnel(i)+anexm(jj)*zps**(jj-1)
            xte(i) =xte(i) +atexm(jj)*zps**(jj-1)
          END DO
        END IF
!.... low limit
        IF (xte(i).lt.10.0d0) THEN
          xte(i)=10.0d0
        END IF

        IF (xnel(i).lt.1.0d11) THEN
          xnel(i)=1.0d11
        END IF

        xnion(i)=xnel(i)

      END DO
!
      DO i=0,loop
        zpsi(i)=psi(i)/psi(loop)
      END DO

      CALL dspln(zvol,psif,mmf,vol(0),zpsi(0),loop+1,ww,0)
      CALL dspln(zvol1,zr2,mmf,vol(0),zpsi(0),loop+1,ww,0)

      DO i=1,mfm1
        psivol(i)=(zvol(i+1)-zvol(i))*1.d6
        sivol(i)=(zvol1(i+1)-zvol1(i))*1.d6
      END DO

      volume=vol(loop)*1.d6
!!      vol1=zvol1(mmf)               ! 200710
!
!     write(6,600) volume,(psivol(i),i=1,mfm1)
 600  FORMAT(1h /,' ******** volume ****',1pe12.3,/' ***** psivol***', &
     &     /(1p10e12.3))
!     write(6,601) vol1,(sivol(i),i=1,mfm1)
 601  FORMAT(1h /,' ******** volume ****',1pe12.3,/' *****  sivol***', &
     &     /(1p10e12.3))

      RETURN
      END SUBROUTINE


!#####################################################################
      SUBROUTINE hrcal
!#####################################################################

!     Calculates the ion deposition from the neutral beam injection.

!--------------------------------------------------------------------c

  USE hfrmod , ONLY : amiss, aplost, bmdep, bmtot, dtime, hbeam, hener, hr2, hrbeam, &
    hrr1, hrr2, ibouce, iecomp, inters, ipacc, ipt, iptt, lbeam, lostpp, mb, mfm1,   &
    mlost, mmf, mnengy, nbeam, nengy, newpar, noplos, npart, npskp, npzone, orlost,  &
    passed, pplost, rpos, vbeam, vx0, vy0, vz0, xpos, ypos, zpos
  IMPLICIT NONE
  INTEGER(4):: i, ii, j, icount, ictp, iflag, ipar, ipl, iploti, jzone
  REAL(8):: wgt


      ipt   = min(500,npart)
      ipl   = max(1,npart/ipt)
      ii    = 0
      ictp  = 0
      iploti= 0
      ipacc = 0

      DO i=1,3
        iptt(i) = 0
      END DO

      DO i=1,mnengy
        bmtot(i) = 0.d0
        bmdep(i) = 0.d0
        orlost(i)= 0.d0
        aplost(i)= 0.d0
        amiss(i) = 0.d0
        passed(i)= 0.d0
        pplost(i)= 0.d0

        DO j=1,mfm1
          hrbeam(i,j) = 0.d0
          hr2(i,j)    = 0.d0
          hrr2(i,j)   = 0.d0
          hrr1(i,j)   = 0.d0
          hbeam(j,i,1)= 0.d0
          hbeam(j,i,2)= 0.d0
        END DO
      END DO

      CALL crsec4
!
      noplos = 0


      DO nbeam=1,mb
        wgt=1.
!
        DO 200 i=1,3

          IF(iecomp(nbeam,i).eq.0) GO TO 200
          nengy=lbeam(nbeam,i)
          WRITE(7,600) i,nengy,hener(nengy),vbeam(nengy)
 600      FORMAT(1h ,' ***** sub. hrcal ,i,nengy,hener,vbeam***',/, 2i5,1p2e13.3)

          newpar=1

          DO 20 ipar=1,iecomp(nbeam,i)
            ipacc=ipacc+1
!x          write(7,601) ipacc
 601        FORMAT(1h //,' ******** ipacc *********',i10//)
            icount=ipar-1
            IF(mod(icount,npskp).eq.0) newpar=1


            CALL inject

!x          write(7,602) ipacc,mlost,inters,istart
 602        FORMAT(1h //,' ****after inject,ipacc,mlost,inters**',4i7/)

            WRITE(21,603) xpos, ypos, zpos, vx0, vy0, vz0
 603        FORMAT(6d24.17)


            IF(mod(iploti,ipl).eq.0.and.ii.lt.500) THEN
              ii=ii+1
!              xposd(ii)=xpos                 !! 200710 nec
!              yposd(ii)=ypos
!              zposd(ii)=zpos
!              rposd(ii)=rpos
!              iptss(ii)=i
              iptt(i)=ii
              ipt=ii
            ENDIF

 11         CONTINUE
            iploti=iploti+1
!
            newpar=0
            bmtot(nengy)=bmtot(nengy)+wgt
            IF(mlost.ne.0) GO TO 95
            IF(inters.eq.0) GO TO 85
            IF(lostpp.ne.0) GO TO 500
            IF(npzone.ge.mmf) GO TO 75
            IF(ibouce.ne.1) GO TO 30
!
            ictp=ictp+1
!
            CALL bounce(iflag,ictp)
            iflag=1
            IF(iflag.ne.1) GO TO 60
            DO jzone=1,mmf
              hbeam(jzone,nengy,nbeam)=hbeam(jzone,nengy,nbeam)+dtime(jzone)
              hrbeam(nengy,npzone)=hrbeam(nengy,npzone)+dtime(jzone)
              hrr2(nengy,jzone)=hrr2(nengy,jzone)+dtime(jzone)
            END DO
            GO TO 50
 30         hrbeam(nengy,npzone)=hrbeam(nengy,npzone)+wgt
            hbeam(npzone,nengy,nbeam)=hbeam(npzone,nengy,nbeam)+wgt
 50         hrr1(nengy,npzone)=hrr1(nengy,npzone)+wgt
            bmdep(nengy)=bmdep(nengy)+wgt
            GO TO 20
 60         orlost(nengy)=orlost(nengy)+1.
            WRITE(6,610) ictp,iflag
 610        FORMAT(1h /,' ***** particle lost***,ictp,iflag',2i10)
10000       FORMAT(1h ,'iflag=',i4)
            GO TO 20
 75         passed(nengy) = passed(nengy)+1.d0
            GO TO 20
 85         amiss(nengy)  = amiss(nengy)+1.d0
            GO TO 20
 95         aplost(nengy) = aplost(nengy)+1.d0
            GO TO 20
 500        pplost(nengy) = pplost(nengy)+1.d0
 20       CONTINUE
 200    CONTINUE
      END DO

      RETURN
      END SUBROUTINE


!#####################################################################
      SUBROUTINE logint(x,y)
!#####################################################################

!     interpolates y(x) quadratically and logarithmically
!     calculation of neutralizing efficiency?

!--------------------------------------------------------------------c

  IMPLICIT NONE
  REAL(8),INTENT(IN) :: x
  REAL(8),INTENT(OUT):: y
  INTEGER(4):: i0, i0end, im, ip, mdat, mdatm
  REAL(8):: d0, d0m, d0p, dm, dm0, dmp, dp, dp0, dpm, fac0, facm, facp, ylog, ylog0,&
    ylogm, ylogp
  REAL(8):: xdat(15)=(/4.d0,6.d0,8.d0,10.d0,20.d0,30.d0,40.d0,60.d0, &
                    80.d0,100.d0,200.d0,300.d0,400.d0,600.d0,800.d0/)
  REAL(8):: ydat(15)=(/8.95d-01,8.75d-01,8.70d-01,8.65d-01,8.20d-01, &
                       7.25d-01,6.25d-01,4.40d-01,2.90d-01,1.90d-01, &
                       2.40d-02,5.25d-03,1.20d-03,1.60d-04,5.40d-05/)

      mdat=15
      mdatm=mdat-1
      DO i0=2,mdatm
        i0end=i0
        IF(xdat(i0).ge.x) GO TO 11
      END DO
 11   i0=i0end
      IF(i0.gt.mdatm) i0=mdatm
      im=i0-1
      ip=i0+1
      ylogp=dlog(ydat(ip))
      ylog0=dlog(ydat(i0))
      ylogm=dlog(ydat(im))
      dm=x-xdat(im)
      d0=x-xdat(i0)
      dp=x-xdat(ip)
      d0m=xdat(i0)-xdat(im)
      dp0=xdat(ip)-xdat(i0)
      dpm=xdat(ip)-xdat(im)
      dm0=-d0m
      d0p=-dp0
      dmp=-dpm
      facm=d0*dp/(dm0*dmp)
      fac0=dm*dp/(d0m*d0p)
      facp=dm*d0/(dpm*dp0)
      ylog=facm*ylogm+fac0*ylog0+facp*ylogp
      y=exp(ylog)
!cc
      RETURN
      END SUBROUTINE


!#####################################################################
      SUBROUTINE crsec4
!#####################################################################
!                                                       1999/2/17
!     Calculates ionization cross sections from fitted results of
!     R. K. Janev, et al., Nuclear Fusion 29 (1989) 2125.

!                                    code by S. Murakami
!--------------------------------------------------------------------c

  USE hfrmod, ONLY : awb, crtc, crtfe, crthe, crthe, crto, hener, mfm1, mfm1,      &
    mnengy, sgxn, sgxne, sgxni, sgxnmi, vbeam, xnel, xnion, xte, zf
  USE hfrmod, ONLY : allocate_crsec4, deallocate_crsec4
  IMPLICIT NONE
  INTEGER(4):: i, ii, j, if7, ia, ja, ka
  REAL(8):: algee, algni, algte, aloge, alogt, awbi, czft, czqc, czqfe, czqhe, czqo,&
    delsig, e, expo, sgxnm, sigcx, sigee, sigi, sigino, sigz, sss1, sszc, sszfe,    &
    sszhe, sszo, veli
  REAL(8):: sgvxne(61), aaa(2,3,2), bbbhe(3,2,2), bbbc(3,2,2), bbbo(3,2,2),  bbbfe(3,2,2)
  REAL(8):: cfione(7)=(/-3.173850d+01,1.143818d+01,-3.833998d0,7.046692d-01, &
                        -7.431486d-02,4.153749d-03,-9.486967d-05/)
  REAL(8):: cfionp(7)=(/-4.203309d+01,3.557321d0,-1.045134d0,3.139238d-01, &
                        -7.454475d-02,8.459113d-03,-3.495444d-04/)

      save
!-------------------------------------------------------c
!     aaa for proton

      aaa(1,1,1)=  4.40d0
      aaa(1,1,2)= -2.49d-2
      aaa(1,2,1)=  7.46d-2
      aaa(1,2,2)=  2.27d-3
      aaa(1,3,1)=  3.16d-3
      aaa(1,3,2)= -2.78d-5
      aaa(2,1,1)=  2.30d-1
      aaa(2,1,2)= -1.15d-2
      aaa(2,2,1)= -2.55d-3
      aaa(2,2,2)= -6.20d-4
      aaa(2,3,1)=  1.32d-3
      aaa(2,3,2)=  3.38d-5

!-------------------------------------------------------c

!-------------------------------------------------------c
!    bbb for He impurity

      bbbhe(1,1,1)= -2.36d0
      bbbhe(1,1,2)=  1.85d-1
      bbbhe(1,2,1)= -2.50d-1
      bbbhe(1,2,2)= -3.81d-2
      bbbhe(2,1,1)=  8.49d-1
      bbbhe(2,1,2)= -4.78d-2
      bbbhe(2,2,1)=  6.77d-2
      bbbhe(2,2,2)=  1.05d-2
      bbbhe(3,1,1)= -5.88d-2
      bbbhe(3,1,2)=  4.34d-3
      bbbhe(3,2,1)= -4.48d-3
      bbbhe(3,2,2)= -6.76d-4

!-------------------------------------------------------c
!    bbb for C impurity

      bbbc(1,1,1)= -1.49d0
      bbbc(1,1,2)= -1.54d-2
      bbbc(1,2,1)= -1.19d-1
      bbbc(1,2,2)= -1.50d-2
      bbbc(2,1,1)=  5.18d-1
      bbbc(2,1,2)=  7.18d-3
      bbbc(2,2,1)=  2.92d-2
      bbbc(2,2,2)=  3.66d-3
      bbbc(3,1,1)= -3.36d-2
      bbbc(3,1,2)=  3.41d-4
      bbbc(3,2,1)= -1.79d-3
      bbbc(3,2,2)= -2.04d-4

!-------------------------------------------------------c
!    bbb for O impurity

      bbbo(1,1,1)= -1.41d0
      bbbo(1,1,2)= -4.08d-4
      bbbo(1,2,1)= -1.08d-1
      bbbo(1,2,2)= -1.38d-2
      bbbo(2,1,1)=  4.77d-1
      bbbo(2,1,2)=  1.57d-3
      bbbo(2,2,1)=  2.59d-2
      bbbo(2,2,2)=  3.33d-3
      bbbo(3,1,1)= -3.05d-2
      bbbo(3,1,2)=  7.35d-4
      bbbo(3,2,1)= -1.57d-3
      bbbo(3,2,2)= -1.86d-4

!-------------------------------------------------------c
!    bbb for Fe impurity

      bbbfe(1,1,1)= -1.03d0
      bbbfe(1,1,2)=  1.06d-1
      bbbfe(1,2,1)= -5.58d-2
      bbbfe(1,2,2)= -3.72d-3
      bbbfe(2,1,1)=  3.22d-1
      bbbfe(2,1,2)= -3.75d-2
      bbbfe(2,2,1)=  1.24d-2
      bbbfe(2,2,2)=  8.61d-4
      bbbfe(3,1,1)= -1.87d-2
      bbbfe(3,1,2)=  3.53d-3
      bbbfe(3,2,1)= -7.43d-4
      bbbfe(3,2,2)= -5.12d-5

!-------------------------------------------------------c

!------------ electron ionization ---------------+
      DO i=1,mfm1
        alogt=0.d0
        IF(xte(i).gt.1.d0)   alogt = log(xte(i))
        IF(xte(i).gt.1.d+05) alogt = 11.51d0
        expo=cfione(7)
        DO ii=1,6
          if7=7-ii
          expo=expo*alogt+cfione(if7)
        END DO
        sgvxne(i)=exp(expo)*xnel(i)
      END DO

!--- charge exchange and ion ionization -----------+

      CALL allocate_crsec4

      awbi=1.d0/awb
      DO j=1,mnengy
      veli=1.d0/vbeam(j)
      DO i=1,mfm1

        e=awbi*hener(j)
        aloge=dlog10(e)
        sigcx=0.6937d-14*(1.d0-0.155d0*aloge)**2/(1.d0+0.1112d-14*e**3.3d0)

!        aloge=aloge*2.302585093d0-6.907755279d0

        aloge=dlog(e*1.d-3)

        IF(aloge.gt.-2.30258d0) THEN

          expo=cfionp(7)
          DO ii=1,6
            if7=7-ii
            expo=expo*aloge+cfionp(if7)
          END DO
!
          sigino=exp(expo)
!
          sss1=0.0d0

          DO ia=1,2
            DO ja=1,3
              DO ka=1,2
                algni = (dlog(xnion(i)/1.d13))**(ja-1)
                algte = (dlog(xte(i)*1.d-3))**(ka-1)
                algee = (aloge)**(ia-1)
                sss1  = sss1 +aaa(ia,ja,ka)*algee*algni*algte

!     print 1112, ia,ja,ka,aaa(ia,ja,ka),algni,algte,algee

              END DO
            END DO
          END DO

          sszhe=0.0d0
          sszc =0.0d0
          sszo =0.0d0
          sszfe=0.0d0
!
          DO ia=1,3
            DO ja=1,2
              DO ka=1,2
                algni = (dlog(xnion(i)/1.d13))**(ja-1)
                algte = (dlog(xte(i)*1.d-3))**(ka-1)
                algee = (aloge)**(ia-1)
                sszhe = sszhe +bbbhe(ia,ja,ka)*algee*algni*algte
                sszc  = sszc  +bbbc(ia,ja,ka)*algee*algni*algte
                sszo  = sszo  +bbbo(ia,ja,ka)*algee*algni*algte
                sszfe = sszfe +bbbfe(ia,ja,ka)*algee*algni*algte
!
!     print *, ia,ja,ka,sssz,bbb(ia,ja,ka),algee,algni,algte
              END DO
            END DO
          END DO
!
!           print *, sszhe, sszc, sszo, sszfe
!
          sigi= exp(sss1)/(e*1.d-3)*1.d-16
!
          IF (zf.eq.1.d0) THEN
            sigz = 1.d0
          ELSE
!
            czqhe = 2.d0
            czqc  = 6.d0
            czqo  = 8.d0
            czqfe = 20.d0
!
            czft = ( (czqhe-1.d0)*czqhe*crthe +(czqc -1.d0)*czqc *crtc &
                    +(czqo -1.d0)*czqo *crto  +(czqfe-1.d0)*czqfe*crtfe )/(zf-1.d0)
!
            sigz= 1.d0     &
                   + (crthe/czft)*(czqhe-1.d0)*czqhe*sszhe &
                   + (crtc /czft)*(czqc -1.d0)*czqc *sszc &
                   + (crto /czft)*(czqo -1.d0)*czqo *sszo &
                   + (crtfe/czft)*(czqfe-1.d0)*czqfe*sszfe
!
          END IF
!
!
        ELSE
          sigi=0.0d0
        END IF
!
!     c       sgxncx(j,i)=sigcx*xnion(i)
!
        sgxni(j,i) =sigi*xnion(i)
        sgxne(j,i) =sgvxne(i)*veli
!     c       sgxn(j,i)  =sgxncx(j,i)+sgxni(j,i)+sgxne(j,i)
        sgxn(j,i)  =sgxni(j,i)*sigz
!
        sigee   =sgvxne(i)*veli/xnion(i)
        delsig = sigi/(sigino+sigcx+sigee)

        PRINT 1111, i,delsig,sigz,sigi,sigino, sigcx,sigee,xte(i),xnion(i)

 1111   FORMAT(' i,dl,sigz,sigi,sigino,xte,xni=',i4,8e12.3)
 1112   FORMAT(' ijk, aaa, nTE=',3i4,4e12.3)
!
        END DO
      END DO

      CALL deallocate_crsec4

!
!     calculate the maximum cross section
!     in ordert to determine the time step

      DO j=1,mnengy
        sgxnm=sgxn(j,1)
        DO i=2,mfm1
          IF(sgxnm.lt.sgxn(j,i)) sgxnm=sgxn(j,i)
        END DO
        sgxnmi(j)=1.d0/sgxnm/5.d0
      END DO
!cc
      RETURN
      END SUBROUTINE

!#####################################################################
      SUBROUTINE inject
!#####################################################################

!     carries out injection from particle generation to ionization

!--------------------------------------------------------------------c

  USE hfrmod, ONLY : dpsii, er1, er1mod, er2, inters, lostpp, mfm1, mlost, mmf,   &
    nengy, newpar, noplos, npzone, phi, phis, prpos, psicut, psimin, psis, pzpos, &
    r1chab, rmax, rmj0, rpos, rr0, sgxn, sgxnmi, thets, vbeam, vx0, vx00, vy0,    &
    vy00, vz0, vz00, x00, xpos, y00, ypos, z00, zmax, zmin, zpos
  IMPLICIT NONE
  INTEGER(4):: i, iin, isol, iwbnd
  REAL(8):: aa, arg, bb, cc, cn4, dfac, dwbnd, dwbnd1, ps, rminer, sn4, sphi, sqr,&
    tbstp, thet, tin, tout, tstp, tt, wbnd, wbnd1, xpos1, xx, ypos1, yy, zpos1, rr,&
    zz, edge(4), timsol(6)
  REAL(8):: ran3
  INTEGER(4):: idum=-12999


      lostpp = 0
      ps     = 2.d0

      IF(newpar.eq.0) GO TO 100
!
      CALL sorspt
      CALL rotat1

      IF(mlost.ne.0) RETURN

      edge(1)=r1chab
      edge(2)=rmax
      edge(3)=zmin
      edge(4)=zmax
      isol=0

      aa=vx0**2+vy0**2
      IF(aa.eq.0.d0) GO TO 90
      bb=2.d0*(vx0*x00+vy0*y00)

      DO 10 i=1,2
        cc=x00**2+y00**2-edge(i)**2
        arg=bb**2-4.*aa*cc
        IF(arg) 10,20,30
 30     sqr=sqrt(arg)
        tt=(-bb-sqr)/(2.d0*aa)
        zz=z00+vz0*tt
        IF(zz.lt.zmin) GO TO 40
        IF(zz.gt.zmax) GO TO 40
        isol=isol+1
        timsol(isol)=tt
 40     tt=(-bb+sqr)/(2.d0*aa)
        zz=z00+vz0*tt
        IF(zz.lt.zmin) GO TO 10
        IF(zz.gt.zmax) GO TO 10
        isol=isol+1
        timsol(isol)=tt
        GO TO 10
 20     tt=-bb/(2.d0*aa)
        zz=z00+vz0*tt
        IF(zz.lt.zmin) GO TO 10
        IF(zz.gt.zmax) GO TO 10
        isol=isol+1
        timsol(isol)=tt
 10   CONTINUE
!cc
 90   DO 50 i=3,4
        IF(vz0.eq.0.d0) GO TO 50
        tt=(edge(i)-z00)/vz0
        xx=x00+vx0*tt
        yy=y00+vy0*tt
        rr=sqrt(xx**2+yy**2)
        IF(rr.lt.r1chab) GO TO 50
        IF(rr.gt.rmax) GO TO 50
        isol=isol+1
        timsol(isol)=tt
 50   CONTINUE
!cc
      IF(isol.gt.1) GO TO 70
      inters=0
      RETURN
!cc
 70   CONTINUE
      tin=timsol(1)
      iin=1
      DO 60 i=2,isol
        IF(timsol(i).ge.tin) GO TO 60
        iin=i
        tin=timsol(i)
 60   CONTINUE
      tout=1.e+25
      DO 80 i=1,isol
        IF(i.eq.iin) GO TO 80
        IF(timsol(i).lt.tout) tout=timsol(i)
 80   CONTINUE
      x00=x00+vx0*tin
      y00=y00+vy0*tin
      z00=z00+vz0*tin
      inters=1
!cc
 100  CONTINUE
!!      istart=0               !200710
      IF(mlost.ne.0) RETURN
      IF(inters.eq.0) RETURN
      xpos=x00
      ypos=y00
      zpos=z00
      tt=tin

!     xi = xpos
!     yi = ypos
!     xf = xi+vx0*(tout-tin)
!     yf = yi+vy0*(tout-tin)
!     grada  = tan(phipp)
!     gradb  = (yf-yi)/(xf-xi)
!     xs = (yi-gradb*xi)/(grada-gradb)
!     ys = grada*xs
!     rs = sqrt(xs**2+ys**2)


      iwbnd = 0
!
 110  CONTINUE

!     dfac=-dlog(rnfl(1))
      dfac=-dlog(ran3(idum))
      tstp=dfac*sgxnmi(nengy)/vbeam(nengy)
      IF(tt+tstp.ge.tout) tstp=tout-tt
      tt=tt+tstp
      xpos=xpos+vx0*tstp
      ypos=ypos+vy0*tstp
      zpos=zpos+vz0*tstp
      rpos=sqrt(xpos**2+ypos**2)

!-----------------------------------------c
!     wall boundary   (helical)
!-----------------------------------------c

      cn4 = (xpos**4 -6.d0*xpos**2*ypos**2 +ypos**4)/ rpos**4
      sn4 = (4.d0*xpos*ypos*(xpos**2 -ypos**2))     / rpos**4

      rminer = rpos -rmj0

      wbnd = (rminer*cn4 -zpos*sn4)**2/(er1 +er1mod*cn4)**2 &
            +(rminer*sn4 +zpos*cn4)**2/ er2**2

      IF ((wbnd.gt.1.d0).and.(iwbnd.eq.1)) THEN

        tbstp = 2.d-9

 700    CONTINUE
        xpos1 = xpos
        ypos1 = ypos
        zpos1 = zpos
        wbnd1 = wbnd
!
        xpos  = xpos1 -vx0*tbstp
        ypos  = ypos1 -vy0*tbstp
        zpos  = zpos1 -vz0*tbstp
        rpos  = sqrt(xpos**2+ypos**2)
!
        cn4 = (xpos**4 -6.d0*xpos**2*ypos**2 +ypos**4)/ rpos**4
        sn4 = (4.d0*xpos*ypos*(xpos**2 -ypos**2))     / rpos**4
!
        rminer = rpos -rmj0
!
        wbnd = (rminer*cn4 -zpos*sn4)**2/(er1 +er1mod*cn4)**2 &
              +(rminer*sn4 +zpos*cn4)**2/ er2**2
!
        IF (wbnd.gt.1.d0) GO TO 700
!
        dwbnd1 = wbnd1 -1.d0
        dwbnd  = 1.d0  -wbnd
!
        xpos   = (xpos1*dwbnd +xpos*dwbnd1)/(dwbnd1 +dwbnd)
        ypos   = (ypos1*dwbnd +ypos*dwbnd1)/(dwbnd1 +dwbnd)
        zpos   = (zpos1*dwbnd +zpos*dwbnd1)/(dwbnd1 +dwbnd)
!
        npzone = 9999
        noplos=noplos+1
        prpos(noplos)=rpos
        pzpos(noplos)=zpos

!cx   write(7,604) tt,tstp,tout,rpos,zpos,dfac
!
        RETURN
!
      END IF
!
      IF ((wbnd.lt.1.d0).and.(iwbnd.eq.0)) THEN
        iwbnd = 1
      END IF
!
!------------------------------------------c

 500  CONTINUE
!x    write(7,604) tt,tstp,tout,rpos,zpos,dfac
 604  FORMAT(' ** tt,tstp,tout,rpos,zpos,dfac**',1p6e12.3)
      IF(tt.lt.tout) GO TO 120
      npzone=mmf+5
      RETURN

 120  CONTINUE

      phi=atan2(ypos,xpos)
!     write(7,600) xpos,ypos,rpos,rpos-r0*100.,zpos,ps
 600  FORMAT(1h ' **** sub.inject**',1p6e12.2)
      rr=1.e-2*rpos
      zz=1.e-2*zpos
!x    write(6,999) rr,zz
 999  FORMAT(//1x,' ***rr & zz***',1p2e12.3)
      CALL findps(rr,phi,zz,ps,thet,sphi)

      psis=ps
      thets=thet
      phis=sphi
      vx00=vx0
      vy00=vy0
      vz00=vz0
      rr0=rr

      IF(ps.lt.psimin.or.ps.gt.psicut) GO TO 110
      npzone=(ps-psimin)*dpsii+1.
      IF(npzone.gt.mfm1) npzone=mfm1
!!      istart=1              ! 200710
      IF(ran3(idum).gt.sgxn(nengy,npzone)*sgxnmi(nengy)) GO TO 110
!c      kdep=kdep+1
!x    write(7,601) ipacc,kdep,ps,xpos,ypos,rpos-100.*r0,zpos
 601  FORMAT(1h //,' *** particle deposited**',2i7,1p5e13.3/)
      RETURN
      END SUBROUTINE

!#####################################################################
      SUBROUTINE rotat1
!#####################################################################
!     3 rotations, 2 translations:
!     rotation 1: y axis to horizontal about x
!     translation 1: origin to pivot pt
!     rotation 2: x axis to horizontal about y
!     rotation 3: x axis through pivot pt and torus center
!                 about z
!     translation 2: origin to torus center

  USE hfrmod, ONLY : aphght, apwdth, blengt, cang1, cang2, cang3, cang4, mlost,    &
    nbeam, npshap, phi0, rmjpvt, rpivot, sang1, sang2, sang3, sang4, vx0, vy0, vz0,&
    x00, y00, z00
  IMPLICIT NONE
  REAL(8):: temp, tpvt, zcos, zphi, zsin
!cc
      IF(sang4(nbeam).eq.0.) GO TO 10
      zcos=cang4(nbeam)
      zsin=sang4(nbeam)
      temp=y00*zcos-z00*zsin
      z00=y00*zsin+z00*zcos
      y00=temp
      temp=vy0*zcos-vz0*zsin
      vz0=vy0*zsin+vz0*zcos
      vy0=temp
   10 CONTINUE
!cc
      tpvt=-blengt(nbeam)/vx0
      x00=0.
      y00=y00+vy0*tpvt
      z00=z00+vz0*tpvt
      mlost=0
      IF(apwdth(nbeam).gt.1.e+03) GO TO 20
      IF(npshap(nbeam).ne.1) GO TO 21
      IF(y00**2+z00**2.gt.(0.5*apwdth(nbeam))**2) mlost=1
      GO TO 20
 21   IF(abs(z00).gt.0.5*aphght(nbeam)) mlost=1
      IF(abs(y00).gt.0.5*apwdth(nbeam)) mlost=1
 20   CONTINUE
      IF(mlost.ne.0d0) RETURN
!cc
      IF(sang2(nbeam).eq.0.d0) GO TO 30
      zcos=cang2(nbeam)
      zsin=sang2(nbeam)
      temp=x00*zcos+z00*zsin
      z00=-x00*zsin+z00*zcos
      x00=temp
      temp=vx0*zcos+vz0*zsin
      vz0=-vx0*zsin+vz0*zcos
      vx0=temp
 30   CONTINUE
!cc
      IF(sang3(nbeam).eq.0.d0) GO TO 40
      zcos=cang3(nbeam)
      zsin=sang3(nbeam)
      temp=x00*zcos+y00*zsin
      y00=-x00*zsin+y00*zcos
      x00=temp
      temp=vx0*zcos+vy0*zsin
      vy0=-vx0*zsin+vy0*zcos
      vx0=temp
 40   CONTINUE
!cc
      x00=x00+rmjpvt(nbeam)+rpivot(nbeam)*cang1(nbeam)
      z00=z00-rpivot(nbeam)*sang1(nbeam)
      zphi=phi0(nbeam)*0.017453293
      zcos=cos(zphi)
      zsin=sin(zphi)
      temp=x00*zcos-y00*zsin
      y00=x00*zsin+y00*zcos
      x00=temp
      temp=vx0*zcos-vy0*zsin
      vy0=vx0*zsin+vy0*zcos
      vx0=temp
!cc
      RETURN
      END SUBROUTINE


!#####################################################################
      SUBROUTINE sorspt
!#####################################################################
!cc
!cc   generates a particle at the injector surface with
!cc   coordinates and velocities x00,y00,z00,vx0,vy0,vz0
!cc

  USE hfrmod, ONLY : bdvghz, bdvgvt, bhght, bhzfoc, bvtfoc, bwdth, nbeam, nengy,  &
    nhght, nshape, nwdth, vbeam, vx0, vy0, vz0, x00, y00, z00
  IMPLICIT NONE
  REAL(8):: dvy, dvz, vdx, vdy, vdz, vsqrt, zsqu
  INTEGER(4):: idum=-12879
  REAL(8):: pio180=0.017453293, thy=0.d0, thz=0.d0
  REAL(8):: ran3

!cc
      x00=0.
      IF(nshape(nbeam).ne.1) GO TO 20
!cc
      IF(nwdth(nbeam).ne.1) GO TO 11
   12 y00=ran3(idum)-0.5
      z00=ran3(idum)-0.5
      zsqu=y00**2+z00**2
      IF(zsqu.gt.0.25) GO TO 12
      y00=y00*bwdth(nbeam)
      z00=z00*bwdth(nbeam)
      GO TO 13
   11 CONTINUE
   13 CONTINUE
      GO TO 30
!cc
   20 CONTINUE
      IF(nwdth(nbeam).ne.1) GO TO 21
      y00=bwdth(nbeam)*(ran3(idum)-0.5)
      GO TO 22
   21 CONTINUE
   22 CONTINUE
      IF(nhght(nbeam).ne.1) GO TO  23
      z00=bhght(nbeam)*(ran3(idum)-0.5)
      GO TO 30
   23 CONTINUE
!cc
   30 CONTINUE
      vdx=-1.
      vdy=-y00/bhzfoc(nbeam)
      vdz=-z00/bvtfoc(nbeam)
      vsqrt=1./sqrt(vdx**2+vdy**2+vdz**2)
      vx0=vbeam(nengy)*vdx*vsqrt
      vy0=vbeam(nengy)*vdy*vsqrt
      vz0=vbeam(nengy)*vdz*vsqrt
!cc
      dvz=pio180*bdvgvt(nbeam)
      dvy=pio180*bdvghz(nbeam)
      CALL gausf(dvz,dvy,0.d0,0.d0,thz,thy)
      vz0=vz0-thz*vx0
      vy0=vy0-thy*vx0
!cc
      RETURN
      END SUBROUTINE

!#####################################################################
      SUBROUTINE gausf(zs1,zs2,za1,za2,z1,z2)
!#####################################################################
!cc
!cc   samples from two gaussians with centers at za1, za2
!cc   and widths zs1, zs2
!cc
  IMPLICIT NONE
  REAL(8),INTENT(IN) :: zs1, zs2, za1, za2
  REAL(8),INTENT(OUT):: z1, z2
  INTEGER(4):: idum=-12345
  REAL(8):: zcos, zsin, zsq, zt
  REAL(8):: ran3


   10 CONTINUE
      z1=2.*ran3(idum)-1.
      z2=2.*ran3(idum)-1.
      zsq=z1*z1+z2*z2
      IF(zsq.gt.1..or.zsq.le.0.) GO TO 10
      zsq=1./zsq
      zsin=2.*z1*z2*zsq
      zcos=(z1*z1-z2*z2)*zsq
      zt=sqrt(-2.*dlog(ran3(idum)))
      z1=zs1*zt*zcos+za1
      z2=zs2*zt*zsin+za2
!cc
      RETURN
      END SUBROUTINE
!      function rnfl(i)
!      implicit real*8 (a-h,o-z)
!      data ix/1234567/
!      rmod=2147483648.d0
!      rlamda=65539.
!      if(ix) 2,4,6
!   2  ix=-ix
!      go to 6
!   4  ix=3907
!   6  w=ix
!      w=rlamda*w
!      if(w-rmod) 20,10,10
! 10   icover=w/rmod
!      w=w-float(icover)*rmod
!  20  ix=w
!      rnfl=w/rmod
!      return
!      end


!-----------------------------------------------------------------------
      REAL(8) FUNCTION gasdev (idum)
!-----------------------------------------------------------------------

!           returns a normally distributed deviate with zero mean and
!           unit variance, using ran3(idum) as the source of uniform
!           deviates

  IMPLICIT NONE
  INTEGER(4),INTENT(IN):: idum
  INTEGER(4):: iset=0
  REAL(8):: fac, gset, r, v1, v2
  REAL(8):: ran3            ! function
!      data iset/0/

      IF (iset.eq.0) THEN
  100   v1 = 2.*ran3(idum)-1.
        v2 = 2.*ran3(idum)-1.
        r  = v1**2+v2**2
        IF(r.ge.1.) GOTO 100
        fac    = sqrt(-2.*log(r)/r)
        gset   = v1*fac
        gasdev = v2*fac
        iset   = 1
      ELSE
        gasdev = gset
        iset   = 0
      ENDIF

      RETURN
      END function


!#####################################################################
       SUBROUTINE output(y,time,kk,k1)
!#####################################################################

  USE hfrmod, ONLY : canfi, canps, canth, clamda, cpsie, een00, energ, hdid, imt, &
    iplot1, kount, loop, psi, r0, ttim1, tx1, ty1, tz1, xxx, yyp, zzz
  USE hfrmod, ONLY : allocate_output1, allocate_output2,                          &
                     deallocate_output1, deallocate_output2
  IMPLICIT NONE
  INTEGER(4),INTENT(IN):: kk, k1
  REAL(8),INTENT(IN):: y(*), time
  INTEGER(4):: ik
  REAL(8):: cfi, cps, cth, een, pp, reen, rr, zz, x(3)

         GO TO (100,120,200,300),kk

!    orbit plot data

 100   CONTINUE

      CALL allocate_output1

      DO ik=1,kount
        cps=0.5*(yyp(2,ik)**2+yyp(4,ik)**2)
        cth=datan2(yyp(4,ik),yyp(2,ik))
        cfi=yyp(1,ik)
        x(1)=cps/psi(loop)
        x(2)=cth
        x(3)=cfi
        CALL rzpfun(x,rr,zz,pp)
!     for deposing particle position information into arries tx1 and
!     tz1 which will be used to plot orbit trajectory latter on.
        iplot1=iplot1+1
        IF(iplot1.gt.imt) GO TO 110
        ttim1(iplot1)=time
        tx1(iplot1)=rr*r0*100.
        ty1(iplot1)=pp
        tz1(iplot1)=zz*r0*100.
      END DO

 110  CONTINUE
      CALL deallocate_output1
      RETURN
 120  CONTINUE

 130  CONTINUE

      CALL deallocate_output1
      RETURN

  200  CONTINUE
!      compute (x,z) of particles from thier (theta,psi,fi) by
!      interpolation.
!      and
!  sampling data of canonical variabls and position variables for
!  printing output.

      CALL allocate_output2

      x(1)=canps(k1)/psi(loop)
      x(2)=canth(k1)
      x(3)=canfi(k1)
      CALL rzpfun(x,rr,zz,pp)
      xxx(k1)=rr*r0*100.0
      zzz(k1)=zz*r0*100.0

      een=een00-energ(k1)
      reen=een/een00
      WRITE(7,704) time,y(1),y(2),clamda(k1),y(3), &
                   canps(k1)/cpsie,rr,zz,hdid,reen,energ(k1)
704    FORMAT(1p10e11.2,e14.5)

      CALL deallocate_output2
      RETURN

300   CONTINUE

      RETURN
      END SUBROUTINE


!#####################################################################
       SUBROUTINE dimles
!#####################################################################
!  the purpose of this subroutine is to make all field quantities
!  dimensionless and compute some constants needed.

  USE hfrmod, ONLY : ai, ambip, amf, b0, bnm, cpsie, dbr, dbr2, gg, loop, modmax,  &
    omega0, psi, qmasb, qmasp, qmc, r0, rnm, zb, znm
  IMPLICIT NONE
  INTEGER(4):: i, j, k
  REAL(8):: chq=1.6021d-19, qmp=1.67252d-27


      WRITE(6,600) r0,b0
 600  FORMAT(1h //,'   r0,b0=',1p2e14.4)
      qmasp=amf*qmp
      dbr=b0*r0
      dbr2=dbr*r0
      qmc=qmasp/(zb*chq)
      omega0=b0/qmc
      DO k=1,modmax
        DO i=0,loop
          bnm(k,i)=bnm(k,i)*b0
        END DO
      END DO
      DO i=0,loop
        ai(i)=ai(i)*b0
        gg(i)=gg(i)*b0
        psi(i)=psi(i)*b0
!c    psip(i)=psip(i)*b0
      END DO
      DO k=1,modmax
        DO i=0,loop
          rnm(k,i)=rnm(k,i)/r0
          znm(k,i)=znm(k,i)/r0
        END DO
      END DO
      DO j=0,loop
        psi(j)=psi(j)/dbr2
        ai(j)=ai(j)/dbr
        gg(j)=gg(j)/dbr
      END DO
      ambip=(ambip*1.6021d-19)/(qmasb*(omega0*r0)**2)
      cpsie=psi(loop)
      RETURN
      END SUBROUTINE


!#####################################################################
      SUBROUTINE roumup(y6,k1)
!#####################################################################

  USE hfrmod, ONLY : agg, aii, amu, canpfi, canpth, clamda, rou, sb, vpa
  IMPLICIT NONE
  INTEGER(4),INTENT(IN):: k1
  REAL(8),INTENT(IN):: y6
  REAL(8):: psip

      psip=agg*rou(k1)-canpfi(k1)
      rou(k1)=vpa(k1)*clamda(k1)/sb
!  amu(*) is magnetic moment.
      amu(k1)=vpa(k1)**2*(1.-clamda(k1)**2)/(2.*sb)
      canpfi(k1)=agg*rou(k1)-psip
      canpth(k1)=aii*rou(k1)+y6
!x    energ(k1)=0.5*(rou(k1)*sb)**2+amu(k1)*sb
      RETURN
      END SUBROUTINE


!#####################################################################
      SUBROUTINE dspln(y,x,n,y0,x0,n0,ww,iw)
!#####################################################################

  IMPLICIT NONE
  INTEGER(4),INTENT(IN):: n, n0, iw
  REAL(8),INTENT(IN) :: x(*), y0(*), x0(*)
  REAL(8),INTENT(OUT):: y(*)
  REAL(8),INTENT(INOUT):: ww(4,*)
  INTEGER(4):: i, ii, i0, ibcbeg, ibcend
  REAL(8):: xx, z

!    dimension y(*),x(*),y0(*),x0(*),ww(4,*)

      ibcbeg=0
      ibcend=0
      DO i=1,n0
        ww(1,i)=y0(i)
      END DO
      CALL cubspl(x0,ww,n0,ibcbeg,ibcend)
      DO i=1,n0
        ww(3,i)=ww(3,i)/2.0
        ww(4,i)=ww(4,i)/6.d0
      END DO

      IF(n.eq.0) RETURN

      DO i=1,n
        xx=x(i)
        DO i0=1,n0-1
          ii=i0
          IF(x(i).ge.x0(i0).and.x(i).lt.x0(i0+1)) GO TO 220
        END DO
        IF(x(i).ge.x0(n0)) ii=n0
        IF(x(i).le.x0(1)) ii=1
 220  CONTINUE
        z=x(i)-x0(ii)
        y(i)=ww(1,ii)+z*(ww(2,ii)+z*(ww(3,ii) +z*(ww(4,ii))))
      END DO
      RETURN
      END SUBROUTINE


!#####################################################################
      SUBROUTINE dspln1(y,x,n,y0,x0,n0,w,i)
!#####################################################################

  IMPLICIT NONE
  INTEGER(4),INTENT(IN):: n, n0, i
  REAL(8),INTENT(IN) :: x(*), y0(*), x0(*)
  REAL(8),INTENT(OUT):: y(*),w(4,*)
  INTEGER(4),PARAMETER :: ndim=202
  INTEGER(4) :: j, kp, m, n1, n2, ip(ndim)
  REAL(8) :: xmin, xmax, xx, z, dx(ndim),dy(ndim),y1(ndim)

 999  FORMAT(//' $$$$$$spln$$$$$ n0 (=',i4,') is larger than', &
     &       ' allowed dimension size (=',i4,')')

      IF(n0.gt.ndim) THEN
        PRINT 999,n0,ndim
        STOP                              999
      ELSE IF(i.eq.0) THEN
        n2= n0-2
        n1= n0-1
        DO m=1,n1
          dx(m)=x0(m+1)-x0(m)
          dy(m)=(y0(m+1)-y0(m))/dx(m)
        END DO
        y1(1)=(-dx(1)*dy(2)+(dx(2)+2.*dx(1))*dy(1))/(dx(2)+dx(1))
        DO m= 2,n1,1
          y1(m)=(dx(m)*dy(m-1)+dx(m-1)*dy(m))/(dx(m-1)+dx(m))
        END DO
        y1(n0)=((2*dx(n1)+dx(n2))*dy(n1)-dx(n1)*dy(n2))/(dx(n1)+dx(n2))
        DO m=1,n1,1
          w(1,m)=y0(m)
          w(2,m)=y1(m)
          w(3,m)=(3.*dy(m)-y1(m+1)-2.*y1(m))/dx(m)
          w(4,m)=(-2.*dy(m)+y1(m+1)+y1(m))/(dx(m)*dx(m))
        END DO
        w(1,n0)=y0(n0)
        w(2,n0)=0.d0
        w(3,n0)=0.d0
        w(4,n0)=0.d0
      ENDIF
      IF(n.ne.0) THEN
        IF(x0(n0).gt.x0(1)) THEN
          j=2
          DO m=1,n
            xx=x(m)
            DO while(xx.ge.x0(j).and.j.lt.n0)
              j=j+1
            END DO
            ip(m)=j
          END DO
          xmin=x0(1)
          xmax=x0(n0)
          DO m=1,n
            kp=min(ip(m),n0)-1
            xx=min(x(m),xmax)
            xx=max(xx,xmin)
            z=max(xx-x0(kp),0.d0)
            y(m)=w(1,kp)+z*(w(2,kp)+z*(w(3,kp)+z*w(4,kp)))
          END DO
        ELSE
          j=n0
          DO m=n,1,-1
            xx=x(m)
            DO while(xx.ge.x0(j).and.j.gt.2)
              j=j-1
            END DO
            ip(m)=j
          END DO
          xmax=x0(1)
          xmin=x0(n0)
          DO m=1,n
            kp=min(ip(m),n0)-1
            xx=min(x(m),xmax)
            xx=max(xx,xmin)
            z=xx-x0(kp)
            y(m)=w(1,kp)+z*(w(2,kp)+z*(w(3,kp)+z*w(4,kp)))
          END DO
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE

!#####################################################################
      SUBROUTINE bounce(iflg,kk1)
!#####################################################################

  USE hfrmod, ONLY : agg, ambip, amu, b0, canfi, canpfi, canps, canth, clamda, cmax, &
    cpsie, dplot, dpsii, dt, dt0, dtime, dtx, dxout, dxsav, een0, energ, eps, hdid,  &
    ipartp, iplot1, irunge, k1, kmax, lxy, mmf, nmax, nn, nplotp, nprntp, qmc, r0,   &
    rou, sb, tlim, vpa
  IMPLICIT NONE
  INTEGER(4),INTENT(IN) :: kk1
  INTEGER(4),INTENT(OUT):: iflg
  INTEGER(4):: i, iorb, nbad, nok, nst
  REAL(8):: h1, hmin, hnext, smamu, smeng, smlmd, smphi, smpsi, smroh, smtht, timeob, &
    x1, x2, xout, xplot, zdt, zsum, y(nmax), ymax(nmax)

!ppp
!cc      print * ,' bounce is called'

      k1=kk1

      dtx=dt0
      cmax=cpsie
      DO i=1,mmf
        dtime(i)=0.0d0
      END DO

      kmax=nplotp/nprntp
      dxsav=tlim/nplotp
      dxout=tlim/nprntp
      dplot=tlim/nplotp
!ppp
!cc      print * ,' kmax is calculated.'

      DO i=1,nmax
        ymax(i)=1.d0
      END DO
!!         iplot0=0       !200710
      iplot1=0
      timeob=0.d0
!ppp
!cc      print * ,' timeob is calculated.'

      CALL initc(k1)
!ppp
!      print * ,' initc is ended'

!      write(6,600) k1,psis/cpsie,clamda(k1),canth(k1)
 600  FORMAT(1h //,' *** start, k1 =',i10,' psi =',1pe14.4, '  lamda=', e14.4,' canth=',e14.4)
      IF(mod(k1,ipartp).eq.0) WRITE(7,703) k1
 703  FORMAT(1h ///,'     particle number =',i5,//, &
             7x,'time',7x,'y(1)',7x,'y(2)',7x,'clam',7x,'y(3)',8x, &
             'psi', 10x,'x',10x,'z',9x,'dt',4x,'denergy',6x,'energy'/)
!!!         theta0=canth(k1)    !! 200710
      xout=dxout
      xplot=dplot
!
 111  CALL derivb(canps(k1),canth(k1),canfi(k1))
!ppp
!cc      print * ,' parab is ended'
      CALL parab(canps(k1))
!ppp
!cc      print * ,' parab is ended'

!     compute the initial rou(parallel),magnetic moment, momentum and
!     energy of particles.
      rou(k1)=vpa(k1)*clamda(k1)/sb
      canpfi(k1)=agg*rou(k1)
      CALL roumup(canps(k1),k1)
!ppp
!      print * ,' roumup is ended'

      energ(k1)=0.5*(rou(k1)*sb)**2+amu(k1)*sb +ambip*(1.-canps(k1)/cpsie)
      een0=energ(k1)
!!          engy00=een0    !! 200710


!x    write(7,7000) k1,y(1),y(2),y(3),y(4),y(5),y(6)
 7000 FORMAT('  k1 =',i5,1p6e11.2)
      hnext=dt0
      hmin=dt0/1.e5
      h1=dt/32.0
!!         kflag=1  !! 200710
!
      smpsi = canps(k1)/cpsie
      smtht = canth(k1)
      smphi = canfi(k1)
      smroh = rou(k1)*r0
!     smroh = rou(k1)*0.957909841d+8*r0
!     smamu = amu(k1)*0.957909841d+8*r0**2*b0
      smamu = amu(k1)*r0**2*b0/(qmc)
!     smeng = energ(k1)*0.957909841d+8*r0**2*b0**2
      smeng = energ(k1)*r0**2*b0**2/(qmc)

      smlmd = clamda(k1)
!
!        write(9,7001) k1,canps(k1),canth(k1),canfi(k1),
!    &        rou(k1),amu(k1),energ(k1),sb,cpsie
         WRITE(9,7001) k1,smpsi,smtht,smphi,smroh,smamu,smeng,smlmd
 7001    FORMAT(i4,8d24.8)
 7002    FORMAT(i4,7d24.6)

!ppppppppp
      RETURN


!     orbit
!
!     call runge-kutta code:
 200  CONTINUE
!
!c    if(canps(k1).le.psibm*cpsie) then

      lxy=.true.
      y(1)=canfi(k1)
      y(2)=dsqrt(2.*canps(k1))*dcos(canth(k1))
      y(3)=rou(k1)
      y(4)=dsqrt(2.*canps(k1))*dsin(canth(k1))

      dt=hnext
      dt=dmin1(dt,dt0/(50.*((canps(k1)/cpsie)**4+0.02)))

 6001 FORMAT(1h //,' too small timestep,canps(k1),y(3),dt  =', i5,1p3e14.4,'  quit'//)

      x1=timeob
      x2=timeob+dxout
      h1=dmin1(dt0,dt)
      iflg=1
      CALL odeint(y,4,x1,x2,eps,h1,hdid,hmin,nok,nbad, iflg,nst,irunge)
      timeob=x2
      dt=h1
      hnext=h1
!!      engy00=engy11   ! 200710

      canfi(k1)=y(1)
      rou(k1)=y(3)

      IF(.not.lxy) THEN
        canth(k1)=y(2)
        canps(k1)=y(4)
      ELSE
        canth(k1)=datan2(y(4),y(2))
        canps(k1)=0.5*(y(2)**2+y(4)**2)
      ENDIF
!
      CALL derivb(canps(k1),canth(k1),canfi(k1))
      energ(k1)=0.5*(y(3)*sb)**2+amu(k1)*sb  +ambip*(1.-canps(k1)/cpsie)

      IF(iflg.le.0) GO TO 1000

      zdt=dt
      IF(mod(k1,ipartp).eq.0) THEN
        CALL output(y,timeob,1,k1)
      ENDIF
      IF(mod(k1,ipartp).eq.0) THEN
      CALL output(y,timeob,3,k1)
!  sampling data of canonical variabls and position variables for
!  printing output.
      xout=xout+dxout
      ENDIF

      iorb=canps(k1)/cpsie*dpsii+1.1
      dtime(iorb)=dtime(iorb)+zdt
!x     if(dabs(y(4)-theta0).ge.4.*cpi) then
!x     dtime(iorb)=dtime(iorb)-zdt*dabs(y(4)-theta0)/dabs(y(4)-y40)

!x     goto 1000
!x     endif

      IF(timeob.ge.tlim) GO TO 1000
      GO TO 200
1000  CONTINUE
      IF(iflg.gt.0) THEN
      zsum=0.0
      DO 300 i=1,mmf
        zsum=zsum+dtime(i)
 300  CONTINUE
      DO 310 i=1,mmf
        dtime(i)=dtime(i)/zsum
 310  CONTINUE
      ENDIF

      IF(mod(k1,ipartp).eq.0) CALL output(y,timeob,4,k1)
!      write(7,601) k1,iflg,nst,timeob,hdid,canps(k1)/cpsie,energ(k1)
601   FORMAT(1h //,'** end**  k1,iflg,nst,time,hdid,psi,energ**' ,3i5,1p4e13.3)

10000 CONTINUE

       RETURN
       END SUBROUTINE


!#####################################################################
       SUBROUTINE initc(k1)
!#####################################################################
  USE hfrmod, ONLY :agg, aiott, ambip, canfi, canps, canth, clamda, cpsie, dbr, een00,&
    iflag, loop, phi, phis, psi, psis, qmc, r0, rr0, sb, thets, vpa, vx00, vy00, vz00
  IMPLICIT NONE
  INTEGER(4),INTENT(IN):: k1
  REAL(8) :: dvp, dvp0, dvr, dvr0, dvz, dvz0, psin, rrr0, vph, vr, vsq, zdv, zvp, &
    x(3), fvec(3), fjac(3,3)
  INTEGER(4):: nfj=3, ldfjac=3

!ppp
!ccc      print * ,' initc is started.'
      psis=psis*psi(loop)
      psin=psis
      rrr0=rr0/r0
      vsq=sqrt(vx00**2+vy00**2+vz00**2)
      vr=vx00*cos(phi)+vy00*sin(phi)
      vph=-vx00*sin(phi)+vy00*cos(phi)

      CALL parab(psin)
!ppp
!cc      print * ,' parab is end.'
      CALL intarz(psin)
!ppp
!cc      print * ,' intarz is end.'
      CALL derivb(psin,thets,phis)
!ppp
!ccc      print * ,' derivb is end.'

!!      psn=psis/psi(loop)   ! 200710
!!      psil=psi(loop)       ! 200710
      iflag=2
      x(1)=psis
      x(2)=thets
      x(3)=phis

      CALL fcn(nfj,x,fvec,fjac,ldfjac,iflag)
!ppp
!cc      print * ,' fcn is end.'
!
      dvr0=(fjac(1,3)+aiott*fjac(1,2))*sb/abs(agg)
      dvz0=(fjac(2,3)+aiott*fjac(2,2))*sb/abs(agg)
      dvp0=(fjac(3,3)+aiott*fjac(3,2))*sb/abs(agg)*rrr0

!x     call dnorm(rrr0,psn,thets,phis,psil,dvr,dvz,dvp,sgg)

      zdv=dsqrt(dvr0**2+dvz0**2+dvp0**2)
      dvr=dvr0/zdv
      dvz=dvz0/zdv
      dvp=dvp0/zdv
      zvp=(vr*dvr+vz00*dvz+vph*dvp)
      vpa(k1)=vsq*1.d-2*qmc/dbr

!  here vpa(k1) is dimensionless initial velosity.

      canps(k1)=psis
      canth(k1)=thets
      clamda(k1)=zvp/vsq
      canfi(k1)=phis
      een00=0.5*vpa(k1)**2+ambip*(1.-canps(k1)/cpsie)


!ppp
!      write(6,600) k1,psis/cpsie,vsq,clamda(k1)
!      write(6,611) zdv,dvr0,dvz0,dvp0

 600  FORMAT(1h ,' ** sub.initc.,k1,psis,vsq,clamda',i5,1p3e11.2)
 611  FORMAT(' *** zdv,dvr,dvz,dvp**',1p4e12.3)

      RETURN
      END SUBROUTINE


!#####################################################################
      SUBROUTINE cubspl( tau, c, n, ibcbeg, ibcend )
!#####################################################################

!     ************************* input **************************
!     n = number of data points. assumed to be .ge. 2.
!     (tau(i), c(1,i), i=1,...,n) = absissae and ordinates of the
!        data points. tau is assumed to be strictly increasing.
!     ibcbeg, ibcend = boundary condition indicators, and
!     c(2,1), c(2,n) = boundary condition information. specifically,
!        ibcbeg = 0  means no boundary condition at tau(1) is given.
!           in this case, the not-a-knot condition is used, i.e. the
!           jump in the third derivative acroos tau(2) is forced to
!           zero, thus the first and the second cubic polynomial pieces
!           are made to coincide.)
!        ibcbeg = 1  means that the slope at tau(1) is made to equal
!           c(2,1), supplied by input.
!        ibcbeg = 2  means that the second derivative at tau(1) is
!           made to equal c(2,1) supplied by input.
!        ibcend = 0, 1, or 2 has analogous meaning concerning the
!           boundary condition at tau(n), with the additional infor-
!           mation taken from c(2,n).
!     ***********************   output  *************************
!     c(j,i), j=1,...,4; i=1,...,l (= n-1) = the polynomial coefficients
!        of the cubic interpolation spline with interior knots (or
!        joints) tau(2), ..., tau(n-1). precisely, in the interval
!        interval (tau(i), tau(i+1)), the spline f is given by
!           f(x) = c(1,i)+h*(c(2,i)+h*(c(3,i)+h*c(4,i)/3.)/2.)
!        where h = x - tau(i).

  IMPLICIT NONE
  INTEGER(4),INTENT(IN):: n, ibcbeg, ibcend
  REAL(8),INTENT(IN)   :: tau(*)
  REAL(8),INTENT(INOUT):: c(4,*)
  INTEGER(4):: i, j, l, m
  REAL(8):: divdf1, divdf3, dtau, g

!      dimension c(4,*),tau(*)
!******* a tridiagonal linear system for the unknown slopes s(i) of
!  f  at tau(i), i=1,...,n, is generated and then solved by gauss elim-
!  ination, with s(i) ending up in c(2,i), all i.
!     c(3,.) and c(4,.) are used initially for temporary storage.
      l=n-1
!ompute first differences of tau sequence and store in c(3,.). also,
!ompute first duvided difference of data and store in c(4,.).
      DO m=2,n
        c(3,m)=tau(m) - tau(m-1)
        c(4,m)=(c(1,m)-c(1,m-1))/c(3,m)
      END DO
!onstructfirst equation from the boundary condition, of the form
!     c(4,1)*s(1) + c(3,1)*s(2) = c(2,1)
      IF(ibcbeg-1)                      11,15,16
 11   IF(n .gt. 2)                      GO TO 12
!     no condition at left end n = 2.
      c(4,1) = 1.
      c(3,1) = 1.
      c(2,1)= 2.*c(4,2)
      GO TO 25
!     not-a-knot condition at left end and n .gt. 2.
 12   c(4,1) = c(3,3)
      c(3,1) = c(3,2) + c(3,3)
      c(2,1) =((c(3,2)+2.*c(3,1))*c(4,2)*c(3,3)+c(3,2)**2*c(4,3))/c(3,1)
      GO TO 19
!     slope prescribed at left end.
 15   c(4,1)= 1.
      c(3,1)= 0.
      GO TO 18
!     second derivative prescribed at left end.
 16   c(4,1)= 2.
      c(3,1)= 1.
      c(2,1)= 3.d0*c(4,2)-c(3,2)/2.*c(2,1)
 18   IF(n .eq. 2)                      GO TO 25
!     if there are interior knots, generate the corresp. equations and car-
!     ry out the forward pass of gauss elimination, after which the m-th
!     equations reads   c(4,m)*s(m) + c(3,m)*s(m+1) = c(2,m).
 19   DO m=2,l
         g = -c(3,m+1)/c(4,m-1)
         c(2,m) = g*c(2,m-1)+3.d0*(c(3,m)*c(4,m+1)+c(3,m+1)*c(4,m))
       c(4,m) = g*c(3,m-1)+2.d0*(c(3,m)+c(3,m+1))
      END DO
!onstruct last equation from the second boundary condition, of the form
!     (-g*c(4,n-1))*s(n-1) + c(4,n)*s(n) = c(2,n)
!     if slope is prescribed at right end, one can go directly to back-
!     substitution, since c array happens to be set uo just right for it
!     at this point.
      IF (ibcend-1)                     21,30,24
 21   IF (n .eq. 3 .and. ibcbeg .eq. 0) GO TO 22
!     not-a-knot and n .ge.3, and either n.gt.3 or  also not-a-knot at
!     left end point.
      g = c(3,n-1) + c(3,n)
      c(2,n) = ((c(3,n)+2.*g)*c(4,n)*c(3,n-1) &
     &     + c(3,n)**2*(c(1,n-1)-c(1,n-2))/c(3,n-1))/g
      g = -g/c(4,n-1)
      c(4,n) = c(3,n-1)
      GO TO 29
!     either (n=3 and not-a-knot also at left) or (n=2 and not-a-
!     knot at left end point).
 22   c(2,n) = 2.*c(4,n)
      c(4,n) = 1.
      GO TO 28
!     second derivative prescribed at right end point.
 24   c(2,n) = 3.d0*c(4,n) + c(3,n)/2.*c(2,n)
      c(4,n) = 2.
      GO TO 28
!
 25   IF(ibcend-1)                      26,30,24
 26   IF(ibcbeg .gt. 0)                 GO TO 22
!     not-a-knot at right endpoint and at left endpoint and n = 2.
      c(2,n)= c(4,n)
      GO TO 30

  28  g=-1./c(4,n-1)
!ompleteforward pass of gauss elimination.
 29   c(4,n)= g*c(3,n-1) + c(4,n)
      c(2,n)= (g*c(2,n-1) + c(2,n))/c(4,n)
!arry back substitution
 30   DO j=l,1,-1
         c(2,j)= (c(2,j)-c(3,j)*c(2,j+1))/c(4,j)
      END DO
!******generate cubic coefficients in each interval, i.e., the deriv.s
!     at its left endpoint, from value at its endpoints.
      DO i=2,n
        dtau= c(3,i)
        divdf1= (c(1,i)-c(1,i-1))/dtau
        divdf3= c(2,i-1) + c(2,i)-2.*divdf1
        c(3,i-1)= 2.*(divdf1 - c(2,i-1) -divdf3)/dtau
        c(4,i-1)= (divdf3/dtau)*(6.d0/dtau)
      END DO
      RETURN
      END SUBROUTINE


!!##################################################################
!      REAL(8) function ran3(idum)
!!##################################################################
!
!  IMPLICIT NONE
!  INTEGER(4),INTENT(INOUT):: idum
!  INTEGER(4),PARAMETER:: mseed=161803398, mbig=1000000000, mz=0
!  REAL(8),PARAMETER:: fac=1.d-9
!  INTEGER(4):: i, ii, inext, inextp, k, mj, mk, ma(55)
!  INTEGER(4):: iff=0
!!!!     parameter (mbig=4000000.,mseed=1618033.,mz=0.,fac=2.5e-7)
!!      parameter (mbig=1000000000,mseed=161803398,mz=0,fac=1.e-9)
!!      dimension ma(55)
!!      data iff /0/
!
!
!      IF(idum.lt.0.or.iff.eq.0)THEN
!        iff=1
!        mj=mseed-iabs(idum)
!        mj=mod(mj,mbig)
!        ma(55)=mj
!        mk=1
!        DO i=1,54
!          ii=mod(21*i,55)
!          ma(ii)=mk
!          mk=mj-mk
!          IF(mk.lt.mz)mk=mk+mbig
!          mj=ma(ii)
!        END DO
!        DO k=1,4
!          DO i=1,55
!            ma(i)=ma(i)-ma(1+mod(i+30,55))
!            IF(ma(i).lt.mz)ma(i)=ma(i)+mbig
!          END DO
!        END DO
!        inext=0
!        inextp=31
!        idum=1
!      ENDIF
!      inext=inext+1
!      IF(inext.eq.56)inext=1
!      inextp=inextp+1
!      IF(inextp.eq.56)inextp=1
!      mj=ma(inext)-ma(inextp)
!      IF(mj.lt.mz)mj=mj+mbig
!      ma(inext)=mj
!      ran3=mj*fac
!      RETURN
!      END FUNCTION
!
!
!#####################################################################
        SUBROUTINE ppout
!#####################################################################

!.................................................................

! output of positions (r,z) of fast ions hitting the first wall

!.................................................................

  USE hfrmod , ONLY : noplos, prpos, pzpos
  IMPLICIT NONE
  INTEGER(4):: n
  REAL(8) :: area, rmax, rmin, zmax, zmin

      IF(noplos.lt.1) THEN
        write(6,*) 'from ppout: noplos=', noplos
      WRITE(6,900)
      RETURN
      ENDIF


      WRITE(6,1000)
      DO n=1,noplos
        WRITE(6,1010) n,prpos(n),pzpos(n)
      END DO


      rmax=prpos(1)
      DO n=1,noplos
        rmax=dmax1(prpos(n),rmax)
      END DO
      rmin=prpos(1)
      DO n=1,noplos
        rmin=dmin1(prpos(n),rmin)
      END DO
      zmax=pzpos(1)
      DO n=1,noplos
        zmax=dmax1(pzpos(n),zmax)
      END DO
      zmin=pzpos(1)
      DO n=1,noplos
        zmin=dmin1(pzpos(n),zmin)
      END DO

      area=(rmax-rmin)*(zmax-zmin)

      WRITE(6,1020) rmax,rmin,zmax,zmin,area

!
      RETURN
!
!....................................................................
!         FORMAT statements

 900  FORMAT(1h1,///,'*****************************************', &
     &             /,'** no fast ions hitting the first wall **', &
     &             /,'*****************************************')
 1000 FORMAT(1h1,'position of fast ions hitting the first wall', &
     &     //1h ,'  no.       r(cm)           z(cm)'/)
 1010 FORMAT(i7,1p2e14.4)
 1020 FORMAT(//1h ,'area hitted by fast ions (cm2)',//, &
     &     1h ,10x, '     rmax =',1pe12.4,'     rmin =',e12.4,/, &
     &     1h ,10x, '     zmax =',1pe12.4,'     zmin =',e12.4,//, &
     &     1h ,10x, 'area(cm2) =',1pe12.4)
!....................................................................

      END SUBROUTINE

!!#####################################################################
!      SUBROUTINE fcnv(x,n,fvec)
!!#####################################################################
!
!  USE hfrmod, ONLY : rr, zz, pphi
!  IMPLICIT NONE
!  INTEGER(4),INTENT(IN) :: n
!  REAL(8),INTENT(INOUT) :: x(n)
!  REAL(8),INTENT(OUT):: fvec(n)
!  REAL(8):: r,z,p
!!      include 'include/field'
!!      common /cpsino/ psino(0:100),xmin
!!      common /comsol/ rr,zz,pphi
!!      dimension x(n),fvec(n)
!
!      IF(x(1).lt.0.0) THEN
!         x(1)=-x(1)
!         x(2)=x(2)+3.1415926535897
!      ENDIF
!!
!!     x1=x(1)
!!     x2=x(2)
!!     x3=x(3)
!
!      CALL rzpfun(x,r,z,p)
!
!      fvec(1)=r-rr
!      fvec(2)=z-zz
!      fvec(3)=p-pphi
!
!      RETURN
!      END SUBROUTINE
!
!
!!#####################################################################
!      SUBROUTINE fcnj(x,n,fjac)
!!#####################################################################
!
!  USE hfrmod, ONLY : abmnum, cp, cr, cz, dp, dr, dz, ep, er, ez, loop, mnumbr,    &
!                      modmax, nnumbr, pnmn, psino, rnmn, znmn
!  IMPLICIT NONE
!  INTEGER(4),INTENT(IN):: n
!  REAL(8),INTENT(INOUT):: x(n)
!  REAL(8),INTENT(OUT)  :: fjac(n,n)
!  INTEGER(4):: i, ii, j, k, ldfjac
!  REAL(8):: cosf, dpharm, drharm, dxabm, dzharm, pharm, rharm, sinf, x1, x2, x3,   &
!            xabm, xx, zharm
!!      include 'include/field'
!!      common /cpsino/ psino(0:100),xmin
!!      common /comsol/ rr,zz,pphi
!!      dimension x(n)
!!      dimension fjac(n,n)
!
!      ldfjac = n
!
!      IF(x(1).lt.0.0) THEN
!         x(1)=-x(1)
!         x(2)=x(2)+3.1415926535897
!      ENDIF
!
!      x1=x(1)
!      x2=x(2)
!      x3=x(3)
!
!      IF(x1.gt.1.0) GO TO 500
!      IF(x1.eq.0.0) THEN
!         ii=0
!         xx=0.
!         GO TO 50
!      ENDIF
!
!      DO i=1,loop
!         IF(x1.gt.psino(i-1).and.x1.le.psino(i)) THEN
!            ii=i-1
!            xx=x1-psino(i-1)
!            GO TO 50
!         ENDIF
!      END DO
!
! 50   CONTINUE
!
!      DO i=1,ldfjac
!         DO j=1,n
!            fjac(i,j)=0.
!         END DO
!      END DO
!
!      DO k=1,modmax
!         xabm=x1**abmnum(k)
!         dxabm=x1**(abmnum(k)-1)*abmnum(k)
!         cosf=dcos(mnumbr(k)*x2-nnumbr(k)*x3)
!         sinf=dsin(mnumbr(k)*x2-nnumbr(k)*x3)
!         rharm=rnmn(k,ii)+(cr(k,ii)+(dr(k,ii)+er(k,ii)*xx)*xx)*xx
!         zharm=znmn(k,ii)+(cz(k,ii)+(dz(k,ii)+ez(k,ii)*xx)*xx)*xx
!         pharm=pnmn(k,ii)+(cp(k,ii)+(dp(k,ii)+ep(k,ii)*xx)*xx)*xx
!         drharm=cr(k,ii)+(2.d0*dr(k,ii)+3.d0*er(k,ii)*xx)*xx
!         dzharm=cz(k,ii)+(2.d0*dz(k,ii)+3.d0*ez(k,ii)*xx)*xx
!         dpharm=cp(k,ii)+(2.d0*dp(k,ii)+3.d0*ep(k,ii)*xx)*xx
!         fjac(1,1)=fjac(1,1)+(drharm*xabm+rharm*dxabm)*cosf
!         fjac(2,1)=fjac(2,1)+(dzharm*xabm+zharm*dxabm)*sinf
!         fjac(3,1)=fjac(3,1)+(dpharm*xabm+pharm*dxabm)*sinf
!         fjac(1,2)=fjac(1,2)-mnumbr(k)*rharm*xabm*sinf
!         fjac(1,3)=fjac(1,3)+nnumbr(k)*rharm*xabm*sinf
!         fjac(2,2)=fjac(2,2)+mnumbr(k)*zharm*xabm*cosf
!         fjac(2,3)=fjac(2,3)-nnumbr(k)*zharm*xabm*cosf
!         fjac(3,2)=fjac(3,2)+mnumbr(k)*pharm*xabm*cosf
!         fjac(3,3)=fjac(3,3)-nnumbr(k)*pharm*xabm*cosf
!      END DO
!      fjac(3,3)=1.+fjac(3,3)
!      RETURN
!
! 500  CONTINUE
! 600  FORMAT(1h //,' **** x **** out of range****',1p3e14.4)
! 800  CONTINUE
!      DO i=1,ldfjac
!         DO j=1,n
!            fjac(i,j)=0.
!         END DO
!      END DO
!!
!      DO k=1,modmax
!         cosf=dcos(mnumbr(k)*x2-nnumbr(k)*x3)
!         sinf=dsin(mnumbr(k)*x2-nnumbr(k)*x3)
!         IF(k.eq.1) GO TO 910
!         fjac(1,1)=fjac(1,1)+rnmn(k,loop)*cosf
!         fjac(1,2)=fjac(1,2)-mnumbr(k)*rnmn(k,loop)*sinf*x1
!         fjac(1,3)=fjac(1,3)+nnumbr(k)*rnmn(k,loop)*sinf*x1
! 910     fjac(2,1)=fjac(2,1)+znmn(k,loop)*sinf
!         fjac(2,2)=fjac(2,2)+mnumbr(k)*znmn(k,loop)*cosf*x1
!         fjac(2,3)=fjac(2,3)-nnumbr(k)*znmn(k,loop)*cosf*x1
!         fjac(3,2)=fjac(3,2)+mnumbr(k)*pnmn(k,loop)*cosf
!         fjac(3,3)=fjac(3,3)-nnumbr(k)*pnmn(k,loop)*cosf
!      END DO
!
!      fjac(3,3)=fjac(3,3)+1.0
!      RETURN
!      END SUBROUTINE

END MODULE
