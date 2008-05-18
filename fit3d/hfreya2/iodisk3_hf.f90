
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  SUBROUTINE iodisk_hfr(ihfmax,ihfmn0,ihfnhl,ihloop,ihnnm,ihmnm,hfpsi &
            ,hfgg,hfai,hfaiot,hfvol,hfbbnm,hfb2nm,hfrnm,hfznm,hfpnm)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     programme for interface vmec(newboz)
!
!                                         v3.0
!                                         modified by S. Murakami
!
!                                         original version
!                                              by n.nakajima (93.12)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!     iodisk ; read input data from newboz
!
!     nsd          number of vmec grid points
!     mbzmax       number of boozer theta harmonics desired
!     nbzmin       number of negative boozer zeta harmonics desired
!     nbzmax       number of positive boozer zeta harmonics desired
!
!     bbozh : mod b spectrum in boozer coordinates
!
!           b(i,t,z) =     sum bbozh(m,i) cos(-t*mboz(m) + z*nboz(m) )
!
!     note   bco(i,m) = bbozh(m,i)
!
!     transformation from boozer to cylindrical coordinates
!
!         ( psi, theta, zeta ) -> ( r, phi, z )
!
!           r(i,t,z) =     sum rbozh(m,i) cos(-t*mboz(m) + z*nboz(m) )
!           p(i,t,z) = z + sum pbozh(m,i) sin(-t*mboz(m) + z*nboz(m) )
!           z(i,t,z) =     sum zbozh(m,i) sin(-t*mboz(m) + z*nboz(m) )
!
!     surface quantities ( after renormalized )
!
!           psi   : toroidal flux within a flux surface / ( 2*pi )
!           eot   : rotational transform
!           cui   : toroidal current within a flux surface
!           cug   : poloidal current without a flux surface

  USE hfrmod, ONLY : kmsh, mmx, ibhar, nbzmax, nbzmin, nsd, bb0, alpha, nmboz, mboz,    &
    nboz, psi2, nfp, eot, cui, cug,bbozh, rbozh, zbozh, pbozh, bco, cm, cn
  IMPLICIT NONE

  INTEGER(4),INTENT(OUT) :: ihfmax, ihfmn0, ihfnhl, ihloop, ihnnm(100), ihmnm(100)
  REAL(8)   ,INTENT(OUT) :: hfpsi(0:100), hfgg(0:100), hfai(0:100),  hfaiot(0:100),      &
    hfvol(0:100), hfbbnm(100,0:100), hfb2nm(100,0:100), hfrnm(100,0:100),                &
    hfznm(100,0:100), hfpnm(100,0:100)
  INTEGER(4), PARAMETER :: ndiskc = 10
  INTEGER(4):: i, igt, ihf, ii, iie, iis, ijf, ijf2, imttl, is, j, kk, m, nout, nprt, nsdm1
  INTEGER(4):: icont(ibhar), msel(ibhar,nsd)
  REAL(8) :: asel, bamax, cmu0, cnorm, critc, pi, pi2, rmaj0
  REAL(8) :: bmin(nsd), bmax(nsd), bsel(ibhar), psel(ibhar),rsel(ibhar),zsel(ibhar),     &
    pco(kmsh+1,mmx+1), rco(kmsh+1,mmx+1), zco(kmsh+1,mmx+1)
  LOGICAL :: lcheck = .false.
  NAMELIST /newbz/ rmaj0, bb0, alpha, critc

      rmaj0 = 3.9
        bb0 = 3.0
      alpha = 8.0
      nfp   = 10

      READ(5,newbz)
      WRITE(6,newbz)

!.... const

      pi     = dacos(-1.0d0)
      pi2    = pi*2.0d0
      cmu0   = 2.0d0*pi2*1.0d-7

!      READ(ndiskc,7001) nmboz, nsdm1, nfp
      READ(ndiskc) nmboz, nsdm1, nfp
 7001 FORMAT(3i4)
      PRINT 6001,  nmboz, nsdm1, nfp

      IF( nmboz - 1 .lt. mmx ) THEN
        PRINT 9001, nmboz - 1, mmx
        STOP
      ELSE
        PRINT 9003, nmboz - 1, mmx
      endif
      IF( nsdm1 .ne. kmsh    ) THEN
        PRINT 9002, nsdm1, kmsh
        STOP
      ELSE
        PRINT 9004, nsdm1, kmsh
      endif

!      do 1000, i = 2, nsd
!c         READ (ndiskc,7002) psi(i), eot(i), cui(i), cug(i)
!         READ (ndiskc) psi(i), eot(i), cui(i), cug(i)
! 1000 continue

      READ(ndiskc) (psi2(i), i = 2, nsd)
      READ(ndiskc) (eot(i), i = 2, nsd)

      READ(ndiskc) (cui(i), i = 2, nsd)
      READ(ndiskc) (cug(i), i = 2, nsd)


! 7002 FORMAT(4d25.17)

      DO m = 1, nmboz
        READ(ndiskc)  mboz(m), nboz(m)
        READ(ndiskc) (bbozh(m,i), i = 2, nsd)
      END DO

      DO m = 1, nmboz
        READ(ndiskc)  mboz(m), nboz(m)
        READ(ndiskc) (rbozh(m,i), i = 2, nsd)
      END DO

      DO m = 1, nmboz
        READ(ndiskc)  mboz(m), nboz(m)
        READ(ndiskc) (zbozh(m,i), i = 2, nsd)
      END DO

      DO m = 1, nmboz
        READ(ndiskc)  mboz(m), nboz(m)
        READ(ndiskc) (pbozh(m,i), i = 2, nsd)
      END DO

!$$$      do 101, m = 1, nmboz
!$$$c
!$$$c         READ (ndiskc,7003) mboz(m), nboz(m)
!$$$         READ (ndiskc) mboz(m), nboz(m)
!$$$         READ(NDISKC) (BBOZH(M,I), I = 2, NSD)
!$$$ 101 continue
!$$$c
!$$$c         do 1001, i = 2, nsd
!$$$c            READ (ndiskc,7004) bbozh(m,i), rbozh(m,i)
!$$$c     &                        ,zbozh(m,i), pbozh(m,i)
!$$$c            READ (ndiskc) bbozh(m,i), rbozh(m,i)
!$$$c     &                        ,zbozh(m,i), pbozh(m,i)
!$$$c 1002    continue

 7003 FORMAT(2i4)
 7004 FORMAT(4d25.17)

      PRINT 6002
      PRINT 6003, (i, psi2(i), eot(i), cui(i), cug(i), i = 2, nsd)

      IF( lcheck ) THEN
        nout  = 10
        nprt  = nsd/nout
        IF( mod(nsd,nout) .ne. 0 ) nprt = nprt + 1
        IF( mod(nsd,nout) .ne. 0 ) nprt = nprt + 1

        DO ii=1,nprt
          iis  =(ii-1)*nout+1
          iie  = ii   *nout
          IF( ii .eq. nprt ) iie  = nsd
          PRINT 6011,         ( kk , kk = iis, iie )

          DO j = 1, nmboz
            PRINT 6012, j, mboz(j), nboz(j), ( bbozh(j,kk), kk = iis, iie )
          END DO
        END DO
      ENDIF

!.... from half mesh to integer mesh

      DO i = 2, nsd-1
        psi2(i)   = ( psi2(i) + psi2(i+1) )/2.0d0
        eot(i)   = ( eot(i) + eot(i+1) )/2.0d0
        cui(i)   = ( cui(i) + cui(i+1) )/2.0d0
        cug(i)   = ( cug(i) + cug(i+1) )/2.0d0
        DO j = 1, nmboz
          bbozh(j,i) = ( bbozh(j,i) + bbozh(j,i+1) )/2.0d0
          rbozh(j,i) = ( rbozh(j,i) + rbozh(j,i+1) )/2.0d0
          zbozh(j,i) = ( zbozh(j,i) + zbozh(j,i+1) )/2.0d0
          pbozh(j,i) = ( pbozh(j,i) + pbozh(j,i+1) )/2.0d0
        END DO
      END DO

      PRINT *, ' - End 301 -'
!+++++ extrapolation of the values at the magnetic axis

      DO j = 1, nmboz
        IF( mboz(j) .eq. 0 ) THEN
          bbozh(j,1) = 3.0d0*bbozh(j,2) - 3.0d0*bbozh(j,3) + bbozh(j,4)
          rbozh(j,1) = 3.0d0*rbozh(j,2) - 3.0d0*rbozh(j,3) + rbozh(j,4)
          zbozh(j,1) = 3.0d0*zbozh(j,2) - 3.0d0*zbozh(j,3) + zbozh(j,4)
          pbozh(j,1) = 3.0d0*pbozh(j,2) - 3.0d0*pbozh(j,3) + pbozh(j,4)
        ELSE
          bbozh(j,1) = 0.0d0
          rbozh(j,1) = 0.0d0
          zbozh(j,1) = 0.0d0
          pbozh(j,1) = 0.0d0
        ENDIF
      END DO
      psi2 (  1) = 0.0d0
      eot  (  1) = 3.0d0*eot  (  2) - 3.0d0*eot  (  3) + eot  (  4)
      cui  (  1) = 3.0d0*cui  (  2) - 3.0d0*cui  (  3) + cui  (  4)
      cug  (  1) = 3.0d0*cug  (  2) - 3.0d0*cug  (  3) + cug  (  4)

      PRINT *, ' - End 302 -'
!+++++ extrapolation of the values at the boundary

      DO j = 1, nmboz
        bbozh(j,nsd) = 3.0d0*bbozh(j,nsd-1) - 3.0d0*bbozh(j,nsd-2) + bbozh(j,nsd-3)
        rbozh(j,nsd) = 3.0d0*rbozh(j,nsd-1) - 3.0d0*rbozh(j,nsd-2) + rbozh(j,nsd-3)
        zbozh(j,nsd) = 3.0d0*zbozh(j,nsd-1) - 3.0d0*zbozh(j,nsd-2) + zbozh(j,nsd-3)
        pbozh(j,nsd) = 3.0d0*pbozh(j,nsd-1) - 3.0d0*pbozh(j,nsd-2) + pbozh(j,nsd-3)
      END DO
!!      if( psic(nsd-1) .gt. 0.0d0 ) then
!!         psisgn   =  1.0d0
!!      else
!!         psisgn   = -1.0d0
!!      endif
!!!     psi  (nsd) = 0.5d0*psisgn
      psi2 (nsd) = 3.0d0*psi2(nsd-1) - 3.0d0*psi2(nsd-2) + psi2(nsd-3)
      eot  (nsd) = 3.0d0*eot(nsd-1)  - 3.0d0*eot(nsd-2)  + eot(nsd-3)
      cui  (nsd) = 3.0d0*cui(nsd-1)  - 3.0d0*cui(nsd-2)  + cui(nsd-3)
      cug  (nsd) = 3.0d0*cug(nsd-1)  - 3.0d0*cug(nsd-2)  + cug(nsd-3)
!
!$$$  c.... normalization to vmec calculation
!$$$  c
      cnorm  = rmaj0*bb0/cug(nsd)
!
!     cnorm = -1.d0
      DO i = 1, nsd
        cui(i)   = cui(i)*cnorm
        cug(i)   = cug(i)*cnorm
        psi2(i)  = psi2(i)*cnorm
        DO j = 1, nmboz
          bbozh(j,i) = bbozh(j,i)*cnorm
        END DO
      END DO
!$$$  c
!
      PRINT *, ' - End 303 -'
!..... assuming if = 1 correponds to mboz = 0 and nboz = 0.


!..... initialization of labels of fourier modes for each surface

      DO is = 1, nsd
        msel( 1,is) = 1
      END DO

      DO ijf = 2, nmboz

        bsel(ijf)=-1.0d0
        rsel(ijf)=-1.0d0
        zsel(ijf)=-1.0d0
        psel(ijf)=-1.0d0

        DO is = 1, nsd
          msel(ijf,is) = 0
        END DO
      END DO

!..... initialization of labels of foutier modes for all surfaces

!!        ifsum     = 0

      PRINT *, ' - End 502 -'

!      do 503 is = 1, nsd
!        ifsum     = ifsum + 1
!  503 continue

!        ifsel( 1) = ifsum

      DO ijf = 2, nmboz
!        ifsel(ijf) = 0
        icont(ijf) = 0
      END DO

      PRINT *, ' - End 504 -'
!..... selection of modes for each surface

      DO is = 1, nsd

        bmax(is)  = -10.0
        bmin(is)  =  10.0

        DO ijf = 2, nmboz
          bmax(is)  = dmax1( bbozh(ijf,is), bmax(is) )
          bmin(is)  = dmin1( bbozh(ijf,is), bmin(is) )
        END DO
!
        bamax = dmax1( dabs(bmax(is)), dabs(bmin(is)) )
!
        DO ijf = 2, nmboz
          asel = dabs(bbozh(ijf,is))/bamax
          bsel(ijf) = dmax1(bsel(ijf),asel)
        END DO
      END DO
!
      PRINT *, ' - End 601 -'
!.....selection of modes for all surfaces

      DO ijf = 2, nmboz
        igt = 0
        DO ijf2 = 2, nmboz
          IF (bsel(ijf2).gt.bsel(ijf)) THEN
             igt = igt+1
          END IF
        END DO
!
        IF (igt.lt.mmx) THEN
          icont(ijf) = 1
        END IF
      END DO

      imttl   = 1

      DO is = 1, nsd
        bco( is,1) = bbozh( 1,is)
        rco( is,1) = rbozh( 1,is)
        pco( is,1) = pbozh( 1,is)
        zco( is,1) = zbozh( 1,is)
      END DO

      PRINT *, ' - End 700 -'


      DO ijf = 2, nmboz
        IF (icont(ijf) .ne. 0 ) THEN
          imttl     = imttl + 1
          DO is = 1, nsd
            bco(is,imttl) = bbozh(ijf,is)
            rco(is,imttl) = rbozh(ijf,is)
            pco(is,imttl) = pbozh(ijf,is)
            zco(is,imttl) = zbozh(ijf,is)
          END DO
          cm (imttl)    = mboz(ijf)
          cn (imttl)    = nboz(ijf)
          ihmnm(imttl)  = mboz(ijf)
          ihnnm(imttl)  = nboz(ijf)
        END IF
      END DO


      PRINT 6031, critc, imttl

      DO is = 1, nsd
        PRINT 6032, is, icont(is), bmin(is), bmax(is)
      END DO

      PRINT *
      PRINT *,' selected foutier spectrum of b '
      PRINT *

      nout  = 10
      nprt  = nsd/nout
      IF( mod(nsd,nout) .ne. 0 ) nprt = nprt + 1

      DO ii=1,nprt
        iis  =(ii-1)*nout+1
        iie  = ii   *nout
        IF( ii .eq. nprt ) iie  = nsd
        PRINT 6033,         ( kk , kk = iis, iie )

        DO j = 1, imttl
          PRINT 6034, j, idint(cm(j)), idint(cn(j)), ( bco(kk,j), kk = iis, iie )
        END DO
      END DO



      PRINT 6004
      PRINT 6003, (i, psi2(i), eot(i), cui(i)*pi2/(1.0d6*cmu0), &
                   cug(i)*pi2/(1.0d6*cmu0), i = 1, nsd)

      nout  = 10
      nprt  = nsd/nout
      IF( mod(nsd,nout) .ne. 0 ) nprt = nprt + 1

      DO ii=1,nprt
        iis  =(ii-1)*nout+1
        iie  = ii   *nout
        IF( ii .eq. nprt ) iie  = nsd
        PRINT 6013,         ( kk , kk = iis, iie )

        DO j = 1, mmx + 1
          PRINT 6012, j, idint(cm(j)), idint(cn(j)), ( bco(kk,j), kk = iis, iie )
        END DO
      END DO

      IF( lcheck ) THEN

      DO ii=1,nprt
        iis  =(ii-1)*nout+1
        iie  = ii   *nout
        IF( ii .eq. nprt ) iie  = nsd
        PRINT 6015,         ( kk , kk = iis, iie )

        DO j = 1, mmx + 1
          PRINT 6014, j, idint(cm(j)), idint(cn(j)), ( rco(kk,j), kk = iis, iie )
        END DO
      END DO

        DO ii=1,nprt
          iis  =(ii-1)*nout+1
          iie  = ii   *nout
          IF( ii .eq. nprt ) iie  = nsd
          PRINT 6017,         ( kk , kk = iis, iie )

          DO j = 1, mmx + 1
            PRINT 6016, j, idint(cm(j)), idint(cn(j)), ( zco(kk,j), kk = iis, iie )
          END DO
        END DO

        DO ii=1,nprt
          iis  =(ii-1)*nout+1
          iie  = ii   *nout
          IF( ii .eq. nprt ) iie  = nsd
          PRINT 6019, ( kk , kk = iis, iie )

          DO j = 1, mmx + 1
            PRINT 6018, j, idint(cm(j)), idint(cn(j)), ( pco(kk,j), kk = iis, iie )
          END DO
        END DO

       END IF

!------------------------------------------c
!    data conversion for HFREYA
!------------------------------------------c

      ihfmax = mmx+1
      ihfmn0 = 1
      ihfnhl  =10
      ihloop = kmsh

      DO i=1,nsd
        ihf = i-1

        hfpsi(ihf)  = psi2(i)
        hfgg(ihf)   = cug(i)
        hfai(ihf)   = cui(i)
        hfaiot(ihf) = eot(i)
        hfvol(ihf)  = 4.d0*pi*pi*rmaj0*psi2(i)/bb0

        DO j=1, mmx+1
          hfbbnm(j,ihf) =  bco(i,j)
          hfb2nm(j,ihf) =  bco(i,j)**2
          hfrnm(j,ihf)  =  rco(i,j)
          hfznm(j,ihf)  = -zco(i,j)
          hfpnm(j,ihf)  = -pco(i,j)
        END DO
      END DO

!      CALL deallocate_iodisk

 6031 FORMAT(//10x,'critical value of spectrum = ',1pd12.5 &
              /10x,'the number of modes needed for all surfaces = ',i6 &
             //10x,'is',10x,'num. modes',10x,'bmin',10x,'bmax'/)
 6032 FORMAT(  10x,i2,15x,i5,2x,1pd12.5,2x,d12.5)
 6033 FORMAT(//2x,' j',2x,'(  m,  n)',10(5x,'b   _',i2)/)
 6034 FORMAT(  2x,i2,2x,'(',i3,',',i3,')',10(1x,1pd11.4))
 6021 FORMAT(//4x,'zeta /pi = ',1pd12.5 &
             //4x,'theta/pi',13x,'b',13x,'r',13x,'z',8x,'phi/pi'/)
 6022 FORMAT(5(2x,1pd12.5))

 6001 FORMAT(//9x,'nmboz  = ',i4,' : nsd - 1 = ',i4,' : nfp    =',i4//)
 6002 FORMAT(//4x,'i',11x,'psi',10x,'iota',13x,'i',13x,'g'/)
 6003 FORMAT(/(2x,i3,4(2x,1pd12.5)))
 6004 FORMAT(//2x,130('+')//10x,'dimensional physical quantities' &
             //4x,'i',3x,'psi(t m**2)', &
              10x,'iota',4x,'2*pi*i(ma)',4x,'2*pi*g(ma)'/)
 6011 FORMAT(//2x,' j',2x,'( m,  n)',10(5x,'b   _',i2)/)
 6013 FORMAT(//2x,' j',2x,'( m,  n)',10(5x,'b(t)_',i2)/)
 6012 FORMAT(  2x,i2,2x,'(',i2,',',i3,')',10(1x,1pd11.4))
 6015 FORMAT(//2x,' j',2x,'( m,  n)',10(5x,'R(t)_',i2)/)
 6014 FORMAT(  2x,i2,2x,'(',i2,',',i3,')',10(1x,1pd11.4))
 6017 FORMAT(//2x,' j',2x,'( m,  n)',10(5x,'Z(t)_',i2)/)
 6016 FORMAT(  2x,i2,2x,'(',i2,',',i3,')',10(1x,1pd11.4))
 6019 FORMAT(//2x,' j',2x,'( m,  n)',10(5x,'P(t)_',i2)/)
 6018 FORMAT(  2x,i2,2x,'(',i2,',',i3,')',10(1x,1pd11.4))
 9001 FORMAT(//10x,'nmboz - 1 = ',i4,' <  mmx  = ',i4//)
 9003 FORMAT(//10x,'nmboz - 1 = ',i4,' >= mmx  = ',i4//)
 9002 FORMAT(//10x,'nsd   - 1 = ',i4,' =/ kmsh = ',i4//)
 9004 FORMAT(//10x,'nsd   - 1 = ',i4,' =  kmsh = ',i4//)

      RETURN
      END
