!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE iodisk
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     programme for interface vmec(newboz)
!                                         v4.1
!                                         modified by S. Murakami
!                                         original version
!                                              by n.nakajima (93.12)
! 04.05.27
!  The B<0 case with a B>0 file is added.
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

  USE mcnmod, ONLY : alpha, bb0, bbozh, bco, cm, cn, cug, cui, eot, ibhar, ispsi,  &
    kmsh, mboz, mmn, mmx, nboz, nfp, nmboz, nmn, nsd, pbozh, psi, rbozh, zbozh ,   &
    allocate_iodisk, deallocate_iodisk
  IMPLICIT NONE
  INTEGER(4),PARAMETER:: ndiskc=10
  INTEGER(4)::i, ifsum, igt, ii, iie, iis, ijf, ijf2, imttl, is, j, kk, m, nout,    &
    nprt, nsdm1, icont(ibhar), msel(ibhar,nsd)

  REAL(8):: asel, bamax, cmu0, cnorm, critc, pi, pi2, psisgn, rmaj0,                &
    bmax(nsd), bmin(nsd), bsel(ibhar), psel(ibhar), rsel(ibhar), zsel(ibhar)
  LOGICAL :: lcheck=.false.

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

!      read(ndiskc,7001) nmboz, nsdm1, nfp
      READ(ndiskc) nmboz, nsdm1, nfp
 7001 FORMAT(3i4)
      WRITE(6,6001)  nmboz, nsdm1, nfp

      IF( nmboz - 1 .lt. mmx ) THEN
        WRITE(6,9001) nmboz - 1, mmx
        STOP
      ELSE
        WRITE(6,9003) nmboz - 1, mmx
      ENDIF
      IF( nsdm1 .ne. kmsh    ) THEN
        WRITE(6,9002) nsdm1, kmsh
        STOP
      ELSE
        WRITE(6,9004) nsdm1, kmsh
      ENDIF

      CALL allocate_iodisk


!      do 1000, i = 2, nsd
!c         read (ndiskc,7002) psi(i), eot(i), cui(i), cug(i)
!         read (ndiskc) psi(i), eot(i), cui(i), cug(i)
! 1000 continue

      READ(ndiskc) (psi(i), i = 2, nsd)
      READ(ndiskc) (eot(i), i = 2, nsd)

      READ(ndiskc) (cui(i), i = 2, nsd)
      READ(ndiskc) (cug(i), i = 2, nsd)

! 7002 format(4d25.17)

!--------------------------------c
!    psi_a < 0 case
!--------------------------------c

      ispsi=0
      IF (psi(nsd).lt.0.d0) THEN
        DO i = 2, nsd
          psi(i) =  -psi(i)
          cui(i) =  -cui(i)
          cug(i) =  -cug(i)
      END DO
        ispsi=1
        WRITE(6,*) '*** psi_a < 0. '
        WRITE(6,*) '*** psi, I, g are replaced by -psi, -I, -g.'
      END IF

!--------------------------------c
!    B < 0 case
!--------------------------------c

      IF (bb0.lt.0.d0) THEN
        ispsi=2
        bb0=-bb0
        WRITE(6,*) '*** BB0 < 0. '
        WRITE(6,*) '*** Reverse filed case is assumed.'
      END IF

!--------------------------------c


      DO M = 1, nmboz
        READ(ndiskc)  mboz(m), nboz(m)
        READ(ndiskc) (bbozh(m,i), i = 2, nsd)
      END DO

      DO M = 1, nmboz
        READ(ndiskc)  mboz(m), nboz(m)
        READ(ndiskc) (rbozh(m,i), i = 2, nsd)
      END DO

      DO M = 1, nmboz
        READ(ndiskc)  mboz(m), nboz(m)
        READ(ndiskc) (zbozh(m,i), i = 2, nsd)
      END DO

      DO M = 1, nmboz
        READ(ndiskc)  mboz(m), nboz(m)
        READ(ndiskc) (pbozh(m,i), i = 2, nsd)
      END DO

!$$$      do 101, m = 1, nmboz
!$$$c
!$$$c         read (ndiskc,7003) mboz(m), nboz(m)
!$$$         read (ndiskc) mboz(m), nboz(m)
!$$$         READ(NDISKC) (BBOZH(M,I), I = 2, NSD)
!$$$ 101 continue
!$$$c
!$$$c         do 1001, i = 2, nsd
!$$$c            read (ndiskc,7004) bbozh(m,i), rbozh(m,i)
!$$$c     &                        ,zbozh(m,i), pbozh(m,i)
!$$$c            read (ndiskc) bbozh(m,i), rbozh(m,i)
!$$$c     &                        ,zbozh(m,i), pbozh(m,i)
!$$$c 1002    continue

 7003 FORMAT(2i4)
 7004 FORMAT(4d25.17)

!--------------------------------c
!    psi_a < 0 case
!--------------------------------c
      IF (ispsi.eq.1) THEN
        DO m = 1, nmboz
          DO i = 2, nsd
            zbozh(m,i) = - zbozh(m,i)
            pbozh(m,i) = - pbozh(m,i)
          END DO
        END DO
          WRITE(6,*) '*** zbozh, pbozh are replaced by -zbozh, -pbozh.'
      ENDIF
!--------------------------------c

      WRITE(6,6002)
      WRITE(6,6003) (i, psi(i), eot(i), cui(i), cug(i), i = 2, nsd)

      IF( lcheck ) THEN
        nout  = 10
        nprt  = nsd/nout
        IF( mod(nsd,nout) .ne. 0 ) nprt = nprt + 1
        IF( mod(nsd,nout) .ne. 0 ) nprt = nprt + 1

        DO ii=1,nprt
          iis  =(ii-1)*nout+1
          iie  = ii   *nout
          IF( ii .eq. nprt ) iie  = nsd
          WRITE(6,6011)         ( kk , kk = iis, iie )

          DO j = 1, nmboz
            WRITE(6,6012) j, mboz(j), nboz(j), ( bbozh(j,kk), kk = iis, iie )
          END DO
        END DO
      ENDIF

!.... from half mesh to integer mesh

      DO i = 2, nsd-1
        psi(i)   = ( psi(i) + psi(i+1) )/2.0d0
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

      WRITE(6,*) ' - End 301 -'
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
        psi  (  1) = 1.0d-18
        eot  (  1) = 3.0d0*eot  (  2) - 3.0d0*eot  (  3) + eot  (  4)
        cui  (  1) = 3.0d0*cui  (  2) - 3.0d0*cui  (  3) + cui  (  4)
        cug  (  1) = 3.0d0*cug  (  2) - 3.0d0*cug  (  3) + cug  (  4)

      WRITE(6,*) ' - End 302 -'
!+++++ extrapolation of the values at the boundary

      DO j = 1, nmboz
        bbozh(j,nsd) = 3.0d0*bbozh(j,nsd-1) - 3.0d0*bbozh(j,nsd-2) + bbozh(j,nsd-3)
        rbozh(j,nsd) = 3.0d0*rbozh(j,nsd-1) - 3.0d0*rbozh(j,nsd-2) + rbozh(j,nsd-3)
        zbozh(j,nsd) = 3.0d0*zbozh(j,nsd-1) - 3.0d0*zbozh(j,nsd-2) + zbozh(j,nsd-3)
        pbozh(j,nsd) = 3.0d0*pbozh(j,nsd-1) - 3.0d0*pbozh(j,nsd-2) + pbozh(j,nsd-3)
      END DO
        IF( psi(nsd-1) .gt. 0.0d0 ) THEN
          psisgn   =  1.0d0
        ELSE
          psisgn   = -1.0d0
        ENDIF
!       psi  (nsd) = 0.5d0*psisgn
        psi  (nsd) = 3.0d0*psi(nsd-1) - 3.0d0*psi(nsd-2) + psi(nsd-3)
        eot  (nsd) = 3.0d0*eot(nsd-1) - 3.0d0*eot(nsd-2) + eot(nsd-3)
        cui  (nsd) = 3.0d0*cui(nsd-1) - 3.0d0*cui(nsd-2) + cui(nsd-3)
        cug  (nsd) = 3.0d0*cug(nsd-1) - 3.0d0*cug(nsd-2) + cug(nsd-3)

!$$$c.... normalization to vmec calculation
!$$$c
      cnorm  = rmaj0*bb0/cug(nsd)

!       cnorm = -1.d0
      DO i = 1, nsd
        cui(i)   = cui(i)*cnorm
        cug(i)   = cug(i)*cnorm
        psi(i)   = psi(i)*cnorm
        DO j = 1, nmboz
          bbozh(j,i) = bbozh(j,i)*cnorm
        END DO
      END DO
  304 CONTINUE
!$$$c

      WRITE(6,*) ' - End 303 -'
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

        ifsum     = 0

      WRITE(6,*) ' - End 502 -'

!      do 503 is = 1, nsd
!        ifsum     = ifsum + 1
!  503 continue

!        ifsel( 1) = ifsum

      DO ijf = 2, nmboz
!        ifsel(ijf) = 0
        icont(ijf) = 0
      END DO

      WRITE(6,*) ' - End 504 -'
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
      WRITE(6,*) ' - End 601 -'
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
      END DO

      WRITE(6,*) ' - End 700 -'


      DO ijf = 2, nmboz
        IF (icont(ijf) .ne. 0 ) THEN
          imttl     = imttl + 1
          DO is = 1, nsd
            bco(is,imttl) = bbozh(ijf,is)
          END DO

          cm (imttl)    = mboz(ijf)
          cn (imttl)    = nboz(ijf)
          mmn(imttl)  = mboz(ijf)
          nmn(imttl)  = nboz(ijf)
        END IF
      END DO


      WRITE(6,6031) critc, imttl

      DO is = 1, nsd
        WRITE(6,6032) is, icont(is), bmin(is), bmax(is)
      END DO

      WRITE(6,*)
      WRITE(6,*) ' selected foutier spectrum of b '
      WRITE(6,*)

      nout  = 10
      nprt  = nsd/nout
      IF( mod(nsd,nout) .ne. 0 ) nprt = nprt + 1

      DO ii=1,nprt
        iis  =(ii-1)*nout+1
        iie  = ii   *nout
        IF( ii .eq. nprt ) iie  = nsd
        WRITE(6,6033)         ( kk , kk = iis, iie )

        DO j = 1, imttl
          WRITE(6,6034) j, idint(cm(j)), idint(cn(j)), ( bco(kk,j), kk = iis, iie )
        END DO
      END DO



      WRITE(6,6004)
      WRITE(6,6003) (i, psi(i), eot(i), cui(i)*pi2/(1.0d6*cmu0), &
                   cug(i)*pi2/(1.0d6*cmu0), i = 1, nsd)

      nout  = 10
      nprt  = nsd/nout
      IF( mod(nsd,nout) .ne. 0 ) nprt = nprt + 1

      DO ii=1,nprt
        iis  =(ii-1)*nout+1
        iie  = ii   *nout
        IF( ii .eq. nprt ) iie  = nsd
        WRITE(6,6013) ( kk , kk = iis, iie )

        DO j = 1, mmx + 1
          WRITE(6,6012) j, idint(cm(j)), idint(cn(j)), ( bco(kk,j), kk = iis, iie )
        END DO
      END DO


 6031 FORMAT(//10x,'critical value of spectrum = ',1pd12.5 &
     &        /10x,'the number of modes needed for all surfaces = ',i6 &
     &       //10x,'is',10x,'num. modes',10x,'bmin',10x,'bmax'/)
 6032 FORMAT(  10x,i2,15x,i5,2x,1pd12.5,2x,d12.5)
 6033 FORMAT(//2x,' j',2x,'(  m,  n)',10(5x,'b   _',i2)/)
 6034 FORMAT(  2x,i2,2x,'(',i3,',',i3,')',10(1x,1pd11.4))
 6021 FORMAT(//4x,'zeta /pi = ',1pd12.5 &
     &       //4x,'theta/pi',13x,'b',13x,'r',13x,'z',8x,'phi/pi'/)
 6022 FORMAT(5(2x,1pd12.5))

 6001 FORMAT(//9x,'nmboz  = ',i4,' : nsd - 1 = ',i4,' : nfp    =',i4//)
 6002 FORMAT(//4x,'i',11x,'psi',10x,'iota',13x,'i',13x,'g'/)
 6003 FORMAT(/(2x,i3,4(2x,1pd12.5)))
 6004 FORMAT(//2x,130('+')//10x,'dimensional physical quantities' &
     &       //4x,'i',3x,'psi(t m**2)', &
     &        10x,'iota',4x,'2*pi*i(ma)',4x,'2*pi*g(ma)'/)
 6011 FORMAT(//2x,' j',2x,'( m,  n)',10(5x,'b   _',i2)/)
 6013 FORMAT(//2x,' j',2x,'( m,  n)',10(5x,'b(t)_',i2)/)
 6012 FORMAT(  2x,i2,2x,'(',i2,',',i3,')',10(1x,1pd11.4))
 9001 FORMAT(//10x,'nmboz - 1 = ',i4,' <  mmx  = ',i4//)
 9003 FORMAT(//10x,'nmboz - 1 = ',i4,' >= mmx  = ',i4//)
 9002 FORMAT(//10x,'nsd   - 1 = ',i4,' =/ kmsh = ',i4//)
 9004 FORMAT(//10x,'nsd   - 1 = ',i4,' =  kmsh = ',i4//)


      CALL deallocate_iodisk

      RETURN
      END
