MODULE  hfrmod
  IMPLICIT NONE
  INTEGER(4), PARAMETER ::                                                          &
    mm=100, nn=100, imsh=50, imax=3, jmax=25, jmsh=1000, kmsh=60, mmx=98,           &
    nsd=61, ipatc=5001, mmf=60, mfm1=mmf-1, imt=5001,  nmax=4,                      &
    nbzmin=-15, nbzmax=15, mbzmax=200,  ibhar=(nbzmax-nbzmin+1)*mbzmax

  INTEGER(4)  ::                                                                    &
    ibouce, ic05pc, icmax, icoef, icut, icut0, ideb, iflag, inters, ipacc, ipartp,  &
    iplot1, ipt, irunge, jobno, k1, kmax, kount, loop, lostpp, lpi2, mb, mlost, mn0,&
    mnengy, modmax, modmxn, nbeam, nengy, newpar, nfp, nhis, nmboz, noplos, npart,  &
    nplotp, nprntp, npskp, nptf, npzone
  INTEGER(4)  ::                                                                    &
    iecomp(2,3),inumpt(2),iptt(3),lbeam(2,3),nhght(2),npshap(2),nshape(2),nwdth(2)
  INTEGER(4), ALLOCATABLE ::                                           &
    mnumbr(:), nnumbr(:), & ! (mm),
    mboz(:), nboz(:)        ! (ibhar)


  REAL(8) ::                                                                          &
    a0, agg, aggps, aii, aiips, aiott, alpha, ambip, amf, ane00, anwall, ate00,       &
    atewal, ati00, atiwal, awb, b0, bb0, bdfi, bdps, bdth, bpartt, cmax, cnel1,       &
    cnel2, cnel3, cpsie, crtc, crtfe, crthe, crto, cte1, cte2, dbr, dbr2, dplot,      &
    dps, dpsii, dt, dt0, dtx, dxout, dxsav, een0, een00, eps, er1, er1mod, er2,       &
    hdid, hrfac, omega0, phi, phipp, phis, pphi, psibm, psicut, psimin, psis, qmasb,  &
    qmasp, qmc, r0, r00, r1chab, r2chab, rbore, rmax, rmin, rmj0, rpos, rpp, rr, rr0, &
    sb, thets, tlim, volume, vx0, vx00, vy0, vy00, vz0, vz00, x00, xmin, xpos, y00,   &
    ypos, z00, zb, zchab, zf, zmax, zmin, zpos, zz
  REAL(8) ::                                                                           &
    amiss(3),  anexm(10), angle1(2), angle2(2), angle3(2), angle4(2), aphght(2),       &
    aplost(3), apwdth(2), atexm(10), atixm(10), bamp(2), bdvghz(2), bdvgvt(2),bhght(2),&
    bhzfoc(2), blengt(2), bmdep(3), bmtot(3), bvtfoc(2), bwdth(2), cang1(2), cang2(2), &
    cang3(2), cang4(2), ebeam(2), ebkev(2),  hener(6), hrlose(3), orlost(3), passed(3),&
    phi0(2), pplost(3),rmjpvt(2),  rpivot(2), sang1(2), sang2(2), sang3(2), sang4(2),  &
    sgxnmi(3), vbeam(3), beff(2,3), bpart(2,3), efract(2,3),  hbeff(3,2)

  REAL(8), ALLOCATABLE ::                                                         &
    ai(:), aiota(:), gg(:), psi(:), vol(:), vprime(:),                            & ! (0:nn)
    psino(:),                                                                     & ! (0:100)
    canth(:), vpa(:), canfi(:), canps(:), canpth(:), canpfi(:), clamda(:),        & ! (ipatc)
    amu(:), rou(:), energ(:), xxx(:), zzz(:),                                     & ! (ipatc)
    abmnum(:), arnm(:), aznm(:), apnm(:), abnm(:), dabnm(:),                      & ! (mm)
    ttim1(:), tx1(:), ty1(:), tz1(:),                                             & ! (imt)
    cm(:), cn(:),                                                                 & ! (mmx+1)
    psi2(:), cug(:), cui(:), eot(:),                                              & ! (kmsh+1)
    xte(:), xnel(:), xnion(:), psif(:),                                           & ! (61)
    dtime(:), psivol(:), sivol(:),                                                & ! (101)
    xp(:),                                                                        & ! (200)
    prpos(:),pzpos(:),                                                            & ! (50000)
    bnmn(:,:), cb(:,:), cp(:,:), cr(:,:), cz(:,:), db(:,:), dp(:,:), dr(:,:),     & ! (mm,0:nn),
    dz(:,:), eb(:,:), znmn(:,:), ep(:,:), er(:,:), ez(:,:), pnmn(:,:), rnmn(:,:), & ! (mm,0:nn),
    bnm(:,:), rnm(:,:), znm(:,:), pnm(:,:),                                       & ! (mm,0:nn)
    b2inm(:,:),                                                                   & ! (mm,0:nn)
    hr2(:,:), hrbeam(:,:), hrr1(:,:), hrr2(:,:),                                  & ! (3,61),
    sgxn(:,:), sgxne(:,:), sgxni(:,:),                                            & ! (3,61)
    yyp(:,:), yp(:,:),                                                            & ! (10,200),
    bco(:,:),                                                                     & ! (kmsh+1,mmx+1)
    bbozh(:,:), pbozh(:,:), rbozh(:,:), zbozh(:,:),                               & ! (ibhar,nsd)
    hbeam(:,:,:),hdep(:,:,:)                                                        ! (61,3,2)

  LOGICAL :: lmonte, lorbp, lxy

!*****************************************

  CONTAINS

  SUBROUTINE allocate_through
    ALLOCATE(mnumbr(mm),nnumbr(mm))
    ALLOCATE(abmnum(mm),ai(0:nn), aiota(0:nn), gg(0:nn), psi(0:nn), vol(0:nn))
    ALLOCATE(canth(ipatc),psino(0:100), dtime(101))
    ALLOCATE(bnmn(mm,0:nn), cb(mm,0:nn), cp(mm,0:nn), cr(mm,0:nn), cz(mm,0:nn))
    ALLOCATE(db(mm,0:nn), dp(mm,0:nn), dr(mm,0:nn), dz(mm,0:nn), eb(mm,0:nn), znmn(mm,0:nn))
    ALLOCATE(ep(mm,0:nn), er(mm,0:nn), ez(mm,0:nn), pnmn(mm,0:nn), rnmn(mm,0:nn))
    ALLOCATE(pbozh(ibhar,nsd),rbozh(ibhar,nsd))
    ALLOCATE(hr2(3,61), hrbeam(3,61), hrr1(3,61), hrr2(3,61), yyp(10,200))
    ALLOCATE(hbeam(61,3,2),hdep(61,3,2))
    RETURN
  END SUBROUTINE
  SUBROUTINE deallocate_through
    DEALLOCATE(mnumbr,nnumbr)
    DEALLOCATE(abmnum, ai, aiota, gg, psi, vol, canth,psino, dtime)
    DEALLOCATE(bnmn, cb, cp, cr, cz, db, dp, dr, dz, eb, znmn, ep, er, ez, pnmn, rnmn)
    DEALLOCATE(pbozh,rbozh,hr2, hrbeam, hrr1, hrr2, yyp)
    DEALLOCATE(hbeam,hdep)
    RETURN
  END SUBROUTINE

  SUBROUTINE allocate_a
    ALLOCATE(bnm(mm,0:nn), rnm(mm,0:nn), znm(mm,0:nn), pnm(mm,0:nn))
    ALLOCATE(vprime(0:nn))
    CALL allocate_iodisk
    RETURN
  END SUBROUTINE

  SUBROUTINE deallocate_a
    DEALLOCATE(bnm, rnm, znm, pnm, vprime)
    CALL deallocate_iodisk
    RETURN
  END SUBROUTINE


  SUBROUTINE allocate_b
    ALLOCATE(vpa(ipatc), canfi(ipatc), canps(ipatc), canpth(ipatc), canpfi(ipatc))
    ALLOCATE(clamda(ipatc), amu(ipatc), rou(ipatc))
    ALLOCATE(energ(ipatc), arnm(mm), aznm(mm), apnm(mm))
    ALLOCATE(abnm(mm), dabnm(mm))
    ALLOCATE(xte(61), xnel(61), xnion(61))
    ALLOCATE(sgxn(3,61), prpos(50000),pzpos(50000), psif(61))
    ALLOCATE(psivol(101), sivol(101), xp(200), yp(10,200))
    RETURN
  END SUBROUTINE

  SUBROUTINE deallocate_b
    DEALLOCATE(vpa, canfi, canps, canpth, canpfi, clamda, amu, rou, energ)
    DEALLOCATE(arnm, aznm, apnm, abnm, dabnm, xte, xnel, xnion, sgxn)
    DEALLOCATE(prpos,pzpos, psif, psivol, sivol, xp, yp)
    RETURN
  END SUBROUTINE


  SUBROUTINE allocate_iodisk
    ALLOCATE(mboz(ibhar), nboz(ibhar))
    ALLOCATE(cm(mmx+1), cn(mmx+1), psi2(kmsh+1), cug(kmsh+1), cui(kmsh+1), eot(kmsh+1))
    ALLOCATE(bco(kmsh+1,mmx+1), b2inm(mm,0:nn), bbozh(ibhar,nsd), zbozh(ibhar,nsd))
    RETURN
  END SUBROUTINE

  SUBROUTINE deallocate_iodisk
    DEALLOCATE(mboz, nboz, cm, cn, psi2, cug, cui, eot, bco, b2inm, bbozh, zbozh)
    RETURN
  END SUBROUTINE

  SUBROUTINE allocate_crsec4
    ALLOCATE(sgxni(3,61), sgxne(3,61))
    RETURN
  END SUBROUTINE

  SUBROUTINE deallocate_crsec4
    DEALLOCATE(sgxni, sgxne)
    RETURN
  END SUBROUTINE

  SUBROUTINE allocate_output1
    ALLOCATE(ttim1(imt), tx1(imt), ty1(imt), tz1(imt))
    RETURN
  END SUBROUTINE

  SUBROUTINE deallocate_output1
    DEALLOCATE(ttim1, tx1, ty1, tz1)
    RETURN
  END SUBROUTINE

  SUBROUTINE allocate_output2
    ALLOCATE(xxx(ipatc), zzz(ipatc))
    RETURN
  END SUBROUTINE

  SUBROUTINE deallocate_output2
    DEALLOCATE(xxx, zzz)
    RETURN
  END SUBROUTINE

END MODULE  hfrmod
