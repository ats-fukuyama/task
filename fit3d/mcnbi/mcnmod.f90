MODULE  mcnmod
  IMPLICIT NONE
  INTEGER(4), PARAMETER ::                                                      &
    ntest=1000, kmsh=60, mmx=99, nsd=61, lplot=2000, iend=100, imsh=50,         &
    nkptl=10, mplot=20000, nplot=1000, mm=20, mm2=ntest/mm,                     &
    nbzmin=-15, nbzmax=15, mbzmax=200, ibhar=(nbzmax-nbzmin+1)*mbzmax,          &
    mvmax=32, mvmax2=2*mvmax-1, mrmax=30, mrmax3=mrmax*3,                       &
    mtmax=50, mtmax3=mtmax*3, mpmax=100,mpmax3=mpmax*3, mxcx=99,                &
    mxmax=100,mymax=100, mymax3=mymax*3
!    ntest=1000, kmsh=100, mmx=49, nsd=101, lplot=2000, iend=100, imsh=50,

  INTEGER(4)  ::                                                                    &
    ich, icoll, ifout, iit, ioptn, iout, ips1, ispsi, istrt, iswrd, it, it0, itmax, &
    ix, iy, jit, jobno, jout, jplot, kout, loss, maxnpl, nfp, nmboz, nout, npl, nsw,&
    ntot


  INTEGER(4), ALLOCATABLE  ::                                                 &
    isw(:), isswt(:), ips(:),           & !(ntest)
    mmn(:), nmn(:),                     & ! (mmx+1),
    mboz(:), nboz(:),                   & ! (ibhar)
    npt(:),nloss(:),                    & ! (lplot)
    iran(:)                               ! (500)

  REAL(8) ::                                                                        &
    a, adiv, ai, alpha, ane, anecgs, anem1, anem2, anew, ani, anicgs, at, b, b0,    &
    bb0, bres, cecgs, ceps0, cfb, cfe0, cfi0, charge, ckb, clight,                  &
    clog, cme, cmi, cmp, cmyu0, cpi, cpi2, crpe, db, dbfi, dbth, ddt, dt, e0int,    &
    einit, el, eloss, em, emin, enorm, eot0, eota, epsa, epsh, erf, esp0, g, pecrh, &
    pna0, pna1, pne0, pne1, psia, r, rtptcl, sgm, taus, tchg, te, teev, teevw, tem1,&
    tem2, tew, ti, tiev, tievw, tim1, tim2, time, tiw, tm2, tmass, ve, vi, vnorm

  REAL(8),  ALLOCATABLE  ::                                                    &
    adinv(:),e0ev(:),ekin0(:),eng(:),fai0(:),pitch0(:), ram(:),rswt(:),sr0(:), &
    swt(:),thta0(:),vel(:),                                                    &  ! (ntest)
    c1et(:),c1g(:),c1i(:),c2et(:),c2g(:),c2i(:), c3et(:),c3g(:),c3i(:),cug(:), &
    cui(:),eot(:), psi(:),                                                     & ! (kmsh+1)
    etlos(:),                                                                  & !(lplot),
    tdep(:), tdepe(:), tdepi(:),                                               & !(mm2),
    cm(:), cn(:),                                                              & !(mmx+1),
    vel0(:),                                                                   & ! (nkptl),
    y(:),                                                                      & ! (4),
    bco(:,:), c1bf(:,:), c2bf(:,:), c3bf(:,:),                                 & !(kmsh+1,mmx+1),
    depepa(:,:), depppa(:,:),                                                  & ! (mpmax3,mm2),
    depera(:,:), deppra(:,:), depvha(:,:), depvra(:,:),                        & ! (mrmax3,mm2),
    depeta(:,:), deppta(:,:),                                                  & !(mtmax3,mm2),
    wmxt(:,:),                                                                 & ! (mrmax3,mm2),
    yc(:,:),                                                                   & !(ntest,4),
    depexa(:,:,:),deppxa(:,:,:),                                               & ! (mxmax,mymax3,mm2),
    depvva(:,:,:,:),                                                           & !(mvmax2,mvmax,mrmax3,mm2)
    fmx0(:,:,:)                                                                  ! (mvmax2,mvmax,mrmax)

  REAL(8), ALLOCATABLE  ::                                             &
    bbozh(:,:), pbozh(:,:), rbozh(:,:), zbozh(:,:),                    & ! (ibhar,nsd)
    denn0(:), engn0(:),                                                & ! (mxcx)
    rgc(:,:), zgc(:,:), pgc(:,:), pang(:,:), tang(:,:),                & ! (mplot,nkptl)
    fs(:),                                                             & ! (iend+1)
    depe(:), depi(:),                                                  & ! (lplot)
    dep1(:), depe1(:), prf1(:), aveng1(:), avens1(:), dtime1(:),       &
    flrf(:), fdrf(:), fderf(:), fld(:),                                & !(nplot)
    echeck(:), bf(:),                                                  & ! (ntest)
    eeer(:), eiir(:), eler(:), elir(:),                                & ! (imsh+1)
    yeer(:), yiir(:)                                                     ! (imsh)

  REAL(8), ALLOCATABLE  ::                                             &
    eir(:,:,:), eer(:,:,:), dir(:,:,:), der(:,:,:),                    & ! (imsh+1,3,mm2)
    dei(:), dee(:),                                                    & ! (mm2)
    f(:), fl(:), ffp(:), ffm(:),                                       &
    dep(:), dtime(:), prf(:) , aveng(:), avengs(:),                    & ! (lplot)
    es(:), vx(:), vy(:), sr1(:),                                       & ! (ntest)
    tmhst(:), vhhst(:), enghst(:), enghsp(:), enghsm(:),               &
    vhhstp(:),vhhstm(:),depths(:), depehs(:), depihs(:)                  !  (100000)

!     ******************************************************
  CONTAINS

  SUBROUTINE allocate_through
    ALLOCATE(isw(ntest), mmn(mmx+1), nmn(mmx+1), isswt(ntest), ips(ntest))
    ALLOCATE(adinv(ntest),e0ev(ntest),ekin0(ntest),eng(ntest),fai0(ntest),pitch0(ntest))
    ALLOCATE(ram(ntest),rswt(ntest),sr0(ntest),swt(ntest),thta0(ntest),vel(ntest))
    ALLOCATE(c1et(kmsh+1),c1g(kmsh+1),c1i(kmsh+1),c2et(kmsh+1),c2g(kmsh+1),c2i(kmsh+1))
    ALLOCATE(c3et(kmsh+1),c3g(kmsh+1),c3i(kmsh+1),cug(kmsh+1), cui(kmsh+1),eot(kmsh+1))
    ALLOCATE(psi(kmsh+1),etlos(lplot),  tdep(mm2),  tdepe(mm2), tdepi(mm2))
    ALLOCATE(cm(mmx+1), cn(mmx+1), vel0(nkptl), y(4))
    ALLOCATE(bco(kmsh+1,mmx+1), c1bf(kmsh+1,mmx+1), c2bf(kmsh+1,mmx+1), c3bf(kmsh+1,mmx+1))
    ALLOCATE(depepa(mpmax3,mm2), depera(mrmax3,mm2), depeta(mtmax3,mm2), depppa(mpmax3,mm2))
    ALLOCATE(deppra(mrmax3,mm2), deppta(mtmax3,mm2), depvha(mrmax3,mm2), depvra(mrmax3,mm2))
    ALLOCATE(wmxt(mrmax3,mm2), yc(ntest,4))
    ALLOCATE(depexa(mxmax,mymax3,mm2),deppxa(mxmax,mymax3,mm2),depvva(mvmax2,mvmax,mrmax3,mm2))
    ALLOCATE(fmx0(mvmax2,mvmax,mrmax))
    RETURN
  END SUBROUTINE allocate_through
  SUBROUTINE deallocate_through
    DEALLOCATE(isw, mmn, nmn, isswt, ips)
    DEALLOCATE(adinv,e0ev,ekin0,eng,fai0,pitch0,ram,rswt,sr0,swt,thta0,vel)
    DEALLOCATE(c1et,c1g,c1i,c2et,c2g,c2i,c3et,c3g,c3i,cug, cui,eot)
    DEALLOCATE(psi,etlos,  tdep,  tdepe, tdepi, cm, cn, vel0, y)
    DEALLOCATE(bco, c1bf, c2bf, c3bf, depepa, depera, depeta, depppa, deppra, deppta)
    DEALLOCATE(depvha, depvra, wmxt, yc)
    DEALLOCATE(depexa,deppxa,depvva, fmx0)
    RETURN
  END SUBROUTINE deallocate_through

  SUBROUTINE allocate_b
    CALL allocate_restrt2
    ALLOCATE(denn0(mxcx), engn0(mxcx))
    RETURN
  END SUBROUTINE allocate_b
  SUBROUTINE deallocate_b
    DEALLOCATE(denn0, engn0)
    RETURN
  END SUBROUTINE deallocate_b

  SUBROUTINE allocate_iodisk
    ALLOCATE(mboz(nmboz), nboz(nmboz))
    ALLOCATE(bbozh(nmboz,nsd), pbozh(nmboz,nsd), rbozh(nmboz,nsd), zbozh(nmboz,nsd))
    RETURN
  END SUBROUTINE allocate_iodisk

  SUBROUTINE deallocate_iodisk
    DEALLOCATE(mboz, nboz,bbozh, pbozh, rbozh, zbozh)
    RETURN
  END SUBROUTINE deallocate_iodisk


  SUBROUTINE allocate_restrt1
    CALL allocate_restrt2
    IF(.not. ALLOCATED(rgc ))   ALLOCATE(rgc(mplot,nkptl))
    IF(.not. ALLOCATED(zgc ))   ALLOCATE(zgc(mplot,nkptl))
    IF(.not. ALLOCATED(pgc ))   ALLOCATE(pgc(mplot,nkptl))
    IF(.not. ALLOCATED(pang ))  ALLOCATE(pang(mplot,nkptl))
    IF(.not. ALLOCATED(tang ))  ALLOCATE(tang(mplot,nkptl))
    IF(.not. ALLOCATED(fs ))    ALLOCATE(fs(iend+1))
    IF(.not. ALLOCATED(depe ))  ALLOCATE(depe(lplot))
    IF(.not. ALLOCATED(depi ))  ALLOCATE(depi(lplot))
    IF(.not. ALLOCATED(dep1 ))  ALLOCATE(dep1(nplot))
    IF(.not. ALLOCATED(depe1 )) ALLOCATE(depe1(nplot))
    IF(.not. ALLOCATED(prf1 ))  ALLOCATE(prf1(nplot))
    IF(.not. ALLOCATED(aveng1)) ALLOCATE(aveng1(nplot))
    IF(.not. ALLOCATED(avens1)) ALLOCATE(avens1(nplot))
    IF(.not. ALLOCATED(dtime1)) ALLOCATE(dtime1(nplot))
    IF(.not. ALLOCATED(flrf ))  ALLOCATE(flrf(nplot))
    IF(.not. ALLOCATED(fdrf ))  ALLOCATE(fdrf(nplot))
    IF(.not. ALLOCATED(fderf )) ALLOCATE(fderf(nplot))
    IF(.not. ALLOCATED(fld ))   ALLOCATE(fld(nplot))
    IF(.not. ALLOCATED(echeck)) ALLOCATE(echeck(ntest))
    IF(.not. ALLOCATED(bf ))    ALLOCATE(bf(ntest))
    IF(.not. ALLOCATED(eeer ))  ALLOCATE(eeer(imsh+1))
    IF(.not. ALLOCATED(eiir ))  ALLOCATE(eiir(imsh+1))
    IF(.not. ALLOCATED(eler ))  ALLOCATE(eler(imsh+1))
    IF(.not. ALLOCATED(elir ))  ALLOCATE(elir(imsh+1))
    IF(.not. ALLOCATED(yeer ))  ALLOCATE(yeer(imsh))
    IF(.not. ALLOCATED(yiir ))  ALLOCATE(yiir(imsh))
    RETURN
  END SUBROUTINE allocate_restrt1

  SUBROUTINE deallocate_restrt1
    DEALLOCATE(rgc, zgc, pgc, echeck, fs, depe, depi, dep1, depe1, prf1)
    DEALLOCATE(aveng1, avens1, dtime1, flrf, fdrf, fderf, fld, bf, eeer,  eiir)
    DEALLOCATE(eler, elir, yeer, yiir, pang, tang)
    RETURN
  END SUBROUTINE deallocate_restrt1

  SUBROUTINE allocate_restrt2
    IF(.not. ALLOCATED(eir ))    ALLOCATE(eir(imsh+1,3,mm2) )
    IF(.not. ALLOCATED(eer ))    ALLOCATE(eer(imsh+1,3,mm2) )
    IF(.not. ALLOCATED(dir ))    ALLOCATE(dir(imsh+1,3,mm2) )
    IF(.not. ALLOCATED(der ))    ALLOCATE(der(imsh+1,3,mm2) )
    IF(.not. ALLOCATED(dei ))    ALLOCATE(dei(mm2) )
    IF(.not. ALLOCATED(dee ))    ALLOCATE(dee(mm2) )
    IF(.not. ALLOCATED(f ))      ALLOCATE(f(iend+1)  )
    IF(.not. ALLOCATED(fl ))     ALLOCATE(fl(iend+1) )
    IF(.not. ALLOCATED(ffp ))    ALLOCATE(ffp(iend+1) )
    IF(.not. ALLOCATED(ffm ))    ALLOCATE(ffm(iend+1) )
    IF(.not. ALLOCATED(dep ))    ALLOCATE(dep(lplot) )
    IF(.not. ALLOCATED(dtime ))  ALLOCATE(dtime(lplot) )
    IF(.not. ALLOCATED(prf ))    ALLOCATE(prf(lplot) )
    IF(.not. ALLOCATED(aveng ))  ALLOCATE(aveng(lplot) )
    IF(.not. ALLOCATED(avengs )) ALLOCATE(avengs(lplot))
    IF(.not. ALLOCATED(npt ))    ALLOCATE(npt(lplot) )
    IF(.not. ALLOCATED(nloss ))  ALLOCATE(nloss(lplot) )
    IF(.not. ALLOCATED(es ))     ALLOCATE(es(ntest)  )
    IF(.not. ALLOCATED(vx ))     ALLOCATE(vx(ntest)  )
    IF(.not. ALLOCATED(vy ))     ALLOCATE(vy(ntest)  )
    IF(.not. ALLOCATED(iran ))   ALLOCATE(iran(500) )
    IF(.not. ALLOCATED(sr1 ))    ALLOCATE(sr1(ntest) )
    IF(.not. ALLOCATED(tmhst ))  ALLOCATE(tmhst(100000) )
    IF(.not. ALLOCATED(vhhst ))  ALLOCATE(vhhst(100000) )
    IF(.not. ALLOCATED(enghst )) ALLOCATE(enghst(100000))
    IF(.not. ALLOCATED(enghsp )) ALLOCATE(enghsp(100000))
    IF(.not. ALLOCATED(enghsm )) ALLOCATE(enghsm(100000))
    IF(.not. ALLOCATED(vhhstp )) ALLOCATE(vhhstp(100000))
    IF(.not. ALLOCATED(vhhstm )) ALLOCATE(vhhstm(100000))
    IF(.not. ALLOCATED(depths )) ALLOCATE(depths(100000))
    IF(.not. ALLOCATED(depehs )) ALLOCATE(depehs(100000))
    IF(.not. ALLOCATED(depihs )) ALLOCATE(depihs(100000))
    RETURN
  END SUBROUTINE allocate_restrt2

  SUBROUTINE deallocate_restrt2
    DEALLOCATE(eir, eer, dir, der, dei, dee, f, ffp, ffm, dep, es, vx, vy, iran, sr1)
    DEALLOCATE(tmhst, vhhst, enghst, enghsp, enghsm, vhhstp, vhhstm)
    RETURN
  END SUBROUTINE deallocate_restrt2
  SUBROUTINE deallocate_restrt22
    DEALLOCATE(depths, depehs, depihs)
    DEALLOCATE(fl, npt, aveng, avengs, dtime, nloss, prf)
    RETURN
  END SUBROUTINE deallocate_restrt22
!     ******************************************************

END MODULE  mcnmod
