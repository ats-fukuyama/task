!     equcomm.f

      module aaa_mod
      implicit none
!      public
!
!=======================================================================
!     #SIZ
!=======================================================================
      integer   isymdm
      integer   irdm,   izdm,   ivdm
      integer   irzdm,  isrzdm
      integer   icvdm,  nsfix
      integer   icvdm1, icvdm2
      integer   itdm
      integer   indm,   inzdm, indmz
      integer   icvmm2, icvmm6, nsfx6
!
!<<ud-symmetry>>
      parameter(isymdm= 1)
!<<equilibrium grid>>
!     parameter (irdm=513,izdm=513,ivdm=201)
      parameter (irdm=257,izdm=257,ivdm=101)
!
      parameter(irzdm=irdm*izdm)
!      parameter(isrzdm=2*(irdm+izdm))
      parameter(isrzdm=3*(irdm+izdm))
!<<vaccuum field>>
      parameter(icvdm= 18, nsfix= 30)
!
      parameter(icvdm1=icvdm+1)
      parameter(icvdm2=icvdm+2)
      parameter(icvmm2=icvdm-2)
      parameter(icvmm6=icvdm-6)
      parameter(nsfx6=nsfix-6)
!-----------------------------------------------------------------------
!<<transport grid>>
      parameter(itdm=101)
!<<particle species>>
      parameter(indm=  2)
      parameter(inzdm= 2)
!
      parameter(indmz=indm+inzdm)
!=======================================================================
!     #PCT --- PHYSICS CONSTANTS  FOR  LIB-2D.FORT77
!=======================================================================
      real*8    cnpi, cnpi2, cnpis, cnmu, cnmu0
      real*8    cnec, cnme,  cnmp,  gam,  gam2
      real*8    denmin, temmin, premin
!
      parameter(cnpi=3.141592654d0,cnpi2=2.d0*cnpi,cnpis=cnpi*cnpi)
      parameter(cnmu=1.25664d-06,cnmu0=1.25664d0)
      parameter(cnec=1.60219d-19)
      parameter(cnme=9.10953d-31)
      parameter(cnmp=1.67265d-27)
      parameter(gam =5.d0/3.d0,gam2=gam-2.d0)
!-----
      parameter(denmin=1.d+17,temmin=10.d0,premin=cnec*denmin*temmin)
!=======================================================================
!     #STP --- STEP CONTROL & OTHER PARAMETERS  FOR  LIB-2D.FORT77
!=======================================================================
      character*70 title
      character*8  xdate
      character*8  device,devnam
!-----
      integer   iudsym,judsym,ivp,step,out,ft05,ft06
!
      end module aaa_mod


      module cnt_mod
      implicit none
      public
!
!=======================================================================
!     #CNT --- #COMDEC  FOR  MEUDAS.FORT77
!=======================================================================
      integer   ieqfm,ieqrd,jeqrd,ieqwt,jeqwt
      integer   ireset,jreset
      integer   ieqout(10)
      real*8    cput
!
!      common/dsk0/ieqfm,ieqrd,jeqrd,ieqwt,jeqwt
!      common/res0/ireset,jreset
!      common/opt0/ieqout
!      common/opt1/cput
!=======================================================================
!     #CNT --- #COMDEC  FOR  LIB-2D.FORT77
!=======================================================================
      integer   iread,jread,iwrite,jwrite,itdsk,jtdsk
      real*8    cputim,cput0
      integer   maxcpu,mcpu
      real*8    time,tmax
      integer   itime,itmax,jtime,jtmax,ktime,ktmax &
               ,ioucnt,joucnt &
               ,itrout(10) &
               ,itrana(10)
!
!<<DISK CONTROL>>
!      common/dsk0x/iread,jread,iwrite,jwrite,itdsk,jtdsk
!<<CPUTIME CONTROL>>
!      common/cpu0/cputim,cput0
!     >           ,maxcpu,mcpu
!<<TIME CONTROL>>
!     >           ,time,tmax
!     >           ,itime,itmax,jtime,jtmax,ktime,ktmax
!<<OUTPUT CONTROL>>
!     >           ,ioucnt,joucnt
!     >           ,itrout
!<<OUTPUT OF POWER BALANCE>>
!     >           ,itrana
!
      end module cnt_mod


      module com_mod
      use aaa_mod
      implicit none
      public
!
!=======================================================================
!     #COM --- #COMDEC  FOR  LIB-2D.FORT77
!=======================================================================
      real*8    ww1(ivdm,4),ww2(ivdm,4)
!      common/www0/ww1,ww2
      end module com_mod


      module eqt_mod
      use aaa_mod
      implicit none
      public
!
!=======================================================================
!     #EQT --- EQUILIBRIUM QUANTITIES AT RO-GRID
!=======================================================================
      integer   nt,ntm
      real*8    gzt(itdm),mut(itdm),nut(itdm) &
               ,hit(itdm),hdt(itdm),sit(itdm),sdt(itdm) &
               ,vlt(itdm),art(itdm) &
               ,ckt(itdm),sst(itdm),aat(itdm),rrt(itdm) &
               ,pdt(itdm),fdt(itdm),cpdt(itdm) &
               ,hih(itdm),vro(itdm),vrh(itdm),sro(itdm),srh(itdm) &
               ,hdh(itdm),aah(itdm),ckh(itdm) &
               ,rph(itdm),rth(itdm),elh(itdm),dlh(itdm) &
               ,bbh(itdm),bih(itdm),ssh(itdm),vlh(itdm) &
               ,fth(itdm),eph(itdm) &
               ,bteqh(itdm),brbpih(itdm) &
               ,r2b2h(itdm)
      integer   nnvlv,mmvlv
      real*8    xxvlv(ivdm),yyvlv(ivdm) &
               ,wwck(ivdm,4),wwaa(ivdm,4) &
               ,wwss(ivdm,4),wwrr(ivdm,4)
!
!      common/eqt0/nt,ntm,gzt,mut,nut,hit,hdt,sit,sdt
!     >           ,vlt,art,ckt,sst,aat,rrt,pdt,fdt,cpdt
!     >           ,hih,vro,vrh,sro,srh,hdh,aah,ckh
!     >           ,rph,rth,elh,dlh,bbh,bih,fth,eph
!     >           ,bteqh,brbpih,r2b2h
!      common/eqt1/nnvlv,mmvlv,xxvlv,yyvlv
!     >           ,wwck,wwaa,wwss,wwrr
      end module eqt_mod


      module equ_mod
      use aaa_mod
      implicit none
      public
!
!=======================================================================
!     #EQU --- #COMDEC  FOR  MEUDAS.FORT77
!=======================================================================
      real*8    psi(irzdm),rcu(irzdm),rbp(irzdm)
      real*8    btv,bts,cpl,tcu
      real*8    raxis,zaxis,saxis,qaxis
      integer   iaxis,iraxis,izaxis,laxis
      real*8    csu(isrzdm),rsu(isrzdm),zsu(isrzdm)
      real*8    rsumax,rsumin,zsumax,zsumin,rzumax,rzumin
      integer   nsumax,nsu
      real*8    eeqmax,eodmax,eadmax
      integer   ieq,ieqmax,iod,iodmax,iad,iadmax
      real*8    bav,baw,bavmax,bavmin
      real*8    error,erch(13),save(ivdm,13)
      real*8    esetup,asol,ssol,esol,qsol
      integer   msetup,nsetup
      real*8    dsep,psep,rspmx,rspmn
      integer   isep,iisep,irsep,izsep
      real*8    bets,beta,betj,betsc,betac,betjc
      real*8    ccpl,ttcu,qsurf,ttpr,ell,trg
      real*8    zzlp,zzli
      real*8    zcsu(isrzdm),zrsu(isrzdm),zzsu(isrzdm)
      integer   nzsu,mzsu
      real*8    zrsmin,zrsmax,zzsmin,zzsmax,rzsmin,rzsmax
      real*8    zsdw,zarv,zvlv,zaav,zckv,zrrv,zssv
      real*8    zrpv,zrtv,zelv,zdlv,zbbv,zbiv
      real*8    q95,qqj,el95,prfac
!
      real*8    cpuequ
      real*8    erreq
      real*8    errod
      real*8    errad
      real*8    errcu
      integer   ieqerr(10)
      real*8    psloop,rloop,zloop
!
      end module equ_mod


      module eqv_mod
      use aaa_mod
      implicit none
      public
!
!=======================================================================
!     #EQV --- #COMDEC  FOR  MEUDAS.FORT77
!=======================================================================
      integer   nv,nvm
      real*8    gzv(ivdm),muv(ivdm),nuv(ivdm)
      real*8    hiv(ivdm),hdv(ivdm),siv(ivdm),sdv(ivdm)
      real*8    siw(ivdm),sdw(ivdm)
      real*8    vlv(ivdm),ckv(ivdm),ssv(ivdm),aav(ivdm),rrv(ivdm)
      real*8    pds(ivdm),fds(ivdm)
      real*8    prv(ivdm),rbv(ivdm),qqv(ivdm)
      real*8    arv(ivdm),rho(ivdm)
      real*8    r2b2(ivdm)
!=======================================================================
      real*8    bbv(ivdm),biv(ivdm)
      real*8    rpv(ivdm),rtv(ivdm),elv(ivdm),dlv(ivdm)
      real*8    ftr(ivdm),epv(ivdm)
      real*8    sha(ivdm),cuv(ivdm),qdv(ivdm)
      integer   nnhit,mmhit
      real*8    xxhit(itdm),yyhit(itdm),zzhit(itdm) &
               ,wwmu(itdm,4),wwnu(itdm,4)
!
      end module eqv_mod


      module geo_mod
      use aaa_mod
      implicit none
      public
!
!=======================================================================
!     #GEO
!=======================================================================
      integer   nr,nz,nrz,nrm,nzm,nzh,nrzh
      real*8    rg(irdm),zg(izdm),rg2(irdm)
      real*8    dr,dz,drz,dr2i,dz2i,ddri,ddzi
      real*8    sf0,sf1,sf2,sf3,sf4
      integer   nsr,nsz,nsrm,nszm,nsr2,nsz2,nsrz,m
      real*8    csrz(irdm)
      integer   irzbnd,jrzbnd
      real*8    rbnd(isrzdm),zbnd(isrzdm)
      integer   lrzbnd(isrzdm)
!
!=======================================================================
      real*8    rcnt
      integer   nnr,nnz,nnrm,nnzm,nnrz
!
      end module geo_mod


      module imp_mod
      use aaa_mod
      implicit none
      public
!
!=======================================================================
!     #IMP --- TRANSPORT PARAMETERS
!=======================================================================
      real*8    fzmas,fzchg
      real*8    zef(itdm),zefnc,zefsp
      real*8    ccze(10)
!
!      common/imp0/fzmas,fzchg
!      common/imp1/zef,zefnc,zefsp
!      common/imp2/ccze
      end module imp_mod


      module par_mod
      use aaa_mod
      implicit none
      public
!=======================================================================
!     #PAR --- #COMDEC  FOR  MEUDAS.FORT77
!=======================================================================
      real*8    rmaj,rpla,zpla,elip,trig,rwmn,rwmx,zwmx,yh,yd
      real*8    elipup,trigup
      real*8    btol,tcur,qaxi,qsur,xxli
      real*8    cd(10),cp(10)
      integer   icd(10),icp(10)
!
      end module par_mod


      module r2d_mod
      use aaa_mod
      implicit none
      public
!
!=======================================================================
!     #R2D --- 2D-INTERFAC
!=======================================================================
      integer   nrr2d,nzr2d,nrzr2d,mrz2d
      real*8    drr2d,dzr2d,rsr2d,zsr2d &
               ,rho2d(irzdm)
      integer   nsur2d,msur2d
      real*8    rsur2d(isrzdm),zsur2d(isrzdm)
      integer   nscl2d,mscl2d
      real*8    rscl2d(isrzdm),zscl2d(isrzdm)
      integer   ionr2d,jonr2d
      real*8    fmsr2d(0:indmz),chgr2d(0:indmz)
      integer   ntr2d,ntr2dm,ntr2db,ntr2dc
      real*8    der2d(itdm,0:indmz),per2d(itdm,0:indmz) &
               ,ter2d(itdm,0:indmz),qir2d(itdm) &
               ,zef2d(itdm) &
               ,vlr2d(itdm)
      character sper2d(0:indmz)*2
!
!      common/r2d0/nrr2d,nzr2d,nrzr2d,mrz2d
!     >           ,drr2d,dzr2d,rsr2d,zsr2d
!     >           ,rho2d
!     >           ,nsur2d,msur2d,rsur2d,zsur2d
!     >           ,nscl2d,mscl2d,rscl2d,zscl2d
!     >           ,ionr2d,jonr2d,fmsr2d,chgr2d
!     >           ,ntr2d,ntr2dm,ntr2db,ntr2dc
!     >           ,der2d,per2d,ter2d,qir2d,zef2d,vlr2d
!      common/r2d1/sper2d
      end module r2d_mod


      module trn_mod
      use aaa_mod
      implicit none
      public
!
!=======================================================================
!     #TRN --- TRANSPORT PARAMETERS
!=======================================================================
      integer   itropt(10),itdopt(0:indmz),itpopt(0:indmz)
      integer   nion,nzion,mion
      real*8    fmass(0:indmz),fchrg(0:indmz)
      character specie(0:indmz)*2
      integer   nro,nrom,nroblk,nroscl,nroblm,nroscm
      real*8    ro(itdm),roh(itdm),rov(itdm),rovh(itdm) &
               ,roscl
      real*8    den(itdm,0:indmz),pre(itdm,0:indmz),qi(itdm) &
               ,tem(itdm,0:indmz),qq(itdm)
      real*8    den0(itdm,0:indmz),pre0(itdm,0:indmz),qi0(itdm) &
               ,tem0(itdm,0:indmz),qq0(itdm)
      real*8    denc(itdm,0:indmz),prec(itdm,0:indmz),qic(itdm) &
               ,temc(itdm,0:indmz),qqc(itdm)
      real*8    denh(itdm,0:indmz),preh(itdm,0:indmz),qih(itdm) &
               ,temh(itdm,0:indmz),qqh(itdm)
      real*8    avden(0:indmz),avpre(0:indmz),avtem(0:indmz) &
               ,totws,wsenr(0:indmz),wsnbi,wsalf,totpow &
               ,ttjoul,ttnbie,ttnbii &
               ,cplt,betst,betat,betjt
      real*8    taue,taup,taup0 &
               ,timold,tprold,wsold,dotw &
               ,cputrn
      real*8    dtime,dtimeh,dtmax,dtmin,dtout,tout
      integer   itrmax,jtrmax
      real*8    etrmax,trdvmx,trdvmn
      integer   itrdbg(10)
!
!      common/trn00/itropt,itdopt,itpopt
!      common/trn0/nion,nzion,mion
!      common/trn1/fmass,fchrg
!      common/trn2/specie
!      common/trn3/nro,nrom,nroblk,nroscl,nroblm,nroscm
!      common/trn4/ro,roh,rov,rovh
!     >           ,roscl
!      common/trn5/den,pre,qi,tem,qq
!      common/trn6/den0,pre0,qi0,tem0,qq0
!      common/trn7/denc,prec,qic,temc,qqc
!      common/trn8/denh,preh,qih,temh,qqh
!      common/trn9/avden,avpre,avtem
!     >           ,totws,wsenr,wsnbi,wsalf,totpow
!     >           ,ttjoul,ttnbie,ttnbii
!     >           ,cplt,betst,betat,betjt
!      common/trn10/taue,taup,taup0
!     >            ,timold,tprold,wsold,dotw
!     >            ,cputrn
!-----
!      common/trt0/dtime,dtimeh,dtmax,dtmin,dtout,tout
!     >           ,itrmax,jtrmax,etrmax,trdvmx,trdvmn
!     >           ,itrdbg
      end module trn_mod


      module vac_mod
      use aaa_mod
      implicit none
      public
!
!=======================================================================
!     #VAC --- #COMDEC  FOR  MEUDAS.FORT77
!=======================================================================
      real*8    rvac(0:nsfix),zvac(0:nsfix)
      real*8    cvac(0:icvdm),cvacst(0:icvdm),cvacwg(0:icvdm), &
                cvact(0:icvdm)
      integer   ivac(0:icvdm),jvac(0:icvdm),msfx,nivac
      real*8    svac(irzdm,icvdm) &
               ,sbvac(isrzdm,icvdm),smvac(irdm,icvdm)
      integer   isvac,jsvac &
               ,icvac(icvdm1),ncvac(icvdm1) &
               ,jcvac(100,icvdm1)
      real*8    cjvac(100,icvdm1) &
               ,rcvac(100,icvdm1),zcvac(100,icvdm1),ccvac(100,icvdm1) &
               ,dcvac(0:9)
      real*8    rcoil(100,icvdm),zcoil(100,icvdm),ccoil(100,icvdm)
      integer   ncoil(0:icvdm)
      character fpfile*20
      integer   ifloop,jfloop
      real*8    efloop &
               ,rfloop(100),zfloop(100),pfloop(100),wfloop(100) &
               ,wsfbi
      integer   iadjpt,jadjpt
      real*8    rrloop(100),zzloop(100),pploop(100)
      integer   ivtab,ivgrp
      real*8    rlimt(100),zlimt(100),slimt(100)
      integer   ilimt,jlimt
      integer   ivbnd,jvbnd
      real*8    rvbnd(50),zvbnd(50)
      end module vac_mod
