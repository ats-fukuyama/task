module tx_graphic
  implicit none
  private
  integer(4), parameter :: NGRM=20, NGTM=5000, NGVM=5000, &
       &                   NGYRM=352, NGYTM=77, NGYVM=70, &
       &                   NGPRM=31, NGPTM=9, NGPVM=18
  type graphic_data
     real(4), dimension(:)    , allocatable :: gnrm
     real(4), dimension(:,:,:), allocatable :: v
  end type graphic_data
  type(graphic_data) :: GY, GYT

  real(4) :: GXM, GYM, GYS
  integer(4) :: MODEG, MODEGL, NGR, NGT, NGVV, NGRSTP, NGTSTP, NGVSTP, NP
  real(4), dimension(:),     allocatable :: GX
  real(4), dimension(:,:,:), allocatable :: GQY!, GY, GYT
  real(4), dimension(0:NGRM) :: GT
  real(4), dimension(0:NGTM) :: GTX
  real(4), dimension(0:NGVM) :: GVX
  real(4), dimension(1:NGYRM) :: gDIV
  real(4), dimension(0:NGTM,1:NGYTM) :: GTY
  real(4), dimension(0:NGVM,1:NGYVM) :: GVY
  public :: TXGOUT, TX_GRAPH_SAVE, TXSTGR, TXSTGT, TXSTGV, TXSTGQ, allocate_txgraf, &
       &    deallocate_txgraf, &
       &    MODEG, MODEGL, NGR, NGT, NGVV, NGYRM, NGYTM, NGYVM, NGRM, NGRSTP, NGTSTP, NGVSTP, &
       &    GY, GQY, GYT, GT, GTX, GVX, GTY, GVY, gDIV

#ifndef nonGSAF
  interface TXWPS
     module procedure TXWPSD
     module procedure TXWPSI
     module procedure TXWPSS
  end interface
#endif

contains

  !***************************************************************
  !
  !   Allocate and deallocate arrays
  !
  !***************************************************************
  
  subroutine allocate_txgraf(ier, icont_in)
    use tx_commons, only : NRMAX, NCM, NQMAX

    integer(4), intent(out) :: ier
    integer(4), intent(in), optional :: icont_in
    integer(4) :: N
    integer(4), dimension(1:3) :: ierl

    ierl(:) = 0

    ! allocation check
    if(allocated(GX)) then
       if(present(icont_in) .and. icont_in /= 0) then
          call deallocate_txgraf
       else
          return
       end if
    end if

    N = NRMAX

    allocate(GX(0:N), GQY(0:N,1:NCM,1:NQMAX)                    , stat = ierl(1))
    allocate(GY%gnrm(1:NGYRM), GYT%gnrm(1:NGYRM), source=1.0    , stat = ierl(2)) 
    allocate(GY%v(0:N,0:NGRM,1:NGYRM), GYT%v(0:N,0:NGTM,1:NGYRM), stat = ierl(3))
    ier = sum(ierl)

    ! All the memories allocated above are clear if some errors occur.
    if(ier /= 0) then
       write(6,*) "XX Allocation error in allocate_txgraf"
       call deallocate_txgraf
    end if

  end subroutine allocate_txgraf

  subroutine deallocate_txgraf

    if(allocated(GX)) then
       deallocate(GX,GQY)
       deallocate(GY%gnrm,GY%v,GYT%gnrm,GYT%v)
    end if

  end subroutine deallocate_txgraf
  
  !***************************************************************
  !
  !   GRAPHIC Command loop
  !
  !***************************************************************
#ifdef nonGSAF
  subroutine txgout
  end subroutine txgout
#else
  subroutine TXGOUT
    use libbes, only : BESINX
    use libgrf, only : GRD1D
    use tx_commons, only : T_TX, TPRE, NQMAX, NRMAX, PI, ieqread!, &
!         &                 rhob, RA, RR, NRA, thrp, kappa, IRPIN, DltRPn, NTCOIL, DltRP_mid, DltRP, epst, rho 
    use tx_interface, only : TXGRUR, TOUPPER, TXGLOD
    use tx_ripple, only : ripple

    integer(4) :: MODE, NGPR, NGPT, NGPV, NGYR, NQ, NQL, NGF, NGFMAX, I, IST, NGRT, NG, IER, J, NGTL
    real(4), dimension(:,:,:), allocatable :: GYL
    character(len=5) :: STR, STR2
    character(len=1) :: KID1, KID2

    integer(4), parameter :: NXMAX = 51, NYMAX = 51, NTHMAX=51
    integer(4) :: NR, NGTDO, IND, IFNT
    integer(4), dimension(:), allocatable :: NGPRL
    real(4) :: GMAX, GMIN, GZ
    real(4), dimension(:), allocatable :: GPXY_IN, GMAXA, GMINA
    real(4), dimension(:,:), allocatable :: GYL2
    real(8), dimension(:), allocatable :: FX
    real(8), dimension(:,:), allocatable :: FY
    real(8) :: TL
    character(len=60) :: STRL
    character(len=60), dimension(:), allocatable :: STRA
!    integer(4) :: NGULEN
!    integer(4) :: NX, NY, NTH, IPAT, IPRD, NSTEP
!    integer(4), dimension(:,:,:), allocatable :: KA
!    real(4) :: DR, DZ, AL, ZORG, ZSTEP, DltRPnL, GXMAX
!    real(4) :: GSRMIN,GSRMAX,GSRSTP,GSZMIN,GSZMAX,GSZSTP
!    real(4), dimension(1:4) :: GPXY
!    real(4), dimension(:), allocatable :: RRL, ZZL
!    real(4), dimension(:,:), allocatable :: VAL, GRPL
!    real(8) :: rhol, theta, dtheta, kappal, thetab
!    character(len=17) :: KOUT

    !     *** MENU ***

    MODE = MODEGL
    OUTER : do
       write(6,*) '# select : Rn: Tn: Un: Vn: A,B: C: E: Dn: Fn: Gn: Hn: S,L,M:file'
       write(6,*) '           W(R,T,V)n:write I:init X:exit'
       read(5,'(A5)',IOSTAT=IST) STR
       if (IST > 0) then
          write(6,*) '### ERROR : Invalid Command : ', STR
          cycle
       else if (IST < 0) then
          return
       end if

       KID1 = STR(1:1)
       KID2 = STR(2:2)
       call TOUPPER(KID1)
       call TOUPPER(KID2)

       select case(KID1)
       case('S')
          call TXGSAV

       case('L')
          call TXGLOD(IER)
          if(IER /= 0) cycle OUTER
          call TXPRFG

       case('I')
          ! *** Initialization ***
          T_TX = 0.D0
          TPRE = 0.D0
          NGT = -1
          NGVV = -1
          call TXSTGT(real(T_TX))
          call TXSTGV(real(T_TX))

       case('R')
          ! *** Time evolution of radial profiles ***
          select case(KID2)
          case('A')
             call TXGRFR(-1,MODE)
          case('B')
             call TXGRFR(-4,MODE)
          case('C')
             call TXGRFR(-7,MODE)
          case('D')
             call TXGRFR(-13,MODE)
          case('E')
             call TXGRFR(-17,MODE)
          case DEFAULT
             read(STR(2:5),*,IOSTAT=IST) NGPR
             if (IST /= 0) then
                write(6,*) '### ERROR : Invalid Command : ', STR
                cycle OUTER
             end if
             if (NGPR >= 0 .and. NGPR <= NGPRM) then
                call TXGRFR(NGPR,MODE)
             end if
          end select

       case('V')
          ! *** Time evolution of variables at a certain position ***
          select case(KID2)
          case('A')
             do NGPV = 1, NGPVM
                call TXGRFV(NGPV,MODE)
             end do
          case('B')
             do NGPV = 1, 6
                call TXGRFV(NGPV,MODE)
             end do
             do NGPV = 9, 10
                call TXGRFV(NGPV,MODE)
             end do

          case DEFAULT
             read(STR(2:5),*,IOSTAT=IST) NGPV
             if (IST < 0) then
                write(6,*) '### ERROR : Invalid Command : ', STR
                cycle OUTER
             end if
             if (NGPV >= 1 .and. NGPV <= NGPVM) then
                call TXGRFV(NGPV,MODE)
             end if
          end select

       case('T')
          ! *** Time evolution of global plasma quantities ***
          select case(KID2)
          case('A')
             do NGPT = 1, NGPTM
                call TXGRFT(NGPT,MODE)
             end do

          case('B')
             do NGPT = 1, 6
                call TXGRFT(NGPT,MODE)
             end do
             do NGPT = 9, 10
                call TXGRFT(NGPT,MODE)
             end do

          case DEFAULT
             read(STR(2:5),*,IOSTAT=IST) NGPT
             if (IST < 0) then
                write(6,*) '### ERROR : Invalid Command : ', STR
                cycle OUTER
             end if
             if (NGPT >= 1 .and. NGPT <= NGPTM) then
                call TXGRFT(NGPT,MODE)
             end if
          end select

       case('U')
          ! *** Radial profiles of the balance among the terms in each equation ***
          if(T_TX == 0.D0) cycle OUTER
          if (KID2 == 'A') then
             do NQ = 1, NQMAX, 4
                call PAGES
                do NQL=NQ,min(NQ+3,NQMAX)
                   call TXGRFQ(NQL,mod(NQL-1,4)+1)
                end do
                call PAGEE
             end do

          else
             read(STR(2:5),*,IOSTAT=IST) NQ
             if (IST < 0) then
                write(6,*) '### ERROR : Invalid Command : ', STR
                cycle OUTER
             end if
             if (NQ >= 1 .and. NQ <= NQMAX) then
                call PAGES
!!$              do NQL=NQ,min(NQ+3,NQMAX)
!!$                 call TXGRFQ(NQL,NQL-NQ+1)
!!$              end do
                call TXGRFQ(NQ,5)
                call PAGEE
             end if
          end if

       case('A')
          ! *** Sequential display of radial profiles (version A) ***
          do NGPR = 1, NGPRM
             call TXGRFR(NGPR,MODE)
          end do

       case('B')
          ! *** Sequential display of radial profiles (version B) ***
          do NGPR = 1, 6
             call TXGRFR(NGPR,MODE)
          end do
          do NGPR = 9, 10
             call TXGRFR(NGPR,MODE)
          end do

       case('C')
          ! *** Comparison of neoclassical characteristics ***
          if (KID2 == 'A') then
             call TXGRCPA(MODE)
          else
             call TXGRCP(MODE)
          end if

       case('D')
          ! *** Time evolution of radial profiles in 3D phase space ***
          read(STR(2:5),*,IOSTAT=IST) NGYR
          if (IST /= 0) then
             write(6,*) '### ERROR : Invalid Command : ', STR
             cycle OUTER
          end if
          if (NGYR >= 0 .and. NGYR <= NGYRM) then
             call TXGRUR(GX,GTX,GYT%v(0:NRMAX,0:NGT,NGYR),NRMAX,NGT,NGTM)!,STR,KV,MODE)
          end if

       case('F')
          ! *** Radial profile of a certain variable ***
          read(STR(2:5),*,IOSTAT=IST) NGYR
          if (IST /= 0) then
             write(6,*) '### ERROR : Invalid Command : ', STR
             cycle OUTER
          end if
          if (NGYR >= 1 .and. NGYR <= NGYRM) then
             allocate(FX(0:NRMAX),FY(0:NRMAX,1))
             do NR=0,NRMAX
                FX(NR)=real(GX(NR),8)
                FY(NR,1)=real(GYT%v(NR,NGT,NGYR),8)
             end do
             STRL='@profile'//STR(2:5)//'@'
             call PAGES
             call GRD1D(0,FX,FY,NRMAX+1,NRMAX+1,1,STRL,0)
             call PAGEE
             deallocate(FX,FY)
          end if

       case('G')
          ! *** Time evolution of radial profile of a certain variable ***
          read(STR(2:5),*,IOSTAT=IST) NGYR
          if (IST /= 0) then
             write(6,*) '### ERROR : Invalid Command : ', STR
             cycle OUTER
          end if
          if (NGYR >= 1 .and. NGYR <= NGYRM) then
             allocate(FX(0:NRMAX),FY(0:NRMAX,0:NGT))
             do NR=0,NRMAX
                FX(NR)=real(GX(NR),8)
             end do
             do NGTDO=0,NGT
             do NR=0,NRMAX
                FY(NR,NGTDO)=real(GYT%v(NR,NGTDO,NGYR),8)
             end do
             end do
             STRL='@profile'//STR(2:5)//'@'
             call PAGES
             call GRD1D(0,FX,FY,NRMAX+1,NRMAX+1,NGT+1,STRL,0)
             call PAGEE
             deallocate(FX,FY)
          end if

       case('E')
          if( ieqread >= 2) call psi_out_gnuplot

!!$          ! *** Contour of ripple amplitude ***
!!$          call PAGES
!!$          call INQFNT(IFNT)
!!$          call SETFNT(32)
!!$
!!$          if(IRPIN == 0) then
!!$             DltRPnL = real(DltRPn) * 100.0
!!$
!!$             allocate(RRL(1:NXMAX),ZZL(1:NYMAX),VAL(1:NXMAX,1:NYMAX),KA(1:8,1:NXMAX,1:NYMAX))
!!$             DR = (4.5 - 2.0) / (NXMAX - 1)
!!$             DZ =  1.5        / (NYMAX - 1)
!!$             do NX = 1, NXMAX
!!$                RRL(NX) = 2.0 + DR * (NX - 1)
!!$                do NY = 1, NYMAX
!!$                   ZZL(NY) = DZ * (NY - 1)
!!$                   AL = SQRT(((RRL(NX) - 2.4)**2 + ZZL(NY)**2) * (2.4 / RRL(NX)))
!!$                   VAL(NX,NY) = real(DltRPnL * BESINX(0,NTCOIL/2.4D0*AL))
!!$                end do
!!$             end do
!!$
!!$             call GQSCAL(2.0,4.5,GSRMIN,GSRMAX,GSRSTP)
!!$             call GQSCAL(0.0,1.5,GSZMIN,GSZMAX,GSZSTP)
!!$
!!$             call GDEFIN(3.0,18.0,1.5,10.5,2.0,4.5,0.0,1.5)
!!$             call GFRAME
!!$
!!$             call GSCALE(GSRMIN,GSRSTP,0.0,0.0,0.1,9)
!!$             call GVALUE(GSRMIN,GSRSTP*2,0.0,0.0,NGULEN(GSRSTP))
!!$             call GSCALE(0.0,0.0,0.0,GSZSTP,0.1,9)
!!$             call GVALUE(0.0,0.0,0.0,GSZSTP*2,NGULEN(GSZSTP*2))
!!$
!!$             IPAT = 1
!!$             IPRD = 0
!!$
!!$             ZORG  = 0.0003
!!$             ZSTEP = 0.0002
!!$             NSTEP = 3
!!$             call CONTP2(VAL,RRL,ZZL,NXMAX,NXMAX,NYMAX,ZORG,ZSTEP,NSTEP,IPRD,IPAT*1,KA)
!!$
!!$             ZORG  = 0.001
!!$             ZSTEP = 0.002
!!$             NSTEP = 4
!!$             call CONTP2(VAL,RRL,ZZL,NXMAX,NXMAX,NYMAX,ZORG,ZSTEP,NSTEP,IPRD,IPAT*2,KA)
!!$
!!$             ZORG  = 0.01
!!$             ZSTEP = 0.02
!!$             NSTEP = 4
!!$             call CONTP2(VAL,RRL,ZZL,NXMAX,NXMAX,NYMAX,ZORG,ZSTEP,NSTEP,IPRD,IPAT*3,KA)
!!$
!!$             ZORG  = 0.0
!!$             ZSTEP = 0.1
!!$             NSTEP = 2
!!$             call CONTP2(VAL,RRL,ZZL,NXMAX,NXMAX,NYMAX,ZORG,ZSTEP,NSTEP,IPRD,IPAT*4,KA)
!!$
!!$             ZORG  = 0.0
!!$             ZSTEP = 0.25
!!$             NSTEP = 4
!!$             call CONTP2(VAL,RRL,ZZL,NXMAX,NXMAX,NYMAX,ZORG,ZSTEP,NSTEP,IPRD,IPAT*4,KA)
!!$
!!$             ZORG  = 1.0
!!$             ZSTEP = 0.5
!!$             NSTEP = 3
!!$             call CONTP2(VAL,RRL,ZZL,NXMAX,NXMAX,NYMAX,ZORG,ZSTEP,NSTEP,IPRD,IPAT*5,KA)
!!$
!!$             deallocate(RRL,ZZL,VAL,KA)
!!$
!!$             allocate(RRL(1:NTHMAX),ZZL(1:NTHMAX))
!!$             dtheta = PI / (NTHMAX - 1)
!!$             theta = 0.d0
!!$             do I = 1, 2
!!$                if(I == 1) then
!!$                   kappal = 1.d0
!!$                else
!!$                   kappal = kappa
!!$                end if
!!$                do NTH = 1, NTHMAX
!!$                   theta = (NTH - 1) * dtheta
!!$                   RRL(NTH) = real(RR + RA * cos(theta))
!!$                   ZZL(NTH) = real(kappal * RA * sin(theta))
!!$                   if(ABS(ZZL(NTH)) < EPSILON(1.0)) ZZL(NTH) = 0.0
!!$                end do
!!$                call LINES2D(RRL,ZZL,NTHMAX)
!!$             end do
!!$             deallocate(RRL,ZZL)
!!$
!!$             allocate(RRL(1:2*NRMAX),ZZL(1:2*NRMAX))
!!$             RRL(1:NRMAX) = RR * (1.d0 + epst(NRMAX:1:-1) * COS(thrp(1:NRMAX)))
!!$             ZZL(1:NRMAX) = kappa * RR * epst(NRMAX:1:-1) * SIN(thrp(1:NRMAX))
!!$             RRL(NRMAX+1:2*NRMAX) = RR * (1.d0 + epst(1:NRMAX) * COS(thrp(NRMAX+1:2*NRMAX)))
!!$             ZZL(NRMAX+1:2*NRMAX) = kappa * RR * epst(1:NRMAX) * SIN(thrp(NRMAX+1:2*NRMAX))
!!$             call LINES2D(RRL,ZZL,2*NRMAX)
!!$             deallocate(RRL,ZZL)
!!$          end if
!!$
!!$          !  *** 1D graphic ***
!!$
!!$          if (MODEG == 2) then
!!$             IND = 9
!!$          else
!!$             IND = 0
!!$          end if
!!$
!!$          GPXY(1) =  3.0
!!$          GPXY(2) = 11.5
!!$          GPXY(3) = 11.5
!!$          GPXY(4) = 17.5
!!$          GXMAX = real(rhob)
!!$          allocate(GRPL(0:NRMAX,1:2))
!!$          do NR = 0, NRMAX
!!$             if(IRPIN == 0) then
!!$                GRPL(NR,1) = real(ripple(rho(NR),0.D0,1.D0))
!!$             else
!!$                GRPL(NR,1) = real(DltRP_mid(NR))
!!$             end if
!!$          end do
!!$          write(KOUT,'(F6.4)') GRPL(NRA,1)
!!$          KOUT = '$#d$#$-a$=='//KOUT(1:6)
!!$          call GTEXT(GPXY(1)+5.5,GPXY(4)-1.5,KOUT,len_trim(KOUT),0)
!!$          STRL = '@DltRP at mid-plane(r)@'
!!$          call TXGRAF(GPXY,GX,GRPL,NRMAX+1,NRMAX+1,1,0.0,GXMAX,STRL,0.26,MODE,IND,0)
!!$
!!$          GPXY(1) = 13.5
!!$          GPXY(2) = 22.0
!!$          GPXY(3) = 11.5
!!$          GPXY(4) = 17.5
!!$          GXMAX = real(rhob)
!!$          thetab = 0.5D0 * PI
!!$          do NR = 0, NRMAX
!!$             rhol = rho(NR) * (1.D0 + (kappa - 1.D0) * sin(thetab))
!!$             if(IRPIN == 0) then
!!$                GRPL(NR,1) = real(ripple(rhol,thetab,1.D0))
!!$                GRPL(NR,2) = real(ripple(rho(NR),thetab,1.D0))
!!$                NG = 2
!!$             else
!!$                GRPL(NR,1) = real(DltRP(NR))
!!$                NG = 1
!!$             end if
!!$          end do
!!$          if(IRPIN == 0) then
!!$             STRL = '@DltRP at tip point with and w/o elongation(r)@'
!!$          else
!!$             STRL = '@DltRP at tip point(r)@'
!!$          end if
!!$          call TXGRAF(GPXY,GX,GRPL,NRMAX+1,NRMAX+1,NG,0.0,GXMAX,STRL,0.26,MODE,IND,0)
!!$          deallocate(GRPL)
!!$
!!$          call SETFNT(IFNT)
!!$          call PAGEE

       case('W')
          ! *** Write out numerical values of a certain variable ***
          NG = 0
          read(STR(3:5),*,IOSTAT=IST) NG
          if (IST /= 0 .and. KID2 /= 'A') then
             write(6,*) '### ERROR : Invalid Command : ', STR
             cycle OUTER
          end if
          call write_console(NG,KID2)

       case('M')
          ! *** Compare the already-saved graphics (Max.6) ***
          do
             if(KID2 == 'T') then
                write(6,*) '## Number of time (Max = 6):'
             else
                write(6,*) '## Number of files (Max = 6):'
             end if
             read(5,*,IOSTAT=IST) NGFMAX
             if (IST > 0) then
                cycle
             else if (IST < 0) then
                cycle OUTER
             else
                exit
             end if
          end do
          if(T_TX == 0.D0) then
             NGRT = 1
          else
             NGRT = 0
             allocate(GYL(0:NRMAX,0:6,1:NGYRM))
             ! Current graphic data temporarily saved in GYL
             GYL(0:NRMAX,0,1:NGYRM) = GY%v(0:NRMAX,NGR,1:NGYRM)
          end if
          NGR=0
          if(KID2 == 'T') then
             call TXGLOD(IER)
             if(IER /= 0) cycle OUTER
             do NGF=1,NGFMAX
                LOOP_NGF: do
                   write(6,'(2(A,I1))') '## INPUT time: ', NGF,' / ',NGFMAX
                   read(5,*,IOSTAT=IST) TL
                   if (IST > 0) then
                      cycle
                   else if (IST < 0) then
                      cycle OUTER
                   else
                      call return_NGT_from_T(NGTL,TL,IER)
                      if(IER == 0) then
                         exit LOOP_NGF
                      else
                         cycle LOOP_NGF
                      end if
                   end if
                end do LOOP_NGF
                if(allocated(GYL) .eqv. .false.) allocate(GYL(0:NRMAX,0:6,1:NGYRM))
                GYL(0:NRMAX,NGF-NGRT,1:NGYRM) = GYT%v(0:NRMAX,NGTL,1:NGYRM)
             end do
          else
             do NGF=1,NGFMAX
                call TXGLOD(IER)
                if(IER /= 0) cycle OUTER
                if(allocated(GYL) .eqv. .false.) allocate(GYL(0:NRMAX,0:6,1:NGYRM))
                GYL(0:NRMAX,NGF-NGRT,1:NGYRM) = GY%v(0:NRMAX,NGR,1:NGYRM)
             end do
          end if
          call TXPRFG
          NGR = NGFMAX-NGRT
          ! Graphic data restored to GY
          GY%v(0:NRMAX,0:NGR,1:NGYRM) = GYL(0:NRMAX,0:NGR,1:NGYRM)
          deallocate(GYL)
          do 
             do
                write(6,*) '## INPUT GRAPH NUMBER: A,B,C,RA,RB,RC,RD,RE,Rn,X:exit'
                read(5,'(A5)',IOSTAT=IST) STR2
                if (IST > 0) then
                   cycle
                else if (IST < 0) then
                   cycle OUTER
                else
                   exit
                end if
             end do
             KID1 = STR2(1:1)
             KID2 = STR2(2:2)
             call TOUPPER(KID1)
             call TOUPPER(KID2)
             select case(KID1)
             case('A')
                !     Correspond to GMA
                do I = 1, NGPRM
                   call TXGRFR(I,MODE)
                end do
                !     Correspond to GMB
             case('B')
                do I = 1, 6
                   call TXGRFR(I,MODE)
                end do
                do I = 9, 10
                   call TXGRFR(I,MODE)
                end do
                !     Correspond to GMC
             case('C')
                do
                   write(6,*) '## HOW MANY GRAPHs ?: 4 or 6'
                   read(5,'(I1)',IOSTAT=IST) NGPR
                   if (IST > 0) then
                      cycle
                   else if (IST < 0) then
                      cycle OUTER
                   else
                      if(NGPR /= 4 .and. NGPR /= 6) cycle
                      exit
                   end if
                end do
                allocate(NGPRL(1:NGPR),GYL2(0:NRMAX,0:NGR),GPXY_IN(1:4))
                write(6,'(A,I1,A)') '## CHOOSE ',NGPR, ' GRAPHs: '
                do
                   read(5,*,IOSTAT=IST) (NGPRL(I), I=1,NGPR)
                   if (IST > 0) then
                      cycle
                   else if (IST < 0) then
                      cycle OUTER
                   else
                      exit
                   end if
                   if (minval(NGPRL) <= 0 .or. maxval(NGPRL) > NGYRM) then
                      write(6,*) 'XX Invalid variable number! Please start over again.'
                      cycle
                   end if
                end do

                call PAGES
                call SETCHS(0.3, 0.0)
                call SETLIN(0, 1, 7)
                if (MODEG == 2) then
                   IND = 9
                else
                   IND = 0
                end if
                do I = 1, NGPR
                   write(STRL,*) NGPRL(I)
                   STRL = "@"//trim(adjustl(STRL))//"@"
                   call APPROPGY(MODEG, GY%v(:,:,NGPRL(I)), GYL2, size(GYL2), STRL, GY%gnrm(NGPRL(I)))
                   J = I - 1
                   if(NGPR == 4) then
                      call TXGRFRX(J, GX, GYL2, NRMAX, NGR, STRL, MODE, IND)
                   else
                      GPXY_IN(1) =  1.8  + 8.45 * mod(J,3)
                      GPXY_IN(2) =  8.45 + 8.45 * mod(J,3)
                      GPXY_IN(3) = 10.55 - 7.55 * real(J/3)
                      GPXY_IN(4) = 15.1  - 7.55 * real(J/3)
                      call TXGRFRX(J, GX, GYL2, NRMAX, NGR, STRL, MODE, IND, GPXY_IN=GPXY_IN)
                   end if
                end do
                call PAGEE

                deallocate(NGPRL,GYL2,GPXY_IN)
             case('R')
                select case(KID2)
                case('A')
                   call TXGRFR(-1,MODE)
                case('B')
                   call TXGRFR(-4,MODE)
                case('C')
                   call TXGRFR(-7,MODE)
                case('D')
                   call TXGRFR(-13,MODE)
                case('E')
                   call TXGRFR(-17,MODE)
                case DEFAULT
                   read(STR2(2:5),'(I4)',IOSTAT=IST) NGPR
                   if (IST < 0) cycle OUTER
                   if      (NGPR == 0) then
                      cycle OUTER
                   else if (NGPR >= 0 .and. NGPR <= NGPRM) then
                      call TXGRFR(NGPR,MODE)
                   end if
                end select
!!$             case('T')
!!$                select case(KID2)
!!$                case('A')
!!$                   do NGPT = 1, NGPTM
!!$                      call TXGRFT(NGPT,MODE)
!!$                   end do
!!$
!!$                case('B')
!!$                   do NGPT = 1, 6
!!$                      call TXGRFT(NGPT,MODE)
!!$                   end do
!!$                   do NGPT = 9, 10
!!$                      call TXGRFT(NGPT,MODE)
!!$                   end do
!!$
!!$                case DEFAULT
!!$                   read(STR2(2:5),'(I4)',IOSTAT=IST) NGPT
!!$                   if (IST < 0) then
!!$                      write(6,*) '### ERROR : Invalid Command : ', STR
!!$                      cycle
!!$                   end if
!!$                   if (NGPT >= 1 .and. NGPT <= NGPTM) then
!!$                      call TXGRFT(NGPT,MODE)
!!$                   end if
!!$                end select
!!$             case('V')
!!$                select case(KID2)
!!$                case('A')
!!$                   do NGPV = 1, NGPVM
!!$                      call TXGRFV(NGPV,MODE)
!!$                   end do
!!$                case('B')
!!$                   do NGPV = 1, 6
!!$                      call TXGRFV(NGPV,MODE)
!!$                   end do
!!$                   do NGPV = 9, 10
!!$                      call TXGRFV(NGPV,MODE)
!!$                   end do
!!$
!!$                case DEFAULT
!!$                   read(STR2(2:5),'(I4)',IOSTAT=IST) NGPV
!!$                   if (IST < 0) then
!!$                      write(6,*) '### ERROR : Invalid Command : ', STR
!!$                      cycle
!!$                   end if
!!$                   if (NGPV >= 1 .and. NGPV <= NGPVM) then
!!$                      call TXGRFV(NGPV,MODE)
!!$                   end if
!!$                end select
             case('X')
                exit
             case DEFAULT
                write(6,*) 'XX UNKNOWN GRAPHIC COMMAND'
             end select
          end do

       case('H')
          ! *** Animation ***
          select case(KID2)
          case('A')
             call TXGRFRA(-1)
          case('B')
             call TXGRFRA(-3)
          case('C')
             call TXGRFRA(-6)
          case('D')
             call TXGRFRA(-10)
          case('N')
             do
                write(6,*) '## HOW MANY GRAPHs ?: 4 or 6'
                read(5,'(I1)',IOSTAT=IST) NGPR
                if (IST > 0) then
                   cycle
                else if (IST < 0) then
                   cycle OUTER
                else
                   if(NGPR /= 4 .and. NGPR /= 6) cycle
                   exit
                end if
             end do
             allocate(NGPRL(1:NGPR),GYL(0:NRMAX,0:NGT,1:NGPR),STRA(1:NGPR),GMAXA(1:NGPR),GMINA(1:NGPR))
             write(6,'(A,I1,A)') '## CHOOSE ',NGPR, ' GRAPHs: '
             do
                read(5,*,IOSTAT=IST) (NGPRL(I), I=1,NGPR)
                if (minval(NGPRL) <= 0 .or. maxval(NGPRL) > NGYRM) then
                   write(6,*) 'XX Invalid variable number! Please start over again.'
                   cycle
                end if
                if (IST > 0) then
                   cycle
                else if (IST < 0) then
                   cycle OUTER
                else
                   exit
                end if
             end do

!             call GSTITL('//') ! Eliminate header eternally
             call PAGES
             call SETCHS(0.3, 0.0)
             call SETLIN(0, 1, 7)
             call INQFNT(IFNT)
             call SETFNT(44)

             if (MODEG == 2) then
                IND = 9
             else
                IND = 0
             end if
             if(NGPR == 4) then
                J = 5 ; GZ = 8.5
             else
                J = 6 ; GZ = 7.0
             end if
             do I = 1, NGPR
                write(STRA(I),*) NGPRL(I)
                STRA(I) = "@"//trim(adjustl(STRA(I)))//"@"
                call APPROPGY(MODEG, GYT%v(:,:,NGPRL(I)), GYL(:,:,I), size(GYL,1)*size(GYL,2), &
                     &        STRA(I), GYT%gnrm(NGPRL(I)), GMAX=GMAXA(I), GMIN=GMINA(I))
             end do
             do NG = 0, NGT
                call animes
                call gtextx(12.5,GZ,'@T=@',0)
                call gnumbr(13.1,GZ,GTX(NG),3,0)
                do I = 1, NGPR
                   call TXGRFRS(I-1, GX, GYL(0:NRMAX,NG:NG,I), NRMAX, 1, STRA(I), 0, IND, 0, J, &
                        &       'ANIME', GMAXA(I), GMINA(I))
                end do
                call animee
             end do
             call SETFNT(IFNT)
             call PAGEE
             deallocate(NGPRL,GYL,STRA,GMAXA,GMINA)

          case DEFAULT
             read(STR(2:5),*,IOSTAT=IST) NGYR
             if (IST /= 0) then
                write(6,*) '### ERROR : Invalid Command : ', STR
                cycle OUTER
             end if
             if (NGYR >= 0 .and. NGYR <= NGYRM) then
                IND = 0
                if (MODEG == 2) IND = 9
                allocate(GYL2(0:NRMAX,0:NGT))
                STRL = '@profile'//STR(2:5)//'@'
                call APPROPGY(MODEG, GYT%v(:,:,NGYR), GYL2, size(GYL2), STRL, GYT%gnrm(NGYR), &
                     &        GMAX=GMAX, GMIN=GMIN)

                call pages

                call SETLIN(0, 1, 7)
                call INQFNT(IFNT)
                call SETFNT(44)

                do ng = 0, ngt
                   call animes
                   call gtextx(12.5,17.7,'@T=@',0)
                   call gnumbr(13.1,17.7,GTX(NG),3,0)
                   call TXGRFRS(0, GX, GYL2(0:NRMAX,NG:NG), NRMAX, 1, STRL, 0, IND, 0, 3, &
                        &       'ANIME', GMAX, GMIN)
                   call animee
                end do

                call SETFNT(IFNT)
                call pagee
                deallocate(GYL2)
             end if
          end select

       case('X')
          ! *** Exit ***
          return

       case DEFAULT
          write(6,*) 'XX UNKNOWN GRAPHIC COMMAND'
       end select
    end do OUTER

  end subroutine TXGOUT
#endif
  !***************************************************************
  !
  !   Save graphic data
  !
  !***************************************************************

  subroutine TX_GRAPH_SAVE

    use tx_commons, only : T_TX

    !  Define radial coordinate for graph

    call TXPRFG

    !  Store center or edge values of variables for showing time-evolution graph

    call TXSTGT(real(T_TX))

    !  Store global quantities for showing time-evolution graph

    call TXSTGV(real(T_TX))

    !  Store profile data for showing graph

    call TXSTGR(NGR,GT,GY,NGRM)

    !  Store balance profile data for showing graph

    if(T_TX /= 0.d0) call TXSTGQ

  end subroutine TX_GRAPH_SAVE

  !***************************************************************
  !
  !   Initialize graphic axis
  !
  !***************************************************************

  subroutine TXPRFG

    use tx_commons, only : NRMAX, RHO

    !  GX(NR) : Integer

    GX(0:NRMAX) = real(RHO(0:NRMAX))

  end subroutine TXPRFG

  !***************************************************************
  !
  !   Store GY
  !
  !***************************************************************

  subroutine TXSTGR(NG,GTL,GYL,NGM)

    use tx_commons

    integer(4), intent(inout) :: NG
    integer(4), intent(in) :: NGM
    real(4), dimension(0:NGM), intent(out), optional :: GTL
    type(graphic_data) :: GYL

    integer(4) :: NR, MDLNEOL, i, j, k

    MDLNEOL = mod(MDLNEO,10)

    if(present(GTL)) then
       if (NG < NGM) NG = NG + 1

       GTL(NG) = real(T_TX)
    end if

    do NR = 0, NRMAX

       !  *** Dependent variables ***************************************

       GYL%v(NR,NG,  1) = real(X(NR,LQm1)) ; GYL%gnrm(  1) = gkilo
       GYL%v(NR,NG,  2) = real(X(NR,LQm2)) ; GYL%gnrm(  2) = 1.e-3
       GYL%v(NR,NG,  3) = real(X(NR,LQm3))
       GYL%v(NR,NG,  4) = real(X(NR,LQm4))
       GYL%v(NR,NG,  5) = real(X(NR,LQm5)) ; GYL%gnrm(  5) = gkilo
       GYL%v(NR,NG,  6) = real(X(NR,LQe1))
       GYL%v(NR,NG,  7) = real(X(NR,LQe2))
       GYL%v(NR,NG,  8) = real(X(NR,LQe3)) ; GYL%gnrm(  8) = gkilo
       GYL%v(NR,NG,  9) = real(X(NR,LQe4)) ; GYL%gnrm(  9) = gkilo
       GYL%v(NR,NG, 10) = real(X(NR,LQe5))
       GYL%v(NR,NG, 11) = real(X(NR,LQe6)) ; GYL%gnrm( 11) = gkilo
       GYL%v(NR,NG, 12) = real(X(NR,LQe7)) ; GYL%gnrm( 12) = gkilo
       GYL%v(NR,NG, 13) = real(X(NR,LQe8)) ; GYL%gnrm( 13) = gkilo
       GYL%v(NR,NG, 14) = real(X(NR,LQi1))
       GYL%v(NR,NG, 15) = real(X(NR,LQi2))
       GYL%v(NR,NG, 16) = real(X(NR,LQi3)) ; GYL%gnrm( 16) = gkilo
       GYL%v(NR,NG, 17) = real(X(NR,LQi4)) ; GYL%gnrm( 17) = gkilo
       GYL%v(NR,NG, 18) = real(X(NR,LQi5))
       GYL%v(NR,NG, 19) = real(X(NR,LQi6)) ; GYL%gnrm( 19) = gkilo
       GYL%v(NR,NG, 20) = real(X(NR,LQi7)) ; GYL%gnrm( 20) = gkilo
       GYL%v(NR,NG, 21) = real(X(NR,LQi8)) ; GYL%gnrm( 21) = gkilo
       GYL%v(NR,NG, 22) = real(X(NR,LQz1))
       GYL%v(NR,NG, 23) = real(X(NR,LQz2))
       GYL%v(NR,NG, 24) = real(X(NR,LQz3)) ; GYL%gnrm( 24) = gkilo
       GYL%v(NR,NG, 25) = real(X(NR,LQz4)) ; GYL%gnrm( 25) = gkilo
       GYL%v(NR,NG, 26) = real(X(NR,LQz5))
       GYL%v(NR,NG, 27) = real(X(NR,LQz6)) ; GYL%gnrm( 27) = gkilo
       GYL%v(NR,NG, 28) = real(X(NR,LQz7)) ; GYL%gnrm( 28) = gkilo
       GYL%v(NR,NG, 29) = real(X(NR,LQz8)) ; GYL%gnrm( 29) = gkilo
       GYL%v(NR,NG, 30) = real(X(NR,LQb1)) ; GYL%gnrm( 30) = 1.d-2
       GYL%v(NR,NG, 31) = real(X(NR,LQb2))
       GYL%v(NR,NG, 32) = real(X(NR,LQb3)) ; GYL%gnrm( 32) = gkilo
       GYL%v(NR,NG, 33) = real(X(NR,LQb4)) ; GYL%gnrm( 33) = gkilo
       GYL%v(NR,NG, 34) = real(X(NR,LQb7)) ; GYL%gnrm( 34) = gkilo
       GYL%v(NR,NG, 35) = real(X(NR,LQb8)) ; GYL%gnrm( 35) = gkilo
       GYL%v(NR,NG, 36) = real(X(NR,LQn1)) ; GYL%gnrm( 36) = 1.e-4
       GYL%v(NR,NG, 37) = real(X(NR,LQn2)) ; GYL%gnrm( 37) = 1.e-6
       GYL%v(NR,NG, 38) = real(X(NR,LQn3)) ; GYL%gnrm( 38) = 1.e-6
       GYL%v(NR,NG, 39) = real(X(NR,LQnz)) ; GYL%gnrm( 39) = 1.e-4
       GYL%v(NR,NG, 40) = real(X(NR,LQr1)) ; GYL%gnrm( 40) = 1.e-4

       ! 41-51: spare for future additions

       !  *** Fluid moment variables ************************************

       do i = 1, NSM
          j = (i - 1) * 13
          GYL%v(NR,NG, 52+j) = real(Var(NR,i)%n*1.d20)      ; GYL%gnrm(52+j) = 1.e20
          GYL%v(NR,NG, 53+j) = real(Var(NR,i)%T)
          GYL%v(NR,NG, 54+j) = real(Var(NR,i)%p*1.d20*rKeV) ; GYL%gnrm(54+j) = gmega
          GYL%v(NR,NG, 55+j) = real(Var(NR,i)%UrV)
          GYL%v(NR,NG, 56+j) = real(Var(NR,i)%BUpar)        ; GYL%gnrm(56+j) = gkilo
          GYL%v(NR,NG, 57+j) = real(Var(NR,i)%Bqpar)        ; GYL%gnrm(57+j) = gkilo
          GYL%v(NR,NG, 58+j) = real(Var(NR,i)%Uthhat)       ; GYL%gnrm(58+j) = gkilo
          GYL%v(NR,NG, 59+j) = real(Var(NR,i)%UphR)         ; GYL%gnrm(59+j) = gkilo
          GYL%v(NR,NG, 60+j) = real(Var(NR,i)%RUph)         ; GYL%gnrm(60+j) = gkilo
          GYL%v(NR,NG, 61+j) = real(Var(NR,i)%Ur)
          GYL%v(NR,NG, 62+j) = real(Var(NR,i)%Uth)          ; GYL%gnrm(62+j) = gkilo
          GYL%v(NR,NG, 63+j) = real(Var(NR,i)%Uph)          ; GYL%gnrm(63+j) = gkilo
          GYL%v(NR,NG, 64+j) = real(Var(NR,i)%BV1)          ; GYL%gnrm(64+j) = gkilo
       end do
       GYL%gnrm( 78) = 1.e18
       GYL%gnrm( 80) = gkilo

       !  *** Electromagnetic quantities ********************************

       GYL%v(NR,NG, 91) = real(ErV(NR))    ; GYL%gnrm( 91) = gkilo
       GYL%v(NR,NG, 92) = real(BEpol(NR))  ; GYL%gnrm( 92) = 1.e-3
       GYL%v(NR,NG, 93) = real(Etor(NR))
       GYL%v(NR,NG, 94) = real(BthV(NR))
       GYL%v(NR,NG, 95) = real(fipol(NR))
       GYL%v(NR,NG, 96) = real(BEpara(NR))
       GYL%v(NR,NG, 97) = real(Q(NR))
       GYL%v(NR,NG, 98) = real(S(NR))

       GYL%v(NR,NG, 99) = real(AJ(NR))
       GYL%v(NR,NG,100) = real(AJOH(NR))
       GYL%v(NR,NG,101) = real(AJBS(NR))
       GYL%v(NR,NG,102) = real(AJNB(NR))
       GYL%v(NR,NG,103) = real(AJRF(NR))
       GYL%v(NR,NG,104) = real(BJPARA(NR))
       GYL%v(NR,NG,105) = real(BJOH(NR))
       GYL%v(NR,NG,106) = real(BJBS(NR))
       GYL%v(NR,NG,107) = real(BJNB(NR))
       GYL%v(NR,NG,108) = real(BJRF(NR))
       GYL%gnrm(99:108) = gmega
       GYL%v(NR,NG,109) = real(aee*achg(1)*Var(NR,1)%n*Var(NR,1)%UphR*1.d20/ait(NR))
       GYL%gnrm(109) = gmega
       GYL%v(NR,NG,110) = real(aee*achg(2)*Var(NR,2)%n*Var(NR,2)%UphR*1.d20/ait(NR))
       GYL%v(NR,NG,111) = real(aee*achg(3)*Var(NR,3)%n*Var(NR,3)%UphR*1.d20/ait(NR))
       GYL%v(NR,NG,112) = real(aee*achgb*PNbV(NR)*(aat(NR)*fipol(NR)/bbt(NR)*BUbparV(NR)) &
            &                 *1.d20/ait(NR))
       GYL%gnrm(110:112) = gkilo

       ! *** Rho vs Psi *************************************************

       GYL%v(NR,NG,113) = real((Psiv(NR)-Psiv(0))/(Psiv(NRA)-Psiv(0)))

       ! *** Beam particles *********************************************

       GYL%v(NR,NG,114) = real(PNbV(NR) * 1.d20)   ; GYL%gnrm(114) = 1.e18
       GYL%v(NR,NG,115) = real(UbrVV(NR))
       GYL%v(NR,NG,116) = real(BUbparV(NR))
       GYL%v(NR,NG,117) = real(RUbphV(NR))
       GYL%v(NR,NG,118) = real(PTbV(NR))
       GYL%v(NR,NG,119) = real(UbphVR(NR))
       GYL%v(NR,NG,120) = real(UbphV(NR))
       GYL%v(NR,NG,121) = real(BVbdiag(NR))
       GYL%gnrm(116:121) = gkilo
       GYL%v(NR,NG,122) = real(PNbrpV(NR) * 1.d20) ; GYL%gnrm(122) = 1.e15
       GYL%v(NR,NG,123) = real(PNbrpV(NR) * PNbVinv(NR))
       GYL%v(NR,NG,124) = real(UbrV(NR))

       ! *** Neutrals ***************************************************
       
       GYL%v(NR,NG,125) = real(PN01V(NR)*1.d20) ; GYL%gnrm(125) = 1.e15
       GYL%v(NR,NG,126) = real(PN02V(NR)*1.d20) ; GYL%gnrm(126) = 1.e13
       GYL%v(NR,NG,127) = real(PN03V(NR)*1.d20) ; GYL%gnrm(127) = 1.e13
       GYL%v(NR,NG,128) = real(PN0zV(NR)*1.d20) ; GYL%gnrm(128) = 1.e15

       ! *** Perpendicular flows ****************************************

       GYL%v(NR,NG,129) = real(Var(NR,1)%RUph/fipol(NR)-Var(NR,1)%BUpar/bbt(NR))
       GYL%gnrm(129) = gkilo
       GYL%v(NR,NG,130) = real(Var(NR,2)%RUph/fipol(NR)-Var(NR,2)%BUpar/bbt(NR))
       GYL%gnrm(130) = gkilo
       GYL%v(NR,NG,131) = real(Var(NR,3)%RUph/fipol(NR)-Var(NR,3)%BUpar/bbt(NR))
       GYL%gnrm(131) = gkilo

       ! *** Diamagnetic flows ******************************************

       do i = 1, NSM
          j = (i - 1) * 2
          GYL%v(NR,NG,132+j) = real(BVsdiag(NR,i,1)) ; GYL%gnrm(132+j) = gkilo
          GYL%v(NR,NG,133+j) = real(BVsdiag(NR,i,2)) ; GYL%gnrm(133+j) = gkilo
       end do

       !  *** Quasineutrality and ambipolarity **************************

       GYL%v(NR,NG,138) = real(qneut(NR) * 1.d20) ; GYL%gnrm(138) = 1.e14 ! quasineutrality
       GYL%v(NR,NG,139) = real(aee*1.d20*sum(achg(:)*Var(NR,:)%n*Var(NR,:)%Ur))
       GYL%v(NR,NG,140) = real(aee*sum(achg(:)*Var(NR,:)%n*Var(NR,:)%UrV)*1.d20) ! sum_s e_s<Gamma_s.grad V>

       !  *** Diffusion coefficients ************************************

       GYL%v(NR,NG,141) = real(sum(Dfs(NR,:)))
       GYL%v(NR,NG,142) = real(Deff(NR,1))
       GYL%v(NR,NG,143) = real(Deff(NR,2))
       GYL%v(NR,NG,144) = real(Deff(NR,3))

       do i = 1, NSM
          j = (i - 1) * 4
          GYL%v(NR,NG,145+j) = real(Chis(NR,i))
          GYL%v(NR,NG,146+j) = real(rMus(NR,i))
          !     ** Effective neoclassical thermal diffusivity
          GYL%v(NR,NG,147+j) = real(ChiNCp(NR,i)+ChiNCt(NR,i))
          !     ** Total thermal diffusivity
          GYL%v(NR,NG,148+j) = real(Chis(NR,i)+ChiNCp(NR,i)+ChiNCt(NR,i))
       end do

       GYL%v(NR,NG,157) = real(D01(NR)) ; GYL%gnrm(157) = gkilo
       GYL%v(NR,NG,158) = real(D02(NR)) ; GYL%gnrm(158) = gkilo
       GYL%v(NR,NG,159) = real(D03(NR)) ; GYL%gnrm(159) = gkilo
       GYL%v(NR,NG,160) = real(D0z(NR)) ; GYL%gnrm(160) = gkilo

       ! 160: spare for future additions

       ! *** Pinch velocities *******************************************

       GYL%v(NR,NG,161) = real(VWpch(NR))
       GYL%v(NR,NG,162) = real(Vmps(NR,1))
       GYL%v(NR,NG,163) = real(Vmps(NR,2))
       GYL%v(NR,NG,164) = real(Vmps(NR,3))
       GYL%v(NR,NG,165) = real(Vhps(NR,1))
       GYL%v(NR,NG,166) = real(Vhps(NR,2))
       GYL%v(NR,NG,167) = real(Vhps(NR,3))
       
       ! *** Drag force coefficients ************************************

       GYL%v(NR,NG,168) = real(rNuL(NR))
       GYL%v(NR,NG,169) = real(rNuLB(NR))
       GYL%v(NR,NG,170) = real(rNuiCX(NR))
       GYL%v(NR,NG,171) = real(rNuION(NR))
       GYL%v(NR,NG,172) = real(rNuOL(NR))
       GYL%v(NR,NG,173) = real(rNu0b(NR))
       GYL%v(NR,NG,174) = real(rNuB(NR))
       GYL%v(NR,NG,175) = real(rNuTei(NR))
       GYL%v(NR,NG,176) = real(rNuTez(NR))
       GYL%v(NR,NG,177) = real(rNuTiz(NR))
       do i = 1, NSM
          j = (i - 1) * 3
          GYL%v(NR,NG,178+j) = real(rNu0s(NR,i))
          GYL%v(NR,NG,179+j) = real(rNuLTs(NR,i))
          GYL%v(NR,NG,180+j) = real(rNuAss(NR,i))
       end do
       GYL%gnrm(179) = gkilo

       ! 187-190: spare for future additions

       ! *** Other coefficients *****************************************
       
       GYL%v(NR,NG,191) = real(WPM(NR))
       !     ** ExB shearing rate
       GYL%v(NR,NG,192) = real(wexb(NR)) ; GYL%gnrm(192) = gmega
       !     ** CDBM
       GYL%v(NR,NG,193) = real(rG1h2(NR))
       GYL%v(NR,NG,194) = real(FCDBM(NR))
       GYL%v(NR,NG,195) = real(Alpha(NR))
       GYL%v(NR,NG,196) = real(rKappa(NR))
       !     ** Ripple loss part
       GYL%v(NR,NG,197) = real(rNubrp1(NR))
       GYL%v(NR,NG,198) = real(DltRP(NR))
       GYL%v(NR,NG,199) = real(Ubrp(NR))
       GYL%v(NR,NG,200) = real(Dbrp(NR))
       GYL%v(NR,NG,201) = real(rNubL(NR))
       !     ** Effective charge
       GYL%v(NR,NG,202) = real(Zeff(NR))

       ! 203-207: spare for future additions

       ! *** Neoclassical friction forces and viscosities ***************

       !      ** thermal - thermal
       do i = 1, NSM
          j = (i - 1) * 4
          GYL%v(NR,NG,208+j) = real(xmu(NR,i,1,1)) ; GYL%gnrm(208+j) = 1.e-6
          GYL%v(NR,NG,209+j) = real(xmu(NR,i,1,2)) ; GYL%gnrm(209+j) = 1.e-6
          GYL%v(NR,NG,210+j) = real(xmu(NR,i,2,1)) ; GYL%gnrm(210+j) = 1.e-6
          GYL%v(NR,NG,211+j) = real(xmu(NR,i,2,2)) ; GYL%gnrm(211+j) = 1.e-6
       end do

       do i = 1, NSM
          do j = 1, NSM
             k = (j - 1) * 4 + (i - 1) * NSM * 4
             GYL%v(NR,NG,220+k) = real(lab(NR,i,j,1,1)) ; GYL%gnrm(220+j) = gkilo
             GYL%v(NR,NG,221+k) = real(lab(NR,i,j,1,2)) ; GYL%gnrm(221+j) = gkilo
             GYL%v(NR,NG,222+k) = real(lab(NR,i,j,2,1)) ; GYL%gnrm(222+j) = gkilo
             GYL%v(NR,NG,223+k) = real(lab(NR,i,j,2,2)) ; GYL%gnrm(223+j) = gkilo
          end do
       end do
       GYL%gnrm(232:255) = 1.0

       !      ** beam - thermal
       GYL%v(NR,NG,256) = real(xmuf(NR,1))     ! fast-ion viscosity
       GYL%gnrm(256) = 1.e20 
       GYL%v(NR,NG,257) = real(laf(NR,1,1,1))  ! ef 11
       GYL%v(NR,NG,258) = real(laf(NR,1,2,1))  ! ef 21
       GYL%v(NR,NG,259) = real(lfb(NR,1,1,1))  ! fe 11, not normalized
       GYL%v(NR,NG,260) = real(laf(NR,2,1,1))  ! if 11
       GYL%v(NR,NG,261) = real(lfb(NR,2,1,1))  ! fi 11, not normalized
       GYL%v(NR,NG,262) = real(laf(NR,3,1,1))  ! zf 11
       GYL%v(NR,NG,263) = real(lfb(NR,3,1,1))  ! fz 11, not normalized
       GYL%v(NR,NG,264) = real(lff(NR,1,1))    ! ff 11, not normalized

       ! 265-267: spare for future additions

       ! *** Sources ****************************************************

       !      ** Particles
       GYL%v(NR,NG,268) = real(SNB(NR))
       GYL%v(NR,NG,269) = real(SNBPDi(NR))
       GYL%v(NR,NG,270) = real(SNBTGi(NR))
       GYL%v(NR,NG,271) = real(SIE(NR))  ; GYL%gnrm(271) = 1.e20
       GYL%v(NR,NG,272) = real(SCX(NR))  ; GYL%gnrm(272) = gmega
       GYL%v(NR,NG,273) = real(SiLC(NR))
       !      ** Heat
       GYL%v(NR,NG,274) = real(PNBe(NR))   ; GYL%gnrm(274) = gmega
       GYL%v(NR,NG,275) = real(PNBi(NR))   ; GYL%gnrm(275) = gmega
       GYL%v(NR,NG,276) = real(PNBz(NR))   ; GYL%gnrm(276) = gmega
       GYL%v(NR,NG,277) = real(PRFe(NR))   ; GYL%gnrm(277) = gmega
       GYL%v(NR,NG,278) = real(PRFi(NR))   ; GYL%gnrm(278) = gmega
       GYL%v(NR,NG,279) = real(PRFz(NR))   ; GYL%gnrm(279) = gmega
       GYL%v(NR,NG,280) = real(PALFe(NR))  ; GYL%gnrm(280) = gmega
       GYL%v(NR,NG,281) = real(PALFi(NR))  ; GYL%gnrm(281) = gmega
       GYL%v(NR,NG,282) = real(PALFz(NR))  ; GYL%gnrm(282) = gmega
       GYL%v(NR,NG,283) = real(POHs(NR,1)) ; GYL%gnrm(283) = gmega
       GYL%v(NR,NG,284) = real(POHs(NR,2))
       GYL%v(NR,NG,285) = real(POHs(NR,3))
       GYL%v(NR,NG,286) = real(PEQei(NR))  ; GYL%gnrm(286) = gmega
       GYL%v(NR,NG,287) = real(PEQez(NR))  ; GYL%gnrm(287) = gmega
       GYL%v(NR,NG,288) = real(PEQiz(NR))  ; GYL%gnrm(288) = gmega
       GYL%v(NR,NG,289) = real(PIE(NR))
       GYL%v(NR,NG,290) = real(PCX(NR))    ; GYL%gnrm(290) = gkilo
       GYL%v(NR,NG,291) = real(PBr(NR))
       GYL%v(NR,NG,292) = real(PNB(NR))    ; GYL%gnrm(292) = gmega
       GYL%v(NR,NG,293) = real(PNBPD(NR)/(Eb*rKeV*1.d20))
       GYL%v(NR,NG,294) = real(PNBTG(NR)/(Eb*rKeV*1.d20))
       GYL%v(NR,NG,295) = real(Vbpara(NR)) ; GYL%gnrm(292) = gkilo
       GYL%v(NR,NG,296) = real(Var(NR,2)%Uthhat*BnablaPi(NR,2))
       GYL%v(NR,NG,297) = real(Var(NR,3)%Uthhat*BnablaPi(NR,3))
       !      ** Torque
       !       External toroidal NBI torque density
       GYL%v(NR,NG,298) = real(BSmb(NR)*fipol(NR)/bbt(NR)*1.d20)
       !       Generated toroidal torque density
       GYL%v(NR,NG,299) = real(sdt(NR)*aee*1.d20*sum(achg(:)*Var(NR,:)%n*Var(NR,:)%UrV)) ! <j.grad psi>
       !       Additional torque
       GYL%v(NR,NG,300) = real(Tqt(NR))
       ! sum_s e_s<Gamma_s.grad V>
       GYL%v(NR,NG,301) = real(aee*( sum(achg(:)*Var(NR,:)%n*Var(NR,:)%UrV) &
            &                       +achgb*PNbV(NR)*UbrVV(NR))*1.d20)
       GYL%v(NR,NG,302) = real(aee*achgb*PNbV(NR)*UbrVV(NR)*1.d20) ! e_b<Gamma_b.grad V>

       ! 303-308: spare for future additions

       ! *** Equilibrium quantities (metrics) ***************************

       GYL%v(NR,NG,309) = real(vlt(NR))
       GYL%v(NR,NG,310) = real(vro(NR))
       GYL%v(NR,NG,311) = real(sdt(NR))
       GYL%v(NR,NG,312) = real(epst(NR))
       GYL%v(NR,NG,313) = real(aat(NR))
       GYL%v(NR,NG,314) = real(rrt(NR))
       GYL%v(NR,NG,315) = real(ckt(NR))
       GYL%v(NR,NG,316) = real(suft(NR))
       GYL%v(NR,NG,317) = real(sst(NR))
       GYL%v(NR,NG,318) = real(vro(NR))
       GYL%v(NR,NG,319) = real(vlt(NR))
       GYL%v(NR,NG,320) = real(art(NR))
       GYL%v(NR,NG,321) = real(ait(NR))
       GYL%v(NR,NG,322) = real(bit(NR))
       GYL%v(NR,NG,323) = real(bbrt(NR))
       GYL%v(NR,NG,324) = real(elip(NR))
       GYL%v(NR,NG,325) = real(trig(NR))
       GYL%v(NR,NG,326) = real(sqrt(bbt(NR)))
       GYL%v(NR,NG,327) = real(ft(NR))
       GYL%v(NR,NG,328) = real(UgV(NR))
       GYL%v(NR,NG,329) = real(qhatsq(NR))

       ! 330-337: spare for future additions


       ! *** Particle fluxes ******************************************** 

       if(NR == 0) then
          GYL%v(NR,NG,338:341) = 0.0
       else
          !      ** Particle flux <Gamma_s . grad rho>
          GYL%v(NR,NG,338) = real(Var(NR,1)%n*Var(NR,1)%UrV/vro(NR))*1.e20
          GYL%v(NR,NG,339) = real(Var(NR,2)%n*Var(NR,2)%UrV/vro(NR))*1.e20
          GYL%v(NR,NG,340) = real(Var(NR,3)%n*Var(NR,3)%UrV/vro(NR))*1.e20
          !      ** Particle flux <Gamma_b . grad rho>
          GYL%v(NR,NG,341) = real(PNbV(NR)*UbrVV(NR)/vro(NR))*1.e20
          GYL%gnrm(338:341) = 1.e20
       end if

       ! *** Quantities computed by neoclassical transport solver ******* 

       do i = 1, NSM
          j = (i - 1) * 2
          ! Neo. solver: parallel flows
          GYL%v(NR,NG,342+j) = real(BusparNCL(NR,i,1)) ; GYL%gnrm(342+j) = gkilo
          ! Neo. solver: parallel heat flows
          GYL%v(NR,NG,343+j) = real(BusparNCL(NR,i,2)) ; GYL%gnrm(343+j) = gkilo
       end do
       GYL%v(NR,NG,348) = real(gflux(NR,1,MDLNEOL))
       GYL%v(NR,NG,349) = real(gflux(NR,2,MDLNEOL))
       GYL%v(NR,NG,350) = real(gflux(NR,3,MDLNEOL))

       ! *** Ionization and recombination *******************************

       GYL%v(NR,NG,351) = real(SiVa6A(NR))
       GYL%v(NR,NG,352) = real(SiVsefA(NR))

    end do

  end subroutine TXSTGR

  !***************************************************************
  !
  !   Store GVY
  !
  !***************************************************************

  subroutine TXSTGV(GTIME)

    use tx_commons, only : NRMAX, NRA, NRC, aee, PI, rKeV, achg, achgb &
         &               , Var, qneut, PNbV, UbrVV, PNbrpV, ErV, BthV, Etor, UbphV &
         &               , PN01V, PN02V, PN03V, BEpol, BUbparV, Q, Dfs &
         &               , rG1h2, FCDBM, S, Alpha, rKappa, rNuION &
         &               , Chis, PIE, PCX, SIE, PBr, Deff &
         &               , AJ, ait, aat, fipol, bbt, BUbparV, qhatsq
    use tx_interface, only : rLINEAVE
    use tx_core_module, only : sub_intg_vol

    real(4), INTENT(IN) :: GTIME
    real(8) :: PNESUM1, PNESUM2
    real(8), dimension(:), allocatable :: PNeION

    if (NGVV < NGVM) NGVV=NGVV+1

    GVX(NGVV) = GTIME

    GVY(NGVV,1)  = real(qneut(0) * 1.d20) ! quasineutrality
    GVY(NGVV,2)  = real(aee*sum(achg(:)*Var(NRC,:)%n*Var(NRC,:)%UrV)*1.d20) ! ambipolarity: e_s<Gamma_s.grad V>
    GVY(NGVV,3)  = GVY(NGVV,2) + real(achgb*aee*PNbV(NRC)*UbrVV(NRC)*1.d20)

    GVY(NGVV,4)  = real(Var(  0,1)%n * 1.d20)
    GVY(NGVV,5)  = real(Var(NRC,1)%Ur)
    GVY(NGVV,6)  = real(Var(NRC,1)%Uth)
    GVY(NGVV,7)  = real(Var(  0,1)%Uph)
    GVY(NGVV,8)  = real(Var(  0,1)%T)
    GVY(NGVV,9)  = real(Var(NRC,1)%Uthhat)
    GVY(NGVV,10) = real(Var(NRC,1)%RUph/fipol(NRC)-Var(NRC,1)%BUpar/bbt(NRC))
    GVY(NGVV,11) = real(Var(NRC,1)%BUpar)
    GVY(NGVV,12) = real(Var(0,1)%p * 1.d20 * rKeV)

    GVY(NGVV,13) = real(Var(  0,2)%n * 1.d20)
    GVY(NGVV,14) = real(Var(NRC,2)%Ur)
    GVY(NGVV,15) = real(Var(NRC,2)%Uth)
    GVY(NGVV,16) = real(Var(NRC,2)%Uph)
    GVY(NGVV,17) = real(Var(  0,2)%T)
    GVY(NGVV,18) = real(Var(NRC,2)%Uthhat)
    GVY(NGVV,19) = real(Var(NRC,2)%RUph/fipol(NRC)-Var(NRC,2)%BUpar/bbt(NRC))
    GVY(NGVV,20) = real(Var(NRC,2)%BUpar)
    GVY(NGVV,21) = real(Var(  0,2)%p * 1.d20 * rKeV)

    GVY(NGVV,22) = real(Var(  0,3)%n * 1.d20)
    GVY(NGVV,23) = real(Var(NRC,3)%Ur)
    GVY(NGVV,24) = real(Var(NRC,3)%Uth)
    GVY(NGVV,25) = real(Var(NRC,3)%Uph)
    GVY(NGVV,26) = real(Var(  0,3)%T)
    GVY(NGVV,27) = real(Var(NRC,3)%Uthhat)
    GVY(NGVV,28) = real(Var(NRC,3)%RUph/fipol(NRC)-Var(NRC,3)%BUpar/bbt(NRC))
    GVY(NGVV,29) = real(Var(NRC,3)%BUpar)
    GVY(NGVV,30) = real(Var(  0,3)%p * 1.d20 * rKeV)

    GVY(NGVV,31) = real(ErV(NRC))
    GVY(NGVV,32) = real(Etor(NRMAX))
    GVY(NGVV,33) = real(BEpol(NRC))
    GVY(NGVV,34) = real(BthV(NRA))
    GVY(NGVV,35) = real(fipol(0))
    GVY(NGVV,36) = real(Q(0))

    GVY(NGVV,37) = real(AJ(0))
    GVY(NGVV,38) = real(aee*achg(1)*Var(0,1)%n*Var(0,1)%UphR*1.d20/ait(0))
    GVY(NGVV,39) = real(aee*achg(2)*Var(0,2)%n*Var(0,2)%UphR*1.d20/ait(0))
    GVY(NGVV,40) = real(aee*achg(3)*Var(0,3)%n*Var(0,3)%UphR*1.d20/ait(0))
    GVY(NGVV,41) = real(aee*achgb*PNbV(0)*(aat(0)*fipol(0)/bbt(0)*BUbparV(0))*1.d20/ait(0))

    GVY(NGVV,42) = real(PNbV(0) * 1.d20)
    GVY(NGVV,43) = real(UbphV(0))
    GVY(NGVV,44) = real(BUbparV(NRC))
    GVY(NGVV,45) = real(PNbrpV(0) * 1.d20)

    GVY(NGVV,46) = real( PN01V(NRA) * 1.d20)
    GVY(NGVV,47) = real( PN02V(0)   * 1.d20)
    GVY(NGVV,48) = real( PN03V(0)   * 1.d20)
    GVY(NGVV,49) = real((PN01V(NRA) + PN02V(NRA) + PN03V(NRA)) * 1.d20)

    GVY(NGVV,50) = real(sum(Dfs(NRC,:)))
    GVY(NGVV,51) = real(Deff(NRC,1))
    GVY(NGVV,52) = real(Chis(NRC,1))
    GVY(NGVV,53) = real(Chis(NRC,2))
    GVY(NGVV,54) = real(Chis(NRC,3))

    GVY(NGVV,55) = real(rG1h2(NRC))
    GVY(NGVV,56) = real(FCDBM(NRC))
    GVY(NGVV,57) = real(S(NRC))
    GVY(NGVV,58) = real(Alpha(NRC))
    GVY(NGVV,59) = real(rKappa(NRC))

    GVY(NGVV,60) = real(PIE(NRC))
    GVY(NGVV,61) = real(PCX(NRC))
    GVY(NGVV,62) = real(SIE(NRC))
    GVY(NGVV,63) = real(PBr(NRC))

    GVY(NGVV,64)  = real(rLINEAVE(0.d0))
    GVY(NGVV,65)  = real(rLINEAVE(0.24d0))
    GVY(NGVV,66)  = real(rLINEAVE(0.6d0))

    call sub_intg_vol(Var(0:NRMAX,1)%n,NRA,PNESUM1)
    allocate(PNeION,mold=rNuION)
    PNeION(0:NRMAX) = Var(0:NRMAX,1)%n*rNuION(0:NRMAX)
    call sub_intg_vol(PNeION,NRA,PNESUM2)
    deallocate(PNeION)
    GVY(NGVV,67) = real(PNESUM1)
    GVY(NGVV,68) = real(PNESUM2)
    if(NGVV == 0.or.ABS(PNESUM2) <= 0.D0) then
       GVY(NGVV,69) = 0.0
    else
       GVY(NGVV,69) = real(PNESUM1/PNESUM2)
    end if

    GVY(NGVV,70) = 1.0+2.0*real(qhatsq(NRC))

  end subroutine TXSTGV

  !***************************************************************
  !
  !   Store GTY
  !
  !***************************************************************

  subroutine TXSTGT(GTIME)

    use tx_commons, only : TS0, TSAV, PINT, POHT, PNBT, PRFT, PRFTe, PRFTi,PNFT, &
         &              AJT, AJOHT, AJNBT, AJBST, POUT, PCXT, PIET, QF, ANS0, &
         &              ANSAV, WPT, WBULKT, WST, TAUE1, TAUE2, TAUEP, BETAA, &
         &              BETA0, BETAPA, BETAP0, VLOOP, ALI, Q, RQ1, ANF0, ANFAV, &
         &              VOLAVN, TAUP, TAUPA, Gamma_a, TNBcol, TTqt, PnumN0, &
         &              CPsi, VPoynt, VLOOPpol, CEjima, PoyntS, PoyntI, PoyntR, &
         &              totmnRV, totjxB
    real(4), INTENT(IN) :: GTIME

    if (NGT < NGTM) NGT=NGT+1

    GTX(NGT) = GTIME

    GTY(NGT,1)  = real(TS0(1))
    GTY(NGT,2)  = real(TS0(2))
    GTY(NGT,3)  = real(TS0(3))
    GTY(NGT,4)  = real(TSAV(1))
    GTY(NGT,5)  = real(TSAV(2))
    GTY(NGT,6)  = real(TSAV(3))
    GTY(NGT,7)  = real(ANS0(1))
    GTY(NGT,8)  = real(ANS0(2))
    GTY(NGT,9)  = real(ANS0(3))
    GTY(NGT,10) = real(ANSAV(1))
    GTY(NGT,11) = real(ANSAV(2))
    GTY(NGT,12) = real(ANSAV(3))

    GTY(NGT,13) = real(PINT)
    GTY(NGT,14) = real(POHT)
    GTY(NGT,15) = real(PNBT)
    GTY(NGT,16) = real(PRFT)
    GTY(NGT,17) = real(PNFT)
    GTY(NGT,18) = real(POUT)
    GTY(NGT,19) = real(PCXT)
    GTY(NGT,20) = real(PIET)
    GTY(NGT,21) = real(PRFTe)
    GTY(NGT,22) = real(PRFTi)
    GTY(NGT,23) = real(QF)
    !  GTY(NGT,24) = real(PRLT)
    !  GTY(NGT,25) = real(PCONT)

    GTY(NGT,26) = real(AJT)
    GTY(NGT,27) = real(AJOHT)
    GTY(NGT,28) = real(AJNBT)
    GTY(NGT,29) = real(AJBST)
    GTY(NGT,30) = real(AJOHT+AJBST+AJNBT)

    GTY(NGT,31) = real(WPT)
    GTY(NGT,32) = real(WBULKT)
    GTY(NGT,33) = real(WST(1))
    GTY(NGT,34) = real(WST(2))
    GTY(NGT,35) = real(WST(3))

    GTY(NGT,36) = real(TAUE1)
    GTY(NGT,37) = real(TAUE2)
    GTY(NGT,38) = real(TAUEP)
    GTY(NGT,39) = real(BETAA) * 100.0
    GTY(NGT,40) = real(BETA0) * 100.0
    GTY(NGT,41) = real(BETAPA)
    GTY(NGT,42) = real(BETAP0)
    GTY(NGT,43) = real(VLOOP) ! = VPoynt(0)
    GTY(NGT,44) = real(ALI)
    GTY(NGT,45) = real(Q(0))
    GTY(NGT,46) = real(RQ1)

    GTY(NGT,47) = real(ANF0(1))  * 100.0
    GTY(NGT,48) = real(ANFAV(1)) * 100.0

    GTY(NGT,49) = real(VOLAVN)
    ! Initial value becomes too big because almost no neutrals exist.
    if(NGT == 0) then
       GTY(NGT,50) = 0.0
       GTY(NGT,51) = 0.0
       GTY(NGT,52) = 0.0
    else
       GTY(NGT,50) = real(TAUP)
       GTY(NGT,51) = real(TAUPA)
       GTY(NGT,52) = real(Gamma_a) * 1.E-20
    end if
    GTY(NGT,53) = real(TNBcol)
    GTY(NGT,54) = real(TTqt)
    GTY(NGT,55) = real(totmnRV)  ! Total toroidal angular momentum
    GTY(NGT,56) = real(totjxB)  ! Total <j.grad psi> torque

    GTY(NGT,57) = real(PnumN0(0))
    GTY(NGT,58) = real(PnumN0(1))
    GTY(NGT,59) = real(PnumN0(2))
    GTY(NGT,60) = real(PnumN0(3))

    GTY(NGT,61) = real(VLOOPpol)
    GTY(NGT,62) = real(CPsi(0))   ! Poloidal flux related to Vloop
    GTY(NGT,63) = real(CPsi(1))   ! Poloidal flux (inductive)
    GTY(NGT,64) = real(CPsi(2))   ! Poloidal flux (resistive)
    GTY(NGT,65) = real(CPsi(3))   ! Poloidal flux related to Vloop_pol
    GTY(NGT,66) = real(VPoynt(1)) ! Voltage (inductive)
    GTY(NGT,67) = real(VPoynt(2)) ! Voltage (resistive)
    GTY(NGT,68) = real(VPoynt(3)) ! Voltage (poloidal current flux)
    GTY(NGT,69) = real(CEjima)    ! Ejima coefficient
    GTY(NGT,70) = real(PoyntS(1)) ! Poynting flux: psidot
    GTY(NGT,71) = real(PoyntS(2)) ! Poynting flux: psitdot
    GTY(NGT,72) = real(PoyntI(1)) ! Poynting flux: inductive, electric
    GTY(NGT,73) = real(PoyntI(2)) ! Poynting flux: inductive, magnetic
    GTY(NGT,74) = real(PoyntR(1)) ! Poynting flux: resistive, toroidal
    GTY(NGT,75) = real(PoyntR(2)) ! Poynting flux: resistive, poloidal
    GTY(NGT,76) = real(PoyntI(3)) ! Poynting flux: inductive, magnetic, poloidal
    GTY(NGT,77) = real(PoyntI(4)) ! Poynting flux: inductive, magnetic, toroidal

    ! Store data for 3D graphics
    call TXSTGR(NGT,GYL=GYT,NGM=NGTM)

  end subroutine TXSTGT

  !***************************************************************
  !
  !   Store GQY
  !
  !***************************************************************

  subroutine TXSTGQ

    use tx_commons, only : NRMAX, NQMAX, NLCMAX, NLC, NLCR, ALC, BLC, CLC, PLC, X!, t_tx,lqi3,lqi4,lqi2
    integer(4) :: NR, NC, NQ, NC1

    do NQ = 1, NQMAX
       do NC = 1, NLCMAX(NQ)
          NR = 0
          NC1 = NLCR(NC,NQ,0)
          if(NC1 == 0) then
             GQY(NR,NC,NQ) = real(PLC(NR,NC,NQ))
          else
             GQY(NR,NC,NQ) = real(  BLC(NR,NC,NQ) * X(NR  ,NC1) &
                  &               + ALC(NR,NC,NQ) * X(NR+1,NC1) &
                  &               + PLC(NR,NC,NQ))
          end if

          NC1 = NLC(NC,NQ)
          do NR = 1, NRMAX - 1
             if(NC1 == 0) then
                GQY(NR,NC,NQ) = real(PLC(NR,NC,NQ))
             else
                GQY(NR,NC,NQ) = real(  CLC(NR,NC,NQ) * X(NR-1,NC1) &
                     &               + BLC(NR,NC,NQ) * X(NR  ,NC1) &
                     &               + ALC(NR,NC,NQ) * X(NR+1,NC1) &
                     &               + PLC(NR,NC,NQ))
             end if
          end do

          NR = NRMAX
          NC1 = NLCR(NC,NQ,1)
          if(NC1 == 0) then
             GQY(NR,NC,NQ) = real(PLC(NR,NC,NQ))
          else
             GQY(NR,NC,NQ) = real(  CLC(NR,NC,NQ) * X(NR-1,NC1) &
                  &               + BLC(NR,NC,NQ) * X(NR  ,NC1) &
                  &               + PLC(NR,NC,NQ))
          end if
       end do
    end do
!    write(6,*) real(t_tx),gqy(38,3,lqi2),gqy(38,4,lqi2),gqy(38,5,lqi2),gqy(38,6,lqi2)
!lqi3    write(6,*) real(t_tx),gqy(38,5,lqi3),gqy(38,6,lqi3),sum(gqy(38,13:16,lqi3))
!lqi4    write(6,*) real(t_tx),gqy(38,4,lqi4),gqy(38,5,lqi4),sum(gqy(38,6:7,lqi4)),sum(gqy(38,12:15,lqi4))

  end subroutine TXSTGQ

#ifndef nonGSAF
  !***************************************************************
  !
  !   Time evolution of radial profile
  !
  !***************************************************************

  subroutine TXGRFR(NGYRIN,MODE)

    use tx_commons, only : NRMAX, DT, rho, NEMAX, NRA, hv, vlt, NRC, amas, amb, amp
    integer(4), INTENT(IN) :: MODE
    integer(4), INTENT(IN) :: NGYRIN
    integer(4) :: IND, NG, NR, NGYR, NE, IFNT, NRMAXL, NMAX, k
    real(4), DIMENSION(:,:), allocatable :: GYL, GYL2
    character(len=60) :: STR
    character(len=1) :: KSTR
    character(len=3) :: Kend,KLABEL
    real(4), dimension(4) :: GPX, GPY
    real(4) :: GPXL, FACT, GYMAX, GYMIN

    NGYR = NGYRIN

    if (NGR <= -1) then
       write(6,*) 'G', NGYR, ' has no data'
       return
    end if

    if (MODEG == 2) then
       IND = 9
    else
       IND = 0
    end if

    allocate(GYL(0:NRMAX,0:NGR))
    NMAX = size(GYL)

    do

    call PAGES
    call SETCHS(0.3, 0.0)
    call SETLIN(0, 1, 7)
    call INQFNT(IFNT)

    call MOVE(2.0,17.7)
    call TEXT('[G', 2)
    call NUMBI(NGYR, '(I3)', 3)
    call TEXT(']  ', 3)
    call TEXT('FROM', 4)
    call NUMBR(GT(0), '(ES9.2)', 9)
    call TEXT(' TO', 3)
    call NUMBR(GT(NGR), '(ES9.2)', 9)
    call TEXT('  DT =', 6)
    call NUMBD(DT, '(ES9.2)', 9)
    call TEXT('  NGRSTP = ', 11)
    call NUMBI(NGRSTP,'(I4)',4)

    select case(NGYR)
    case(0)
       ! Draw frame of rho
       GPX(1) =  2.5 ; GPY(1) = 10.5
       GPX(2) =  2.5 ; GPY(2) = 16.5
       GPX(3) = 24.5 ; GPY(3) = 16.5
       GPX(4) = 24.5 ; GPY(4) = 10.5
       call LINESC(GPX,GPY,4)
       ! Draw lines
       FACT = (GPX(4) - GPX(1)) / rho(NRMAX)
       GPXL = GPX(1)
       do NE = 1, NEMAX-1
          GPXL = GPXL + real((rho(NE)-rho(NE-1))) * FACT
          if(NE == NRA) call SETLIN(-1,-1,6)
          call LINE1(GPXL,10.5,GPXL,16.5)
          if(NE == NRA) call SETLIN(-1,-1,7)
       end do
       write(KSTR,'(I1)') 0
       write(Kend,'(I3)') NRMAX
       KLABEL='$#r'
       call SETCHS(0.4,0.0)
       call GTEXT(GPX(1),GPY(1)-0.5,KSTR,1,2)
       call GTEXT(GPX(4),GPY(4)-0.5,Kend,3,2)
       call SETFNT(33)
       call SETCHS(0.6,0.0)
       call GTEXT(GPX(1)-1.0,0.5*(GPY(1)+GPY(2)),KLABEL,3,2)
       call SETFNT(IFNT)

       ! Draw frame of vv
       GPX(1) =  2.5 ; GPY(1) =  2.0
       GPX(2) =  2.5 ; GPY(2) =  8.0
       GPX(3) = 24.5 ; GPY(3) =  8.0
       GPX(4) = 24.5 ; GPY(4) =  2.0
       call LINESC(GPX,GPY,4)
       ! Draw lines
       FACT = (GPX(4) - GPX(1)) / vlt(NRMAX)
       GPXL = GPX(1)
       do NE = 1, NEMAX-1
          GPXL = GPXL + real(hv(NE)) * FACT
          if(NE == NRA) call SETLIN(-1,-1,6)
          call LINE1(GPXL,2.0,GPXL,8.0)
          if(NE == NRA) call SETLIN(-1,-1,7)
       end do
       write(KSTR,'(I1)') 0
       write(Kend,'(I3)') NRMAX
       KLABEL='V'
       call SETCHS(0.4,0.0)
       call GTEXT(GPX(1),GPY(1)-0.5,KSTR,1,2)
       call GTEXT(GPX(4),GPY(4)-0.5,Kend,3,2)
       call SETFNT(33)
       call SETCHS(0.6,0.0)
       call GTEXT(GPX(1)-1.0,0.5*(GPY(1)+GPY(2)),KLABEL,1,2)
       call SETFNT(IFNT)
       call SETCHS(0.3,0.0)

    case(1)
       k =  52 ; STR = '@n$-e$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       ! NOTE: GY%v(:,:,1) means GY%v(0:NRMAX,0:NGR,1). NMAX includes the array rank (size) information.
       call TXGRFRX(0, GX, GYL, NRMAX, NGR, STR, MODE, IND)!, GYMAX=0.2)

       k = 138 ; STR = '@Z*n$-i$=+Z*n$-b$=+Z*n$-brp$=-n$-e$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(1, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       k =  91 ; STR = '@E$-r$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(2, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       call TXWPGR

    case(2)
       k =  61 ; STR = '@u$-er$=@'
       call TXGRFRX(0, GX, GY%v(0,0,k), NRMAX, NGR, STR, MODE, IND)!,GYMAX=0.4,GYMIN=0.0)

       k =  62 ; STR = '@u$-e$#q$#$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(1, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       k =  63 ; STR = '@u$-e$#f$#$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(2, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       k =  58 ; STR = '@uhat$-e$#q$#$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(3, GX, GYL, NRMAX, NGR, STR, MODE, IND)

    case(3)
       k =  74 ; STR = '@u$-ir$=@'
       call TXGRFRX(0, GX, GY%v(0,0,k), NRMAX, NGR, STR, MODE, IND)

       k =  75 ; STR = '@u$-i$#q$#$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(1, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       k =  76 ; STR = '@u$-i$#f$#$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(2, GX, GYL, NRMAX, NGR, STR, MODE, IND)!, GYMAX=200.0, GYMIN=-100.0)

       k =  71 ; STR = '@uhat$-i$#q$#$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(3, GX, GYL, NRMAX, NGR, STR, MODE, IND)

    case(4)
       k =  97 ; STR = '@q@'
       call TXGRFRX(0,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k =  93 ; STR = '@E$-tor$=@'
       call TXGRFRX(1,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k =  99 ; STR = '@j$-$#f$#$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(2,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k =  94 ; STR = '@B$-$#q$#$=@'
       call TXGRFRX(3,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

    case(5)
       k = 114 ; STR = '@n$-b$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(0,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k = 120 ; STR = '@u$-b$#f$#$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(1,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k = 122 ; STR = '@Ripple n$-b$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(2,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k = 116 ; STR = '@Bu$-b//$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(3,GX,GYL,NRMAX,NGR,STR,MODE,IND)

!       k = 123 ; STR = '@Ripple n$-b$=/n$-b$=@'
!       call APPROPGY(MODEG, GY%v(0:,:,k), GYL, NMAX, STR, GY%gnrm(k))
!       call TXGRFRX(3,GX,GYL,NRMAX,NGR,STR,MODE,IND)

    case(6)
       k = 109 ; STR = '@j$-e$#f$#$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(0,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k = 110 ; STR = '@j$-i$#f$#$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(1,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k = 111 ; STR = '@j$-z$#f$#$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(2,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k = 139 ; STR = '@j$-r$=@'
       call TXGRFRX(3,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

    case(7)
       k =  56 ; STR = '@Bu$-e//$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(0,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k = 129 ; STR = '@u$-e$#$/136$#$=/B@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(1,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k =  57 ; STR = '@Bq$-e//$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(2,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k = 104 ; STR = '@Bj$-//$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(3,GX,GYL,NRMAX,NGR,STR,MODE,IND)

    case(8)
       k =  69 ; STR = '@Bu$-i//$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(0,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k = 130 ; STR = '@u$-i$#$/136$#$=/B@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(1,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k =  70 ; STR = '@Bq$-i//$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(2,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k =  96 ; STR = '@BE$-//$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(3,GX,GYL,NRMAX,NGR,STR,MODE,IND)

    case(9)
       k = 142 ; STR = '@D$-eff$=@'
       GYMAX = maxval(GY%v(0:NRA,0:NGR,k))
       call TXGRFRX(0,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND,GYMAX=GYMAX)

       k = 141 ; STR = '@sum(D$-s$=)@'
       call TXGRFRX(1,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 145 ; STR = '@$#c$#$-tbe$=@'
       call TXGRFRX(2,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 147 ; STR = '@$#c$#$-NCi$=@'
       call TXGRFRX(3,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

    case(10)
       k =  53 ; STR = '@T$-e$=@'
       call TXGRFRX(0,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k =  66 ; STR = '@T$-i$=@'
       call TXGRFRX(1,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)!,GYMAX=6.0)

       k =  54 ; STR = '@p$-e$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(2,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k =  67 ; STR = '@p$-i$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(3,GX,GYL,NRMAX,NGR,STR,MODE,IND)

    case(11)
       k =  98 ; STR = '@s@'
       call TXGRFRX(0,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 195 ; STR = '@$#a$#@'
       call TXGRFRX(1,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 196 ; STR = '@$#k$#@'
       call TXGRFRX(2,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 194 ; STR = '@F$-CDBM$=@'
       call TXGRFRX(3,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

    case(12)
       k =  92 ; STR = '@BE$-pol$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(0,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k =  95 ; STR = '@I@'
       call TXGRFRX(1,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 328 ; STR = '@U$-g$=V@'
       call TXGRFRX(2,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 329 ; STR = '@qhat$+2$=@'
       call TXGRFRX(3,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

    case(13)
       k = 125 ; STR = '@SLOW N$-0$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(0,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k = 126 ; STR = '@THERMAL N$-0$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(1,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k = 127 ; STR = '@HALO N$-0$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(2,GX,GYL,NRMAX,NGR,STR,MODE,IND)

                 STR = '@TOTAL N$-0$=@'
       do NG = 0, NGR
          do NR = 0, NRMAX
             GYL(NR,NG) = GLOG(real(sum(GY%v(NR,NG,125:127)),8),1.D0,1.D23)
          end do
       end do
       call TXGRFRX(3,GX,GYL,NRMAX,NGR,STR,MODE,IND,ILOGIN=1)

    case(14)
       k = 272 ; STR = '@SCX@'
       do NG = 0, NGR
          do NR = 0, NRMAX
             GYL(NR,NG) = GLOG(real(GY%v(NR,NG,k),8),1.D10,1.D23)
          end do
       end do
       call TXGRFRX(0,GX,GYL,NRMAX,NGR,STR,MODE,IND,GYMIN=18.0,ILOGIN=1)

       k = 271 ; STR = '@SIE@'
       do NG = 0, NGR
          do NR = 0, NRMAX
             GYL(NR,NG) = GLOG(real(GY%v(NR,NG,k),8),1.D10,1.D23)
          end do
       end do
       call TXGRFRX(1,GX,GYL,NRMAX,NGR,STR,MODE,IND,GYMIN=18.0,ILOGIN=1)

       k = 289 ; STR = '@PIE@'
       call TXGRFRX(2,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 290 ; STR = '@PCX@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(3,GX,GYL,NRMAX,NGR,STR,MODE,IND)

    case(15)
       k = 291 ; STR = '@PBr@'
       call TXGRFRX(0,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 283 ; STR = '@POHe@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(1,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k = 286 ; STR = '@PEQei@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(2,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k = 284 ; STR = '@POHi@'
       call TXGRFRX(3,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

    case(16)
       k = 268 ; STR = '@SNB@'
       call TXGRFRX(0,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 274 ; STR = '@PNBe@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(1,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k = 292 ; STR = '@NBI deposition@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(2,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k = 275 ; STR = '@PNBi@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(3,GX,GYL,NRMAX,NGR,STR,MODE,IND)

    case(17)
       k = 293 ; STR = '@SNB perp@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(0,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k = 294 ; STR = '@SNB tang@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(1,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k = 269 ; STR = '@SNB perp ion@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(2,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k = 270 ; STR = '@SNB tang ion@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(3,GX,GYL,NRMAX,NGR,STR,MODE,IND)

    case(18)
!       STR = '@{j$/264$#$/321y$#} [N/m$+2$=]@'
       k = 299 ; STR = '@{j.grad $#y$#} [N/m$+2$=]@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(0,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k = 298 ; STR = '@NBI tor. coll. torque density [N/m$+2$=]@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(1,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k = 280 ; STR = '@PALFe@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(2,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k = 281 ; STR = '@PALFi@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(3,GX,GYL,NRMAX,NGR,STR,MODE,IND)

    case(19)
       do NR = 0, NRMAX
          if(GX(NR) >= 0.95) exit
       end do
       NR = NR - 1

       NRMAXL = NRMAX-NR
       allocate(GYL2(0:NRMAXL,0:NGR))

       k =  52 ; STR = '@n$-e$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       GYL2(0:NRMAXL,0:NGR) = GYL(NR:NRMAX,0:NGR)
       GYMAX = maxval(GYL2(0:NRMAXL,0:NGR))
       call TXGRFRX(0, GX(NR:NRMAX), GYL2, NRMAXL, NGR, STR, MODE, IND, &
            &       GX(NR), GYMAX)

       k =  53 ; STR = '@T$-e$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       GYL2(0:NRMAXL,0:NGR) = GYL(NR:NRMAX,0:NGR)
       GYMAX = maxval(GYL2(0:NRMAXL,0:NGR))
       call TXGRFRX(1, GX(NR:NRMAX), GYL2, NRMAXL, NGR, STR, MODE, IND, &
            &       GX(NR), GYMAX)

       k =  65 ; STR = '@n$-i$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       GYL2(0:NRMAXL,0:NGR) = GYL(NR:NRMAX,0:NGR)
       GYMAX = maxval(GYL2(0:NRMAXL,0:NGR))
       call TXGRFRX(2, GX(NR:NRMAX), GYL2, NRMAXL, NGR, STR, MODE, IND, &
            &       GX(NR), GYMAX)

       k =  66 ; STR = '@T$-i$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       GYL2(0:NRMAXL,0:NGR) = GYL(NR:NRMAX,0:NGR)
       GYMAX = maxval(GYL2(0:NRMAXL,0:NGR))
       call TXGRFRX(3, GX(NR:NRMAX), GYL2, NRMAXL, NGR, STR, MODE, IND, &
            &       GX(NR), GYMAX)

       deallocate(GYL2)

    case(20)
!       k =  55 ; STR = '@u$-e$=$+V$=@'
!       call TXGRFRX(0,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

!       k =  68 ; STR = '@u$-i$=$+V$=@'
!       call TXGRFRX(1,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 338 ; STR = '@$#G$#$-e$=.grad $#r$#@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
!       call TXGRFRX(2,GX,GYL,NRMAX,NGR,STR,MODE,IND,GYMIN=0.0,GYMAX=0.05)
       call TXGRFRX(0,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k = 339 ; STR = '@$#G$#$-i$=.grad $#r$#@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(1,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k = 340 ; STR = '@$#G$#$-z$=.grad $#r$#@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(2,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k = 140 ; STR = '@$#S$#e$-s$=$#G$#$-s$=.grad V w/ beam@'
       call TXGRFRX(3, GX, GY%v(0,0,k), NRMAX, NGR, STR, MODE, IND)

    case(21)
       k = 301 ; STR = '@$#S$#e$-s$=$#G$#$-s$=.grad V@'
       call TXGRFRX(0, GX, GY%v(0,0,k), NRMAX, NGR, STR, MODE, IND)

       k = 140 ; STR = '@$#S$#e$-s$=$#G$#$-s$=.grad V w/ beam@'
       call TXGRFRX(1, GX, GY%v(0,0,k), NRMAX, NGR, STR, MODE, IND)

       k = 341 ; STR = '@$#G$#$-b$=.grad $#r$#@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(2, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       k = 302 ; STR = '@e$-b$=$#G$#$-b$=.grad V@'
       call TXGRFRX(3, GX, GY%v(0,0,k), NRMAX, NGR, STR, MODE, IND)

    case(22)
       k = 132 ; STR = '@BVediag 1@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(0,GX,GYL        ,NRMAX,NGR,STR,MODE,IND)

       k = 134 ; STR = '@BVidiag 1@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(1,GX,GYL        ,NRMAX,NGR,STR,MODE,IND)

       k = 133 ; STR = '@BVediag 2@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(2,GX,GYL        ,NRMAX,NGR,STR,MODE,IND)

       k = 135 ; STR = '@BVidiag 2@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(3,GX,GYL        ,NRMAX,NGR,STR,MODE,IND)

    case(23)
!       k = 161 ; STR = '@VWpch@'
!       call TXGRFRX(0,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 146 ; STR = '@rMue@'
       call TXGRFRX( 0,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 150 ; STR = '@rMui@'
       call TXGRFRX( 1,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 162 ; STR = '@Vmpe@'
       call TXGRFRX( 2,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 163 ; STR = '@Vmpi@'
       call TXGRFRX( 3,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

    case(24)
       k = 192 ; STR = '@$#w$#$-ExB$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(0,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k = 193 ; STR = '@G$-1$=h$+2$=@'
       call TXGRFRX(1,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 148 ; STR = '@$#c$#$-e$=@'
       call TXGRFRX(2,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 152 ; STR = '@$#c$#$-i$=@'
       call TXGRFRX(3,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

    case(25)
       k = 277 ; STR = '@PRFe@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(0,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k = 278 ; STR = '@PRFi@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(1,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k = 296 ; STR = '@Viscous heating@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(2,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k = 300 ; STR = '@Additional Torque@'
       call TXGRFRX(3,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

    case(26)
!!$       k = 138 ; STR = '@sum(Z*n$-j$=)+Z*n$-b$=-n$-e$= (closeup)@'
!!$       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
!!$       GYMIN = minval(GYL(0:NRA-5,0:NGR))
!!$       GYMAX = maxval(GYL(0:NRA-5,0:NGR))
!!$       call TXGRFRX(0, GX, GYL, NRMAX, NGR, STR, MODE, IND, GYMAX=GYMAX, GYMIN=GYMIN)
       k = 202 ; STR = '@Zeff@'
       call TXGRFRX(0,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k =  52 ; STR = '@n$-e$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(1, GX, GYL, NRMAX, NGR, STR, MODE, IND)!, GYMAX=0.2)

       k =  65 ; STR = '@n$-i$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(2, GX, GYL, NRMAX, NGR, STR, MODE, IND)!, GYMAX=0.2)

       k =  78 ; STR = '@n$-z$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(3, GX, GYL, NRMAX, NGR, STR, MODE, IND)!, GYMAX=0.2)

    case(27)
       k =  87 ; STR = '@u$-zr$=@'
       call TXGRFRX(0, GX, GY%v(0,0,k), NRMAX, NGR, STR, MODE, IND)!,GYMAX=0.4,GYMIN=0.0)

       k =  88 ; STR = '@u$-z$#q$#$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(1, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       k =  89 ; STR = '@u$-z$#f$#$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(2, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       k =  84 ; STR = '@uhat$-z$#q$#$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(3, GX, GYL, NRMAX, NGR, STR, MODE, IND)

    case(28)
       k =  82 ; STR = '@Bu$-z//$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(0,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k = 131 ; STR = '@u$-z$#$/136$#$=/B@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(1,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k =  83 ; STR = '@Bq$-z//$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(2,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       call TXWPGR

    case(29)
       k =  79 ; STR = '@T$-z$=@'
       call TXGRFRX(0,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 136 ; STR = '@BVzdiag 1@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(1,GX,GYL        ,NRMAX,NGR,STR,MODE,IND)

       k =  80 ; STR = '@p$-z$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(2,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k = 137 ; STR = '@BVzdiag 2@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(3,GX,GYL        ,NRMAX,NGR,STR,MODE,IND)

    case(30)
       k = 142 ; STR = '@D$-eff,e$=@'
       GYMAX = maxval(GY%v(0:NRA,0:NGR,k))
       GYMIN = minval(GY%v(0:NRA,0:NGR,k))
       call TXGRFRX(0,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND,GYMAX=GYMAX,GYMIN=GYMIN)

       k = 141 ; STR = '@sum(D$-s$=)@'
       call TXGRFRX(1,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 143 ; STR = '@D$-eff,i$=@'
!       GYMAX = maxval(GY%v(0:NRA,0:NGR,k))
!       GYMIN = minval(GY%v(0:NRA,0:NGR,k))
       call TXGRFRX(2,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND,GYMAX=GYMAX,GYMIN=GYMIN)

       k = 144 ; STR = '@D$-eff,z$=@'
!       GYMAX = maxval(GY%v(0:NRC,0:NGR,k))
!       GYMIN = minval(GY%v(0:NRC,0:NGR,k))
       call TXGRFRX(3,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND,GYMAX=GYMAX,GYMIN=GYMIN)

    case(31)
       k =  78 ; STR = '@n$-z$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(0, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       k = 128 ; STR = '@n$-0z$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX(1, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       k = 351 ; STR = '@A$-6$=@'
       do NG = 0, NGR
          do NR = 0, NRMAX
             GYL(NR,NG) = GLOG(real(GY%v(NR,NG,k),8),1.D-30,1.D-10)
          end do
       end do
       call TXGRFRX(2,GX,GYL,NRMAX,NGR,STR,MODE,IND,ILOGIN=1)

       k = 352 ; STR = '@F$-5$=S$-5$=@'
       do NG = 0, NGR
          do NR = 0, NRMAX
             GYL(NR,NG) = GLOG(real(GY%v(NR,NG,k),8),1.D-30,1.D-10)
          end do
       end do
       call TXGRFRX(3,GX,GYL,NRMAX,NGR,STR,MODE,IND,ILOGIN=1)

    ! ------------------------------------------------------------------

    case(-1)
       k =  91 ; STR = '@E$-r$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 0, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       k = 138 ; STR = '@Z*n$-i$=+Z*n$-b$=-n$-e$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 1, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       k = 283 ; STR = '@$#S$#e$-s$=$#G$#$-s$=.grad V@'
       call TXGRFRXS( 2, GX, GY%v(0,0,k), NRMAX, NGR, STR, MODE, IND)

       k = 140 ; STR = '@$#S$#e$-s$=$#G$#$-s$=.grad V w/ beam@'
       call TXGRFRXS( 3, GX, GY%v(0,0,k), NRMAX, NGR, STR, MODE, IND)

       k =  52 ; STR = '@n$-e$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 4,GX, GYL, NRMAX, NGR, STR, MODE, IND)

       k =  61 ; STR = '@u$-er$=@'
       call TXGRFRXS( 5, GX, GY%v(0,0,k), NRMAX, NGR, STR, MODE, IND)

       k =  56 ; STR = '@Bu$-e//$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 6, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       k =  63 ; STR = '@u$-e$#f$#$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 7, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       k =  65 ; STR = '@n$-i$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 8,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k =  74 ; STR = '@u$-ir$=@'
       call TXGRFRXS( 9, GX, GY%v(0,0,k), NRMAX, NGR, STR, MODE, IND)

       k =  69 ; STR = '@Bu$-i//$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS(10, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       k =  76 ; STR = '@u$-i$#f$#$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS(11, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       k =  78 ; STR = '@n$-z$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS(12,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k =  87 ; STR = '@u$-zr$=@'
       call TXGRFRXS(13, GX, GY%v(0,0,k), NRMAX, NGR, STR, MODE, IND)

       k =  82 ; STR = '@Bu$-z//$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS(14, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       k =  89 ; STR = '@u$-z$#f$#$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS(15, GX, GYL, NRMAX, NGR, STR, MODE, IND)

    case(-2)
       k =  53 ; STR = '@T$-e$=@'
       call TXGRFRXS( 0,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k =  57 ; STR = '@Bq$-e//$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 1, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       k =  92 ; STR = '@BE$-pol$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 2,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k =  93 ; STR = '@E$-tor$=@'
       call TXGRFRXS( 3,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k =  66 ; STR = '@T$-i$=@'
       call TXGRFRXS( 4,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k =  70 ; STR = '@Bq$-i//$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 5, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       k =  94 ; STR = '@B$-$#q$#$=@'
       call TXGRFRXS( 6,GX,GY%v(0,0,k),NRMAX  ,NGR,STR,MODE,IND)

       k =  95 ; STR = '@I@'
       call TXGRFRXS( 7,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k =  79 ; STR = '@T$-z$=@'
       call TXGRFRXS( 8,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k =  83 ; STR = '@Bq$-z//$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 9, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       k = 327 ; STR = '@f$-t$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS(10,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k = 102 ; STR = '@j$-NB$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS(11,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k = 114 ; STR = '@n$-b$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS(12,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k = 124 ; STR = '@u$-br$=@'
       call TXGRFRXS(13,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 116 ; STR = '@Bu$-b//$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS(14,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k = 120 ; STR = '@u$-b$#f$#$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS(15,GX,GYL,NRMAX  ,NGR,STR,MODE,IND)     

    case(-3)
       k =  97 ; STR = '@q@'
       call TXGRFRXS( 0,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k =  98 ; STR = '@S@'
       call TXGRFRXS( 1,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 100 ; STR = '@j$-OH$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 2,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k = 101 ; STR = '@j$-BS$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 3,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k = 128 ; STR = '@N$-0z$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 4,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k = 195 ; STR = '@$#a$#@'
       call TXGRFRXS( 5,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 194 ; STR = '@F$-CDBM$=@'
       call TXGRFRXS( 6,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)     

       k = 193 ; STR = '@G$-1$=h$+2$=@'
       call TXGRFRXS( 7,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

                 STR = '@TOTAL N$-0$=@'
       do NG = 0, NGR
          do NR = 0, NRMAX
             GYL(NR,NG) = GLOG(real(sum(GY%v(NR,NG,125:127)),8),1.D0,1.D23)
          end do
       end do
       call TXGRFRXS( 8,GX,GYL,NRMAX,NGR,STR,MODE,IND,ILOGIN=1)

       k = 125 ; STR = '@SLOW N$-0$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 9,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k = 126 ; STR = '@THERMAL N$-0$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS(10,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k = 127 ; STR = '@HALO N$-0$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS(11,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       k =   9 ; STR = '@L$-e$=@'
       call APPROPGY(MODEG, real(amas(1)*amp,4)*GY%v(:,:,k)*1e20, GYL, NMAX, STR, 1e-3)
       call TXGRFRXS(12,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k =  17 ; STR = '@L$-i$=@'
       call APPROPGY(MODEG, real(amas(2)*amp,4)*GY%v(:,:,k)*1e20, GYL, NMAX, STR, 1e-3)
       call TXGRFRXS(13,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k =  25 ; STR = '@L$-z$=@'
       call APPROPGY(MODEG, real(amas(3)*amp,4)*GY%v(:,:,k)*1e20, GYL, NMAX, STR, 1e-3)
       call TXGRFRXS(14,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k =  33 ; STR = '@L$-b$=@'
       call APPROPGY(MODEG, real(amb    *amp,4)*GY%v(:,:,k)*1e20, GYL, NMAX, STR, 1e-3)
       call TXGRFRXS(15,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

    case(-4)
       k =   1 ; STR = '@LQm1@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 0,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k =   2 ; STR = '@LQm2@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 1,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k =   3 ; STR = '@LQm3@'
       call TXGRFRXS( 2,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k =   4 ; STR = '@LQm4@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 3,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k =   5 ; STR = '@LQm5@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 4,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k =   6 ; STR = '@LQe1@'
       call TXGRFRXS( 8,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k =   7 ; STR = '@LQe2@'
       call TXGRFRXS( 9,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k =   8 ; STR = '@LQe3@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS(10,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k =   9 ; STR = '@LQe4@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS(11,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k =  14 ; STR = '@LQi1@'
       call TXGRFRXS(12,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k =  15 ; STR = '@LQi2@'
       call TXGRFRXS(13,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k =  16 ; STR = '@LQi3@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS(14,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k =  17 ; STR = '@LQi4@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS(15,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

    case(-5)
       k =  22 ; STR = '@LQz1@'
       call TXGRFRXS( 0,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k =  23 ; STR = '@LQz2@'
       call TXGRFRXS( 1,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k =  24 ; STR = '@LQz3@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 2,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k =  25 ; STR = '@LQz4@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 3,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k =  10 ; STR = '@LQe5@'
       call TXGRFRXS( 4,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k =  11 ; STR = '@LQe6@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 5,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k =  12 ; STR = '@LQe7@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 6,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k =  13 ; STR = '@LQe8@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 7,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)     

       k =  18 ; STR = '@LQi5@'
       call TXGRFRXS( 8,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k =  19 ; STR = '@LQi6@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 9,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k =  20 ; STR = '@LQi7@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS(10,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k =  21 ; STR = '@LQi8@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS(11,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k =  26 ; STR = '@LQz5@'
       call TXGRFRXS(12,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k =  27 ; STR = '@LQz6@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS(13,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k =  28 ; STR = '@LQz7@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS(14,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k =  29 ; STR = '@LQz8@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS(15,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

    case(-6)
       k =  36 ; STR = '@LQn1@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 0,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)     

       k =  37 ; STR = '@LQn2@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 1,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k =  38 ; STR = '@LQn3@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 2,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k =  39 ; STR = '@LQnz@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 3,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k =  30 ; STR = '@LQb1@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 4,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k =  31 ; STR = '@LQb2@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 5,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k =  32 ; STR = '@LQb3@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 6,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k =  33 ; STR = '@LQb4@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 7,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k =  34 ; STR = '@LQb7@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 8,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k =  35 ; STR = '@LQb8@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 9,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k =  40 ; STR = '@LQr1@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS(10,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

    case(-7)
       k = 146 ; STR = '@rMue@'
       call TXGRFRXS( 0,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 150 ; STR = '@rMui@'
       call TXGRFRXS( 1,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 154 ; STR = '@rMuz@'
       call TXGRFRXS( 2,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 157 ; STR = '@D01@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 3,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 145 ; STR = '@$#c$#$-e$=@'
       call TXGRFRXS( 4,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 149 ; STR = '@$#c$#$-i$=@'
       call TXGRFRXS( 5,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 153 ; STR = '@$#c$#$-z$=@'
       call TXGRFRXS( 6,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 158 ; STR = '@D02@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 7,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 178 ; STR = '@rNu0e@'
       call TXGRFRXS( 8,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 179 ; STR = '@rNu0i@'
       call TXGRFRXS( 9,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 184 ; STR = '@rNu0z@'
       call TXGRFRXS(10,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 159 ; STR = '@D03@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS(11,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 165 ; STR = '@rNuLTe@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS(12,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 182 ; STR = '@rNuLTi@'
       call TXGRFRXS(13,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 185 ; STR = '@rNuLTz@'
       call TXGRFRXS(14,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 164 ; STR = '@rNuL@'
       call TXGRFRXS(15,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

    case(-8)
       k = 175 ; STR = '@rNuTei@'
       call TXGRFRXS( 0,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 176 ; STR = '@rNuTez@'
       call TXGRFRXS( 1,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 177 ; STR = '@rNuTiz@'
       call TXGRFRXS( 2,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 170 ; STR = '@rNuiCX@'
       call TXGRFRXS( 3,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 180 ; STR = '@rNuAsse@'
       do NG = 0, NGR
          do NR = 0, NRMAX
             GYL(NR,NG) = GLOG(real(GY%v(NR,NG,k),8),1.d-3,1.d3)
          end do
       end do
       call TXGRFRXS( 4,GX,GYL,NRMAX,NGR,STR,MODE,IND,ILOGIN=1)

       k = 183 ; STR = '@rNuAssi@'
       do NG = 0, NGR
          do NR = 0, NRMAX
             GYL(NR,NG) = GLOG(real(GY%v(NR,NG,k),8),1.d-3,1.d3)
          end do
       end do
       call TXGRFRXS( 5,GX,GYL,NRMAX,NGR,STR,MODE,IND,ILOGIN=1)

       k = 186 ; STR = '@rNuAssz@'
       do NG = 0, NGR
          do NR = 0, NRMAX
             GYL(NR,NG) = GLOG(real(GY%v(NR,NG,k),8),1.d-3,1.d3)
          end do
       end do
       call TXGRFRXS( 6,GX,GYL,NRMAX,NGR,STR,MODE,IND,ILOGIN=1)

       k = 171 ; STR = '@rNuION@'
       call TXGRFRXS( 7,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)     

       k = 172 ; STR = '@rNuOL@'
       call TXGRFRXS( 8,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 273 ; STR = '@SiLC@'
       call TXGRFRXS( 9,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 191 ; STR = '@WPM@'
       call TXGRFRXS(10,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 173 ; STR = '@rNu0b@'
       call TXGRFRXS(11,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 295 ; STR = '@Vbpara@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS(12,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 174 ; STR = '@rNuB@'
       call TXGRFRXS(13,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 201 ; STR = '@rNubL@'
       call TXGRFRXS(14,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 169 ; STR = '@rNuLB@'
       call TXGRFRXS(15,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

!!$       k = 197 ; STR = '@rNubrp1@'
!!$       call TXGRFRXS( 4,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)
!!$
!!$       k = 198 ; STR = '@DltRP@'
!!$       call TXGRFRXS( 5,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

!!$       k = 199 ; STR = '@Ubrp@'
!!$       call TXGRFRXS(12,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)
!!$
!!$       k = 200 ; STR = '@Dbrp@'
!!$       call TXGRFRXS(13,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

    case(-9)
       k = 220 ; STR = '@l$-11$=$+ee$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 0,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 221 ; STR = '@l$-12$=$+ee$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 1,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 222 ; STR = '@l$-21$=$+ee$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 2,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 223 ; STR = '@l$-22$=$+ee$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 3,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 224 ; STR = '@l$-11$=$+ei$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 4,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 225 ; STR = '@l$-12$=$+ei$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 5,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 226 ; STR = '@l$-21$=$+ei$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 6,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 227 ; STR = '@l$-22$=$+ei$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 7,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 232 ; STR = '@l$-11$=$+ie$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 8,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 233 ; STR = '@l$-12$=$+ie$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 9,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 234 ; STR = '@l$-21$=$+ie$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS(10,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 235 ; STR = '@l$-22$=$+ie$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS(11,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 236 ; STR = '@l$-11$=$+ii$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS(12,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 237 ; STR = '@l$-12$=$+ii$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS(13,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 238 ; STR = '@l$-21$=$+ii$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS(14,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 239 ; STR = '@l$-22$=$+ii$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS(15,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

    case(-10)
       k = 228 ; STR = '@l$-11$=$+ez$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 0,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 229 ; STR = '@l$-12$=$+ez$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 1,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 230 ; STR = '@l$-21$=$+ez$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 2,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 231 ; STR = '@l$-22$=$+ez$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 3,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 244 ; STR = '@l$-11$=$+ze$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 4,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 245 ; STR = '@l$-12$=$+ze$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 5,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 246 ; STR = '@l$-21$=$+ze$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 6,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 247 ; STR = '@l$-22$=$+ze$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 7,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 240 ; STR = '@l$-11$=$+iz$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 8,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 241 ; STR = '@l$-12$=$+iz$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 9,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 242 ; STR = '@l$-21$=$+iz$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS(10,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 243 ; STR = '@l$-22$=$+iz$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS(11,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 248 ; STR = '@l$-11$=$+zi$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS(12,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 249 ; STR = '@l$-12$=$+zi$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS(13,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 250 ; STR = '@l$-21$=$+zi$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS(14,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 251 ; STR = '@l$-22$=$+zi$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS(15,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

    case(-11)
       k = 252 ; STR = '@l$-11$=$+zz$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 0,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 253 ; STR = '@l$-12$=$+zz$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 1,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 254 ; STR = '@l$-21$=$+zz$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 2,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 255 ; STR = '@l$-22$=$+zz$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 3,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 257 ; STR = '@l$-11$=$+ef$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 4,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 258 ; STR = '@l$-21$=$+ef$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 5,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 259 ; STR = '@l$-11$=$+fe$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 6,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 260 ; STR = '@l$-11$=$+if$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 7,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 261 ; STR = '@l$-11$=$+fi$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 8,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 262 ; STR = '@l$-11$=$+zf$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 9,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 263 ; STR = '@l$-11$=$+fz$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS(10,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 264 ; STR = '@l$-11$=$+ff$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS(11,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

    case(-12)
       k = 208 ; STR = '@xmu$-e11$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 0,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 209 ; STR = '@xmu$-e12$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 1,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 210 ; STR = '@xmu$-e21$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 2,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 211 ; STR = '@xmu$-e22$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 3,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 212 ; STR = '@xmu$-i11$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 4,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 213 ; STR = '@xmu$-i12$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 5,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 214 ; STR = '@xmu$-i21$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 6,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 215 ; STR = '@xmu$-i22$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 7,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 216 ; STR = '@xmu$-z11$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 8,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 217 ; STR = '@xmu$-z12$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS( 9,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 218 ; STR = '@xmu$-z21$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS(10,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 219 ; STR = '@xmu$-z22$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS(11,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 256 ; STR = '@xmu$-f$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRXS(12,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

    case(-13)
       k =  56 ; STR = '@Bu$-e//$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX( 0,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k =  57 ; STR = '@Bqhat$-e//$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX( 1,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 342 ; STR = '@Bu$-e//$= Neo. solver@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
!       call TXGRFRX( 2,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)
       call TXGRFRX( 2,GX(1:NRMAX),GYL(1:NRMAX,:),NRMAX-1,NGR,STR,MODE,IND)

       k = 343 ; STR = '@Bqhat$-e//$= Neo. solver@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX( 3,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

    case(-14)
       k =  69 ; STR = '@Bu$-i//$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX( 0,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k =  70 ; STR = '@Bqhat$-i//$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX( 1,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 344 ; STR = '@Bu$-i//$= Neo. solver@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
!       call TXGRFRX( 3,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)
       call TXGRFRX( 2,GX(1:NRMAX),GYL(1:NRMAX,:),NRMAX-1,NGR,STR,MODE,IND)

       k = 345 ; STR = '@Bqhat$-i//$= Neo. solver@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX( 3,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

    case(-15)
       k =  82 ; STR = '@Bu$-z//$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX( 0,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k =  83 ; STR = '@Bqhat$-z//$=@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX( 1,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       k = 346 ; STR = '@Bu$-z//$= Neo. solver@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX( 2,GX(1:NRMAX),GYL(1:NRMAX,:),NRMAX-1,NGR,STR,MODE,IND)

       k = 347 ; STR = '@Bqhat$-z//$= Neo. solver@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX( 3,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

    case(-16)
       k = 281 ; STR = '@ambipolarity@'
       call APPROPGY(MODEG, GY%v(:,:,k), GYL, NMAX, STR, GY%gnrm(k))
       call TXGRFRX( 0,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

    case(-17)
       k = 312 ; STR = '@epst@'
       call TXGRFRXS( 0,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 313 ; STR = '@aat@'
       call TXGRFRXS( 1,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 314 ; STR = '@rrt@'
       call TXGRFRXS( 2,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 315 ; STR = '@ckt@'
       call TXGRFRXS( 3,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 316 ; STR = '@suft@'
       call TXGRFRXS( 4,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 317 ; STR = '@sst@'
       call TXGRFRXS( 5,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 318 ; STR = '@vro@'
       call TXGRFRXS( 6,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 319 ; STR = '@vlt@'
       call TXGRFRXS( 7,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 320 ; STR = '@art@'
       call TXGRFRXS( 8,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 321 ; STR = '@ait@'
       call TXGRFRXS( 9,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 322 ; STR = '@bit@'
       call TXGRFRXS(10,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 323 ; STR = '@bbrt@'
       call TXGRFRXS(11,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 324 ; STR = '@elip@'
       call TXGRFRXS(12,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 325 ; STR = '@trig@'
       call TXGRFRXS(13,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k =  95 ; STR = '@fipol@'
       call TXGRFRXS(14,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

       k = 311 ; STR = '@sdt@'
       call TXGRFRXS(15,GX,GY%v(0,0,k),NRMAX,NGR,STR,MODE,IND)

    case DEFAULT
       write(6,*) 'Unknown NGYR: NGYR = ',NGYR
    end select

    call PAGEE

    select case(NGYR)
    case(-1)
       NGYR = -2  ; cycle
    case(-2)
       NGYR = -3  ; cycle
    case(-3)
       NGYR =  0  ; exit
    case(-4)
       NGYR = -5  ; cycle
    case(-5)
       NGYR = -6  ; cycle
    case(-6)
       NGYR =  0  ; exit
    case(-7)
       NGYR = -8  ; cycle
    case(-8)
       NGYR = -9  ; cycle
    case(-9)
       NGYR = -10 ; cycle
    case(-10)
       NGYR = -11 ; cycle
    case(-11)
       NGYR = -12 ; cycle
    case(-12)
       NGYR =  0  ; exit
    case(-13)
       NGYR = -14 ; cycle
    case(-14)
       NGYR = -15 ; cycle
    case(-15)
       NGYR = -16 ; cycle
    case(-16)
       NGYR =  0  ; exit
    case(-17)
       NGYR =  0  ; exit
    case DEFAULT
       exit
    end select

    end do

    deallocate(GYL)

  end subroutine TXGRFR

  !***************************************************************
  !
  !   Animation of radial profile evolution
  !
  !***************************************************************

  subroutine TXGRFRA(NGYRIN)

    use tx_commons, only : NRMAX
    integer(4), INTENT(IN) :: NGYRIN
    integer(4) :: IND, NG, NGYR, IFNT, I, j, k, NMAX
    real(4), dimension(:),     allocatable :: GMAX, GMIN
    real(4), dimension(:,:,:), allocatable :: GYL
    character(len=60), dimension(:), allocatable :: STR

    NGYR = NGYRIN

    if (NGT <= -1) then
       write(6,*) 'G', NGYR, ' has no data'
       return
    end if

    if (MODEG == 2) then
       IND = 9
    else
       IND = 0
    end if

    allocate(GYL(0:NRMAX,0:NGT,0:15), STR(0:15), GMAX(0:15), GMIN(0:15))
    NMAX = size(GYL,1)*size(GYL,2)

    do

    call PAGES
    call SETCHS(0.3, 0.0)
    call SETLIN(0, 1, 7)
    call INQFNT(IFNT)
    call SETFNT(44)

    select case(NGYR)
    case(-1)
       k =  91 ; j = 0
       STR(j) = '@E$-r$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k),  GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  92 ; j = j + 1
       STR(j) = '@BE$-pol$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  93 ; j = j + 1
       STR(j) = '@E$-tor$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  94 ; j = j + 1
       STR(j) = '@B$-$#q$#$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 138 ; j = j + 1
       STR(j) = '@Z*n$-i$=+Z*n$-b$=-n$-e$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k),  GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  53 ; j = j + 1
       STR(j) = '@T$-e$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  66 ; j = j + 1
       STR(j) = '@T$-i$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  95 ; j = j + 1
       STR(j) = '@I@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  52 ; j = j + 1
       STR(j) = '@n$-e$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k),  GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  61 ; j = j + 1
       STR(j) = '@u$-er$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k),  GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  62 ; j = j + 1
       STR(j) = '@u$-e$#q$#$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k),  GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  63 ; j = j + 1
       STR(j) = '@u$-e$#f$#$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k),  GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  65 ; j = j + 1
       STR(j) = '@n$-i$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  74 ; j = j + 1
       STR(j) = '@u$-ir$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k),  GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  75 ; j = j + 1
       STR(j) = '@u$-i$#q$#$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k),  GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  76 ; j = j + 1
       STR(j) = '@u$-i$#f$#$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k),  GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       do NG = 0, NGT
          call animes
          call gtextx(12.5,17.7,'@T=@',0)
          call gnumbr(13.1,17.7,GTX(NG),3,0)
          do I = 0, j
             call TXGRFRS(I, GX, GYL(0:NRMAX,NG:NG,I), NRMAX, 1, STR(I), 0, IND, 0, 1, &
                  &       'ANIME', GMAX(I), GMIN(I))
          end do
          call animee
       end do

    case(-2)
       j = 0
       k = 114 ; STR(j) = '@n$-b$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 116 ; j = j + 1
       STR(j) = '@Bu$-b//$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 120 ; j = j + 1
       STR(j) = '@u$-b$#f$#$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 122 ; j = j + 1
       STR(j) = '@Ripple n$-b$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

                 j = j + 1
       STR(j) = '@TOTAL N$-0$=@'
       call APPROPGY(MODEG, GYT%v(:,:,125)+GYT%v(:,:,126)+GYT%v(:,:,127), &
            &        GYL(:,:,j), NMAX, STR(j), GYT%gnrm(125), GMAX=GMAX(j), GMIN=GMIN(j))

       k = 125 ; j = j + 1
       STR(j) = '@SLOW N$-0$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 126 ; j = j + 1
       STR(j) = '@THERMAL N$-0$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 127 ; j = j + 1
       STR(j) = '@NBI N$-0$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  97 ; j = j + 1
       STR(j) = '@q@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  99 ; j = j + 1
       STR(j) = '@j$-$#f$#$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 109 ; j = j + 1
       STR(j) = '@j$-e$#f$#$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 110 ; j = j + 1
       STR(j) = '@j$-i$#f$#$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  98 ; j = j + 1
       STR(j) = '@S@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 195 ; j = j + 1
       STR(j) = '@$#a$#@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 194 ; j = j + 1
       STR(j) = '@F$-CDBM$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 193 ; j = j + 1
       STR(j) = '@G$-1$=h$+2$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       do NG = 0, NGT
          call animes
          call gtextx(12.5,17.7,'@T=@',0)
          call gnumbr(13.1,17.7,GTX(NG),3,0)
          do I = 0, j
             call TXGRFRS(I, GX, GYL(0:NRMAX,NG:NG,I), NRMAX, 1, STR(I), 0, IND, 0, 1, &
                  &       'ANIME', GMAX(I), GMIN(I))
          end do
          call animee
       end do

    case(-3)
       k =   1 ; j = 0
       STR(j) = '@LQm1@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =   2 ; j = j + 1
       STR(j) = '@LQm2@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =   3 ; j = j + 1
       STR(j) = '@LQm3@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =   4 ; j = j + 1
       STR(j) = '@LQm4@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =   5 ; j = j + 1
       STR(j) = '@LQm5@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =   6 ; j = 8
       STR(j) = '@LQe1@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =   7 ; j = j + 1
       STR(j) = '@LQe2@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =   8 ; j = j + 1
       STR(j) = '@LQe3@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =   9 ; j = j + 1
       STR(j) = '@LQe4@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  14 ; j = j + 1
       STR(j) = '@LQi1@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  15 ; j = j + 1
       STR(j) = '@LQi2@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  16 ; j = j + 1
       STR(j) = '@LQi3@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  17 ; j = j + 1
       STR(j) = '@LQi4@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       do NG = 0, NGT
          call animes
          call gtextx(12.5,17.7,'@T=@',0)
          call gnumbr(13.1,17.7,GTX(NG),3,0)
          do I = 0, 4
             call TXGRFRS(I, GX, GYL(0:NRMAX,NG:NG,I), NRMAX, 1, STR(I), 0, IND, 0, 1, &
                  &       'ANIME', GMAX(I), GMIN(I))
          end do
          do I = 8, j
             call TXGRFRS(I, GX, GYL(0:NRMAX,NG:NG,I), NRMAX, 1, STR(I), 0, IND, 0, 1, &
                  &       'ANIME', GMAX(I), GMIN(I))
          end do
          call animee
       end do

    case(-4)
       k =  22 ; j = 0
       STR(j) = '@LQz1@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  23 ; j = j + 1
       STR(j) = '@LQz2@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  24 ; j = j + 1
       STR(j) = '@LQz3@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  25 ; j = j + 1
       STR(j) = '@LQz4@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  10 ; j = j + 1
       STR(j) = '@LQe5@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  11 ; j = j + 1
       STR(j) = '@LQe6@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  12 ; j = j + 1
       STR(j) = '@LQe7@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  13 ; j = j + 1
       STR(j) = '@LQe8@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  18 ; j = j + 1
       STR(j) = '@LQi5@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  19 ; j = j + 1
       STR(j) = '@LQi6@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  20 ; j = j + 1
       STR(j) = '@LQi7@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  21 ; j = j + 1
       STR(j) = '@LQi8@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  26 ; j = j + 1
       STR(j) = '@LQz5@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  27 ; j = j + 1
       STR(j) = '@LQz6@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  28 ; j = j + 1
       STR(j) = '@LQz7@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  29 ; j = j + 1
       STR(j) = '@LQz8@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       do NG = 0, NGT
          call animes
          call gtextx(12.5,17.7,'@T=@',0)
          call gnumbr(13.1,17.7,GTX(NG),3,0)
          do I = 0, 15
             call TXGRFRS(I, GX, GYL(:,NG:NG,I), NRMAX, 1, STR(I), 0, IND, 0, 1, &
                  &       'ANIME', GMAX(I), GMIN(I))
          end do
          call animee
       end do

    case(-5)
       k =  36 ; j = 0
       STR(j) = '@LQn1@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  37 ; j = j + 1
       STR(j) = '@LQn2@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  38 ; j = j + 1
       STR(j) = '@LQn3@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  30 ; j = 4
       STR(j) = '@LQb1@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  31 ; j = j + 1
       STR(j) = '@LQb2@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  32 ; j = j + 1
       STR(j) = '@LQb3@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  33 ; j = j + 1
       STR(j) = '@LQb4@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  34 ; j = j + 1
       STR(j) = '@LQb7@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  35 ; j = j + 1
       STR(j) = '@LQb8'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  40 ; j = j + 1
       STR(j) = '@LQr1@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       do NG = 0, NGT
          call animes
          call gtextx(12.5,17.7,'@T=@',0)
          call gnumbr(13.1,17.7,GTX(NG),3,0)
          do I = 0, 2
             call TXGRFRS(I, GX, GYL(:,NG:NG,I), NRMAX, 1, STR(I), 0, IND, 0, 1, &
                  &       'ANIME', GMAX(I), GMIN(I))
          end do
          do I = 4, 10
             call TXGRFRS(I, GX, GYL(:,NG:NG,I), NRMAX, 1, STR(I), 0, IND, 0, 1, &
                  &       'ANIME', GMAX(I), GMIN(I))
          end do
          call animee
       end do

    case(-6)
       k = 146 ; j = 0
       STR(j) = '@rMue@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 150 ; j = j + 1
       STR(j) = '@rMui@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 145 ; j = j + 1
       STR(j) = '@$#c$#$-e$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 149 ; j = j + 1
       STR(j) = '@$#c$#$-i$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 168 ; j = j + 1
       STR(j) = '@rNuL@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 165 ; j = j + 1
       STR(j) = '@rNuLTe@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 182 ; j = j + 1
       STR(j) = '@rNuLTi@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 172 ; j = j + 1
       STR(j) = '@rNuOL@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 170 ; j = j + 1
       STR(j) = '@rNuiCX@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 171 ; j = j + 1
       STR(j) = '@rNuION@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 164 ; j = j + 1
       STR(j) = '@rNu0e@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 181 ; j = j + 1
       STR(j) = '@rNu0i@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 175 ; j = j + 1
       STR(j) = '@rNuTei@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 157 ; j = j + 1
       STR(j) = '@D01@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 158 ; j = j + 1
       STR(j) = '@D02@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 159 ; j = j + 1
       STR(j) = '@D03@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

!!$       k = 277 ; j = j + 1
!!$       STR(j) = '@Vbpara@'
!!$       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
!!$            &        GMAX=GMAX(j), GMIN=GMIN(j))

       do NG = 0, NGT
          call animes
          call gtextx(12.5,17.7,'@T=@',0)
          call gnumbr(13.1,17.7,GTX(NG),3,0)
          do I = 0, j
             call TXGRFRS(I, GX, GYL(:,NG:NG,I), NRMAX, 1, STR(I), 0, IND, 0, 1, &
                  &       'ANIME', GMAX(I), GMIN(I))
          end do
          call animee
       end do

    case(-7)
       k = 208 ; j = 0
       STR(j) = '@xmu$-e11$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 209 ; j = j + 1
       STR(j) = '@xmu$-e12$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 210 ; j = j + 1
       STR(j) = '@xmu$-e21$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 211 ; j = j + 1
       STR(j) = '@xmu$-e22$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 212 ; j = j + 1
       STR(j) = '@xmu$-i11$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 213 ; j = j + 1
       STR(j) = '@xmu$-i12$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 214 ; j = j + 1
       STR(j) = '@xmu$-i21$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 215 ; j = j + 1
       STR(j) = '@xmu$-i22$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 216 ; j = j + 1
       STR(j) = '@WPM@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 217 ; j = j + 1
       STR(j) = '@rNu0b@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 277 ; j = j + 1
       STR(j) = '@Vbpara@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 174 ; j = j + 1
       STR(j) = '@rNuB@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 201 ; j = j + 1
       STR(j) = '@rNubL@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 169 ; j = j + 1
       STR(j) = '@rNuLB@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       do NG = 0, NGT
          call animes
          call gtextx(12.5,17.7,'@T=@',0)
          call gnumbr(13.1,17.7,GTX(NG),3,0)
          do I = 0, j
             call TXGRFRS(I, GX, GYL(:,NG:NG,I), NRMAX, 1, STR(I), 0, IND, 0, 1, &
                  &       'ANIME', GMAX(I), GMIN(I))
          end do
          call animee
       end do

    case(-8)
       k = 220 ; j = 0
       STR(j) = '@l$-11$=$+ee$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 221 ; j = j + 1
       STR(j) = '@l$-12$=$+ee$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 222 ; j = j + 1
       STR(j) = '@l$-21$=$+ee$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 223 ; j = j + 1
       STR(j) = '@l$-22$=$+ee$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 224 ; j = j + 1
       STR(j) = '@l$-11$=$+ei$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 225 ; j = j + 1
       STR(j) = '@l$-12$=$+ei$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 226 ; j = j + 1
       STR(j) = '@l$-21$=$+ei$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 227 ; j = j + 1
       STR(j) = '@l$-22$=$+ei$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 232 ; j = j + 1
       STR(j) = '@l$-11$=$+ie$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 233 ; j = j + 1
       STR(j) = '@l$-12$=$+ie$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 234 ; j = j + 1
       STR(j) = '@l$-21$=$+ie$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 235 ; j = j + 1
       STR(j) = '@l$-22$=$+ie$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 236 ; j = j + 1
       STR(j) = '@l$-11$=$+ii$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 237 ; j = j + 1
       STR(j) = '@l$-12$=$+ii$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 238 ; j = j + 1
       STR(j) = '@l$-21$=$+ii$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 239 ; j = j + 1
       STR(j) = '@l$-22$=$+ii$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       do NG = 0, NGT
          call animes
          call gtextx(12.5,17.7,'@T=@',0)
          call gnumbr(13.1,17.7,GTX(NG),3,0)
          do I = 0, j
             call TXGRFRS(I, GX, GYL(:,NG:NG,I), NRMAX, 1, STR(I), 0, IND, 0, 1, &
                  &       'ANIME', GMAX(I), GMIN(I))
          end do
          call animee
       end do

    case(-9)
       k = 257 ; j = 0
       STR(j) = '@l$-11$=$+ef$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 258 ; j = j + 1
       STR(j) = '@l$-21$=$+ef$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 259 ; j = j + 1
       STR(j) = '@l$-11$=$+fe$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 260 ; j = j + 1
       STR(j) = '@l$-11$=$+if$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 261 ; j = j + 1
       STR(j) = '@l$-11$=$+fi$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 246 ; j = j + 1
       STR(j) = '@l$-11$=$+ff$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       do NG = 0, NGT
          call animes
          call gtextx(12.5,17.7,'@T=@',0)
          call gnumbr(13.1,17.7,GTX(NG),3,0)
          do I = 0, j
             call TXGRFRS(I, GX, GYL(:,NG:NG,I), NRMAX, 1, STR(I), 0, IND, 0, 1, &
                  &       'ANIME', GMAX(I), GMIN(I))
          end do
          call animee
       end do

    case(-10)
       k =  91 ; j = 0
       STR(j) = '@E$-r$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  52 ; j = j + 1
       STR(j) = '@n$-e$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  65 ; j = j + 1
       STR(j) = '@n$-i$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k = 114 ; j = j + 1
       STR(j) = '@n$-b$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  92 ; j = j + 1
       STR(j) = '@BE$-pol$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  53 ; j = j + 1
       STR(j) = '@T$-e$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  66 ; j = j + 1
       STR(j) = '@T$-i$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  94 ; j = j + 1
       STR(j) = '@B$-$#q$#$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  61 ; j = j + 1
       STR(j) = '@u$-er$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  62 ; j = j + 1
       STR(j) = '@u$-e$#q$#$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  63 ; j = j + 1
       STR(j) = '@u$-e$#f$#$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  97 ; j = j + 1
       STR(j) = '@q@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j))!, GMIN(j) = 0.0)

       k = 139 ; j = j + 1
       STR(j) = '@j$-r$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  75 ; j = j + 1
       STR(j) = '@u$-i$#q$#$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  76 ; j = j + 1
       STR(j) = '@u$-i$#f$#$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       k =  99 ; j = j + 1
       STR(j) = '@j$-$#f$#$=@'
       call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,j), NMAX, STR(j), GYT%gnrm(k), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       do NG = 0, NGT
          call animes
          call gtextx(12.5,17.7,'@T=@',0)
          call gnumbr(13.1,17.7,GTX(NG),3,0)
          do I = 0, 15
             call TXGRFRS(I, GX, GYL(:,NG:NG,I), NRMAX, 1, STR(I), 0, IND, 0, 1, &
                  &       'ANIME', GMAX(I), GMIN(I))
          end do
          call animee
       end do

    case DEFAULT
       write(6,*) 'Unknown NGYR: NGYR = ',NGYR
    end select

    call SETFNT(IFNT)
    call PAGEE

    select case(NGYR)
    case(-1)
       NGYR = -2  ; cycle
    case(-2)
       NGYR =  0  ; exit
    case(-3)
       NGYR = -4  ; cycle
    case(-4)
       NGYR = -5  ; cycle
    case(-5)
       NGYR =  0  ; exit
    case(-6)
       NGYR = -7  ; cycle
    case(-7)
       NGYR = -8  ; cycle
    case(-8)
       NGYR = -9  ; cycle
    case(-9)
       NGYR =  0  ; exit
    case(-10)
       NGYR =  0  ; exit
    case DEFAULT
       exit
    end select

    end do

    deallocate(GYL,STR,GMAX,GMIN)

  end subroutine TXGRFRA

  !***************************************************************
  !
  !   Comparison with radial profiles at one slice time
  !
  !***************************************************************

  subroutine TXGRCP(MODE)

    use tx_commons, only : NRMAX, DT, NRA, UsthNCL, BJPARA, BEpara, BJBSvar, ETAvar, Var, &
         &                 BVsdiag, bbt, QsthNCL, MDLNEO, sdt, gflux, aee, achg, drhodr
    use tx_variables, only : TXCALC
    integer(4), intent(in) :: MODE
    character(len=100) :: STR
    integer(4) :: IND, IFNT, NR, inum, ipage, MDLNEOstore, MDLNEOL, i
    real(4), dimension(:,:), allocatable :: GYL, GYL2

    if (NGR <= -1) then
       write(6,*) 'G', NGR, ' has no data'
       return
    end if

    if (MODEG == 2) then
       IND = 9
    else
       IND = 0
    end if

    allocate(GYL(0:NRMAX,1:4),GYL2(0:NRMAX,1:4))

    ! *** 1st page **

    ipage = 1
    call PAGES
    call SETCHS(0.3, 0.0)
    call SETLIN(0, 1, 7)
    call INQFNT(IFNT)

    call MOVE(2.0,17.7)
    call TEXT('[G', 2)
    call NUMBI(ipage, '(I2)', 2)
    call TEXT(']  ', 3)
    call TEXT('FROM', 4)
    call NUMBR(GT(0), '(ES9.2)', 9)
    call TEXT(' TO', 3)
    call NUMBR(GT(NGR), '(ES9.2)', 9)
    call TEXT('  DT =', 6)
    call NUMBD(DT, '(ES9.2)', 9)
    call TEXT('  NGRSTP = ', 11)
    call NUMBI(NGRSTP,'(I4)',4)

!!    call SETFNT(44)
!    call gtextx(5.5,9.4,'@Please compare ETA and BJBS in the limit of Zeff = 1.@',0)
!!    call SETFNT(-1)

    MDLNEOstore = MDLNEO
    MDLNEO = 10 + MDLNEO
    call TXCALC(1)
    MDLNEO = MDLNEOstore

    ! Resistivity

    inum = 4
    do NR = 0, NRMAX
       GYL(NR,1) = GLOG(ETAvar(NR,1),1.D-10,1.D0) ! MI
       GYL(NR,2) = GLOG(ETAvar(NR,2),1.D-10,1.D0) ! NCLASS
       GYL(NR,3) = GLOG(ETAvar(NR,3),1.D-10,1.D0) ! Sauter
       GYL(NR,4) = GLOG(ETAvar(NR,0),1.D-10,1.D0) ! TX
!       GYL(NR,3) = GLOG(ETAS(NR),1.D-10,1.D0)
!       write(6,'(I3,5es15.7)') nr,etavar(nr,0:4)
    end do

    STR = '@LOG: ETA MI, NCLASS, Sauter, TX@'
    call TXGRFRS(0, GX(1:NRA), GYL(1:NRA,:), NRA-1, inum, STR, MODE, IND, 1, 0, 'STATIC')

    ! Bootstrap current

    inum = 3
    GYL(0:NRMAX,1) = real(BJBSvar(0:NRMAX,1)) ! MI
    GYL(0:NRMAX,2) = real(BJBSvar(0:NRMAX,2)) ! NCLASS
    GYL(0:NRMAX,3) = real(BJBSvar(0:NRMAX,3)) ! Sauter

    STR = '@BJ$-BS$= MI, NCLASS, Sauter@'
    call APPROPGY(MODEG, GYL, GYL2, (NRMAX+1)*inum, STR, GY%gnrm(106))
    call TXGRFRS(1, GX(1:NRA), GYL2(1:NRA,1:inum), NRA-1, inum, STR, MODE, IND, 0, 0, 'STATIC')

    ! Total current

    inum = 4
    GYL(0:NRMAX,1) = real(BJPARA(0:NRMAX))                                      ! TX
    GYL(0:NRMAX,2) = real(BEpara(0:NRMAX)/ETAvar(0:NRMAX,1)+BJBSvar(0:NRMAX,1)) ! MI
    GYL(0:NRMAX,3) = real(BEpara(0:NRMAX)/ETAvar(0:NRMAX,2)+BJBSvar(0:NRMAX,2)) ! NCLASS
    GYL(0:NRMAX,4) = real(BEpara(0:NRMAX)/ETAvar(0:NRMAX,3)+BJBSvar(0:NRMAX,3)) ! Sauter

    STR = '@BJ// TX, MI, NCLASS, SAUTER@'
    call APPROPGY(MODEG, GYL, GYL2, (NRMAX+1)*inum, STR, GY%gnrm(104))
    call TXGRFRS(2, GX, GYL2, NRMAX, inum, STR, MODE, IND, 0, 0, 'STATIC')

    ! Comparison of u_r

    inum = 3
    GYL(0:NRMAX,1) = GY%v(0:NRMAX,NGR,61) ! u_er
    GYL(0:NRMAX,2) = GY%v(0:NRMAX,NGR,74) ! u_ir
    GYL(0:NRMAX,3) = GY%v(0:NRMAX,NGR,87) ! u_zr

    STR = '@u$-er$=,u$-ir$=,u$-zr$=@'
    call APPROPGY(MODEG, GYL, GYL2, (NRMAX+1)*inum, STR, GY%gnrm(61))
    call TXGRFRS(3, GX, GYL2(:,1:inum), NRMAX, inum, STR, MODE, IND, 0, 0, 'STATIC' &
         &     , GMAX=min(maxval(GYL(:,1:3)),2.0), GMIN=max(minval(GYL(:,1:3)),-2.0))

    call PAGEE

    ! *** 2nd page **

    ipage = 2
    call PAGES
    call SETCHS(0.3, 0.0)
    call SETLIN(0, 1, 7)
    call INQFNT(IFNT)

    call MOVE(2.0,17.7)
    call TEXT('[G', 2)
    call NUMBI(ipage, '(I2)', 2)
    call TEXT(']  ', 3)
    call TEXT('FROM', 4)
    call NUMBR(GT(0), '(ES9.2)', 9)
    call TEXT(' TO', 3)
    call NUMBR(GT(NGR), '(ES9.2)', 9)
    call TEXT('  DT =', 6)
    call NUMBD(DT, '(ES9.2)', 9)
    call TEXT('  NGRSTP = ', 11)
    call NUMBI(NGRSTP,'(I4)',4)

    ! Particle flux  <Gamma_s . grad psi>

    MDLNEOL = mod(MDLNEO,10)
    if( MDLNEOL == 3 ) MDLNEOL = 2

    inum = 3
    i = 1 ! electrons
    GYL(0:NRMAX,1) = real(Var(0:NRMAX,i)%n*Var(0:NRMAX,i)%UrV*sdt(0:NRMAX))*1.e20 ! TX
    GYL(0:NRMAX,2) = real(gflux(0:NRMAX,i,MDLNEOL))                               ! MI
    GYL(0:NRMAX,3) = real(gflux(0:NRMAX,i,3))                                     ! NCLASS

    STR = '@$#G$#$-e$=$+psi$= [10$+20$=T/m$+2$=s$+1$=] TX, MI, NCLASS@'
    call APPROPGY(MODEG, GYL, GYL2, (NRMAX+1)*inum, STR, GY%gnrm(338))
    call TXGRFRS(0, GX, GYL2, NRMAX, inum, STR, MODE, IND, 0, 0, 'STATIC')

    i = 2 ! ions
    GYL(0:NRMAX,1) = real(Var(0:NRMAX,i)%n*Var(0:NRMAX,i)%UrV*sdt(0:NRMAX))*1.e20 ! TX
    GYL(0:NRMAX,2) = real(gflux(0:NRMAX,i,MDLNEOL))                               ! MI
    GYL(0:NRMAX,3) = real(gflux(0:NRMAX,i,3))                                     ! NCLASS

    STR = '@$#G$#$-i$=$+psi$= [10$+20$=T/m$+2$=s$+1$=] TX, MI, NCLASS@'
    call APPROPGY(MODEG, GYL, GYL2, (NRMAX+1)*inum, STR, GY%gnrm(339))
    call TXGRFRS(1, GX, GYL2, NRMAX, inum, STR, MODE, IND, 0, 0, 'STATIC')

    i = 3 ! impurities
    GYL(0:NRMAX,1) = real(Var(0:NRMAX,i)%n*Var(0:NRMAX,i)%UrV*sdt(0:NRMAX))*1.e20 ! TX
    GYL(0:NRMAX,2) = real(gflux(0:NRMAX,i,MDLNEOL))                               ! MI
    GYL(0:NRMAX,3) = real(gflux(0:NRMAX,i,3))                                     ! NCLASS

    STR = '@$#G$#$-z$=$+psi$= [10$+20$=T/m$+2$=s$+1$=] TX, MI, NCLASS@'
    call APPROPGY(MODEG, GYL, GYL2, (NRMAX+1)*inum, STR, GY%gnrm(340))
    call TXGRFRS(2, GX, GYL2, NRMAX, inum, STR, MODE, IND, 0, 0, 'STATIC')

    ! Comparison of e_s<Gamma_s.grad rho>

    inum = 3
    do i = 1, inum
       GYL(0,      i) = 0.0
       GYL(1:NRMAX,i) = real(achg(i) * aee / drhodr(1:NRMAX)) * GY%v(1:NRMAX,NGR,338+i-1) ! e_s<Gamma_s.grad r>
    end do

    STR = '@e$-e$=$#G$#$-e$=.grad r,e$-i$=$#G$#$-i$=.grad r,e$-z$=$#G$#$-z$=.grad r@'
    call TXGRFRS(3, GX, GYL(:,1:inum), NRMAX, inum, STR, MODE, IND, 0, 0, 'STATIC')

    call PAGEE

    ! *** 3rd page **

    ipage = 3
    call PAGES
    call SETCHS(0.3, 0.0)
    call SETLIN(0, 1, 7)
    call INQFNT(IFNT)

    call MOVE(2.0,17.7)
    call TEXT('[G', 2)
    call NUMBI(ipage, '(I2)', 2)
    call TEXT(']  ', 3)
    call TEXT('FROM', 4)
    call NUMBR(GT(0), '(ES9.2)', 9)
    call TEXT(' TO', 3)
    call NUMBR(GT(NGR), '(ES9.2)', 9)
    call TEXT('  DT =', 6)
    call NUMBD(DT, '(ES9.2)', 9)
    call TEXT('  NGRSTP = ', 11)
    call NUMBI(NGRSTP,'(I4)',4)

    ! Parallel particle and heat flows

    inum = 2
    GYL(0:NRMAX,1) = GY%v(0:NRMAX,NGR, 56) ! TX
    GYL(0:NRMAX,2) = GY%v(0:NRMAX,NGR,342) ! Neo. solver

    STR = '@Bu$-e//$= TX, Neo. solver@'
    call APPROPGY(MODEG, GYL, GYL2, (NRMAX+1)*inum, STR, GY%gnrm(56))
    call TXGRFRS(0, GX(1:NRMAX), GYL2(1:NRMAX,1:inum), NRMAX-1, inum, STR, MODE, IND, 0, 0, 'STATIC')
!    call TXGRFRS(0, GX, GYL2, NRMAX, inum, STR, MODE, IND, 0, 0, 'STATIC')

    GYL(0:NRMAX,1) = GY%v(0:NRMAX,NGR, 69) ! TX
    GYL(0:NRMAX,2) = GY%v(0:NRMAX,NGR,344) ! Neo. solver

    STR = '@Bu$-i//$= TX, Neo. solver@'
    call APPROPGY(MODEG, GYL, GYL2, (NRMAX+1)*inum, STR, GY%gnrm(69))
    call TXGRFRS(1, GX(1:NRMAX), GYL2(1:NRMAX,1:inum), NRMAX-1, inum, STR, MODE, IND, 0, 0, 'STATIC')
!    call TXGRFRS(1, GX, GYL2, NRMAX, inum, STR, MODE, IND, 0, 0, 'STATIC')

    GYL(0:NRMAX,1) = GY%v(0:NRMAX,NGR, 57) ! TX
    GYL(0:NRMAX,2) = GY%v(0:NRMAX,NGR,343) ! Neo. solver

    STR = '@Bqhat$-e//$= TX, Neo. solver@'
    call APPROPGY(MODEG, GYL, GYL2, (NRMAX+1)*inum, STR, GY%gnrm(57))
    call TXGRFRS(2, GX, GYL2, NRMAX, inum, STR, MODE, IND, 0, 0, 'STATIC')

    GYL(0:NRMAX,1) = GY%v(0:NRMAX,NGR, 70) ! TX
    GYL(0:NRMAX,2) = GY%v(0:NRMAX,NGR,345) ! Neo. solver

    STR = '@Bqhat$-i//$= TX, Neo. solver@'
    call APPROPGY(MODEG, GYL, GYL2, (NRMAX+1)*inum, STR, GY%gnrm(70))
    call TXGRFRS(3, GX, GYL2, NRMAX, inum, STR, MODE, IND, 0, 0, 'STATIC')

    call PAGEE

    ! *** 4th page **

    ipage = 4
    call PAGES
    call SETCHS(0.3, 0.0)
    call SETLIN(0, 1, 7)
    call INQFNT(IFNT)

    call MOVE(2.0,17.7)
    call TEXT('[G', 2)
    call NUMBI(ipage, '(I2)', 2)
    call TEXT(']  ', 3)
    call TEXT('FROM', 4)
    call NUMBR(GT(0), '(ES9.2)', 9)
    call TEXT(' TO', 3)
    call NUMBR(GT(NGR), '(ES9.2)', 9)
    call TEXT('  DT =', 6)
    call NUMBD(DT, '(ES9.2)', 9)
    call TEXT('  NGRSTP = ', 11)
    call NUMBI(NGRSTP,'(I4)',4)

    ! Parallel particle flows

    inum = 2
    GYL(0:NRMAX,1) = GY%v(0:NRMAX,NGR, 82) ! TX
    GYL(0:NRMAX,2) = GY%v(0:NRMAX,NGR,346) ! Neo. solver

    STR = '@Bu$-z//$= TX, Neo. solver@'
    call APPROPGY(MODEG, GYL, GYL2, (NRMAX+1)*inum, STR, GY%gnrm(82))
    call TXGRFRS(0, GX(1:NRMAX), GYL2(1:NRMAX,1:inum), NRMAX-1, inum, STR, MODE, IND, 0, 0, 'STATIC')
!    call TXGRFRS(0, GX, GYL2, NRMAX, inum, STR, MODE, IND, 0, 0, 'STATIC')

    ! Parallel heat flows

    inum = 2
    GYL(0:NRMAX,1) = GY%v(0:NRMAX,NGR, 83) ! TX
    GYL(0:NRMAX,2) = GY%v(0:NRMAX,NGR,347) ! Neo. solver

    STR = '@Bqhat$-z//$= TX, Neo. solver@'
    call APPROPGY(MODEG, GYL, GYL2, (NRMAX+1)*inum, STR, GY%gnrm(83))
    call TXGRFRS(2, GX, GYL2, NRMAX, inum, STR, MODE, IND, 0, 0, 'STATIC')

    ! Poloidal rotations

    inum = 4
    i = 1
    GYL(0:NRMAX,1) = real(Var(0:NRMAX,i)%Uthhat) ! TX
    GYL(0:NRMAX,3) = real(UsthNCL(0:NRMAX,i,1))  ! Neo. solver, from p' T'
    GYL(0:NRMAX,4) = real(UsthNCL(0:NRMAX,i,2))  ! Neo. solver, from <E.B>
    GYL(0:NRMAX,2) = GYL(0:NRMAX,3) + GYL(0:NRMAX,4) ! Neo. solver, total

    STR = '@uhat$-e$#q$#$= TX, Neo solver@'
    call APPROPGY(MODEG, GYL, GYL2, (NRMAX+1)*inum, STR, GY%gnrm(58))
    call TXGRFRS(1, GX(1:NRA), GYL2(1:NRA,:), NRA-1, 4, STR, MODE, IND, 0, 0, 'STATIC')

    i = 2
    GYL(0:NRMAX,1) = real(Var(0:NRMAX,i)%Uthhat)
    GYL(0:NRMAX,3) = real(UsthNCL(0:NRMAX,i,1))
    GYL(0:NRMAX,4) = real(UsthNCL(0:NRMAX,i,2))
    GYL(0:NRMAX,2) = GYL(0:NRMAX,3) + GYL(0:NRMAX,4)

    STR = '@uhat$-i$#q$#$= TX, Neo solver@'
    call APPROPGY(MODEG, GYL, GYL2, (NRMAX+1)*inum, STR, GY%gnrm(71))
    call TXGRFRS(3, GX(1:NRA), GYL2(1:NRA,:), NRA-1, 4, STR, MODE, IND, 0, 0, 'STATIC')

    call PAGEE

    ! *** 5th page **

    ipage = 5
    call PAGES
    call SETCHS(0.3, 0.0)
    call SETLIN(0, 1, 7)
    call INQFNT(IFNT)

    call MOVE(2.0,17.7)
    call TEXT('[G', 2)
    call NUMBI(ipage, '(I2)', 2)
    call TEXT(']  ', 3)
    call TEXT('FROM', 4)
    call NUMBR(GT(0), '(ES9.2)', 9)
    call TEXT(' TO', 3)
    call NUMBR(GT(NGR), '(ES9.2)', 9)
    call TEXT('  DT =', 6)
    call NUMBD(DT, '(ES9.2)', 9)
    call TEXT('  NGRSTP = ', 11)
    call NUMBI(NGRSTP,'(I4)',4)

    ! Poloidal rotations

    inum = 4
    i = 3
    GYL(0:NRMAX,1) = real(Var(0:NRMAX,i)%Uthhat) ! TX
    GYL(0:NRMAX,3) = real(UsthNCL(0:NRMAX,i,1))  ! Neo. solver, from p' T'
    GYL(0:NRMAX,4) = real(UsthNCL(0:NRMAX,i,2))  ! Neo. solver, from <E.B>
    GYL(0:NRMAX,2) = GYL(0:NRMAX,3) + GYL(0:NRMAX,4) ! Neo. solver, total

    STR = '@uhat$-z$#q$#$= TX, Neo solver@'
    call APPROPGY(MODEG, GYL, GYL2, (NRMAX+1)*inum, STR, GY%gnrm(58))
    call TXGRFRS(0, GX(1:NRA), GYL2(1:NRA,:), NRA-1, 4, STR, MODE, IND, 0, 0, 'STATIC')

    ! Poloidal heat rotations

    i = 1
    GYL(0:NRMAX,1) = real((Var(0:NRMAX,i)%Bqpar - BVsdiag(0:NRMAX,i,2))/bbt(0:NRMAX))
    GYL(0:NRMAX,3) = real(QsthNCL(0:NRMAX,i,1))
    GYL(0:NRMAX,4) = real(QsthNCL(0:NRMAX,i,2))
    GYL(0:NRMAX,2) = GYL(0:NRMAX,3) + GYL(0:NRMAX,4)

    STR = '@qhat$-e$#q$#$= TX, Neo solver@'
    call APPROPGY(MODEG, GYL, GYL2, (NRMAX+1)*inum, STR, GY%gnrm(71))
    call TXGRFRS(1, GX(1:NRA), GYL2(1:NRA,:), NRA-1, 4, STR, MODE, IND, 0, 0, 'STATIC')

    i = 2
    GYL(0:NRMAX,1) = real((Var(0:NRMAX,i)%Bqpar - BVsdiag(0:NRMAX,i,2))/bbt(0:NRMAX))
    GYL(0:NRMAX,3) = real(QsthNCL(0:NRMAX,i,1))
    GYL(0:NRMAX,4) = real(QsthNCL(0:NRMAX,i,2))
    GYL(0:NRMAX,2) = GYL(0:NRMAX,3) + GYL(0:NRMAX,4)

    STR = '@qhat$-i$#q$#$= TX, Neo solver@'
    call APPROPGY(MODEG, GYL, GYL2, (NRMAX+1)*inum, STR, GY%gnrm(71))
    call TXGRFRS(2, GX(1:NRA), GYL2(1:NRA,:), NRA-1, 4, STR, MODE, IND, 0, 0, 'STATIC')

    i = 3
    GYL(0:NRMAX,1) = real((Var(0:NRMAX,i)%Bqpar - BVsdiag(0:NRMAX,i,2))/bbt(0:NRMAX))
    GYL(0:NRMAX,3) = real(QsthNCL(0:NRMAX,i,1))
    GYL(0:NRMAX,4) = real(QsthNCL(0:NRMAX,i,2))
    GYL(0:NRMAX,2) = GYL(0:NRMAX,3) + GYL(0:NRMAX,4)

    STR = '@qhat$-z$#q$#$= TX, Neo solver@'
    call APPROPGY(MODEG, GYL, GYL2, (NRMAX+1)*inum, STR, GY%gnrm(71))
    call TXGRFRS(3, GX(1:NRA), GYL2(1:NRA,:), NRA-1, 4, STR, MODE, IND, 0, 0, 'STATIC')

    call PAGEE

    ! *** 6th page **

    ipage = 6
    call PAGES
    call SETCHS(0.3, 0.0)
    call SETLIN(0, 1, 7)
    call INQFNT(IFNT)

    call MOVE(2.0,17.7)
    call TEXT('[G', 2)
    call NUMBI(ipage, '(I2)', 2)
    call TEXT(']  ', 3)
    call TEXT('FROM', 4)
    call NUMBR(GT(0), '(ES9.2)', 9)
    call TEXT(' TO', 3)
    call NUMBR(GT(NGR), '(ES9.2)', 9)
    call TEXT('  DT =', 6)
    call NUMBD(DT, '(ES9.2)', 9)
    call TEXT('  NGRSTP = ', 11)
    call NUMBI(NGRSTP,'(I4)',4)

    ! Ambipolarity & Self-adjointness of the linearlized collision operator

    inum = 3
    i = 1
    do NR = 0, NRMAX
       GYL(NR,1) = real(sum(achg(:)*Var(NR,:)%n*Var(NR,:)%UrV)*sdt(NR)*aee*1.d20) ! TX, same as 299
       GYL(NR,2) = real(sum(achg(:)*gflux(NR,:,MDLNEOL)*aee))                     ! MI
       GYL(NR,3) = real(sum(achg(:)*gflux(NR,:,3)*aee))                           ! NCLASS
    end do

    STR = '@sum e$-s$=$#G$#$-s$=$+psi$= [Pa] TX@'
    call TXGRFRS(0, GX(1:NRA), GYL(1:NRA,1:), NRA-1, 1, STR, MODE, IND, 0, 0, 'STATIC')

    STR = '@sum e$-s$=$#G$#$-s$=$+psi$= [Pa] MI@'
    call TXGRFRS(1, GX(1:NRA), GYL(1:NRA,2:), NRA-1, 1, STR, MODE, IND, 0, 0, 'STATIC')

    STR = '@sum e$-s$=$#G$#$-s$=$+psi$= [Pa] NCLASS@'
    call TXGRFRS(2, GX(1:NRA), GYL(1:NRA,3:), NRA-1, 1, STR, MODE, IND, 0, 0, 'STATIC')

    call PAGEE

    deallocate(GYL,GYL2)

  end subroutine TXGRCP

  !***************************************************************
  !
  !   Animation of comparison with radial profiles
  !
  !***************************************************************

  subroutine TXGRCPA(MODE)

    use tx_commons, only : NRMAX, DT, achg, aee, drhodr
    integer(4), intent(in) :: MODE
    integer(4) :: IND, IFNT, NG, I, N, j, NMAX, k, NR!, inum
    real(4), dimension(:),     allocatable :: GMAX, GMIN
    real(4), dimension(:,:,:), allocatable :: GYL
    character(len=60), dimension(:), allocatable :: STR

    if (NGR <= -1) then
       write(6,*) 'G', NGR, ' has no data'
       return
    end if

    if (MODEG == 2) then
       IND = 9
    else
       IND = 0
    end if

    N = NRMAX
    allocate(GYL(0:N,0:NGT,0:9), STR(0:3), GMAX(0:9), GMIN(0:9))
    NMAX = size(GYL,1)*size(GYL,2)
    STR = ' '

    ! *** Correspond to 3rd page of TXGRCP ***

    call PAGES
    call SETCHS(0.3, 0.0)
    call SETLIN(0, 1, 7)
    call INQFNT(IFNT)

    call MOVE(2.0,17.7)
!    call TEXT('[G', 2)
!    call NUMBI(ipage, '(I2)', 2)
!    call TEXT(']  ', 3)
!    call TEXT('FROM', 4)
!    call NUMBR(GT(0), '(ES9.2)', 9)
!    call TEXT(' TO', 3)
!    call NUMBR(GT(NGR), '(ES9.2)', 9)
    call TEXT('  DT =', 6)
    call NUMBD(DT, '(ES9.2)', 9)
!    call TEXT('  NGRSTP = ', 11)
!    call NUMBI(NGRSTP,'(I4)',4)

    ! Total current

    j = 0
    k =  56
    call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,0), NMAX, STR(j), GYT%gnrm(k), &
         &        GMAX=GMAX(0), GMIN=GMIN(0)) ! TX
    k = 342 ; STR(j) = '@Bu$-e//$= TX, Neo. solver@'
    call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,1), NMAX, STR(j), GYT%gnrm(k), &
         &        GMAX=GMAX(1), GMIN=GMIN(1)) ! Neo. solver
    GMAX(j) = max(GMAX(0),GMAX(1))
    GMIN(j) = min(GMIN(0),GMIN(1))

    j = j + 1
    k =  69
    call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,2), NMAX, STR(j), GYT%gnrm(k), &
         &        GMAX=GMAX(2), GMIN=GMIN(2)) ! TX
    k = 344 ; STR(1) = '@Bu$-i//$= TX, Neo. solver@'
    call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,3), NMAX, STR(j), GYT%gnrm(k), &
         &        GMAX=GMAX(3), GMIN=GMIN(3)) ! Neo. solver
!    GMAX(j) = max(GMAX(2),GMAX(3))
!    GMIN(j) = min(GMIN(2),GMIN(3))
    GMAX(j) = min(GMAX(2),GMAX(3))
    GMIN(j) = max(GMIN(2),GMIN(3))

    j = j + 1
    k =  57
    call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,4), NMAX, STR(j), GYT%gnrm(k), &
         &        GMAX=GMAX(4), GMIN=GMIN(4)) ! TX
    k = 343 ; STR(2) = '@Bqhat$-e//$= TX, Neo. solver@'
    call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,5), NMAX, STR(j), GYT%gnrm(k), &
         &        GMAX=GMAX(5), GMIN=GMIN(5)) ! Neo. solver
    GMAX(j) = max(GMAX(4),GMAX(5))
    GMIN(j) = min(GMIN(4),GMIN(5))

    j = j + 1
    k = 70
    call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,6), NMAX, STR(j), GYT%gnrm(k), &
         &        GMAX=GMAX(6), GMIN=GMIN(6)) ! TX
    k = 345 ; STR(j) = '@Bqhat$-i//$= TX, Neo. solver@'
    call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,7), NMAX, STR(j), GYT%gnrm(k), &
         &        GMAX=GMAX(7), GMIN=GMIN(7)) ! Neo. solver
    GMAX(j) = max(GMAX(6),GMAX(7))
    GMIN(j) = min(GMIN(6),GMIN(7))

    call SETFNT(44)
    do NG = 0, NGT
       call animes
       call SETLIN(-1, -1, 7)
       call gtextx(12.5,17.7,'@T=@',0)
       call gnumbr(13.1,17.7,GTX(NG),3,0)
       do I = 0, 1
          call TXGRFRS(I, GX(1:N), GYL(1:N,NG,2*I:2*I+1), N-1, 2, STR(I), MODE, IND, 0, 0, &
               &       'ANIME', GMAX(I), GMIN(I))
       end do
       do I = 2, 3
          call TXGRFRS(I, GX(0:N), GYL(0:N,NG,2*I:2*I+1), N  , 2, STR(I), MODE, IND, 0, 0, &
               &       'ANIME', GMAX(I), GMIN(I))
       end do
       call animee
    end do

!!$    inum = 2
!!$    GYL(0:NRMAX,1) = GY%v(0:NRMAX,NGR, 69) ! TX
!!$    GYL(0:NRMAX,2) = GY%v(0:NRMAX,NGR,344) ! Neo. solver
!!$
!!$    NMAX = (NRMAX+1)*inum
!!$
!!$    STR = '@Bu$-i//$= TX, Neo. solver@'
!!$    call APPROPGY(MODEG, GYL, GYL2, NMAX, STR, GY%gnrm(69))
!!$    call TXGRFRS(1, GX(1:NRMAX), GYL2(1:NRMAX,1:inum), NRMAX-1, inum, STR, MODE, IND, 0, 0, 'STATIC')
!!$
!!$    GYL(0:NRMAX,1) = GY%v(0:NRMAX,NGR, 57) ! TX
!!$    GYL(0:NRMAX,2) = GY%v(0:NRMAX,NGR,343) ! Neo. solver
!!$
!!$    STR = '@Bqhat$-e//$= TX, Neo. solver@'
!!$    call APPROPGY(MODEG, GYL, GYL2, NMAX, STR, GY%gnrm(57))
!!$    call TXGRFRS(2, GX, GYL2, NRMAX, inum, STR, MODE, IND, 0, 0, 'STATIC')
!!$
!!$    GYL(0:NRMAX,1) = GY%v(0:NRMAX,NGR,70) ! TX
!!$    GYL(0:NRMAX,2) = GY%v(0:NRMAX,NGR,345) ! Neo. solver
!!$
!!$    STR = '@Bqhat$-i//$= TX, Neo. solver@'
!!$    call APPROPGY(MODEG, GYL, GYL2, NMAX, STR, GY%gnrm(70))
!!$    call TXGRFRS(3, GX, GYL2, NRMAX, inum, STR, MODE, IND, 0, 0, 'STATIC')

    call PAGEE

    ! *** Correspond to 4th page of TXGRCP ***

    call PAGES
    call SETCHS(0.3, 0.0)
    call SETLIN(0, 1, 7)
    call INQFNT(IFNT)

    call MOVE(2.0,17.7)
!    call TEXT('[G', 2)
!    call NUMBI(ipage, '(I2)', 2)
!    call TEXT(']  ', 3)
!    call TEXT('FROM', 4)
!    call NUMBR(GT(0), '(ES9.2)', 9)
!    call TEXT(' TO', 3)
!    call NUMBR(GT(NGR), '(ES9.2)', 9)
    call TEXT('  DT =', 6)
    call NUMBD(DT, '(ES9.2)', 9)
!    call TEXT('  NGRSTP = ', 11)
!    call NUMBI(NGRSTP,'(I4)',4)

    ! Total current

    j = 0
    k =  82
    call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,0), NMAX, STR(j), GYT%gnrm(k), &
         &        GMAX=GMAX(0), GMIN=GMIN(0)) ! TX
    k = 346 ; STR(j) = '@Bu$-z//$= TX, Neo. solver@'
    call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,1), NMAX, STR(j), GYT%gnrm(k), &
         &        GMAX=GMAX(1), GMIN=GMIN(1)) ! Neo. solver
    GMAX(j) = max(GMAX(0),GMAX(1))
    GMIN(j) = min(GMIN(0),GMIN(1))

!!$    j = j + 1
!!$    k = 338
!!$    call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,2), NMAX, STR(j), GYT%gnrm(k), &
!!$         &        GMAX=GMAX(2), GMIN=GMIN(2))
!!$    k = 339
!!$    call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,3), NMAX, STR(j), GYT%gnrm(k), &
!!$         &        GMAX=GMAX(3), GMIN=GMIN(3))
!!$    k = 340 ; STR(j) = '@e$-e$=$#G$#$-e$=.grad r,e$-i$=$#G$#$-i$=.grad r,e$-z$=$#G$#$-z$=.grad r@'
!!$    call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,4), NMAX, STR(j), GYT%gnrm(k), &
!!$         &        GMAX=GMAX(4), GMIN=GMIN(4))
!!$    do i = 1, 3
!!$       do NR = 1, NRMAX
!!$          GYL(NR,:,i+1) = GYL(NR,:,i+1) * real(achg(i) * aee / drhodr(NR))
!!$       end do
!!$       GMAX(i+1) = GMAX(i+1) * real(achg(i) * aee / maxval(drhodr(:)))
!!$       GMIN(i+1) = GMIN(i+1) * real(achg(i) * aee / minval(drhodr(:)))
!!$    end do
!!$    GMAX(j) = min(GMAX(2),GMAX(3),GMAX(4))
!!$    GMIN(j) = max(GMIN(2),GMIN(3),GMIN(4))
    j = j + 1
    STR(j) = '@e$-e$=$#G$#$-e$=.grad r,e$-i$=$#G$#$-i$=.grad r,e$-z$=$#G$#$-z$=.grad r@'
    GYL(0,:,2:4) = 0.0
    do NR = 1, NRMAX
       k = 338
       GYL(NR,:,2) = GYT%v(NR,:,k) * real(achg(1) * aee / drhodr(NR))
       k = 339
       GYL(NR,:,3) = GYT%v(NR,:,k) * real(achg(2) * aee / drhodr(NR))
       k = 340
       GYL(NR,:,4) = GYT%v(NR,:,k) * real(achg(3) * aee / drhodr(NR))
    end do
    do i = 2, 4
       GMAX(i) = maxval(GYL(:,:,i))
       GMIN(i) = minval(GYL(:,:,i))
    end do
    GMAX(j) = maxval(GMAX(2:4))
    GMIN(j) = minval(GMIN(2:4))

    j = j + 1
    k =  83
    call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,5), NMAX, STR(j), GYT%gnrm(k), &
         &        GMAX=GMAX(5), GMIN=GMIN(5)) ! TX
    k = 347 ; STR(j) = '@Bqhat$-z//$= TX, Neo. solver@'
    call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,6), NMAX, STR(j), GYT%gnrm(k), &
         &        GMAX=GMAX(6), GMIN=GMIN(6)) ! Neo. solver
    GMAX(j) = max(GMAX(5),GMAX(6))
    GMIN(j) = min(GMIN(5),GMIN(6))

    j = j + 1
    k =  61
    call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,7), NMAX, STR(j), GYT%gnrm(k), &
         &        GMAX=GMAX(7), GMIN=GMIN(7))
    k =  74
    call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,8), NMAX, STR(j), GYT%gnrm(k), &
         &        GMAX=GMAX(8), GMIN=GMIN(8))
    k =  87 ; STR(j) = '@u$-er$=,u$-ir$=,u$-zr$=@'
    call APPROPGY(MODEG, GYT%v(:,:,k), GYL(:,:,9), NMAX, STR(j), GYT%gnrm(k), &
         &        GMAX=GMAX(9), GMIN=GMIN(9))
    GMAX(j) = max(GMAX(7),GMAX(8),GMAX(9))
    GMIN(j) = min(GMIN(7),GMIN(8),GMIN(9))

    call SETFNT(44)
    do NG = 0, NGT
       call animes
       call SETLIN(-1, -1, 7)
       call gtextx(12.5,17.7,'@T=@',0)
       call gnumbr(13.1,17.7,GTX(NG),3,0)

       call TXGRFRS(0, GX(1:N), GYL(1:N,NG,0:1), N-1, 2, STR(0), MODE, IND, 0, 0, &
            &       'ANIME', GMAX(0), GMIN(0))
       call TXGRFRS(1, GX(0:N), GYL(0:N,NG,2:4), N  , 3, STR(1), MODE, IND, 0, 0, &
            &       'ANIME', GMAX(1), GMIN(1))
       call TXGRFRS(2, GX(0:N), GYL(0:N,NG,5:6), N  , 2, STR(2), MODE, IND, 0, 0, &
            &       'ANIME', GMAX(2), GMIN(2))
       call TXGRFRS(3, GX(0:N), GYL(0:N,NG,7:9), N  , 3, STR(3), MODE, IND, 0, 0, &
            &       'ANIME', GMAX(3), GMIN(3))
       call animee
    end do

    call PAGEE

    deallocate(GYL,STR,GMAX,GMIN)

  end subroutine TXGRCPA

  !***************************************************************
  !
  !   SLAVE ROUTINE TO Draw graph of GX, GY
  !
  !***************************************************************

  subroutine TXGRFRX(K, GXL, GYL, NRMAX, NGMAX, STR, MODE, IND, &
       &             GXMIN, GYMAX, GYMIN, ILOGIN, GPXY_IN)

    use tx_commons, only : rhob
    integer(4), intent(in) :: K, NRMAX, NGMAX, MODE, IND
    real(4), dimension(:), intent(in) :: GXL
    real(4), dimension(0:NRMAX,1:NGMAX+1), intent(in) :: GYL
    real(4), intent(in), optional :: GXMIN, GYMAX, GYMIN
    integer(4), intent(in), optional :: ILOGIN
    real(4), dimension(:), intent(in), optional :: GPXY_IN
    character(len=*), intent(in) :: STR
    integer(4) :: ILOG, IPRES
    real(4) :: GXMAX, GXMINL
    real(4), dimension(4) :: GPXY

    if(present(GPXY_IN)) then
       GPXY(1:4) = GPXY_IN(1:4)
    else
       GPXY(1) =  3.0 + 12.5 * mod(K,2)
       GPXY(2) = 12.5 + 12.5 * mod(K,2)
       GPXY(3) = 10.5 -  8.5 * real(K/2)
       GPXY(4) = 17.0 -  8.5 * real(K/2)
       ! square
!!$       GPXY(1) =  3.0 + 12.5 * mod(K,2)
!!$       GPXY(2) = 10.4286 + 12.5 * mod(K,2)
!!$       GPXY(3) = 10.5 -  8.5 * real(K/2)
!!$       GPXY(4) = 17.0 -  8.5 * real(K/2)
    end if

    GXMAX = real(rhob)

    if(present(GXMIN)) then
       GXMINL = GXMIN
    else
       GXMINL = 0.0
    end if
    if(present(ILOGIN)) then
       ! Semi-Log scale
       ILOG = ILOGIN
    else
       ! Normal scale
       ILOG = 0
    end if

    IPRES = 0
    if(present(GYMAX)) IPRES = IPRES + 1
    if(present(GYMIN)) IPRES = IPRES + 2

    if(IPRES == 0) then
       call TXGRAF(GPXY, GXL, GYL, NRMAX+1, NRMAX+1, NGMAX+1, &
            &            GXMINL, GXMAX, STR, 0.3, MODE, IND, ILOG)
    else if(IPRES == 1) then
       call TXGRAF(GPXY, GXL, GYL, NRMAX+1, NRMAX+1, NGMAX+1, &
            &            GXMINL, GXMAX, STR, 0.3, MODE, IND, ILOG, GYMAX_IN=GYMAX)
    else if(IPRES == 2) then
       call TXGRAF(GPXY, GXL, GYL, NRMAX+1, NRMAX+1, NGMAX+1, &
            &            GXMINL, GXMAX, STR, 0.3, MODE, IND, ILOG, GYMIN_IN=GYMIN)
    else
       call TXGRAF(GPXY, GXL, GYL, NRMAX+1, NRMAX+1, NGMAX+1, &
            &            GXMINL, GXMAX, STR, 0.3, MODE, IND, ILOG, GYMAX_IN=GYMAX, &
            &                                                      GYMIN_IN=GYMIN)
    end if

  end subroutine TXGRFRX

  !***************************************************************
  !
  !   SLAVE ROUTINE TO Draw graph of GX, GY (small version)
  !
  !***************************************************************

  subroutine TXGRFRXS(K, GXL, GYL, NRMAX, NGMAX, STR, MODE, IND, GYMAX, ILOGIN)

    use tx_commons, only : rhob
    integer(4), INTENT(IN) :: K, NRMAX, NGMAX, MODE, IND
    real(4), DIMENSION(:), INTENT(IN) :: GXL
    real(4), DIMENSION(0:NRMAX,1:NGMAX+1), INTENT(IN) :: GYL
    character(len=*), INTENT(IN) :: STR
    real(4), intent(in), optional :: GYMAX
    integer(4), intent(in), optional :: ILOGIN
    integer(4) :: ILOG
    real(4) :: GXMAX
    real(4), DIMENSION(4) :: GPXY

    GPXY(1) =   2.0 + 6.1  * mod(K,4)
    GPXY(2) =   6.7 + 6.1  * mod(K,4)
    GPXY(3) = 13.75 - 4.25 * real(K/4)
    GPXY(4) = 17.0  - 4.25 * real(K/4)
    GXMAX = real(rhob)

    if(present(ILOGIN)) then
       ILOG = ILOGIN
    else
       ILOG = 0
    end if

    if(present(GYMAX)) then
       call TXGRAF(GPXY, GXL, GYL, NRMAX+1, NRMAX+1, NGMAX+1, &
            &            0.0, GXMAX, STR, 0.26, MODE, IND, ILOG, GYMAX)
    else
       call TXGRAF(GPXY, GXL, GYL, NRMAX+1, NRMAX+1, NGMAX+1, &
            &            0.0, GXMAX, STR, 0.26, MODE, IND, ILOG)
    end if

  end subroutine TXGRFRXS

  !***************************************************************
  !
  !   SLAVE ROUTINE TO Draw graph of GX, GY at one slice time
  !
  !     All arguments are input
  !
  !***************************************************************

  subroutine TXGRFRS(K, GXL, GYL, NXMAX, NGMAX, STR, MODE, IND, ILOGIN, ISIZE, Kindex, &
       &             GMAX, GMIN)
    ! Forth argument, NXMAX, is not always "NRMAX" defined as the number of grid points.

    use tx_commons, only : rhob
    integer(4), INTENT(IN) :: K, NXMAX, NGMAX, MODE, IND, ISIZE
    real(4), intent(in), optional :: GMAX, GMIN
    real(4), DIMENSION(:), INTENT(IN) :: GXL
    real(4), DIMENSION(:,:), INTENT(IN) :: GYL
    character(len=*), INTENT(IN) :: STR, Kindex
    integer(4), intent(in), optional :: ILOGIN
    integer(4) :: ILOG
    real(4) :: GXMAX, FNTSIZE
    real(4), DIMENSION(4) :: GPXY

    FNTSIZE = 0.3
    if(ISIZE == 1) then
       ! Small
       GPXY(1) =   2.0 + 6.1  * mod(K,4)
       GPXY(2) =   6.7 + 6.1  * mod(K,4)
       GPXY(3) = 13.75 - 4.25 * real(K/4)
       GPXY(4) = 17.0  - 4.25 * real(K/4)
       FNTSIZE = 0.26
    else if(ISIZE == 2) then
       ! Square
       GPXY(1) =  3.0 + 12.5 * mod(K,2)
!?       GPXY(2) = 10.4286 + 12.5 * mod(K,2)
       GPXY(2) = 12.5 + 12.5 * mod(K,2)
       GPXY(3) = 10.5 -  8.5 * real(K/2)
       GPXY(4) = 17.0 -  8.5 * real(K/2)
    else if(ISIZE == 3) then
       ! One graph per page (K is no longer valid.)
       GPXY(1) =  3.0
       GPXY(2) = 23.0
       GPXY(3) =  2.0
       GPXY(4) = 17.0
    else if(ISIZE == 4) then
       ! Six graphs per page
       GPXY(1) =  1.8  + 8.45 * mod(K,3)
       GPXY(2) =  8.45 + 8.45 * mod(K,3)
       GPXY(3) = 10.55 - 7.55 * real(K/3)
       GPXY(4) = 15.1  - 7.55 * real(K/3)
    else if(ISIZE == 5) then
       ! Standard (lower set)
       GPXY(1) =  3.0 + 12.5 * mod(K,2)
       GPXY(2) = 12.5 + 12.5 * mod(K,2)
       GPXY(3) =  9.5 -  8.5 * real(K/2)
       GPXY(4) = 16.0 -  8.5 * real(K/2)
    else if(ISIZE == 6) then
       ! Six graphs per page (lower set)
       GPXY(1) =  1.8  + 8.45 * mod(K,3)
       GPXY(2) =  8.45 + 8.45 * mod(K,3)
       GPXY(3) =  8.55 - 7.55 * real(K/3)
       GPXY(4) = 13.1  - 7.55 * real(K/3)
    else
       ! Standard
       GPXY(1) =  3.0 + 12.5 * mod(K,2)
       GPXY(2) = 12.5 + 12.5 * mod(K,2)
       GPXY(3) = 10.5 -  8.5 * real(K/2)
       GPXY(4) = 17.0 -  8.5 * real(K/2)
    end if
    GXMAX = real(rhob)

    if(present(ILOGIN)) then
       ILOG = ILOGIN
    else
       ILOG = 0
    end if

    if(Kindex == 'STATIC') then
       call TXGRAF(GPXY, GXL, GYL(1:NXMAX+1,1:NGMAX), NXMAX+1, NXMAX+1, NGMAX, &
            &            0.0, GXMAX, STR, FNTSIZE, MODE, IND, ILOG, &
            &            GYMAX_IN=GMAX, GYMIN_IN=GMIN)
    else if(Kindex == 'ANIME') then
       call TXGRAF(GPXY, GXL, GYL, NXMAX+1, NXMAX+1, NGMAX, &
            &            0.0, GXMAX, STR, FNTSIZE, MODE, IND, ILOG, &
            &            GYMAX_IN=GMAX, GYMIN_IN=GMIN)!, Kindex=Kindex)
    else
       call TXGRAF(GPXY, GXL, GYL, NXMAX+1, NXMAX+1, NGMAX, &
            &            0.0, GXMAX, STR, FNTSIZE, MODE, IND, ILOG, &
            &            GYMAX_IN=GMAX, GYMIN_IN=GMIN)
    end if

  end subroutine TXGRFRS

  !***************************************************************
  !
  !   Draw graph of GVX, GVY
  !
  !***************************************************************

  subroutine TXGRFV(NGYV,MODE)

    use tx_commons, only : DT, gkilo, gmega
    integer(4), intent(IN) :: NGYV, MODE
    integer(4) :: IND, NG
    character(len=60) ::  STR
    real(4), dimension(:,:), allocatable :: GVYL

    if (NGVV <= 1) then
       write(6,*) 'GV', NGYV, ' has no data'
       return
    end if

    if (MODEG == 2) then
       IND = 9
    else
       IND = 0
    end if

    allocate(GVYL(0:NGVV,1:4))

    call PAGES
    call SETCHS(0.3, 0.0)
    call SETLIN(0, 1, 7)

    call MOVE(2.0,17.7)
    call TEXT('[GV', 3)
    call NUMBI(NGYV, '(I2)', 2)
    call TEXT(']  ', 3)
    call TEXT('FROM', 4)
    call NUMBR(GVX(0), '(ES9.2)', 9)
    call TEXT(' TO', 3)
    call NUMBR(GVX(NGVV), '(ES9.2)', 9)
    call TEXT('  DT =', 6)
    call NUMBD(DT, '(ES9.2)', 9)
    call TEXT('  NGVSTP = ', 11)
    call NUMBI(NGVSTP,'(I4)',4)

    select case(NGYV)
    case(1)
       STR = '@n$-e$=(0),LAVE(0),(0.24),(0.6)@'
       do NG = 0, NGVV
          GVYL(NG,1)   = GVY(NG,1)
          GVYL(NG,2:4) = GVY(NG,64:66)
       end do
       call APPROPGY(MODEG, GVYL, GVYL, 4*(NGVV+1), STR, 1e20)
       call TXGREVM(0, GVX, GVYL, NGVV, 4, STR, MODE, IND)

       STR = '@Z*n$-i$=+Z*n$-b$=-n$-e$=(0)@'
       call APPROPGY(MODEG, GVY(:,1), GVYL, NGVV+1, STR, 1e14)
       call TXGREVM(1, GVX, GVYL, NGVV, 1, STR, MODE, IND)

       STR = '@E$-r$=(a/2)@'
       call APPROPGY(MODEG, GVY(:,31), GVYL, NGVV+1, STR, gkilo)
       call TXGREVM(2, GVX, GVYL, NGVV, 1, STR, MODE, IND)

       call TXWPGR

    case(2)
       STR = '@u$-er$=(a/2)@'
       call TXGREVM(0, GVX, GVY(:,5), NGVV, 1, STR, MODE, IND)

       STR = '@u$-e$#q$#$=(a/2)@'
       call APPROPGY(MODEG, GVY(:,6), GVYL, NGVV+1, STR, gkilo)
       call TXGREVM(1, GVX, GVYL, NGVV, 1, STR, MODE, IND)

       STR = '@u$-e$#f$#$=(0)@'
       call APPROPGY(MODEG, GVY(:,7), GVYL, NGVV+1, STR, gkilo)
       call TXGREVM(2, GVX, GVYL,  NGVV, 1, STR, MODE, IND)

       STR = '@j$-r$=(a/2)@'
       call TXGREVM(3, GVX, GVY(:,3), NGVV, 1, STR, MODE, IND)

       !call TXWPGR

    case(3)
       STR = '@u$-ir$=(a/2)@'
       call TXGREVM(0, GVX, GVY(:,14), NGVV, 1, STR, MODE, IND)

       STR = '@u$-i$#q$#$=(a/2)@'
       call APPROPGY(MODEG, GVY(:,15), GVYL, NGVV+1, STR, gkilo)
       call TXGREVM(1, GVX, GVYL, NGVV, 1, STR, MODE, IND)

       STR = '@u$-i$#f$#$=(a/2)@'
       call APPROPGY(MODEG, GVY(:,16), GVYL, NGVV+1, STR, gkilo)
       call TXGREVM(2, GVX, GVYL, NGVV, 1, STR, MODE, IND)

       call TXWPGR

    case(4)
       STR = '@u$-zr$=(a/2)@'
       call TXGREVM(0, GVX, GVY(:,23), NGVV, 1, STR, MODE, IND)

       STR = '@u$-z$#q$#$=(a/2)@'
       call APPROPGY(MODEG, GVY(:,24), GVYL, NGVV+1, STR, gkilo)
       call TXGREVM(1, GVX, GVYL, NGVV, 1, STR, MODE, IND)

       STR = '@u$-z$#f$#$=(a/2)@'
       call APPROPGY(MODEG, GVY(:,25), GVYL, NGVV+1, STR, gkilo)
       call TXGREVM(2, GVX, GVYL, NGVV, 1, STR, MODE, IND)

       call TXWPGR

    case(5)
       STR = '@q(0)@'
       call TXGREVM(0, GVX, GVY(:,36), NGVV, 1, STR, MODE, IND)

       STR = '@E$-tor$=(b)@'
       call TXGREVM(1, GVX, GVY(:,32), NGVV, 1, STR, MODE, IND)

       STR = '@j$-$#f$#$=(0)@'
       call APPROPGY(MODEG, GVY(:,37), GVYL, NGVV+1, STR, gmega)
       call TXGREVM(2, GVX, GVYL, NGVV, 1, STR, MODE, IND)

       STR = '@B$-$#q$#$=(a)@'
       call TXGREVM(3, GVX, GVY(:,34), NGVV, 1, STR, MODE, IND)


    case(6)
       STR = '@n$-b$=(0)@'
       call APPROPGY(MODEG, GVY(:,42), GVYL, NGVV+1, STR, 1e18)
       call TXGREVM(0, GVX, GVYL, NGVV, 1, STR, MODE, IND)

       STR = '@u$-b$#f$#$=(0)@'
       call APPROPGY(MODEG, GVY(:,43), GVYL, NGVV+1, STR, gkilo)
       call TXGREVM(1, GVX, GVYL, NGVV, 1, STR, MODE, IND)

       STR = '@Ripple n$-b$=(0)@'
       call APPROPGY(MODEG, GVY(:,45), GVYL, NGVV+1, STR, gkilo)
       call TXGREVM(2, GVX, GVYL, NGVV, 1, STR, MODE, IND)

       STR = '@Bu$-b//$=(a/2)@'
       call APPROPGY(MODEG, GVY(:,44), GVYL, NGVV+1, STR, gkilo)
       call TXGREVM(3, GVX, GVYL, NGVV, 1, STR, MODE, IND)

    case(7)
       STR = '@j$-e$#f$#$=(0)@'
       call APPROPGY(MODEG, GVY(:,38), GVYL, NGVV+1, STR, gmega)
       call TXGREVM(0, GVX, GVYL, NGVV, 1, STR, MODE, IND)

       STR = '@j$-i$#f$#$=(0)@'
       call APPROPGY(MODEG, GVY(:,39), GVYL, NGVV+1, STR, gkilo)
       call TXGREVM(1, GVX, GVYL, NGVV, 1, STR, MODE, IND)

       STR = '@j$-z$#f$#$=(0)@'
       call APPROPGY(MODEG, GVY(:,40), GVYL, NGVV+1, STR, gkilo)
       call TXGREVM(2, GVX, GVYL, NGVV, 1, STR, MODE, IND)

       STR = '@j$-b$#f$#$=(0)@'
       call APPROPGY(MODEG, GVY(:,41), GVYL, NGVV+1, STR, gkilo)
       call TXGREVM(3, GVX, GVYL, NGVV, 1, STR, MODE, IND)

    case(8)
       STR = '@u$-er$=(a/2)@'
       call TXGREVM(0, GVX, GVY(:, 5), NGVV, 1, STR, MODE, IND)

       STR = '@u$-e$#$/136$#$=(a/2)@'
       call APPROPGY(MODEG, GVY(:,10), GVYL, NGVV+1, STR, gkilo)
       call TXGREVM(1, GVX, GVYL, NGVV, 1, STR, MODE, IND)

       STR = '@Bu$-e//$=(a/2)@'
       call APPROPGY(MODEG, GVY(:,11), GVYL, NGVV+1, STR, gkilo)
       call TXGREVM(2, GVX, GVYL, NGVV, 1, STR, MODE, IND)

       call TXWPGR

    case(9)
       STR = '@u$-ir$=(a/2)@'
       call TXGREVM(0, GVX, GVY(:,14), NGVV, 1, STR, MODE, IND)

       STR = '@u$-i$#$/136$#$=(a/2)@'
       call APPROPGY(MODEG, GVY(:,19), GVYL, NGVV+1, STR, gkilo)
       call TXGREVM(1, GVX, GVYL, NGVV, 1, STR, MODE, IND)

       STR = '@Bu$-i//$=(a/2)@'
       call APPROPGY(MODEG, GVY(:,20), GVYL, NGVV+1, STR, gkilo)
       call TXGREVM(2, GVX, GVYL, NGVV, 1, STR, MODE, IND)

       call TXWPGR

    case(10)
       STR = '@u$-zr$=(a/2)@'
       call TXGREVM(0, GVX, GVY(:,23), NGVV, 1, STR, MODE, IND)

       STR = '@u$-z$#$/136$#$=(a/2)@'
       call APPROPGY(MODEG, GVY(:,28), GVYL, NGVV+1, STR, gkilo)
       call TXGREVM(1, GVX, GVYL, NGVV, 1, STR, MODE, IND)

       STR = '@Bu$-z//$=(a/2)@'
       call APPROPGY(MODEG, GVY(:,29), GVYL, NGVV+1, STR, gkilo)
       call TXGREVM(2, GVX, GVYL, NGVV, 1, STR, MODE, IND)

       call TXWPGR

    case(11)
       STR = '@D$-eff$=(a/2)@'
       call TXGREVM(0, GVX, GVY(:,51), NGVV, 1, STR, MODE, IND)

       STR = '@D$-e$=(a/2)+D$-i$=(a/2)+D$-z$=(a/2)@'
       call TXGREVM(1, GVX, GVY(:,50), NGVV, 1, STR, MODE, IND)

       STR = '@$#c$#$-tbe$=(a/2)@'
       call TXGREVM(2, GVX, GVY(:,52), NGVV, 1, STR, MODE, IND)

       STR = '@$#c$#$-tbi$=(a/2)@'
       call TXGREVM(3, GVX, GVY(:,53), NGVV, 1, STR, MODE, IND)

!!$       STR = '@G$-1$=h$+2$=(a/2)@'
!!$       call TXGREVM(1, GVX, GVY(:,55), NGVV, 1, STR, MODE, IND)
!!$
!!$       STR = '@F$-CDBM$=(a/2)@'
!!$       call TXGREVM(2, GVX, GVY(:,56), NGVV, 1, STR, MODE, IND)

    case(12)
       STR = '@T$-e$=(0)@'
       call TXGREVM(0, GVX, GVY(:,8), NGVV, 1, STR, MODE, IND)

       STR = '@T$-i$=(0)@'
       call TXGREVM(1, GVX, GVY(:,17), NGVV, 1, STR, MODE, IND)

       STR = '@T$-z$=(0)@'
       call TXGREVM(2, GVX, GVY(:,26), NGVV, 1, STR, MODE, IND)

       call TXWPGR

    case(13)
       STR = '@p$-e$=(0)@'
       call TXGREVM(0, GVX, GVY(:,12), NGVV, 1, STR, MODE, IND)

       STR = '@p$-i$=(0)@'
       call TXGREVM(1, GVX, GVY(:,21), NGVV, 1, STR, MODE, IND)

       STR = '@p$-z$=(0)@'
       call TXGREVM(2, GVX, GVY(:,30), NGVV, 1, STR, MODE, IND)

       call TXWPGR

    case(14)
       STR = '@s(a/2)@'
       call TXGREVM(0, GVX, GVY(:,57), NGVV, 1, STR, MODE, IND)

       STR = '@$#a$#(a/2)@'
       call TXGREVM(1, GVX, GVY(:,58), NGVV, 1, STR, MODE, IND)

       STR = '@$#k$#(a/2)@'
       call TXGREVM(2, GVX, GVY(:,59), NGVV, 1, STR, MODE, IND)

       call TXWPGR

    case(15)
       STR = '@BE$-pol$=(a/2)@'
       call TXGREVM(0, GVX, GVY(:,33), NGVV, 1, STR, MODE, IND)

       STR = '@I(0)@'
       call TXGREVM(1, GVX, GVY(:,35), NGVV, 1, STR, MODE, IND)

       STR = '@uhat$-i$#q$#$=(a/2)@'
       call APPROPGY(MODEG, GVY(:,9), GVYL, NGVV+1, STR, gkilo)
       call TXGREVM(2, GVX, GVYL, NGVV, 1, STR, MODE, IND)

       STR = '@1+2qhat$+2$=@'
       call TXGREVM(3, GVX, GVY(:,70), NGVV, 1, STR, MODE, IND)

    case(16)
       STR = '@SLOW N$-0$=(a)@'
       call APPROPGY(MODEG, GVY(:,46), GVYL, NGVV+1, STR, 1e15)
       call TXGREVM(0, GVX, GVYL, NGVV, 1, STR, MODE, IND)

       STR = '@THERMAL N$-0$=(0)@'
       call APPROPGY(MODEG, GVY(:,47), GVYL, NGVV+1, STR, 1e13)
       call TXGREVM(1, GVX, GVYL, NGVV, 1, STR, MODE, IND)

       STR = '@HALO N$-0$=@'
       call APPROPGY(MODEG, GVY(:,48), GVYL, NGVV+1, STR, 1e13)
       call TXGREVM(2, GVX, GVYL, NGVV, 1, STR, MODE, IND)

       STR = '@TOTAL N$-0$=(a)@'
       call APPROPGY(MODEG, GVY(:,49), GVYL, NGVV+1, STR, 1e16)
       call TXGREVM(3, GVX, GVYL, NGVV, 1, STR, MODE, IND)

    case(17)
       STR = '@PIE@'
       call TXGREVM(0, GVX, GVY(:,60), NGVV, 1, STR, MODE, IND)

       STR = '@PCX@'
       call APPROPGY(MODEG, GVY(:,61), GVYL, NGVV+1, STR, gkilo)
       call TXGREVM(1, GVX, GVYL, NGVV, 1, STR, MODE, IND)

       STR = '@SIE@'
       call APPROPGY(MODEG, GVY(:,62), GVYL, NGVV+1, STR, 1e20)
       call TXGREVM(2, GVX, GVYL, NGVV, 1, STR, MODE, IND)

       STR = '@PBr@'
       call TXGREVM(3, GVX, GVY(:,63), NGVV, 1, STR, MODE, IND)

    case(18)
       STR = '@n$-e$=:total@'
       call TXGREVM(0, GVX, GVY(:,67), NGVV, 1, STR, MODE, IND)

       STR = '@$#g$#:total@'
       call TXGREVM(1, GVX, GVY(:,68), NGVV, 1, STR, MODE, IND)

       STR = '@$#t$#$-p$=@'
       call TXGREVM(2, GVX, GVY(:,69), NGVV, 1, STR, MODE, IND)

       call TXWPGR

    case DEFAULT
       write(6,*) 'Unknown NGYV: NGYV = ',NGYV
    end select

    call PAGEE

    deallocate(GVYL)

  end subroutine TXGRFV

  !***************************************************************
  !
  !   Draw graph of GTX, GTY
  !
  !***************************************************************

  subroutine TXGRFT(NGYT,MODE)

    use tx_commons, only : DT
    integer(4), intent(IN) :: NGYT, MODE
    integer(4) :: IND, i, NG
    real(4), dimension(:,:), allocatable :: GTYL
    character(len=60) ::  STR
    real(4) :: GYMIN

    if (NGT <= 1) then
       write(6,*) 'GT', NGYT, ' has no data'
       return
    end if

    if (MODEG == 2) then
       IND = 9
    else
       IND = 0
    end if

    call PAGES
    call SETCHS(0.3, 0.0)
    call SETLIN(0, 1, 7)

    call MOVE(2.0,17.7)
    call TEXT('[GT', 3)
    call NUMBI(NGYT, '(I2)', 2)
    call TEXT(']  ', 3)
    call TEXT('FROM', 4)
    call NUMBR(GTX(0), '(ES9.2)', 9)
    call TEXT(' TO', 3)
    call NUMBR(GTX(NGT), '(ES9.2)', 9)
    call TEXT('  DT =', 6)
    call NUMBD(DT, '(ES9.2)', 9)
    call TEXT('  NGTSTP = ', 11)
    call NUMBI(NGTSTP,'(I4)',4)

    select case(NGYT)
    case(1)
       allocate(GTYL(0:NGT,3))
       do NG = 0, NGT
          do i = 1, 3
             GTYL(NG,i) = GTY(NG,7+(i-1)) / GTY(NG,10+(i-1))
          end do
       end do
       STR = '@Ne0,Ni0,Nz0,{Ne},{Ni},{Nz} [10$+20$=/m$+3$=] vs t@'
       call TXGREVM(0, GTX, GTY(0:NGT,7:), NGT, 6, STR, MODE, IND)

       STR = '@Te0,Ti0,Tz0,{Te},{Ti},{Tz} [keV] vs t@'
       call TXGREVM(1, GTX, GTY(0:NGT, 1:), NGT, 6, STR, MODE, IND)

       STR = '@Ne0/{Ne},Ni0/{Ni},Nz0/{Nz} vs t@'
       call TXGREVM(2, GTX, GTYL, NGT, 3, STR, MODE, IND)

       do NG = 0, NGT
          GTYL(NG,1) = GTY(NG,10)
          GTYL(NG,2) = GTY(NG,11)
          GTYL(NG,3) = GTY(NG,12) * 100
       end do
       STR = '@{Ne},{Ni},{Nz}*100 vs t@'
       call TXGREVM(3, GTX, GTYL, NGT, 3, STR, MODE, IND)

       deallocate(GTYL)

    case(2)
       STR = '@IP,IOH,INB,IBS [MA] vs t@'
       call TXGREVM(0, GTX, GTY(0:NGT,26:), NGT, 4, STR, MODE, IND)

       STR = '@PIN,POH,PNB,PRF,PNF [MW] vs t@'
       call TXGREVM(1, GTX, GTY(0:NGT,13:), NGT, 5, STR, MODE, IND)

       STR = '@QF vs t@'
       call TXGREVM(2, GTX, GTY(0:NGT,23), NGT, 1, STR, MODE, IND)

       STR = '@POUT,PCX,PIE [MW] vs t@'
       call TXGREVM(3, GTX, GTY(0:NGT,18:), NGT, 3, STR, MODE, IND)

    case(3)
       STR = '@WF,WB,We,Wi,Wz [MJ] vs t@'
       call TXGREVM(0, GTX, GTY(0:NGT,31:), NGT, 5, STR, MODE, IND)

       STR = '@TAUE1,TAUE2,TAUEP [s] vs t@'
       call TXGREVM(1, GTX, GTY(0:NGT,36:), NGT, 3, STR, MODE, IND)

       STR = '@BETAa,BETA0 [%] vs t@'
       call TXGREVM(2, GTX, GTY(0:NGT,39:), NGT, 2, STR, MODE, IND)

       STR = '@BETAPa,BETAP0 vs t@'
       call TXGREVM(3, GTX, GTY(0:NGT,41:), NGT, 2, STR, MODE, IND)

    case(4)
       allocate(GTYL(0:NGT,4))
       do NG = 0, NGT
          GTYL(NG,1) = GTY(NG,43)
          GTYL(NG,2) = sum(GTY(NG,66:68))
          GTYL(NG,3) = sum(GTY(NG,66:67))
          GTYL(NG,4) = GTY(NG,61)
       end do
       STR = '@VLOOP, VPoyntRHS, VPoynt ind+res, VLOOPpol [V]@'
       GYMIN = min(minval(GTYL(:,1)),minval(GTYL(:,4)))
       call TXGREVM(0, GTX, GTYL, NGT, 4, STR, MODE, IND, GYMIN=GYMIN)

       STR = '@ALI@'
       call TXGREVM(1, GTX, GTY(0:NGT,44), NGT, 1, STR, MODE, IND)

       STR = '@Q(0)@'
       call TXGREVM(2, GTX, GTY(0:NGT,45), NGT, 1, STR, MODE, IND)

       STR = '@RQ1 [m]@'
       call TXGREVM(3, GTX, GTY(0:NGT,46), NGT, 1, STR, MODE, IND)

       deallocate(GTYL)

    case(5)
       allocate(GTYL(0:NGT,5))
       do NG = 0, NGT
          GTYL(NG,1:4) = GTY(NG,62:65)
          GTYL(NG,5) = sum(GTY(NG,63:65))
       end do
       STR = '@TAUP, TAUPA [s] vs t@'
       call TXGREVM(0, GTX, GTY(0:NGT,50:), NGT, 2, STR, MODE, IND)

       STR = '@Gamma$-a$= [10$+20$=/m$+2$=s$+1$=] vs t@'
       call TXGREVM(1, GTX, GTY(0:NGT,52), NGT, 1, STR, MODE, IND)

       STR = '@integrated {L} [Nms] vs t@'
       call TXGREVM(2, GTX, GTY(0:NGT,55), NGT, 1, STR, MODE, IND)

       STR = '@volt-sec balance [Wb] vs t@'
       call TXGREVM(3, GTX, GTYL, NGT, 5, STR, MODE, IND)

       deallocate(GTYL)

    case(6)
       STR = '@Te0,Ti0,Tz0,{Te},{Ti},{Tz} [keV] vs t@'
       call TXGREVS(0, GTX, GTY(0:NGT, 1:), NGT, 6, STR, IND)

       STR = '@IP,IOH,INB,IBS,ITOT [MA] vs t@'
       call TXGREVS(2, GTX, GTY(0:NGT,26:), NGT, 5, STR, IND)

       STR = '@PIN,POH,PNB,PRF,PNF [MW] vs t@'
       call TXGREVS(4, GTX, GTY(0:NGT,13:), NGT, 5, STR, IND)

       STR = '@POUT,PCX,PIE [MW] vs t@'
       call TXGREVS(6, GTX, GTY(0:NGT,18:), NGT, 3, STR, IND)

       STR = '@QF vs t@'
       call TXGREVS(1, GTX, GTY(0:NGT,23), NGT, 1, STR, IND)

       STR = '@{Te},{Ti},{Tz} [keV] vs t@'
       call TXGREVS(3, GTX, GTY(0:NGT, 4:), NGT, 3, STR, IND)

       STR = '@Te0,Ti0,Tz0 [keV] vs t@'
       call TXGREVS(5, GTX, GTY(0:NGT, 1:), NGT, 3, STR, IND)

       STR = '@Ne0,Ni0,Nz0,{Ne},{Ni},{Nz} [10$+20$=/m$+3$=] vs t@'
       call TXGREVS(7, GTX, GTY(0:NGT, 7:), NGT, 6, STR, IND)

    case(7)
       allocate(GTYL(0:NGT,4))
       do NG = 0, NGT
          GTYL(NG,1) = GTY(NG,43)
          GTYL(NG,2) = sum(GTY(NG,66:68))
          GTYL(NG,3) = sum(GTY(NG,66:67))
          GTYL(NG,4) = GTY(NG,61)
       end do
       STR = '@WF,WB,We,Wi,Wz [MJ] vs t@'
       call TXGREVS(0, GTX, GTY(0:NGT,31:), NGT, 5, STR, IND)

       STR = '@BETAa,BETA0 [%] vs t@'
       call TXGREVS(2, GTX, GTY(0:NGT,39:), NGT, 2, STR, IND)

       STR = '@TAUE1,TAUE2,TAUEP [s] vs t@'
       call TXGREVS(4, GTX, GTY(0:NGT,36:), NGT, 3, STR, IND)

       STR = '@BETAPa,BETAP0 vs t@'
       call TXGREVS(6, GTX, GTY(0:NGT,41:), NGT, 2, STR, IND)

       STR = '@VLOOP, VPoyntRHS, VPoynt ind+res, VLOOPpol [V]@'
       GYMIN = min(minval(GTYL(0:NGT,1)),minval(GTYL(0:NGT,4)))
       call TXGREVS(1, GTX, GTYL, NGT, 4, STR, IND, GYMIN=GYMIN)

       STR = '@Q(0)@'
       call TXGREVS(3, GTX, GTY(:,45), NGT, 1, STR, IND)

       STR = '@ALI@'
       call TXGREVS(5, GTX, GTY(:,44), NGT, 1, STR, IND)

       STR = '@RQ1 [m]@'
       call TXGREVS(7, GTX, GTY(:,46), NGT, 1, STR, IND)

       deallocate(GTYL)

    case(8)
       STR = '@Nb0,{Nb} [10$+18$=/m$+3$=] vs t@'
       call TXGREVS(0, GTX, GTY(0:NGT,47:), NGT, 2, STR, IND)

       STR = '@neutrality [/m$+3$=] vs t@'
       call TXGREVS(1, GTX, GTY(:,49), NGT, 1, STR, IND)

       STR = '@TAUP, TAUPA [s] vs t@'
       call TXGREVS(2, GTX, GTY(0:NGT,50:), NGT, 2, STR, IND)

       STR = '@Gamma$-a$= [10$+20$=/m$+2$=s$+1$=] vs t@'
       call TXGREVS(3, GTX, GTY(:,52), NGT, 1, STR, IND)

       STR = '@TNBcol, TTqt [Nm] vs t@'
       call TXGREVS(4, GTX, GTY(0:NGT,53:), NGT, 2, STR, IND)

       STR = '@CEjima vs t@'
       call TXGREVS(5, GTX, GTY(:,69), NGT, 1, STR, IND, GYMAX=1.0, GYMIN=0.0)

       STR = '@Num. of Neutrals vs t@'
       call TXGREVS(6, GTX, GTY(:,57), NGT, 1, STR, IND)

       STR = '@VPoynt ind,res,pol [Wb] vs t@'
       call TXGREVS(7, GTX, GTY(0:NGT,66:), NGT, 3, STR, IND)

    case(9)
       allocate(GTYL(0:NGT,4))
       do NG = 0, NGT
          GTYL(NG,1) = GTY(NG,10)
          GTYL(NG,2) = GTY(NG,11)
          GTYL(NG,3) = GTY(NG,12)
          GTYL(NG,4) = GTY(NG,48) * 1.e-2
       end do
       STR = '@{Ne},{Ni},{Nz},{Nb} [10$+20$=/m$+3$=] vs t@'
       call TXGREVM(0, GTX, GTYL, NGT, 4, STR, MODE, IND)

       STR = '@integrated {j.grad psi} [Nm] vs t@'
       call TXGREVM(1, GTX, GTY(:,56), NGT, 1, STR, MODE, IND)

       STR = '@neutrality [/m$+3$=] vs t@'
       call TXGREVM(2, GTX, GTY(:,49), NGT, 1, STR, MODE, IND)

       STR = '@Num. of Neutrals vs t@'
       call TXGREVM(3, GTX, GTY(0:NGT,57:), NGT, 4, STR, MODE, IND)

       deallocate(GTYL)

    case DEFAULT
       write(6,*) 'Unknown NGYT: NGYT = ',NGYT
    end select

    call PAGEE

  end subroutine TXGRFT

  !***************************************************************
  !
  !   draw medium-sized graph of time history of quantities 
  !     TXGREVM: TX GRaph EVolution Medium
  !              ^^ ^^    ^^        ^
  !     Note: 1st actual argument of GYIN can be written as ":" if NG=1.
  !           However, it must be given as, say, "0:NGT" if NG>1.
  !***************************************************************

  subroutine TXGREVM(K, GXL, GYIN, NMAX, NG, STR, MODE, IND, &
       &             GYMAX, GYMIN)

    integer(4), intent(in) :: K, NMAX, NG, MODE, IND
    real(4), dimension(NMAX+1), intent(in) :: GXL
    real(4), dimension((NMAX+1)*NG), intent(in) :: GYIN
    real(4), intent(in), optional :: GYMAX, GYMIN
    character(len=*), intent(in) ::  STR
    integer(4) :: IPRES, NMAXX
    real(4), dimension(4) :: GPXY
    real(4), dimension(:,:), allocatable :: GYL
   
    NMAXX = NMAX + 1
    allocate(GYL(NMAXX,NG))
    GYL = reshape(GYIN, [NMAXX,NG])

    GPXY(1) =  3.0 + 12.5 * mod(K,2)
    GPXY(2) = 12.5 + 12.5 * mod(K,2)
    GPXY(3) = 10.5 -  8.5 * real(K/2)
    GPXY(4) = 17.0 -  8.5 * real(K/2)

    IPRES = 0
    if(present(GYMAX)) IPRES = IPRES + 1
    if(present(GYMIN)) IPRES = IPRES + 2

    if(IPRES == 0) then
       call TXGRAF(GPXY, GXL, GYL, NMAXX, NMAXX, NG, &
            &            GXL(1), GXL(NMAXX), STR, 0.3, MODE, IND, 0)
    else if(IPRES == 1) then
       call TXGRAF(GPXY, GXL, GYL, NMAXX, NMAXX, NG, &
            &            GXL(1), GXL(NMAXX), STR, 0.3, MODE, IND, 0, GYMAX_IN=GYMAX)
    else if(IPRES == 2) then
       call TXGRAF(GPXY, GXL, GYL, NMAXX, NMAXX, NG, &
            &            GXL(1), GXL(NMAXX), STR, 0.3, MODE, IND, 0, GYMIN_IN=GYMIN)
    else
       call TXGRAF(GPXY, GXL, GYL, NMAXX, NMAXX, NG, &
            &            GXL(1), GXL(NMAXX), STR, 0.3, MODE, IND, 0, GYMAX_IN=GYMAX, &
            &                                                        GYMIN_IN=GYMIN)
    end if

    deallocate(GYL)

  end subroutine TXGREVM

  !***************************************************************
  !
  !   draw small-sized graph of time history of quantities 
  !     TXGREVS: TX GRaph EVolution Small
  !              ^^ ^^    ^^        ^
  !     Note: 1st actual argument of GYIN can be written as ":" if NG=1.
  !           However, it must be given as, say, "0:NGT" if NG>1.
  !***************************************************************

  subroutine TXGREVS(K, GXL, GYIN, NMAX, NG, STR, IND, &
       &             GYMAX, GYMIN)

    integer(4), intent(in) :: K, NMAX, NG, IND
    real(4), dimension(NMAX+1), intent(in) :: GXL
    real(4), dimension((NMAX+1)*NG), intent(in) :: GYIN
    real(4), intent(in), optional :: GYMAX, GYMIN
    character(len=*), intent(in) ::  STR
    integer(4) :: IPRES, NMAXX
    real(4), dimension(4) :: GPXY
    real(4), dimension(:,:), allocatable :: GYL
   
    NMAXX = NMAX + 1
    allocate(GYL(NMAXX,NG))
    GYL = reshape(GYIN, [NMAXX,NG])

    GPXY(1) =  3.0 + 12.0 * mod(K,2)
    GPXY(2) = 12.0 + 12.0 * mod(K,2)
    GPXY(3) = 14.0 -  4.3 * real(K/2)
    GPXY(4) = 17.0 -  4.3 * real(K/2)

    IPRES = 0
    if(present(GYMAX)) IPRES = IPRES + 1
    if(present(GYMIN)) IPRES = IPRES + 2

    if(IPRES == 0) then
       call TXGRAF(GPXY, GXL, GYL, NMAXX, NMAXX, NG, &
            &            GXL(1), GXL(NMAXX), STR, 0.3, 1, IND, 0)
    else if(IPRES == 1) then
       call TXGRAF(GPXY, GXL, GYL, NMAXX, NMAXX, NG, &
            &            GXL(1), GXL(NMAXX), STR, 0.3, 1, IND, 0, GYMAX_IN=GYMAX)
    else if(IPRES == 2) then
       call TXGRAF(GPXY, GXL, GYL, NMAXX, NMAXX, NG, &
            &            GXL(1), GXL(NMAXX), STR, 0.3, 1, IND, 0, GYMIN_IN=GYMIN)
    else
       call TXGRAF(GPXY, GXL, GYL, NMAXX, NMAXX, NG, &
            &            GXL(1), GXL(NMAXX), STR, 0.3, 1, IND, 0, GYMAX_IN=GYMAX, &
            &                                                     GYMIN_IN=GYMIN)
    end if

    deallocate(GYL)

  end subroutine TXGREVS

  !***************************************************************
  !
  !   Graph of each TERM
  !
  !***************************************************************

  subroutine TXGRFQ(NQ,ID)

    use tx_commons, only : NRMAX, NQM, NCM, NLCMAX, rhob, NRC, NRA
    use tx_interface, only : APTOS

    integer(4), INTENT(IN) :: NQ, ID
    integer(4) :: NC, NSTR, IND
    real(4) :: GXMAX
    real(4), DIMENSION(0:NRMAX,1:NCM) :: GQYL
    real(4), DIMENSION(1:4) :: GPXY
    real(4), DIMENSION(1:4,1:5) :: GPXYA
    character(len=80), DIMENSION(1:NQM) :: STRGQ
    character(len=80) :: STR

    !        Title
    DATA STRGQ /'$#f$#$',"$#y$#$-t$='","$#y$#'",'$#y$#','$#y$#$-t$=',  &
         &      'n$-e$=','n$-e$=u$-e$=$+V$=', &
         &      'Bu$-e//$=','n$-e$={Ru$-e$#f$#$=}','n$-e$=T$-e$=', &
         &      'Bq$-e//$=', 'n$-e$={u$-e$#f$#$=/R}', 'BV$-1e$=', &
         &      'n$-i$=','n$-i$=u$-i$=$+V$=', &
         &      'Bu$-i//$=','n$-i$={Ru$-i$#f$#$=}','n$-i$=T$-i$=', &
         &      'Bq$-i//$=', 'n$-i$={u$-i$#f$#$=/R}', 'BV$-1i$=', &
         &      'n$-z$=','n$-z$=u$-z$=$+V$=', &
         &      'Bu$-z//$=','n$-z$={Ru$-z$#f$#$=}','n$-z$=T$-z$=', &
         &      'Bq$-z//$=', 'n$-z$={u$-z$#f$#$=/R}', 'BV$-1z$=', &
         &      'n$-b$=','n$-b$=u$-b$=$+V$=','Bu$-b//$=', &
         &      'n$-b$={Ru$-b$#f$#$=}','n$-b$={u$-b$#f$#$=/R}', 'BV$-1b$=', &
         &      'Slow n$-0$=', 'Thermal n$-0$=', 'Halo n$-0$=', 'n$-0z$=', &
         &      'Ripple n$-b$='/

    !        Graph coordinate
    DATA GPXYA/ 2.5, 9.5, 10.5,17.0, &
         &      2.5, 9.5,  1.5, 8.0, &
         &      15.0,22.0,10.5,17.0, &
         &      15.0,22.0, 1.5, 8.0, &
         &      2.5 ,22.0, 1.5,17.5/

    GPXY(1) = GPXYA(1,ID)
    GPXY(2) = GPXYA(2,ID)
    GPXY(3) = GPXYA(3,ID)
    GPXY(4) = GPXYA(4,ID)

    NSTR = 0
    call APTOS(STR,NSTR,NQ)
    call APTOS(STR,NSTR, ': ',2)
    call APTOS(STR,NSTR, STRGQ(NQ), len_trim(STRGQ(NQ)))

    GQYL(0:NRMAX,1:NLCMAX(NQ)) = GQY(0:NRMAX,1:NLCMAX(NQ),NQ)

    if (MODEG == 2) then
       IND = 9
    else
       IND = 0
    end if
    GXMAX = real(rhob)
    call TXGRAF(GPXY, GX, GQYL, NRMAX+1, NRMAX+1, NLCMAX(NQ), &
         &      0.0, GXMAX, '@'//STR(1:NSTR)//'@', 0.3, 2, IND, 0)
!         &      0.0, GXMAX, '@'//STR(1:NSTR)//'@', 0.3, 4, IND, 0)
    write(6,'(A3,I2,X,A4,6X,A3,3(11X,A7))') 'EQ=',NQ,'term','sum','NRC','NRA','NRMAX'
    do nc=1,nlcmax(nq)
       write(6,'(5X,I3,4(1ES15.5,3X))') nc, sum(gqyl(0:NRMAX,nc)) &
            & ,gqyl(NRC,nc), gqyl(NRA,nc), gqyl(NRMAX,nc)
    end do

  end subroutine TXGRFQ

  !***************************************************************
  !
  !   Write parameter to graphic screen
  !
  !***************************************************************

  subroutine TXWPGR

    use tx_commons, only : SLID, PNBCD, BB, rIp, FSDFIX, FSANOM, FSBOHM, FSPCL, &
         &              PROFD, PROFC, FSCX, FSRP, FSLC, FSNC, FSLP, FSLPB, FSION, FSD02, &
         &              PNBHT1, PNBHT2, PNBHP, PNBHex, PRFHe, PRFHi, Vb, Dfs0, rMus0, &
         &              Chis0, PTe0, PTea, PTi0, PTia, V0, rGamm0, rGASPF, ZEFF0
    integer(4) :: IFNT

    call INQFNT(IFNT)
!    call SETFNT(32)
    call SETFNT(0)
    call SETCHS(0.3, 0.0)
    call SETLIN(0, 1, 7)

    GXM = 13.0 + 1.0
    GYM =  8.5 + 0.2
    GYS =  0.5
    NP  = 0

    call TXWPS('+'//SLID//'+')
    call TXWPS('@PNBCD @', PNBCD)
!!!    call TXWPS('@NRMAX @', NRMAX)

!    NP = NP + 1
    call TXWPS('@BB    @', BB)
    call TXWPS('@rIp   @', rIp)
    call TXWPS('@FSDFX1@', FSDFIX(1))
    call TXWPS('@FANOM3@', FSANOM(3))
    call TXWPS('@FSBOHM@', FSBOHM)
    call TXWPS('@FSPCL1@', FSPCL(1))
    call TXWPS('@FSPCL2@', FSPCL(2))
    call TXWPS('@FSPCL3@', FSPCL(3))
    call TXWPS('@PROFD @', PROFD)
    call TXWPS('@PROFC @', PROFC)
    call TXWPS('@FSCX  @', FSCX)
    call TXWPS('@FSRP  @', FSRP)
    call TXWPS('@FSLC  @', FSLC)
    call TXWPS('@FSNC1 @', FSNC(1))
    call TXWPS('@FSNC2 @', FSNC(2))
    call TXWPS('@FSLP  @', FSLP)
    call TXWPS('@FSLPB @', FSLPB)
    call TXWPS('@FSION @', FSION)
    call TXWPS('@FSD02 @', FSD02)

    GXM = GXM + 0.35 * 17
    NP = 0
!!!    call TXWPS('@PNBH  @', PNBH)
    call TXWPS('@PNBHT @', PNBHT1+PNBHT2)
    if(PNBHex == 0.D0) then
       call TXWPS('@PNBHP @', PNBHP)
    else
       call TXWPS('@PNBHex@', PNBHex)
    end if
    call TXWPS('@PRFH  @', PRFHe+PRFHi)
    call TXWPS('@Vb    @', Vb)
    call TXWPS('@De0   @', Dfs0(1))
!!!    call TXWPS('@Di0   @', Dfs0(2))
    call TXWPS('@rMue0 @', rMus0(1))
    call TXWPS('@rMui0 @', rMus0(2))
    call TXWPS('@Chie0 @', Chis0(1))
    call TXWPS('@Chii0 @', Chis0(2))
!!!    call TXWPS('@WPM0  @', WPM0)
    call TXWPS('@PTe0  @', PTe0)
    call TXWPS('@PTea  @', PTea)
    call TXWPS('@PTi0  @', PTi0)
    call TXWPS('@PTia  @', PTia)
    call TXWPS('@V0    @', V0)
    call TXWPS('@rGamm0@', rGamm0)
    call TXWPS('@rGASPF@', rGASPF)
    call TXWPS('@ZEFF0 @', ZEFF0)

    call SETFNT(IFNT)

  end subroutine TXWPGR

  !***************************************************************
  !
  !   WPgr's Sub : write Double
  !
  !***************************************************************

  subroutine TXWPSD(STR, VAL)

    real(8), INTENT(IN) :: VAL
    character(len=*), INTENT(IN) :: STR

    call MOVE(GXM, GYM - GYS * NP)
    call TEXTX(STR)
    call NUMBD(VAL,'(ES9.2)', 9)
    NP = NP + 1

  end subroutine TXWPSD

  !***************************************************************
  !
  !   WPgr's Sub : write Integer
  !
  !***************************************************************

  subroutine TXWPSI(STR, IVAL)

    integer(4), INTENT(IN) :: IVAL
    character(len=*), INTENT(IN) :: STR

    call MOVE(GXM, GYM - GYS * NP)
    call TEXTX(STR)
    call TEXT(' = ', 3)
    call NUMBI(IVAL,'(I6)',6)
    NP = NP + 1

  end subroutine TXWPSI

  !***************************************************************
  !
  !   WPgr's Sub : write Strings
  !
  !***************************************************************

  subroutine TXWPSS(STR)

    character(len=*), INTENT(IN) :: STR

    call MOVE(GXM, GYM - GYS * NP)
    call TEXTX(STR)
    NP = NP + 1

  end subroutine TXWPSS

  !***************************************************************
  !
  !   Draw graph
  !
  !***************************************************************

  subroutine TXGRAF(GPXY, GX, GY, NXM, NXMAX, NGMAX, &
       &                  GXMIN, GXMAX, STR, FONT, MODE, IND, ILOG, &
       &                  GYMAX_IN, GYMIN_IN)!, Kindex)

    integer(4), INTENT(IN) :: NXM, NXMAX, NGMAX, MODE, IND
    real(4), INTENT(IN) :: GXMIN, GXMAX, FONT
    real(4), DIMENSION(4), INTENT(IN) :: GPXY
    real(4), DIMENSION(:), INTENT(IN) :: GX
    real(4), DIMENSION(:,:), INTENT(IN) :: GY
    integer(4), intent(in) :: ILOG
    real(4), INTENT(IN), OPTIONAL :: GYMAX_IN, GYMIN_IN
    character(len=*), INTENT(IN) :: STR
!    character(len=*), INTENT(IN), optional :: Kindex

    integer(4) :: IFNT, NGV, NGULEN, ICL, IPAT, IMRK, ISTEP, NG
    real(4) :: GX1, GX2, GY1, GY2, gSLEN, GSXMIN, GSXMAX, GXSTEP, &
         &        GYMIN, GYMAX, GSYMIN, GSYMAX, GYSTEP, GYORG,  &
         &        GMRK, GCHH, GXL, GYL
    integer(4), DIMENSION(0:4) :: NLTYPE
    real(4), dimension(:,:), allocatable :: GYAR
    DATA NLTYPE/0,2,3,4,6/

    if (MODE < 0 .or. MODE > 4) then
       write(6,*) '### ERROR(TXGRAF) : MODE = ', MODE
       return
    end if

    GX1=GPXY(1)
    GX2=GPXY(2)
    GY1=GPXY(3)
    GY2=GPXY(4)

    call INQFNT(IFNT)

    if (IND == 0) then
       gSLEN = 1.0
    else
       gSLEN = 0.2
    end if

    call SETFNT(32)
    call SETCHS(FONT, 0.0)
    call SETLIN(0, 0, 7)
    call GTEXTX(GX1,GY2+0.2,STR,0)

    call SETFNT(44)
    call SETCHS(FONT, 0.0)
    call SETLIN(0, 0, 7)
    call SETLNW(0.035)

    allocate(GYAR(size(GY,1),size(GY,2)))
    GYAR = GY
!    if(present(GYMAX_IN)) where(GYAR > GYMAX_IN) GYAR = GYMAX_IN
!    if(present(GYMIN_IN)) where(GYAR < GYMIN_IN) GYAR = GYMIN_IN

    call GQSCAL(GXMIN, GXMAX, GSXMIN, GSXMAX, GXSTEP)
    GSXMIN = GXMIN
    GSXMAX = GXMAX
    call GMNMX2(GYAR,NXM,1,NXMAX,1,1,NGMAX,1,GYMIN,GYMAX)
    if(ILOG == 0) then
       if (GYMAX > 0.0) then
          if (GYMIN > 0.0) GYMIN=0.0
       else
          GYMAX=0.0
       end if
    end if
!    if(present(Kindex)) then
       if(present(GYMAX_IN)) GYMAX=GYMAX_IN
       if(present(GYMIN_IN)) GYMIN=GYMIN_IN
!    end if
    call GQSCAL(GYMIN, GYMAX, GSYMIN, GSYMAX, GYSTEP)
    if (GSYMIN > GYMIN) GSYMIN = GSYMIN - GYSTEP
    if (GSYMAX < GYMAX) GSYMAX = GSYMAX + GYSTEP
    if(GYMIN * GYMAX <= 0.0) then
       GYORG = 0.0
    else
       GYORG = GSYMIN
    endIF

    call GDEFIN(GX1, GX2, GY1, GY2, GSXMIN, GSXMAX, GSYMIN, GSYMAX)
    call GFRAME
    call SETLNW( 0.017)
    call GSCALE(GSXMIN, GXSTEP, 0.0, 0.0, gSLEN, IND)
!!$    if (GXSTEP < 0.01) then
!!$       NGV=NGULEN(GXSTEP*5)
!!$       if (NGV < 0) then
!!$          NGV=NGV-3200
!!$       else
!!$          NGV=NGV+3200
!!$       end if
!!$       call GVALUE(GSXMIN, GXSTEP*5, 0.0, 0.0, NGV)
       NGV=NGULEN(GXSTEP*2)
       call GVALUE(GSXMIN, GXSTEP*2, 0.0, 0.0, NGV)
!!$    else
!!$       NGV=NGULEN(GXSTEP*2)
!!$       if (NGV < 0) then
!!$          NGV=NGV-3200
!!$       else
!!$          NGV=NGV+3200
!!$       end if
!!$       call GVALUE(GSXMIN, GXSTEP*2, 0.0, 0.0, NGV)
!!$    end if

    ! Semi-Log or not

    if(ILOG == 0) then
       call GSCALE(0.0, 0.0, GYORG, GYSTEP, gSLEN, IND)
       call SETLNW(-0.017)
       if (GSYMIN < 0.0 .and. GSYMAX > 0.0) &
            &     call GSCALE(0.0, 0.0, 0.0, GSYMAX-GSYMIN,  2.0, 0)
       call GVALUE(0.0,0.0,GYORG,GYSTEP*2,NGULEN(GYSTEP*2))
    else
       call GScall(0.0, 0.0, GYORG, 4, gSLEN, IND)
       call SETLNW(-0.017)
       if (GSYMIN < 0.0 .and. GSYMAX > 0.0) &
            &     call GScall(0.0, 0.0, 0.0, 3,  2.0, 0)
       call GVALUL(0.0,0.0,GYORG,1,NGULEN(GYSTEP*2))
    end if

    ! MODE = 0: Change Line Color (Last Color Fixed)

    select case(MODE)
    case (0)
       do NG = 1, NGMAX
          ICL = 7 - mod(NGMAX - NG, 5)
          call SETLIN(0, 1, ICL)
          call GPLOTP(GX, GYAR(1,NG), 1, NXMAX, 1, 0, 0, 0)
       end do

       ! MODE = 1: Change Line Color and Style
       ! MODE = 2: Change Line Color and Style (With Legend)

    case (1:2)
       if (MODE == 1) then
          do NG = 1, NGMAX
             ICL  = 7 - mod(NG-1, 5)
             IPAT = NLTYPE(mod(NG-1, 5))
             call SETLIN(0, 1, ICL)
             call GPLOTP(GX, GYAR(1,NG), 1, NXMAX, 1, 0, 0, IPAT)
          end do
       else
          do NG = 1, NGMAX
             ICL = 7 - mod(NG - 1, 5)
             ISTEP = NXMAX / 10
             IPAT  = (NG - 1) / 5
             call SETLIN(0, 1, ICL)
             call GPLOTP(GX, GYAR(1,NG), 1, NXMAX, 1, 0, ISTEP, IPAT)
          end do
       end if
       ! Legend
       if (MODE == 2) then
          call SETCHS(0.25,0.0)
          GCHH = 0.25
          GXL = GX2 + GCHH
          GYL = GY2 - GCHH
          do NG = 1, NGMAX
             ICL = 7 - mod(NG - 1, 5)
             call SETLIN(0, 1, ICL)
             IPAT = (NG - 1) / 5
             call MOVEPT(GXL + GCHH * 3.0, GYL + GCHH / 2.0, IPAT)
             call DRAWPT(GXL + GCHH * 7.0, GYL + GCHH / 2.0)
             call MOVE(GXL, GYL)
             call NUMBI(NG, '(I2)', 2)
             GYL = GYL - GCHH * 2.0
          end do
       end if
       ! End of Legend

       ! MODE = 3: Change Line Color, Style and Mark
       ! MODE = 4: Change Line Color, Style and Mark (With Legend)

    case (3:4)
       IMRK = 0
       GMRK = 0.3
       call SETMKS(IMRK, GMRK)
       do NG = 1, NGMAX
          ICL = 7 - mod(NG - 1, 5)
          IMRK  = mod(NG - 1, 5) + 1
          ISTEP = NXMAX / 10
          IPAT  = (NG - 1) / 5
          call SETLIN(0, 1, ICL)
          call GPLOTP(GX, GYAR(1,NG), 1, NXMAX, 1, -IMRK, ISTEP, IPAT)
       end do
       ! Legend
       if (MODE == 4) then
          call SETCHS(0.25,0.0)
          GCHH = 0.25
          GXL = GX2 + GCHH
          GYL = GY2 - GCHH
          do NG = 1, NGMAX
             ICL = 7 - mod(NG - 1, 5)
             call SETLIN(0, 1, ICL)
             IPAT = (NG - 1) / 5
             call MOVEPT(GXL + GCHH * 3.0, GYL + GCHH / 2.0, IPAT)
             call DRAWPT(GXL + GCHH * 7.0, GYL + GCHH / 2.0)
             call MOVE(GXL, GYL)
             call NUMBI(NG, '(I2)', 2)
             IMRK = mod(NG - 1, 5) + 1
             if (IMRK /= 0) then
                call SETMKS(IMRK, GMRK)
                call MARK(GXL + GCHH * 5.0, GYL + GCHH / 2.0)
             end if
             GYL = GYL - GCHH * 2.0
          end do
       end if
       ! End of Legend
       IMRK = 0
       GMRK = 0.2
       call SETMKS(IMRK, GMRK)
    end select

    call SETFNT(IFNT)

    deallocate(GYAR)

  end subroutine TXGRAF

  !***************************************************************
  !
  !   Appropriate GY and STR for graphic routine
  !
  !     Note: NMAX is the array size of GOUT, i.e., NMAX=size(GOUT)
  !
  !***************************************************************

  subroutine APPROPGY(MODE, GIN, GOUT, NMAX, STR, gDIV, GMAX, GMIN)
    use tx_interface, only : APTOS

    integer(4), intent(in) :: MODE, NMAX
    real(4), intent(in) :: gDIV
    real(4), intent(out), optional :: GMAX, GMIN
    real(4), dimension(NMAX),intent(in)  :: GIN
    real(4), dimension(NMAX),intent(out) :: GOUT
    character(len=*), intent(inout) :: STR

    integer(4) :: NSTR, POSAT
    real(4) :: gDIVL

    gDIVL = 1.0
    if(MODE /= 0) gDIVL = gDIV

    ! Append gDIV to string for showing multiplication factor
    if (gDIVL /= 1.0) then
       POSAT = index(STR,'@',.true.)
       if (POSAT /= 0) STR = STR(1:POSAT-1)
       NSTR = len_trim(STR)+1
       call APTOS(STR, NSTR, ' [', 2)
       call APTOS(STR, NSTR, gDIVL, 'E0')
       call APTOS(STR, NSTR, ']@', 2)
    end if

    GOUT(:) = GIN(:) / gDIVL

    if(present(GMAX)) GMAX = maxval(GOUT)
    if(present(GMIN)) GMIN = minval(GOUT)

  end subroutine APPROPGY

  !***********************************************************
  !
  !   Ceiling function for LOG10 plot
  !
  !***********************************************************

  real(4) function GLOG(X,XMIN,XMAX)

    use libplog, only : plog
    implicit none
    real(4) :: GUCLIP
    real(8), intent(in) :: X, XMIN, XMAX

    GLOG = GUCLIP(PLOG(X,XMIN,XMAX))

  end function GLOG
#endif
  !***********************************************************
  !
  !   write RAW VALUE ONTO CONSOLE
  !
  !***********************************************************

  subroutine write_console(n,char)

    use tx_commons, only : NRMAX, rpt
    use tx_interface, only : dfdx
    integer(4), intent(in) :: n
    character(len=1), intent(in) :: char
    character(len=12) :: cfmt
    integer(4) :: i, j, iout = 111
    real(8), dimension(:), allocatable :: dVdr

    if (ngr <= -1) then
       write(6,*) 'No data.'
       return
    end if

    select case(char)
    case('R')
       if (n >= 0 .and. n <= NGYRM) then
          write(6,'(A1,6X,A3,11X,A5)') "#","rho","Value"
          do i = 0, nrmax
             write(6,*) gx(i),gy%v(i,ngr,n)
          end do
       end if
    case('T')
       if (n >= 0 .and. n <= NGYTM) then
          write(6,'(A1,6X,A1,13X,A5)') "#","T","Value"
          do i = 0, ngt
             write(6,*) gtx(i),gty(i,n)
          end do
       end if
    case('V')
       if (n >= 0 .and. n <= NGYVM) then
          write(6,'(A1,6X,A1,13X,A5)') "#","T","Value"
          do i = 0, ngvv
             write(6,*) gvx(i),gvy(i,n)
          end do
       end if
    case('D') ! R-derivative of a radial profile
       if (n >= 0 .and. n <= NGYRM) then
          allocate(dVdr,mold=rpt)
          dVdr(0:nrmax) = dfdx(rpt,real(gy%v(0:nrmax,ngr,n),8),nrmax,0)
          write(6,'(A1,6X,A3,11X,A5)') "#","rho","dValue/dR"
          do i = 0, nrmax
             write(6,*) gx(i),real(dVdr(i))
          end do
          deallocate(dVdr)
       end if
    case('A') ! All radial profiles
       cfmt = '(1PxxxE15.7)'
       write(cfmt(4:6),'(I3)') NGYRM+1

       do i = 0, NRMAX
          write(iout,cfmt) GX(i), (GY%v(i,NGR,j), j = 1, NGYRM)
       end do
    case default
       write(6,*) '### ERROR : Invalid Command : '
    end select

  end subroutine write_console

  !***********************************************************
  !
  !   Return NGT corresponding to T
  !
  !***********************************************************

  subroutine return_NGT_from_T(NGT,T,IERR)

    real(8), intent(in) :: T
    integer(4), intent(out) :: NGT, IERR
    integer(4) :: n

    IERR = 0
    do n = 0, NGTM
       if(abs(GTX(n) - real(T)) < epsilon(1.d0)) then
          NGT = n
          return
       end if
    end do

    IERR = 1

  end subroutine return_NGT_from_T

  !***********************************************************
  !
  !   Psi contour output for gnuplot
  !
  !***********************************************************

  subroutine psi_out_gnuplot

    use equ_params, only : psi, rg, zg, nsr, nsz
    integer(4) :: i, j, ist, iopsi = 21, iogpl = 22
    character(len=9)  :: file_psi = 'psi_out.d'
    character(len=11) :: file_gpl = 'psi_out.gpl' &
         &              ,file_out = 'psi_out.eps'

    ! Create psi contour data

    open(iopsi,file=file_psi,iostat=ist,status='replace',form='formatted')
    do j = 1, nsz
       do i = 1, nsr
          write(iopsi,*) rg(i), zg(j), psi(i+(j-1)*nsr)
       end do
       write(iopsi,*) 
    end do

    ! Create gnuplot script

    open(iogpl,file=file_gpl,iostat=ist,status='replace',form='formatted')
    write(iogpl,*) 'set term postscript enh color "Helvetica" 16'
    write(iogpl,*) 'set output "',file_out,'"'
    write(iogpl,*) 'unset surface'
    write(iogpl,*) 'unset key'
    write(iogpl,*) 'unset ztics'
    write(iogpl,*) 'set ytics offset 1,0'
    write(iogpl,*) 'set view 0,0'
    write(iogpl,*) 'set view equal xy' ! N.B. 'set size ratio -1' for 2D
    write(iogpl,*) 'set style textbox opaque noborder'
    write(iogpl,*) 'set contour base'
    write(iogpl,*) 'set cntrparam levels auto 20'
    write(iogpl,*) 'set cntrlabel onecolor'
    write(iogpl,*) 'set cntrlabel format "%8.3g" font ",7" start 5 interval 100'
!!$    write(iogpl,*) 'splot "',file_psi,'" using 1:2:3 t "psi" w l,\'
!!$    write(iogpl,*) '"" using 1:2:3 w labels boxed'
    write(iogpl,*) 'splot "',file_psi,'" using 1:2:3 t "psi" w l'
    write(iogpl,*) 'exit'

    ! Call gnuplot

    call system("gnuplot "//file_gpl)
    call system("gv "//file_out)

    close(iopsi,status='delete')
    close(iogpl,status='delete')

  end subroutine psi_out_gnuplot

end module tx_graphic
