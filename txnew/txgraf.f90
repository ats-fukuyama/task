!     $Id: txgraf.f90,v 1.83 2011/01/20 07:06:17 honda Exp $
module tx_graphic
  implicit none
  private
  integer(4), parameter :: NGRM=20, NGTM=5000, NGVM=5000, &
       &                   NGYRM=200, NGYTM=53, NGYVM=58, &
       &                   NGPRM=24, NGPTM=8, NGPVM=15
  real(4) :: GXM, GYM, GYS
  integer(4) :: MODEG, MODEGL, NGR, NGT, NGVV, NGRSTP, NGTSTP, NGVSTP, NP
  real(4), dimension(:),     allocatable :: GX
  real(4), dimension(:,:,:), allocatable :: GY, GQY, GYT
  real(4), dimension(0:NGRM) :: GT
  real(4), dimension(0:NGTM) :: GTX
  real(4), dimension(0:NGVM) :: GVX
  real(4), dimension(1:NGYRM) :: gDIV
  real(4), dimension(0:NGTM,1:NGYTM) :: GTY
  real(4), dimension(0:NGVM,1:NGYVM) :: GVY
  public :: TXGOUT, TX_GRAPH_SAVE, TXSTGR, TXSTGT, TXSTGV, TXSTGQ, allocate_txgraf, deallocate_txgraf, &
       &    MODEG, MODEGL, NGR, NGT, NGVV, NGYRM, NGYTM, NGYVM, NGRM, NGRSTP, NGTSTP, NGVSTP, &
       &    GY, GQY, GYT, GT, GTX, GVX, gDIV, GTY, GVY

  interface TXWPS
     module procedure TXWPSD
     module procedure TXWPSI
     module procedure TXWPSS
  end interface

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
    integer(4), dimension(1:2) :: ierl

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

    allocate(GX(0:N), GY(0:N,0:NGRM,1:NGYRM), GQY(0:N,1:NCM,1:NQMAX),      stat = ierl(1))
    allocate(GYT(0:N,0:NGTM,1:NGYRM),                                      stat = ierl(2))
    ier = sum(ierl)

    ! All the memories allocated above are clear if some errors occur.
    if(ier /= 0) then
       write(6,*) "XX Allocation error in allocate_txgraf"
       call deallocate_txgraf
    end if

  end subroutine allocate_txgraf

  subroutine deallocate_txgraf

    deallocate(GX, GY, GQY)
    deallocate(GYT)

  end subroutine deallocate_txgraf
  
  !***************************************************************
  !
  !   GRAPHIC Command loop
  !
  !***************************************************************
  
  SUBROUTINE TXGOUT
    use libbes, only : BESIN
    use libgrf, only : GRD1D
    use tx_commons, only : T_TX, TPRE, NQMAX, NRMAX, rhob, RA, PI, RR, NRA, thrp, kappa, &
         &                 IRPIN, DltRPn, NTCOIL, DltRP_mid, DltRP, epst, rho
    use tx_interface, only : TXGRUR, TOUPPER, TXGLOD
    use tx_ripple, only : ripple

    INTEGER(4) :: MODE, NGPR, NGPT, NGPV, NGYR, NQ, NQL, NGF, NGFMAX, I, IST, NGRT, NG, IER, J, NGTL
!    real(4), dimension(0:NRMAX,0:5,1:NGYRM) :: GYL
    real(4), dimension(:,:,:), allocatable :: GYL
    character(len=5) :: STR, STR2
    character(len=1) :: KID1, KID2

    integer(4), parameter :: NXMAX = 51, NYMAX = 51, NTHMAX=51
    integer(4) :: NX, NY, NR, NTH, NGTDO, IPAT, IPRD, NSTEP, IND, IFNT
    integer(4), dimension(:,:,:), allocatable :: KA
    integer(4), dimension(:), allocatable :: NGPRL
    real(4) :: DR, DZ, AL, ZORG, ZSTEP, DltRPnL, GXMAX, GMAX, GMIN, GZ
    real(4), dimension(1:4) :: GPXY
    real(4), dimension(:), allocatable :: RRL, ZZL, GPXY_IN, GMAXA, GMINA
    real(4), dimension(:,:), allocatable :: VAL, GRPL, GYL2
    real(8), dimension(:), allocatable :: FX
    real(8), dimension(:,:), allocatable :: FY
    real(4) :: GSRMIN,GSRMAX,GSRSTP,GSZMIN,GSZMAX,GSZSTP
    real(8) :: rhol, theta, dtheta, kappal, thetab, TL
    character(len=50) :: STRL
    character(len=17) :: KOUT
    character(len=50), dimension(:), allocatable :: STRA
    integer(4) :: NGULEN

    !     *** MENU ***

    MODE = MODEGL
    OUTER : DO
       WRITE(6,*) '# SELECT : Rn: Tn: Un: Vn: A,B: C: E: Dn: Fn: Gn: Hn: S,L,M:file'
       WRITE(6,*) '           W(R,T,V)n:write I:init X:exit'
       READ(5,'(A5)',IOSTAT=IST) STR
       IF (IST > 0) THEN
          WRITE(6,*) '### ERROR : Invalid Command : ', STR
          CYCLE
       ELSEIF (IST < 0) THEN
          RETURN
       END IF

       KID1 = STR(1:1)
       KID2 = STR(2:2)
       CALL TOUPPER(KID1)
       CALL TOUPPER(KID2)

       SELECT CASE(KID1)
       CASE('S')
          CALL TXGSAV

       CASE('L')
          CALL TXGLOD(IER)
          IF(IER /= 0) CYCLE OUTER
          CALL TXPRFG

       CASE('I')
          ! *** Initialization ***
          T_TX = 0.D0
          TPRE = 0.D0
          NGT = -1
          NGVV = -1
          CALL TXSTGT(REAL(T_TX))
          CALL TXSTGV(REAL(T_TX))

       CASE('R')
          ! *** Time evolution of radial profiles ***
          SELECT CASE(KID2)
          CASE('A')
             CALL TXGRFR(-1,MODE)
          CASE('B')
             CALL TXGRFR(-3,MODE)
          CASE('C')
             CALL TXGRFR(-5,MODE)
          CASE('D')
             CALL TXGRFR(-9,MODE)
          CASE DEFAULT
             READ(STR(2:5),*,IOSTAT=IST) NGPR
             IF (IST /= 0) THEN
                WRITE(6,*) '### ERROR : Invalid Command : ', STR
                CYCLE OUTER
             END IF
             IF (NGPR >= 0 .AND. NGPR <= NGPRM) THEN
                CALL TXGRFR(NGPR,MODE)
             END IF
          END SELECT

       CASE('V')
          ! *** Time evolution of variables at a certain position ***
          SELECT CASE(KID2)
          CASE('A')
             DO NGPV = 1, NGPVM
                CALL TXGRFV(NGPV,MODE)
             END DO
          CASE('B')
             DO NGPV = 1, 6
                CALL TXGRFV(NGPV,MODE)
             END DO
             DO NGPV = 9, 10
                CALL TXGRFV(NGPV,MODE)
             END DO

          CASE DEFAULT
             READ(STR(2:5),*,IOSTAT=IST) NGPV
             IF (IST < 0) THEN
                WRITE(6,*) '### ERROR : Invalid Command : ', STR
                CYCLE OUTER
             END IF
             IF (NGPV >= 1 .AND. NGPV <= NGPVM) THEN
                CALL TXGRFV(NGPV,MODE)
             END IF
          END SELECT

       CASE('T')
          ! *** Time evolution of global plasma quantities ***
          SELECT CASE(KID2)
          CASE('A')
             DO NGPT = 1, NGPTM
                CALL TXGRFT(NGPT,MODE)
             END DO

          CASE('B')
             DO NGPT = 1, 6
                CALL TXGRFT(NGPT,MODE)
             END DO
             DO NGPT = 9, 10
                CALL TXGRFT(NGPT,MODE)
             END DO

          CASE DEFAULT
             READ(STR(2:5),*,IOSTAT=IST) NGPT
             IF (IST < 0) THEN
                WRITE(6,*) '### ERROR : Invalid Command : ', STR
                CYCLE OUTER
             END IF
             IF (NGPT >= 1 .AND. NGPT <= NGPTM) THEN
                CALL TXGRFT(NGPT,MODE)
             END IF
          END SELECT

       CASE('U')
          ! *** Radial profiles of the balance among the terms in each equation ***
          IF(T_TX == 0.D0) CYCLE OUTER
          IF (KID2 == 'A') THEN
             DO NQ = 1, NQMAX, 4
                CALL PAGES
                DO NQL=NQ,MIN(NQ+3,NQMAX)
                   CALL TXGRFQ(NQL,MOD(NQL-1,4)+1)
                END DO
                CALL PAGEE
             END DO

          ELSE
             READ(STR(2:5),*,IOSTAT=IST) NQ
             IF (IST < 0) THEN
                WRITE(6,*) '### ERROR : Invalid Command : ', STR
                CYCLE OUTER
             END IF
             IF (NQ >= 1 .AND. NQ <= NQMAX) THEN
                CALL PAGES
!!$              DO NQL=NQ,MIN(NQ+3,NQMAX)
!!$                 CALL TXGRFQ(NQL,NQL-NQ+1)
!!$              END DO
                CALL TXGRFQ(NQ,5)
                CALL PAGEE
             END IF
          END IF

       CASE('A')
          ! *** Sequential display of radial profiles (version A) ***
          DO NGPR = 1, NGPRM
             CALL TXGRFR(NGPR,MODE)
          END DO

       CASE('B')
          ! *** Sequential display of radial profiles (version B) ***
          DO NGPR = 1, 6
             CALL TXGRFR(NGPR,MODE)
          END DO
          DO NGPR = 9, 10
             CALL TXGRFR(NGPR,MODE)
          END DO

       CASE('C')
          ! *** Comparison of neoclassical characteristics ***
          IF (KID2 == 'A') THEN
             CALL TXGRCPA(MODE)
          ELSE
             CALL TXGRCP(MODE)
          END IF

       CASE('D')
          ! *** Time evolution of radial profiles in 3D phase space ***
          READ(STR(2:5),*,IOSTAT=IST) NGYR
          IF (IST /= 0) THEN
             WRITE(6,*) '### ERROR : Invalid Command : ', STR
             CYCLE OUTER
          END IF
          IF (NGYR >= 0 .AND. NGYR <= NGYRM) THEN
             CALL TXGRUR(GX,GTX,GYT(0:NRMAX,0:NGT,NGYR),NRMAX,NGT,NGTM)!,STR,KV,MODE)
          END IF

       CASE('F')
          ! *** Radial profile of a certain variable ***
          READ(STR(2:5),*,IOSTAT=IST) NGYR
          IF (IST /= 0) THEN
             WRITE(6,*) '### ERROR : Invalid Command : ', STR
             CYCLE OUTER
          END IF
          IF (NGYR >= 1 .AND. NGYR <= NGYRM) THEN
             ALLOCATE(FX(0:NRMAX),FY(0:NRMAX,1))
             DO NR=0,NRMAX
                FX(NR)=DBLE(GX(NR))
                FY(NR,1)=DBLE(GYT(NR,NGT,NGYR))
             ENDDO
             STRL='@profile'//STR(2:5)//'@'
             CALL PAGES
             CALL GRD1D(0,FX,FY,NRMAX+1,NRMAX+1,1,STRL,0)
             CALL PAGEE
             DEALLOCATE(FX,FY)
          END IF

       CASE('G')
          ! *** Time evolution of radial profile of a certain variable ***
          READ(STR(2:5),*,IOSTAT=IST) NGYR
          IF (IST /= 0) THEN
             WRITE(6,*) '### ERROR : Invalid Command : ', STR
             CYCLE OUTER
          END IF
          IF (NGYR >= 1 .AND. NGYR <= NGYRM) THEN
             ALLOCATE(FX(0:NRMAX),FY(0:NRMAX,0:NGT))
             DO NR=0,NRMAX
                FX(NR)=DBLE(GX(NR))
             ENDDO
             DO NGTDO=0,NGT
             DO NR=0,NRMAX
                FY(NR,NGTDO)=DBLE(GYT(NR,NGTDO,NGYR))
             ENDDO
             ENDDO
             STRL='@profile'//STR(2:5)//'@'
             CALL PAGES
             CALL GRD1D(0,FX,FY,NRMAX+1,NRMAX+1,NGT+1,STRL,0)
             CALL PAGEE
             DEALLOCATE(FX,FY)
          END IF

       CASE('E')
          ! *** Contour of ripple amplitude ***
          CALL PAGES
          CALL INQFNT(IFNT)
          CALL SETFNT(32)

          if(IRPIN == 0) then
             DltRPnL = REAL(DltRPn) * 100.0

             allocate(RRL(1:NXMAX),ZZL(1:NYMAX),VAL(1:NXMAX,1:NYMAX),KA(1:8,1:NXMAX,1:NYMAX))
             DR = (4.5 - 2.0) / (NXMAX - 1)
             DZ =  1.5        / (NYMAX - 1)
             DO NX = 1, NXMAX
                RRL(NX) = 2.0 + DR * (NX - 1)
                DO NY = 1, NYMAX
                   ZZL(NY) = DZ * (NY - 1)
                   AL = SQRT(((RRL(NX) - 2.4)**2 + ZZL(NY)**2) * (2.4 / RRL(NX)))
                   VAL(NX,NY) = REAL(DltRPnL * BESIN(0,NTCOIL/2.4D0*AL))
                END DO
             END DO

             CALL GQSCAL(2.0,4.5,GSRMIN,GSRMAX,GSRSTP)
             CALL GQSCAL(0.0,1.5,GSZMIN,GSZMAX,GSZSTP)

             CALL GDEFIN(3.0,18.0,1.5,10.5,2.0,4.5,0.0,1.5)
             CALL GFRAME

             CALL GSCALE(GSRMIN,GSRSTP,0.0,0.0,0.1,9)
             CALL GVALUE(GSRMIN,GSRSTP*2,0.0,0.0,NGULEN(GSRSTP))
             CALL GSCALE(0.0,0.0,0.0,GSZSTP,0.1,9)
             CALL GVALUE(0.0,0.0,0.0,GSZSTP*2,NGULEN(GSZSTP*2))

             IPAT = 1
             IPRD = 0

             ZORG  = 0.0003
             ZSTEP = 0.0002
             NSTEP = 3
             CALL CONTP2(VAL,RRL,ZZL,NXMAX,NXMAX,NYMAX,ZORG,ZSTEP,NSTEP,IPRD,IPAT*1,KA)

             ZORG  = 0.001
             ZSTEP = 0.002
             NSTEP = 4
             CALL CONTP2(VAL,RRL,ZZL,NXMAX,NXMAX,NYMAX,ZORG,ZSTEP,NSTEP,IPRD,IPAT*2,KA)

             ZORG  = 0.01
             ZSTEP = 0.02
             NSTEP = 4
             CALL CONTP2(VAL,RRL,ZZL,NXMAX,NXMAX,NYMAX,ZORG,ZSTEP,NSTEP,IPRD,IPAT*3,KA)

             ZORG  = 0.0
             ZSTEP = 0.1
             NSTEP = 2
             CALL CONTP2(VAL,RRL,ZZL,NXMAX,NXMAX,NYMAX,ZORG,ZSTEP,NSTEP,IPRD,IPAT*4,KA)

             ZORG  = 0.0
             ZSTEP = 0.25
             NSTEP = 4
             CALL CONTP2(VAL,RRL,ZZL,NXMAX,NXMAX,NYMAX,ZORG,ZSTEP,NSTEP,IPRD,IPAT*4,KA)

             ZORG  = 1.0
             ZSTEP = 0.5
             NSTEP = 3
             CALL CONTP2(VAL,RRL,ZZL,NXMAX,NXMAX,NYMAX,ZORG,ZSTEP,NSTEP,IPRD,IPAT*5,KA)

             deallocate(RRL,ZZL,VAL,KA)

             allocate(RRL(1:NTHMAX),ZZL(1:NTHMAX))
             dtheta = PI / (NTHMAX - 1)
             theta = 0.d0
             DO I = 1, 2
                IF(I == 1) THEN
                   kappal = 1.d0
                ELSE
                   kappal = kappa
                END IF
                DO NTH = 1, NTHMAX
                   theta = (NTH - 1) * dtheta
                   RRL(NTH) = real(RR + RA * cos(theta))
                   ZZL(NTH) = real(kappal * RA * sin(theta))
                   IF(ABS(ZZL(NTH)) < EPSILON(1.0)) ZZL(NTH) = 0.0
                END DO
                CALL LINES2D(RRL,ZZL,NTHMAX)
             END DO
             deallocate(RRL,ZZL)

             allocate(RRL(1:2*NRMAX),ZZL(1:2*NRMAX))
             RRL(1:NRMAX) = RR * (1.d0 + epst(NRMAX:1:-1) * COS(thrp(1:NRMAX)))
             ZZL(1:NRMAX) = kappa * RR * epst(NRMAX:1:-1) * SIN(thrp(1:NRMAX))
             RRL(NRMAX+1:2*NRMAX) = RR * (1.d0 + epst(1:NRMAX) * COS(thrp(NRMAX+1:2*NRMAX)))
             ZZL(NRMAX+1:2*NRMAX) = kappa * RR * epst(1:NRMAX) * SIN(thrp(NRMAX+1:2*NRMAX))
             CALL LINES2D(RRL,ZZL,2*NRMAX)
             deallocate(RRL,ZZL)
          end if

          !  *** 1D graphic ***

          IF (MODEG == 2) THEN
             IND = 9
          ELSE
             IND = 0
          END IF

          GPXY(1) =  3.0
          GPXY(2) = 11.5
          GPXY(3) = 11.5
          GPXY(4) = 17.5
          GXMAX = REAL(rhob)
          allocate(GRPL(0:NRMAX,1:2))
          DO NR = 0, NRMAX
             if(IRPIN == 0) then
                GRPL(NR,1) = real(ripple(rho(NR),0.D0,1.D0))
             else
                GRPL(NR,1) = real(DltRP_mid(NR))
             end if
          END DO
          WRITE(KOUT,'(F6.4)') GRPL(NRA,1)
          KOUT = '$#d$#$-a$=='//KOUT(1:6)
          CALL GTEXT(GPXY(1)+5.5,GPXY(4)-1.5,KOUT,len_trim(KOUT),0)
          STRL = '@DltRP at mid-plane(r)@'
          CALL TXGRAF(GPXY,GX,GRPL,NRMAX+1,NRMAX+1,1,0.0,GXMAX,STRL,0.26,MODE,IND,0)

          GPXY(1) = 13.5
          GPXY(2) = 22.0
          GPXY(3) = 11.5
          GPXY(4) = 17.5
          GXMAX = REAL(rhob)
          thetab = 0.5D0 * PI
          DO NR = 0, NRMAX
             rhol = rho(NR) * (1.D0 + (kappa - 1.D0) * sin(thetab))
             if(IRPIN == 0) then
                GRPL(NR,1) = real(ripple(rhol,thetab,1.D0))
                GRPL(NR,2) = real(ripple(rho(NR),thetab,1.D0))
                NG = 2
             else
                GRPL(NR,1) = real(DltRP(NR))
                NG = 1
             end if
          END DO
          if(IRPIN == 0) then
             STRL = '@DltRP at tip point with and w/o elongation(r)@'
          else
             STRL = '@DltRP at tip point(r)@'
          end if
          CALL TXGRAF(GPXY,GX,GRPL,NRMAX+1,NRMAX+1,NG,0.0,GXMAX,STRL,0.26,MODE,IND,0)
          deallocate(GRPL)

          CALL SETFNT(IFNT)
          CALL PAGEE

       CASE('W')
          ! *** Write out numerical values of a certain variable ***
          READ(STR(3:5),*,IOSTAT=IST) NG
          IF (IST /= 0) THEN
             WRITE(6,*) '### ERROR : Invalid Command : ', STR
             CYCLE OUTER
          END IF
          CALL write_console(NG,KID2)

       CASE('M')
          ! *** Compare the already-saved graphics (Max.5) ***
          DO
             IF(KID2 == 'T') THEN
                WRITE(6,*) '## Number of time (Max = 5):'
             ELSE
                WRITE(6,*) '## Number of files (Max = 5):'
             END IF
             READ(5,*,IOSTAT=IST) NGFMAX
             IF (IST > 0) THEN
                CYCLE
             ELSEIF (IST < 0) THEN
                CYCLE OUTER
             ELSE
                EXIT
             END IF
          END DO
          IF(T_TX == 0.D0) THEN
             NGRT = 1
          ELSE
             NGRT = 0
             allocate(GYL(0:NRMAX,0:5,1:NGYRM))
             ! Current graphic data temporarily saved in GYL
             GYL(0:NRMAX,0,1:NGYRM) = GY(0:NRMAX,NGR,1:NGYRM)
          END IF
          NGR=0
          IF(KID2 == 'T') THEN
             CALL TXGLOD(IER)
             IF(IER /= 0) CYCLE OUTER
             DO NGF=1,NGFMAX
                LOOP_NGF: DO
                   WRITE(6,'(2(A,I1))') '## INPUT time: ', NGF,' / ',NGFMAX
                   READ(5,*,IOSTAT=IST) TL
                   IF (IST > 0) THEN
                      CYCLE
                   ELSEIF (IST < 0) THEN
                      CYCLE OUTER
                   ELSE
                      CALL return_NGT_from_T(NGTL,TL,IER)
                      IF(IER == 0) THEN
                         EXIT LOOP_NGF
                      ELSE
                         CYCLE LOOP_NGF
                      END IF
                   END IF
                END DO LOOP_NGF
                if(allocated(GYL) .eqv. .false.) allocate(GYL(0:NRMAX,0:5,1:NGYRM))
                GYL(0:NRMAX,NGF-NGRT,1:NGYRM) = GYT(0:NRMAX,NGTL,1:NGYRM)
             END DO
          ELSE
             DO NGF=1,NGFMAX
                CALL TXGLOD(IER)
                IF(IER /= 0) CYCLE OUTER
                if(allocated(GYL) .eqv. .false.) allocate(GYL(0:NRMAX,0:5,1:NGYRM))
                GYL(0:NRMAX,NGF-NGRT,1:NGYRM) = GY(0:NRMAX,NGR,1:NGYRM)
             END DO
          END IF
          CALL TXPRFG
          NGR = NGFMAX-NGRT
          ! Graphic data restored to GY
          GY(0:NRMAX,0:NGR,1:NGYRM) = GYL(0:NRMAX,0:NGR,1:NGYRM)
          deallocate(GYL)
          DO 
             DO
                WRITE(6,*) '## INPUT GRAPH NUMBER: A,B,C,RA,RB,RC,RD,Rn,X:exit'
                READ(5,'(A5)',IOSTAT=IST) STR2
                IF (IST > 0) THEN
                   CYCLE
                ELSEIF (IST < 0) THEN
                   CYCLE OUTER
                ELSE
                   EXIT
                END IF
             END DO
             KID1 = STR2(1:1)
             KID2 = STR2(2:2)
             CALL TOUPPER(KID1)
             CALL TOUPPER(KID2)
             SELECT CASE(KID1)
             CASE('A')
                !     Correspond to GMA
                DO I = 1, NGPRM
                   CALL TXGRFR(I,MODE)
                END DO
                !     Correspond to GMB
             CASE('B')
                DO I = 1, 6
                   CALL TXGRFR(I,MODE)
                END DO
                DO I = 9, 10
                   CALL TXGRFR(I,MODE)
                END DO
                !     Correspond to GMC
             CASE('C')
                DO
                   WRITE(6,*) '## HOW MANY GRAPHs ?: 4 or 6'
                   READ(5,'(I1)',IOSTAT=IST) NGPR
                   IF (IST > 0) THEN
                      CYCLE
                   ELSEIF (IST < 0) THEN
                      CYCLE OUTER
                   ELSE
                      IF(NGPR /= 4 .and. NGPR /= 6) CYCLE
                      EXIT
                   END IF
                END DO
                allocate(NGPRL(1:NGPR),GYL2(0:NRMAX,0:NGR),GPXY_IN(1:4))
                WRITE(6,'(A,I1,A)') '## CHOOSE ',NGPR, ' GRAPHs: '
                DO
                   READ(5,*,IOSTAT=IST) (NGPRL(I), I=1,NGPR)
                   IF (MINVAL(NGPRL) <= 0 .OR. MAXVAL(NGPRL) > NGYRM) THEN
                      WRITE(6,*) 'XX Invalid variable number! Please start over again.'
                      CYCLE
                   END IF
                   IF (IST > 0) THEN
                      CYCLE
                   ELSEIF (IST < 0) THEN
                      CYCLE OUTER
                   ELSE
                      EXIT
                   END IF
                END DO

                CALL PAGES
                CALL SETCHS(0.3, 0.0)
                CALL SETLIN(0, 1, 7)
                IF (MODEG == 2) THEN
                   IND = 9
                ELSE
                   IND = 0
                END IF
                DO I = 1, NGPR
                   WRITE(STRL,'(I)') NGPRL(I)
                   STRL = "@"//trim(adjustl(STRL))//"@"
                   CALL APPROPGY(MODEG, GY(0,0,NGPRL(I)), GYL2, STRL, NRMAX, NGR, gDIV(NGPRL(I)))
                   J = I - 1
                   IF(NGPR == 4) THEN
                      CALL TXGRFRX(J, GX, GYL2, NRMAX, NGR, STRL, MODE, IND)
                   ELSE
                      GPXY_IN(1) =  1.8  + 8.45 * MOD(J,3)
                      GPXY_IN(2) =  8.45 + 8.45 * MOD(J,3)
                      GPXY_IN(3) = 10.55 - 7.55 * REAL(J/3)
                      GPXY_IN(4) = 15.1  - 7.55 * REAL(J/3)
                      CALL TXGRFRX(J, GX, GYL2, NRMAX, NGR, STRL, MODE, IND, GPXY_IN=GPXY_IN)
                   END IF
                END DO
                CALL PAGEE

                deallocate(NGPRL,GYL2,GPXY_IN)
             CASE('R')
                SELECT CASE(KID2)
                CASE('A')
                   CALL TXGRFR(-1,MODE)
                CASE('B')
                   CALL TXGRFR(-3,MODE)
                CASE('C')
                   CALL TXGRFR(-5,MODE)
                CASE DEFAULT
                   READ(STR2(2:5),'(I4)',IOSTAT=IST) NGPR
                   IF (IST < 0) CYCLE OUTER
                   IF      (NGPR == 0) THEN
                      CYCLE OUTER
                   ELSE IF (NGPR >= 0 .AND. NGPR <= NGPRM) THEN
                      CALL TXGRFR(NGPR,MODE)
                   END IF
                END SELECT
!!$             CASE('T')
!!$                SELECT CASE(KID2)
!!$                CASE('A')
!!$                   DO NGPT = 1, NGPTM
!!$                      CALL TXGRFT(NGPT,MODE)
!!$                   END DO
!!$
!!$                CASE('B')
!!$                   DO NGPT = 1, 6
!!$                      CALL TXGRFT(NGPT,MODE)
!!$                   END DO
!!$                   DO NGPT = 9, 10
!!$                      CALL TXGRFT(NGPT,MODE)
!!$                   END DO
!!$
!!$                CASE DEFAULT
!!$                   READ(STR2(2:5),'(I4)',IOSTAT=IST) NGPT
!!$                   IF (IST < 0) THEN
!!$                      WRITE(6,*) '### ERROR : Invalid Command : ', STR
!!$                      CYCLE
!!$                   END IF
!!$                   IF (NGPT >= 1 .AND. NGPT <= NGPTM) THEN
!!$                      CALL TXGRFT(NGPT,MODE)
!!$                   END IF
!!$                END SELECT
!!$             CASE('V')
!!$                SELECT CASE(KID2)
!!$                CASE('A')
!!$                   DO NGPV = 1, NGPVM
!!$                      CALL TXGRFV(NGPV,MODE)
!!$                   END DO
!!$                CASE('B')
!!$                   DO NGPV = 1, 6
!!$                      CALL TXGRFV(NGPV,MODE)
!!$                   END DO
!!$                   DO NGPV = 9, 10
!!$                      CALL TXGRFV(NGPV,MODE)
!!$                   END DO
!!$
!!$                CASE DEFAULT
!!$                   READ(STR2(2:5),'(I4)',IOSTAT=IST) NGPV
!!$                   IF (IST < 0) THEN
!!$                      WRITE(6,*) '### ERROR : Invalid Command : ', STR
!!$                      CYCLE
!!$                   END IF
!!$                   IF (NGPV >= 1 .AND. NGPV <= NGPVM) THEN
!!$                      CALL TXGRFV(NGPV,MODE)
!!$                   END IF
!!$                END SELECT
             CASE('X')
                EXIT
             CASE DEFAULT
                WRITE(6,*) 'XX UNKNOWN GRAPHIC COMMAND'
             END SELECT
          END DO

       CASE('H')
          ! *** Animation ***
          SELECT CASE(KID2)
          CASE('A')
             CALL TXGRFRA(-1)
          CASE('B')
             CALL TXGRFRA(-3)
          CASE('C')
             CALL TXGRFRA(-5)
          CASE('D')
             CALL TXGRFRA(-9)
          CASE('N')
             DO
                WRITE(6,*) '## HOW MANY GRAPHs ?: 4 or 6'
                READ(5,'(I1)',IOSTAT=IST) NGPR
                IF (IST > 0) THEN
                   CYCLE
                ELSEIF (IST < 0) THEN
                   CYCLE OUTER
                ELSE
                   IF(NGPR /= 4 .and. NGPR /= 6) CYCLE
                   EXIT
                END IF
             END DO
             allocate(NGPRL(1:NGPR),GYL(0:NRMAX,0:NGT,1:NGPR),STRA(1:NGPR),GMAXA(1:NGPR),GMINA(1:NGPR))
             WRITE(6,'(A,I1,A)') '## CHOOSE ',NGPR, ' GRAPHs: '
             DO
                READ(5,*,IOSTAT=IST) (NGPRL(I), I=1,NGPR)
                IF (MINVAL(NGPRL) <= 0 .OR. MAXVAL(NGPRL) > NGYRM) THEN
                   WRITE(6,*) 'XX Invalid variable number! Please start over again.'
                   CYCLE
                END IF
                IF (IST > 0) THEN
                   CYCLE
                ELSEIF (IST < 0) THEN
                   CYCLE OUTER
                ELSE
                   EXIT
                END IF
             END DO

!             CALL GSTITL('//') ! Eliminate header eternally
             CALL PAGES
             CALL SETCHS(0.3, 0.0)
             CALL SETLIN(0, 1, 7)
             CALL INQFNT(IFNT)
             CALL SETFNT(44)

             IF (MODEG == 2) THEN
                IND = 9
             ELSE
                IND = 0
             END IF
             IF(NGPR == 4) THEN
                J = 5 ; GZ = 8.5
             ELSE
                J = 6 ; GZ = 7.0
             END IF
             DO I = 1, NGPR
                WRITE(STRA(I),'(I)') NGPRL(I)
                STRA(I) = "@"//trim(adjustl(STRA(I)))//"@"
                CALL APPROPGY(MODEG, GYT(0,0,NGPRL(I)), GYL(0,0,I), STRA(I), NRMAX, NGT, &
                     &        gDIV(NGPRL(I)), GMAX=GMAXA(I), GMIN=GMINA(I))
             END DO
             DO NG = 0, NGT
                call animes
                call gtextx(12.5,GZ,'@T=@',0)
                call gnumbr(13.1,GZ,GTX(NG),3,0)
                DO I = 1, NGPR
                   CALL TXGRFRS(I-1, GX, GYL(0:NRMAX,NG:NG,I), NRMAX, 1, STRA(I), 0, IND, 0, J, &
                        &       'ANIME', GMAXA(I), GMINA(I))
                END DO
                call animee
             END DO
             CALL SETFNT(IFNT)
             CALL PAGEE
             deallocate(NGPRL,GYL,STRA,GMAXA,GMINA)

          CASE DEFAULT
             READ(STR(2:5),*,IOSTAT=IST) NGYR
             IF (IST /= 0) THEN
                WRITE(6,*) '### ERROR : Invalid Command : ', STR
                CYCLE OUTER
             END IF
             IF (NGYR >= 0 .AND. NGYR <= NGYRM) THEN
                IND = 0
                IF (MODEG == 2) IND = 9
                allocate(GYL2(0:NRMAX,0:NGT))
                STRL = '@profile'//STR(2:5)//'@'
                CALL APPROPGY(MODEG, GYT(0,0,NGYR), GYL2, STRL, NRMAX, NGT, gDIV(NGYR), &
                     &        GMAX=GMAX, GMIN=GMIN)

                call pages

                CALL SETLIN(0, 1, 7)
                CALL INQFNT(IFNT)
                CALL SETFNT(44)

                do ng = 0, ngt
                   call animes
                   call gtextx(12.5,17.7,'@T=@',0)
                   call gnumbr(13.1,17.7,GTX(NG),3,0)
                   CALL TXGRFRS(0, GX, GYL2(0:NRMAX,NG:NG), NRMAX, 1, STRL, 0, IND, 0, 3, &
                        &       'ANIME', GMAX, GMIN)
                   call animee
                end do

                CALL SETFNT(IFNT)
                call pagee
                deallocate(GYL2)
             END IF
          END SELECT

       CASE('X')
          ! *** Exit ***
          RETURN

       CASE DEFAULT
          WRITE(6,*) 'XX UNKNOWN GRAPHIC COMMAND'
       END SELECT
    END DO OUTER

    RETURN
  END SUBROUTINE TXGOUT

  !***************************************************************
  !
  !   Save graphic data
  !
  !***************************************************************

  SUBROUTINE TX_GRAPH_SAVE

    use tx_commons, only : T_TX

    !  Define radial coordinate for graph

    CALL TXPRFG

    !  Store center or edge values of variables for showing time-evolution graph

    CALL TXSTGT(REAL(T_TX))

    !  Store global quantities for showing time-evolution graph

    CALL TXSTGV(REAL(T_TX))

    !  Store profile data for showing graph

    CALL TXSTGR(NGR,GT,GY,NGRM)

    !  Store balance profile data for showing graph

    IF(T_TX /= 0.D0) CALL TXSTGQ

  END SUBROUTINE TX_GRAPH_SAVE

  !***************************************************************
  !
  !   Initialize graphic axis
  !
  !***************************************************************

  SUBROUTINE TXPRFG

    use tx_commons, only : NRMAX, RHO

    !  GX(NR) : Integer

    GX(0:NRMAX) = REAL(RHO(0:NRMAX))

    RETURN
  END SUBROUTINE TXPRFG

  !***************************************************************
  !
  !   Store GY
  !
  !***************************************************************

  SUBROUTINE TXSTGR(NG,GTL,GYL,NGM)

    use tx_commons, NXM => NRMAX

    integer(4), intent(inout) :: NG
    integer(4), intent(in) :: NGM
    real(4), dimension(0:NGM), intent(out), optional :: GTL
    real(4), dimension(0:NXM,0:NGM,1:NGYRM), intent(out) :: GYL

    INTEGER(4) :: NX

    if(present(GTL)) then
       IF (NG < NGM) NG = NG + 1

       GTL(NG) = REAL(T_TX)
    end IF

    do NX = 0, NXM

       GYL(NX,NG,1)  = REAL(Var(NX,1)%n*1.D20)
       GYL(NX,NG,2)  = REAL((achg(2)*Var(NX,2)%n+achgb*PNbV(NX) &
            &               +achgb*rip_rat(NX)*PNbrpV(NX)-Var(NX,1)%n)* 1.D20)
       GYL(NX,NG,3)  = REAL(Var(NX,1)%Ur)
       GYL(NX,NG,4)  = REAL(Var(NX,1)%Uth)
       GYL(NX,NG,5)  = REAL(Var(NX,1)%Uph)
       GYL(NX,NG,6)  = REAL(Var(NX,2)%Ur)
       GYL(NX,NG,7)  = REAL(Var(NX,2)%Uth)
       GYL(NX,NG,8)  = REAL(Var(NX,2)%Uph)
       GYL(NX,NG,9)  = REAL(ErV(NX))
       GYL(NX,NG,10) = REAL(BthV(NX))
       GYL(NX,NG,11) = REAL(Etor(NX))
       GYL(NX,NG,12) = REAL(PNbV(NX) * 1.D20)
       GYL(NX,NG,13) = REAL(UbphV(NX))
       GYL(NX,NG,14) = REAL(Var(NX,1)%T)
       GYL(NX,NG,15) = REAL(Var(NX,2)%T)
       GYL(NX,NG,16) = REAL((PN01V(NX)+PN02V(NX)+PN03V(NX))*1.D20)
       GYL(NX,NG,17) = REAL(BEpol(NX))
       GYL(NX,NG,18) = real(fipol(NX))
       GYL(NX,NG,19) = REAL(BUbparV(NX))
       GYL(NX,NG,20) = REAL(Q(NX))

       GYL(NX,NG,21) = REAL(AJ(NX))
       GYL(NX,NG,22) = REAL(aee*achg(1)*Var(NX,1)%n*Var(NX,1)%UphR*1.d20/d_rrr(NX))
       GYL(NX,NG,23) = REAL(aee*achg(2)*Var(NX,2)%n*Var(NX,2)%UphR*1.d20/d_rrr(NX))
       GYL(NX,NG,24) = REAL(aee*achgb*PNbV(NX)*(aat(NX)*fipol(NX)/bbt(NX)*BUbparV(NX))*1.d20/d_rrr(NX))

       GYL(NX,NG,25) = real(Var(NX,1)%RUph/fipol(NX)-Var(NX,1)%BUpar/bbt(NX))
       GYL(NX,NG,26) = REAL(Var(NX,1)%BUpar)
       GYL(NX,NG,27) = real(Var(NX,2)%RUph/fipol(NX)-Var(NX,2)%BUpar/bbt(NX))
       GYL(NX,NG,28) = REAL(Var(NX,2)%BUpar)
       GYL(NX,NG,29) = REAL(Di(NX)+De(NX))
!!$       IF (rG1h2(NX) > 3.D0) THEN
!!$          GYL(NX,NG,30) = 3.0
!!$       ELSE
!!$          GYL(NX,NG,30) = REAL(rG1h2(NX))
!!$       END IF
       GYL(NX,NG,30) = REAL(rG1h2(NX))
       GYL(NX,NG,31) = REAL(FCDBM(NX))
       GYL(NX,NG,32) = REAL(S(NX))
       GYL(NX,NG,33) = REAL(Alpha(NX))
       GYL(NX,NG,34) = REAL(rKappa(NX))
       GYL(NX,NG,35) = REAL(PN01V(NX)*1.D20)
       GYL(NX,NG,36) = REAL(PN02V(NX)*1.D20)
       GYL(NX,NG,37) = REAL(Var(NX,2)%n*1.D20)

       !  Raw data

       GYL(NX,NG,38) = real(X(NX,LQm1))
       GYL(NX,NG,39) = real(X(NX,LQm2))
       GYL(NX,NG,40) = real(X(NX,LQm3))
       GYL(NX,NG,41) = real(X(NX,LQm4))
       GYL(NX,NG,42) = real(X(NX,LQm5))
       GYL(NX,NG,43) = real(X(NX,LQe1))
       GYL(NX,NG,44) = real(X(NX,LQe2))
       GYL(NX,NG,45) = real(X(NX,LQe3))
       GYL(NX,NG,46) = real(X(NX,LQe4))
       GYL(NX,NG,47) = real(X(NX,LQe5))
       GYL(NX,NG,48) = real(X(NX,LQe6))
       GYL(NX,NG,49) = real(X(NX,LQe7))
       GYL(NX,NG,50) = real(X(NX,LQi1))
       GYL(NX,NG,51) = real(X(NX,LQi2))
       GYL(NX,NG,52) = real(X(NX,LQi3))
       GYL(NX,NG,53) = real(X(NX,LQi4))
       GYL(NX,NG,54) = real(X(NX,LQi5))
       GYL(NX,NG,55) = real(X(NX,LQi6))
       GYL(NX,NG,56) = real(X(NX,LQi7))
       GYL(NX,NG,57) = real(X(NX,LQb1))
       GYL(NX,NG,58) = real(X(NX,LQb3))
       GYL(NX,NG,59) = real(X(NX,LQn1))
       GYL(NX,NG,60) = real(X(NX,LQn2))
       GYL(NX,NG,61) = real(X(NX,LQn3))
       GYL(NX,NG,62) = real(X(NX,LQr1))

       !  Coefficients

       GYL(NX,NG,63) = real(rNuL(NX))
!       GYL(NX,NG,64) = 
!       GYL(NX,NG,65) = 
       GYL(NX,NG,66) = REAL(WPM(NX))
       GYL(NX,NG,67) = REAL(rNuTei(NX))
       GYL(NX,NG,68) = REAL(rNu0e(NX))
       GYL(NX,NG,69) = REAL(rNu0i(NX))
       GYL(NX,NG,70) = REAL(Chie(NX))
       GYL(NX,NG,71) = REAL(Chii(NX))
       GYL(NX,NG,72) = REAL(D01(NX))
       GYL(NX,NG,73) = REAL(D02(NX))
       GYL(NX,NG,74) = real(rMue(NX))
       GYL(NX,NG,75) = real(rMui(NX))
!       GYL(NX,NG,76) = REAL(rNueHL(NX))
!       GYL(NX,NG,77) = REAL(rNuiHL(NX))

       GYL(NX,NG,78) = REAL(rNuiCX(NX))
       GYL(NX,NG,79) = REAL(rNuB(NX))
       GYL(NX,NG,80) = REAL(rNuION(NX))
       GYL(NX,NG,81) = REAL(SiLC(NX))
       GYL(NX,NG,82) = REAL(rNuOL(NX))
       GYL(NX,NG,83) = REAL(Deff(NX))
       GYL(NX,NG,84) = REAL(SNB(NX))
       GYL(NX,NG,85) = REAL(PRFe(NX))
       GYL(NX,NG,86) = REAL(PRFi(NX))
       GYL(NX,NG,87) = REAL(rNu0b(NX))

       GYL(NX,NG,88) = REAL(PIE(NX))
       GYL(NX,NG,89) = REAL(PCX(NX))
       GYL(NX,NG,90) = REAL(SIE(NX))
       GYL(NX,NG,91) = REAL(PBr(NX))

       GYL(NX,NG,92) = REAL(rNuLTe(NX))
       GYL(NX,NG,93) = REAL(rNuLTi(NX))
       GYL(NX,NG,94) = REAL(rNuAse(NX))
       GYL(NX,NG,95) = REAL(rNuASi(NX))

       GYL(NX,NG,96) = REAL(PNBe(NX))
       GYL(NX,NG,97) = REAL(PNBi(NX))
       GYL(NX,NG,98) = REAL(POHe(NX))
       GYL(NX,NG,99) = REAL(POHi(NX))
       GYL(NX,NG,100) = REAL(PEQe(NX))
       GYL(NX,NG,101) = REAL(PEQi(NX))

       GYL(NX,NG,102) = REAL(Var(NX,1)%p*1.D20*rKeV)
       GYL(NX,NG,103) = REAL(Var(NX,2)%p*1.D20*rKeV)
       GYL(NX,NG,104) = REAL(PNB(NX))
       GYL(NX,NG,105) = REAL(PNBPD(NX)/(Eb*rKeV*1.D20))
       GYL(NX,NG,106) = REAL(PNBTG(NX)/(Eb*rKeV*1.D20))
       GYL(NX,NG,107) = REAL(SNBPDi(NX))
       GYL(NX,NG,108) = REAL(SNBTGi(NX))
       GYL(NX,NG,109) = real(Vbpara(NX))

       ! *** Ripple loss part ******************************************

       GYL(NX,NG,110) = real(PNbrpV(NX) * 1.D20) 
       GYL(NX,NG,111) = REAL(rNubrp1(NX))
       GYL(NX,NG,112) = REAL(DltRP(NX))
       GYL(NX,NG,113) = REAL(Ubrp(NX))
       GYL(NX,NG,114) = REAL(Dbrp(NX))
       IF(PNbV(NX) == 0.D0) THEN
          GYL(NX,NG,115) = 0.D0
       ELSE
          GYL(NX,NG,115) = REAL(PNbrpV(NX)/PNbV(NX))
       END IF
!       GYL(NX,NG,115) = REAL(PNbrpLV(NX)*1.D20)

       GYL(NX,NG,116) = REAL(rNuLB(NX))
       GYL(NX,NG,117) = real(aee*1.D20*sum(achg(:)*Var(NX,:)%n*Var(NX,:)%Ur))
       ! ***************************************************************

       GYL(NX,NG,118) = REAL(rNubL(NX))
       ! External toroidal NBI torque density
       GYL(NX,NG,119) = real(BSmb(NX)*fipol(NX)/bbt(NX)*1.D20)
       ! Generated toroidal torque density
       GYL(NX,NG,120) = real(sdt(NX)*aee*1.d20*sum(achg(:)*Var(NX,:)%n*Var(NX,:)%UrV)) ! <j.grad psi>

       GYL(NX,NG,121) = REAL(BJPARA(NX))
       GYL(NX,NG,122) = REAL(BEpara(NX))
       GYL(NX,NG,123) = REAL(PALFe(NX))
       GYL(NX,NG,124) = REAL(PALFi(NX))

       ! *** Particle diffusion due to magnetic braidng ***AF (2008-06-08)
       GYL(NX,NG,125) = REAL(DMAG(NX))
       GYL(NX,NG,126) = REAL(DMAGe(NX))
       GYL(NX,NG,127) = REAL(DMAGi(NX))

       ! *** Rho vs Psi *************************************************
       GYL(NX,NG,128) = real(Psiv(NX)-Psiv(0))/(Psiv(NRA)-Psiv(0))

       ! *** Effective neoclassical thermal diffusivity *****************
       GYL(NX,NG,129) = REAL(ChiNCpe(NX)+ChiNCte(NX))
       GYL(NX,NG,130) = REAL(ChiNCpi(NX)+ChiNCti(NX))

       ! *** Additional torque *********************************************
       GYL(NX,NG,131) = REAL(Tqt(NX))

       ! *** Turbulent pinch velocity *********************************************
       GYL(NX,NG,132) = REAL(VWpch(NX))

       GYL(NX,NG,133) = REAL(SCX(NX))

       ! *** Total thermal diffusivity *****************
       GYL(NX,NG,134) = REAL(Chie(NX)+ChiNCpe(NX)+ChiNCte(NX))
       GYL(NX,NG,135) = REAL(Chii(NX)+ChiNCpi(NX)+ChiNCti(NX))

       ! *** Particle flux <Gamma_s . grad rho>
       if(NX == 0) then
          GYL(NX,NG,136) = 0.0
          GYL(NX,NG,137) = 0.0
       else
          GYL(NX,NG,136) = real(Var(NX,1)%n*Var(NX,1)%UrV/vro(NX))*1.e20
          GYL(NX,NG,137) = real(Var(NX,2)%n*Var(NX,2)%UrV/vro(NX))*1.e20
       end if

       ! *** Halo neutral ***
       GYL(NX,NG,138) = REAL(PN03V(NX)*1.D20)
       GYL(NX,NG,139) = REAL(D03(NX))

       ! *** CDIM mode 09/07/13 miki_m ***
       GYL(NX,NG,140) = REAL(FCDIM(NX))

       ! *** ExB shearing rate ***
       GYL(NX,NG,141) = REAL(wexb(NX))

       GYL(NX,NG,142) = real(BVsdiag(NX,1,1)) ! elec. particle
       GYL(NX,NG,143) = real(BVsdiag(NX,2,1)) ! ion   particle
       GYL(NX,NG,144) = real(BVsdiag(NX,1,2)) ! elec. heat
       GYL(NX,NG,145) = real(BVsdiag(NX,2,2)) ! ion   heat

       GYL(NX,NG,146) = real(xmu(NX,1,1,1))
       GYL(NX,NG,147) = real(xmu(NX,1,1,2))
       GYL(NX,NG,148) = real(xmu(NX,1,2,1))
       GYL(NX,NG,149) = real(xmu(NX,1,2,2))
       GYL(NX,NG,150) = real(xmu(NX,2,1,1))
       GYL(NX,NG,151) = real(xmu(NX,2,1,2))
       GYL(NX,NG,152) = real(xmu(NX,2,2,1))
       GYL(NX,NG,153) = real(xmu(NX,2,2,2))

       GYL(NX,NG,154) = real(Var(NX,1)%Uthhat) ! elec.
       GYL(NX,NG,155) = real(Var(NX,2)%Uthhat) ! ion
       GYL(NX,NG,156) = real(Var(NX,1)%UrV)    ! elec.
       GYL(NX,NG,157) = real(Var(NX,2)%UrV)    ! ion
!       GYL(NX,NG,156) = real(Var(NX,1)%UrV-UgV(NX))    ! elec.
!       GYL(NX,NG,157) = real(Var(NX,2)%UrV-UgV(NX))    ! ion
       GYL(NX,NG,158) = real(Var(NX,1)%Bqpar)  ! parallel ele. heat flow
       GYL(NX,NG,159) = real(Var(NX,2)%Bqpar)  ! parallel ion  heat flow
       GYL(NX,NG,160) = real(aee*sum(achg(:)*Var(NX,:)%n*Var(NX,:)%UrV)*1.d20) ! e_s<Gamma_s.grad V>
       GYL(NX,NG,161) = real(sdt(NX))

       GYL(NX,NG,162) = real(UsparNCL(NX,1,1)) ! Neo. solver: parallel ele. flow
       GYL(NX,NG,163) = real(UsparNCL(NX,2,1)) ! Neo. solver: parallel ion  flow
       GYL(NX,NG,164) = real(UsparNCL(NX,1,2)) ! Neo. solver: parallel ele. heat flow
       GYL(NX,NG,165) = real(UsparNCL(NX,2,2)) ! Neo. solver: parallel ion  heat flow

       GYL(NX,NG,166) = real(lab(NX,1,1,1,1)) ! ee 11
       GYL(NX,NG,167) = real(lab(NX,1,1,1,2)) ! ee 12
       GYL(NX,NG,168) = real(lab(NX,1,1,2,1)) ! ee 21
       GYL(NX,NG,169) = real(lab(NX,1,1,2,2)) ! ee 22
       GYL(NX,NG,170) = real(lab(NX,1,2,1,1)) ! ei 11
       GYL(NX,NG,171) = real(lab(NX,1,2,1,2)) ! ei 12
       GYL(NX,NG,172) = real(lab(NX,1,2,2,1)) ! ei 21
       GYL(NX,NG,173) = real(lab(NX,1,2,2,2)) ! ei 22
       GYL(NX,NG,174) = real(lab(NX,2,1,1,1)) ! ie 11
       GYL(NX,NG,175) = real(lab(NX,2,1,1,2)) ! ie 12
       GYL(NX,NG,176) = real(lab(NX,2,1,2,1)) ! ie 21
       GYL(NX,NG,177) = real(lab(NX,2,1,2,2)) ! ie 22
       GYL(NX,NG,178) = real(lab(NX,2,2,1,1)) ! ii 11
       GYL(NX,NG,179) = real(lab(NX,2,2,1,2)) ! ii 12
       GYL(NX,NG,180) = real(lab(NX,2,2,2,1)) ! ii 21
       GYL(NX,NG,181) = real(lab(NX,2,2,2,2)) ! ii 22

       GYL(NX,NG,182) = real(laf(NX,1,1,1))   ! ef 11
       GYL(NX,NG,183) = real(laf(NX,1,2,1))   ! ef 21
       GYL(NX,NG,184) = real(lfb(NX,1,1,1))   ! fe 11, not normalized
       GYL(NX,NG,185) = real(laf(NX,2,1,1))   ! if 11
       GYL(NX,NG,186) = real(lfb(NX,2,1,1))   ! fi 11, not normalized
       GYL(NX,NG,187) = real(lff(NX,1,1))     ! ff 11, not normalized
       GYL(NX,NG,188) = real(xmuf(NX,1))      ! fast-ion viscosity

       GYL(NX,NG,189) = real(Vmps(NX,1))
       GYL(NX,NG,190) = real(Vmps(NX,2))

       GYL(NX,NG,191) = real(UgV(NX))
       GYL(NX,NG,192) = real(qhatsq(NX))

       GYL(NX,NG,193) = real(sqrt(bbt(NX)))
       GYL(NX,NG,194) = real(ft(NX))
       GYL(NX,NG,195) = real(vro(NX))

       GYL(NX,NG,196) = real(X(NX,LQb4))
       GYL(NX,NG,197) = real(X(NX,LQb7))
       GYL(NX,NG,198) = real(X(NX,LQb2))
       GYL(NX,NG,199) = real(aee*(sum(achg(:)*Var(NX,:)%n*Var(NX,:)%UrV)+achgb*PNbV(NX)*UbrVV(NX))*1.d20) ! e_s<Gamma_s.grad V>
       if(NX == 0) then
          GYL(NX,NG,200) = 0.0
       else
          GYL(NX,NG,200) = real(PNbV(NX)*UbrVV(NX)/vro(NX))*1.e20
       end if

    end do

    RETURN
  END SUBROUTINE TXSTGR

  !***************************************************************
  !
  !   Store GVY
  !
  !***************************************************************

  SUBROUTINE TXSTGV(GTIME)

    use tx_commons, only : AEE, PI, achg, achgb, PNbV, &
         &              rip_rat, PNbrpV, NRC, &
         &              ErV, BthV, Etor, NRMAX, UbphV, &
         &              PN01V, PN02V, PN03V, BEpol, BUbparV, Q, Di, De, &
         &              rG1h2, FCDBM, S, Alpha, rKappa, NRA, rNuION, &
         &              Chie,  Chii, PIE, PCX, SIE, PBr, Deff, rKeV, FCDIM, &
         &              AJ, d_rrr, aat, fipol, bbt, BUbparV, qhatsq, Var
    use tx_interface, only : rLINEAVE, sub_intg_vol

    REAL(4), INTENT(IN) :: GTIME
    integer(4) :: NRL
    REAL(8) :: PNESUM1, PNESUM2
    real(8), dimension(:), allocatable :: PNeION

    IF (NGVV < NGVM) NGVV=NGVV+1

    GVX(NGVV) = GTIME

    GVY(NGVV,1)  = REAL(Var(0,1)%n * 1.D20)
    GVY(NGVV,2)  = REAL((achg(2) * Var(0,2)%n + achgb * PNbV(0) + achgb * rip_rat(0) * PNbrpV(0) &
         &                       - Var(0,1)%n) * 1.D20)
    GVY(NGVV,3)  = REAL(Var(NRC,1)%Ur)
    GVY(NGVV,4)  = REAL(Var(NRC,1)%Uth)
    GVY(NGVV,5)  = REAL(Var(0,1)%Uph)
    GVY(NGVV,6)  = REAL(Var(NRC,2)%Ur)
    GVY(NGVV,7)  = REAL(Var(NRC,2)%Uth)
    GVY(NGVV,8)  = REAL(Var(NRC,2)%Uph)
    GVY(NGVV,9)  = REAL(ErV(NRC))
!!    GVY(NGVV,9)  = REAL(ErV(29)) ! rho=0.8 when nrmax=60
    GVY(NGVV,10) = REAL(BthV(NRA))
    GVY(NGVV,11) = REAL(Etor(NRMAX))
    GVY(NGVV,12) = REAL(PNbV(0) * 1.D20)
    GVY(NGVV,13) = REAL(UbphV(0))
    GVY(NGVV,14) = REAL(Var(0,1)%T)
    GVY(NGVV,15) = REAL(Var(0,2)%T)
    GVY(NGVV,16) = REAL((PN01V(NRA) + PN02V(NRA) + PN03V(NRA)) * 1.D20)
    GVY(NGVV,17) = REAL(BEpol(NRC))
    GVY(NGVV,18) = REAL(fipol(0))
    GVY(NGVV,19) = REAL(BUbparV(NRC))
    GVY(NGVV,20) = REAL(Q(0))

    GVY(NGVV,21) = REAL(AJ(0))
    GVY(NGVV,22) = REAL(aee*achg(1)*Var(0,1)%n*Var(0,1)%UphR*1.d20/d_rrr(0))
    GVY(NGVV,23) = REAL(aee*achg(2)*Var(0,2)%n*Var(0,2)%UphR*1.d20/d_rrr(0))
    GVY(NGVV,24) = REAL(aee*achgb*PNbV(0)*(aat(0)*fipol(0)/bbt(0)*BUbparV(0))*1.d20/d_rrr(0))

    GVY(NGVV,25) = real(Var(NRC,1)%RUph/fipol(NRC)-Var(NRC,1)%BUpar/bbt(NRC))
    GVY(NGVV,26) = real(Var(NRC,1)%BUpar)
    GVY(NGVV,27) = real(Var(NRC,2)%RUph/fipol(NRC)-Var(NRC,2)%BUpar/bbt(NRC))
    GVY(NGVV,28) = real(Var(NRC,2)%BUpar)
    GVY(NGVV,29) = REAL(Di(NRC)+De(NRC))
    GVY(NGVV,30) = REAL(rG1h2(NRC))
    GVY(NGVV,31) = REAL(FCDBM(NRC))
    GVY(NGVV,32) = REAL(S(NRC))
    GVY(NGVV,33) = REAL(Alpha(NRC))
    GVY(NGVV,34) = REAL(rKappa(NRC))
    GVY(NGVV,35) = REAL(PN01V(NRA) * 1.D20)
    GVY(NGVV,36) = REAL(PN02V(0) * 1.D20)
    GVY(NGVV,37) = REAL(Var(0,2)%n * 1.D20)

    GVY(NGVV,38)  = REAL(rLINEAVE(0.D0))
    GVY(NGVV,39)  = REAL(rLINEAVE(0.24D0))
    GVY(NGVV,40)  = REAL(rLINEAVE(0.6D0))

    call sub_intg_vol(Var(0:NRMAX,1)%n,NRA,PNESUM1)
    allocate(PNeION(0:NRMAX))
    PNeION(0:NRMAX) = Var(0:NRMAX,1)%n*rNuION(0:NRMAX)
    call sub_intg_vol(PNeION,NRA,PNESUM2)
    deallocate(PNeION)

    GVY(NGVV,41) = REAL(PNESUM1)
    GVY(NGVV,42) = REAL(PNESUM2)
    IF(NGVV == 0.OR.ABS(PNESUM2) <= 0.D0) THEN
       GVY(NGVV,43) = 0.0
    ELSE
       GVY(NGVV,43) = REAL(PNESUM1/PNESUM2)
    END IF

    GVY(NGVV,44) = REAL(Chie(NRC))
    GVY(NGVV,45) = REAL(Chii(NRC))
    GVY(NGVV,46) = REAL(PIE(NRC))
    GVY(NGVV,47) = REAL(PCX(NRC))
    GVY(NGVV,48) = REAL(SIE(NRC))
    GVY(NGVV,49) = REAL(PBr(NRC))

    GVY(NGVV,50) = REAL(PNbrpV(0) * 1.D20)
    GVY(NGVV,51) = REAL(PN03V(0) * 1.D20)
    GVY(NGVV,52) = REAL(Deff(NRC))
    GVY(NGVV,53) = REAL(Var(0,1)%p * 1.D20 * rKeV)
    GVY(NGVV,54) = REAL(Var(0,2)%p * 1.D20 * rKeV)

    GVY(NGVV,55) = REAL(FCDIM(NRC)) !09/07/13 miki_m
    GVY(NGVV,56) = REAL( achg(1)*AEE*Var(NRC,1)%n*Var(NRC,1)%Ur &
         &              +achg(2)*AEE*Var(NRC,2)%n*Var(NRC,2)%Ur)*1.D20

    NRL = NRC
!    NRL = 29
    GVY(NGVV,57) = real(Var(NRL,2)%Uthhat)
    GVY(NGVV,58) = 1.0+2.0*real(qhatsq(NRL))

    RETURN
  END SUBROUTINE TXSTGV

  !***************************************************************
  !
  !   Store GTY
  !
  !***************************************************************

  SUBROUTINE TXSTGT(GTIME)

    use tx_commons, only : TS0, TSAV, PINT, POHT, PNBT, PRFT, PRFTe, PRFTi,PNFT, &
         &              AJT, AJOHT, AJNBT, AJBST, POUT, PCXT, PIET, QF, ANS0, &
         &              ANSAV, WPT, WBULKT, WST, TAUE1, TAUE2, TAUEP, BETAA, &
         &              BETA0, BETAPA, BETAP0, VLOOP, ALI, Q, RQ1, ANF0, ANFAV, &
         &              VOLAVN, TAUP, TAUPA, Gamma_a, TNBcol, TTqt
    REAL(4), INTENT(IN) :: GTIME

    IF (NGT < NGTM) NGT=NGT+1

    GTX(NGT) = GTIME

    GTY(NGT,1)  = REAL(TS0(1))
    GTY(NGT,2)  = REAL(TS0(2))
    GTY(NGT,3)  = REAL(TSAV(1))
    GTY(NGT,4)  = REAL(TSAV(2))
    GTY(NGT,5)  = REAL(PINT)
    GTY(NGT,6)  = REAL(POHT)
    GTY(NGT,7)  = REAL(PNBT)
    GTY(NGT,8)  = REAL(PRFT)
    GTY(NGT,9)  = REAL(PNFT)

    GTY(NGT,10) = REAL(AJT)
    GTY(NGT,11) = REAL(AJOHT)
    GTY(NGT,12) = REAL(AJNBT)
    GTY(NGT,13) = REAL(AJBST)
    GTY(NGT,14) = REAL(AJOHT+AJBST+AJNBT)
    GTY(NGT,15) = REAL(POUT)
    GTY(NGT,16) = REAL(PCXT)
    GTY(NGT,17) = REAL(PIET)
    !  GTY(NGT,18) = REAL(PRLT)
    !  GTY(NGT,19) = REAL(PCONT)
    GTY(NGT,20) = REAL(QF)
    GTY(NGT,21) = REAL(TS0(1))
    GTY(NGT,22) = REAL(TS0(2))
    GTY(NGT,23) = REAL(TSAV(1))
    GTY(NGT,24) = REAL(TSAV(2))
    GTY(NGT,25) = REAL(ANS0(1))
    GTY(NGT,26) = REAL(ANS0(2))
    GTY(NGT,27) = REAL(ANSAV(1))
    GTY(NGT,28) = REAL(ANSAV(2))
    GTY(NGT,29) = REAL(WPT)
    GTY(NGT,30) = REAL(WBULKT)

    GTY(NGT,31) = REAL(WST(1))
    GTY(NGT,32) = REAL(WST(2))

    GTY(NGT,33) = REAL(TAUE1)
    GTY(NGT,34) = REAL(TAUE2)
    GTY(NGT,35) = REAL(TAUEP)
    GTY(NGT,36) = REAL(BETAA) * 100.0
    GTY(NGT,37) = REAL(BETA0) * 100.0
    GTY(NGT,38) = REAL(BETAPA)
    GTY(NGT,39) = REAL(BETAP0)
    GTY(NGT,40) = REAL(VLOOP)
    GTY(NGT,41) = REAL(ALI)
    GTY(NGT,42) = REAL(Q(0))
    GTY(NGT,43) = REAL(RQ1)

    GTY(NGT,44) = REAL(ANF0(1))  * 100.0
    GTY(NGT,45) = REAL(ANFAV(1)) * 100.0

    GTY(NGT,46) = REAL(VOLAVN)
    GTY(NGT,47) = REAL(PRFTe)
    GTY(NGT,48) = REAL(PRFTi)
    ! Initial value becomes too big because almost no neutrals exist.
    IF(NGT == 0) THEN
       GTY(NGT,49) = 0.0
       GTY(NGT,50) = 0.0
       GTY(NGT,51) = 0.0
    ELSE
       GTY(NGT,49) = REAL(TAUP)
       GTY(NGT,50) = REAL(TAUPA)
       GTY(NGT,51) = REAL(Gamma_a) * 1.E-20
    END IF
    GTY(NGT,52) = real(TNBcol)
    GTY(NGT,53) = real(TTqt)

    ! Store data for 3D graphics
    CALL TXSTGR(NGT,GYL=GYT,NGM=NGTM)

    RETURN
  END SUBROUTINE TXSTGT

  !***************************************************************
  !
  !   Store GQY
  !
  !***************************************************************

  SUBROUTINE TXSTGQ

    use tx_commons, only : NRMAX, NQMAX, NLCMAX, NLC, NLCR, ALC, BLC, CLC, PLC, X!, t_tx,lqi3,lqi4,lqi2
    integer(4) :: NR, NC, NQ, NC1

    DO NQ = 1, NQMAX
       DO NC = 1, NLCMAX(NQ)
          NR = 0
          NC1 = NLCR(NC,NQ,0)
          IF(NC1 == 0) THEN
             GQY(NR,NC,NQ) = REAL(PLC(NR,NC,NQ))
          ELSE
             GQY(NR,NC,NQ) = REAL(  BLC(NR,NC,NQ) * X(NR  ,NC1) &
                  &               + ALC(NR,NC,NQ) * X(NR+1,NC1) &
                  &               + PLC(NR,NC,NQ))
          END IF

          NC1 = NLC(NC,NQ)
          DO NR = 1, NRMAX - 1
             IF(NC1 == 0) THEN
                GQY(NR,NC,NQ) = REAL(PLC(NR,NC,NQ))
             ELSE
                GQY(NR,NC,NQ) = REAL(  CLC(NR,NC,NQ) * X(NR-1,NC1) &
                     &               + BLC(NR,NC,NQ) * X(NR  ,NC1) &
                     &               + ALC(NR,NC,NQ) * X(NR+1,NC1) &
                     &               + PLC(NR,NC,NQ))
             END IF
          END DO

          NR = NRMAX
          NC1 = NLCR(NC,NQ,1)
          IF(NC1 == 0) THEN
             GQY(NR,NC,NQ) = REAL(PLC(NR,NC,NQ))
          ELSE
             GQY(NR,NC,NQ) = REAL(  CLC(NR,NC,NQ) * X(NR-1,NC1) &
                  &               + BLC(NR,NC,NQ) * X(NR  ,NC1) &
                  &               + PLC(NR,NC,NQ))
          END IF
       END DO
    END DO
!    write(6,*) real(t_tx),gqy(38,3,lqi2),gqy(38,4,lqi2),gqy(38,5,lqi2),gqy(38,6,lqi2)
!lqi3    write(6,*) real(t_tx),gqy(38,5,lqi3),gqy(38,6,lqi3),sum(gqy(38,13:16,lqi3))
!lqi4    write(6,*) real(t_tx),gqy(38,4,lqi4),gqy(38,5,lqi4),sum(gqy(38,6:7,lqi4)),sum(gqy(38,12:15,lqi4))

  END SUBROUTINE TXSTGQ

  !***************************************************************
  !
  !   Time evolution of radial profile
  !
  !***************************************************************

  SUBROUTINE TXGRFR(NGYRIN,MODE)

    use tx_commons, only : NRMAX, DT, rho, NEMAX, NRA, hv, vlt, pi
    INTEGER(4), INTENT(IN) :: MODE
    INTEGER(4), INTENT(IN) :: NGYRIN
    INTEGER(4) :: IND, NG, NR, NGYR, NE, IFNT, NRMAXL
    REAL(4), DIMENSION(:,:), allocatable :: GYL, GYL2
    character(len=50) :: STR
    character(len=1) :: KSTR
    character(len=3) :: KEND,KLABEL
    real(4), dimension(4) :: GPX, GPY
    real(4) :: GPXL, FACT, GYMAX

    NGYR = NGYRIN

    IF (NGR <= -1) THEN
       WRITE(6,*) 'G', NGYR, ' has no data'
       RETURN
    END IF

    IF (MODEG == 2) THEN
       IND = 9
    ELSE
       IND = 0
    END IF

    allocate(GYL(0:NRMAX,0:NGR))

    DO

    CALL PAGES
    CALL SETCHS(0.3, 0.0)
    CALL SETLIN(0, 1, 7)
    CALL INQFNT(IFNT)

    CALL MOVE(2.0,17.7)
    CALL TEXT('[G', 2)
    CALL NUMBI(NGYR, '(I2)', 2)
    CALL TEXT(']  ', 3)
    CALL TEXT('FROM', 4)
    CALL NUMBR(GT(0), '(1PE9.2)', 9)
    CALL TEXT(' TO', 3)
    CALL NUMBR(GT(NGR), '(1PE9.2)', 9)
    CALL TEXT('  DT =', 6)
    CALL NUMBD(DT, '(1PD9.2)', 9)
    CALL TEXT('  NGRSTP = ', 11)
    CALL NUMBI(NGRSTP,'(I4)',4)

    SELECT CASE(NGYR)
    CASE(0)
       ! Draw frame of rho
       GPX(1) =  2.5 ; GPY(1) = 10.5
       GPX(2) =  2.5 ; GPY(2) = 16.5
       GPX(3) = 24.5 ; GPY(3) = 16.5
       GPX(4) = 24.5 ; GPY(4) = 10.5
       CALL LINESC(GPX,GPY,4)
       ! Draw lines
       FACT = (GPX(4) - GPX(1)) / rho(NRMAX)
       GPXL = GPX(1)
       DO NE = 1, NEMAX-1
          GPXL = GPXL + REAL((rho(NE)-rho(NE-1))) * FACT
          IF(NE == NRA) CALL SETLIN(-1,-1,6)
          CALL LINE1(GPXL,10.5,GPXL,16.5)
          IF(NE == NRA) CALL SETLIN(-1,-1,7)
       END DO
       WRITE(KSTR,'(I1)') 0
       WRITE(KEND,'(I3)') NRMAX
       KLABEL='$#r'
       CALL SETCHS(0.4,0.0)
       CALL GTEXT(GPX(1),GPY(1)-0.5,KSTR,1,2)
       CALL GTEXT(GPX(4),GPY(4)-0.5,KEND,3,2)
       CALL SETFNT(33)
       CALL SETCHS(0.6,0.0)
       CALL GTEXT(GPX(1)-1.0,0.5*(GPY(1)+GPY(2)),KLABEL,3,2)
       CALL SETFNT(IFNT)

       ! Draw frame of vv
       GPX(1) =  2.5 ; GPY(1) =  2.0
       GPX(2) =  2.5 ; GPY(2) =  8.0
       GPX(3) = 24.5 ; GPY(3) =  8.0
       GPX(4) = 24.5 ; GPY(4) =  2.0
       CALL LINESC(GPX,GPY,4)
       ! Draw lines
       FACT = (GPX(4) - GPX(1)) / vlt(NRMAX)
       GPXL = GPX(1)
       DO NE = 1, NEMAX-1
          GPXL = GPXL + REAL(hv(NE)) * FACT
          IF(NE == NRA) CALL SETLIN(-1,-1,6)
          CALL LINE1(GPXL,2.0,GPXL,8.0)
          IF(NE == NRA) CALL SETLIN(-1,-1,7)
       END DO
       WRITE(KSTR,'(I1)') 0
       WRITE(KEND,'(I3)') NRMAX
       KLABEL='V'
       CALL SETCHS(0.4,0.0)
       CALL GTEXT(GPX(1),GPY(1)-0.5,KSTR,1,2)
       CALL GTEXT(GPX(4),GPY(4)-0.5,KEND,3,2)
       CALL SETFNT(33)
       CALL SETCHS(0.6,0.0)
       CALL GTEXT(GPX(1)-1.0,0.5*(GPY(1)+GPY(2)),KLABEL,1,2)
       CALL SETFNT(IFNT)
       CALL SETCHS(0.3,0.0)

    CASE(1)
       STR = '@n$-e$=@'
       CALL APPROPGY(MODEG, GY(0,0,1), GYL, STR, NRMAX, NGR, gDIV(1))
       CALL TXGRFRX(0, GX, GYL, NRMAX, NGR, STR, MODE, IND)!, GYMAX=0.2)

       STR = '@Z*n$-i$=+Z*n$-b$=+Z*n$-brp$=-n$-e$=@'
       CALL APPROPGY(MODEG, GY(0,0,2), GYL, STR, NRMAX, NGR, gDIV(2))
       CALL TXGRFRX(1, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       STR = '@E$-r$=@'
       CALL APPROPGY(MODEG, GY(0,0,9), GYL, STR, NRMAX, NGR, gDIV(9))
       CALL TXGRFRX(2, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       CALL TXWPGR

    CASE(2)
       STR = '@u$-er$=@'
       CALL TXGRFRX(0, GX, GY(0,0,3), NRMAX, NGR, STR, MODE, IND)!,GYMAX=0.4,GYMIN=0.0)

       STR = '@u$-e$#q$#$=@'
       CALL APPROPGY(MODEG, GY(0,0,4), GYL, STR, NRMAX, NGR, gDIV(4))
       CALL TXGRFRX(1, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       STR = '@u$-e$#f$#$=@'
       CALL APPROPGY(MODEG, GY(0,0,5), GYL, STR, NRMAX, NGR, gDIV(5))
       CALL TXGRFRX(2, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       STR = '@uhat$-e$#q$#$=@'
       CALL APPROPGY(MODEG, GY(0,0,154), GYL, STR, NRMAX, NGR, gDIV(154))
       CALL TXGRFRX(3, GX, GYL, NRMAX, NGR, STR, MODE, IND)

!       CALL TXWPGR

    CASE(3)
       STR = '@u$-ir$=@'
       CALL TXGRFRX(0, GX, GY(0,0,6), NRMAX, NGR, STR, MODE, IND)

       STR = '@u$-i$#q$#$=@'
       CALL APPROPGY(MODEG, GY(0,0,7), GYL, STR, NRMAX, NGR, gDIV(7))
       CALL TXGRFRX(1, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       STR = '@u$-i$#f$#$=@'
       CALL APPROPGY(MODEG, GY(0,0,8), GYL, STR, NRMAX, NGR, gDIV(8))
       CALL TXGRFRX(2, GX, GYL, NRMAX, NGR, STR, MODE, IND)!, GYMAX=200.0, GYMIN=-100.0)

       STR = '@uhat$-i$#q$#$=@'
       CALL APPROPGY(MODEG, GY(0,0,155), GYL, STR, NRMAX, NGR, gDIV(155))
       CALL TXGRFRX(3, GX, GYL, NRMAX, NGR, STR, MODE, IND)

!       CALL TXWPGR

    CASE(4)
       STR = '@q@'
       CALL TXGRFRX(0,GX,GY(0,0,20),NRMAX,NGR,STR,MODE,IND)

       STR = '@E$-tor$=@'
       CALL TXGRFRX(1,GX,GY(0,0,11),NRMAX,NGR,STR,MODE,IND)

       STR = '@j$-$#f$#$=@'
       CALL APPROPGY(MODEG, GY(0,0,21), GYL, STR, NRMAX, NGR, gDIV(21))
       CALL TXGRFRX(2,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@B$-$#q$#$=@'
       CALL TXGRFRX(3,GX,GY(0,0,10),NRMAX,NGR,STR,MODE,IND)

!       CALL TXWPGR

    CASE(5)
       STR = '@n$-b$=@'
       CALL APPROPGY(MODEG, GY(0,0,12), GYL, STR, NRMAX, NGR, gDIV(12))
       CALL TXGRFRX(0,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@u$-b$#f$#$=@'
       CALL APPROPGY(MODEG, GY(0,0,13), GYL, STR, NRMAX, NGR, gDIV(13))
       CALL TXGRFRX(1,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@Ripple n$-b$=@'
       CALL APPROPGY(MODEG, GY(0,0,110), GYL, STR, NRMAX, NGR, gDIV(110))
       CALL TXGRFRX(2,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@Bu$-b//$=@'
       CALL APPROPGY(MODEG, GY(0,0,19), GYL, STR, NRMAX, NGR, gDIV(19))
       CALL TXGRFRX(3,GX,GYL,NRMAX,NGR,STR,MODE,IND)

!       STR = '@Ripple n$-b$=/n$-b$=@'
!!       STR = '@PNbrpLV@'
!       CALL APPROPGY(MODEG, GY(0,0,115), GYL, STR, NRMAX, NGR, gDIV(115))
!       CALL TXGRFRX(3,GX,GYL,NRMAX,NGR,STR,MODE,IND)

!       CALL TXWPGR

    CASE(6)
       STR = '@j$-e$#f$#$=@'
       CALL APPROPGY(MODEG, GY(0,0,22), GYL, STR, NRMAX, NGR, gDIV(22))
       CALL TXGRFRX(0,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@j$-i$#f$#$=@'
       CALL APPROPGY(MODEG, GY(0,0,23), GYL, STR, NRMAX, NGR, gDIV(23))
       CALL TXGRFRX(1,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@j$-b$#f$#$=@'
       CALL APPROPGY(MODEG, GY(0,0,24), GYL, STR, NRMAX, NGR, gDIV(24))
       CALL TXGRFRX(2,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@j$-r$=@'
       CALL TXGRFRX(3,GX,GY(0,0,117),NRMAX,NGR,STR,MODE,IND)

!!$       write(6,*) 'Er=', GY(21,NGR,9)
!!$       write(6,*) 'jr=', GY(21,NGR,117)

!       CALL TXWPGR

    CASE(7)
       STR = '@Bu$-e//$=@'
       CALL APPROPGY(MODEG, GY(0,0,26), GYL, STR, NRMAX, NGR, gDIV(26))
       CALL TXGRFRX(0,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@u$-e$#$/136$#$=/B@'
       CALL APPROPGY(MODEG, GY(0,0,25), GYL, STR, NRMAX, NGR, gDIV(25))
       CALL TXGRFRX(1,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@Bq$-e//$=@'
       CALL APPROPGY(MODEG, GY(0,0,158), GYL, STR, NRMAX, NGR, gDIV(159))
       CALL TXGRFRX(2,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@Bj$-//$=@'
       CALL APPROPGY(MODEG, GY(0,0,121), GYL, STR, NRMAX, NGR, gDIV(121))
       CALL TXGRFRX(3,GX,GYL,NRMAX,NGR,STR,MODE,IND)

!       CALL TXWPGR

    CASE(8)
       STR = '@Bu$-i//$=@'
       CALL APPROPGY(MODEG, GY(0,0,28), GYL, STR, NRMAX, NGR, gDIV(28))
       CALL TXGRFRX(0,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@u$-i$#$/136$#$=/B@'
       CALL APPROPGY(MODEG, GY(0,0,27), GYL, STR, NRMAX, NGR, gDIV(27))
       CALL TXGRFRX(1,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@Bq$-i//$=@'
       CALL APPROPGY(MODEG, GY(0,0,159), GYL, STR, NRMAX, NGR, gDIV(159))
       CALL TXGRFRX(2,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@BE$-//$=@'
       CALL APPROPGY(MODEG, GY(0,0,122), GYL, STR, NRMAX, NGR, gDIV(122))
       CALL TXGRFRX(3,GX,GYL,NRMAX,NGR,STR,MODE,IND)

!       CALL TXWPGR

    CASE(9)
       STR = '@D$-i,eff$=@'
!       CALL TXGRFRX(0,GX,GY(0,0,83),NRMAX,NGR,STR,MODE,IND,GYMAX=1.0)
       CALL TXGRFRX(0,GX,GY(0,0,83),NRMAX,NGR,STR,MODE,IND)

       STR = '@D$-i$=+D$-e$=@'
       CALL TXGRFRX(1,GX,GY(0,0,29),NRMAX,NGR,STR,MODE,IND)

       STR = '@$#c$#$-tbe$=@'
       CALL TXGRFRX(2,GX,GY(0,0,70),NRMAX,NGR,STR,MODE,IND)

       STR = '@$#c$#$-NCi$=@'
       CALL TXGRFRX(3,GX,GY(0,0,130),NRMAX,NGR,STR,MODE,IND)

!       CALL TXWPGR

    CASE(10)
       STR = '@T$-e$=@'
       CALL TXGRFRX(0,GX,GY(0,0,14),NRMAX,NGR,STR,MODE,IND)

       STR = '@T$-i$=@'
       CALL TXGRFRX(1,GX,GY(0,0,15),NRMAX,NGR,STR,MODE,IND)!,GYMAX=6.0)

       STR = '@p$-e$=@'
       CALL APPROPGY(MODEG, GY(0,0,102), GYL, STR, NRMAX, NGR, gDIV(102))
       CALL TXGRFRX(2,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@p$-i$=@'
       CALL APPROPGY(MODEG, GY(0,0,103), GYL, STR, NRMAX, NGR, gDIV(103))
       CALL TXGRFRX(3,GX,GYL,NRMAX,NGR,STR,MODE,IND)

!       CALL TXWPGR

    CASE(11)
       STR = '@s@'
       CALL TXGRFRX(0,GX,GY(0,0,32),NRMAX,NGR,STR,MODE,IND)

       STR = '@$#a$#@'
       CALL TXGRFRX(1,GX,GY(0,0,33),NRMAX,NGR,STR,MODE,IND)

       STR = '@$#k$#@'
       CALL TXGRFRX(2,GX,GY(0,0,34),NRMAX,NGR,STR,MODE,IND)

       STR = '@F$-CDBM$=@'
       CALL TXGRFRX(3,GX,GY(0,0,31),NRMAX,NGR,STR,MODE,IND)

!       CALL TXWPGR

    CASE(12)
       STR = '@BE$-pol$=@'
       CALL APPROPGY(MODEG, GY(0,0,17), GYL, STR, NRMAX, NGR, gDIV(17))
       CALL TXGRFRX(0,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@I@'
       CALL TXGRFRX(1,GX,GY(0,0,18),NRMAX,NGR,STR,MODE,IND)

       STR = '@U$-g$=V@'
       CALL TXGRFRX(2,GX,GY(0,0,191),NRMAX,NGR,STR,MODE,IND)

       STR = '@qhat$+2$=@'
       CALL TXGRFRX(3,GX,GY(0,0,192),NRMAX,NGR,STR,MODE,IND)

!       CALL TXWPGR

    CASE(13)
       STR = '@SLOW N$-0$=@'
       CALL APPROPGY(MODEG, GY(0,0,35), GYL, STR, NRMAX, NGR, gDIV(35))
       CALL TXGRFRX(0,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@THERMAL N$-0$=@'
       CALL APPROPGY(MODEG, GY(0,0,36), GYL, STR, NRMAX, NGR, gDIV(36))
       CALL TXGRFRX(1,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@HALO N$-0$=@'
       CALL APPROPGY(MODEG, GY(0,0,138), GYL, STR, NRMAX, NGR, gDIV(138))
       CALL TXGRFRX(2,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@TOTAL N$-0$=@'
!       CALL APPROPGY(MODEG, GY(0,0,16), GYL, STR, NRMAX, NGR, gDIV(16))
!       CALL TXGRFRX(3,GX,GYL,NRMAX,NGR,STR,MODE,IND)
       DO NG = 0, NGR
          DO NR = 0, NRMAX
             GYL(NR,NG) = GLOG(DBLE(GY(NR,NG,16)),1.D0,1.D23)
          END DO
       END DO
       CALL TXGRFRX(3,GX,GYL,NRMAX,NGR,STR,MODE,IND,ILOGIN=1)

!       CALL TXWPGR

    CASE(14)
       STR = '@SCX@'
       DO NG = 0, NGR
          DO NR = 0, NRMAX
             GYL(NR,NG) = GLOG(DBLE(GY(NR,NG,133)),1.D10,1.D23)
          END DO
       END DO
       CALL TXGRFRX(0,GX,GYL,NRMAX,NGR,STR,MODE,IND,GYMIN=18.0,ILOGIN=1)

       STR = '@SIE@'
       DO NG = 0, NGR
          DO NR = 0, NRMAX
             GYL(NR,NG) = GLOG(DBLE(GY(NR,NG,90)),1.D10,1.D23)
          END DO
       END DO
       CALL TXGRFRX(1,GX,GYL,NRMAX,NGR,STR,MODE,IND,GYMIN=18.0,ILOGIN=1)

       STR = '@PIE@'
       CALL TXGRFRX(2,GX,GY(0,0,88),NRMAX,NGR,STR,MODE,IND)

       STR = '@PCX@'
       CALL APPROPGY(MODEG, GY(0,0,89), GYL, STR, NRMAX, NGR, gDIV(89))
       CALL TXGRFRX(3,GX,GYL,NRMAX,NGR,STR,MODE,IND)

!       CALL TXWPGR

    CASE(15)
       STR = '@PBr@'
       CALL TXGRFRX(0,GX,GY(0,0,91),NRMAX,NGR,STR,MODE,IND)

       STR = '@POHe@'
       CALL APPROPGY(MODEG, GY(0,0,98), GYL, STR, NRMAX, NGR, gDIV(98))
       CALL TXGRFRX(1,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@PEQe, -PEQi@'
       CALL APPROPGY(MODEG, GY(0,0,100), GYL, STR, NRMAX, NGR, gDIV(100))
       CALL TXGRFRX(2,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@POHi@'
       CALL TXGRFRX(3,GX,GY(0,0,99),NRMAX,NGR,STR,MODE,IND)

    CASE(16)
       STR = '@SNB@'
       CALL TXGRFRX(0,GX,GY(0,0,84),NRMAX,NGR,STR,MODE,IND)

       STR = '@PNBe@'
       CALL APPROPGY(MODEG, GY(0,0,96), GYL, STR, NRMAX, NGR, gDIV(96))
       CALL TXGRFRX(1,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@NBI deposition@'
       CALL APPROPGY(MODEG, GY(0,0,104), GYL, STR, NRMAX, NGR, gDIV(104))
       CALL TXGRFRX(2,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@PNBi@'
       CALL APPROPGY(MODEG, GY(0,0,97), GYL, STR, NRMAX, NGR, gDIV(97))
       CALL TXGRFRX(3,GX,GYL,NRMAX,NGR,STR,MODE,IND)

!       CALL TXWPGR

    CASE(17)
       STR = '@SNB perp@'
       CALL APPROPGY(MODEG, GY(0,0,105), GYL, STR, NRMAX, NGR, gDIV(105))
       CALL TXGRFRX(0,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@SNB tang@'
       CALL APPROPGY(MODEG, GY(0,0,106), GYL, STR, NRMAX, NGR, gDIV(106))
       CALL TXGRFRX(1,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@SNB perp ion@'
       CALL APPROPGY(MODEG, GY(0,0,107), GYL, STR, NRMAX, NGR, gDIV(107))
       CALL TXGRFRX(2,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@SNB tang ion@'
       CALL APPROPGY(MODEG, GY(0,0,108), GYL, STR, NRMAX, NGR, gDIV(108))
       CALL TXGRFRX(3,GX,GYL,NRMAX,NGR,STR,MODE,IND)

    CASE(18)
       STR = '@Generated toroidal torque density@'
       CALL APPROPGY(MODEG, GY(0,0,120), GYL, STR, NRMAX, NGR, gDIV(120))
       CALL TXGRFRX(0,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@NBI tor. coll. torque density@'
       CALL APPROPGY(MODEG, GY(0,0,119), GYL, STR, NRMAX, NGR, gDIV(119))
       CALL TXGRFRX(1,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@PALFe@'
       CALL APPROPGY(MODEG, GY(0,0,123), GYL, STR, NRMAX, NGR, gDIV(123))
       CALL TXGRFRX(2,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@PALFi@'
       CALL APPROPGY(MODEG, GY(0,0,124), GYL, STR, NRMAX, NGR, gDIV(124))
       CALL TXGRFRX(3,GX,GYL,NRMAX,NGR,STR,MODE,IND)

!       CALL TXWPGR

    CASE(19)
       DO NR = 0, NRMAX
          IF(GX(NR) >= 0.95) EXIT
       END DO
       NR = NR - 1

       NRMAXL = NRMAX-NR
       allocate(GYL2(0:NRMAXL,0:NGR))

       STR = '@n$-e$=@'
       CALL APPROPGY(MODEG, GY(0,0,1), GYL, STR, NRMAX, NGR, gDIV(1))
       GYL2(0:NRMAXL,0:NGR) = GYL(NR:NRMAX,0:NGR)
       GYMAX = MAXVAL(GYL2(0:NRMAXL,0:NGR))
       CALL TXGRFRX(0, GX(NR:NRMAX), GYL2, NRMAXL, NGR, STR, MODE, IND, &
            &       GX(NR), GYMAX)

       STR = '@T$-e$=@'
       CALL APPROPGY(MODEG, GY(0,0,14), GYL, STR, NRMAX, NGR, gDIV(14))
       GYL2(0:NRMAXL,0:NGR) = GYL(NR:NRMAX,0:NGR)
       GYMAX = MAXVAL(GYL2(0:NRMAXL,0:NGR))
       CALL TXGRFRX(1, GX(NR:NRMAX), GYL2, NRMAXL, NGR, STR, MODE, IND, &
            &       GX(NR), GYMAX)

       STR = '@n$-i$=@'
       CALL APPROPGY(MODEG, GY(0,0,37), GYL, STR, NRMAX, NGR, gDIV(37))
       GYL2(0:NRMAXL,0:NGR) = GYL(NR:NRMAX,0:NGR)
       GYMAX = MAXVAL(GYL2(0:NRMAXL,0:NGR))
       CALL TXGRFRX(2, GX(NR:NRMAX), GYL2, NRMAXL, NGR, STR, MODE, IND, &
            &       GX(NR), GYMAX)

       STR = '@T$-i$=@'
       CALL APPROPGY(MODEG, GY(0,0,15), GYL, STR, NRMAX, NGR, gDIV(15))
       GYL2(0:NRMAXL,0:NGR) = GYL(NR:NRMAX,0:NGR)
       GYMAX = MAXVAL(GYL2(0:NRMAXL,0:NGR))
       CALL TXGRFRX(3, GX(NR:NRMAX), GYL2, NRMAXL, NGR, STR, MODE, IND, &
            &       GX(NR), GYMAX)

       deallocate(GYL2)

    CASE(20)
       STR = '@u$-e$=$+V$=@'
       CALL TXGRFRX(0,GX,GY(0,0,156),NRMAX,NGR,STR,MODE,IND)

       STR = '@u$-i$=$+V$=@'
       CALL TXGRFRX(1,GX,GY(0,0,157),NRMAX,NGR,STR,MODE,IND)

!       STR = '@$#G$#$-e$=.$#$/321r$#@'
       STR = '@$#G$#$-e$=.grad $#r$#@'
       CALL APPROPGY(MODEG, GY(0,0,136), GYL, STR, NRMAX, NGR, gDIV(136))
!       CALL TXGRFRX(2,GX,GYL,NRMAX,NGR,STR,MODE,IND,GYMIN=0.0,GYMAX=0.05)
       CALL TXGRFRX(2,GX,GYL,NRMAX,NGR,STR,MODE,IND)

!       STR = '@$#G$#$-i$=.$#$/321r$#@'
       STR = '@$#G$#$-i$=.grad $#r$#@'
       CALL APPROPGY(MODEG, GY(0,0,137), GYL, STR, NRMAX, NGR, gDIV(137))
       CALL TXGRFRX(3,GX,GYL,NRMAX,NGR,STR,MODE,IND)

    CASE(21)
       STR = '@$#S$#e$-s$=$#G$#$-s$=.grad V@'
       CALL TXGRFRX(0,GX,GY(0,0,160),NRMAX,NGR,STR,MODE,IND)

       STR = '@Additional Torque@'
       CALL TXGRFRX(1,GX,GY(0,0,131),NRMAX,NGR,STR,MODE,IND)

       STR = '@$#S$#e$-s$=$#G$#$-s$=.grad V w/ beam@'
       CALL TXGRFRX(2,GX,GY(0,0,199),NRMAX,NGR,STR,MODE,IND)

       STR = '@$#G$#$-b$=.grad $#r$#@'
       CALL APPROPGY(MODEG, GY(0,0,200), GYL, STR, NRMAX, NGR, gDIV(200))
       CALL TXGRFRX(3,GX,GYL,NRMAX,NGR,STR,MODE,IND)

!       STR = '@F$-CDIM$=@'
!       CALL TXGRFRX(3,GX,GY(0,0,140),NRMAX,NGR,STR,MODE,IND)

!       CALL TXWPGR

    CASE(22)
       STR = '@BVediag 1@'
       CALL APPROPGY(MODEG, GY(0,0,142), GYL, STR, NRMAX, NGR, gDIV(142))
       CALL TXGRFRX(0,GX,GYL        ,NRMAX,NGR,STR,MODE,IND)

       STR = '@BVidiag 1@'
       CALL APPROPGY(MODEG, GY(0,0,143), GYL, STR, NRMAX, NGR, gDIV(143))
       CALL TXGRFRX(1,GX,GYL        ,NRMAX,NGR,STR,MODE,IND)

       STR = '@BVediag 2@'
       CALL APPROPGY(MODEG, GY(0,0,144), GYL, STR, NRMAX, NGR, gDIV(144))
       CALL TXGRFRX(2,GX,GYL        ,NRMAX,NGR,STR,MODE,IND)

       STR = '@BVidiag 2@'
       CALL APPROPGY(MODEG, GY(0,0,145), GYL, STR, NRMAX, NGR, gDIV(145))
       CALL TXGRFRX(3,GX,GYL        ,NRMAX,NGR,STR,MODE,IND)

    CASE(23)
!       STR = '@VWpch@'
!       CALL TXGRFRX(0,GX,GY(0,0,132),NRMAX,NGR,STR,MODE,IND)

       STR = '@rMue@'
       CALL TXGRFRX( 0,GX,GY(0,0,74),NRMAX,NGR,STR,MODE,IND)

       STR = '@rMui@'
       CALL TXGRFRX( 1,GX,GY(0,0,75),NRMAX,NGR,STR,MODE,IND)

       STR = '@Vmpe@'
       CALL TXGRFRX( 2,GX,GY(0,0,189),NRMAX,NGR,STR,MODE,IND)

       STR = '@Vmpi@'
       CALL TXGRFRX( 3,GX,GY(0,0,190),NRMAX,NGR,STR,MODE,IND)

    CASE(24)
       STR = '@PRFe@'
       CALL APPROPGY(MODEG, GY(0,0,85), GYL, STR, NRMAX, NGR, gDIV(85))
       CALL TXGRFRX(0,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@PRFi@'
       CALL APPROPGY(MODEG, GY(0,0,86), GYL, STR, NRMAX, NGR, gDIV(86))
       CALL TXGRFRX(1,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@$#w$#$-ExB$=@'
       CALL APPROPGY(MODEG, GY(0,0,141), GYL, STR, NRMAX, NGR, gDIV(141))
       CALL TXGRFRX(2,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@$#c$#$-i$=@'
       CALL TXGRFRX(3,GX,GY(0,0,135),NRMAX,NGR,STR,MODE,IND)

    CASE(-1)
       STR = '@E$-r$=@'
       CALL APPROPGY(MODEG, GY(0,0,9), GYL, STR, NRMAX, NGR, gDIV(9))
       CALL TXGRFRXS(0, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       STR = '@BE$-pol$=@'
       CALL APPROPGY(MODEG, GY(0,0,17), GYL, STR, NRMAX, NGR, gDIV(17))
       CALL TXGRFRXS(1,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@E$-tor$=@'
       CALL TXGRFRXS(2,GX,GY(0,0,11),NRMAX,NGR,STR,MODE,IND)

       STR = '@B$-$#q$#$=@'
       CALL TXGRFRXS(3,GX,GY(0,0,10),NRMAX  ,NGR,STR,MODE,IND)

       STR = '@Z*n$-i$=+Z*n$-b$=-n$-e$=@'
       CALL APPROPGY(MODEG, GY(0,0,2), GYL, STR, NRMAX, NGR, gDIV(2))
       CALL TXGRFRXS(4, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       STR = '@T$-e$=@'
       CALL TXGRFRXS(5,GX,GY(0,0,14),NRMAX,NGR,STR,MODE,IND)

       STR = '@T$-i$=@'
       CALL TXGRFRXS(6,GX,GY(0,0,15),NRMAX,NGR,STR,MODE,IND)

       STR = '@I@'
       CALL TXGRFRXS(7,GX,GY(0,0,18),NRMAX,NGR,STR,MODE,IND)

       STR = '@n$-e$=@'
       CALL APPROPGY(MODEG, GY(0,0,1), GYL, STR, NRMAX, NGR, gDIV(1))
       CALL TXGRFRXS(8,GX, GYL, NRMAX, NGR, STR, MODE, IND)

       STR = '@u$-er$=@'
       CALL TXGRFRXS(9, GX, GY(0,0,3), NRMAX, NGR, STR, MODE, IND)

!!$       STR = '@u$-e$#q$#$=@'
!!$       CALL APPROPGY(MODEG, GY(0,0,4), GYL, STR, NRMAX, NGR, gDIV(4))
!!$       CALL TXGRFRXS(10, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       STR = '@Bu$-e//$=@'
       CALL APPROPGY(MODEG, GY(0,0,26), GYL, STR, NRMAX, NGR, gDIV(26))
       CALL TXGRFRXS(10, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       STR = '@u$-e$#f$#$=@'
       CALL APPROPGY(MODEG, GY(0,0,5), GYL, STR, NRMAX, NGR, gDIV(5))
       CALL TXGRFRXS(11, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       STR = '@n$-i$=@'
       CALL APPROPGY(MODEG, GY(0,0,37), GYL, STR, NRMAX, NGR, gDIV(37))
       CALL TXGRFRXS(12,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@u$-ir$=@'
       CALL TXGRFRXS(13, GX, GY(0,0,6), NRMAX, NGR, STR, MODE, IND)

!!$       STR = '@u$-i$#q$#$=@'
!!$       CALL APPROPGY(MODEG, GY(0,0,7), GYL, STR, NRMAX, NGR, gDIV(7))
!!$       CALL TXGRFRXS(14, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       STR = '@Bu$-i//$=@'
       CALL APPROPGY(MODEG, GY(0,0,28), GYL, STR, NRMAX, NGR, gDIV(28))
       CALL TXGRFRXS(14, GX, GYL, NRMAX, NGR, STR, MODE, IND)

       STR = '@u$-i$#f$#$=@'
       CALL APPROPGY(MODEG, GY(0,0,8), GYL, STR, NRMAX, NGR, gDIV(8))
       CALL TXGRFRXS(15, GX, GYL, NRMAX, NGR, STR, MODE, IND)

    CASE(-2)
       STR = '@n$-b$=@'
       CALL APPROPGY(MODEG, GY(0,0,12), GYL, STR, NRMAX, NGR, gDIV(12))
       CALL TXGRFRXS(0,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@Bu$-b//$=@'
       CALL APPROPGY(MODEG, GY(0,0,19), GYL, STR, NRMAX, NGR, gDIV(13))
       CALL TXGRFRXS(1,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@u$-b$#f$#$=@'
       CALL APPROPGY(MODEG, GY(0,0,13), GYL, STR, NRMAX, NGR, gDIV(13))
       CALL TXGRFRXS(2,GX,GYL,NRMAX  ,NGR,STR,MODE,IND)     

       STR = '@Ripple n$-b$=@'
       CALL APPROPGY(MODEG, GY(0,0,109), GYL, STR, NRMAX, NGR, gDIV(109))
       CALL TXGRFRXS(3,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@TOTAL N$-0$=@'
!       CALL APPROPGY(MODEG, GY(0,0,16), GYL, STR, NRMAX, NGR, gDIV(16))
!       CALL TXGRFRXS(4,GX,GYL,NRMAX,NGR,STR,MODE,IND)
       DO NG = 0, NGR
          DO NR = 0, NRMAX
             GYL(NR,NG) = GLOG(DBLE(GY(NR,NG,16)),1.D0,1.D23)
          END DO
       END DO
       CALL TXGRFRXS(4,GX,GYL,NRMAX,NGR,STR,MODE,IND,ILOGIN=1)

       STR = '@SLOW N$-0$=@'
       CALL APPROPGY(MODEG, GY(0,0,35), GYL, STR, NRMAX, NGR, gDIV(35))
       CALL TXGRFRXS(5,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@THERMAL N$-0$=@'
       CALL APPROPGY(MODEG, GY(0,0,36), GYL, STR, NRMAX, NGR, gDIV(36))
       CALL TXGRFRXS(6,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@HALO N$-0$=@'
       CALL APPROPGY(MODEG, GY(0,0,138), GYL, STR, NRMAX, NGR, gDIV(138))
       CALL TXGRFRXS(7,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@q@'
       CALL TXGRFRXS(8,GX,GY(0,0,20),NRMAX,NGR,STR,MODE,IND)

       STR = '@j$-$#f$#$=@'
       CALL APPROPGY(MODEG, GY(0,0,21), GYL, STR, NRMAX, NGR, gDIV(21))
       CALL TXGRFRXS(9,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@j$-e$#f$#$=@'
       CALL APPROPGY(MODEG, GY(0,0,22), GYL, STR, NRMAX, NGR, gDIV(22))
       CALL TXGRFRXS(10,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@j$-i$#f$#$=@'
       CALL APPROPGY(MODEG, GY(0,0,23), GYL, STR, NRMAX, NGR, gDIV(23))
       CALL TXGRFRXS(11,GX,GYL,NRMAX,NGR,STR,MODE,IND)

       STR = '@S@'
       CALL TXGRFRXS(12,GX,GY(0,0,32),NRMAX,NGR,STR,MODE,IND)

       STR = '@$#a$#@'
       CALL TXGRFRXS(13,GX,GY(0,0,33),NRMAX,NGR,STR,MODE,IND)

       STR = '@F$-CDBM$=@'
       CALL TXGRFRXS(14,GX,GY(0,0,31),NRMAX,NGR,STR,MODE,IND)     

       STR = '@G$-1$=h$+2$=@'
       CALL TXGRFRXS(15,GX,GY(0,0,30),NRMAX,NGR,STR,MODE,IND)

    CASE(-3)
       STR = '@LQm1@'
       CALL APPROPGY(MODEG, GY(0,0,38), GYL, STR, NRMAX, NGR, gDIV(38))
       CALL TXGRFRXS( 0,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@LQm2@'
       CALL APPROPGY(MODEG, GY(0,0,39), GYL, STR, NRMAX, NGR, gDIV(39))
       CALL TXGRFRXS( 1,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@LQm3@'
       CALL TXGRFRXS( 2,GX,GY(0,0,40),NRMAX,NGR,STR,MODE,IND)

       STR = '@LQm4@'
       CALL APPROPGY(MODEG, GY(0,0,41), GYL, STR, NRMAX, NGR, gDIV(41))
       CALL TXGRFRXS( 3,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@LQm5@'
       CALL APPROPGY(MODEG, GY(0,0,42), GYL, STR, NRMAX, NGR, gDIV(42))
       CALL TXGRFRXS( 4,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@LQe1@'
       CALL TXGRFRXS( 8,GX,GY(0,0,43),NRMAX,NGR,STR,MODE,IND)

       STR = '@LQe2@'
       CALL TXGRFRXS( 9,GX,GY(0,0,44),NRMAX,NGR,STR,MODE,IND)

       STR = '@LQe3@'
       CALL APPROPGY(MODEG, GY(0,0,45), GYL, STR, NRMAX, NGR, gDIV(45))
       CALL TXGRFRXS(10,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@LQe4@'
       CALL APPROPGY(MODEG, GY(0,0,46), GYL, STR, NRMAX, NGR, gDIV(46))
       CALL TXGRFRXS(11,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@LQi1@'
       CALL TXGRFRXS(12,GX,GY(0,0,50),NRMAX,NGR,STR,MODE,IND)

       STR = '@LQi2@'
       CALL TXGRFRXS(13,GX,GY(0,0,51),NRMAX,NGR,STR,MODE,IND)

       STR = '@LQi3@'
       CALL APPROPGY(MODEG, GY(0,0,52), GYL, STR, NRMAX, NGR, gDIV(52))
       CALL TXGRFRXS(14,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@LQi4@'
       CALL APPROPGY(MODEG, GY(0,0,53), GYL, STR, NRMAX, NGR, gDIV(53))
       CALL TXGRFRXS(15,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

    CASE(-4)
       STR = '@LQe5@'
       CALL TXGRFRXS( 0,GX,GY(0,0,47),NRMAX,NGR,STR,MODE,IND)

       STR = '@LQe6@'
       CALL APPROPGY(MODEG, GY(0,0,48), GYL, STR, NRMAX, NGR, gDIV(48))
       CALL TXGRFRXS( 1,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@LQe7@'
       CALL APPROPGY(MODEG, GY(0,0,49), GYL, STR, NRMAX, NGR, gDIV(49))
       CALL TXGRFRXS( 2,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@LQn1@'
       CALL APPROPGY(MODEG, GY(0,0,59), GYL, STR, NRMAX, NGR, gDIV(59))
       CALL TXGRFRXS( 3,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)     

       STR = '@LQi5@'
       CALL TXGRFRXS( 4,GX,GY(0,0,54),NRMAX,NGR,STR,MODE,IND)

       STR = '@LQi6@'
       CALL APPROPGY(MODEG, GY(0,0,55), GYL, STR, NRMAX, NGR, gDIV(55))
       CALL TXGRFRXS( 5,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@LQi7@'
       CALL APPROPGY(MODEG, GY(0,0,56), GYL, STR, NRMAX, NGR, gDIV(56))
       CALL TXGRFRXS( 6,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@LQn2@'
       CALL APPROPGY(MODEG, GY(0,0,60), GYL, STR, NRMAX, NGR, gDIV(60))
       CALL TXGRFRXS( 7,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@LQb1@'
       CALL APPROPGY(MODEG, GY(0,0,57), GYL, STR, NRMAX, NGR, gDIV(57))
       CALL TXGRFRXS( 8,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@LQb2@'
       CALL APPROPGY(MODEG, GY(0,0,198), GYL, STR, NRMAX, NGR, gDIV(198))
       CALL TXGRFRXS( 9,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@LQb3@'
       CALL APPROPGY(MODEG, GY(0,0,58), GYL, STR, NRMAX, NGR, gDIV(58))
       CALL TXGRFRXS(10,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@LQn3@'
       CALL APPROPGY(MODEG, GY(0,0,61), GYL, STR, NRMAX, NGR, gDIV(61))
       CALL TXGRFRXS(11,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@LQb4@'
       CALL APPROPGY(MODEG, GY(0,0,196), GYL, STR, NRMAX, NGR, gDIV(196))
       CALL TXGRFRXS(12,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@LQb7@'
       CALL APPROPGY(MODEG, GY(0,0,197), GYL, STR, NRMAX, NGR, gDIV(197))
       CALL TXGRFRXS(13,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@LQr1@'
       CALL APPROPGY(MODEG, GY(0,0,62), GYL, STR, NRMAX, NGR, gDIV(62))
       CALL TXGRFRXS(14,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

    CASE(-5)
       STR = '@rMue@'
       CALL TXGRFRXS( 0,GX,GY(0,0,74),NRMAX,NGR,STR,MODE,IND)

       STR = '@rMui@'
       CALL TXGRFRXS( 1,GX,GY(0,0,75),NRMAX,NGR,STR,MODE,IND)

       STR = '@$#c$#$-e$=@'
       CALL TXGRFRXS( 2,GX,GY(0,0,70),NRMAX,NGR,STR,MODE,IND)

       STR = '@$#c$#$-i$=@'
       CALL TXGRFRXS( 3,GX,GY(0,0,71),NRMAX,NGR,STR,MODE,IND)

       STR = '@rNuL@'
       CALL TXGRFRXS( 4,GX,GY(0,0,63),NRMAX,NGR,STR,MODE,IND)

       STR = '@rNuLTe@'
       CALL APPROPGY(MODEG, GY(0,0,92), GYL, STR, NRMAX, NGR, gDIV(92))
       CALL TXGRFRXS( 5,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@rNuLTi@'
       CALL TXGRFRXS( 6,GX,GY(0,0,93),NRMAX,NGR,STR,MODE,IND)

       STR = '@rNuOL@'
       CALL TXGRFRXS( 7,GX,GY(0,0,82),NRMAX,NGR,STR,MODE,IND)

       STR = '@rNuiCX@'
       CALL TXGRFRXS( 8,GX,GY(0,0,78),NRMAX,NGR,STR,MODE,IND)

       STR = '@rNuION@'
       CALL TXGRFRXS( 9,GX,GY(0,0,80),NRMAX,NGR,STR,MODE,IND)     

       STR = '@rNu0e@'
       CALL TXGRFRXS(10,GX,GY(0,0,68),NRMAX,NGR,STR,MODE,IND)

       STR = '@rNu0i@'
       CALL TXGRFRXS(11,GX,GY(0,0,69),NRMAX,NGR,STR,MODE,IND)

       STR = '@rNuTei@'
       CALL TXGRFRXS(12,GX,GY(0,0,67),NRMAX,NGR,STR,MODE,IND)

       STR = '@D01@'
       CALL APPROPGY(MODEG, GY(0,0,72), GYL, STR, NRMAX, NGR, gDIV(72))
       CALL TXGRFRXS(13,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@D02@'
       CALL APPROPGY(MODEG, GY(0,0,73), GYL, STR, NRMAX, NGR, gDIV(73))
       CALL TXGRFRXS(14,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@D03@'
       CALL APPROPGY(MODEG, GY(0,0,139), GYL, STR, NRMAX, NGR, gDIV(139))
       CALL TXGRFRXS(15,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

    CASE(-6)
       STR = '@xmu$-e11$=@'
       CALL APPROPGY(MODEG, GY(0,0,146), GYL, STR, NRMAX, NGR, gDIV(146))
       CALL TXGRFRXS( 0,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@xmu$-e12$=@'
       CALL APPROPGY(MODEG, GY(0,0,147), GYL, STR, NRMAX, NGR, gDIV(147))
       CALL TXGRFRXS( 1,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@xmu$-e21$=@'
       CALL APPROPGY(MODEG, GY(0,0,148), GYL, STR, NRMAX, NGR, gDIV(148))
       CALL TXGRFRXS( 2,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@xmu$-e22$=@'
       CALL APPROPGY(MODEG, GY(0,0,149), GYL, STR, NRMAX, NGR, gDIV(149))
       CALL TXGRFRXS( 3,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@xmu$-i11$=@'
       CALL APPROPGY(MODEG, GY(0,0,150), GYL, STR, NRMAX, NGR, gDIV(150))
       CALL TXGRFRXS( 4,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@xmu$-i12$=@'
       CALL APPROPGY(MODEG, GY(0,0,151), GYL, STR, NRMAX, NGR, gDIV(151))
       CALL TXGRFRXS( 5,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@xmu$-i21$=@'
       CALL APPROPGY(MODEG, GY(0,0,152), GYL, STR, NRMAX, NGR, gDIV(152))
       CALL TXGRFRXS( 6,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@xmu$-i22$=@'
       CALL APPROPGY(MODEG, GY(0,0,153), GYL, STR, NRMAX, NGR, gDIV(153))
       CALL TXGRFRXS( 7,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

!!$       STR = '@FWthe@'
!!$       CALL APPROPGY(MODEG, GY(0,0,64), GYL, STR, NRMAX, NGR, gDIV(64))
!!$       CALL TXGRFRXS( 8,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@WPM@'
       CALL TXGRFRXS( 9,GX,GY(0,0,66),NRMAX,NGR,STR,MODE,IND)

       STR = '@rNu0b@'
       CALL TXGRFRXS(10,GX,GY(0,0,87),NRMAX,NGR,STR,MODE,IND)

       STR = '@Vbpara@'
       CALL APPROPGY(MODEG, GY(0,0,109), GYL, STR, NRMAX, NGR, gDIV(109))
       CALL TXGRFRXS(11,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@rNuB@'
       CALL TXGRFRXS(12,GX,GY(0,0,79),NRMAX,NGR,STR,MODE,IND)

       STR = '@rNubL@'
       CALL TXGRFRXS(13,GX,GY(0,0,118),NRMAX,NGR,STR,MODE,IND)

       STR = '@rNuLB@'
       CALL TXGRFRXS(14,GX,GY(0,0,116),NRMAX,NGR,STR,MODE,IND)

!!$       STR = '@SiLC@'
!!$       CALL TXGRFRXS(11,GX,GY(0,0,81),NRMAX,NGR,STR,MODE,IND)

!!$       STR = '@rNubrp1@'
!!$       CALL TXGRFRXS( 4,GX,GY(0,0,111),NRMAX,NGR,STR,MODE,IND)
!!$
!!$       STR = '@DltRP@'
!!$       CALL TXGRFRXS( 5,GX,GY(0,0,112),NRMAX,NGR,STR,MODE,IND)

!!$       STR = '@Ubrp@'
!!$       CALL TXGRFRXS(12,GX,GY(0,0,113),NRMAX,NGR,STR,MODE,IND)
!!$
!!$       STR = '@Dbrp@'
!!$       CALL TXGRFRXS(13,GX,GY(0,0,114),NRMAX,NGR,STR,MODE,IND)

    CASE(-7)
       STR = '@l$-11$=$+ee$=@'
       CALL APPROPGY(MODEG, GY(0,0,166), GYL, STR, NRMAX, NGR, gDIV(166))
       CALL TXGRFRXS( 0,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@l$-12$=$+ee$=@'
       CALL APPROPGY(MODEG, GY(0,0,167), GYL, STR, NRMAX, NGR, gDIV(167))
       CALL TXGRFRXS( 1,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@l$-21$=$+ee$=@'
       CALL APPROPGY(MODEG, GY(0,0,168), GYL, STR, NRMAX, NGR, gDIV(168))
       CALL TXGRFRXS( 2,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@l$-22$=$+ee$=@'
       CALL APPROPGY(MODEG, GY(0,0,169), GYL, STR, NRMAX, NGR, gDIV(169))
       CALL TXGRFRXS( 3,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@l$-11$=$+ei$=@'
       CALL APPROPGY(MODEG, GY(0,0,170), GYL, STR, NRMAX, NGR, gDIV(170))
       CALL TXGRFRXS( 4,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@l$-12$=$+ei$=@'
       CALL APPROPGY(MODEG, GY(0,0,171), GYL, STR, NRMAX, NGR, gDIV(171))
       CALL TXGRFRXS( 5,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@l$-21$=$+ei$=@'
       CALL APPROPGY(MODEG, GY(0,0,172), GYL, STR, NRMAX, NGR, gDIV(172))
       CALL TXGRFRXS( 6,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@l$-22$=$+ei$=@'
       CALL APPROPGY(MODEG, GY(0,0,173), GYL, STR, NRMAX, NGR, gDIV(173))
       CALL TXGRFRXS( 7,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@l$-11$=$+ie$=@'
       CALL APPROPGY(MODEG, GY(0,0,174), GYL, STR, NRMAX, NGR, gDIV(174))
       CALL TXGRFRXS( 8,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@l$-12$=$+ie$=@'
       CALL APPROPGY(MODEG, GY(0,0,175), GYL, STR, NRMAX, NGR, gDIV(175))
       CALL TXGRFRXS( 9,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@l$-21$=$+ie$=@'
       CALL APPROPGY(MODEG, GY(0,0,176), GYL, STR, NRMAX, NGR, gDIV(176))
       CALL TXGRFRXS(10,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@l$-22$=$+ie$=@'
       CALL APPROPGY(MODEG, GY(0,0,177), GYL, STR, NRMAX, NGR, gDIV(177))
       CALL TXGRFRXS(11,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@l$-11$=$+ii$=@'
       CALL APPROPGY(MODEG, GY(0,0,178), GYL, STR, NRMAX, NGR, gDIV(178))
       CALL TXGRFRXS(12,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@l$-12$=$+ii$=@'
       CALL APPROPGY(MODEG, GY(0,0,179), GYL, STR, NRMAX, NGR, gDIV(179))
       CALL TXGRFRXS(13,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@l$-21$=$+ii$=@'
       CALL APPROPGY(MODEG, GY(0,0,180), GYL, STR, NRMAX, NGR, gDIV(180))
       CALL TXGRFRXS(14,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@l$-22$=$+ii$=@'
       CALL APPROPGY(MODEG, GY(0,0,181), GYL, STR, NRMAX, NGR, gDIV(181))
       CALL TXGRFRXS(15,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

    CASE(-8)
       STR = '@l$-11$=$+ef$=@'
       CALL APPROPGY(MODEG, GY(0,0,182), GYL, STR, NRMAX, NGR, gDIV(182))
       CALL TXGRFRXS( 0,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@l$-21$=$+ef$=@'
       CALL APPROPGY(MODEG, GY(0,0,183), GYL, STR, NRMAX, NGR, gDIV(183))
       CALL TXGRFRXS( 1,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@l$-11$=$+fe$=@'
       CALL APPROPGY(MODEG, GY(0,0,184), GYL, STR, NRMAX, NGR, gDIV(184))
       CALL TXGRFRXS( 2,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@l$-11$=$+if$=@'
       CALL APPROPGY(MODEG, GY(0,0,185), GYL, STR, NRMAX, NGR, gDIV(185))
       CALL TXGRFRXS( 3,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@l$-11$=$+fi$=@'
       CALL APPROPGY(MODEG, GY(0,0,186), GYL, STR, NRMAX, NGR, gDIV(186))
       CALL TXGRFRXS( 4,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@l$-11$=$+ff$=@'
       CALL APPROPGY(MODEG, GY(0,0,187), GYL, STR, NRMAX, NGR, gDIV(187))
       CALL TXGRFRXS( 5,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@xmu$-f$=@'
       CALL APPROPGY(MODEG, GY(0,0,188), GYL, STR, NRMAX, NGR, gDIV(188))
       CALL TXGRFRXS( 6,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

    CASE(-9)
       STR = '@Bu$-e//$=@'
       CALL APPROPGY(MODEG, GY(0,0,26), GYL, STR, NRMAX, NGR, gDIV(26))
       CALL TXGRFRX( 0,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@Bu$-i//$=@'
       CALL APPROPGY(MODEG, GY(0,0,28), GYL, STR, NRMAX, NGR, gDIV(28))
       CALL TXGRFRX( 1,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@Bu$-e//$= Neo. solver@'
       CALL APPROPGY(MODEG, GY(0,0,162), GYL, STR, NRMAX, NGR, gDIV(162))
!       CALL TXGRFRX( 2,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)
       CALL TXGRFRX( 2,GX(1:NRMAX),GYL(1:NRMAX,:),NRMAX-1,NGR,STR,MODE,IND)

       STR = '@Bu$-i//$= Neo. solver@'
       CALL APPROPGY(MODEG, GY(0,0,163), GYL, STR, NRMAX, NGR, gDIV(163))
!       CALL TXGRFRX( 3,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)
       CALL TXGRFRX( 3,GX(1:NRMAX),GYL(1:NRMAX,:),NRMAX-1,NGR,STR,MODE,IND)

    CASE(-10)
       STR = '@Bqhat$-e//$=@'
       CALL APPROPGY(MODEG, GY(0,0,158), GYL, STR, NRMAX, NGR, gDIV(158))
       CALL TXGRFRX( 0,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@Bqhat$-i//$=@'
       CALL APPROPGY(MODEG, GY(0,0,159), GYL, STR, NRMAX, NGR, gDIV(159))
       CALL TXGRFRX( 1,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@Bqhat$-e//$= Neo. solver@'
       CALL APPROPGY(MODEG, GY(0,0,164), GYL, STR, NRMAX, NGR, gDIV(164))
       CALL TXGRFRX( 2,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

       STR = '@Bqhat$-i//$= Neo. solver@'
       CALL APPROPGY(MODEG, GY(0,0,165), GYL, STR, NRMAX, NGR, gDIV(165))
       CALL TXGRFRX( 3,GX,GYL       ,NRMAX,NGR,STR,MODE,IND)

    CASE DEFAULT
       WRITE(6,*) 'Unknown NGYR: NGYR = ',NGYR
    END SELECT

    CALL PAGEE

    SELECT CASE(NGYR)
    CASE(-1)
       NGYR = -2  ; CYCLE
    CASE(-2)
       NGYR =  0  ; EXIT
    CASE(-3)
       NGYR = -4  ; CYCLE
    CASE(-4)
       NGYR =  0  ; EXIT
    CASE(-5)
       NGYR = -6  ; CYCLE
    CASE(-6)
       NGYR = -7  ; CYCLE
    CASE(-7)
       NGYR = -8  ; CYCLE
    CASE(-8)
       NGYR =  0  ; EXIT
    CASE(-9)
       NGYR = -10 ; CYCLE
    CASE(-10)
       NGYR =  0  ; EXIT
    CASE DEFAULT
       EXIT
    END SELECT

    END DO

    deallocate(GYL)

    RETURN
  END SUBROUTINE TXGRFR

  !***************************************************************
  !
  !   Animation of radial profile evolution
  !
  !***************************************************************

  SUBROUTINE TXGRFRA(NGYRIN)

    use tx_commons, only : NRMAX
    INTEGER(4), INTENT(IN) :: NGYRIN
    INTEGER(4) :: IND, NG, NGYR, IFNT, I, j
    real(4), dimension(:),     allocatable :: GMAX, GMIN
    real(4), dimension(:,:,:), allocatable :: GYL
    character(len=50), dimension(:), allocatable :: STR

    NGYR = NGYRIN

    IF (NGT <= -1) THEN
       WRITE(6,*) 'G', NGYR, ' has no data'
       RETURN
    END IF

    IF (MODEG == 2) THEN
       IND = 9
    ELSE
       IND = 0
    END IF

    allocate(GYL(0:NRMAX,0:NGT,0:15), STR(0:15), GMAX(0:15), GMIN(0:15))

    DO

    CALL PAGES
    CALL SETCHS(0.3, 0.0)
    CALL SETLIN(0, 1, 7)
    CALL INQFNT(IFNT)
    CALL SETFNT(44)

    SELECT CASE(NGYR)
    CASE(-1)
       j = 0
       STR(j) = '@E$-r$=@'
       CALL APPROPGY(MODEG, GYT(0,0,9),  GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(9), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@BE$-pol$=@'
       CALL APPROPGY(MODEG, GYT(0,0,17), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(17), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@E$-tor$=@'
       CALL APPROPGY(MODEG, GYT(0,0,11), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(11), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@B$-$#q$#$=@'
       CALL APPROPGY(MODEG, GYT(0,0,10), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(10), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@Z*n$-i$=+Z*n$-b$=-n$-e$=@'
       CALL APPROPGY(MODEG, GYT(0,0,2),  GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(2), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@T$-e$=@'
       CALL APPROPGY(MODEG, GYT(0,0,14), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(14), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@T$-i$=@'
       CALL APPROPGY(MODEG, GYT(0,0,15), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(15), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@I@'
       CALL APPROPGY(MODEG, GYT(0,0,18), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(18), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@n$-e$=@'
       CALL APPROPGY(MODEG, GYT(0,0,1),  GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(1), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@u$-er$=@'
       CALL APPROPGY(MODEG, GYT(0,0,3),  GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(3), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@u$-e$#q$#$=@'
       CALL APPROPGY(MODEG, GYT(0,0,4),  GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(4), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@u$-e$#f$#$=@'
       CALL APPROPGY(MODEG, GYT(0,0,5),  GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(5), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@n$-i$=@'
       CALL APPROPGY(MODEG, GYT(0,0,37), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(37), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@u$-ir$=@'
       CALL APPROPGY(MODEG, GYT(0,0,6),  GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(6), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@u$-i$#q$#$=@'
       CALL APPROPGY(MODEG, GYT(0,0,7),  GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(7), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@u$-i$#f$#$=@'
       CALL APPROPGY(MODEG, GYT(0,0,8),  GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(8), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       DO NG = 0, NGT
          call animes
          call gtextx(12.5,17.7,'@T=@',0)
          call gnumbr(13.1,17.7,GTX(NG),3,0)
          DO I = 0, j
             CALL TXGRFRS(I, GX, GYL(0:NRMAX,NG:NG,I), NRMAX, 1, STR(I), 0, IND, 0, 1, &
                  &       'ANIME', GMAX(I), GMIN(I))
          END DO
          call animee
       END DO

    CASE(-2)
       j = 0
       STR(j) = '@n$-b$=@'
       CALL APPROPGY(MODEG, GYT(0,0,12), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(12), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@u$-b$#q$#$=@'
       CALL APPROPGY(MODEG, GYT(0,0,19), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(19), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@u$-b$#f$#$=@'
       CALL APPROPGY(MODEG, GYT(0,0,13), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(13), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@Ripple n$-b$=@'
       CALL APPROPGY(MODEG, GYT(0,0,109), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(109), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@TOTAL N$-0$=@'
       CALL APPROPGY(MODEG, GYT(0,0,16), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(16), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@SLOW N$-0$=@'
       CALL APPROPGY(MODEG, GYT(0,0,35), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(35), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@THERMAL N$-0$=@'
       CALL APPROPGY(MODEG, GYT(0,0,36), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(36), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@NBI N$-0$=@'
       CALL APPROPGY(MODEG, GYT(0,0,138), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(138), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@q@'
       CALL APPROPGY(MODEG, GYT(0,0,20), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(20), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@j$-$#f$#$=@'
       CALL APPROPGY(MODEG, GYT(0,0,21), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(21), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@j$-e$#f$#$=@'
       CALL APPROPGY(MODEG, GYT(0,0,22), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(22), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@j$-i$#f$#$=@'
       CALL APPROPGY(MODEG, GYT(0,0,23), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(23), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@S@'
       CALL APPROPGY(MODEG, GYT(0,0,32), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(32), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@$#a$#@'
       CALL APPROPGY(MODEG, GYT(0,0,33), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(33), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@F$-CDBM$=@'
       CALL APPROPGY(MODEG, GYT(0,0,31), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(31), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@G$-1$=h$+2$=@'
       CALL APPROPGY(MODEG, GYT(0,0,30), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(30), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       DO NG = 0, NGT
          call animes
          call gtextx(12.5,17.7,'@T=@',0)
          call gnumbr(13.1,17.7,GTX(NG),3,0)
          DO I = 0, j
             CALL TXGRFRS(I, GX, GYL(0:NRMAX,NG:NG,I), NRMAX, 1, STR(I), 0, IND, 0, 1, &
                  &       'ANIME', GMAX(I), GMIN(I))
          END DO
          call animee
       END DO

    CASE(-3)
       j = 0
       STR(j) = '@LQm1@'
       CALL APPROPGY(MODEG, GYT(0,0,38), GYL(0,0,0), STR(j), NRMAX, NGT, gDIV(38), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@LQm2@'
       CALL APPROPGY(MODEG, GYT(0,0,39), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(39), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@LQm3@'
       CALL APPROPGY(MODEG, GYT(0,0,40), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(40), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@LQm4@'
       CALL APPROPGY(MODEG, GYT(0,0,41), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(41), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@LQm5@'
       CALL APPROPGY(MODEG, GYT(0,0,42), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(42), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = 8
       STR(j) = '@LQe1@'
       CALL APPROPGY(MODEG, GYT(0,0,43), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(43), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@LQe2@'
       CALL APPROPGY(MODEG, GYT(0,0,44), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(44), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@LQe3@'
       CALL APPROPGY(MODEG, GYT(0,0,45), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(45), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@LQe4@'
       CALL APPROPGY(MODEG, GYT(0,0,46), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(46), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@LQi1@'
       CALL APPROPGY(MODEG, GYT(0,0,50), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(50), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@LQi2@'
       CALL APPROPGY(MODEG, GYT(0,0,51), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(51), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@LQi3@'
       CALL APPROPGY(MODEG, GYT(0,0,52), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(52), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@LQi4@'
       CALL APPROPGY(MODEG, GYT(0,0,53), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(53), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       DO NG = 0, NGT
          call animes
          call gtextx(12.5,17.7,'@T=@',0)
          call gnumbr(13.1,17.7,GTX(NG),3,0)
          DO I = 0, 4
             CALL TXGRFRS(I, GX, GYL(0:NRMAX,NG:NG,I), NRMAX, 1, STR(I), 0, IND, 0, 1, &
                  &       'ANIME', GMAX(I), GMIN(I))
          END DO
          DO I = 8, j
             CALL TXGRFRS(I, GX, GYL(0:NRMAX,NG:NG,I), NRMAX, 1, STR(I), 0, IND, 0, 1, &
                  &       'ANIME', GMAX(I), GMIN(I))
          END DO
          call animee
       END DO

    CASE(-4)
       j = 0
       STR(j) = '@LQe5@'
       CALL APPROPGY(MODEG, GYT(0,0,47), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(47), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@LQe6@'
       CALL APPROPGY(MODEG, GYT(0,0,48), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(48), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@LQe7@'
       CALL APPROPGY(MODEG, GYT(0,0,49), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(49), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = 4
       STR(j) = '@LQi5@'
       CALL APPROPGY(MODEG, GYT(0,0,54), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(54), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@LQi6@'
       CALL APPROPGY(MODEG, GYT(0,0,55), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(55), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@LQi7@'
       CALL APPROPGY(MODEG, GYT(0,0,56), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(56), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = 8
       STR(j) = '@LQb1@'
       CALL APPROPGY(MODEG, GYT(0,0,57), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(57), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@LQb3@'
       CALL APPROPGY(MODEG, GYT(0,0,58), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(58), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@LQr1@'
       CALL APPROPGY(MODEG, GYT(0,0,62), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(62), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = 12
       STR(j) = '@LQn1@'
       CALL APPROPGY(MODEG, GYT(0,0,59), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(59), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@LQn2@'
       CALL APPROPGY(MODEG, GYT(0,0,60), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(60), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@LQn3@'
       CALL APPROPGY(MODEG, GYT(0,0,61), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(61), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       DO NG = 0, NGT
          call animes
          call gtextx(12.5,17.7,'@T=@',0)
          call gnumbr(13.1,17.7,GTX(NG),3,0)
          DO I = 0, 2
             CALL TXGRFRS(I, GX, GYL(:,NG:NG,I), NRMAX, 1, STR(I), 0, IND, 0, 1, &
                  &       'ANIME', GMAX(I), GMIN(I))
          END DO
          DO I = 4, 6
             CALL TXGRFRS(I, GX, GYL(:,NG:NG,I), NRMAX, 1, STR(I), 0, IND, 0, 1, &
                  &       'ANIME', GMAX(I), GMIN(I))
          END DO
          DO I = 8, 10
             CALL TXGRFRS(I, GX, GYL(:,NG:NG,I), NRMAX, 1, STR(I), 0, IND, 0, 1, &
                  &       'ANIME', GMAX(I), GMIN(I))
          END DO
          DO I = 12, 14
             CALL TXGRFRS(I, GX, GYL(:,NG:NG,I), NRMAX, 1, STR(I), 0, IND, 0, 1, &
                  &       'ANIME', GMAX(I), GMIN(I))
          END DO
          call animee
       END DO

    CASE(-5)
       j=0
       STR(j) = '@rMue@'
       CALL APPROPGY(MODEG, GYT(0,0,74), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(74), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@rMui@'
       CALL APPROPGY(MODEG, GYT(0,0,75), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(75), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@$#c$#$-e$=@'
       CALL APPROPGY(MODEG, GYT(0,0,70), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(70), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@$#c$#$-i$=@'
       CALL APPROPGY(MODEG, GYT(0,0,71), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(71), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@rNuL@'
       CALL APPROPGY(MODEG, GYT(0,0,63), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(63), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@rNuLTe@'
       CALL APPROPGY(MODEG, GYT(0,0,92), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(92), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@rNuLTi@'
       CALL APPROPGY(MODEG, GYT(0,0,93), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(93), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@rNuOL@'
       CALL APPROPGY(MODEG, GYT(0,0,82), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(82), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@rNuiCX@'
       CALL APPROPGY(MODEG, GYT(0,0,78), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(78), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@rNuION@'
       CALL APPROPGY(MODEG, GYT(0,0,80), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(80), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@rNu0e@'
       CALL APPROPGY(MODEG, GYT(0,0,68), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(68), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@rNu0i@'
       CALL APPROPGY(MODEG, GYT(0,0,69), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(69), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@rNuTei@'
       CALL APPROPGY(MODEG, GYT(0,0,67), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(67), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@D01@'
       CALL APPROPGY(MODEG, GYT(0,0,72), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(72), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@D02@'
       CALL APPROPGY(MODEG, GYT(0,0,73), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(73), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@D03@'
       CALL APPROPGY(MODEG, GYT(0,0,139), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(139), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

!!$       STR(15) = '@Vbpara@'
!!$       CALL APPROPGY(MODEG, GYT(0,0,109), GYL(0,0,15), STR(15), NRMAX, NGT, gDIV(109), &
!!$            &        GMAX=GMAX(15), GMIN=GMIN(15))

       DO NG = 0, NGT
          call animes
          call gtextx(12.5,17.7,'@T=@',0)
          call gnumbr(13.1,17.7,GTX(NG),3,0)
          DO I = 0, j
             CALL TXGRFRS(I, GX, GYL(:,NG:NG,I), NRMAX, 1, STR(I), 0, IND, 0, 1, &
                  &       'ANIME', GMAX(I), GMIN(I))
          END DO
          call animee
       END DO

    CASE(-6)
       j = 0
       STR(0) = '@xmu$-e11$=@'
       CALL APPROPGY(MODEG, GYT(0,0,146), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(146), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@xmu$-e12$=@'
       CALL APPROPGY(MODEG, GYT(0,0,147), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(147), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@xmu$-e21$=@'
       CALL APPROPGY(MODEG, GYT(0,0,148), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(148), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@xmu$-e22$=@'
       CALL APPROPGY(MODEG, GYT(0,0,149), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(149), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@xmu$-i11$=@'
       CALL APPROPGY(MODEG, GYT(0,0,150), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(150), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@xmu$-i12$=@'
       CALL APPROPGY(MODEG, GYT(0,0,151), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(151), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@xmu$-i21$=@'
       CALL APPROPGY(MODEG, GYT(0,0,152), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(152), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@xmu$-i22$=@'
       CALL APPROPGY(MODEG, GYT(0,0,153), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(153), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

!!$       j = j + 1
!!$       STR(j) = '@FWthe@'
!!$       CALL APPROPGY(MODEG, GYT(0,0,64), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(64), &
!!$            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@WPM@'
       CALL APPROPGY(MODEG, GYT(0,0,66), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(66), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@rNu0b@'
       CALL APPROPGY(MODEG, GYT(0,0,87), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(87), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@Vbpara@'
       CALL APPROPGY(MODEG, GYT(0,0,109), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(109), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@rNuB@'
       CALL APPROPGY(MODEG, GYT(0,0,79), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(79), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@rNubL@'
       CALL APPROPGY(MODEG, GYT(0,0,118), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(118), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@rNuLB@'
       CALL APPROPGY(MODEG, GYT(0,0,116), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(116), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

!!$       STR(15) = '@rNuLB@'
!!$       CALL APPROPGY(MODEG, GYT(0,0,116), GYL(0,0,15), STR(15), NRMAX, NGT, gDIV(116), &
!!$            &        GMAX=GMAX(15), GMIN=GMIN(15))

       DO NG = 0, NGT
          call animes
          call gtextx(12.5,17.7,'@T=@',0)
          call gnumbr(13.1,17.7,GTX(NG),3,0)
          DO I = 0, j
             CALL TXGRFRS(I, GX, GYL(:,NG:NG,I), NRMAX, 1, STR(I), 0, IND, 0, 1, &
                  &       'ANIME', GMAX(I), GMIN(I))
          END DO
          call animee
       END DO

    CASE(-7)
       j = 0
       STR(0) = '@l$-11$=$+ee$=$=@'
       CALL APPROPGY(MODEG, GYT(0,0,166), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(166), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@l$-12$=$+ee$=$=@'
       CALL APPROPGY(MODEG, GYT(0,0,167), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(167), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@l$-21$=$+ee$=$=@'
       CALL APPROPGY(MODEG, GYT(0,0,168), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(168), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@l$-22$=$+ee$=@'
       CALL APPROPGY(MODEG, GYT(0,0,169), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(169), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@l$-11$=$+ei$=@'
       CALL APPROPGY(MODEG, GYT(0,0,170), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(170), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@l$-12$=$+ei$=@'
       CALL APPROPGY(MODEG, GYT(0,0,171), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(171), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@l$-21$=$+ei$=@'
       CALL APPROPGY(MODEG, GYT(0,0,172), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(172), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@l$-22$=$+ei$=@'
       CALL APPROPGY(MODEG, GYT(0,0,173), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(173), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@l$-11$=$+ie$=@'
       CALL APPROPGY(MODEG, GYT(0,0,174), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(174), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@l$-12$=$+ie$=@'
       CALL APPROPGY(MODEG, GYT(0,0,175), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(175), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@l$-21$=$+ie$=@'
       CALL APPROPGY(MODEG, GYT(0,0,176), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(176), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@l$-22$=$+ie$=@'
       CALL APPROPGY(MODEG, GYT(0,0,177), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(177), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@l$-11$=$+ii$=@'
       CALL APPROPGY(MODEG, GYT(0,0,178), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(178), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@l$-12$=$+ii$=@'
       CALL APPROPGY(MODEG, GYT(0,0,179), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(179), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@l$-21$=$+ii$=@'
       CALL APPROPGY(MODEG, GYT(0,0,180), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(180), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@l$-22$=$+ii$=@'
       CALL APPROPGY(MODEG, GYT(0,0,181), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(181), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       DO NG = 0, NGT
          call animes
          call gtextx(12.5,17.7,'@T=@',0)
          call gnumbr(13.1,17.7,GTX(NG),3,0)
          DO I = 0, j
             CALL TXGRFRS(I, GX, GYL(:,NG:NG,I), NRMAX, 1, STR(I), 0, IND, 0, 1, &
                  &       'ANIME', GMAX(I), GMIN(I))
          END DO
          call animee
       END DO

    CASE(-8)
       j = 0
       STR(0) = '@l$-11$=$+ef$=@'
       CALL APPROPGY(MODEG, GYT(0,0,182), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(182), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@l$-21$=$+ef$=@'
       CALL APPROPGY(MODEG, GYT(0,0,183), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(183), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@l$-11$=$+fe$=@'
       CALL APPROPGY(MODEG, GYT(0,0,184), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(184), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@l$-11$=$+if$=@'
       CALL APPROPGY(MODEG, GYT(0,0,185), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(185), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@l$-11$=$+fi$=@'
       CALL APPROPGY(MODEG, GYT(0,0,186), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(186), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       j = j + 1
       STR(j) = '@l$-11$=$+ff$=@'
       CALL APPROPGY(MODEG, GYT(0,0,187), GYL(0,0,j), STR(j), NRMAX, NGT, gDIV(187), &
            &        GMAX=GMAX(j), GMIN=GMIN(j))

       DO NG = 0, NGT
          call animes
          call gtextx(12.5,17.7,'@T=@',0)
          call gnumbr(13.1,17.7,GTX(NG),3,0)
          DO I = 0, j
             CALL TXGRFRS(I, GX, GYL(:,NG:NG,I), NRMAX, 1, STR(I), 0, IND, 0, 1, &
                  &       'ANIME', GMAX(I), GMIN(I))
          END DO
          call animee
       END DO

    CASE(-9)
       STR(0) = '@E$-r$=@'
       CALL APPROPGY(MODEG, GYT(0,0,9),  GYL(0,0,0), STR(0), NRMAX, NGT, gDIV(9), &
            &        GMAX=GMAX(0), GMIN=GMIN(0))

       STR(1) = '@n$-e$=@'
       CALL APPROPGY(MODEG, GYT(0,0,1),  GYL(0,0,1), STR(1), NRMAX, NGT, gDIV(1), &
            &        GMAX=GMAX(1), GMIN=GMIN(1))

       STR(2) = '@n$-i$=@'
       CALL APPROPGY(MODEG, GYT(0,0,37), GYL(0,0,2), STR(2), NRMAX, NGT, gDIV(37), &
            &        GMAX=GMAX(2), GMIN=GMIN(2))

       STR(3) = '@n$-b$=@'
       CALL APPROPGY(MODEG, GYT(0,0,12), GYL(0,0,3), STR(3), NRMAX, NGT, gDIV(12), &
            &        GMAX=GMAX(3), GMIN=GMIN(3))

       STR(4) = '@E$-$#f$#$=@'
       CALL APPROPGY(MODEG, GYT(0,0,11), GYL(0,0,4), STR(4), NRMAX, NGT, gDIV(11), &
            &        GMAX=GMAX(4), GMIN=GMIN(4))

       STR(5) = '@T$-e$=@'
       CALL APPROPGY(MODEG, GYT(0,0,14), GYL(0,0,5), STR(5), NRMAX, NGT, gDIV(14), &
            &        GMAX=GMAX(5), GMIN=GMIN(5))

       STR(6) = '@T$-i$=@'
       CALL APPROPGY(MODEG, GYT(0,0,15), GYL(0,0,6), STR(6), NRMAX, NGT, gDIV(15), &
            &        GMAX=GMAX(6), GMIN=GMIN(6))

       STR(7) = '@B$-$#q$#$=@'
       CALL APPROPGY(MODEG, GYT(0,0,10), GYL(0,0,7), STR(7), NRMAX, NGT, gDIV(10), &
            &        GMAX=GMAX(7), GMIN=GMIN(7))

       STR(8) = '@u$-er$=@'
       CALL APPROPGY(MODEG, GYT(0,0,3),  GYL(0,0,8), STR(8), NRMAX, NGT, gDIV(3), &
            &        GMAX=GMAX(8), GMIN=GMIN(8))

       STR(9) = '@u$-e$#q$#$=@'
       CALL APPROPGY(MODEG, GYT(0,0,4),  GYL(0,0,9), STR(9), NRMAX, NGT, gDIV(4), &
            &        GMAX=GMAX(9), GMIN=GMIN(9))

       STR(10) = '@u$-e$#f$#$=@'
       CALL APPROPGY(MODEG, GYT(0,0,5),  GYL(0,0,10), STR(10), NRMAX, NGT, gDIV(5), &
            &        GMAX=GMAX(10), GMIN=GMIN(10))

       STR(11) = '@q@'
       CALL APPROPGY(MODEG, GYT(0,0,20), GYL(0,0,11), STR(11), NRMAX, NGT, gDIV(20), &
            &        GMAX=GMAX(11))!, GMIN(11) = 0.0)

       STR(12) = '@j$-r$=@'
       CALL APPROPGY(MODEG, GYT(0,0,117),  GYL(0,0,12), STR(12), NRMAX, NGT, gDIV(117), &
            &        GMAX=GMAX(12), GMIN=GMIN(12))

       STR(13) = '@u$-i$#q$#$=@'
       CALL APPROPGY(MODEG, GYT(0,0,7),  GYL(0,0,13), STR(13), NRMAX, NGT, gDIV(7), &
            &        GMAX=GMAX(13), GMIN=GMIN(13))

       STR(14) = '@u$-i$#f$#$=@'
       CALL APPROPGY(MODEG, GYT(0,0,8),  GYL(0,0,14), STR(14), NRMAX, NGT, gDIV(8), &
            &        GMAX=GMAX(14), GMIN=GMIN(14))

       STR(15) = '@j$-$#f$#$=@'
       CALL APPROPGY(MODEG, GYT(0,0,21), GYL(0,0,15), STR(15), NRMAX, NGT, gDIV(21), &
            &        GMAX=GMAX(15), GMIN=GMIN(15))

       DO NG = 0, NGT
          call animes
          call gtextx(12.5,17.7,'@T=@',0)
          call gnumbr(13.1,17.7,GTX(NG),3,0)
          DO I = 0, 15
             CALL TXGRFRS(I, GX, GYL(:,NG:NG,I), NRMAX, 1, STR(I), 0, IND, 0, 1, &
                  &       'ANIME', GMAX(I), GMIN(I))
          END DO
          call animee
       END DO

    CASE DEFAULT
       WRITE(6,*) 'Unknown NGYR: NGYR = ',NGYR
    END SELECT

    CALL SETFNT(IFNT)
    CALL PAGEE

    SELECT CASE(NGYR)
    CASE(-1)
       NGYR = -2  ; CYCLE
    CASE(-2)
       NGYR =  0  ; EXIT
    CASE(-3)
       NGYR = -4  ; CYCLE
    CASE(-4)
       NGYR =  0  ; EXIT
    CASE(-5)
       NGYR = -6  ; CYCLE
    CASE(-6)
       NGYR = -7  ; CYCLE
    CASE(-7)
       NGYR = -8  ; CYCLE
    CASE(-8)
       NGYR =  0  ; EXIT
    CASE(-9)
       NGYR =  0  ; EXIT
    CASE DEFAULT
       EXIT
    END SELECT

    END DO

    deallocate(GYL,STR,GMAX,GMIN)

    RETURN
  END SUBROUTINE TXGRFRA

  !***************************************************************
  !
  !   Comparison with radial profiles at one slice time
  !
  !***************************************************************

  SUBROUTINE TXGRCP(MODE)

    use tx_commons, only : NRMAX, DT, NRA, UsthNCL, BJPARA, BEpara, BJBSvar, ETAvar, Var
    integer(4), intent(in) :: MODE
    character(len=50) :: STR
    integer(4) :: IND, IFNT, NR, inum, ipage
    real(4), dimension(:,:), allocatable :: GYL, GYL2

    IF (NGR <= -1) THEN
       WRITE(6,*) 'G', NGR, ' has no data'
       RETURN
    END IF

    IF (MODEG == 2) THEN
       IND = 9
    ELSE
       IND = 0
    END IF

    allocate(GYL(0:NRMAX,1:4),GYL2(0:NRMAX,1:4))

    ! *** 1st page **

    ipage = 1
    CALL PAGES
    CALL SETCHS(0.3, 0.0)
    CALL SETLIN(0, 1, 7)
    CALL INQFNT(IFNT)

    CALL MOVE(2.0,17.7)
    CALL TEXT('[G', 2)
    CALL NUMBI(ipage, '(I2)', 2)
    CALL TEXT(']  ', 3)
    CALL TEXT('FROM', 4)
    CALL NUMBR(GT(0), '(1PE9.2)', 9)
    CALL TEXT(' TO', 3)
    CALL NUMBR(GT(NGR), '(1PE9.2)', 9)
    CALL TEXT('  DT =', 6)
    CALL NUMBD(DT, '(1PD9.2)', 9)
    CALL TEXT('  NGRSTP = ', 11)
    CALL NUMBI(NGRSTP,'(I4)',4)

!!    CALL SETFNT(44)
!    call gtextx(5.5,9.4,'@Please compare ETA and BJBS in the limit of Zeff = 1.@',0)
!!    CALL SETFNT(-1)

    ! Resistivity

    inum = 4
    DO NR = 0, NRMAX
       GYL(NR,1) = GLOG(ETAvar(NR,2),1.D-10,1.D0) ! NCLASS
       GYL(NR,2) = GLOG(ETAvar(NR,3),1.D-10,1.D0) ! Sauter
       GYL(NR,3) = GLOG(ETAvar(NR,5),1.D-10,1.D0) ! MI
       GYL(NR,4) = GLOG(ETAvar(NR,1),1.D-10,1.D0) ! TX, tentative
!       GYL(NR,3) = GLOG(ETAS(NR),1.D-10,1.D0)
    END DO

    STR = '@LOG: ETA NCLASS, SAUTER, MI, TX@'
    CALL TXGRFRS(0, GX(1:NRA), GYL(1:NRA,:), NRA-1, inum, STR, MODE, IND, 1, 0, 'STATIC')

    ! Bootstrap current

    inum = 3
    GYL(0:NRMAX,1) = REAL(BJBSvar(0:NRMAX,2)) ! NCLASS
    GYL(0,1) = 0.0
    GYL(0:NRMAX,2) = REAL(BJBSvar(0:NRMAX,3)) ! Sauter
    GYL(0:NRMAX,3) = REAL(BJBSvar(0:NRMAX,5)) ! MI

    STR = '@BJ$-BS$= NCLASS, SAUTER, MI@'
    CALL APPROPGY(MODEG, GYL, GYL2, STR, NRMAX, inum-1, gDIV(22))
    CALL TXGRFRS(1, GX(1:NRA), GYL2(1:NRA,:), NRA-1, inum, STR, MODE, IND, 0, 0, 'STATIC')

    ! Poloidal rotations

    GYL(0:NRMAX,1) = REAL(Var(0:NRMAX,1)%Uthhat) ! TX
    GYL(0:NRMAX,3) = REAL(UsthNCL(0:NRMAX,1,1))  ! Neo. solver, from p' T'
    GYL(0:NRMAX,4) = REAL(UsthNCL(0:NRMAX,1,2))  ! Neo. solver, from <E.B>
    GYL(0:NRMAX,2) = GYL(0:NRMAX,3) + GYL(0:NRMAX,4) ! Neo. solver, total

    STR = '@Ueth@'
    CALL APPROPGY(MODEG, GYL, GYL2, STR, NRMAX, 4-1, gDIV(4))
    CALL TXGRFRS(2, GX(1:NRA), GYL2(1:NRA,:), NRA-1, 4, STR, MODE, IND, 0, 0, 'STATIC')

    GYL(0:NRMAX,1) = REAL(Var(0:NRMAX,2)%Uthhat)
    GYL(0:NRMAX,3) = REAL(UsthNCL(0:NRMAX,2,1))
    GYL(0:NRMAX,4) = REAL(UsthNCL(0:NRMAX,2,2))
    GYL(0:NRMAX,2) = GYL(0:NRMAX,3) + GYL(0:NRMAX,4)

    STR = '@Uith@'
    CALL APPROPGY(MODEG, GYL, GYL2, STR, NRMAX, 4-1, gDIV(7))
    CALL TXGRFRS(3, GX(1:NRA), GYL2(1:NRA,:), NRA-1, 4, STR, MODE, IND, 0, 0, 'STATIC')

    CALL PAGEE

    ! *** 2nd page **

    ipage = 2
    CALL PAGES
    CALL SETCHS(0.3, 0.0)
    CALL SETLIN(0, 1, 7)
    CALL INQFNT(IFNT)

    CALL MOVE(2.0,17.7)
    CALL TEXT('[G', 2)
    CALL NUMBI(ipage, '(I2)', 2)
    CALL TEXT(']  ', 3)
    CALL TEXT('FROM', 4)
    CALL NUMBR(GT(0), '(1PE9.2)', 9)
    CALL TEXT(' TO', 3)
    CALL NUMBR(GT(NGR), '(1PE9.2)', 9)
    CALL TEXT('  DT =', 6)
    CALL NUMBD(DT, '(1PD9.2)', 9)
    CALL TEXT('  NGRSTP = ', 11)
    CALL NUMBI(NGRSTP,'(I4)',4)

    ! Total current

    inum = 4
    GYL(0:NRMAX,1) = real(BJPARA(0:NRMAX))                              ! TX
    GYL(0:NRMAX,2) = real(BEpara(0:NRMAX)/ETAvar(0:NRMAX,2)+BJBSvar(0:NRMAX,2)) ! NCLASS
    GYL(0:NRMAX,3) = real(BEpara(0:NRMAX)/ETAvar(0:NRMAX,3)+BJBSvar(0:NRMAX,3)) ! Sauter
    GYL(0:NRMAX,4) = real(BEpara(0:NRMAX)/ETAvar(0:NRMAX,5)+BJBSvar(0:NRMAX,5)) ! MI

    STR = '@BJ// TX, NCLASS, SAUTER, MI@'
    CALL APPROPGY(MODEG, GYL, GYL2, STR, NRMAX, inum-1, gDIV(22))
    CALL TXGRFRS(0, GX, GYL2, NRMAX, inum, STR, MODE, IND, 0, 0, 'STATIC')

    CALL PAGEE

    ! *** 3rd page **

    ipage = 3
    CALL PAGES
    CALL SETCHS(0.3, 0.0)
    CALL SETLIN(0, 1, 7)
    CALL INQFNT(IFNT)

    CALL MOVE(2.0,17.7)
    CALL TEXT('[G', 2)
    CALL NUMBI(ipage, '(I2)', 2)
    CALL TEXT(']  ', 3)
    CALL TEXT('FROM', 4)
    CALL NUMBR(GT(0), '(1PE9.2)', 9)
    CALL TEXT(' TO', 3)
    CALL NUMBR(GT(NGR), '(1PE9.2)', 9)
    CALL TEXT('  DT =', 6)
    CALL NUMBD(DT, '(1PD9.2)', 9)
    CALL TEXT('  NGRSTP = ', 11)
    CALL NUMBI(NGRSTP,'(I4)',4)

    ! Total current

    inum = 2
    GYL(0:NRMAX,1) = GY(0:NRMAX,NGR,26)  ! TX
    GYL(0:NRMAX,2) = GY(0:NRMAX,NGR,162) ! Neo. solver

    STR = '@Bu$-e//$= TX, Neo. solver@'
    CALL APPROPGY(MODEG, GYL, GYL2, STR, NRMAX, inum-1, gDIV(26))
    CALL TXGRFRS(0, GX(1:NRMAX), GYL2(1:NRMAX,:), NRMAX-1, inum, STR, MODE, IND, 0, 0, 'STATIC')
!    CALL TXGRFRS(0, GX, GYL2, NRMAX, inum, STR, MODE, IND, 0, 0, 'STATIC')

    GYL(0:NRMAX,1) = GY(0:NRMAX,NGR,28)  ! TX
    GYL(0:NRMAX,2) = GY(0:NRMAX,NGR,163) ! Neo. solver

    STR = '@Bu$-i//$= TX, Neo. solver@'
    CALL APPROPGY(MODEG, GYL, GYL2, STR, NRMAX, inum-1, gDIV(28))
    CALL TXGRFRS(1, GX(1:NRMAX), GYL2(1:NRMAX,:), NRMAX-1, inum, STR, MODE, IND, 0, 0, 'STATIC')
!    CALL TXGRFRS(1, GX, GYL2, NRMAX, inum, STR, MODE, IND, 0, 0, 'STATIC')

    GYL(0:NRMAX,1) = GY(0:NRMAX,NGR,158) ! TX
    GYL(0:NRMAX,2) = GY(0:NRMAX,NGR,164) ! Neo. solver

    STR = '@Bqhat$-e//$= TX, Neo. solver@'
    CALL APPROPGY(MODEG, GYL, GYL2, STR, NRMAX, inum-1, gDIV(158))
    CALL TXGRFRS(2, GX, GYL2, NRMAX, inum, STR, MODE, IND, 0, 0, 'STATIC')

    GYL(0:NRMAX,1) = GY(0:NRMAX,NGR,159) ! TX
    GYL(0:NRMAX,2) = GY(0:NRMAX,NGR,165) ! Neo. solver

    STR = '@Bqhat$-i//$= TX, Neo. solver@'
    CALL APPROPGY(MODEG, GYL, GYL2, STR, NRMAX, inum-1, gDIV(164))
    CALL TXGRFRS(3, GX, GYL2, NRMAX, inum, STR, MODE, IND, 0, 0, 'STATIC')

    CALL PAGEE

    deallocate(GYL,GYL2)

  END SUBROUTINE TXGRCP

  !***************************************************************
  !
  !   Animation of comparison with radial profiles
  !
  !***************************************************************

  SUBROUTINE TXGRCPA(MODE)

    use tx_commons, only : NRMAX, DT
    integer(4), intent(in) :: MODE
    integer(4) :: IND, IFNT, NG, I, N, j
    real(4), dimension(:),     allocatable :: GMAX, GMIN
    real(4), dimension(:,:,:), allocatable :: GYL
    character(len=50), dimension(:), allocatable :: STR

    IF (NGR <= -1) THEN
       WRITE(6,*) 'G', NGR, ' has no data'
       RETURN
    END IF

    IF (MODEG == 2) THEN
       IND = 9
    ELSE
       IND = 0
    END IF

    N = NRMAX
    allocate(GYL(0:N,0:NGT,0:7), STR(0:3), GMAX(0:7), GMIN(0:7))
    STR = ' '

    ! *** Correspond to 3rd page of TXGRCP ***

    CALL PAGES
    CALL SETCHS(0.3, 0.0)
    CALL SETLIN(0, 1, 7)
    CALL INQFNT(IFNT)

    CALL MOVE(2.0,17.7)
!    CALL TEXT('[G', 2)
!    CALL NUMBI(ipage, '(I2)', 2)
!    CALL TEXT(']  ', 3)
!    CALL TEXT('FROM', 4)
!    CALL NUMBR(GT(0), '(1PE9.2)', 9)
!    CALL TEXT(' TO', 3)
!    CALL NUMBR(GT(NGR), '(1PE9.2)', 9)
    CALL TEXT('  DT =', 6)
    CALL NUMBD(DT, '(1PD9.2)', 9)
!    CALL TEXT('  NGRSTP = ', 11)
!    CALL NUMBI(NGRSTP,'(I4)',4)

    ! Total current

    j = 0
    CALL APPROPGY(MODEG, GYT(0,0, 26), GYL(0,0,0), STR(j), N, NGT, gDIV( 26), &
         &        GMAX=GMAX(0), GMIN=GMIN(0)) ! TX
    STR(j) = '@Bu$-e//$= TX, Neo. solver@'
    CALL APPROPGY(MODEG, GYT(0,0,162), GYL(0,0,1), STR(j), N, NGT, gDIV(162), &
         &        GMAX=GMAX(1), GMIN=GMIN(1)) ! Neo. solver
    GMAX(j) = max(GMAX(0),GMAX(1))
    GMIN(j) = min(GMIN(0),GMIN(1))

    j = j + 1
    CALL APPROPGY(MODEG, GYT(0,0, 28), GYL(0,0,2), STR(j), N, NGT, gDIV( 28), &
         &        GMAX=GMAX(2), GMIN=GMIN(2)) ! TX
    STR(1) = '@Bu$-i//$= TX, Neo. solver@'
    CALL APPROPGY(MODEG, GYT(0,0,163), GYL(0,0,3), STR(j), N, NGT, gDIV(163), &
         &        GMAX=GMAX(3), GMIN=GMIN(3)) ! Neo. solver
!    GMAX(j) = max(GMAX(2),GMAX(3))
!    GMIN(j) = min(GMIN(2),GMIN(3))
    GMAX(j) = min(GMAX(2),GMAX(3))
    GMIN(j) = max(GMIN(2),GMIN(3))

    j = j + 1
    CALL APPROPGY(MODEG, GYT(0,0,158), GYL(0,0,4), STR(j), N, NGT, gDIV(158), &
         &        GMAX=GMAX(4), GMIN=GMIN(4)) ! TX
    STR(2) = '@Bqhat$-e//$= TX, Neo. solver@'
    CALL APPROPGY(MODEG, GYT(0,0,164), GYL(0,0,5), STR(j), N, NGT, gDIV(164), &
         &        GMAX=GMAX(5), GMIN=GMIN(5)) ! Neo. solver
    GMAX(j) = max(GMAX(4),GMAX(5))
    GMIN(j) = min(GMIN(4),GMIN(5))

    j = j + 1
    CALL APPROPGY(MODEG, GYT(0,0,159), GYL(0,0,6), STR(j), N, NGT, gDIV(159), &
         &        GMAX=GMAX(6), GMIN=GMIN(6)) ! TX
    STR(j) = '@Bqhat$-i//$= TX, Neo. solver@'
    CALL APPROPGY(MODEG, GYT(0,0,165), GYL(0,0,7), STR(j), N, NGT, gDIV(165), &
         &        GMAX=GMAX(7), GMIN=GMIN(7)) ! Neo. solver
    GMAX(j) = max(GMAX(6),GMAX(7))
    GMIN(j) = min(GMIN(6),GMIN(7))

    CALL SETFNT(44)
    do NG = 0, NGT
       call animes
       CALL SETLIN(-1, -1, 7)
       call gtextx(12.5,17.7,'@T=@',0)
       call gnumbr(13.1,17.7,GTX(NG),3,0)
       do I = 0, 1
          CALL TXGRFRS(I, GX(1:N), GYL(1:N,NG,2*I:2*I+1), N-1, 2, STR(I), MODE, IND, 0, 0, &
               &       'ANIME', GMAX(I), GMIN(I))
       end do
       do I = 2, 3
          CALL TXGRFRS(I, GX(0:N), GYL(0:N,NG,2*I:2*I+1), N  , 2, STR(I), MODE, IND, 0, 0, &
               &       'ANIME', GMAX(I), GMIN(I))
       end do
       call animee
    end do

!!$    GYL(0:NRMAX,1) = GY(0:NRMAX,NGR,28)  ! TX
!!$    GYL(0:NRMAX,2) = GY(0:NRMAX,NGR,163) ! Neo. solver
!!$
!!$    STR = '@Bu$-i//$= TX, Neo. solver@'
!!$    CALL APPROPGY(MODEG, GYL, GYL2, STR, NRMAX, inum-1, gDIV(28))
!!$    CALL TXGRFRS(1, GX(1:NRMAX), GYL2(1:NRMAX,:), NRMAX-1, inum, STR, MODE, IND, 0, 0, 'STATIC')
!!$
!!$    GYL(0:NRMAX,1) = GY(0:NRMAX,NGR,158) ! TX
!!$    GYL(0:NRMAX,2) = GY(0:NRMAX,NGR,164) ! Neo. solver
!!$
!!$    STR = '@Bqhat$-e//$= TX, Neo. solver@'
!!$    CALL APPROPGY(MODEG, GYL, GYL2, STR, NRMAX, inum-1, gDIV(158))
!!$    CALL TXGRFRS(2, GX, GYL2, NRMAX, inum, STR, MODE, IND, 0, 0, 'STATIC')
!!$
!!$    GYL(0:NRMAX,1) = GY(0:NRMAX,NGR,159) ! TX
!!$    GYL(0:NRMAX,2) = GY(0:NRMAX,NGR,165) ! Neo. solver
!!$
!!$    STR = '@Bqhat$-i//$= TX, Neo. solver@'
!!$    CALL APPROPGY(MODEG, GYL, GYL2, STR, NRMAX, inum-1, gDIV(164))
!!$    CALL TXGRFRS(3, GX, GYL2, NRMAX, inum, STR, MODE, IND, 0, 0, 'STATIC')

    CALL PAGEE

    deallocate(GYL,STR,GMAX,GMIN)

  END SUBROUTINE TXGRCPA

  !***************************************************************
  !
  !   SLAVE ROUTINE TO Write graph of GX, GY
  !
  !***************************************************************

  SUBROUTINE TXGRFRX(K, GXL, GYL, NRMAX, NGMAX, STR, MODE, IND, &
       &             GXMIN, GYMAX, GYMIN, ILOGIN, GPXY_IN)

    use tx_commons, only : rhob
    INTEGER(4), INTENT(IN) :: K, NRMAX, NGMAX, MODE, IND
    REAL(4), DIMENSION(:), INTENT(IN) :: GXL
    REAL(4), DIMENSION(0:NRMAX,1:NGMAX+1), INTENT(IN) :: GYL
    real(4), intent(in), optional :: GXMIN, GYMAX, GYMIN
    integer(4), intent(in), optional :: ILOGIN
    real(4), dimension(:), intent(in), optional :: GPXY_IN
    character(len=*), INTENT(IN) :: STR
    integer(4) :: ILOG, IPRES
    REAL(4) :: GXMAX, GXMINL
    REAL(4), DIMENSION(4) :: GPXY

    if(present(GPXY_IN)) then
       GPXY(1:4) = GPXY_IN(1:4)
    else
       GPXY(1) =  3.0 + 12.5 * MOD(K,2)
       GPXY(2) = 12.5 + 12.5 * MOD(K,2)
       GPXY(3) = 10.5 -  8.5 * REAL(K/2)
       GPXY(4) = 17.0 -  8.5 * REAL(K/2)
       ! square
!!$       GPXY(1) =  3.0 + 12.5 * MOD(K,2)
!!$       GPXY(2) = 10.4286 + 12.5 * MOD(K,2)
!!$       GPXY(3) = 10.5 -  8.5 * REAL(K/2)
!!$       GPXY(4) = 17.0 -  8.5 * REAL(K/2)
    end if

    GXMAX = REAL(rhob)

    IF(PRESENT(GXMIN)) THEN
       GXMINL = GXMIN
    ELSE
       GXMINL = 0.0
    END IF
    IF(PRESENT(ILOGIN)) THEN
       ! Semi-Log scale
       ILOG = ILOGIN
    ELSE
       ! Normal scale
       ILOG = 0
    END IF

    IPRES = 0
    IF(PRESENT(GYMAX)) IPRES = IPRES + 1
    IF(PRESENT(GYMIN)) IPRES = IPRES + 2

    IF(IPRES == 0) THEN
       CALL TXGRAF(GPXY, GXL, GYL, NRMAX+1, NRMAX+1, NGMAX+1, &
            &            GXMINL, GXMAX, STR, 0.3, MODE, IND, ILOG)
    ELSEIF(IPRES == 1) THEN
       CALL TXGRAF(GPXY, GXL, GYL, NRMAX+1, NRMAX+1, NGMAX+1, &
            &            GXMINL, GXMAX, STR, 0.3, MODE, IND, ILOG, GYMAX_IN=GYMAX)
    ELSEIF(IPRES == 2) THEN
       CALL TXGRAF(GPXY, GXL, GYL, NRMAX+1, NRMAX+1, NGMAX+1, &
            &            GXMINL, GXMAX, STR, 0.3, MODE, IND, ILOG, GYMIN_IN=GYMIN)
    ELSE
       CALL TXGRAF(GPXY, GXL, GYL, NRMAX+1, NRMAX+1, NGMAX+1, &
            &            GXMINL, GXMAX, STR, 0.3, MODE, IND, ILOG, GYMAX_IN=GYMAX, &
            &                                                      GYMIN_IN=GYMIN)
    END IF

    RETURN
  END SUBROUTINE TXGRFRX

  !***************************************************************
  !
  !   SLAVE ROUTINE TO Write graph of GX, GY (small version)
  !
  !***************************************************************

  SUBROUTINE TXGRFRXS(K, GXL, GYL, NRMAX, NGMAX, STR, MODE, IND, GYMAX, ILOGIN)

    use tx_commons, only : rhob
    INTEGER(4), INTENT(IN) :: K, NRMAX, NGMAX, MODE, IND
    REAL(4), DIMENSION(:), INTENT(IN) :: GXL
    REAL(4), DIMENSION(0:NRMAX,1:NGMAX+1), INTENT(IN) :: GYL
    character(len=*), INTENT(IN) :: STR
    real(4), intent(in), optional :: GYMAX
    integer(4), intent(in), optional :: ILOGIN
    integer(4) :: ILOG
    REAL(4) :: GXMAX
    REAL(4), DIMENSION(4) :: GPXY

    GPXY(1) =   2.0 + 6.1  * MOD(K,4)
    GPXY(2) =   6.7 + 6.1  * MOD(K,4)
    GPXY(3) = 13.75 - 4.25 * REAL(K/4)
    GPXY(4) = 17.0  - 4.25 * REAL(K/4)
    GXMAX = REAL(rhob)

    IF(PRESENT(ILOGIN)) THEN
       ILOG = ILOGIN
    ELSE
       ILOG = 0
    END IF

    IF(PRESENT(GYMAX)) THEN
       CALL TXGRAF(GPXY, GXL, GYL, NRMAX+1, NRMAX+1, NGMAX+1, &
            &            0.0, GXMAX, STR, 0.26, MODE, IND, ILOG, GYMAX)
    ELSE
       CALL TXGRAF(GPXY, GXL, GYL, NRMAX+1, NRMAX+1, NGMAX+1, &
            &            0.0, GXMAX, STR, 0.26, MODE, IND, ILOG)
    END IF

    RETURN
  END SUBROUTINE TXGRFRXS

  !***************************************************************
  !
  !   SLAVE ROUTINE TO Write graph of GX, GY at one slice time
  !
  !     All arguments are input
  !
  !***************************************************************

  SUBROUTINE TXGRFRS(K, GXL, GYL, NXMAX, NGMAX, STR, MODE, IND, ILOGIN, ISIZE, KINDEX, &
       &             GMAX, GMIN)
    ! Forth argument, NXMAX, is not always "NRMAX" defined as the number of grid points.

    use tx_commons, only : rhob
    INTEGER(4), INTENT(IN) :: K, NXMAX, NGMAX, MODE, IND, ISIZE
    real(4), intent(in), optional :: GMAX, GMIN
    REAL(4), DIMENSION(:), INTENT(IN) :: GXL
    REAL(4), DIMENSION(:,:), INTENT(IN) :: GYL
    character(len=*), INTENT(IN) :: STR, KINDEX
    integer(4), intent(in), optional :: ILOGIN
    integer(4) :: ILOG
    REAL(4) :: GXMAX, FNTSIZE
    REAL(4), DIMENSION(4) :: GPXY

    FNTSIZE = 0.3
    IF(ISIZE == 1) THEN
       ! Small
       GPXY(1) =   2.0 + 6.1  * MOD(K,4)
       GPXY(2) =   6.7 + 6.1  * MOD(K,4)
       GPXY(3) = 13.75 - 4.25 * REAL(K/4)
       GPXY(4) = 17.0  - 4.25 * REAL(K/4)
       FNTSIZE = 0.26
    ELSE IF(ISIZE == 2) THEN
       ! Square
       GPXY(1) =  3.0 + 12.5 * MOD(K,2)
       GPXY(2) = 10.4286 + 12.5 * MOD(K,2)
       GPXY(3) = 10.5 -  8.5 * REAL(K/2)
       GPXY(4) = 17.0 -  8.5 * REAL(K/2)
    ELSE IF(ISIZE == 3) THEN
       ! One graph per page (K is no longer valid.)
       GPXY(1) =  3.0
       GPXY(2) = 23.0
       GPXY(3) =  2.0
       GPXY(4) = 17.0
    ELSE IF(ISIZE == 4) THEN
       ! Six graphs per page
       GPXY(1) =  1.8  + 8.45 * MOD(K,3)
       GPXY(2) =  8.45 + 8.45 * MOD(K,3)
       GPXY(3) = 10.55 - 7.55 * REAL(K/3)
       GPXY(4) = 15.1  - 7.55 * REAL(K/3)
    ELSE IF(ISIZE == 5) THEN
       ! Standard (lower set)
       GPXY(1) =  3.0 + 12.5 * MOD(K,2)
       GPXY(2) = 12.5 + 12.5 * MOD(K,2)
       GPXY(3) =  9.5 -  8.5 * REAL(K/2)
       GPXY(4) = 16.0 -  8.5 * REAL(K/2)
    ELSE IF(ISIZE == 6) THEN
       ! Six graphs per page (lower set)
       GPXY(1) =  1.8  + 8.45 * MOD(K,3)
       GPXY(2) =  8.45 + 8.45 * MOD(K,3)
       GPXY(3) =  8.55 - 7.55 * REAL(K/3)
       GPXY(4) = 13.1  - 7.55 * REAL(K/3)
    ELSE
       ! Standard
       GPXY(1) =  3.0 + 12.5 * MOD(K,2)
       GPXY(2) = 12.5 + 12.5 * MOD(K,2)
       GPXY(3) = 10.5 -  8.5 * REAL(K/2)
       GPXY(4) = 17.0 -  8.5 * REAL(K/2)
    END IF
    GXMAX = REAL(rhob)

    IF(PRESENT(ILOGIN)) THEN
       ILOG = ILOGIN
    ELSE
       ILOG = 0
    END IF

    IF(KINDEX == 'STATIC') THEN
       CALL TXGRAF(GPXY, GXL, GYL(1:NXMAX+1,1:NGMAX), NXMAX+1, NXMAX+1, NGMAX, &
            &            0.0, GXMAX, STR, FNTSIZE, MODE, IND, ILOG)
    ELSE IF(KINDEX == 'ANIME') THEN
       CALL TXGRAF(GPXY, GXL, GYL, NXMAX+1, NXMAX+1, NGMAX, &
            &            0.0, GXMAX, STR, FNTSIZE, MODE, IND, ILOG, &
            &            GYMAX_IN=GMAX, GYMIN_IN=GMIN)!, KINDEX=KINDEX)
    ELSE
       CALL TXGRAF(GPXY, GXL, GYL, NXMAX+1, NXMAX+1, NGMAX, &
            &            0.0, GXMAX, STR, FNTSIZE, MODE, IND, ILOG)
    END IF

    RETURN
  END SUBROUTINE TXGRFRS

  !***************************************************************
  !
  !   Write graph of GVX, GVY
  !
  !***************************************************************

  SUBROUTINE TXGRFV(NGYV,MODE)

    use tx_commons, only : DT
    INTEGER(4), INTENT(IN) :: NGYV, MODE
    INTEGER(4) :: IND
    character(len=50) ::  STR
    REAL(4), DIMENSION(0:NGVM,1:4) :: GVYL

    IF (NGVV <= 1) THEN
       WRITE(6,*) 'GV', NGYV, ' has no data'
       RETURN
    END IF

    IF (MODEG == 2) THEN
       IND = 9
    ELSE
       IND = 0
    END IF

    CALL PAGES
    CALL SETCHS(0.3, 0.0)
    CALL SETLIN(0, 1, 7)

    CALL MOVE(2.0,17.7)
    CALL TEXT('[GV', 3)
    CALL NUMBI(NGYV, '(I2)', 2)
    CALL TEXT(']  ', 3)
    CALL TEXT('FROM', 4)
    CALL NUMBR(GVX(0), '(1PE9.2)', 9)
    CALL TEXT(' TO', 3)
    CALL NUMBR(GVX(NGVV), '(1PE9.2)', 9)
    CALL TEXT('  DT =', 6)
    CALL NUMBD(DT, '(1PD9.2)', 9)
    CALL TEXT('  NGVSTP = ', 11)
    CALL NUMBI(NGVSTP,'(I4)',4)

    SELECT CASE(NGYV)
    CASE(1)
       STR = '@n$-e$=(0),LAVE(0),(0.24),(0.6)@'
       CALL APPROPGY(MODEG, GVY(0,38), GVYL, STR, NGVV, 3, gDIV(1), GVY(0, 1))
       CALL TXGRFVX(0, GVX, GVYL, NGVM, NGVV, 4, STR, MODE, IND)

       STR = '@Z*n$-i$=+Z*n$-b$=-n$-e$=(0)@'
       CALL APPROPGY(MODEG, GVY(0, 2), GVYL, STR, NGVV, 1, gDIV(2))
       CALL TXGRFVX(1, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@E$-r$=(a/2)@'
       CALL APPROPGY(MODEG, GVY(0, 9), GVYL, STR, NGVV, 1, gDIV(9))
       CALL TXGRFVX(2, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

       CALL TXWPGR

    CASE(2)
       STR = '@u$-er$=(a/2)@'
       CALL TXGRFVX(0, GVX, GVY(0, 3), NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@u$-e$#q$#$=(a/2)@'
       CALL APPROPGY(MODEG, GVY(0, 4), GVYL, STR, NGVV, 1, gDIV(4))
       CALL TXGRFVX(1, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@u$-e$#f$#$=(0)@'
       CALL APPROPGY(MODEG, GVY(0, 5), GVYL, STR, NGVV, 1, gDIV(5))
       CALL TXGRFVX(2, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

       CALL TXWPGR

    CASE(3)
       STR = '@u$-ir$=(a/2)@'
       CALL TXGRFVX(0, GVX, GVY(0, 6), NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@u$-i$#q$#$=(a/2)@'
       CALL APPROPGY(MODEG, GVY(0, 7), GVYL, STR, NGVV, 1, gDIV(7))
       CALL TXGRFVX(1, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@u$-i$#f$#$=(a/2)@'
       CALL APPROPGY(MODEG, GVY(0, 8), GVYL, STR, NGVV, 1, gDIV(8))
       CALL TXGRFVX(2, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

       CALL TXWPGR

    CASE(4)
       STR = '@q(0)@'
       CALL TXGRFVX(0, GVX, GVY(0,20), NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@E$-tor$=(b)@'
       CALL TXGRFVX(1, GVX, GVY(0,11), NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@j$-$#f$#$=(0)@'
       CALL APPROPGY(MODEG, GVY(0,21), GVYL, STR, NGVV, 1, gDIV(21))
       CALL TXGRFVX(2, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@B$-$#q$#$=(a)@'
       CALL TXGRFVX(3, GVX, GVY(0,10), NGVM, NGVV, 1, STR, MODE, IND)

!       CALL TXWPGR

    CASE(5)
       STR = '@n$-b$=(0)@'
       CALL APPROPGY(MODEG, GVY(0,12), GVYL, STR, NGVV, 1, gDIV(12))
       CALL TXGRFVX(0, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@u$-b$#f$#$=(0)@'
       CALL APPROPGY(MODEG, GVY(0,13), GVYL, STR, NGVV, 1, gDIV(13))
       CALL TXGRFVX(1, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@Ripple n$-b$=(0)@'
       CALL APPROPGY(MODEG, GVY(0,50), GVYL, STR, NGVV, 1, gDIV(109))
       CALL TXGRFVX(2, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@Bu$-b//$=(a/2)@'
       CALL APPROPGY(MODEG, GVY(0,19), GVYL, STR, NGVV, 1, gDIV(19))
       CALL TXGRFVX(3, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

!       CALL TXWPGR

    CASE(6)
       STR = '@j$-e$#f$#$=(0)@'
       CALL APPROPGY(MODEG, GVY(0,22), GVYL, STR, NGVV, 1, gDIV(22))
       CALL TXGRFVX(0, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@j$-i$#f$#$=(0)@'
       CALL APPROPGY(MODEG, GVY(0,23), GVYL, STR, NGVV, 1, gDIV(23))
       CALL TXGRFVX(1, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@j$-b$#f$#$=(0)@'
       CALL APPROPGY(MODEG, GVY(0,24), GVYL, STR, NGVV, 1, gDIV(24))
       CALL TXGRFVX(2, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@j$-r$=(a/2)@'
       CALL TXGRFVX(3, GVX, GVY(0,56), NGVM, NGVV, 1, STR, MODE, IND)

!       CALL TXWPGR

    CASE(7)
       STR = '@u$-er$=(a/2)@'
       CALL TXGRFVX(0, GVX, GVY(0, 3), NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@u$-e$#$/136$#$=(a/2)@'
       CALL APPROPGY(MODEG, GVY(0,25), GVYL, STR, NGVV, 1, gDIV(25))
       CALL TXGRFVX(1, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@Bu$-e//$=(a/2)@'
       CALL APPROPGY(MODEG, GVY(0,26), GVYL, STR, NGVV, 1, gDIV(26))
       CALL TXGRFVX(2, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

       CALL TXWPGR

    CASE(8)
       STR = '@u$-ir$=(a/2)@'
       CALL TXGRFVX(0, GVX, GVY(0, 6), NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@u$-i$#$/136$#$=(a/2)@'
       CALL APPROPGY(MODEG, GVY(0,27), GVYL, STR, NGVV, 1, gDIV(27))
       CALL TXGRFVX(1, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@Bu$-i//$=(a/2)@'
       CALL APPROPGY(MODEG, GVY(0,28), GVYL, STR, NGVV, 1, gDIV(28))
       CALL TXGRFVX(2, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

       CALL TXWPGR

    CASE(9)
       STR = '@D$-i,eff$=(a/2)@'
       CALL TXGRFVX(0, GVX, GVY(0,52), NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@D$-i$=(a/2)+D$-e$=(a/2)@'
       CALL TXGRFVX(1, GVX, GVY(0,29), NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@$#c$#$-tbe$=(a/2)@'
       CALL TXGRFVX(2, GVX, GVY(0,44), NGVM, NGVV, 1, STR, MODE, IND)

!!$       STR = '@G$-1$=h$+2$=(a/2)@'
!!$       CALL TXGRFVX(1, GVX, GVY(0,30), NGVM, NGVV, 1, STR, MODE, IND)
!!$
!!$       STR = '@F$-CDBM$=(a/2)@'
!!$       CALL TXGRFVX(2, GVX, GVY(0,31), NGVM, NGVV, 1, STR, MODE, IND)

       CALL TXWPGR

    CASE(10)
       STR = '@T$-e$=(0)@'
       CALL TXGRFVX(0, GVX, GVY(0,14), NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@T$-i$=(0)@'
       CALL TXGRFVX(1, GVX, GVY(0,14), NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@p$-e$=(0)@'
       CALL TXGRFVX(2, GVX, GVY(0,53), NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@p$-i$=(0)@'
       CALL TXGRFVX(3, GVX, GVY(0,54), NGVM, NGVV, 1, STR, MODE, IND)

!       CALL TXWPGR

    CASE(11)
       STR = '@s(a/2)@'
       CALL TXGRFVX(0, GVX, GVY(0,32), NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@$#a$#(a/2)@'
       CALL TXGRFVX(1, GVX, GVY(0,33), NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@$#k$#(a/2)@'
       CALL TXGRFVX(2, GVX, GVY(0,34), NGVM, NGVV, 1, STR, MODE, IND)

       CALL TXWPGR

    CASE(12)
       STR = '@BE$-pol$=(a/2)@'
       CALL TXGRFVX(0, GVX, GVY(0,17), NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@I(0)@'
       CALL TXGRFVX(1, GVX, GVY(0,18), NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@uhat$-i$#q$#$=(a/2)@'
       CALL TXGRFVX(2, GVX, GVY(0,57), NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@1+2qhat$+2$=@'
       CALL TXGRFVX(3, GVX, GVY(0,58), NGVM, NGVV, 1, STR, MODE, IND)

!       CALL TXWPGR

    CASE(13)
       STR = '@SLOW N$-0$=(a)@'
       CALL APPROPGY(MODEG, GVY(0,35), GVYL, STR, NGVV, 1, gDIV(35))
       CALL TXGRFVX(0, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@THERMAL N$-0$=(0)@'
       CALL APPROPGY(MODEG, GVY(0,36), GVYL, STR, NGVV, 1, gDIV(36))
       CALL TXGRFVX(1, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@HALO N$-0$=@'
       CALL APPROPGY(MODEG, GVY(0,50), GVYL, STR, NGVV, 1, gDIV(138))
       CALL TXGRFVX(2, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@TOTAL N$-0$=(a)@'
       CALL APPROPGY(MODEG, GVY(0,16), GVYL, STR, NGVV, 1, gDIV(16))
       CALL TXGRFVX(3, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

!       CALL TXWPGR

    CASE(14)
       STR = '@PIE@'
       CALL TXGRFVX(0, GVX, GVY(0,46), NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@PCX@'
       CALL APPROPGY(MODEG, GVY(0,47), GVYL, STR, NGVV, 1, gDIV(89))
       CALL TXGRFVX(1, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@SIE@'
       CALL APPROPGY(MODEG, GVY(0,48), GVYL, STR, NGVV, 1, gDIV(90))
       CALL TXGRFVX(2, GVX, GVYL, NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@PBr@'
       CALL TXGRFVX(3, GVX, GVY(0,49), NGVM, NGVV, 1, STR, MODE, IND)
       !CALL TXWPGR

    CASE(15)
       STR = '@n$-e$=:total@'
       CALL TXGRFVX(0, GVX, GVY(0,41), NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@$#g$#:total@'
       CALL TXGRFVX(1, GVX, GVY(0,42), NGVM, NGVV, 1, STR, MODE, IND)

       STR = '@$#t$#$-p$=@'
       CALL TXGRFVX(2, GVX, GVY(0,43), NGVM, NGVV, 1, STR, MODE, IND)

       CALL TXWPGR

    CASE DEFAULT
       WRITE(6,*) 'Unknown NGYV: NGYV = ',NGYV
    END SELECT

    CALL PAGEE

    RETURN
  END SUBROUTINE TXGRFV

  !***************************************************************
  !
  !   Write graph of GTX, GTY
  !
  !***************************************************************

  SUBROUTINE TXGRFT(NGYT,MODE)

    use tx_commons, only : DT
    INTEGER(4), INTENT(IN) :: NGYT, MODE
    INTEGER(4) :: IND
    character(len=50) ::  STR

    IF (NGT <= 1) THEN
       WRITE(6,*) 'GT', NGYT, ' has no data'
       RETURN
    END IF

    IF (MODEG == 2) THEN
       IND = 9
    ELSE
       IND = 0
    END IF

    CALL PAGES
    CALL SETCHS(0.3, 0.0)
    CALL SETLIN(0, 1, 7)

    CALL MOVE(2.0,17.7)
    CALL TEXT('[GT', 3)
    CALL NUMBI(NGYT, '(I2)', 2)
    CALL TEXT(']  ', 3)
    CALL TEXT('FROM', 4)
    CALL NUMBR(GTX(0), '(1PE9.2)', 9)
    CALL TEXT(' TO', 3)
    CALL NUMBR(GTX(NGT), '(1PE9.2)', 9)
    CALL TEXT('  DT =', 6)
    CALL NUMBD(DT, '(1PD9.2)', 9)
    CALL TEXT('  NGTSTP = ', 11)
    CALL NUMBI(NGTSTP,'(I4)',4)

    SELECT CASE(NGYT)
    CASE(1)
       STR = '@Te,Ti,<Te>,<Ti> [keV] vs t@'
       CALL TXGRFVX(0, GTX, GTY(0,1), NGTM, NGT, 4, STR, MODE, IND)

       STR = '@PIN,POH,PNB,PRF,PNF [MW] vs t@'
       CALL TXGRFVX(1, GTX, GTY(0,5), NGTM, NGT, 5, STR, MODE, IND)

       STR = '@IP,IOH,INB,IOH+INB [MA] vs t@'
       CALL TXGRFVX(2, GTX, GTY(0,10), NGTM, NGT, 4, STR, MODE, IND)

       STR = '@POUT,PCX,PIE,PRL,PCON [MW] vs t@'
       CALL TXGRFVX(3, GTX, GTY(0,15), NGTM, NGT, 5, STR, MODE, IND)

    CASE(2)
       STR = '@QF vs t@'
       CALL TXGRFVX(0, GTX, GTY(0,20), NGTM, NGT, 1, STR, MODE, IND)

       STR = '@Te0,Ti0 [keV] vs t@'
       CALL TXGRFVX(1, GTX, GTY(0,21), NGTM, NGT, 2, STR, MODE, IND)

       STR = '@<Te>,<Ti> [keV] vs t@'
       CALL TXGRFVX(2, GTX, GTY(0,23), NGTM, NGT, 2, STR, MODE, IND)

       STR = '@Ne0,Ni0,<Ne>,<Ni> [10$+20$=/m$+3$=] vs t@'
       CALL TXGRFVX(3, GTX, GTY(0,25), NGTM, NGT, 4, STR, MODE, IND)

    CASE(3)
       STR = '@WF,WB,Wi,We [MJ] vs t@'
       CALL TXGRFVX(0, GTX, GTY(0,29), NGTM, NGT, 4, STR, MODE, IND)

       STR = '@TAUE1,TAUE2,TAUEP [s] vs t@'
       CALL TXGRFVX(1, GTX, GTY(0,33), NGTM, NGT, 3, STR, MODE, IND)

       STR = '@BETAa,BETA0 [%] vs t@'
       CALL TXGRFVX(2, GTX, GTY(0,36), NGTM, NGT, 2, STR, MODE, IND)

       STR = '@BETAPa,BETAP0 vs t@'
       CALL TXGRFVX(3, GTX, GTY(0,38), NGTM, NGT, 2, STR, MODE, IND)

    CASE(4)
       STR = '@VLOOP [V]@'
       CALL TXGRFVX(0, GTX, GTY(0,40), NGTM, NGT, 1, STR, MODE, IND)

       STR = '@ALI@'
       CALL TXGRFVX(1, GTX, GTY(0,41), NGTM, NGT, 1, STR, MODE, IND)

       STR = '@Q(0)@'
       CALL TXGRFVX(2, GTX, GTY(0,42), NGTM, NGT, 1, STR, MODE, IND)

       STR = '@RQ1 [m]@'
       CALL TXGRFVX(3, GTX, GTY(0,43), NGTM, NGT, 1, STR, MODE, IND)

    CASE(5)
       STR = '@TAUP, TAUPA [s] vs t@'
       CALL TXGRFVX(0, GTX, GTY(0,49), NGTM, NGT, 2, STR, MODE, IND)

       STR = '@Gamma_a [10$+20$=/m$+2$=s$+1$=] vs t@'
       CALL TXGRFVX(1, GTX, GTY(0,51), NGTM, NGT, 1, STR, MODE, IND)

    CASE(6)
       STR = '@Te,Ti,<Te>,<Ti> [keV] vs t@'
       CALL TXGRFTX(0, GTX, GTY(0, 1), NGTM, NGT, 4, STR, IND)

       STR = '@IP,IOH,INB,IBS,ITOT [MA] vs t@'
       CALL TXGRFTX(2, GTX, GTY(0,10), NGTM, NGT, 5, STR, IND)

       STR = '@PIN,POH,PNB,PRF,PNF [MW] vs t@'
       CALL TXGRFTX(4, GTX, GTY(0, 5), NGTM, NGT, 5, STR, IND)

       STR = '@POUT,PCX,PIE,PRL,PCON [MW] vs t@'
       CALL TXGRFTX(6, GTX, GTY(0,15), NGTM, NGT, 5, STR, IND)

       STR = '@QF vs t@'
       CALL TXGRFTX(1, GTX, GTY(0,20), NGTM, NGT, 1, STR, IND)

       STR = '@<Te>,<Ti> [keV] vs t@'
       CALL TXGRFTX(3, GTX, GTY(0,23), NGTM, NGT, 2, STR, IND)

       STR = '@Te0,Ti0 [keV] vs t@'
       CALL TXGRFTX(5, GTX, GTY(0,21), NGTM, NGT, 2, STR, IND)

       STR = '@Ne0,Ni0,<Ne>,<Ni> [10$+20$=/m$+3$=] vs t@'
       CALL TXGRFTX(7, GTX, GTY(0,25), NGTM, NGT, 4, STR, IND)

    CASE(7)
       STR = '@WF,WB,Wi,We [MJ] vs t@'
       CALL TXGRFTX(0, GTX, GTY(0,29), NGTM, NGT, 4, STR, IND)

       STR = '@BETAa,BETA0 [%] vs t@'
       CALL TXGRFTX(2, GTX, GTY(0,36), NGTM, NGT, 2, STR, IND)

       STR = '@TAUE1,TAUE2,TAUEP [s] vs t@'
       CALL TXGRFTX(4, GTX, GTY(0,33), NGTM, NGT, 3, STR, IND)

       STR = '@BETAPa,BETAP0 vs t@'
       CALL TXGRFTX(6, GTX, GTY(0,38), NGTM, NGT, 2, STR, IND)

       STR = '@VLOOP [V]@'
       CALL TXGRFTX(1, GTX, GTY(0,40), NGTM, NGT, 1, STR, IND)

       STR = '@Q(0)@'
       CALL TXGRFTX(3, GTX, GTY(0,42), NGTM, NGT, 1, STR, IND)

       STR = '@ALI@'
       CALL TXGRFTX(5, GTX, GTY(0,41), NGTM, NGT, 1, STR, IND)

       STR = '@RQ1 [m]@'
       CALL TXGRFTX(7, GTX, GTY(0,43), NGTM, NGT, 1, STR, IND)

    CASE(8)
       STR = '@Nb0,<Nb> [10$+18$=/m$+3$=] vs t@'
       CALL TXGRFTX(0, GTX, GTY(0,44), NGTM, NGT, 2, STR, IND)

       STR = '@<neutrality> [/m$+3$=] vs t@'
       CALL TXGRFTX(1, GTX, GTY(0,46), NGTM, NGT, 1, STR, IND)

       STR = '@TAUP, TAUPA [s] vs t@'
       CALL TXGRFTX(2, GTX, GTY(0,49), NGTM, NGT, 2, STR, IND)

       STR = '@Gamma_a [10$+20$=/m$+2$=s$+1$=] vs t@'
       CALL TXGRFTX(3, GTX, GTY(0,51), NGTM, NGT, 1, STR, IND)

       STR = '@TNBcol [Nm] vs t@'
       CALL TXGRFTX(4, GTX, GTY(0,52), NGTM, NGT, 1, STR, IND)

       STR = '@TTqt [Nm] vs t@'
       CALL TXGRFTX(5, GTX, GTY(0,53), NGTM, NGT, 1, STR, IND)

    CASE DEFAULT
       WRITE(6,*) 'Unknown NGYT: NGYT = ',NGYT
    END SELECT

    CALL PAGEE

    RETURN
  END SUBROUTINE TXGRFT

  !***************************************************************
  !
  !   SLAVE ROUTINE TO Write graph of GVX, GVY
  !
  !***************************************************************

  SUBROUTINE TXGRFVX(K, GTXL, GTYL, NGTM, NGTL, NG, STR, MODE, IND, GYMIN)

    INTEGER(4), INTENT(IN) :: K, NGTM, NGTL, NG
    REAL(4), DIMENSION(0:NGTM), INTENT(IN) :: GTXL
    REAL(4), DIMENSION(0:NGTM,1:NG), INTENT(IN) :: GTYL
    real(4), intent(in), optional :: GYMIN
    character(len=*), INTENT(IN) ::  STR
    INTEGER(4) :: MODE, IND
    REAL(4), DIMENSION(4) :: GPXY

    GPXY(1) =  3.0 + 12.5 * MOD(K,2)
    GPXY(2) = 12.5 + 12.5 * MOD(K,2)
    GPXY(3) = 10.5 -  8.5 * REAL(K/2)
    GPXY(4) = 17.0 -  8.5 * REAL(K/2)
    IF(PRESENT(GYMIN)) THEN
       CALL TXGRAF(GPXY, GTXL, GTYL, NGTM+1, NGTL+1, NG, &
            &            GTXL(0), GTXL(NGTL), STR, 0.3, MODE, IND, 0, GYMIN_IN=GYMIN)
    ELSE
       CALL TXGRAF(GPXY, GTXL, GTYL, NGTM+1, NGTL+1, NG, &
            &            GTXL(0), GTXL(NGTL), STR, 0.3, MODE, IND, 0)
    END IF

    RETURN
  END SUBROUTINE TXGRFVX

  !***************************************************************
  !
  !   SLAVE ROUTINE TO Write graph of GTX, GTY
  !
  !***************************************************************

  SUBROUTINE TXGRFTX(K, GTXL, GTYL, NGTM, NGTL, NG, STR, IND)

    INTEGER(4), INTENT(IN) :: K, NGTM, NGTL, NG
    REAL(4), DIMENSION(0:NGTM), INTENT(IN) :: GTXL
    REAL(4), DIMENSION(0:NGTM,1:NG), INTENT(IN) :: GTYL
    character(len=*), INTENT(IN) ::  STR
    INTEGER(4) :: IND
    REAL(4), DIMENSION(4) :: GPXY

    GPXY(1) =  3.0 + 12.0 * MOD(K,2)
    GPXY(2) = 12.0 + 12.0 * MOD(K,2)
    GPXY(3) = 14.0 -  4.3 * REAL(K/2)
    GPXY(4) = 17.0 -  4.3 * REAL(K/2)
    CALL TXGRAF(GPXY, GTXL, GTYL, NGTM+1, NGTL+1, NG, &
         &            GTXL(0), GTXL(NGTL), STR, 0.3, 1, IND, 0)

    RETURN
  END SUBROUTINE TXGRFTX

  !***************************************************************
  !
  !   Graph of each TERM
  !
  !***************************************************************

  SUBROUTINE TXGRFQ(NQ,ID)

    use tx_commons, only : NRMAX, NQM, NCM, NLCMAX, rhob, NRC
    use tx_interface, only : APTOS

    INTEGER(4), INTENT(IN) :: NQ, ID
    INTEGER(4) :: NC, NSTR, IND
    REAL(4) :: GXMAX
    REAL(4), DIMENSION(0:NRMAX,1:NCM) :: GQYL
    REAL(4), DIMENSION(1:4) :: GPXY
    REAL(4), DIMENSION(1:4,1:5) :: GPXYA
    character(len=80), DIMENSION(1:NQM) :: STRGQ
    character(len=80) :: STR

    !        Title
    DATA STRGQ /'$#f$#$',"$#y$#$-t$='","$#y$#'",'$#y$#','$#y$#$-t$=',  &
         &      'n$-e$=','n$-e$=u$-e$=$+V$=', &
         &      'Bu$-e//$=','n$-e$=<Ru$-e$#f$#$=>','n$-e$=T$-e$=', &
         &      'Bq$-e//$=', 'n$-e$=<u$-e$#f$#$=/R>', &
         &      'n$-i$=','n$-i$=u$-i$=$+V$=', &
         &      'Bu$-i//$=','n$-i$=<Ru$-i$#f$#$=>','n$-i$=T$-i$=', &
         &      'Bq$-i//$=', 'n$-i$=<u$-i$#f$#$=/R>',&
         &      'n$-b$=','n$-b$=u$-b$=$+V$=','Bu$-b//$=', &
         &      'n$-b$=<Ru$-b$#f$#$=>','n$-b$=<u$-b$#f$#$=/R>', &
         &      'Slow n$-0$=', 'Thermal n$-0$=', 'Halo n$-0$=', 'Ripple n$-b$='/

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
    CALL APTOS(STR,NSTR,NQ)
    CALL APTOS(STR,NSTR, ': ',2)
    CALL APTOS(STR,NSTR, STRGQ(NQ), LEN_TRIM(STRGQ(NQ)))

    GQYL(0:NRMAX,1:NLCMAX(NQ)) = GQY(0:NRMAX,1:NLCMAX(NQ),NQ)

    IF (MODEG == 2) THEN
       IND = 9
    ELSE
       IND = 0
    END IF
    GXMAX = REAL(rhob)
    CALL TXGRAF(GPXY, GX, GQYL, NRMAX+1, NRMAX+1, NLCMAX(NQ), &
         &      0.0, GXMAX, '@'//STR(1:NSTR)//'@', 0.3, 2, IND, 0)
!         &      0.0, GXMAX, '@'//STR(1:NSTR)//'@', 0.3, 4, IND, 0)
    write(6,'(A3,I2,X,A4,6X,A3,11X,A7)') "EQ=",NQ,"term","sum","nrc"
    do nc=1,nlcmax(nq)
       write(6,'(5X,I3,1P1E15.5,3X,1E15.5)') nc,sum(gqyl(0:nrmax,nc)),gqyl(nrc,nc)
    end do

    RETURN
  END SUBROUTINE TXGRFQ

  !***************************************************************
  !
  !   Write parameter to graphic screen
  !
  !***************************************************************

  SUBROUTINE TXWPGR

    use tx_commons, only : SLID, PNBCD, BB, rIp, FSDFIX, FSANOM, FSBOHM, FSPCLD, FSPCLC, &
         &              PROFD, PROFC, FSCX, FSRP, FSLC, FSNC, FSLP, FSION, FSD02, &
         &              PNBHT1, PNBHT2, PNBHP, PNBHex, PRFHe, PRFHi, Vb, De0, rMue0, rMui0, &
         &              Chie0, Chii0, PTe0, PTea, PTi0, PTia, V0, rGamm0, rGASPF, Zeff
    INTEGER(4) :: IFNT

    CALL INQFNT(IFNT)
!    CALL SETFNT(32)
    CALL SETFNT(0)
    CALL SETCHS(0.3, 0.0)
    CALL SETLIN(0, 1, 7)

    GXM = 13.0 + 1.0
    GYM =  8.5 + 0.2
    GYS =  0.5
    NP  = 0

    CALL TXWPS('+'//SLID//'+')
    CALL TXWPS('@PNBCD @', PNBCD)
!!!    CALL TXWPS('@NRMAX @', NRMAX)

!    NP = NP + 1
    CALL TXWPS('@BB    @', BB)
    CALL TXWPS('@rIp   @', rIp)
    CALL TXWPS('@FSDFX1@', FSDFIX(1))
    CALL TXWPS('@FANOM3@', FSANOM(3))
    CALL TXWPS('@FSBOHM@', FSBOHM)
    CALL TXWPS('@FSPCLD@', FSPCLD)
    CALL TXWPS('@FSPCLC@', FSPCLC)
    CALL TXWPS('@PROFD @', PROFD)
    CALL TXWPS('@PROFC @', PROFC)
    CALL TXWPS('@FSCX  @', FSCX)
    CALL TXWPS('@FSRP  @', FSRP)
    CALL TXWPS('@FSLC  @', FSLC)
    CALL TXWPS('@FSNC  @', FSNC)
    CALL TXWPS('@FSLP  @', FSLP)
    CALL TXWPS('@FSION @', FSION)
    CALL TXWPS('@FSD02 @', FSD02)

    GXM = GXM + 0.35 * 17
    NP = 0
!!!    CALL TXWPS('@PNBH  @', PNBH)
    CALL TXWPS('@PNBHT @', PNBHT1+PNBHT2)
    IF(PNBHex == 0.D0) THEN
       CALL TXWPS('@PNBHP @', PNBHP)
    ELSE
       CALL TXWPS('@PNBHex@', PNBHex)
    END IF
    CALL TXWPS('@PRFH  @', PRFHe+PRFHi)
    CALL TXWPS('@Vb    @', Vb)
    CALL TXWPS('@De0   @', De0)
!!!    CALL TXWPS('@Di0   @', Di0)
    CALL TXWPS('@rMue0 @', rMue0)
    CALL TXWPS('@rMui0 @', rMui0)
    CALL TXWPS('@Chie0 @', Chie0)
    CALL TXWPS('@Chii0 @', Chii0)
!!!    CALL TXWPS('@WPM0  @', WPM0)
    CALL TXWPS('@PTe0  @', PTe0)
    CALL TXWPS('@PTea  @', PTea)
    CALL TXWPS('@PTi0  @', PTi0)
    CALL TXWPS('@PTia  @', PTia)
    CALL TXWPS('@V0    @', V0)
    CALL TXWPS('@rGamm0@', rGamm0)
    CALL TXWPS('@rGASPF@', rGASPF)
    CALL TXWPS('@Zeff  @', Zeff)

    CALL SETFNT(IFNT)
    RETURN
  END SUBROUTINE TXWPGR

  !***************************************************************
  !
  !   WPgr's Sub : write Double
  !
  !***************************************************************

  SUBROUTINE TXWPSD(STR, VAL)

    REAL(8), INTENT(IN) :: VAL
    character(len=*), INTENT(IN) :: STR

    CALL MOVE(GXM, GYM - GYS * NP)
    CALL TEXTX(STR)
    CALL NUMBD(VAL,'(1PD9.2)', 9)
    NP = NP + 1

    RETURN
  END SUBROUTINE TXWPSD

  !***************************************************************
  !
  !   WPgr's Sub : write Integer
  !
  !***************************************************************

  SUBROUTINE TXWPSI(STR, IVAL)

    INTEGER(4), INTENT(IN) :: IVAL
    character(len=*), INTENT(IN) :: STR

    CALL MOVE(GXM, GYM - GYS * NP)
    CALL TEXTX(STR)
    CALL TEXT(' = ', 3)
    CALL NUMBI(IVAL,'(I6)',6)
    NP = NP + 1

    RETURN
  END SUBROUTINE TXWPSI

  !***************************************************************
  !
  !   WPgr's Sub : write Strings
  !
  !***************************************************************

  SUBROUTINE TXWPSS(STR)

    character(len=*), INTENT(IN) :: STR

    CALL MOVE(GXM, GYM - GYS * NP)
    CALL TEXTX(STR)
    NP = NP + 1

    RETURN
  END SUBROUTINE TXWPSS

  !***************************************************************
  !
  !   Draw graph
  !
  !***************************************************************

  SUBROUTINE TXGRAF(GPXY, GX, GY, NXM, NXMAX, NGMAX, &
       &                  GXMIN, GXMAX, STR, FONT, MODE, IND, ILOG, &
       &                  GYMAX_IN, GYMIN_IN)!, KINDEX)

    INTEGER(4), INTENT(IN) :: NXM, NXMAX, NGMAX, MODE, IND
    REAL(4), INTENT(IN) :: GXMIN, GXMAX, FONT
    REAL(4), DIMENSION(4), INTENT(IN) :: GPXY
    REAL(4), DIMENSION(:), INTENT(IN) :: GX
    REAL(4), DIMENSION(:,:), INTENT(IN) :: GY
    integer(4), intent(in) :: ILOG
    REAL(4), INTENT(IN), OPTIONAL :: GYMAX_IN, GYMIN_IN
    character(len=*), INTENT(IN) :: STR
!    character(len=*), INTENT(IN), optional :: KINDEX

    INTEGER(4) :: IFNT, NGV, NGULEN, ICL, IPAT, IMRK, ISTEP, NG
    REAL(4) :: GX1, GX2, GY1, GY2, gSLEN, GSXMIN, GSXMAX, GXSTEP, &
         &        GYMIN, GYMAX, GSYMIN, GSYMAX, GYSTEP, GYORG,  &
         &        GMRK, GCHH, GXL, GYL
    INTEGER(4), DIMENSION(0:4) :: NLTYPE
    real(4), dimension(:,:), allocatable :: GYAR
    DATA NLTYPE/0,2,3,4,6/

    IF (MODE < 0 .OR. MODE > 4) THEN
       WRITE(6,*) '### ERROR(TXGRAF) : MODE = ', MODE
       RETURN
    END IF

    GX1=GPXY(1)
    GX2=GPXY(2)
    GY1=GPXY(3)
    GY2=GPXY(4)

    CALL INQFNT(IFNT)

    IF (IND == 0) THEN
       gSLEN = 1.0
    ELSE
       gSLEN = 0.2
    END IF

    CALL SETFNT(32)
    CALL SETCHS(FONT, 0.0)
    CALL SETLIN(0, 0, 7)
    CALL GTEXTX(GX1,GY2+0.2,STR,0)

    CALL SETFNT(44)
    CALL SETCHS(FONT, 0.0)
    CALL SETLIN(0, 0, 7)
    CALL SETLNW(0.035)

    allocate(GYAR(size(GY,1),size(GY,2)))
    GYAR = GY
!    IF(PRESENT(GYMAX_IN)) where(GYAR > GYMAX_IN) GYAR = GYMAX_IN
!    IF(PRESENT(GYMIN_IN)) where(GYAR < GYMIN_IN) GYAR = GYMIN_IN

    CALL GQSCAL(GXMIN, GXMAX, GSXMIN, GSXMAX, GXSTEP)
    GSXMIN = GXMIN
    GSXMAX = GXMAX
    CALL GMNMX2(GYAR,NXM,1,NXMAX,1,1,NGMAX,1,GYMIN,GYMAX)
    IF(ILOG == 0) THEN
       IF (GYMAX > 0.0) THEN
          IF (GYMIN > 0.0) GYMIN=0.0
       ELSE
          GYMAX=0.0
       END IF
    END IF
!    IF(PRESENT(KINDEX)) THEN
       IF(PRESENT(GYMAX_IN)) GYMAX=GYMAX_IN
       IF(PRESENT(GYMIN_IN)) GYMIN=GYMIN_IN
!    END IF
    CALL GQSCAL(GYMIN, GYMAX, GSYMIN, GSYMAX, GYSTEP)
    IF (GSYMIN > GYMIN) GSYMIN = GSYMIN - GYSTEP
    IF (GSYMAX < GYMAX) GSYMAX = GSYMAX + GYSTEP
    IF(GYMIN * GYMAX <= 0.0) THEN
       GYORG = 0.0
    ELSE
       GYORG = GSYMIN
    ENDIF

    CALL GDEFIN(GX1, GX2, GY1, GY2, GSXMIN, GSXMAX, GSYMIN, GSYMAX)
    CALL GFRAME
    CALL SETLNW( 0.017)
    CALL GSCALE(GSXMIN, GXSTEP, 0.0, 0.0, gSLEN, IND)
!!$    IF (GXSTEP < 0.01) THEN
!!$       NGV=NGULEN(GXSTEP*5)
!!$       IF (NGV < 0) THEN
!!$          NGV=NGV-3200
!!$       ELSE
!!$          NGV=NGV+3200
!!$       END IF
!!$       CALL GVALUE(GSXMIN, GXSTEP*5, 0.0, 0.0, NGV)
       NGV=NGULEN(GXSTEP*2)
       CALL GVALUE(GSXMIN, GXSTEP*2, 0.0, 0.0, NGV)
!!$    ELSE
!!$       NGV=NGULEN(GXSTEP*2)
!!$       IF (NGV < 0) THEN
!!$          NGV=NGV-3200
!!$       ELSE
!!$          NGV=NGV+3200
!!$       END IF
!!$       CALL GVALUE(GSXMIN, GXSTEP*2, 0.0, 0.0, NGV)
!!$    END IF

    ! Semi-Log or not

    IF(ILOG == 0) THEN
       CALL GSCALE(0.0, 0.0, GYORG, GYSTEP, gSLEN, IND)
       CALL SETLNW(-0.017)
       IF (GSYMIN < 0.0 .AND. GSYMAX > 0.0) &
            &     CALL GSCALE(0.0, 0.0, 0.0, GSYMAX-GSYMIN,  2.0, 0)
       CALL GVALUE(0.0,0.0,GYORG,GYSTEP*2,NGULEN(GYSTEP*2))
    ELSE
       CALL GSCALL(0.0, 0.0, GYORG, 4, gSLEN, IND)
       CALL SETLNW(-0.017)
       IF (GSYMIN < 0.0 .AND. GSYMAX > 0.0) &
            &     CALL GSCALL(0.0, 0.0, 0.0, 3,  2.0, 0)
       CALL GVALUL(0.0,0.0,GYORG,1,NGULEN(GYSTEP*2))
    END IF

    ! MODE = 0: Change Line Color (Last Color Fixed)

    SELECT CASE(MODE)
    CASE (0)
       DO NG = 1, NGMAX
          ICL = 7 - MOD(NGMAX - NG, 5)
          CALL SETLIN(0, 1, ICL)
          CALL GPLOTP(GX, GYAR(1,NG), 1, NXMAX, 1, 0, 0, 0)
       END DO

       ! MODE = 1: Change Line Color and Style
       ! MODE = 2: Change Line Color and Style (With Legend)

    CASE (1:2)
       IF (MODE == 1) THEN
          DO NG = 1, NGMAX
             ICL  = 7 - MOD(NG-1, 5)
             IPAT = NLTYPE(MOD(NG-1, 5))
             CALL SETLIN(0, 1, ICL)
             CALL GPLOTP(GX, GYAR(1,NG), 1, NXMAX, 1, 0, 0, IPAT)
          END DO
       ELSE
          DO NG = 1, NGMAX
             ICL = 7 - MOD(NG - 1, 5)
             ISTEP = NXMAX / 10
             IPAT  = (NG - 1) / 5
             CALL SETLIN(0, 1, ICL)
             CALL GPLOTP(GX, GYAR(1,NG), 1, NXMAX, 1, 0, ISTEP, IPAT)
          END DO
       END IF
       ! Legend
       IF (MODE == 2) THEN
          CALL SETCHS(0.25,0.0)
          GCHH = 0.25
          GXL = GX2 + GCHH
          GYL = GY2 - GCHH
          DO NG = 1, NGMAX
             ICL = 7 - MOD(NG - 1, 5)
             CALL SETLIN(0, 1, ICL)
             IPAT = (NG - 1) / 5
             CALL MOVEPT(GXL + GCHH * 3.0, GYL + GCHH / 2.0, IPAT)
             CALL DRAWPT(GXL + GCHH * 7.0, GYL + GCHH / 2.0)
             CALL MOVE(GXL, GYL)
             CALL NUMBI(NG, '(I2)', 2)
             GYL = GYL - GCHH * 2.0
          END DO
       END IF
       ! End of Legend

       ! MODE = 3: Change Line Color, Style and Mark
       ! MODE = 4: Change Line Color, Style and Mark (With Legend)

    CASE (3:4)
       IMRK = 0
       GMRK = 0.3
       CALL SETMKS(IMRK, GMRK)
       DO NG = 1, NGMAX
          ICL = 7 - MOD(NG - 1, 5)
          IMRK  = MOD(NG - 1, 5) + 1
          ISTEP = NXMAX / 10
          IPAT  = (NG - 1) / 5
          CALL SETLIN(0, 1, ICL)
          CALL GPLOTP(GX, GYAR(1,NG), 1, NXMAX, 1, -IMRK, ISTEP, IPAT)
       END DO
       ! Legend
       IF (MODE == 4) THEN
          CALL SETCHS(0.25,0.0)
          GCHH = 0.25
          GXL = GX2 + GCHH
          GYL = GY2 - GCHH
          DO NG = 1, NGMAX
             ICL = 7 - MOD(NG - 1, 5)
             CALL SETLIN(0, 1, ICL)
             IPAT = (NG - 1) / 5
             CALL MOVEPT(GXL + GCHH * 3.0, GYL + GCHH / 2.0, IPAT)
             CALL DRAWPT(GXL + GCHH * 7.0, GYL + GCHH / 2.0)
             CALL MOVE(GXL, GYL)
             CALL NUMBI(NG, '(I2)', 2)
             IMRK = MOD(NG - 1, 5) + 1
             IF (IMRK /= 0) THEN
                CALL SETMKS(IMRK, GMRK)
                CALL MARK(GXL + GCHH * 5.0, GYL + GCHH / 2.0)
             END IF
             GYL = GYL - GCHH * 2.0
          END DO
       END IF
       ! End of Legend
       IMRK = 0
       GMRK = 0.2
       CALL SETMKS(IMRK, GMRK)
    END SELECT

    CALL SETFNT(IFNT)

    deallocate(GYAR)

    RETURN
  END SUBROUTINE TXGRAF

  !***************************************************************
  !
  !   Appropriate GY and STR for graphic routine
  !
  !***************************************************************

  SUBROUTINE APPROPGY(MODE, GIN, GOUT, STR, NXMAX, NYMAX, gDIV, GIN1, GMAX, GMIN)
    use tx_interface, only : APTOS

    INTEGER(4), INTENT(IN) :: MODE, NXMAX, NYMAX
    REAL(4), INTENT(IN) :: gDIV
    real(4), intent(out), optional :: GMAX, GMIN
    REAL(4), DIMENSION(0:NXMAX,0:NYMAX),INTENT(IN)  :: GIN
    REAL(4), DIMENSION(0:NXMAX,0:NYMAX),INTENT(OUT) :: GOUT
    REAL(4), DIMENSION(0:NXMAX),INTENT(IN), optional :: GIN1
    character(len=*), INTENT(INOUT) :: STR

    INTEGER(4) :: NSTR, POSAT
    REAL(4) :: gDIVL

    gDIVL = 1.0
    IF(MODE /= 0) gDIVL = gDIV

    ! Append gDIV to string for showing multiplication factor
    IF (gDIVL /= 1.0) THEN
       POSAT = INDEX(STR,'@',.TRUE.)
       IF (POSAT /= 0) STR = STR(1:POSAT-1)
       NSTR = LEN_TRIM(STR)+1
       CALL APTOS(STR, NSTR, ' [', 2)
       CALL APTOS(STR, NSTR, gDIVL, 'E0')
       CALL APTOS(STR, NSTR, ']@', 2)
    END IF

    if(present(GIN1)) then
       GOUT(0:NXMAX,0) = GIN1(0:NXMAX) / gDIVL
       GOUT(0:NXMAX,1:NYMAX) = GIN(0:NXMAX,0:NYMAX-1) / gDIVL
    else
       GOUT(0:NXMAX,0:NYMAX) = GIN(0:NXMAX,0:NYMAX) / gDIVL
    end if

    if(present(GMAX)) GMAX = maxval(GOUT)
    if(present(GMIN)) GMIN = minval(GOUT)

    RETURN
  END SUBROUTINE APPROPGY

  !***********************************************************
  !
  !   Ceiling FUNCTION for LOG10 plot
  !
  !***********************************************************

  REAL(4) FUNCTION GLOG(X,XMIN,XMAX)

    implicit none
    real(4) :: GUCLIP
    real(8), intent(in) :: X, XMIN, XMAX
    real(8) :: PLOG

    GLOG = GUCLIP(PLOG(X,XMIN,XMAX))

    RETURN
  END FUNCTION GLOG

  !***********************************************************
  !
  !   WRITE RAW VALUE ONTO CONSOLE
  !
  !***********************************************************

  subroutine write_console(n,char)

    use tx_commons, only : NRMAX, R
    use tx_interface, only : dfdx
    integer(4), intent(in) :: n
    character(len=1), intent(in) :: char
    integer(4) :: i
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
             write(6,*) gx(i),gy(i,ngr,n)
          end do
       end IF
    case('T')
       if (n >= 0 .and. n <= NGYTM) then
          write(6,'(A1,6X,A1,13X,A5)') "#","T","Value"
          do i = 0, ngt
             write(6,*) gtx(i),gty(i,n)
          end do
       end IF
    case('V')
       if (n >= 0 .and. n <= NGYVM) then
          write(6,'(A1,6X,A1,13X,A5)') "#","T","Value"
          do i = 0, ngvv
             write(6,*) gvx(i),gvy(i,n)
          end do
       end IF
    case('D') ! R-derivative of a radial profile
       if (n >= 0 .and. n <= NGYRM) then
          allocate(dVdr(0:nrmax))
          dVdr(0:nrmax) = dfdx(R,dble(gy(0:nrmax,ngr,n)),nrmax,0)
          write(6,'(A1,6X,A3,11X,A5)') "#","rho","dValue/dR"
          do i = 0, nrmax
             write(6,*) gx(i),real(dVdr(i))
          end do
          deallocate(dVdr)
       end IF
    case default
       WRITE(6,*) '### ERROR : Invalid Command : '
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

end module tx_graphic
