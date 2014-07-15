  subroutine GSOPEN
  end subroutine GSOPEN

  subroutine GSCLOS
  end subroutine GSCLOS

  subroutine GRMODE

  end subroutine GRMODE

  subroutine CPLOT3D1(i, r)
    implicit none
    integer :: i
    external :: r
  end subroutine CPLOT3D1

  subroutine R2W2B(r, r3)
    real(4) :: r
    real(4), dimension(3) :: r3
  end subroutine R2W2B

  subroutine W2G2B(r, r3)
    real(4) :: r
    real(4), dimension(3) :: r3
  end subroutine W2G2B

  subroutine R2Y2W(r, r3)
    real(4) :: r
    real(4), dimension(3) :: r3
  end subroutine R2Y2W

  subroutine R2G2B
  end subroutine R2G2B

  subroutine WRGBW
  end subroutine WRGBW

  subroutine GTEXT(r1, r2, c, i)
    character(len=80) :: c
    real(8) :: r1, r2
    integer :: i
  end subroutine GTEXT

  subroutine TEXT(ktext, nchar)
    character(len=256) :: ktext
    integer :: nchar
  end subroutine TEXT

  integer(4) function NGULEN(gs)
    real(4) :: gs
    NGULEN = 0
  end function NGULEN

  subroutine GPLOTP(r1, r2, i1, i2, i3, i4, i5, i6)
    real(4) :: r1, r2
    integer :: i1, i2, i3, i4, i5, i6
  end subroutine GPLOTP

  subroutine GAXIS3D(i)
    integer :: i
  end subroutine GAXIS3D

  subroutine GDEFIN3D(GX1,GX2,GY1,GY2,GXL,GYL,GZL)
    real(4) :: GX1, GX2, GY1, GY2, GXL, GYL, GZL
  end subroutine GDEFIN3D

  subroutine GDEFIN(GX1,GX2,GY1,GY2,GSXMIN,GSXMAX,GSYMIN,GSYMAX)
    real(4) :: GX1, GX2, GY1, GY2, GSXMIN, GSXMAX, GSYMIN, GSYMAX
  end subroutine GDEFIN

  subroutine GQSCAL(GXMIN,GXMAX,GSXMIN,GSXMAX,GSTEPX)
    real(4) :: GXMIN, GXMAX, GSXMIN, GSXMAX, GSTEPX
  end subroutine GQSCAL

  subroutine GPLOTPE(GXE,GEUN,GEUP,i1,NXEMAX,i2,r)
    integer :: i1, i2, nxemax
    real(4) :: r
    real(4), dimension(NXEMAX) :: GXE, GEUN, GEUP
    
  end subroutine GPLOTPE

  subroutine GSCALE3DX(r1, r2, r3, i)
    real(4) :: r1, r2, r3
    integer :: i
  end subroutine GSCALE3DX

  subroutine GSCALE3DY(r1, r2, r3, i)
    real(4) :: r1, r2, r3
    integer :: i
  end subroutine GSCALE3DY

  subroutine GSCALE3DZ(r1, r2, r3, i)
    real(4) :: r1, r2, r3
    integer :: i
  end subroutine GSCALE3DZ

  SUBROUTINE GVALUL(XORG,NXSTEP,YORG,NYSTEP,IND)
    real(4) :: XORG, YORG
    integer :: NXSTEP, NYSTEP, IND
  end SUBROUTINE GVALUL

  SUBROUTINE GVALUE(XORG,XSTEP,YORG,YSTEP,IND)
    real(4) :: XORG, YORG, XTEP, YSTEP
    integer :: IND
  end SUBROUTINE GVALUE

  subroutine GUFLSH

  end subroutine GUFLSH

  subroutine PERSE1(GZ,NXM,NXMAX,NYMAX,GZMIN,GZMAX,IXY,IND,XL,YL,ZL,A,B,C,D,E)
    real(4) :: GZ, GZMIN, GZMAX
    integer :: NXM, NXMAX, NYMAX, IXY, IND
    real(8) :: XL, YL, ZL, A, B, C, D, E
  end subroutine PERSE1

  subroutine NUMBR(r, kform, nchar)
    real(4) :: r
    character(len=16) :: kform
    integer :: nchar
  end subroutine NUMBR

  subroutine NUMBD(r, kform, nchar)
    real(8) :: r
    character(len=16) :: kform
    integer :: nchar
  end subroutine NUMBD

  subroutine CHMODE
  end subroutine CHMODE

  subroutine PAGES
  end subroutine PAGES

  subroutine PAGEE
  end subroutine PAGEE

  subroutine SETRGB(r1, r2, r3)
    real(4) :: r1, r2, r3
  end subroutine SETRGB

  subroutine GUCPTL(k)
    character(len=1) :: k
  end subroutine GUCPTL

  subroutine SETLNW(r)
    real(4) :: r
  end subroutine SETLNW

  subroutine GMNMX1(GX, i1, NXMAX, i2, GXMIN, GXMAX)
    real(4), dimension(*) :: GX
    integer :: i1, nxmax, i2
    real(4) :: gxmin, gxmax
  end subroutine GMNMX1

  subroutine GMNMX2(GY,NXM,i1,NXMAX,i2,i3,NGMAX,i4,GYMIN,GYMAX)
    integer :: nxm, nxmax, ngmax, i1, i2, i3, i4
    real(4), dimension(nxm, ngmax) :: gy
    real(4) :: gymin, gymax
  end subroutine GMNMX2

  subroutine SETLIN(i1, i2, i3)
    integer :: i1, i2, i3
  end subroutine SETLIN

  subroutine  BLACK(GFACT,r3)
    real(4) :: gfact
    real(4), dimension(3) :: r3
  end subroutine BLACK
    
  subroutine SETFNT(i)
    integer :: i
  end subroutine SETFNT

  real(4) function GUCLIP(r)
    real(4) :: r
    guclip = 0.
  end function GUCLIP

  subroutine SETCHS(r1, r2)
    real(4) :: r1, r2
  end subroutine SETCHS

  subroutine GSCALE(r1,r2,r3,r4,r5,i)
    real(4) :: r1, r2, r3, r4, r5
    integer :: i
  end subroutine GSCALE

  subroutine MOVE(r1, r2)
    real(4) :: r1, r2
  end subroutine MOVE

  subroutine GUTIME(T)
    real(4) :: T
  end subroutine GUTIME

  subroutine GFRAME
  end subroutine GFRAME

  SUBROUTINE CONTP1(Z,NXA,NXMAX,NYMAX, &
       ZORG,ZSTEP,NSTEP,IPRD,IPAT,KA)
    integer :: nxa, nxmax, nymax, nstep, iprd, ipat
    real(4), dimension(nxa, *) :: z
    real(4), dimension(2,*) :: ka
    real(4) :: zorg, zstep
  end SUBROUTINE CONTP1

  SUBROUTINE CONTP0(Z,X,Y,NXA,NXMAX,NYMAX, &
       ZORG,ZSTEP,NSTEP,IPRD,IPAT,KA,SUB)
    integer :: nxa, nxmax, nymax, nstep, iprd, iapt
    real(4) :: z(nxa, *), x(*), y(*)
    integer :: ka(2,*)
    real(4) :: zorg, zstep
    external sub
  end SUBROUTINE CONTP0

  SUBROUTINE CONTR0(Z,X,Y,NXA,NX,NY,ZORG,ZSTEP,NSTEP,IPRD,SUB)
    integer :: nxa, nx, ny, nstep, iprd
    real(4) :: z(nxa,*), x(*), y(*)
    real(4) :: zorg, zstep
    external :: sub
  end SUBROUTINE CONTR0

  SUBROUTINE CONTQ0(Z,X,Y,NXA,NXMAX,NYMAX, &
       ZORG,ZSTEP,NSTEP,IPRD,IPAT,KA,SUB)
    integer :: nxa, nxmax, nymax, nstep, iprd, ipat
    real(4) :: z(nxa,*), x(*), y(*)
    integer :: ka(2,*)
    real(4) :: zorg, zstep
    external :: sub
  end SUBROUTINE CONTQ0

  SUBROUTINE CONTQ1(Z,NXA,NXMAX,NYMAX, &
       ZORG,ZSTEP,NSTEP,IPRD,IPAT,KA)
    integer :: nxa, nxmax, nymax, nstep, iprd, ipat
    real(4) :: z(nxa,*)
    integer :: ka(2,*)
    real(4) :: zorg, zstep
  end SUBROUTINE CONTQ1

  SUBROUTINE CONTG1(Z,NXA,NXMAX,NYMAX, &
       ZL,RGB,ILN,WLN,NSTEP,ISPL,IPRD,KA)
    integer :: nxa, nxmax, nymax, nstep, ispl, iprd
    real(4) :: z(nxa, nymax), zl(nstep), rgb(3,nstep)
    integer :: iln(nstep), wln(nstep), ka(2,nxmax*nymax)
  end SUBROUTINE CONTG1

  SUBROUTINE CONTF1(Z,NXA,NXMAX,NYMAX,ZL,RGB,NSTEP,IPRD)
    integer :: nxa, nxmax, nymax, nstep, iprd
    real(4) :: z(nxa, nymax), zl(nstep), rgb(3,0:nstep)
  end SUBROUTINE CONTF1

  SUBROUTINE GSCALL(XORG,NXSTEP,YORG,NYSTEP,SLENG,IND)
    integer :: nxstep, nystep, ind
    real(4) :: xorg, yorg, sleng
  end SUBROUTINE GSCALL

  SUBROUTINE GDATA3D1(Z,NXA,NX1,NY1,XMIN1,XMAX1,YMIN1,YMAX1, &
       ZMIN1,ZMAX1)
    integer :: nxa, nx1, ny1
    real(4) :: z(nxa,*)
    real(4) :: xmin1, xmax1, ymin1, ymax1
  end SUBROUTINE GDATA3D1

  SUBROUTINE GVIEW3D(PHI1,THETA1,RADIUS1,R1,IAXIS1,OX1,OY1,OZ1)
    real(4) :: phi1, theta1, radius1, r1, ox1, oy1, oz1
    integer :: iaxis1
  end SUBROUTINE GVIEW3D

  subroutine DRAWPT(x, y)
    real(4) :: x, y
  end subroutine DRAWPT

  SUBROUTINE GUDATE(NDY,NDM,NDD,NTH,NTM,NTS)
    integer :: ndy, ndm, ndd, nth, ntm, nts
  end SUBROUTINE GUDATE

  SUBROUTINE GVALUE3DX(XORG,DPX,IPOSV,IND)
    real(4) :: xorg, dpx
    integer :: iposv, ind
  end SUBROUTINE GVALUE3DX

  SUBROUTINE GVALUE3DY(XORG,DPX,IPOSV,IND)
    real(4) :: xorg, dpx
    integer :: iposv, ind
  end SUBROUTINE GVALUE3DY

  SUBROUTINE GVALUE3DZ(XORG,DPX,IPOSV,IND)
    real(4) :: xorg, dpx
    integer :: iposv, ind
  end SUBROUTINE GVALUE3DZ

  SUBROUTINE MOVEPT(X,Y,IPAT)
    real(4) :: x, y
    integer :: ipat
  end SUBROUTINE MOVEPT

  SUBROUTINE INQGDEFIN(PXMIN,PXMAX,PYMIN,PYMAX, &
       GXMIN,GXMAX,GYMIN,GYMAX)
    real(4) :: pxmin, pxmax, pymin, pymax, gxmin, gxmax, gymin, gymax
  end SUBROUTINE INQGDEFIN

  subroutine EQGOUT(MODE)
    integer :: MODE
  end subroutine EQGOUT
