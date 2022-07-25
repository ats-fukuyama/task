!**************************************************************
!
!   GRAPHIC 3D : UNIVERSAL ROUTINE
!
!**************************************************************

subroutine TXGRUR(GX,GTX,GYL,NRMAX,NGT,NGTM)!,STR,KV,INQ)

  implicit none
  integer(4),                  intent(in) :: NRMAX, NGT, NGTM!, INQ
  real(4), dimension(0:NRMAX), intent(in) :: GX
  real(4), dimension(0:NGT),   intent(in) :: GTX
  real(4), dimension(0:NRMAX,0:NGTM), intent(in) :: GYL
!  CHARACTER(LEN=80),intent(IN):: STR, KV
  integer(4) :: INQ
  character(LEN=80) :: STR, KV
  real(4) :: GX1, GX2, GY1, GY2

  STR = '@3D Graphic@'
  KV  = ''
  INQ = 0

  GX1=3.0
  GX2=18.0
  GY1=2.0
  GY2=17.0

  call PAGES
  call TXGR3D(GX1,GX2,GY1,GY2,GX,GTX,GYL,NRMAX+1,NRMAX+1,NGT,STR,KV,2+INQ)
  call PAGEE

  
end subroutine TXGRUR

!***********************************************************
!
!   SUBPROGRAM FOR 3D PROFILE FOR GSAF
!
!***********************************************************

subroutine TXGR3D(GX1,GX2,GY1,GY2,GX,GY,GZ,NXM,NXMAX,NYMAX,STR,KV,MODE)

  implicit none
  real(4),        intent(IN):: GX1,GX2,GY1,GY2
  integer(4),     intent(IN):: NXM,NXMAX,NYMAX,MODE
  real(4),dimension(NXMAX),      intent(IN):: GX
  real(4),dimension(NYMAX),      intent(IN):: GY
  real(4),dimension(NXM,NYMAX),  intent(IN):: GZ
  character(LEN=80),             intent(IN):: STR,KV
  integer(4) :: I, NGULEN
  real(4)    :: GOX, GOY, GOZ, GPHI, GRADIUS, GSTEPX, GSTEPY, GSTEPZ, GSXMAX, GSXMIN, &
       &        GSYMAX, GSYMIN, GSZMAX, GSZMIN, GTHETA, GXL, GXMAX, GXMIN, GXORG, &
       &        GYL, GYMAX, GYMIN, GYORG, GZL, GZMAX, GZMIN, GZVAL
  character(LEN=80) :: KT
  character(LEN=1)  :: KDL
  external R2W2B,W2G2B,R2Y2W,WRGBW

  call SETCHS(0.3,0.0)
  call SETFNT(32)
  call SETLNW(0.035)
  call SETRGB(0.0,0.0,0.0)
  call SETLIN(-1,-1,7)
  KDL=STR(1:1)
  I=2
1 if(STR(I:I) == KDL .or. I == 80) go to 2
  KT(I-1:I-1)=STR(I:I)
  I=I+1
  go to 1

2 call MOVE(GX1,GY2+0.2)
  call TEXT(KT,I-2)

  call GMNMX2(GZ,NXM,1,NXMAX,1,1,NYMAX,1,GZMIN,GZMAX)
  GZVAL=0.5*(abs(GZMIN)+abs(GZMAX))
  if((GZVAL <= 1.D-6) .or. (abs(GZMAX-GZMIN)/GZVAL < 1.D-6)) then
     write(6,*) GZMIN,GZMAX
     flush(6)
     return
  end if
  !      write(6,*) GZMIN,GZMAX
  !      flush(6)

  if(mod(MODE,2) == 0) then
     if(GZMIN >= 0.0) then
        GZMIN=0.0
     else if(GZMAX <= 0.0) then
        GZMAX=0.0
     end if
  end if

  call GMNMX1(GX,1,NXMAX,1,GXMIN,GXMAX)
  call GMNMX1(GY,1,NYMAX,1,GYMIN,GYMAX)

  if(mod(MODE/4,2) == 1) then
     call CHMODE
     write(6,*) '## TXGR : XMIN,XMAX,YMIN,YMAX = ',GXMIN,GXMAX,GYMIN,GYMAX
     read(5,*) GXMIN,GXMAX,GYMIN,GYMAX
     call GRMODE
  end if

  call GQSCAL(GXMIN,GXMAX,GSXMIN,GSXMAX,GSTEPX)
  call GQSCAL(GYMIN,GYMAX,GSYMIN,GSYMAX,GSTEPY)
  call GQSCAL(GZMIN,GZMAX,GSZMIN,GSZMAX,GSTEPZ)
  !      write(6,*) GSZMIN,GSZMAX,GSTEPZ
  !      flush(6)

  if(GXMIN*GXMAX < 0.0) then
     GXORG=0.0
  else
     GXORG=GSXMIN
  end if
  if(GYMIN*GYMAX < 0.0) then
     GYORG=0.0
  else
     GYORG=GSYMIN
  end if
  if(GZMIN*GZMAX < 0.0) then
     GSZMAX=max(abs(GSZMIN),abs(GSZMAX))
     GSZMIN=-GSZMAX
  end if

  GXL=10.0*1.5
  GYL=20.0*1.5
  GZL=10.0*1.5
  call GDEFIN3D(GX1,GX2,GY1,GY2,GXL,GYL,GZL)

  ! viewpoint
  GPHI=-65.0
  GTHETA=45.0
  GRADIUS=35.0
!  GPHI=-60.0
!  GTHETA=65.0
!  GRADIUS=100.0

  GOX=0.5*(GSXMIN+GSXMAX)
  GOY=0.5*(GYMIN+GYMAX)
  GOZ=0.5*(GSZMIN+GSZMAX)
  call GVIEW3D(GPHI,GTHETA,GRADIUS,1.0,1,GOX,GOY,GOZ)
  call GDATA3D2(GZ,GX,GY,NXM,NXMAX,NYMAX,GSZMIN,GSZMAX)
  call SETCHS(0.2,0.0)
  call SETLIN(0,0,7)

  call GSCALE3DX(GSXMIN,GSTEPX,0.3,0)
  call GSCALE3DY(GSYMIN,GSTEPY,0.3,0)
  call GSCALE3DZ(GSZMIN,GSTEPZ,0.3,0)
  call GVALUE3DX(GSXMIN,GSTEPX,1,NGULEN(GSTEPX))
  call GVALUE3DY(GSYMIN,GSTEPY,1,NGULEN(GSTEPY))
  call GVALUE3DZ(GSZMIN,GSTEPZ,2,NGULEN(GSTEPZ))

  !      call GTEXTX(GX1-0.3,0.5*(GY1+GY2),
  !     &              '@TIME (sec)@',
  !     &              2)
  !      call GTEXTX(0.5*(GX1+GX2),GY1-0.3,
  !     &              '@R@',
  !     &              2)
  !      call GTEXTX(GX1,GY2+0.1,
  !     &              KV,
  !     &              0)

  if(GZMIN*GZMAX < 0.0) then
     call CPLOT3D1(7,R2W2B)
  else if(GZMIN < 0.0) then
     call CPLOT3D1(7,W2G2B)
  else
     call CPLOT3D1(7,R2Y2W)
  end if

  !      call CONTQ3D1(GZMIN,0.1*(GZMAX-GZMIN),11,0,0,KA,R2W2B,2)
  !      call PERSE3D(3,1)
  call GAXIS3D(0)
  !      call GDrawBack3D(0.5, 0.5, 0.5)

  call SETLIN(0,0,7)
  
end subroutine TXGR3D
