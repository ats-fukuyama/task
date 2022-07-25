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
!  character(LEN=80),INTENT(IN):: STR, KV
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
  call GSGLENABLELIGHTING
  call TXGR3D(GX1,GX2,GY1,GY2,GX,GTX,GYL,NRMAX+1,NRMAX+1,NGT+1,STR,KV,2+INQ)
  call PAGEE
  
end subroutine TXGRUR

!***********************************************************
!
!   SUBPROGRAM FOR 3D PROFILE
!
!***********************************************************

subroutine TXGR3D(GX1,GX2,GY1,GY2,GX,GY,GZ,NXM,NXMAX,NYMAX,STR,KV,MODE)

  implicit none
  real(4),    intent(IN) :: GX1, GX2, GY1, GY2
  integer(4), intent(IN) :: NXM, NXMAX, NYMAX, MODE
  real(4), dimension(NXMAX),     intent(IN) :: GX
  real(4), dimension(NYMAX),     intent(IN) :: GY
  real(4), dimension(NXM,NYMAX), intent(IN) :: GZ
  integer(4) :: I
  real(4) :: GXMIN, GXMAX, GYMIN, GYMAX, GZMIN, GZMAX, GSXMIN, GSXMAX, GSYMIN, GSYMAX, &
       &     GSZMIN, GSZMAX, GSTEPX, GSTEPY, GSTEPZ, GXL, GYL, GZL, GPHI, GTHETA, &
       &     GRADIUS, GOX, GOY, GOZ
  character(LEN=80) :: STR, KV
  character(LEN=80) :: KT
  character(LEN=1)  :: KDL
  external R2G2B_Gat0

  call SETFNT(32)
  call SETCHS(0.3,0.0)
  call SETLIN(0,0,7)
  KDL=STR(1:1)
  I=2
1 if(STR(I:I) == KDL .or. I == 80) go to 2
  KT(I-1:I-1)=STR(I:I)
  I=I+1
  go to 1

2 call MOVE(GX1,GY2+0.2)
  call TEXT(KT,I-2)

  call GMNMX2(GZ,NXM,1,NXMAX,1,1,NYMAX,1,GZMIN,GZMAX)
  if(abs(GZMAX-GZMIN) < 1.D-6) then
     GZMIN=GZMIN-0.999E-6
     GZMAX=GZMAX+1.000E-6
  end if

  if(mod(MODE,2) == 0) then
     if(GZMIN >= 0.0) then
        GZMIN=0.0
     elseif(GZMAX <= 0.0) then
        GZMAX=0.0
     end if
  end if

  call GMNMX1(GX,1,NXMAX,1,GXMIN,GXMAX)
  call GMNMX1(GY,1,NYMAX,1,GYMIN,GYMAX)
  if(abs(GXMAX-GXMIN) < 1.D-6) then
     GXMIN=GXMIN-0.999E-6
     GXMAX=GXMAX+1.000E-6
  end if

  !      GXMIN=GX(1)
  !      GXMAX=GX(NXMAX)

  call GQSCAL(GXMIN,GXMAX,GSXMIN,GSXMAX,GSTEPX)
  call GQSCAL(GYMIN,GYMAX,GSYMIN,GSYMAX,GSTEPY)
  call GQSCAL(GZMIN,GZMAX,GSZMIN,GSZMAX,GSTEPZ)
  ! Correct the misestimation by GQSCAL
  if(GZMAX > GSZMAX) GSZMAX = GSZMAX + GSTEPZ
  if(GZMIN < GSZMIN) GSZMIN = GSZMIN - GSTEPZ

  ! Origin of the vertical axis is forced to set zero
  !   if either all the values are positive or negative.
  if(mod(MODE,2) == 0) then
     if(GZMIN >= 0.0) then
        GSZMIN=0.0
     elseif(GZMAX <= 0.0) then
        GSZMAX=0.0
     end if
  end if
  GZMIN=GSZMIN
  GZMAX=GSZMAX
  if(mod(MODE/4,2) == 1) then
     call CHMODE
     write(6,*) '## TXGR : XMIN,XMAX,YMIN,YMAX = ',GXMIN,GXMAX,GYMIN,GYMAX
     read(5,*) GXMIN,GXMAX,GYMIN,GYMAX
     call GRMODE
  end if
  call GQSCAL(GXMIN,GXMAX,GSXMIN,GSXMAX,GSTEPX)
  call GQSCAL(GYMIN,GYMAX,GSYMIN,GSYMAX,GSTEPY)
  call GQSCAL(GZMIN,GZMAX,GSZMIN,GSZMAX,GSTEPZ)

  !      if(GXMIN*GXMAX.LE.0.0) then
  !         GXORG=0.0
  !      else
  !         GXORG=GSXMIN
  !      end if
  !      if(GYMIN*GYMAX.LE.0.0) then
  !         GYORG=0.0
  !      else
  !         GYORG=GSYMIN
  !      end if

  call GDEFIN(GX1,GX2,GY1,GY2,GSXMIN,GSXMAX,GSYMIN,GSYMAX)
  GXL=10.0*1.5
  GYL=20.0*1.5
  GZL=10.0*1.5
  call GDEFIN3D(GXL,GYL,GZL,GSXMIN,GSXMAX,GYMIN,GYMAX,GSZMIN,GSZMAX)

  ! viewpoint
  GPHI=-65.0
  GTHETA=45.0
  GRADIUS=35.0

  GOX=0.5*(GSXMIN+GSXMAX)
  GOY=0.5*(GYMIN+GYMAX)
  GOZ=0.5*(GSZMIN+GSZMAX)
  call GVIEW3D(GPHI,GTHETA,GRADIUS,GOX,GOY,GOZ)
  call SETCHS(0.3,0.0)
  call SETLIN(0,0,4)

  !      call GSCALE3DX(GSXMIN,GSTEPX,0.3,0)
  !      call GSCALE3DY(GT(1),GSTEPT,0.3,0)
  !      call GSCALE3DZ(GSYMIN,GSTEPY,0.3,10)
  !      call GVALUE3DX(GSXMIN,GSTEPX,1,1)
  !      call GVALUE3DY(GT(1),GSTEPT,1,1)
  !      call GVALUE3DZ(GSYMIN,GSTEPY,11,-2)
  call GSCALE3DX(GSXMIN,GSTEPX,0.3,0)
  call GSCALE3DY(GSYMIN,GSTEPY,0.3,0)
  call GSCALE3DZ(GSZMIN,GSTEPZ,0.3,10)
  call GVALUE3DX(GSXMIN,GSTEPX,1,1)
  call GVALUE3DY(GSYMIN,GSTEPY,1,3)
  call GVALUE3DZ(GSZMIN,GSTEPZ,11,-2)

  call Set3DTextBaseLine(0.0, 1.0, 0.0, -1.0, 0.0, 0.0)
  call GTEXTX3D(GSXMAX+0.15*(GSXMAX-GSXMIN),0.5*(GY(1)+GY(NYMAX)), &
       &              GSZMIN,'@TIME (sec)@',2)
  call Set3DTextBaseLine(0.0, 1.0, 0.0, -1.0, 0.0, 0.0)
  call GTEXTX3D(0.5*(GSXMIN+GSXMAX),GY(1)+0.1*(GY(1)-GY(NYMAX)), &
       &              GSZMIN,'@R@',2)
  call Set3DTextBaseLine(0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
  call GTEXTX3D(GSXMIN,GY(1)+0.05*(GY(1)-GY(NYMAX)), &
       &              GSZMAX+0.1*(GSZMAX-GSZMIN),KV,2)

  call PERS3D2(GZ,GX,GY,NXM,NXMAX,NYMAX,-27,R2G2B_Gat0)
  call GAxis3D(0)
  call GDrawBack3D(0.5, 0.5, 0.5)

  call SETLIN(0,0,4)
  
end subroutine TXGR3D
