!     $Id: wfgsub.f90,v 1.19 2012/03/05 06:29:02 maruyama Exp $

!     ****** CALCULATE RANGE OF WINDOW ******

SUBROUTINE WFGWIN(NW,NWMAX,PXMIN,PXMAX,PYMIN,PYMAX)

  use wfcomm
  implicit none
  integer,intent(in) :: NW,NWMAX
  integer :: NWW,NWYMAX,MIN,NWX,NWY
  real(rkind),intent(out) :: PXMIN,PXMAX,PYMIN,PYMAX
  real(rkind) :: DXLEN,DYLEN,DRATIO,PXLEN,PYLEN

  DXLEN= RNDMAX-RNDMIN
  DYLEN= ZNDMAX-ZNDMIN
  DRATIO=DYLEN/DXLEN
  
  IF(NWXMAX.EQ.0) THEN
     IF(DRATIO.GT.1.25D0) THEN
        NWW=6
     ELSEIF(DRATIO.GT.0.5D0) THEN
        IF(NWMAX.LE.3) then
           NWW=3
        ELSE IF(NWMAX.LE.4) THEN
           NWW=2
        ELSE IF(NWMAX.LE.6) THEN
           NWW=3
        ELSE IF(NWMAX.LE.8) THEN
           NWW=4
        ELSE
           NWW=5
        ENDIF
        NWW=4
     ELSE
        NWW=2
     ENDIF
  ELSE
     NWW=NWXMAX
  ENDIF
  
  PXMIN=0.0D0
  PXMAX=25.6D0
  PYMIN=0.0D0
  PYMAX=14.0D0
  
  NWYMAX=(NWMAX-1)/NWW+1
  PXLEN=(PXMAX-PXMIN)/MIN(NWW,NWMAX)
  PYLEN=(PYMAX-PYMIN)/NWYMAX
  NWX=MOD(NW-1,NWW)+1
  PXMIN=PXMIN+(NWX-1)*PXLEN
  PXMAX=PXMIN+PXLEN
  NWY=NWYMAX-(NW-1)/NWW
  PYMIN=PYMIN+(NWY-1)*PYLEN
  PYMAX=PYMIN+PYLEN

  RETURN
END SUBROUTINE WFGWIN

!     ****** PLOT POWER ABSORPTION ******

SUBROUTINE PWRPLOT(NS)

  use wfcomm
  USE libgrf
  implicit none

  integer,intent(in) :: NS
  integer    :: IE,NGX,NGY
  real(rkind)    :: X,Y
  integer    :: I,J,K,IN
  real(rkind)    :: PABS
  real(rkind)    :: RW,WGT(3)
  complex(rkind) :: DTENS(NSM,3,3,3)
  complex(rkind) :: CIWE,CTENS(3,3),CER,CEP,CEZ,JP(3)

  ! JP: plasma current
  ! --- initialize ---
  
  RW=2.D0*PI*RF*1.D6
  CIWE=CII*RW*EPS0

  CALL WFGMESH
  DO NGY=1,NGYMAX
     Y=G2Y(NGY)
     DO NGX=1,NGXMAX
        X=G2X(NGX)
        IE=IEGZ(NGX,NGY)
        IF(IE.EQ.0) THEN
           GZ(NGX,NGY)=0.0
        ELSE
           ! --- calculate conductivity tensor ---

           CALL DTENSR(IE,DTENS)
           call WFWGT(IE,X,Y,WGT)

           CTENS=(0.d0,0.d0)

           do J=1,3
              do I=1,3
                 do IN=1,3
                    CTENS(I,J)=CTENS(I,J)-WGT(IN)*CIWE*DTENS(NS,IN,I,J)
                 end do
              end do
           end do

           ! --- calculate power absorption ---
           call FIELDCR(IE,X,Y,CESD,CER)
           call FIELDCP(IE,X,Y,CEND,CEP)
           call FIELDCZ(IE,X,Y,CESD,CEZ)

           JP  =(0.d0,0.d0)
           PABS=0.d0

           do K=1,3
              JP(K)= CTENS(K,1)*CER &
                    +CTENS(K,2)*CEP &
                    +CTENS(K,3)*CEZ
           end do

           PABS=0.5d0*real(JP(1)*conjg(CER)&
                          +JP(2)*conjg(CEP)&
                          +JP(3)*conjg(CEZ))
           IF(ABS(PABS).LE.1.D-20) PABS=0.D0
           GZ(NGX,NGY)=gdclip(PABS)

!        write(6,'(A,3I8,1P3E12.4)') 'NGX,NGY,IE,X,Y,PABS=',NGX,NGY,IE,X,Y,PABS
        ENDIF
     END DO
  ENDDO
     
  RETURN
END SUBROUTINE PWRPLOT

!    ***** PLOT DENSITY PROFILE *****

SUBROUTINE NPLOT

  use wfcomm
  USE libgrf
  implicit none
  
  integer :: NGX,NGY,NE
  real(rkind) :: DX,DY,X,Y
  real(rkind) :: RN(NSM),RTPR(NSM),RTPP(NSM),RZCL(NSM)
  
  NE=0

  DY=(ZNDMAX-ZNDMIN)/(NGYMAX-1)
  DX=(RNDMAX-RNDMIN)/(NGXMAX-1)
  DO NGX=1,NGXMAX
     G2X(NGX)=gdclip(RNDMIN+DX*(NGX-1))
  ENDDO
  DO NGY=1,NGYMAX
     G2Y(NGY)=gdclip(ZNDMIN+DY*(NGY-1))
  ENDDO

  DO NGY=1,NGYMAX
     Y=ZNDMIN+DY*(NGY-1)
     DO NGX=1,NGXMAX
        X=RNDMIN+DX*(NGX-1)
        CALL FEP(X,Y,NE)
        IF(NE.EQ.0) THEN
           GZ(NGX,NGY)=0.0
        ELSE
           CALL WFSDEN(X,Y,RN,RTPR,RTPP,RZCL)
           GZ(NGX,NGY)=gdclip(RN(1)*1.d20)
        end IF
     end DO
  end DO

  return
END SUBROUTINE NPLOT

!     ****** PLOT 2D STRUCTURE OF FIELD AND POWER ******

SUBROUTINE WFGPPC(NW,NWMAX,KWD)
  
  use wfcomm
  USE libgrf
  implicit none
  integer,intent(in) :: NW,NWMAX
  integer :: NGX,NGY,ISTEP
  real :: GXMIN,GXMAX,GYMIN,GYMAX,GPYMIN,GPYMAX
  real :: GPXMIN,GPXMAX,GSXMIN,GSXMAX,GXSCAL,GSYMIN,GYSCAL
  real :: GXORG,GYORG,GZMIN,GZMAX,GQZMIN,GQZMAX,GZSCAL,GZDEL,GZORG,GXPOS,GYPOS
  real(rkind) :: XMIN,XMAX,YMIN,YMAX,PXMIN,PXMAX,PYMIN,PYMAX
  real(rkind) :: DXLEN,DYLEN,DRATIO,PRATIO
  real(rkind) :: PXMID,PYLEN,PXLEN,PYMID,DX,DY
  real(rkind) :: n_para
  character,intent(in):: KWD*(NCHM)

  integer :: NS,NSDO
  real(rkind) :: X,Y,WC(NSMAX),WP(NSMAX),WW
  real(rkind) :: BABS,AL(3),RTPR(NSM),RTPP(NSM),RZCL(NSM),RN(NSM)
  REAL :: GTCO_MAX,GTCR_MAX,GTHR_MAX,GTRC_MAX,GTLC_MAX
  INTEGER :: NBSD,NSD,NN1,NN2
  REAL :: R1,Z1,R2,Z2

  REAL(rkind),DIMENSION(:,:),ALLOCATABLE:: TCO,TCR,THR,TRC,TLC
  real :: GUCLIP
  integer,DIMENSION (:,:,:),ALLOCATABLE :: KA
  real,DIMENSION (:),ALLOCATABLE :: GAX,GAY
  real,DIMENSION (:,:),ALLOCATABLE :: GTCO,GTCR,GTHR,GTRC,GTLC

  IF(NWMAX.LE.0) RETURN

  ALLOCATE(KA(4,NGXMAX,NGYMAX))
  ALLOCATE(GAX(NGXMAX),GAY(NGYMAX))
  ALLOCATE(GTCO(NGXMAX,NGYMAX),GTCR(NGXMAX,NGYMAX),GTHR(NGXMAX,NGYMAX))
  ALLOCATE(GTRC(NGXMAX,NGYMAX),GTLC(NGXMAX,NGYMAX))
  ALLOCATE(TCO(NGXMAX,NGYMAX),TCR(NGXMAX,NGYMAX),THR(NGXMAX,NGYMAX))
  ALLOCATE(TRC(NGXMAX,NGYMAX),TLC(NGXMAX,NGYMAX))
  
  ! --- initialize ---

  CALL WFGWIN(NW,NWMAX,PXMIN,PXMAX,PYMIN,PYMAX)
  
  XMIN = RNDMIN
  XMAX = RNDMAX
  YMIN = ZNDMIN
  YMAX = ZNDMAX  
  DXLEN= XMAX-XMIN
  DYLEN= YMAX-YMIN
  DRATIO=DYLEN/DXLEN
  PRATIO=(PYMAX-PYMIN-1.5D0)/(PXMAX-PXMIN-0.6D0)
  GXMIN=gdclip(XMIN)
  GXMAX=gdclip(XMAX)
  GYMIN=gdclip(YMIN)
  GYMAX=gdclip(YMAX)
  IF(DRATIO.GT.PRATIO) THEN
     GPYMIN=gdclip(PYMIN)+1.5
     GPYMAX=gdclip(PYMAX)
     PXMID=0.5D0*(PXMIN+PXMAX)
     PYLEN=PYMAX-PYMIN-1.5D0
     PXLEN=PYLEN/DRATIO
     GPXMIN=gdclip(PXMID-0.5D0*PXLEN)
     GPXMAX=gdclip(PXMID+0.5D0*PXLEN)
  ELSE
     GPXMIN=gdclip(PXMIN)+0.3
     GPXMAX=gdclip(PXMAX)-0.3
     PYMID=0.5D0*(PYMIN+PYMAX+1.5D0)
     PXLEN=PXMAX-PXMIN-0.6D0
     PYLEN=PXLEN*DRATIO
     GPYMIN=gdclip(PYMID-0.5D0*PYLEN)
     GPYMAX=gdclip(PYMID+0.5D0*PYLEN)
  ENDIF

  ! --- set GAX,GAY ---
  
  DX=(XMAX-XMIN)/(NGXMAX-1)
  DO NGX=1,NGXMAX
     GAX(NGX)=gdclip(XMIN+(NGX-1)*DX)
  ENDDO
  DY=(YMAX-YMIN)/(NGYMAX-1)
  DO NGY=1,NGYMAX
     GAY(NGY)=gdclip(YMIN+(NGY-1)*DY)
  ENDDO

  ! --- scaling X,Y ---
  
  CALL GQSCAL(GXMIN,GXMAX,GSXMIN,GSXMAX,GXSCAL)
  CALL GQSCAL(GYMIN,GYMAX,GSYMIN,GSXMAX,GYSCAL)
  IF(GXMIN*GXMAX.LT.0.0) THEN
     GXORG=0.0
  ELSE
     GXORG=GSXMIN+GXSCAL
  ENDIF
  IF(GYMIN*GYMAX.LT.0.0) THEN
     GYORG=0.0
  ELSE
     GYORG=GSYMIN+GYSCAL
  ENDIF

  ! --- draw cutoff and resonance layer ---
  if (KWD(1:1).eq.'P') then

     WW=2*PI*RF*1.0d6
     if     (KWD(2:2).eq.'1') then
        NS=1
     elseif (KWD(2:2).eq.'2') then
        NS=2
     elseif (KWD(2:2).eq.'3') then
        NS=3
     end if

!!! the following should be corrected by including poloidal magnetic field
     SELECT CASE(MODELG)
     CASE(0,12)
        n_para=0.D0
     CASE(1)
        n_para=NPH*VC/(RA*WW)
     CASE(2:6)
        n_para=NPH*VC/(RR*WW)
     END SELECT

     GTCO_MAX=0.D0
     GTCR_MAX=0.D0
     GTHR_MAX=0.D0
     GTRC_MAX=0.D0
     GTLC_MAX=0.D0

     DO NGY=1,NGYMAX
        Y=ZNDMIN+DY*(NGY-1)
        DO NGX=1,NGXMAX
           X=RNDMIN+DX*(NGX-1)

           CALL WFSMAG(X,Y,BABS,AL)
           CALL WFSDEN(X,Y,RN,RTPR,RTPP,RZCL)

           DO NSDO=1,NSMAX
              WP(NSDO)=PZ(NSDO)*PZ(NSDO)*AEE*AEE*RN(NSDO)*1.D20 &
                            /(PA(NSDO)*AMP*EPS0*WW*WW)
              WC(NSDO)=PZ(NSDO)*AEE*BABS/(PA(NSDO)*AMP*WW)
           end DO

           TCR(NGX,NGY)=WC(1)**2
           TCO(NGX,NGY)=1.D0
           THR(NGX,NGY)=1.D0
           TRC(NGX,NGY)=1.D0
           TLC(NGX,NGY)=1.D0
           do NSDO=1,NSMAX
              TCO(NGX,NGY)=TCO(NGX,NGY)-WP(NSDO)
              THR(NGX,NGY)=THR(NGX,NGY)-WP(NSDO)/(1.D0-WC(NSDO)**2)
              TRC(NGX,NGY)=TRC(NGX,NGY)-WP(NSDO)/(1.D0+WC(NSDO))
              TLC(NGX,NGY)=TLC(NGX,NGY)-WP(NSDO)/(1.D0-WC(NSDO))
           end do

           GTCR(NGX,NGY)=GUCLIP(1.D0-TCR(NGX,NGY))
           GTCO(NGX,NGY)=GUCLIP(TCO(NGX,NGY))
           GTHR(NGX,NGY)=GUCLIP(THR(NGX,NGY))
           GTRC(NGX,NGY)=GUCLIP(n_para**2-TRC(NGX,NGY))
           GTLC(NGX,NGY)=GUCLIP(n_para**2-TLC(NGX,NGY))
           GTCR_MAX=MAX(GTCR_MAX,ABS(GTCR(NGX,NGY)))
           GTCO_MAX=MAX(GTCO_MAX,ABS(GTCO(NGX,NGY)))
           GTHR_MAX=MAX(GTHR_MAX,ABS(GTHR(NGX,NGY)))
           GTRC_MAX=MAX(GTRC_MAX,ABS(GTRC(NGX,NGY)))
           GTLC_MAX=MAX(GTLC_MAX,ABS(GTLC(NGX,NGY)))

        END DO
     END DO

     CALL GDEFIN(GPXMIN,GPXMAX,GPYMIN,GPYMAX,GXMIN,GXMAX,GYMIN,GYMAX)
     CALL SETCHS(0.3,0.0)
     CALL SETLNW(0.03)
     
     !  ****** DRAW CYCROTRON REASONANCE SURFACE (yellow, dot-dashed) ******
     
     CALL SETRGB(1.0,0.5,0.0)
     CALL CONTP2(GTCR,GAX,GAY,NGXMAX,NGXMAX,NGYMAX,0.0,0.5*GTCR_MAX,1,0,4,KA)
     
     !   ****** DRAW PLASMA CUT OFF SURFACE (black, two-dot-dashed) ******
     
     CALL SETRGB(0.0,0.0,0.0)
     CALL CONTP2(GTCO,GAX,GAY,NGXMAX,NGXMAX,NGYMAX,0.0,0.5*GTCO_MAX,1,0,6,KA)
     
     !     ****** DRAW RIGHT CUT OFF SURFACE (blue, two-dot-dashed) ******

     CALL SETRGB(0.0,0.0,1.0)
     CALL CONTP2(GTRC,GAX,GAY,NGXMAX,NGXMAX,NGYMAX,0.0,0.5*GTRC_MAX,1,0,6,KA)
     
     !     ****** DRAW LEFT CUT OFF SURFACE (green, two-dots-dashed) ******
     
     CALL SETRGB(0.0,1.0,0.0)
     CALL CONTP2(GTLC,GAX,GAY,NGXMAX,NGXMAX,NGYMAX,0.0,0.5*GTLC_MAX,1,0,6,KA)
     
     !     ****** DRAW HYBRID RESONANCE SURFACE (purple, long-dashed) ******
     
     CALL SETRGB(1.0,0.0,1.0)
     CALL CONTP2(GTHR,GAX,GAY,NGXMAX,NGXMAX,NGYMAX,0.0,0.5*GTHR_MAX,1,0,3,KA)
     
  end if

  ! --- smoozing Z ---
  do NGY=1,NGYMAX
     do NGX=1,NGXMAX
        GZ_temp(NGX,NGY)=GZ(NGX,NGY)
     end do
  end do

  do NGY=1,NGYMAX
     do NGX=1,NGXMAX
        GZ(NGX,NGY)=0.0
        if(NGX.ne.1.and.NGY.ne.1) &
             GZ(NGX,NGY)=GZ(NGX,NGY)+GZ_temp(NGX-1,NGY-1)*0.0625
        if(NGY.ne.1) &
             GZ(NGX,NGY)=GZ(NGX,NGY)+GZ_temp(NGX  ,NGY-1)*0.125
        if(NGX.ne.NGXMAX.and.NGY.ne.1) &
             GZ(NGX,NGY)=GZ(NGX,NGY)+GZ_temp(NGX+1,NGY-1)*0.0625
        if(NGX.ne.1) &
             GZ(NGX,NGY)=GZ(NGX,NGY)+GZ_temp(NGX-1,NGY  )*0.125
        GZ(NGX,NGY)=GZ(NGX,NGY)+GZ_temp(NGX  ,NGY  )*0.25
        if(NGX.ne.NGXMAX) &
             GZ(NGX,NGY)=GZ(NGX,NGY)+GZ_temp(NGX+1,NGY  )*0.125
        if(NGX.ne.1.and.NGY.ne.NGYMAX) &
             GZ(NGX,NGY)=GZ(NGX,NGY)+GZ_temp(NGX-1,NGY+1)*0.0625
        if(NGY.ne.NGYMAX) &
             GZ(NGX,NGY)=GZ(NGX,NGY)+GZ_temp(NGX  ,NGY+1)*0.125
        if(NGX.ne.NGXMAX.and.NGY.ne.NGYMAX) &
             GZ(NGX,NGY)=GZ(NGX,NGY)+GZ_temp(NGX+1,NGY+1)*0.0625
     end do
  end do

  ! --- scaling Z & PLOT---

  CALL GMNMX2(GZ,NGXMAX,1,NGXMAX,1,1,NGYMAX,1,GZMIN,GZMAX)
  CALL GQSCAL(GZMIN,GZMAX,GQZMIN,GQZMAX,GZSCAL)
  GZDEL=REAL(GFACTOR)*GZSCAL
  IF(GZDEL.EQ.0.0) GOTO 1000
  ISTEP=INT((GZMAX-GZMIN)/GZDEL)
  
  CALL GDEFIN(GPXMIN,GPXMAX,GPYMIN,GPYMAX,GXMIN,GXMAX,GYMIN,GYMAX)
  CALL SETCHS(0.3,0.0)
  CALL SETLNW(0.03)
  
  CALL SETLIN(0,0,7)
  CALL GFRAME
  CALL GSCALE(GXORG,GXSCAL,0.0,0.0,0.1,9)
  CALL GSCALE(0.0,0.0,GYORG,GYSCAL,0.1,9)

  CALL SETLIN(0,0,7)
  IF(GZMIN*GZMAX.GT.0.0) THEN
!     GZORG=GQZMIN-GZDEL
     GZORG=0.5*GZDEL
  ELSE
     GZORG=0.5*GZDEL
  ENDIF

  CALL SETLIN(0,0,6)
  CALL CONTP2(GZ,GAX,GAY,NGXMAX,NGXMAX,NGYMAX, GZORG, GZDEL,ISTEP,0,0,KA)
  CALL SETLIN(0,0,5)
  CALL CONTP2(GZ,GAX,GAY,NGXMAX,NGXMAX,NGYMAX,-GZORG,-GZDEL,ISTEP,0,0,KA)

! ----- PRINT MAX,MIN,STP ----- 

1000 CONTINUE
  CALL SETLIN(0,2,7)
  DO NBSD=1,NBSID
     NSD=NSDBS(NBSD)
     NN1=NDSID(1,NSD)
     NN2=NDSID(2,NSD)
     R1=guclip(RNODE(NN1))
     Z1=guclip(ZNODE(NN1))
     R2=guclip(RNODE(NN2))
     Z2=guclip(ZNODE(NN2))
     CALL MOVE2D(R1,Z1)
     CALL DRAW2D(R2,Z2)
  END DO

  CALL SETLIN(0,0,7)
  CALL SETCHS(0.2,0.0)
  GXPOS=0.5*(GPXMIN+GPXMAX-17*0.2)
  GYPOS=GPYMIN-0.4
  CALL MOVE(GXPOS,GYPOS)
  CALL TEXT(KWD(1:3),3)
  GYPOS=GYPOS-0.3
  CALL MOVE(GXPOS,GYPOS)
  CALL TEXT('MAX=',4)
  CALL NUMBR(GZMAX,'(1PE9.2)',9)
  GYPOS=GYPOS-0.3
  CALL MOVE(GXPOS,GYPOS)
  CALL TEXT('MIN=',4)
  CALL NUMBR(GZMIN,'(1PE9.2)',9)
  GYPOS=GYPOS-0.3
  CALL MOVE(GXPOS,GYPOS)
  CALL TEXT('STP=',4)
  CALL NUMBR(GZDEL,'(1PE9.2)',9)

  DEALLOCATE(KA)
  DEALLOCATE(GAX,GAY)
  DEALLOCATE(GTCO,GTCR,GTHR,GTRC,GTLC)
  DEALLOCATE(TCO,TCR,THR,TRC,TLC)
  RETURN
END SUBROUTINE WFGPPC

!     ****** PLOT 2D STRUCTURE OF FIELD AND POWER ******

SUBROUTINE WFGPFC(NW,NWMAX,KWD)

  use wfcomm
  USE libgrf
  implicit none
  integer,parameter :: NSTEPM=101
  integer,parameter :: NRGBA=5
  integer,parameter :: NRGBB=7
  integer,intent(in):: NWMAX,NW
  integer :: NGX,NGY,ISTEP,I
  real(rkind) :: PXMIN,PXMAX,PYMIN,PYMAX,XMIN,XMAX,YMIN,YMAX
  real(rkind) :: DXLEN,DYLEN,DRATIO,PRATIO
  real(rkind) :: PXMID,PYLEN,PXLEN,PYMID,DX,DY
  real :: GXMIN,GXMAX,GYMIN,GYMAX,GPYMIN,GPYMAX,GPXMIN,GPXMAX
  real :: GSXMAX,GSXMIN,GXSCAL,GSYMIN,GYSCAL,GXORG,GYORG,GZMIN,GZMAX,GZA
  real :: GDZ,GFACT,GXPOS,GYPOS
  real :: GAX(NGXMAX),GAY(NGYMAX),GZL(NSTEPM),GRGBL(3,0:NSTEPM)
  real :: GRGBA(3,NRGBA),GLA(NRGBA),GRGBB(3,NRGBB),GLB(NRGBB)
  character,intent(in) :: KWD*(NCHM)

  DATA GRGBA/0.0,0.0,1.0,&
       &     0.0,1.0,1.0,&
       &     1.0,1.0,1.0,&
       &     1.0,1.0,0.0,&
       &     1.0,0.0,0.0/
  DATA GLA/0.0,0.40,0.5,0.60,1.0/
  DATA GRGBB/1.0,1.0,1.0,&
       &     0.0,1.0,1.0,&
       &     0.0,0.0,1.0,&
       &     0.0,1.0,0.0,&
       &     1.0,0.0,0.0,&
       &     1.0,1.0,0.0,&
       &     1.0,1.0,1.0/
  DATA GLB/0.0,0.15,0.3,0.5,0.7,0.85,1.0/
  
  IF(NWMAX.LE.0) RETURN
  
  CALL WFGWIN(NW,NWMAX,PXMIN,PXMAX,PYMIN,PYMAX)
  
  XMIN = RNDMIN
  XMAX = RNDMAX
  YMIN = ZNDMIN
  YMAX = ZNDMAX
  
  DXLEN= XMAX-XMIN
  DYLEN= YMAX-YMIN
  DRATIO=DYLEN/DXLEN
  PRATIO=(PYMAX-PYMIN-1.5D0)/(PXMAX-PXMIN-0.6D0)
  GXMIN=gdclip(XMIN)
  GXMAX=gdclip(XMAX)
  GYMIN=gdclip(YMIN)
  GYMAX=gdclip(YMAX)
  IF(DRATIO.GT.PRATIO) THEN
     GPYMIN=gdclip(PYMIN)+1.5
     GPYMAX=gdclip(PYMAX)
     PXMID=0.5D0*(PXMIN+PXMAX)
     PYLEN=PYMAX-PYMIN-1.5D0
     PXLEN=PYLEN/DRATIO
     GPXMIN=gdclip(PXMID-0.5D0*PXLEN)
     GPXMAX=gdclip(PXMID+0.5D0*PXLEN)
  ELSE
     GPXMIN=gdclip(PXMIN)+0.3
     GPXMAX=gdclip(PXMAX)-0.3
     PYMID=0.5D0*(PYMIN+PYMAX+1.5D0)
     PXLEN=PXMAX-PXMIN-0.6D0
     PYLEN=PXLEN*DRATIO
     GPYMIN=gdclip(PYMID-0.5D0*PYLEN)
     GPYMAX=gdclip(PYMID+0.5D0*PYLEN)
  ENDIF
  
  DX=(XMAX-XMIN)/(NGXMAX-1)
  DO NGX=1,NGXMAX
     GAX(NGX)=gdclip(XMIN+(NGX-1)*DX)
  ENDDO
  DY=(YMAX-YMIN)/(NGYMAX-1)
  DO NGY=1,NGYMAX
     GAY(NGY)=gdclip(YMIN+(NGY-1)*DY)
  ENDDO
  
  CALL GQSCAL(GXMIN,GXMAX,GSXMIN,GSXMAX,GXSCAL)
  CALL GQSCAL(GYMIN,GYMAX,GSYMIN,GSXMAX,GYSCAL)
  IF(GXMIN*GXMAX.LT.0.0) THEN
     GXORG=0.0
  ELSE
     GXORG=GSXMIN+GXSCAL
  ENDIF
  IF(GYMIN*GYMAX.LT.0.0) THEN
     GYORG=0.0
  ELSE
     GYORG=GSYMIN+GYSCAL
  ENDIF
  
  ISTEP=50
  
  CALL GMNMX2(GZ,NGXMAX,1,NGXMAX,1,1,NGYMAX,1,GZMIN,GZMAX)
  IF(GZMIN*GZMAX.LT.0.D0) THEN
     GZA=MAX(ABS(GZMAX),ABS(GZMIN))
     GDZ=2*GZA/ISTEP
     DO I=1,ISTEP
        GZL(I)=GDZ*(I-0.5)-GZA
     ENDDO
     
     DO I=0,ISTEP
        GFACT=REAL(I)/REAL(ISTEP)
        CALL GUSRGB(GFACT,GRGBL(1,I),NRGBA,GLA,GRGBA)
     ENDDO
  ELSE
     GZA=GZMIN
     GDZ=(GZMAX-GZMIN)/ISTEP
     DO I=1,ISTEP
        GZL(I)=GDZ*I+GZA
     ENDDO
     DO I=0,ISTEP
        GFACT=REAL(I)/REAL(ISTEP)
        CALL GUSRGB(GFACT,GRGBL(1,I),NRGBB,GLB,GRGBB)
     ENDDO
  ENDIF
  
  CALL GDEFIN(GPXMIN,GPXMAX,GPYMIN,GPYMAX,GXMIN,GXMAX,GYMIN,GYMAX)
  CALL SETCHS(0.3,0.0)
  CALL SETLNW(0.03)
  
  CALL CONTF2(GZ,GAX,GAY,NGXMAX,NGXMAX,NGYMAX,GZL,GRGBL,ISTEP,0)
  CALL SETLIN(0,0,7)
  CALL GFRAME
  CALL GSCALE(GXORG,GXSCAL,0.0,0.0,0.1,9)
  CALL GSCALE(0.0,0.0,GYORG,GYSCAL,0.1,9)
  
  CALL RGBBAR(GPXMAX+0.2,GPXMAX+0.5,GPYMIN,GPYMAX,GRGBL,ISTEP+1,1)
  
!      IF(KID.EQ.'Z') THEN
!         CALL SETVEW(GPXMIN,GPXMAX,GPYMIN,GPYMAX,&
!     &               GXMIN,GXMAX,GYMIN,GYMAX)

!         CALL SETLIN(0,0,7)
!         CALL WFGBDY

!         CALL SETLIN(0,0,4)
!         CALL WFGPLA

!         CALL SETLIN(0,0,4)
!         CALL WFGANT

!         CALL OFFVEW
!      ENDIF

  CALL SETLIN(0,0,7)
  CALL SETCHS(0.2,0.0)
  GXPOS=0.5*(GPXMIN+GPXMAX-17*0.2)
  GYPOS=GPYMIN-0.4
  CALL MOVE(GXPOS,GYPOS)
  CALL TEXT(KWD(1:3),3)
  IF(KWD(4:4).EQ.'X') THEN
     CALL TEXT('(YZ) ',5)
  ELSE IF(KWD(4:4).EQ.'Y') THEN
     CALL TEXT('(XZ) ',5)
  ELSE IF(KWD(4:4).EQ.'Z') THEN
     CALL TEXT('(XY) ',5)
  ENDIF
  CALL TEXT(KWD(4:4),1)
  CALL TEXT('=',1)
  CALL TEXT(KWD(5:NCHM),NCHM-4)
  GYPOS=GYPOS-0.3
  CALL MOVE(GXPOS,GYPOS)
  CALL TEXT('MAX=',4)
  CALL NUMBR(GZMAX,'(1PE9.2)',9)
  GYPOS=GYPOS-0.3
  CALL MOVE(GXPOS,GYPOS)
  CALL TEXT('MIN=',4)
  CALL NUMBR(GZMIN,'(1PE9.2)',9)
  GYPOS=GYPOS-0.3
  CALL MOVE(GXPOS,GYPOS)
  CALL TEXT('TOP=',4)
  CALL NUMBR(GZA,'(1PE9.2)',9)
  RETURN
END SUBROUTINE WFGPFC

!     ****** PLOT 2D STRUCTURE OF FIELD AND POWER ******

SUBROUTINE WFGPBC(NW,NWMAX,KWD)
  
  use wfcomm
  USE libgrf
  implicit none
  integer,intent(in) :: NWMAX,NW
  integer,DIMENSION(:,:,:),ALLOCATABLE :: KA
  integer :: NGX,NGY
  real :: GAX(NGXMAX),GAY(NGYMAX),GXMIN,GXMAX,GYMIN,GYMAX
  real :: GPYMIN,GPYMAX,GPXMIN,GPXMAX
  real(rkind) :: PXMIN,PXMAX,PYMIN,PYMAX,XMIN,XMAX,YMIN,YMAX
  real(rkind) :: DXLEN,DYLEN,DRATIO,PRATIO
  real(rkind) :: PXMID,PYLEN,PXLEN,PYMID,DX,DY
  real :: GSXMIN,GSXMAX,GXSCAL,GSYMIN,GYSCAL,GXORG,GYORG,GZMIN,GZMAX,GSZMIN
  real :: GSZMAX,GZSCAL,GZORG,GXL,GYL,GZL,GPHI,GTHETA,GRADIUS,GXPOS,GYPOS,GZA
  EXTERNAL R2W2B
  character,intent(in) :: KWD*(NCHM)

  IF(NWMAX.LE.0) RETURN
  
  ALLOCATE(KA(8,NGXMAX,NGYMAX))
  CALL WFGWIN(NW,NWMAX,PXMIN,PXMAX,PYMIN,PYMAX)
  
  XMIN = RNDMIN
  XMAX = RNDMAX
  YMIN = ZNDMIN
  YMAX = ZNDMAX
  
  DXLEN= XMAX-XMIN
  DYLEN= YMAX-YMIN
  DRATIO=DYLEN/DXLEN
  PRATIO=(PYMAX-PYMIN-1.5D0)/(PXMAX-PXMIN-0.6D0)
  GXMIN=gdclip(XMIN)
  GXMAX=gdclip(XMAX)
  GYMIN=gdclip(YMIN)
  GYMAX=gdclip(YMAX)
  IF(DRATIO.GT.PRATIO) THEN
     GPYMIN=gdclip(PYMIN)+1.5
     GPYMAX=gdclip(PYMAX)
     PXMID=0.5D0*(PXMIN+PXMAX)
     PYLEN=PYMAX-PYMIN-1.5D0
     PXLEN=PYLEN/DRATIO
     GPXMIN=gdclip(PXMID-0.5D0*PXLEN)
     GPXMAX=gdclip(PXMID+0.5D0*PXLEN)
  ELSE
     GPXMIN=gdclip(PXMIN)+0.3
     GPXMAX=gdclip(PXMAX)-0.3
     PYMID=0.5D0*(PYMIN+PYMAX+1.5D0)
     PXLEN=PXMAX-PXMIN-0.6D0
     PYLEN=PXLEN*DRATIO
     GPYMIN=gdclip(PYMID-0.5D0*PYLEN)
     GPYMAX=gdclip(PYMID+0.5D0*PYLEN)
  ENDIF
  
  DX=(XMAX-XMIN)/(NGXMAX-1)
  DO NGX=1,NGXMAX
     GAX(NGX)=gdclip(XMIN+(NGX-1)*DX)
  ENDDO
  DY=(YMAX-YMIN)/(NGYMAX-1)
  DO NGY=1,NGYMAX
     GAY(NGY)=gdclip(YMIN+(NGY-1)*DY)
  ENDDO
  
  CALL GQSCAL(GXMIN,GXMAX,GSXMIN,GSXMAX,GXSCAL)
  CALL GQSCAL(GYMIN,GYMAX,GSYMIN,GSXMAX,GYSCAL)
  IF(GXMIN*GXMAX.LT.0.0) THEN
     GXORG=0.0
  ELSE
     GXORG=(GSXMIN/(2*GXSCAL)+1)*2*GXSCAL
  ENDIF
  IF(GYMIN*GYMAX.LT.0.0) THEN
     GYORG=0.0
  ELSE
     GYORG=(GSYMIN/(2*GYSCAL)+1)*2*GYSCAL
  ENDIF
  
  CALL GMNMX2(GZ,NGXMAX,1,NGXMAX,1,1,NGYMAX,1,GZMIN,GZMAX)
  GZMAX=MAX(ABS(GZMIN),ABS(GZMAX))
  GZMIN=-GZMAX
  CALL GQSCAL(GZMIN,GZMAX,GSZMIN,GSZMAX,GZSCAL)
  GZORG=0.0
  
  CALL SETLIN(0,0,7)
  GXL=(GPXMAX-GPXMIN)
  GYL=GXL
  GZL=0.5*(GPYMAX-GPYMIN)
  IF(NGRAPH.EQ.3) THEN
     GPHI=-80.0
  ELSEIF(NGRAPH.EQ.4) THEN
     GPHI=190.0
  ELSEIF(NGRAPH.EQ.5) THEN
     GPHI=100.0
  ELSEIF(NGRAPH.EQ.6) THEN
     GPHI=10.0
  ENDIF
  GTHETA=40.0
  GRADIUS=100.0
  
  CALL GDEFIN3D(GPXMIN,GPXMAX,GPYMIN,GPYMAX,GXL,GYL,GZL)
  CALL GVIEW3D(GPHI,GTHETA,GRADIUS,0.9,1,0.5*(GXMIN+GXMAX),0.5*(GYMIN+GYMAX),0.0)
  CALL GDATA3D1(GZ,NGXMAX,NGXMAX,NGYMAX,GXMIN,GXMAX,GYMIN,GYMAX,GZMIN,GZMAX)
  
  CALL CPLOT3D1(1,R2W2B)
  CALL CONTQ3D1(GZMIN,0.1*(GZMAX-GZMIN),11,0,0,KA,R2W2B,0)
  
  CALL GAXIS3D(0)
  CALL GSCALE3DX(GXORG,GXSCAL,0.2,2)
  CALL GSCALE3DY(GYORG,GYSCAL,0.2,2)
  CALL GSCALE3DZ(GZORG,GZSCAL,0.2,2)
  CALL SETCHS(0.2,0.0)
  CALL GVALUE3DX(GXORG,2.*GXSCAL,-6,2)
  CALL GVALUE3DY(GYORG,2.*GYSCAL,-6,2)
  CALL GVALUE3DZ(GZORG,2.*GZSCAL,-2,0)
  
  CALL SETLIN(0,0,7)
  CALL SETCHS(0.2,0.0)
  GXPOS=0.5*(GPXMIN+GPXMAX-17*0.2)
  GYPOS=GPYMIN-0.4
  CALL MOVE(GXPOS,GYPOS)
  CALL TEXT(KWD(1:3),3)
!  IF(KWD(4:4).EQ.'X') THEN
!     CALL TEXT('(YZ) ',5)
!  ELSE IF(KWD(4:4).EQ.'Y') THEN
!     CALL TEXT('(XZ) ',5)
!  ELSE IF(KWD(4:4).EQ.'Z') THEN
  CALL TEXT('(XY) ',5)
!  ENDIF
  CALL TEXT(KWD(4:4),1)
  CALL TEXT('=',1)
  CALL TEXT(KWD(5:NCHM),NCHM-4)
  GYPOS=GYPOS-0.3
  CALL MOVE(GXPOS,GYPOS)
  CALL TEXT('MAX=',4)
  CALL NUMBR(GZMAX,'(1PE9.2)',9)
  GYPOS=GYPOS-0.3
  CALL MOVE(GXPOS,GYPOS)
  CALL TEXT('MIN=',4)
  CALL NUMBR(GZMIN,'(1PE9.2)',9)
  GYPOS=GYPOS-0.3
  CALL MOVE(GXPOS,GYPOS)
  CALL TEXT('TOP=',4)
  CALL NUMBR(GZA,'(1PE9.2)',9)

  DEALLOCATE(KA)

  RETURN
END SUBROUTINE WFGPBC

!     ****** DRAW 1D PROFILE ******

SUBROUTINE WFGPFR(NW,NWMAX,KWD)
  
  use wfcomm
  USE libgrf
  implicit none
  integer,intent(in) :: NWMAX,NW
  integer :: IPAT(3),NGMAX,NG
  real(rkind) :: PXMIN,PYMIN,PXMAX,PYMAX,PRATIO,PXLEN,PYLEN,PXMID
  real :: GPXMIN
  real :: GXMIN1,GXMAX1,GYMIN1,GYMAX1,GXMIN,GXMAX,GXSTEP,GYMIN,GYMAX,GYSTEP
  real :: GXORG,GYORG,GPXMAX,GPYMIN,GPYMAX,GXPOS,GYPOS
  character,intent(in) :: KWD*(NCHM)
  DATA IPAT/0,2,4/
  
  IF(NWMAX.LE.0) RETURN
  
  CALL WFGWIN(NW,NWMAX,PXMIN,PXMAX,PYMIN,PYMAX)

  !     ------ DO NOT ALLOW LANDSCAPE GRAPH ------
  
  PRATIO=(PYMAX-PYMIN-1.5D0)/(PXMAX-PXMIN-2.5D0)
  IF(PRATIO.LT.1.D0) THEN
     PYLEN=PYMAX-PYMIN-1.5D0
     PXLEN=PYLEN+2.5D0
     PXMID=0.5D0*(PXMIN+PXMAX)
     PXMIN=PXMID-0.5D0*PXLEN
     PXMAX=PXMID+0.5D0*PXLEN
  ENDIF

  GPXMIN=gdclip(PXMIN)+2.2
  GPXMAX=gdclip(PXMAX)-0.3
  GPYMIN=gdclip(PYMIN)+1.5
  GPYMAX=gdclip(PYMAX)
  
  IF(KWD(1:1).EQ.'E')THEN!.OR.&
!  &  KWD(1:1).EQ.'D'.OR.&
!  &  KWD(1:1).EQ.'B'.OR.&
!  &  KWD(1:1).EQ.'A'    ) THEN
     NGMAX=3
  ELSE
     NGMAX=1
  ENDIF
  
  CALL GMNMX1(GX,1,NGVMAX,1,GXMIN1,GXMAX1)
  CALL GMNMX2(GV,NGVMAX,1,NGVMAX,1,1,NGMAX,1,GYMIN1,GYMAX1)
  CALL GQSCAL(GXMIN1,GXMAX1,GXMIN,GXMAX,GXSTEP)
  CALL GQSCAL(GYMIN1,GYMAX1,GYMIN,GYMAX,GYSTEP)
  IF(GXMIN*GXMAX.LT.0.0) THEN
     GXORG=0.0
  ELSE
     GXORG=GXMIN
  ENDIF
  IF(GYMIN*GYMAX.LT.0.0) THEN
     GYORG=0.0
  ELSE
     GYORG=GYMIN
  ENDIF
  
  CALL GDEFIN(GPXMIN,GPXMAX,GPYMIN,GPYMAX,&
       &      GXMIN ,GXMAX ,GYMIN ,GYMAX  )
  
  CALL SETLIN(0,0,4)
  CALL SETCHS(0.25,0.0)
  CALL GFRAME
  CALL GSCALE(GXORG,GXSTEP,0.0,0.0,0.2,9)
  CALL GSCALE(0.0,0.0,GYORG,GYSTEP,0.2,9)
  IF(GYMIN*GYMAX.LT.0.0) THEN
     CALL GSCALE(0.0,0.0,GYORG,100*GYSTEP,0.0,0)
  ENDIF
  CALL GVALUE(GXORG,4*GXSTEP,0.0,0.0,NGSLEN(4*GXSTEP))
  CALL GVALUE(0.0,0.0,GYORG,2*GYSTEP,NGSLEN(2*GYSTEP))
  
  DO NG=1,NGMAX
     if(NG.eq.1) CALL SETLIN(0,0,7) !Real part
     if(NG.eq.2) CALL SETLIN(0,0,5) !Imag part
     if(NG.eq.3) CALL SETLIN(0,0,6) !Amplitude
     CALL GPLOTP(GX,GV(1,NG),1,NGVMAX,1,0,0,IPAT(NG))
  ENDDO
  
  CALL SETLIN(0,0,7)
  CALL SETCHS(0.2,0.0)
  GXPOS=0.5*(GPXMIN+GPXMAX-20*0.2)
  GYPOS=GPYMIN-1.0
  CALL MOVE(GXPOS,GYPOS)
  CALL TEXT(KWD(1:2),2)
  IF(KWD(3:3).EQ.'X') THEN
     CALL TEXT('(X): Y=',8)
  ELSE IF(KWD(3:3).EQ.'Y') THEN
     CALL TEXT('(Y): X=',8)
  ENDIF
  CALL TEXT(KWD(4:NCHM),NCHM-3)
  RETURN
END SUBROUTINE WFGPFR

!     ****** DRAW PARAMETER ON GRAPHIC SCREEN ******

SUBROUTINE wf_gdraw_parm
  
  use wfcomm
  implicit none
  integer :: NA,NB,L,NS!,NK,NM
  real(rkind) :: REST(NAM),REAT(NAM),WW,RNZ
  real :: GXMIN,GYMAX,GRCHH,GDX,GDY,GXL,GYL
  real(rkind) :: SRFR(NMDM,NBM),SRFI(NMDM,NBM),SRFL(NMDM,NBM)
  
  DO NA=1,NAMAX
     REST(NA)=DBLE(CIMP(NA))
     REAT(NA)=AIMAG(CIMP(NA))
  ENDDO
  DO NB=1,NBMAX
     DO L=1,NMBDY(NB)
        SRFR(L,NB)=DBLE(CRFL(L,NB))
        SRFI(L,NB)=AIMAG(CRFL(L,NB))
        SRFL(L,NB)=ABS(CRFL(L,NB))**2
     ENDDO
  ENDDO
  
  GXMIN=0.0
  GYMAX=18.2
  GRCHH=0.30
  CALL SETCHS(GRCHH,0.0)
  CALL SETLNW(0.03)
  GDX=15.*GRCHH
  GDY=-1.5*GRCHH
  GXL=GXMIN
  GYL=GYMAX
  
  GXL=GXL+GRCHH
  GYL=GYL+GDY
  CALL MOVE(GXL,GYL)
  CALL TEXT('RF  =',5)
  IF(RF.GE.1.D5) THEN
     CALL NUMBD(RF  ,'(F7.0)',7)
  ELSEIF(RF.GE.1.D4) THEN
     CALL NUMBD(RF  ,'(F7.1)',7)
  ELSEIF(RF.GE.1.D3) THEN
     CALL NUMBD(RF  ,'(F7.2)',7)
  ELSE
     CALL NUMBD(RF  ,'(F7.3)',7)
  ENDIF
  
  GXL=GXL+GDX
  CALL MOVE(GXL,GYL)
  CALL TEXT('BB  =',5)
  CALL NUMBD(BB  ,'(F7.3)',7)
  
  GXL=GXL+GDX
  CALL MOVE(GXL,GYL)
  SELECT CASE(MODELG)
  CASE(0)
     CALL TEXT('RKZ=',4)
     CALL NUMBD(RKZ,'(F7.2)',7)
  CASE(1,2)
     CALL TEXT('NPH=',4)
     CALL NUMBI(NPH,'(I3)',3)
  CASE(12)
     WW=2*PI*RF*1.0d6
     RNZ=RKZ*VC/WW
     CALL TEXT('RNZ=',4)
     CALL NUMBD(RNZ,'(F7.3)',7)
  END SELECT
  
  GXL=GXMIN
  GXL=GXL+GRCHH
  GYL=GYL+GDY
  CALL MOVE(GXL,GYL)
  CALL TEXT('NNMAX=',6)
  CALL NUMBI(NNMAX,'(I8)',8)
  
  GXL=GXL+GDX
  CALL MOVE(GXL,GYL)
  CALL TEXT('R,Z MAX=',8)
  CALL NUMBD(RNDMAX,'(F7.3)',7)
  CALL NUMBD(ZNDMAX,'(F7.3)',7)
  
  GXL=GXMIN
  GXL=GXL+GRCHH
  GYL=GYL+GDY
  CALL MOVE(GXL,GYL)
  CALL TEXT('NEMAX=',6)
  CALL NUMBI(NEMAX,'(I8)',8)
  
  GXL=GXL+GDX
  CALL MOVE(GXL,GYL)
  CALL TEXT('R,Z MIN=',8)
  CALL NUMBD(RNDMIN,'(F7.3)',7)
  CALL NUMBD(ZNDMIN,'(F7.3)',7)

!      GXL=GXL+GDX
!      CALL MOVE(GXL,GYL)
!      CALL TEXT('P=',2)
!      CALL NUMBD(TSPWR,'(1PE10.3)',10)

  GXL=GXMIN
!  GYL=GYL+GDY
  
!  CALL MOVE(GXL,GYL)
!  CALL TEXT(' NK',3)
!  CALL TEXT(' NM',3)
!  CALL TEXT(' PABS     ',10)
!  CALL TEXT(' NK',3)
!  CALL TEXT(' NM',3)
!  CALL TEXT(' PABS     ',10)
  
!  I=0
!  DO NK=1,NKMAX
!     NM=NMKA(NK)
!     IF(MOD(I,2).EQ.0) THEN
!        GXL=GXMIN
!        GYL=GYL+GDY
!     ELSE
!        GXL=GXMIN+16.*GRCHH
!     ENDIF
!     I=I+1
!     CALL MOVE(GXL,GYL)
!     CALL NUMBI(NK,'(I3)',3)
!     CALL NUMBI(NM,'(I3)',3)
!     CALL NUMBD(PABSK(NK),'(1PE10.2)',10)
!  ENDDO
  
  IF(NSMAX.GT.0) THEN
     GXL=GXMIN
     GYL=GYL+GDY
     CALL MOVE(GXL,GYL)
     CALL TEXT(' NS',3)
     CALL TEXT('    PA    ',10)
     CALL TEXT('    PZ    ',10)
     CALL TEXT('    PN    ',10)
     CALL TEXT('   PZCL   ',10)
     CALL TEXT('   PABS   ',10)
     
     DO NS=1,NSMAX
        GXL=GXMIN
        GYL=GYL+GDY
        CALL MOVE(GXL,GYL)
        CALL NUMBI(NS,'(I3)',3)
        CALL NUMBD(PA(NS),   '(1PE10.3)',10)
        CALL NUMBD(PZ(NS),   '(1PE10.3)',10)
        IF(MODELG.EQ.0.OR.MODELB.EQ.2) THEN
           CALL NUMBD(pn_corner(1,NS),'(1PE10.3)',10)
        ELSE
           CALL NUMBD(PN(NS),   '(1PE10.3)',10)
        ENDIF
        CALL NUMBD(PZCL(NS) ,'(1PE10.3)',10)
        CALL NUMBD(PABST(NS),'(1PE10.3)',10)
     ENDDO
  ENDIF
  
  GXL=GXMIN+45.*GRCHH
  GYL=GYMAX+GDY
  IF(NAMAX.GT.0) THEN
     CALL MOVE(GXL,GYL)
     CALL TEXT('IJ',2)
     CALL TEXT('  AJ ',5)
     CALL TEXT(' PHASE ',7)
     CALL TEXT('     R     ',11)
     CALL TEXT('     X     ',11)
     
     DO NA=1,NAMAX
        GXL=GXMIN+45.*GRCHH
        GYL=GYL+GDY
        CALL MOVE(GXL,GYL)
        CALL NUMBI(NA,'(I2)',2)
        CALL NUMBD(AJ(NA),'(F5.1)',5)
        CALL NUMBD(APH(NA),'(F7.1)',7)
        CALL NUMBD(REST(NA),'(1PE11.3)',11)
        CALL NUMBD(REAT(NA),'(1PE11.3)',11)
     ENDDO
  ELSEIF(NBMAX.GT.0) THEN
     CALL MOVE(GXL,GYL)
     CALL TEXT('NB',2)
     CALL TEXT(' MD',3)
     CALL TEXT('  Real   ',10)
     CALL TEXT('  Imag   ',10)
     CALL TEXT('  REFL     ',11)
     
     DO NB=1,NBMAX
        IF(KABDY(NB).GE.8) THEN
           DO L=1,NMBDY(NB)
              GXL=GXMIN+45.*GRCHH
              GYL=GYL+GDY
              CALL MOVE(GXL,GYL)
              CALL NUMBI(NB,'(I2)',2)
              CALL NUMBI(L,'(I3)',3)
              CALL NUMBD(SRFR(L,NB),'(1PE10.2)',10)
              CALL NUMBD(SRFI(L,NB),'(1PE10.2)',10)
              CALL NUMBD(SRFL(L,NB),'(1PE11.3)',11)
           ENDDO
        ENDIF
     ENDDO
  ENDIF
  
  RETURN
END SUBROUTINE wf_gdraw_parm

!     ****** Draw Vessel Boundary ******

SUBROUTINE wf_gdraw_wall

  use wfcomm
  RETURN
END SUBROUTINE wf_gdraw_wall

!     ****** Draw Plasma Boundary ******

SUBROUTINE wf_gdraw_plasma
  
  use wfcomm
  USE libgrf
  implicit none
  integer :: NPMAX,I
  real(rkind) :: DTH,THETA
  real :: GXL,GYL
  
  NPMAX=100
  DTH=2.D0*PI/NPMAX
  GXL=gdclip(RA+RR)
  GYL=0.0
  CALL MOVE2D(GXL,GYL)
  DO I=1,NPMAX
     THETA=DTH*I
     GXL=gdclip(RA*COS(THETA)+RR)
     GYL=gdclip(RA*SIN(THETA))
     CALL DRAW2D(GXL,GYL)
  END DO
  RETURN
END SUBROUTINE wf_gdraw_plasma

!     ****** Draw Antenna ******

SUBROUTINE wf_gr_antenna

  use wfcomm
  USE libgrf
  implicit none
  integer :: NTEMP
  real(rkind) :: YMIN,YMAX,XMIN,XMAX
  real(rkind) :: DXLEN,DYLEN,DRATIO,XMID,XLEN,YMID,YLEN
  REAL(rkind):: xnode_min,xnode_max,ynode_min,ynode_max

  xnode_min=RNDMIN
  xnode_max=RNDMAX
  ynode_min=ZNDMIN
  ynode_max=ZNDMAX

  DXLEN= xnode_max-xnode_min
  DYLEN= ynode_max-ynode_min
  DRATIO=DYLEN/DXLEN

  IF(DRATIO.GT.0.75D0) THEN
     YMIN=ynode_min-0.02D0*DYLEN
     YMAX=ynode_max+0.02D0*DYLEN
     XMID=0.5D0*(xnode_min+xnode_max)
     XLEN=DYLEN/0.75D0
     XMIN=XMID-0.5D0*XLEN
     XMAX=XMID+0.5D0*XLEN
  ELSE
     XMIN=xnode_min-0.02D0*DXLEN
     XMAX=xnode_max+0.02D0*DXLEN
     YMID=0.5D0*(ynode_min+ynode_max)
     YLEN=DXLEN*0.75D0
     YMIN=YMID-0.5D0*YLEN
     YMAX=YMID+0.5D0*YLEN
  ENDIF
  
  CALL PAGES
  CALL wf_gdraw_parm_ant
  CALL grd2d_frame_start(0,XMIN,XMAX,YMIN,YMAX, &
       XSCALE_ZERO=0,YSCALE_ZERO=0)
  CALL SETLIN(0,0,4)
  IF(NDRAWA.LE.1) THEN
     CALL wf_gdraw_wall
  ELSE
     NTEMP=NDRAWD
     NDRAWD=NDRAWA-1
     CALL wf_gdraw_element
     NDRAWD=NTEMP
  ENDIF
  
  CALL SETLIN(0,0,4)
  CALL wf_gdraw_plasma

  CALL SETLIN(0,0,6)
  CALL wf_gdraw_antenna
  CALL grd2d_frame_end
  
  CALL PAGEE
  RETURN
END SUBROUTINE wf_gr_antenna

SUBROUTINE wf_gdraw_antenna

  use wfcomm
  USE libgrf
  implicit none
  integer :: NA,I
  real :: GXL,GYL

  IF(NDRAWA.EQ.0) THEN
     DO NA=1,NAMAX
        GXL=gdclip(RJ0(1,NA))
        GYL=gdclip(ZJ0(1,NA))
        CALL MOVE2D(GXL,GYL)
        DO I=2,JNUM0(NA)
           GXL=gdclip(RJ0(I,NA))
           GYL=gdclip(ZJ0(I,NA))
           CALL DRAW2D(GXL,GYL)
           WRITE(6,*) I,GXL,GYL
        END DO
     END DO
  ELSE
     DO NA=1,NAMAX
        GXL=gdclip(RJ(1,NA))
        GYL=gdclip(ZJ(1,NA))
        CALL MOVE2D(GXL,GYL)
        DO I=2,JNUM(NA)
           GXL=gdclip(RJ(I,NA))
           GYL=gdclip(ZJ(I,NA))
           CALL DRAW2D(GXL,GYL)
           WRITE(6,*) I,GXL,GYL
        END DO
     END DO
  ENDIF
  RETURN
END SUBROUTINE wf_gdraw_antenna

!     ****** Draw Element Data ******

SUBROUTINE wf_gdraw_element

  use wfcomm
  USE libgrf
  implicit none
  integer :: IE,IN1,IN2,IN3,IEL,IN,INL
  real :: GX1,GX2,GX3,GY1,GY2,GY3,GXC,GYC
  
  CALL SETCHR(0.2,0.15,0.2,0.,-30.)
  
  DO IE=1,NEMAX
     IN1=NDELM(1,IE)
     IN2=NDELM(2,IE)
     IN3=NDELM(3,IE)
     GX1=gdclip(RNODE(IN1))
     GY1=gdclip(ZNODE(IN1))
     GX2=gdclip(RNODE(IN2))
     GY2=gdclip(ZNODE(IN2))
     GX3=gdclip(RNODE(IN3))
     GY3=gdclip(ZNODE(IN3))
     IF (GX1.GT.GX2) THEN
        CALL MOVE2D(GX2,GY2)
        CALL DRAW2D(GX1,GY1)
     ELSE
        CALL MOVE2D(GX1,GY1)
        CALL DRAW2D(GX2,GY2)
     END IF
     IF (GX2.GT.GX3) THEN
        CALL MOVE2D(GX3,GY3)
        CALL DRAW2D(GX2,GY2)
     ELSE
        CALL MOVE2D(GX2,GY2)
        CALL DRAW2D(GX3,GY3)
     END IF
     IF (GX3.GT.GX1) THEN
        CALL MOVE2D(GX1,GY1)
        CALL DRAW2D(GX3,GY3)
     ELSE
        CALL MOVE2D(GX3,GY3)
        CALL DRAW2D(GX1,GY1)
     END IF
     
     IF(NDRAWD.GE.2) THEN
        GXC=(GX1+GX2+GX3)/3.
        GYC=(GY1+GY2+GY3)/3.
        IEL=IE
        CALL GNUMBI2D(GXC,GYC,IEL,2)
     ENDIF
  ENDDO
  
  IF(NDRAWD.GE.3) THEN
     CALL SETCHS(0.2,0.)
     DO IN=1,NNMAX
        GX1=gdclip(RNODE(IN))
        GY1=gdclip(ZNODE(IN))
        INL=IN
        CALL GNUMBI2D(GX1,GY1,INL,0)
     END DO
  ENDIF
  
  RETURN
END SUBROUTINE wf_gdraw_element

!     ****** Draw Element Data ******

SUBROUTINE WFGNAS!(ID)
  
  use wfcomm
  USE libgrf
  implicit none
!  integer :: ID
  real(rkind) :: PXMIN,PXMAX,PYMIN,PYMAX
  real :: GPXMIN,GPXMAX,GPYMIN,GPYMAX,GYMIN,GYMAX,GXMIN,GXMAX
  real(rkind) :: DXLEN,DYLEN,DRATIO,PRATIO,XMID,XLEN,YMID,YLEN
  
  CALL WFGWIN(1,1,PXMIN,PXMAX,PYMIN,PYMAX)
  GPXMIN=gdclip(PXMIN)
  GPXMAX=gdclip(PXMAX)
  GPYMIN=gdclip(PYMIN)+1.0
  GPYMAX=gdclip(PYMAX)+1.0
  
  DXLEN= RNDMAX-RNDMIN+0.5D0*(ZNDMAX-ZNDMIN)
  !!!! TO BE MODIFIED?
  DYLEN=               0.5D0*(ZNDMAX-ZNDMIN)
  DRATIO=DYLEN/DXLEN
  PRATIO=(PYMAX-PYMIN)/(PXMAX-PXMIN)
  IF(DRATIO.GT.PRATIO) THEN
     GYMIN=gdclip(0.5D0*ZNDMIN)
     GYMAX=gdclip(0.5D0*ZNDMAX)
     XMID=0.5D0*(RNDMIN+RNDMAX+0.5D0*(ZNDMAX+ZNDMIN))
     XLEN=DYLEN/PRATIO
     GXMIN=gdclip(XMID-0.5D0*XLEN)
     GXMAX=gdclip(XMID+0.5D0*XLEN)
  ELSE
     GXMIN=gdclip(RNDMIN+0.5D0*ZNDMIN)
     GXMAX=gdclip(RNDMAX+0.5D0*ZNDMIN)
     YMID=0.5D0*(0.5D0*(ZNDMAX+ZNDMIN))
     YLEN=DXLEN*PRATIO
     GYMIN=gdclip(YMID-0.5D0*YLEN)
     GYMAX=gdclip(YMID+0.5D0*YLEN)
  ENDIF
  
  CALL PAGES
  CALL wf_gdraw_element
  CALL SETVEW(GPXMIN,GPXMAX,GPYMIN,GPYMAX,&
       &      GXMIN ,GXMAX ,GYMIN ,GYMAX )
  !     &            GXMIN,GXMAX,GYMIN,GYMAX)
  !     &     -0.03,-0.01,-0.005,0.005)
  !     &     -0.03,0.03,-0.02,0.02)
  
  CALL SETLIN(0,0,7)
!  CALL WFGELM3(ID)
  
  CALL PAGEE
  RETURN
END SUBROUTINE WFGNAS

!     ****** Draw Element Data ******

SUBROUTINE wf_gr_element

  use wfcomm
  USE libgrf
  implicit none
  real(rkind) :: YMIN,YMAX,XMIN,XMAX
  real(rkind) :: DXLEN,DYLEN,DRATIO,XMID,XLEN,YMID,YLEN
  REAL(rkind):: xnode_min,xnode_max,ynode_min,ynode_max

  xnode_min=RNDMIN
  xnode_max=RNDMAX
  ynode_min=ZNDMIN
  ynode_max=ZNDMAX

  
  DXLEN= xnode_max-xnode_min
  DYLEN= ynode_max-ynode_min
  DRATIO=DYLEN/DXLEN

  IF(DRATIO.GT.0.75D0) THEN
     YMIN=ynode_min
     YMAX=ynode_max
     XMID=0.5D0*(xnode_min+xnode_max)
     XLEN=DYLEN/0.75D0
     XMIN=XMID-0.5D0*XLEN
     XMAX=XMID+0.5D0*XLEN
  ELSE
     XMIN=xnode_min
     XMAX=xnode_max
     YMID=0.5D0*(ynode_min+ynode_max)
     YLEN=DXLEN*0.75D0
     YMIN=YMID-0.5D0*YLEN
     YMAX=YMID+0.5D0*YLEN
  ENDIF
  WRITE(6,'(A,5ES12.4)') 'x,y,R:',XMIN,XMAX,YMIN,YMAX
  
  CALL PAGES
  CALL wf_gdraw_parm_elm
  CALL grd2d_frame_start(0,XMIN,XMAX,YMIN,YMAX, &
       XSCALE_ZERO=0,YSCALE_ZERO=0)

  IF(NDRAWD.EQ.0) THEN
     CALL wf_gdraw_wall
  ELSE
     CALL wf_gdraw_element
  ENDIF
  
  CALL grd2d_frame_end
  CALL PAGEE
  RETURN
END SUBROUTINE wf_gr_element

!     ****** Draw Element Paramters ******

SUBROUTINE wf_gdraw_parm_elm

  use wfcomm
  implicit none
  real :: GXMIN,GYMAX,GDY,GXL,GYL
  
  GXMIN=4.0
  GYMAX=17.8
  CALL SETCHS(0.25,0.)
  GDY=0.3
  
  GXL=GXMIN
  GYL=GYMAX
  CALL MOVE(GXL,GYL)
  CALL TEXT('node_max=',9)
  CALL NUMBI(NNMAX,'(I8)',8)
  
  GYL=GYL-GDY
  CALL MOVE(GXL,GYL)
  CALL TEXT('nelm_max=',9)
  CALL NUMBI(NEMAX,'(I8)',8)
  
  GYL=GYL-GDY
  CALL MOVE(GXL,GYL)
  CALL TEXT('mtx_len= ',9)
  CALL NUMBI(MLEN,'(I8)',8)
  RETURN
END SUBROUTINE wf_gdraw_parm_elm

!     ****** Draw Antenna Paramters ******

SUBROUTINE wf_gdraw_parm_ant

  use wfcomm
  implicit none
  integer :: NA
  real :: GXMIN,GYMAX,GDY,GXL,GYL

  GXMIN=4.0
  GYMAX=17.8
  CALL SETCHS(0.25,0.)
  GDY=0.3
  
  GXL=GXMIN
  GYL=GYMAX
  CALL MOVE(GXL,GYL)
  CALL TEXT('JNUM0=',6)
  GXL=GXL+1.5
  DO NA=1,NAMAX
     CALL MOVE(GXL,GYL)
     CALL NUMBI(JNUM0(NA),'(I4)',4)
     GXL=GXL+1.0
  END DO
  GXMIN=4.0
  GYMAX=17.8
  CALL SETCHS(0.25,0.)
  GDY=0.3
  
  GXL=GXMIN
  GYL=GYL-GDY
  CALL MOVE(GXL,GYL)
  CALL TEXT('JNUM =',6)
  GXL=GXL+1.5
  DO NA=1,NAMAX
     CALL MOVE(GXL,GYL)
     CALL NUMBI(JNUM(NA),'(I4)',4)
     GXL=GXL+1.0
  END DO
  RETURN
END SUBROUTINE wf_gdraw_parm_ant

!     ****** Draw Waveguide ******

SUBROUTINE wf_gr_waveguide

  use wfcomm
  USE libgrf
  implicit none
  integer :: NTEMP
  real(rkind) :: YMIN,YMAX,XMIN,XMAX
  real(rkind) :: DXLEN,DYLEN,DRATIO,XMID,XLEN,YMID,YLEN
  REAL(rkind):: xnode_min,xnode_max,ynode_min,ynode_max

  xnode_min=RNDMIN
  xnode_max=RNDMAX
  ynode_min=ZNDMIN
  ynode_max=ZNDMAX

  DXLEN= xnode_max-xnode_min
  DYLEN= ynode_max-ynode_min
  DRATIO=DYLEN/DXLEN

  IF(DRATIO.GT.0.75D0) THEN
     YMIN=ynode_min-0.02D0*DYLEN
     YMAX=ynode_max+0.02D0*DYLEN
     XMID=0.5D0*(xnode_min+xnode_max)
     XLEN=DYLEN/0.75D0
     XMIN=XMID-0.5D0*XLEN
     XMAX=XMID+0.5D0*XLEN
  ELSE
     XMIN=xnode_min-0.02D0*DXLEN
     XMAX=xnode_max+0.02D0*DXLEN
     YMID=0.5D0*(ynode_min+ynode_max)
     YLEN=DXLEN*0.75D0
     YMIN=YMID-0.5D0*YLEN
     YMAX=YMID+0.5D0*YLEN
  ENDIF
  
!  CALL PAGES
  CALL wf_gdraw_parm_waveguide
  CALL grd2d_frame_start(0,XMIN,XMAX,YMIN,YMAX, &
       XSCALE_ZERO=0,YSCALE_ZERO=0)
  
  CALL SETLIN(0,0,4)
  IF(NDRAWA.LE.1) THEN
     CALL wf_gdraw_wall
  ELSE
     NTEMP=NDRAWD
     NDRAWD=NDRAWA-1
     CALL wf_gdraw_element
     NDRAWD=NTEMP
  ENDIF
  
  CALL SETLIN(0,0,5)
  CALL wf_gdraw_plasma
  
  CALL SETLIN(0,0,6)
  CALL wf_gdraw_waveguide
  CALL grd2d_frame_end
  
!  CALL PAGEE
  RETURN
END SUBROUTINE wf_gr_waveguide

!     ****** Draw Waveguide ******

SUBROUTINE wf_gdraw_waveguide

  use wfcomm
  USE libgrf
  implicit none
  real :: GXL,GYL
  REAL(rkind):: xwg_min,xwg_max,ywg_min,ywg_max

  xwg_min=R1WG
  xwg_max=R2WG
  ywg_min=Z1WG
  ywg_max=Z2WG
  
  GXL=gdclip(xwg_min)
  GYL=gdclip(ywg_min)
  CALL MOVE2D(GXL,GYL)
  GXL=gdclip(xwg_max)
  CALL DRAW2D(GXL,GYL)
  GYL=gdclip(ywg_max)
  CALL DRAW2D(GXL,GYL)
  GXL=gdclip(xwg_min)
  CALL DRAW2D(GXL,GYL)
  GYL=gdclip(ywg_min)
  CALL DRAW2D(GXL,GYL)
  RETURN
END SUBROUTINE wf_gdraw_waveguide

!     ****** Draw Waveguide Paramters ******

SUBROUTINE wf_gdraw_parm_waveguide

  use wfcomm
  implicit none
  real :: GXMIN,GYMAX,GDY,GXL,GYL

  GXMIN=4.0
  GYMAX=17.8
  CALL SETCHS(0.25,0.)
  GDY=0.3
  
  GXL=GXMIN
  GYL=GYMAX
  CALL MOVE(GXL,GYL)
  CALL TEXT('RF   =',6)
  CALL NUMBD(RF,'(ES10.3)',10)
  
  GYL=GYL-GDY
  CALL MOVE(GXL,GYL)
  CALL TEXT('NPH  =',6)
  CALL NUMBI(NPH,'(I10)',10)
  
END SUBROUTINE wf_gdraw_parm_waveguide

SUBROUTINE WFGMESH
  USE wfcomm
  USE libgrf
  IMPLICIT NONE
  INTEGER,SAVE:: NGXMAX_SAVE=0,NGYMAX_SAVE=0
  REAL(rkind),SAVE:: RNDMIN_SAVE=0.D0,RNDMAX_SAVE=0.D0
  REAL(rkind),SAVE:: ZNDMIN_SAVE=0.D0,ZNDMAX_SAVE=0.D0
  REAL(rkind):: DX,DY,X,Y
  INTEGER:: NGX,NGY,IE

  IF(NGXMAX.NE.NGXMAX_SAVE.OR. &
     NGYMAX.NE.NGYMAX_SAVE.OR. &
     RNDMIN.NE.RNDMIN_SAVE.OR. &
     RNDMAX.NE.RNDMAX_SAVE.OR. &
     ZNDMIN.NE.ZNDMIN_SAVE.OR. &
     ZNDMAX.NE.ZNDMAX_SAVE) THEN
     IE=0
     DY=(ZNDMAX-ZNDMIN)/(NGYMAX-1)
     DX=(RNDMAX-RNDMIN)/(NGXMAX-1)
     DO NGX=1,NGXMAX
        G2X(NGX)=gdclip(RNDMIN+DX*(NGX-1))
     ENDDO
     DO NGY=1,NGYMAX
        G2Y(NGY)=gdclip(ZNDMIN+DY*(NGY-1))
     ENDDO
     DO NGY=1,NGYMAX
        Y=ZNDMIN+DY*(NGY-1)
        DO NGX=1,NGXMAX
           X=RNDMIN+DX*(NGX-1)
           CALL FEP(X,Y,IE)
           IEGZ(NGX,NGY)=IE
        END DO
     END DO
     NGXMAX_SAVE=NGXMAX
     NGYMAX_SAVE=NGYMAX
     RNDMIN_SAVE=RNDMIN
     RNDMAX_SAVE=RNDMAX
     ZNDMIN_SAVE=ZNDMIN
     ZNDMAX_SAVE=ZNDMAX
  END IF
  RETURN
END SUBROUTINE WFGMESH
