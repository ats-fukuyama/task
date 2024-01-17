! wfgout.f90

MODULE wfgout
  
  PRIVATE
  PUBLIC wf_gout
  PUBLIC wf_gpwr
  PUBLIC wf_gden
  PUBLIC wf_g1d_plot
  PUBLIC wf_g1d_binary
  PUBLIC wf_g1d_text
  PUBLIC wf_g2d_contour
  PUBLIC wf_g2d_paint
  PUBLIC wf_g2d_birdeye
  PUBLIC wf_g2d_binary
  PUBLIC wf_g2d_text
  PUBLIC wf_g2d_vector
  
CONTAINS

!     *********  /TASKX/WFT/GOUT  *********
!     $Id: wfgout.f90,v 1.24 2012/02/11 01:23:26 maruyama Exp $
!
!       GRAPHIC DATA PROCESSING PROGRAM
!             FOR FEM COMPUTATION
!
!     ***********************************

SUBROUTINE wf_gout

  use wfcomm
  USE wfgsub
  USE wfsub
  USE libchar
  implicit none

  integer :: NL,NWD,NCH,NWMAX,NW,idraw_parm
  CHARACTER KLINE*80,KWORD*(NCHM),KWD*(NCHM),KID*1,KTAIL*7
  CHARACTER KG1*1,KG2*1
  DIMENSION KWORD(NWDM)

  call wf_win_allocate
  CALL wf_set_node_range

  if(NDRAWV.eq.1)  call wf_g2d_vector

  IF(xnode_max.EQ.xnode_min) THEN
     WRITE(6,*) 'XX NO DATA IS LOADED FOR GRAPHICS'
     GOTO 9000
  ENDIF
  
1 WRITE(6,*) '# INPUT : E,X/Y/Z,R/I/A'
  WRITE(6,*) '          Xyz/Yxz/Zxy/A9 0-9 V,0-9, L'
  WRITE(6,*) '          P,1/2  N  A  W  X=EXIT'
  CALL GUFLSH
  READ(5,'(A80)',ERR=1,END=9000) KLINE
  NWXMAX=0
  idraw_parm=0
  
9 NL=0
  NWD=0
  NCH=0
10 IF(NL.GE.80) GOTO 20
  NL=NL+1
  KID=KLINE(NL:NL)
  CALL toupper(KID)
  IF(KID.NE.' ') THEN
     IF(NCH.LT.NCHM) NCH=NCH+1
     KWD(NCH:NCH)=KID
  ELSE
     IF(NCH.NE.0) THEN
        IF(NWD.LT.NWDM) NWD=NWD+1
        KWORD(NWD)=KWD(1:NCH)
        NCH=0
     ENDIF
  ENDIF
  GOTO 10
  
20 CONTINUE
  NWMAX=NWD
  KWD=KWORD(1)
  KG1=KWD(1:1)
  KG2=KWD(2:2)
  IF(KG1.EQ.'M') THEN
     DO NCH=3,9
        IF(KWD(NCH:NCH).EQ.' ') GOTO 30
     ENDDO
     NCH=NCH-1
30   KTAIL=KWD(3:NCH-1)//' '
!     WRITE(6,*) 'KTAIL = /',KTAIL,'/'
     
     IF(KG2.EQ.'E') THEN
        if(NDRAWE.eq.0) then
           KLINE='ERR'//KTAIL//'ERI'//KTAIL//&
                &'EZR'//KTAIL//'EZI'//KTAIL//&
                &'EPR'//KTAIL//'EPI'//KTAIL//&
                &'P1C'//KTAIL//'P2C'//KTAIL
        elseif(NDRAWE.eq.1) then
           KLINE='ERR'//KTAIL//'ERI'//KTAIL//&
                &'ETR'//KTAIL//'ETI'//KTAIL//&
                &'EPR'//KTAIL//'EPI'//KTAIL//&
                &'P1C'//KTAIL//'P2C'//KTAIL
        end if
        GOTO 9
     ELSE
        WRITE(6,*) 'XX UNKNOWN KID2:',KG2
        GOTO 1
     ENDIF
  ELSEIF(KWD(1:1).EQ.'X') THEN
     GOTO 9000
  ENDIF
  
  IF(NWMAX.EQ.0) GOTO 1

  IF(NGRAPH.GE.1) CALL PAGES
  DO NW=1,NWMAX
     KWD=KWORD(NW)
     WRITE(6,*) 'KWD=',KWD(1:4)
     KID=KWD(1:1)
     IF(    KID.EQ.'E') THEN
        KID=KWD(2:2)
        IF(    KID.EQ.'R') THEN
           if(NDRAWE.eq.0) CALL wf_ctog_side(1,KWD)  ! E_R    (R,Z,psi)
           if(NDRAWE.eq.1) CALL wf_ctog_side(4,KWD)  ! E_rho  (rho,theta,-psi)
        ELSEIF(KID.EQ.'P') THEN
           if(NDRAWE.eq.0) CALL wf_ctog_node(2,KWD)  ! E_psi
           if(NDRAWE.eq.1) CALL wf_ctog_node(5,KWD)  ! -E_psi
        ELSEIF(KID.EQ.'Z') THEN
           CALL wf_ctog_side(3,KWD)                  ! E_Z
        ELSEIF(KID.eq.'T') THEN
           CALL wf_ctog_side(6,KWD)                  ! E_theta
        ELSE
           WRITE(6,*) 'XX UNKNOWN KID2:',KID
           GOTO 1000
        ENDIF
        idraw_parm=1

     ELSE IF(KID.eq.'P') THEN
        KID=KWD(2:2)
        IF(KID.eq.'1') THEN
           call wf_gpwr(1)
        ELSE IF(KID.eq.'2') THEN
           call wf_gpwr(2)
        ELSE IF(KID.eq.'3') THEN
           call wf_gpwr(3)
        ELSE
           WRITE(6,*) 'XX UNKNOWN KID2:',KID
           GOTO 1000
        ENDIF
        idraw_parm=1

     ELSE IF(KID.EQ.'N') THEN
        CALL wf_gden
     ELSE IF(KID.EQ.'A') THEN
        CALL wf_gr_antenna
     ELSE IF(KID.EQ.'W') THEN
        CALL wf_gr_waveguide
     ELSE
        WRITE(6,*) 'XX UNKNOWN KID1:',KID
        GOTO 1000
     ENDIF
     KID=KWD(3:3)
     IF(    KID.EQ.'R'.OR.&
       &    KID.EQ.'I'.OR.&
       &    KID.EQ.'A'.OR.&
       &    KID.EQ.'C') THEN
        SELECT CASE(NGRAPH)
        CASE(-1)
           CALL wf_g2d_binary(KWD)
        CASE(0)
           CALL wf_g2d_text(KWD)
        CASE(1)
           CALL wf_g2d_contour(NW,NWMAX,KWD)
        CASE(2)
           CALL wf_g2d_paint(NW,NWMAX,KWD)
        CASE(3:6)
           CALL wf_g2d_birdeye(NW,NWMAX,KWD)
        END SELECT
     ELSEIF(KID.EQ.'X'.OR.&
          & KID.EQ.'Y') THEN
        SELECT CASE(NGRAPH)
        CASE(-1)
           CALL wf_g1d_binary(KWD)
        CASE(0)
           CALL wf_g1d_text(KWD)
        CASE(1:6)
           CALL wf_g1d_plot(NW,NWMAX,KWD)
        END SELECT
     ELSE
        WRITE(6,*) 'XX UNKNOWN KID3:',KID
     ENDIF
  END DO

1000 continue
  IF(NGRAPH.GE.1) THEN
     IF(idraw_parm.EQ.1) CALL wf_gdraw_parm
     CALL PAGEE
  ENDIF
  GOTO 1     
  
9000 RETURN
END SUBROUTINE wf_gout

!     ****** PLOT POWER ABSORPTION ******

SUBROUTINE wf_gpwr(NS)

  use wfcomm
  USE wfwave
  USE wfsolv
  USE wfsub,ONLY: wf_fieldcr,wf_fieldcz,wf_fieldcp,wf_set_weight
  USE wfgsub
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

  CALL wf_gen_elm_mesh
  DO NGY=1,NGYMAX
     Y=G2Y(NGY)
     DO NGX=1,NGXMAX
        X=G2X(NGX)
        IE=IEGZ(NGX,NGY)
        IF(IE.EQ.0) THEN
           GZ(NGX,NGY)=0.0
        ELSE
           ! --- calculate conductivity tensor ---

           CALL wf_dtensr(IE,DTENS)
           call wf_set_weight(IE,X,Y,WGT)

           CTENS=(0.d0,0.d0)

           do J=1,3
              do I=1,3
                 do IN=1,3
                    CTENS(I,J)=CTENS(I,J)-WGT(IN)*CIWE*DTENS(NS,IN,I,J)
                 end do
              end do
           end do

           ! --- calculate power absorption ---
           call wf_fieldcr(IE,X,Y,CESD_nseg,CER)
           call wf_fieldcp(IE,X,Y,CEND_node,CEP)
           call wf_fieldcz(IE,X,Y,CESD_nseg,CEZ)

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
END SUBROUTINE wf_gpwr

!    ***** PLOT DENSITY PROFILE *****

SUBROUTINE wf_gden

  use wfcomm
  USE wfprof
  USE wfindex
  USE libgrf
  implicit none
  
  integer :: NGX,NGY,NE
  real(rkind) :: DX,DY,X,Y
  real(rkind) :: RN(NSM),RTPR(NSM),RTPP(NSM),RZCL(NSM)
  
  NE=0

  DY=(ynode_max-ynode_min)/(NGYMAX-1)
  DX=(xnode_max-xnode_min)/(NGXMAX-1)
  DO NGX=1,NGXMAX
     G2X(NGX)=gdclip(xnode_min+DX*(NGX-1))
  ENDDO
  DO NGY=1,NGYMAX
     G2Y(NGY)=gdclip(ynode_min+DY*(NGY-1))
  ENDDO

  DO NGY=1,NGYMAX
     Y=ynode_min+DY*(NGY-1)
     DO NGX=1,NGXMAX
        X=xnode_min+DX*(NGX-1)
        CALL wf_fep(X,Y,NE)
        IF(NE.EQ.0) THEN
           GZ(NGX,NGY)=0.0
        ELSE
           CALL wf_sden(X,Y,RN,RTPR,RTPP,RZCL)
           GZ(NGX,NGY)=gdclip(RN(1)*1.d20)
        end IF
     end DO
  end DO

  return
END SUBROUTINE wf_gden

!     ****** DRAW 1D PROFILE ******

SUBROUTINE wf_g1d_plot(NW,NWMAX,KWD)
  
  use wfcomm
  USE wfgsub
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
  
  CALL wf_gwin_range(NW,NWMAX,PXMIN,PXMAX,PYMIN,PYMAX)

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
     CALL TEXT('(X): Y=',9)
  ELSE IF(KWD(3:3).EQ.'Y') THEN
     CALL TEXT('(Y): X=',9)
  ENDIF
  CALL TEXT(KWD(4:NCHM),NCHM-3)
  RETURN
END SUBROUTINE wf_g1d_plot

!     ****** WRITE 1D PROFILE IN BINARY FILE ******

SUBROUTINE wf_g1d_binary(KWD)

  use wfcomm
  implicit none
  integer :: NGMAX,NFD,NGV,NG
  character,intent(in) :: KWD*(NCHM)

  IF(KWD(1:1).EQ.'E'.OR.&
 &   KWD(1:1).EQ.'D'.OR.&
 &   KWD(1:1).EQ.'B'.OR.&
 &   KWD(1:1).EQ.'A') THEN
     NGMAX=3
  ELSE
     NGMAX=1
  ENDIF

  NFD=23
  WRITE(NFD) KWD
  WRITE(NFD) 1
  WRITE(NFD) NGVMAX,NGMAX
  WRITE(NFD) (GX(NGV),NGV=1,NGVMAX)
  WRITE(NFD) ((GV(NGV,NG),NGV=1,NGVMAX),NG=1,NGMAX)

  RETURN
END SUBROUTINE wf_g1d_binary

!     ****** WRITE 1D PROFILE IN TEXT FILE ******

SUBROUTINE wf_g1d_text(KWD)

  use wfcomm
  implicit none
  integer :: NGMAX,NFD,NGV,NG
  character,intent(in) :: KWD*(NCHM)
  
  IF(KWD(1:1).EQ.'E'.OR.&
 &   KWD(1:1).EQ.'D'.OR.&
 &   KWD(1:1).EQ.'B'.OR.&
 &   KWD(1:1).EQ.'A') THEN
     NGMAX=3
  ELSE
     NGMAX=1
  ENDIF
  
  NFD=22
  WRITE(NFD,'(A79)') KWD
  WRITE(NFD,'(I8)') 1
  WRITE(NFD,'(2I8)') NGVMAX,NGMAX
  WRITE(NFD,'(1P5E15.7)') (GX(NGV),NGV=1,NGVMAX)
  WRITE(NFD,'(1P5E15.7)') ((GV(NGV,NG),NGV=1,NGVMAX),NG=1,NGMAX)
  
  RETURN
END SUBROUTINE wf_g1d_text

!     ****** PLOT 2D STRUCTURE OF FIELD AND POWER ******

SUBROUTINE wf_g2d_contour(NW,NWMAX,KWD)
  
  use wfcomm
  USE wfprof
  USE wfgsub
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
  INTEGER :: nbdy,nseg,node1,node2
  REAL :: X1,Y1,X2,Y2

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

  CALL wf_gwin_range(NW,NWMAX,PXMIN,PXMAX,PYMIN,PYMAX)
  
  XMIN = xnode_min
  XMAX = xnode_max
  YMIN = ynode_min
  YMAX = ynode_max  
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
        Y=ynode_min+DY*(NGY-1)
        DO NGX=1,NGXMAX
           X=xnode_min+DX*(NGX-1)

           CALL wf_smag(X,Y,BABS,AL)
           CALL wf_sden(X,Y,RN,RTPR,RTPP,RZCL)

           DO NSDO=1,NSMAX
              WP(NSDO)=PZ(NSDO)*PZ(NSDO)*AEE*AEE*RN(NSDO)*1.D20 &
                            /(PA(NSDO)*AMP*EPS0*WW*WW)
              WC(NSDO)=PZ(NSDO)*AEE*BABS/(PA(NSDO)*AMP*WW)
           end DO

           TCR(NGX,NGY)=WC(NS)**2
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
  GZDEL=REAL(gaspect)*GZSCAL
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
  DO nbdy=1,nbdy_max
     nseg=nseg_nbdy(nbdy)
     node1=node_nseg(1,nseg)
     node2=node_nseg(2,nseg)
     X1=guclip(xnode(node1))
     Y1=guclip(ynode(node1))
     Y2=guclip(xnode(node2))
     Y2=guclip(ynode(node2))
     CALL MOVE2D(X1,Y1)
     CALL DRAW2D(X2,Y2)
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
END SUBROUTINE wf_g2d_contour

!     ****** PLOT 2D STRUCTURE OF FIELD AND POWER ******

SUBROUTINE wf_g2d_paint(NW,NWMAX,KWD)

  use wfcomm
  USE wfgsub
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
  
  CALL wf_gwin_range(NW,NWMAX,PXMIN,PXMAX,PYMIN,PYMAX)
  
  XMIN = xnode_min
  XMAX = xnode_max
  YMIN = ynode_min
  YMAX = ynode_max
  
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
END SUBROUTINE wf_g2d_paint

!     ****** PLOT 2D STRUCTURE OF FIELD AND POWER ******

SUBROUTINE wf_g2d_birdeye(NW,NWMAX,KWD)
  
  use wfcomm
  USE wfgsub
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
  CALL wf_gwin_range(NW,NWMAX,PXMIN,PXMAX,PYMIN,PYMAX)
  
  XMIN = xnode_min
  XMAX = xnode_max
  YMIN = ynode_min
  YMAX = ynode_max
  
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
END SUBROUTINE wf_g2d_birdeye

!     ****** WRITE 2D PROFILE IN BINARY FILE ******

SUBROUTINE wf_g2d_binary(KWD)

  use wfcomm
  implicit none
  integer :: NFD,NGX,NGY
  character,intent(in) :: KWD*(NCHM)

  NFD=23
  WRITE(NFD) KWD
  WRITE(NFD) 2
  WRITE(NFD) NGXMAX,NGYMAX
  WRITE(NFD) (G2X(NGX),NGX=1,NGXMAX)
  WRITE(NFD) (G2Y(NGY),NGY=1,NGYMAX)
  WRITE(NFD) ((GZ(NGX,NGY),NGX=1,NGXMAX),NGY=1,NGYMAX)

  RETURN
END SUBROUTINE wf_g2d_binary

!     ****** WRITE 2D PROFILE IN TEXT FILE ******

SUBROUTINE wf_g2d_text(KWD)

  use wfcomm
  implicit none
  integer :: NFD,NGX,NGY
  character,intent(in) :: KWD*(NCHM)

  NFD=22
  WRITE(NFD,'(A79)') KWD
  WRITE(NFD,'(I8)') 2
  WRITE(NFD,'(2I8)') NGXMAX,NGYMAX
  WRITE(NFD,'(1P5E15.7)') (G2X(NGX),NGX=1,NGXMAX)
  WRITE(NFD,'(1P5E15.7)') (G2Y(NGY),NGY=1,NGYMAX)
  WRITE(NFD,'(1P5E15.7)') ((GZ(NGX,NGY),NGX=1,NGXMAX),NGY=1,NGYMAX)

  RETURN
END SUBROUTINE wf_g2d_text

! ---- output text data (2D vector field) ----
! This subroutine is created in order to see waveguide eigenmode
! as vector field.
! Output file includes only Er and Ez, not E_phi.

SUBROUTINE wf_g2d_vector

  use wfcomm
  USE wfindex
  USE wfsub
  USE libgrf
  implicit none

  integer :: nelm,NGX,NGY
  complex(rkind):: CE
  real(rkind):: DX,DY,X,Y
  real(rkind),dimension(:,:),ALLOCATABLE::GZ_x,GZ_y

  allocate(GZ_x(NGXMAX,NGYMAX),GZ_y(NGXMAX,NGYMAX))
  nelm=0

  DY=(ynode_max-ynode_min)/(NGYMAX-1)
  DX=(xnode_max-xnode_min)/(NGXMAX-1)
  DO NGX=1,NGXMAX
     G2X(NGX)=gdclip(xnode_min+DX*(NGX-1))
  ENDDO
  DO NGY=1,NGYMAX
     G2Y(NGY)=gdclip(ynode_min+DY*(NGY-1))
  ENDDO
  DO NGY=1,NGYMAX
     Y=ynode_min+DY*(NGY-1)
     DO NGX=1,NGXMAX
        X=xnode_min+DX*(NGX-1)
        CALL wf_fep(X,Y,nelm)
        IF(nelm.EQ.0) THEN
           GZ_x(NGX,NGY)=0.0
           GZ_y(NGX,NGY)=0.0
        ELSE
           CALL wf_fieldcr(nelm,X,Y,CESD_nseg,CE)
           GZ_x(NGX,NGY)=gdclip(AIMAG(CE))
           CALL wf_fieldcz(nelm,X,Y,CESD_nseg,CE)
           GZ_y(NGX,NGY)=gdclip(AIMAG(CE))
       ENDIF
     END DO
  ENDDO
  
  open(100,file="vfield")
  write(100,'(7X,A2,15X,A2,15X,A2,15X,A2)') " x"," y","Ex","Ey"
  DO NGY=1,NGYMAX
     Y=ynode_max+DY*(NGY-1)
     DO NGX=1,NGXMAX
        X=xnode_max+DX*(NGX-1)
        write(100,'(E16.8,1X,E16.8,1X,E16.8,1X,E16.8)') &
             X,Y,GZ_x(NGX,NGY),GZ_y(NGX,NGY)
     end DO
  end DO
  close(100)

  deallocate(GZ_x,GZ_y)
  return
end subroutine wf_g2d_vector

END MODULE wfgout
