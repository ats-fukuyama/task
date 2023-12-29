! wfgsub.f90

MODULE wfgsub

  PRIVATE
  PUBLIC wf_ctog_side
  PUBLIC wf_ctog_node
  PUBLIC wf_gwin_range
  PUBLIC wf_gdraw_parm
  PUBLIC wf_gr_element
  PUBLIC wf_gr_antenna
  PUBLIC wf_gr_waveguide
  PUBLIC wf_gen_elm_mesh
  
CONTAINS

!     ****** EXTRACT FROM COMPLEX DATA ******

SUBROUTINE wf_ctog_side(ID,KWD)

  use wfcomm
  USE wfsub
  use feminterpolate
  USE libgrf
  implicit none
  integer,intent(in) :: ID
  integer :: IE,NGX,NGY,NGV
  real(rkind) :: DX,DY,X,Y
  real(rkind) :: XPOS,YPOS
  complex(rkind) :: CE
  character,intent(in) :: KWD*(NCHM)
  CHARACTER  :: KID*1

  real(rkind) :: theta
  complex(rkind)::CE_R,CE_Z

  KID=KWD(3:3)
  IE=0

! --- EXTRACT DATA FOR 2D PLOT ---  

  if(KID.eq.'R'.or.&
     KID.eq.'I'.or.&
     KID.eq.'C') then
     CALL wf_gen_elm_mesh
     DO NGY=1,NGYMAX
        Y=G2Y(NGY)
        DO NGX=1,NGXMAX
           X=G2X(NGX)
           !           IE=IEGZ(NGX,NGY)
           CALL fem_find_nelm_for_xy(x,y,ie)
           IF(IE.EQ.0) THEN
              GZ(NGX,NGY)=0.0
           ELSE
              IF(ID.eq.1) THEN
                 CALL wf_fieldcr(IE,X,Y,CESD,CE)
              ELSE IF(ID.eq.3) THEN
                 CALL wf_fieldcz(IE,X,Y,CESD,CE)
              ELSEIF(ID.eq.4.or.ID.eq.6) THEN
                 theta=datan2(Y,(X-RR))
                 CALL wf_fieldcr(IE,X,Y,CESD,CE)
                 CE_R=CE
                 CALL wf_fieldcz(IE,X,Y,CESD,CE)
                 CE_Z=CE
                 if(ID.eq.4) then
                    CE=-CE_R*cos(theta)-CE_Z*sin(theta)
                 elseif(ID.eq.6) then
                    CE= CE_R*sin(theta)-CE_Z*cos(theta)
                 end if
              ENDIF
              IF(ABS(CE).LT.1.D-12) CE=(0.D0,0.D0)
!              IF(Y.GT.0.3D0.AND.Y.LT.0.5D0) THEN
!              N1=nseg_nside_nelm(1,IE)
!              IF(N1.GT.0) THEN
!                 CE1=CESD(N1)
!              ELSE
!                 CE1=-CESD(-N1)
!              END IF
!              N2=nseg_nside_nelm(2,IE)
!              IF(N2.GT.0) THEN
!                 CE2=CESD(N2)
!              ELSE
!                 CE2=-CESD(-N2)
!              END IF
!              N3=nseg_nside_nelm(3,IE)
!              IF(N3.GT.0) THEN
!                 CE3=CESD(N3)
!              ELSE
!                 CE3=-CESD(-N3)
!              END IF
             
!              IF(ABS(CE).NE.0.D0) THEN
!                 WRITE(6,'(A,I12,1P5E12.4)') 'CESD:',IE,X,Y,CE
!                 WRITE(6,'(A,1P6E12.4)')     '     ', CE1,CE2,CE3
!                 WRITE(6,'(A,1P6E12.4)')     '     ', &
!                      xnode(node_nside_nelm(1,IE)),ynode(node_nside_nelm(1,IE)), &
!                      xnode(node_nside_nelm(2,IE)),ynode(node_nside_nelm(2,IE)), &
!                      xnode(node_nside_nelm(3,IE)),ynode(node_nside_nelm(3,IE))
!              END IF
!              END IF
                 
              ! -----------
              IF(KID.EQ.'R') THEN
                 GZ(NGX,NGY)=gdclip(REAL(CE))
              ELSE IF(KID.EQ.'I') THEN
                 GZ(NGX,NGY)=gdclip(AIMAG(CE))
              ELSE IF(KID.EQ.'A') THEN
                 GZ(NGX,NGY)=gdclip(ABS(CE))
              ENDIF
           ENDIF
        END DO
     ENDDO

! --- EXTRACT DATA FOR 1D X PLOT ---

  elseif(KID.eq.'X') then
     READ(KWD(4:NCHM),*,ERR=9000) YPOS
     DX=(RNDMAX-RNDMIN)/(NGVMAX-1)
     DO NGV=1,NGVMAX
        X=RNDMIN+DX*(NGV-1)
        CALL fem_find_nelm_for_xy(x,ypos,ie)
        IF(IE.EQ.0) THEN
           GX(NGV)=gdclip(X)
           GV(NGV,1)=0.0
           GV(NGV,2)=0.0
           GV(NGV,3)=0.0
        ELSE
           IF(ID.eq.1) THEN
              CALL wf_fieldcr(IE,X,YPOS,CESD,CE)
           ELSE IF(ID.eq.3) THEN
              CALL wf_fieldcz(IE,X,YPOS,CESD,CE)
           ELSEIF(ID.eq.4.or.ID.eq.6) THEN
              theta=datan2(Y,(X-RR))
              CALL wf_fieldcr(IE,X,Y,CESD,CE)
              CE_R=CE
              CALL wf_fieldcz(IE,X,Y,CESD,CE)
              CE_Z=CE
              if(ID.eq.4) then
                 CE=-CE_R*cos(theta)-CE_Z*sin(theta)
              elseif(ID.eq.6) then
                 CE= CE_R*sin(theta)-CE_Z*cos(theta)
              end if
           ENDIF
           IF(ABS(CE).LT.1.D-12) CE=(0.D0,0.D0)
           GX(NGV)=gdclip(X)
           GV(NGV,1)=gdclip(REAL(CE))
           GV(NGV,2)=gdclip(AIMAG(CE))
           GV(NGV,3)=gdclip(ABS(CE))
        ENDIF
!        IF(ABS(CE).NE.0.D0) THEN
!           WRITE(6,'(A,I12,1P5E12.4)') 'CESD:',IE,X,YPOS,CE
!           WRITE(6,'(A,1P6E12.4)')     '     ', &
!                CESD(ABS(nseg_nside_nelm(1,IE))),CESD(ABS(nseg_nside_nelm(2,IE))), &
!                CESD(ABS(nseg_nside_nelm(3,IE)))
!           WRITE(6,'(A,1P6E12.4)')     '     ', &
!                xnode(node_nside_nelm(1,IE)),ynode(node_nside_nelm(1,IE)), &
!                xnode(node_nside_nelm(2,IE)),ynode(node_nside_nelm(2,IE)), &
!                xnode(node_nside_nelm(3,IE)),ynode(node_nside_nelm(3,IE))
!        END IF
     ENDDO

  ! --- EXTRACT DATA FOR 1D Y PLOT ---

  else if(KID.eq.'Y') then
     READ(KWD(4:NCHM),*,ERR=9000) XPOS
     DY=(ZNDMAX-ZNDMIN)/(NGVMAX-1)
     DO NGV=1,NGVMAX
        Y=ZNDMIN+DY*(NGV-1)
        IF(Y.GT.0.4D0.AND.Y.LT.0.41D0) THEN
           idebug=-1
        ELSE
           idebug=0
        END IF
        CALL fem_find_nelm_for_xy(xpos,y,ie)
        IF(IE.EQ.0) THEN
           GX(NGV)=gdclip(Y)
           GV(NGV,1)=0.0
           GV(NGV,2)=0.0
           GV(NGV,3)=0.0
        ELSE
           IF(ID.eq.1) THEN
              CALL wf_fieldcr(IE,XPOS,Y,CESD,CE)
           ELSE IF(ID.eq.3) THEN
              CALL wf_fieldcz(IE,XPOS,Y,CESD,CE)
           ELSEIF(ID.eq.4.or.ID.eq.6) THEN
              theta=datan2(Y,(X-RR))
              CALL wf_fieldcr(IE,X,Y,CESD,CE)
              CE_R=CE
              CALL wf_fieldcz(IE,X,Y,CESD,CE)
              CE_Z=CE
              if(ID.eq.4) then
                 CE=-CE_R*cos(theta)-CE_Z*sin(theta)
              elseif(ID.eq.6) then
                 CE= CE_R*sin(theta)-CE_Z*cos(theta)
              end if
           ENDIF
!              IF(idebug.EQ.-1) THEN
!              N1=nseg_nside_nelm(1,IE)
!              IF(N1.GT.0) THEN
!                 CE1=CESD(N1)
!              ELSE
!                 CE1=-CESD(-N1)
!              END IF
!              N2=nseg_nside_nelm(2,IE)
!              IF(N2.GT.0) THEN
!                 CE2=CESD(N2)
!              ELSE
!                 CE2=-CESD(-N2)
!              END IF
!              N3=nseg_nside_nelm(3,IE)
!              IF(N3.GT.0) THEN
!                 CE3=CESD(N3)
!              ELSE
!                 CE3=-CESD(-N3)
!              END IF
             
!              IF(ABS(CE).NE.0.D0) THEN
!                 WRITE(6,'(A,I12,1P5E12.4)') 'CESD:',IE,XPOS,Y,CE
!                 WRITE(6,'(A,1P6E12.4)')     '     ', CE1,CE2,CE3
!                 WRITE(6,'(A,1P6E12.4)')     '     ', &
!                      xnode(node_nside_nelm(1,IE)),ynode(node_nside_nelm(1,IE)), &
!                      xnode(node_nside_nelm(2,IE)),ynode(node_nside_nelm(2,IE)), &
!                      xnode(node_nside_nelm(3,IE)),ynode(node_nside_nelm(3,IE))
!              END IF
!           END IF
           IF(ABS(CE).LT.1.D-12) CE=(0.D0,0.D0)
           GX(NGV)=gdclip(Y)
           GV(NGV,1)=gdclip(REAL(CE))
           GV(NGV,2)=gdclip(AIMAG(CE))
           GV(NGV,3)=gdclip(ABS(CE))
        ENDIF
     ENDDO
  end if
  
9000 CONTINUE
  RETURN
END SUBROUTINE wf_ctog_side

!     ****** EXTRACT FROM COMPLEX DATA ******

SUBROUTINE wf_ctog_node(ID,KWD)

  use wfcomm
  USE wfsub
  USE libgrf
  use feminterpolate
  implicit none
  integer,intent(in) :: ID
  integer :: IE,NGX,NGY,NGV
  real(rkind) :: DX,DY,X,Y
  real(rkind) :: XPOS,YPOS
  complex(rkind) :: CE
  character,intent(in) :: KWD*(NCHM)
  CHARACTER  :: KID*1

  KID=KWD(3:3)
  IE=0

! --- EXTRACT DATA FOR 2D PLOT ---  

  if(KID.eq.'R'.or.&
     KID.eq.'I'.or.&
     KID.eq.'C') then
     CALL wf_gen_elm_mesh
     DO NGY=1,NGYMAX
        Y=G2Y(NGY)
        DO NGX=1,NGXMAX
           X=G2X(NGX)
           !           IE=IEGZ(NGX,NGY)
           CALL fem_find_nelm_for_xy(x,y,ie)
           IF(IE.EQ.0) THEN
              GZ(NGX,NGY)=0.0
           ELSE
              IF(ID.eq.2) THEN
                 CALL wf_fieldcp(IE,X,Y,CEND,CE)
              ! -----------
              ELSEIF(ID.eq.5) THEN
                 CALL wf_fieldcp(IE,X,Y,CEND,CE)
                 CE=-CE
              ENDIF
              IF(ABS(CE).LT.1.D-12) CE=(0.D0,0.D0)
!              IF(ABS(CE).NE.0.D0) THEN
!                 WRITE(6,'(A,I12,1P5E12.4)') 'CEND:',IE,X,Y,CE
!                 WRITE(6,'(A,1P6E12.4)')     '     ', &
!                      CEND(node_nside_nelm(1,IE)),CEND(node_nside_nelm(2,IE)),CEND(node_nside_nelm(3,IE))
!                 WRITE(6,'(A,1P6E12.4)')     '     ', &
!                      xnode(node_nside_nelm(1,IE)),ynode(node_nside_nelm(1,IE)), &
!                      xnode(node_nside_nelm(2,IE)),ynode(node_nside_nelm(2,IE)), &
!                      xnode(node_nside_nelm(3,IE)),ynode(node_nside_nelm(3,IE))
!              END IF
                 
              ! -----------
              IF(KID.EQ.'R') THEN
                 GZ(NGX,NGY)=gdclip(REAL(CE))
              ELSE IF(KID.EQ.'I') THEN
                 GZ(NGX,NGY)=gdclip(AIMAG(CE))
              ELSE IF(KID.EQ.'A') THEN
                 GZ(NGX,NGY)=gdclip(ABS(CE))
              ENDIF
           ENDIF
        END DO
     ENDDO

! --- EXTRACT DATA FOR 1D PLOT ---

  elseif(KID.eq.'X') then
     READ(KWD(4:NCHM),*,ERR=9000) YPOS
     DX=(RNDMAX-RNDMIN)/(NGVMAX-1)
     DO NGV=1,NGVMAX
        X=RNDMIN+DX*(NGV-1)
        CALL fem_find_nelm_for_xy(x,ypos,ie)
        IF(IE.EQ.0) THEN
           GX(NGV)=gdclip(X)
           GV(NGV,1)=0.0
           GV(NGV,2)=0.0
           GV(NGV,3)=0.0
        ELSE
           IF(ID.eq.2) THEN
              CALL wf_fieldcp(IE,X,YPOS,CEND,CE)
           ELSE IF(ID.eq.5) THEN
              CALL wf_fieldcp(IE,X,YPOS,CEND,CE)
              CE=-CE
           ENDIF
           IF(ABS(CE).LT.1.D-12) CE=(0.D0,0.D0)
           GX(NGV)=gdclip(X)
           GV(NGV,1)=gdclip(REAL(CE))
           GV(NGV,2)=gdclip(AIMAG(CE))
           GV(NGV,3)=gdclip(ABS(CE))
        ENDIF
!        IF(ABS(CE).NE.0.D0) THEN
!           WRITE(6,'(A,I12,1P5E12.4)') 'CEND:',IE,X,Y,CE
!           WRITE(6,'(A,1P6E12.4)')     '     ', &
!                CEND(node_nside_nelm(1,IE)),CEND(node_nside_nelm(2,IE)),CEND(node_nside_nelm(3,IE))
!           WRITE(6,'(A,1P6E12.4)')     '     ', &
!                xnode(node_nside_nelm(1,IE)),ynode(node_nside_nelm(1,IE)), &
!                xnode(node_nside_nelm(2,IE)),ynode(node_nside_nelm(2,IE)), &
!                xnode(node_nside_nelm(3,IE)),ynode(node_nside_nelm(3,IE))
!        END IF
     ENDDO
  else if(KID.eq.'Y') then
     READ(KWD(4:NCHM),*,ERR=9000) XPOS
     DY=(ZNDMAX-ZNDMIN)/(NGVMAX-1)
     DO NGV=1,NGVMAX
        Y=ZNDMIN+DY*(NGV-1)
        CALL fem_find_nelm_for_xy(xpos,y,ie)
        IF(IE.EQ.0) THEN
           GX(NGV)=gdclip(Y)
           GV(NGV,1)=0.0
           GV(NGV,2)=0.0
           GV(NGV,3)=0.0
        ELSE
           IF(ID.eq.2) THEN
              CALL wf_fieldcp(IE,XPOS,Y,CEND,CE)
           ELSE IF(ID.eq.5) THEN
              CALL wf_fieldcz(IE,XPOS,Y,CEND,CE)
              CE=-CE
           ENDIF
           IF(ABS(CE).LT.1.D-12) CE=(0.D0,0.D0)
           GX(NGV)=gdclip(Y)
           GV(NGV,1)=gdclip(REAL(CE))
           GV(NGV,2)=gdclip(AIMAG(CE))
           GV(NGV,3)=gdclip(ABS(CE))
        ENDIF
!        IF(ABS(CE).NE.0.D0) THEN
!           WRITE(6,'(A,I12,1P5E12.4)') 'CEND:',IE,X,Y,CE
!           WRITE(6,'(A,1P6E12.4)')     '     ', &
!                CEND(node_nside_nelm(1,IE)),CEND(node_nside_nelm(2,IE)),CEND(node_nside_nelm(3,IE))
!           WRITE(6,'(A,1P6E12.4)')     '     ', &
!                xnode(node_nside_nelm(1,IE)),ynode(node_nside_nelm(1,IE)), &
!                xnode(node_nside_nelm(2,IE)),ynode(node_nside_nelm(2,IE)), &
!                xnode(node_nside_nelm(3,IE)),ynode(node_nside_nelm(3,IE))
!        END IF
     ENDDO
  end if
  
9000 CONTINUE
  RETURN
END SUBROUTINE wf_ctog_node

!     ****** CALCULATE RANGE OF WINDOW ******

SUBROUTINE wf_gwin_range(NW,NWMAX,PXMIN,PXMAX,PYMIN,PYMAX)

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
END SUBROUTINE wf_gwin_range


!     ****** DRAW PARAMETER ON GRAPHIC SCREEN ******

SUBROUTINE wf_gdraw_parm
  
  use wfcomm
  implicit none
  integer :: nant,NS
  real(rkind) :: REST(NAM),REAT(NAM),WW,RNZ
  real :: GXMIN,GYMAX,GRCHH,GDX,GDY,GXL,GYL
  
  DO nant=1,nant_max
     REST(nant)=DBLE(CIMP(nant))
     REAT(nant)=AIMAG(CIMP(nant))
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
  CALL TEXT('node_max=',6)
  CALL NUMBI(node_max,'(I8)',8)
  
  GXL=GXL+GDX
  CALL MOVE(GXL,GYL)
  CALL TEXT('R,Z MAX=',8)
  CALL NUMBD(RNDMAX,'(F7.3)',7)
  CALL NUMBD(ZNDMAX,'(F7.3)',7)
  
  GXL=GXMIN
  GXL=GXL+GRCHH
  GYL=GYL+GDY
  CALL MOVE(GXL,GYL)
  CALL TEXT('nelm_max=',9)
  CALL NUMBI(nelm_max,'(I8)',8)
  
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
  GYL=GYL+GDY
  
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
  IF(nant_max.GT.0) THEN
     CALL MOVE(GXL,GYL)
     CALL TEXT('IJ',2)
     CALL TEXT('  AJ ',5)
     CALL TEXT(' PHASE ',7)
     CALL TEXT('     R     ',11)
     CALL TEXT('     X     ',11)
     
     DO nant=1,nant_max
        GXL=GXMIN+45.*GRCHH
        GYL=GYL+GDY
        CALL MOVE(GXL,GYL)
        CALL NUMBI(nant,'(I2)',2)
        CALL NUMBD(AJ(nant),'(F5.1)',5)
        CALL NUMBD(APH(nant),'(F7.1)',7)
        CALL NUMBD(REST(nant),'(1PE11.3)',11)
        CALL NUMBD(REAT(nant),'(1PE11.3)',11)
     ENDDO
  ENDIF
  
  RETURN
END SUBROUTINE wf_gdraw_parm

!     ****** Draw Vessel wall ******

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

!     ****** Draw Antenna Path ******

SUBROUTINE wf_gdraw_antenna

  use wfcomm
  USE libgrf
  implicit none
  integer :: nant,npoint
  real :: GXL,GYL

  IF(NDRAWA.EQ.0) THEN
     DO nant=1,nant_max
        GXL=gdclip(RJ0(1,nant))
        GYL=gdclip(ZJ0(1,nant))
        CALL MOVE2D(GXL,GYL)
        DO npoint=2,JNUM0(nant)
           GXL=gdclip(RJ0(npoint,nant))
           GYL=gdclip(ZJ0(npoint,nant))
           CALL DRAW2D(GXL,GYL)
        END DO
     END DO
  ELSE
     DO nant=1,nant_max
        GXL=gdclip(RJ(1,nant))
        GYL=gdclip(ZJ(1,nant))
        CALL MOVE2D(GXL,GYL)
        DO npoint=2,JNUM(nant)
           GXL=gdclip(RJ(npoint,nant))
           GYL=gdclip(ZJ(npoint,nant))
           CALL DRAW2D(GXL,GYL)
        END DO
     END DO
  ENDIF
  RETURN
END SUBROUTINE wf_gdraw_antenna

!     ****** Draw Waveguide ******

SUBROUTINE wf_gdraw_waveguide

  use wfcomm
  USE libgrf
  implicit none
  real :: GXL,GYL

  GXL=gdclip(x1wg)
  GYL=gdclip(y1wg)
  CALL MOVE2D(GXL,GYL)
  GXL=gdclip(x2wg)
  CALL DRAW2D(GXL,GYL)
  GYL=gdclip(y2wg)
  CALL DRAW2D(GXL,GYL)
  GXL=gdclip(x1wg)
  CALL DRAW2D(GXL,GYL)
  GYL=gdclip(y1wg)
  CALL DRAW2D(GXL,GYL)
  RETURN
END SUBROUTINE wf_gdraw_waveguide

!     ****** Draw Element Data ******

SUBROUTINE wf_gdraw_element

  use wfcomm
  USE libgrf
  implicit none
  integer :: IE,IN1,IN2,IN3,IEL,IN,INL
  real :: GX1,GX2,GX3,GY1,GY2,GY3,GXC,GYC
  
  CALL SETCHR(0.2,0.15,0.2,0.,-30.)

  DO IE=1,nelm_max
     IN1=node_nside_nelm(1,IE)
     IN2=node_nside_nelm(2,IE)
     IN3=node_nside_nelm(3,IE)
     GX1=gdclip(xnode(IN1))
     GY1=gdclip(ynode(IN1))
     GX2=gdclip(xnode(IN2))
     GY2=gdclip(ynode(IN2))
     GX3=gdclip(xnode(IN3))
     GY3=gdclip(ynode(IN3))
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
        CALL GNUMBI(GXC,GYC,IEL,2)
     ENDIF
  ENDDO
  
  IF(NDRAWD.GE.3) THEN
     CALL SETCHS(0.2,0.)
     DO IN=1,node_max
        GX1=gdclip(xnode(IN))
        GY1=gdclip(ynode(IN))
        INL=IN
        CALL GNUMBI(GX1,GY1,INL,0)
     END DO
  ENDIF
  
  RETURN
END SUBROUTINE wf_gdraw_element

!     ****** Draw Element Data ******

SUBROUTINE wf_gr_element

  use wfcomm
  USE libgrf
  implicit none
  real(rkind) :: YMIN,YMAX,XMIN,XMAX
  real(rkind) :: DXLEN,DYLEN,DRATIO,XMID,XLEN,YMID,YLEN
  
  DXLEN= RNDMAX-RNDMIN
  DYLEN= ZNDMAX-ZNDMIN
  DRATIO=DYLEN/DXLEN

  IF(DRATIO.GT.0.75D0) THEN
     YMIN=ZNDMIN
     YMAX=ZNDMAX
     XMID=0.5D0*(RNDMIN+RNDMAX)
     XLEN=DYLEN/0.75D0
     XMIN=XMID-0.5D0*XLEN
     XMAX=XMID+0.5D0*XLEN
  ELSE
     XMIN=RNDMIN
     XMAX=RNDMAX
     YMID=0.5D0*(ZNDMIN+ZNDMAX)
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
  CALL TEXT('node_max=',6)
  CALL NUMBI(node_max,'(I8)',8)
  
  GYL=GYL-GDY
  CALL MOVE(GXL,GYL)
  CALL TEXT('nelm_max=',9)
  CALL NUMBI(nelm_max,'(I8)',8)
  
  GXL= 8.0
  GYL=17.8
  CALL MOVE(GXL,GYL)
  CALL TEXT('met_len=',8)
  CALL NUMBI(mtx_len,'(I8)',8)
  RETURN
END SUBROUTINE wf_gdraw_parm_elm

!     ****** Draw Antenna ******

SUBROUTINE wf_gr_antenna

  use wfcomm
  USE libgrf
  implicit none
  integer :: NTEMP
  real(rkind) :: YMIN,YMAX,XMIN,XMAX
  real(rkind) :: DXLEN,DYLEN,DRATIO,XMID,XLEN,YMID,YLEN

  DXLEN= RNDMAX-RNDMIN
  DYLEN= ZNDMAX-ZNDMIN
  DRATIO=DYLEN/DXLEN

  IF(DRATIO.GT.0.75D0) THEN
     YMIN=ZNDMIN
     YMAX=ZNDMAX
     XMID=0.5D0*(RNDMIN+RNDMAX)
     XLEN=DYLEN/0.75D0
     XMIN=XMID-0.5D0*XLEN
     XMAX=XMID+0.5D0*XLEN
  ELSE
     XMIN=RNDMIN
     XMAX=RNDMAX
     YMID=0.5D0*(ZNDMIN+ZNDMAX)
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
  
  CALL SETLIN(0,0,5)
  CALL wf_gdraw_plasma
  
  CALL SETLIN(0,0,6)
  CALL wf_gdraw_antenna
  CALL grd2d_frame_end
  
  CALL PAGEE
  RETURN
END SUBROUTINE wf_gr_antenna

!     ****** Draw Waveguide ******

SUBROUTINE wf_gr_waveguide

  use wfcomm
  USE libgrf
  implicit none
  integer :: NTEMP
  real(rkind) :: YMIN,YMAX,XMIN,XMAX
  real(rkind) :: DXLEN,DYLEN,DRATIO,XMID,XLEN,YMID,YLEN

  DXLEN= RNDMAX-RNDMIN
  DYLEN= ZNDMAX-ZNDMIN
  DRATIO=DYLEN/DXLEN

  IF(DRATIO.GT.0.75D0) THEN
     YMIN=ZNDMIN
     YMAX=ZNDMAX
     XMID=0.5D0*(RNDMIN+RNDMAX)
     XLEN=DYLEN/0.75D0
     XMIN=XMID-0.5D0*XLEN
     XMAX=XMID+0.5D0*XLEN
  ELSE
     XMIN=RNDMIN
     XMAX=RNDMAX
     YMID=0.5D0*(ZNDMIN+ZNDMAX)
     YLEN=DXLEN*0.75D0
     YMIN=YMID-0.5D0*YLEN
     YMAX=YMID+0.5D0*YLEN
  ENDIF
  
!  CALL PAGES
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
  
  CALL SETLIN(0,0,5)
  CALL wf_gdraw_plasma
  
  CALL SETLIN(0,0,6)
  CALL wf_gdraw_waveguide
  CALL grd2d_frame_end
  
!  CALL PAGEE
  RETURN
END SUBROUTINE wf_gr_waveguide

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
  CALL TEXT('node_max=',6)
  CALL NUMBI(node_max,'(I8)',8)
  
  GYL=GYL-GDY
  CALL MOVE(GXL,GYL)
  CALL TEXT('nelm_max=',9)
  CALL NUMBI(nelm_max,'(I8)',8)
  
  GXL= 8.0
  GYL=17.8
  CALL MOVE(GXL,GYL)
  CALL TEXT('JNUM =',6)
  GXL=GXL+6*0.3
  DO NA=1,nant_max
     CALL MOVE(GXL,GYL)
     CALL NUMBI(JNUM(NA),'(I5)',5)
     GYL=GYL-GDY
  END DO
  RETURN
END SUBROUTINE wf_gdraw_parm_ant

SUBROUTINE wf_gen_elm_mesh
  USE wfcomm
  USE wfindex
  USE libgrf
  IMPLICIT NONE
  INTEGER,SAVE:: NGXMAX_mesh_SAVE=0,NGYMAX_mesh_SAVE=0
  REAL(rkind),SAVE:: RNDMIN_mesh_SAVE=0.D0,RNDMAX_mesh_SAVE=0.D0
  REAL(rkind),SAVE:: ZNDMIN_mesh_SAVE=0.D0,ZNDMAX_mesh_SAVE=0.D0
  REAL(rkind):: DX,DY,X,Y
  INTEGER:: NGX,NGY,IE

  IF(NGXMAX.NE.NGXMAX_mesh_SAVE.OR. &
     NGYMAX.NE.NGYMAX_mesh_SAVE.OR. &
     RNDMIN.NE.RNDMIN_mesh_SAVE.OR. &
     RNDMAX.NE.RNDMAX_mesh_SAVE.OR. &
     ZNDMIN.NE.ZNDMIN_mesh_SAVE.OR. &
     ZNDMAX.NE.ZNDMAX_mesh_SAVE) THEN
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
           CALL wf_fep(X,Y,IE)
           IEGZ(NGX,NGY)=IE
        END DO
     END DO
     NGXMAX_mesh_SAVE=NGXMAX
     NGYMAX_mesh_SAVE=NGYMAX
     RNDMIN_mesh_SAVE=RNDMIN
     RNDMAX_mesh_SAVE=RNDMAX
     ZNDMIN_mesh_SAVE=ZNDMIN
     ZNDMAX_mesh_SAVE=ZNDMAX
  END IF
  RETURN
END SUBROUTINE wf_gen_elm_mesh

end MODULE wfgsub
