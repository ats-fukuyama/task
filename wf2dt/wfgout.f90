!     *********  /TASKX/WFT/GOUT  *********
!     $Id: wfgout.f90,v 1.24 2012/02/11 01:23:26 maruyama Exp $
!
!       GRAPHIC DATA PROCESSING PROGRAM
!             FOR FEM COMPUTATION
!
!     ***********************************

SUBROUTINE WFGOUT

  use wfcomm
  USE libchar
  implicit none

  integer :: NL,NWD,NCH,NWMAX,NW
  CHARACTER KLINE*80,KWORD*(NCHM),KWD*(NCHM),KID*1,KTAIL*7
  CHARACTER KG1*1,KG2*1
  DIMENSION KWORD(NWDM)

  call wfwin_allocate
  CALL WFSLIM

  if(NDRAWV.eq.1)  call wfgout_2d_vector

  IF(RNDMAX.EQ.RNDMIN) THEN
     WRITE(6,*) 'XX NO DATA IS LOADED FOR GRAPHICS'
     GOTO 9000
  ENDIF
  
1 WRITE(6,*) '# INPUT : E/B/D,X/Y/Z/+/-/P/F,[R/I/A][Xx/Yy/Zz/B6]/'
  WRITE(6,*) '          Xyz/Yxz/Zxy/A9 0-9 V,0-9, L'
  WRITE(6,*) '          P/F,1/2,CXx/CYy/CZz/Xyz/Yxz/Zxy  X=EXIT'
  CALL GUFLSH
  READ(5,'(A80)',ERR=1,END=9000) KLINE
  NWXMAX=0
  
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
  
20 NWMAX=NWD
21 KWD=KWORD(1)
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
           if(NDRAWE.eq.0) CALL WFCTOGSD(1,KWD)  ! E_R    (R,Z,psi)
           if(NDRAWE.eq.1) CALL WFCTOGSD(4,KWD)  ! E_rho  (rho,theta,-psi)
        ELSEIF(KID.EQ.'P') THEN
           if(NDRAWE.eq.0) CALL WFCTOGND(2,KWD)  ! E_psi
           if(NDRAWE.eq.1) CALL WFCTOGND(5,KWD)  ! -E_psi
        ELSEIF(KID.EQ.'Z') THEN
           CALL WFCTOGSD(3,KWD)                  ! E_Z
        ELSEIF(KID.eq.'T') THEN
           CALL WFCTOGSD(6,KWD)                  ! E_theta
        ELSE
           WRITE(6,*) 'XX UNKNOWN KID2:',KID
           GOTO 1000
        ENDIF

     ELSE IF(KID.eq.'P') THEN
        KID=KWD(2:2)
        IF(KID.eq.'1') THEN
           call PWRPLOT(1)
        ELSE IF(KID.eq.'2') THEN
           call PWRPLOT(2)
        ELSE IF(KID.eq.'3') THEN
           call PWRPLOT(3)
        ELSE
           WRITE(6,*) 'XX UNKNOWN KID2:',KID
           GOTO 1000
        ENDIF
     ELSE IF(KID.EQ.'N') THEN
        CALL NPLOT
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
           CALL WFGBFC(KWD)
        CASE(0)
           CALL WFGWFC(KWD)
        CASE(1)
           CALL WFGPPC(NW,NWMAX,KWD)
        CASE(2)
           CALL WFGPFC(NW,NWMAX,KWD)
        CASE(3:6)
           CALL WFGPBC(NW,NWMAX,KWD)
        END SELECT
     ELSEIF(KID.EQ.'X'.OR.&
          & KID.EQ.'Y') THEN
        SELECT CASE(NGRAPH)
        CASE(-1)
           CALL WFGBFR(KWD)
        CASE(0)
           CALL WFGWFR(KWD)
        CASE(1:6)
           CALL WFGPFR(NW,NWMAX,KWD)
        END SELECT
     ELSE
        WRITE(6,*) 'XX UNKNOWN KID3:',KID
     ENDIF
  END DO

1000 continue
  IF(NGRAPH.GE.1) THEN
     CALL WFGPRM
     CALL PAGEE
  ENDIF
  GOTO 1     
  
9000 RETURN
END SUBROUTINE WFGOUT

!     ****** EXTRACT FROM COMPLEX DATA ******

SUBROUTINE WFCTOGSD(ID,KWD)

  use wfcomm
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
  INTEGER:: N1,N2,N3
  COMPLEX(rkind) :: CE1,CE2,CE3

  real(rkind) :: theta
  complex(rkind)::CE_R,CE_Z

  KID=KWD(3:3)
  IE=0

! --- EXTRACT DATA FOR 2D PLOT ---  

  if(KID.eq.'R'.or.&
     KID.eq.'I'.or.&
     KID.eq.'C') then
     CALL WFGMESH
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
                 CALL FIELDCR(IE,X,Y,CESD,CE)
              ELSE IF(ID.eq.3) THEN
                 CALL FIELDCZ(IE,X,Y,CESD,CE)
              ELSEIF(ID.eq.4.or.ID.eq.6) THEN
                 theta=datan2(Y,(X-RR))
                 CALL FIELDCR(IE,X,Y,CESD,CE)
                 CE_R=CE
                 CALL FIELDCZ(IE,X,Y,CESD,CE)
                 CE_Z=CE
                 if(ID.eq.4) then
                    CE=-CE_R*cos(theta)-CE_Z*sin(theta)
                 elseif(ID.eq.6) then
                    CE= CE_R*sin(theta)-CE_Z*cos(theta)
                 end if
              ENDIF
              IF(ABS(CE).LT.1.D-12) CE=(0.D0,0.D0)
!              IF(Y.GT.0.3D0.AND.Y.LT.0.5D0) THEN
!              N1=NSDELM(1,IE)
!              IF(N1.GT.0) THEN
!                 CE1=CESD(N1)
!              ELSE
!                 CE1=-CESD(-N1)
!              END IF
!              N2=NSDELM(2,IE)
!              IF(N2.GT.0) THEN
!                 CE2=CESD(N2)
!              ELSE
!                 CE2=-CESD(-N2)
!              END IF
!              N3=NSDELM(3,IE)
!              IF(N3.GT.0) THEN
!                 CE3=CESD(N3)
!              ELSE
!                 CE3=-CESD(-N3)
!              END IF
             
!              IF(ABS(CE).NE.0.D0) THEN
!                 WRITE(6,'(A,I12,1P5E12.4)') 'CESD:',IE,X,Y,CE
!                 WRITE(6,'(A,1P6E12.4)')     '     ', CE1,CE2,CE3
!                 WRITE(6,'(A,1P6E12.4)')     '     ', &
!                      RNODE(NDELM(1,IE)),ZNODE(NDELM(1,IE)), &
!                      RNODE(NDELM(2,IE)),ZNODE(NDELM(2,IE)), &
!                      RNODE(NDELM(3,IE)),ZNODE(NDELM(3,IE))
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
!        CALL FEP(X,YPOS,IE)
        IF(IE.EQ.0) THEN
           GX(NGV)=gdclip(X)
           GV(NGV,1)=0.0
           GV(NGV,2)=0.0
           GV(NGV,3)=0.0
        ELSE
           IF(ID.eq.1) THEN
              CALL FIELDCR(IE,X,YPOS,CESD,CE)
           ELSE IF(ID.eq.3) THEN
              CALL FIELDCZ(IE,X,YPOS,CESD,CE)
           ELSEIF(ID.eq.4.or.ID.eq.6) THEN
              theta=datan2(Y,(X-RR))
              CALL FIELDCR(IE,X,Y,CESD,CE)
              CE_R=CE
              CALL FIELDCZ(IE,X,Y,CESD,CE)
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
!                CESD(ABS(NSDELM(1,IE))),CESD(ABS(NSDELM(2,IE))), &
!                CESD(ABS(NSDELM(3,IE)))
!           WRITE(6,'(A,1P6E12.4)')     '     ', &
!                RNODE(NDELM(1,IE)),ZNODE(NDELM(1,IE)), &
!                RNODE(NDELM(2,IE)),ZNODE(NDELM(2,IE)), &
!                RNODE(NDELM(3,IE)),ZNODE(NDELM(3,IE))
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
!        CALL FEP(XPOS,Y,IE)
        IF(IE.EQ.0) THEN
           GX(NGV)=gdclip(Y)
           GV(NGV,1)=0.0
           GV(NGV,2)=0.0
           GV(NGV,3)=0.0
        ELSE
           IF(ID.eq.1) THEN
              CALL FIELDCR(IE,XPOS,Y,CESD,CE)
           ELSE IF(ID.eq.3) THEN
              CALL FIELDCZ(IE,XPOS,Y,CESD,CE)
           ELSEIF(ID.eq.4.or.ID.eq.6) THEN
              theta=datan2(Y,(X-RR))
              CALL FIELDCR(IE,X,Y,CESD,CE)
              CE_R=CE
              CALL FIELDCZ(IE,X,Y,CESD,CE)
              CE_Z=CE
              if(ID.eq.4) then
                 CE=-CE_R*cos(theta)-CE_Z*sin(theta)
              elseif(ID.eq.6) then
                 CE= CE_R*sin(theta)-CE_Z*cos(theta)
              end if
           ENDIF
!              IF(idebug.EQ.-1) THEN
!              N1=NSDELM(1,IE)
!              IF(N1.GT.0) THEN
!                 CE1=CESD(N1)
!              ELSE
!                 CE1=-CESD(-N1)
!              END IF
!              N2=NSDELM(2,IE)
!              IF(N2.GT.0) THEN
!                 CE2=CESD(N2)
!              ELSE
!                 CE2=-CESD(-N2)
!              END IF
!              N3=NSDELM(3,IE)
!              IF(N3.GT.0) THEN
!                 CE3=CESD(N3)
!              ELSE
!                 CE3=-CESD(-N3)
!              END IF
             
!              IF(ABS(CE).NE.0.D0) THEN
!                 WRITE(6,'(A,I12,1P5E12.4)') 'CESD:',IE,XPOS,Y,CE
!                 WRITE(6,'(A,1P6E12.4)')     '     ', CE1,CE2,CE3
!                 WRITE(6,'(A,1P6E12.4)')     '     ', &
!                      RNODE(NDELM(1,IE)),ZNODE(NDELM(1,IE)), &
!                      RNODE(NDELM(2,IE)),ZNODE(NDELM(2,IE)), &
!                      RNODE(NDELM(3,IE)),ZNODE(NDELM(3,IE))
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
END SUBROUTINE WFCTOGSD

!     ****** EXTRACT FROM COMPLEX DATA ******

SUBROUTINE WFCTOGND(ID,KWD)

  use wfcomm
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

  real(rkind) :: theta
  complex(rkind)::CE_R,CE_Z

  KID=KWD(3:3)
  IE=0

! --- EXTRACT DATA FOR 2D PLOT ---  

  if(KID.eq.'R'.or.&
     KID.eq.'I'.or.&
     KID.eq.'C') then
     CALL WFGMESH
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
                 CALL FIELDCP(IE,X,Y,CEND,CE)
              ! -----------
              ELSEIF(ID.eq.5) THEN
                 CALL FIELDCP(IE,X,Y,CEND,CE)
                 CE=-CE
              ENDIF
              IF(ABS(CE).LT.1.D-12) CE=(0.D0,0.D0)
!              IF(ABS(CE).NE.0.D0) THEN
!                 WRITE(6,'(A,I12,1P5E12.4)') 'CEND:',IE,X,Y,CE
!                 WRITE(6,'(A,1P6E12.4)')     '     ', &
!                      CEND(NDELM(1,IE)),CEND(NDELM(2,IE)),CEND(NDELM(3,IE))
!                 WRITE(6,'(A,1P6E12.4)')     '     ', &
!                      RNODE(NDELM(1,IE)),ZNODE(NDELM(1,IE)), &
!                      RNODE(NDELM(2,IE)),ZNODE(NDELM(2,IE)), &
!                      RNODE(NDELM(3,IE)),ZNODE(NDELM(3,IE))
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
!        CALL FEP(X,YPOS,IE)
        IF(IE.EQ.0) THEN
           GX(NGV)=gdclip(X)
           GV(NGV,1)=0.0
           GV(NGV,2)=0.0
           GV(NGV,3)=0.0
        ELSE
           IF(ID.eq.2) THEN
              CALL FIELDCP(IE,X,YPOS,CEND,CE)
           ELSE IF(ID.eq.5) THEN
              CALL FIELDCP(IE,X,YPOS,CEND,CE)
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
!                CEND(NDELM(1,IE)),CEND(NDELM(2,IE)),CEND(NDELM(3,IE))
!           WRITE(6,'(A,1P6E12.4)')     '     ', &
!                RNODE(NDELM(1,IE)),ZNODE(NDELM(1,IE)), &
!                RNODE(NDELM(2,IE)),ZNODE(NDELM(2,IE)), &
!                RNODE(NDELM(3,IE)),ZNODE(NDELM(3,IE))
!        END IF
     ENDDO
  else if(KID.eq.'Y') then
     READ(KWD(4:NCHM),*,ERR=9000) XPOS
     DY=(ZNDMAX-ZNDMIN)/(NGVMAX-1)
     DO NGV=1,NGVMAX
        Y=ZNDMIN+DY*(NGV-1)
        CALL fem_find_nelm_for_xy(xpos,y,ie)
!        CALL FEP(XPOS,Y,IE)
        IF(IE.EQ.0) THEN
           GX(NGV)=gdclip(Y)
           GV(NGV,1)=0.0
           GV(NGV,2)=0.0
           GV(NGV,3)=0.0
        ELSE
           IF(ID.eq.2) THEN
              CALL FIELDCP(IE,XPOS,Y,CEND,CE)
           ELSE IF(ID.eq.5) THEN
              CALL FIELDCZ(IE,XPOS,Y,CEND,CE)
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
!                CEND(NDELM(1,IE)),CEND(NDELM(2,IE)),CEND(NDELM(3,IE))
!           WRITE(6,'(A,1P6E12.4)')     '     ', &
!                RNODE(NDELM(1,IE)),ZNODE(NDELM(1,IE)), &
!                RNODE(NDELM(2,IE)),ZNODE(NDELM(2,IE)), &
!                RNODE(NDELM(3,IE)),ZNODE(NDELM(3,IE))
!        END IF
     ENDDO
  end if
  
9000 CONTINUE
  RETURN
END SUBROUTINE WFCTOGND

!     ****** WRITE 2D PROFILE IN BINARY FILE ******

SUBROUTINE WFGBFC(KWD)

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
END SUBROUTINE WFGBFC

!     ****** WRITE 2D PROFILE IN TEXT FILE ******

SUBROUTINE WFGWFC(KWD)

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
END SUBROUTINE WFGWFC

!     ****** WRITE 1D PROFILE IN BINARY FILE ******

SUBROUTINE WFGBFR(KWD)

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
END SUBROUTINE WFGBFR

!     ****** WRITE 1D PROFILE IN TEXT FILE ******

SUBROUTINE WFGWFR(KWD)

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
END SUBROUTINE WFGWFR

! ---- output text data (2D vector field) ----
! This subroutine is created in order to see waveguide eigenmode as vector field.
! Output file includes only Er and Ez, not E_phi.

subroutine wfgout_2d_vector

  use wfcomm
  USE libgrf
  implicit none

  integer :: IE,NGX,NGY
  complex(rkind):: CE
  real(rkind):: DX,DY,X,Y
  real(rkind),dimension(:,:),ALLOCATABLE::GZ_r,GZ_z

  allocate(GZ_r(NGXMAX,NGYMAX),GZ_z(NGXMAX,NGYMAX))
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
        IF(IE.EQ.0) THEN
           GZ_r(NGX,NGY)=0.0
           GZ_z(NGX,NGY)=0.0
        ELSE
           CALL FIELDCR(IE,X,Y,CESD,CE)
           GZ_r(NGX,NGY)=gdclip(AIMAG(CE))
           CALL FIELDCZ(IE,X,Y,CESD,CE)
           GZ_z(NGX,NGY)=gdclip(AIMAG(CE))
       ENDIF
     END DO
  ENDDO
  
  open(100,file="vfield")
  write(100,'(7X,A2,15X,A2,15X,A2,15X,A2)') " r"," z","Er","Ez"
  DO NGY=1,NGYMAX
     Y=ZNDMIN+DY*(NGY-1)
     DO NGX=1,NGXMAX
        X=RNDMIN+DX*(NGX-1)
        write(100,'(E16.8,1X,E16.8,1X,E16.8,1X,E16.8)') &
             X,Y,GZ_r(NGX,NGY),GZ_z(NGX,NGY)
     end DO
  end DO
  close(100)

  deallocate(GZ_r,GZ_z)
  return
end subroutine wfgout_2d_vector
