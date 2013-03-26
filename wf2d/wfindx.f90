!     $Id$

!     ***** SORT INDEX FUNCTION ******

FUNCTION FINDEX(X,Y,Z)

  use wfcomm
  implicit none
  real(8) :: FINDEX,XN,X,YN,Y,ZN,Z

  XN=(X-XNDMIN)/(XNDMAX-XNDMIN)
  YN=(Y-YNDMIN)/(YNDMAX-YNDMIN)
  ZN=(Z-ZNDMIN)/(ZNDMAX-ZNDMIN)
  IF(MODELS.EQ.0) THEN
     FINDEX=1.001D0*XN+0.999D0*YN+ZN
  ELSEIF(MODELS.EQ.1) THEN
     FINDEX=1.D3*XN+1.001D0*YN+0.999D0*ZN
  ELSEIF(MODELS.EQ.2) THEN
     FINDEX=0.999D0*XN+1.D3*YN+1.001D0*ZN
  ELSEIF(MODELS.EQ.3) THEN
     FINDEX=1.001D0*XN+0.999D0*YN+1.D3*ZN
  ELSEIF(MODELS.EQ.4) THEN
     FINDEX=XN+1.001D3*YN+0.999D3*ZN
  ELSEIF(MODELS.EQ.5) THEN
     FINDEX=0.999D3*XN+YN+1.001D3*ZN
  ELSEIF(MODELS.EQ.6) THEN
     FINDEX=1.001D3*XN+0.999D3*YN+ZN
  ENDIF
  RETURN
END FUNCTION FINDEX

!     ***** SORT ELEMENTS BY SINDEX *****

SUBROUTINE WFINDX
  
  use wfcomm
  implicit none
  integer :: NE,IN,NN
  real(8) :: XC,YC,ZC,FINDEX
  EXTERNAL WFSRTS,WFSRTX

  CALL WFVLIM
  call wfsrt_allocate
  
  DO NE=1,NEMAX
     XC=0.D0
     YC=0.D0
     ZC=0.D0
     DO IN=1,4
        NN=NDELM(IN,NE)
        XC=XC+XND(NN)
        YC=YC+YND(NN)
        ZC=ZC+ZND(NN)
     ENDDO
     SINDEX(NE)=FINDEX(0.25D0*XC,0.25D0*YC,0.25D0*ZC)
  ENDDO
  
  DO NE=1,NEMAX
     IVELM(NE)=NE
  ENDDO

!  DO NE=1,NEMAX
!     WRITE(6,'(A,2I8,1P,E12.4)') & 
!          &        'NE,IV,SINDEX=',NE,IVELM(NE),SINDEX(NE)
!  ENDDO

  CALL WFSORT(NEMAX,SINDEX,WFSRTS,WFSRTX)

  DO NE=1,NEMAX
     IWELM(NE)=0
  ENDDO
  DO NE=1,NEMAX
     IWELM(IVELM(NE))=NE
  ENDDO
!  DO NE=1,NEMAX
!     WRITE(6,'(A,3I8,1P,E12.4)') &
!          &        'NE,IV,IW,SINDEX=',NE,IVELM(NE),IWELM(NE),SINDEX(NE)
!  ENDDO
  DO NE=1,NEMAX
     if (nrank.eq.0) then
        IF(IWELM(NE).EQ.0) WRITE(6,*) 'XXXX IWELM UNDEFINED'
     end if
  ENDDO
  
  CALL WFVLIM
  CALL WFVELM
  
  RETURN
END SUBROUTINE WFINDX

!     ***** SET FEP DATA *****

SUBROUTINE WFFEPI

  use wfcomm
  implicit none
  integer :: NE,IN,NN
  real(8) :: SINDXL(4),XNDL(4),YNDL(4),ZNDL(4),FINDEX
  real(8) :: SMAX,SMIN

  DO NE=1,NEMAX
     DO IN=1,4
        NN=NDELM(IN,NE)
        SINDXL(IN)=FINDEX(XND(NN),YND(NN),ZND(NN))
        XNDL(IN)=XND(NN)
        YNDL(IN)=YND(NN)
        ZNDL(IN)=ZND(NN)
     ENDDO
     SINDEX_MIN(NE)=MIN(SINDXL(1),SINDXL(2),SINDXL(3),SINDXL(4))
     SINDEX_MAX(NE)=MAX(SINDXL(1),SINDXL(2),SINDXL(3),SINDXL(4))
     XEMIN(NE)=MIN(XNDL(1),XNDL(2),XNDL(3),XNDL(4))
     XEMAX(NE)=MAX(XNDL(1),XNDL(2),XNDL(3),XNDL(4))
     YEMIN(NE)=MIN(YNDL(1),YNDL(2),YNDL(3),YNDL(4))
     YEMAX(NE)=MAX(YNDL(1),YNDL(2),YNDL(3),YNDL(4))
     ZEMIN(NE)=MIN(ZNDL(1),ZNDL(2),ZNDL(3),ZNDL(4))
     ZEMAX(NE)=MAX(ZNDL(1),ZNDL(2),ZNDL(3),ZNDL(4))
  ENDDO
  
  SMAX=SINDEX_MAX(1)
  DO NE=1,NEMAX
     IF(SINDEX_MAX(NE).GT.SMAX) SMAX=SINDEX_MAX(NE)
     SINDEX_MAX(NE)=SMAX
  ENDDO
  SMIN=SINDEX_MIN(NEMAX)
  DO NE=NEMAX,1,-1
     IF(SINDEX_MIN(NE).LT.SMIN) SMIN=SINDEX_MIN(NE)
     SINDEX_MIN(NE)=SMIN
  ENDDO
  
  RETURN
END SUBROUTINE WFFEPI

!     ******* FIND ELEMENT INCLUDING NODES N1,N2,N3 *******

SUBROUTINE EFINDK(IES,N1,N2,N3,IE)

  use wfcomm
  implicit none
  integer :: IES,NELMIN,NELMAX,I,J,K,ND1,N1,L,ND2,N2,M,ND3,N3,IE
  real(8) :: SIDX1,FINDEX,SIDX2,SIDX3,SIDX_MIN,SIDX_MAX

  
  IF(IES.LT.0.OR.IES.GT.NEMAX) GOTO 9000
  
  SIDX1=FINDEX(XND(N1),YND(N1),ZND(N1))
  SIDX2=FINDEX(XND(N2),YND(N2),ZND(N2))
  SIDX3=FINDEX(XND(N3),YND(N3),ZND(N3))
  SIDX_MIN=MIN(SIDX1,SIDX2,SIDX3)
  SIDX_MAX=MAX(SIDX1,SIDX2,SIDX3)
  
  IF(SIDX_MAX.GT.SINDEX_MAX(NEMAX)) THEN
     NELMIN=NEMAX
10   IF(SINDEX_MAX(NELMIN-1).EQ.SINDEX_MAX(NEMAX)) THEN
        NELMIN=NELMIN-1
        IF(NELMIN.GT.1) GOTO 10
     ENDIF
  ELSEIF(SIDX_MAX.LT.SINDEX_MAX(1)) THEN
     NELMIN=1
  ELSE
     CALL WFLCAT_MAX(SINDEX_MAX,NEMAX,SIDX_MAX,NELMIN)
  ENDIF
  
  IF(SIDX_MIN.LT.SINDEX_MIN(1)) THEN
     NELMAX=1
20   IF(SINDEX_MIN(NELMAX+1).EQ.SINDEX_MIN(1)) THEN
        NELMAX=NELMAX+1
        IF(NELMAX.LT.NEMAX) GOTO 20
     ENDIF
  ELSEIF(SIDX_MIN.GT.SINDEX_MIN(NEMAX)) THEN
     NELMAX=NEMAX
  ELSE
     CALL WFLCAT_MIN(SINDEX_MIN,NEMAX,SIDX_MIN,NELMAX)
  ENDIF
  
  IF(IES.LT.NELMIN.OR.IES.GT.NELMAX) THEN
     if (nrank.eq.0) WRITE(6,*) 'XX NELMIN,IES,NELMAX=',NELMIN,IES,NELMAX
!         WRITE(6,'(A,1P3E12.4)') &
!     &              'MIN:',SINDEX_MIN(1),SIDX_MIN,SINDEX_MIN(NEMAX)
!         WRITE(6,'(A,1P3E12.4)') &
!     &              'MAX:',SINDEX_MAX(1),SIDX_MAX,SINDEX_MAX(NEMAX)
!         WRITE(6,'(A,I5,1P3E12.4)') 'N1,X,Y,Z=', &
!     &                               N1,XND(N1),YND(N1),ZND(N1)
!         WRITE(6,'(A,I5,1P3E12.4)') 'N2,X,Y,Z=', &
!     &                               N2,XND(N2),YND(N2),ZND(N2)
!         WRITE(6,'(A,I5,1P3E12.4)') 'N3,X,Y,Z=', &
!     &                               N3,XND(N3),YND(N3),ZND(N3)
     IES=(NELMIN+NELMAX)/2
  ENDIF
  
  DO I=1,MAX(NELMAX-IES,IES-NELMIN)
     DO J=1,2
        IF(J.EQ.1) THEN
           IE=IES+I
        ELSE
           IE=IES-I
        ENDIF
        IF(IE.GE.1.AND.IE.LE.NEMAX) THEN
           DO K=1,4
              ND1=NDELM(K,IE)
              IF(ND1.EQ.N1) THEN
                 DO L=1,4
                    ND2=NDELM(L,IE)
                    IF(ND2.EQ.N2) THEN
                       DO M=1,4
                          ND3=NDELM(M,IE)
                          IF(ND3.EQ.N3) RETURN
                       ENDDO
                    ENDIF
                 ENDDO
              ENDIF
           ENDDO
        ENDIF
     ENDDO
  ENDDO
  
  IE=0
  RETURN
  
9000 IE=0
  if (nrank.eq.0) WRITE(6,*) 'XX EFINDK: 9000'
  RETURN
  
! 9001 IE=0
!      WRITE(6,*) 'XX EFINDK: 9001'
!      WRITE(6,'(1P3E12.4)') SINDEX_MAX(1),SIDX_MAX,SINDEX_MAX(NEMAX)
!      RETURN

! 9002 IE=0
!      WRITE(6,*) 'XX EFINDK: 9002'
!      WRITE(6,'(1P3E12.4)') SINDEX_MIN(1),SIDX_MIN,SINDEX_MIN(NEMAX)
!      RETURN
END SUBROUTINE EFINDK

!     ******* FIND SUFACE INCLUDING NODES N1,N2 *******

SUBROUTINE EFINDS(NSF,N1,N2,ISF)

  use wfcomm
  implicit none
  integer :: NSF,I,J,K,ND1,N1,L,ND2,N2,ISF

  IF(NSF.LT.0.OR.NSF.GT.NSFMAX) GOTO 9000
  DO I=1,MAX(NSFMAX-NSF,NSF)
     DO J=1,2
        IF(J.EQ.1) THEN
           ISF=NSF+I
        ELSE
           ISF=NSF-I
        ENDIF
        IF(ISF.GE.1.AND.ISF.LE.NSFMAX) THEN
           DO K=1,3
              ND1=NDSRF(K,ISF)
              IF(ND1.EQ.N1) THEN
                 DO L=1,3
                    ND2=NDSRF(L,ISF)
                    IF(ND2.EQ.N2) RETURN
                 ENDDO
              ENDIF
           ENDDO
        ENDIF
     ENDDO
  ENDDO
  
  ISF=0
  RETURN
  
9000 ISF=0
  RETURN
END SUBROUTINE EFINDS

!     ******* FIND ELEMENT INCLUDING POINT (X,Y) *******

SUBROUTINE FEP(X,Y,Z,IE)

  use wfcomm
  implicit none
  integer :: IE,ICOUNT,INMIN,NN1,NN2,NN3,NELMAX,NELMIN,IES,I,IDELT,J
  real(8) :: EPS,X,Y,Z,WGT(4),WGTMIN,XC,YC,ZC,RES,SIDX,FINDEX
  DATA EPS/1.D-12/
  
  IF(IE.NE.0) THEN
     ICOUNT=0
     
1000 CONTINUE
     ICOUNT=ICOUNT+1
     IF(ICOUNT.GT.NEMAX/4)  GOTO 2000
     CALL WFWGT(X,Y,Z,IE,WGT)
     WGTMIN=WGT(1)
     INMIN=1
     IF(WGT(2).LT.WGTMIN) THEN
        WGTMIN=WGT(2)
        INMIN=2
     ENDIF
     IF(WGT(3).LT.WGTMIN) THEN
        WGTMIN=WGT(3)
        INMIN=3
     ENDIF
     IF(WGT(4).LT.WGTMIN) THEN
        WGTMIN=WGT(4)
        INMIN=4
     ENDIF
     
     IF(IDEBUG.NE.0) THEN
        NN1=NDELM(1,IE)
        NN2=NDELM(2,IE)
        NN3=NDELM(3,IE)
        XC=(XND(NN1)+XND(NN2)+XND(NN3))/3.D0
        YC=(YND(NN1)+YND(NN2)+YND(NN3))/3.D0
        ZC=(ZND(NN1)+ZND(NN2)+ZND(NN3))/3.D0
        RES=SQRT((X-XC)**2+(Y-YC)**2+(Z-ZC)**2)
        if (nrank.eq.0) WRITE(6,'(A,I5,1P5E12.4)') 'IE,WGT:', &
             &              IE,WGT(1),WGT(2),WGT(3),WGT(4),RES
     ENDIF
     
     IF(WGTMIN.GE.-EPS) RETURN
     
     IE=KNELM(INMIN,IE)
     IF(IE.LE.0) GOTO 2000
     GOTO 1000
  ENDIF
  
2000 SIDX=FINDEX(X,Y,Z)
  if (nrank.eq.0) then
     IF(IDEBUG.NE.0) WRITE(6,'(A,1P3E12.4)') 'SIDX:MIN=',&
          &     SIDX,SINDEX_MIN(1),SINDEX_MIN(NEMAX)
     IF(IDEBUG.NE.0) WRITE(6,'(A,1P3E12.4)') 'SIDX:MAX=',&
          &     SIDX,SINDEX_MAX(1),SINDEX_MAX(NEMAX)
  end  if
  IF(SIDX.GE.SINDEX_MIN(1)) THEN
     IF(SIDX.GT.SINDEX_MIN(NEMAX)) THEN
        NELMAX=NEMAX
     ELSE
        CALL WFLCAT_MIN(SINDEX_MIN,NEMAX,SIDX,NELMAX)
     ENDIF
     IF(SIDX.LE.SINDEX_MAX(NEMAX)) THEN
        IF(SIDX.LT.SINDEX_MAX(1)) THEN
           NELMIN=1
        ELSE
           CALL WFLCAT_MIN(SINDEX_MAX,NEMAX,SIDX,NELMIN)
        ENDIF
        
        IES=(NELMIN+NELMAX)/2
        if (nrank.eq.0) then
           IF(IDEBUG.NE.0) WRITE(6,*) 'IE,NELMIN,NELMAX=',IE,NELMIN,NELMAX
        end if
        DO I=0,MAX(NELMAX-IES,IES-NELMIN)+1
           IDELT=I
           DO J=1,2
              IDELT=-IDELT
              IE=IES+IDELT
              IF(IE.GE.1.AND.IE.LE.NEMAX) THEN
!               WRITE(6,'(A,1P3E12.4)') 'X:',XEMIN(IE),X,XEMAX(IE)
!               WRITE(6,'(A,1P3E12.4)') 'Y:',YEMIN(IE),Y,YEMAX(IE)
!               WRITE(6,'(A,1P3E12.4)') 'Z:',ZEMIN(IE),Z,ZEMAX(IE)
                 CALL WFWGT(X,Y,Z,IE,WGT)
                 if (nrank.eq.0) then
                    IF(IDEBUG.NE.0) WRITE(6,'(A,I8,1PE12.4)')&
                         &                 'IE,WGTMIN=',IE,WGTMIN
                 end if
                 IF(X.GE.XEMIN(IE).AND.X.LE.XEMAX(IE).AND.&
                      &            Y.GE.YEMIN(IE).AND.Y.LE.YEMAX(IE).AND.&
                      &            Z.GE.ZEMIN(IE).AND.Z.LE.ZEMAX(IE)) THEN
                    CALL WFWGT(X,Y,Z,IE,WGT)
                    WGTMIN=MIN(WGT(1),WGT(2),WGT(3),WGT(4))
                    if (nrank.eq.0) then
                       IF(IDEBUG.NE.0) WRITE(6,'(A,I8,1PE12.4)')&
                            &                 'IE,WGTMIN=',IE,WGTMIN
                    end if
                    IF(WGTMIN.GE.-EPS) RETURN
                 ENDIF
              ENDIF
           ENDDO
        ENDDO
     ENDIF
  ENDIF
  IE=0
  RETURN
END SUBROUTINE FEP
