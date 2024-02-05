!     $Id: wfindx.f90,v 1.9 2011/11/06 07:51:15 fukuyama Exp $

!     ***** SORT INDEX FUNCTION ******

FUNCTION FINDEX(R,Z)

  use wfcomm
  implicit none
  real(rkind) :: FINDEX,RN,R,ZN,Z

  RN=(R-RNDMIN)/(RNDMAX-RNDMIN)
  ZN=(Z-ZNDMIN)/(ZNDMAX-ZNDMIN)
  FINDEX=0.999D3*RN+ZN

  RETURN
END FUNCTION FINDEX

!     ***** SORT ELEMENTS BY SINDEX *****

SUBROUTINE WFINDX
  
  use wfcomm
  implicit none
  integer :: NE,IN,NN
  real(rkind) :: RC,ZC,FINDEX
  EXTERNAL WFSRTS,WFSRTX

  CALL WFSLIM
  call wfsrt_allocate
  DO NE=1,nelm_max
     RC=0.D0
     ZC=0.D0
     DO IN=1,3
        NN=NDELM(IN,NE)
        RC=RC+RNODE(NN)
        ZC=ZC+ZNODE(NN)
     ENDDO
     SINDEX(NE)=FINDEX(0.333D0*RC,0.333D0*ZC)
  ENDDO
  
  DO NE=1,nelm_max
     IVELM(NE)=NE
     IDELM(NE)=0
     KAELM(NE)=0
  ENDDO

!  DO NE=1,nelm_max
!     WRITE(6,'(A,2I8,1P,E12.4)') & 
!          &        'NE,IV,SINDEX=',NE,IVELM(NE),SINDEX(NE)
!  ENDDO

  CALL WFSORT(nelm_max,SINDEX,WFSRTS,WFSRTX)

  DO NE=1,nelm_max
     IWELM(NE)=0
  ENDDO
  DO NE=1,nelm_max
     IWELM(IVELM(NE))=NE
  ENDDO

!  DO NE=1,nelm_max
!     WRITE(6,'(A,3I8,1P,E12.4)') &
!          &        'NE,IV,IW,SINDEX=',NE,IVELM(NE),IWELM(NE),SINDEX(NE)
!  ENDDO

  DO NE=1,nelm_max
     if (nrank.eq.0) then
        IF(IWELM(NE).EQ.0) WRITE(6,*) 'XXXX IWELM UNDEFINED'
     end if
  ENDDO
  
  CALL WFSLIM
  CALL WFSELM

  RETURN
END SUBROUTINE WFINDX

!     ***** SET FEP DATA *****

SUBROUTINE WFFEPI

  use wfcomm
  implicit none
  integer :: NE,IN,NN
  real(rkind) :: SINDXL(3),RNDL(3),ZNDL(3),FINDEX
  real(rkind) :: SMAX,SMIN

  DO NE=1,nelm_max
     DO IN=1,3
        NN=NDELM(IN,NE)
        SINDXL(IN)=FINDEX(RNODE(NN),ZNODE(NN))
        RNDL(IN)=RNODE(NN)
        ZNDL(IN)=ZNODE(NN)
     ENDDO
     SINDEX_MIN(NE)=MIN(SINDXL(1),SINDXL(2),SINDXL(3))
     SINDEX_MAX(NE)=MAX(SINDXL(1),SINDXL(2),SINDXL(3))
     REMIN(NE)=MIN(RNDL(1),RNDL(2),RNDL(3))
     REMAX(NE)=MAX(RNDL(1),RNDL(2),RNDL(3))
     ZEMIN(NE)=MIN(ZNDL(1),ZNDL(2),ZNDL(3))
     ZEMAX(NE)=MAX(ZNDL(1),ZNDL(2),ZNDL(3))
  ENDDO
  
  SMAX=SINDEX_MAX(1)
  DO NE=1,nelm_max
     IF(SINDEX_MAX(NE).GT.SMAX) SMAX=SINDEX_MAX(NE)
     SINDEX_MAX(NE)=SMAX
  ENDDO
  SMIN=SINDEX_MIN(nelm_max)
  DO NE=nelm_max,1,-1
     IF(SINDEX_MIN(NE).LT.SMIN) SMIN=SINDEX_MIN(NE)
     SINDEX_MIN(NE)=SMIN
  ENDDO
  
  RETURN
END SUBROUTINE WFFEPI

!     ******* FIND ELEMENT INCLUDING NODES N1,N2 *******

SUBROUTINE EFINDS(IES,N1,N2,IE)

  use wfcomm
  implicit none
  integer,intent(in):: IES,N1,N2
  integer,intent(out):: IE
  integer :: I,J,K,ND1,L,ND2

  IF(IES.LT.0.OR.IES.GT.nelm_max) GOTO 9000
  DO I=1,MAX(nelm_max-IES,IES)
     DO J=1,2
        IF(J.EQ.1) THEN
           IE=IES+I
        ELSE
           IE=IES-I
        ENDIF
        IF(IE.GE.1.AND.IE.LE.nelm_max) THEN
           DO K=1,3
              ND1=NDELM(K,IE)
              IF(ND1.EQ.N1) THEN
                 DO L=1,3
                    ND2=NDELM(L,IE)
                    IF(ND2.EQ.N2) RETURN
                 ENDDO
              ENDIF
           ENDDO
        ENDIF
     ENDDO
  ENDDO
  
  IE=0
  RETURN
  
9000 IE=0
  RETURN
END SUBROUTINE EFINDS

!     ******* FIND ELEMENT INCLUDING POINT (R,Z) *******

SUBROUTINE FEP(R,Z,IE)

  use wfcomm
  implicit none
  real(rkind),intent(in) :: R,Z
  integer,intent(inout):: IE
  integer :: ICOUNT,INMIN,NELMAX,NELMIN,IES,I,IDELT,J
  real(rkind),parameter :: EPS = 1.d-12
  real(rkind) :: WGT(3),WGTMIN,SIDX,FINDEX
  INTEGER:: IN,NN
  
  IF(IE.NE.0) THEN   ! with initial guess

!    move to a neighboring element colser to the point

     ICOUNT=0
     
1000 CONTINUE
     ICOUNT=ICOUNT+1
     IF(ICOUNT.GT.nelm_max/3)  GOTO 2000 ! Could not find, try another method

     CALL WFWGT(IE,R,Z,WGT)

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

     IF(IDEBUG.EQ.2) THEN
        if (nrank.eq.0) THEN
           WRITE(6,'(A,I8,1P5E12.4)') 'IE,WGT,R,Z:', &
                        IE,WGT(1),WGT(2),WGT(3),R,Z
           DO IN=1,3
              NN=NDELM(IN,IE)
              WRITE(6,'(A,3I8,1P2E12.4)') 'IE,IN,NN,R,Z,:', &
                        IE,IN,NN,RNODE(NN),ZNODE(NN)
           END DO
        END if
     ENDIF
     
     IF(WGTMIN.GE.-EPS) RETURN
     
     IE=KNELM(INMIN,IE)
     IF(IE.LE.0) GOTO 2000 ! out of region
     GOTO 1000
  ENDIF

! without initial guess or could not find by sequential scheme
! find the element using the index of elements

2000 SIDX=FINDEX(R,Z)
  IF(SIDX.GE.SINDEX_MIN(1)) THEN
     IF(SIDX.GT.SINDEX_MIN(nelm_max)) THEN
        NELMAX=nelm_max
     ELSE
        CALL WFLCAT_MIN(SINDEX_MIN,nelm_max,SIDX,NELMAX)
     ENDIF
     IF(SIDX.LE.SINDEX_MAX(nelm_max)) THEN
        IF(SIDX.LT.SINDEX_MAX(1)) THEN
           NELMIN=1
        ELSE
           CALL WFLCAT_MIN(SINDEX_MAX,nelm_max,SIDX,NELMIN)
        ENDIF
        
        IES=(NELMIN+NELMAX)/2
        DO I=0,MAX(NELMAX-IES,IES-NELMIN)+1
           IDELT=I
           DO J=1,2
              IDELT=-IDELT
              IE=IES+IDELT
              IF(IE.GE.1.AND.IE.LE.nelm_max) THEN
                 IF(R.GE.REMIN(IE).AND.R.LE.REMAX(IE).AND.&
                  & Z.GE.ZEMIN(IE).AND.Z.LE.ZEMAX(IE)) THEN
                    CALL WFWGT(IE,R,Z,WGT)
                    WGTMIN=MIN(WGT(1),WGT(2),WGT(3))
                    if (nrank.eq.0) then
                       IF(IDEBUG.EQ.2) THEN
                          WRITE(6,'(A,I8,1P3E12.4)') &
                                     'IE,R,Z,WGTMIN=', &
                                      IE,R,Z,WGTMIN
                       end if
                    END if
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
