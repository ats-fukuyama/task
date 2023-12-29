! wfindex.f90

MODULE wfindex

  PRIVATE
  PUBLIC wf_findex
  PUBLIC wf_index
  PUBLIC wf_find_elm
  PUBLIC wf_set_fep
  PUBLIC wf_fep
  PUBLIC wf_set_bdy

CONTAINS

!     ***** SORT INDEX FUNCTION ******

FUNCTION wf_findex(R,Z)

  use wfcomm
  implicit none
  real(rkind) :: wf_findex,RN,R,ZN,Z

  RN=(R-RNDMIN)/(RNDMAX-RNDMIN)
  ZN=(Z-ZNDMIN)/(ZNDMAX-ZNDMIN)
  wf_findex=0.999D3*RN+ZN

  RETURN
END FUNCTION wf_findex

!     ***** SORT ELEMENTS BY SINDEX *****

SUBROUTINE wf_index
  
  use wfcomm
  USE wfsub
  USE wfsort
  implicit none
  integer :: NE,IN,NN
  real(rkind) :: RC,ZC

  CALL wf_set_node_range
  call wf_sort_allocate
  DO NE=1,nelm_max
     RC=0.D0
     ZC=0.D0
     DO IN=1,3
        NN=node_nside_nelm(IN,NE)
        RC=RC+xnode(NN)
        ZC=ZC+ynode(NN)
     ENDDO
     SINDEX(NE)=wf_findex(0.333D0*RC,0.333D0*ZC)
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

  CALL wf_sort(nelm_max,SINDEX,wf_sort_subst,wf_sort_exchange)

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
  
  CALL wf_set_node_range
  CALL wf_set_elm_area

  RETURN
END SUBROUTINE wf_index

!     ***** SET FEP DATA *****

SUBROUTINE wf_set_fep

  use wfcomm
  implicit none
  integer :: NE,IN,NN
  real(rkind) :: SINDXL(3),RNDL(3),ZNDL(3)
  real(rkind) :: SMAX,SMIN

  DO NE=1,nelm_max
     DO IN=1,3
        NN=node_nside_nelm(IN,NE)
        SINDXL(IN)=wf_findex(xnode(NN),ynode(NN))
        RNDL(IN)=xnode(NN)
        ZNDL(IN)=ynode(NN)
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
END SUBROUTINE wf_set_fep

!     ******* FIND ELEMENT INCLUDING NODES N1,N2 *******

SUBROUTINE wf_find_elm(IES,N1,N2,IE)

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
              ND1=node_nside_nelm(K,IE)
              IF(ND1.EQ.N1) THEN
                 DO L=1,3
                    ND2=node_nside_nelm(L,IE)
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
END SUBROUTINE wf_find_elm

!     ******* FIND ELEMENT INCLUDING POINT (R,Z) *******

SUBROUTINE wf_fep(R,Z,IE)

  use wfcomm
  USE wfsort
  USE wfsub
  implicit none
  real(rkind),intent(in) :: R,Z
  integer,intent(inout):: IE
  integer :: ICOUNT,INMIN,NELMAX,NELMIN,IES,I,IDELT,J
  real(rkind),parameter :: EPS = 1.d-12
  real(rkind) :: WGT(3),WGTMIN,SIDX
  INTEGER:: IN,NN
  
  IF(IE.NE.0) THEN   ! with initial guess

!    move to a neighboring element colser to the point

     ICOUNT=0
     
1000 CONTINUE
     ICOUNT=ICOUNT+1
     IF(ICOUNT.GT.nelm_max/3)  GOTO 2000 ! Could not find, try another method

     CALL wf_set_weight(IE,R,Z,WGT)

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

     IF(IDEBUG.EQ.4) THEN
        if (nrank.eq.0) THEN
           WRITE(6,'(A,I8,1P5E12.4)') 'IE,WGT,R,Z:', &
                        IE,WGT(1),WGT(2),WGT(3),R,Z
           DO IN=1,3
              NN=node_nside_nelm(IN,IE)
              WRITE(6,'(A,3I8,1P2E12.4)') 'IE,IN,NN,R,Z,:', &
                        IE,IN,NN,xnode(NN),ynode(NN)
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

2000 CONTINUE
  SIDX=wf_findex(R,Z)
  IF(SIDX.GE.SINDEX_MIN(1)) THEN
     IF(SIDX.GT.SINDEX_MIN(nelm_max)) THEN
        NELMAX=nelm_max
     ELSE
        CALL wf_eval_index_min(SINDEX_MIN,nelm_max,SIDX,NELMAX)
     ENDIF
     IF(SIDX.LE.SINDEX_MAX(nelm_max)) THEN
        IF(SIDX.LT.SINDEX_MAX(1)) THEN
           NELMIN=1
        ELSE
           CALL wf_eval_index_min(SINDEX_MAX,nelm_max,SIDX,NELMIN)
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
                    CALL wf_set_weight(IE,R,Z,WGT)
                    WGTMIN=MIN(WGT(1),WGT(2),WGT(3))
                    if (nrank.eq.0) then
                       IF(IDEBUG.EQ.4) THEN
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
END SUBROUTINE wf_fep

!     ****** Set Boundary Attribute for Side and Node ******

SUBROUTINE wf_set_bdy(IERR)

  use wfcomm
  implicit none
  integer,intent(out) :: IERR
  integer :: ISD,NSD,NE
  integer :: NN1,NN2,NN3,NSD1,NSD2,NSD3
  integer :: NEL,NSDL

! ----- SET nseg_max & KNELM -----
! nseg_bdy_max :: Number of Boundary Side
! KNELM(ISD,NE) :: THE ELEMENT WHICH SHARE THE SIDE (ISD,NE)
!                  IF KNELM.EQ.0, SIDE IS BOUNDARY

  if(nrank.eq.0) WRITE(6,*) '------- wf_set_bdy set nseg_max & KNELM start ---'

!  WRITE(6,'(A,2I8)') 'nelm_max,node_max=',nelm_max,node_max
  nseg_bdy_max=0
  DO NE=1,nelm_max
     NN1=node_nside_nelm(1,NE)
     NN2=node_nside_nelm(2,NE)
     NN3=node_nside_nelm(3,NE)
     CALL wf_find_elm(NE,NN1,NN2,KNELM(1,NE))
     CALL wf_find_elm(NE,NN2,NN3,KNELM(2,NE))
     CALL wf_find_elm(NE,NN3,NN1,KNELM(3,NE))
     IF(KNELM(1,NE).EQ.0) nseg_bdy_max=nseg_bdy_max+1
     IF(KNELM(2,NE).EQ.0) nseg_bdy_max=nseg_bdy_max+1
     IF(KNELM(3,NE).EQ.0) nseg_bdy_max=nseg_bdy_max+1
!     IF(nrank.EQ.0.AND.MOD(NE-1,100).EQ.0) &
!          WRITE(6,'(A,7I8)') 'wf_set_bdy NE,NNs,NEs:', &
!          NE,NN1,NN2,NN3,KNELM(1,NE),KNELM(2,NE),KNELM(3,NE)
  ENDDO
  nseg_max=(3*nelm_max-nseg_bdy_max)/2+nseg_bdy_max

  call wf_nelm_allocate
  call wf_nseg_allocate

!  ----- SET nseg_nside_nelm -----
!  KANOD :: IF NODE IS ON THE BOUNDARY, 1

  if(nrank.eq.0) WRITE(6,*) '------- wf_set_bdy set nseg_nside_nelm start ---'

  NSD=0

  DO NE=1,nelm_max
     NN1=node_nside_nelm(1,NE)
     NN2=node_nside_nelm(2,NE)
     NN3=node_nside_nelm(3,NE)

! --- CREATE FOR THE FIRST TIME, SIDE IS NEW ---

     IF(KNELM(1,NE).GE.NE) THEN
        KANOD(NN1)=0
        KANOD(NN2)=0
        NSD=NSD+1
        nseg_nside_nelm(1,NE)=NSD
        IF(NSD.LE.nseg_max) THEN
           INSID(NSD)=1
           NESID(NSD)=NE
           node_nseg(1,NSD)=NN1
           node_nseg(2,NSD)=NN2
        ENDIF
     ENDIF
     IF(KNELM(2,NE).GE.NE) THEN
        KANOD(NN2)=0
        KANOD(NN3)=0
        NSD=NSD+1
        nseg_nside_nelm(2,NE)=NSD
        IF(NSD.LE.nseg_max) THEN
           INSID(NSD)=2
           NESID(NSD)=NE
           node_nseg(1,NSD)=NN2
           node_nseg(2,NSD)=NN3
        ENDIF
     ENDIF
     IF(KNELM(3,NE).GE.NE) THEN
        KANOD(NN3)=0
        KANOD(NN1)=0
        NSD=NSD+1
        nseg_nside_nelm(3,NE)=NSD
        IF(NSD.LE.nseg_max) THEN
           INSID(NSD)=3
           NESID(NSD)=NE
           node_nseg(1,NSD)=NN3
           node_nseg(2,NSD)=NN1
        ENDIF
     ENDIF

! --- ALREADY CREATED ---

     IF(KNELM(1,NE).LT.NE.AND.KNELM(1,NE).GT.0) THEN
        NEL=KNELM(1,NE)
        DO ISD=1,3
           NSDL=ABS(nseg_nside_nelm(ISD,NEL))
           IF (NSDL.NE.0) THEN
              IF  ((node_nseg(1,NSDL).EQ.NN1).AND.&
                  &(node_nseg(2,NSDL).EQ.NN2)) THEN
                 nseg_nside_nelm(1,NE)=NSDL
              ELSEIF((node_nseg(1,NSDL).EQ.NN2).AND.&
                   & (node_nseg(2,NSDL).EQ.NN1)) THEN
                 nseg_nside_nelm(1,NE)=-NSDL
              END IF
           ENDIF
        ENDDO
     ENDIF
     IF(KNELM(2,NE).LT.NE.AND.KNELM(2,NE).GT.0) THEN
        NEL=KNELM(2,NE)
        DO ISD=1,3
           NSDL=ABS(nseg_nside_nelm(ISD,NEL))
           IF (NSDL.NE.0) THEN
              IF  ((node_nseg(1,NSDL).EQ.NN2).AND.&
                  &(node_nseg(2,NSDL).EQ.NN3)) THEN
                 nseg_nside_nelm(2,NE)=NSDL
              ELSEIF((node_nseg(1,NSDL).EQ.NN3).AND.&
                   & (node_nseg(2,NSDL).EQ.NN2)) THEN
                 nseg_nside_nelm(2,NE)=-NSDL
              ENDIF
           ENDIF
        ENDDO
     ENDIF
     IF(KNELM(3,NE).LT.NE.AND.KNELM(3,NE).GT.0) THEN
        NEL=KNELM(3,NE)
        DO ISD=1,3
           NSDL=ABS(nseg_nside_nelm(ISD,NEL))
           IF (NSDL.NE.0) THEN
              IF  ((node_nseg(1,NSDL).EQ.NN3).AND.&
                  &(node_nseg(2,NSDL).EQ.NN1)) THEN
                 nseg_nside_nelm(3,NE)=NSDL
              ELSEIF((node_nseg(1,NSDL).EQ.NN1).AND.&
                   & (node_nseg(2,NSDL).EQ.NN3)) THEN
                 nseg_nside_nelm(3,NE)=-NSDL
              ENDIF
           END IF
        ENDDO
     ENDIF

! --- SIDE IS ON THE BOUNDARY, SIDE IS NEW ---

     IF(KNELM(1,NE).EQ.0) THEN
        KANOD(NN1)=1
        KANOD(NN2)=1
        NSD=NSD+1
        nseg_nside_nelm(1,NE)=NSD
        IF(NSD.LE.nseg_max) THEN
           INSID(NSD)=1
           NESID(NSD)=NE
           node_nseg(1,NSD)=NN1
           node_nseg(2,NSD)=NN2
        ENDIF
     ENDIF
     IF(KNELM(2,NE).EQ.0) THEN
        KANOD(NN2)=1
        KANOD(NN3)=1
        NSD=NSD+1
        nseg_nside_nelm(2,NE)=NSD
        IF(NSD.LE.nseg_max) THEN
           INSID(NSD)=2
           NESID(NSD)=NE
           node_nseg(1,NSD)=NN2
           node_nseg(2,NSD)=NN3
        ENDIF
     ENDIF
     IF(KNELM(3,NE).EQ.0) THEN
        KANOD(NN3)=1
        KANOD(NN1)=1
        NSD=NSD+1
        nseg_nside_nelm(3,NE)=NSD
        IF(NSD.LE.nseg_max) THEN
           INSID(NSD)=3
           NESID(NSD)=NE
           node_nseg(1,NSD)=NN3
           node_nseg(2,NSD)=NN1
        ENDIF
     ENDIF
  ENDDO

  
!  DO NSD=1,nseg_max
!     write(*,*) "NESID,INSID,NSD",NESID(NSD),INSID(NSD),NSD
!  ENDDO

!  DO NE=1,nelm_max
!     DO I=1,3
!        write(*,*) "NE,I,nseg_nside_nelm",NE,I,nseg_nside_nelm(I,NE)
!     end DO
!  end DO

!    DO NE=1,nelm_max
!     DO I=1,3
!        write(*,*) "NE,I,KNELM",NE,I,KNELM(I,NE)
!     end DO
!  end DO

 ! ----- SET KASID -----

  if(nrank.eq.0) WRITE(6,*) '------- wf_set_bdy set KASID start ---'

  DO NSD=1,nseg_max
     KASID(NSD)=0
  ENDDO

  DO NE=1,nelm_max
     NSD1=ABS(nseg_nside_nelm(1,NE))
     NSD2=ABS(nseg_nside_nelm(2,NE))
     NSD3=ABS(nseg_nside_nelm(3,NE))
     IF(KNELM(1,NE).EQ.0) KASID(NSD1)=1
     IF(KNELM(2,NE).EQ.0) KASID(NSD2)=1
     IF(KNELM(3,NE).EQ.0) KASID(NSD3)=1
  ENDDO
  
  IERR=0
  RETURN
END SUBROUTINE wf_set_bdy

END MODULE wfindex
