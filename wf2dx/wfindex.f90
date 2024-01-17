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

FUNCTION wf_findex(x,y)

  use wfcomm
  implicit none
  real(rkind),INTENT(IN) :: x,y
  real(rkind) :: wf_findex,xn,yn

  xn=(x-xnode_min)/(xnode_max-xnode_min)
  yn=(y-ynode_min)/(ynode_max-ynode_min)
  wf_findex=0.999D3*xn+yn

  RETURN
END FUNCTION wf_findex

!     ***** SORT ELEMENTS BY SINDEX *****

SUBROUTINE wf_index
  
  use wfcomm
  USE wfsub
  USE wfsort
  implicit none
  integer :: nelm,nside,node
  real(rkind) :: xc,yc

  CALL wf_set_node_range

  DO nelm=1,nelm_max
     xc=0.D0
     yc=0.D0
     DO nside=1,3
        node=node_nside_nelm(nside,nelm)
        xc=xc+xnode(node)
        yc=yc+ynode(node)
     ENDDO
     sindex_nelm(nelm)=wf_findex(xc/3.D0,yc/3.D0)
  ENDDO
  
  DO nelm=1,nelm_max
     iv_nelm(nelm)=nelm
     id_nelm(nelm)=0
     nmed_nelm(nelm)=0
  ENDDO

!  DO nelm=1,nelm_max
!     WRITE(6,'(A,2I8,1P,E12.4)') & 
!          &        'nelm,IV,SINDEX=',nelm,IVELM(nelm),SINDEX(nelm)
!  ENDDO

  CALL wf_sort(nelm_max,sindex_nelm,wf_sort_subst,wf_sort_exchange)

  DO nelm=1,nelm_max
     iw_nelm(nelm)=0
  ENDDO
  DO nelm=1,nelm_max
     iw_nelm(iv_nelm(nelm))=nelm
  ENDDO

!  DO nelm=1,nelm_max
!     WRITE(6,'(A,3I8,1P,E12.4)') &
!          &        'nelm,IV,IW,SINDEX=',nelm,IVELM(nelm),IWELM(nelm),SINDEX(nelm)
!  ENDDO

  DO nelm=1,nelm_max
     IF (nrank.EQ.0) THEN
        IF(iw_nelm(nelm).EQ.0) WRITE(6,*) 'XX wf_index: iw_nelm undefined'
     END IF
  ENDDO
  
  CALL wf_set_node_range
  CALL wf_set_elm_area

  RETURN
END SUBROUTINE wf_index

!     ***** SET FEP DATA *****

SUBROUTINE wf_set_fep

  use wfcomm
  implicit none
  integer :: nelm,nside,node
  real(rkind) :: sindexL(3),xnodeL(3),ynodeL(3)
  real(rkind) :: smax,smin

  DO nelm=1,nelm_max
     DO nside=1,3
        node=node_nside_nelm(nside,nelm)
        sindexL(nside)=wf_findex(xnode(node),ynode(node))
        xnodeL(nside)=xnode(node)
        ynodeL(nside)=ynode(node)
     ENDDO
     sindex_min_nelm(nelm)=MIN(sindexL(1),sindexL(2),sindexL(3))
     sindex_max_nelm(nelm)=MAX(sindexL(1),sindexL(2),sindexL(3))
     xmin_nelm(nelm)=MIN(xnodeL(1),xnodeL(2),xnodeL(3))
     xmax_nelm(nelm)=MAX(xnodeL(1),xnodeL(2),xnodeL(3))
     ymin_nelm(nelm)=MIN(ynodeL(1),ynodeL(2),ynodeL(3))
     ymax_nelm(nelm)=MAX(ynodeL(1),ynodeL(2),ynodeL(3))
  ENDDO
  
  SMAX=sindex_max_nelm(1)
  DO nelm=1,nelm_max
     IF(sindex_max_nelm(nelm).GT.SMAX) SMAX=sindex_max_nelm(nelm)
     sindex_max_nelm(nelm)=SMAX
  ENDDO
  SMIN=sindex_min_nelm(nelm_max)
  DO nelm=nelm_max,1,-1
     IF(sindex_min_nelm(nelm).LT.SMIN) SMIN=sindex_min_nelm(nelm)
     sindex_min_nelm(nelm)=SMIN
  ENDDO
  
  RETURN
END SUBROUTINE wf_set_fep

!     ******* FIND ELEMENT INCLUDING NODES N1,N2 *******

SUBROUTINE wf_find_elm(nelm_s,node1,node2,nelm)

  use wfcomm
  implicit none
  integer,intent(in):: nelm_s,node1,node2
  integer,intent(out):: nelm
  integer :: I,J,K,L,node1L,node2L

  IF(nelm_s.LT.0.OR.nelm_s.GT.nelm_max) GOTO 9000
  DO I=1,MAX(nelm_max-nelm_s,nelm_s)
     DO J=1,2
        IF(J.EQ.1) THEN
           nelm=nelm_s+I
        ELSE
           nelm=nelm_s-I
        ENDIF
        IF(nelm.GE.1.AND.nelm.LE.nelm_max) THEN
           DO K=1,3
              node1L=node_nside_nelm(K,nelm)
              IF(node1L.EQ.node1) THEN
                 DO L=1,3
                    node2L=node_nside_nelm(L,nelm)
                    IF(node2L.EQ.node2) RETURN
                 ENDDO
              ENDIF
           ENDDO
        ENDIF
     ENDDO
  ENDDO
  
  nelm=0
  RETURN
  
9000 nelm=0
  RETURN
END SUBROUTINE wf_find_elm

!     ******* FIND ELEMENT INCLUDING POINT (x,y) *******

SUBROUTINE wf_fep(x,y,nelm)

  use wfcomm
  USE wfsort
  USE wfsub
  implicit none
  real(rkind),intent(in) :: x,y
  integer,intent(inout):: nelm
  integer :: ICOUNT,INMIN,NELMAX,NELMIN,nelmS,I,IDELT,J
  real(rkind),parameter :: EPS = 1.d-12
  real(rkind) :: WGT(3),WGTMIN,SIDX
  INTEGER:: IN,NN
  
  IF(nelm.NE.0) THEN   ! with initial guess

!    move to a neighboring element colser to the point

     ICOUNT=0
     
1000 CONTINUE
     ICOUNT=ICOUNT+1
     IF(ICOUNT.GT.nelm_max/3)  GOTO 2000 ! Could not find, try another method

     CALL wf_set_weight(nelm,x,y,WGT)

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
           WRITE(6,'(A,I8,1P5E12.4)') 'nelm,WGT,x,y:', &
                        nelm,WGT(1),WGT(2),WGT(3),x,y
           DO IN=1,3
              NN=node_nside_nelm(IN,nelm)
              WRITE(6,'(A,3I8,1P2E12.4)') 'nelm,IN,NN,x,y,:', &
                        nelm,IN,NN,xnode(NN),ynode(NN)
           END DO
        END if
     ENDIF
     
     IF(WGTMIN.GE.-EPS) RETURN
     
     nelm=nelm1_nside_nelm(INMIN,nelm)
     IF(nelm.LE.0) GOTO 2000 ! out of region
     GOTO 1000
  ENDIF

! without initial guess or could not find by sequential scheme
! find the element using the index of elements

2000 CONTINUE
  SIDX=wf_findex(x,y)
  IF(SIDX.GE.sindex_min_nelm(1)) THEN
     IF(SIDX.GT.sindex_min_nelm(nelm_max)) THEN
        NELMAX=nelm_max
     ELSE
        CALL wf_eval_index_min(sindex_min_nelm,nelm_max,SIDX,NELMAX)
     ENDIF
     IF(SIDX.LE.sindex_max_nelm(nelm_max)) THEN
        IF(SIDX.LT.sindex_max_nelm(1)) THEN
           NELMIN=1
        ELSE
           CALL wf_eval_index_min(sindex_max_nelm,nelm_max,SIDX,NELMIN)
        ENDIF
        
        nelmS=(NELMIN+NELMAX)/2
        DO I=0,MAX(NELMAX-nelmS,nelmS-NELMIN)+1
           IDELT=I
           DO J=1,2
              IDELT=-IDELT
              nelm=nelmS+IDELT
              IF(nelm.GE.1.AND.nelm.LE.nelm_max) THEN
                 IF(x.GE.xmin_nelm(nelm).AND.x.LE.xmax_nelm(nelm).AND.&
                    y.GE.ymin_nelm(nelm).AND.y.LE.ymax_nelm(nelm)) THEN
                    CALL wf_set_weight(nelm,x,y,WGT)
                    WGTMIN=MIN(WGT(1),WGT(2),WGT(3))
                    if (nrank.eq.0) then
                       IF(IDEBUG.EQ.4) THEN
                          WRITE(6,'(A,I8,1P3E12.4)') &
                                     'nelm,x,y,WGTMIN=', &
                                      nelm,x,y,WGTMIN
                       end if
                    END if
                    IF(WGTMIN.GE.-EPS) RETURN
                 ENDIF
              ENDIF
           ENDDO
        ENDDO
     ENDIF
  ENDIF
  nelm=0
  RETURN
END SUBROUTINE wf_fep

!     ****** Set Boundary Attribute for Side and Node ******

SUBROUTINE wf_set_bdy(IERR)

  use wfcomm
  USE femcomm
  implicit none
  integer,intent(out) :: IERR
  integer :: nside,nseg,nelm,nbdy,node
  integer :: node1,node2,node3,nseg1,nseg2,nseg3
  integer :: nelmL,nsegL

! ----- SET nseg_max & nelm1_nside_nelm -----
! nbdy_max :: Number of Boundary Side
! nelm1_nside_nelm(nside,NE) :: THE ELEMENT WHICH SHARE THE SIDE (nside,NE)
!                  IF nelm1_nside_nelm.EQ.0, SIDE IS BOUNDARY

  if(nrank.eq.0) WRITE(6,*) '------- wf_set_bdy set nseg_max & nelm1_nside_nelm start ---'

!  WRITE(6,'(A,2I8)') 'nelm_max,node_max=',nelm_max,node_max

  nbdy=0
  DO nelm=1,nelm_max
     node1=node_nside_nelm(1,nelm)
     node2=node_nside_nelm(2,nelm)
     node3=node_nside_nelm(3,nelm)
     CALL wf_find_elm(nelm,node1,node2,nelm1_nside_nelm(1,nelm))
     CALL wf_find_elm(nelm,node2,node3,nelm1_nside_nelm(2,nelm))
     CALL wf_find_elm(nelm,node3,node1,nelm1_nside_nelm(3,nelm))
     IF(nelm1_nside_nelm(1,nelm).EQ.0) nbdy=nbdy+1
     IF(nelm1_nside_nelm(2,nelm).EQ.0) nbdy=nbdy+1
     IF(nelm1_nside_nelm(3,nelm).EQ.0) nbdy=nbdy+1
  ENDDO
  nbdy_max=nbdy
  nseg_max=(3*nelm_max-nbdy_max)/2+nbdy_max

  call fem_mesh_allocate

!  ----- SET nseg_nside_nelm -----
!  mode_node :: IF NODE IS ON THE BOUNDARY, 1

  if(nrank.eq.0) WRITE(6,*) '------- wf_set_bdy set nseg_nside_nelm start ---'

  nseg=0

  DO node=1,node_max
     mode_node(node)=0
  END DO
  
  DO nelm=1,nelm_max
     node1=node_nside_nelm(1,nelm)
     node2=node_nside_nelm(2,nelm)
     node3=node_nside_nelm(3,nelm)

! --- CREATE nseg FOR THE FIRST TIME, SIDE IS NEW ---

     IF(nelm1_nside_nelm(1,nelm).GE.nelm) THEN  ! nelm over nside 1
        nseg=nseg+1
        nseg_nside_nelm(1,nelm)=nseg
        IF(nseg.LE.nseg_max) THEN
           node_nseg(1,nseg)=node1
           node_nseg(2,nseg)=node2
        ENDIF
     ENDIF
     IF(nelm1_nside_nelm(2,nelm).GE.nelm) THEN  ! nelm over nside 2
        nseg=nseg+1
        nseg_nside_nelm(2,nelm)=nseg
        IF(nseg.LE.nseg_max) THEN
           node_nseg(1,nseg)=node2
           node_nseg(2,nseg)=node3
        ENDIF
     ENDIF
     IF(nelm1_nside_nelm(3,nelm).GE.nelm) THEN  ! nelm over nside 3
        mode_node(node3)=0
        mode_node(node1)=0
        nseg=nseg+1
        nseg_nside_nelm(3,nelm)=nseg
        IF(nseg.LE.nseg_max) THEN
           node_nseg(1,nseg)=node3
           node_nseg(2,nseg)=node1
        ENDIF
     ENDIF

! --- Adjust the sign of nseg_nside_nelm ---

     nelmL=nelm1_nside_nelm(1,nelm)
     IF(nelmL.GT.0.AND.nelmL.LT.nelm) THEN
        DO nside=1,3
           nsegL=ABS(nseg_nside_nelm(nside,nelmL))
           IF (nsegL.NE.0) THEN
              IF    ((node_nseg(1,nsegL).EQ.node1).AND.&
                     (node_nseg(2,nsegL).EQ.node2)) THEN
                 nseg_nside_nelm(1,nelm)=nsegL
              ELSEIF((node_nseg(1,nsegL).EQ.node2).AND.&
                     (node_nseg(2,nsegL).EQ.node1)) THEN
                 nseg_nside_nelm(1,nelm)=-nsegL
              END IF
           ENDIF
        ENDDO
     ENDIF
     nelmL=nelm1_nside_nelm(2,nelm)
     IF(nelmL.GT.0.AND.nelmL.LT.nelm) THEN
        DO nside=1,3
           nsegL=ABS(nseg_nside_nelm(nside,nelmL))
           IF (nsegL.NE.0) THEN
              IF    ((node_nseg(1,nsegL).EQ.node2).AND.&
                     (node_nseg(2,nsegL).EQ.node3)) THEN
                 nseg_nside_nelm(2,nelm)=nsegL
              ELSEIF((node_nseg(1,nsegL).EQ.node3).AND.&
                     (node_nseg(2,nsegL).EQ.node2)) THEN
                 nseg_nside_nelm(2,nelm)=-nsegL
              ENDIF
           ENDIF
        ENDDO
     ENDIF
     nelmL=nelm1_nside_nelm(3,nelm)
     IF(nelmL.LT.nelm.AND.nelmL.GT.0) THEN
        DO nside=1,3
           nsegL=ABS(nseg_nside_nelm(nside,nelmL))
           IF (nsegL.NE.0) THEN
              IF    ((node_nseg(1,nsegL).EQ.node3).AND.&
                     (node_nseg(2,nsegL).EQ.node1)) THEN
                 nseg_nside_nelm(3,nelm)=nsegL
              ELSEIF((node_nseg(1,nsegL).EQ.node1).AND.&
                     (node_nseg(2,nsegL).EQ.node3)) THEN
                 nseg_nside_nelm(3,nelm)=-nsegL
              ENDIF
           END IF
        ENDDO
     ENDIF

! --- SIDE IS ON THE BOUNDARY, SIDE IS NEW ---

     IF(nelm1_nside_nelm(1,nelm).EQ.0) THEN
        mode_node(node1)=1
        mode_node(node2)=1
        nseg=nseg+1
        nseg_nside_nelm(1,nelm)=nseg
        IF(nseg.LE.nseg_max) THEN
           node_nseg(1,nseg)=node1
           node_nseg(2,nseg)=node2
        ENDIF
     ENDIF
     IF(nelm1_nside_nelm(2,nelm).EQ.0) THEN
        mode_node(node2)=1
        mode_node(node3)=1
        nseg=nseg+1
        nseg_nside_nelm(2,nelm)=nseg
        IF(nseg.LE.nseg_max) THEN
           node_nseg(1,nseg)=node2
           node_nseg(2,nseg)=node3
        ENDIF
     ENDIF
     IF(nelm1_nside_nelm(3,nelm).EQ.0) THEN
        mode_node(node3)=1
        mode_node(node1)=1
        nseg=nseg+1
        nseg_nside_nelm(3,nelm)=nseg
        IF(nseg.LE.nseg_max) THEN
           node_nseg(1,nseg)=node3
           node_nseg(2,nseg)=node1
        ENDIF
     ENDIF
  ENDDO

  nseg_max=nseg
  
!  DO nseg=1,nseg_max
!     write(*,*) "nelmSID,INSID,nseg",nelmSID(nseg),INSID(nseg),nseg
!  ENDDO

!  DO nelm=1,nelm_max
!     DO I=1,3
!        write(*,*) "nelm,I,nseg_nside_nelm",nelm,I,nseg_nside_nelm(I,nelm)
!     end DO
!  end DO

!    DO nelm=1,nelm_max
!     DO I=1,3
!        write(*,*) "nelm,I,nelm1_nside_nelm",nelm,I,nelm1_nside_nelm(I,nelm)
!     end DO
!  end DO

 ! ----- SET mode_nseg -----

  if(nrank.eq.0) WRITE(6,*) '------- wf_set_bdy set mode_nseg start ---'

  DO nseg=1,nseg_max
     mode_nseg(nseg)=0
  ENDDO

  DO nelm=1,nelm_max
     nseg1=ABS(nseg_nside_nelm(1,nelm))
     nseg2=ABS(nseg_nside_nelm(2,nelm))
     nseg3=ABS(nseg_nside_nelm(3,nelm))
     IF(nelm1_nside_nelm(1,nelm).EQ.0) mode_nseg(nseg1)=1
     IF(nelm1_nside_nelm(2,nelm).EQ.0) mode_nseg(nseg2)=1
     IF(nelm1_nside_nelm(3,nelm).EQ.0) mode_nseg(nseg3)=1
  ENDDO
  
  IERR=0
  RETURN
END SUBROUTINE wf_set_bdy

END MODULE wfindex
