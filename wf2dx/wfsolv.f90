! wfsolv.f90

MODULE wfsolv

  PRIVATE
  PUBLIC wf_defmlen
  PUBLIC wf_cvsolv
  PUBLIC wf_dtensr
  PUBLIC wf_mutensr
  

CONTAINS
  
!    ******* SET MLEN *******

SUBROUTINE wf_defmlen

  use wfcomm
  implicit none

!  nbdy_max=0
!  do NSD=1,nseg_max
!     if(KASID(NSD).eq.1) nbdy_max=nbdy_max+1
!  end do

!  nbdy_max=0
!  do NN=1,node_max
!     if(mode_node(NN).eq.1) nbdy_max=nbdy_max+1
!  end do

!  mtx_len=nseg_max+node_max-nbdy_max-nbdy_max
  
  mtx_len=nseg_max+node_max-2*nbdy_max

  call wf_solve_allocate

  RETURN
END SUBROUTINE wf_defmlen

!     ****** SOLV MATRIX EQUATION *****

SUBROUTINE wf_cvsolv

  use wfcomm
  use libmpi
  use libmtx
  use libqsort
  implicit none
  integer :: nelm,node,nseg,nvar,nvar_max,nbdy
  INTEGER:: node1,node2,nside
  integer :: I,J,KK,LL
  integer :: nsegI,nsegJ,nodeI,nodeJ,nvarI,nvarJ
  integer :: itype
  integer :: its
  integer :: JMIN,JMAX,MILEN,MJLEN
  integer :: NNZ,NNZMAX,NNZME      !Number of Non-Zero Matrix Element
  integer,dimension(:),ALLOCATABLE :: NEFLAG
  integer :: ORIENTJ,ORIENTI
  real :: cputime1,cputime2
  real(rkind) :: x,y,val
  complex(rkind):: CEB
  complex(rkind),dimension(:),ALLOCATABLE :: CRVP,CEQP
  integer,dimension(:),ALLOCATABLE :: NSEQ
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: value_sort
  INTEGER,DIMENSION(:),ALLOCATABLE:: nvar_sort
  INTEGER,DIMENSION(:),ALLOCATABLE:: ntype_nvar,npos_nvar
  INTEGER:: IX,IY,nvar_old

  ! ----- initialize ------
  
  do nvar=1,mtx_len
     CSV(nvar) =(0.d0,0.d0)
  end do

  ! --- decide istart,iend ---

  IF(nrank.EQ.0) write(6,*) 'mtx_len=',mtx_len
  call mtxc_setup(mtx_len,istart,iend,0)
  call mtxc_cleanup

  ! ----- set NV_NSD & NV_NN -----

  ALLOCATE(ntype_nvar(node_max+nseg_max))
  ALLOCATE(npos_nvar(node_max+nseg_max))

  nvar=0
  DO nseg=1,nseg_max
     IF(mode_nseg(nseg).EQ.1) THEN
        nvar_nseg(nseg)=0
     ELSE
        nvar=nvar+1
        nvar_nseg(nseg)=nvar
        ntype_nvar(nvar)=2
        npos_nvar(nvar)=nseg
     END IF
  END DO
  DO node=1,node_max
     IF(mode_node(node).EQ.1) THEN
        nvar_node(node)=0
     ELSE
        nvar=nvar+1
        nvar_node(node)=nvar
        ntype_nvar(nvar)=1
        npos_nvar(nvar)=node
     END IF
  END DO
  nvar_max=nvar

  ALLOCATE(nvar_sort(nvar_max),value_sort(nvar_max))

  DO nvar=1,nvar_max
     IF(ntype_nvar(nvar).EQ.1) THEN
        node=npos_nvar(nvar)
        x=xnode(node)
        y=ynode(node)
     ELSE
        nseg=npos_nvar(nvar)
        node1=node_nseg(1,nseg)
        node2=node_nseg(2,nseg)
        x=0.5D0*(xnode(node1)+xnode(node2))
        y=0.5D0*(ynode(node1)+ynode(node2))
     END IF
     val=sort_weight_x*x+sort_weight_y*y
     nvar_sort(nvar)=nvar
     value_sort(nvar)=val
  END DO

!  CALL qsort_dl(value_sort,nv_sort)
  CALL qsort_di(value_sort,nvar_sort)

  DO nvar=1,nvar_max
     nvar_old=nvar_sort(nvar)
     IF(ntype_nvar(nvar_old).EQ.1) THEN
        node=npos_nvar(nvar_old)
        nvar_node(node)=nvar
!        X=xnode(node)
!        Y=ynode(node)
!        VAL=sort_weight_x*x+sort_weight_y*y
!        write(21,'(I10,A,I10,1P3E12.4)') nvar,' node ',node,x,y,val
     ELSE
        nseg=npos_nvar(nvar_old)
        nvar_nseg(nseg)=nvar
!        x=xcenter_nseg(nseg)
!        y_ysenter_nseg(nseg)
!        VAL=sort_weight_x*x+sort_weight_y*y
!        write(21,'(I10,A,I10,1P3E12.4)') nvar,' nseg ',nseg,x,y,val
     END IF
  END DO
     
  DEALLOCATE(nvar_sort,value_sort,ntype_nvar,npos_nvar)

!  do NSD=1,nseg_max
!     write(*,*) NSD,KASID(NSD),NV_NSD(NSD)
!  end do
!  do NN=1,node_max
!     write(*,*) NN,mode_node(NN),NV_NN(NN)
!  end do

  ! ----- set NEFLAG ------

  ALLOCATE(NEFLAG(nelm_max))
  DO nelm=1,nelm_max
     NEFLAG(nelm)=0
  END DO

  DO nelm=1,nelm_max
     DO nside=1,3
        nseg=ABS(nseg_nside_nelm(nside,nelm))
        nvar=nvar_nseg(nseg)
        IF(nvar.GE.istart.AND.&
           nvar.LE.iend ) then
           NEFLAG(nelm)=1
        END IF
     END DO

     IF(NEFLAG(nelm).eq.0) then
        DO nside=1,3
           node=node_nside_nelm(nside,nelm)
           nvar=nvar_node(node)
           IF(nvar.GE.istart.AND.&
              nvar.LE.iend) THEN
              NEFLAG(nelm)=1
           END IF
        END DO
     END IF

  END DO

  ! ----- set MJLEN -----

  JMIN=mtx_len
  JMAX=0
  DO nelm=1,nelm_max

     IF(NEFLAG(nelm).EQ.0) GOTO 8100

     LL=0
     DO J=1,6
        IF(J.GE.1.AND.J.LE.3) THEN
           nsegJ=ABS(nseg_nside_nelm(J,nelm))
           nvarJ=nvar_nseg(nsegJ)
        else
           nodeJ=node_nside_nelm(J-3,nelm)
           nvarJ=nvar_node(nodeJ)
        end if
        if (nvarJ.eq.0) goto 8110
        LL=nvarJ

        KK=0
        DO I=1,6
           IF(I.GE.1.AND.I.LE.3) THEN
              nsegI=ABS(nseg_nside_nelm(I,nelm))
              nvarI=nvar_nseg(nsegI)
           ELSE
              nodeI=node_nside_nelm(I-3,nelm)
              nvarI=nvar_node(nodeI)
           END IF
           IF(nvarI.EQ.0) GOTO 8120
           KK=nvarI

           IF((KK.GE.istart).AND.&
              (KK.LE.iend  )) THEN
              IF(LL.LT.JMIN) JMIN=LL
              IF(LL.GT.JMAX) JMAX=LL
           END IF

8120       CONTINUE
        ENDDO
8110    CONTINUE
     ENDDO
     
8100 CONTINUE
  END do

! ------ Count non-zero component -------

  NNZ=0

  DO nelm=1,nelm_max
     IF(NEFLAG(nelm).NE.0) THEN

        LL=0
        DO J=1,6
           IF(J.GE.1.AND.J.LE.3) THEN
              nsegJ=nseg_nside_nelm(J,nelm)
              IF(nsegJ.LT.0) THEN
                 nsegJ=-nsegJ
              END IF
              nvarJ =nvar_nseg(nsegJ)
           ELSE
              nodeJ=node_nside_nelm(J-3,nelm)
              nvarJ=nvar_node(nodeJ)
           END IF
           LL=nvarJ

           IF ((LL.GE.JMIN).AND.(LL.LE.JMAX)) THEN

              KK=0
              DO I=1,6
                 IF(I.GE.1.AND.I.LE.3) THEN
                    nsegI=nseg_nside_nelm(I,nelm)
                    IF(nsegI.LT.0) THEN
                       nsegI=-nsegI
                    END IF
                    nvarI =nvar_nseg(nsegI)
                 ELSE
                    nodeI=node_nside_nelm(I-3,nelm)
                    nvarI=nvar_node(nodeI)
                 END IF
                 KK=nvarI

                 IF((KK.GE.istart).AND.&
                      (KK.LE.iend  )) THEN
                    NNZ=NNZ+1
                 END IF
              END DO
           END IF
        ENDDO
     END IF
  END DO

  NNZMAX=NNZ
  MILEN=iend-istart+1
  MJLEN=JMAX-JMIN+1

  ! ----- set CEQP,CRVP -----

  allocate(CEQP(NNZMAX),NSEQ(NNZMAX),CRVP(MILEN))
  DO NNZ=1,NNZMAX
     CEQP(NNZ)=(0.d0,0.d0)
     NSEQ(NNZ)=0  ! (i-istart)*MJLEN+j-jmin
  END DO
  DO I=1,MILEN
     CRVP(I)=(0.d0,0.d0)
  END DO

! ------ set grobal matrix -------

  NNZ=0

  DO nelm=1,nelm_max
     IF(NEFLAG(nelm).EQ.0) GOTO 8000
     CALL wf_cmcalc(nelm)
     
!    === ASSEMBLY ===
!    If KK (or LL) is out of assigned range,  
!      CVSOLV do not save the matrix element. 

!    --- inside of the boundary ---
     LL=0
     DO J=1,6
        ORIENTJ=1
        IF(J.GE.1.AND.J.LE.3) THEN
           nsegJ=nseg_nside_nelm(J,nelm)
           IF(nsegJ.LT.0) THEN
              nsegJ=-nsegJ
              ORIENTJ=-1
           END IF
           nvarJ =nvar_nseg(nsegJ)
        ELSE
           nodeJ=node_nside_nelm(J-3,nelm)
           nvarJ=nvar_node(nodeJ)
        end if
        LL=nvarJ

        IF ((LL.GE.JMIN).AND.(LL.LE.JMAX)) THEN

           KK=0
           DO I=1,6
              ORIENTI=1
              IF(I.GE.1.AND.I.LE.3) THEN
                 nsegI=nseg_nside_nelm(I,nelm)
                 IF(nsegI.LT.0) THEN
                    nsegI=-nsegI
                    ORIENTI=-1
                 END IF
                 nvarI =nvar_nseg(nsegI)
              ELSE
                 node=node_nside_nelm(I-3,nelm)
                 nvarI=nvar_node(nodeI)
              END IF
              KK=nvarI
              IF((KK.GE.istart).AND.&
                 (KK.LE.iend  )) THEN
                 IF(ABS(CM(I,J)).NE.0.d0) THEN 
                    NNZ=NNZ+1
                    IX=KK-istart
                    IY=MJLEN
                    NSEQ(NNZ)=IX*IY+LL-jmin
                    CEQP(NNZ)=ORIENTJ*ORIENTI*CM(I,J)
                 END if
              end if
           END DO
        END if
     ENDDO

!    --- Contribution from the boundary electric field ---

     LL=0
     DO J=1,6
        ORIENTJ=1
        IF(J.GE.1.AND.J.LE.3) THEN
           nsegJ=nseg_nside_nelm(J,nelm)
           IF(nsegJ.LT.0) THEN
              nsegJ=-nsegJ
              ORIENTJ=-1
           END IF
           nbdy =nbdy_nseg(nsegJ)
           IF(nbdy.NE.0) CEB=CESD_nbdy(nbdy)
        ELSE
           nodeJ=node_nside_nelm(J-3,nelm)
           nbdy=nbdy_node(nodeJ)
           IF(nbdy.NE.0) CEB=CEND_nbdy(nbdy)
        END IF

        IF(nbdy.NE.0) THEN
           IF(ABS(CEB).GT.0.D0) THEN
              KK=0
              DO I=1,6
                 ORIENTI=1
                 IF(I.GE.1.AND.I.LE.3) THEN
                    nsegI=nseg_nside_nelm(I,nelm)
                    IF(nsegI.LT.0) THEN
                       nsegI=-nsegI
                       ORIENTI=-1
                    END IF
                    nvarI=nvar_nseg(nsegI)
                 ELSE
                    nodeI=node_nside_nelm(I-3,nelm)
                    nvarI=nvar_node(nodeI)
                 END IF
                 KK=nvarI
                 IF((KK.GE.istart).AND.&
                      (KK.LE.iend  )) THEN
                    CRVP(KK-istart+1)=CRVP(KK-istart+1) &
                         -ORIENTI*ORIENTJ*CM(I,J)*CEB
                 END IF
              END DO
           END IF
        END IF
     ENDDO

!    --- Contribution from the antenna current ---

     KK=0
     DO I=1,6
        ORIENTI=1
        IF(I.GE.1.AND.I.LE.3) THEN
           nsegI=nseg_nside_nelm(I,nelm)
           IF(nsegI.LT.0) THEN
              nsegI=-nsegI
              ORIENTI=-1
           END IF
           nvarI =nvar_nseg(nsegI)
        ELSE
           nodeI=node_nside_nelm(I-3,nelm)
           nvarI=nvar_node(nodeI)
        END IF
        KK=nvarI
        IF((KK.GE.istart).AND.&
           (KK.LE.iend  )) THEN
           CRVP(KK-istart+1)=CRVP(KK-istart+1)+ORIENTI*CVTOT(I,nelm)
        END IF
     ENDDO

8000 CONTINUE
  END DO

  IF(nrank.EQ.0) WRITE(6,*) 'wfsolv: sort started'
  CALL qsort_ic(NSEQ,CEQP)
  !  CALL qsort_lc(NSEQ,CEQP)
  
  IF(nrank.EQ.0) WRITE(6,*) 'wfsolv: reduction started'
  NNZME=1
  DO NNZ=2,NNZMAX
     IF(NSEQ(NNZ).EQ.NSEQ(NNZME)) THEN
        CEQP(NNZME)=CEQP(NNZME)+CEQP(NNZ)
     ELSE
        NNZME=NNZME+1
        NSEQ(NNZME)=NSEQ(NNZ)
        CEQP(NNZME)=CEQP(NNZ)
     END IF
  END DO

  IF(nrank.EQ.0) WRITE(6,'(A77)') &
  '      nrank     istart       iend      MILEN      MJLEN     NNZMAX      NNZME'
  CALL mtx_barrier
  WRITE(6,'(7I11)') nrank,istart,iend,iend-istart+1, &
                    JMAX-JMIN+1,NNZMAX,NNZME

  ! ----- initialize for parallel computing -----

  itype = 0
  CALL mtxc_setup(mtx_len,istart,iend,nzmax=NNZME)

  DO NNZ=1,NNZME
     IF(ABS(CEQP(NNZ)).GT.0.D0) THEN
        i=NSEQ(NNZ)/MJLEN
        j=NSEQ(NNZ)-i*MJLEN
        IF(i.LT.0.OR.i.GT.iend-istart) THEN
           WRITE(6,'(A/6I12)') 'NNZ,NSEQ(NNZ),i,istart,i+istart,iend=', &
                               NNZ,NSEQ(NNZ),i,istart,i+istart,iend
           STOP
        END IF
        CALL mtxc_set_matrix(i+istart,j+jmin,CEQP(NNZ))
     END IF
  END DO

  DO i=istart,iend
     CALL mtxc_set_source(i,CRVP(i-istart+1))
  END DO

  CALL GUTIME(cputime1)

  CALL mtxc_solve(itype,tolerance,its)
  !zmumps always return "its = 0"
  IF(nrank.EQ.0) WRITE(6,*) 'Iteration Number=',its

  CALL GUTIME(cputime2)

  CALL mtxc_gather_vector(CSV)

  DEALLOCATE(CEQP,NSEQ,CRVP)
  DEALLOCATE(NEFLAG)
  CALL mtxc_cleanup
  RETURN
END SUBROUTINE wf_cvsolv

!     ****** LOCAL ELEMENT MATRIX  ******

  SUBROUTINE wf_cmcalc(NE)

  USE wfcomm
  IMPLICIT NONE
  INTEGER,INTENT(IN):: NE
  INTEGER:: I,J

  ! --- initialize CM ---

  DO J=1,6
     DO I=1,6
        CM(I,J)=(0.D0,0.D0)
     END DO
  END DO

  CALL CMCALCV(NE)

  CALL CMCALCS(NE)

  CALL CMCALCP(NE)

  RETURN
  END SUBROUTINE wf_cmcalc


!     ****** LOCAL ELEMENT MATRIX : volume integral ******

  SUBROUTINE cmcalcv(nelm)

    use wfcomm
    USE wfsub
  implicit none
  integer,intent(in) :: nelm
  integer :: I,J,K,M,N,nside,nseg
  complex(rkind) :: CM1(3,3)
  real(rkind) :: S,L(3)
  real(rkind) :: A(3),B(3),C(3),AW(3),BW(3),CW(3)
  real(rkind) :: x(3),y(3)

  ! --- initialize ---

  S=area_nelm(nelm)
  CALL wf_set_abc(nelm,A,B,C)
  CALL wf_set_node(nelm,x,y)

  do nside=1,3
     nseg=ABS(nseg_nside_nelm(nside,nelm))
     IF(model_wf.EQ.0) THEN
        L(nside)=len_nseg(nseg)
     ELSE
        IF(nseg_nside_nelm(nside,nelm).GT.0.D0) THEN
           L(nside)=len_nseg(nseg)
        ELSE
           L(nside)=-len_nseg(nseg)
        END IF
     END IF
  end do

  DO nside=1,3
     M=nside
     N=nside+1
     if(N.gt.3) N=N-3
     AW(nside)=L(nside)*(A(M)*B(N)-A(N)*B(M))
     BW(nside)=L(nside)*(B(M)*C(N)-B(N)*C(M))
     CW(nside)=L(nside)*(C(M)*A(N)-C(N)*A(M))
  end do

  ! ----- rotErotF term -----

  ! --- E1F1 ---

  do J=1,3
     do I=1,3
        CM1(I,J)=(0.d0,0.d0)
     end do
  end do

  SELECT CASE(MODELG)
  CASE(0,12)
     do K=1,3 
        do J=1,3
           do I=1,3
              CM1(I,J)=CM1(I,J) &
                      +(RKZ**2)*RR &
                       *( AW(I)*AW(J)-(AW(I)*BW(J)+AW(J)*BW(I))*y(K) &
                         +CW(I)*CW(J)-(BW(I)*CW(J)+BW(J)*CW(I))*x(K) &
                         +BW(I)*BW(J)*(x(K)**2+y(K)**2)) &
                       *S*AIF1(K) &
                      +4.d0*BW(I)*BW(J)*RR*S*AIF1(K)
           end do
        end do
     end do
  CASE(1:10,13)
     do K=1,3 
        do J=1,3
           do I=1,3
              CM1(I,J)=CM1(I,J) &
                      +(real(NPH)**2)/x(K) &
                       *( AW(I)*AW(J)-(AW(I)*BW(J)+AW(J)*BW(I))*y(K) &
                         +CW(I)*CW(J)-(BW(I)*CW(J)+BW(J)*CW(I))*x(K) &
                         +BW(I)*BW(J)*(x(K)**2+y(K)**2)) &
                       *S*AIF1(K) &
                      +4.d0*BW(I)*BW(J)*x(K)*S*AIF1(K)
           end do
        end do
     end do
  END SELECT
  do J=1,3
     do I=1,3
        CM(I,J)=CM(I,J)+CM1(I,J)
     end do
  end do

  ! --- E1F2 --- 

  do J=1,3
     do I=1,3
        CM1(I,J)=(0.d0,0.d0)
     end do
  end do

  SELECT CASE(MODELG)
  CASE(0,12)
     do K=1,3
        do J=1,3
           do I=1,3
              CM1(I,J)=CM1(I,J) &
                      +(CII*RKZ*RR) &
                       *(-B(I) &
                          *(AW(J)-BW(J)*y(K)) &
                         +C(I) &
                          *(CW(J)-BW(J)*x(K))) &
                       *S*AIF1(K)
           end do
        end do
     end do
  CASE(1:10,13)
     do K=1,3
        do J=1,3
           do I=1,3
              CM1(I,J)=CM1(I,J) &
!
                       +(CII*real(NPH)) &
                       *(-B(I) &
                          *(AW(J)-BW(J)*y(K)) &
                         +C(I) &
                          *(CW(J)-BW(J)*x(K))) &
                       *S*AIF1(K)&
!
                      +(CII*real(NPH))&
                       *(-(AW(J)-BW(J)*y(K))/x(K))&
                       *S*AIF2(I,K)
           end do
        end do
     end do
  END SELECT

  do J=1,3
     do I=1,3
        CM(I+3,J)=CM(I+3,J)+CM1(I,J)
     end do
  end do

  ! --- E2F1 ---

  do J=1,3
     do I=1,3
        CM1(I,J)=(0.d0,0.d0)
     end do
  end do

  SELECT CASE(MODELG)
  CASE(0,12)
     do K=1,3
        do J=1,3
           do I=1,3
              CM1(I,J)=CM1(I,J) &
                      -(CII*RKZ*RR) &
                       *(-B(J) &
                          *(AW(I)-BW(I)*y(K)) &
                         +C(J)&
                          *(CW(I)-BW(I)*x(K))) &
                        *S*AIF1(K)
           end do
        end do
     end do
  CASE(1:10,13)
     do K=1,3
        do J=1,3
           do I=1,3
              CM1(I,J)=CM1(I,J) &
!
                      -(CII*real(NPH)) &
                       *(-B(J) &
                          *(AW(I)-BW(I)*y(K)) &
                         +C(J)&
                          *(CW(I)-BW(I)*x(K))) &
                        *S*AIF1(K) &
!
                      -(CII*real(NPH)) &
                       *(-(AW(I)-BW(I)*y(K))/x(K)) &
                        *S*AIF2(J,K)
           end do
        end do
     end do
  END SELECT
  do J=1,3
     do I=1,3
        CM(I,J+3)=CM(I,J+3)+CM1(I,J)
     end do
  end do

  ! --- E2F2 ---

  do J=1,3
     do I=1,3
        CM1(I,J)=(0.d0,0.d0)
     end do
  end do

  SELECT CASE(MODELG)
  CASE(0,12)
     do K=1,3
        do J=1,3
           do I=1,3
              CM1(I,J)=CM1(I,J) &
                      +(B(I)*B(J)+C(I)*C(J))*RR*S*AIF1(K)
!                      +B(J)*S*AIF2(I,K) &
!                      +B(I)*S*AIF2(J,K) &
!                      +1.D0/RR*S*AIF3(I,J,K)
           end do
        end do
     end do
  CASE(1:10,13)
     do K=1,3
        do J=1,3
           do I=1,3
              CM1(I,J)=CM1(I,J) &
                      +(B(I)*B(J)+C(I)*C(J))*x(K)*S*AIF1(K) &
                      +B(J)*S*AIF2(I,K) &
                      +B(I)*S*AIF2(J,K) &
                      +1.D0/(x(K))*S*AIF3(I,J,K)
           end do
        end do
     end do
  END SELECT
  do J=1,3
     do I=1,3
        CM(I+3,J+3)=CM(I+3,J+3)+CM1(I,J)
     end do
  end do

  RETURN
  END SUBROUTINE cmcalcv

!     ****** LOCAL ELEMENT MATRIX : Surface integral ******

  SUBROUTINE cmcalcs(nelm)

    use wfcomm
    USE wfsub
  implicit none
  integer,intent(in) :: nelm
  integer :: I,J,K,I1,J1,K1,ID,nside,nseg,IND(2),node1,node2
  real(rkind) :: S,L(3),xn,yn
  real(rkind) :: A(3),B(3),C(3),AW(3),BW(3),CW(3)
  real(rkind) :: x(3),y(3)
  COMPLEX(rkind):: CTEMP

  ! --- if no boundary side, return ---

  ID=0
  DO nside=1,3
     nseg=ABS(nseg_nside_nelm(nside,nelm))
     IF(mode_nseg(nseg).EQ.1) ID=1
  END DO
  IF(ID.EQ.0) RETURN

  ! --- initialize ---

  S=area_nelm(nelm)
  call wf_set_abc(nelm,A,B,C)
  call wf_set_node(nelm,x,y)

  DO nside=1,3
     nseg=ABS(nseg_nside_nelm(nside,nelm))
     IF(mode_nseg(nseg).EQ.1) THEN
        IF(model_wf.EQ.0) THEN
           L(nside)=len_nseg(nseg)
        ELSE
           IF(nseg_nside_nelm(nside,nelm).GT.0.D0) THEN
              L(nside)=len_nseg(nseg)
           ELSE
              L(nside)=-len_nseg(nseg)
           END IF
        END IF

        node1=nside
        node2=nside+1
        IF(node2.GT.3) node2=node2-3
        AW(nside)=L(nside)*(A(node1)*B(node2)-A(node2)*B(node1))
        BW(nside)=L(nside)*(B(node1)*C(node2)-B(node2)*C(node1))
        CW(nside)=L(nside)*(C(node1)*A(node2)-C(node2)*A(node1))
     END IF
  END DO

  DO nside=1,3
     nseg=ABS(nseg_nside_nelm(nside,nelm))
     IF(mode_nseg(nseg).EQ.1) THEN

        IND(1)=nside
        IND(2)=nside+1
        IF(IND(2).GT.3) IND(2)=IND(2)-3
        xn= (y(IND(2))-y(IND(1)))/L(nside)
        yn=-(x(IND(2))-x(IND(1)))/L(nside)

  ! ----- rotErotF term -----

  ! --- E1F1 ---

        I=nside
        J=nside
        SELECT CASE(MODELG)
        CASE(0,12)
           DO K1=1,2
              K=IND(K1)
              CTEMP=CM(I,J)
              CM(I,J)=CM(I,J) &
                      -2.D0*BW(J)*(yn*(AW(I)-BW(I)*y(K)) &
                                  +xn*(CW(I)-BW(I)*x(K)))*RR*L(nside)*AIE1(K)
!                    WRITE(21,'(A,I10,3I5,1P4E12.4)') &
!                         'CM:',NE,I,J,K,CM(I,J),CM(I,J)-CTEMP
!                    WRITE(21,'(5X,1P5E12.4)') &
!                         BW(J),ZN,AW(I),BW(I),Z(K)
!                    WRITE(21,'(5X,1P5E12.4)') &
!                         RN,CW(I),R(K),RR,L(ISD),AIE1(K)
           END DO
        CASE(1:10,13)
           DO K1=1,2
              K=IND(K1)
              CM(I,J)=CM(I,J) &
                      -2.D0*BW(J)*(yn*(AW(I)-BW(I)*y(K)) &
                                  +xn*(CW(I)-BW(I)*x(K)))*x(K)*L(nside)*AIE1(K)
           END DO
        END SELECT

  ! --- E2F1 --- 

        I=nside
        DO J1=1,2
           J=IND(J1)
           SELECT CASE(MODELG)
           CASE(0,12)
              DO K1=1,2
                 K=IND(K1)
                 CTEMP=CM(I,J+3)
                 CM(I,J+3)=CM(I,J+3) &
                          +CII*RKZ*RR &
                           *(-yn*(-CW(J)+BW(J)*x(K)) &
                             -xn*( AW(J)-BW(J)*y(K)))*L(nside)*AIE2(J,K)
                    WRITE(21,'(A,I10,3I5,1P4E12.4)') &
                         'CM:',nelm,I,J+3,K,CM(I,J+3),CM(I,J+3)-CTEMP
              END DO
           CASE(1:10,13)
              DO K1=1,2
                 K=IND(K1)
                 CM(I,J+3)=CM(I,J+3) &
                          +CII*NPH &
                           *(-yn*(-CW(J)+BW(J)*x(K)) &
                             -xn*( AW(J)-BW(J)*y(K)))*L(nside)*AIE2(J,K)
              END DO
           END SELECT
        END DO

  ! --- E2F2 ---

        DO J1=1,2
           J=IND(J1)
           DO I1=1,2
              I=IND(I1)
              SELECT CASE(MODELG)
              CASE(0,12)
                 DO K1=1,2
                    K=IND(K1)
                    CTEMP=CM(I+3,J+3)
                    CM(I+3,J+3)=CM(I+3,J+3) &
                               +(xn*B(J)+yn*C(J))*RR &
                                *L(nside)*AIE3(I,J,K)
                    WRITE(21,'(A,I10,3I5,1P4E12.4)') &
                         'CM:',nelm,I+3,J+3,K,CM(I+3,J+3),CM(I+3,J+3)-CTEMP
                 END DO
              CASE(1:10,13)
                 DO K1=1,2
                    K=IND(K1)
                    CM(I+3,J+3)=CM(I+3,J+3) &
                               +((xn*B(J)+yn*C(J))*x(K) &
                                 +xn*(A(J)+B(J)*x(K)+C(J)*y(K))) &
                                *L(nside)*AIE3(I,J,K)
                 END DO
              END SELECT
           END DO
        END DO
     END IF
  END DO

  RETURN
  END SUBROUTINE cmcalcs

!     ****** LOCAL ELEMENT MATRIX : Surface integral ******

  SUBROUTINE cmcalcp(nelm)

    use wfcomm
    USE wfsub
  implicit none
  integer,intent(in) :: nelm
  integer :: I,J,K,NS,IN,JJ,II
  real(rkind) :: RW,WC,WC2
  real(rkind) :: S
  real(rkind) :: A(3),B(3),C(3)
  real(rkind) :: x(3),y(3),MU(3,3,6)
  COMPLEX(rkind):: CM2(6,6)
  complex(rkind) :: DTENS(NSM,3,3,3)
  complex(rkind) :: DTENST(3,3,3)

  ! --- initialize ---

  RW=2.D0*PI*RF*1.D6
  WC=RW/VC
  WC2=WC**2

  S=area_nelm(nelm)  
  call wf_set_abc(nelm,A,B,C)
  call wf_set_node(nelm,x,y)

  ! ----- dielectric tensor term -----

  call wf_dtensr(nelm,DTENS)
  call wf_mutensr(nelm,MU)

  DO J=1,6
     DO I=1,6
        CM2(I,J)=(0.D0,0.D0)
     END DO
  END DO

  do J=1,3
     do I=1,3
        do IN=1,3
           if(I.eq.J) then
              DTENST(IN,I,J)=(1.d0,0.d0)
           else
              DTENST(IN,I,J)=(0.d0,0.d0)
           end if
        end do
     end do
  end do

  ! --- assemble dielectric tensor ---

  do J=1,3
     do I=1,3
        do IN=1,3
           do NS=1,NSMAX
              DTENST(IN,I,J)=DTENST(IN,I,J)+DTENS(NS,IN,I,J)
           end do
        end do
     end do
  end do

  SELECT CASE(MODELG)
  CASE(0,12)
     do JJ=1,6
        do II=1,6
           do K=1,3
              do J=1,3
                 do I=1,3
                    CM2(II,JJ)= CM2(II,JJ)&
                                +((MU(I,1,II)*DTENST(J,1,1)&
                                  +MU(I,2,II)*DTENST(J,2,1)&
                                  +MU(I,3,II)*DTENST(J,3,1))&
                                  *MU(K,1,JJ)&
                                 +(MU(I,1,II)*DTENST(J,1,2)&
                                  +MU(I,2,II)*DTENST(J,2,2)&
                                  +MU(I,3,II)*DTENST(J,3,2))&
                                  *MU(K,2,JJ)&
                                 +(MU(I,1,II)*DTENST(J,1,3)&
                                  +MU(I,2,II)*DTENST(J,2,3)&
                                  +MU(I,3,II)*DTENST(J,3,3))&
                                  *MU(K,3,JJ))&
                                 *RR*S*AIF3(I,J,K)
                 end do
              end do
           end do
        end do
     end do
  CASE(1:10,13)
     do JJ=1,6
        do II=1,6
           do K=1,3
              do J=1,3
                 do I=1,3
                    CM2(II,JJ)= CM2(II,JJ)&
                                +((MU(I,1,II)*DTENST(J,1,1)&
                                  +MU(I,2,II)*DTENST(J,2,1)&
                                  +MU(I,3,II)*DTENST(J,3,1))&
                                  *MU(K,1,JJ)&
                                 +(MU(I,1,II)*DTENST(J,1,2)&
                                  +MU(I,2,II)*DTENST(J,2,2)&
                                  +MU(I,3,II)*DTENST(J,3,2))&
                                  *MU(K,2,JJ)&
                                 +(MU(I,1,II)*DTENST(J,1,3)&
                                  +MU(I,2,II)*DTENST(J,2,3)&
                                  +MU(I,3,II)*DTENST(J,3,3))&
                                  *MU(K,3,JJ))&
                                 *x(J)*S*AIF3(I,J,K)
                 end do
              end do
           end do
        end do
     end do
  END SELECT

  do J=1,6
     do I=1,6
        CM(I,J)=CM(I,J)-WC2*CM2(I,J)
     end do
  end do

!  write(*,*) 'NE=',NE
!  do I=1,6
!     write(6,'(6(A1,E9.3,A1,E9.3,A1))')&
!          "(",real(CM(I,1)),',',aimag(CM(I,1)),")",&
!          "(",real(CM(I,2)),',',aimag(CM(I,2)),")",&
!          "(",real(CM(I,3)),',',aimag(CM(I,3)),")",&
!          "(",real(CM(I,4)),',',aimag(CM(I,4)),")",&
!          "(",real(CM(I,5)),',',aimag(CM(I,5)),")",&
!          "(",real(CM(I,6)),',',aimag(CM(I,6)),")"
!  end do

  return
END SUBROUTINE cmcalcp

!     ****** DIELECTRIC TENSOR ******

  subroutine wf_dtensr(nelm,DTENS)

    use wfcomm
    USE wfprof
  implicit none
  integer,intent(in) :: nelm
  integer :: IN,I,J,ID
  complex(rkind),intent(out):: DTENS(NSM,3,3,3)
  integer    :: NS,node
  real(rkind)    :: x,y,WW,WP(NSM),WC(NSM),BABS,AL(3),RN(NSM),RTPR(NSM)
  real(rkind)    :: RTPP(NSM),RZCL(NSM),FP,FR,FZ,dx,dy,F
  complex(rkind) :: CWP,CWC,CDT0,CDX0,CDP0,CDT,CDP,CDX,CDAMP
  complex(rkind) :: CRR,CRP,CRZ,CPR,CPP,CPZ,CZR,CZP,CZZ

  ! ----- initialize -----  

  do J=1,3
     do I=1,3
        do IN=1,3
           do NS=1,NSMAX
              DTENS(NS,IN,I,J)=(0.d0,0.d0)
           end do
        end do
     end do
  end do

  WW=2.D0*PI*RF*1.D6

  DO NS=1,NSMAX
     WP(NS)=PZ(NS)*PZ(NS)*AEE*AEE*1.D20/(PA(NS)*AMP*EPS0*WW*WW)
     WC(NS)=PZ(NS)*AEE/(PA(NS)*AMP*WW)
  ENDDO

  ! ----- collisional cold plasma model -----

  do IN=1,3
     
     node=node_nside_nelm(IN,nelm)
     x=xnode(node)
     y=ynode(node)
     
     CALL wf_smag(x,y,BABS,AL)
     FR=AL(1)
     FP=AL(2)
     FZ=AL(3)
     
     CALL wf_sden(x,y,RN,RTPR,RTPP,RZCL)

     do NS=1,NSMAX
        
        CWP = WP(NS)*RN(NS)/(1.D0+CII*RZCL(NS))
        CWC = WC(NS)*BABS  /(1.D0+CII*RZCL(NS))
        CDT0= CWP/(1.D0-CWC**2)
        CDX0= CII*CWP*CWC/(1.D0-CWC**2)
        CDP0= CWP
        
        CDT=CDT0       
        CDP=CDP0-CDT0
        CDX=CDX0

        CRR= CDT   +CDP*FR*FR
        CRP= CDX*FZ+CDP*FR*FP
        CRZ=-CDX*FP+CDP*FR*FZ
        CPR=-CDX*FZ+CDP*FP*FR
        CPP= CDT   +CDP*FP*FP
        CPZ= CDX*FR+CDP*FP*FZ
        CZR= CDX*FP+CDP*FZ*FR
        CZP=-CDX*FR+CDP*FZ*FP
        CZZ= CDT   +CDP*FZ*FZ
        
        DTENS(NS,IN,1,1)=DTENS(NS,IN,1,1)-CRR
        DTENS(NS,IN,1,2)=DTENS(NS,IN,1,2)-CRP
        DTENS(NS,IN,1,3)=DTENS(NS,IN,1,3)-CRZ
        DTENS(NS,IN,2,1)=DTENS(NS,IN,2,1)-CPR
        DTENS(NS,IN,2,2)=DTENS(NS,IN,2,2)-CPP
        DTENS(NS,IN,2,3)=DTENS(NS,IN,2,3)-CPZ
        DTENS(NS,IN,3,1)=DTENS(NS,IN,3,1)-CZR
        DTENS(NS,IN,3,2)=DTENS(NS,IN,3,2)-CZP
        DTENS(NS,IN,3,3)=DTENS(NS,IN,3,3)-CZZ

     END do


     IF(width_damp.GT.0.D0) THEN
        CDAMP=CII*PZCL(NSMAX)
        F=factor_damp
        IF(x-xnode_min.LT.width_damp) THEN
           ID=1
           IF(model_damp.EQ.1.AND. &
              y.GT.ydamp_min.AND. &
              y.LT.ydamp_max) ID=0
           IF(ID.EQ.1) THEN
              dx=x-xnode_min
         !    DTENS(NSMAX,IN,1,1)=DTENS(NSMAX,IN,1,1) &
         !         +F*(width_damp-dx)/(dx-CDAMP)
              DTENS(NSMAX,IN,2,2)=DTENS(NSMAX,IN,2,2) &
                   +F*(width_damp-dx)/(dx-CDAMP)
              DTENS(NSMAX,IN,3,3)=DTENS(NSMAX,IN,3,3) &
                   +F*(width_damp-Dx)/(dx-CDAMP)
           END IF
        END IF
        IF(xnode_max-x.LT.width_damp) THEN
           ID=1
           IF(model_damp.EQ.2.AND. &
              y.GT.ydamp_min.AND. &
              y.LT.ydamp_max) ID=0
           IF(ID.EQ.1) THEN
              Dx=xnode_max-x
   !          DTENS(NSMAX,IN,1,1)=DTENS(NSMAX,IN,1,1) &
   !               +F*(width_damp-Dx)/(Dx-CDAMP)
              DTENS(NSMAX,IN,2,2)=DTENS(NSMAX,IN,2,2) &
                   +F*(width_damp-dx)/(dx-CDAMP)
              DTENS(NSMAX,IN,3,3)=DTENS(NSMAX,IN,3,3) &
                   +F*(width_damp-dx)/(dx-CDAMP)
           END IF
        END IF
        IF(y-ynode_min.LT.width_damp) THEN
           ID=1
           IF(model_damp.EQ.3.AND. &
              x.GT.xdamp_min.AND.&
              x.LT.xdamp_max) ID=0
           IF(ID.EQ.1) THEN
              dy=y-ynode_min
              DTENS(NSMAX,IN,1,1)=DTENS(NSMAX,IN,1,1) &
                   +F*(width_damp-dy)/(dy-CDAMP)
              DTENS(NSMAX,IN,2,2)=DTENS(NSMAX,IN,2,2) &
                   +F*(width_damp-dy)/(dy-CDAMP)
            ! DTENS(NSMAX,IN,3,3)=DTENS(NSMAX,IN,3,3) &
            !      +F*(width_damp-dy)/(dy-CDAMP)
           END IF
        END IF
        IF(ynode_max-y.LT.width_damp) THEN
           ID=1
           IF(model_damp.EQ.4.AND. &
              x.GT.xdamp_min.AND. &
              x.LT.xdamp_max) ID=0
           IF(ID.EQ.1) THEN
              dy=ynode_max-y
              DTENS(NSMAX,IN,1,1)=DTENS(NSMAX,IN,1,1) &
                   +F*(width_damp-dy)/(dy-CDAMP)
              DTENS(NSMAX,IN,2,2)=DTENS(NSMAX,IN,2,2) &
                   +F*(width_damp-dy)/(dy-CDAMP)
            ! DTENS(NSMAX,IN,3,3)=DTENS(NSMAX,IN,3,3) &
            !      +F*(width_damp-dy)/(dy-CDAMP)
           END IF
        END IF
     END IF
  end do

  return
  end subroutine wf_dtensr

!     ***** INTERPOLATION TENSOR *****

  subroutine wf_mutensr(nelm,MU)

    use wfcomm
    USE wfsub
  implicit none
  integer,intent(in) :: nelm
  integer :: nside,nseg,I,J,K
  real(rkind),intent(out)::MU(3,3,6)
  real(rkind) :: A(3),B(3),C(3),L(3)

  do nside=1,3
     nseg=ABS(nseg_nside_nelm(nseg,nelm))
     IF(model_wf.EQ.0) THEN
        L(nside)=len_nseg(nseg)
     ELSE
        IF(nseg_nside_nelm(nside,nelm).GT.0.D0) THEN
           L(nside)=len_nseg(nseg)
        ELSE
           L(nside)=-len_nseg(nseg)
        END IF
     END IF
  end do

  call wf_set_abc(nelm,A,B,C)

  do K=1,6
     do J=1,3
        do I=1,3
           MU(I,J,K)=0.d0
        end do
     end do
  end do

  MU(1,1,1)= L(1)*B(2)
  MU(1,3,1)= L(1)*C(2)
  MU(1,2,4)= 1.D0
  MU(1,1,3)=-L(3)*B(3)
  MU(1,3,3)=-L(3)*C(3)

  MU(2,1,1)=-L(1)*B(1)
  MU(2,3,1)=-L(1)*C(1)
  MU(2,2,5)= 1.D0
  MU(2,1,2)= L(2)*B(3)
  MU(2,3,2)= L(2)*C(3)

  MU(3,1,2)=-L(2)*B(2)
  MU(3,3,2)=-L(2)*C(2)
  MU(3,2,6)= 1.D0
  MU(3,1,3)= L(3)*B(1)
  MU(3,3,3)= L(3)*C(1)

  return
  end subroutine wf_mutensr
END MODULE wfsolv
