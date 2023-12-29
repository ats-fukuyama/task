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
  integer :: NSD,NN

  nseg_bdy_max=0
  do NSD=1,nseg_max
     if(KASID(NSD).eq.1) nseg_bdy_max=nseg_bdy_max+1
  end do

  node_bdy_max=0
  do NN=1,node_max
     if(KANOD(NN).eq.1) node_bdy_max=node_bdy_max+1
  end do

  mtx_len=nseg_max+node_max-nseg_bdy_max-node_bdy_max

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
  integer :: ISD,NSD,nnd,nnd1,nnd2
  integer :: NE,NN,nv,nvmax
  integer :: I,J,KK,LL
  integer :: JNSD,JNN,INSD,INN
  integer :: IN,INV,JNV,KB
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
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: VAL_SORT
  INTEGER,DIMENSION(:),ALLOCATABLE:: NV_SORT
  INTEGER,DIMENSION(:),ALLOCATABLE:: ntyp_nv,nnsd_nv
  INTEGER:: IX,IY,nv_old

  ! ----- initialize ------
  
  do NV=1,mtx_len
     CSV(NV) =(0.d0,0.d0)
  end do

  ! --- decide istart,iend ---

  IF(nrank.EQ.0) write(6,*) 'mtx_len=',mtx_len
  call mtxc_setup(mtx_len,istart,iend,0)
  call mtxc_cleanup

  ! ----- set NV_NSD & NV_NN -----

  ALLOCATE(ntyp_nv(node_max+nseg_max),nnsd_nv(node_max+nseg_max))

  NV=0
  do NSD=1,nseg_max
     if(KASID(NSD).eq.1) then
        NVNSD(NSD)=0
     else
        NV=NV+1
        NVNSD(NSD)=NV
        NTYP_NV(NV)=2
        NNSD_NV(NV)=NSD
     end if
  end do
  do NN=1,node_max
     if(KANOD(NN).eq.1) then
        NVNN(NN)=0
     else
        NV=NV+1
        NVNN(NN)=NV
        NTYP_NV(NV)=1
        NNSD_NV(NV)=NN
     end if
  end do
  NVMAX=NV

  ALLOCATE(NV_SORT(NVMAX),VAL_SORT(NVMAX))

  DO NV=1,NVMAX
     IF(NTYP_NV(NV).EQ.1) THEN
        NN=NNSD_NV(NV)
        X=xnode(NN)
        Y=ynode(NN)
        VAL=sort_weight_x*X+sort_weight_y*Y
     ELSE
        NSD=NNSD_NV(NV)
        NND1=node_nseg(1,NSD)
        NND2=node_nseg(2,NSD)
        X=0.5D0*(xnode(NND1)+xnode(NND2))
        Y=0.5D0*(ynode(NND1)+ynode(NND2))
     END IF
     VAL=sort_weight_x*X+sort_weight_y*Y
     nv_sort(NV)=NV
     val_sort(NV)=VAL
  END DO

!  CALL qsort_dl(val_sort,nv_sort)
  CALL qsort_di(val_sort,nv_sort)

  DO nv=1,nvmax
     nv_old=nv_sort(nv)
     IF(ntyp_nv(nv_old).EQ.1) THEN
        nnd=nnsd_nv(nv_old)
        nvnn(nnd)=nv
!        X=xnode(nnd)
!        Y=ynode(nnd)
!        VAL=sort_weight_x+sort_weight_y*Y
!        write(21,'(I10,A,I10,1P3E12.4)') nv,' node ',nnd,X,Y,VAL
     ELSE
        nsd=nnsd_nv(nv_old)
        nvnsd(nsd)=nv
!        NND1=node_nseg(1,NSD)
!        NND2=node_nseg(2,NSD)
!        X=0.5D0*(xnode(NND1)+xnode(NND2))
!        Y=0.5D0*(ynode(NND1)+ynode(NND2))
!        VAL=sort_weight_x*X+sort_weight_y*Y
!        write(21,'(I10,A,I10,1P3E12.4)') nv,' side ',nsd,X,Y,VAL
     END IF
  END DO
     
  DEALLOCATE(NV_SORT,VAL_SORT,ntyp_nv,nnsd_nv)

!  do NSD=1,nseg_max
!     write(*,*) NSD,KASID(NSD),NV_NSD(NSD)
!  end do
!  do NN=1,node_max
!     write(*,*) NN,KANOD(NN),NV_NN(NN)
!  end do

  ! ----- set NEFLAG ------

  allocate(NEFLAG(nelm_max))
  do NE=1,nelm_max
     NEFLAG(NE)=0
  end do

  do NE=1,nelm_max
     do ISD=1,3
        NSD=ABS(nseg_nside_nelm(ISD,NE))
        NV=NVNSD(NSD)
        if(NV.ge.istart.and.&
           NV.le.iend ) then
           NEFLAG(NE)=1
        end if
     end do

     if(NEFLAG(NE).eq.0) then
        do IN=1,3
           NN=node_nside_nelm(IN,NE)
           NV=NVNN(NN)
           if(NV.ge.istart.and.&
              NV.le.iend) then
              NEFLAG(NE)=1
           end if
        end do
     end if

  end do

  ! ----- set MJLEN -----

  JMIN=mtx_len
  JMAX=0
  do NE=1,nelm_max

     if(NEFLAG(NE).eq.0) goto 8100

     LL=0
     DO J=1,6
        if(J.ge.1.and.J.le.3) then
           JNSD=ABS(nseg_nside_nelm(J,NE))
           JNV =NVNSD(JNSD)
        else
           JNN=node_nside_nelm(J-3,NE)
           JNV=NVNN(JNN)
        end if
        if (JNV.eq.0) goto 8110
        LL=JNV

        KK=0
        DO I=1,6
           if(I.ge.1.and.I.le.3) then
              INSD=ABS(nseg_nside_nelm(I,NE))
              INV=NVNSD(INSD)
           else
              INN=node_nside_nelm(I-3,NE)
              INV=NVNN(INN)
           end if
           if(INV.eq.0) goto 8120
           KK=INV

           if((KK.ge.istart).and.&
              (KK.le.iend  )) then
              if(LL.lt.JMIN) JMIN=LL
              if(LL.gt.JMAX) JMAX=LL
           end if

8120       continue
        ENDDO
8110    continue
     ENDDO
     
8100 continue
  end do

! ------ Count non-zero component -------

  NNZ=0

  DO NE=1,nelm_max
     IF(NEFLAG(NE).NE.0) THEN

        LL=0
        DO J=1,6
           IF(J.ge.1.and.J.le.3) then
              JNSD=nseg_nside_nelm(J,NE)
              if(JNSD.lt.0) then
                 JNSD=-JNSD
              end if
              JNV =NVNSD(JNSD)
           else
              JNN=node_nside_nelm(J-3,NE)
              JNV=NVNN(JNN)
           END IF
           LL=JNV

           if ((LL.GE.JMIN).AND.(LL.LE.JMAX)) THEN

              KK=0
              DO I=1,6
                 if(I.ge.1.and.I.le.3) then
                    INSD=nseg_nside_nelm(I,NE)
                    if(INSD.lt.0) then
                       INSD=-INSD
                    end if
                    INV =NVNSD(INSD)
                 else
                    INN=node_nside_nelm(I-3,NE)
                    INV=NVNN(INN)
                 end if
                 KK=INV

                 if((KK.ge.istart).and.&
                      (KK.le.iend  )) then
                    NNZ=NNZ+1
                 end if
              END DO
           END if
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

  do NE=1,nelm_max
     if(NEFLAG(NE).eq.0) goto 8000
     CALL wf_cmcalc(NE)
     
!    === ASSEMBLY ===
!    If KK (or LL) is out of assigned range,  
!      CVSOLV do not save the matrix element. 

!    --- inside of the boundary ---
     LL=0
     DO J=1,6
        ORIENTJ=1
        if(J.ge.1.and.J.le.3) then
           JNSD=nseg_nside_nelm(J,NE)
           if(JNSD.lt.0) then
              JNSD=-JNSD
              ORIENTJ=-1
           end if
           JNV =NVNSD(JNSD)
        else
           JNN=node_nside_nelm(J-3,NE)
           JNV=NVNN(JNN)
        end if
        LL=JNV

        if ((LL.GE.JMIN).AND.(LL.LE.JMAX)) THEN

           KK=0
           DO I=1,6
              ORIENTI=1
              if(I.ge.1.and.I.le.3) then
                 INSD=nseg_nside_nelm(I,NE)
                 if(INSD.lt.0) then
                    INSD=-INSD
                    ORIENTI=-1
                 end if
                 INV =NVNSD(INSD)
              else
                 INN=node_nside_nelm(I-3,NE)
                 INV=NVNN(INN)
              end if
              KK=INV
              if((KK.ge.istart).and.&
                   (KK.le.iend  )) then
                 if(abs(CM(I,J)).ne.0.d0) THEN 
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
        if(J.ge.1.and.J.le.3) then
           JNSD=nseg_nside_nelm(J,NE)
           if(JNSD.lt.0) then
              JNSD=-JNSD
              ORIENTJ=-1
           end if
           KB =KBSID(JNSD)
           IF(KB.NE.0) CEB=CEBSD(KB)
        else
           JNN=node_nside_nelm(J-3,NE)
           KB=KBNOD(JNN)
           IF(KB.NE.0) CEB=CEBND(KB)
        end if

        IF(KB.NE.0) THEN
           IF(ABS(CEB).GT.0.D0) THEN
              KK=0
              DO I=1,6
                 ORIENTI=1
                 if(I.ge.1.and.I.le.3) then
                    INSD=nseg_nside_nelm(I,NE)
                    if(INSD.lt.0) then
                       INSD=-INSD
                       ORIENTI=-1
                    end if
                    INV =NVNSD(INSD)
                 else
                    INN=node_nside_nelm(I-3,NE)
                    INV=NVNN(INN)
                 end if
                 KK=INV
                 if((KK.ge.istart).and.&
                      (KK.le.iend  )) then
                    CRVP(KK-istart+1)=CRVP(KK-istart+1) &
                         -ORIENTI*ORIENTJ*CM(I,J)*CEB
                 end if
              END DO
           END IF
        END IF
     ENDDO

!    --- Contribution from the antenna current ---

     KK=0
     DO I=1,6
        ORIENTI=1
        if(I.ge.1.and.I.le.3) then
           INSD=nseg_nside_nelm(I,NE)
           if(INSD.lt.0) then
              INSD=-INSD
              ORIENTI=-1
           end if
           INV =NVNSD(INSD)
        else
           INN=node_nside_nelm(I-3,NE)
           INV=NVNN(INN)
        end if
        KK=INV
        if((KK.ge.istart).and.&
           (KK.le.iend  )) then
           CRVP(KK-istart+1)=CRVP(KK-istart+1)+ORIENTI*CVTOT(I,NE)
        end if
     ENDDO

8000 continue
  end do

  IF(nrank.EQ.0) write(6,*) 'wfsolv: sort started'
  CALL qsort_ic(NSEQ,CEQP)
  !  CALL qsort_lc(NSEQ,CEQP)
  
  IF(nrank.EQ.0) write(6,*) 'wfsolv: reduction started'
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

  if(nrank.eq.0) write(6,'(A77)') &
  '      nrank     istart       iend      MILEN      MJLEN     NNZMAX      NNZME'
  call mtx_barrier
  write(6,'(7I11)') nrank,istart,iend,iend-istart+1, &
                    JMAX-JMIN+1,NNZMAX,NNZME

  ! ----- initialize for parallel computing -----

  itype = 0
  call mtxc_setup(mtx_len,istart,iend,nzmax=NNZME)

  do NNZ=1,NNZME
     IF(ABS(CEQP(NNZ)).GT.0.D0) THEN
        i=NSEQ(NNZ)/MJLEN
        j=NSEQ(NNZ)-i*MJLEN
        if(i.lt.0.OR.i.GT.iend-istart) then
           WRITE(6,'(A/6I12)') 'NNZ,NSEQ(NNZ),i,istart,i+istart,iend=', &
                               NNZ,NSEQ(NNZ),i,istart,i+istart,iend
           STOP
        END if
        call mtxc_set_matrix(i+istart,j+jmin,CEQP(NNZ))
     END IF
  end do

  do i=istart,iend
     call mtxc_set_source(i,CRVP(i-istart+1))
  end do

  call GUTIME(cputime1)

  call mtxc_solve(itype,tolerance,its)
  !zmumps always return "its = 0"
  if(nrank.eq.0) write(6,*) 'Iteration Number=',its

  call GUTIME(cputime2)

  call mtxc_gather_vector(CSV)

  deallocate(CEQP,NSEQ,CRVP)
  deallocate(NEFLAG)
  call mtxc_cleanup
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

  SUBROUTINE cmcalcv(ne)

    use wfcomm
    USE wfsub
  implicit none
  integer,intent(in) :: NE 
  integer :: I,J,K,M,N,ISD,NSD
  complex(rkind) :: CM1(3,3)
  real(rkind) :: S,L(3)
  real(rkind) :: A(3),B(3),C(3),AW(3),BW(3),CW(3)
  real(rkind) :: R(3),Z(3)

  ! --- initialize ---

  S=SELM(NE)
  call wf_set_abc(NE,A,B,C)
  call wf_set_node(NE,R,Z)

  do ISD=1,3
     NSD=ABS(nseg_nside_nelm(ISD,NE))
     IF(MODELWF.EQ.0) THEN
        L(ISD)=LSID(NSD)
     ELSE
        IF(nseg_nside_nelm(ISD,NE).GT.0.D0) THEN
           L(ISD)=LSID(NSD)
        ELSE
           L(ISD)=-LSID(NSD)
        END IF
     END IF
  end do

  do ISD=1,3
     M=ISD
     N=ISD+1
     if(N.gt.3) N=N-3
     AW(ISD)=L(ISD)*(A(M)*B(N)-A(N)*B(M))
     BW(ISD)=L(ISD)*(B(M)*C(N)-B(N)*C(M))
     CW(ISD)=L(ISD)*(C(M)*A(N)-C(N)*A(M))
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
                       *( AW(I)*AW(J)-(AW(I)*BW(J)+AW(J)*BW(I))*Z(K) &
                         +CW(I)*CW(J)-(BW(I)*CW(J)+BW(J)*CW(I))*R(K) &
                         +BW(I)*BW(J)*(R(K)**2+Z(K)**2)) &
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
                      +(real(NPH)**2)/R(K) &
                       *( AW(I)*AW(J)-(AW(I)*BW(J)+AW(J)*BW(I))*Z(K) &
                         +CW(I)*CW(J)-(BW(I)*CW(J)+BW(J)*CW(I))*R(K) &
                         +BW(I)*BW(J)*(R(K)**2+Z(K)**2)) &
                       *S*AIF1(K) &
                      +4.d0*BW(I)*BW(J)*R(K)*S*AIF1(K)
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
                          *(AW(J)-BW(J)*Z(K)) &
                         +C(I) &
                          *(CW(J)-BW(J)*R(K))) &
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
                          *(AW(J)-BW(J)*Z(K)) &
                         +C(I) &
                          *(CW(J)-BW(J)*R(K))) &
                       *S*AIF1(K)&
!
                      +(CII*real(NPH))&
                       *(-(AW(J)-BW(J)*Z(K))/R(K))&
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
                          *(AW(I)-BW(I)*Z(K)) &
                         +C(J)&
                          *(CW(I)-BW(I)*R(K))) &
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
                          *(AW(I)-BW(I)*Z(K)) &
                         +C(J)&
                          *(CW(I)-BW(I)*R(K))) &
                        *S*AIF1(K) &
!
                      -(CII*real(NPH)) &
                       *(-(AW(I)-BW(I)*Z(K))/R(K)) &
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
                      +(B(I)*B(J)+C(I)*C(J))*R(K)*S*AIF1(K) &
                      +B(J)*S*AIF2(I,K) &
                      +B(I)*S*AIF2(J,K) &
                      +1.D0/(R(K))*S*AIF3(I,J,K)
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

  SUBROUTINE cmcalcs(ne)

    use wfcomm
    USE wfsub
  implicit none
  integer,intent(in) :: NE 
  integer :: I,J,K,I1,J1,K1,ID,ISD,NSD,IND(2),ND1,ND2
  real(rkind) :: S,L(3),RN,ZN
  real(rkind) :: A(3),B(3),C(3),AW(3),BW(3),CW(3)
  real(rkind) :: R(3),Z(3)
  COMPLEX(rkind):: CTEMP

  ! --- if no boundary side, return ---

  ID=0
  DO ISD=1,3
     NSD=ABS(nseg_nside_nelm(ISD,NE))
     IF(KASID(NSD).EQ.1) ID=1
  END DO
  IF(ID.EQ.0) RETURN

  ! --- initialize ---

  S=SELM(NE)
  call wf_set_abc(NE,A,B,C)
  call wf_set_node(NE,R,Z)

  DO ISD=1,3
     NSD=ABS(nseg_nside_nelm(ISD,NE))
     IF(KASID(NSD).EQ.1) THEN
        IF(MODELWF.EQ.0) THEN
           L(ISD)=LSID(NSD)
        ELSE
           IF(nseg_nside_nelm(ISD,NE).GT.0.D0) THEN
              L(ISD)=LSID(NSD)
           ELSE
              L(ISD)=-LSID(NSD)
           END IF
        END IF

        ND1=ISD
        ND2=ISD+1
        IF(ND2.GT.3) ND2=ND2-3
        AW(ISD)=L(ISD)*(A(ND1)*B(ND2)-A(ND2)*B(ND1))
        BW(ISD)=L(ISD)*(B(ND1)*C(ND2)-B(ND2)*C(ND1))
        CW(ISD)=L(ISD)*(C(ND1)*A(ND2)-C(ND2)*A(ND1))
     END IF
  END DO

  DO ISD=1,3
     NSD=ABS(nseg_nside_nelm(ISD,NE))
     IF(KASID(NSD).EQ.1) THEN

        IND(1)=ISD
        IND(2)=ISD+1
        IF(IND(2).GT.3) IND(2)=IND(2)-3
        RN= (Z(IND(2))-Z(IND(1)))/L(ISD)
        ZN=-(R(IND(2))-R(IND(1)))/L(ISD)

  ! ----- rotErotF term -----

  ! --- E1F1 ---

        I=ISD
        J=ISD
        SELECT CASE(MODELG)
        CASE(0,12)
           DO K1=1,2
              K=IND(K1)
              CTEMP=CM(I,J)
              CM(I,J)=CM(I,J) &
                      -2.D0*BW(J)*(ZN*(AW(I)-BW(I)*Z(K)) &
                                  +RN*(CW(I)-BW(I)*R(K)))*RR*L(ISD)*AIE1(K)
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
                      -2.D0*BW(J)*(ZN*(AW(I)-BW(I)*Z(K)) &
                                  +RN*(CW(I)-BW(I)*R(K)))*R(K)*L(ISD)*AIE1(K)
           END DO
        END SELECT

  ! --- E2F1 --- 

        I=ISD
        DO J1=1,2
           J=IND(J1)
           SELECT CASE(MODELG)
           CASE(0,12)
              DO K1=1,2
                 K=IND(K1)
                 CTEMP=CM(I,J+3)
                 CM(I,J+3)=CM(I,J+3) &
                          +CII*RKZ*RR &
                           *(-ZN*(-CW(J)+BW(J)*R(K)) &
                             -RN*( AW(J)-BW(J)*Z(K)))*L(ISD)*AIE2(J,K)
                    WRITE(21,'(A,I10,3I5,1P4E12.4)') &
                         'CM:',NE,I,J+3,K,CM(I,J+3),CM(I,J+3)-CTEMP
              END DO
           CASE(1:10,13)
              DO K1=1,2
                 K=IND(K1)
                 CM(I,J+3)=CM(I,J+3) &
                          +CII*NPH &
                           *(-ZN*(-CW(J)+BW(J)*R(K)) &
                             -RN*( AW(J)-BW(J)*Z(K)))*L(ISD)*AIE2(J,K)
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
                               +(RN*B(J)+ZN*C(J))*RR &
                                *L(ISD)*AIE3(I,J,K)
                    WRITE(21,'(A,I10,3I5,1P4E12.4)') &
                         'CM:',NE,I+3,J+3,K,CM(I+3,J+3),CM(I+3,J+3)-CTEMP
                 END DO
              CASE(1:10,13)
                 DO K1=1,2
                    K=IND(K1)
                    CM(I+3,J+3)=CM(I+3,J+3) &
                               +((RN*B(J)+ZN*C(J))*R(K) &
                                 +RN*(A(J)+B(J)*R(K)+C(J)*Z(K))) &
                                *L(ISD)*AIE3(I,J,K)
                 END DO
              END SELECT
           END DO
        END DO
     END IF
  END DO

  RETURN
  END SUBROUTINE cmcalcs

!     ****** LOCAL ELEMENT MATRIX : Surface integral ******

  SUBROUTINE cmcalcp(ne)

    use wfcomm
    USE wfsub
  implicit none
  integer,intent(in) :: NE 
  integer :: I,J,K,NS,IN,JJ,II
  real(rkind) :: RW,WC,WC2
  real(rkind) :: S
  real(rkind) :: A(3),B(3),C(3)
  real(rkind) :: R(3),Z(3),MU(3,3,6)
  COMPLEX(rkind):: CM2(6,6)
  complex(rkind) :: DTENS(NSM,3,3,3)
  complex(rkind) :: DTENST(3,3,3)

  ! --- initialize ---

  RW=2.D0*PI*RF*1.D6
  WC=RW/VC
  WC2=WC**2

  S=SELM(NE)  
  call wf_set_abc(NE,A,B,C)
  call wf_set_node(NE,R,Z)

  ! ----- dielectric tensor term -----

  call wf_dtensr(NE,DTENS)
  call wf_mutensr(NE,MU)

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
                                 *R(J)*S*AIF3(I,J,K)
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

  subroutine wf_dtensr(NE,DTENS)

    use wfcomm
    USE wfprof
  implicit none
  integer,intent(in) :: NE
  integer :: IN,I,J,ID
  complex(rkind),intent(out):: DTENS(NSM,3,3,3)
  integer    :: NS,NN
  real(rkind)    :: R,Z,WW,WP(NSM),WC(NSM),BABS,AL(3),RN(NSM),RTPR(NSM)
  real(rkind)    :: RTPP(NSM),RZCL(NSM),FP,FR,FZ,DR,DZ,F
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
!     WRITE(6,'(I8,1P2E12.4)') NS,PA(NS),PZ(NS)
     WP(NS)=PZ(NS)*PZ(NS)*AEE*AEE*1.D20/(PA(NS)*AMP*EPS0*WW*WW)
     WC(NS)=PZ(NS)*AEE/(PA(NS)*AMP*WW)
  ENDDO

  ! ----- collisional cold plasma model -----

  do IN=1,3
     
     NN=node_nside_nelm(IN,NE)
     R=xnode(NN)
     Z=ynode(NN)
     
     CALL wf_smag(R,Z,BABS,AL)
     FR=AL(1)
     FP=AL(2)
     FZ=AL(3)
     
     CALL wf_sden(R,Z,RN,RTPR,RTPP,RZCL)

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


     IF(WDAMP.GT.0.D0) THEN
        CDAMP=CII*PZCL(NSMAX)
        F=FDAMP
        IF(R-BDRMIN.LT.WDAMP) THEN
           ID=1
           IF(MDAMP.EQ.1.AND. &
              Z.GT.ZDAMP_MIN.AND.Z.LT.ZDAMP_MAX) ID=0
           IF(ID.EQ.1) THEN
              DR=R-BDRMIN
!              DTENS(NSMAX,IN,1,1)=DTENS(NSMAX,IN,1,1)+F*(WDAMP-DR)/(DR-CDAMP)
              DTENS(NSMAX,IN,2,2)=DTENS(NSMAX,IN,2,2)+F*(WDAMP-DR)/(DR-CDAMP)
              DTENS(NSMAX,IN,3,3)=DTENS(NSMAX,IN,3,3)+F*(WDAMP-DR)/(DR-CDAMP)
           END IF
        END IF
        IF(BDRMAX-R.LT.WDAMP) THEN
           ID=1
           IF(MDAMP.EQ.2.AND. &
              Z.GT.ZDAMP_MIN.AND.Z.LT.ZDAMP_MAX) ID=0
           IF(ID.EQ.1) THEN
              DR=BDRMAX-R
!              DTENS(NSMAX,IN,1,1)=DTENS(NSMAX,IN,1,1)+F*(WDAMP-DR)/(DR-CDAMP)
              DTENS(NSMAX,IN,2,2)=DTENS(NSMAX,IN,2,2)+F*(WDAMP-DR)/(DR-CDAMP)
              DTENS(NSMAX,IN,3,3)=DTENS(NSMAX,IN,3,3)+F*(WDAMP-DR)/(DR-CDAMP)
           END IF
        END IF
        IF(Z-BDZMIN.LT.WDAMP) THEN
           ID=1
           IF(MDAMP.EQ.3.AND. &
              R.GT.RDAMP_MIN.AND.R.LT.RDAMP_MAX) ID=0
           IF(ID.EQ.1) THEN
              DZ=Z-BDZMIN
              DTENS(NSMAX,IN,1,1)=DTENS(NSMAX,IN,1,1)+F*(WDAMP-DZ)/(DZ-CDAMP)
              DTENS(NSMAX,IN,2,2)=DTENS(NSMAX,IN,2,2)+F*(WDAMP-DZ)/(DZ-CDAMP)
!              DTENS(NSMAX,IN,3,3)=DTENS(NSMAX,IN,3,3)+F*(WDAMP-DZ)/(DZ-CDAMP)
           END IF
        END IF
        IF(BDZMAX-Z.LT.WDAMP) THEN
           ID=1
           IF(MDAMP.EQ.4.AND. &
              R.GT.RDAMP_MIN.AND.R.LT.RDAMP_MAX) ID=0
           IF(ID.EQ.1) THEN
              DZ=BDZMAX-Z
              DTENS(NSMAX,IN,1,1)=DTENS(NSMAX,IN,1,1)+F*(WDAMP-DZ)/(DZ-CDAMP)
              DTENS(NSMAX,IN,2,2)=DTENS(NSMAX,IN,2,2)+F*(WDAMP-DZ)/(DZ-CDAMP)
!              DTENS(NSMAX,IN,3,3)=DTENS(NSMAX,IN,3,3)+F*(WDAMP-DZ)/(DZ-CDAMP)
           END IF
        END IF
     END IF
  end do

  return
  end subroutine wf_dtensr

!     ***** INTERPOLATION TENSOR *****

  subroutine wf_mutensr(NE,MU)

    use wfcomm
    USE wfsub
  implicit none
  integer,intent(in) :: NE
  integer :: ISD,NSD,I,J,K
  real(rkind),intent(out)::MU(3,3,6)
  real(rkind) :: A(3),B(3),C(3),L(3)

  do ISD=1,3
     NSD=ABS(nseg_nside_nelm(ISD,NE))
     IF(MODELWF.EQ.0) THEN
        L(ISD)=LSID(NSD)
     ELSE
        IF(nseg_nside_nelm(ISD,NE).GT.0.D0) THEN
           L(ISD)=LSID(NSD)
        ELSE
           L(ISD)=-LSID(NSD)
        END IF
     END IF
  end do

  call wf_set_abc(NE,A,B,C)

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
