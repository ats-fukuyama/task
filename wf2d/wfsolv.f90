!     $Id: wfsolv.f90,v 1.19 2012/03/05 06:29:02 maruyama Exp $

!    ******* SET MLEN *******

SUBROUTINE DEFMLEN

  use wfcomm
  implicit none
  integer :: NSD,NN

  NBSID=0
  do NSD=1,NSDMAX
     if(KASID(NSD).eq.1) NBSID=NBSID+1
  end do

  NBNOD=0
  do NN=1,node_max
     if(KANOD(NN).eq.1) NBNOD=NBNOD+1
  end do

  MLEN=NSDMAX+node_max-NBSID-NBNOD

  call wfslv_allocate

  RETURN
END SUBROUTINE DEFMLEN

!     ****** SOLV MATRIX EQUATION *****

SUBROUTINE CVSOLV

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
  integer(long),dimension(:),ALLOCATABLE :: NSEQ
  REAL(rkind),DIMENSION(:),ALLOCATABLE:: VAL_SORT
  INTEGER(long),DIMENSION(:),ALLOCATABLE:: NV_SORT
  INTEGER,DIMENSION(:),ALLOCATABLE:: ntyp_nv,nnsd_nv
  INTEGER(long):: IX,IY,nv_old

  ! ----- initialize ------
  
  do NV=1,MLEN
     CSV(NV) =(0.d0,0.d0)
  end do

  ! --- decide istart,iend ---

  IF(nrank.EQ.0) write(6,*) 'MLEN=',MLEN
  call mtxc_setup(MLEN,istart,iend,0)
  call mtxc_cleanup

  ! ----- set NV_NSD & NV_NN -----

  ALLOCATE(ntyp_nv(node_max+nsdmax),nnsd_nv(node_max+nsdmax))

  NV=0
  do NSD=1,NSDMAX
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
        X=RNODE(NN)
        Y=ZNODE(NN)
        VAL=sort_weight_x*X+sort_weight_y*Y
     ELSE
        NSD=NNSD_NV(NV)
        NND1=NDSID(1,NSD)
        NND2=NDSID(2,NSD)
        X=0.5D0*(RNODE(NND1)+RNODE(NND2))
        Y=0.5D0*(ZNODE(NND1)+ZNODE(NND2))
     END IF
     VAL=sort_weight_x*X+sort_weight_y*Y
     nv_sort(NV)=NV
     val_sort(NV)=VAL
  END DO

  CALL qsort_dl(val_sort,nv_sort)

  DO nv=1,nvmax
     nv_old=nv_sort(nv)
     IF(ntyp_nv(nv_old).EQ.1) THEN
        nnd=nnsd_nv(nv_old)
        nvnn(nnd)=nv
!        X=RNODE(nnd)
!        Y=ZNODE(nnd)
!        VAL=sort_weight_x+sort_weight_y*Y
!        write(21,'(I10,A,I10,1P3E12.4)') nv,' node ',nnd,X,Y,VAL
     ELSE
        nsd=nnsd_nv(nv_old)
        nvnsd(nsd)=nv
!        NND1=NDSID(1,NSD)
!        NND2=NDSID(2,NSD)
!        X=0.5D0*(RNODE(NND1)+RNODE(NND2))
!        Y=0.5D0*(ZNODE(NND1)+ZNODE(NND2))
!        VAL=sort_weight_x*X+sort_weight_y*Y
!        write(21,'(I10,A,I10,1P3E12.4)') nv,' side ',nsd,X,Y,VAL
     END IF
  END DO
     
  DEALLOCATE(NV_SORT,VAL_SORT,ntyp_nv,nnsd_nv)

!  do NSD=1,NSDMAX
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
        NSD=ABS(NSDELM(ISD,NE))
        NV=NVNSD(NSD)
        if(NV.ge.istart.and.&
           NV.le.iend ) then
           NEFLAG(NE)=1
        end if
     end do

     if(NEFLAG(NE).eq.0) then
        do IN=1,3
           NN=NDELM(IN,NE)
           NV=NVNN(NN)
           if(NV.ge.istart.and.&
              NV.le.iend) then
              NEFLAG(NE)=1
           end if
        end do
     end if

  end do

  ! ----- set MJLEN -----

  JMIN=MLEN
  JMAX=0
  do NE=1,nelm_max

     if(NEFLAG(NE).eq.0) goto 8100

     LL=0
     DO J=1,6
        if(J.ge.1.and.J.le.3) then
           JNSD=ABS(NSDELM(J,NE))
           JNV =NVNSD(JNSD)
        else
           JNN=NDELM(J-3,NE)
           JNV=NVNN(JNN)
        end if
        if (JNV.eq.0) goto 8110
        LL=JNV

        KK=0
        DO I=1,6
           if(I.ge.1.and.I.le.3) then
              INSD=ABS(NSDELM(I,NE))
              INV=NVNSD(INSD)
           else
              INN=NDELM(I-3,NE)
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
              JNSD=NSDELM(J,NE)
              if(JNSD.lt.0) then
                 JNSD=-JNSD
              end if
              JNV =NVNSD(JNSD)
           else
              JNN=NDELM(J-3,NE)
              JNV=NVNN(JNN)
           END IF
           LL=JNV

           if ((LL.GE.JMIN).AND.(LL.LE.JMAX)) THEN

              KK=0
              DO I=1,6
                 if(I.ge.1.and.I.le.3) then
                    INSD=NSDELM(I,NE)
                    if(INSD.lt.0) then
                       INSD=-INSD
                    end if
                    INV =NVNSD(INSD)
                 else
                    INN=NDELM(I-3,NE)
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
     CALL CMCALC(NE)
     
!    === ASSEMBLY ===
!    If KK (or LL) is out of assigned range,  
!      CVSOLV do not save the matrix element. 

!    --- inside of the boundary ---
     LL=0
     DO J=1,6
        ORIENTJ=1
        if(J.ge.1.and.J.le.3) then
           JNSD=NSDELM(J,NE)
           if(JNSD.lt.0) then
              JNSD=-JNSD
              ORIENTJ=-1
           end if
           JNV =NVNSD(JNSD)
        else
           JNN=NDELM(J-3,NE)
           JNV=NVNN(JNN)
        end if
        LL=JNV

        if ((LL.GE.JMIN).AND.(LL.LE.JMAX)) THEN

           KK=0
           DO I=1,6
              ORIENTI=1
              if(I.ge.1.and.I.le.3) then
                 INSD=NSDELM(I,NE)
                 if(INSD.lt.0) then
                    INSD=-INSD
                    ORIENTI=-1
                 end if
                 INV =NVNSD(INSD)
              else
                 INN=NDELM(I-3,NE)
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
           JNSD=NSDELM(J,NE)
           if(JNSD.lt.0) then
              JNSD=-JNSD
              ORIENTJ=-1
           end if
           KB =KBSID(JNSD)
           IF(KB.NE.0) CEB=CEBSD(KB)
        else
           JNN=NDELM(J-3,NE)
           KB=KBNOD(JNN)
           IF(KB.NE.0) CEB=CEBND(KB)
        end if

        IF(KB.NE.0) THEN
           IF(ABS(CEB).GT.0.D0) THEN
              KK=0
              DO I=1,6
                 ORIENTI=1
                 if(I.ge.1.and.I.le.3) then
                    INSD=NSDELM(I,NE)
                    if(INSD.lt.0) then
                       INSD=-INSD
                       ORIENTI=-1
                    end if
                    INV =NVNSD(INSD)
                 else
                    INN=NDELM(I-3,NE)
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
           INSD=NSDELM(I,NE)
           if(INSD.lt.0) then
              INSD=-INSD
              ORIENTI=-1
           end if
           INV =NVNSD(INSD)
        else
           INN=NDELM(I-3,NE)
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
  CALL qsort_lc(NSEQ,CEQP)
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
  call mtxc_setup(MLEN,istart,iend,nzmax=NNZME)

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
END SUBROUTINE CVSOLV
