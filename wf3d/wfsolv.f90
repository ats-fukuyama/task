!     $Id$


!     ******* SET KANOD ARRAY *******

SUBROUTINE DEFBND(IERR)

  use wfcomm
  implicit none
  integer :: IERR,IL,NSD,KA,K,NB,MAX,NE,IN,NN

! ----- caluculate MLEN -------

! -- Count number of variables to describe a node --

  IERR=0
  IL=0
  DO NSD=1,NSDMAX
     KA=KASID(NSD)
     IF(KA.EQ.0) THEN
        K=1
     ELSEIF(KA.EQ.1) THEN
        K=0
     ELSEIF(KA.EQ.4) THEN
        K=1
     ELSEIF(KA.LT.0) THEN
        K=0
     ELSE
        WRITE(6,*) 'XX DEFBND: UNDEFINED KANOD=',KA
        IERR=2
     ENDIF
     IL=IL+K
  ENDDO

! -- WG boundary block requires a node --

  DO NB=1,NBMAX
     IF(KABDY(NB).GE.8) THEN
        K=NMBDY(NB)
        IL=IL+K
     ENDIF
  ENDDO

  MLEN=IL

  call wfslv_allocate

! ------------------------------------------------

!     Set number of variables to describe a node

  IERR=0
  IL=1
  DO NSD=1,NSDMAX
     KA=KASID(NSD)
     IF(KA.EQ.0) THEN
        K=1
     ELSEIF(KA.EQ.1) THEN
        K=0
     ELSEIF(KA.EQ.4) THEN
        K=1
     ELSEIF(KA.LT.0) THEN
        K=0
     ELSE
        WRITE(6,*) 'XX DEFBND: UNDEFINED KANOD=',KA
        IERR=2
     ENDIF
     INLEN(NSD)=K
     IMLEN(NSD)=IL
     IL=IL+K
!     WRITE(6,'(A,5I5)') 'NSD,KA,K,IL=',NSD,KA,K,IL
  ENDDO

!     WG boundary block requires a node

  NSD=NSDMAX
  NMDMAX=1
  DO NB=1,NBMAX
     IF(KABDY(NB).GE.8) THEN
        K=NMBDY(NB)
        NMDMAX=MAX(NMDMAX,K)
        NSD=NSD+1
        NDBDY(NB)=NSD
        INLEN(NSD)=K
        IMLEN(NSD)=IL
        IL=IL+K
     ENDIF
  ENDDO
  
  NNBMAX=NSD
  MLEN=IL-1
  IMLEN(NSD+1)=IL
  
!     Set NBELM to identify an element with WG boundary
!         NDELM for additional node for the element

  DO NE=1,NEMAX
     NBELM(NE)=0
  ENDDO
  DO NE=NEMAX,1,-1
     DO IN=1,4
        NN=NDELM(IN,NE)
        KA=KANOD(NN)
        IF(KA.LT.0) THEN
           NB=-KA
           IF(KABDY(NB).GT.8) THEN
              NSDELM(7,NE)=NDBDY(NB)
              IF(NBELM(NE).EQ.0) THEN
                 NBELM(NE)=NB
              ELSE
                 IF(NBELM(NE).NE.NB) THEN
                    WRITE(6,*) 'XX Element faces two boundary block'
                    WRITE(6,*) '   NE=',NE
                    IERR=3
                 ENDIF
              ENDIF
           ENDIF
        ENDIF
     ENDDO
  ENDDO
! DO NB=1,NBMAX
!     NN=NDBDY(NB)
!     WRITE(6,'(A,6I5)') 'NB,KABDY,NMBDY,NDBDY,INLEN,IMLEN=',&
!     &        NB,KABDY(NB),NMBDY(NB),NN,INLEN(NN),IMLEN(NN)
!  ENDDO
!  WRITE(6,*) 'MLEN,MLENM =',MLEN,MLENM
  RETURN
  
!9000 WRITE(6,*) 'XX DEFBND ERROR : MLEN,MLENM =',MLEN,MLENM
!  IERR=900
  RETURN
END SUBROUTINE DEFBND

! ****************************
SUBROUTINE CVSOLV(IERR)

  use wfcomm
  use libmtxc
  implicit none
  integer :: I,ISDMAX,ISD
  integer :: NSD,NE,NVMAX,J,KC,M
  integer :: K,L,II,JJ,LL,IE,NV,IV
  integer :: KK,IN,IERR,INSD,JNSD
  integer :: isize,itype
  integer :: its
  integer :: MLENP
  integer :: JMIN,JMAX,MILEN,MJLEN
  integer :: NNZME!Number of Non-Zero Matrix Element
  integer,dimension(:),pointer :: NEFLAG
  real(8) :: tolerance
  real(8),dimension(1) :: ddata
  real(4) :: cputime1,cputime2
  complex(8),dimension(:),pointer   :: CRVP
  complex(8),dimension(:,:),pointer :: CEQP

! ----- initialize for parallel computing -----

!  if (nrank.eq.0) then
!     write(*,*) "## INPUT: tolerance"
!     read (*,*) tolerance
!     write(*,'(A,1P,D16.4)') " tolerance =", tolerance
     ddata(1)=tolerance
!  end if
!  
!  call mtx_broadcast_real8(ddata,1)

  tolerance=EPSWF
  itype = 0

  call mtx_setup(MLEN,istart,iend,MLEN)

! ----- initialize ------

  CSV =(0.d0,0.d0)

  DO NE=1,NEMAX
     IF(NBELM(NE).EQ.0) THEN
        ISDMAX=6
     ELSE
        ISDMAX=7
     ENDIF
     DO ISD=1,ISDMAX
        ISDELM(ISD,NE)=ABS(NSDELM(ISD,NE))
     ENDDO
  ENDDO

! ------ create NV-NSD table -----

  NV=0
  do NSD=1,NSDMAX
     do IV=1,INLEN(NSD)
        NV=NV+1
        NSDNV(NV)=NSD
     end do
  end do

! ----- set NEFLAG ------

  allocate(NEFLAG(NEMAX))
  NEFLAG=0

  do NE=1,NEMAX
     IF(NBELM(NE).EQ.0) THEN
        ISDMAX=6
     ELSE
        ISDMAX=7
     ENDIF
     do ISD=1,ISDMAX
        NSD=ISDELM(ISD,NE)
        if(NSD.ge.NSDNV(istart).and.&
           NSD.le.NSDNV(iend) )then  
           NEFLAG(NE)=1
        end if
     end do
  end do

  ! --- decide MJLEN ---

  JMIN=MLEN
  JMAX=0
  do NE=1,NEMAX
     if(NEFLAG(NE).eq.0) goto 8100
     IF(NBELM(NE).EQ.0) THEN
        ISDMAX=6
     ELSE
        ISDMAX=7
     ENDIF
     do ISD=1,ISDMAX
        LL=0
        DO J=1,ISDMAX
           JNSD=ABS(ISDELM(J,NE))
           DO JJ=1,INLEN(JNSD)
              LL=IMLEN(JNSD)+JJ-1
              KK=0
              DO I=1,ISDMAX
                 INSD=ABS(ISDELM(I,NE))
                 DO II=1,INLEN(INSD)
                    KK=IMLEN(INSD)+II-1
                    if((KK.ge.istart).and.&
                       (KK.le.iend  )) then
                          if(LL.lt.JMIN) JMIN=LL
                          if(LL.gt.JMAX) JMAX=LL
                    end if
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     end do
8100 continue
  end do

  MJLEN=JMAX-JMIN+1
  MILEN=iend-istart+1

  allocate(CEQP(MILEN,MJLEN),CRVP(MILEN))
  CEQP=(0.d0,0.d0)
  CRVP=(0.d0,0.d0)

! ------ set grobal matrix -------

  NNZME=0
  do NE=1,NEMAX
     if(NEFLAG(NE).eq.0) goto 8000
     IF(NBELM(NE).EQ.0) THEN
        ISDMAX=6
     ELSE
        ISDMAX=7
     ENDIF
     CALL CMCALC(NE)     
     
     !    ----- ASSEMBLY -----
     !     --- if KK are out of assigned range,  
     !           CVSOLV do not save the matrix element. ---

     LL=0
     DO J=1,ISDMAX
        JNSD=ABS(ISDELM(J,NE))
        DO JJ=1,INLEN(JNSD)
           LL=IMLEN(JNSD)+JJ-1
           if(LL.gt.JMAX) goto 8200
           KK=0
           DO I=1,ISDMAX
              INSD=ABS(ISDELM(I,NE))
              DO II=1,INLEN(INSD)
                 KK=IMLEN(INSD)+II-1
                 if((KK.ge.istart).and.&
                    (KK.le.iend  )) then
                      CEQP(KK-istart+1,LL+1-JMIN)=CEQP(KK-istart+1,LL+1-JMIN)+CM(II,JJ,I,J)
                      if(CM(II,JJ,I,J).ne.(0.d0,0.d0)) NNZME=NNZME+1
                 end if
              ENDDO
           ENDDO
8200       continue
        ENDDO
     ENDDO
     
     KK=0
     DO I=1,ISDMAX
        INSD=ABS(ISDELM(I,NE))
        DO II=1,INLEN(INSD)
           KK=IMLEN(INSD)+II-1
           if((KK.ge.istart).and.&
              (KK.le.iend  )) then
                CRVP(KK-istart+1)=CRVP(KK-istart+1)+CV(II,I)
           end if
        ENDDO
     ENDDO
8000 continue
  end do

  if(nrank.eq.0) write(6,'(56A)') "  nrank   IMIN   IMAX  MILEN   JMIN   JMAX  MJLEN  NNZME"
  call mtx_barrier
  write(6,'(8I7)') nrank,istart,iend,iend-istart+1,JMIN,JMAX,jMAX-JMIN+1,NNZME

  do j=1,MJLEN
     do i=istart,iend
        if (abs(CEQP(i-istart+1,j)).ne.0.d0) &
             call mtx_set_matrix(i,JMIN-1+j,CEQP(i-istart+1,j))
     end do
  end do

  do i=istart,iend
     call mtx_set_source(i,CRVP(i-istart+1))
  end do

  call GUTIME(cputime1)
  
  call mtx_solve(itype,tolerance,its)
  !zmumps always return "its = 0"
  if(nrank.eq.0) write(6,*) 'Iteration Number=',its
  
  call GUTIME(cputime2)
  write(*,*) nrank,cputime2-cputime1

  call mtx_gather_vector(CSV)

  deallocate(CEQP,CRVP)
  deallocate(NEFLAG)
  RETURN
END SUBROUTINE CVSOLV
    
