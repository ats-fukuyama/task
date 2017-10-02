!     $Id: fpwmin.f90,v 1.4 2013/01/14 16:48:26 fukuyama Exp $
!
!     ***** READ AND INTERPOLATE wmdata *****
!
      MODULE fpwmin

      use fpcomm

      integer:: NWTHMAX,NWPHMAX,NWRMAX,NPH0W,NTH0W
      real(rkind):: RAW,RRW,BBW,RFWR,RFWI,RABS
      complex(rkind),dimension(:,:,:,:),POINTER:: CEWV 
                                                  ! (3,NWTHM,NWPHM,NWRM)
      real(rkind),dimension(:),POINTER ::  RWSPL  ! (NWRM)
      real(rkind),dimension(:),POINTER ::  THWSPL ! (NWTHM)
      real(rkind),dimension(:),POINTER ::  PHWSPL ! (NWPHM)
      complex(rkind),dimension(:),POINTER :: CEWL ! (NWPHM)
      complex(rkind),dimension(:,:),POINTER :: CEWL2,CEWX2,CEWY2,CEWXY2
                                                  ! (NWTHMP,NWRM)
      complex(rkind),dimension(:,:,:),POINTER :: CEWL3,CEWX3,CEWY3,CEWZ3, &
                                                 CEWXY3,CEWYZ3,CEWZX3,CEWXYZ3
                                                  ! (NWTHMP,NWPHMP,NWRM)
      complex(rkind),dimension(:,:,:,:,:),POINTER :: UCEW2 
                                                  ! (4,4,NWTHM,NWRM,3)
      complex(rkind),dimension(:,:,:,:,:,:),POINTER :: UCEW3
                                                  ! (4,4,NWTHM,NEPHM,NWRM,3)

      complex(rkind),dimension(:),POINTER:: CFFT  ! (NWTHM)
      real(rkind),dimension(:),POINTER:: RFFT     ! (NWTHM)
      integer,dimension(:),POINTER:: LFFT         ! (NWTHM) 

      contains

!--------------------------------------------

      SUBROUTINE fp_wm_read(IERR)

      USE libmpi
      USE libmtx
      IMPLICIT NONE
      integer:: IERR,NSMAX1,I,NWR,NWTH,NWPH

      IF(nrank.EQ.0) THEN
         CALL FROPEN(21,KNAMWM,0,MODEFR,'WM',IERR)
         IF(IERR.NE.0) THEN
            WRITE(6,*) 'XX FPWMREAD: FROPEN: IERR=',IERR
            RETURN
         ENDIF

         REWIND(21)
         READ(21) NWTHMAX,NWPHMAX,NWRMAX,NSMAX1
         WRITE(6,*) 'NWTHMAX,NWPHMAX,NWRMAX=',NWTHMAX,NWPHMAX,NWRMAX

         IF(ASSOCIATED(CEWV)) deallocate(CEWV)
         allocate(CEWV(3,NWTHMAX,NWPHMAX,NWRMAX))

         READ(21) RAW,RRW,BBW,RFWR,RFWI,NPH0W,NTH0W

         READ(21) ((((CEWV(I,NWTH,NWPH,NWR),I=1,3),NWTH=1,NWTHMAX), &
                      NWPH=1,NWPHMAX), NWR=1,NWRMAX)
         CLOSE(21)
      ENDIF

      CALL fp_wm_broadcast

      IF(ASSOCIATED(RWSPL)) deallocate(RWSPL)
      IF(ASSOCIATED(THWSPL)) deallocate(THWSPL)
      IF(ASSOCIATED(PHWSPL)) deallocate(PHWSPL)
      IF(ASSOCIATED(CEWL)) deallocate(CEWL)
      IF(ASSOCIATED(UCEW2)) deallocate(UCEW2)
      IF(ASSOCIATED(CEWL2)) deallocate(CEWL2)
      IF(ASSOCIATED(CEWX2)) deallocate(CEWX2)
      IF(ASSOCIATED(CEWY2)) deallocate(CEWY2)
      IF(ASSOCIATED(CEWXY2)) deallocate(CEWXY2)
      IF(ASSOCIATED(UCEW3)) deallocate(UCEW3)
      IF(ASSOCIATED(CEWL3)) deallocate(CEWL3)
      IF(ASSOCIATED(CEWX3)) deallocate(CEWX3)
      IF(ASSOCIATED(CEWY3)) deallocate(CEWY3)
      IF(ASSOCIATED(CEWZ3)) deallocate(CEWZ3)
      IF(ASSOCIATED(CEWXY3)) deallocate(CEWXY3)
      IF(ASSOCIATED(CEWYZ3)) deallocate(CEWYZ3)
      IF(ASSOCIATED(CEWZX3)) deallocate(CEWZX3)
      IF(ASSOCIATED(CEWXYZ3)) deallocate(CEWXYZ3)

      allocate(RWSPL(NWRMAX))
      allocate(THWSPL(NWTHMAX+1))
      allocate(PHWSPL(NWPHMAX+1))
      allocate(CEWL(NWTHMAX))
      IF(NWPHMAX.EQ.1) THEN
         allocate(UCEW2(4,4,NWTHMAX+1,NWRMAX,3))
         allocate(CEWL2(NWTHMAX+1,NWRMAX))
         allocate(CEWX2(NWTHMAX+1,NWRMAX))
         allocate(CEWY2(NWTHMAX+1,NWRMAX))
         allocate(CEWXY2(NWTHMAX+1,NWRMAX))
      ELSE
         allocate(UCEW3(4,4,NWTHMAX+1,NWPHMAX+1,NWRMAX,3))
         allocate(CEWL3(NWTHMAX+1,NWPHMAX+1,NWRMAX))
         allocate(CEWX3(NWTHMAX+1,NWPHMAX+1,NWRMAX))
         allocate(CEWY3(NWTHMAX+1,NWPHMAX+1,NWRMAX))
         allocate(CEWZ3(NWTHMAX+1,NWPHMAX+1,NWRMAX))
         allocate(CEWXY3(NWTHMAX+1,NWPHMAX+1,NWRMAX))
         allocate(CEWYZ3(NWTHMAX+1,NWPHMAX+1,NWRMAX))
         allocate(CEWZX3(NWTHMAX+1,NWPHMAX+1,NWRMAX))
         allocate(CEWXYZ3(NWTHMAX+1,NWPHMAX+1,NWRMAX))
      ENDIF

      IF(NWPHMAX.EQ.1) THEN
         DO NWR=1,NWRMAX
            DO I=1,3
               DO NWTH=1,NWTHMAX
                  CEWL(NWTH)=FACT_WM*CEWV(I,NWTH,1,NWR)
               ENDDO
               CALL FPFFT(CEWL,NWTHMAX,1)
               DO NWTH=1,NWTHMAX
                  CEWV(I,NWTH,1,NWR)=CEWL(NWTH)
               ENDDO
            ENDDO
         ENDDO
      ELSE
         WRITE(6,*) 'XX fpwmin: NWPHMAX.NE.1 not supprted yet.'
         STOP
      ENDIF

      DO NWR=1,NWRMAX
         RWSPL(NWR)=DBLE(NWR-1)/DBLE(NWRMAX-1)
      ENDDO
      DO NWTH=1,NWTHMAX+1
         THWSPL(NWTH)=DBLE(NWTH-1)*2.D0*PI/DBLE(NWTHMAX)
      ENDDO
      DO NWPH=1,NWPHMAX+1
         PHWSPL(NWPH)=DBLE(NWPH-1)*2.D0*PI/DBLE(NWPHMAX)
      ENDDO

      IF(NWPHMAX.EQ.1) THEN
         NWPH=1
         DO I=1,3
            DO NWR=1,NWRMAX
               DO NWTH=1,NWTHMAX
                  CEWL2(NWTH,NWR)=CEWV(I,NWTH,NWPH,NWR)
               ENDDO
                  CEWL2(NWTHMAX+1,NWR)=CEWV(I,1,NWPH,NWR)
            ENDDO
            CALL CSPL2D(THWSPL,RWSPL,CEWL2,CEWX2,CEWY2,CEWXY2, &
                        UCEW2(1,1,1,1,I), &
                        NWTHMAX+1,NWTHMAX+1,NWRMAX,4,0,IERR)
            IF(IERR.NE.0) THEN
               WRITE(6,*) 'XX FPWMREAD: CSPL2D: IERR=',IERR
               IERR=201
               RETURN
            ENDIF
         ENDDO
      ELSE
         DO I=1,3
            DO NWR=1,NWRMAX
               DO NWTH=1,NWTHMAX
                  DO NWPH=1,NWPHMAX
                     CEWL3(NWTH,NWPH,NWR)=CEWV(I,NWTH,NWPH,NWR)
                  ENDDO
                  CEWL3(NWTH,NWPHMAX+1,NWR)=CEWV(I,NWTH,1,NWR)
               ENDDO
               DO NWPH=1,NWPHMAX
                  CEWL3(NWTHMAX+1,NWPH,NWR)=CEWV(I,1,NWPH,NWR)
               ENDDO
               CEWL3(NWTHMAX+1,NWPHMAX+1,NWR)=CEWV(I,1,1,NWR)
            ENDDO
            CALL CSPL3D(THWSPL,PHWSPL,RWSPL,CEWL3, &
                        CEWX3,CEWY3,CEWZ3,CEWXY3,CEWYZ3,CEWZX3,CEWXYZ3, &
                        UCEW3(1,1,1,1,1,I), &
                        NWTHMAX+1,NWPHMAX+1,NWTHMAX+1,NWPHMAX+1,NWRMAX, &
                        4,4,0,IERR)
            IF(IERR.NE.0) THEN
               WRITE(6,*) 'XX FPWMREAD: CSPL2D: IERR=',IERR
               IERR=201
               RETURN
            ENDIF
         ENDDO
      ENDIF

      RETURN
      END SUBROUTINE fp_wm_read

!     ***** BROADCAST WM DATA *****

      SUBROUTINE fp_wm_broadcast

      USE libmpi
      USE libmtx
      USE fpcomm
      IMPLICIT NONE
      INTEGER,DIMENSION(5):: idata
      REAL(8),DIMENSION(5):: ddata
      COMPLEX(8),DIMENSION(:),POINTER:: temp
      INTEGER:: nr1,md1,nd1,n

      IF(nrank.eq.0) THEN
         idata(1)=NWTHMAX
         idata(2)=NWPHMAX
         idata(3)=NWRMAX
         idata(4)=NPH0W
         idata(5)=NTH0W
         ddata(1)=RAW
         ddata(2)=RRW
         ddata(3)=BBW
         ddata(4)=RFWR
         ddata(5)=RFWI
      ENDIF
      CALL mtx_broadcast_integer(idata,5)
      CALL mtx_broadcast_real8(ddata,5)
      NWTHMAX =idata(1)
      NWPHMAX =idata(2)
      NWRMAX  =idata(3)
      NPH0W   =idata(4)
      NTH0W   =idata(5)
      RAW  =ddata(1)
      RRW  =ddata(2)
      BBW  =ddata(3)
      RFWR =ddata(4)
      RFWI =ddata(5)

      IF(nrank.NE.0) THEN
         IF(ASSOCIATED(CEWV)) DEALLOCATE(CEWV)
         ALLOCATE(CEWV(3,NWTHMAX,NWPHMAX,NWRMAX))
      ENDIF

      ALLOCATE(TEMP(3*NWTHMAX*NWPHMAX*NWRMAX))

      IF(nrank.eq.0) THEN
         DO NR1=1,NWRMAX
            DO ND1=1,NWPHMAX
               DO MD1=1,NWTHMAX
                  N=3*(MD1-1)+3*NWTHMAX*(ND1-1)+3*NWTHMAX*NWPHMAX*(NR1-1)
                  TEMP(N+1)=CEWV(1,MD1,ND1,NR1)
                  TEMP(N+2)=CEWV(2,MD1,ND1,NR1)
                  TEMP(N+3)=CEWV(3,MD1,ND1,NR1)
               ENDDO
            ENDDO
         ENDDO
      ENDIF
      CALL mtx_broadcast_complex8(TEMP,3*NWTHMAX*NWPHMAX*NWRMAX)
      IF(nrank.ne.0) THEN
         DO NR1=1,NWRMAX
            DO ND1=1,NWPHMAX
               DO MD1=1,NWTHMAX
                  N=3*(MD1-1)+3*NWTHMAX*(ND1-1)+3*NWTHMAX*NWPHMAX*(NR1-1)
                  CEWV(1,MD1,ND1,NR1)=TEMP(N+1)
                  CEWV(2,MD1,ND1,NR1)=TEMP(N+2)
                  CEWV(3,MD1,ND1,NR1)=TEMP(N+3)
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      DEALLOCATE(TEMP)

      END SUBROUTINE fp_wm_broadcast

!     ***** check INTERPOLATE wmdata *****

      SUBROUTINE FPWMCHEK

      IMPLICIT NONE
      integer:: IERR
      COMPLEX(8):: CEWR1,CEWTH1,CEWPH1,CKWR1,CKWTH1,CKWPH1
      real(8):: RL, THL, PHL, RFWM
      DATA RL,THL,PHL/0.D0,0.D0,0.D0/

 1010 CONTINUE
         WRITE(6,'(A,1P2E12.4)') 'Input RL,THL,PHL:',RL,THL,PHL
         READ(5,*,END=9000,ERR=1010) RL,THL,PHL
         IF(RL.EQ.0.D0) GOTO 9000
         CALL FPWMGET(RL,THL,PHL,RFWM,CEWR1,CEWTH1,CEWPH1, &
                                      CKWR1,CKWTH1,CKWPH1,IERR)
         WRITE(6,'(2I10)') NPH0W,NTH0W
         WRITE(6,'(1P2E12.4)') RRW,RAW
         WRITE(6,'(1P6E12.4)') CEWR1,CEWTH1,CEWPH1,CKWR1,CKWTH1,CKWPH1
         GOTO 1010

 9000    RETURN
         END SUBROUTINE FPWMCHEK
!
!     ***** INTERPOLATE wmdata *****
!
      SUBROUTINE FPWMGET(RL,THL,PHL,RFWM,CEWR1,CEWTH1,CEWPH1, &
                                         CKWR1,CKWTH1,CKWPH1,IERR)

      IMPLICIT NONE
      REAL(8),INTENT(IN):: RL,THL,PHL
      REAL(8),INTENT(OUT):: RFWM
      COMPLEX(8),INTENT(OUT):: CEWR1,CEWTH1,CEWPH1,CKWR1,CKWTH1,CKWPH1
      INTEGER,INTENT(OUT):: IERR
      COMPLEX(8):: CEWDTH,CEWDPH,CEWDR

      IERR=0
      RFWM=RFWR

      IF(NWPHMAX.EQ.1) THEN
         CALL CSPL2DD(THL,RL,CEWR1,CEWDTH,CEWDR,THWSPL,RWSPL, &
                      UCEW2(1,1,1,1,1),NWTHMAX+1,NWTHMAX+1,NWRMAX,IERR)
      ELSE
         CALL CSPL3DD(THL,PHL,RL,CEWR1,CEWDTH,CEWDPH,CEWDR, &
                      THWSPL,PHWSPL,RWSPL,UCEW3(1,1,1,1,1,1), &
                      NWTHMAX+1,NWPHMAX+1,NWTHMAX+1,NWPHMAX+1,NWRMAX,IERR)
      ENDIF
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX FPWMGET: 1: CSPL2DD: IERR=',IERR
         IERR=1
         RETURN
      ENDIF
!      write(*,'(1P8E10.2)') THL, RL, CEWR1, CEWDTH, CEWDR
      CKWR1=-CI*CEWDR/(CEWR1*RAW)

      IF(NWPHMAX.EQ.1) THEN
         CALL CSPL2DD(THL,RL,CEWTH1,CEWDTH,CEWDR,THWSPL,RWSPL, &
                      UCEW2(1,1,1,1,2),NWTHMAX+1,NWTHMAX+1,NWRMAX,IERR)
      ELSE
         CALL CSPL3DD(THL,PHL,RL,CEWTH1,CEWDTH,CEWDPH,CEWDR, &
                      THWSPL,PHWSPL,RWSPL,UCEW3(1,1,1,1,1,2), &
                      NWTHMAX+1,NWPHMAX+1,NWTHMAX+1,NWPHMAX+1,NWRMAX,IERR)
      ENDIF
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX FPWMGET: 2: CSPL2DD: IERR=',IERR
         IERR=2
         RETURN
      ENDIF
!      write(*,'(1P8E10.2)') THL, RL, CEWTH1, CEWDTH, CEWDR

      IF(RL.LE.0.D0) THEN
         CKWTH1=0.D0
      ELSE
         CKWTH1=-CI*CEWDTH/(CEWTH1*RL*RAW)
      ENDIF

      IF(NWPHMAX.EQ.1) THEN
         CALL CSPL2DD(THL,RL,CEWPH1,CEWDTH,CEWDR,THWSPL,RWSPL, &
                      UCEW2(1,1,1,1,3),NWTHMAX+1,NWTHMAX+1,NWRMAX,IERR)
         CKWPH1=NPH0W/RRW
      ELSE
         CALL CSPL3DD(THL,PHL,RL,CEWPH1,CEWDTH,CEWDPH,CEWDR, &
                      THWSPL,PHWSPL,RWSPL,UCEW3(1,1,1,1,1,3), &
                      NWTHMAX+1,NWPHMAX+1,NWTHMAX+1,NWPHMAX+1,NWRMAX,IERR)
         CKWPH1=-CI*CEWDPH/(CEWPH1*RRW)
      ENDIF
      IF(IERR.NE.0) THEN
         WRITE(6,*) 'XX FPWMGET: 3: CSPL2DD: IERR=',IERR
         RETURN
      ENDIF

      RETURN
      END SUBROUTINE FPWMGET
!
!     ****** INTERFACE FOR FFT ******
!
      SUBROUTINE FPFFT(CA,N,KEY)

      USE libfft,ONLY: FFT2L
      IMPLICIT NONE
      integer:: N, KEY, IND, IX, I
      integer,save:: NS=0
      COMPLEX(8),DIMENSION(N):: CA
!      complex(8),dimension(N):: CFFT ! (NWTHM) 
!      real(8),dimension(N):: RFFT ! (NWTHM) 
!      integer,dimension(N):: LFFT ! (NWTHM)

      IF(N.NE.1) THEN
         IF(N.EQ.NS) THEN
            IND=0
         ELSE
            IF(NS.NE.0) deallocate(CFFT,RFFT,LFFT)
            allocate(CFFT(N),RFFT(N),LFFT(N))               
            IND=1
            NS=N
         ENDIF
         IF(KEY.EQ.0) THEN
            CALL FFT2L(CA,CFFT,RFFT,LFFT,N,IND,KEY)
            DO I=1,N
               IX=I+N/2-1
               IF(IX.GT.N) IX=IX-N
               CA(IX)=CFFT(I)
            ENDDO
         ELSE
            DO I=1,N
               IX=I+N/2-1
               IF(IX.GT.N) IX=IX-N
               CFFT(I)=CA(IX)
            ENDDO
            CALL FFT2L(CFFT,CA,RFFT,LFFT,N,IND,KEY)
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE FPFFT

!--------------------------

      end module fpwmin
