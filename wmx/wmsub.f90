! wmsub.f90

MODULE wmsub

  PRIVATE
  PUBLIC wmsubc,wmsube,wmsube_f,wmsubg,wmsubg_f,wmsubf,wmsubf_f,wmxfft

CONTAINS
  
!     ****** 2D FOURIER TRANSFORM ******

  SUBROUTINE WMSUBC(CF1)

    USE wmcomm,ONLY: rkind,nthmax,nhhmax,ldsiz,kdsiz
    IMPLICIT NONE
    COMPLEX(rkind),INTENT(INOUT):: CF1(nthmax,nhhmax)
    COMPLEX(rkind),ALLOCATABLE:: CF2(:,:),CFM(:),CFN(:)
    INTEGER:: nth,nhh,ldx,kdx

    ALLOCATE(CF2(nthmax,nhhmax),CFM(nthmax),CFN(nhhmax))
    
    DO NHH=1,NHHMAX
       DO NTH=1,NTHMAX
          CFM(NTH)=CF1(NTH,NHH)
       ENDDO
       CALL WMXFFT(CFM,NTHMAX,0)
       DO LDX=1,LDSIZ
          CF2(LDX,NHH)=CFM(LDX)
       ENDDO
    ENDDO

    DO LDX=1,LDSIZ
       DO NHH=1,NHHMAX
          CFN(NHH)=CF2(LDX,NHH)
       ENDDO
       CALL WMXFFT(CFN,NHHMAX,0)
       DO KDX=1,KDSIZ
          CF1(LDX,KDX)=CFN(KDX)
       ENDDO
    ENDDO
    RETURN
  END SUBROUTINE WMSUBC

!     ****** 2D FOURIER TRANSFORM ******

  SUBROUTINE WMSUBE(CF1,CF2)

    USE wmcomm,ONLY: rkind,nthmax,nhhmax,mdsiz,ndsiz
    IMPLICIT NONE
    COMPLEX(rkind),INTENT(IN):: CF1(nthmax,nhhmax)
    COMPLEX(rkind),INTENT(OUT):: CF2(nthmax,nhhmax)
    COMPLEX(rkind),ALLOCATABLE:: CFM(:),CFN(:)
    INTEGER:: nth,nhh,mdx,ndx

    ALLOCATE(CFM(nthmax),CFN(nhhmax))

    DO NDX=1,NDSIZ
       DO MDX=1,MDSIZ
          CFM(MDX)=CF1(MDX,NDX)
       ENDDO
       CALL WMXFFT(CFM,NTHMAX,1)
       DO NTH=1,NTHMAX
          CF2(NTH,NDX)=CFM(NTH)
       ENDDO
    ENDDO

    DO NTH=1,NTHMAX
       DO NDX=1,NDSIZ
          CFN(NDX)=CF2(NTH,NDX)
       ENDDO
       CALL WMXFFT(CFN,NHHMAX,1)
       DO NHH=1,NHHMAX
          CF2(NTH,NHH)=CFN(NHH)
       ENDDO
    ENDDO
    RETURN
  END SUBROUTINE WMSUBE

!     ****** 2D FOURIER TRANSFORM ******

  SUBROUTINE WMSUBE_F(CF1,CF2)

    USE wmcomm,ONLY: rkind,nthmax_f,nhhmax_f,mdsiz_f,ndsiz_f
    IMPLICIT NONE
    COMPLEX(rkind),INTENT(IN):: CF1(nthmax_f,nhhmax_f)
    COMPLEX(rkind),INTENT(OUT):: CF2(nthmax_f,nhhmax_f)
    COMPLEX(rkind),ALLOCATABLE:: CFM(:),CFN(:)
    INTEGER:: nth,nhh,mdx,ndx

    ALLOCATE(CFM(nthmax_f),CFN(nhhmax_f))

    DO NDX=1,NDSIZ_F
       DO MDX=1,MDSIZ_F
          CFM(MDX)=CF1(MDX,NDX)
       ENDDO
       CALL WMXFFT(CFM,NTHMAX_F,1)
       DO NTH=1,NTHMAX_F
          CF2(NTH,NDX)=CFM(NTH)
       ENDDO
    ENDDO

    DO NTH=1,NTHMAX_F
       DO NDX=1,NDSIZ_F
          CFN(NDX)=CF2(NTH,NDX)
       ENDDO
       CALL WMXFFT(CFN,NHHMAX_F,1)
       DO NHH=1,NHHMAX_F
          CF2(NTH,NHH)=CFN(NHH)
       ENDDO
    ENDDO
    RETURN
  END SUBROUTINE WMSUBE_F

!     ****** 2D FOURIER TRANSFORM ******

  SUBROUTINE WMSUBG(RF1,RF2,CF)

    USE wmcomm
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: RF1(nthmax,nhhmax),RF2(nthmax,nhhmax)
    COMPLEX(rkind),INTENT(OUT):: CF(nthmax,nhhmax)
    COMPLEX(rkind):: CFM(nthmax),CFN(nhhmax)
    INTEGER:: NHH,NTH,LDX,KDX

    DO NHH=1,NHHMAX
       DO NTH=1,NTHMAX
          CFM(NTH)=RF1(NTH,NHH)/RF2(NTH,NHH)
       ENDDO
       CALL WMXFFT(CFM,NTHMAX,0)
       IF (LDSIZ == 1)THEN
          LDX=LDSIZ
          CF(LDX,NHH)=CFM(LDX)
       ELSE
          DO LDX=1,LDSIZ
             CF(LDX,NHH)=CFM(LDX)
          ENDDO
       ENDIF
    ENDDO

    DO LDX=1,LDSIZ
       DO NHH=1,NHHMAX
          CFN(NHH)=CF(LDX,NHH)
       ENDDO
       CALL WMXFFT(CFN,NHHMAX,0)
       IF (KDSIZ == 1)THEN
          KDX=KDSIZ
          CF(LDX,KDX)=CFN(KDX)
       ELSE
          DO KDX=1,KDSIZ
             CF(LDX,KDX)=CFN(KDX)
          ENDDO
       ENDIF
    ENDDO
    RETURN
  END SUBROUTINE WMSUBG

!     ****** 2D FOURIER TRANSFORM ******

  SUBROUTINE WMSUBG_F(RF1,RF2,CF)

    USE wmcomm
    IMPLICIT NONE
    REAL(rkind),INTENT(IN):: RF1(nthmax_f,nhhmax_f),RF2(nthmax_f,nhhmax_f)
    COMPLEX(rkind),INTENT(OUT):: CF(nthmax_f,nhhmax_f)
    COMPLEX(rkind):: CFM(nthmax_f),CFN(nhhmax_f)
    INTEGER:: NHH,NTH,LDX,KDX

    DO NHH=1,NHHMAX_F
       DO NTH=1,NTHMAX_F
          CFM(NTH)=RF1(NTH,NHH)/RF2(NTH,NHH)
       ENDDO
       CALL WMXFFT(CFM,NTHMAX_F,0)
       IF (LDSIZ_F == 1)THEN
          LDX=LDSIZ_F
          CF(LDX,NHH)=CFM(LDX)
       ELSE
          DO LDX=1,LDSIZ_F
             CF(LDX,NHH)=CFM(LDX)
          ENDDO
       ENDIF
    ENDDO

    DO LDX=1,LDSIZ_F
       DO NHH=1,NHHMAX_F
          CFN(NHH)=CF(LDX,NHH)
       ENDDO
       CALL WMXFFT(CFN,NHHMAX_F,0)
       IF (KDSIZ_F == 1)THEN
          KDX=KDSIZ_F
          CF(LDX,KDX)=CFN(KDX)
       ELSE
          DO KDX=1,KDSIZ_F
             CF(LDX,KDX)=CFN(KDX)
          ENDDO
       ENDIF
    ENDDO
    RETURN
  END SUBROUTINE WMSUBG_F

!     ****** 2D FOURIER TRANSFORM ******

  SUBROUTINE WMSUBF(CF1,CF2)

    USE wmcomm
    IMPLICIT NONE
    COMPLEX(rkind),INTENT(IN):: CF1(nthmax,nhhmax)
    COMPLEX(rkind),INTENT(OUT):: CF2(nthmax,nhhmax)
    COMPLEX(rkind):: CFM(nthmax),CFN(nhhmax)
    INTEGER:: NHH,NTH,LDX,KDX

    DO NHH=1,NHHMAX
       DO NTH=1,NTHMAX
          CFM(NTH)=CF1(NTH,NHH)
       ENDDO
       CALL WMXFFT(CFM,NTHMAX,0)
       IF (LDSIZ == 1)THEN
          LDX=LDSIZ
          CF2(LDX,NHH)=CFM(LDX)
       ELSE
          DO LDX=1,LDSIZ
             CF2(LDX,NHH)=CFM(LDX)
          ENDDO
       ENDIF
    ENDDO

    DO LDX=1,LDSIZ
       DO NHH=1,NHHMAX
          CFN(NHH)=CF2(LDX,NHH)
       ENDDO
       CALL WMXFFT(CFN,NHHMAX,0)
       IF (KDSIZ == 1)THEN
          KDX=KDSIZ
          CF2(LDX,KDX)=CFN(KDX)
       ELSE
          DO KDX=1,KDSIZ
             CF2(LDX,KDX)=CFN(KDX)
          ENDDO
       ENDIF
    ENDDO
    RETURN
  END SUBROUTINE WMSUBF

!     ****** 2D FOURIER TRANSFORM 2 ******

  SUBROUTINE WMSUBF_F(CF1,CF2)

    USE wmcomm
    IMPLICIT NONE
    COMPLEX(rkind),INTENT(IN):: CF1(nthmax_f,nhhmax_f)
    COMPLEX(rkind),INTENT(OUT):: CF2(nthmax_f,nhhmax_f)
    COMPLEX(rkind):: CFM(nthmax_f),CFN(nhhmax_f)
    INTEGER:: NHH,NTH,LDX,KDX

    DO NHH=1,NHHMAX_F
       DO NTH=1,NTHMAX_F
          CFM(NTH)=CF1(NTH,NHH)
       ENDDO
       CALL WMXFFT(CFM,NTHMAX_F,0)
       IF (LDSIZ_F == 1)THEN
          LDX=LDSIZ_F
          CF2(LDX,NHH)=CFM(LDX)
       ELSE
          DO LDX=1,LDSIZ_F
             CF2(LDX,NHH)=CFM(LDX)
          ENDDO
       ENDIF
    ENDDO

    DO LDX=1,LDSIZ_F
       DO NHH=1,NHHMAX_F
          CFN(NHH)=CF2(LDX,NHH)
       ENDDO
       CALL WMXFFT(CFN,NHHMAX_F,0)
       IF (KDSIZ_F == 1)THEN
          KDX=KDSIZ_F
          CF2(LDX,KDX)=CFN(KDX)
       ELSE
          DO KDX=1,KDSIZ_F
             CF2(LDX,KDX)=CFN(KDX)
          ENDDO
       ENDIF
    ENDDO
    RETURN
  END SUBROUTINE WMSUBF_F

!     ****** INTERFACE FOR FFT ******

  SUBROUTINE WMXFFT(CA,N,KEY)

    USE wmcomm
    USE libfft,ONLY: FFT2L
    IMPLICIT NONE
    INTEGER,INTENT(IN):: N,KEY
    COMPLEX(rkind),INTENT(INOUT):: CA(N)
    INTEGER,SAVE:: NS=0
    INTEGER:: IND,I,IX
!
    IF(N.NE.1) THEN
       IF(N.EQ.NS) THEN
          IND=0
       ELSE
          IF(ALLOCATED(CFFT)) DEALLOCATE(CFFT,RFFT,LFFT)
          ALLOCATE(CFFT(N*2),RFFT(N*2),LFFT(N*2))
          IND=1
          NS=N
       END IF
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
  END SUBROUTINE WMXFFT
END MODULE wmsub
