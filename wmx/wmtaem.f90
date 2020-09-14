! wmtaem.f90

MODULE wmtaem

  PRIVATE
  PUBLIC wm_taem

CONTAINS

!   ***********************************************
!   **     Analysis of Alfven Eigenmode          **
!   ***********************************************

  SUBROUTINE wm_taem

    USE wmcomm
    USE plprof
    USE pllocal
    USE wmprof
    USE wmsub
    IMPLICIT NONE
    COMPLEX(rkind),ALLOCATABLE:: CWALF(:,:),CWALFK(:,:)
    REAL(rkind),ALLOCATABLE:: RNX(:),RTPRX(:),RTPPX(:),RUX(:)
    REAL(rkind),ALLOCATABLE:: DM(:),EM(:),FM(:,:),QM(:,:),W(:),WORK(:)
    REAL(4),ALLOCATABLE:: GX(:),GY(:,:)
    INTEGER,ALLOCATABLE:: NGY(:),IBLOCK(:),ISPLIT(:),IWORK(:)
    INTEGER:: MMM,INMAX,NRMAXP,MMMAX,MHMAX,NR,NS,NHH,NTH
    INTEGER:: ND,NDX,NN,MD,MDX,MM,NW,NW0,ND1,NDX1,NN1,MD1,MDX1,MM1,NW1,NY
    INTEGER:: IMIN,IMAX,NSPLIT,INFO1,INFO2
    REAL(rkind):: RHOM,RIOTAL,BSUPTH,BSUPPH,VALF,RHON,RKPR1,RKPR2
    REAL(rkind):: EPSEG,VU,VL


    MMM=nthmax*nhhmax
    ALLOCATE(CWALF(nthmax,nhhmax),CWALFK(nthmax,nhhmax))
    ALLOCATE(NGY(nrmax+1))
    ALLOCATE(RNX(nsmax),RTPRX(nsmax),RTPPX(nsmax),RUX(nsmax))
    ALLOCATE(DM(MMM),EM(MMM),FM(MMM,MMM),QM(1,MMM),W(MMM))
    ALLOCATE(IBLOCK(MMM),ISPLIT(MMM))
    ALLOCATE(GX(nrmax+1),GY(nrmax+1,MMM))
    ALLOCATE(WORK(5*MMM),IWORK(3*MMM))

    INMAX=nrmax+1
    NRMAXP=NRMAX

    MMMAX=MAX(1,NDSIZ-1)*MAX(1,MDSIZ-1)
    MHMAX=MAX(1,NDSIZ/2)*MAX(1,MDSIZ/2)-1

    DO NR=1,NRMAX
       CALL WMCDEN(NR,RNX,RTPRX,RTPPX,RUX)
       RHON=XRHO(NR)
       CALL PL_PROF_OLD(RHON)
       RHOM=0.D0
       DO NS=2,NSMAX
          IF(MODELP(NS).LT.0) THEN
             RHOM=RHOM+RNX(NS)*1.D20*PA(NS)*AMP
          ELSE
             RHOM=RHOM+RN(NS)*1.D20*PA(NS)*AMP
          ENDIF
       ENDDO
       IF(RHOM.LE.0.D0) THEN
          NRMAXP=NR-1
          GOTO 1000
       ENDIF

!         RIOTA_AV=0.D0
!         BABS_AV=0.D0
!         DO NHH=1,NHHMAX
!            DO NTH=1,NTHMAX
!               CALL WMCMAG(NR,NTH,NHH,BABS,BSUPTH,BSUPPH)
!               RIOTA_AV=RIOTA_AV+BSUPPH/BSUPTH
!               BABS_AV=BABS_AV+BABS
!            ENDDO
!         ENDDO
!         RIOTA_AV=RIOTA_AV/(NTHMAX*NHHMAX)
!         BABS_AV=BABS_AV/(NTHMAX*NHHMAX)

!         WRITE(6,'(I5,1P3E12.4)') NR,BABS_AV,RIOTA_AV,QPS(NR)

       RIOTAL=QPS(NR)

       DO NHH=1,NHHMAX
          DO NTH=1,NTHMAX
             CALL WMCMAG(NR,NTH,NHH,BABS,BSUPTH,BSUPPH)
             VALF=BABS/SQRT(RMU0*RHOM)
             CWALF(NTH,NHH)=RPS(NTH,NR)**2/VALF**2
          ENDDO
       ENDDO

       CALL wmsubf(CWALF,CWALFK)

       DO ND=NDMIN,MAX(NDMIN,NDMAX-1)
          NDX=ND-NDMIN+1
          NN=NPH0+NHC*ND
          DO MD=MDMIN,MAX(MDMIN,MDMAX-1)
             MDX=MD-MDMIN+1
             MM=NTH0+MD
               
             RKPR1=MM+NN*RIOTAL

             NW=MAX(1,MDSIZ-1)*(NDX-1)+(MDX-1)+1
             NW0=MAX(1,MDSIZ-1)*(-NDMIN)+(-MDMIN)+1
             DO ND1=NDMIN,MAX(NDMIN,NDMAX-1)
                NDX1=ND1-NDMIN+1
                NN1=NN+NHC*ND1
                DO MD1=MDMIN,MAX(MDMIN,MDMAX-1)
                   MDX1=MD1-MDMIN+1
                   MM1=MM+MD1

                   RKPR2=MM1+NN1*RIOTAL

                   NW1=MAX(1,MDSIZ-1)*(NDX1-1)+(MDX1-1)+1
                   NW1=NW1-NW0+1
                   IF(NW1.GE.1.AND.NW1.LE.MHMAX+1) THEN
                      FM(NW1,NW)=DBLE(CWALFK(MDX1,NDX1))/(RKPR1*RKPR2) &
                           *RIOTAL**2
                   ENDIF
                ENDDO
             ENDDO

!            DM(NWW)=DBLE(CWALFK(0-MDMIN+1,0-NDMIN+1))
!            EM(NWW)=DBLE(CWALFK(1-MDMIN+1,0-NDMIN+1))
          ENDDO
       ENDDO

       CALL LAPACK_DSBTRD('N','L',MMMAX,MHMAX, &
                          FM,MMM,DM,EM,QM,1,WORK,INFO1)

       EPSEG=0.D0
       VU=1.D0/(2.D0*PI*1.D6*WAEMIN)**2
       VL=1.D0/(2.D0*PI*1.D6*WAEMAX)**2
!         VU=(2.D0*PI*1.D6*OME)**2
!         VL=(2.D0*PI*1.D6*OMS)**2

       CALL LAPACK_DSTEBZ('V','B',MMMAX,VL,VU,IMIN,IMAX,EPSEG, &
                          DM,EM,INMAX,NSPLIT,W, &
                          IBLOCK,ISPLIT,WORK,IWORK,INFO2)
!         WRITE(6,*) 'NR,INFO1,INFO2,INMAX=',NR,INFO1,INFO2,INMAX
!         WRITE(6,601) (W(IN),IN=1,INMAX)
!  601    FORMAT(1P5E14.6)

       GX(NR) =GUCLIP(XRHO(NR))
       NGY(NR)=INMAX
       DO NY=1,INMAX
          GY(NR,NY)=GUCLIP(1.D0/(2*PI*1.D6*SQRT(W(NY))))
!               WRITE(6,602) NR,XRHO(NR),GY(NR,NY)
!  602          FORMAT(I3,1P2E13.5)
       ENDDO
    ENDDO

1000 CONTINUE
    CALL GSGRAF(GX,GY,NGY,NRMAXP,nrmax+1, &
                0.0,REAL(RB/RA),REAL(WAEMIN),REAL(WAEMAX))

    RETURN
  END SUBROUTINE wm_taem

!
!**********************************************************************

  SUBROUTINE GSGRAF(GX,GY,NGY,NXMAX,NXM,GXMIN1,GXMAX1,GYMIN1,GYMAX1)

    IMPLICIT NONE
    INTEGER,INTENT(IN):: NXMAX,NXM
    REAL(4),INTENT(IN):: GX(NXM),GY(NXM,*)
    INTEGER,INTENT(IN):: NGY(NXM)
    REAL(4),INTENT(IN):: GXMIN1,GXMAX1,GYMIN1,GYMAX1
    INTERFACE
       FUNCTION GUCLIP(X)
         USE bpsd_kinds
         REAL(rkind):: X
         REAL(sp):: GUCLIP
       END FUNCTION GUCLIP
       FUNCTION NGULEN(X)
         USE bpsd_kinds
         REAL(sp):: X
         INTEGER:: NGULEN
       END FUNCTION NGULEN
    END INTERFACE

    
    REAL(4):: GXMIN,GXMAX,GXSTEP,GYMIN,GYMAX,GYSTEP
    INTEGER:: NX,NY
    
    CALL GQSCAL(GXMIN1,GXMAX1,GXMIN,GXMAX,GXSTEP)
    CALL GQSCAL(GYMIN1,GYMAX1,GYMIN,GYMAX,GYSTEP)
    CALL PAGES
    CALL SETCHS(0.5,0.0)
    CALL GDEFIN(3.0,16.0,3.0,17.0,GXMIN,GXMAX,GYMIN,GYMAX)
    CALL SETLIN(0,0,7)
    CALL GFRAME
    CALL GSCALE(GXMIN,GXSTEP,GYMIN,GYSTEP,0.2,9)
    CALL GVALUE(GXMIN,2*GXSTEP,0.0,0.0,NGULEN(2*GXSTEP))
    CALL GVALUE(0.0,0.0,GYMIN,2*GYSTEP,NGULEN(2*GYSTEP))
    CALL SETMKS(1,0.1)
    DO NX=1,NXMAX
!         WRITE(6,*) NX,NGY(NX)
       DO NY=1,NGY(NX)
!            WRITE(6,*) NY,GX(NX),GY(NX,NY)
          CALL GPLOTP(GX(NX),GY(NX,NY),1,1,1,2,1,0)
       ENDDO
    ENDDO
    CALL PAGEE
    RETURN
  END SUBROUTINE GSGRAF
END MODULE wmtaem
