C     $Id$
C   ***********************************************
C   **     Analysis of Alfven Eigenmode          **
C   ***********************************************
C
      SUBROUTINE WMTAE
C
      USE plprof,ONLY:pl_prof_old
      USE pllocal
      INCLUDE 'wmcomm.inc'
C
      PARAMETER(MMM=MDM*NDM)
      DIMENSION CWALF(MDM,NDM),CWALFK(MDM,NDM)
      DIMENSION GX(NRM),GY(NRM,MMM)
      DIMENSION NGY(NRM)
C
      DIMENSION DM(MMM),EM(MMM),FM(MMM,MMM),QM(1,MMM),W(MMM)
      DIMENSION IBLOCK(MMM),ISPLIT(MMM)
      DIMENSION WORK(5*MMM),IWORK(3*MMM)
C
      INMAX=NRM 
      NRMAXP=NRMAX
C
      MMMAX=MAX(1,NDSIZ-1)*MAX(1,MDSIZ-1)
      MHMAX=MAX(1,NDSIZ/2)*MAX(1,MDSIZ/2)-1
C
      DO NR=1,NRMAX
         RHON=XRHO(NR)
         CALL PL_PROF_OLD(RHON)
         RHOM=0.D0
         DO NS=2,NSMAX
            RHOM=RHOM+RN(NS)*1.D20*PA(NS)*AMP
         ENDDO
         IF(RHOM.LE.0.D0) THEN
            NRMAXP=NR-1
            GOTO 1000
         ENDIF
C
         RIOTAL=QPS(NR)
C
         rho=xrho(NR)
         DO NHH=1,NHHMAX
            ph=(NHH-1)*2.D0*PI/NHHMAX
         DO NTH=1,NTHMAX
            th=(NTH-1)*2.D0*PI/NTHMAX
            CALL wmfem_magnetic(rho,th,ph,babs,bsupr,bsupth,bsupph)
            VALF=BABS/SQRT(RMU0*RHOM)
C            CWALF(NTH,NHH)=BABS**2/(VALF**2*BSUPTH**2)
            CWALF(NTH,NHH)=RPS(NTH,NR)**2/VALF**2
         ENDDO
         ENDDO
C
         CALL WMSUBFX(CWALF,CWALFK,NTHMAX,NHHMAX)
C
         DO ND=NDMIN,MAX(NDMIN,NDMAX-1)
            NDX=ND-NDMIN+1
            NN=NPH0+NHC*ND
         DO MD=MDMIN,MAX(MDMIN,MDMAX-1)
            MDX=MD-MDMIN+1
            MM=NTH0+MD
               
            RKPR1=MM+NN*RIOTAL
C
            NW=MAX(1,MDSIZ-1)*(NDX-1)+(MDX-1)+1
            NW0=MAX(1,MDSIZ-1)*(-NDMIN)+(-MDMIN)+1
            DO ND1=NDMIN,MAX(NDMIN,NDMAX-1)
               NDX1=ND1-NDMIN+1
               NN1=NN+NHC*ND1
            DO MD1=MDMIN,MAX(MDMIN,MDMAX-1)
               MDX1=MD1-MDMIN+1
               MM1=MM+MD1
C
               RKPR2=MM1+NN1*RIOTAL
C
               NW1=MAX(1,MDSIZ-1)*(NDX1-1)+(MDX1-1)+1
               NW1=NW1-NW0+1
               IF(NW1.GE.1.AND.NW1.LE.MHMAX+1) THEN
                  FM(NW1,NW)=DBLE(CWALFK(MDX1,NDX1))/(RKPR1*RKPR2)
     &                      *RIOTAL**2
               ENDIF
            ENDDO
            ENDDO
C
C            DM(NWW)=DBLE(CWALFK(0-MDMIN+1,0-NDMIN+1))
C            EM(NWW)=DBLE(CWALFK(1-MDMIN+1,0-NDMIN+1))
         ENDDO
         ENDDO
C
         CALL LAPACK_DSBTRD('N','L',MMMAX,MHMAX,
     &                      FM,MMM,DM,EM,QM,1,WORK,INFO1)
C
         EPSEG=0.D0
         VU=1.D0/(2.D0*PI*1.D6*WAEMIN)**2
         VL=1.D0/(2.D0*PI*1.D6*WAEMAX)**2
C         VU=(2.D0*PI*1.D6*OME)**2
C         VL=(2.D0*PI*1.D6*OMS)**2
C
         CALL LAPACK_DSTEBZ('V','B',MMMAX,VL,VU,IMIN,IMAX,EPSEG,
     &               DM,EM,INMAX,NSPLIT,W,
     &               IBLOCK,ISPLIT,WORK,IWORK,INFO2)
C         WRITE(6,*) 'NR,INFO1,INFO2,INMAX=',NR,INFO1,INFO2,INMAX
C         WRITE(6,601) (W(IN),IN=1,INMAX)
C  601    FORMAT(1P5E14.6)
C
         GX(NR) =GUCLIP(XRHO(NR))
         NGY(NR)=INMAX
         DO NY=1,INMAX
            GY(NR,NY)=GUCLIP(1.D0/(2*PI*1.D6*SQRT(W(NY))))
C               WRITE(6,602) NR,XRHO(NR),GY(NR,NY)
C  602          FORMAT(I3,1P2E13.5)
         ENDDO
      ENDDO
C
 1000 CALL GSGRAF(GX,GY,NGY,NRMAXP,NRM,
     &            0.0,REAL(RB/RA),REAL(WAEMIN),REAL(WAEMAX))
C
      RETURN
      END
C
C
C**********************************************************************
C
      SUBROUTINE GSGRAF(GX,GY,NGY,NXMAX,NXM,
     &                  GXMIN1,GXMAX1,GYMIN1,GYMAX1)
C
      DIMENSION GX(NXM),GY(NXM,*),NGY(NXM)
C
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
C         WRITE(6,*) NX,NGY(NX)
         DO NY=1,NGY(NX)
C            WRITE(6,*) NY,GX(NX),GY(NX,NY)
            CALL GPLOTP(GX(NX),GY(NX,NY),1,1,1,2,1,0)
         ENDDO
      ENDDO
      CALL PAGEE
      RETURN
      END
