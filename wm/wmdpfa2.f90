
! ******************************************************
!                       WMDPFA2      local definition 
!      dielectric tensor used by arbital istribution function 
!      without relativistic effects

! ******************************************************

SUBROUTINE WMDPFA2(CW,AM,RHOL,RKPR,VTA,CPM1,CPM2,CQM1,CQM2,CRM1,CRM2)

!  USE libbes,ONLY: bessjn
  USE plcomm
  USE pllocal
  USE bpsd_constants,ONLY : CI,PI
!  USE dpfpin,ONLY: dpfmfl 
!  INCLUDE '../dp/dpcomm.inc'

  IMPLICIT NONE
  COMPLEX(8),INTENT(IN) :: CW
  REAL(8),INTENT(IN) :: RKPR,RHOL,AM,VTA
  COMPLEX(8),INTENT(OUT) :: CPM1,CPM2,CQM1,CQM2,CRM1,CRM2
  REAL(8) :: p,v
  REAL(8) :: p0,v0,PTA,Norml
  COMPLEX(8) :: CX
  REAL(8) :: AP,AQ,AR
!  REAL(8) :: DIFFRHO,DIFFRHOMIN
!  INTEGER :: NRFA2
  REAL(8) :: DELP1,DELP2
  COMPLEX(8) :: funcD
  INTEGER :: NP,NPMAX,NTH,NTHMAX,NR,NRMAX
  REAL(8) :: PMAX,DELP,RMAX,DELR,DELTH
  REAL(8) :: FM,SUMFM,SUMFMN
  REAL(8) :: PM,THM,RM
  REAL(8) :: PG1,PG2
  REAL(8),ALLOCATABLE :: DFR(:,:),DFP(:,:)

 NPMAX=100
 NTHMAX=32
 NRMAX=32
 PMAX=5.D-19  !not normalized
 RMAX=1.D0
 DELP=PMAX/NPMAX
 DELTH=PI/NTHMAX
 DELR=RMAX/NRMAX

 PTA=AM*VTA
 Norml=1.D0/(SQRT(PI)*PI*VTA*VTA*VTA)
!-- NR in dp from rho of wm
!   DIFFRHOMIN=1.2D0
!   RM=0.D0
!    DO NR=1,NRMAX
!      RM=(NR-0.5D0)*DELR
!      DIFFRHO=ABS(RHOL-RM)
!     IF(DIFFRHO.LT.DIFFRHOMIN) THEN
!      DIFFRHOMIN=DIFFRHO
!      NRFA2=NR
!     ENDIF
!    END DO 
!--

   ALLOCATE (DFR(1:NPMAX,1:NTHMAX))
   ALLOCATE (DFP(1:NPMAX-1,1:NTHMAX))

!-- df/dp,df/dtheta with Test function Maxwellian

      DO NTH=1,NTHMAX
      DO NP=1,NPMAX-1
          PM=(NP-0.5D0)*DELP
          THM=(NTH-0.5D0)*DELTH
          DFP(NP,NTH)=Norml*(FM(PM+DELP,THM,PTA)-FM(PM,THM,PTA))/DELP
      END DO
      END DO

      DO NTH=1,NTHMAX
      DO NP=1,NPMAX
          PM=(NP-0.5D0)*DELP
          THM=(NTH-0.5D0)*DELTH
          DFR(NP,NTH)=Norml*FM(PM,THM,PTA)
      END DO
      END DO
       SUMFM=0.D0
      DO NTH=1,NTHMAX
      DO NP=1,NPMAX-1
          PM=(NP-0.5D0)*DELP
          THM=(NTH-0.5D0)*DELTH
      WRITE(6,'(2I5,1P6E12.4)') NP,NTH,Norml*FM(PM,THM,PTA) !DFP(NP,NTH),DFR(NP,NTH)
          SUMFM=SUMFM+PM*PM*SIN(THM)*FM(PM,THM,PTA)
      ENDDO
      ENDDO
          SUMFM=2*PI*SUMFM*DELP*DELTH/(AM*AM*AM)
          SUMFMN=Norml*SUMFM        
!      WRITE(6,'(1P6E12.4)')  SUMFM,Norml,SUMFMN
!--

!*************SUM1********************
      CX=ABS(RKPR)/(CW*AM)
      CPM1=(0.D0,0.D0)
      CPM2=(0.D0,0.D0)
      CQM1=(0.D0,0.D0)
      CQM2=(0.D0,0.D0)
      CRM1=(0.D0,0.D0)
      CRM2=(0.D0,0.D0)

     DO NTH=1,NTHMAX
        THM=(NTH-0.5D0)*DELTH
        IF(ABS(COS(THM)).LT.1.D-21) THEN
         p0=5.D-19
        ELSE
         p0=1.D0/(REAL(CX)*COS(THM))
        ENDIF
         v0=p0/AM
         AP=SIN(THM)*((1.D0+COS(THM)*COS(THM))**2)
         AQ=SIN(THM)*COS(THM)*(1.D0+COS(THM)*COS(THM))
         AR=SIN(THM)*COS(THM)*COS(THM)
     DO NP=1,NPMAX-1
         p=(NP-0.5D0)*DELP
         v=p/AM
         funcD=1.D0-CX*p*COS(THM)
         PG1=NP*DELP
         PG2=(NP+1)*DELP
!!----INCLUDE SINGULAR POINT
      IF((PG1.LT.p0).AND.(PG2.GT.p0)) THEN
        DELP1=p0-PG1
        DELP2=PG2-p0

!-----PRICIPAL VALUE
        CPM1=CPM1 &
        +(DFP(NP,NTH)*(v**5)*AP/(1.D0-CX*COS(THM)*(PG1+0.5D0*DELP1)))*DELP1*DELTH &
        +(DFP(NP,NTH)*(v**5)*AP/(1.D0-CX*COS(THM)*(PG2-0.5D0*DELP2)))*DELP2*DELTH
        CPM2=CPM2 &
        +(DFR(NP,NTH)*(v**6)*AP/(1.D0-CX*COS(THM)*(PG1+0.5D0*DELP1)))*DELP1*DELTH &
        +(DFR(NP,NTH)*(v**6)*AP/(1.D0-CX*COS(THM)*(PG2-0.5D0*DELP2)))*DELP2*DELTH

        CQM1=CQM1 &
        +(DFP(NP,NTH)*(v**4)*AQ/(1.D0-CX*COS(THM)*(PG1+0.5D0*DELP1)))*DELP1*DELTH &
        +(DFP(NP,NTH)*(v**4)*AQ/(1.D0-CX*COS(THM)*(PG2-0.5D0*DELP2)))*DELP2*DELTH
        CQM2=CQM2 &
        +(DFR(NP,NTH)*(v**5)*AQ/(1.D0-CX*COS(THM)*(PG1+0.5D0*DELP1)))*DELP1*DELTH &
        +(DFR(NP,NTH)*(v**5)*AQ/(1.D0-CX*COS(THM)*(PG2-0.5D0*DELP2)))*DELP2*DELTH

        CRM1=CRM1 &
        +(DFP(NP,NTH)*(v**3)*AR/(1.D0-CX*COS(THM)*(PG1+0.5D0*DELP1)))*DELP1*DELTH &
        +(DFP(NP,NTH)*(v**3)*AR/(1.D0-CX*COS(THM)*(PG2-0.5D0*DELP2)))*DELP2*DELTH
        CRM2=CRM2 &
        +(DFR(NP,NTH)*(v**4)*AR/(1.D0-CX*COS(THM)*(PG1+0.5D0*DELP1)))*DELP1*DELTH &
        +(DFR(NP,NTH)*(v**4)*AR/(1.D0-CX*COS(THM)*(PG2-0.5D0*DELP2)))*DELP2*DELTH
!------END PRINCIPAL VALUE 

!-----SINGULAR POINT
        IF(AIMAG(CW).GT.0.D0) THEN
         CPM1=CPM1
         CPM2=CPM2
         CQM1=CQM1
         CQM2=CQM2
         CRM1=CRM1
         CRM2=CRM2
        ELSEIF(AIMAG(CW).EQ.0.D0) THEN  
         CPM1=CPM1-CI*PI*(DFP(NP,NTH)*(v0**5)*AP)*p0 *DELTH
         CPM2=CPM2-CI*PI*(DFR(NP,NTH)*(v0**6)*AP)*p0 *DELTH
         CQM1=CQM1-CI*PI*(DFP(NP,NTH)*(v0**4)*AQ)*p0 *DELTH
         CQM2=CQM2-CI*PI*(DFR(NP,NTH)*(v0**5)*AQ)*p0 *DELTH
         CRM1=CRM1-CI*PI*(DFP(NP,NTH)*(v0**3)*AR)*p0 *DELTH
         CRM2=CRM2-CI*PI*(DFR(NP,NTH)*(v0**4)*AR)*p0 *DELTH
      ELSE  
!(AIMAG(CW).LT.0.D0)
         CPM1=CPM1-2.D0*CI*PI*(DFP(NP,NTH)*(v0**5)*AP)*p0 *DELTH
         CPM2=CPM2-2.D0*CI*PI*(DFR(NP,NTH)*(v0**6)*AP)*p0 *DELTH
         CQM1=CQM1-2.D0*CI*PI*(DFP(NP,NTH)*(v0**4)*AQ)*p0 *DELTH
         CQM2=CQM2-2.D0*CI*PI*(DFR(NP,NTH)*(v0**5)*AQ)*p0 *DELTH
         CRM1=CRM1-2.D0*CI*PI*(DFP(NP,NTH)*(v0**3)*AR)*p0 *DELTH
         CRM2=CRM2-2.D0*CI*PI*(DFR(NP,NTH)*(v0**4)*AR)*p0 *DELTH
        ENDIF
!-----END SINGULAR POINT

!!-----END INCLUDE SINGULAR POINT
      ELSE
!-----WITHOUT SINGULAR POINT
        CPM1=CPM1+(DFP(NP,NTH)*(v**5)*AP/funcD)*DELP*DELTH
        CPM2=CPM2+(DFR(NP,NTH)*(v**6)*AP/funcD)*DELP*DELTH

        CQM1=CQM1+(DFP(NP,NTH)*(v**4)*AQ/funcD)*DELP*DELTH
        CQM2=CQM2+(DFR(NP,NTH)*(v**5)*AQ/funcD)*DELP*DELTH

        CRM1=CRM1+(DFP(NP,NTH)*(v**3)*AR/funcD)*DELP*DELTH
        CRM2=CRM2+(DFR(NP,NTH)*(v**4)*AR/funcD)*DELP*DELTH
      ENDIF
!-----END WITHOUT SINGULAR POINT
     END DO
     END DO

  DEALLOCATE(DFR)
  DEALLOCATE(DFP)
!     CPM1=(0.D0,0.D0)
!     CPM2=(0.D0,0.D0)
!     CQM1=(0.D0,0.D0)
!     CQM2=(0.D0,0.D0)
!     CRM1=(0.D0,0.D0)
!     CRM2=(0.D0,0.D0)  

    RETURN
END SUBROUTINE WMDPFA2                    

FUNCTION FM(PM,THM,PTA)
 IMPLICIT NONE
 REAL(8),INTENT(IN) :: PM,THM,PTA
 REAL(8) :: FM
 REAL(8) :: ET

 ET=PM*PM/(PTA*PTA)
 FM=EXP(-ET)
END FUNCTION FM
