
! ******************************************************
!                       WMDPFA2              
!      dielectric tensor used by arbital istribution function 
!      without relativistic effects

! ******************************************************

SUBROUTINE WMDPFA2(NS,CW,RHOL,RKPR,VTA,CPM1,CPM2,CQM1,CQM2,CRM1,CRM2)

  USE plcomm
  USE pllocal
  USE bpsd_constants,ONLY : CI,PI,AMP
! USE dpfpin,ONLY: dpfmfl 
  INCLUDE '../dp/dpcomm.inc'

!  IMPLICIT NONE
  INTEGER,INTENT(IN) :: NS 
  COMPLEX(8),INTENT(IN) :: CW
  REAL(8),INTENT(IN) :: RKPR,RHOL,VTA
  COMPLEX(8),INTENT(OUT) :: CPM1,CPM2,CQM1,CQM2,CRM1,CRM2
  REAL(8) :: p,v
  REAL(8) :: p0,v0,PTA
  REAL(8) :: AM,Norml
  COMPLEX(8) :: CX
  REAL(8) :: AP,AQ,AR
  REAL(8) :: DIFFRHO,DIFFRHOMIN,DELR
  INTEGER :: NRFA2
  REAL(8) :: DELP1,DELP2
  COMPLEX(8) :: funcD

  REAL(8),ALLOCATABLE :: DFR(:,:)

!-- NR in dp from rho of wm
   DIFFRHOMIN=1.2D0
    DO NR=1,NRMAX
      DIFFRHO=ABS(RHOL-RM(NR))
     IF(DIFFRHO.LT.DIFFRHOMIN) THEN
      DIFFRHOMIN=DIFFRHO
      NRFA2=NR
!     ELSE
!      RETURN 
     ENDIF
    END DO
!--
     NPMAX=128
   ALLOCATE (DFR(1:NPMAX,1:NTHMAX))
     
     DELTH=PI/NTHMAX
     AM=AMP*PA(NS)
     PTA=VTA*AM
     CX=ABS(RKPR)/(CW*AM)
     CPM1=(0.D0,0.D0)
     CPM2=(0.D0,0.D0)
     CQM1=(0.D0,0.D0)
     CQM2=(0.D0,0.D0)
     CRM1=(0.D0,0.D0)
     CRM2=(0.D0,0.D0)
!   WRITE(6,'(1P6E12.4)') PM(NPMAX,NS)
!!-- df/dp,df/dtheta with FP
!      DO NTH=1,NTHMAX
!      DO NP=1,NPMAX-1
!         DFP(NP,NTH) = (FP(NRFA2,NP+1,NTH) - FP(NRFA2,NP,NTH))/DELP(NS)
!      END DO
!      END DO
!        DELR=RM(NRFA2+1)-RM(NRFA2)
!      DO NP=1,NPMAX
!      DO NTH=1,NTHMAX
!         DFR(NP,NTH) = (FP(NRFA2+1,NP,NTH) - FP(NRFA2,NP,NTH))/DELR
!      END DO
!      END DO
!!--

!-- df/dp,df/dtheta with Test function Maxwellian
      ID=0
      Norml=1.D0/(VTA*VTA*VTA)
     CALL DPFMFL(NS,ID)
      DO NTH=1,NTHMAX
      DO NP=1,NPMAX-1
          DFP(NP,NTH)=Norml*(FM(NP+1,NTH)-FM(NP,NTH))/(DELP(NS)*PTA)
!         WRITE(6,'(1P6E12.4)') DFP(NP,NTH),FM(NP,NTH)
      END DO
      END DO

      DO NTH=1,NTHMAX
      DO NP=1,NPMAX
          DFR(NP,NTH)=Norml*FM(NP,NTH)
!        WRITE(6,'(1P6E12.4)') FM(NP,NTH),PM(NP,NS),REAL(NP)
      END DO
      END DO
!         WRITE(6,'(1P6E12.4)') TPP,FACTOR

!--
     DO NTH=1,NTHMAX
        IF(ABS(TCSM(NTH)).LT.1.D-21) THEN
         p0=7.D0*PTA
        ELSE
         p0=1.D0/(REAL(CX)*TCSM(NTH))
        ENDIF
         v0=p0/AM
         AP=TSNM(NTH)*((1.D0+TCSM(NTH)*TCSM(NTH))**2)
         AQ=TSNM(NTH)*TCSM(NTH)*(1.D0+TCSM(NTH)*TCSM(NTH))
         AR=TSNM(NTH)*TCSM(NTH)*TCSM(NTH)
     DO NP=1,NPMAX-1
         p=PM(NP,NS)*PTA
         v=p/AM
!    WRITE(6,'(1P6E12.4)') p0,v0,p,v
         funcD=1.D0-CX*p*TCSM(NTH)
!----INCLUDE SINGULAR POINT
      IF(((PG(NP,NS)*PTA).LT.p0).AND.((PG(NP+1,NS)*PTA).GT.p0)) THEN
        DELP1=p0-PG(NP,NS)*PTA
        DELP2=PG(NP+1,NS)*PTA-p0
!    WRITE(6,'(1P6E12.4)') p0,v0,REAL(NP),REAL(NTH)
!-----PRICIPAL VALUE
        CPM1=CPM1   &
        +(DFP(NP,NTH)*(v**5)*AP/(1.D0-CX*TCSM(NTH)*(PG(NP,NS)*PTA+0.5D0*DELP1)))*DELP1*DELTH &
        +(DFP(NP,NTH)*(v**5)*AP/(1.D0-CX*TCSM(NTH)*(PG(NP+1,NS)*PTA-0.5D0*DELP2)))*DELP2*DELTH
        CPM2=CPM2   &
        +(DFR(NP,NTH)*(v**6)*AP/(1.D0-CX*TCSM(NTH)*(PG(NP,NS)*PTA+0.5D0*DELP1)))*DELP1*DELTH &
        +(DFR(NP,NTH)*(v**6)*AP/(1.D0-CX*TCSM(NTH)*(PG(NP+1,NS)*PTA-0.5D0*DELP2)))*DELP2*DELTH

        CQM1=CQM1   &
        +(DFP(NP,NTH)*(v**4)*AQ/(1.D0-CX*TCSM(NTH)*(PG(NP,NS)*PTA+0.5D0*DELP1)))*DELP1*DELTH &
        +(DFP(NP,NTH)*(v**4)*AQ/(1.D0-CX*TCSM(NTH)*(PG(NP+1,NS)*PTA-0.5D0*DELP2)))*DELP2*DELTH
        CQM2=CQM2   &
        +(DFR(NP,NTH)*(v**5)*AQ/(1.D0-CX*TCSM(NTH)*(PG(NP,NS)*PTA+0.5D0*DELP1)))*DELP1*DELTH &
        +(DFR(NP,NTH)*(v**5)*AQ/(1.D0-CX*TCSM(NTH)*(PG(NP+1,NS)*PTA-0.5D0*DELP2)))*DELP2*DELTH

        CRM1=CRM1   &
        +(DFP(NP,NTH)*(v**3)*AR/(1.D0-CX*TCSM(NTH)*(PG(NP,NS)*PTA+0.5D0*DELP1)))*DELP1*DELTH &
        +(DFP(NP,NTH)*(v**3)*AR/(1.D0-CX*TCSM(NTH)*(PG(NP+1,NS)*PTA-0.5D0*DELP2)))*DELP2*DELTH
        CRM2=CRM2   &
        +(DFR(NP,NTH)*(v**4)*AR/(1.D0-CX*TCSM(NTH)*(PG(NP,NS)*PTA+0.5D0*DELP1)))*DELP1*DELTH &
        +(DFR(NP,NTH)*(v**4)*AR/(1.D0-CX*TCSM(NTH)*(PG(NP+1,NS)*PTA-0.5D0*DELP2)))*DELP2*DELTH
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
         CPM1=CPM1-CI*PI*(DFP(NP,NTH)*(v0**5)*AP)*p0*DELTH
         CPM2=CPM2-CI*PI*(DFR(NP,NTH)*(v0**6)*AP)*p0*DELTH
         CQM1=CQM1-CI*PI*(DFP(NP,NTH)*(v0**4)*AQ)*p0*DELTH
         CQM2=CQM2-CI*PI*(DFR(NP,NTH)*(v0**5)*AQ)*p0*DELTH
         CRM1=CRM1-CI*PI*(DFP(NP,NTH)*(v0**3)*AR)*p0*DELTH
         CRM2=CRM2-CI*PI*(DFR(NP,NTH)*(v0**4)*AR)*p0*DELTH
        ELSE  !(AIMAG(CW).LT.0.D0)
         CPM1=CPM1-2.D0*CI*PI*(DFP(NP,NTH)*(v0**5)*AP)*p0*DELTH
         CPM2=CPM2-2.D0*CI*PI*(DFR(NP,NTH)*(v0**6)*AP)*p0*DELTH
         CQM1=CQM1-2.D0*CI*PI*(DFP(NP,NTH)*(v0**4)*AQ)*p0*DELTH
         CQM2=CQM2-2.D0*CI*PI*(DFR(NP,NTH)*(v0**5)*AQ)*p0*DELTH
         CRM1=CRM1-2.D0*CI*PI*(DFP(NP,NTH)*(v0**3)*AR)*p0*DELTH
         CRM2=CRM2-2.D0*CI*PI*(DFR(NP,NTH)*(v0**4)*AR)*p0*DELTH
        ENDIF
!-----END SINGULAR POINT

!-----END INCLUDE SINGULAR POINT
      ELSE
!-----WITHOUT SINGULAR POINT

        CPM1=CPM1+(DFP(NP,NTH)*(v**5)*AP/funcD)*DELP(NS)*PTA*DELTH
        CPM2=CPM2+(DFR(NP,NTH)*(v**6)*AP/funcD)*DELP(NS)*PTA*DELTH

        CQM1=CQM1+(DFP(NP,NTH)*(v**4)*AQ/funcD)*DELP(NS)*PTA*DELTH
        CQM2=CQM2+(DFR(NP,NTH)*(v**5)*AQ/funcD)*DELP(NS)*PTA*DELTH

        CRM1=CRM1+(DFP(NP,NTH)*(v**3)*AR/funcD)*DELP(NS)*PTA*DELTH
        CRM2=CRM2+(DFR(NP,NTH)*(v**4)*AR/funcD)*DELP(NS)*PTA*DELTH
!         WRITE(6,'(1P6E12.4)') DFP(NP,NTH),DFR(NP,NTH)
      ENDIF
!-----END WITHOUT SINGULAR POINT
     END DO
     END DO

  DEALLOCATE(DFR)
     
!     CPM1=(0.D0,0.D0)
!     CPM2=(0.D0,0.D0)
!     CQM1=(0.D0,0.D0)
!     CQM2=(0.D0,0.D0)
!     CRM1=(0.D0,0.D0)
!     CRM2=(0.D0,0.D0)  

!    RETURN
END SUBROUTINE WMDPFA2 
