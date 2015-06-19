
! ******************************************************
!                       WMDPFA2              
!      dielectric tensor used by arbital istribution function 
!      without relativistic effects

! ******************************************************

SUBROUTINE WMDPFA2(NS,CW,RHOL,RKPR,CPM1,CPM2,CQM1,CQM2,CRM1,CRM2)

!  USE libbes,ONLY: bessjn
  USE plcomm
  USE pllocal
  USE bpsd_constants,ONLY : CI,PI,AMP
! USE dpfpin,ONLY: dpfmfl 
  INCLUDE '../dp/dpcomm.inc'

!  IMPLICIT NONE
  INTEGER,INTENT(IN) :: NS 
  COMPLEX(8),INTENT(IN) :: CW
  REAL(8),INTENT(IN) :: RKPR,RHOL
  COMPLEX(8),INTENT(OUT) :: CPM1,CPM2,CQM1,CQM2,CRM1,CRM2
  REAL(8) :: p,v
  REAL(8) :: p0,v0
  REAL(8) :: AM
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

   ALLOCATE (DFR(1:NPMAX,1:NTHMAX))
   DELTH=PI/NTHMAX

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
     CALL DPFMFL(NS,ID)
      DO NTH=1,NTHMAX
      DO NP=1,NPMAX-1
          DFP(NP,NTH)=(FM(NP+1,NTH) - FM(NP,NTH))/DELP(NS)
      END DO
      END DO

      DO NTH=1,NTHMAX
      DO NP=1,NPMAX
          DFR(NP,NTH)=FM(NP,NTH)
      END DO
      END DO
!--

!***********DGP1,DGP2,DGT1,DGT2************
!
!      DO NP=1,NPMAX
!      DO NTH=1,NTHMAX
!         DGP1(NP,NTH)=-DKPRW*TCSM(NTH)
!         DGP2(NP,NTH)=-DKPRW*TCSG(NTH)
!         DGT1(NP,NTH)= DKPRW*PG(NP,NS)*TSNM(NTH)
!         DGT2(NP,NTH)= DKPRW*PM(NP,NS)*TSNG(NTH)
!      END DO
!      END DO
!--
!*************SUM1********************
      AM=AMP*PA(NS)
      CX=ABS(RKPR)/(CW*AM)
      CPM1=(0.D0,0.D0)
      CPM2=(0.D0,0.D0)
      CQM1=(0.D0,0.D0)
      CQM2=(0.D0,0.D0)
      CRM1=(0.D0,0.D0)
      CRM2=(0.D0,0.D0)

     DO NTH=1,NTHMAX
        IF(ABS(TCSM(NTH)).LT.1.D-21) THEN
         p0=3.D1
        ELSE
         p0=REAL(1.D0/(CX*TCSM(NTH)))
        ENDIF
         v0=p0/AM
         AP=TSNM(NTH)*((1.D0+TCSM(NTH)*TCSM(NTH))**2)
         AQ=TSNM(NTH)*TCSM(NTH)*(1.D0+TCSM(NTH)*TCSM(NTH))
         AR=TSNM(NTH)*TCSM(NTH)*TCSM(NTH)
     DO NP=1,NPMAX-1
         p=PM(NP,NS)
         v=p/AM
         funcD=1.D0-CX*p*TCSM(NTH)
!----INCLUDE SINGULAR POINT
      IF((PG(NP,NS).LT.p0).AND.(PG(NP+1,NS).GT.p0)) THEN
        DELP1=p0-PG(NP,NS)
        DELP2=PG(NP+1,NS)-p0  
!!--detailed integral
!        DELDELP1=(p0-PG(NP,NS))/5.D0
!        DELDELP2=(PG(NP+1,NS)-p0)/5.D0
!  DO i=1,5
!      CPM1=CPM1   & 
!     +(DFP(NP,NTH)*(v**5)*AP/(1.D0-CX*TCSM(NTH)*(PG(NP,NS)+(i-1)*DELDELP1)))*DELDELP1*DELTH &
!     +(DFP(NP,NTH)*(v**5)*AP/(1.D0-CX*TCSM(NTH)*(PG(NP+1,NS)-(i-1)*DELDELP2)))*DELDELP2*DELTH
!          :
!          :
!  ENDDO
!!---detailed integral END

!-----PRICIPAL VALUE
        CPM1=CPM1   &
        +(DFP(NP,NTH)*(v**5)*AP/(1.D0-CX*TCSM(NTH)*(PG(NP,NS)+0.5D0*DELP1)))*DELP1*DELTH &
        +(DFP(NP,NTH)*(v**5)*AP/(1.D0-CX*TCSM(NTH)*(PG(NP+1,NS)-0.5D0*DELP2)))*DELP2*DELTH
        CPM2=CPM2   &
        +(DFR(NP,NTH)*(v**6)*AP/(1.D0-CX*TCSM(NTH)*(PG(NP,NS)+0.5D0*DELP1)))*DELP1*DELTH &
        +(DFR(NP,NTH)*(v**6)*AP/(1.D0-CX*TCSM(NTH)*(PG(NP+1,NS)-0.5D0*DELP2)))*DELP2*DELTH

        CQM1=CQM1   &
        +(DFP(NP,NTH)*(v**4)*AQ/(1.D0-CX*TCSM(NTH)*(PG(NP,NS)+0.5D0*DELP1)))*DELP1*DELTH &
        +(DFP(NP,NTH)*(v**4)*AQ/(1.D0-CX*TCSM(NTH)*(PG(NP+1,NS)-0.5D0*DELP2)))*DELP2*DELTH
        CQM2=CQM2   &
        +(DFR(NP,NTH)*(v**5)*AQ/(1.D0-CX*TCSM(NTH)*(PG(NP,NS)+0.5D0*DELP1)))*DELP1*DELTH &
        +(DFR(NP,NTH)*(v**5)*AQ/(1.D0-CX*TCSM(NTH)*(PG(NP+1,NS)-0.5D0*DELP2)))*DELP2*DELTH

        CRM1=CRM1   &
        +(DFP(NP,NTH)*(v**3)*AR/(1.D0-CX*TCSM(NTH)*(PG(NP,NS)+0.5D0*DELP1)))*DELP1*DELTH &
        +(DFP(NP,NTH)*(v**3)*AR/(1.D0-CX*TCSM(NTH)*(PG(NP+1,NS)-0.5D0*DELP2)))*DELP2*DELTH
        CRM2=CRM2   &
        +(DFR(NP,NTH)*(v**4)*AR/(1.D0-CX*TCSM(NTH)*(PG(NP,NS)+0.5D0*DELP1)))*DELP1*DELTH &
        +(DFR(NP,NTH)*(v**4)*AR/(1.D0-CX*TCSM(NTH)*(PG(NP+1,NS)-0.5D0*DELP2)))*DELP2*DELTH 
!------PRINCIPAL VALUE END

!-----SINGULAR POINT
        IF(AIMAG(CW).GT.0.D0) THEN
         CPM1=CPM1
         CPM2=CPM2
         CQM1=CQM1
         CQM2=CQM2
         CRM1=CRM1
         CRM2=CRM2
        ELSEIF(AIMAG(CW).EQ.0.D0) THEN
         CPM1=CPM1+CI*PI*(DFP(NP,NTH)*(v0**5)*AP)*p0
         CPM2=CPM2+CI*PI*(DFR(NP,NTH)*(v0**6)*AP)*p0
         CQM1=CQM1+CI*PI*(DFP(NP,NTH)*(v0**4)*AQ)*p0
         CQM2=CQM2+CI*PI*(DFR(NP,NTH)*(v0**5)*AQ)*p0
         CRM1=CRM1+CI*PI*(DFP(NP,NTH)*(v0**3)*AR)*p0
         CRM2=CRM2+CI*PI*(DFR(NP,NTH)*(v0**4)*AR)*p0
        ELSE  !(AIMAG(CW).LT.0.D0)
         CPM1=CPM1+2.D0*CI*PI*(DFP(NP,NTH)*(v0**5)*AP)*p0
         CPM2=CPM2+2.D0*CI*PI*(DFR(NP,NTH)*(v0**6)*AP)*p0
         CQM1=CQM1+2.D0*CI*PI*(DFP(NP,NTH)*(v0**4)*AQ)*p0
         CQM2=CQM2+2.D0*CI*PI*(DFR(NP,NTH)*(v0**5)*AQ)*p0
         CRM1=CRM1+2.D0*CI*PI*(DFP(NP,NTH)*(v0**3)*AR)*p0
         CRM2=CRM2+2.D0*CI*PI*(DFR(NP,NTH)*(v0**4)*AR)*p0
        ENDIF
!-----SINGULAR POINT END

!-----INCLUDE SINGULAR POINT END   
      ELSE 
!-----WITHOUT SINGULAR POINT
    
        CPM1=CPM1+(DFP(NP,NTH)*(v**5)*AP/funcD)*DELP(NS)*DELTH
        CPM2=CPM2+(DFR(NP,NTH)*(v**6)*AP/funcD)*DELP(NS)*DELTH

        CQM1=CQM1+(DFP(NP,NTH)*(v**4)*AQ/funcD)*DELP(NS)*DELTH   
        CQM2=CQM2+(DFR(NP,NTH)*(v**5)*AQ/funcD)*DELP(NS)*DELTH

        CRM1=CRM1+(DFP(NP,NTH)*(v**3)*AR/funcD)*DELP(NS)*DELTH
        CRM2=CRM2+(DFR(NP,NTH)*(v**4)*AR/funcD)*DELP(NS)*DELTH
      ENDIF
!-----WITHOUT SINGULAR POINT
     END DO
     END DO

  DEALLOCATE(DFR)
     
!     CPM1=(0.D0,0.D0)
!     CPM2=(0.D0,0.D0)
!     CQM1=(0.D0,0.D0)
!     CQM2=(0.D0,0.D0)
!     CRM1=(0.D0,0.D0)
!     CRM2=(0.D0,0.D0)  

    RETURN
END SUBROUTINE WMDPFA2                    
