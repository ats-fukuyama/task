MODULE wmfa
 IMPLICIT NONE
 REAL(8) ::paN
 REAL(8) :: DFPax,DFRax,fax
 REAL(8) ::RHO0,p0N,v0N ! **N -> normalized by Thermal (pN=p/PTA(=Pthermal)) 
 REAL(8) :: A0
 INTEGER :: la,n,s,NTHv
 COMPLEX(8) :: Cp0N,Cv0N,I2Cp0N
 REAL(8),PARAMETER :: PmaxN=7.D0

END MODULE wmfa
 
! ******************************************************
!                       WMDPFAA              
!      dielectric tensor used by arbitral distribution function 
!      without relativistic effects
! ******************************************************

SUBROUTINE WMDPFAA(CW,RHOWM,RKPR,AE2N0,CPM1,CPM2,CQM1,CQM2,CRM1,CRM2)
  USE plcomm
  USE pllocal
  USE bpsd_constants,ONLY : CI,PI,AMP,AEE
  USE wmfa
  INCLUDE '../dp/dpcomm.inc'

!  IMPLICIT NONE
!  INTEGER,INTENT(IN) :: NRWM
  COMPLEX(8),INTENT(IN) :: CW
  REAL(8),INTENT(IN) :: RHOWM,RKPR
  REAL(8),INTENT(OUT) ::AE2N0
  COMPLEX(8),INTENT(OUT) :: CPM1,CPM2,CQM1,CQM2,CRM1,CRM2
  REAL(8) :: p,v
  REAL(8) :: p0,v0,PTA
  REAL(8) :: RTA,VTA
  REAL(8) :: AM,Norml
  REAL(8) :: RHOM_MIN,RHOM_MAX
  COMPLEX(8) :: CX,Cp0
  REAL(8) :: AP,AQ,AR
  REAL(8) :: fa1,DFPa1,DFRa1
  INTEGER :: NS,ID
  COMPLEX(8) :: CIPL1,CIPR1,CIPL2,CIPR2,CIP1,CIP2
  COMPLEX(8) :: CIQL1,CIQR1,CIQL2,CIQR2,CIQ1,CIQ2
  COMPLEX(8) :: CIRL1,CIRR1,CIRL2,CIRR2,CIR1,CIR2

     NS=3
     RHO0=RHOWM
     RHOM_MIN=RHOa0(1)
     RHOM_MAX=RHOa0(NRMAXFP+2)
     AM=AMFP(NS)
     AE=AEFP(NS)
     AE2N0=AE*AE*RNFP0(NS)
     RTA=RTFP0(NS)*AEE*1.D3
     VTA=SQRT(2.D0*RTA/AM)
     PTA=VTA*AM
     CX=ABS(RKPR)/(CW*AM)
     CPM1=(0.D0,0.D0)
     CPM2=(0.D0,0.D0)
     CQM1=(0.D0,0.D0)
     CQM2=(0.D0,0.D0)
     CRM1=(0.D0,0.D0)
     CRM2=(0.D0,0.D0)

! IF((RHO0.GE.RHOM_MIN).AND.(RHO0.LE.RHOM_MAX)) THEN

!     DO NR=1,NRMAX+2   NRMAX => WM's NRMAX -> US(x,x,x,NR>20,x) don't exist  
!         WRITE(6,'(1P4E12.4)') US(1,1,NPMAX,NR,1)
!     ENDDO

!  DO NTH=1,NTHMAX
!  DO NP=1,NPMAX+2
!      NR=2
!      WRITE(6,'(2I5,1P4E12.4)') NP,NTH,fa0(NP,NR,NTH)
!  ENDDO
!  ENDDO
!RETURN 

!--
     DO NTH=1,NTHMAX
        p0N=1.D0/(PTA*REAL(CX)*TCSM(NTH))
        Cp0N=1.D0/(PTA*CX*TCSM(NTH))

        v0N=p0N*PTA/(AM*VTA)
        Cv0N=Cp0N*PTA/(AM*VTA)
        NTHv=NTH
!        WRITE(6,'(1I5,1P4E12.4)') NTHv,Cp0N,p0N
!        WRITE(6,'(A,I5,1P5E12.4)') 'NTH,Cp0N=',NTH,Cp0N,CX,TCSM(NTH)

        AP=TSNM(NTH)*((1.D0+TCSM(NTH)*TCSM(NTH))**2)
        AQ=TSNM(NTH)*TCSM(NTH)*(1.D0+TCSM(NTH)*TCSM(NTH))
        AR=TSNM(NTH)*TCSM(NTH)*TCSM(NTH)
!-----PRICIPAL VALUE
!       IF((p0N.GE.0.D0).AND.(p0N.LT.1.D-3)) THEN
!         CALL PVINT(1,2,CIPR1)
!          CPM1=CPM1 + AP*CIPR1
!         CALL PVINT(2,2,CIPR2)
!          CPM2=CPM2 + AP*CIPR2

!         CALL PVINT(3,2,CIQR1)
!          CQM1=CQM1 + AQ*CIQR1
!         CALL PVINT(4,2,CIQR2)
!          CQM2=CQM2 + AQ*CIQR2

!         CALL PVINT(5,2,CIRR1)
!          CRM1=CRM1 + AR*CIRR1
!         CALL PVINT(6,2,CIRR2)
!          CRM2=CRM2 + AR*CIRR2
!       ElSEIF((p0N.LT.0.D0).OR.(p0N.GE.PmaxN)) THEN

        IF((p0N.LT.0.D0).OR.(p0N.GE.PmaxN)) THEN
         CALL PVINT(1,3,CIP1)
          CPM1=CPM1 + AP*CIP1
         CALL PVINT(2,3,CIP2)
          CPM2=CPM2 + AP*CIP2

         CALL PVINT(3,3,CIQ1)
          CQM1=CQM1 + AQ*CIQ1
         CALL PVINT(4,3,CIQ2)
          CQM2=CQM2 + AQ*CIQ2

         CALL PVINT(5,3,CIR1)
          CRM1=CRM1 + AR*CIR1
         CALL PVINT(6,3,CIR2)
          CRM2=CRM2 + AR*CIR2
!         WRITE(6,'(1P8E12.4)') CIP1,CIP2
       ELSE
         CALL PVINT(1,1,CIPL1)
         CALL PVINT(1,2,CIPR1)
          CPM1=CPM1 + AP*(CIPL1+CIPR1)
         CALL PVINT(2,1,CIPL2)
         CALL PVINT(2,2,CIPR2)
          CPM2=CPM2 + AP*(CIPL2+CIPR2)

         CALL PVINT(3,1,CIQL1)
         CALL PVINT(3,2,CIQR1)
          CQM1=CQM1 + AQ*(CIQL1+CIQR1)
         CALL PVINT(4,1,CIQL2)
         CALL PVINT(4,2,CIQR2)
          CQM2=CQM2 + AQ*(CIQL2+CIQR2)

         CALL PVINT(5,1,CIRL1)
         CALL PVINT(5,2,CIRR1)
          CRM1=CRM1 + AR*(CIRL1+CIRR1)
         CALL PVINT(6,1,CIRL2)
         CALL PVINT(6,2,CIRR2)
          CRM2=CRM2 + AR*(CIRL2+CIRR2)
!         WRITE(6,'(1P8E12.4)') CIRL1,CIRR1,CIRL2,CIRR2
       ENDIF
!------END PRINCIPAL VALUE

!-----SINGULAR POINT
        IF((p0N.GT.0.D0).AND.(p0N.LT.PmaxN)) THEN
         CALL SPL2DD(p0N,RHOWM,fa1,DFPa1,DFRa1,PMa0,RHOa0,&
                    US(1:4,1:4,1:NPM,1:NRM,NTHv),NPM,NPMAX+2,NRMAXFP+2,IERR)
!          WRITE(6,'(3I5)') 2,NRMAX,NRMAXFP
!        WRITE(6,'(1I5,1P6E12.4)') NTHv,p0N,RHO0,fa1,DFPa1,DFRa1
         IF(AIMAG(CW).GT.0.D0) THEN
          CPM1=CPM1
          CPM2=CPM2
          CQM1=CQM1
          CQM2=CQM2
          CRM1=CRM1
          CRM2=CRM2
         ELSEIF(AIMAG(CW).EQ.0.D0) THEN
          CPM1=CPM1-CI*PI*(DFPa1*((v0N)**5)*AP)*p0N !!????
          CPM2=CPM2-CI*PI*(DFRa1*((v0N)**6)*AP)*p0N
          CQM1=CQM1-CI*PI*(DFPa1*((v0N)**4)*AQ)*p0N
          CQM2=CQM2-CI*PI*(DFRa1*((v0N)**5)*AQ)*p0N
          CRM1=CRM1-CI*PI*(DFPa1*((v0N)**3)*AR)*p0N
          CRM2=CRM2-CI*PI*(DFRa1*((v0N)**4)*AR)*p0N
         ELSE  !(AIMAG(CW).LT.0.D0)
          CPM1=CPM1-2.D0*CI*PI*(DFPa1*((Cv0N)**5)*AP)*Cp0N
          CPM2=CPM2-2.D0*CI*PI*(DFRa1*((Cv0N)**6)*AP)*Cp0N
          CQM1=CQM1-2.D0*CI*PI*(DFPa1*((Cv0N)**4)*AQ)*Cp0N
          CQM2=CQM2-2.D0*CI*PI*(DFRa1*((Cv0N)**5)*AQ)*Cp0N
          CRM1=CRM1-2.D0*CI*PI*(DFPa1*((Cv0N)**3)*AR)*Cp0N
          CRM2=CRM2-2.D0*CI*PI*(DFRa1*((Cv0N)**4)*AR)*Cp0N
         ENDIF
!         WRITE(6,'(1P8E12.4)') v0N,p0N 
        ENDIF
!-----END SINGULAR POINT
     END DO

     CPM1=CPM1*DELTH*(VTA**5)   ! dp=p_0*dp_n -> df/dp=(df/dp_n)*dp_n/dp=dfp_n/p_0
     CPM2=CPM2*DELTH*(VTA**6)*PTA
     CQM1=CQM1*DELTH*(VTA**4)
     CQM2=CQM2*DELTH*(VTA**5)*PTA
     CRM1=CRM1*DELTH*(VTA**3)
     CRM2=CRM2*DELTH*(VTA**4)*PTA
! ELSE
!  RETURN
! ENDIF 


!  WRITE(6,'(1P6E12.4)') CPM1,CPM2
!     CPM1=(0.D0,0.D0)
!     CPM2=(0.D0,0.D0)
!     CQM1=(0.D0,0.D0)
!     CQM2=(0.D0,0.D0)
!     CRM1=(0.D0,0.D0)
!     CRM2=(0.D0,0.D0)  
!    RETURN
END SUBROUTINE WMDPFAA


SUBROUTINE PVINT(j,l,CINT)
 USE plcomm
 USE pllocal
 USE bpsd_constants,ONLY : CI,PI
 USE wmfa
 USE libde,ONLY : DEFTC
 INCLUDE '../dp/dpcomm.inc'
! IMPLICIT NONE
 INTEGER,INTENT(IN) :: j  ! j-> (df/dp or df/dr) and power 
 INTEGER,INTENT(IN) :: l  ! l=1-> boundary (0,p0N), l=2-> (p0N,pmaxN), l=3-> (0,pmaxN)
 COMPLEX(8),INTENT(OUT) :: CINT ! Integral
! COMPLEX(8), :: CFUNC1,CFUNC2,CFUNC3,CFUNC4,CFUNC5,CFUNC6
 REAL(8) :: x,xm,xp 
 REAL(8) :: H0,EPS,ES ! --> task/lib/libde.f90  
 CHARACTER(80):: LINE
 INTEGER :: ILST
  la=l
  H0=0.5
  EPS=1.D-4 !!!??
  ILST=0

   WRITE(LINE,'(3I5,1P5E12.4)') j,l,NTHv,p0N,RHO0
!   WRITE(6,'(3I5,1P5E12.4)') j,l,NTHv,p0N,RHO0,x,xm,xp
!  IF(RHO0.EQ.2.58D-1) ILST=1
!ILST=1
   IF(mod(j,2).EQ.1) THEN
    s=1
    n=5-(j-1)/2
   ELSE
    s=2
    n=6-(j-2)/2
   ENDIF
   I2Cp0N=2*CI*AIMAG(Cp0N)
    
   CALL DEFTC(CINT,ES,H0,EPS,ILST,CFUNCS,KID='CFUNCS'//TRIM(LINE))

 CONTAINS

 FUNCTION CFUNCS(x,xm,xp)
  USE plcomm
  USE pllocal
  USE wmfa
  INCLUDE '../dp/dpcomm.inc'
  REAL(8),INTENT(IN) :: x,xm,xp
  REAL(8) :: delp0
  REAL(8) :: DFS
  COMPLEX(8) :: CFUNCS

!  WRITE(6,'(1I5,1P5E12.4)') NTHv,p0N,RHO0,x,xm,xp
   IF(la.EQ.1) THEN
    paN=0.5D0*(1.D0+x)*p0N*(1.D0+1.D-10) 
    A0=0.5D0*p0N
   ELSEIF(la.EQ.2) THEN
    paN=0.5D0*(PmaxN+p0N*(1.D0+1.D-10))+0.5D0*(PmaxN-p0N*(1.D0+1.D-10))*x 
    A0=0.5D0*(PmaxN+p0N)
   ELSEIF(la.EQ.3) THEN
    paN=0.5D0*(1.D0+x)*PmaxN
    A0=0.5D0*PmaxN
   ELSE
    RETURN
   ENDIF
    CALL SPL2DD(paN,RHO0,fax,DFPax,DFRax,PMa0,RHOa0,&
         US(1:4,1:4,1:NPM,1:NRM,NTHv),NPM,NPMAX+2,NRMAXFP+2,IERR)
!   WRITE(6,'(1I5,1P6E12.4)') la,p0N,paN,RHO0,fax,DFPax,DFRax
   IF(s.EQ.1) THEN
    DFS=DFPax
   ELSE
    DFS=DFRax
   ENDIF

     delp0=ABS(p0N-paN)/p0N
   IF(delp0.LT.1.D-2) THEN
     IF(la.EQ.1) THEN
      CFUNCS=A0*DFS*(paN**n)*2*Cp0N/(I2Cp0N+p0N*xm+1.D-10*p0N*xp)
     ELSEIF(la.EQ.2) THEN
      CFUNCS=A0*DFS*(paN**n)*2*Cp0N/(I2Cp0N+p0N*xp-PmaxN*xp-1.D-10*p0N*xm)
     ELSE
      CFUNCS=A0*DFS*(paN**n)*2*Cp0N/(I2Cp0N+p0N*xm+1.D-10*p0N*xp)
     ENDIF
   ELSE
      CFUNCS=A0*DFS*(paN**n)/(1.D0-paN/Cp0N)
   ENDIF
! WRITE(6,'(1I5,1P6E12.4)') NTHv,p0N,x,xm,xp,CFUNCS
! WRITE(6,'(1P5E12.4)') A0,DFS,paN,Cp0N
 END FUNCTION CFUNCS 

END SUBROUTINE PVINT
