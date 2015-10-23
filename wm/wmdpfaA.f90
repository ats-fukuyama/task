MODULE wmfa
 IMPLICIT NONE
 REAL(8) ::paN
! REAL(8) :: x,xm,xp
 REAL(8) :: DFPax,DFRax,fax
 REAL(8) ::RHO0,p0N
 REAL(8) :: A0
 INTEGER :: la,NTHv
 COMPLEX(8) :: Cp0N
 REAL(8),PARAMETER :: PmaxN=7.D0

END MODULE wmfa
 
! ******************************************************
!                       WMDPFAA              
!      dielectric tensor used by arbitral distribution function 
!      without relativistic effects
! ******************************************************

SUBROUTINE WMDPFAA(CW,RHOWM,RKPR,CPM1,CPM2,CQM1,CQM2,CRM1,CRM2)
  USE plcomm
  USE pllocal
  USE bpsd_constants,ONLY : CI,PI,AMP
  USE wmfa
  INCLUDE '../dp/dpcomm.inc'

!  IMPLICIT NONE
!  INTEGER,INTENT(IN) :: NRWM
  COMPLEX(8),INTENT(IN) :: CW
  REAL(8),INTENT(IN) :: RHOWM,RKPR
  COMPLEX(8),INTENT(OUT) :: CPM1,CPM2,CQM1,CQM2,CRM1,CRM2
  REAL(8) :: p,v
  REAL(8) :: p0,v0,PTA
  REAL(8) :: RTA,VTA
  REAL(8) :: v0N ! **N -> normalized by Thermal (pN=p/PTA(=Pthermal))
  REAL(8) :: AM,Norml
  COMPLEX(8) :: CX,Cp0,Cv0N
  REAL(8) :: AP,AQ,AR
  REAL(8) :: fa1,DFPa1,DFRa1
  INTEGER :: NS,ID
  COMPLEX(8) :: CIPL1,CIPR1,CIPL2,CIPR2,CIP1,CIP2
  COMPLEX(8) :: CIQL1,CIQR1,CIQL2,CIQR2,CIQ1,CIQ2
  COMPLEX(8) :: CIRL1,CIRR1,CIRL2,CIRR2,CIR1,CIR2

     NS=3
     RHO0=RHOWM
     AM=AMP*PA(NS)
     RTA=RTFP0(NS)
!     WRITE(6,'(1P6E12.4)') RTA
! RETURN 
     VTA=SQRT(2.D0*RTA/AM)
     PTA=VTA*AM
     CX=ABS(RKPR)/(CW*AM)
     CPM1=(0.D0,0.D0)
     CPM2=(0.D0,0.D0)
     CQM1=(0.D0,0.D0)
     CQM2=(0.D0,0.D0)
     CRM1=(0.D0,0.D0)
     CRM2=(0.D0,0.D0)

!--
     DO NTH=1,NTHMAX
!       IF(ABS(TCSM(NTH)).LT.1.D-21) THEN
!        p0=1.D2*PTA
!        Cp0=1.D2*PTA+1.D21*AIMAG(1.D0/CX)  !!??
!       ELSE
        p0=1.D0/(REAL(CX)*TCSM(NTH))
        Cp0=1.D0/(CX*TCSM(NTH))
!       ENDIF
        v0=p0/AM
        v0N=v0/VTA
        p0N=p0/PTA
        Cp0N=Cp0/PTA
        Cv0N=Cp0N/AM
        NTHv=NTH

        AP=TSNM(NTH)*((1.D0+TCSM(NTH)*TCSM(NTH))**2)
        AQ=TSNM(NTH)*TCSM(NTH)*(1.D0+TCSM(NTH)*TCSM(NTH))
        AR=TSNM(NTH)*TCSM(NTH)*TCSM(NTH)
  
!-----PRICIPAL VALUE
       IF(p0N.LT.PmaxN) THEN
         CALL PVINT(1,1,Cp0N,RHOWM,CIPL1)
         CALL PVINT(1,2,Cp0N,RHOWM,CIPR1)
          CPM1=CPM1 + AP*(CIPL1+CIPR1)
         CALL PVINT(2,1,Cp0N,RHOWM,CIPL2)
         CALL PVINT(2,2,Cp0N,RHOWM,CIPR2)
          CPM2=CPM2 + AP*(CIPL2+CIPR2)

         CALL PVINT(3,1,Cp0N,RHOWM,CIQL1)
         CALL PVINT(3,2,Cp0N,RHOWM,CIQR1)
          CQM1=CQM1 + AQ*(CIQL1+CIQR1)
         CALL PVINT(4,1,Cp0N,RHOWM,CIQL2)
         CALL PVINT(4,2,Cp0N,RHOWM,CIQR2)
          CQM2=CQM2 + AQ*(CIQL2+CIQR2)

         CALL PVINT(5,1,Cp0N,RHOWM,CIRL1)
         CALL PVINT(5,2,Cp0N,RHOWM,CIRR1)
          CRM1=CRM1 + AR*(CIRL1+CIRR1)
         CALL PVINT(6,1,Cp0N,RHOWM,CIRL2)
         CALL PVINT(6,2,Cp0N,RHOWM,CIRR2)
          CRM2=CRM2 + AR*(CIRL2+CIRR2)
       ELSE
         CALL PVINT(1,3,Cp0N,RHOWM,CIP1)
          CPM1=CPM1 + AP*CIP1
         CALL PVINT(2,3,Cp0N,RHOWM,CIP2)
          CPM1=CPM2 + AP*CIP2

         CALL PVINT(3,3,Cp0N,RHOWM,CIQ1)
          CQM1=CQM1 + AQ*CIQ1
         CALL PVINT(4,3,Cp0N,RHOWM,CIQ2)
          CQM2=CQM2 + AQ*CIQ2

         CALL PVINT(5,3,Cp0N,RHOWM,CIR1)
          CRM1=CRM1 + AR*CIR1
         CALL PVINT(6,3,Cp0N,RHOWM,CIR2)
          CRM2=CRM2 + AR*CIR2
       ENDIF
!------END PRINCIPAL VALUE

!-----SINGULAR POINT
        CALL SPL2DD(p0N,RHOWM,fa1,DFPa1,DFRa1,PMa0,RHOa0,&
                    US(1:4,1:4,1:NPM,1:NRM,NTHv),NPM,NPMAX,NRMAX,IERR)
        IF(AIMAG(CW).GT.0.D0) THEN
         CPM1=CPM1
         CPM2=CPM2
         CQM1=CQM1
         CQM2=CQM2
         CRM1=CRM1
         CRM2=CRM2
        ELSEIF(AIMAG(CW).EQ.0.D0) THEN
         CPM1=CPM1-CI*PI*(DFPa1*(v0N**5)*AP)*p0N !!????
         CPM2=CPM2-CI*PI*(DFRa1*(v0N**6)*AP)*p0N
         CQM1=CQM1-CI*PI*(DFPa1*(v0N**4)*AQ)*p0N
         CQM2=CQM2-CI*PI*(DFRa1*(v0N**5)*AQ)*p0N
         CRM1=CRM1-CI*PI*(DFPa1*(v0N**3)*AR)*p0N
         CRM2=CRM2-CI*PI*(DFRa1*(v0N**4)*AR)*p0N
        ELSE  !(AIMAG(CW).LT.0.D0)
         CPM1=CPM1-2.D0*CI*PI*(DFPa1*(Cv0N**5)*AP)*Cp0N
         CPM2=CPM2-2.D0*CI*PI*(DFRa1*(Cv0N**6)*AP)*Cp0N
         CQM1=CQM1-2.D0*CI*PI*(DFPa1*(Cv0N**4)*AQ)*Cp0N
         CQM2=CQM2-2.D0*CI*PI*(DFRa1*(Cv0N**5)*AQ)*Cp0N
         CRM1=CRM1-2.D0*CI*PI*(DFPa1*(Cv0N**3)*AR)*Cp0N
         CRM2=CRM2-2.D0*CI*PI*(DFRa1*(Cv0N**4)*AR)*Cp0N
        ENDIF
!-----END SINGULAR POINT
     END DO

     CPM1=CPM1*DELTH 
     CPM2=CPM2*DELTH
     CQM1=CQM1*DELTH
     CQM2=CQM2*DELTH
     CRM1=CRM1*DELTH
     CRM2=CRM2*DELTH

!  WRITE(6,'(1P6E12.4)') CPM1,CPM2
!     CPM1=(0.D0,0.D0)
!     CPM2=(0.D0,0.D0)
!     CQM1=(0.D0,0.D0)
!     CQM2=(0.D0,0.D0)
!     CRM1=(0.D0,0.D0)
!     CRM2=(0.D0,0.D0)  
!    RETURN
END SUBROUTINE WMDPFAA


SUBROUTINE PVINT(j,l,Cp0a,RHOWM,CINT)
 USE plcomm
 USE pllocal
 USE bpsd_constants,ONLY : CI,PI,AMP
 USE wmfa
 USE libde,ONLY : DEFTC
 INCLUDE '../dp/dpcomm.inc'
! IMPLICIT NONE
 INTEGER,INTENT(IN) :: j  ! j-> (df/dp or df/dr) and power 
 INTEGER,INTENT(IN) :: l  ! l=1-> boundary (0,p0N), l=2-> (p0N,pmaxN), l=3-> (0,pmaxN)
 REAL(8),INTENT(IN) :: RHOWM
 COMPLEX(8),INTENT(IN) :: Cp0a
 COMPLEX(8),INTENT(OUT) :: CINT ! Integral
! COMPLEX(8), :: CFUNC1,CFUNC2,CFUNC3,CFUNC4,CFUNC5,CFUNC6
 REAL(8) :: x,xm,xp 
 REAL(8) :: H0,EPS,ES ! --> task/lib/libde.f90  
 INTEGER :: ILST
  la=l
  H0=0.5
  EPS=1.D-8
  ILST=0 
 
  IF(j.EQ.1) THEN
   CALL DEFTC(CINT,ES,H0,EPS,ILST,CFUNC1,KID='CFUNC1')
  ELSEIF(j.EQ.2) THEN
   CALL DEFTC(CINT,ES,H0,EPS,ILST,CFUNC2,KID='CFUNC2')
  ELSEIF(j.EQ.3) THEN
   CALL DEFTC(CINT,ES,H0,EPS,ILST,CFUNC3,KID='CFUNC3')
  ELSEIF(j.EQ.4) THEN
   CALL DEFTC(CINT,ES,H0,EPS,ILST,CFUNC4,KID='CFUNC4')
  ELSEIF(j.EQ.5) THEN
   CALL DEFTC(CINT,ES,H0,EPS,ILST,CFUNC5,KID='CFUNC5')
  ELSEIF(j.EQ.6) THEN
   CALL DEFTC(CINT,ES,H0,EPS,ILST,CFUNC6,KID='CFUNC6')
  ELSE
   RETURN
  ENDIF

 CONTAINS

  SUBROUTINE PXOUT(lo,xo,paNo,Ao)
   USE wmfa
   IMPLICIT NONE
   INTEGER,INTENT(IN) :: lo ! l=1-> boundary (0,p0N), l=2-> (p0N,pmaxN), l=3-> (0,pmaxN)
   REAL(8),INTENT(IN) :: xo
   REAL(8),INTENT(OUT) :: paNo !paN=p(x)
   REAL(8),INTENT(OUT) :: Ao ! A -> factor dp=Adx
    
   IF(lo.EQ.1) THEN
    paNo=0.5D0*(1.D0+xo)*p0N    
    Ao=0.5D0*p0N
   ELSEIF(lo.EQ.2) THEN
    paNo=0.5D0*(PmaxN+p0N)+0.5D0*(PmaxN-p0N)*xo
    Ao=0.5D0*(PmaxN+p0N)
   ELSEIF(lo.EQ.3) THEN
    paNo=0.5D0*(1.D0+xo)*PmaxN
    Ao=0.5D0*PmaxN
   ELSE
    RETURN
   ENDIF
  END SUBROUTINE PXOUT

  FUNCTION CFUNC1(x,xm,xp) 
    USE plcomm
    USE pllocal
    USE wmfa
    INCLUDE '../dp/dpcomm.inc'
    REAL(8),INTENT(IN) :: x,xm,xp
    COMPLEX(8) :: CFUNC1  
    CALL PXOUT(la,x,paN,A0)
    CALL SPL2DD(paN,RHO0,fax,DFPax,DFRax,PMa0,RHOa0,&
         US(1:4,1:4,1:NPM,1:NRM,NTHv),NPM,NPMAX,NRMAX,IERR)
    CFUNC1=A0*DFPax*(paN**5)/(1.D0-paN/Cp0N)
  END FUNCTION CFUNC1

  FUNCTION CFUNC2(x,xm,xp)
    USE plcomm
    USE pllocal
    USE wmfa
    INCLUDE '../dp/dpcomm.inc'
    REAL(8),INTENT(IN) :: x,xm,xp
    COMPLEX(8) :: CFUNC2
    CALL PXOUT(la,x,paN,A0)
    CALL SPL2DD(paN,RHO0,fax,DFPax,DFRax,PMa0,RHOa0,&
         US(1:4,1:4,1:NPM,1:NRM,NTHv),NPM,NPMAX,NRMAX,IERR)
    CFUNC2=A0*DFRax*(paN**6)/(1.D0-paN/Cp0N)
  END FUNCTION CFUNC2

  FUNCTION CFUNC3(x,xm,xp)
    USE plcomm
    USE pllocal
    USE wmfa
    INCLUDE '../dp/dpcomm.inc'
    REAL(8),INTENT(IN) :: x,xm,xp
    COMPLEX(8) :: CFUNC3
    CALL PXOUT(la,x,paN,A0)
    CALL SPL2DD(paN,RHO0,fax,DFPax,DFRax,PMa0,RHOa0,&
         US(1:4,1:4,1:NPM,1:NRM,NTHv),NPM,NPMAX,NRMAX,IERR)
    CFUNC3=A0*DFPax*(paN**4)/(1.D0-paN/Cp0N)
  END FUNCTION CFUNC3

  FUNCTION CFUNC4(x,xm,xp)
    USE plcomm
    USE pllocal
    USE wmfa
    INCLUDE '../dp/dpcomm.inc'
    REAL(8),INTENT(IN) :: x,xm,xp
    COMPLEX(8) :: CFUNC4
    CALL PXOUT(la,x,paN,A0)
    CALL SPL2DD(paN,RHO0,fax,DFPax,DFRax,PMa0,RHOa0,&
         US(1:4,1:4,1:NPM,1:NRM,NTHv),NPM,NPMAX,NRMAX,IERR)
    CFUNC4=A0*DFRax*(paN**5)/(1.D0-paN/Cp0N)
  END FUNCTION CFUNC4

  FUNCTION CFUNC5(x,xm,xp)
    USE plcomm
    USE pllocal
    USE wmfa
    INCLUDE '../dp/dpcomm.inc'
    REAL(8),INTENT(IN) :: x,xm,xp
    COMPLEX(8) :: CFUNC5
    CALL PXOUT(la,x,paN,A0)
    CALL SPL2DD(paN,RHO0,fax,DFPax,DFRax,PMa0,RHOa0,&
         US(1:4,1:4,1:NPM,1:NRM,NTHv),NPM,NPMAX,NRMAX,IERR)
    CFUNC5=A0*DFPax*(paN**3)/(1.D0-paN/Cp0N)
  END FUNCTION CFUNC5

  FUNCTION CFUNC6(x,xm,xp)
    USE plcomm
    USE pllocal
    USE wmfa
    INCLUDE '../dp/dpcomm.inc'
    REAL(8),INTENT(IN) :: x,xm,xp
    COMPLEX(8) :: CFUNC6
    CALL PXOUT(la,x,paN,A0)
    CALL SPL2DD(paN,RHO0,fax,DFPax,DFRax,PMa0,RHOa0,&
         US(1:4,1:4,1:NPM,1:NRM,NTHv),NPM,NPMAX,NRMAX,IERR)
    CFUNC6=A0*DFRax*(paN**4)/(1.D0-paN/Cp0N)
  END FUNCTION CFUNC6
  
END SUBROUTINE PVINT
