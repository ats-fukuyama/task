MODULE wmfa
 IMPLICIT NONE
 REAL(8) ::paN,paNL,paNR  ! **L => left of pole, **R => right of pole
 REAL(8) :: DFPax,DFRax,fax
 REAL(8) :: DFPaxR,DFRaxR,faxR,DFPaxL,DFRaxL,faxL
 REAL(8) ::RHO0,p0N,v0N ! **N -> normalized by Thermal (pN=p/PTA(=Pthermal)) 
 REAL(8) :: A0
 REAL(8) :: pintN,LBi1,RBi1,LBi2,RBi2,LBi3,RBi3
 REAL(8) :: LBint,RBint
 INTEGER :: la,n,s,NTHv
 COMPLEX(8) :: Cp0N,Cv0N,ImCp0N
 COMPLEX(8) :: CFfa
 REAL(8) :: PmaxN
! REAL(8),PARAMETER :: PmaxN=5.D0
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
  USE libgrf, ONLY : grd1d,grd2d
  USE dpcomm

!  IMPLICIT NONE
!  INTEGER,INTENT(IN) :: NRWM
  COMPLEX(8),INTENT(IN) :: CW
  REAL(8),INTENT(IN) :: RHOWM,RKPR
  REAL(8),INTENT(OUT) ::AE2N0
  COMPLEX(8),INTENT(OUT) :: CPM1,CPM2,CQM1,CQM2,CRM1,CRM2
  REAL(8) :: p,v
  REAL(8) :: p0,v0,PTA
  REAL(8) :: RTA,VTA
  REAL(8) :: AM,Norml,fpr1,fprnp
  REAL(8) :: RHOM_MIN,RHOM_MAX
  COMPLEX(8) :: CX,Cp0
  REAL(8) :: AP,AQ,AR
  REAL(8) :: fa1,DFPa1,DFRa1
  INTEGER :: NS,ID,NRWM
  COMPLEX(8) :: CIPL1,CIPR1,CIPL2,CIPR2,CIP1,CIP2
  COMPLEX(8) :: CIQL1,CIQR1,CIQL2,CIQR2,CIQ1,CIQ2
  COMPLEX(8) :: CIRL1,CIRR1,CIRL2,CIRR2,CIR1,CIR2
!
  REAL(8) :: xg(1:120),yg(1:120),zg(1:120,1:120)
  REAL(8) :: SUMfa!(1:NRMAXFP+2) !!!!
!
     NS=3
     PmaxN=5.D0
     RHO0=RHOWM
     RHOM_MIN=RHOa0(1)
     RHOM_MAX=RHOa0(NRMAXFP+2)
     AM=AMFP(NS)
     AE=AEFP(NS)
     AE2N0=AE*AE*RNFP0(NS)
     RTA=RTFP0(3)*AEE*1.D3 !  T_thermal difinition ????
     VTA=SQRT(RTA/AM) ! VTA=SQRT(2.D0*RTA/AM) ???
     PTA=VTA*AM
     CFfa=CW*1.D-6/(2.D0*PI) !???
     CX=ABS(RKPR)/(CW*AM)
     CPM1=(0.D0,0.D0)
     CPM2=(0.D0,0.D0)
     CQM1=(0.D0,0.D0)
     CQM2=(0.D0,0.D0)
     CRM1=(0.D0,0.D0)
     CRM2=(0.D0,0.D0)

!  WRITE(6,'(1P5E12.4)') CFfa!RTA,VTA,PTA,DELTH

! distribution function from FP --------------------------
! DO NNP=1,100
!  xg(NNP)=(400.D0+REAL(NNP)-0.5D0)*5.D0/500.D0 !PMa0(NNP+1)
! DO NNR=1,100
!  yg(NNR)=(400.D0+REAL(NNR)-0.5D0)/500.D0 !RHOa0(NNR+1)
!     CALL SPL2DD(xg(NNP),yg(NNR),fa1,DFPa1,DFRa1,PMa0,RHOa0,&
!                 US(1:4,1:4,1:NPM,1:NRM,16),NPM,NPMAX+2,NRMAXFP+2,IERR) 
!!  yg(NNP)=fa1
!  zg(NNP,NNR)=fa1 !fa0(NNP+1,NNR+1,16) 
!  WRITE(6,'(1P5E12.4)') xg(NNP),yg(NNR),zg(NNP,NNR)
! ENDDO
! ENDDO
 
! CALL PAGES
!!  CALL GRD1D(0,xg,yg,120,120,1,'fa1')
! CALL GRD2D(0,xg(1:120),yg(1:120),zg(1:120,1:120),120,120,120,'fa1',MODE_2d=12)
! CALL PAGEE 
 
! RETURN
! ---------------------------
!  IF(RHOWM.GT.9.D-1) THEN
!   DO NR=1,NRMAXFP+2
!      IF(RHOWM-RHOa0(NR).LE.0.D0) THEN
!        NRWM=NR-1
!        EXIT
!      ENDIF
!   ENDDO
!   DO NTH=1,NTHMAXFP
!      fpr1=fpr1+fa0(1,NRWM,NTH)*TSNM(NTH)*DELTH
!   ENDDO

!   DO NR=1,NRMAXFP+2
!      SUMfa=0.D0
!   DO NP=1,NPMAX+1
!   DO NTH=1,NTHMAXFP
!       fprnp=0.D0
!       fprnp=0.5D0*(fa0(NP,NRWM,NTH)+fa0(NP+1,NRWM,NTH))*TSNM(NTH)*DELTH
!      SUMfa=SUMfa&
!      +0.5D0*(PMa0(NP+1)**2+PMa0(NP)**2)*fprnp
!      WRITE(6,'(3I5,1P5E12.4)') NR,NP,NTH,fa0(NP+1,NR,NTH)
!   ENDDO
!      IF((fprnp.LT.fpr1*1.D-8).OR.(fprnp.LT.0.D0)) THEN 
!        PmaxN=(NP-0.5D0)*DELP(3)
!        IF(PmaxN.LT.1.D0) RETURN
!        EXIT
!      ENDIF
!   ENDDO   
!     SUMfa=RNFP0(3)*SUMfa*2*PI*DELP(3) ! *VTA**3 ? not need? 
!     WRITE(6,'(1I5,1P5E12.4)') NR,SUMfa
!   ENDDO

!  ENDIF

! RETURN                       !!----

!--
     DO NTH=1,NTHMAXFP
        p0N=1.D0/(PTA*REAL(CX)*TCSM(NTH))
        Cp0N=1.D0/(PTA*CX*TCSM(NTH))

        v0N=p0N*PTA/(AM*VTA)
        Cv0N=Cp0N*PTA/(AM*VTA)
        NTHv=NTH
!        WRITE(6,'(A,1I5,1P5E12.4)') 'NTH,Cp0N=',NTH,Cp0N,CX,TCSM(NTH)

        AP=TSNM(NTH)*((1.D0+TCSM(NTH)*TCSM(NTH))**2)
        AQ=TSNM(NTH)*TCSM(NTH)*(1.D0+TCSM(NTH)*TCSM(NTH))
        AR=TSNM(NTH)*TCSM(NTH)*TCSM(NTH)
!-----PRICIPAL VALUE

        IF((p0N.LT.0.5D0*DELP(3)).OR.(p0N.GE.PmaxN-0.5D0*DELP(3))) THEN
           LBi3=DELP(3)
           RBi3=PmaxN-DELP(3)
         CALL PVINT(1,3,CIP1)   ! 3 => INTEGRAL 0.14 < pN < 4.86
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
        IF(p0N.LE.0.5D0*PmaxN) THEN 
           ! Integral 1:(0,p0N)+(p0N,2p0N) 2:(2p0N,PmaxN) 
           LBi1=p0N
           RBi1=2*p0N
           LBi2=2*p0N
           RBi2=PmaxN !-5.D-2
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

        ELSE
           ! Integral 1:(PmaxN-2*(PmaxN-p0N),p0N)+(p0N,PmaxN) 2:(0,Pmax-2*(PmaxN-p0N)) 
           LBi1=p0N
           RBi1=PmaxN
           LBi2=0.D0
           RBi2=PmaxN-2*(PmaxN-p0N)
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
        ENDIF
!         WRITE(6,'(1P8E12.4)') CIRL1,CIRR1,CIRL2,CIRR2
       ENDIF
!------END PRINCIPAL VALUE

!-----SINGULAR POINT
        IF((p0N.GT.0.D0).AND.(p0N.LT.PmaxN)) THEN
         CALL SPL2DD(p0N,RHOWM,fa1,DFPa1,DFRa1,PMa0,RHOa0,&
                    US(1:4,1:4,1:NPM,1:NRM,NTHv),NPM,NPMAX+2,NRMAXFP+2,IERR)
!        WRITE(6,'(1I5,1P6E12.4)') NTHv,p0N,RHO0,fa1,DFPa1,DFRa1
         IF(AIMAG(Cp0N).GT.0.D0) THEN
          CPM1=CPM1
          CPM2=CPM2
          CQM1=CQM1
          CQM2=CQM2
          CRM1=CRM1
          CRM2=CRM2
         ELSEIF(AIMAG(Cp0N).EQ.0.D0) THEN
          CPM1=CPM1-CI*PI*(DFPa1*((v0N)**5)*AP)*p0N !!????
          CPM2=CPM2-CI*PI*(DFRa1*((v0N)**6)*AP)*p0N
          CQM1=CQM1-CI*PI*(DFPa1*((v0N)**4)*AQ)*p0N
          CQM2=CQM2-CI*PI*(DFRa1*((v0N)**5)*AQ)*p0N
          CRM1=CRM1-CI*PI*(DFPa1*((v0N)**3)*AR)*p0N
          CRM2=CRM2-CI*PI*(DFRa1*((v0N)**4)*AR)*p0N
         ELSE 
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

     CPM1=CPM1*DELTH*(VTA**5)/(VTA**3)   ! dp=p_0*dp_n -> df/dp=(df/dp_n)*dp_n/dp=dfp_n/p_0
     CPM2=CPM2*DELTH*(VTA**6)*PTA/(VTA**3) ! ---> here may be wrong
     CQM1=CQM1*DELTH*(VTA**4)/(VTA**3)
     CQM2=CQM2*DELTH*(VTA**5)*PTA/(VTA**3)
     CRM1=CRM1*DELTH*(VTA**3)/(VTA**3)
     CRM2=CRM2*DELTH*(VTA**4)*PTA/(VTA**3)

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
! USE libde,ONLY : DEFTC
 USE libgrf
 USE dpcomm
! IMPLICIT NONE
 INTEGER,INTENT(IN) :: j  ! j-> (df/dp or df/dr) and power 
 INTEGER,INTENT(IN) :: l  ! l -> Integral boundary 
 COMPLEX(8),INTENT(OUT) :: CINT ! Integral
 REAL(8) :: x,xm,xp 
 REAL(8) :: H0,EPS,ES ! --> task/lib/libde.f90  
 CHARACTER(80):: LINE
 INTEGER :: ILST
 INTEGER:: m
 REAL(8):: xg(0:1000),yg(0:1000)
  H0=0.5
  EPS=1.D-6 !!!??
  ILST=0

!   WRITE(6,'(1I5,1P5E12.4)') NTHv, RHO0, CWfa
   WRITE(LINE,'(1I5,1P5E12.4)') NTHv,RHO0,CFfa,Cp0N
!   WRITE(LINE,'(1P4E12.4)') j.l,Cp0N,LBint,RBint   
!   IF(RHO0.GE.9.84D-1) ILST=1
   la=l
   IF(mod(j,2).EQ.1) THEN
    s=1
    n=5-(j-1)/2
   ELSE
    s=2
    n=6-(j-2)/2
   ENDIF

   IF(l.EQ.1) THEN
     LBint=LBi1
     RBint=RBi1
     ImCp0N=Cp0N-p0N

     IF(s.EQ.1) THEN
      CALL DEFTCFA(CINT,ES,H0,EPS,ILST,CFUNCpldp,KID='CFUNCpldp'//TRIM(LINE))
     ELSE
      CALL DEFTCFA(CINT,ES,H0,EPS,ILST,CFUNCpldr,KID='CFUNCpldr'//TRIM(LINE))
     ENDIF

!      IF(ILST.EQ.1) THEN
!      WRITE(6,'(1I5,1P5E12.4)') NTHv, RHO0, CWfa, LBint, RBint
!      do m=0,1000
!        xg(m)=-9.999999999D-1+dble(m)/2000
!        yg(m)=CFUNCpldp(xg(m),1.D0-xg(m),1.D0+xg(m))
!      end do
!      CALL pages
!      CALL GRD1D(0,xg,yg,1001,1001,1,'CFUNCpldp')
!      CALL pagee
!      ENDIF

   ELSEIF(l.EQ.2) THEN
     LBint=LBi2
     RBint=RBi2

     IF(s.EQ.1) THEN
      CALL DEFTCFA(CINT,ES,H0,EPS,ILST,CFUNCdp,KID='CFUNCdp'//TRIM(LINE))
     ELSE
      CALL DEFTCFA(CINT,ES,H0,EPS,ILST,CFUNCdr,KID='CFUNCdr'//TRIM(LINE))
     ENDIF

!     IF(ILST.EQ.1) THEN
!      WRITE(6,'(1I5,1P5E12.4)') NTHv, RHO0, CWfa, LBint, RBint
!      do m=0,1000
!        xg(m)=-9.999999999D-1+dble(m)/2000
!        yg(m)=CFUNCdp(xg(m),1.D0-xg(m),1.D0+xg(m))
!      end do
!      CALL pages
!      CALL GRD1D(0,xg,yg,1001,1001,1,'CFUNCdp')
!      CALL pagee
!     ENDIF

   ELSEIF(l.EQ.3) THEN
!      ILST=1
     LBint=LBi3
     RBint=RBi3

     IF(s.EQ.1) THEN
      CALL DEFTCFA(CINT,ES,H0,EPS,ILST,CFUNCdp,KID='CFUNCdp'//TRIM(LINE))
     ELSE
      CALL DEFTCFA(CINT,ES,H0,EPS,ILST,CFUNCdr,KID='CFUNCdr'//TRIM(LINE))
     ENDIF

!     IF(ILST.GT.1) THEN
!      WRITE(6,'(1I5,1P5E12.4)') NTHv, RHO0, CWfa, LBint, RBint
!      do m=0,1000
!        xg(m)=-9.999999999D-1+dble(m)/2000
!        yg(m)=CFUNCdp(xg(m),1.D0-xg(m),1.D0+xg(m))
!      end do
!      CALL pages
!      CALL GRD1D(0,xg,yg,1001,1001,1,'CFUNCdp')
!      CALL pagee
!     ENDIF

   ELSE
    RETURN
   ENDIF

!===Check convergence
!     IF(ILST.NE.0) THEN       
!      WRITE(6,'(1I5,1P5E12.4)') NTHv, RHO0, CWfa, LBint, RBint
!      do m=0,1000
!        xg(m)=-9.999999999D-1+dble(m)/500
!        yg(m)=CFUNCdr(xg(m),1.D0-xg(m),1.D0+xg(m))
!      end do
!      CALL pages
!      CALL GRD1D(0,xg,yg,1001,1001,1,'CFUNCdr')
!      CALL pagee
!     ENDIF
!===Check convergence

 CONTAINS

 FUNCTION CFUNCdp(x,xm,xp)
  USE plcomm
  USE pllocal
  USE wmfa
  USE dpcomm
  REAL(8),INTENT(IN) :: x,xm,xp
  COMPLEX(8) :: CFUNCdp

    paN=0.5D0*(RBint+LBint)+0.5D0*(RBint-LBint)*x
    A0 =0.5D0*(RBint-LBint)

    CALL SPL2DD(paN,RHO0,fax,DFPax,DFRax,PMa0,RHOa0,&
         US(1:4,1:4,1:NPM,1:NRM,NTHv),NPM,NPMAX+2,NRMAXFP+2,IERR)

      CFUNCdp=A0*DFPax*(paN**n)/(1-paN/Cp0N)
 END FUNCTION CFUNCdp

 FUNCTION CFUNCdr(x,xm,xp)
  USE plcomm
  USE pllocal
  USE wmfa
  USE dpcomm
  REAL(8),INTENT(IN) :: x,xm,xp
  COMPLEX(8) :: CFUNCdr

    paN=0.5D0*(RBint+LBint)+0.5D0*(RBint-LBint)*x
    A0 =0.5D0*(RBint-LBint)

    CALL SPL2DD(paN,RHO0,fax,DFPax,DFRax,PMa0,RHOa0,&
         US(1:4,1:4,1:NPM,1:NRM,NTHv),NPM,NPMAX+2,NRMAXFP+2,IERR)

      CFUNCdr=A0*DFRax*(paN**n)/(1-paN/Cp0N)
 END FUNCTION CFUNCdr

 FUNCTION CFUNCpldp(x,xm,xp)
  USE plcomm
  USE pllocal
  USE wmfa
  use dpcomm
  REAL(8),INTENT(IN) :: x,xm,xp
  REAL(8) :: FL,FR
  COMPLEX(8) :: CFUNCpldp

    paNL=LBint-0.5D0*(RBint-LBint)*xp
    paNR=LBint+0.5D0*(RBint-LBint)*xp 
    A0  =0.5D0*(RBint-LBint)

    CALL SPL2DD(paNL,RHO0,faxL,DFPaxL,DFRaxL,PMa0,RHOa0,&
         US(1:4,1:4,1:NPM,1:NRM,NTHv),NPM,NPMAX+2,NRMAXFP+2,IERR)
    CALL SPL2DD(paNR,RHO0,faxR,DFPaxR,DFRaxR,PMa0,RHOa0,&
         US(1:4,1:4,1:NPM,1:NRM,NTHv),NPM,NPMAX+2,NRMAXFP+2,IERR)  

    FL=(paNL**n)*DFPaxL
    FR=(paNR**n)*DFPaxR

      CFUNCpldp=A0*Cp0N*(FL/(ImCp0N+A0*xp)+FR/(ImCp0N-A0*xp))
!      CFUNCpldp=A0*Cp0N*(ImCp0N*(FR+FL)+A0*xp*(FR-FL))/(ImCp0N**2-(A0*xp)**2)
! WRITE(6,'(2I5,1P10E12.4)') NTHv,la,x,xm,xp,paNL,paNR,p0N
! WRITE(6,'(1P10E12.4)') A0,FL,FR,Cp0N,ImCp0N,CFUNCpldp 
 END FUNCTION CFUNCpldp

 FUNCTION CFUNCpldr(x,xm,xp)
  USE plcomm
  USE pllocal
  USE wmfa
  USE dpcomm
  REAL(8),INTENT(IN) :: x,xm,xp
  REAL(8) :: FL,FR
  COMPLEX(8) :: CFUNCpldr

    paNL=LBint-0.5D0*(RBint-LBint)*xp
    paNR=LBint+0.5D0*(RBint-LBint)*xp 
    A0  =0.5D0*(RBint-LBint)

    CALL SPL2DD(paNL,RHO0,faxL,DFPaxL,DFRaxL,PMa0,RHOa0,&
         US(1:4,1:4,1:NPM,1:NRM,NTHv),NPM,NPMAX+2,NRMAXFP+2,IERR)
    CALL SPL2DD(paNR,RHO0,faxR,DFPaxR,DFRaxR,PMa0,RHOa0,&
         US(1:4,1:4,1:NPM,1:NRM,NTHv),NPM,NPMAX+2,NRMAXFP+2,IERR)  

    FL=(paNL**n)*DFRaxL
    FR=(paNR**n)*DFRaxR

      CFUNCpldr=A0*Cp0N*(ImCp0N*(FR+FL)+A0*xp*(FR-FL))/(ImCp0N**2-(A0*xp)**2)

 END FUNCTION CFUNCpldr

!END SUBROUTINE PVINT

!     ************************************************************
!        FINITE INTEGRAL BY DOUBLE-EXPONENTIAL FORMULA
!                    (-1.D0, +1.D0)
!         INTEGRAND SHOULD BE DEFINED BY FUNC(X,1-X,1+X)
!     ************************************************************
!
!        For integral with respect to y from a to b,
!            variable transformation should be
!                 x = (2*y-a-b)/(b-a)
!
!                 near y ~ a   y-a=(b-a)*(1+x)/2
!                 near y ~ b   b-y=(b-a)*(1-x)/2
!                 otherwize      y=(b-a)*x/2+(a+b)/2

 SUBROUTINE DEFTCFA(CS,ES,H0,EPS,ILST,CFUNC,KID)
    IMPLICIT NONE
    COMPLEX(8),INTENT(OUT):: CS   ! Integral
    REAL(8),INTENT(OUT):: ES   ! Estimated error 
    REAL(8),INTENT(IN)::  H0   ! Initial step size
    REAL(8),INTENT(IN)::  EPS  ! Convergence thrshold
    INTEGER,INTENT(INOUT)::  ILST ! print out control: 0 for no print out
    INTERFACE
       FUNCTION CFUNC(X,XM,XP)
         COMPLEX(8):: CFUNC
         REAL(8),INTENT(IN):: X,XM,XP
       END FUNCTION CFUNC
    END INTERFACE
    CHARACTER(LEN=*),INTENT(IN),OPTIONAL:: KID   ! function identifier string
    REAL(8),PARAMETER:: HP=1.5707963267948966192D0

    COMPLEX(8):: CSI,CSP,CT
    REAL(8):: EPS1,H,X,ATP,ATM,HN,HC,HS,CC,XM,XP,AT
    INTEGER:: NP,NM,NPMIN,NMMIN,IND,NPD,NMD,NMAX

      EPS1=EPS**0.75D0
      H=H0
      X=0.D0
      CSI=HP*CFUNC(X,1.D0-X,1.D0+X)
      CS=H*CSI
      CSP=0.D0
      NP=0
      NM=0
      NPMIN=1
      NMMIN=1

    1 IND=0
      ATP=ABS(CSI)
      ATM=ATP
      NPD=2
      IF(NP.EQ.0) NPD=1
      NMD=2
      IF(NM.EQ.0) NMD=1

   10 IF(IND.NE.1) THEN
         IF(NP.EQ.NPMIN+2) NPD=1
         NP=NP+NPD
         HN=DBLE(NP)*H
         HC=HP*H*COSH(HN)
         HS=HP*SINH(-HN)
         X=TANH(HS)
         CC=1.D0/COSH(HS)
         XM=CC*EXP(-HS)
         XP=CC*EXP( HS)
         CT=HC*CFUNC(X,XM,XP)*CC*CC
         CS=CS+CT
         AT=ATP
         ATP=ABS(CT)/H
         IF(NP.GE.NPMIN) THEN
            IF(AT+ATP.LE.EPS1*ABS(CS)) THEN
               IF(IND.EQ.-1) GO TO 100
               IND=1
            ENDIF
         ENDIF
      ENDIF

      IF(IND.NE.-1) THEN
         IF(NM.EQ.NMMIN+2) NMD=1
         NM=NM+NMD
         HN=DBLE(NM)*H
         HC=HP*H*COSH(HN)
         HS=HP*SINH( HN)
         X=TANH(HS)
         CC=1.D0/COSH(HS)
         XM=CC*EXP(-HS)
         XP=CC*EXP( HS)
         CT=HC*CFUNC(X,XM,XP)*CC*CC
         CS=CS+CT
         AT=ATM
         ATM=ABS(CT)/H
         IF(NM.GE.NMMIN) THEN
            IF(AT+ATM.LE.EPS1*ABS(CS)) THEN
               IF(IND.EQ.1) GO TO 100
               IND=-1
            ENDIF
         ENDIF
      ENDIF
      GO TO 10

  100 ES=ABS(CS-CSP)
      IF(ILST.NE.0) THEN
         IF(H.GE.H0) WRITE(6,601) H,NP,NM,CS
         IF(H.LT.H0) WRITE(6,602) H,NP,NM,CS,ES
      ENDIF
      CSP=CS
      IF(ES.LE.EPS1*ABS(CS)) GO TO 200
      NMAX=MAX0(NP,NM)
      IF(NMAX.GT.1000) THEN
         WRITE(6,603) TRIM(KID)
         IF(ABS(CS).LT.1.D-20) THEN
            CS=(0.D0,0.D0)
            GO TO 200            
         ELSE
            GO TO 9999
         ENDIF
      ENDIF
      H=0.5D0*H
      CS=0.5D0*CS
      NPMIN=NP*2-1
      NMMIN=NM*2-1
      NP=-1
      NM=-1
      GO TO 1

  200 RETURN

  501 FORMAT(A1)
  601 FORMAT(1H ,1PD13.5,2I8,1PD24.15)
  602 FORMAT(1H ,1PD13.5,2I8,1PD24.15,1PD14.5)
  603 FORMAT(1H ,'XX DEFT: NMAX EXCEEDS 1000: FUNC=',A)
 9999 STOP
 END SUBROUTINE DEFTCFA

END SUBROUTINE PVINT
