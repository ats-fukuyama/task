C     $Id$
C
C   ************************************************  
C   **            Mesh Definition                 **
C   ************************************************
C
      SUBROUTINE EQMESH
C      
      INCLUDE 'eqcomc.h'
C
      EXTERNAL EQFBND
C
      DSG=1.D0/NSGMAX
      DTG=2.D0*PI/NTGMAX 
C
      DO NSG=1,NSGMAX
         SIGM(NSG)=DSG*(NSG-0.5D0)
      ENDDO
      DO NSG=1,NSGMAX+1
         SIGG(NSG)=DSG*(NSG-1.D0)
      ENDDO
      DO NTG=1,NTGMAX
         THGM(NTG)=DTG*(NTG-0.5D0)
      ENDDO
      DO NTG=1,NTGMAX+1
         THGG(NTG)=DTG*(NTG-1.D0)
      ENDDO
      RETURN
      END
C
C   ************************************************ 
C   **                Initial Psi                 **
C   ************************************************
C
      SUBROUTINE EQPSIN
C
      INCLUDE 'eqcomc.h'
C
      RAXIS=RR
      ZAXIS=0.D0
      PSI0=-0.5D0*RMU0*RIP*1.D6*RR
      DO NSG=1,NSGMAX
      DO NTG=1,NTGMAX
          PSI(NTG,NSG)=PSI0*(1-SIGM(NSG)*SIGM(NSG))
          DELPSI(NTG,NSG)=0.D0
      ENDDO
      ENDDO
      RETURN
      END
C
C   ************************************************
C   **               Iteration Loop               **
C   ************************************************
C
      SUBROUTINE EQLOOP(IERR)
C
      INCLUDE 'eqcomc.h'
C
      CALL EQDEFB
      DO NLOOP=1,20
         CALL EQBAND
         CALL EQRHSV(IERR)
         IF(IERR.NE.0) GOTO 200
         CALL EQSOLV
         SUM0=0.D0
         SUM1=0.D0
         DO NSG=1,NSGMAX
         DO NTG=1,NTGMAX
            SUM0=SUM0+PSI(NTG,NSG)**2
            SUM1=SUM1+DELPSI(NTG,NSG)**2
         ENDDO
         ENDDO
         SUM=SQRT(SUM1/SUM0)
         WRITE(6,'(A,1P4E14.6)')
     &        'SUM,RAXIS,ZAXIS,PSI0=',SUM,RAXIS,ZAXIS,PSI0
         IF(SUM.LT.EPSEQ)GOTO 100
      ENDDO
  100 CONTINUE 
      IERR=0
      RETURN
C
  200 CONTINUE
      IERR=200
      RETURN
      END
C
C   ************************************************  
C   **          Boundary Definition               **
C   ************************************************
C
      SUBROUTINE EQDEFB
C      
      INCLUDE 'eqcomc.h'
C
      DIMENSION DRHOM(NTGM),DRHOG(NTGMP)
      EXTERNAL EQFBND
C
      EPSZ=1.D-8
C
      DO NTG=1,NTGMAX
         ZBRF=TAN(THGM(NTG))
         THDASH=ZBRENT(EQFBND,THGM(NTG)-1.0D0,THGM(NTG)+1.0D0,EPSZ)
         RHOM(NTG)=RA*SQRT(COS(THDASH+RDLT*SIN(THDASH))**2
     &                    +RKAP**2*SIN(THDASH)**2)
      ENDDO
      DRHOM(1)=(RHOM(2)-RHOM(NTGMAX))/(2*DTG)
      DO NTG=2,NTGMAX-1
         DRHOM(NTG)=(RHOM(NTG+1)-RHOM(NTG-1))/(2*DTG)
      ENDDO
      DRHOM(NTGMAX)=(RHOM(1)-RHOM(NTGMAX-1))/(2*DTG)
C
      DO NTG=1,NTGMAX
         ZBRF=TAN(THGG(NTG))
         THDASH=ZBRENT(EQFBND,THGG(NTG)-1.0D0,THGG(NTG)+1.0D0,EPSZ)
         RHOG(NTG)=RA*SQRT(COS(THDASH+RDLT*SIN(THDASH))**2
     &                    +RKAP**2*SIN(THDASH)**2)
      ENDDO
      RHOG(NTGMAX+1)=RHOG(1)
      DRHOG(1)=(RHOG(2)-RHOG(NTGMAX))/(2*DTG)
      DO NTG=2,NTGMAX-1
         DRHOG(NTG)=(RHOG(NTG+1)-RHOG(NTG-1))/(2*DTG)
      ENDDO
      DRHOG(NTGMAX)=(RHOG(1)-RHOG(NTGMAX-1))/(2*DTG)
      DRHOG(NTGMAX+1)=DRHOG(1)
C
      DO NSG=1,NSGMAX
      DO NTG=1,NTGMAX+1
         RMG(NSG,NTG)=RR+SIGM(NSG)*RHOG(NTG)*COS(THGG(NTG))
      ENDDO
      ENDDO
      DO NSG=1,NSGMAX+1
      DO NTG=1,NTGMAX
         RGM(NSG,NTG)=RR+SIGG(NSG)*RHOM(NTG)*COS(THGM(NTG))       
      ENDDO
      ENDDO
C
      DO NSG=1,NSGMAX+1
      DO NTG=1,NTGMAX
         AA(NSG,NTG)=SIGG(NSG)/RGM(NSG,NTG)
     &              +(DRHOM(NTG)*DRHOM(NTG)*SIGG(NSG))
     &              /(RGM(NSG,NTG)*RHOM(NTG)*RHOM(NTG))  
         AB(NSG,NTG)=-DRHOM(NTG)/(RGM(NSG,NTG)*RHOM(NTG))
      ENDDO
      ENDDO
      DO NSG=1,NSGMAX
      DO NTG=1,NTGMAX+1
         AC(NSG,NTG)=-DRHOG(NTG)/(RMG(NSG,NTG)*RHOG(NTG))
         AD(NSG,NTG)=1/(RMG(NSG,NTG)*SIGM(NSG))
      ENDDO
      ENDDO
      RETURN
      END
C
C   ***********************************************
C   **            Matrix calculation             **
C   ***********************************************
C
      SUBROUTINE EQBAND
C
      INCLUDE 'eqcomc.h'
C
      MMAX=NSGMAX*NTGMAX
      NBND=2*NTGMAX
      DO N=1,MMAX
      DO M=1,2*NBND-1
          Q(M,N)=0.D0
      ENDDO
      ENDDO
C
      DO NSG=1,NSGMAX
      DO NTG=1,NTGMAX
         I=(NSG-1)*NTGMAX+NTG
         IF(NTG.EQ.1)THEN
            Q(NBND         -1,I)= (AB(NSG,NTG)+AC(NSG,NTG))
     &                           /(4*DSG*DTG)
            Q(NBND  +NTGMAX-1,I)= (AB(NSG,NTG)-AB(NSG+1,NTG))
     &                           /(4*DSG*DTG)
     &                           +AD(NSG,NTG)/(DTG*DTG)
            Q(NBND+2*NTGMAX-1,I)=-(AB(NSG+1,NTG)+AC(NSG,NTG))
     &                           /(4*DSG*DTG)
         ELSE
            Q(NBND  -NTGMAX-1,I)= (AB(NSG,NTG)+AC(NSG,NTG))
     &                           /(4*DSG*DTG)
            Q(NBND         -1,I)=(AB(NSG,NTG)-AB(NSG+1,NTG))
     &                           /(4*DSG*DTG)
     &                           +AD(NSG,NTG)/(DTG*DTG)
            Q(NBND  +NTGMAX-1,I)=-(AB(NSG+1,NTG)+AC(NSG,NTG))
     &                           /(4*DSG*DTG)
         ENDIF
         Q(NBND-NTGMAX,I)=  AA(NSG,NTG)/(DSG*DSG)
     &                   -(AC(NSG,NTG+1)-AC(NSG,NTG))/(4*DSG*DTG)
         Q(NBND       ,I)=-(AA(NSG+1,NTG)+AA(NSG,NTG))/(DSG*DSG)
     &                   -(AD(NSG,NTG+1)+AD(NSG,NTG))/(DTG*DTG)
         Q(NBND+NTGMAX,I)=  AA(NSG+1,NTG)/(DSG*DSG)
     &                   +(AC(NSG,NTG+1)-AC(NSG,NTG))/(4*DSG*DTG)
         IF(NTG.EQ.NTGMAX)THEN	
            Q(NBND-2*NTGMAX+1,I)=-(AB(NSG,NTG)+AC(NSG,NTG+1))
     &                            /(4*DSG*DTG)
            Q(NBND-  NTGMAX+1,I)= (AB(NSG+1,NTG)-AB(NSG,NTG))
     &                            /(4*DSG*DTG)
     &                          +  AD(NSG,NTG+1)/(DTG*DTG)
            Q(NBND         +1,I)= (AB(NSG+1,NTG)+AC(NSG,NTG+1))
     &                            /(4*DSG*DTG)
         ELSE
            Q(NBND-NTGMAX+1,I)=-(AB(NSG,NTG)+AC(NSG,NTG+1))/(4*DSG*DTG)
            Q(NBND       +1,I)= (AB(NSG+1,NTG)-AB(NSG,NTG))
     &                          /(4*DSG*DTG)
     &                        +  AD(NSG,NTG+1)/(DTG*DTG)
            Q(NBND+NTGMAX+1,I)= (AB(NSG+1,NTG)+AC(NSG,NTG+1))
     &                          /(4*DSG*DTG)
         ENDIF
      ENDDO
      ENDDO
C
      DO I=1,NTGMAX
      DO J=1,NTGMAX/2
         Q(NBND+J-I,I)=Q(NBND+J-I,I)
     &                +Q(NBND+J-I-NTGMAX/2,I)
         Q(NBND+J-I-NTGMAX/2,I)=0.D0
         Q(NBND+J-I+NTGMAX/2,I)=Q(NBND+J-I+NTGMAX/2,I)
     &                         +Q(NBND+J-I-NTGMAX,I)
         Q(NBND+J-I-NTGMAX,I)=0.D0
      ENDDO
      ENDDO
      DO I=1,NTGMAX
      DO J=1,NTGMAX
         Q(NBND+J-I,I+(NSGMAX-1)*NTGMAX)
     &                   =Q(NBND+J-I,I+(NSGMAX-1)*NTGMAX)
     &                   -Q(NBND+J-I+NTGMAX,I+(NSGMAX-1)*NTGMAX)
         Q(NBND+J-I+NTGMAX,I+(NSGMAX-1)*NTGMAX)=0.D0
      ENDDO
      ENDDO
      RETURN
      END
C
C   ************************************************ 
C   **                  RHS vector                **
C   ************************************************
C
      SUBROUTINE EQRHSV(IERR)
C
      INCLUDE 'eqcomc.h'
C
      DIMENSION PSIRG(NRGM,NZGM),PSIZG(NRGM,NZGM),PSIRZG(NRGM,NZGM)
      EXTERNAL EQPSID
C
      IERR=0
      CALL EQTORZ
      CALL SPL2D(RG,ZG,PSIRZ,PSIRG,PSIZG,PSIRZG,URZ,
     &           NRGM,NRGMAX,NZGMAX,0,0,IERR)
C
      DELT=1.D-8
      EPS=1.D-8
      ILMAX=20
      LIST=0
      RINIT=RAXIS
      ZINIT=ZAXIS
      CALL NEWTN(EQPSID,RINIT,ZINIT,RAXIS,ZAXIS,
     &            DELT,EPS,ILMAX,LIST,IER)
      IF(IER.NE.0) THEN
         WRITE(6,'(A,I5)') 'XX EQRHSV: NEWTN ERROR: IER=',IER
         IERR=101
         RETURN
      ENDIF
      IF(RAXIS.LE.RR+RB) THEN
         PSI0=PSIG(RAXIS,ZAXIS)
      ELSE
         WRITE(6,'(A)') 'XX EQRHSV: AXIS OUT OF PLASMA:'
         IERR=102
         RETURN
      ENDIF
      PSIMIN=PSI(1,1)
      PSIMAX=PSI(1,1)
      DO NSG=1,NSGMAX
      DO NTG=1,NTGMAX
         IF(PSI(NTG,NSG).LT.PSIMIN) PSIMIN=PSI(NTG,NSG)
         IF(PSI(NTG,NSG).GT.PSIMAX) PSIMAX=PSI(NTG,NSG)
      ENDDO
      ENDDO
      PSIMIN=PSIMIN/PSI0
      PSIMAX=PSIMAX/PSI0
      IF(MAX(ABS(PSIMIN),ABS(PSIMAX)).GT.3.D0*ABS(PSI0)) THEN
         WRITE(6,'(A)') 'XX EQRHSV: PSI OUT OF RANGE:'
         WRITE(6,'(A,1P3E12.4)') 
     &        '  PSIMIN,PSIMAX,PSI0=',PSIMIN,PSIMAX,PSI0
         IERR=103
         RETURN
      ENDIF
         
C     ----- positive current density, jp.gt.0-----
      RRC=RR-RA
C     ----- quasi-symmetric current density, jp:anti-symmetric -----
C      RRC=RR
C
      FJP=0.D0
      FJT=0.D0
      DO NSG=1,NSGMAX
      DO NTG=1,NTGMAX
         PSIN=PSI(NTG,NSG)/PSI0
         PP(NTG,NSG)=PPSI(PSIN)
         RMM(NTG,NSG)=RR+SIGM(NSG)*RHOM(NTG)*COS(THGM(NTG))
         HJP1(NTG,NSG)=RMM(NTG,NSG)*DPPSI(PSIN)
         HJT1(NTG,NSG)=HJPSI(PSIN)
         HJP2(NTG,NSG)=(1.D0-RRC**2/RMM(NTG,NSG)**2)*HJP1(NTG,NSG)
         HJT2(NTG,NSG)=(RRC/RMM(NTG,NSG))*HJT1(NTG,NSG)
         DVOL=SIGM(NSG)*RHOM(NTG)*RHOM(NTG)*DSG*DTG
         FJP=FJP+HJP2(NTG,NSG)*DVOL
         FJT=FJT+HJT2(NTG,NSG)*DVOL
C         IF(FJT.EQ.0.D0) WRITE(6,'(A,2I5,1P4E12.4)') 
C     &        'NTG,NSG,PSIN,PSI0,HJP1,HJT1=',
C     &         NTG,NSG,PSIN,PSI0,HJP1(NTG,NSG),HJT1(NTG,NSG)

      ENDDO
      ENDDO
C
C      WRITE(6,'(A,1P3E12.4)') 'RIP,FJP,FJT=',RIP,FJP,FJT
      TJ=(-RIP*1.D6-FJP)/FJT
      DO NSG=1,NSGMAX
      DO NTG=1,NTGMAX
         HJT(NTG,NSG)=HJP2(NTG,NSG)+TJ*HJT2(NTG,NSG)
         PSIN=PSI(NTG,NSG)/PSI0
         TT(NTG,NSG)=SQRT(BB**2*RR**2
     &                  +2.D0*RMU0*RRC
     &                  *(TJ*HJPSID(PSIN)-RRC*PPSI(PSIN)))
      ENDDO
      ENDDO
      RETURN
      END
C
C   ************************************************
C   **              Matrix Solver                 **
C   ************************************************
C
      SUBROUTINE EQSOLV
C
      INCLUDE 'eqcomc.h'
C
      DIMENSION FJT(MLM),PSIOLD(NTGM,NSGM)
C
      DO NSG=1,NSGMAX
      DO NTG=1,NTGMAX
         PSIOLD(NTG,NSG)=PSI(NTG,NSG)
      ENDDO
      ENDDO
C
      DO NSG=1,NSGMAX
         I=(NSG-1)*NTGMAX
      DO NTG=1,NTGMAX
         FJT(I+NTG)=-RMU0*HJT(NTG,NSG)*SIGM(NSG)*RHOM(NTG)*RHOM(NTG)
      ENDDO
      ENDDO
C
      CALL BANDRD(Q,FJT,NTGMAX*NSGMAX,4*NTGMAX-1,MWM,IERR)
         IF(IERR.NE.0) THEN
            WRITE(6,*) 'XX EQSOLV: BANDRD ERROR: IERR = ',IERR
         ENDIF
C
      DO NSG=1,NSGMAX
         I=(NSG-1)*NTGMAX
      DO NTG=1,NTGMAX
         PSI(NTG,NSG)=FJT(I+NTG)
      ENDDO
      ENDDO
      DO NSG=1,NSGMAX
      DO NTG=1,NTGMAX
         DELPSI(NTG,NSG)=PSI(NTG,NSG)-PSIOLD(NTG,NSG)
      ENDDO
      ENDDO
      RETURN
      END
C
C   ************************************************
C   **          sigma,theta to R,Z                **
C   ************************************************
C
      SUBROUTINE EQTORZ
C
      INCLUDE 'eqcomc.h'
      EXTERNAL EQFBND
C
      RMIN= RR-RB
      RMAX= RR+RB
      ZMIN=-RKAP*RB
      ZMAX= RKAP*RB
C
      SIG1=1.D0
      SIG2=0.95D0
C
      DTRG=(RMAX-RMIN)/(NRGMAX-1)
      DTZG=(ZMAX-ZMIN)/(NZGMAX-1)
      EPSZ=1.D-8
C
      DO NRG=1,NRGMAX
         RG(NRG)=RMIN+DTRG*(NRG-1)
      ENDDO
      DO NZG=1,NZGMAX
         ZG(NZG)=ZMIN+DTZG*(NZG-1)
      ENDDO
C
      CALL EQSETF
C
      DO NRG=1,NRGMAX
      DO NZG=1,NZGMAX
         THL=ATAN2(ZG(NZG),RG(NRG)-RR)
         IF(THL.LT.0.D0) THL=THL+2.D0*PI
         ZBRF=TAN(THL)
         THDASH=ZBRENT(EQFBND,THL-1.0D0,THL+1.0D0,EPSZ)
         RHOL=RA*SQRT(COS(THDASH+RDLT*SIN(THDASH))**2
     &               +RKAP**2*SIN(THDASH)**2)
         SIGL=SQRT((RG(NRG)-RR)**2+ZG(NZG)**2)/RHOL
         IF(SIGL.LT.1.D0) THEN
            PSIRZ(NRG,NZG)=PSIF(SIGL,THL)
         ELSE
            PSI1=PSIF(SIG1,THL)
            PSI2=PSIF(SIG2,THL)
            PSIRZ(NRG,NZG)=PSI2+(PSI1-PSI2)*(SIGL-SIG2)/(SIG1-SIG2)
         ENDIF
      ENDDO
      ENDDO
      CALL EQCALP
      RETURN
      END
C
C   ************************************************
C   **          CALCULATE pp,tt                   **
C   ************************************************
C
      SUBROUTINE EQCALP
C
      INCLUDE 'eqcomc.h'
C
      DPS=PSI0/(NPSMAX-1)
      DO NPS=1,NPSMAX
         PSIPS(NPS)=DPS*(NPS-1)
         PSIN=PSIPS(NPS)/PSI0
         PPPS(NPS)=PPSI(PSIN)
         TTPS(NPS)=SQRT(BB**2*RR**2
     &                  +2.D0*RMU0*RRC
     &                   *(TJ*HJPSID(PSIN)-RRC*PPSI(PSIN)))
      ENDDO
      RETURN
      END
C
C   ************************************************
C   **                                            **
C   ************************************************
C
      SUBROUTINE EQSETF
C
      INCLUDE 'eqcomc.h'
C
      DIMENSION PSISX(NTGPM,NSGPM),PSITX(NTGPM,NSGPM)
      DIMENSION PSISTX(NTGPM,NSGPM)
      NSGPMAX=NSGMAX+2
      NTGPMAX=NTGMAX+2
C
      SIGMX(1)=0.D0
      DO NSG=1,NSGMAX
         SIGMX(NSG+1)=SIGM(NSG)
      ENDDo
      SIGMX(NSGMAX+2)=1.D0
C
      THGMX(1)=0.D0
      DO NTG=1,NTGMAX
         THGMX(NTG+1)=THGM(NTG)
      ENDDO
      THGMX(NTGMAX+2)=2.D0*PI
C
      SUM=0.D0
      DO NTG=1,NTGMAX
         SUM=SUM+(9*PSI(NTG,1)-PSI(NTG,2))/8.D0
      ENDDO
      PSI1=SUM/NTGMAX
      DO NTG=1,NTGMAX+2
         PSIX(NTG,1)=PSI1
      ENDDO
C
      DO NSG=1,NSGMAX
      DO NTG=1,NTGMAX
         PSIX(NTG+1,NSG+1)=PSI(NTG,NSG)
      ENDDO
      ENDDO
C
      DO NSG=1,NSGMAX
         PSIX(       1,NSG+1)=(9*PSI(     1,NSG)-PSI(       2,NSG))
     &                        /16.D0
     &                       +(9*PSI(NTGMAX,NSG)-PSI(NTGMAX-1,NSG))
     &                        /16.D0
         PSIX(NTGMAX+2,NSG+1)=PSIX(1,NSG+1)
      ENDDO
C
      DO NTG=1,NTGMAX+2
         PSIX(NTG,NSGMAX+2)=0.D0
      ENDDO
C
      DO NSGP=1,NSGMAX+2
         PSITX(       1,NSGP)=(PSIX(2,NSGP)-PSIX(NTGMAX+1,NSGP))
     &                       /(2.D0*PI+THGMX(2)-THGMX(NTGMAX+1))
         PSITX(NTGMAX+2,NSGP)=PSITX(       1,NSGP)
      ENDDO
      PSISTX(       1,       1)=(PSITX(1,2)-PSITX(1,1))
     &                         /(SIGMX(  2)-SIGMX(  1))
      PSISTX(NTGMAX+2,       1)=PSISTX(1,1)
      PSISTX(       1,NSGMAX+2)=(PSITX(1,NSGMAX+2)-PSITX(1,NSGMAX+1))
     &                         /(SIGMX(  NSGMAX+2)-SIGMX(  NSGMAX+1))
      PSISTX(NTGMAX+2,NSGMAX+2)=PSISTX(1,NSGMAX+2)
C
      CALL SPL2D(THGMX,SIGMX,PSIX,PSITX,PSISX,PSISTX,U,
     &           NTGPM,NTGPMAX,NSGPMAX,3,0,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX EQSETF: SPL2D ERROR : IERR=',IERR
      RETURN
      END
C
C   ************************************************
C   **                                            **
C   ************************************************
C
      FUNCTION PSIF(RSIG,RTHG)
C
      INCLUDE 'eqcomc.h'
C
      CALL SPL2DF(RTHG,RSIG,PSIL,THGMX,SIGMX,U,
     &            NTGPM,NTGPMAX,NSGPMAX,IERR)
      IF(IERR.NE.0) WRITE(6,*) 'XX PSIF: SPL2DF ERROR : IERR=',IERR
      PSIF=PSIL
      RETURN
      END
C
C   ************************************************
C
      SUBROUTINE EQSAVE
C
      INCLUDE 'eqcomc.h'
C
      CHARACTER*32 KNAM
C      CHARACTER*1 KID
      LOGICAL LEX
C
    1 WRITE(6,*) '#EQ> INPUT : SAVE FILE NAME : ',KNAMEQ
      READ(5,'(A32)',ERR=1,END=900) KNAM
      IF(KNAM(1:2).NE.'/ ') KNAMEQ=KNAM
      IF(KNAM(1:2).EQ.'  ') GOTO 900
C
      INQUIRE(FILE=KNAMEQ,EXIST=LEX,ERR=1)
      IF(LEX) THEN
         OPEN(21,FILE=KNAMEQ,IOSTAT=IST,STATUS='OLD',ERR=10,
     &        FORM='UNFORMATTED')
         WRITE(6,*) '# OLD FILE (',KNAMEQ,') IS ASSIGNED FOR OUTPUT.'
         GOTO 30
   10    WRITE(6,*) 'XX OLD FILE OPEN ERROR : IOSTAT = ',IST
         GOTO 1
      ELSE
         OPEN(21,FILE=KNAMEQ,IOSTAT=IST,STATUS='NEW',ERR=20,
     &        FORM='UNFORMATTED')
         WRITE(6,*) '# NEW FILE (',KNAMEQ,') IS CREATED FOR OUTPUT.'
         GOTO 30
   20    WRITE(6,*) 'XX NEW FILE OPEN ERROR : IOSTAT = ',IST
         GOTO 1
      ENDIF
C
   30 WRITE(21) RR,BB,RIP
      WRITE(21) NRGMAX,NZGMAX
      WRITE(21) (RG(NRG),NRG=1,NRGMAX)
      WRITE(21) (ZG(NZG),NZG=1,NZGMAX)
      WRITE(21) ((PSIRZ(NRG,NZG),NRG=1,NRGMAX),NZG=1,NZGMAX)
      WRITE(21) NPSMAX
      WRITE(21) (PSIPS(NPS),NPS=1,NPSMAX)
      WRITE(21) (PPPS(NPS),NPS=1,NPSMAX)
      WRITE(21) (TTPS(NPS),NPS=1,NPSMAX)
C
      WRITE(21) NSGMAX,NTGMAX
      WRITE(21) RA,RKAP,RDLT,RB
      WRITE(21) PJ0,PJ1,PJ2,PROFJ0,PROFJ1,PROFJ2
      WRITE(21) PP0,PP1,PP2,PROFP0,PROFP1,PROFP2
      WRITE(21) PROFR0,PROFR1,PROFR2
      CLOSE(21)
C
      WRITE(6,*) '# DATA WAS SUCCESSFULLY SAVED TO THE FILE.'
C
  900 RETURN
      END
