C     $Id$
C
      SUBROUTINE PLDATA_CLEAR
C
      INCLUDE 'plcomd.inc'
C
      DO NS=1,NSMAXPL
         DO NR=1,NRMAXPL
            RNPL(NR,NS)=0.D0
            RTPL(NR,NS)=3.D-5
            RUPL(NR,NS)=0.D0
         ENDDO
      ENDDO
      DO NWR=1,NWRM
         KIDWR(NWR)=' '
         DO NS=1,NSMAXPL
            DO NR=1,NRMAXPL
               PWRPL(NR,NS,NWR) =0.D0
            ENDDO
         ENDDO
         DO NR=1,NRMAXPL
            AJWRPL(NR,NWR)=0.D0
         ENDDO
      ENDDO
      DO NWM=1,NWMM
         KIDWM(NWM)=' '
         DO NS=1,NSMAXPL
            DO NR=1,NRMAXPL
               PWMPL(NR,NS,NWM) =0.D0
            ENDDO
         ENDDO
         DO NR=1,NRMAXPL
            AJWMPL(NR,NWM)=0.D0
         ENDDO
      ENDDO
      NWRMAX=0
      NWMMAX=0
      NSPL=0
      RETURN
      END
C
      SUBROUTINE PLDATA_SETN(NRMAXPL1,NSMAXPL1)
C
      INCLUDE 'plcomd.inc'
C
      NRMAXPL=NRMAXPL1
      NSMAXPL=NSMAXPL1
      RETURN
      END
C
      SUBROUTINE PLDATA_GETN(NRMAXPL1,NSMAXPL1)
C
      INCLUDE 'plcomd.inc'
C
      NRMAXPL1=NRMAXPL
      NSMAXPL1=NSMAXPL
      RETURN
      END
C
      SUBROUTINE PLDATA_GETNW(NWRMAX1,NWMMAX1)
C
      INCLUDE 'plcomd.inc'
C
      NWRMAX1=NWRMAX
      NWMMAX1=NWMMAX
      RETURN
      END
C
      SUBROUTINE PLDATA_SETR(RGPL1,RMPL1)
C
      INCLUDE 'plcomd.inc'
      DIMENSION RGPL1(NRM),RMPL1(NRM)
C
      DO NR=1,NRMAXPL
         RGPL(NR)=RGPL1(NR)
         RMPL(NR)=RMPL1(NR)
      ENDDO
      RETURN
      END
C
      SUBROUTINE PLDATA_GETR(RGPL1,RMPL1)
C
      INCLUDE 'plcomd.inc'
      DIMENSION RGPL1(NRM),RMPL1(NRM)
C
      DO NR=1,NRMAXPL
         RGPL1(NR)=RGPL(NR)
         RMPL1(NR)=RMPL(NR)
      ENDDO
      RETURN
      END
C
      SUBROUTINE PLDATA_SETP(NRM1,RNPL1,RTPL1,RUPL1)
C
      INCLUDE 'plcomd.inc'
      DIMENSION RNPL1(NRM1,NSM),RTPL1(NRM1,NSM),RUPL1(NRM1,NSM)
C
C      WRITE(6,*) NRM,NRM1
C
      DO NS=1,NSMAXPL
      DO NR=1,NRMAXPL
         RNPL(NR,NS)=RNPL1(NR,NS)
         RTPL(NR,NS)=RTPL1(NR,NS)
         RUPL(NR,NS)=RUPL1(NR,NS)
      ENDDO
      ENDDO
C      WRITE(6,'(1P4E12.4)') RNPL(1,1),RTPL(1,1),RNPL(1,2),RTPL(1,2)
      NSPL=0
      RETURN
      END
C
      SUBROUTINE PLDATA_GETP(RNPL1,RTPL1,RUPL1)
C
      INCLUDE 'plcomd.inc'
      DIMENSION RNPL1(NRM,NSM),RTPL1(NRM,NSM),RUPL1(NRM,NSM)
C
      DO NS=1,NSMAXPL
      DO NR=1,NRMAXPL
         RNPL1(NR,NS)=RNPL(NR,NS)
         RTPL1(NR,NS)=RTPL(NR,NS)
         RUPL1(NR,NS)=RUPL(NR,NS)
      ENDDO
      ENDDO
      RETURN
      END
C
      SUBROUTINE PLDATA_GETPL(RHOL,RNPL1,RTPL1,RUPL1)
C
      INCLUDE 'plcomd.inc'
      DIMENSION RNPL1(NSM),RTPL1(NSM),RUPL1(NSM)
C
      IF(NSPL.EQ.0) CALL PLDATA_SPL
C
      RHOLL=RHOL
      IF(RHOLL.GE.1.D0) RHOLL=1.D0
C
      DO NS=1,NSMAXPL
         CALL SPL1DF(RHOLL,RNPL1(NS),URPL,URNPL(1,1,NS),
     &               NRMAXPL+2,IERR)
         IF(IERR.NE.0) 
     &        WRITE(6,*) 'XX PLDATA_GETPL/RN: SPL1DF ERROR: IERR=',IERR
         CALL SPL1DF(RHOLL,RTPL1(NS),URPL,URTPL(1,1,NS),
     &               NRMAXPL+2,IERR)
         IF(IERR.NE.0) 
     &        WRITE(6,*) 'XX PLDATA_GETPL/RT: SPL1DF ERROR: IERR=',IERR
         CALL SPL1DF(RHOLL,RUPL1(NS),URPL,URUPL(1,1,NS),
     &               NRMAXPL+2,IERR)
         IF(IERR.NE.0) 
     &        WRITE(6,*) 'XX PLDATA_GETP/RUL: SPL1DF ERROR: IERR=',IERR
      ENDDO
      RETURN
      END
C
      SUBROUTINE PLDATA_SPL
C
      INCLUDE 'plcomd.inc'
      DIMENSION F(NRMPP),DF(NRMPP)
C
      DO NR=1,NRMAXPL
         URPL(NR+1)=RMPL(NR)
      ENDDO
      URPL(1)=0.D0
      URPL(NRMAXPL+2)=1.D0
C
      DO NS=1,NSMAXPL
         DO NR=1,NRMAXPL
            F(NR+1)=RNPL(NR,NS)
         ENDDO
         F(1)=(9.D0*F(2)-F(3))/8.D0
         F(NRMAXPL+2)=(3.D0*F(NRMAXPL+1)-F(NRMAXPL))/2.D0
         DF(1)=0.D0         
         CALL SPL1D(URPL,F,DF,URNPL(1,1,NS),NRMAXPL+2,1,IERR)
         IF(IERR.NE.0) 
     &        WRITE(6,*) 'XX PLDATA_SPL/RN: SPL1D ERROR: IERR=',IERR
C
         DO NR=1,NRMAXPL
            F(NR+1)=RTPL(NR,NS)
         ENDDO
         F(1)=(9.D0*F(2)-F(3))/8.D0
         F(NRMAXPL+2)=(3.D0*F(NRMAXPL+1)-F(NRMAXPL))/2.D0
         DF(1)=0.D0         
         CALL SPL1D(URPL,F,DF,URTPL(1,1,NS),NRMAXPL+2,1,IERR)
         IF(IERR.NE.0) 
     &        WRITE(6,*) 'XX PLDATA_SPL/RT: SPL1D ERROR: IERR=',IERR
C
         DO NR=1,NRMAXPL
            F(NR+1)=RUPL(NR,NS)
         ENDDO
         F(1)=(9.D0*F(2)-F(3))/8.D0
         F(NRMAXPL+2)=(3.D0*F(NRMAXPL+1)-F(NRMAXPL))/2.D0
         DF(1)=0.D0         
         CALL SPL1D(URPL,F,DF,URUPL(1,1,NS),NRMAXPL+2,1,IERR)
         IF(IERR.NE.0) 
     &        WRITE(6,*) 'XX PLDATA_SPL/RU: SPL1D ERROR: IERR=',IERR
C
      ENDDO
      NSPL=1
      RETURN
      END
C
      SUBROUTINE PLDATA_SETWR(NWR,KID,PWRPL1,AJWRPL1)
C
      INCLUDE 'plcomd.inc'
      DIMENSION PWRPL1(NRM,NSM),AJWRPL1(NRM)
      CHARACTER KID*80
C
      IF(NWR.GT.NWRM) THEN
         WRITE(6,*) 'XX PLDATA_SETWR: NWR.GT.NWRM'
         STOP
      ENDIF
      IF(NWR.GT.NWRMAX) NWRMAX=NWR
      KIDWR(NWR)=KID
C
      DO NS=1,NSMAXPL
      DO NR=1,NRMAXPL
         PWRPL(NR,NS,NWR)=PWRPL1(NR,NS)
      ENDDO
      ENDDO
      DO NR=1,NRMAXPL
         AJWRPL(NR,NWR)=AJWRPL1(NR)
      ENDDO
      RETURN
      END
C
      SUBROUTINE PLDATA_GETWR(NWR,KID,PWRPL1,AJWRPL1)
C
      INCLUDE 'plcomd.inc'
      DIMENSION PWRPL1(NRM,NSM),AJWRPL1(NRM)
      CHARACTER KID*80
C
      IF(NWR.GT.NWRM) THEN
         WRITE(6,*) 'XX PLDATA_GETWR: NWR.GT.NWRM'
         RETURN
      ENDIF
C
      KID=KIDWR(NWR)
C
      DO NS=1,NSMAXPL
      DO NR=1,NRMAXPL
         PWRPL1(NR,NS)=PWRPL(NR,NS,NWR)
      ENDDO
      ENDDO
      DO NR=1,NRMAXPL
         AJWRPL1(NR)=AJWRPL(NR,NWR)
      ENDDO
      RETURN
      END
C
      SUBROUTINE PLDATA_SETWM(NWM,KID,PWMPL1,AJWMPL1)
C
      INCLUDE 'plcomd.inc'
      DIMENSION PWMPL1(NRM,NSM),AJWMPL1(NRM)
      CHARACTER KID*80
C
      IF(NWM.GT.NWMM) THEN
         WRITE(6,*) 'XX PLDATA_SETWM: NWM.GT.NWMM'
         STOP
      ENDIF
      IF(NWM.GT.NWMMAX) NWMMAX=NWM
C
      KIDWM(NWM)=KID
C
      DO NS=1,NSMAXPL
      DO NR=1,NRMAXPL
         PWMPL(NR,NS,NWM)=PWMPL1(NR,NS)
      ENDDO
      ENDDO
      DO NR=1,NRMAXPL
         AJWMPL(NR,NWM)=AJWMPL1(NR)
      ENDDO
      RETURN
      END
C
      SUBROUTINE PLDATA_GETWM(NWM,KID,PWMPL1,AJWMPL1)
C
      INCLUDE 'plcomd.inc'
      DIMENSION PWMPL1(NRM,NSM),AJWMPL1(NRM)
      CHARACTER KID*80
C
      IF(NWM.GT.NWMM) THEN
         WRITE(6,*) 'XX PLDATA_GETWW: NWM.GT.NWMM'
         RETURN
      ENDIF
C
      KID=KIDWM(NWM)
C
      DO NS=1,NSMAXPL
      DO NR=1,NRMAXPL
         PWMPL1(NR,NS)=PWMPL(NR,NS,NWM)
      ENDDO
      ENDDO
      DO NR=1,NRMAXPL
         AJWMPL1(NR)=AJWMPL(NR,NWM)
      ENDDO
      RETURN
      END
