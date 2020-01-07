C     $Id$
C
C     ***** SAVE TASK/EQ DATA *****
C
      SUBROUTINE EQSAVE
C
      USE libfio
      INCLUDE '../eq/eqcomc.inc'
C
      CALL FWOPEN(21,KNAMEQ,0,MODEFW,'EQ',IERR)
      IF(IERR.NE.0) RETURN
C
      REWIND(21)
      WRITE(21) RR,BB,RIP
      WRITE(21) NRGMAX,NZGMAX
      WRITE(21) (RG(NRG),NRG=1,NRGMAX)
      WRITE(21) (ZG(NZG),NZG=1,NZGMAX)
      WRITE(21) ((PSIRZ(NRG,NZG),NRG=1,NRGMAX),NZG=1,NZGMAX)
      WRITE(21) NPSMAX
      WRITE(21) (PSIPS(NPS),NPS=1,NPSMAX)
      WRITE(21) (PPPS(NPS),NPS=1,NPSMAX)
      WRITE(21) (TTPS(NPS),NPS=1,NPSMAX)
      WRITE(21) (TEPS(NPS),NPS=1,NPSMAX)
      WRITE(21) (OMPS(NPS),NPS=1,NPSMAX)
C
      WRITE(21) NSGMAX,NTGMAX,NUGMAX,NRMAX,NTHMAX,NSUMAX,NRVMAX,NTVMAX
      WRITE(21) ((PSI(NSG,NTG),NSG=1,NSGMAX),NTG=1,NTGMAX)
      WRITE(21) ((DELPSI(NSG,NTG),NSG=1,NSGMAX),NTG=1,NTGMAX)
      WRITE(21) ((HJT(NSG,NTG),NSG=1,NSGMAX),NTG=1,NTGMAX)
      WRITE(21) RAXIS,ZAXIS,PSITA,PSIPA,PSI0
      WRITE(21) (PSIPNV(NRV),NRV=1,NRVMAX)
      WRITE(21) (PSIPV(NRV),NRV=1,NRVMAX)
      WRITE(21) (PSITV(NRV),NRV=1,NRVMAX)
      WRITE(21) (QPV(NRV),NRV=1,NRVMAX)
      WRITE(21) (TTV(NRV),NRV=1,NRVMAX)
      WRITE(21) RA,RKAP,RDLT,RB,FRBIN
      WRITE(21) PJ0,PJ1,PJ2,PROFJ0,PROFJ1,PROFJ2
      WRITE(21) PP0,PP1,PP2,PROFP0,PROFP1,PROFP2
      WRITE(21) PT0,PT1,PT2,PROFTP0,PROFTP1,PROFTP2
      WRITE(21) PV0,PV1,PV2,PROFV0,PROFV1,PROFV2
      WRITE(21) PROFR0,PROFR1,PROFR2
C      WRITE(21) PTS,PN0,HM
      WRITE(21) ((HJTRZ(NRG,NZG),NRG=1,NRGMAX),NZG=1,NZGMAX)
      CLOSE(21)
C
C      WRITE(6,*) 'HJTRZ=',HJTRZ(10,10)
C
      WRITE(6,*) '# DATA WAS SUCCESSFULLY SAVED TO THE FILE.'
C
      RETURN
      END
C
C     ***** LOAD EQUILIBRIUM DATA *****
C
      SUBROUTINE EQLOAD(MODELG1,KNAMEQ1,IERR)
C
      USE eq_bpsd_mod
      INCLUDE '../eq/eqcomc.inc'
      CHARACTER KNAMEQ1*80
      INTEGER ierr
C
      MODELG=MODELG1
      KNAMEQ=KNAMEQ1
      CALL EQREAD(IERR)
      RETURN
      END
C
C     ***** LOAD EQUILIBRIUM DATA *****
C
      SUBROUTINE EQREAD(IERR)
C
      USE eqread_mod
      INCLUDE '../eq/eqcomc.inc'
C
      IF(MODELG.EQ.3.OR.MODELG.EQ.9) THEN
         CALL EQRTSK(IERR)
      ELSEIF(MODELG.EQ.5) THEN
         CALL EQDSKR(IERR)
         CALL EQCALQ(IERR)
      ELSEIF(MODELG.EQ.8) THEN
         CALL EQJAEAR(IERR)
      ELSEIF(MODELG.EQ.15) THEN
         CALL EQDSK
      ELSE
         WRITE(6,*) 'XX EQLOAD: UNKNOWN MODELG: MODELG=',MODELG
      ENDIF
C
      RETURN
      END
C
C     ***** LOAD TASK/EQ DATA *****
C
      SUBROUTINE EQRTSK(IERR)
C
      USE libfio
      INCLUDE '../eq/eqcomc.inc'
      DIMENSION DERIV(NRVM)
C
      CALL FROPEN(21,KNAMEQ,0,MODEFR,'EQ',IERR)
      IF(IERR.NE.0) RETURN
C
      READ(21) RR,BB,RIP
      READ(21) NRGMAX,NZGMAX
      READ(21) (RG(NRG),NRG=1,NRGMAX)
      READ(21) (ZG(NZG),NZG=1,NZGMAX)
      READ(21) ((PSIRZ(NRG,NZG),NRG=1,NRGMAX),NZG=1,NZGMAX)
      READ(21) NPSMAX
      READ(21) (PSIPS(NPS),NPS=1,NPSMAX)
      READ(21) (PPPS(NPS),NPS=1,NPSMAX)
      READ(21) (TTPS(NPS),NPS=1,NPSMAX)
      READ(21) (TEPS(NPS),NPS=1,NPSMAX)
      READ(21) (OMPS(NPS),NPS=1,NPSMAX)
C
      READ(21) NSGMAX,NTGMAX,NUGMAX,NRMAX,NTHMAX,NSUMAX,NRVMAX,NTVMAX
      READ(21) ((PSI(NTG,NSG),NTG=1,NTGMAX),NSG=1,NSGMAX)
      READ(21) ((DELPSI(NTG,NSG),NTG=1,NTGMAX),NSG=1,NSGMAX)
      READ(21) ((HJT(NTG,NSG),NTG=1,NTGMAX),NSG=1,NSGMAX)
      READ(21) RAXIS,ZAXIS,PSITA,PSIPA,PSI0
      READ(21) (PSIPNV(NRV),NRV=1,NRVMAX)
      READ(21) (PSIPV(NRV),NRV=1,NRVMAX)
      READ(21) (PSITV(NRV),NRV=1,NRVMAX)
      READ(21) (QPV(NRV),NRV=1,NRVMAX)
      READ(21) (TTV(NRV),NRV=1,NRVMAX)
      READ(21) RA,RKAP,RDLT,RB,FRBIN
      READ(21) PJ0,PJ1,PJ2,PROFJ0,PROFJ1,PROFJ2
      READ(21) PP0,PP1,PP2,PROFP0,PROFP1,PROFP2
      READ(21) PT0,PT1,PT2,PROFT0,PROFT1,PROFT2
      READ(21) PV0,PV1,PV2,PROFV0,PROFV1,PROFV2
      READ(21) PROFR0,PROFR1,PROFR2
C      READ(21) PTS,PN0,HM
      READ(21,ERR=1000) ((HJTRZ(NRG,NZG),NRG=1,NRGMAX),NZG=1,NZGMAX)
      GOTO 1001
 1000 CONTINUE
         DO NZG=1,NZGMAX
         DO NRG=1,NRGMAX
            HJTRZ(NRG,NZG)=0.D0
         ENDDO
         ENDDO
 1001 CONTINUE
C      
      CLOSE(21)
      RIPX=RIP
C
      IF(MODELG.EQ.9) THEN
         CALL EQMESH
         CALL SPL1D(PSIPNV,PSITV,DERIV,UPSITV,NRVMAX,0,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for PSITV: IERR=',IERR
         CALL SPL1D(PSIPNV,QPV,DERIV,UQPV,NRVMAX,0,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for QPV: IERR=',IERR
         CALL SPL1D(PSIPNV,TTV,DERIV,UTTV,NRVMAX,0,IERR)
         IF(IERR.NE.0) WRITE(6,*) 'XX SPL1D for TTV: IERR=',IERR
         CALL EQDEFB
C         DO NSG=1,NSGMAX
C            WRITE(6,'(A,I5)') 'PSI NSG=',NSG
C            WRITE(6,'(1P5E12.4)') (PSI(NTG,NSG),NTG=1,NTGMAX)
C         ENDDO
C         DO NSG=1,NSGMAX
C            WRITE(6,'(A,I5)') 'HJT NSG=',NSG
C            WRITE(6,'(1P5E12.4)') (HJT(NTG,NSG),NTG=1,NTGMAX)
C         ENDDO
C         WRITE(6,'(A,I5,1P4E12.4)')
C     &        ('NV:',NV,PSIPNV(NV),PSITV(NV),QPV(NV),TTV(NV),
C     &         NV=1,NRVMAX)
      ENDIF
C     
C     WRITE(6,*) 'HJTRZ=',HJTRZ(10,10)
C
      RETURN
      END
C
C     ***** SAVE METRICS *****
C
      SUBROUTINE EQMETRIC(IERR)
C
      USE libfio
      INCLUDE '../eq/eqcomq.inc'
C
      character KNAMET*80
      data KNAMET /'eq_metric.dat'/ 
C
      nmetric=21
      CALL FWOPEN(nmetric,KNAMET,1,MODEFW,'EQ',IERR)
      IF(IERR.NE.0) RETURN
C
      REWIND(nmetric)
      WRITE(nmetric,'(A,2X,A,6X,A,8X,A,3X,A)') '#','rho_tor','dV/drho',
     &     '<1/R^2>',"<|grad rho|^2/R^2>"
      DO NR = 1, NRPMAX
         CALL SPL1DF(FNPSIP(RHOT(NR)),DAT1,PSIP,UDVDRHO ,NRMAX,IERR)
         WRITE(nmetric,'(1X,F10.7,1P3E15.7)') RHOT(NR),DAT1,
     &        fnavir2(rhot(nr)),fnavgrr2(nr)
      ENDDO
      WRITE(nmetric,'(80X)')
      WRITE(nmetric,'(A,2X,A,3X,A,3X,A,5X,A,9X,A)') '#','rho_tor',
     &     "<|grad rho|>","<|grad rho|^2>","<B^2>","<1/B^2>"
      WRITE(nmetric,'(1X,0PF10.7,1P4E15.7)') (RHOT(NR),
     &     fnavgr(rhot(nr)),fnavgr2(rhot(nr)),
     &     fnavbb2(rhot(nr)),fnavib2(rhot(nr)),NR=1,NRPMAX)
      CLOSE(nmetric)
C
      WRITE(6,*) '# METRIC DATA WAS SUCCESSFULLY SAVED TO "',
     &           KNAMET(1:13),'".'
C
      RETURN
      END
C
C     ***** READ RIPPLE CONTOUR DATA FROM OFMC ****
C
      subroutine read_rppl(ierr)
C
      USE libfio
      INCLUDE '../eq/eqcomq.inc'
C
      character kfile*20, kline*130
C
      kfile='ripple.profile'
      nrppl=21
      CALL FROPEN(nrppl,kfile,1,MODEFR,'EQ',IERR)
      IF(IERR.NE.0) RETURN
c
      rewind(nrppl)
c
c     *** R-coordinates ***
      idx = 0
      do
         if(idx == 0) then
            read(nrppl,'(A130)',iostat=ist) kline
            if(index(kline,"R-coordinate") /= 0) then ! detect the start position of the data chunk
               idx = 1
            end if
            cycle
         end if
c
         read(nrppl,'(1x,10e13.5)') (Rrp(i),i=1,NRrpM)
         exit
      end do
c
c     *** Z-coordinates ***
      idx = 0
      do
         if(idx == 0) then
            read(nrppl,'(A130)',iostat=ist) kline
            if(index(kline,"Z-coordinate") /= 0) then ! detect the start position of the data chunk
               idx = 1
            end if
            cycle
         end if
c
         read(nrppl,'(1x,10e13.5)') (Zrp(i),i=1,NZrpM)
         exit
      end do
c
c     *** Ripple contour ***
      idx = 0
      j   = 0
      do
         if(idx == 0) then
            read(nrppl,'(A130)',iostat=ist) kline
            if(index(kline,"at") /= 0) then ! detect the start position of the data chunk
               idx = 1
               j = j + 1
            end if
            cycle
         end if
c
         read(nrppl,'(1x,10e13.5)') (RpplRZ(i,j),i=1,NRrpM)
         idx = 0
         if(j == NZrpM) then
            exit
         else
            cycle
         end if
      end do
c
      close(nrppl)
c
      return
      end

      SUBROUTINE draw_cross(x,y,len)
      IMPLICIT NONE
      REAL(8),INTENT(IN):: x,y,len
      INTERFACE
         FUNCTION GUCLIP(x)
            REAL(8),INTENT(IN):: x
            REAL(4):: GUCLIP
         END FUNCTION GUCLIP
      END INTERFACE

      CALL MOVE2D(GUCLIP(x-0.5D0*len),GUCLIP(y))
      CALL DRAW2D(GUCLIP(x+0.5D0*len),GUCLIP(y))
      CALL MOVE2D(GUCLIP(x),GUCLIP(y-0.5D0*len))
      CALL DRAW2D(GUCLIP(x),GUCLIP(y+0.5D0*len))
      RETURN
      END
