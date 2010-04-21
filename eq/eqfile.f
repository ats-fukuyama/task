C     $Id$
C
C     ***** SAVE TASK/EQ DATA *****
C
      SUBROUTINE EQSAVE
C
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
      WRITE(21) RA,RKAP,RDLT,RB
      WRITE(21) PJ0,PJ1,PJ2,PROFJ0,PROFJ1,PROFJ2
      WRITE(21) PP0,PP1,PP2,PROFP0,PROFP1,PROFP2
      WRITE(21) PT0,PT1,PT2,PROFT0,PROFT1,PROFT2
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
      INCLUDE '../eq/eqcomc.inc'
      CHARACTER KNAMEQ1*80
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
      INCLUDE '../eq/eqcomc.inc'
C
      IF(MODELG.EQ.3.OR.MODELG.EQ.9) THEN
         CALL EQRTSK(IERR)
      ELSEIF(MODELG.EQ.5) THEN
         CALL EQDSKR(IERR)
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
      READ(21) RA,RKAP,RDLT,RB
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
C     ***** READ EQDSK FORMAT FILE *****
C
      SUBROUTINE EQDSKR(IERR)
C
      INCLUDE '../eq/eqcomc.inc'
C
      character case(6)*10
      dimension rlim(NSUM),zlim(NSUM),pressw(NPSM),pwprim(NPSM),
     &          dmion(NSUM),rhovn(NSUM),ajtor(NPSM)
C
      neqdsk=21
      CALL FROPEN(neqdsk,KNAMEQ,1,MODEFR,'EQ',IERR)
      IF(IERR.NE.0) RETURN
c
      REWIND(neqdsk)
      read (neqdsk,2000) (case(i),i=1,6),idum,NRGMAX,NZGMAX
      NPSMAX=NRGMAX
      read (neqdsk,2020) rdim,zdim,Rctr,rleft,zmid
CChonda      RA=0.5D0*rdim
CChonda      RKAP=zdim/rdim
CChonda      RDLT=0.D0
CChonda      RB=1.1D0*RA
      DR=rdim/(NRGMAX-1)
      DZ=zdim/(NZGMAX-1)
      DO NRG=1,NRGMAX
         RG(NRG)=rleft+DR*(NRG-1)
      ENDDO
      DO NZG=1,NZGMAX
         ZG(NZG)=zmid-0.5D0*zdim+DZ*(NZG-1)
      ENDDO
      read (neqdsk,2020) RAXIS,ZAXIS,PSI0,PSIA,Bctr
C
      PSI0=2.D0*PI*PSI0
      PSIA=2.D0*PI*PSIA
C
      DPS=(PSIA-PSI0)/(NPSMAX-1)
      PSI0=PSI0-PSIA
      PSIPA=-PSI0
      DO NPS=1,NPSMAX
         PSIPS(NPS)=DPS*(NPS-1)
      ENDDO
      read (neqdsk,2020) RIP,simag,xdum,rmaxis,xdum
      RIP=RIP/1.D6
      read (neqdsk,2020) zmaxis,xdum,sibry,xdum,xdum
      read (neqdsk,2020) (TTPS(i),i=1,NPSMAX)
      read (neqdsk,2020) (PPPS(i),i=1,NPSMAX)
      read (neqdsk,2020) (TTDTTPS(i),i=1,NPSMAX)
      read (neqdsk,2020) (DPPPS(i),i=1,NPSMAX)
      read (neqdsk,2020) ((PSIRZ(i,j),i=1,NRGMAX),j=1,NZGMAX)
      read (neqdsk,2020) (QQPS(i),i=1,NPSMAX)
      read (neqdsk,2022) NSUMAX,limitr
      read (neqdsk,2020) (RSU(i),ZSU(i),i=1,NSUMAX)
      read (neqdsk,2020) (rlim(i),zlim(i),i=1,limitr)
C
      RSUMAX = Rctr
      RSUMIN = Rctr
      ZSUMAX = 0.d0
      ZSUMIN = 0.d0
      R_ZSUMAX = 0.d0
      R_ZSUMIN = 0.d0
      do i = 1, nsumax
         if(RSU(i) .GT. RSUMAX) RSUMAX = RSU(i)
         if(RSU(i) .LT. RSUMIN) RSUMIN = RSU(i)
         if(ZSU(i) .GT. ZSUMAX) then
            ZSUMAX   = ZSU(i)
            R_ZSUMAX = RSU(i)
         end if
         if(ZSU(i) .LT. ZSUMIN) then
            ZSUMIN   = ZSU(i)
            R_ZSUMIN = RSU(i)
         end if
      enddo
C *** The following variable defined in Tokamaks 3rd, Sec. 14.14 ***
      RR   = 0.5d0 * (RSUMAX - RSUMIN) + RSUMIN
      RA   = RR - RSUMIN
      !==  RB: wall minor radius  ======================
      !    Multiplication factor 1.1 is tentatively set.
      RB   = 1.1d0 * RA
!      RB   = 1.2d0 * RA
      !=================================================
      RKAP = (ZSUMAX - ZSUMIN) / (RSUMAX - RSUMIN)

!  ---- corrected on 2010/01/18 for negative triangularity ----
!      RDLT = 0.5d0 * (ABS(RR-R_ZSUMIN) + ABS(RR-R_ZSUMAX)) / RA

      RDLT = 0.5d0 * ((RR-R_ZSUMIN) + (RR-R_ZSUMAX)) / RA
      BB   = Rctr * Bctr / RR
      RIPX = RIP
C
C      GOTO 1000
C      kvtor=0
C      rvtor=0
C      nmass=0
C      read (neqdsk,2024,end=1000) kvtor,rvtor,nmass
C      WRITE(6,*) kvtor,rvtor,nmass
C      if (kvtor.gt.0) then
C         read (neqdsk,2020) (pressw(i),i=1,NPSMAX)
C         read (neqdsk,2020) (pwprim(i),i=1,NPSMAX)
C      endif
C      if (nmass.gt.0) then
C         read (neqdsk,2020) (dmion(i),i=1,NPSMAX)
C      endif
C      read (neqdsk,2020,end=1000) (rhovn(i),i=1,NPSMAX)
C 1000 CONTINUE
C
      REWIND(neqdsk)
      CLOSE(neqdsk)
!      write (6,'(I5,1PE12.4)') (i,QQPS(i),i=1,NPSMAX)
C
C      WRITE(6,'(1P3E12.4)') RR,BB,RIP
C      WRITE(6,'(1P4E12.4)') RAXIS,ZAXIS,PSI0,PSIA
C      WRITE(6,'(1P4E12.4)') RG(1),RG(2),RG(NRGMAX-1),RG(NRGMAX)
C      WRITE(6,'(1P4E12.4)') ZG(1),ZG(2),ZG(NZGMAX-1),ZG(NZGMAX)
C      WRITE(6,'(1P4E12.4)') PSIPS(1),PSIPS(2),
C     &                      PSIPS(NPSMAX-1),PSIPS(NPSMAX)
C
      DO NZG=1,NZGMAX
         DO NRG=1,NRGMAX
            PSIRZ(NRG,NZG)=2.D0*PI*PSIRZ(NRG,NZG)-PSIA
         ENDDO
      ENDDO
      DO i=1,NPSMAX
         TTPS(i)   =2.D0*PI*TTPS(i)
         TTDTTPS(i)=4.D0*PI**2*TTDTTPS(i)
         DTTPS(i)  =TTDTTPS(i)/TTPS(i)
Chonda         write(6,*) PSIPS(i),QQPS(i)
      ENDDO
C
      DO NZG=1,NZGMAX
      DO NRG=1,NRGMAX
         HJTRZ(NRG,NZG)=0.D0
      ENDDO
      ENDDO
C
C     ** Simplified check for Toroidal current and parallel current **
C
c$$$      DO i=1,NPSMAX
c$$$         write(6,*) PSIPS(i),
c$$$     &              -RR*DPPPS(i)-TTDTTPS(i)/(4.D0*PI**2*RR*RMU0),
c$$$     &              (-TTPS(i)*DPPPS(i)/BB-DTTPS(i)*BB/RMU0)/(2.D0*PI)
c$$$      ENDDO
C
      return
c     
 2000 format (6a8,3i4)
 2020 format (5e16.9)
 2022 format (2i5)
c 2024 format (i5,e16.9,i5)
       end
C
C     ***** SAVE METRICS *****
C
      SUBROUTINE EQMETRIC(IERR)
C
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
         WRITE(nmetric,'(1X,F10.8,1P3E15.7)') RHOT(NR),DAT1,
     &        fnavir2(rhot(nr)),fnavgrr2(nr)
      ENDDO
      WRITE(nmetric,'(80X)')
      WRITE(nmetric,'(A,2X,A,3X,A,3X,A,5X,A,9X,A)') '#','rho_tor',
     &     "<|grad rho|>","<|grad rho|^2>","<B^2>","<1/B^2>"
      WRITE(nmetric,'(1X,0PF10.8,1P4E15.7)') (RHOT(NR),
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
