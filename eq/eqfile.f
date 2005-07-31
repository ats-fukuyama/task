C     $Id$
C
C     ***** SAVE TASK/EQ DATA *****
C
      SUBROUTINE EQSAVE
C
      INCLUDE '../eq/eqcomc.inc'
C
      CALL FWOPEN(21,KNAMEQ,0,1,'EQ',IERR)
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
      WRITE(21) NSGMAX,NTGMAX
      WRITE(21) RA,RKAP,RDLT,RB
      WRITE(21) PJ0,PJ1,PJ2,PROFJ0,PROFJ1,PROFJ2
      WRITE(21) PP0,PP1,PP2,PROFP0,PROFP1,PROFP2
      WRITE(21) PT0,PT1,PT2,PROFT0,PROFT1,PROFT2
      WRITE(21) PV0,PV1,PV2,PROFV0,PROFV1,PROFV2
      WRITE(21) PROFR0,PROFR1,PROFR2
      WRITE(21) PTS,PN0,HM
      CLOSE(21)
C
      WRITE(6,*) '# DATA WAS SUCCESSFULLY SAVED TO THE FILE.'
C
      RETURN
      END
C
C     ***** LOAD EQUILIBRIUM DATA *****
C
      SUBROUTINE EQLOAD(MODELG,KNAMEQ1,IERR)
C
      INCLUDE '../eq/eqcomc.inc'
      CHARACTER KNAMEQ1*80
C
      PI=2.D0*ASIN(1.D0)
      RMU0=4.D0*PI*1.D-7
      KNAMEQ=KNAMEQ1
C
      IF(MODELG.EQ.3) THEN
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
C
      CALL FROPEN(21,KNAMEQ,0,0,'EQ',IERR)
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
      READ(21) NSGMAX,NTGMAX
      READ(21) RA,RKAP,RDLT,RB
      READ(21) PJ0,PJ1,PJ2,PROFJ0,PROFJ1,PROFJ2
      READ(21) PP0,PP1,PP2,PROFP0,PROFP1,PROFP2
      READ(21) PT0,PT1,PT2,PROFT0,PROFT1,PROFT2
      READ(21) PV0,PV1,PV2,PROFV0,PROFV1,PROFV2
      READ(21) PROFR0,PROFR1,PROFR2
C     READ(21) PTS,PN0,HM
      CLOSE(21)
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
     &          dmion(NSUM),rhovn(NSUM)
C
      neqdsk=21
      CALL FROPEN(neqdsk,KNAMEQ,1,0,'EQ',IERR)
      IF(IERR.NE.0) RETURN
c
      REWIND(neqdsk)
      read (neqdsk,2000) (case(i),i=1,6),idum,NRGMAX,NZGMAX
      NPSMAX=NRGMAX
      read (neqdsk,2020) rdim,zdim,RR,rleft,zmid
      RA=0.5D0*rdim
      RKAP=zdim/rdim
      RDLT=0.D0
      RB=1.1D0*RA
      DR=rdim/(NRGMAX-1)
      DZ=zdim/(NZGMAX-1)
      DO NRG=1,NRGMAX
         RG(NRG)=rleft+DR*(NRG-1)
      ENDDO
      DO NZG=1,NZGMAX
         ZG(NZG)=zmid-0.5D0*zdim+DZ*(NZG-1)
      ENDDO
      read (neqdsk,2020) RAXIS,ZAXIS,PSI0,PSIS,BB
C
      PSI0=2.D0*PI*PSI0
      PSIS=2.D0*PI*PSIS
C
      DPS=(PSI0-PSIS)/(NPSMAX-1)
      PSI0=PSI0-PSIS
      PSIPA=-PSI0
      DO NPS=1,NPSMAX
         PSIPS(NPS)=PSI0-DPS*(NPS-1)
      ENDDO
      read (neqdsk,2020) RIP,simag,xdum,rmaxis,xdum
      RIP=RIP/1.D6
      read (neqdsk,2020) zmaxis,xdum,sibry,xdum,xdum
      read (neqdsk,2020) (TTPS(i),i=1,NPSMAX)
      read (neqdsk,2020) (PPPS(i),i=1,NPSMAX)
      read (neqdsk,2020) (DTTPS(i),i=1,NPSMAX)
      read (neqdsk,2020) (DPPPS(i),i=1,NPSMAX)
      read (neqdsk,2020) ((PSIRZ(i,j),i=1,NRGMAX),j=1,NZGMAX)
      read (neqdsk,2020) (QQPS(i),i=1,NPSMAX)
      read (neqdsk,2022) NSUMAX,limitr
      read (neqdsk,2020) (RSU(i),ZSU(i),i=1,NSUMAX)
      read (neqdsk,2020) (rlim(i),zlim(i),i=1,limitr)
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
C
C      WRITE(6,'(1P3E12.4)') RR,BB,RIP
C      WRITE(6,'(1P4E12.4)') RAXIS,ZAXIS,PSI0,PSIS
C      WRITE(6,'(1P4E12.4)') RG(1),RG(2),RG(NRGMAX-1),RG(NRGMAX)
C      WRITE(6,'(1P4E12.4)') ZG(1),ZG(2),ZG(NZGMAX-1),ZG(NZGMAX)
C      WRITE(6,'(1P4E12.4)') PSIPS(1),PSIPS(2),
C     &                      PSIPS(NPSMAX-1),PSIPS(NPSMAX)
C
      DO NZG=1,NZGMAX
         DO NRG=1,NRGMAX
            PSIRZ(NRG,NZG)=2.D0*PI*PSIRZ(NRG,NZG)-PSIS
         ENDDO
      ENDDO
      DO i=1,NPSMAX
         TTPS(i) =2.D0*PI*TTPS(i)
         DTTPS(i)=2.D0*PI*DTTPS(i)
      ENDDO
      return
c     
 2000 format (6a8,3i4)
 2020 format (5e16.9)
 2022 format (2i5)
 2024 format (i5,e16.9,i5)
       end
