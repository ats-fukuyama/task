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
C for negative Ip and negative BB
      IF(RIP.LT.0.D0) RIP=-RIP
      IF(Bctr.LT.0.D0) Bctr=-Bctr
      IF(TTPS(1).LT.0.D0) THEN
         DO NPS=1,NPSMAX
            TTPS(NPS)=-TTPS(NPS)
         ENDDO
      ENDIF

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
