C     $Id$
C
C     ***** READ EQDSK FORMAT FILE *****
C
      SUBROUTINE EQDSKR(IERR)
C
      INCLUDE '../eq/eqcomq.inc'
C
      character case(6)*10
      dimension rlim(NSUM),zlim(NSUM),pressw(NPSM),pwprim(NPSM),
     &          dmion(NSUM),rhovn(NSUM)
      CHARACTER KLINE*80
C
      LOGICAL LEX
C
      IERR=0
      neqdsk=21
      IF(KNAMEQ(1:2).EQ.'  ') GOTO 9000
      INQUIRE(FILE=KNAMEQ,EXIST=LEX,ERR=9000)
      IF(LEX) THEN
         OPEN(neqdsk,FILE=KNAMEQ,IOSTAT=IST,STATUS='OLD',ERR=10,
     &        FORM='FORMATTED')
         WRITE(6,*) '# OLD FILE (',KNAMEQ,') IS ASSIGNED FOR INPUT.'
         GOTO 30
   10    WRITE(6,*) 'XX OLD FILE OPEN ERROR : IOSTAT = ',IST
         GOTO 9000
      ELSE
         WRITE(6,*) 'XX FILE (',KNAMEQ,') NOT FOUND'
         GOTO 9000
      ENDIF
c
   30 read (neqdsk,2000) (case(i),i=1,6),idum,NRGMAX,NZGMAX
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
      read (neqdsk,2020) RAXIS,ZAXIS,SAXIS,SBNDY,BB
      DPS=(SBNDY-SAXIS)/(NPSMAX-1)
      DO NPS=1,NPSMAX
         PSIPS(NPS)=SAXIS+DPS*(NPS-1)
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
      read (neqdsk,2024) kvtor,rvtor,nmass
      if (kvtor.gt.0) then
         read (neqdsk,2020) (pressw(i),i=1,NPSMAX)
         read (neqdsk,2020) (pwprim(i),i=1,NPSMAX)
      endif
      if (nmass.gt.0) then
         read (neqdsk,2020) (dmion(i),i=1,NPSMAX)
      endif
      read (neqdsk,2020) (rhovn(i),i=1,NPSMAX)
C
C      WRITE(6,'(1P3E12.4)') RR,BB,RIP
C      WRITE(6,'(1P4E12.4)') RAXIS,ZAXIS,SAXIS,SBNDY
C      WRITE(6,'(1P4E12.4)') RG(1),RG(2),RG(NRGMAX-1),RG(NRGMAX)
C      WRITE(6,'(1P4E12.4)') ZG(1),ZG(2),ZG(NZGMAX-1),ZG(NZGMAX)
C      WRITE(6,'(1P4E12.4)') PSIPS(1),PSIPS(2),
C     &                      PSIPS(NPSMAX-1),PSIPS(NPSMAX)
      return
c     
 2000 format (6a8,3i4)
 2020 format (5e16.9)
 2022 format (2i5)
 2024 format (i5,e16.9,i5)
C
 9000 IERR=1
      RETURN
       end
