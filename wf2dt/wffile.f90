!     $Id: wffile.f90,v 1.1 2011/07/19 07:16:11 maruyama Exp $

!     ******* OUTPUT ELEMENT DATA *******

SUBROUTINE WFWELM(ID)

  use wfcomm
  implicit none
  integer :: IST,ID,I,J,NM,NB
  CHARACTER KNAME*32
  LOGICAL LEX
  
1 WRITE(6,*) '## ENTER ELEMENT FILE NAME: ',KFNAME
  READ(5,'(A32)',ERR=1,END=9000) KNAME
  IF(KNAME(1:2).NE.'/ ') KFNAME=KNAME
  INQUIRE(FILE=KFNAME,EXIST=LEX)
  IF(LEX) THEN
     OPEN(25,FILE=KFNAME,IOSTAT=IST,STATUS='OLD',ERR=10,&
          &        FORM='UNFORMATTED')
     WRITE(6,*) '## OLD FILE (',KFNAME,') IS ASSIGNED FOR OUTPUT.'
     GOTO 30
10   WRITE(6,*) 'XX OLD FILE OPEN ERROR : IOSTAT = ',IST
     GOTO 1
  ELSE
     OPEN(25,FILE=KFNAME,IOSTAT=IST,STATUS='NEW',ERR=20,&
          &        FORM='UNFORMATTED')
     WRITE(6,*) '## NEW FILE (',KFNAME,') IS CREATED FOR OUTPUT.'
     GOTO 30
20   WRITE(6,*) 'XX NEW FILE OPEN ERROR : IOSTAT = ',IST
     GOTO 1
  ENDIF
  
30 REWIND 25
  IF(ID.EQ.0) THEN
     WRITE(25) 'PAF-ELM0-V01'
     WRITE(25) NNMAX,NEMAX
     WRITE(25) (XND(I),YND(I),ZND(I),I=1,NNMAX)
     WRITE(25) ((NDELM(J,I),J=1,3),I=1,NEMAX)
  ELSEIF(ID.EQ.1) THEN
     WRITE(25) 'PAF-ELM1-V01'
     WRITE(25) NNMAX,NEMAX,NMMAX,NBMAX
     WRITE(25) (XND(I),YND(I),ZND(I),I=1,NNMAX)
     WRITE(25) (KANOD(I),I=1,NNMAX)
     WRITE(25) ((NDELM(J,I),J=1,3),I=1,NEMAX)
     WRITE(25) ((KNELM(J,I),J=1,3),I=1,NEMAX)
     WRITE(25) (KAELM(I),I=1,NEMAX)
     WRITE(25) (EPSDM(NM),AMUDM(NM),SIGDM(NM),NM=1,NMMAX)
     WRITE(25) (KABDY(NB),PHIBDY(NB),RESBDY(NB),&
          &              PWRBDY(NB),PHABDY(NB),          &
          &              XGBDY(NB),YGBDY(NB),ZGBDY(NB),  &
          &              (XNBDY(I,NB),YNBDY(I,NB),       &
          &               ZNBDY(I,NB),I=1,3),            &
          &              XPBDY(NB),YPBDY(NB),ZPBDY(NB),  &
          &              SZBDY(1,NB),SZBDY(2,NB),NB=1,NBMAX)
  ENDIF
  CLOSE(25)
  
  WRITE(6,*) '## ELEMENT DATA SAVED IN FILE: ',KFNAME

9000 RETURN
END SUBROUTINE WFWELM

!     ******* INPUT ELEMENT DATA *******

SUBROUTINE WFRELM(ID)

  use wfcomm
  implicit none
  integer :: IST,I,J,NE,NN,NM,NB,ID,idata(2)
  LOGICAL LEX
  CHARACTER KID*12

  INQUIRE(FILE=KFNAME,EXIST=LEX)
  IF(LEX) THEN
     OPEN(25,FILE=KFNAME,IOSTAT=IST,STATUS='OLD',ERR=10,&
          &        FORM='UNFORMATTED')
     if(nrank.eq.0) WRITE(6,*) '## FILE (',KFNAME,') IS ASSIGNED FOR ELM INPUT'
     GOTO 30
10   if(nrank.eq.0) WRITE(6,*) 'XX OLD FILE OPEN ERROR : IOSTAT = ',IST
     GOTO 9000
  ELSE
     if(nrank.eq.0) WRITE(6,*) 'XX FILE (',KFNAME,') IS NOT FOUND'
     GOTO 9000
  ENDIF

30 READ(25,ERR=9100,END=9000) KID

  IF(KID.EQ.'PAF-ELM0-V01') THEN
     ID=0
     
     call wfelm_allocate
     READ(25,ERR=9100,END=9200) NNMAX,NEMAX
     call wfelm_allocate
     
     READ(25,ERR=9100,END=9200) (XND(I),YND(I),ZND(I),I=1,NNMAX)
     READ(25,ERR=9100,END=9200) ((NDELM(J,I),J=1,3),I=1,NEMAX)

     CALL WFINDX
     CALL WFFEPI
     
     NKMAX=1
     DO NE=1,NEMAX
        KAELM(NE)=1
     ENDDO
     NMKA(1)=0
     NMMAX=0
     
     NBMAX=0
     DO NN=1,NNMAX
        KANOD(NN)=0
     ENDDO

  ELSEIF(KID.EQ.'PAF-ELM1-V01') THEN
     ID=1
     READ(25,ERR=9100,END=9200) NNMAX,NEMAX,NMMAX,NBMAX
     call wfelm_allocate
     
     READ(25,ERR=9100,END=9200) (XND(I),YND(I),ZND(I),I=1,NNMAX)
     READ(25,ERR=9100,END=9200) (KANOD(I),I=1,NNMAX)
     READ(25,ERR=9100,END=9200) ((NDELM(J,I),J=1,3),I=1,NEMAX)
     READ(25,ERR=9100,END=9200) ((KNELM(J,I),J=1,3),I=1,NEMAX)
     READ(25,ERR=9100,END=9200) (KAELM(I),I=1,NEMAX)
     READ(25,ERR=9100,END=9200) (EPSDM(NM),AMUDM(NM),SIGDM(NM),&
          &                              NM=1,NMMAX)
     READ(25,ERR=9100,END=9200) (KABDY(NB),PHIBDY(NB),RESBDY(NB),&
          &                               PWRBDY(NB),PHABDY(NB),&
          &                               XGBDY(NB),YGBDY(NB),ZGBDY(NB),&
          &                               (XNBDY(I,NB),YNBDY(I,NB),&
          &                                ZNBDY(I,NB),I=1,3),&
          &                               XPBDY(NB),YPBDY(NB),ZPBDY(NB),&
          &                               SZBDY(1,NB),SZBDY(2,NB),NB=1,NBMAX)
  ELSE
     ID=9999
     if(nrank.eq.0) WRITE(6,*) 'XX INVALID ELEMENT DATA : KID = ',KID
     CLOSE(25)
     GOTO 9000
  ENDIF
  CLOSE(25)
     
9000 RETURN
  
9100 if(nrank.eq.0) WRITE(6,*) 'XX ERROR IN READING FILE: ',KFNAME
  GOTO 9000
  
9200 if(nrank.eq.0) WRITE(6,*) 'XX UNEXPECTED END OF FILE: ',KFNAME
  GOTO 9000

END SUBROUTINE WFRELM

!     ******* OUTPUT ANTENNA DATA *******

SUBROUTINE WFWANT

  use wfcomm
  implicit none
  integer :: IST,NA,I
  CHARACTER KNAME*32
  LOGICAL LEX

1 WRITE(6,*) '## ENTER ANTENNA FILE NAME: ',KFNAMA
  READ(5,'(A32)',ERR=1,END=9000) KNAME
  IF(KNAME(1:2).NE.'/ ') KFNAMA=KNAME
  INQUIRE(FILE=KFNAMA,EXIST=LEX)
  IF(LEX) THEN
     OPEN(24,FILE=KFNAMA,IOSTAT=IST,STATUS='OLD',ERR=10,&
          &        FORM='FORMATTED')
     WRITE(6,*) '## OLD FILE (',KFNAMA,') IS ASSIGNED FOR OUTPUT.'
     GOTO 30
10   WRITE(6,*) 'XX OLD FILE OPEN ERROR : IOSTAT = ',IST
     GOTO 1
  ELSE
     OPEN(24,FILE=KFNAMA,IOSTAT=IST,STATUS='NEW',ERR=20,&
          &        FORM='FORMATTED')
     WRITE(6,*) '## NEW FILE (',KFNAMA,') IS CREATED FOR OUTPUT.'
     GOTO 30
20   WRITE(6,*) 'XX NEW FILE OPEN ERROR : IOSTAT = ',IST
     GOTO 1
  ENDIF
  
30 REWIND 24
  WRITE(24,200) NAMAX
200 FORMAT(I6)
  DO NA=1,NAMAX
     WRITE(24,200) JNUM0(NA)
     WRITE(24,300) (XJ0(I,NA),YJ0(I,NA),ZJ0(I,NA),I=1,JNUM0(NA))
300  FORMAT(3E23.15)
  END DO
  CLOSE(24)
  
  WRITE(6,*) '## ANTENNA DATA SAVED IN FILE: ',KFNAMA
  
9000 RETURN
END SUBROUTINE WFWANT

!     ******* INPUT ANTENNA DATA FROM FILE 24 *******

SUBROUTINE WFRANT

  use wfcomm
  implicit none
  integer :: IST,NA,I,IERR
  LOGICAL LEX
  
  INQUIRE(FILE=KFNAMA,EXIST=LEX)
  IF(LEX) THEN
     OPEN(24,FILE=KFNAMA,IOSTAT=IST,STATUS='OLD',ERR=10,&
          &        FORM='FORMATTED')
     if(nrank.eq.0) WRITE(6,*) '## FILE (',KFNAMA,') IS ASSIGNED FOR ANT INPUT'
     GOTO 30
10   if(nrank.eq.0) WRITE(6,*) 'XX OLD FILE OPEN ERROR : IOSTAT = ',IST
     GOTO 9000
  ELSE
     if(nrank.eq.0) WRITE(6,*) 'XX FILE (',KFNAMA,') IS NOT FOUND'
     GOTO 9000
  ENDIF
  
30 READ(24,150,ERR=9100,END=9200) NAMAX
150 FORMAT(I6)
  IF(NAMAX.LT.1.OR.NAMAX.GT.NAM) THEN
     if(nrank.eq.0) WRITE(6,*) 'XX INVALID ANTENNA DATA: NAMAX,NAM = ',NAMAX,NAM
     GOTO 9000
  ENDIF
  
  DO NA=1,NAMAX
     READ(24,150,ERR=9100,END=9200) JNUM0(NA)
     IF(JNUM0(NA).GT.NJM.OR.JNUM0(NA).LT.1) THEN
        if(nrank.eq.0) WRITE(6,*) 'XX INVALID ANTENNA DATA: NA,JNUM0,NJM =',&
             &                 NA,JNUM0(NA),NJM
        GOTO 9000
     ENDIF
     
     READ(24,250,ERR=9100,END=9200) &
          &       (XJ0(I,NA),YJ0(I,NA),ZJ0(I,NA),I=1,JNUM0(NA))
250  FORMAT(3E23.15)
  END DO
  CLOSE(24)
  
  CALL MODANT(IERR)
  
9000 RETURN
  
9100 if(nrank.eq.0) WRITE(6,*) 'XX ERROR IN READING FILE: ',KFNAMA
  GOTO 9000
  
9200 if(nrank.eq.0) WRITE(6,*) 'XX UNEXPECTED END OF FILE: ',KFNAMA
  GOTO 9000
END SUBROUTINE WFRANT

!     ******* OUTPUT FIELD DATA *******

SUBROUTINE WFWFLD

  use wfcomm
  implicit none
  integer :: IST,M
  CHARACTER KNAME*32
  LOGICAL LEX

  IF(NFOPEN.EQ.0) THEN
        
1    WRITE(6,*) '## ENTER FIELD FILE NAME: ',KFNAMF
     READ(5,'(A32)',ERR=1,END=9000) KNAME
     IF(KNAME(1:2).NE.'/ ') KFNAMF=KNAME
     INQUIRE(FILE=KFNAMF,EXIST=LEX)
     IF(LEX) THEN
        OPEN(26,FILE=KFNAMF,IOSTAT=IST,STATUS='OLD',ERR=10,&
             &        FORM='UNFORMATTED')
        WRITE(6,*) '## OLD FILE (',KFNAMF,') IS ASSIGNED FOR OUTPUT.'
        GOTO 30
10      WRITE(6,*) 'XX OLD FILE OPEN ERROR : IOSTAT = ',IST
        GOTO 1
     ELSE
        OPEN(26,FILE=KFNAMF,IOSTAT=IST,STATUS='NEW',ERR=20,&
             &        FORM='UNFORMATTED')
        WRITE(6,*) '## NEW FILE (',KFNAMF,') IS CREATED FOR OUTPUT.'
        GOTO 30
20      WRITE(6,*) 'XX NEW FILE OPEN ERROR : IOSTAT = ',IST
        GOTO 1
     ENDIF
30   NFOPEN=1
     
  ENDIF
  
  WRITE(26) MLEN
  WRITE(26) (CSV(M),M=1,MLEN)
  
  WRITE(6,*) '## FIELD DATA SAVED IN FILE: ',KFNAMF
  
9000 RETURN
END SUBROUTINE WFWFLD

!     ****** INPUT FIELD DATA ******

SUBROUTINE WFRFLD

  use wfcomm
  implicit none
  integer :: IST,IERR,M
  LOGICAL LEX
  
  IF(NFOPEN.EQ.0) THEN
     
     INQUIRE(FILE=KFNAMF,EXIST=LEX)
     IF(LEX) THEN
        OPEN(26,FILE=KFNAMF,IOSTAT=IST,STATUS='OLD',ERR=10,&
             &        FORM='UNFORMATTED')
        if(nrank.eq.0) WRITE(6,*) '## FILE (',KFNAMF,') IS ASSIGNED FOR FLD INPUT'
        GOTO 30
10      if(nrank.eq.0) WRITE(6,*) 'XX OLD FILE OPEN ERROR : IOSTAT = ',IST
        GOTO 9000
     ELSE
        if(nrank.eq.0) WRITE(6,*) 'XX FILE (',KFNAMF,') IS NOT FOUND'
        GOTO 9000
     ENDIF
30   NFOPEN=1
     
  ENDIF
  
  READ(26,ERR=9100,END=9200) MLEN
  READ(26,ERR=9100,END=9200) (CSV(M),M=1,MLEN)
  
  if(nrank.eq.0) WRITE(6,*) '--- WFWPRE started ---'
  CALL WFWPRE(IERR)
  IF(IERR.NE.0) GOTO 9000
  
  CALL CALFLD
  CALL PWRABS
  CALL PWRRAD
  CALL TERMEP
  CALL WFCALB
  CALL LPEFLD
  
9000 RETURN
  
9100 if(nrank.eq.0) WRITE(6,*) 'XX ERROR IN READING FILE: ',KFNAMF
  GOTO 9000
  
9200 if(nrank.eq.0) WRITE(6,*) 'XX UNEXPECTED END OF FILE: ',KFNAMF
  GOTO 9000

END SUBROUTINE WFRFLD
