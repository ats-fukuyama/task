!     $Id$

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
     WRITE(25) ((NDELM(J,I),J=1,4),I=1,NEMAX)
  ELSEIF(ID.EQ.1) THEN
     WRITE(25) 'PAF-ELM1-V01'
     WRITE(25) NNMAX,NEMAX,NMMAX,NBMAX
     WRITE(25) (XND(I),YND(I),ZND(I),I=1,NNMAX)
     WRITE(25) (KANOD(I),I=1,NNMAX)
     WRITE(25) ((NDELM(J,I),J=1,4),I=1,NEMAX)
     WRITE(25) ((KNELM(J,I),J=1,4),I=1,NEMAX)
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
     READ(25,ERR=9100,END=9200) ((NDELM(J,I),J=1,4),I=1,NEMAX)

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
     READ(25,ERR=9100,END=9200) ((NDELM(J,I),J=1,4),I=1,NEMAX)
     READ(25,ERR=9100,END=9200) ((KNELM(J,I),J=1,4),I=1,NEMAX)
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

!
! ----- Add. By YOKOYAMA 01/03/2013 ----
!
!     G10磁力線座標の読み込み
!     Rcc=18cm (0,0,18cm)を出発点とする磁力線
!     主プラズマと周辺プラズマの境界線
!     セントラル中央に存在するリミタ(半径18cm)で切られる磁力管
!
!     FLZ: Z座標
!     FLX: ZX平面上における，Rcc=18cmの磁力線到達位置(X座標)
!     FLY: ZX平面上における，Rcc=18cmの磁力線到達位置(Y座標)
!     XYR: Y/X
!
!     ****** G10 B-FIELD LINE DATA ******
!
      SUBROUTINE G10BFL
!
      use wfcomm
      CHARACTER KGFN*15, KLINE1*30, KLINE2*30

      KGFN = 'G10BFLD_ZXY.txt'
!
      OPEN(UNIT=7,FILE=KGFN,STATUS='OLD',ACCESS='SEQUENTIAL',&
     &     FORM='FORMATTED',IOSTAT=IST,ERR=500)
!
         WRITE(6,100) KGFN
 100     FORMAT('## FILE (',1H ,A15,1H ,') IS ASSIGNED')
!
!        NGFLIN: ファイルに記述されているデータの行数
         READ(7,*,ERR=600) KLINE1, NGFLIN
         READ(7,*,ERR=600) KLINE2
         DO J=1,NGFLIN
            READ(7,*,ERR=600)  FLZ(J),FLX(J),FLY(J),XYR(J)
!            WRITE(6,'(4F9.4)') FLZ(J),FLX(J),FLY(J),XYR(J)
         ENDDO
 200     GOTO 900
!
 500     WRITE(6,*)'!! FILE OPEN ERROR : IOSTAT = ',IST
         STOP
 600     WRITE(6,*)'!! FILE READING ERROR'
         STOP
!
 900  CLOSE(7,STATUS='KEEP')
      RETURN
END SUBROUTINE G10BFL
!
!
!
!     磁力管断面内の波動電磁界を書き出す
!
!     ****** FIELD DATA OUT ******
!
      SUBROUTINE XYPROF
!
      use wfcomm
      CHARACTER KFNCEP*10
      DIMENSION ZPOINT(NNMAX),ESUM(2),BSUM(2)

!     格納用のファイルを開く
!     open files for saving
      KFNCEP = 'XYPROF.txt'
      OPEN(UNIT=8,FILE=KFNCEP,STATUS='OLD',ACCESS='SEQUENTIAL',&
     &     FORM='FORMATTED',IOSTAT=IST,ERR=500)
!
      WRITE(6,100) KFNCEP
 100  FORMAT('## FILE (',1H ,A10,1H ,') IS ASSIGNED')
      WRITE(8,*) '   NNOD   COUNT  ZPT(NOD)'//&
     &           '    ZPT(BFL)    E+       '//&
     &           '   E-          B+          B-'
!
!     節点のZ座標を抽出
!     extract z-coordinate of node
      NUMZND = 0
      DO NN=1,NNMAX
         IF((XND(NN).EQ.0.D0).and.&
     &      (YND(NN).EQ.0.D0)) THEN
            NUMZND = NUMZND + 1
            ZPOINT(NUMZND) = ZND(NN)
         ENDIF
      ENDDO
!
!     各Z座標において，磁力管断面上の点を探す
!     In each points, find the points on the flux tube cross section 
      NUMNOD = 0
      DO NZ=1,NUMZND
!        与えられたZ座標に最も近いZ座標を持つ磁力管断面を探す
!        FLZ,FLX,FLY ->  'cm' unit
!        ZPT,XPT,YPT ->  'm' unit
         DFZMIN = 1.D2
         DO J=1,NGFLIN
            REF = FLZ(J)/1.D2
            DFZ = ABS(ZPOINT(NZ)-REF)
            IF(DFZ.LT.DFZMIN) THEN
               DFZMIN = DFZ
               XPT = FLX(J)/1.D2
               YPT = FLY(J)/1.D2
               ZPT = FLZ(J)/1.D2
            ENDIF
         ENDDO
!
!        磁力管断面内の波動電磁界振幅の総和を計算する
!        総和は，磁力管断面内の節点数（加算の回数）で規格化する
         ESUM(1) = 0.D0
         ESUM(2) = 0.D0
         BSUM(1) = 0.D0
         BSUM(2) = 0.D0
         NCOUNT  = 0
         DO NN=1,NNMAX
            IF(ZND(NN).EQ.ZPOINT(NZ)) THEN
               NUMNOD = NUMNOD + 1
!              楕円の方程式より
!              According to the Elliptic equation
               ELLIPL  =  XND(NN)*XND(NN)/(XPT*XPT)&
     &                  + YND(NN)*YND(NN)/(YPT*YPT)
               IF(ELLIPL.LE.1.D0) THEN
                  NCOUNT = NCOUNT + 1
                  ESUM(1) = ESUM(1) + ABS(CEP(1,NN))
                  ESUM(2) = ESUM(2) + ABS(CEP(2,NN))
                  BSUM(1) = BSUM(1) + ABS(CBP(1,NN))
                  BSUM(2) = BSUM(2) + ABS(CBP(2,NN))
               ENDIF
            ENDIF
         ENDDO
!        加算の回数で規格化
!        Normalized by the number of times of addition
         EPOS = ESUM(1)/NCOUNT
         ENEG = ESUM(2)/NCOUNT
         BPOS = BSUM(1)/NCOUNT
         BNEG = BSUM(2)/NCOUNT
         WRITE(8,'(2I8,1P6E12.4)',ERR=600) &
     &      NUMNOD,NCOUNT,ZPOINT(NZ),&
     &      ZPT,EPOS,ENEG,BPOS,BNEG
      ENDDO
      GOTO 900
!
 500  WRITE(6,*)'!! FILE OPEN ERROR : IOSTAT = ',IST
      STOP
 600  WRITE(6,*)'!! FILE WRITING ERROR'
      STOP
!
 900  CLOSE(8,STATUS='KEEP')
      RETURN
END SUBROUTINE XYPROF
!
!
!     各パラメータのZ軸上における分布を書き出す
!     ****** FIELD DATA OUT ******
!
      SUBROUTINE AXPROF
!
      use wfcomm
      CHARACTER KFNCEP*10
      DIMENSION RNAX(2),AEP(2),ABP(2)

!     格納用のファイルを開く
!     opening file for saving
      KFNCEP = 'AXPROF.txt'
      OPEN(UNIT=8,FILE=KFNCEP,STATUS='OLD',ACCESS='SEQUENTIAL',&
     &     FORM='FORMATTED',IOSTAT=IST,ERR=500)
!
!      WRITE(6,100) KFNCEP
! 100  FORMAT('## FILE (',1H ,A10,1H ,') IS ASSIGNED')
!      WRITE(8,'(A111)',ERR=600)
!     &   ' ZPT(NOD)    BABS       '//
!     &   ' DEN.        ZPT(BFL)   '//
!     &   ' XPT(BFL)    YPT(BFL)   '//
!     &   ' E+          E-         '//
!     &   ' B+          B-'
!
      DO NN=1,NNMAX
!        節点のZ座標を抽出
         IF((XND(NN).EQ.0.D0).AND.&
     &      (YND(NN).EQ.0.D0)) THEN
!
!           与えられたZ座標に最も近いZ座標を持つ磁力管断面を探す
!           FLZ,FLX,FLYは cm 単位
!           ZPT,XPT,YPTは  m 単位
            DFZMIN = 1.D2
            DO J=1,NGFLIN
               REF = FLZ(J)/1.D2
               DFZ = ABS(ZND(NN)-REF)
               IF(DFZ.LT.DFZMIN) THEN
                  DFZMIN = DFZ
                  XPT = FLX(J)/1.D2
                  YPT = FLY(J)/1.D2
                  ZPT = FLZ(J)/1.D2 
              ENDIF
            ENDDO
!
!           軸上の磁場強度
            CALL WFSMAG(NN,BABSAX,DUMMY1)
!
!           軸上のプラズマ密度
            CALL WFSDEN(NN,RNAX,DUMMY2,DUMMY3,DUMMY4)
!
!           E+ E- B+ B- の絶対値
            AEP(1) = ABS(CEP(1,NN))
            AEP(2) = ABS(CEP(2,NN))
            ABP(1) = ABS(CBP(1,NN))
            ABP(2) = ABS(CBP(2,NN))
!
!           ファイルへの書き込み
!            WRITE(8,'(1P10(E12.4))',ERR=600) 
!     &         ZND(NN),BABSAX,RNAX(1),
!     &         ZPT,XPT,YPT,
            WRITE(8,200,ERR=600) ZND(NN),CEF(1,NN), &
                                 AEP(1),AEP(2), &
                                 ABP(1),ABP(2)
 200        FORMAT(1P7(E12.4,','))
         ENDIF
      ENDDO
      GOTO 900
!
 500  WRITE(6,*)'!! FILE OPEN ERROR : IOSTAT = ',IST
      STOP
 600  WRITE(6,*)'!! FILE WRITING ERROR : ',KFNCEP
      STOP
!
 900  CLOSE(8,STATUS='KEEP')
      RETURN
END SUBROUTINE AXPROF
!
!
!     ****** LOADING DATA OUT ******
!
      SUBROUTINE ANTIMP
!
      use wfcomm
      CHARACTER KFNCEP*10
      DIMENSION RNAX(2),AEP(2),ABP(2)

!     格納用のファイルを開く
!     opening file for saving
      KFNCEP = 'ANTIMP.txt'
      OPEN(UNIT=9,FILE=KFNCEP,STATUS='OLD',ACCESS='SEQUENTIAL',&
     &     FORM='FORMATTED',IOSTAT=IST,ERR=500)
!
!        ファイルへの書き込み
         WRITE(9,200,ERR=600) RF,ZANT,PABST,&
     &                        CLOAD(1),CLOAD(2),&
     &                        CLOAD(3),CLOAD(4),&
     &                        CLOAD(5),CLOAD(6),&
     &                        CLOAD(7),CLOAD(8),&
     &                        CLOAD(9),CLOAD(10)
 200     FORMAT(1P12(E12.4,','))
      GOTO 900
!
 500  WRITE(6,*)'!! FILE OPEN ERROR : IOSTAT = ',IST
      STOP
 600  WRITE(6,*)'!! FILE WRITING ERROR : ',KFNCEP
      STOP
!
 900  CLOSE(9,STATUS='KEEP')
      RETURN
END SUBROUTINE ANTIMP
!
!
!      節点上の波動電磁界を書き出す
!
!     ****** FIELD DATA OUT ******
!
      SUBROUTINE FLDOUT
!
      use wfcomm
      CHARACTER KFNCEP*10
      DIMENSION ZPOINT(NNMAX),ESUM(2),BSUM(2),RZCL(NSM),WC(NSM)
      DIMENSION CDUMMY(3,3,NSM)

!     格納用のファイルを開く
!     opening file for saving
      KFNCEP = 'FLDOUT.txt'
      OPEN(UNIT=23,FILE=KFNCEP,STATUS='OLD',ACCESS='SEQUENTIAL',&
     &     FORM='FORMATTED',IOSTAT=IST,ERR=500)

      WRITE(6,100) KFNCEP
 100  FORMAT('## FILE (',1H ,A10,1H ,') IS ASSIGNED')

!      WW=2.D0*PI*RF*1.D6

      DO NN=1,NNMAX
         CALL WFSMAG(NN,BFIELD,DUMMY1)
         CALL WFSDEN(NN,YDEN,DUMMY2,DUMMY3,RZCL)
         CALL DTENSR(NN,CDUMMY)
!         DO NS=1,NSMAX
!            WC(NS) = PZ(NS)*AEE/(PA(NS)*AMP*WW)
!            PCRT   = RZCL(NS)
!            PBKG   = 0.D0
!            WCI    = WC(NS)*WW*BFIELD
!            WDIF   = (WW/WCI-1.D0)*(WW/WCI-1.D0)
!            RZCL(NS) = PCRT*EXP(-500.D0*WDIF)+PBKG
!         ENDDO

         WRITE(23,'(I8,1P33E12.4)',ERR=600) &
!           N, X,      Y,      Z,      BABS     ,DENSITY
     &      NN,XND(NN),YND(NN),ZND(NN),BFIELD,YDEN(1),&
!           Ex,       Ey,       Ez
     &      CEF(1,NN),CEF(2,NN),CEF(3,NN),&
!           Bx,       By,       Bz
     &      CBF(1,NN),CBF(2,NN),CBF(3,NN),&
!           E-r,       E-theta,
     &      CERT(1,NN),CERT(2,NN), &
!           B-r,       B-theta,
     &      CBRT(1,NN),CBRT(2,NN), &
!           E+,       E-
     &      CEP(1,NN),CEP(2,NN),&
!           Electron,    Proton
     &      PABSSN(NN,1),PABSSN(NN,2),&
!           Electron,Proton
     &      RZCO(1,NN),RZCO(2,NN)
      ENDDO
!
      GOTO 900
 500  WRITE(6,*)'!! FILE OPEN ERROR : IOSTAT = ',IST
      STOP
 600  WRITE(6,*)'!! FILE WRITING ERROR'
      STOP
!
 900  CLOSE(23,STATUS='KEEP')
      RETURN
END SUBROUTINE FLDOUT
