C     $Id$
C  
C     ***********************************************************
C
C           GRAPHIC 3D : CONTROL ROUTINE
C
C     ***********************************************************
C
      SUBROUTINE TRGRD0(K2,K3,INQ)
C
      INCLUDE 'trcomm.inc'
      CHARACTER K2*1,K3*3,KK*4
      CHARACTER STRL*80,KVL*80
C
      READ(K2,'(I1)',ERR=600) I2
      READ(K3,'(I1)') I3
      IF (K3.EQ.' ') THEN
         NMB=I2
         IF (K2.EQ.' ') THEN
            CALL VIEW3DLIST
            GOTO 100
         ENDIF
      ELSE
         NMB=I2*10+I3
      ENDIF
      GOTO 200
C
 600  KK=K2//K3
      CALL CHNBPR(KK,NMB,IERR)
      IF (IERR.EQ.1) GOTO 100
C
 200  CALL TRGRTD(STRL,KVL,NMB)
      IF(NMB.GE.61) NMB=60
      CALL TRGRUR(G3D(1,1,NMB),STRL,KVL,INQ)
C
 100  RETURN
      END
C
C     **************************************************************
C
C           GRAPHIC 3D : GRAPH TITLE DATA
C
C     **************************************************************
C
      SUBROUTINE TRGRTD(STRL,KVL,NMB)
C
      INCLUDE 'trcomm.inc'
      CHARACTER STRL*80,KVL*80
      CHARACTER STR0(60)*80,KV0(60)*80
C
      DATA STR0/'@TE [keV] vs t@','@TD [keV] vs t@',
     &          '@TT [keV] vs t@','@TA [keV] vs t@',
     &          '@NE [10^20/m^3] vs t@','@ND [10^20/m^3] vs t@',
     &          '@NT [10^20/m^3] vs t@','@NA [10^20/m^3] vs t@',
     &          '@AJ [A/m^2] vs t@','@AJOH [A/m^2] vs t@',
     &          '@AJNB [A/m^2] vs t@','@AJRF [A/m^2] vs t@',
     &          '@AJBS [A/m^2] vs t@',
     &          '@PIN [W/m^3] vs t@','@POH [W/m^3] vs t@',
     &          '@PNB [W/m^3] vs t@','@PNF [W/m^3] vs t@',
     &          '@PRFE [W/m^3] vs t@','@PRFD [W/m^3] vs t@',
     &          '@PRFT [W/m^3] vs t@','@PRFA [W/m^3] vs t@',
     &          '@PRL [W/m^3] vs t@','@PCX [W/m^3] vs t@',
     &          '@PIE [W/m^3] vs t@',
     &          '@PEXE [W/m^3] vs t@','@PEXI [W/m^3] vs t@',
     &          '@QP vs t@','@EZOH [V/m] vs t@',
     &          '@BETA vs t@','@BETAP vs t@','@VLOOP [V] vs t@',
     &          '@ETA [Ohm*m] vs t@','@ZEFF vs t@',
     &          '@AKE [m^2/s] vs t@','@AKD [m^2/s] vs t@',
     &          '@PRECE [W/m^3] vs t@','@PRLHE [W/m^3] vs t@',
     &          '@PRICE [W/m^3] vs t@','@PRECI [W/m^3] vs t@',
     &          '@PRLHI [W/m^3] vs t@','@PRICI [W/m^3] vs t@',
     &          '@AJEC [A/m^2] vs t@','@AJLH [A/m^2] vs t@',
     &          '@AJIC [A/m^2] vs t@',
     &          '@NFAST [10^20/m^3] vs t@','@NIMP [10^20/m^3] vs t@',
     &          '@BPOL [T] vs t@','@PSI [Wb] vs t@',
     &          '@RMAJOR [m] vs t@','@RMINOR [m] vs t@',
     &          '@VOLUME [m^3] sv t@','@KAPPAR vs t@',
     &          '@DELTAR@',
     &          '@GRHO1 vs t@','@GRHO2 vs t@',
     &          '@AKDWE vs t@','@AKDWI vs t@',
     &          '@PE [MPa] vs t@','@PI [MPa] vs t@',
     &          '@VTOR [m/s] vs t@'/
C
      DATA KV0 /'@TE@','@TD@','@TT@','@TA@',
     &          '@NE@','@ND@','@NT@','@NA@',
     &          '@AJ@','@AJOH@','@AJNB@','@AJRF@','@AJBS@',
     &          '@PIN@','@POH@','@PNB@','@PNF@',
     &          '@PRFE@','@PRFD@','@PRFT@','@PRFA@',
     &          '@PRL@','@PCX@','@PIE@','@PEXE@','@PEXI@',
     &          '@QP@','@EZOH@','@BETA@','@BETAP@','@VLOOP@',
     &          '@ETA@','@ZEFF@','@AKE@','@AKD@',
     &          '@PRECE@','@PRLHE@','@PRICE@',
     &          '@PRECI@','@PRLHI@','@PRICI@',
     &          '@AJEC@','@AJLH@','@AJIC@',
     &          '@NFAST@','@NIMP@','@BPOL@','@PSI@',
     &          '@RMAJOR@','@RMINOR@','@VOLUME@','@KAPPAR@',
     &          '@DELTAR@','@GRHO1@','@GRHO2@',
     &          '@AKTBE@','@AKTBI@','@PE@','@PI@','@VTOR@'/
      SAVE STR0, KV0
C
      STRL=STR0(NMB)
      KVL =KV0 (NMB)
C
      RETURN
      END
C
C     **************************************************************
C
C           GRAPHIC 3D : PARAMETER STRING FOR NP
C
C     **************************************************************
C
      SUBROUTINE GETVPL(KVPL,NP)
C
      CHARACTER KVPL*4,KVP(60)*4
C
      DATA KVP /'TE  ','TD  ','TT  ','TA  ',
     &          'NE  ','ND  ','NT  ','NA  ',
     &          'AJ  ','AJOH','AJNB','AJRF','AJBS',
     &          'PIN ','POH ','PNB ','PNF ',
     &          'PRFE','PRFD','PRFT','PRFA',
     &          'PRL ','PCX ','PIE ','PEXE','PEXI',
     &          'QP  ','EZOH','BETA','BETP','VLOP',
     &          'ETA ','ZEFF','AKE ','AKD ',
     &          'PREE','PRLE','PRIE','PREI','PRLI','PRII',
     &          'AJEC','AJLH','AJIC',
     &          'NFST','NIMP','BP  ','PSI',
     &          'RMJ ','RMN ','VOL ','LKAP','DLT ','GRH1','GRH2',
     &          'AKTE','AKTI','PE','PI','VTOR'/
C      SAVE KVP
C
      KVPL=KVP(NP)
      RETURN
      END
C
C     **************************************************************
C
C           GRAPHIC 3D : CHECK A NUMBER CORRESPONDING A PARAMETER
C
C     **************************************************************
C
      SUBROUTINE CHNBPR(KK,NMB,IERR)
C
      CHARACTER KK*4,KVPL*4
C
      IERR=0
      DO NP=1,60
         CALL GETVPL(KVPL,NP)
         IF (KK.EQ.KVPL) THEN
            NMB=NP
            GOTO 1000
         ENDIF
      ENDDO
      IERR=1
C
 1000 RETURN
      END
C
C     **************************************************************
C
C           GRAPHIC 3D : VIEW 3D GRAPHICS LIST 
C
C     **************************************************************
C
      SUBROUTINE VIEW3DLIST
C
      CHARACTER KVP(60)*4
C
      DO I=1,60
         CALL GETVPL(KVP(I),I)
      ENDDO
      WRITE(6,700)
      WRITE(6,710) (KVP(I),I= 1, 5)
      WRITE(6,720) (KVP(I),I= 6,10)
      WRITE(6,730) (KVP(I),I=11,15)
      WRITE(6,740) (KVP(I),I=16,20)
      WRITE(6,750) (KVP(I),I=21,25)
      WRITE(6,760) (KVP(I),I=26,30)
      WRITE(6,770) (KVP(I),I=31,35)
      WRITE(6,780) (KVP(I),I=36,40)
      WRITE(6,790) (KVP(I),I=41,45)
      WRITE(6,800) (KVP(I),I=46,50)
      WRITE(6,810) (KVP(I),I=51,55)
      WRITE(6,820) (KVP(I),I=56,60)
C
 700  FORMAT(' ',
     &       '# 3D GRAPHICS LIST (EACH ONE IS CONSIST OF 4 LETTERS)')
 710  FORMAT(' ',' 1- 5: ',4(A4,', '),A4)
 720  FORMAT(' ',' 6-10: ',4(A4,', '),A4)
 730  FORMAT(' ','11-15: ',4(A4,', '),A4)
 740  FORMAT(' ','16-20: ',4(A4,', '),A4)
 750  FORMAT(' ','21-25: ',4(A4,', '),A4)
 760  FORMAT(' ','26-30: ',4(A4,', '),A4)
 770  FORMAT(' ','31-35: ',4(A4,', '),A4)
 780  FORMAT(' ','36-40: ',4(A4,', '),A4)
 790  FORMAT(' ','41-45: ',4(A4,', '),A4)
 800  FORMAT(' ','46-50: ',4(A4,', '),A4)
 810  FORMAT(' ','51-55: ',4(A4,', '),A4)
 820  FORMAT(' ','56-60: ',4(A4,', '),A4)
C
      RETURN
      END
