C  
C     ***********************************************************
C
C           GRAPHIC 3D : CONTROL ROUTINE
C
C     ***********************************************************
C
      SUBROUTINE TRGRD0(K2,K3,INQ)
C
      INCLUDE 'trcomm.h'
      CHARACTER K2*1,K3*3,KK*4
      CHARACTER STR(88)*80,KV(88)*80
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
 200  CALL TRGRTD(STR,KV,NMB)
      CALL TRGRUR(G3D(1,1,NMB),STR(NMB),KV(NMB),INQ)
C
 100  RETURN
      END
C
C     **************************************************************
C
C           GRAPHIC 3D : UNIVERSAL ROUTINE
C
C     **************************************************************
C
      SUBROUTINE TRGRUR(GVD,STR,KV,INQ)
C
      INCLUDE 'trcomm.h'
      CHARACTER STR*80,KV*80
      DIMENSION GVD(NRM,NTM)
C
      GX1=3.0
      GX2=18.0
      GY1=2.0
      GY2=17.0
C
      CALL PAGES
      CALL TRGR3D(GX1,GX2,GY1,GY2,GRM,GT,GVD,NRM,NRMAX,NGT,
     &            STR,KV,2+INQ)
      CALL PAGEE
C
      RETURN
      END
C
C     **************************************************************
C
C           GRAPHIC 3D : GRAPH TITLE DATA
C
C     **************************************************************
C
      SUBROUTINE TRGRTD(STR,KV,NMB)
C
      INCLUDE 'trcomm.h'
      CHARACTER STR(88)*80,KV(88)*80
      CHARACTER STR0(88)*80,KV0(88)*80
C
      DATA STR0/'@TE [keV] vs t@','@TD [keV] vs t@',
     &          '@TT [keV] vs t@','@TA [keV] vs t@',
     &          '@NE [10^20/m^3] vs t@','@ND [10^20/m^3] vs t@',
     &          '@NT [10^20/m^3] vs t@','@NA [10^20/m^3] vs t@',
     &          '@AJ [MA/m^2] vs t@','@AJOH [MA/m^2] vs t@',
     &          '@AJNB [MA/m^2] vs t@','@AJRF [MA/m^2] vs t@',
     &          '@AJBS [MA/m^2] vs t@',
     &          '@PIN [MW/m^3] vs t@','@POH [MW/m^3] vs t@',
     &          '@PNB [MW/m^3] vs t@','@PNF [MW/m^3] vs t@',
     &          '@PRFE [MW/m^3] vs t@','@PRFD [MW/m^3] vs t@',
     &          '@PRFT [MW/m^3] vs t@','@PRFA [MW/m^3] vs t@',
     &          '@PRL [MW/m^3] vs t@','@PCX [MW/m^3] vs t@',
     &          '@PIE [MW/m^3] vs t@',
     &          '@QP vs t@','@EZOH [V/m] vs t@',
     &          '@BETA vs t@','@BETAP vs t@','@VLOOP [V] vs t@'/
C
      DATA KV0 /'@TE@','@TD@','@TT@','@TA@',
     &          '@NE@','@ND@','@NT@','@NA@',
     &          '@AJ@','@AJOH@','@AJNB@','@AJRF@','@AJBS@',
     &          '@PIN@','@POH@','@PNB@','@PNF@',
     &          '@PRFE@','@PRFD@','@PRFT@','@PRFA@',
     &          '@PRL@','@PCX@','@PIE@',
     &          '@QP@','@EZOH@','@BETA@','@BETAP@','@VLOOP@'/
      SAVE STR0, KV0
C
      STR(NMB)=STR0(NMB)
      KV (NMB)=KV0 (NMB)
C
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
      CHARACTER KK*4,KVP(88)*4
      COMMON KVP
C
      DATA KVP /'TE  ','TD  ','TT  ','TA  ',
     &          'NE  ','ND  ','NT  ','NA  ',
     &          'AJ  ','AJOH','AJNB','AJRF','AJBS',
     &          'PIN ','POH ','PNB ','PNF ',
     &          'PRFE','PRFD','PRFT','PRFA',
     &          'PRL ','PCX ','PIE ',
     &          'QP  ','EZOH','BETA','BETP','VLOP'/
C      SAVE KVP
C
      IERR=0
      DO NP=1,88
         IF (KK.EQ.KVP(NP)) THEN
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
      CHARACTER KK*4,KVP(88)*4
      COMMON KVP
C
      WRITE(6,700)
      WRITE(6,710) (KVP(I),I= 1, 5)
      WRITE(6,720) (KVP(I),I= 6,10)
      WRITE(6,730) (KVP(I),I=11,15)
      WRITE(6,740) (KVP(I),I=16,20)
      WRITE(6,750) (KVP(I),I=21,25)
      WRITE(6,760) (KVP(I),I=26,29)
C
 700  FORMAT(' ',
     &       '# 3D GRAPHICS LIST (EACH ONE IS CONSIST OF 4 LETTERS)')
 710  FORMAT(' ',' 1- 5: ',4(A4,', '),A4,)
 720  FORMAT(' ',' 6-10: ',4(A4,', '),A4,)
 730  FORMAT(' ','11-15: ',4(A4,', '),A4,)
 740  FORMAT(' ','16-20: ',4(A4,', '),A4,)
 750  FORMAT(' ','21-25: ',4(A4,', '),A4,)
 760  FORMAT(' ','26-29: ',3(A4,', '),A4,)
C
      RETURN
      END
