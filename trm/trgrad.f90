!     ***********************************************************

!           GRAPHIC 3D : CONTROL ROUTINE

!     ***********************************************************

      SUBROUTINE TRGRD0(K2,K3,INQ)

      USE TRCOMM, ONLY : NT, NGTSTP, NCRTM
      IMPLICIT NONE
      CHARACTER(LEN=1), INTENT(IN):: K2
      CHARACTER(LEN=3), INTENT(IN):: K3
      INTEGER,       INTENT(IN):: INQ
      INTEGER :: I2, I3, IERR, NMB
      CHARACTER(LEN=4) :: KK
      CHARACTER(LEN=80):: KVL, STRL

      READ(K2,'(I1)',ERR=600) I2
      READ(K3,'(I1)') I3
      IF (K3.EQ.' ') THEN
         NMB=I2
         IF (K2.EQ.' ') THEN
            CALL VIEWRTLIST(NCRTM)
            GOTO 100
         ENDIF
      ELSE
         NMB=I2*10+I3
      ENDIF
      IF(NT.LT.NGTSTP) RETURN
      GOTO 200

 600  KK=K2//K3
      CALL CHNBPR(KK,NMB,IERR)
      IF (IERR.EQ.1) GOTO 100

 200  CALL TRGRTD(STRL,KVL,NMB)
      IF(NMB.GE.68) NMB=67
      CALL TRGRUR(NMB,STRL,KVL,INQ)

 100  RETURN
      END SUBROUTINE TRGRD0

!     **************************************************************

!           GRAPHIC 3D : GRAPH TITLE DATA

!     **************************************************************

      SUBROUTINE TRGRTD(STRL,KVL,NMB)

      IMPLICIT NONE
      CHARACTER(LEN=80),INTENT(OUT):: STRL, KVL
      INTEGER,       INTENT(IN) :: NMB
      CHARACTER(LEN=80), DIMENSION(67), SAVE:: STR0, KV0
      DATA STR0/ &
           '@TE [keV] vs t@','@TD [keV] vs t@', &
           '@TT [keV] vs t@','@TA [keV] vs t@', &
           '@NE [10^20/m^3] vs t@','@ND [10^20/m^3] vs t@', &
           '@NT [10^20/m^3] vs t@','@NA [10^20/m^3] vs t@', &
           '@AJ [A/m^2] vs t@','@AJOH [A/m^2] vs t@', &
           '@AJNB [A/m^2] vs t@','@AJRF [A/m^2] vs t@', &
           '@AJBS [A/m^2] vs t@','@PIN [W/m^3] vs t@', &
           '@POH [W/m^3] vs t@','@PNB [W/m^3] vs t@', &
           '@PNF [W/m^3] vs t@', '@PRFE [W/m^3] vs t@', &
           '@PRFD [W/m^3] vs t@','@PRFT [W/m^3] vs t@', &
           '@PRFA [W/m^3] vs t@', '@PRL [W/m^3] vs t@', &
           '@PCX [W/m^3] vs t@','@PIE [W/m^3] vs t@', &
           '@PEXE [W/m^3] vs t@','@PEXI [W/m^3] vs t@',&
           '@QP vs t@','@EZOH [V/m] vs t@', &
           '@BETA vs t@','@BETAP vs t@', &
           '@VLOOP [V] vs t@', '@ETA [Ohm*m] vs t@', &
           '@ZEFF vs t@','@AKE [m^2/s] vs t@', &
           '@AKD [m^2/s] vs t@', '@PRECE [W/m^3] vs t@', &
           '@PRLHE [W/m^3] vs t@','@PRICE [W/m^3] vs t@', &
           '@PRECI [W/m^3] vs t@', '@PRLHI [W/m^3] vs t@', &
           '@PRICI [W/m^3] vs t@','@AJEC [A/m^2] vs t@', &
           '@AJLH [A/m^2] vs t@', '@AJIC [A/m^2] vs t@', &
           '@NFAST [10^20/m^3] vs t@','@NIMP [10^20/m^3] vs t@', &
           '@BPOL [T] vs t@','@PSI [Wb] vs t@', &
           '@RMAJOR [m] vs t@','@RMINOR [m] vs t@', &
           '@VOLUME [m^3] sv t@','@KAPPAR vs t@', &
           '@DELTAR@','@GRHO1 vs t@', &
           '@GRHO2 vs t@','@AKDWE vs t@', &
           '@AKDWI vs t@','@PE [MPa] vs t@', &
           '@PI [MPa] vs t@','@VTOR [m/s] vs t@', &
           '@VPOL [m/s] vs t@','@S-ALPHA vs t@', &
           '@ER [V/m] vs t@', '@S vs t@', &
           '@ALPHA vs t@','@G vs t@', &
           '@IOTA vs t@'/

      DATA KV0 / &
           '@TE@','@TD@','@TT@','@TA@','@NE@', &
           '@ND@','@NT@','@NA@','@AJ@','@AJOH@', &
           '@AJNB@','@AJRF@','@AJBS@','@PIN@','@POH@', &
           '@PNB@','@PNF@','@PRFE@','@PRFD@','@PRFT@', &
           '@PRFA@','@PRL@','@PCX@','@PIE@','@PEXE@', &
           '@PEXI@','@QP@','@EZOH@','@BETA@','@BETAP@', &
           '@VLOOP@','@ETA@','@ZEFF@','@AKE@','@AKD@', &
           '@PRECE@','@PRLHE@','@PRICE@','@PRECI@','@PRLHI@', &
           '@PRICI@','@AJEC@','@AJLH@','@AJIC@','@NFAST@', &
           '@NIMP@','@BPOL@','@PSI@','@RMAJOR@','@RMINOR@', &
           '@VOLUME@','@KAPPAR@','@DELTAR@','@GRHO1@','@GRHO2@', &
           '@AKTBE@','@AKTBI@','@PE@','@PI@','@VTOR@', &
           '@VPOL@','@SALPHA@','@ER@','@S@','@ALPHA@', &
           '@G@','@IOTA@'/

      STRL=STR0(NMB)
      KVL =KV0 (NMB)

      RETURN
      END SUBROUTINE TRGRTD

!     **************************************************************

!           GRAPHIC 3D : PARAMETER STRING FOR NP

!     **************************************************************

      SUBROUTINE GETKRT(NP,KRT)

      IMPLICIT NONE
      INTEGER,       INTENT(IN) :: NP
      CHARACTER(LEN=4), INTENT(OUT):: KRT
      CHARACTER(LEN=4),DIMENSION(67)::KRTDATA
      DATA KRTDATA / &
           'TE  ','TD  ','TT  ','TA  ','NE  ', &
           'ND  ','NT  ','NA  ','AJ  ','AJOH', &
           'AJNB','AJRF','AJBS','PIN ','POH ', &
           'PNB ','PNF ','PRFE','PRFD','PRFT', &
           'PRFA','PRL ','PCX ','PIE ','PEXE', &
           'PEXI','QP  ','EZOH','BETA','BETP', &
           'VLOP','ETA ','ZEFF','AKE ','AKD ', &
           'PREE','PRLE','PRIE','PREI','PRLI', &
           'PRII','AJEC','AJLH','AJIC','NFST', &
           'NIMP','BP  ','PSI ','RMJ ','RMN ', &
           'VOL ','LKAP','DLT ','GRH1','GRH2', &
           'AKTE','AKTI','PE  ','PI  ','VTOR', &
           'VPOL','SALF','ER  ','S   ','ALFA', &
           'G   ','IOTA'/

      KRT=KRTDATA(NP)
      RETURN
      END SUBROUTINE GETKRT

!     **************************************************************

!           PARAMETER NAME STRING FOR GT

!     **************************************************************

      SUBROUTINE GETKGT(NP,KGT)

      USE TRCOMM,ONLY:NCTM
      IMPLICIT NONE
      INTEGER,       INTENT(IN) :: NP
      CHARACTER(LEN=4), INTENT(OUT):: KGT
      CHARACTER(LEN=4),DIMENSION(NCTM)::KGTDATA
      DATA KGTDATA / &
           'NE  ','ND  ','NT  ','NA  ','<NE>', &
           '<ND>','<NT>','<NA>','TE  ','TD  ', &
           'TD  ','TA  ','<TE>','<TD>','<TT>', &
           '<TA>','WE  ','WD  ','WT  ','WA  ', &
           'NB  ','NF  ','<NB>','<NF>','TB  ', &
           'TF  ','<TB>','<TF>','WB  ','WF  ', &
           'WBUL','WTAL','WTOT','IP  ','IOH ', &
           'INB ','IRF ','IBS ','PINT','POHT', &
           'PNBT','PRFE','PRFD','PRFT','PRFA', &
           'PNF ','PBIN','PBCE','PBCD','PBCT', &
           'PBCA','PFIN','PFCE','PFCD','PFCT', &
           'PFCA','POUT','PCXT','PIET','PRDT', &
           'PLE ','PLD ','PLT ','PLA ','SINT', &
           'SIET','SNBT','SNFT','SOUT','SLE ', &
           'SLD ','SLT ','SLA ','VLOP','ALI ', &
           'RQ1 ','Q0  ','WDOT','TE1 ','TE2 ', &
           'TE89','BTP0','BTPA','BT0 ','BTA ', &
           'ZEFF','QF  ','RIP ','PRFE','PRFI', &
           'PEE ','PECI','PLHE','PLHI','PICE', &
           'PICI','RR  ','RA  ','BB  ','RKAP', &
           'ITOT','TE98','H98Y','ANLE','ANLD', &
           'ANLT','ANLA','PRBT','PRCT','PRLT'/

      KGT=KGTDATA(NP)
      RETURN
      END SUBROUTINE GETKGT

!     **************************************************************

!           GRAPHIC 3D : CHECK A NUMBER CORRESPONDING A PARAMETER

!     **************************************************************

      SUBROUTINE CHNBPR(KK,NMB,IERR)

      USE TRCOMM,ONLY: NCRTM
      IMPLICIT NONE
      INTEGER,      INTENT(OUT):: NMB,IERR
      CHARACTER(LEN=4),INTENT(IN) :: KK
      INTEGER      :: NP
      CHARACTER(LEN=4):: KRT

      IERR=0
      DO NP=1,NCRTM
         CALL GETKRT(NP,KRT)
         IF (KK.EQ.KRT) THEN
            NMB=NP
            GOTO 1000
         ENDIF
      ENDDO
      IERR=1

 1000 RETURN
      END SUBROUTINE CHNBPR

!     **************************************************************

!           GRAPHIC 3D : VIEW 3D GRAPHICS LIST

!     **************************************************************

      SUBROUTINE VIEWRTLIST(NMAX)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: NMAX
      INTEGER :: I,L
      CHARACTER(LEN=4),DIMENSION(NMAX):: KRT

      DO I=1,NMAX
         CALL GETKRT(I,KRT(I))
      ENDDO
      WRITE(6,700)
      DO L=1,NMAX,10
         WRITE(6,710)  L,MIN(L+9,NMAX),(KRT(I),I= L,MIN(L+9,NMAX))
      ENDDO

 700  FORMAT(' ','# RT VARIABLE LIST')
 710  FORMAT(' ',I2,'-',I2,': ',10(A4:', '))

      RETURN
      END SUBROUTINE VIEWRTLIST

!     **************************************************************

!           GT VARIABLE LIST

!     **************************************************************

      SUBROUTINE VIEWGTLIST(NMAX)

      IMPLICIT NONE
      INTEGER,INTENT(IN):: NMAX
      INTEGER :: I,L
      CHARACTER(LEN=4),DIMENSION(NMAX):: KGT

      DO I=1,NMAX
         CALL GETKGT(I,KGT(I))
      ENDDO
      WRITE(6,700)
      DO L=1,NMAX,10
         WRITE(6,710)  L,MIN(L+9,NMAX),(KGT(I),I= L,MIN(L+9,NMAX))
      ENDDO

 700  FORMAT(' ','# GT VARIABLE LIST')
 710  FORMAT(' ',I2,'-',I2,': ',10(A4:', '))

      RETURN
      END SUBROUTINE VIEWGTLIST
