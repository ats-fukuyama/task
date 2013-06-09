C     $Id$
C
C     ******* SET KANOD ARRAY *******
C
      SUBROUTINE DEFBND(IERR)
C
      INCLUDE 'wfcomm.inc'
C
C     Set number of variables to describe a node
C
      IERR=0
      IL=1
      DO NSD=1,NSDMAX
         KA=KASID(NSD)
         IF(KA.EQ.0) THEN
            K=1
         ELSEIF(KA.EQ.1) THEN
            K=0
         ELSEIF(KA.EQ.4) THEN
            K=1
         ELSEIF(KA.LT.0) THEN
            K=0
         ELSE
            WRITE(6,*) 'XX DEFBND: UNDEFINED KANOD=',KA
            IERR=2
         ENDIF
         INLEN(NSD)=K
         IMLEN(NSD)=IL
         IL=IL+K
C         WRITE(6,'(A,5I5)') 'NSD,KA,K,IL=',NSD,KA,K,IL
      ENDDO
C
C     WG boundary block requires a node
C
      NSD=NSDMAX
      NMDMAX=1
      DO NB=1,NBMAX
         IF(KABDY(NB).GE.8) THEN
            K=NMBDY(NB)
            NMDMAX=MAX(NMDMAX,K)
            NSD=NSD+1
            NDBDY(NB)=NSD
            INLEN(NSD)=K
            IMLEN(NSD)=IL
            IL=IL+K
         ENDIF
      ENDDO
C
      NNBMAX=NSD
      MLEN=IL-1
      IMLEN(NSD+1)=IL
      IF(MLEN.GT.MLENM) GOTO 9000
C
C     Set NBELM to identify an element with WG boundary
C         NDELM for additional node for the element
C
      DO NE=1,NEMAX
         NBELM(NE)=0
      ENDDO
      DO NE=NEMAX,1,-1
         DO IN=1,4
            NN=NDELM(IN,NE)
            KA=KANOD(NN)
            IF(KA.LT.0) THEN
              NB=-KA
              IF(KABDY(NB).GT.8) THEN
                 NSDELM(7,NE)=NDBDY(NB)
                 IF(NBELM(NE).EQ.0) THEN
                    NBELM(NE)=NB
                 ELSE
                    IF(NBELM(NE).NE.NB) THEN
                       WRITE(6,*) 'XX Element faces two boundary block'
                       WRITE(6,*) '   NE=',NE
                       IERR=3
                    ENDIF
                 ENDIF
              ENDIF
           ENDIF
        ENDDO
      ENDDO
C      DO NB=1,NBMAX
C         NN=NDBDY(NB)
C         WRITE(6,'(A,6I5)') 'NB,KABDY,NMBDY,NDBDY,INLEN,IMLEN=',
C     &        NB,KABDY(NB),NMBDY(NB),NN,INLEN(NN),IMLEN(NN)
C      ENDDO
C      WRITE(6,*) 'MLEN,MLENM =',MLEN,MLENM
      RETURN
C
 9000 WRITE(6,*) 'XX DEFBND ERROR : MLEN,MLENM =',MLEN,MLENM
      IERR=900
      RETURN
      END
C
C     ******* DEF MBND *******
C
      SUBROUTINE DEFMBN(IERR)
C
      INCLUDE 'wfcomm.inc'
C
C     FIND LAST APPEAREANCE OF EACH NODE
C
      DO NE=1,NEMAX
         IF(NBELM(NE).EQ.0) THEN
            ISDMAX=6
         ELSE
            ISDMAX=7
         ENDIF
         DO ISD=1,ISDMAX
            ISDELM(ISD,NE)=ABS(NSDELM(ISD,NE))
         ENDDO
      ENDDO
      DO NSD=1,NNBMAX
         NFLG(NSD)=0
      ENDDO
      DO M=1,MLEN
         LHED(M)=0
      ENDDO
C
      DO NE=NEMAX,1,-1
         IF(NBELM(NE).EQ.0) THEN
            ISDMAX=6
         ELSE
            ISDMAX=7
         ENDIF
         DO ISD=1,ISDMAX
            NSD=ABS(ISDELM(ISD,NE))
            IF(NFLG(NSD).EQ.0) THEN
               NFLG(NSD)=1
               ISDELM(ISD,NE)=-ISDELM(ISD,NE)
            ENDIF
         ENDDO
      ENDDO
C
C     PREFRONT
C
      MBND=0
      NELL=0
      LCOL=0
C
 1000 NELL=NELL+1
C
      IF(NBELM(NELL).EQ.0) THEN
         ISDMAX=6
      ELSE
         ISDMAX=7
      ENDIF
C      WRITE(6,*) 'NE,NB=',NELL,NBELM(NELL)
C
      KC=0
      DO J=1,ISDMAX
         NSD=ISDELM(J,NELL)
         M=ABS(NSD)
         K=IMLEN(M)
         DO L=1,INLEN(M)
            KC=KC+1
            II=K+L-1
            IF(NSD.LT.0) II=-II
            NODEK(KC)=II
         ENDDO
      ENDDO
C
C     SET UP HEADING VECTORS
C
      DO LK=1,KC
         NODE1=NODEK(LK)
         DO L=1,LCOL
            IF(ABS(NODE1).EQ.ABS(LHED(L))) THEN
               LHED(L)=NODE1
               GOTO 60
            ENDIF
         ENDDO
         LCOL=LCOL+1
         IF(LCOL.GT.MBND) MBND=LCOL
         IF(MBND.GT.MLENM) GOTO 9200
         LHED(LCOL)=NODE1
   60    CONTINUE
      ENDDO
C
      L=1
   65 IF(LHED(L).LT.0) THEN
         DO LL=L+1,LCOL
            LHED(LL-1)=LHED(LL)
         ENDDO
         LCOL=LCOL-1
      ELSE
         L=L+1
      ENDIF
      IF(L.LE.LCOL) GOTO 65
      IF(NELL.LT.NEMAX) GOTO 1000
      IF(LCOL.GT.0) GOTO 9100
      IF(MBND.GT.MBNDM) GOTO 9000
      IERR=0
      RETURN
C
 9000 IERR=9000
      WRITE(6,900) MBND,MBNDM
  900 FORMAT(' ','XX DEFMBN ERROR: MBND.GT.MBNDM: MBND,MBNDM=',2I8)
      RETURN
C
 9100 IERR=9100
      WRITE(6,*) 'XX DEFMBN ERROR: LCOL.GT.0: LCOL=',LCOL
      RETURN
C
 9200 IERR=9200
      WRITE(6,920) MBND,MLENM
  920 FORMAT(' ','XX DEFMBN ERROR : MBND.GT.MLENM : MBND,MLENM=',2I8)
      RETURN
      END
C
C     ******* FRONTAL ELIMINATION WITH LEAST DIAGONAL PIVOTING *******
C
      SUBROUTINE CVSOLV(IERR)
C
      INCLUDE 'wfcomm.inc'
      COMMON /WFBUF1/ CBUF(MBUFM),MHED(MBUFM)
      COMMON /WFBUF2/ CM(NMDM,NMDM,7,7),CV(NMDM,7)
C
      DATA ND1/23/
C
      MTRACK=MBUFM*(4+16)
      OPEN(ND1,FILE=KFNAMB,FORM='UNFORMATTED',
     &         ACCESS='DIRECT',RECL=MTRACK)
      IBUFF=0
C
C     CLEAR BUFFER
C
      DO MBUF=1,MBUFM
         CBUF(MBUF)=(0.D0,0.D0)
      ENDDO
      DO MBUF=1,MBUFM
         MHED(MBUF)=0
      ENDDO
      DO I=1,MLEN
         CRV(I)=0.D0
         CSV(I)=0.D0
      ENDDO
C
C     FIND LAST APPEAREANCE OF EACH NODE
C
      DO NE=1,NEMAX
         IF(NBELM(NE).EQ.0) THEN
            ISDMAX=6
         ELSE
            ISDMAX=7
         ENDIF
         DO ISD=1,ISDMAX
            ISDELM(ISD,NE)=ABS(NSDELM(ISD,NE))
         ENDDO
      ENDDO
      DO NSD=1,NNBMAX
         NFLG(NSD)=0
      ENDDO
C
      DO NE=NEMAX,1,-1
         IF(NBELM(NE).EQ.0) THEN
            ISDMAX=6
         ELSE
            ISDMAX=7
         ENDIF
         DO ISD=1,ISDMAX
            NSD=ABS(ISDELM(ISD,NE))
            IF(NFLG(NSD).EQ.0) THEN
               NFLG(NSD)=1
               ISDELM(ISD,NE)=-ISDELM(ISD,NE)
            ENDIF
         ENDDO
      ENDDO
C
C     PREFRONT
C
      NMAX=MBNDM
      NELL=0
      LCOL=0
      ILEN=0
      IPOS=0
      DO I=1,NMAX
      DO J=1,NMAX
         CEQ(J,I)=0.D0
      ENDDO
      ENDDO
C
 1000 NELL=NELL+1
      IF(MOD(NELL,10000).EQ.0) WRITE(6,'(A,I10)') '--- NE = ',NELL
C
      IF(NBELM(NELL).EQ.0) THEN
         ISDMAX=6
      ELSE
         ISDMAX=7
      ENDIF
      CALL CMCALC(NELL,CM,CV)
C      WRITE(6,*) 'NE,NB=',NELL,NBELM(NELL)
C
      KC=0
      DO J=1,ISDMAX
         NSD=ISDELM(J,NELL)
         M=ABS(NSD)
         K=IMLEN(M)
         DO L=1,INLEN(M)
            KC=KC+1
            II=K+L-1
            IF(NSD.LT.0) II=-II
            NODEK(KC)=II
         ENDDO
      ENDDO
C
C     SET UP HEADING VECTORS
C
      DO LK=1,KC
         NODE1=NODEK(LK)
         DO L=1,LCOL
            IF(ABS(NODE1).EQ.ABS(LHED(L))) THEN
               LDEST(LK)=L
               LHED(LDEST(LK))=NODE1
               GOTO 60
            ENDIF
         ENDDO
         LCOL=LCOL+1
         IF(LCOL.GT.NMAX) GOTO 9000
         LDEST(LK)=LCOL
         LHED(LCOL)=NODE1
   60    CONTINUE
      ENDDO
C
C     ASSEMBLY
C
      L=0
      DO J=1,ISDMAX
         DO JJ=1,INLEN(ABS(ISDELM(J,NELL)))
            L=L+1
            LL=LDEST(L)
            K=0
            DO I=1,ISDMAX
               DO II=1,INLEN(ABS(ISDELM(I,NELL)))
                  K=K+1
                  KK=LDEST(K)
                  CEQ(KK,LL)=CEQ(KK,LL)+CM(II,JJ,I,J)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      K=0
      DO I=1,ISDMAX
         DO II=1,INLEN(ABS(ISDELM(I,NELL)))
            K=K+1
            KK=LDEST(K)
            LCO=ABS(LHED(KK))
            CRV(LCO)=CRV(LCO)+CV(II,I)
         ENDDO
      ENDDO
C
C     FIND OUT WHICH MATRIX ELEMENTS ARE FULLY SUMMED
C
 2000 LC=0
      DO L=1,LCOL
         IF(LHED(L).LT.0) THEN
            LC=LC+1
            LPIV(LC)=L
         ENDIF
      ENDDO
C
      IF(LC.LE.0) GOTO 1000
C
C     SEARCH FOR ABSOLUTE PIVOT
C
      CPIVOT=0.D0
      DO L=1,LC
         LPIVC=LPIV(L)
         CPIVA=CEQ(LPIVC,LPIVC)
         IF(CDABS(CPIVA).GE.CDABS(CPIVOT)) THEN
            CPIVOT=CPIVA
            LPIVCO=LPIVC
         ENDIF
      ENDDO
C
C     NORMALIZE PIVOTAL ROW
C
      LCO=ABS(LHED(LPIVCO))
      IF(CDABS(CPIVOT).LT.1.D-15) THEN
         WRITE(6,601)
         WRITE(6,'(A,3I8)') 'NELL,LCO,LCOL=',NELL,LCO,LCOL
            IE=NELL
            WRITE(6,'(A,5I8)') '  NE,KN=',IE,
     &               KNELM(1,IE),KNELM(2,IE),KNELM(3,IE),KNELM(4,IE)
            WRITE(6,'(A,5I8)') '  KA,ND=',KAELM(IE),
     &               NDELM(1,IE),NDELM(2,IE),NDELM(3,IE),NDELM(4,IE)
            DO I=1,4
               IN=NDELM(I,IE)
               RND=SQRT(XND(IN)**2+YND(IN)**2)
               WRITE(6,'(A,I5,1P4E12.4,I3,I5)') 'IN,X,Y,Z,R,KA,IM =',
     &           IN,XND(IN),YND(IN),ZND(IN),RND,KANOD(IN),IMLEN(IN)
            ENDDO
         STOP
      ENDIF
      DO L=1,LCOL
         CQQ(L)=CEQ(LPIVCO,L)/CPIVOT
      ENDDO
      CRHS=CRV(LCO)/CPIVOT
      CRV(LCO)=CRHS
C
C     ELIMINATE THEN DELETE PIVOTAL ROW AND COLUMN
C
      DO K=1,LPIVCO-1
         CFAC=CEQ(K,LPIVCO)
         KRW=ABS(LHED(K))
         CRV(KRW)=CRV(KRW)-CFAC*CRHS
         DO L=1,LPIVCO-1
            CEQ(K,L)=CEQ(K,L)-CFAC*CQQ(L)
         ENDDO
         DO L=LPIVCO+1,LCOL
            CEQ(K,L-1)=CEQ(K,L)-CFAC*CQQ(L)
         ENDDO
      ENDDO
      DO  K=LPIVCO+1,LCOL
         CFAC=CEQ(K,LPIVCO)
         KRW=ABS(LHED(K))
         CRV(KRW)=CRV(KRW)-CFAC*CRHS
         DO L=1,LPIVCO-1
            CEQ(K-1,L)=CEQ(K,L)-CFAC*CQQ(L)
         ENDDO
         DO L=LPIVCO+1,LCOL
            CEQ(K-1,L-1)=CEQ(K,L)-CFAC*CQQ(L)
         ENDDO
      ENDDO
      DO L=LPIVCO+1,LCOL
         LHED(L-1)=LHED(L)
         CQQ (L-1)=CQQ (L)
      ENDDO
      DO K=1,LCOL
         CEQ(K,LCOL)=0.D0
      ENDDO
      DO L=1,LCOL-1
         CEQ(LCOL,L)=0.D0
      ENDDO
      LCOL=LCOL-1
C
C     WRITE PIVOTAL EQUATION ON DISC
C
      IF(IPOS+LCOL.GT.MBUFM) THEN
         IBUFF=IBUFF+1
         WRITE(ND1,REC=IBUFF) MHED,CBUF
         IPOS=0
      ENDIF
      ILEN=ILEN+1
      MLCO(ILEN)=LCO
      MCOL(ILEN)=LCOL
      MPOS(ILEN)=IPOS
      DO L=1,LCOL
         MHED(IPOS+L)=ABS(LHED(L))
         CBUF(IPOS+L)=CQQ(L)
      ENDDO
      IPOS=IPOS+LCOL
C
C     DETERMINE WHETHER TO ASSEMBLE OR BACKSUBSTITUTE
C
      IF(LCOL.GT.1) GOTO 2000
      LCO=ABS(LHED(1))
      CPIVOT=CEQ(1,1)
      IF(ABS(CPIVOT).LT.1.D-15) THEN
         WRITE(6,601)
         WRITE(6,'(A,2I8)') 'NELL,LCO=',NELL,LCO
         STOP
      ENDIF
      CSV(LCO)=CRV(LCO)/CPIVOT
      IF(IBUFF.NE.0) WRITE(6,*) '## BUFFER MAX = ',IBUFF
      WRITE(6,*) '## CVSOLV : MBUFMAX = ',IPOS
      WRITE(6,'(A,1PE12.4)') '   CPIVOT = ',ABS(CPIVOT)
C
C     BACK SUBSTITUTION
C
      DO ILEN=MLEN-1,1,-1
         IF(IPOS.EQ.0) THEN
            READ(ND1,REC=IBUFF) MHED,CBUF
            IBUFF=IBUFF-1
         ENDIF
         LCO =MLCO(ILEN)
         IPOS=MPOS(ILEN)
         CGASH=0.D0
         DO L=1,MCOL(ILEN)
            CGASH=CGASH-CBUF(IPOS+L)*CSV(MHED(IPOS+L))
         ENDDO
         IF(LCO.EQ.0) WRITE(6,*) 'XX LCO.EQ.0: IPOS,ILNE=',IPOS,ILEN
         CSV(LCO)=CRV(LCO)+CGASH
      ENDDO
      CLOSE(ND1)
      IERR=0
      RETURN
C
 9000 IERR=2
      WRITE(6,602) IERR
      WRITE(6,'(A,3I8)') 'XX LCOL,NMAX,NELL=',LCOL,NMAX,NELL
      RETURN
C
  601 FORMAT(' ','## FRONT WARNING : SINGULAR OR ILL CONDITIONED')
  602 FORMAT(' ','## FRONT ERROR : ERR =',I5/
     &       ' ','       LCOL EXCEEDS NMAX')
      END
