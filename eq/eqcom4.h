C
      PARAMETER (NTRM=200)
C
      COMMON /TREQL1/ PSIRHO(NTRM+2),UPSIRHO(4,NTRM+2)
      COMMON /TREQL2/ PSITR(NTRM),PSITRG(NTRM+1),PSITRX(NTRM+2)
      COMMON /TREQL2/ RHOTR(NTRM),RHOTRG(NTRM+1),RHOTRX(NTRM+2)
C
      COMMON /TREQG1/ UPPSI(4,NTRM+2),UJPSI(4,NTRM+2),UJPSI0(NTRM+2)
      COMMON /TREQG2/ UVTPSI(4,NTRM+2),UTPSI(4,NTRM+2)
C
      COMMON /TREQP1/ NTRMAX
