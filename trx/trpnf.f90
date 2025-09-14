! trpnf.f90

MODULE trpnf

  PRIVATE
  PUBLIC tr_prep_pnf
  PUBLIC tr_pnf
  PRIVATE tr_nf_dt
  PRIVATE tr_nf_dd1
  PRIVATE tr_nf_dd2
  PRIVATE tr_nf_dd3
  PRIVATE tr_nf_dhe31
  PRIVATE tr_nf_dhe32
  PRIVATE tr_nf_tt
  PRIVATE tr_nf_the31
  PRIVATE tr_nf_the32
  PRIVATE tr_nf_the33
  PRIVATE tr_nf_the34
  PRIVATE tr_nf_the35
  PRIVATE tr_nf_the36

CONTAINS

  SUBROUTINE tr_prep_pnf

    USE trcomm
    USE libnf
    IMPLICIT NONE
    INTEGER:: nnf,id_nf

    ! --- initialize fusion cross section and reaction rate ---
    
    DO nnf=1,nnfmax
       id_nf=id_nf_nnf(nnf)
       ns1_nnf(nnf)=ns1_idnf(id_nf)
       ns2_nnf(nnf)=ns2_idnf(id_nf)
       nsp_nnf(nnf)=nsp_idnf(id_nf)
       wgt_nnf(nnf)=wgt_idnf(id_nf)
       eng_nnf(nnf)=eng_idnf(id_nf)
       enn_nnf(nnf)=enn_idnf(id_nf)
    END DO
    
  END SUBROUTINE tr_prep_pnf

  ! *** calculate fusion power ***

  SUBROUTINE tr_pnf

    USE trcomm
    USE libnf
    USE trlib
    IMPLICIT NONE
    REAL(rkind):: ANE,TE,P1,VC3,VCR,WF,VF,TAUS,HYF
    REAL(rkind):: PN1,PN2,PT1,RATE_NF,SNF
    REAL(rkind):: wgt,eng,enn
    INTEGER:: nnf,nr,id_nf,ns1,ns2,nsp,ns

    SNF_NSNNFNR(1:NSMAX,1:NNFMAX,1:NRMAX)=0.D0   ! particle source
    PNFCL_NSNNFNR(1:NSMAX,1:NNFMAX,1:NRMAX)=0.D0 ! collisional transfer in
    SNFNN_NNFNR(1:NNFMAX,1:NRMAX)=0.D0  ! neutron number
    PNFNN_NNFNR(1:NNFMAX,1:NRMAX)=0.D0  ! neutron power
    
    DO nnf=1,nnfmax
       id_nf=id_nf_nnf(nnf)
       ns1=ns1_nnf(nnf)
       ns2=ns2_nnf(nnf)
       wgt=wgt_nnf(nnf)
       nsp=nsp_nnf(nnf)
       eng=eng_nnf(nnf)
       enn=enn_nnf(nnf)
       DO NR=1,NRMAX
          PN1=RN(NR,ns1)
          PN2=RN(NR,ns2)
          PT1=RT(NR,ns1)
          RATE_NF=sigmav_nf(id_nf,PT1)
          SNF=wgt*PN1*PN2*1.D20*RATE_NF
          SNF_NSNNFNR(ns1,nnf,nr)=SNF_NSNNFNR(ns1,nnf,nr)-SNF
          SNF_NSNNFNR(ns2,nnf,nr)=SNF_NSNNFNR(ns2,nnf,nr)-SNF
          SNF_NSNNFNR(nsp,nnf,nr)=SNF_NSNNFNR(nsp,nnf,nr)+SNF
          PNF_NSNNFNR(nsp,nnf,nr)=PNF_NSNNFNR(nsp,nnf,nr)+eng*SNF
          IF(enn.GT.0.D0) THEN
             SNFNN_NNFNR(nnf,nr)=SNFNN_NNFNR(nnf,nr)+SNF
             PNFNN_NNFNR(nnf,nr)=PNFNN_NNFNR(nnf,nr)+enn*SNF
          END IF
       END DO
    END DO

    DO NR=1,NRMAX
       ANE= RN(NR,NS_e)
       TE = RT(NR,NS_e)
       P1   = 3.D0*SQRT(0.5D0*PI)*AME/ANE *(ABS(TE)*RKEV/AME)**1.5D0
       VC3=0.D0
       DO NS=1,NSMAX
          IF(PZ(NS).GT.0.D0) &    ! sum over ions
               VC3=VC3+P1*RN(NR,NS)*PZ(NS)**2/(PA(NS)*AMP)
       END DO
       VCR  = VC3**(1.D0/3.D0)
       DO nnf=1,nnfmax
          nsp=nsp_nnf(nnf)
          WF = RW(NR,NNBMAX+NNF)
          VF =SQRT(2.D0*eng_nnf(nnf)*RKEV/(PA(ns)*AMP))
          HYF=HY(VF/VCR)
          TAUS = 0.2D0*PA(ns)*ABS(TE)**1.5D0 &
               /(PZ(ns)**2*ANE*COULOG(1,ns,ANE,TE))
          TAUF(NNF,NR)= 0.5D0*TAUS*(1.D0-HYF)
       END DO
    END DO
          
    ! --- following variables are used in trcalc at every step ---
    
    DO NR=1,NRMAX
       DO NS=1,NSMAX
          SNF_NSNR(NS,NR)=SUM(SNF_NSNNFNR(NS,1:NNFMAX,NR))
          PNFCL_NSNR(NS,NR)=SUM(PNFCL_NSNNFNR(NS,1:NNFMAX,NR))
       END DO
    END DO

    RETURN
  END SUBROUTINE tr_pnf

END MODULE trpnf
