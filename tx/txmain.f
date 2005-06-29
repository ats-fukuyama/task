C     $Id$
C
C        ***********************************************************
C
C                        TOKAMAK TRANSPORT SIMULATION
C            INCLUDING RADIAL ELECTRIC FIELD AND PLASMA ROTATION
C
C             DEVELOPED BY Y. FUJI AND A. FUKUYAMA
C
C                     DEPARTMENT OF NUCLEAR ENGINEERING
C                      GRADUATE SCHOOL OF ENGINEERING
C                              KYOTO UNIVERSITY
C                              KYOTO 606-8501
C
C        ***********************************************************
C
C        X1  = Er                                       integer
C        X2  = Etheta                                   integer
C        X3  = Ephi                                half integer
C        X4  = Btheta                                   integer
C        X5  = Bphi                                half integer
C        X6  = Ne             X11  = Ni            half integer
C        X7  = Ne * Uer       X12  = Ni * Uir           integer
C        X8  = Ne * UeTheta   X13  = Ni * UiTheta       integer
C        X9  = Ne * UePhi     X14  = Ni * UiPhi    half integer
C        X10 = Ne * Te        X15  = Ni * Ti       half integer
C        X16 = Nb                                  half integer
C        X17 = Nb * UbTheta                             integer
C        X18 = Nb * UbPhi                          half integer
C        X19 = Slow n01                            half integer
C        X20 = Fast n02                            half integer
C
C        G16 = Total n0                            half integer
C        G20 = Q                                        integer
C        G21 = Jphi                                     integer
C        G22 = JePhi     G23 = JiPhi                    integer
C        G24 = JbPhi                                    integer
C        G25 = UePerp    G26 = UePara                   integer
C        G27 = UiPerp    G28 = UiPara                   integer
C        G29 = D    G30 = G1h2    G31 = MuI        half integer
C        G32 = S    G33 = Alpha   G34 = rKappa     half integer
C        G35 = Slow N0   G36 = Fast N0             half integer
C
      INCLUDE 'txcomm.inc'
      CHARACTER KPNAME*80
C
      CALL GSOPEN
      OPEN(7,STATUS='SCRATCH',FORM='FORMATTED')
C
C     ***** Version ID *****
C     SLID is used to identify data file.
      SLID = 'tx300.0'
      WRITE(6,*) '######## /TASK/TX V3.00.00 04/11/03 ########'
C
      CALL TXINIT
      KPNAME='txparm'
      CALL TXPARF(KPNAME)
C
      CALL TXMENU
C
      CLOSE(7)
      CALL GSCLOS
      STOP
      END
