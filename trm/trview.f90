! trview.f90

MODULE trview

  PRIVATE
  PUBLIC tr_view

CONTAINS

!     ***********************************************************

!           VIEW INPUT PARAMETER

!     ***********************************************************

    SUBROUTINE tr_view(ID)

      USE trcomm

      IMPLICIT NONE
      INTEGER,INTENT(IN) :: ID
      INTEGER :: NS,NPSC


      WRITE(6,*) '** TRANSPORT **'
      WRITE(6,602) 'MDLEQB',MDLEQB,'MDLEQN',MDLEQN,'MDLEQT',MDLEQT,'MDLEQU',MDLEQU
      WRITE(6,602) 'MDLEQZ',MDLEQZ,'MDLEQ0',MDLEQ0,'MDLEQE',MDLEQE,'MDLEOI',MDLEOI
      WRITE(6,602) 'NSMAX ',NSMAX, 'NSZMAX',NSZMAX,'NSNMAX',NSNMAX
      WRITE(6,601) 'RR    ',RR,    'RA    ',RA,    'RKAP  ',RKAP,  'RDLT  ',RDLT
      WRITE(6,601) 'RIPS  ',RIPS,  'RIPE  ',RIPE,  'BB    ',BB

      WRITE(6,611)
  611 FORMAT(' ','NS',2X,'PA           PZ      PN(E20)  PNS(E20) ','PT(KEV)  PTS(KEV) PELPAT')
      DO NS=1,NSMAX
         WRITE(6,612) NS,PA(NS),PZ(NS),PN(NS),PNS(NS),PT(NS),PTS(NS),PELPAT(NS)
  612    FORMAT(' ',I2,1PD12.4,0P,F8.3,5F9.4)
      ENDDO

      WRITE(6,601) 'PROFN1',PROFN1,'PROFT1',PROFT1,'PROFU1',PROFU1,'PROFJ1',PROFJ1
      WRITE(6,601) 'PROFN2',PROFN2,'PROFT2',PROFT2,'PROFU2',PROFU2,'PROFJ2',PROFJ2
      WRITE(6,601) 'ALP(1)',ALP(1),'ALP(2)',ALP(2),'ALP(3)',ALP(3),'PBSCD ',PBSCD
      WRITE(6,602) 'MDLKAI',MDLKAI,'MDLETA',MDLETA,'MDLAD ',MDLAD, 'MDLAVK',MDLAVK
      WRITE(6,602) 'MDLJBS',MDLJBS,'MDLKNC',MDLKNC,'MDLTPF',MDLTPF,'MDNCLS',MDNCLS
      WRITE(6,604) 'MDLUF ',MDLUF, 'KUFDEV',KUFDEV,'KUFDCG',KUFDCG,'MDNI  ',MDNI
      WRITE(6,605) 'MDLJQ ',MDLJQ, 'MDLFLX',MDLFLX,'MDTC  ',MDTC,  'RHOA  ',RHOA
      WRITE(6,602) 'MDLWLD',MDLWLD,'MDLER ',MDLER, 'MODELG',MODELG,'NTEQIT',NTEQIT
      WRITE(6,603) 'MDCD05',MDCD05,'CK0   ',CK0,   'CK1   ',CK1
      WRITE(6,603) 'MDEDGE',MDEDGE,'CSPRS ',CSPRS, 'CNN   ',CNN
      WRITE(6,601) 'CNP   ',CNP,   'CNH   ',CNH,   'CDP   ',CDP,   'CDH   ',CDH
      WRITE(6,601) 'AD0   ',AD0,   'CHP   ',CHP,   'CWEB  ',CWEB,  'CALF  ',CALF
      WRITE(6,631)     'model_prof  ',model_prof
      WRITE(6,'(A,A)') 'knam_prof   ',knam_prof
      IF((MDLKAI.GE.1.AND.MDLKAI.LT.10).OR.ID.EQ.1) &
         WRITE(6,601) 'CKALFA',CKALFA,'CKBETA',CKBETA,'CKGUMA',CKGUMA

      IF((MDLKAI.GE.10.AND.MDLKAI.LT.20).OR.ID.EQ.1)  &
         WRITE(6,613) CDW(1),CDW(2),CDW(3),CDW(4),CDW(5),CDW(6),CDW(7),CDW(8)
  613 FORMAT(' ','    AKDW(E) =  ',0PF6.3,' DEDW + ',F6.3,' DIDW'/ &
             ' ','    AKDW(D) =  ',0PF6.3,' DEDW + ',F6.3,' DIDW'/ &
             ' ','    AKDW(T) =  ',0PF6.3,' DEDW + ',F6.3,' DIDW'/ &
             ' ','    AKDW(A) =  ',0PF6.3,' DEDW + ',F6.3,' DIDW')

      WRITE(6,601) 'DT    ',DT,    'EPSLTR',EPSLTR,'TSST  ',TSST,  'TPRST ',TPRST
      WRITE(6,602) 'LMAXTR',LMAXTR,'NRMAX ',NRMAX, 'NTMAX ',NTMAX, 'NTSTEP',NTSTEP
      WRITE(6,602) 'NGRSTP',NGRSTP,'NGTSTP',NGTSTP,'NGPST ',NGPST, 'IZERO ',IZERO
      WRITE(6,602) 'MDLST ',MDLST, 'MDLCD ',MDLCD, 'MDLNF ',MDLNF

      IF(MDLIMP.GT.0) THEN
         WRITE(6,602) 'MDLIMP',MDLIMP
         WRITE(6,601) 'PNC   ',PNC,   'PNFE  ',PNFE
      END IF
!         WRITE(6,601) 'PNNU  ',PNNU,  'PNNUS ',PNNUS

      IF((PNBTOT.GT.0.D0).OR.(ID.EQ.1)) THEN
         WRITE(6,601) 'PNBTOT',PNBTOT,'PNBR0 ',PNBR0,'PNBRW ',PNBRW,'PNBENG',PNBENG
         WRITE(6,603) 'MDLNB ',MDLNB, 'PNBRTG',PNBRTG,'PNBCD ',PNBCD,'PNBVY ',PNBVY
         WRITE(6,623) 'NRNBMAX ',NRNBMAX,'PNBVW   ',PNBVW
      ENDIF

      IF((PECTOT.GT.0.D0).OR.(ID.EQ.1)) THEN
         WRITE(6,601) 'PECTOT',PECTOT, 'PECR0 ',PECR0,'PECRW ',PECRW,'PECTOE',PECTOE
         WRITE(6,603) 'MDLEC ',MDLEC,  'PECNPR',PECNPR,'PECCD ',PECCD
      ENDIF

      IF((PLHTOT.GT.0.D0).OR.(ID.EQ.1)) THEN
         WRITE(6,601) 'PLHTOT',PLHTOT, 'PLHR0 ',PLHR0,'PLHRW ',PLHRW, 'PLHTOE',PLHTOE
         WRITE(6,603) 'MDLLH ',MDLLH,  'PLHNPR',PLHNPR,'PLHCD ',PLHCD
      ENDIF

      IF((PICTOT.GT.0.D0).OR.(ID.EQ.1)) THEN
         WRITE(6,601) 'PICTOT',PICTOT, 'PICR0 ',PICR0, 'PICRW ',PICRW,'PICTOE',PICTOE
         WRITE(6,603) 'MDLIC ',MDLIC,  'PICNPR',PICNPR,'PICCD ',PICCD
      ENDIF

      IF((PELTOT.GT.0.D0).OR.(ID.EQ.1)) THEN
         WRITE(6,601) 'PELTOT',PELTOT,'PELR0 ',PELR0,'PELRW ',PELRW
         WRITE(6,603) 'MDLPEL',MDLPEL,'PELRAD',PELRAD, &
                      'PELVEL',PELVEL,'PELTIM',PELTIM
         WRITE(6,'(2(A24,ES12.4))') &
              'pellet_time_start'      ,pellet_time_start, &
              'pellet_time_interval'   ,pellet_time_interval
         WRITE(6,'(2(A24,I8,4X))') &
              'number_of_pellet_repeat',number_of_pellet_repeat
      ENDIF
      IF(MDLPSC.GT.0) THEN
         WRITE(6,622) 'MDLPSC  ',MDLPSC,'NPSCMAX ',NPSCMAX
         DO NPSC=1,NPSCMAX
            WRITE(6,603) 'NSPSC ',NSPSC(NPSC) ,'PSCTOT',PSCTOT(NPSC), &
                         'PSCR0 ',PSCR0(NPSC) ,'PSCRW ',PSCRW(NPSC)
         END DO
      END IF

      IF((MDLPR.GE.1).OR.(ID.EQ.1)) THEN
         WRITE(6,623) 'MDLPR   ',MDLPR,   'SYNCABS ',SYNCABS, &
                      'SYNCSELF',SYNCSELF
      ENDIF

      WRITE(6,'(A,A)') 'KNAMEQ =',knameq
      WRITE(6,'(A,A)') 'KNAMEQ2=',knameq2
      WRITE(6,'(A,A)') 'KNAMTR =',knamtr
      WRITE(6,'(A,A)') 'KFNLOG =',kfnlog
      WRITE(6,'(A,A)') 'KFNTXT =',kfntxt
      WRITE(6,'(A,A)') 'KFNCVS =',kfncvs

      WRITE(6,601) 'CNB   ',CNB
      RETURN

  601 FORMAT(' ',A6,'=',1PE11.3 :2X,A6,'=',1PE11.3: &
              2X,A6,'=',1PE11.3 :2X,A6,'=',1PE11.3)
  602 FORMAT(' ',A6,'=',I7,4X   :2X,A6,'=',I7,4X  : &
              2X,A6,'=',I7,4X   :2X,A6,'=',I7)
  603 FORMAT(' ',A6,'=',I7,4X   :2X,A6,'=',1PE11.3: &
              2X,A6,'=',1PE11.3 :2X,A6,'=',1PE11.3)
  604 FORMAT(' ',A6,'=',I7,4X   :2X,A6,'=',1X,A6,4X: &
              2X,A6,'=',1X,A6,4X:2X,A6,'=',I7)
  605 FORMAT(' ',A6,'=',I7,4X   :2X,A6,'=',I7,4X  : &
              2X,A6,'=',I7,4X   :2X,A6,'=',1PE11.3)
  622 FORMAT(' ',A8,'=',I5,4X   :2X,A8,'=',I5,4X  : &
              2X,A8,'=',I5,4X   :2X,A8,'=',I5)
  623 FORMAT(' ',A8,'=',I7,4X   :2X,A8,'=',1PE11.3: &
              2X,A8,'=',1PE11.3)
631   FORMAT(' ',A12,'=',ES11.3)
    END SUBROUTINE tr_view
END MODULE trview
