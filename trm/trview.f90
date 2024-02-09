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
      INTEGER :: NS,NNB,NEC,NLH,NIC,NPEL,NPSC,idx

      WRITE(6,*) '** TRANSPORT **'
      WRITE(6,602) 'MDLEQB',MDLEQB,'MDLEQN',MDLEQN, &
                   'MDLEQT',MDLEQT,'MDLEQU',MDLEQU
      WRITE(6,602) 'MDLEQZ',MDLEQZ,'MDLEQ0',MDLEQ0, &
                   'MDLEQE',MDLEQE,'MDLEOI',MDLEOI
      WRITE(6,602) 'NSMAX ',NSMAX, 'NSZMAX',NSZMAX, &
                   'NSNMAX',NSNMAX
      WRITE(6,601) 'RR    ',RR,    'RA    ',RA,     &
                   'RKAP  ',RKAP,  'RDLT  ',RDLT
      WRITE(6,601) 'RIPS  ',RIPS,  'RIPE  ',RIPE,   &
                   'BB    ',BB

      WRITE(6,611)
611   FORMAT(' ','NS',2X, &
             'PA           PZ      PN(E20)  PNS(E20) ','PT(KEV)  PTS(KEV)')
      DO NS=1,NSMAX
         WRITE(6,612) NS,PA(NS),PZ(NS),PN(NS),PNS(NS),PT(NS),PTS(NS)
      ENDDO

      WRITE(6,601) 'PROFN1',PROFN1,'PROFT1',PROFT1,'PROFU1',PROFU1,'PROFJ1',PROFJ1
      WRITE(6,601) 'PROFN2',PROFN2,'PROFT2',PROFT2,'PROFU2',PROFU2,'PROFJ2',PROFJ2
      WRITE(6,601) 'ALP(1)',ALP(1),'ALP(2)',ALP(2),'ALP(3)',ALP(3),'PBSCD ',PBSCD
      WRITE(6,602) 'MDLKAI',MDLKAI,'MDLETA',MDLETA,'MDLAD ',MDLAD, 'MDLAVK',MDLAVK
      WRITE(6,602) 'MDLJBS',MDLJBS,'MDLKNC',MDLKNC,'MDLTPF',MDLTPF,'MDLNCL',MDLNCL
      WRITE(6,604) 'MDLUF ',MDLUF, 'KUFDEV',KUFDEV,'KUFDCG',KUFDCG,'MDLNI ',MDLNI
      WRITE(6,605) 'MDLJQ ',MDLJQ, 'MDLFLX',MDLFLX,'MDLTC ',MDLTC, 'RHOA  ',RHOA
      WRITE(6,602) 'MDLWLD',MDLWLD,'MDLER ',MDLER, 'MODELG',MODELG,'NTEQIT',NTEQIT
      WRITE(6,603) 'MDLCD05',MDLCD05,'CK0   ',CK0,   'CK1   ',CK1
      WRITE(6,603) 'MDLEDGE',MDLEDGE,'CSPRS ',CSPRS, 'CNN   ',CNN
      WRITE(6,601) 'CNP   ',CNP,   'CNH   ',CNH,   'CDP   ',CDP,   'CDH   ',CDH
      WRITE(6,601) 'AD0   ',AD0,   'CHP   ',CHP,   'CWEB  ',CWEB,  'CALF  ',CALF
      WRITE(6,630)     'model_prof  ',model_prof
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

         DO NNB=1,NNBMAX
            IF(PNBIN(NNB).GT.0.D0) THEN
               WRITE(6,632) NNB, &
                    'MDLNB ',MDLNB(nnb), &
                    'PNBIN ',PNBIN(NNB), &
                    'PNBR0 ',PNBR0(NNB), &
                    'PNBRW ',PNBRW(nnb)
            END IF
         END DO
         DO NNB=1,NNBMAX
            IF(PNBIN(NNB).GT.0.D0) THEN
               WRITE(6,633) NNB, &
                    'PNBVY ',PNBVY(nnb), &
                    'PNBVW ',PNBVW(nnb), &
                    'PNBENG',PNBENG(nnb), &
                    'PNBRTG',PNBRTG(nnb)
            END IF
         END DO
         DO NNB=1,NNBMAX
            IF(PNBIN(NNB).GT.0.D0) THEN
               WRITE(6,634) NNB, &
                    'NRNBMAX ',NRNBMAX(nnb), &
                    'NSPNB   ',NSPNB(nnb), &
                    'PNBCD   ',PNBCD(nnb)
            END IF
         END DO

         DO NEC=1,NECMAX
            IF(PECIN(NEC).GT.0.D0) THEN
               WRITE(6,632) NEC, &
                    'MDLEC ',MDLEC(nec), &
                    'PECIN ',PECIN(nec), &
                    'PECR0 ',PECR0(nec), &
                    'PECRW ',PECRW(nec)
            END IF
         END DO
         DO NEC=1,NECMAX
            IF(PECIN(NEC).GT.0.D0) THEN
               WRITE(6,633) NEC, &
                    'PECTOE',PECTOE(nec), &
                    'PECNPR',PECNPR(nec), &
                    'PECCD ',PECCD(nec)
            END IF
         END DO

         DO NLH=1,NLHMAX
            IF(PLHIN(NLH).GT.0.D0) THEN
               WRITE(6,632) NLH, &
                    'MDLLH ',MDLLH(nlh), &
                    'PLHIN ',PLHIN(nlh), &
                    'PLHR0 ',PLHR0(nlh), &
                    'PLHRW ',PLHRW(nlh)
            END IF
         END DO
         DO NLH=1,NLHMAX
            IF(PLHIN(NLH).GT.0.D0) THEN
               WRITE(6,633) NLH, &
                    'PLHCD ',PLHCD(NLH), &
                    'PLHTOE',PLHTOE(NLH), &
                    'PLHNPR',PLHNPR(NLH)
            END IF
         END DO

         DO NIC=1,NICMAX
            IF(PICIN(NIC).GT.0.D0) THEN
               WRITE(6,632) NIC, &
                    'MDLIC ',MDLIC(NIC), &
                    'PICIN ',PICIN(NIC), &
                    'PICR0 ',PICR0(NIC), &
                    'PICRW ',PICRW(NIC)
            END IF
         END DO
         DO NIC=1,NICMAX
            IF(PICIN(NIC).GT.0.D0) THEN
               WRITE(6,633) NIC, &
                    'PICCD ',PICCD(NIC), &
                    'PICTOE',PICTOE(NIC), &
                    'PICNPR',PICNPR(NIC)
            ENDIF
         END DO

      IF(PELIN.GT.0.D0) THEN
         WRITE(6,632) 1, &
              'MDLPEL',MDLPEL, &
              'PELIN ',PELIN, &
              'PELR0 ',PELR0, &
              'PELRW ',PELRW
         WRITE(6,633) 1, &
                 'PELRAD',PELRAD, &
                 'PELVEL',PELVEL, &
                 'PELTIM',PELTIM
         WRITE(6,'(A24,4ES12.4)') &
              'pellet_time_start', &
               pellet_time_start
         WRITE(6,'(A24,4ES12.4)') &
              'pellet_time_interval', &
               pellet_time_interval
         WRITE(6,'(A24,4I12)') &
              'number_of_pellet_repeat', &
               number_of_pellet_repeat
         DO NS=1,NSMAX
            WRITE(6,'(A7,I2,A3,ES12.4)') 'PELPAT(',NS,'): ', &
                 PELPAT(NS)
         END DO
         ENDIF

!      idx=0
!      DO NPEL=1,NPELMAX
!         IF((PELTOT(NPEL).GT.0.D0).OR.(ID.EQ.1)) THEN
!            idx=1
!            WRITE(6,632) NPEL, &
!                'MDLPEL',MDLPEL(NPEL), &
!                 'PELTOT',PELTOT(NPEL), &
!                 'PELR0 ',PELR0(NPEL), &
!                 'PELRW ',PELRW(NPEL)
!            WRITE(6,633) NPEL, &
!                 'PELRAD',PELRAD(NPEL), &
!                 'PELVEL',PELVEL(NPEL), &
!                 'PELTIM',PELTIM(NPEL)
!         END IF
!      END DO
!      IF(idx.EQ.1) THEN
!         WRITE(6,'(A24,4ES12.4)') &
!              'pellet_time_start', &
!              (pellet_time_start(NPEL),NPEL=1,NPELMAX)
!         WRITE(6,'(A24,4ES12.4)') &
!              'pellet_time_interval', &
!              (pellet_time_interval(NPEL),NPEL=1,NPELMAX)
!         WRITE(6,'(A24,4I12)') &
!              'number_of_pellet_repeat', &
!              (number_of_pellet_repeat(NPEL),NPEL=1,NPELMAX)
!         DO NS=1,NSMAX
!            WRITE(6,'(A7,I2,A3)') 'PELPAT(',NS,'): ', &
!                 (PELPAT(NS,NPEL),NPEL=1,NPELMAX)
!         END DO
!      END IF

      DO NPSC=1,NPSCMAX
         IF(MDLPSC(NPSC).GT.0) THEN
            WRITE(6,632) NPSC, &
                 'MDLPSC',MDLPSC(NPSC), &
                 'PSCIN ',PSCIN(NPSC), &
                 'PSCR0 ',PSCR0(NPSC) , &
                 'PSCRW ',PSCRW(NPSC)
            WRITE(6,632) NPSC, &
                 'NSPSC ',NSPSC(NPSC)
         END IF
      END DO
 
      IF(MDLPR.GE.1) THEN
         WRITE(6,623) 'MDLPR   ',MDLPR,   'SYNCABS ',SYNCABS, &
                      'SYNCSELF',SYNCSELF
      ENDIF

      WRITE(6,'(A,A)') 'KNAMEQ =',TRIM(knameq)
      WRITE(6,'(A,A)') 'KNAMEQ2=',TRIM(knameq2)
      WRITE(6,'(A,A)') 'KNAMTR =',TRIM(knamtr)
      WRITE(6,'(A,A)') 'KFNLOG =',TRIM(kfnlog)
      WRITE(6,'(A,A)') 'KFNTXT =',TRIM(kfntxt)
      WRITE(6,'(A,A)') 'KFNCVS =',TRIM(kfncvs)

      RETURN

601   FORMAT(' ',A6,'=',1PE11.3 :2X,A6,'=',1PE11.3: &
              2X,A6,'=',1PE11.3 :2X,A6,'=',1PE11.3)
602   FORMAT(' ',A6,'=',I7,4X   :2X,A6,'=',I7,4X  : &
              2X,A6,'=',I7,4X   :2X,A6,'=',I7)
603   FORMAT(' ',A7,'=',I6,4X   :2X,A6,'=',1PE11.3: &
              2X,A6,'=',1PE11.3 :2X,A6,'=',1PE11.3)
604   FORMAT(' ',A6,'=',I7,4X   :2X,A6,'=',1X,A6,4X: &
              2X,A6,'=',1X,A6,4X:2X,A6,'=',I7)
605   FORMAT(' ',A6,'=',I7,4X   :2X,A6,'=',I7,4X  : &
              2X,A6,'=',I7,4X   :2X,A6,'=',1PE11.3)
612   FORMAT(' ',I2,1PD12.4,0P,F8.3,4F9.4)
622   FORMAT(' ',A8,'=',I5,4X   :2X,A8,'=',I5,4X  : &
              2X,A8,'=',I5,4X   :2X,A8,'=',I5)
623   FORMAT(' ',A8,'=',I7,4X   :2X,A8,'=',1PE11.3: &
              2X,A8,'=',1PE11.3)
630   FORMAT(' ',A12,'=',I12)
631   FORMAT(' ',A12,'=',ES12.4)
632   FORMAT(' ',I2,1X,A6,I8,4X,3(1X,A6,ES12.4))
633   FORMAT(' ',I2,4(1X,A6,ES12.4))
634   FORMAT(' ',I2,2(1X,A6,I8,4X),2(1X,A6,ES12.4))
    END SUBROUTINE tr_view
END MODULE trview
