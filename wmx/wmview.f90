! wmview.f90

MODULE wmview

  PRIVATE
  PUBLIC wm_view

CONTAINS

!     ****** DISPLAY INPUT DATA ******

  SUBROUTINE wm_view
    USE wmcomm_parm
    IMPLICIT NONE
    INTEGER:: NS,NA,i

    SELECT CASE(MODELG)
    CASE(0)
       WRITE(6,*) '## MODELG=0: UNIFORM ##'
    CASE(1)
       WRITE(6,*) '## MODELG=1: CYLINDRICAL ##'
    CASE(2)
       WRITE(6,*) '## MODELG=2: TOROIDAL ##'
    CASE(3)
       WRITE(6,*) '## MODELG=3: TASK/EQ ##'
    CASE(4)
       WRITE(6,*) '## MODELG=4: VMEC ##'
    CASE(5)
       WRITE(6,*) '## MODELG=5: EQDSK ##'
    CASE(6)
       WRITE(6,*) '## MODELG=6: BOOZER ##'
    END SELECT

    SELECT CASE(MODELJ)
    CASE(0)
       WRITE(6,*) '## MODELJ=0: Loop antenna ##'
    CASE(1)
       WRITE(6,*) '## MODELJ=1: Waveguide ##'
    CASE(2)
       WRITE(6,*) '## MODELJ=2: POLOIDAL MODE ##'
    CASE(3)
       WRITE(6,*) '## MODELJ=3: TOROIDAL MODE ##'
    CASE DEFAULT
       WRITE(6,*) '## MODELJ=',MODELJ,': VACUUM EIGEN MODE ##'
    END SELECT

    SELECT CASE(MODELN)
    CASE(7)
       WRITE(6,*) '## MODELN=7: READ PROFILE DATA : WMDPRF ##'
    CASE(8)
       WRITE(6,*) '## MODELN=8: READ PROFILE DATA : WMXPRF ##'
    CASE(9)
       WRITE(6,*) '## MODELN=9: READ PROFILE DATA : TRDATA ##'
    END SELECT

    SELECT CASE(MODELA)
    CASE(1)
       WRITE(6,*) '## MODELA=1: ALPHA PARTICLE EFFECT ##'
    CASE(2)
       WRITE(6,*) '## MODELA=2: ELECTRON BETA EFFECT ##'
    CASE(3)
       WRITE(6,*) '## MODELA=3: ALPHA PARTICLE AND ELECTRON BETA EFFECTS ##'
    CASE(4)
       WRITE(6,*) '## MODELA=4: ALPHA PARTICLE DENSITY CALCULATED ##'
    END SELECT

    WRITE(6,601) 'BB    ',BB    ,'RR    ',RR    , &
                 'RA    ',RA    ,'RB    ',RB
    WRITE(6,601) 'Q0    ',Q0    ,'QA    ',QA    , &
                 'RKAP  ',RKAP  ,'RDLT  ',RDLT
    WRITE(6,601) 'ZEFF  ',ZEFF  ,'PNA   ',PNA   , &
                 'PNAL  ',PNAL  ,'PTA   ',PTA
    WRITE(6,601) 'RHOMIN',RHOMIN,'QMIN  ',QMIN  , &
                 'PRFIN ',PRFIN
    WRITE(6,601) 'RF    ',RF    ,'RFI   ',RFI   , &
         'RD    ',RD
    WRITE(6,604) 'factor_nth  ',factor_nth, &
                 'factor_nhh  ',factor_nhh, &
                 'factor_nph  ',factor_nph
    WRITE(6,602) 'NRMAX ',NRMAX ,'NTHMAX',NTHMAX, &
                 'NHHMAX',NHHMAX,'NPHMAX',NPHMAX
    WRITE(6,602) 'NTH0  ',NTH0  ,'NPH0  ',NPH0  , &
                 'NHC   ',NHC   
    WRITE(6,602) 'MODELG',MODELG,'MODELB',MODELB, &
                 'MODELN',MODELN,'MODELQ',MODELQ
    WRITE(6,602) 'MODELA',MODELA, & 
                 'MODEFR',MODEFR,'MODEFW',MODEFW
    WRITE(6,602) 'MODELM',MODELM, & 
                 'MDLWMK',MDLWMK, &
                 'MDLWMX',MDLWMX
    WRITE(6,604) 'MODEL_PROF  ',MODEL_PROF, &
                 'MODEL_NPROF ',MODEL_NPROF, &
                 'nthmax_g     ',nthmax_g

    WRITE(6,692)
    DO NS=1,NSMAX
       WRITE(6,611) NS,PA(NS),PZ(NS),PN(NS),PNS(NS), &
                       PTPR(NS),PTPP(NS),PTS(NS)
    ENDDO
    DO NS=1,NSMAX
       WRITE(6,612) NS,MODELP(NS),MODELV(NS),NCMIN(NS),NCMAX(NS), &
                       PZCL(NS),PU(NS),PUS(NS),PNITB(NS),PTITB(NS),PUITB(NS)
    ENDDO

    DO NS=1,NSMAX
       WRITE(6,613) NS,PROFN1(NS),PROFN2(NS), &
                       PROFT1(NS),PROFT2(NS), &
                       PROFU1(NS),PROFU2(NS),RHOITB(NS)
    END DO

    WRITE(6,693)
    DO NA=1,NAMAX
       WRITE(6,610) NA,AJ(NA),APH(NA),THJ1(NA),THJ2(NA), &
                                      PHJ1(NA),PHJ2(NA)
       WRITE(6,614)    BETAJ(NA),AEWGT(NA),AEWGZ(NA),ANTANG(NA)
    ENDDO

    DO i=1,idebug_max
       IF(idebuga(i).NE.0) WRITE(6,'(A,I2,A,I6)') 'idebuga(',i,'):',idebuga(i)
    END DO
    IF(idebuga(61).NE.0) WRITE(6,'(A,A)') 'knam_dump=',TRIM(knam_dump)
    RETURN

  601 FORMAT(A6,'=',1PE12.4:1X,A6,'=',1PE12.4:1X, &
             A6,'=',1PE12.4:1X,A6,'=',1PE12.4)
  602 FORMAT(A6,'=',I8,4X  :1X,A6,'=',I8,4X  :1X, &
             A6,'=',I8,4X  :1X,A6,'=',I8)
  604 FORMAT(A12,'=',I8,4X  :1X,A12,'=',I8,4X  :1X, &
             A12,'=',I8)
  610 FORMAT(' ',I1,6(1PE11.3))
  611 FORMAT(' ',I1,7(1PE11.3))
  612 FORMAT(' ',I1,3I3,I2,6(1PE11.3))
  613 FORMAT(' ',I1,7(1PE11.3))
  614 FORMAT(' ',1X,6(1PE11.3))
  692 FORMAT(' ','NS    PA',9X,'PZ',9X,'PN',9X,'PNS', &
                            8X,'PTPR',7X,'PTPP',7X,'PTS'/ &
             ' ','  MP MV ND1 ND2',2X,'PZCL',7X,'PU',9X,'PUS', &
                            8X,'PNITB',6X,'PTITB',6X,'PUITB'/ & 
             ' ','  PROFN1',5X,'PROFN2',5X,'PROFT1',5X,'PROFT2', &
                            5X,'PROFU1',5X,'PROFU2',5X,'RHOITB')
  693 FORMAT(' ','NA    AJ',9X,'APH',8X,'THJ1',7X,'THJ2', &
                            7X,'PHJ1',7X,'PHJ2')
  END SUBROUTINE wm_view
END MODULE wmview
