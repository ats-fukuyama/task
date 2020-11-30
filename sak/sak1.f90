! sak1.f90

MODULE sak1

  PRIVATE
  PUBLIC sak_1

CONTAINS

  SUBROUTINE sak_1

    USE sakcomm
    USE libgrf
    IMPLICIT NONE
    REAL(dp):: rkmin,rkmax,sgmin,sgmax
    REAL(dp):: rkminl,rkmaxl,sgminl,sgmaxl,delrkl,delsgl
    REAL(dp),ALLOCATABLE:: rkl(:),sgl(:),f1(:,:),f2(:,:),f3(:,:),f4(:,:)
    INTEGER:: nrkmax,nsgmax,nrk,nsg
    INTEGER:: nlmax,npmax,nl,np
    REAL(dp),ALLOCATABLE:: line_value(:),line_rgb(:,:)
    INTEGER,ALLOCATABLE:: line_pat(:)
    REAL(dp),ALLOCATABLE:: paint_value(:),paint_rgb(:,:)
    EXTERNAL PAGES,PAGEE

    rkmin=0.01D0
    rkmax=100.D0
    nrkmax=51
    sgmin=0.01D0
    sgmax=100.D0
    nsgmax=51

1   CONTINUE

    WRITE(6,'(A)') '## INPUT rkmin,rkmax,sgmin,sgmax?'
    WRITE(6,'(A,4ES12.4)') '      ',rkmin,rkmax,sgmin,sgmax
    READ(5,*,ERR=1,END=9000) rkmin,rkmax,sgmin,sgmax

    rkminl=LOG10(rkmin)
    rkmaxl=LOG10(rkmax)
    sgminl=LOG10(sgmin)
    sgmaxl=LOG10(sgmax)
    delrkl=(rkmaxl-rkminl)/(nrkmax-1)
    delsgl=(sgmaxl-sgminl)/(nsgmax-1)
    ALLOCATE(rkl(nrkmax),sgl(nsgmax))
    ALLOCATE(f1(nrkmax,nsgmax),f2(nrkmax,nsgmax))
    ALLOCATE(f3(nrkmax,nsgmax),f4(nrkmax,nsgmax))

    nlmax=13
    npmax=13
    ALLOCATE(line_value(nlmax),line_rgb(3,nlmax),line_pat(nlmax))
    ALLOCATE(paint_value(npmax),paint_rgb(3,npmax))
    DO nl=1,nlmax
       line_value(nl)=-3.D0+0.5D0*(nl-1)
       line_pat(nl)=MOD(nl-1,2)
       IF(nl.LT.7) THEN
          line_rgb(1,nl)=0.D0
          line_rgb(2,nl)=0.D0
          line_rgb(3,nl)=1.D0
       ELSE IF(nl.EQ.7) THEN
          line_rgb(1,nl)=0.D0
          line_rgb(2,nl)=0.D0
          line_rgb(3,nl)=0.D0
       ELSE
          line_rgb(1,nl)=1.D0
          line_rgb(2,nl)=0.D0
          line_rgb(3,nl)=0.D0
       END IF
    END DO
    DO np=1,npmax
       paint_value(np)=-5.5D0+(np-1)
       IF(np.LE.7) THEN
          paint_rgb(1,np)=(np-1)/6.D0
          paint_rgb(2,np)=(np-1)/6.D0
          paint_rgb(3,np)=1.0D0
       ELSE
          paint_rgb(1,np)=1.0D0
          paint_rgb(2,np)=(13-np)/6.D0
          paint_rgb(3,np)=(13-np)/6.D0
       END IF
    END DO

    DO nsg=1,nsgmax
       sgl(nsg)=sgminl+delsgl*(nsg-1)
    END DO
    DO nrk=1,nrkmax
       rkl(nrk)=rkminl+delrkl*(nrk-1)
    END DO
    DO nsg=1,nsgmax
       DO nrk=1,nrkmax
          CALL sub_sak11(rkl(nrk),sgl(nsg),f1(nrk,nsg))
          CALL sub_sak12(rkl(nrk),sgl(nsg),f2(nrk,nsg))
          f3(nrk,nsg)=LOG10(MAX(f1(nrk,nsg),1.D-8))
          f4(nrk,nsg)=LOG10(MAX(f2(nrk,nsg),1.D-8))
       END DO
    END DO

    CALL PAGES
    CALL grd2d(1,rkl,sgl,f3,nrkmax,nrkmax,nsgmax, &
         '@gamma OAM SAK(klamdad,sigma)@',3, &
         MODE_2d=1,NLMAX=nlmax,ASPECT=1.D0,NOINFO=1, &
         LINE_VALUE=line_value,LINE_RGB=line_rgb, LINE_PAT=line_pat, &
         SCALE_THICKNESS=0.035D0, &
         XMIN=rkminl,XMAX=rkmaxl,YMIN=sgminl,YMAX=sgmaxl)
    CALL grd2d(2,rkl,sgl,f4,nrkmax,nrkmax,nsgmax, &
         '@gamma OAM JTM(klamdad,sigma)@',3, &
         MODE_2d=1,NLMAX=nlmax,ASPECT=1.D0,NOINFO=1, &
         LINE_VALUE=line_value,LINE_RGB=line_rgb, LINE_PAT=line_pat, &
         SCALE_THICKNESS=0.035D0, &
         XMIN=rkminl,XMAX=rkmaxl,YMIN=sgminl,YMAX=sgmaxl)
!    line_rgb(1,7)=0.D0
!    line_rgb(2,7)=1.D0
!    line_rgb(3,7)=0.D0
!    CALL grd2d(3,rkl,sgl,f3,nrkmax,nrkmax,nsgmax, &
!         '@gamma OAM SAK(klamdad,sigma)@',3, &
!         MODE_2d=3,NLMAX=nlmax,ASPECT=1.D0,NOINFO=1, &
!         LINE_VALUE=line_value,LINE_RGB=line_rgb, LINE_PAT=line_pat, &
!         SCALE_THICKNESS=0.035D0)
!    CALL grd2d(4,rkl,sgl,f4,nrkmax,nrkmax,nsgmax, &
!         '@gamma OAM JTM(klamdad,sigma)@',3, &
!         MODE_2d=3,NLMAX=nlmax,ASPECT=1.D0,NOINFO=1, &
!         LINE_VALUE=line_value,LINE_RGB=line_rgb, LINE_PAT=line_pat, &
!         SCALE_THICKNESS=0.035D0)
    CALL grd2d(3,rkl,sgl,f3,nrkmax,nrkmax,nsgmax, &
         '@gamma OAM SAK(klamdad,sigma)@',3, &
         MODE_2d=3,NLMAX=nlmax,NPMAX=npmax,ASPECT=1.D0,NOINFO=1, &
         LINE_VALUE=line_value,LINE_RGB=line_rgb, LINE_PAT=line_pat, &
         PAINT_VALUE=paint_value,PAINT_RGB=paint_RGB, &
         SCALE_THICKNESS=0.035D0, &
         XMIN=rkminl,XMAX=rkmaxl,YMIN=sgminl,YMAX=sgmaxl)
    CALL grd2d(4,rkl,sgl,f4,nrkmax,nrkmax,nsgmax, &
         '@gamma OAM JTM(klamdad,sigma)@',3, &
         MODE_2d=3,NLMAX=nlmax,NPMAX=npmax,ASPECT=1.D0,NOINFO=1, &
         LINE_VALUE=line_value,LINE_RGB=line_rgb, LINE_PAT=line_pat, &
         PAINT_VALUE=paint_value,PAINT_RGB=paint_RGB, &
         SCALE_THICKNESS=0.035D0, &
         XMIN=rkminl,XMAX=rkmaxl,YMIN=sgminl,YMAX=sgmaxl)
    DEALLOCATE(line_value,line_rgb,line_pat)
    DEALLOCATE(paint_value,paint_rgb)
    CALL PAGEE
    DEALLOCATE(rkl,sgl,f1,f2,f3,f4)
    GO TO 1

9000 CONTINUE
    RETURN
  END SUBROUTINE sak_1

  SUBROUTINE sub_sak11(rkl,sgl,f1)
    USE sakcomm
    IMPLICIT NONE
    REAL(dp),INTENT(IN):: rkl,sgl
    REAL(dp),INTENT(out):: f1
    REAL(DP):: rk,sg

    rk=10.D0**rkl
    sg=10.D0**sgl
    f1=SQRT(Pi/8.D0) &
         *(1.D0/(rk**3*SQRT(1.D0+1.D0/sg**2)**3)) &
         *EXP(-1.D0/(2.D0*rk**2*(1.D0+1.D0/sg**2)))
    RETURN
  END SUBROUTINE sub_sak11

  SUBROUTINE sub_sak12(rkl,sgl,f2)
    USE sakcomm
    IMPLICIT NONE
    REAL(dp),INTENT(IN):: rkl,sgl
    REAL(dp),INTENT(out):: f2
    REAL(DP):: rk,sg

    rk=10.D0**rkl
    sg=10.D0**sgl
    f2=SQRT(Pi/8.D0) &
         *(1.D0/((rk**3)*(1.D0+1.D0/sg**2))) &
         *(EXP(-1.D0/rk**2)*EXP(-1.D0/(2.D0*sg**2)) &
         +sg*EXP(-sg**2/rk**2)*EXP(-sg**2/2.D0))
    RETURN
  END SUBROUTINE sub_sak12
END MODULE sak1
