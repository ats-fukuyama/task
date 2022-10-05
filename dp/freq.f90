! freq.f90

PROGRAM freq
  USE bpsd
  USE libgrf
  IMPLICIT NONE
  REAL(rkind):: pnemin,pnemax,rfmin,rfmax
  REAL(rkind):: LOG_pnemin,LOG_pnemax,LOG_rfmin,LOG_rfmax
  REAL(rkind):: pa,pz,bb
  INTEGER:: nxmax
  
  REAL(rkind):: pne,pni,rf
  REAL(rkind):: pnemin_temp,f1,f2,f3,del_LOG_pne,LOG_pne,PI2
  REAL(rkind):: omega,omega_ce,omega_ci,omega_pe,omega_pi
  REAL(rkind):: omega_UH,omega_LH,omega_CR,omega_CL
  INTEGER:: nx,i
  REAL(rkind),ALLOCATABLE:: gx(:),gf(:,:)
  REAL(rkind),ALLOCATABLE:: x(:),f(:,:)
  REAL(rkind):: LINE_RGB(3,9)
  INTEGER:: LINE_PAT(9)
  INTEGER,PARAMETER:: NLMAX=9
  DATA LINE_RGB/1.0,0.0,0.0, 0.0,0.0,1.0, 1.0,0.0,0.0, 0.0,0.0,1.0, &
                1.0,0.0,0.0, 0.0,0.0,1.0, 1.0,0.0,0.0, 0.0,0.0,1.0, &
                0.0,1.0,0.0/
  DATA LINE_PAT/0,0, 1,1, 2,2, 3,3, 4/

  CALL GSOPEN
  PI2=2.D0*PI
  
!  pnemin=1.D12
!  pnemax=1.D20
!  rfmin=1.D3
!  rfmax=1.D12
!  rf=13.56D6
!  pa=39.D0
!  pz=1.D0
!  bb=0.04D0
!  nxmax=101

  pnemin=1.D18
  pnemax=1.D20
  rfmin=1.D6
  rfmax=1.D11
  rf=13.56D6
  pa=1.D0
  pz=1.D0
  bb=2.5D0
  nxmax=101

1 WRITE(6,*)
  WRITE(6,'(A,1P4E12.4)') &
       'pnemin/max[1/m^3],rfmin/max[MHz]=', &
       pnemin,pnemax,rfmin,rfmax
  WRITE(6,'(A,1P3E12.4,I5)') 'pa,pz,bb=', pa,pz,bb,nxmax
  WRITE(6,'(A)') &
       '## INPUT : pnemin,pnemax,rfmin,rfmax,pa,pz,bb,nxmax (0/ for end)'
  pnemin_temp=pnemin
  READ(5,*,ERR=1,END=9) pnemin_temp,pnemax,rfmin,rfmax,pa,pz,bb,nxmax
  IF(pnemin_temp.LE.0.D0) GOTO 9
  IF(nxmax.LE.0) GOTO 9

  ALLOCATE(gx(nxmax),gf(nxmax,9))
  ALLOCATE(x(nxmax),f(nxmax,9))
  pnemin=pnemin_temp

  omega_ce=AEE*bb/AME
  omega_ci=pz*AEE*bb/(pa*AMP)

  LOG_pnemin=LOG10(pnemin)
  LOG_pnemax=LOG10(pnemax)
  del_LOG_pne=(LOG_pnemax-LOG_pnemin)/(nxmax-1)

  WRITE(6,'(A)') 'f_ce/ci/pe/pi/UH/LH/CR/CL:'
  DO nx=1,nxmax
     LOG_pne=LOG_pnemin+del_LOG_pne*(nx-1)
     pne=10.D0**(LOG_pne)
     pni=pne/pz
     omega_pe=SQRT(pne*AEE**2/(AME*EPS0))
     omega_pi=SQRT(pni*pz**2*AEE**2/(pa*AMP*EPS0))
     f1= omega_ce**2+omega_ci**2+omega_pe**2+omega_pi**2
     f2=(omega_ce**2-omega_ci**2+omega_pe**2-omega_pi**2)**2 &
          +4.D0*omega_pe**2*omega_pi**2
     omega_UH=SQRT(0.5D0*(f1+SQRT(f2)))
!     WRITE(6,'(A,1P4E12.4)') 'f1,f2=',f1,f2,SQRT(f2),f1-SQRT(f2)
     IF(f1-SQRT(f2).LE.0.D0) THEN
        omega_LH=omega_ci
     ELSE
        omega_LH=SQRT(0.5D0*(f1-SQRT(f2)))
     END IF
     f3=(omega_ce+omega_ci)**2+4.D0*(omega_pe**2+omega_pi**2)
     omega_CR=0.5D0*(SQRT(f3)+omega_ce-omega_ci)
     omega_CL=0.5D0*(SQRT(f3)-omega_ce+omega_ci)
     x(nx)=pne
     f(nx,1)=omega_ce/PI2
     f(nx,2)=omega_ci/PI2
     f(nx,3)=omega_pe/PI2
     f(nx,4)=omega_pi/PI2
     f(nx,5)=omega_UH/PI2
     f(nx,6)=omega_LH/PI2
     f(nx,7)=omega_CR/PI2
     f(nx,8)=omega_CL/PI2
     f(nx,9)=rf
     gx(nx)=LOG10(x(nx))
     DO i=1,9
        gf(nx,i)=LOG10(f(nx,i))
     END DO
     WRITE(6,'(ES8.0,8ES9.1)') x(nx),(f(nx,i),i=1,8)
  END DO

  CALL PAGES
  CALL grd1d(0,gx,gf,nxmax,nxmax,9,&
       '@f_ce/ci/pe/pi/UH/LH/CR/CL [Hz] vs n_e [1/m^3]@',3, &
       FMIN=LOG10(rfmin),FMAX=LOG10(rfmax),NLMAX=NLMAX, &
       LINE_RGB=LINE_RGB,LINE_PAT=LINE_PAT)
  CALL PAGEE
  DEALLOCATE(gx,gf)

  GOTO 1

9 CONTINUE
  CALL GSCLOS
  STOP
END PROGRAM freq
     
