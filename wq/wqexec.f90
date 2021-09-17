! wqexec.f90

MODULE wqexec

  PRIVATE
  PUBLIC wq_exec

CONTAINS

  SUBROUTINE wq_exec
    USE wqcomm
    USE wqsolv
    IMPLICIT NONE
    INTEGER:: nx,ny,k,i,j,nt,N,NA,ILL
    REAL(rkind):: t,pabs_tot,factor
    REAL(rkind):: sys_len,y0,y,pulse_width

    t=ttot
    do nt=1,ntmax
       nttot=nttot+1
       ! give Boundary Condition
       SELECT CASE(model_pulse)
       CASE(0)
          factor=1.D0
       CASE(1)
          IF(ttot.LE.pulse_cycle) THEN
             WRITE(21,'(A,2I8,3ES12.4)') 'tot:',nt,nttot,t,ttot,pulse_cycle
             factor=1.D0
          ELSE
             factor=0.D0
          END IF
       CASE(2)
          IF(ttot.LE.6.D0*pulse_cycle) THEN
             factor=EXP(-(t-3.D0*pulse_cycle)**2)
          ELSE
             factor=0.D0
          END IF
       END SELECT

       nx=nxmax
       sys_len=dy*(nymax-1)
       y0=0.5D0*sys_len
       pulse_width=2.D0*wavelength
       do ny=1,nymax
          y=dy*(ny-1)
          if(INMODE.eq.1) then
             EY(nx,ny) = factor/sqrt(2.d0)*exp(-(y-y0)**2/pulse_width**2)
             EZ(nx,ny) = factor/sqrt(2.d0)*exp(-(y-y0)**2/pulse_width**2)
          else if(INMODE.eq.2) then
             EY(nx,ny) = factor*exp(-(y-y0)**2/pulse_width**2)
          else if(INMODE.eq.3) then
             EZ(nx,ny) = factor*exp(-(y-y0)**2/pulse_width**2)
          else
             write(*,*) "ERROR:INMODE needs to be 1or2or3"
          end if
       end do

       !compute next E

       call wq_solv

       t=t+dt
       ttot=ttot+dt

       if(mod(nt,1000).eq.0.OR.nt.EQ.ntmax) THEN
          DO ny=1,nymax
             pabs(    1,ny)=0.D0
             pabs(nxmax,ny)=0.D0
          END DO
          DO nx=1,nxmax
             pabs(nx,    1)=0.D0
             pabs(nx,nymax)=0.D0
          END DO
          
        ! compute Absorbed Power
          do ny=2,nymax-1
             do nx=2,nxmax-1
                pabs(nx,ny)=0.5d0*real( &
                     conjg(EX(nx,ny))*(CD(1,1,nx,ny)*EX(nx,ny) &
                                      +CD(1,2,nx,ny)*EY(nx,ny) &
                                      +CD(1,3,nx,ny)*EZ(nx,ny))&
                    +conjg(EY(nx,ny))*(CD(2,1,nx,ny)*EX(nx,ny) &
                                      +CD(2,2,nx,ny)*EY(nx,ny) &
                                      +CD(2,3,nx,ny)*EZ(nx,ny))&
                    +conjg(EZ(nx,ny))*(CD(3,1,nx,ny)*EX(nx,ny) &
                                      +CD(3,2,nx,ny)*EY(nx,ny) &
                                      +CD(3,3,nx,ny)*EZ(nx,ny)))
             end do
          end do
          pabs_tot=0.d0
          do ny=1,nymax
             do nx=1,nxmax
                pabs_tot=pabs_tot+pabs(nx,ny)
             end do
          end do

          WRITE(6,'(A,I8,2ES12.4)') '## nt,t,pabs:',nt,t,pabs_tot

       end if

       IF(MOD(nttot,ntplot_interval).EQ.0.AND.ntplot.LT.ntplot_max) THEN
          ntplot=ntplot+1
          DO ny=1,nymax
             DO nx=1,nxmax
                EX_save(nx,ny,ntplot)=EX(nx,ny)
                EY_save(nx,ny,ntplot)=EY(nx,ny)
                EZ_save(nx,ny,ntplot)=EZ(nx,ny)
             END DO
          END DO
          WRITE(6,'(A,2I8,ES12.4)') &
               '## Plot saved: ntplot,nttot,ttot=',ntplot,nttot,ttot
       END IF
    end do

    RETURN

  END SUBROUTINE wq_exec
END MODULE wqexec
