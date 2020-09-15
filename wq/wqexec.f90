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
    REAL(rkind):: t,pabs_tot
   
    t=ttot
    do nt=1,ntmax
       ! give Boundary Condition
       nx=nxmax
       do ny=1,nymax
          if(INMODE.eq.1) then
             EY(nx,ny) = 1.0d0/sqrt(2.d0) &
                  *exp(-5.0d-4/wavelength*(((ny-nymax/2))**2))
             EZ(nx,ny) = 1.0d0/sqrt(2.d0) &
                  *exp(-5.0d-4/wavelength*(((ny-nymax/2))**2))
          else if(INMODE.eq.2) then
             EY(nx,ny) = 1.0d0*exp(-5.0d-4/wavelength*(((ny-nymax/2))**2))
          else if(INMODE.eq.3) then
             EZ(nx,ny) = 1.0d0*exp(-5.0d-4/wavelength*(((ny-nymax/2))**2))
          else
             write(*,*) "ERROR:INMODE needs to be 1or2or3"
          end if
       end do

       !compute next E
       call wq_solv
       t=t+dt
       if(mod(nt,1000).eq.0.OR.nt.EQ.ntmax)then
          write(6,'(I4)') (ntmax-nt)/1000

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
       end if
    end do

    ttot=t
    nttot=nttot+ntmax
    WRITE(6,'(A,I6,2ES12.4)') '## nt,t,pabs:',nttot,t,pabs_tot
    RETURN

  END SUBROUTINE wq_exec
END MODULE wqexec
