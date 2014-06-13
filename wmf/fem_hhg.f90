!   $Id$

MODULE libfem_hhg

  PRIVATE
  PUBLIC fem_hhg

!----- calculate coefficint matrix fmd ---
!          with hhg: hermite, hermite, third-order-polnomial

!     h     : Hermite (continuous, derivative-continuous, 2 nodes)
!           : h1=(1-x)^2(1+2x)   dh1=-6x(1-x)
!           : h2= x(1-x)^2       dh2=(1-x)(1-3x)
!           : h3= x^2(3-2x)      dh3= 6x(1-x)
!           : h4=-x^2(1-x)       dh4=x(-2+3x)
!     g     : Polynomial-third order (discontinuous, 1 node)
!           : g1=1               dg1=0
!           : g2=x-1/2           dg2=1
!           : g3=(x-1/2)^2/2     dg3=(x-1/2)
!           : g4=(x-1/2)^3/6     dg4=(x-1/2)^2/2

!     h     : h1=1 h1'=0 at x=0, h1=0 h1'=0 at x=1
!           : h2=0 h2'=1 at x=0, h2=0 h2'=0 at x=1
!           : h3=0 h3'=1 at x=0, h3=1 h3'=0 at x=1
!           : h4=0 h4'=1 at x=0, h4=0 h4'=1 at x=1
!     g     : g1=1 g1'=0 g1''=0 g1'''=0  at x=1/2
!           : g2=0 g2'=1 g2''=0 g2'''=0  at x=1/2
!           : g3=0 g3'=0 g3''=1 g2'''=0  at x=1/2
!           : g4=0 g4'=0 g4''=0 g4'''=1  at x=1/2

CONTAINS
  SUBROUTINE fem_hhg(fmd,drho,fml)
    USE libfem,ONLY: table_hhg,table_initialize_flag,table_initialize
    USE wmfem_comm,ONLY: mbmax,nfcmax
    IMPLICIT NONE
    COMPLEX(8),DIMENSION(4,4,4,nfcmax,nfcmax,4),INTENT(IN) :: fmd
    REAL(8),INTENT(IN):: drho
    COMPLEX(8),DIMENSION(mbmax,mbmax),INTENT(OUT):: fml
    INTEGER:: mr,mc,i,j,k,nf1,nf2,inod,mb1,mb2,ig

    IF(table_initialize_flag.eq.0) THEN
       CALL table_initialize
       table_initialize_flag=1
    ENDIF

    mr=8*nfcmax  ! line interval between radial points
    mc=mbmax     ! diagonal column at the top line

    DO nf2=1,nfcmax
       DO nf1=1,nfcmax
          DO j=1,4
             DO i=1,4
                mb1=8*(nf1-1)+2*(i-1)+1
                mb2=8*(nf2-1)+2*(j-1)+1
                DO inod=1,4
                   ig=2*inod-1
                   fml(mb1     ,mb2     )=fml(mb1     ,mb2     ) &
                        +fmd(i,j,1,nf1,nf2,inod)*table_hhg(1,1,ig)*drho &
                        +fmd(i,j,2,nf1,nf2,inod)*table_hhg(5,1,ig) &
                        +fmd(i,j,3,nf1,nf2,inod)*table_hhg(1,5,ig) &
                        +fmd(i,j,4,nf1,nf2,inod)*table_hhg(5,5,ig)/drho
                   fml(mb1+1   ,mb2     )=fml(mb1+1   ,mb2     ) &
                        +fmd(i,j,1,nf1,nf2,inod)*table_hhg(1,2,ig)*drho**2 &
                        +fmd(i,j,2,nf1,nf2,inod)*table_hhg(5,2,ig)*drho &
                        +fmd(i,j,3,nf1,nf2,inod)*table_hhg(1,6,ig)*drho &
                        +fmd(i,j,4,nf1,nf2,inod)*table_hhg(5,6,ig)
                   fml(mb1+mr  ,mb2     )=fml(mb1+mr  ,mb2     ) &
                        +fmd(i,j,1,nf1,nf2,inod)*table_hhg(1,3,ig)*drho &
                        +fmd(i,j,2,nf1,nf2,inod)*table_hhg(5,3,ig) &
                        +fmd(i,j,3,nf1,nf2,inod)*table_hhg(1,7,ig) &
                        +fmd(i,j,4,nf1,nf2,inod)*table_hhg(5,7,ig)/drho
                   fml(mb1+mr+1,mb2     )=fml(mb1+mr+1,mb2     ) &
                        +fmd(i,j,1,nf1,nf2,inod)*table_hhg(1,4,ig)*drho**2 &
                        +fmd(i,j,2,nf1,nf2,inod)*table_hhg(5,4,ig)*drho &
                        +fmd(i,j,3,nf1,nf2,inod)*table_hhg(1,8,ig)*drho &
                        +fmd(i,j,4,nf1,nf2,inod)*table_hhg(5,8,ig)

                   fml(mb1     ,mb2+1   )=fml(mb1     ,mb2+1   ) &
                        +fmd(i,j,1,nf1,nf2,inod)*table_hhg(2,1,ig)*drho**2 &
                        +fmd(i,j,2,nf1,nf2,inod)*table_hhg(6,1,ig)*drho &
                        +fmd(i,j,3,nf1,nf2,inod)*table_hhg(2,5,ig)*drho &
                        +fmd(i,j,4,nf1,nf2,inod)*table_hhg(6,5,ig)
                   fml(mb1+1   ,mb2+1   )=fml(mb1+1   ,mb2+1   ) &
                        +fmd(i,j,1,nf1,nf2,inod)*table_hhg(2,2,ig)*drho**3 &
                        +fmd(i,j,2,nf1,nf2,inod)*table_hhg(6,2,ig)*drho**2 &
                        +fmd(i,j,3,nf1,nf2,inod)*table_hhg(2,6,ig)*drho**2 &
                        +fmd(i,j,4,nf1,nf2,inod)*table_hhg(6,6,ig)*drho
                   fml(mb1+mr  ,mb2+1   )=fml(mb1+mr  ,mb2+1   ) &
                        +fmd(i,j,1,nf1,nf2,inod)*table_hhg(2,3,ig)*drho**2 &
                        +fmd(i,j,2,nf1,nf2,inod)*table_hhg(6,3,ig)*drho &
                        +fmd(i,j,3,nf1,nf2,inod)*table_hhg(2,7,ig)*drho &
                        +fmd(i,j,4,nf1,nf2,inod)*table_hhg(6,7,ig)
                   fml(mb1+mr+1,mb2+1   )=fml(mb1+mr+1,mb2+1   ) &
                        +fmd(i,j,1,nf1,nf2,inod)*table_hhg(2,4,ig)*drho**3 &
                        +fmd(i,j,2,nf1,nf2,inod)*table_hhg(6,4,ig)*drho**2 &
                        +fmd(i,j,3,nf1,nf2,inod)*table_hhg(2,8,ig)*drho**2 &
                        +fmd(i,j,4,nf1,nf2,inod)*table_hhg(6,8,ig)*drho


                   fml(mb1     ,mb2+mr  )=fml(mb1    ,mb2+mr  ) &
                        +fmd(i,j,1,nf1,nf2,inod)*table_hhg(3,1,ig)*drho &
                        +fmd(i,j,2,nf1,nf2,inod)*table_hhg(7,1,ig) &
                        +fmd(i,j,3,nf1,nf2,inod)*table_hhg(3,5,ig) &
                        +fmd(i,j,4,nf1,nf2,inod)*table_hhg(7,5,ig)/drho
                   fml(mb1+1   ,mb2+mr  )=fml(mb1+1  ,mb2+mr  ) &
                        +fmd(i,j,1,nf1,nf2,inod)*table_hhg(3,2,ig)*drho**2 &
                        +fmd(i,j,2,nf1,nf2,inod)*table_hhg(7,2,ig)*drho &
                        +fmd(i,j,3,nf1,nf2,inod)*table_hhg(3,6,ig)*drho &
                        +fmd(i,j,4,nf1,nf2,inod)*table_hhg(7,6,ig)
                   fml(mb1+mr  ,mb2+mr  )=fml(mb1+mr  ,mb2+mr  ) &
                        +fmd(i,j,1,nf1,nf2,inod)*table_hhg(3,3,ig)*drho &
                        +fmd(i,j,2,nf1,nf2,inod)*table_hhg(7,3,ig) &
                        +fmd(i,j,3,nf1,nf2,inod)*table_hhg(3,7,ig) &
                        +fmd(i,j,4,nf1,nf2,inod)*table_hhg(7,7,ig)/drho
                   fml(mb1+mr+1,mb2+mr  )=fml(mb1+mr+1,mb2+mr  ) &
                        +fmd(i,j,1,nf1,nf2,inod)*table_hhg(3,4,ig)*drho**2 &
                        +fmd(i,j,2,nf1,nf2,inod)*table_hhg(7,4,ig)*drho &
                        +fmd(i,j,3,nf1,nf2,inod)*table_hhg(3,8,ig)*drho &
                        +fmd(i,j,4,nf1,nf2,inod)*table_hhg(7,8,ig)

                   fml(mb1     ,mb2+mr+1)=fml(mb1     ,mb2+mr+1) &
                        +fmd(i,j,1,nf1,nf2,inod)*table_hhg(4,1,ig)*drho**2 &
                        +fmd(i,j,2,nf1,nf2,inod)*table_hhg(8,1,ig)*drho &
                        +fmd(i,j,3,nf1,nf2,inod)*table_hhg(4,5,ig)*drho &
                        +fmd(i,j,4,nf1,nf2,inod)*table_hhg(8,5,ig)
                   fml(mb1+1   ,mb2+mr+1)=fml(mb1+1   ,mb2+mr+1) &
                        +fmd(i,j,1,nf1,nf2,inod)*table_hhg(4,2,ig)*drho**3 &
                        +fmd(i,j,2,nf1,nf2,inod)*table_hhg(8,2,ig)*drho**2 &
                        +fmd(i,j,3,nf1,nf2,inod)*table_hhg(4,6,ig)*drho**2 &
                        +fmd(i,j,4,nf1,nf2,inod)*table_hhg(8,6,ig)*drho
                   fml(mb1+mr  ,mb2+mr+1)=fml(mb1+mr  ,mb2+mr+1) &
                        +fmd(i,j,1,nf1,nf2,inod)*table_hhg(4,3,ig)*drho**2 &
                        +fmd(i,j,2,nf1,nf2,inod)*table_hhg(8,3,ig)*drho &
                        +fmd(i,j,3,nf1,nf2,inod)*table_hhg(4,7,ig)*drho &
                        +fmd(i,j,4,nf1,nf2,inod)*table_hhg(8,7,ig)
                   fml(mb1+mr+1,mb2+mr+1)=fml(mb1+mr+1,mb2+mr+1) &
                        +fmd(i,j,1,nf1,nf2,inod)*table_hhg(4,4,ig)*drho**3 &
                        +fmd(i,j,2,nf1,nf2,inod)*table_hhg(8,4,ig)*drho**2 &
                        +fmd(i,j,3,nf1,nf2,inod)*table_hhg(4,8,ig)*drho**2 &
                        +fmd(i,j,4,nf1,nf2,inod)*table_hhg(8,8,ig)*drho
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    RETURN
  END SUBROUTINE fem_hhg
END MODULE libfem_hhg
