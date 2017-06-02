module eqread_mod

!F_UFMTENDIAN=big:17

contains

!=======================================================================
! Allocate and deallocate arrays associated with equilibria
!=======================================================================

  subroutine alloc_equ(mode)
    use equ_params
    IMPLICIT NONE
    integer(4), intent(in) :: mode

    select case(mode)
    case(1)
       ! *** allocate arrays that are widely used ***
       allocate(rg(irdm),zg(izdm2),psi(irzdm2),rbp(irzdm2))
       allocate(pds(ivdm),fds(ivdm),vlv(ivdm),qqv(ivdm),prv(ivdm))
       allocate(csu(isrzdm),rsu(isrzdm),zsu(isrzdm))
       allocate(hiv(ivdm),siv(ivdm),siw(ivdm),sdw(ivdm),ckv(ivdm),ssv(ivdm),aav(ivdm),rrv(ivdm), &
            &   rbv(ivdm),arv(ivdm),bbv(ivdm),biv(ivdm),r2b2v(ivdm),shv(ivdm),grbm2v(ivdm), &
            &   rov(ivdm),aiv(ivdm),brv(ivdm),epsv(ivdm),elipv(ivdm),trigv(ivdm),ftv(ivdm))
       
    case(2)
       ! *** allocate arrays that are used only when reading an equilibium data ***
       allocate(ieqout(10),ieqerr(10),icp(10),cp(10))
       allocate(ivac(0:nsfix),ncoil(0:nsfix),cvac(0:nsfix),rvac(0:nsfix),zvac(0:nsfix))
       allocate(rcoil(100,icvdm),zcoil(100,icvdm),ccoil(100,icvdm),rlimt(200),zlimt(200))

    case(-1)
       ! *** deallocate arrays that are widely used ***
       if(allocated(rg)) then
          deallocate(rg,zg,psi,rbp)
          deallocate(pds,fds,vlv,qqv,prv)
          deallocate(csu,rsu,zsu)
          deallocate(hiv,siv,siw,sdw,ckv,ssv,aav,rrv, &
               &     rbv,arv,bbv,biv,r2b2v,shv,grbm2v,rov, &
               &     aiv,brv,epsv,elipv,trigv,ftv)
       end if

    case(-2)
       ! *** deallocate arrays that are used only when reading an equilibium data ***
       deallocate(ieqout,ieqerr,icp,cp)
       deallocate(ivac,ncoil,cvac,rvac,zvac)
       deallocate(rcoil,zcoil,ccoil,rlimt,zlimt)

    end select

  end subroutine alloc_equ

!=======================================================================
!             eqdsk               level=32       date=81.12.23
!=======================================================================
  subroutine eqdsk
!=======================================================================
    use bpsd_constants, only : cnpi => PI, rMU0
    use equ_params
    use libgrf

    logical :: lex
    integer(4) :: i, j, n, ir, iz, ist
    real(8) :: rsep, zsep, rrmax, rrmin, zzmax, zzmin, rzmax, rzmin, rpmax, rpmin, betp, zzlam
    real(8) :: dpsi, btv2, dsr, dsz, dxx, dl
    real(8):: xxx(ivdm),dummy(1000)
    integer(4):: idummy(200)
    CHARACTER(LEN=3):: K3
    CHARACTER(LEN=110):: K110
!=======================================================================
    inquire(file=eqfile,exist=lex)
    if( lex .neqv. .true. ) stop 'No equilibrium file.'
    open(ieqrd,file=eqfile,iostat=ist,form='unformatted',access='sequential',action='read')
    if(ist /= 0) stop 'file open error.'
!-----------------------------------------------------------------------

    call alloc_equ(1)
    call alloc_equ(2)

    read(ieqrd,iostat=ist) K3,K110
    write(6,*) K110
    read(ieqrd,iostat=ist)nsr,nsz  &
         &        ,(rg(i),i=1,nsr),(zg(j),j=1,nsz)  &
         &        ,(psi(i),i=1,nsr*nsz)  &
         &        ,nv,(pds(i),i=1,nv),(fds(i),i=1,nv)  &
         &        ,(vlv(i),i=1,nv),(qqv(i),i=1,nv),(prv(i),i=1,nv)  &
         &        ,(xxx(i),i=1,nv)  &
         &        ,btv,ttcu,ttpr,bets,beta,betj &
         &        ,(dummy(i),i=1,6) &
         &        ,(idummy(i),i=1,4) &
         &        ,(dummy(i),i=7,24) &
         &        ,nsu,(rsu(i),i=1,nsu),(zsu(i),i=1,nsu),(csu(i),i=1,nsu) &
         &        ,isep,dsep,rsep,zsep &
         &        ,idummy(5) &
         &        ,idummy(6) &
         &        ,(dummy(i),i=25,29) &
         &        ,(idummy(i),i=7,14) &
         &        ,(dummy(i),i=30,50) &
         &        ,(idummy(i),i=15,19) &
         &        ,(dummy(i),i=51,1000)
!         &        ,dummy(25) &
!         &        ,(idummy(i),i=7,100)
!         &        ,(icp(i),i=1,10),(cp(i),i=1,10)
!         &        ,saxis,raxis,zaxis,ell,trg
!         &        ,nsu,(rsu(i),i=1,nsu),(zsu(i),i=1,nsu),(csu(i),i=1,nsu) &
!         &        ,isep,dsep,rsep,zsep
!         &        ,(ivac(i),i=0,nsfix)  &
!         &        ,(rvac(i),i=0,nsfix),(zvac(i),i=0,nsfix)  &
!         &        ,(cvac(i),i=0,nsfix),(ncoil(i),i=1,nsfix)  &
!         &        ,((rcoil(i,j),i=1,100),j=1,icvdm)  &
!         &        ,((zcoil(i,j),i=1,100),j=1,icvdm)  &
!         &        ,((ccoil(i,j),i=1,100),j=1,icvdm)  &
!         &        ,ilimt,(rlimt(i),i=1,100),(zlimt(i),i=1,100)
    close(ieqrd)
    if( ist > 0 ) then ! error
       write(6,*) 'ist=',ist
       stop '===== dataio/read error ====='
    else if( ist < 0 ) then ! end
       backspace ieqrd
       write(6,"(5x,'dataio/end:backspace')")
       read(*,*)
    end if
!-----------------------------------------------------------------------
    write(6,40)eqfile,ieqrd,jeqrd
40  format(//5x,'diskio:',a19,':','read(',i3,'/',i3,')')

    write(6,'(A,7I10)') 'N:',nsr,nsz,nv,nsu,nsfix,icvdm,ilimt
    write(6,'(A,1P2E12.4)') 'rg: ',rg(1),rg(nsr)
    write(6,'(A,1P2E12.4)') 'zg: ',zg(1),zg(nsz)
    write(6,'(A,1P4E12.4)') 'psi:',psi(1),psi(nsr), &
                                   psi(nsr*(nsz-1)+1),psi(nsr*nsz)
    call pages
    call grd2d(0,rg,zg,psi,nsr,nsr,nsz,'psi',ASPECT=0.D0)
    call pagee
    write(6,'(A,1P2E12.4)') 'pds:',pds(1),pds(nv)
    write(6,'(A,1P2E12.4)') 'fds:',fds(1),fds(nv)
    write(6,'(A,1P2E12.4)') 'vlv:',vlv(1),vlv(nv)
    write(6,'(A,1P2E12.4)') 'qqv:',qqv(1),qqv(nv)
    write(6,'(A,1P2E12.4)') 'prv:',prv(1),prv(nv)
    write(6,'(A,1P2E12.4)') 'xxx:',xxx(1),xxx(nv)
    write(6,'(A,1P6E12.4)') 'btv:',btv,ttcu,ttpr,bets,beta,betj
!    write(6,'(A,5I10)')     'icp:',(icp(i),i=1,5)
!    write(6,'(A,1P5E12.4)') 'cp :',(cp(i),i=1,5)
    write(6,'(A,1P5E12.4)') 'sax:',saxis,raxis,zaxis,ell,trg
    write(6,'(A,1P3E12.4)') 'rsu:',rsu(1),rsu((nsu-1)/2+1),rsu(nsu)
    write(6,'(A,1P3E12.4)') 'zsu:',zsu(1),zsu((nsu-1)/2+1),zsu(nsu)
    write(6,'(A,1P3E12.4)') 'csu:',csu(1),csu((nsu-1)/2+1),csu(nsu)
    write(6,'(A,i10,1P3E12.4)') 'sep:',isep,dsep,rsep,zsep
    write(6,'(5I10)') (idummy(i),i=1,20)
    write(6,'(1P5E12.4)') (dummy(i),i=1,100)
    write(6,*)
    write(6,'(1P5E12.4)') (dummy(i),i=101,200)
    write(6,*)
    write(6,'(1P5E12.4)') (dummy(i),i=201,300)
    
    CALL equ_set_var1(nsr,nsz,nv,nsu,ilimt,btv,saxis,ell,trg)
    CALL equ_set_psi(psi)

   RETURN
 END subroutine eqdsk

!=======================================================================
!             eqrbp               revised on oct.,1986
!=======================================================================
  subroutine eqrbp
!=======================================================================
!     load rbp                                                     jaeri
!=======================================================================
    use equ_params
    implicit none
    integer(4) :: i, ir, iz
    real(8) :: x, y
!-----------------------------------------------------------------------
    i=nrm
    do iz=2,nszm
       i=i+2
       do ir=2,nrm
          i=i+1
          x=dr2i*(psi(i+ 1)-psi(i- 1))
          y=dz2i*(psi(i+nr)-psi(i-nr))
          rbp(i)=sqrt(x*x+y*y)
       end do
    end do
  end subroutine eqrbp

!=======================================================================
!             eqlin               revised        nov.,1986
!=======================================================================
  subroutine eqlin
!=======================================================================
!     line integrals
!=======================================================================
    use bpsd_constants, only : zpi => PI
    use equ_params
    implicit none
    ! intf: num. of division in the direction of lambda for trapped particle fraction
    integer(4), parameter :: intf = 100
    integer(4) :: i,i1,i2,i3,i4,ir,iz,istep,irst,ist,jr,k,ll,lm,lp,n,j
    real(8) :: dpsi,psi0,x,r0,z0,r1,z1,s1,s2,s3,s4,dl &
         &     ,bp0,bl0,ds0,ck0,ss0,vl0,aa0,rr0,bb0,bi0,sh0 &
         &     ,bp1,bl1,ds1,ck1,ss1,vl1,aa1,rr1,bb1,bi1,sh1 &
         &     ,ai0,ai1,br0,br1,bm0,bm1,zzmax,zzmin,rzmax,rzmin &
         &     ,rrmax,rrmin,rmajl,rplal,sdw2
    real(8) :: fintx, hsq, h
    integer(4), dimension(:), allocatable :: nsul
    real(8), dimension(:), allocatable :: bmax,fint,flam,dll,zbl
!=======================================================================
    ieqerr(1)=0
    allocate(bmax(ivdm),fint(0:intf),flam(0:intf))
    allocate(nsul(isrzdm),dll(isrzdm),zbl(isrzdm))
    do j = 0, intf
       flam(j)  = real( j, 8 ) / real( intf, 8 )
    end do
!-----------------------------------------------------------------------
    istep=0
999 nsu=0
    dpsi=saxis/float(nvm)
    do jr=iraxis,nrm
       irst=jr
       ist= (izaxis-1)*nr + jr
       if(sigcu*psi(ist+1) > 0.d0) exit
    end do
!-----------------------------------------------------------------------
    do n=nv,2,-1
!-----------------------------------------------------------------------
!..for double integral to compute trapped particle fraction
       nsul(n)=0
       fint(:) = 0.d0
!-----------------------------------------------------------------------
       psi0=dpsi*dfloat(nv-n)
5      siw(n)=psi0
!-----------------------------------------------------------------------
!..search starting point
       i=ist+1
       ir=irst+1
       iz=izaxis
10   i=i-1
       ir=ir-1
       if(ir < iraxis) go to 800
       if((psi0-psi(i))*(psi0-psi(i+1)) > 0.d0) go to 10
       ist=i
       irst=ir
!-----------------------------------------------------------------------
!..initialize line integrals
       arv(n)=0.d0
       vlv(n)=0.d0
       sdw(n)=0.d0
       ckv(n)=0.d0
       ssv(n)=0.d0
       aav(n)=0.d0
       rrv(n)=0.d0
       bbv(n)=0.d0
       biv(n)=0.d0
       shv(n)=0.d0
       grbm2v(n)=0.d0
       aiv(n)=0.d0
       brv(n)=0.d0
       bmax(n)=0.d0

       x=(psi0-psi(i))/(psi(i+1)-psi(i))
       r1=rg(ir)+x*dr
       z1=zg(iz)
       bp1=(rbp(i)+x*(rbp(i+1)-rbp(i)))/r1
       vl1=r1*r1
       ds1=1.d0/bp1
       ck1=bp1
       ss1=vl1*bp1
       aa1=ds1/vl1
       rr1=vl1/bp1
       bb1=rbv(n)*rbv(n)/(r1*r1)+bp1*bp1
       bi1=ds1/bb1
       bm1=sqrt(bb1)
       br1=bm1*ds1
       bb1=bb1*ds1
       sh1=r1
       ai1=ds1/r1

       rrmax=r1
       rrmin=r1
       zzmax=z1
       zzmin=z1
       rzmax=-1000.d0
       rzmin= 1000.d0

!..trace contour
       k=1
       nsul(n)=1
       if(n == nv) then
          nsu=1
          rsu(1)=r1
          zsu(1)=z1
          csu(1)=bp1*sigcu
       endif
20     r0=r1
       z0=z1
       bp0=bp1
       vl0=vl1
       ds0=ds1
       ck0=ck1
       ss0=ss1
       aa0=aa1
       rr0=rr1
       bi0=bi1
       bb0=bb1
       sh0=sh1
       ai0=ai1
       br0=br1
       bm0=bm1
       i1=i
       i2=i1+1
       i3=i2-nr
       i4=i1-nr
       s1=psi0-psi(i1)
       s2=psi0-psi(i2)
       s3=psi0-psi(i3)
       s4=psi0-psi(i4)
       if(s1*s2 < 0.d0 .and. k /= 1) go to 30
       if(s2*s3 < 0.d0 .and. k /= 2) go to 40
       if(s3*s4 < 0.d0 .and. k /= 3) go to 50
       if(s4*s1 < 0.d0 .and. k /= 4) go to 60
       psi0=psi0+1.d-08*saxis
       go to 5
30     x=s1/(s1-s2)
       r1=rg(ir)+dr*x
       z1=zg(iz)
       bp1=(rbp(i1)+x*(rbp(i2)-rbp(i1)))/r1
       i=i+nr
       iz=iz+1
       k=3
       go to 70
40     x=s2/(s2-s3)
       r1=rg(ir+1)
       z1=zg(iz)-dz*x
       bp1=(rbp(i2)+x*(rbp(i3)-rbp(i2)))/r1
       i=i+1
       ir=ir+1
       k=4
       go to 70
50     x=s4/(s4-s3)
       r1=rg(ir)+dr*x
       z1=zg(iz-1)
       bp1=(rbp(i4)+x*(rbp(i3)-rbp(i4)))/r1
       i=i-nr
       iz=iz-1
       k=1
       go to 70
60     x=s1/(s1-s4)
       r1=rg(ir)
       z1=zg(iz)-dz*x
       bp1=(rbp(i1)+x*(rbp(i4)-rbp(i1)))/r1
       i=i-1
       ir=ir-1
       k=2
!..check
70     if( ir <=  1   ) go to 810
       if( ir >= nrm  ) go to 810
       if( iz <=  1   ) go to 810
       if( iz >= nszm ) go to 810
!..line integrals
       vl1=r1*r1
       ds1=1.d0/bp1
       ck1=bp1
       ss1=vl1*bp1
       aa1=ds1/vl1
       rr1=vl1/bp1
       bb1=rbv(n)*rbv(n)/(r1*r1)+bp1*bp1
       bi1=ds1/bb1
       bm1=sqrt(bb1)
       br1=bm1*ds1
       bb1=bb1*ds1
       sh1=r1
       ai1=ds1/r1
       dl=sqrt((r1-r0)*(r1-r0)+(z1-z0)*(z1-z0))
       arv(n)=arv(n)+(r0+r1)*(z0-z1)*0.5d0
       vlv(n)=vlv(n)+(vl0+vl1)*(z0-z1)*0.5d0
       sdw(n)=sdw(n)+dl*(ds0+ds1)*0.5d0
       ckv(n)=ckv(n)+dl*(ck0+ck1)*0.5d0
       ssv(n)=ssv(n)+dl*(ss0+ss1)*0.5d0
       aav(n)=aav(n)+dl*(aa0+aa1)*0.5d0
       rrv(n)=rrv(n)+dl*(rr0+rr1)*0.5d0
       bbv(n)=bbv(n)+dl*(bb0+bb1)*0.5d0
       biv(n)=biv(n)+dl*(bi0+bi1)*0.5d0
       shv(n)=shv(n)+dl*(sh0+sh1)*0.5d0
       grbm2v(n)=grbm2v(n)+dl*(r0**2*bp0**2*bi0 &
            &                 +r1**2*bp1**2*bi1)*0.5d0
       aiv(n)=aiv(n)+dl*(ai0+ai1)*0.5d0
       brv(n)=brv(n)+dl*(br0+br1)*0.5d0

       dll(nsul(n)) = dl*(ds0+ds1)*0.5d0 ! dlp/Bp
       zbl(nsul(n)) = 0.5d0*(bm0+bm1)    ! B

       bmax(n) = max(bmax(n),zbl(nsul(n)))

!..geometric factors

       if(r1 > rrmax) rrmax=r1 ! rrmax = R_max
       if(r1 < rrmin) rrmin=r1 ! rrmin = R_min
       if(z1 > zzmax) then
          zzmax=z1 ! zzmax = Z_max
          rzmax=r1 ! rzmax = R at Z_max
       end if
       if(z1 < zzmin) then
          zzmin=z1 ! zzmin = Z_min
          rzmin=r1 ! rzmin = R at Z_min
       end if

       nsul(n)=nsul(n)+1
       if(n == nv) then
          nsu=nsu+1
          rsu(nsu)=r1
          zsu(nsu)=z1
          csu(nsu)=bp1*sigcu
       endif
!..
80     if( iudsym == 1 .and. iz > nz) go to 100
       if( i /= ist ) go to 20
!=======================================================================
100    if( iudsym == 1 ) then
          vlv(n)=2.d0*vlv(n)
          arv(n)=2.d0*arv(n)
          sdw(n)=2.d0*sdw(n)
          ckv(n)=2.d0*ckv(n)
          ssv(n)=2.d0*ssv(n)
          aav(n)=2.d0*aav(n)
          rrv(n)=2.d0*rrv(n)
          bbv(n)=2.d0*bbv(n)
          biv(n)=2.d0*biv(n)
          shv(n)=2.d0*shv(n)
          grbm2v(n)=2.d0*grbm2v(n)
          aiv(n)=2.d0*aiv(n)
          brv(n)=2.d0*brv(n)
          if( n == nv ) then
             lp=nsu
             lm=nsu
             do ll=2,nsu-1
                lp=lp+1
                lm=lm-1
                rsu(lp)= rsu(lm)
                zsu(lp)=-zsu(lm)
                csu(lp)= csu(lm)
             end do
             nsu=2*nsu-1
          endif
       endif
!-----
       vlv(n)=zpi*vlv(n)
       sdw(n)=1.d0/(2.d0*zpi*sdw(n))
       ckv(n)=2.d0*zpi*ckv(n)/sdw(n)
       r2b2v(n)=2.d0*zpi*ssv(n)*sdw(n)
       ssv(n)=2.d0*zpi*ssv(n)/sdw(n)
       aav(n)=2.d0*zpi*aav(n)*sdw(n)
       rrv(n)=2.d0*zpi*rrv(n)*sdw(n)
       bbv(n)=2.d0*zpi*bbv(n)*sdw(n)
       biv(n)=2.d0*zpi*biv(n)*sdw(n)
       shv(n)=2.d0*zpi*shv(n)*sdw(n)
       grbm2v(n)=2.d0*zpi*grbm2v(n)*sdw(n)
       aiv(n)=2.d0*zpi*aiv(n)*sdw(n)
       brv(n)=2.d0*zpi*brv(n)*sdw(n)
       sdw(n)=sdw(n)*sigcu

       rmajl=0.5d0*(rrmax+rrmin)
       rplal=0.5d0*(rrmax-rrmin)
       epsv(n) =(rrmax-rrmin)/(rrmax+rrmin)
       elipv(n)=(zzmax-zzmin)/(rrmax-rrmin)
       trigv(n)=(rmajl-0.5d0*(rzmax+rzmin))/rplal

       ! *** trapped particle fraction **********************************
       !
       !   ft = 1 - 0.75 <h^2> int_0^1 flam dflam / <sqrt(1 - flam * h)>
       !
       ! ****************************************************************

       do i = 1, nsul(n)
          h = zbl(i) / bmax(n) ! h
          do j = 0, intf
             if( 1.d0 - flam(j) * h < 0.d0 ) then
                write(6,'(a,1pe26.18,a)') &
                     &        '** WARNING  IN SQRT(DX),DX < 0.0(DX=',1.d0-flam(j)*h,').'
             endif
             fint(j) = fint(j) + sqrt( abs( 1.d0 - flam(j) * h ) ) * dll(i)
          enddo
       enddo

       fintx = 0.d0
       do j = 1, intf
          fintx = fintx + ( flam(j)**2 - flam(j-1)**2 ) / ( fint(j) + fint(j-1) )
       enddo
       hsq     = bbv(n) / bmax(n)**2 ! <h^2>
       fintx   = 0.75d0 * hsq * fintx / (2.d0 * zpi * sdw(n))
       ftv(n) = 1.d0 - fintx

    end do
!..evaluate on the axis
    siw(1)=saxis
    arv(1)=0.d0
    vlv(1)=0.d0
!    sdw(1)=-dpsi*(vlv(3)**2-2.d0*vlv(2)**2)/(vlv(3)*vlv(2)*(vlv(3)-vlv(2)))
    sdw2   = 0.5_8 * ( ( siw(2) - siw(1) ) / vlv(2) &
         &           + ( siw(3) - siw(2) ) / ( vlv(3) - vlv(2) ) )
    sdw(1) = ( vlv(3) * sdw2 - 2.0_8 * vlv(2) / vlv(3) &
         &           * ( siw(3) - siw(1) ) ) / ( vlv(3) - 2.0_8 * vlv(2) )
    ckv(1)=0.d0
    r2b2v(1)=0.d0
    ssv(1)=0.d0
    aav(1)=1.d0/(raxis*raxis)
    rrv(1)=raxis*raxis
    bbv(1)=(rbv(1)/raxis)**2
    biv(1)=1.d0/bbv(1)
    shv(1)=0.d0
    grbm2v(1)=0.d0
    aiv(1)=1.d0/raxis
    brv(1)=rbv(1)/raxis
    nsu = nsu - 1

    epsv(1)=0.d0
    elipv(1)=1.d0
    trigv(1)=0.d0

    ftv(1) =0.d0

    deallocate(bmax,fint,flam,nsul,dll,zbl)
    return
!*
800 stop 'eqlin : no starting point                '
810 stop 'eqlin : out of range                     '
  end subroutine eqlin

end module eqread_mod

