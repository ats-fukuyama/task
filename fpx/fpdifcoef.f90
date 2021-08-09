module fpdiffusivity

  use fpcomm

  implicit none
    real(rkind),dimension(:,:),allocatable::gammap,deff,dndr
contains
!-------------------------------------------------------------------------------

  subroutine fp_deff

    implicit none
    real(rkind),dimension(nrmax,nsamax)::dndt
    real(rkind)::vprimem,vprimeg,gsum
    integer:: nsa,nr

    allocate(gammap(nrmax,nsamax))
    allocate(dndr(nrmax,nsamax))
    allocate(deff(nrmax,nsamax))

    dndt(:,:)=(rns(:,:)-rns_prev(:,:))/delt

    do nsa=1,nsamax
      dndr(1,nsa)=0.d0
      gsum=0.d0
      gammap(1,nsa)=0.d0
      do nr=1,nrmax-1
        vprimem=4.d0*pi**2*rr*(nr-0.5d0)/nrmax*ra
        vprimeg=4.d0*pi**2*rr*nr/nrmax*ra
        dndr(nr+1,nsa)=(rns(nr+1,nsa)-rns(nr,nsa))/delr/ra
        gsum=gsum+(-dndt(nr,nsa)+tps(nr,nsa))*vprimem*delr*ra
        gammap(nr+1,nsa)=gsum/vprimeg
      enddo
    enddo

    do nsa=1,nsamax
       do nr=2,nrmax
         if(dndr(nr,nsa)==0.d0) then
            deff(nr,nsa)=0.d0
          else
            deff(nr,nsa)=-gammap(nr,nsa)/dndr(nr,nsa)
          endif
       enddo
       deff(1,nsa)=(4.d0*deff(2,nsa)-deff(3,nsa))/3.d0
    enddo

  end subroutine fp_deff

!-------------------------------------------------------------------------------

  subroutine fp_chieff
    use libgrf
    use fpwrite
    implicit none
    integer::nr,nsa,t=0
    REAL(rkind)::vprimem,vprimeg,qsum,mj2kev
    REAL(rkind),dimension(nrmax,nsamax)::temps,drhotdt,ps,hf,dtdr,chieff,keff,dpdr,inPower

    mj2kev=1.d-17/aee    !MJ to 10^20keV
    temps(:,:)=rws(:,:)/(1.5d0*rns(:,:))
    drhotdt(:,:)=(rws(:,:)-rws_prev(:,:))/delt
    inPower(:,:)=rspb(:,:)+rspf(:,:)+rsps(:,:)+rspl(:,:)+rsps_cx(:,:)
    ps(:,:)=inPower(:,:)+rpcs(:,:)+rpss(:,:)! &
            !+rpws(:,:)+rpes(:,:)+rlhs(:,:)+rfws(:,:)+recs(:,:)+rpss(:,:) &
            !+rpls(:,:)+rws123(:,:)+rwrs(:,:)+rwms(:,:)

    do nsa=1,nsamax
      dtdr(1,nsa)=0.d0
      qsum=0.d0
      hf(1,nsa)=0.d0
      do nr=1,nrmax-1
        dtdr(nr+1,nsa)=(temps(nr+1,nsa)-temps(nr,nsa))/delr/ra
        dpdr(nr+1,nsa)=(rws(nr+1,nsa)-rws(nr,nsa))/delr/ra
        vprimem=4.d0*pi**2*rr*(nr-0.5d0)/nrmax*ra
        vprimeg=4.d0*pi**2*rr*nr/nrmax*ra
        qsum=qsum+(-drhotdt(nr,nsa)+ps(nr,nsa))*vprimem*delr*ra
        hf(nr+1,nsa)=qsum/vprimeg
      enddo
    enddo
    do nsa=1,nsamax
      do nr=2,nrmax
        if (dtdr(nr,nsa)==0.d0.or.(rns(nr,nsa)-rns(nr-1,nsa))==0.d0.or.dpdr(nr,nsa)==0.d0) then
          keff(nr,nsa)=0.d0               ! effective pressure diffusivity
          chieff(nr,nsa)=0.d0             ! effective temperature diffusivity
        else
          keff(nr,nsa)=-hf(nr,nsa)/(dpdr(nr,nsa))
          chieff(nr,nsa)=-hf(nr,nsa)/((rns(nr,nsa)+rns(nr-1,nsa))/2.d0*dtdr(nr,nsa))
          ! chieff(nr,nsa)=-(hf(nr,nsa)-1.5d0*(temps(nr,nsa)+temps(nr-1,nsa))/2.d0*gammap(nr,nsa)) &
          !                 /((rns(nr,nsa)+rns(nr-1,nsa))/2.d0*dtdr(nr,nsa))
        end if
      enddo
      keff(1,nsa)=(4.d0*keff(2,nsa)-keff(3,nsa))/3.d0
      chieff(1,nsa)=(4.d0*chieff(2,nsa)-chieff(3,nsa))/3.d0
    enddo

    open(1010,file='input.csv')
      do nr=1,nrmax
        write(1010,*)tps(nr,1),',',tps(nr,2),',',inPower(nr,1),',',inPower(nr,2)
      end do
    close(1010)

    call fpcsv(deff,chieff,keff,gammap,hf,temps*mj2kev,ps,rws)

    deallocate(gammap,deff,dndr)

  end subroutine fp_chieff
!-------------------------------------------------------------------------------
end module fpdiffusivity
