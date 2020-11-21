subroutine fowweight(NSA,IERR) ! proposed by Chang and Cooper [30] in Karney

    use fpcomm
    use fowcomm

    implicit none
    integer:: nsa, np, nth, nr, ntha, nthb, ns
    integer:: ierr
    real(8):: dfdth, fvel, dfdp
    
    integer :: nra, nrb, npa, npb
    double precision :: dfdr, delps, factab_p, factab_r, factab_t

    !     +++++ calculation of weigthing (including off-diagonal terms) +++++

    real(8)::epswt=1.d-70

    ns = ns_nsa(nsa)

    do nr = nrstart, nrend
      delps = psimg(nr+1)-psimg(nr)
      do np = npstart, npendwg
        do nth = 1, nthmax
            dfdth = 0.d0
            dfdr  = 0.d0
            if ( np /= 1 ) then
              ntha = min(nth+1,nthmax)
              nthb = max(nth-1,1)
              nra  = min(nr+1,nrmax)
              nrb  = max(nr-1,1)

              if ( np == npmax+1 ) then
                if ( abs(f(nth,np-1,nr)) > epswt ) then
                  dfdth = (f(ntha,np-1,nr)-f(nthb,np-1,nr)) &
                          /(2.d0*pg(np,ns)*delthm_pg(nth,np,nr,nsa)*f(nth,np-1,nr))

                  dfdr  = (f(nth,np-1,nra)-f(nth,np-1,nrb)) &
                          /(2.d0*delps*f(nth,np-1,nr))
                end if
              else

                if ( abs(f(nth,np-1,nr)) > epswt ) then
                  if ( abs(f(nth,np,nr)) > epswt ) then
                    dfdth = (f(ntha,np-1,nr)-f(nthb,np-1,nr)) &
                            /(4.d0*pg(np,ns)*delthm(nth,np,nr,nsa)*f(nth,np-1,nr)) &
                          + (f(ntha,np  ,nr)-f(nthb,np  ,nr)) &
                            /(4.d0*pg(np,ns)*delthm(nth,np,nr,nsa)*f(nth,np  ,nr))

                    dfdr  = (f(nth,np-1,nra)-f(nthb,np-1,nrb)) &
                            /(4.d0*delps*f(nth,np-1,nr)) &
                          + (f(ntha,np  ,nra)-f(nthb,np  ,nrb)) &
                            /(4.d0*delps*f(nth,np  ,nr))
                  else
                    dfdth = (f(ntha,np-1,nr)-f(nthb,np-1,nr)) &
                          /(2.d0*pg(np,ns)*delthm(nth,np,nr,nsa)*f(nth,np-1,nr)) 

                    dfdr  = (f(nth,np-1,nra)-f(nth,np-1,nrb)) &
                          /(2.d0*delps*f(nth,np-1,nr)) 
                  end if
                else
                  if ( abs(f(nth,np,nr)) > epswt ) then
                    dfdth = (f(ntha,np  ,nr)-f(nthb,np  ,nr)) &
                          /(2.d0*pg(np,ns)*delthm(nth,np,nr,nsa)*f(nth,np  ,nr))

                    dfdr  = (f(nth,np  ,nra)-f(nth,np  ,nrb)) &
                          /(2.d0*delps*f(nth,np  ,nr))
                  end if
                end if
              end if
            end if
            fvel = Fppfow(nth,np,nr,nsa)-Dptfow(nth,np,nr,nsa)*dfdth-Dprfow(nth,np,nr,nsa)*dfdr
            weighp(nth,np,nr,nsa) = fowwegh(-delp(ns)*fvel,dpp(nth,np,nr,nsa))
        end do
      end do
    end do

    do nr = nrstart, nrend
      delps = psimg(nr+1)-psimg(nr)
      do np = npstart, npend
        do nth = 1, nthmax+1
          dfdp = 0.d0          
          npa = min(npmax,np+1)
          npb = max(1, np-1)
          factab_p = dble(npa-npb)

          dfdr = 0.d0
          nra = min(nrmax,nr+1)
          nrb = max(1,nr-1)
          factab_r = dble(nra-nrb)
          
          if ( nth == 1 ) then
            dfdp = (f(nth,npa,nr)-f(nth,npb,nr))/(factab_p*delp(ns)*f(nth,np,nr))
          else if ( nth == nthmax+1 ) then
            dfdp = (f(nthmax,npa,nr)-f(nthmax,npb,nr))/(factab_p*delp(ns)*f(nthmax,np,nr))
            dfdr = (f(nthmax,np,nra)-f(nthmax,np,nrb))/(factab_r*delps*f(nthmax,np,nr))
          else
            dfdp = (
              (f(nth-1,npa,nr)-f(nth-1,npb,nr))/(factab_p*delp(ns)*f(nth-1,np,nr))+&
              (f(nth  ,npa,nr)-f(nth  ,npb,nr))/(factab_p*delp(ns)*f(nth  ,np,nr)) &
            )/2.d0

            dfdr = (
              (f(nth-1,np,nra)-f(nth-1,np,nrb))/(facrab_r*delps*f(nth-1,np,nr))+&
              (f(nth  ,np,nra)-f(nth  ,np,nrb))/(facrab_r*delps*f(nth  ,np,nr)) &
            )
          end if
          fvel = Ftpfow(nth,np,nr,nsa)-Dtpfow(nth,np,nr,nsa)*dfdp-Dtrfow(nth,np,nr,nsa)*dfdr
          weight(nth,np,nr,nsa) = fowwegh(-delthm(nth,np,nr,nsa)*fvel,Dttfow(nth,np,nr,nsa))
        end do
      end do
    end do

    do nr = nrstart, nrendwg
      delps = psimg(nr+1)-psimg(nr)
      do np = npstart, npend
        do nth = 1, nthmax
          dfdp = 0.d0          
          npa = min(npmax,np+1)
          npb = max(1, np-1)
          factab_p = dble(npa-npb)

          dfdth = 0.d0
          ntha = min(nth+1,nthmax)
          nthb = max(nth-1,1)
          factab_t = dble(ntha-nthb)

          if ( nr == 1 ) then
            dfdp = (f(nth,npa,1)-f(nth,npb,1))/(factab_p*delp(ns)*f(nth,np,1))
            dfdth= (f(ntha,np,1)-f(nthb,np,1))/(factab_t*delthm(nth,np,nr,nsa)*f(nth,np,1))
          else if ( nr == nrmax+1 ) then
            dfdp = (f(nth,npa,nrmax)-f(nth,npb,nrmax))/(factab_p*delp(ns)*f(nth,np,nrmax))
            dfdth= (f(ntha,np,nrmax)-f(nthb,np,nrmax))/(factab_t*delthm(nth,np,nr,nsa)*f(nth,np,nrmax))
          else
            dfdp = (
              (f(nth,npa,nr-1)-f(nth,npb,nr-1))/(factab_p*delp(ns)*f(nth,np,nr-1))+&
              (f(nth,npa,nr  )-f(nth,npb,nr  ))/(factab_p*delp(ns)*f(nth,np,nr  )) &
            )/2.d0

            dfdth= (
              (f(ntha,np,nr-1)-f(nthb,np,nr-1))/(factab_t*delthm(nth,np,nr-1,nsa)*f(nth,np,nr))+&
              (f(ntha,np,nr  )-f(nthb,np,nr  ))/(factab_t*delthm(nth,np,nr  ,nsa)*f(nth,np,nr)) &
            )/2.d0
          end if
          fvel = Frrfow(nth,np,nr,nsa)-Drpfow(nth,np,nr,nsa)*dfdp-Drtfow(nth,np,nr,nsa)*dfdth
          weighr(nth,np,nr,nsa) = fowwegh(-delps*fvel,Drrfow(nth,np,nr,nsa))
        end do
      end do
    end do

    return
  end subroutine fowweight
