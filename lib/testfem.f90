program testfem

      USE libgrf
      USE libtestfem
      USE fem_calc
      implicit none
!      complex(8),parameter:: ci=(0.d0,1.d0)
      integer:: id,nrmax,npow,nth,nph,nr
      real(8):: rf,angl

      call gsopen

      write(6,*) ' * id=-3: Laplace eq (cylinder): hhg'
      write(6,*) ' * id=-2: Laplace eq (cylinder): hhh'
      write(6,*) ' * id=-1: Laplace eq (slab): hh'
      write(6,*) ' * id= 0: stop'
      write(6,*) ' * id= 1: Maxwell eq (slab) Linear'
      write(6,*) ' * id= 2: Maxwell eq (slab) Linear + edge 0th'
      write(6,*) ' * id= 3: Maxwell eq (slab) Linear + edge 1st'
      write(6,*) ' * id= 4: Maxwell eq (slab) Hermite'
      write(6,*) ' * id= 6: Maxwell eq (cylinder) Linear'
      write(6,*) ' * id= 7: Maxwell eq (cylinder) Linear + edge 0th'
      write(6,*) ' * id= 8: Maxwell eq (cylinder) Linear + edge 1st'
      write(6,*) ' * id= 9: Maxwell eq (cylinder) Hermite'
      write(6,*) ' * id=10: Maxwell eq (cylinder) Hermite+ quadra'
      write(6,*) ' * id=11: Maxwell eq (cylinder) quadra + linear (cont)'
      write(6,*) ' * id=12: Maxwell eq (cylinder) quadra + linear (discont)'
      write(6,*) ' * id=13: Maxwell eq (cylinder) E+- Hermite'
      write(6,*) ' * id=14: Maxwell eq (cylinder) Laplacian E Hermite'
      write(6,*) ' * id=15: Maxwell eq (cylinder) Potential A,phi Hermite+step'
      write(6,*) ' * id=16: Maxwell eq (cylinder) Potential A,phi Hermite+cubic'
      write(6,*) ' * *******************************'

      id=16
      nrmax=11
      npow=1
      nth=0
      nph=0
      rf=1.d0
      angl=0.d0

    1 write(6,'(A,I3,I5,I3,I3,I3,1PE12.4,0PF8.4)') &
              '## input: id,nrmax,npow,nth,nph,rf,angl=', &
                         id,nrmax,npow,nth,nph,rf,angl
      read(5,*,err=1,end=9000) id,nrmax,npow,nth,nph,rf,angl
      if(id.eq.0) goto 9000

!$$$      call fem_exec(id,nrmax,npow,nth,nph,rf,angl)
!$$$      call pages
!$$$      call femgr1dc( 1,rho,cf1,nrmax,'@cf1(rho)@')
!$$$      call femgr1dc( 2,rho,cf2,nrmax,'@cf2(rho)@')
!$$$      call femgr1dc( 3,rho,cf3,nrmax,'@cf3(rho)@')
!$$$      call pagee
!$$$      do nr=1,nrmax
!$$$         write(6,'(I5,1P7E10.2)') nr,rho(nr),cf1(nr),cf2(nr),cf3(nr)
!$$$      enddo

      call pages
      call fem_exec(id,nrmax,npow,nth,nph,rf,angl)
      call femgr1dc( 5,rho,cf1,nrmax,'@cf1(rho)@')
      call femgr1dc( 6,rho,cf2,nrmax,'@cf2(rho)@')
      call femgr1dc( 7,rho,cf3,nrmax,'@cf3(rho)@')
      call fem_exec(id,nrmax,npow,nth+1,nph,rf,angl)
      call femgr1dc( 8,rho,cf1,nrmax,'@cf1(rho)@')
      call femgr1dc( 9,rho,cf2,nrmax,'@cf2(rho)@')
      call femgr1dc(10,rho,cf3,nrmax,'@cf3(rho)@')
      call fem_exec(id,nrmax,npow,nth+2,nph,rf,angl)
      call femgr1dc(11,rho,cf1,nrmax,'@cf1(rho)@')
      call femgr1dc(12,rho,cf2,nrmax,'@cf2(rho)@')
      call femgr1dc(13,rho,cf3,nrmax,'@cf3(rho)@')
      call pagee

      goto 1

 9000 call gsclos
      stop

      contains

      subroutine fem_exec(id,nrmax,npow,nth,nph,rf,angl)

      implicit none
      integer,intent(in):: id,nrmax,npow,nth,nph
      real(8),intent(in):: rf,angl
      real(8):: rho0,rkth,rkph
      integer:: mw,ml,nr,ierr

      select case(id)
         case(-1)
            call fem_calc_l1(nrmax,npow)
         case(-2)
            call fem_calc_l2(nrmax,npow)
         case(-3)
            call fem_calc_l3(nrmax,npow)
         case(1)
            call fem_calc_x1(nrmax,npow,nth,nph,rf,angl)
         case(2)
            call fem_calc_x2(nrmax,npow,nth,nph,rf,angl)
         case(3)
            call fem_calc_x3(nrmax,npow,nth,nph,rf,angl)
         case(4)
            call fem_calc_x4(nrmax,npow,nth,nph,rf,angl)
         case(6)
            call fem_calc_r1(nrmax,npow,nth,nph,rf,angl)
         case(7)
            call fem_calc_r2(nrmax,npow,nth,nph,rf,angl)
         case(8)
            call fem_calc_r3(nrmax,npow,nth,nph,rf,angl)
         case(9)
            call fem_calc_r4(nrmax,npow,nth,nph,rf,angl)
         case(10)
            call fem_calc_r5(nrmax,npow,nth,nph,rf,angl)
         case(11)
            call fem_calc_r6(nrmax,npow,nth,nph,rf,angl)
         case(12)
            call fem_calc_r7(nrmax,npow,nth,nph,rf,angl)
         case(13)
            call fem_calc_r8(nrmax,npow,nth,nph,rf,angl)
         case(14)
            call fem_calc_r9(nrmax,npow,nth,nph,rf,angl)
         case(15)
            call fem_calc_ra(nrmax,npow,nth,nph,rf,angl)
         case(16)
            call fem_calc_16(nrmax,npow,nth,nph,rf,angl)
      end select

!      do ml=1,mlmax
!         write(16,'(1P6E12.4)') (fma(mw,ml),mw=1,mwmax)
!      enddo
!      write(16,'(1P6E12.4)') (fvb(ml),ml=1,mlmax)

      do ml=1,mlmax
         fvx(ml)=fvb(ml)
      enddo
      call bandcd(fma,fvx,mlmax,mwmax,mwmax,ierr)
      if(ierr.ne.0) write(6,*) '# ierr= ',ierr

      if(id.eq.-1.or.id.eq.-2.or.id.eq.-3) then
         do nr=1,nrmax
            cf1(nr)=fvx(2*(nr-1)+1)
            cf2(nr)=fvx(2*(nr-1)+2)
            cf3(nr)=0.d0
         enddo

      else if(id.eq.1.or.id.eq.6) then
         do nr=1,nrmax
            cf1(nr)=fvx(3*(nr-1)+1)
            cf2(nr)=fvx(3*(nr-1)+2)
            cf3(nr)=fvx(3*(nr-1)+3)
         enddo

      else if(id.eq.2.or.id.eq.7) then
         do nr=1,nrmax
            cf1(nr)=fvx(3*(nr-1)+1)
            cf2(nr)=fvx(3*(nr-1)+2)
            cf3(nr)=fvx(3*(nr-1)+3)
         enddo

      else if(id.eq.3.or.id.eq.8) then
         do nr=1,nrmax
            cf1(nr)=fvx(4*(nr-1)+1)
            cf2(nr)=fvx(4*(nr-1)+3)
            cf3(nr)=fvx(4*(nr-1)+4)
         enddo

      else if(id.eq.4.or.id.eq.9.or.id.eq.14) then
         do nr=1,nrmax
            cf1(nr)=fvx(6*(nr-1)+1)
            cf2(nr)=fvx(6*(nr-1)+3)
            cf3(nr)=fvx(6*(nr-1)+5)
         enddo

      else if(id.eq.10) then
         do nr=1,nrmax
            cf1(nr)=fvx(6*(nr-1)+1)
            cf2(nr)=fvx(6*(nr-1)+3)
            cf3(nr)=fvx(6*(nr-1)+5)
         enddo
         cf1(nrmax)=cf1(nrmax-1)

      else if(id.eq.11) then
         do nr=1,nrmax
            cf1(nr)=fvx(5*(nr-1)+1)
            cf2(nr)=fvx(5*(nr-1)+2)
            cf3(nr)=fvx(5*(nr-1)+4)
         enddo
         cf1(nrmax)=cf1(nrmax-1)
      else if(id.eq.12) then
         do nr=1,nrmax-1
            cf1(nr)=fvx(6*(nr-1)+3)
            cf2(nr)=fvx(6*(nr-1)+1)
            cf3(nr)=fvx(6*(nr-1)+2)
         enddo
         nr=nrmax
         cf1(nr)=fvx(6*(nr-1)  )
         cf2(nr)=fvx(6*(nr-1)+1)
         cf3(nr)=fvx(6*(nr-1)+2)
      else if(id.eq.13) then
         do nr=1,nrmax
            cf1(nr)=0.5D0*(fvx(6*(nr-1)+1)+fvx(6*(nr-1)+3))
            cf2(nr)=0.5D0*(fvx(6*(nr-1)+1)-fvx(6*(nr-1)+3))/CI
            cf3(nr)=fvx(6*(nr-1)+5)
         enddo
      else if(id.eq.15 .or. id.eq.16) then
         do nr=1,nrmax
            IF(nr.EQ.1) THEN
               IF(ABS(nth).eq.1) THEN
                  cf2(nr)=ci*fvx(8*(nr-1)+3)+ci*nth*fvx(8*(nr-1)+8)
               ELSE
                  cf2(nr)=0.D0
               ENDIF
            ELSE
               rho0=rho(nr)
               rkth=nth/rho0
               cf2(nr)=ci*fvx(8*(nr-1)+3)-ci*rkth*fvx(8*(nr-1)+7)
            ENDIF
            rkph=nph
            cf1(nr)=ci*fvx(8*(nr-1)+1)-fvx(8*(nr-1)+8)
            cf3(nr)=ci*fvx(8*(nr-1)+5)-ci*rkph*fvx(8*(nr-1)+7)

            cf1(nr)=ci*fvx(8*(nr-1)+1)
            cf2(nr)=ci*fvx(8*(nr-1)+3)
            cf3(nr)=ci*fvx(8*(nr-1)+7)
         enddo

      endif
      return
      end subroutine fem_exec

      include 'testfem-sub.f90'

      end program testfem
