      module eqinit_mod
      character(len=80) knamcoil
      public
      contains
c     
c=======================================================================
      subroutine eqinit
c=======================================================================
c     DEFAULT VALUES FOR JT60
c=======================================================================
      USE plcomm,ONLY: KNAMEQ,MODEFW,pl_allocate_ns
      use aaa_mod
      use par_mod
      use geo_mod
      use equ_mod
      use vac_mod
      use eqv_mod
      use cnt_mod
      implicit none
c
c-----CURRENT PROFILE
c
      integer   i
c
      integer   icp_l(10)
      real*8    cp_l(10)
      data icp_l  /  1,  1,  0,  0,  0,  0,  0,  0,  0,  0/
      data cp_l   /0.1,1.0,1.0, 0., 0., 0., 0., 0., 0., 0.05/
c
c-----PLASMA SHAPE CONTROLL
c
      integer ivac_l(0:icvdm)
      real*8  cvac_l(0:icvdm)
      data ivac_l /  0,  1,  2, -1,-1,-1,-1,icvmm6*-1/
      data cvac_l / 0., icvdm*0./
c
c-----OUTPUT
      integer ieqout_l(10)
      data ieqout_l/5,9*0/
c
c=======================================================================
c     JT-60 UPGRADE COIL DATA
c=======================================================================
c
      integer  ncoil_l(0:icvdm)
      real*8   rcoil_l(100,icvdm),zcoil_l(100,icvdm),ccoil_l(100,icvdm)
c
      data ncoil_l/0,70,78,26,10, 0, 0,icvmm6*0/
c
c********   RCOIL  ********
c
c*** OH-COIL
      data(rcoil_l(i,1),i=1,100)/
     * 5.0605, 4.8708, 4.1943, 3.7172, 3.3038, 3.3038, 2.8903,
     * 2.8735, 2.4803, 2.4803, 2.4775, 2.2353, 2.2353, 2.0803,
     * 2.0803, 2.0803, 2.0803, 1.8903, 1.8903, 1.8903, 1.8903,
     * 1.8903, 1.7803, 1.7803, 1.7803, 1.7803, 1.7803, 1.7803,
     * 1.7053, 1.7053, 1.7053, 1.7053, 1.7053, 1.6903, 1.6903,
     * 1.6903, 1.6903, 1.6903, 1.6903, 1.7013, 1.7013, 1.7013,
     * 1.7013, 1.7013, 1.7803, 1.7803, 1.7803, 1.7803, 1.7803,
     * 1.7803, 1.9543, 1.9543, 1.9543, 1.9543, 1.9543, 2.2053,
     * 2.2053, 2.2053, 2.2053, 2.2053, 2.4793, 2.4793, 2.4793,
     * 2.4793, 3.8213, 3.8213, 3.8213, 3.8213, 3.8063, 5.0430,
     * 30*0.0/
c*** V-COIL
      data(rcoil_l(i,2),i=1,100)/
     * 4.7663, 4.7663, 4.7633, 4.7633, 4.7633, 4.7633, 4.7633,
     * 4.7663, 4.7663, 4.7633, 4.7633, 4.8943, 4.8943, 4.9088,
     * 4.9088, 4.9088, 4.9088, 4.9088, 4.9088, 4.9088, 4.9088,
     * 4.9088, 1.9493, 1.9493, 1.9493, 1.9493, 1.9493, 1.9493,
     * 1.9493, 1.9493, 1.9493, 1.9493, 1.9493, 1.9493, 4.7198,
     * 4.7198, 4.7198, 4.7198, 4.7198, 4.7198, 4.7198, 4.7198,
     * 4.7198, 4.7198, 4.7198, 4.7198, 4.8763, 4.8763, 4.8763,
     * 4.8763, 4.8763, 4.8763, 4.8763, 4.8763, 4.8763, 4.8763,
     * 4.8605, 4.8525,
     * 4.0313, 4.0313, 4.0313, 4.0313, 4.0380, 3.8870, 3.8925,
     * 3.8925, 3.8925, 3.8975, 4.0763, 4.0763, 4.0763, 4.0763,
     * 4.0763, 4.1853, 4.1838, 4.1838, 4.1838, 4.1830, 22*0.0/
c*** H-COIL
      data(rcoil_l(i,3),i=1,100)/
     * 4.6168, 4.2288, 4.2288, 3.7313, 3.0360, 3.0800, 3.1240,
     * 2.4058, 2.2783, 2.2783, 2.1108, 2.1108, 2.0433, 2.0033,
     * 2.2188, 2.2188, 2.2188, 2.2188, 2.4100, 2.4250, 2.4250,
     * 3.8020, 3.7880, 3.9080, 3.9080, 4.5665, 74*0.0/
c***  D-COIL
      data(rcoil_l(i,4),i=1,100)/
     * 2.9008, 2.9008, 3.0758, 3.0758, 3.0758, 3.1943, 3.1943,
     * 3.1943, 3.3593, 3.3693, 90*0.0/
c********* ZCOIL ************
c*** OH-COIL
      data(zcoil_l(i,1),i=1,100)/
     * 0.5755, 0.9690, 1.6150, 1.7990, 1.8620, 1.8240, 1.8105,
     * 1.7595, 1.6224, 1.5650, 1.5077, 1.4605, 1.3995, 1.2995,
     * 1.2465, 1.1935, 1.1405, 1.0228, 0.9764, 0.9300, 0.8836,
     * 0.8372, 0.7334, 0.6880, 0.6427, 0.5973, 0.5520, 0.5067,
     * 0.4068, 0.3584, 0.3100, 0.2616, 0.2132, 0.1134, 0.0680,
     * 0.0227,-0.0227,-0.0680,-0.1134,-0.2374,-0.2858,-0.3342,
     *-0.3826,-0.4310,-0.5068,-0.5521,-0.5974,-0.6427,-0.6880,
     *-0.7333,-0.8372,-0.8836,-0.9300,-0.9764,-1.0228,-1.1612,
     *-1.2076,-1.2540,-1.3004,-1.3468,-1.4230,-1.4810,-1.5390,
     *-1.5970,-1.6280,-1.6665,-1.7050,-1.7435,-1.7820,-0.6070,
     * 30*0.0/
c*** V-COIL
      data(zcoil_l(i,2),i=1,100)/
     * 0.8790, 0.8390, 0.7990, 0.7590, 0.7190, 0.6790, 0.6390,
     * 0.5990, 0.5590, 0.5190, 0.4790, 0.8757, 0.8300, 0.7874,
     * 0.7489, 0.7104, 0.6719, 0.6334, 0.5949, 0.5564, 0.5179,
     * 0.4794, 0.3675, 0.3205, 0.2735, 0.2265, 0.1795, 0.1325,
     *-0.1325,-0.1795,-0.2265,-0.2735,-0.3205,-0.3675,-0.4770,
     *-0.5130,-0.5490,-0.5850,-0.6210,-0.6570,-0.6930,-0.7290,
     *-0.7650,-0.8010,-0.8370,-0.8730,-0.4775,-0.5125,-0.5475,
     *-0.5825,-0.6175,-0.6525,-0.6875,-0.7225,-0.7575,-0.7925,
     *-0.8305,-0.8715,
     * 1.6918, 1.6534, 1.6150, 1.5766, 1.5383, 1.7618, 1.7234,
     * 1.6850, 1.6466, 1.6082,-1.4562,-1.5006,-1.5450,-1.5894,
     *-1.6338,-1.4562,-1.5006,-1.5450,-1.5894,-1.6338, 22*0.0/
c*** H-COIL
      data(zcoil_l(i,3),i=1,100)/
     * 0.7225, 1.5330, 1.4970, 1.7085, 1.8268, 1.8268, 1.8268,
     * 1.4150, 1.2855, 1.2345, 1.0255, 0.9745, 0.6600,-0.6600,
     *-0.9095,-0.9565,-1.0035,-1.0505,-1.2054,-1.2540,-1.3185,
     *-1.4800,-1.5490,-1.4775,-1.5775,-0.7340, 74*0.0/
c***  D-COIL
      data(zcoil_l(i,4),i=1,100)/
     *-1.7350,-1.7950,-1.7310,-1.7850,-1.8390,-1.7310,
     *-1.7850,-1.8390,-1.7445,-1.8155, 90*0.0/
c********* CCOIL ***********
c*** OH-COIL
      data(ccoil_l(i,1),i=1,100)/
     * 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
     * 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
     * 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
     * 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
     * 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
     * 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
     * 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
     * 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
     * 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
     * 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
     * 30*0.0/
c*** V-COIL
      data(ccoil_l(i,2),i=1,100)/
     *-1.0000,-1.0000,-1.0000,-1.0000,-1.0000,-1.0000,-1.0000,
     *-1.0000,-1.0000,-1.0000,-1.0000,-1.0000,-1.0000,-1.0000,
     *-1.0000,-1.0000,-1.0000,-1.0000,-1.0000,-1.0000,-1.0000,
     *-1.0000,-1.0000,-1.0000,-1.0000,-1.0000,-1.0000,-1.0000,
     *-1.0000,-1.0000,-1.0000,-1.0000,-1.0000,-1.0000,-1.0000,
     *-1.0000,-1.0000,-1.0000,-1.0000,-1.0000,-1.0000,-1.0000,
     *-1.0000,-1.0000,-1.0000,-1.0000,-1.0000,-1.0000,-1.0000,
     *-1.0000,-1.0000,-1.0000,-1.0000,-1.0000,-1.0000,-1.0000,
     *-1.0000,-1.0000,
     *-1.0000,-1.0000,-1.0000,-1.0000,-1.0000,-1.0000,-1.0000,
     *-1.0000,-1.0000,-1.0000,-1.0000,-1.0000,-1.0000,-1.0000,
     *-1.0000,-1.0000,-1.0000,-1.0000,-1.0000,-1.0000, 22*0.0/
c*** H-COIL
      data(ccoil_l(i,3),i=1,100)/
     * 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
     * 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,-1.0000,
     *-1.0000,-1.0000,-1.0000,-1.0000,-1.0000,-1.0000,-1.0000,
     *-1.0000,-1.0000,-1.0000,-1.0000,-1.0000, 74*0.0/
c***  D-COIL
      data(ccoil_l(i,4),i=1,100)/
     * 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
     * 1.0000,-1.0000,-1.0000, 90*0.0/
c
c=======================================================================
c     JT-60 LIMITER DATA
c=======================================================================
c
      real*8 rlimt_l(100),zlimt_l(100)
c
      data rlimt_l
     >          / 4.548 , 4.537 , 4.504 , 4.377 , 4.271 , 4.105 ,
     >            3.785 , 3.500 , 3.210 , 2.926 , 2.729 , 2.562 ,
     >            2.428 , 2.335 , 2.281 , 2.262 , 2.250 , 2.273 ,
     >            2.324 , 2.408 , 2.524 , 2.685 , 2.857 , 3.039 ,
     >            3.323 , 3.402 , 3.795 , 4.083 , 4.321 , 4.480 ,
     >            70*0./
      data zlimt_l
     >          / 0.0   , 0.160 , 0.335 , 0.795 , 1.032 , 1.220 ,
     >            1.425 , 1.570 , 1.627 , 1.550 , 1.420 , 1.255 ,
     >            1.061 , 0.845 , 0.581 , 0.291 , 0.0   ,-0.257 ,
     >           -0.510 ,-0.754 ,-0.985 ,-1.211 ,-1.433 ,-1.518 ,
     >           -1.472 ,-1.441 ,-1.219 ,-0.955 ,-0.643 ,-0.335 ,
     >            70*0./
c=======================================================================
c
      device='JT60U'   ! device name
      devnam='JT60U'   ! default device name
      title=' '        ! job title
c
      nr    =irdm      ! number of grid in R direction
      nz    =izdm      ! number of grid in Z direction
      nv    =ivdm      ! number of grid in psi direction
c
      btv   = 13.5     ! toroidal magnetic field multiplied by R: R[m]*B[T]
      tcur  = 2.0      ! plasma current [MA]
c
c-----CURRENT PROFILE  ! see subroutine eqequ and eqpdf
c
      do i=1,10
         icp(i)=icp_l(i)  
         cp(i) =cp_l(i)
      enddo
c----------------------
c
      rwmn  = 1.90     ! minimum of R in calculation box
      rwmx  = 4.20     ! maximum of R in calculation box
      zwmx  = 1.20     ! maximum of Z in calculation box
c
      iudsym= 0        ! updown symmetry
      ivp   = 0        ! ???
c
c----- plasma shape -----     
c
      msfx  = 9        ! number of marker points
c     rvac             ! position R of marker points
c     zvac             ! position Z of marker points
c
c----- optional shape parameters -----
c
      rmaj  = 3.05     ! major radius [m]
      rpla  = 0.9      ! minor radius [m]
      zpla  = 0.0      ! posiion Z of magnetic axis [m]
      elip  = 1.0      ! ellipticity
      trig  = 0.0      ! triangularity
      yh    = 0.5      ! half height
      yd    = 0.995    ! edge height
      elipup= 1.0
      trigup= 0.0
c
      qaxi  = 1.0
      qsur  = 3.0
      xxli  = 1.0
c
c-----PARAMETERS OF EQUILIBRIUM SOLVER
c
      msetup = 20      ! maximum iteration count of initial equilibrium 
      esetup = 1.d-03  ! covergence criterion of initial equilibrium 
c
      ieqmax = 20      ! maximum iteration count of FCT equilibrium 
      eeqmax = 1.d-03  ! covergence criterion of FCT equilibrium 
c
      iodmax = 10      ! maximum iteration count of ODE calculation
      eodmax = 1.d-05  ! covergence criterion of ODE calculation
c
      iadmax =  1      ! maximum iteration count of vaccum adjustment
      eadmax = 1.d-04  ! covergence criterion of vaccum adjustment
c
      bavmax = 0.8     ! FCT convergence parameter (see eqfct and eqrcu)
      bavmin = 0.2     ! FCT convergence parameter (see eqfct and eqrcu)
c
c-----number of points on plasma surface
c
      nsumax = ivdm
c
c-----coil parameters
c
      do i=0,icvdm
         ivac(i)=ivac_l(i)    ! mode of coil setting : 
                              !         negative: fixed
                              !         positive: adjustable
         cvac(i)=cvac_l(i)    ! coil current [AT]
         ncoil(i)=ncoil_l(i)  ! coil turn
         cvacst(i)=0.d0       ! coil current (initial value) [AT]
         cvacwg(i)=0.d0       ! weight of coil current
         cvact(i)=0.d0        ! time derivative of coil current
      enddo
c
c----- separatrix option
c
      isep   = 2
      dsep   = 0.001
c
c-----
c
      ivtab  = 0
      ivgrp  = 0
c
c-----
c
      rloop  = 2.786
      zloop  = 0.952
c
c----- limiter option
c
      ilimt = 0
      ivbnd = 0
c
c----- output control
c
      do i=1,10
         ieqout(i)=ieqout_l(i)
      enddo
c
c----- Initialize parameters -----
c
      isvac  = 0
      ieqrd  = 0
      jeqrd  = 0
      ieqwt  = 0
      jeqwt  = 0
c
c----- coil position
c
      do i=1,100
         rcoil(i,1)=rcoil_l(i,1)
         rcoil(i,2)=rcoil_l(i,2)
         rcoil(i,3)=rcoil_l(i,3)
         rcoil(i,4)=rcoil_l(i,4)
         zcoil(i,1)=zcoil_l(i,1)
         zcoil(i,2)=zcoil_l(i,2)
         zcoil(i,3)=zcoil_l(i,3)
         zcoil(i,4)=zcoil_l(i,4)
         ccoil(i,1)=ccoil_l(i,1)
         ccoil(i,2)=ccoil_l(i,2)
         ccoil(i,3)=ccoil_l(i,3)
         ccoil(i,4)=ccoil_l(i,4)
      enddo
c
c----- limiter position
c
      do i=1,100
         rlimt(i)=rlimt_l(i)
         zlimt(i)=zlimt_l(i)
      enddo
c
c----- coil file name : 'init' for stdin
c
      knamcoil='init'
c
c----- equ file name : 'equfile'
c
      KNAMEQ='equdata'
      MODEFW=5
      CALL pl_allocate_ns    ! necessary to link plcomm with iort
c
      return
      end subroutine eqinit
c
c=======================================================================
c           input parameters
c=======================================================================
      subroutine eqparm(mode,kin,ierr)
c
c     mode=0 : standard namelinst input
c     mode=1 : namelist file input
c     mode=2 : namelist line input
c
c     ierr=0 : normal end
c     ierr=1 : namelist standard input error
c     ierr=2 : namelist file does not exist
c     ierr=3 : namelist file open error
c     ierr=4 : namelist file read error
c     ierr=5 : namelist file abormal end of file
c     ierr=6 : namelist line input error
c     ierr=7 : unknown mode
c     ierr=10x : input parameter out of range
c
      use aaa_mod
      use vac_mod
      USE libkio
      implicit none
      character kin*(*)
      integer mode,ierr
c
    1 call task_parm(mode,'equ',kin,eqnlin,eqplst,ierr)
      if(ierr.ne.0) return
C
C      call eqchek(ierr)
C      if(mode.eq.0.and.ierr.ne.0) goto 1
C      if(ierr.ne.0) ierr=ierr+100
c
      return
      end subroutine eqparm
c
c=======================================================================
      subroutine eqnlin(nid,ist,ierr)
c=======================================================================
c     read namelist &equ
c=======================================================================
      use aaa_mod
      use par_mod
      use geo_mod
      use equ_mod
      use vac_mod
      use eqv_mod
      use cnt_mod
      implicit none
      integer nid,ist,ierr
      real*8    eqfile ! dummy namelist
      real*8    rwal,zwal
      character aname*(10)
      integer n,nn,il,nid1
c=======================================================================
      namelist /equ/ 
     &     device,devnam,title,nv,nr,nz,
     &     btv,tcur,icp,cp,
     &     rwmx,rwmn,zwmx,iudsym,ivp,
     &     msfx,rvac,zvac,
     &     rmaj,rpla,zpla,elip,trig,yh,yd,elipup,trigup,
     &     qaxi,qsur,xxli,
     &     msetup,esetup,ieqmax,eeqmax,iodmax,eodmax,
     &     iadmax,eadmax,bavmax,bavmin,nsumax,
     &     ivac,cvac,ncoil,cvacst,cvacwg,
     &     isep,dsep,ivtab,ivgrp,rloop,zloop,
     &     ilimt,ivbnd,isvac,ieqout,
     &     ieqrd,jeqrd,ieqwt,jeqwt,
     &     eqfile,rwal,zwal,knamcoil,
     &     cvact
c=======================================================================
      if(nid.lt.0) then
         write(-nid,equ,IOSTAT=ist,ERR=9800)
      ELSE
         read(nid,equ,IOSTAT=ist,ERR=9800,END=9900)
      ENDIF
C
      ierr=0
      return
C
 9800 ierr=8
      return
 9900 ierr=9
      return
      end subroutine eqnlin
c
c=======================================================================
      subroutine eqflin(nid,ierr)
c=======================================================================
c     read coil data 
c=======================================================================
      use aaa_mod
      use par_mod
      use geo_mod
      use equ_mod
      use vac_mod
      use eqv_mod
      use cnt_mod
      implicit none
      integer nid,ist,ierr
      real*8    eqfile ! dummy namelist
      real*8    rwal,zwal
      character aname*(10)
      integer n,nn,il,nid1
c=======================================================================
      if(knamcoil.ne.'done') then
         if(nid.ne.5) 
     &        open(nid,FILE=knamcoil,STATUS='OLD',FORM='FORMATTED')
         do n=1,icvdm
            if(ncoil(n).gt.0)then
               read(nid,'(a10)',ERR=9800,END=9900)aname
               do nn=1,ncoil(n)
                  read(nid,'(3d10.3)',ERR=9800,END=9900)
     &                 rcoil(nn,n),zcoil(nn,n),ccoil(nn,n)
               enddo
            endif
         enddo
c-----------------------------------------------------------------------
c-----read limiter positions
c-----------------------------------------------------------------------
         if(ilimt.gt.0)then
            read(nid,'(a10)',ERR=9800,END=9900)aname
            do il=1,ilimt
               read(nid,'(2d10.3)',ERR=9800,END=9900) 
     &              rlimt(il),zlimt(il)
            enddo
         endif
         if(nid.ne.5) close(nid)
         devnam=device
C     endif
         knamcoil='done'
      endif
C
      ierr=0
      return
C
 9800 ierr=8
      return
 9900 ierr=9
      return
      end subroutine eqflin
c
c=======================================================================
      subroutine eqplst
c=======================================================================
c     show variable name of namelist &equ
c=======================================================================
      write(6,*) '   device,devnam,title,nv,nr,nz,'
      write(6,*) '   btv,tcur,icp,cp,qaxi,qsur,xxli,'
      write(6,*) '   rwmx,rwmn,zwmx,iudsym,ivp,msfx,rvac,zvac,'
      write(6,*) '   rmaj,rpla,zpla,elip,trig,yh,yd,elipup,trigup,'
      write(6,*) '   msetup,esetup,ieqmax,eeqmax,iodmax,eodmax,'
      write(6,*) '   iadmax,eadmax,bavmax,bavmin,nsumax,'
      write(6,*) '   ivac,cvac,ncoil,cvacst,cvacwg,'
      write(6,*) '   isep,dsep,ivtab,ivgrp,rloop,zloop,'
      write(6,*) '   ilimt,ivbnd,isvac,ieqout,'
      write(6,*) '   ieqrd,jeqrd,ieqwt,jeqwt,'
      write(6,*) '   eqfile,rwal,zwal'
      return
      end subroutine eqplst
c
c=======================================================================
      subroutine eqchek(ierr)
c=======================================================================
c     check input data
c=======================================================================
      use aaa_mod
      use par_mod
      use geo_mod
      use equ_mod
      use vac_mod
      use eqv_mod
      use cnt_mod
      implicit none
      integer ierr
! local variables
      integer   i, ic, il, k, kstepk, kudsym, ludsym, n, nc, nfix
     >        , nn, nnc,  istop
      real*8    ex, ex2, ex2up, qx, qxup, r1, r2, ra2
     >        , rrmax, rrmin, rwal, rx, rx2, rx2up, rxup, rzmax
     >        , rzmin, y, z1, zd, zd2, zd2up, zd3, zd3up
     >        , zdup, zh, zh2, zh2up, zh3, zh3up, zhup, zwal, zx
     >        , zx2, zx2up, zxup
c======================================================================
      ierr=0
      if(ieqrd.ne.0) then
         ierr=1
         return
      endif
c======================================================================
      if(nr.gt.irdm)nr=irdm
      if(nz.gt.izdm)nz=izdm
      if(nv.gt.ivdm)nv=ivdm
c-----
      if(nr*nz*nv.eq.0)then
      write(ft06,*)'========== data error : nr*nz*nv=0 ====='
      ierr=2
      return
      endif
c-----
      iudsym=iudsym*isymdm
c-----
      do ic=1,icvdm
      if(cvacst(ic).ne.0.)cvac(ic)=cvacst(ic)
      if(cvac(ic).ne.0.)cvacst(ic)=cvac(ic)
      enddo
c-----------------------------------------------------------------------
c-----read vaccume field coils
c-----------------------------------------------------------------------
      if(device.ne.devnam)then
      kudsym=0
      ludsym=0
      do n=1,icvdm
      if(ncoil(n).gt.0)then
      do nn=1,ncoil(n)
      if(zcoil(nn,n).lt.0.)kudsym=1
      if(zcoil(nn,n).gt.0.)ludsym=1
      enddo
      if(kudsym.eq.1.and.ludsym.eq.1)iudsym=0
      endif
      enddo
      endif
c-----------------------------------------------------------------------
c-----set major radius etc. from marker points ( msfx < 0 )
c-----------------------------------------------------------------------
      if(msfx.lt.0) then
      rrmax=-1.d+10
      rzmax=-1.d+10
      rrmin= 1.d+10
      rzmin= 1.d+10
      do n=0,iabs(msfx)
      if(rvac(n).gt.0.and.rvac(n).gt.rrmax)then
      rrmax=rvac(n)
      rzmax=zvac(n)
      endif
      if(rvac(n).gt.0.and.rvac(n).lt.rrmin)then
      rrmin=rvac(n)
      rzmin=zvac(n)
      endif
      enddo
      rmaj = 0.5*(rrmax+rrmin)
      rpla = 0.5*dabs(rrmax-rrmin)
      zpla = 0.5*(rzmax+rzmin)
      endif
c-----------------------------------------------------------------------
c     parameters of Soloveiev equilibrium
c-----------------------------------------------------------------------
      ra2=rmaj**2+rpla**2
      rx=rmaj-trig*rpla
      zx=elip*rpla
      rx2=rx*rx
      zx2=zx*zx
      qx=1.-2.*ra2*(ra2-rx2)/((rx2-ra2)**2+(2.*rmaj*rpla)**2)
      ex2=4.*ra2*zx2/((rx2-ra2)**2+(2.*rmaj*rpla)**2)
      ex=dsqrt(ex2)
c-----
      asol=ra2
      ssol=ex*(rmaj*rpla)**2
      esol=ex
      qsol=qx
c-----
      if(elipup.le.0.)then
      elipup=elip
      trigup=trig
      endif
c-----
      rxup=rmaj-trigup*rpla
      zxup=elipup*rpla
      rx2up=rxup*rxup
      zx2up=zxup*zxup
      qxup=
     >   1.d0-2.d0*ra2*(ra2-rx2up)/((rx2up-ra2)**2+(2.d0*rmaj*rpla)**2)
      ex2up=4.d0*ra2*zx2up/((rx2up-ra2)**2+(2.d0*rmaj*rpla)**2)
c-----------------------------------------------------------------------
c-----17 marker points (msfx > 0)
c-----------------------------------------------------------------------
      if(msfx.gt.0)then
      if(iudsym.eq.0)then
c-----
        rvac(0)=rmaj-rpla
        rvac(1)=rmaj+rpla
        zvac(0)=zpla
        zvac(1)=zpla
c----
        zd=yd*zx
        zd2=zd*zd
        zd3=zd2/ex2
        y=2.d0*dsqrt(((1.d0-qx)*zd3)**2-ra2*zd3+(rmaj*rpla)**2)
        rvac(2)= dsqrt(ra2-2.d0*(1.d0-qx)*zd3+y)
        rvac(3)= dsqrt(ra2-2.d0*(1.d0-qx)*zd3-y)
        zvac(2)=-zd+zpla
        zvac(3)=-zd+zpla
c--
        zdup=yd*zxup
        zd2up=zdup*zdup
        zd3up=zd2up/ex2up
        y=2.d0*dsqrt(((1.d0-qxup)*zd3up)**2-ra2*zd3up+(rmaj*rpla)**2)
        rvac(6)= dsqrt(ra2-2.d0*(1.d0-qxup)*zd3up+y)
        rvac(7)= dsqrt(ra2-2.d0*(1.d0-qxup)*zd3up-y)
        zvac(6)= zdup+zpla
        zvac(7)= zdup+zpla
c----
        zh=yh*zx
        zh2=zh*zh
        zh3=zh2/ex2
        y=2.d0*dsqrt(((1.d0-qx)*zh3)**2-ra2*zh3+(rmaj*rpla)**2)
        rvac(4)= dsqrt(ra2-2.d0*(1.d0-qx)*zh3+y)
        rvac(5)= dsqrt(ra2-2.d0*(1.d0-qx)*zh3-y)
        zvac(4)=-zh+zpla
        zvac(5)=-zh+zpla
c--
        zhup=yh*zxup
        zh2up=zhup*zhup
        zh3up=zh2up/ex2up
        y=2.d0*dsqrt(((1.d0-qxup)*zh3up)**2-ra2*zh3up+(rmaj*rpla)**2)
        rvac(8)= dsqrt(ra2-2.d0*(1.d0-qxup)*zh3up+y)
        rvac(9)= dsqrt(ra2-2.d0*(1.d0-qxup)*zh3up-y)
        zvac(8)= zhup+zpla
        zvac(9)= zhup+zpla
      if(iabs(msfx).gt.9)then
        zh=zx/4.d0
        zh2=zh*zh
        zh3=zh2/ex2
        y=2.d0*dsqrt(((1.d0-qx)*zh3)**2-ra2*zh3+(rmaj*rpla)**2)
        rvac(10)= dsqrt(ra2-2.d0*(1.d0-qx)*zh3+y)
        rvac(11)= dsqrt(ra2-2.d0*(1.d0-qx)*zh3-y)
        zvac(10)=-zh+zpla
        zvac(11)=-zh+zpla
c--
        zhup=zxup/4.d0
        zh2up=zhup*zhup
        zh3up=zh2up/ex2up
        y=2.d0*dsqrt(((1.d0-qxup)*zh3up)**2-ra2*zh3up+(rmaj*rpla)**2)
        rvac(12)= dsqrt(ra2-2.d0*(1.d0-qxup)*zh3up+y)
        rvac(13)= dsqrt(ra2-2.d0*(1.d0-qxup)*zh3up-y)
        zvac(12)= zhup+zpla
        zvac(13)= zhup+zpla
c-----
        zh=zx/4.d0*3.d0
        zh2=zh*zh
        zh3=zh2/ex2
        y=2.d0*dsqrt(((1.d0-qx)*zh3)**2-ra2*zh3+(rmaj*rpla)**2)
        rvac(14)= dsqrt(ra2-2.d0*(1.d0-qx)*zh3+y)
        rvac(15)= dsqrt(ra2-2.d0*(1.d0-qx)*zh3-y)
        zvac(14)=-zh+zpla
        zvac(15)=-zh+zpla
c--
        zhup=zxup/4.d0*3.d0
        zh2up=zhup*zhup
        zh3up=zh2up/ex2up
        y=2.d0*dsqrt(((1.d0-qxup)*zh3up)**2-ra2*zh3up+(rmaj*rpla)**2)
        rvac(16)= dsqrt(ra2-2.d0*(1.d0-qxup)*zh3up+y)
        rvac(17)= dsqrt(ra2-2.d0*(1.d0-qxup)*zh3up-y)
        zvac(16)= zhup+zpla
        zvac(17)= zhup+zpla
        if(msfx.gt.17)msfx=17
      endif
      else
c-----
        rvac(0)=rmaj-rpla
        rvac(1)=rmaj+rpla
        zvac(0)=zpla
        zvac(1)=zpla
c----
        zd=0.99d0*zx
        zd2=zd*zd
        zd3=zd2/ex2
        y=2.d0*dsqrt(((1.d0-qx)*zd3)**2-ra2*zd3+(rmaj*rpla)**2)
        rvac(2)= dsqrt(ra2-2.d0*(1.d0-qx)*zd3+y)
        rvac(3)= dsqrt(ra2-2.d0*(1.d0-qx)*zd3-y)
        zvac(2)=-zd+zpla
        zvac(3)=-zd+zpla
c----
        zh=zx/2.d0
        zh2=zh*zh
        zh3=zh2/ex2
        y=2.d0*dsqrt(((1.d0-qx)*zh3)**2-ra2*zh3+(rmaj*rpla)**2)
        rvac(4)= dsqrt(ra2-2.d0*(1.d0-qx)*zh3+y)
        rvac(5)= dsqrt(ra2-2.d0*(1.d0-qx)*zh3-y)
        zvac(4)=-zh+zpla
        zvac(5)=-zh+zpla
c----
        zh=3.d0*zx/4.d0
        zh2=zh*zh
        zh3=zh2/ex2
        y=2.d0*dsqrt(((1.d0-qx)*zh3)**2-ra2*zh3+(rmaj*rpla)**2)
        rvac(6)= dsqrt(ra2-2.d0*(1.d0-qx)*zh3+y)
        rvac(7)= dsqrt(ra2-2.d0*(1.d0-qx)*zh3-y)
        zvac(6)=-zh+zpla
        zvac(7)=-zh+zpla
c----
        zh=zx/4.d0
        zh2=zh*zh
        zh3=zh2/ex2
        y=2.d0*dsqrt(((1.d0-qx)*zh3)**2-ra2*zh3+(rmaj*rpla)**2)
        rvac(8)= dsqrt(ra2-2.d0*(1.d0-qx)*zh3+y)
        rvac(9)= dsqrt(ra2-2.d0*(1.d0-qx)*zh3-y)
        zvac(8)=-zh+zpla
        zvac(9)=-zh+zpla
        if(msfx.gt.9)msfx=9
c----
      endif
      endif
c-----------------------------------------------------------------------
c     check the size of calculation box
c-----------------------------------------------------------------------
      r1=rmaj-1.2d0*rpla
      r2=rmaj+1.2d0*rpla
      z1=1.2d0*elip*rpla
      if(rwmn.gt.r1)then
       write(ft06,'(2x,a35,f8.3,a1,f8.3)')
     > 'rwmn is too close to plasma : rwmn=',rwmn,'>',r1
        rwmn=r1
        endif
      if(rwmx.lt.r2)then
       write(ft06,'(2x,a35,f8.3,a1,f8.3)')
     > 'rwmx is too close to plasma : rwmx=',rwmx,'>',r2
        rwmx=r2
        endif
      if(zwmx.lt.z1)then
        write(ft06,'(2x,a35,f8.3,a1,f8.3)')
     > 'zwmx is too close to plasma : zwmx=',zwmx,'>',z1
        zwmx=z1
        endif
c=======================================================================
      nivac = -1
      do i = 0, icvdm
      if( ivac(i).ge.0 ) nivac = nivac + 1
      enddo
c------
      if(nivac.gt.iabs(msfx)) then
        write(ft06,*)'### nivac=',nivac,'  msfx=',msfx
        write(ft06,*)'     **** INCONSISTENT INPUT DATA IS FOUND. ****'
        ierr=3
        return
      endif
c-----
      return
      end subroutine eqchek
c
c=======================================================================
      subroutine eqview
c=======================================================================
c     show input variables
c=======================================================================
      use aaa_mod
      use par_mod
      use geo_mod
      use equ_mod
      use vac_mod
      use eqv_mod
      use cnt_mod
      implicit none
! local variables
      integer   i, n, nc, nnc, mm
c=======================================================================
c     print out
c=======================================================================
      write(ft06,'(//)')
      write(ft06,*)'   namelist &equ'
      write(ft06,'(5x,a12,a10,5x,a4,i1)')
     >       'device name=', device,'sym=',iudsym
      write(ft06,*)'  === equilibrium parameters'
      write(ft06,'(5x,a5,i5,5x,a5,i5,5x,a5,i5)')
     >       'nr  =',nr,'nz  =',nz,'nv  =',nv
      write(ft06,'(3x,5(2x,a5,f8.3))')
     >       'rmaj=',rmaj,'rpla=',rpla,'zpla=',zpla
     >      ,'elip=',elip,'trig=',trig
      write(ft06,'(3x,5(2x,a5,f8.3))')
     >       'rwmn=',rwmn,'rwmx=',rwmx,'zwmx=',zwmx
     >      ,'yh  =',yh  ,'yd  =',yd
      write(ft06,'(3x,5(2x,a5,f8.3))')
     >       'btv =',btv ,'btol=',btol,'tcur=',tcur
     >      ,'qaxi=',qaxi,'qsur=',qsur
      write(ft06,*)'  === initial profile controll parameters'
      write(ft06,'(5x,a4,10(i3,4x))')'icp=',(icp(i),i=1,10)
      write(ft06,'(5x,a4,10f7.3)')'cp =',(cp(i),i=1,10)
      write(ft06,*)'  === equilibrium solver'
      write(ft06,'(5x,a7,i5,5x,a7,1pd10.3,2x,a7,i5)')
     >     'msetup=', msetup,'esetup=',esetup,'nsumax=',nsumax
      write(ft06,'(5x,a7,i5,5x,a7,1pd10.3,2x,a7,i5,7x,a7,1pd10.3)')
     >     'ieqmax=',ieqmax,'eeqmax=',eeqmax 
     >    ,'iodmax=',iodmax,'eodmax=',eodmax 
      write(ft06,'(5x,a7,i5,5x,a7,1pd10.3,2x,a7,1pd10.3,2x,a7,1pd10.3)')
     >     'iadmax=',iadmax,'eadmax=',eadmax 
     >    ,'bavmax=',bavmax,'bavmin=',bavmin 
      write(ft06,*)'  === separatrix'
      write(ft06,'(5x,a5,i2,7x,a5,1pd10.3,2x,a5,i2,2x,a6,i2)')
     >   'isep=',isep,'dsep=',dsep
c-----------------------------------------------------------------------
      write(ft06,*)'  === poloidal coil current control'
      do i = 0, icvdm
       write(ft06,'(3x,a7,i2,a8,i2,a7,f7.3,a9,f7.3,a9,f6.3)')
     > '  ivac=',ivac(i),'  ncoil=',ncoil(i),
     > '  cvac=',cvac(i),'  cvacst=',cvacst(i),'  cvacwg=',cvacwg(i)
      enddo
c------
      write(ft06,*)'  === marker points on the plasma surface '
      write(ft06,'(4x,a5,i6,2x,a6,i6)')'msfz=',msfx,'nivac=',nivac 
      write(ft06,'(3x,i2,a7,f6.3,a7,f6.3)')
     >     (i,'  rvac=',rvac(i),'  zvac=',zvac(i),i=0,iabs(msfx))
c------
      nnc=0
      do n=1,icvdm
      if(ncoil(n).gt.0)nnc=1
      enddo
c----
      if(ieqout(5).eq.0)nnc=0
      if(devnam.ne.device.and.nnc.gt.0)then
      write(ft06,*)'  === poloidal coil system'
      do n=1,icvdm
      nc=ncoil(n)
      if(nc.gt.0)then
      write(ft06,'(a8,i2)')'==coil #',n
      do mm=1,nc
      write(ft06,'(i10,3(2x,f10.5))')
     &        mm,rcoil(mm,n),zcoil(mm,n),ccoil(mm,n)
      enddo
      endif
      enddo
      endif
c-----
      return
      end subroutine eqview
c
      end module eqinit_mod
