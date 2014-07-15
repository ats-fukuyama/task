!
!	BOZdata_transfer.f90
!	task3d_beta_s109695_nbi12x45_experimentaldensity
!
!	Created by WAKASA Arimitsu on 11/11/02.
!	Copyright 2011 __MyCompanyName__. All rights reserved.
!

subroutine bozdata_translation_yt2mb(BOZDATA_PATH)
    use task3d_tentative_flags

    implicit none
    integer,parameter :: NSD_MAX=300
    integer,parameter :: NMBOZ_MAX=3000
    character*(*) :: BOZDATA_PATH

!    integer ::NSD=60
    integer, parameter:: ndiskc_mb=31
    integer, parameter:: ndiskc_yt=32
    character(21) ::BOZ_mb="BOZdata_mb_format.dat"
    integer m,i
    integer :: nsd
    integer :: nmboz,nsdm1,nfp
    real(8) :: rmajmy, dummy0
    real(8), dimension(NSD_MAX) :: psi,eot,cui,cug,dumnsd
    integer(4), dimension(NMBOZ_MAX) :: MBOZ, NBOZ
    real(8),dimension(NMBOZ_MAX,NSD_MAX) ::  BBOZH, RBOZH, ZBOZH, PBOZH

! =============== Yokoyama txt ==============================

    open(ndiskc_yt,file=BOZDATA_PATH,status='old')

    read(ndiskc_yt,*) nmboz, nsdm1, nfp
    NSD=nsdm1+1
    if( nmboz .gt. NMBOZ_MAX ) then
        print*,'nmboz=',NMBOZ,'>','NMBOZ_MAX=',NMBOZ_MAX
        stop
    endif
    if( nsdm1+1 .gt. NSD_MAX    ) then
        print*,'nsdm1+1',nsdm1+1,'>','NSD_MAX',NSD_MAX
        stop
    endif    
    
    
    READ(NDISKC_yt,*) (psi(i), i = 2, NSD)
    READ(NDISKC_yt,*) (eot(i), i = 2, NSD)
    READ(NDISKC_yt,*) (cui(i), i = 2, NSD)
    READ(NDISKC_yt,*) (cug(i), i = 2, NSD)

    DO M = 1, NMBOZ
        READ(NDISKC_yt,*)  MBOZ(M), NBOZ(M)
        READ(NDISKC_yt,*) (BBOZH(M,I), I = 2, NSD)
    enddo
    DO M = 1, NMBOZ
        READ(NDISKC_yt,*)  MBOZ(M), NBOZ(M)
        READ(NDISKC_yt,*) (RBOZH(M,I), I = 2, NSD)
    enddo
    DO M = 1, NMBOZ
        READ(NDISKC_yt,*)  MBOZ(M), NBOZ(M)
        READ(NDISKC_yt,*) (ZBOZH(M,I), I = 2, NSD)
    end do
    DO M = 1, NMBOZ
        READ(NDISKC_yt,*)  MBOZ(M), NBOZ(M)
        READ(NDISKC_yt,*) (PBOZH(M,I), I = 2, NSD)
    enddo
    
! =============== Translate from RHS to LHS  ================
    if(need_boz_Translate_flg==1)then
    do i=2,nsd
        cug(i)=-1.d0*cug(i)
        eot(i)=-1.d0*eot(i)
    enddo
    DO M = 1, NMBOZ
!        MBOZ(M)=-1*MBOZ(M)
        NBOZ(M)=-1*NBOZ(M)        
        do i=2,NSD
            PBOZH(M,I)=PBOZH(M,I)*-1.d0
        enddo
    enddo
    endif
    
! =============== Murakami bin ==============================    

    open(ndiskc_mb,file=BOZ_mb,status='replace',FORM='UNFORMATTED', CONVERT='BIG_ENDIAN')

    write(ndiskc_mb) nmboz, nsdm1, nfp
    write(ndiskc_mb) (psi(i), i = 2, NSD)
    write(ndiskc_mb) (eot(i), i = 2, NSD)
!
    write(ndiskc_mb) (cui(i), i = 2, NSD)
    write(ndiskc_mb) (cug(i), i = 2, NSD)

    DO M = 1, NMBOZ
        write(ndiskc_mb)  MBOZ(M), NBOZ(M)
        write(ndiskc_mb) (BBOZH(M,I), I = 2, NSD)
    enddo
    DO M = 1, NMBOZ
        write(ndiskc_mb)  MBOZ(M), NBOZ(M)
        write(ndiskc_mb) (RBOZH(M,I), I = 2, NSD)
    enddo
    DO M = 1, NMBOZ
        write(ndiskc_mb)  MBOZ(M), NBOZ(M)
        write(ndiskc_mb) (ZBOZH(M,I), I = 2, NSD)
    enddo
!
    DO M = 1, NMBOZ
        write(ndiskc_mb)  MBOZ(M), NBOZ(M)
        write(ndiskc_mb) (PBOZH(M,I), I = 2, NSD)
    enddo
    
! ↓↓BOZデータ仕様確認用↓↓
!    close(ndiskc_yt)
!    close(ndiskc_mb)
!    
!    print*,'BOZ_CHECK!!'
!    open(ndiskc_mb,file='boz10.r360q100b050a8020.dat',status='old',FORM='UNFORMATTED', CONVERT='BIG_ENDIAN')
!    open(ndiskc_yt,file='boz10.r360q100b050a8020.txt',status='replace')
!
!    read(ndiskc_mb) nmboz, nsdm1, nfp
!    read(ndiskc_mb) (psi(i), i = 2, NSD)
!    read(ndiskc_mb) (eot(i), i = 2, NSD)
!    read(ndiskc_mb) (cui(i), i = 2, NSD)
!    read(ndiskc_mb) (cug(i), i = 2, NSD)
!
!    DO M = 1, NMBOZ
!        read(ndiskc_mb)  MBOZ(M), NBOZ(M)
!        read(ndiskc_mb) (BBOZH(M,I), I = 2, NSD)
!    enddo
!    DO M = 1, NMBOZ
!        read(ndiskc_mb)  MBOZ(M), NBOZ(M)
!        read(ndiskc_mb) (RBOZH(M,I), I = 2, NSD)
!    enddo
!    DO M = 1, NMBOZ
!        read(ndiskc_mb)  MBOZ(M), NBOZ(M)
!        read(ndiskc_mb) (ZBOZH(M,I), I = 2, NSD)
!    enddo
!    DO M = 1, NMBOZ
!        read(ndiskc_mb)  MBOZ(M), NBOZ(M)
!        read(ndiskc_mb) (PBOZH(M,I), I = 2, NSD)
!    enddo
!    
!    write(ndiskc_yt,*) nmboz, nsdm1, nfp
!    write(ndiskc_yt,*) (psi(i), i = 2, NSD)
!    write(ndiskc_yt,*) (eot(i), i = 2, NSD)
!    write(ndiskc_yt,*) (cui(i), i = 2, NSD)
!    write(ndiskc_yt,*) (cug(i), i = 2, NSD)
!
!    DO M = 1, NMBOZ
!        write(ndiskc_yt,*)  MBOZ(M), NBOZ(M)
!        write(ndiskc_yt,*) (BBOZH(M,I), I = 2, NSD)
!    enddo
!    DO M = 1, NMBOZ
!        write(ndiskc_yt,*)  MBOZ(M), NBOZ(M)
!        write(ndiskc_yt,*) (RBOZH(M,I), I = 2, NSD)
!    enddo
!    DO M = 1, NMBOZ
!        write(ndiskc_yt,*)  MBOZ(M), NBOZ(M)
!        write(ndiskc_yt,*) (ZBOZH(M,I), I = 2, NSD)
!    enddo
!    DO M = 1, NMBOZ
!        write(ndiskc_yt,*)  MBOZ(M), NBOZ(M)
!        write(ndiskc_yt,*) (PBOZH(M,I), I = 2, NSD)
!    enddo      
!          
! ↑↑旧BOZデータ仕様確認用↑↑
    
    close(ndiskc_mb)
    close(ndiskc_yt)

end subroutine bozdata_translation_Yt2Mb
