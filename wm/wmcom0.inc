C
C      PARAMETER (NRM=201)
C      PARAMETER (NDPM=0)
      PARAMETER (NRM=201)
      PARAMETER (NDPM=0)
      PARAMETER (MDPM=5)
      PARAMETER (NAM=4)
      PARAMETER (NSUM=2049)
C      PARAMETER (NGZM=1001)
      PARAMETER (NGZM=401)
      PARAMETER (NTHGM=64)
      PARAMETER (MWGM=8)
      PARAMETER (NPHM=128)
C
      PARAMETER (NDM=2**NDPM)
      PARAMETER (MDM=2**MDPM)
      PARAMETER (NMDM=NDM*MDM)
      PARAMETER (NDMX=3*NDM/2)
      PARAMETER (MDMX=3*MDM/2)
      PARAMETER (NHHM=NDM)
      PARAMETER (NTHM=MDM)
      PARAMETER (LFFTM=2*MDM*(MDPM+1))
      PARAMETER (MSIZM=3*NRM*MDM*NDM+MWGM*NAM,MBNDM=12*MDM*NDM)
      PARAMETER (NCPUMIN=1,NCPUMAX=1)
      PARAMETER (NRMD=NRM/NCPUMIN+1)
      PARAMETER (MSIZD=3*NRMD*MDM*NDM+1)
      PARAMETER (MSIZP=MSIZD+MBNDM)
      PARAMETER (NRGM=NRM+1)
