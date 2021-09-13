       Program read

       IMPLICIT NONE
       INTEGER:: NP, NTH, i, j, k 
       INTEGER,parameter:: NPMAX=150, NTHMAX=50
       INTEGER,dimension(150*50):: TH, P, TH2, P2
       REAL(rkind),dimension(150,50)::FNS1, FNS2
       REAL(rkind)::f


       open(8,file='FNS2_init.dat')
       DO NP=1,NPMAX
          DO NTH=1,NTHMAX   
             k=NTH+(NP-1)*NTHMAX
             read(8,*) TH(k), P(k), FNS1(NTH,NP)
          END DO
!          read(8,*) i,j,f
!          read(8,*) i,j,f
       END DO
       close(8)

       open(8,file='FNS2_t1.dat')
       DO NP=1,NPMAX
          DO NTH=1,NTHMAX   
             k=NTH+(NP-1)*NTHMAX
             read(8,*) TH2(k), P2(k), FNS2(NTH,NP)
          END DO
!          read(8,*) i,j,f
!          read(8,*) i,j,f
       END DO
       close(8)

       open(8,file='dif_FNS.dat')
       DO NP=1,NPMAX
          DO NTH=1,NTHMAX   
             write(8,*) TH(k), P(k), FNS2(NTH,NP) - FNS1(NTH,NP)
          END DO
          WRITE(8,*) " "
          WRITE(8,*) " "
       END DO
       close(8)


       END Program read
