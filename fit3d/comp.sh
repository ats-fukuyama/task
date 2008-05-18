#!/bin/csh


cd ./hfreya2
rm *.o *.mod
sxf90 -Chopt -Wf"-A dbl4 -O extendreorder" -c hfrmod.f90  >& ../hfreya2.comp.log
sxf90 -Chopt -Wf"-A dbl4 -O extendreorder" -c hfreya2z_020502Tnlim_gauss2_rdv01.f90 drive_rdv01.f90 iodisk3_hf.f90 hfrfnc.f90 hfrmod.o hfrmod.mod  -lasl >>& ../hfreya2.comp.log
cd ../

cd ./mcnbi
rm *.o *.mod
sxf90 -Chopt -c -Wf"-init heap=zero" mcnmod.f90  >& ../mcnbi.comp.log
sxf90 -Chopt -c -Wf"-init heap=zero" nbidb_040527pb.f90  rkhn2.f90  iodisk41.f90 dinspl.f90 ranu2.f90 mcnmod.o mcnmod.mod -lasl  >>& ../mcnbi.comp.log
cd ../

cd ./fit
rm *.o
sxf90  -c fit_read20_050323.f90 depsum_070416.f   >& ../fit.comp.log
cd ../

sxf90 -o ../exec/mcnbi_all.lm mcnbi_all.f90  ./hfreya2/*.o  ./hfreya2/*.mod ./mcnbi/*.o ./mcnbi/*.mod  ./fit/*.o -lasl  >& ./main.comp.log


