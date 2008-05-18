#!/bin/csh
 setenv F_PROGINF DETAIL
 setenv F_UFMTIEEE 10
#
 setenv F_FF10 ./indata/boz10.r360q100b004a8020.dat
 setenv F_FF07 ./outdata_pb4/s64773@t03035pb4.out7
 setenv F_FF08 ./outdata_pb4/s64773@t03035pb4.ps
 setenv F_FF09 ./outdata_pb4/s64773@t03035pb4.out9
 setenv F_FF15 ./outdata_pb4/s64773@t03035pb4.out15
 setenv F_FF21 ./outdata_pb4/s64773@t03035pb4.out21
#
 setenv F_FF11 /dev/null
 setenv F_FF12 /dev/null
 setenv F_FF14 /dev/null
 setenv F_FF16 ./indata/aurora.fu08
 setenv F_FF20 ./outdata_pb4/s64773@t03035pb4.out20
 setenv F_FF30 ./outdata_pb4/s64773@t03035pb4.out30
 setenv F_FF40 ./outdata_pb4/s64773@t03035pb4.out40
#
 setenv F_FF60 ./outdata_pb4/fit_s64773@t03035pb4.out10
 setenv F_FF61 ./outdata_pb4/fit_s64773@t03035pb4.out11
 setenv F_FF70 ./outdata_pb4/fit_s64773@t03035pb4.out20
 setenv F_FF71 ./outdata_pb4/fit_s64773@t03035pb4.out21
 setenv F_FF80 ./outdata_pb4/fit_s64773@t03035pb4.out30
 setenv F_FF81 ./outdata_pb4/fit_s64773@t03035pb4.out31
 setenv F_FF90 ./outdata_pb4/fit_s64773@t03035pb4.out40
#

 ln -s ./indata/s64773@t03035pb4.in       hfreya.in
 ln -s ./indata/mcnbi_s64773@t03035pb4.in mcnbi.in
 ln -s ./indata/fit_s64773@t03035pb4.in   fit.in

 time ./mcnbi_all.lm pb4

 mv hfreya.out ./outdata_pb4/s64773@t03035pb4.out6
 mv mcnbi.out  ./outdata_pb4/mcnbi_s64773@t03035pb4.out
 mv fit.out    ./outdata_pb4/fit_s64773@t03035pb4.out

 rm  hfreya.in mcnbi.in fit.in

