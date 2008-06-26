#!/bin/csh
cd /home/murakami/fort/ktestNBI
sxf90 -C hopt -o mcnbi_LHDstd_N1000_040527pb.lm nbidb_040527pb.f  rkhn2.f  iodisk41.f dinspl.f ranu2.f ranu2n.f -lasl
