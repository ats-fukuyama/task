#!/bin/bash

if [ -f testmmm ]; then
   for i in case*
   do
      echo
      echo ====================
      echo "  $i"
      echo ====================
      echo

      cd $i
      cp sample_input input
      ../testmmm
      diff -q output sample_output
      cd ..

   done
else
   echo "ERROR: testmmm executable has not been built!"
fi

