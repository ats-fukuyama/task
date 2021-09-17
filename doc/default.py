#!/usr/bin/env python

# Build PETSc, with macports including mpich and gcc

configure_options = [
  '--with-mpi-dir=/opt/local',
  '--with-shared-libraries=0',
  '--with-cxx-dialect=C++11',
  '--download-mpich=0',
  '--download-hypre=0',
  '--download-fblaslapack=1',
  '--download-spooles=1',
  '--download-superlu=1',
  '--download-metis=1',
  '--download-parmetis=1',
  '--download-superlu_dist=0',
  '--download-blacs=1',
  '--download-scalapack=1',
  '--download-mumps=1'
  ]

if __name__ == '__main__':
  import sys,os
  sys.path.insert(0,os.path.abspath('config'))
  import configure
  configure.petsc_configure(configure_options)
