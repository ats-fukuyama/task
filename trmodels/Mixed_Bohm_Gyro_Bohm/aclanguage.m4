dnl Copyright (C) 1996 John W. Eaton <jwe@bevo.che.wisc.edu>
dnl Copyright (C) 1998 Eleftherios Gkioulekas <lf@amath.washington.edu>
dnl Copyright (C) 1998 Gary Mathers <Gary.Mathers@cern.ch>
dnl Copyright (C) 2002 Alex Pankin <pankin@fusion.physics.lehigh.edu>
dnl 

dnl This program is free software; you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation; either version 2 of the License, or
dnl (at your option) any later version.
dnl 

dnl This program is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.
dnl 
dnl You should have received a copy of the GNU General Public License
dnl along with this program; if not, write to the Free Software 
dnl Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

dnl The following macros are contained within this file:
dnl	GM_PROG_F90
dnl	GM_PROG_F90_VENDOR
dnl	GM_PATH_F90_LIBS
dnl	GM_LANG_FORTRAN90
dnl	AC_LANG_SAVE
dnl	AC_LANG_RESTORE

# The following set of macros will allow you to mix Fortran90 and C in
# a portable manner. This work is based on the autoconf macros written
# by Eleftherios Gkioulekas who in turn modified work by John W. Eaton
# for GNU Octave, which is also distributed under the terms of the GNU
# public license.


# -------------------------------------------------------------------------
# This is the macro that you want to call if you want to use Fortran 90.
# This macro sets F90 equal to a valid Fortran 90 compiler and F90FLAGS
# to a set of flags to pass to that compiler.
# -------------------------------------------------------------------------

AC_DEFUN(GM_PROG_F90,[
  dnl Initialize the following use variables to false
  dnl These variables indicate which compiler we want to use
  dnl
  dnl Now assign F90 with the appropriate Fortran compiler 
  dnl 
  dnl This is to be fixed when I figure out what to do about all
  dnl fortran90 compilers being called f90 ???
  dnl As new names become known e.g. vf90 insert into variable below

  if test -z "$F90"; then
    gm_f90_native_compiler_list="pgf90 lf95 xlf90 f90 f95"
    for f90_comp in $gm_f90_native_compiler_list; do
      AC_PATH_PROG(F90, [$f90_comp],)
      if test -n "$F90"; then
        break
      fi
    done
    dnl Change this line to an error after development is completed
  fi
  dnl replace the warn with the error msg after development.

  if test -z "$F90"; then
    AC_MSG_ERROR(unable to locate f90 compiler)
  fi

  dnl 
  dnl It ensures that an unknown compiler behaves like the NAG. This may 
  dnl need to be changed as I have very little experience with F90 compilers
  dnl and this may not be suitable.
  dnl
  GM_PROG_F90_VENDOR

  dnl By default, compile Fortran with optimization
  dnl F90FLAGS="-O"

  dnl Export F90 and FLAGS to Automake
  FC="$F90"
  FFLAGS="$F90FLAGS"
  AC_SUBST(FC)
  AC_SUBST(F90)
  AC_SUBST(FFLAGS)
  AC_SUBST(F90FLAGS)
  AC_SUBST(F90DFLAGS)
  AC_SUBST(MFLAG)
])

# ---------------------------------------------------------------------------
# THIS macro tests whether the compiler assigned to F90 is the F90
# compiler. If this is the NAG compiler, then set f90_is_nag equal to "true".
# Otherwise, it is set to be an empty string.
# ---------------------------------------------------------------------------

AC_DEFUN(GM_PROG_F90_VENDOR,[
  AC_MSG_CHECKING([the compiler vendor])
  f90_is_nag=no
  f90_is_dec=no
  f90_is_fujitsu=no
  f90_is_lahey=no
  f90_is_pgp=no
  f90_is_compaq=no
  f90_is_sun=no
  f90_is_sgi=no
  f90_is_ibm=no
  f90_is_unknown=no
  AC_TRY_COMMAND(${F90} -version)>&5
  AC_TRY_COMMAND(${F90} --version)>&5
  AC_TRY_COMMAND(${F90} -V) >&5
  dnl AC_TRY_COMMAND(${F90} -version)
  dnl (eval "${F90} -version"; "${F90} --version"; "${F90} -V")>&AC_FD_CC
  if egrep 'NAGWare f90' config.log>/dev/null 2>&1; then
    f90_is_nag=yes
    if (test -z "$F90FLAGS"); then
      F90FLAGS="-O -w=all -mismatch_all"
    fi
    if (test -z "$F90DFLAGS"); then
      F90DFLAGS="-g -g90 -w=all -mismatch_all"
    fi
    MFLAG="-M"
    AC_MSG_RESULT([NAGWare Fortran 90])
  elif egrep 'NAGWare Fortran 95' config.log >/dev/null 2>&1; then
    f90_is_nag=yes
    if (test -z "$F90FLAGS"); then
      F90FLAGS="-O -w=all -mismatch_all"
    fi
    if (test -z "$F90DFLAGS"); then
      F90DFLAGS="-g -g90 -w=all -mismatch_all"
    fi
    MFLAG="-M"
    AC_MSG_RESULT([NAGWare Fortran 95])
  elif egrep 'Lahey/Fujitsu' config.log>/dev/null 2>&1; then
    f90_is_lahey=yes
    if test -z "$F90FLAGS"; then
      F90FLAGS="-O"
    fi
    if test -z "$F90DFLAGS"; then
      F90DFLAGS="-g --chk a,e,s,u"
    fi
    MFLAG="-M"
    AC_MSG_RESULT([Lahey/Fujitsu])
  elif egrep 'MIPSpro' config.log >/dev/null 2>&1; then
    f90_is_sgi=yes
    if (test -z "$F90FLAGS"); then
      F90FLAGS="-O"
    fi
    if (test -z "$F90DFLAGS"); then
      F90DFLAGS="-g -O0 -check_bounds"
    fi
    MFLAG="-I"
    AC_MSG_RESULT([SGI MIPSpro Compiler])
  elif egrep 'DIGITAL Fortran 90' config.log>/dev/null 2>&1; then
    f90_is_dec=yes
    if (test -z "$F90FLAGS"); then
      F90FLAGS="-O -align dcommons"
    fi
    if (test -z "$F90DFLAGS"); then
      F90DFLAGS="-g -O0 -C -align dcommons"
    fi
    MFLAG="-M"
    AC_MSG_RESULT([DEC Compiler])
  elif egrep 'AIX XL Fortran' config.log>/dev/null 2>&1; then
    f90_is_ibm=yes
    if (test -z "$F90FLAGS"); then
      F90FLAGS="-O"
    fi
    if (test -z "$F90DFLAGS"); then
      F90DFLAGS="-g"
    fi
    MFLAG="-M"
    AC_MSG_RESULT([AIX XL Fortran Compiler])
  elif egrep 'WorkShop Compilers' config.log>/dev/null 2>&1; then
    f90_is_sun=yes
    if (test -z "$F90FLAGS"); then
      F90FLAGS="-fast -ftrap=common"
    fi
    if (test -z "$F90DFLAGS"); then
      F90DFLAGS="-dalign -fns -libmil -xlibmopt -g -ftrap=common"
    fi
    MFLAG="-M"
    AC_MSG_RESULT([SUN Fortran-90 Compiler])
  elif egrep 'Fujitsu' config.log>/dev/null 2>&1; then
    f90_is_fujitsu=yes
    if (test -z "$F90FLAGS"); then
      F90FLAGS="-O -f1444,2004,2006,2008 -X9  -Kfast -Kfap -Am -Nallextput"
    fi
    if (test -z "$F90DFLAGS"); then
      F90DFLAGS="-g -f1444,2004,2006,2008 -X9 -Haesu -Kfap -Am -Nallextput"
    fi
    MFLAG="-M"
    AC_MSG_RESULT([Fujitsu])
  elif egrep 'Compaq Fortran' config.log>/dev/null 2>&1; then
    f90_is_compaq=yes
    if (test -z "$F90FLAGS"); then
      F90FLAGS="-O -align dcommons -align sequence -assume no2underscore"
    fi
    if (test -z "$F90DFLAGS"); then
      F90DFLAGS="-g -O0 -C -align dcommons -align sequence -assume no2underscore"
    fi
    MFLAG="-I"
    AC_MSG_RESULT([Compaq Fortran])
  else
    f90_is_unknown=yes
    if (test -z "$F90FLAGS"); then
      F90FLAGS="-O"
    fi
    if (test -z "$F90DFLAGS"); then
      F90DFLAGS="-g -O0 -check_bounds"
    fi
    MFLAG="-M"
    AC_MSG_RESULT(Fortran 90 Compiler is unknown ... defaulting to SGI Compiler)
  fi
])




# --------------------------------------------------------------------------
# See what libraries are used by the Fortran compiler
# Write a minimal program and compiler it with -v. I don't know what to
# do if your compiler doesn't have -v
# The result is returned in the variable F90LIBS which is made
# available in Makefile.am
# ALSO: requires ac_cv_prog_gcc
# --------------------------------------------------------------------------

AC_DEFUN(GM_PATH_F90_LIBS,[
  AC_MSG_CHECKING(for Fortran libraries)
  dnl
  dnl Write a minimal program and compile it with -v. I don't know
  dnl what to do if your compiler doesn't have -v
  dnl
  changequote(, )dnl
  echo "      END" > conftest.f90
   dnl FIXME: import f90flags incase of Qpath option $F90FLAGS
  if test "$f90_is_nag" = "yes"; then
    foutput=`${F90-f90} ${FFLAGS} -ldarg -v -o conftest conftest.f90 2>&1`
  else
    foutput=`${F90-f90} ${FFLAGS} -v -o conftest conftest.f90 2>&1`
  fi
  dnl
  dnl The easiest thing to do for xlf output is to replace all the commas
  dnl with spaces.  Try to only do that if the output is really from xlf,
  dnl since doing that causes problems on other systems.
  dnl
  xlf_p=`echo $foutput | grep xlfentry`
  if test -n "$xlf_p"; then
    foutput=`echo $foutput | sed 's/,/ /g'`
  fi
  dnl
  ld_run_path=`echo $foutput | \
    sed -n -e 's/^.*LD_RUN_PATH *= *\([^ ]*\).*/\1/p'`
  dnl
  dnl We are only supposed to find this on Solaris systems...
  dnl Uh, the run path should be absolute, shouldn't it?
  dnl
  case "$ld_run_path" in
    /*)
      if test "$ac_cv_prog_gcc" = yes; then
        ld_run_path="-Xlinker -R -Xlinker $ld_run_path"
      else
        ld_run_path="-R $ld_run_path"
      fi
    ;;
    *)
      ld_run_path=
    ;;
  esac
  dnl
  f90libs=
  lflags=
  dnl
  dnl If want_arg is set, we know we want the arg to be added to the list,
  dnl so we don't have to examine it.
  dnl
  want_arg=
  dnl
  for arg in $foutput; do
    old_want_arg=$want_arg
    want_arg=
  dnl
  dnl None of the options that take arguments expect the argument to
  dnl start with a -, so pretend we didn't see anything special.
  dnl
    if test -n "$old_want_arg"; then
      case "$arg" in
        -*)
        old_want_arg=
        ;;
      esac
    fi
    case "$old_want_arg" in
      '')
        case $arg in
        /*.a)
          exists=false
          for f in $lflags; do
            if test x$arg = x$f; then
              exists=true
            fi
          done
          if $exists; then
            arg=
          else
            lflags="$lflags $arg"
          fi
        ;;
        -bI:*)
          exists=false
          for f in $lflags; do
            if test x$arg = x$f; then
              exists=true
            fi
          done
          if $exists; then
            arg=
          else
            if test "$ac_cv_prog_gcc" = yes; then
              lflags="$lflags -Xlinker $arg"
            else
              lflags="$lflags $arg"
            fi
          fi
        ;;
        -lang* | -lcrt0.o | -lc | -lgcc)
          arg=
        ;;
        -[lLR])
          want_arg=$arg
          arg=
        ;;
        -[lLR]*)
          exists=false
          for f in $lflags; do
            if test x$arg = x$f; then
              exists=true
            fi
          done
          if $exists; then
            arg=
          else
            case "$arg" in
              -lkernel32)
                case "$canonical_host_type" in
                  *-*-cygwin32)
                  ;;
                  *)
                    lflags="$lflags $arg"
                  ;;
                esac
              ;;
              -lm)
              ;;
              *)
                lflags="$lflags $arg"
              ;;
            esac
          fi
        ;;
        -u)
          want_arg=$arg
          arg=
        ;;
        -Y)
          want_arg=$arg
          arg=
        ;;
        *)
          arg=
        ;;
        esac
      ;;
      -[lLR])
        arg="$old_want_arg $arg"
      ;;
      -u)
        arg="-u $arg"
      ;;
      -Y)
  dnl
  dnl Should probably try to ensure unique directory options here too.
  dnl This probably only applies to Solaris systems, and then will only
  dnl work with gcc...
  dnl
        arg=`echo $arg | sed -e 's%^P,%%'`
        SAVE_IFS=$IFS
        IFS=:
        list=
        for elt in $arg; do
        list="$list -L$elt"
        done
        IFS=$SAVE_IFS
        arg="$list"
      ;;
    esac
  dnl
    if test -n "$arg"; then
      f90libs="$f90libs $arg"
    fi
  done
  if test -n "$ld_run_path"; then
    f90libs_result="$ld_run_path $f90libs"
  else
    f90libs_result="$f90libs"
  fi
  changequote([, ])dnl
  rm -f conftest.f conftest.o conftest
  dnl
  dnl Phew! Done! Now, output the result
  dnl
  F90LIBS="$f90libs_result"
  FLIBS="$F90LIBS"
  AC_MSG_RESULT([$FLIBS])
  AC_SUBST(FLIBS)
])

dnl Do compilation tests using F90 and use extension .f90 for
dnl test programs.  
dnl 
dnl GM_LANG_FORTRAN90()
AC_DEFUN([GM_LANG_FORTRAN90],
[define([AC_LANG], [FORTRAN90])dnl
ac_ext=f90
ac_compile='$F90 $F90LIBS $F90FLAGS -c conftest.$ac_ext 1>&AC_FD_CC'
ac_link='$F90 $F90FLAGS $LDFLAGS -c conftest.$ac_ext -o conftest $F90LIBS $LIBS 1>&AC_FD_CC'
AC_MSG_RESULT(default language set to Fortran 90)
])


dnl Remember the current language {AC_LANG_C},{AC_LANG_CPLUSPLUS}, {AC_LANG_FORTRAN77}
dnl or GM_LANG_FORTRAN90 on a stack.
dnl Does not change which language is current.  Use this macro and
dnl AC_LANG_RESTORE in macros that need to temporarily switch to
dnl a particular language.
dnl AC_LANG_SAVE()
pushdef([AC_LANG_SAVE],
[pushdef([AC_LANG_STACK], AC_LANG)])


dnl Select the language that is saved on the top of the stack, as set by
dnl @code{AC_LANG_SAVE}, and remove it from the stack.  This macro is
dnl equivalent to either AC_LANG_C, AC_LANG_CPLUSPLUS, AC_LANG_FORTRAN77 or
dnl GM_LANG_FORTRAN90 whichever had been run most recently when
dnl AC_LANG_SAVE was last called.
dnl AC_LANG_RESTORE()
pushdef([AC_LANG_RESTORE],
[ifelse(AC_LANG_STACK, [C], [AC_LANG_C],dnl
AC_LANG_STACK, [CPLUSPLUS], [AC_LANG_CPLUSPLUS],dnl
AC_LANG_STACK, [FORTRAN77], [MDL_LANG_FORTRAN77],dnl
AC_LANG_STACK, [FORTRAN90], [GM_LANG_FORTRAN90]),[]popdef([AC_LANG_STACK])])
