************************************************
*   DOCUMENTATION   Written by HONDA Mitsuru   *
************************************************

*** GUIREAD ***

    This program is a GUI-based reader for ufiles and data provided by
    mdsplus.
    It requires GSGL (GSAF/OpenGL) libraries (v0.4.10 or more), and
    mdsplus libraries if one would like to obtain data from mdsplus
    server. Note that mdsplus libraries are not indispensable for
    compliling this program.
    In the first place, this program reads "guiparm" file existing
    the same directory. "Guiparm" file contains information about 
	       DIR (CHAR) : directory where ufiles are stored
			    The standard directory is "./data".
	       DEV (CHAR) : device(tokamak) name (ex. 'jt60u')
	       SHOT (INT) : discharge number (ex. 29728)
	       DIM  (INT) : dimension (1:time traces, 2:radial profiles)
	       ID  (CHAR) : parameter (ex. 'te')
	       STIME (REAL*4) : slice time (ex. 0.0)
			        this option is valid only if radial
                                profiles have time-evoluted data.
				If one'd like to investigate the radial
				profile at an arbitrary time, one should
				input the finite value, otherwise one
				should input zero.
    and one do not need to define all the parameters.
    This file need not be necessary. Appropriate values are
    automatically set as a default.

    If one read 2-D profile ufiles with time evolution, 3-D view graph
    will appear. If one'd like to know how the radial profile is 
    at an arbitrary time within the available range of time, 
    one should input the value in the box "SLICE=" and then push "DRAW",
    so that the radial profile at user-established time will appear. 
    To see original 3-D type graph, please set "SLICE=0.0" and "DRAW".
    One can also rotate, scale and translate the figure and switch the
    lighting with right-clicking or selecting from the menu bar.

*** UFREAD, UFREAD3 ***

    These programs are CUI-based readers for ufiles.
    UFREAD requires GSAF libraries and UFREAD3 requires GSGL
    (GSAF/OpenGL) libraries.
    Firstly this program reads "ufparm" file existing the same 
    directory. "Ufparm" file has information about only
	       DIR (CHAR) : directory where ufiles are stored.
			    The standard directory is "./data".
    This file need not be necessary. Appropriate value is
    automatically set as a default.
    Ctrl+D terminates it.

*** USPLIT ***

    Quote:
    This routine splits an iter database file into 1D or 2D u-files.
    (The dimension is deduced from the file name)
    last revision: 12/94 s.e.attenberger, w.a.houlberg, ornl
    DISCLAIMER: The authors do not assume any liability for
    deficiencies in this software.  Please report any problems to
    Stan Attenberger at attenberger@fed.ornl.gov
       ***   This software is licensed for export under  ***
       ***   general license agreements GTDA and GUS.    ***

    I mean one need not use this program alone.
    I provided the script program named "UFMAKE" in the following,
    so one might as well use UFMAKE instead of usplit.

*** UFMAKE ***

    This program is written by Perl script for splitting ufiles and
    storing them in an appropriate directory.
    It calls USPLIT and some shell program inside it, so it does not
    work well unless USPLIT.
    It assumes that data files are stored in the subdirectory of 
    '~/profiledb/profile_data/'. For instance, it assumes that data
    files of "JET, 53532" (jet_53532_1d.dat or jet_53532_2d.dat)
    are stored in '~/profiledb/profile_data/jet/53532/'.
    I therefore strongly recommend that all data files be stored in this
    type of form, however if one have one's own style of store, please
    edit "$dir_ufile" in UFMAKE.

    Just answering two questions, one can produce separate ufiles
    from data files into "in" subdirectory, which is automatically made
    if there does not exist it.
    If one is in '~/profiledb/profile_data/jet/53532', for instance,
    one need not answer the questions from UFMAKE.
    Just push the enter key twice.
    This script automatically complements appropriate answers
    instead of you.

*** UFCOMB ***

    This program is written by Perl script for combining separate
    ufiles into one file (original "dat" file).
    It executes the opposite process of USPLIT.

