#
# Python preprocessor using Doug's $CODESYSDIR/source/bpython/fpp.py code.
#

# ---- import ----
import os
import sys

#
# transp paths and imports
#

if (os.environ.has_key('CODESYSDIR')):

  CODESYSDIR = os.environ['CODESYSDIR']
  bpython    = os.path.join(CODESYSDIR,"source","bpython")
else:
  bpython    = os.path.join("..","fpreproc")

sys.path.append(bpython)   # transp's python scripts
from fpp import *

# -- parse arguments --
args = sys.argv
if (len(args)<3):
    print """
    usage,

    python fppfile.py  fileIn  fileOut [ -132 -free  options... ]

         fileIn  = input file
         fileOut = output file
         -132    = set line length to 132 characters
         -free   = free format source
         -info   = debug messages
         -keep   = do not delete intermediate file (for debugging)
         -fppcom <fppcom> = change the preprocessor executable
         options = -I and -D preprocessor directives

    """
    sys.exit(1)
    
fileIn  = args[1]
fileOut = args[2]

macros = ""
maxlen = 72
free   = 0
info   = 0
keep   = 0
fppcom = None

xfppcom=0    # flag next argument is fppcom

for x in args[3:]:
  if (xfppcom):
    fppcom = x
    xfppcom=0  
  elif (x == "-132"):
    maxlen = 132
  elif (x == "-free"):
    free = 1
  elif (x == "-keep"):
    keep = 1
  elif (x == "-info"):
    info = 1
  elif (x == "-fppcom"):
    xfppcom = 1
  else:
    macros += " %s"%x


#
# -- preprocess --
# temporary file will be based off python's os.tempnam() function which will use TMPDIR
# along with the filename as prefix
#
fppFile(fileIn, fileOut, macros, free=free, maxlen=maxlen, fileTemp=None,
        fppcom=fppcom, info=info, keep=keep)

