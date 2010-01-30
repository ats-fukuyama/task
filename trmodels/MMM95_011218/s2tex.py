#!/loc/bin/python
# s2tex.py
#
# This python program converts from source (f, py, c, java, h files) 
# to tex files and back again
#
# When converting from file.tex to file.src:
#   All running text (TeX or LaTeX) lines are converted to lines
#     that start with comment_specificator|
#   All lines between \begin{verbatim} and \end{verbatim}
#     are converted to source code (ie, they are copied verbatim)
#   The \begin{verbatim} and \end{verbatim} lines are removed
#
# When converting from file.src to file.tex, the process is reversed.
# Besides two addtional files are produced. The first one is source files
# without tex comment (file.pure.src) and the other one is tex file that
# consists only from tex part (file.pure.tex)
#

import sys                      # load the system module
import string                   # load the string module

# Global Variables
languageSpecific = {"py": (("#"), ()), \
		    "f": (("!", "c", "C"), ()), \
		    "for": (("!", "c", "C"), ()), \
                    "f90": (("!", "c", "C"), ()), \
                    "java": (("//", ""), ()), \
                    "c": (("//", ""), ("/*", "*/")), \
                    "C": (("//", ""), ("/*", "*/")), \
                    "h": (("//", ""), ("/*", "*/")), \
                    "cpp": (("//", ""), ("/*", "*/"))} 
texSpecific      = {"tex": (("%", ""), ())}

def Initialize(fileName) :
  # determine the file basename and suffix extension
  global suffix
  global input
  global output, outputTex, outputSrc
  try :
    baseFileNameLength = string.rfind( fileName, '.')
    baseFileName = fileName[ :baseFileNameLength ]
    suffix = fileName[ baseFileNameLength+1: ]
    print "Input file name:  " + fileName
  except (Exception, e) :  
    print "Invalid input: " + filename
    print "Exception: " + e
    return 0
  try :  
    if languageSpecific.has_key(suffix) :
      input  = open( fileName, 'r' )          # open input  file
      outputFileName = baseFileName + ".tex"
      output = open( outputFileName, 'w' )    # open output file
      output.write(texSpecific["tex"][0][0] + fileName + "\n")
      print "Output file name: " + outputFileName
      outputFNTex = baseFileName + ".pure.tex"
      outputTex = open( outputFNTex, 'w' )
      outputTex.write(2*texSpecific["tex"][0][0] + fileName + "\n")
      print "Output file name: " + outputFNTex + " (tex file without source)"
      outputFNSrc = baseFileName + ".pure." + suffix
      outputSrc = open( outputFNSrc, 'w' )
      outputSrc.write(languageSpecific[suffix][0][0] + fileName + "\n")
      print "Output file name: " + outputFNSrc + " (source file without tex comments)"
      return 1
    elif texSpecific.has_key(suffix) :
      input  = open( fileName, 'r' )          # open input  file
      line = input.readline()                 # read first line
      texcomment = texSpecific[suffix][0][0]
      if (line[0:1] == 2 * texcomment) :
        print fileName + " could not be converted to the source file."
        return 0
      elif line[0] == texcomment :
        list = string.split(line[1:])  
        outputFileName = list[0]
        output = open( outputFileName, 'w' )    # open output file
        print "Output file name: " + outputFileName
        baseFileNameLength = string.rfind(outputFileName, '.')
        suffix = outputFileName[ baseFileNameLength+1: ]
        return 2
      else : 
        print "Invalid structure of input file " + fileName
        return 0
    else :
      print "Unknown extension of the input file " + fileName
      return 0
  except (Exception, e) :
    print "Error: " + fileName
    print e
    return 0
  return 1
  
  
#argv = ['c:\\data\\python\\fbes.f']
#for filename in argv :
for filename in sys.argv[1:]:   # go through command line items
  
  iniFile = Initialize(filename)
  
  # If the file is a source file, convert to a *.tex file

  if ( iniFile == 1 ):              # convert to .tex file
    # linetype = 0 for source code line
    # linetype > 0 for tex line
    linetype = 1

    # start .tex file with the name of the source file
    line = ""
    while 1:
      oldline = line                # save last line
      oldlinetype = linetype
      line = input.readline()       # read one line at a time
      if not line: break            # break out at end of file
      line = string.rstrip(string.expandtabs(line, 8)) + '\012'

      linetype = 0
      for texcomment in languageSpecific[suffix][0] :
        texcomment1 = texcomment + "|"   # for lines with no space after |
        texcomment2 = texcomment1 + " "  # for lines with a space after |
        if ( string.find(line, texcomment1) >= 0 ) :
          linetype = 1
          if (oldlinetype == 0):
            output.write("\\end{verbatim}\n")
          if ( string.find(line, texcomment2) >= 0 ) :
            for substr in string.split(line, texcomment2) :
              output.write(substr)
              outputTex.write(substr)
            break
          else :
            for substr in string.split(line, texcomment1) :
              output.write(substr)
              outputTex.write(substr)
            break
      if (linetype == 0) :
        if (oldlinetype > 0) : 
          output.write("\\begin{verbatim}\n")
        output.write(line)
        outputSrc.write(line)
          
    # insert "\end{verbatim}" if needed
    if ( oldlinetype == 0):
      output.write("\\end{verbatim}\n")

    output.close()
    outputSrc.close()
    outputTex.close()
    input.close()

  # If the file is a *.tex file, convert to a file whose name is given
  # in the comment on the first line of the *.tex file

  elif (iniFile == 2) :
    texcomment = languageSpecific[suffix][0][0] + "| "
    linetype = 1
    line = ""
    # read lines one at a time

    while 1:
      oldline = line                # save last line
      oldlinetype = linetype
      line = input.readline()       # read one line at a time
      if not line: break            # break out at end of file

      if ( string.find (line, "\\begin{verbatim}" ) > -1):
        linetype = 0
      elif ( string.find (line, "\\end{verbatim}" ) > -1):
        linetype = 1
      else:
        if ( linetype > 0 ):
          nBlanks = len(line) - len(string.lstrip(line))
          # newline = (nBlanks-1)*' ' + texcomment + string.lstrip(line)
          newline = texcomment + line
          if (newline[len(newline) - 1] != '\012') :
            newline = newline + '\012'
          output.write(newline)
        else:
          output.write(line)

    output.close()
    input.close()
