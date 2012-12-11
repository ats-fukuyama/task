#/usr/bin/python
#

#
# This is an object version of Doug's prefpp.py and postfpp.py code.  It expects
# the input as list of lines and returns list of lines.
#

import sys
import string
import os
import time

class prefpp:
    """
    #  python version of prefpp (experimental, dmc Jan 2005)
    #
    #  prefpp "encodes" fortran source code so that it can be safely fed through
    #  the gcc C pre-processor, even if it is fixed form source and gcc doesn't
    #  respect whitespace:
    #
    #    (a) comments and quoted literals are protected from gcc by encoding
    #    (b) source code seen by gcc is converted to uppercase, but, original
    #        case is stored as an encoded literal, so that (as much as possible)
    #        the original typesetting of the source can be restored by postfpp
    #        (a 2nd program which decodes that which this prograf prefpp encodes).
    #    (c) cpp lines (which start with "#" are passed through as literals), but,
    #        comments are encoded so that the processed source shows in comments
    #        all cpp directives that were present in the original.
    #
    #  ------- example --------
    ## > cat foon.zork
    #
    ## define MYPROG PROGRAM
    #
    ## > cat foo.f90
    #
    # (strip off the leading "## " to recover the input file)
    #
    ## #include "foon.zork"
    ## 
    ## myprog foo
    ## 
    ## ! a cmt line (Doug's code)
    ## 
    ##   write(6,*) 'this is a MYPROG test'
    ##   call foosub
    ## 
    ##   stop   ! end of prog
    ## 
    ## end myprog foo
    #
    # preprocess the foo.f90 file in python,
    #
    ##f = open("foo.f90",'r')
    ##s = f.readlines()    # grab source in list of lines
    ##f.close()
    ##
    ##pre = prefpp()       # create an instance
    ##   
    ##out = pre(s,free=1)  # preprocess the source code producing a list output
    # 
    # contents of out list,
    #
    ##  ['/&! ! [ ~~P i n c l u d e ~~2 f o o n . z o r k~2 ~ ]&/',
    ## '#include "foon.zork"',
    ## ' /&! ! [ ~~P i n c l u d e ~~2 f o o n . z o r k~2 ~ ]&/',
    ## ' /&b&/', 
    ## ' /&s m y p r o g ~ f o o&/',
    ## 'MYPROG /~FOO',
    ## ' /&b&/ ',
    ## ' /&! ! ~ a ~ c m t ~ l i n e&/',
    ## ' /&b&/ ',
    ## ' /&s ~ ~ w r i t e ( 6 , * ) ~~1 t h i s ~ i s ~ a ~ t e s t~1&/',
    ## ' /~ /~WRITE(6,*) /~ /&1~1 t h i s ~ i s ~ a ~ t e s t~1&/',
    ## ' /&s ~ ~ c a l l ~ f o o s u b&/',
    ## ' /~ /~CALL /~FOOSUB',
    ## ' /&b&/ ',
    ## ' /&s ~ ~ s t o p&/',
    ## ' /~ /~STOP /&! ~ ~ ~ ! ~ e n d ~ o f ~ p r o g&/',
    ## ' /&b&/ ',
    ## ' /&s e n d ~ m y p r o g ~ f o o&/',
    ## 'END /~MYPROG /~FOO']
    #
    # --- start up gcc in a pipe ---
    ##f = os.popen("gcc -E -P -x c foo_pre.out",'r')  # read output from stdout of gcc
    ##gout = f.readlines()
    ##f.close()
    #
    # -- post process --
    ##post = postfpp()              # create an instance of a post processor
    ##fout = post(gout,free=1)      # post process
    #
    #
    # contents of fout list,
    #
    ## ['![ #include "foon.zork" ]',
    ## ' ',
    ## 'PROGRAM foo',
    ## ' ',
    ## '! a cmt line (Doug's code)',
    ## ' ',
    ## '  write(6,*) 'this is a MYPROG test'',
    ## '  call foosub',
    ## ' ',
    ## '  stop   ! end of prog',
    ## ' ',
    ## 'end PROGRAM foo']
    #
    #  note that the "macro" redefining MYPROG as PROGRAM has been applied, and,
    #  the #include reference is commented out, but otherwise the source code is
    #  unchanged.
    #
    #  if gcc were used directly, myprog would not be converted because it is
    #  lowercase, and gcc is case sensitive.  Further, an error would occur in
    #  the comment line with "Doug's code" because being fortran it is not
    #  recognized as a comment by gcc, so it thinks there is an unterminated
    #  quote string in the source.
    #
    #  examples of more practical use -- interchanging "single precision" and
    #  "double precision" BLAS, LAPACK or NAG subroutine calls in code which is
    #  being ported back and forth between CRAY and 32 bit workstation machines.
    #  Case sensitivity makes this tough for gcc alone...
    #
    #  By this method, all gcc functionality is available in a "safe way" for
    #  fortran programmers-- this is more powerful and more standard than the
    #  ".F" file preprocessing supported in various ways by various fortran
    #  vendors.
    #
    """
    def clear(self):
        """
        restore to initial state
        """
        self.ilen = 0
        self.ifnb = 0
        self.line = ""
       
        self.freeform = 0
        self.const_conv = 0
        self.maxlen = 72
        self.infile = ""
        
        self.f77exp = "D"
        
        self.oline = ""
        self.nquote=0
        self.out = []       # processed output
    
    
    def mylen(self):
        """
        get length of line; trailing whitespace was removed.  Blank lines
        give zero length.  For nonblank lines return position of first nonblank
        character.
        """
        self.ilen=len(self.line)
    
        itab = 0  # flag if tab before 1st non-whitespace character
        self.ifnb = 0
        isig = 0
        while self.ifnb < self.ilen and isig == 0:
            if string.find(string.whitespace,self.line[self.ifnb]) <> -1:
                if self.line[self.ifnb] == '\t':
                    itab = 1
                self.ifnb += 1
            else:
                isig=1
            
        return itab

    
    def prefpp_cmt(self,itab):
        """
        determine if line of source text is blank or comment or cpp line
        """
        icmt=self.ilen
        icpp=self.ilen
        if self.ilen == 0:
            iblank=1
            return [iblank,icmt,icpp]  # blank line
        
        # nonblank line...
        
        iblank=0
        if self.line[self.ifnb] == '#':
            isig = 0
            if self.freeform == 0:
                if self.ifnb == 1 and itab == 1: isig = 1
                if self.ifnb == 5 and itab == 0: isig = 1
            if isig == 1:
                self.line = self.line[0:self.ifnb] + ">" + self.line[self.ifnb+1:] 
                self.line = self.line + ' ! prefpp warning: # continuation '
                self.line = self.line + 'marker, col. 6'
                itab = self.mylen()
                icmt=self.ilen
                icpp=self.ilen
            else:
                icpp=self.ifnb
                return [iblank,icmt,icpp]  # cpp line, starts with "#"
            
        if self.freeform == 0 and self.ifnb == 0:
            if self.line[self.ifnb] == 'c' or self.line[self.ifnb] == 'C': icmt=0
            if self.line[self.ifnb] == 'd' or self.line[self.ifnb] == 'D': icmt=0
            if self.line[self.ifnb] == '!': icmt=0
            if self.line[self.ifnb] == '*': icmt=0
        elif self.line[self.ifnb] == '!':
            #  first non-blank is a comment character... so whole line is a comment
            #  unless fixed form f77 code is using "!" as a continuation character
            icmt = 0
            if self.freeform == 0:
                if self.ifnb == 1 and itab == 1: icmt = 0
                if self.ifnb == 5 and itab == 0: icmt = 0
                
        if icmt == 0:
            return [iblank,icmt,icpp]  # entire line is a comment
    
        # at least part of line is fortran non-comment and non-cpp
    
        kquote = 0
        if self.line[self.ifnb] == "'": kquote = 1
        if self.line[self.ifnb] == '"': kquote = 2
    
        ilnbp = self.ifnb
        ilnb = self.ifnb
        for i in range(self.ifnb+1,self.ilen):
            if kquote > 0:
                if kquote == 1 and self.line[i] == "'":  kquote = 0
                if kquote == 2 and self.line[i] == '"':  kquote = 0
                if kquote == 0:
                    ilnbp = i
                    ilnb = i
            else:
                if self.line[i] == "!":
                    #  inline comment found; absorb preceding whitespace and "&"...
                    if self.line[ilnb] == "&":  ilnb = ilnbp
                    icmt = ilnb + 1
                    return [iblank,icmt,icpp]
                elif self.line[i] == "'":  kquote = 1
                elif self.line[i] == '"':  kquote = 2
                elif self.line[i] <> " " and self.line[i] <> "\t":
                    ilnbp = ilnb
                    ilnb = i
    
        # fall through: no comment; flag ampersand though
    
        if self.line[ilnb] == "&":
            icmt = ilnbp + 1
    
        return [iblank,icmt,icpp]
    
    def prefpp_const_conv(self):
        """
        -------------------------------------------------------------
        replace construct CONST(<mantissa>,<exp>) with
         <mantissa>[D|E]<exp>
         
        this implements a floating point expansion macro that fails under some
        versions of gcc cpp...
        """
    
        legal_mm = string.digits + "." + "+" + "-"  # legal mantissa characters
        legal_ee = string.digits + "+" + "-"        # legal exponent characters
    
        ileno = len(self.oline)
        ilast = ileno - 1
        ic = string.find(self.oline,'CONST(')
        while ic > -1:
            if ic + 10 > ileno: break     # not enough room for CONST(m,e) macro
    
            ip1 = ic  + 5                 # first parenthesis
            im1 = ip1 + 1
            ip2 = ic + string.find(self.oline[ic:],')')
            if ip2 < 0: break             # ")" not found.
    
            nump = string.split(self.oline[ip1+1:ip2],',')
            if len(nump) <> 2: break      # m,e form not inside parens
    
            mm = string.strip(nump[0])
            ee = string.strip(nump[1])
            
            if mm == "": break            # non-null mantissa required.
            if ee == "": break            # non-null exponent required.
    
            ichk = 0
            for i in range(len(mm)):
                if string.find(legal_mm,mm[i]) < 0: ichk = ichk + 1
            if ichk > 0: break            # illegal character in mantissa
    
            ichk = 0
            for i in range(len(ee)):
                if string.find(legal_ee,ee[i]) < 0: ichk = ichk + 1
            if ichk > 0: break            # illegal character in exponent
    
            #  final check...
    
            try:
                iee = string.atoi(ee)
            except ValueError:
                ichk = ichk + 1
    
            try:
                zmm = string.atof(mm)
            except ValueError:
                ichk = ichk + 1
    
            if ichk > 0: break            # decode failed
    
            newnum = mm + self.f77exp + ee
    
            if ip2 == ileno:
                self.oline = self.oline[:ic] + " " + newnum
            else:
                self.oline = self.oline[:ic] + " " + newnum + self.oline[ip2+1:]
    
            # update info on oline -- re-enter loop to process more than one CONST
            # per line...
            
            ileno = len(self.oline)
            ilast = ileno - 1
            ic = string.find(self.oline,'CONST(')
            
        return 0
            
    def prefpp_out(self,icont):
        """
        write a non-empty line of output, possibly with continuation
        clear line after written
        """
        if len(self.oline) == 0:
            return 0
    
        if icont ==1:
            self.oline = self.oline + "&-"
    
        if self.const_conv == 1:
            self.prefpp_const_conv()
    
        self.out.append(self.oline)
        self.oline = ""
    
        return 0
    
    def prefpp_cgen(self,gchar,str):
        """
        transcribe a comment...
        gchar = "!" means a true comment; gchar = "s" means a
                transcription of source to preserve case etc.
                on reconstruction
        """
        self.oline = self.oline + " /&" + gchar
        if len(self.oline) > 75: self.prefpp_out(1)
    
        for i in range (0,len(str)):
            if str[i] == '/':
                self.oline = self.oline + '~S'
            elif str[i] == '#':
                self.oline = self.oline + '~P'
            elif str[i] == "'":
                self.oline = self.oline + '~1'
            elif str[i] == '"':
                self.oline = self.oline + '~2'
            elif str[i] == ' ':
                self.oline = self.oline + ' ~'
            elif str[i] == '\t':
                self.oline = self.oline + '~^'
            elif str[i] == '~':
                self.oline = self.oline + '~~'
            else:
                self.oline = self.oline + ' ' + str[i]
    
            if len(self.oline) > 75: self.prefpp_out(1)
    
        self.oline = self.oline + '&/'
        return 0

    def prefpp_sgen(self,icmt):
        """
        transcribe a source line with possible trailing comment
        """
    
        iquote=-1
        if self.nquote > 0: iquote = 0
        
        for i in range(icmt):
            jquote=0
            if self.line[i] == "'": jquote = 1
            if self.line[i] == '"': jquote = 2
            if self.nquote > 0:
                if self.nquote == jquote:
                    # end of quote string; process...
                    self.nquote = 0
                    self.prefpp_cgen(str(jquote),self.line[iquote:i+1])
            else:
                if jquote > 0:
                    # start of quote string
                    self.nquote = jquote
                    iquote = i
                else:
                    if self.line[i] == '/':
                        self.oline = self.oline + ' /&S&/'
                    elif self.line[i] == ' ':
                        self.oline = self.oline + ' /~'
                    elif self.line[i] == '~':
                        self.oline = self.oline + '~~'
                    elif self.line[i] == '\t':
                        self.oline = self.oline + '~^'
                    elif self.line[i] == '#':
                        self.oline = self.oline + '~P'
                    else:
                        achar = self.line[i]
                        achar = string.upper(achar)
                        self.oline = self.oline + achar
    
        # pick up partial quote
        if self.nquote > 0: self.prefpp_cgen(str(self.nquote),self.line[iquote:icmt])
    
        # pick up trailing comment
        if icmt < self.ilen: self.prefpp_cgen('!',self.line[icmt:])
    
        self.prefpp_out(0)
    
        return 0

    
    def prefpp_p1(self):
        """
        process one line of source text
        """
        itab = self.mylen()
        [iblank,icmt,icpp] = self.prefpp_cmt(itab)
    
        if iblank == 1:
            self.oline=" /&b&/ "
        elif icpp < self.ilen:
            # cpp line; check for "fp_const.h"
            if string.find(self.line,'#include') > -1:
                if string.find(self.line,'fp_const.h"') > -1:
                    self.const_conv = 1

            if string.find(self.line,'#define') > -1:
                if string.find(self.line,'CONST') > -1:
                    self.const_conv = 1
            
            self.prefpp_cgen("!","![ " + self.line + " ]")
            self.prefpp_out(0)
            self.oline = self.line
            self.prefpp_out(0)
            self.prefpp_cgen("!","![ " + self.line + " ]")
        else:
            if icmt == 0:
                # pure comment
                self.prefpp_cgen("!",self.line)
            else:
                self.prefpp_cgen("s",self.line[:icmt])
                self.prefpp_out(0)
                self.prefpp_sgen(icmt)
    
        self.prefpp_out(0)
        return 0

    def __call__(self,flines,free=0,f77exp="D",maxlen=72):
        """
        Read in the list of fortran lines and return them prepared for preprocessing as a new list.
          flines = list of fortran code
          free   = nonzero for free form f90
          f77exp = expression to use for exponent
          maxlen = maximum length of code line
        """
        self.clear()
        
        self.freeform = free
        self.f77exp   = f77exp
        self.maxlen   = maxlen
    
        for aline in flines:
            self.line = string.rstrip(aline)
            self.prefpp_p1()

        return self.out


#
#-------------------------------- postfpp --------------------------------------------
#
class postfpp:
    """
    #  python version of postfpp (experimental, dmc Jan 2005) -- see also
    #  prefpp.py ... (transcribed from codesys/source/postfpp.f90)
    #
    #  postfpp "decodes" prefpp/gcc processed fortran source code so that it can
    #  be fed safely to a fortran compiler.
    #
    #  for more information please see the prefpp.py comments...
    """

    def __init__(self):
        self.clear()
        

    def clear(self):
        """
        initialize state
        """
        self.freeform = 0
        self.maxlen = 72
    
        self.line = ""
        self.linep = ""
        self.saveline = ""
    
        self.include_state = 0
    
        self.iline = -1  # current input line
        self.inptr = -1  # current character within input line
        self.out   = []  # output lines
    
    
    def postfpp_dcod(self):
        """
        this procedure decodes a comment or literal string
        that starts at [inptr:]; the decoded string is returned and
        inptr is updated to point past the end in the input line.
        """
        
        ostr = ""
        i = self.inptr -2
    
        while 1:
            i = i + 2
            if i >= len(self.line):
                raise RuntimeError("postfpp.py: postfpp_dcod internal error #1., line->\n"+self.line)
            if self.line[i:i+2] == '&-':
                # continuation-- read next line (which must exist)
                self.iline = self.iline + 1
                if self.iline == len(self.flines):
                    raise RuntimeError("postfpp.py: premature EOF in postfpp_dcod.")
                self.line = string.rstrip(self.flines[self.iline])
                i = -2
            elif self.line[i:i+2] == '&/':
                # end of decode
                self.inptr = i+1  # needs to point to last character decoded...
                break
            elif self.line[i] == ' ':
                if self.line[i+1] == '~':
                    ostr = ostr + " "
                else:
                    ostr = ostr + self.line[i+1]
            elif self.line[i] == '~':
                if self.line[i+1] == 'S':
                    ostr = ostr + "/"
                elif self.line[i+1] == 'P':
                    ostr = ostr + "#"
                elif self.line[i+1] == '1':
                    ostr = ostr + "'"
                elif self.line[i+1] == '2':
                    ostr = ostr + '"'
                elif self.line[i+1] == '~':
                    ostr = ostr + "~"
                elif self.line[i+1] == '^':
                    ostr = ostr + "\t"
                else:
                    raise RuntimeError("postfpp.py: postfpp_dcod internal error #2.")
            else:
                raise RuntimeError("postfpp.py: postfpp_dcod internal error #3.")
    
        return ostr
            
    
    def postfpp_src(self,model,umodel):
        """
        this procedure regenerates a source line, using model information
        as provided.  Trailing comments are also appended.  Regenerated
        source line with trailing comments is returned.  This may involve
        reading through multiple lines of the input file.
        """
    
        #  model -- orig. source line
        #  umodel -- orig. source line converted to uppercase except
        #            for literal strings occurring within the source line
    
    
        self.iline = self.iline + 1
        if self.iline == len(self.flines):
            s = " ?postfpp_src: premature EOF.\n  model source line: "+str(model)
            raise RuntimeError(s)
        self.line = string.rstrip(self.flines[self.iline])
        self.inptr = -1
    
        ostr = ""
        j = 0
        
        cont_quote = 0      # nquote from start of line
        while 1:
            self.inptr = self.inptr + 1
            if self.inptr == len(self.line): break
    
            if self.line[self.inptr] == " ":
                iskip = 0
                if self.inptr + 1 < len(self.line):
                    if self.line[self.inptr+1] == " ":
                        # skip consecutive blanks: presumed cpp artifact
                        iskip=1
                if iskip: continue
    
                if self.line[self.inptr:self.inptr+2] == ' /':
                    if self.line[self.inptr:self.inptr+3] == ' /~':
                        # encoded blank
                        ostr = ostr + " "
                        j = j + 1
                        self.inptr = self.inptr+2
    
                    elif self.line[self.inptr:self.inptr+4] == ' /&1' or self.line[self.inptr:self.inptr+4] == ' /&2':
                        # quoted string literal
                        if self.inptr==0 and len(self.line)>=6:
                            # is this a continuation of a previous line's quote?
                            xx = self.line[self.inptr+4:self.inptr+5]
                            if xx!='~1' and xx!='~2' :
                                if self.line[self.inptr:self.inptr+4] == ' /&1':
                                    cont_quote=1   # continue single quote
                                else:
                                    cont_quote=2   # continue double quote
                        self.inptr = self.inptr + 4
                        ostr = ostr + self.postfpp_dcod()
                        j = len(ostr)
    
                    elif self.line[self.inptr:self.inptr+4] == ' /&!':
                        # trailing comment
                        self.inptr = self.inptr + 4
                        ostr = ostr + self.postfpp_dcod()
                        break   # ostr[0:j] -- source; ostr[j:] -- trailing cmt
    
                    elif self.line[self.inptr:self.inptr+6] == ' /&S&/':
                        # encoded slash
                        ostr = ostr + "/"
                        j = j + 1
                        self.inptr = self.inptr + 5
                        
                    else:
                        s = self.line + "\n at position: "+str(self.inptr)+"\n postfpp_src decode error #1."
                        raise RuntimeError(s)
                else:
                    continue  # extraneous blank, presumed cpp artifacte
    
            elif self.line[self.inptr] == "~":
                if self.line[self.inptr+1] == "~":
                    # encoded tilda
                    ostr = ostr + "~"
                    j = j + 1
                    self.inptr = self.inptr + 1
    
                elif self.line[self.inptr+1] == "^":
                    # encoded tab
                    ostr = ostr + "\t"
                    j = j + 1
                    self.inptr = self.inptr + 1
    
                elif self.line[self.inptr+1] == "P":
                    # encoded #-sign
                    ostr = ostr + "#"
                    j = j + 1
                    self.inptr = self.inptr + 1
    
                else:
                    s = self.line+"\n at position: "+str(self.inptr)+"\n postfpp_src decode error #1."
                    raise RuntimeError(s)
    
            else:
                # ordinary non-blank:
                j = j+1
                ostr = ostr + self.line[self.inptr]
    
        # if uppercase conversion matches reconstructed source, then
        # there was no change, return the original source line...
        
        if ostr[0:j] == umodel:
            if j == len(ostr):
                ostr = model
            else:
                ostr = model + ostr[j:]
    
            return ostr
        
        # if it does not match, then there were cpp changes; restore
        # line to matching case for characters before the first change
        # and after the last change only...
    
##         print ' umodel  = "',umodel,'"'
##         print ' model   = "',model,'"'
##         print ' ostr[:j]= "',ostr[:j],'"'
        
        lhs = -1
    
        while 1:
            lhs = lhs + 1
            if ostr[lhs] <> umodel[lhs]: break
    
        rhs = j
        rhsu = len(umodel)
    
        while 1:
            rhs = rhs - 1
            rhsu = rhsu - 1
            if ostr[rhs] <> umodel[rhsu]: break
    
        nstr = ""
        if lhs > 0:
            nstr = nstr + model[0:lhs]
    
        rhs = rhs+1
        rhsu = rhsu+1
        
        nstr = nstr + ostr[lhs:rhs]
        if rhsu < len(umodel):
            nstr = nstr + model[rhsu:]
    
        # and include trailing comment, if any...
        if j < len(ostr):
            nstr = nstr + ostr[j:]

        if self.freeform or j <= self.maxlen:
            pass
        else:
            # attempt repair of cpp caused source line length error in fixed form
            # source: first try to squeeze out redundant whitespace...
    
            ostr = nstr
            nstr = ""
            
            self.nquote = cont_quote
            ilnb = -1
            itab = 0
            for i in range(len(ostr)):
                iblank = 0
                iskip = 0
                if ostr[i] == '\t': itab = itab + 5
                if ostr[i] <> '\t' and ostr[i] <> ' ': ilnb = i
                if self.nquote > 0:
                    if self.nquote == 1:
                        if ostr[i] == "'":  self.nquote = 0
                    if self.nquote == 2:
                        if ostr[i] == '"':  self.nquote = 0
                else:
                    if ostr[i] == "'":  self.nquote = 1
                    if ostr[i] == '"':  self.nquote = 2
                    if ostr[i] == ' ' or ostr[i] == '\t':
                        iblank = 1
                        if ilnb+itab >= 5 and j > self.maxlen and i < j:
                            if ostr[i+1] == ' ' or ostr[i+1] == '\t':
                                iskip = 1
                                j = j - 1
    
                if iskip:
                    pass
                else:
                    if iblank:
                        nstr = nstr + ' '
                    else:
                        nstr = nstr + ostr[i]
        
        if self.freeform or j <= self.maxlen:
            pass
        else:
            # compression insufficient; now look for blank after non-blanks or
            # start of literal as line break opportunity...
            
            self.nquote = cont_quote
            ilnb = -1
            k = 0
            
            for i in range(self.maxlen):
                if nstr[i] <> " " and nstr[i] <> "\t":  ilnb=i
                if self.nquote > 0:
                    if self.nquote == 1:
                        if nstr[i] == "'":  self.nquote = 0
                    if self.nquote == 2:
                        if nstr[i] == '"':  self.nquote = 0
                else:
                    if nstr[i] == "'":  self.nquote = 1
                    if nstr[i] == '"':  self.nquote = 2
                    if self.nquote > 0 and i > 0:
                        if nstr[i] <> nstr[i-1]:  k = i  # start of literal
                    if nstr[i] == ' ' or nstr[i] == '\t':
                        if ilnb > 9:  k = i   # available blank
                    if nstr[i] == '*' and ilnb > 9:
                        if nstr[i-1] <> '*' and nstr[i+1] <> '*':  k = i # multiply
                    if nstr[i] == '/' and ilnb > 9:
                        if nstr[i-1] <> '/' and nstr[i+1] <> '/':  k = i # divide
    
            if k < 10:
                nstr = nstr + ' !(postfpp: ?line-length error - repair failed.)'
            else:
                self.out.append(nstr[0:k] + '   ! (postfpp: inserted line break)')
                nstr = '     >   ' + nstr[k:]
                
        return nstr
    

    def postfpp_p1(self):
        """
        this procedure processes a single line of prefpp|gcc output
        """
        if self.line == " /&b&/":
            # blank line in original source...
            self.out.append(" ")
    
        elif self.line[0:4] == " /&!":
            # encoded comment
            self.inptr = 4
            ostr = self.postfpp_dcod()
            self.out.append(ostr)
    
        elif self.line[0:4] == " /&s":
            # encoded source, followed by a source line
            self.inptr = 4
            model = self.postfpp_dcod()
            modeu = string.upper(model)
            umodel = ""
            iquote = 0
            self.nquote = 0
            for i in range(len(model)):
                if self.nquote > 0:
                    iendq = 0
                    if self.nquote == 1 and model[i] == "'": iendq = 1
                    if self.nquote == 2 and model[i] == '"': iendq = 1
                    if iendq:
                        umodel = umodel + model[iquote:i+1]
                        iquote = i+1
                        self.nquote = 0
                else:
                    if model[i] == "'": self.nquote = 1
                    if model[i] == '"': self.nquote = 2
                    if self.nquote > 0:
                        umodel = umodel + modeu[iquote:i+1]
                        iquote = i+1
            if iquote < len(model):
                if self.nquote > 0:
                    umodel = umodel + model[iquote:]
                else:
                    umodel = umodel + modeu[iquote:]
    
            ostr = self.postfpp_src(model,umodel)
            self.out.append(ostr)
            
        else:
            if self.include_state:
                self.out.append(self.line)
            else:
                raise RuntimeError("postfpp_p1: internal error #1\n at input line: "+self.line)
    
        return 0

    def __call__(self, flines, free=0,f77exp="D",maxlen=72):
        """
        Read in the list of preprocessed lines and return them as fortran lines in a new list.
          flines = list of preprocessed code
          free   = nonzero for free form f90
          maxlen = maximum length of code line
        """
        self.clear()
        self.freeform = free
        self.maxlen   = maxlen
        self.flines   = flines
        
        while 1:
            if self.saveline <> "":
                self.line = self.saveline
                self.saveline = ""
            else:
                self.iline = self.iline + 1
                if self.iline == len(self.flines): break   # at EOF
                self.line = string.rstrip(self.flines[self.iline])
                if self.line == "": continue          # skip gcc-generated blank lines
            
            ifind = string.find(self.line,' /&! ! [ ~')
            if ifind > -1:
                # set include_state based on comments recording #include stmts
                # when reaching end -- take care of fact that gcc might not
                # keep whitespace, so, the marker could happen in the middle of
                # a line
                for i in range(ifind+10,len(self.line)-1,2):
                    if self.line[i:i+2] <> " ~":
                        if self.line[i:i+2] == "~P":
                            if self.line[i+2:i+16] == ' i n c l u d e':
                                # #include stmt record found...
                                if self.include_state:
                                    if ifind > 0:
                                        self.saveline = self.line[ifind:]
                                        self.line = self.line[:ifind]
                                    else:
                                        self.include_state = 0
                                else:
                                    if ifind == 0:
                                        self.include_state = 1
    
                        # exit after non " ~" pair is found...
                        break
                    
            # this scary looking code prevents duplicated generated cpp
            # comment lines from being copied; also, clear state flag for
            # detected #include files which may insert verbatim lines
            # which do not want postfpp decoding!
    
            icancel = 0
            if self.line == self.linep:
                if self.line[0:10] == ' /&! ! [ ~':
                    for i in range(10,len(self.line)-1,2):
                        if self.line[i:i+2] <> " ~":
                            if self.line[i:i+2] == "~P":
                                icancel = 1
    
                            # exit after non " ~" pair is found...
                            break
    
            self.linep = self.line
            if icancel == 0:
                self.postfpp_p1()
            else:
                while 1:
                    # skip continuations of canceled line
                    if self.line[-2:] == '&-':
                        self.iline = self.iline + 1
                        if self.iline == len(self.flines): break  # at EOF
                        self.line = string.rstrip(self.flines[self.iline])
                    else:
                        break
                if self.iline == len(self.flines): break  # at EOF -- go all the way out.
                            
        # end of processing...
        return self.out



def fppFile(fileIn, fileOut, fppdefs, free=0, maxlen=72, fileTemp=None, fppcom=None, info=0, keep=0):
    """
    Preprocess a file into another file.
      fileIn   = source file
      fileOut  = name of resulting preprocessed file
      fppdefs  = preprocessor directives such as "-D_Linux -I/usr/loca/lib"
      free     = nonzero if this is a free format fortran file
      maxlen   = maximum number of characters in a source line
      fileTemp = name of temporary file or None to invent one, this is deleted after preprocessing
      fppcom   = command for executing the preprocessor or None to use the gcc default
      info     = nonzero for info messages
      keep     = nonzero to keep temporary file
    """
    if (not fppcom):
        fppcom = "gcc -E -P -x c"      # OSX might need -traditional-cpp flag
        if (os.environ.get("WORKSTATION_TYPE","")=="OSX"):
            fppcom += " -xassembler-with-cpp"

    if (info):
        print "   -> fileIn      ", fileIn
        print "   -> fileOut     ", fileOut
        print "   -> free:       ", free
        print "   -> maxlen:     ", maxlen
        print "   -> fppdefs:    ", fppdefs
        print "   -> fppcom:     ", fppcom
        
    if (not fileTemp):
        basename = os.path.basename(fileIn)
        tmpdir   = string.strip(os.environ.get('TMPDIR',""))
        if (len(tmpdir)==0):
            tmpdir = string.strip(os.environ.get('TEMP',""))
        if (len(tmpdir)==0):
            tmpdir = string.strip(os.environ.get('TMP',""))
        if (len(tmpdir)==0):
            tmpdir = "/tmp"
        # create a unique name with basename prefix in TMPDIR or somewhere
        myid = "_%f" % time.time()
        fileTemp = os.path.join(tmpdir,basename+myid)   
        if (info):
            print "   -> temp file:  ", fileTemp
            
    if (os.path.isfile(fileOut)):
        os.unlink(fileOut)             # remove previous preprocessed file

    f = open(fileIn,'r')
    slines = f.readlines()             # grab source in list of lines
    f.close()

    pre = prefpp()
    out = pre(slines, free=free, maxlen=maxlen)
    f = open(fileTemp,'w')
    
    status = 0      # return status of preprocessor
    fcmd   = ""     # command used for preprocessing

    for x in out:
        f.write(x+'\n')
    f.close()

    # -- preprocess through gcc --
    fcmd = "%s %s %s"%(fppcom,fppdefs,fileTemp)
    if (info):
        print "   -> preprocess: ",fcmd
    f = os.popen(fcmd,'r')                   # read output from stdout of gcc
    gout = f.readlines()                     # grab lines of preprocessed output
    #for x in gout:
    #    print x,
    status = f.close()                       # will be None if no error
    if (info):
        print "   -> status:     ",status

    # The following seems to fail occasionally with the message that tempFile does
    # not exist. I believe this could occur whenever the file is cleaned up before
    # preprocessing step is complete. To prevent this the pipe should be closed
    # before cleaning up. Removed the try/finalize block. A. Pletzer April 2009
    
    if (os.path.isfile(fileTemp) and keep==0):
        os.unlink(fileTemp)                  # clean up intermediate file
    

##     try:
##         for x in out:
##             f.write(x+'\n')
##         f.close()
        
##         # -- preprocess through gcc --
##         fcmd = "%s %s %s"%(fppcom,fppdefs,fileTemp)
##         if (info):
##             print "   -> preprocess: ",fcmd
##         f = os.popen(fcmd,'r')                   # read output from stdout of gcc
##         gout = f.readlines()                     # grab lines of preprocessed output
##         #for x in gout:
##         #    print x,
##         status = f.close()                       # will be None if no error
##         if (info):
##             print "   -> status:     ",status
##     finally:
##         if (os.path.isfile(fileTemp) and keep==0):
##             os.unlink(fileTemp)                  # clean up intermediate file
    if (status):
        raise RuntimeError("Error running preprocessor on %s, command = %s"%(fileIn,fcmd))
    
    # -- post-preprocess --
    post = postfpp()
    fout = post(gout, free=free, maxlen=maxlen)
    f = open(fileOut,'w')
    for x in fout:
        f.write(x+'\n')      # write preprocessed source file
    f.close()


#
# ========================== testing ==========================
#
if __name__ == "__main__":
    #
    # run the test example
    #
    import os, pprint

    pp = pprint.PrettyPrinter(indent=4)
    
    #
    # -- create files foo.f90, foon.zork --
    #
    f = open("foo.f90",'w')
    
    f.write("""
#include "foon.zork"

myprog foo

! a cmt line (Doug's code)

  write(6,*) 'this is a MYPROG test'
  call foosub
  call BOWL(3.14159)

  stop   ! end of prog

end myprog foo
""")
    f.close()

    f = open("foon.zork",'w')
    f.write("""
#define MYPROG PROGRAM
""")
    f.close()

    # -- read in source file --
    f = open("foo.f90",'r')
    s = f.readlines()
    f.close()

    print " "
    print "input -> "
    pp.pprint(s)

    # -- preprocess --
    pre = prefpp()

    out = pre(s,free=1)

    print " "
    print "preproceesor out -> "
    pp.pprint(out)

    f = open("foo_pre.out",'w')
    for x in out:
        f.write(x+"\n")
    f.close()

    # --- start up gcc in a pipe to read the output ---
    f = os.popen("gcc -E -P -x c foo_pre.out",'r')
    gout = f.readlines()
    f.close()

    print " "
    print "gcc preprocessor out -> "
    pp.pprint(gout)
    
    # -- postprocess --
    post = postfpp()

    fpost = post(gout,free=1)
    
    print " "
    print "post processing out -> "
    pp.pprint(fpost)
    
    f = open("foo_post.f90",'w')
    for x in fpost:
        f.write(x+"\n")
    f.close()
                                   
    print " ------- fppFile test ---------"
    fppFile("foo.f90","foo_fpp.f90","-I. -DBOWL=dogdish",
            free=1, maxlen=132, info=1)

    f = open("foo_fpp.f90",'r')
    pout = f.readlines()
    f.close()
    for x in pout:
        print x,

                                   
    print " ------- Bad fppFile test should fail ---------"
    caught = 0
    try:
        fppFile("foo.f90","foo_fpp.f90","-DBOWL=dogdish",
                free=1, maxlen=132, info=1)
    except RuntimeError, e:
        caught = 1
        print "caught = ",e

    if (not caught):
        raise RuntimeError("expected 'Bad fppFile test' to fail because there is no -I. argument")
    print "finished ok"

