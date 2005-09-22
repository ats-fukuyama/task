c**********************************************************************
c*****This routine splits an iter database file into 1D or 2D u-files.*
c*****(The dimension is deduced from the file name)                   *
c*****last revision: 12/94 s.e.attenberger, w.a.houlberg, ornl        *
c*****DISCLAIMER: The authors do not assume any liability for         *
c*****deficiencies in this software.  Please report any problems to   *
c*****Stan Attenberger at attenberger@fed.ornl.gov                    *
c*****   ***   This software is licensed for export under  ***        *
c*****   ***   general license agreements GTDA and GUS.    ***        *   
c**********************************************************************
      data nin/22/,nhead/20/,nout/23/
      character dum*50,fname*50,fout*50,foutnm*60
      character line(20)*80,flab*20
      dimension iwrite(20)
10    write(*,*)'  give name of iter database file'
      read(*,'(a50)')fname
c*****Begin to compose output file name(s)
      i1=index(fname,'_')
      dum=fname(i1+1:)
      i2=index(dum,'_')
      if(i1.eq.0 .or. i2.eq.0)then
        write(*,*)'name must have the form dev_shotnum_Nd.dat',
     #            ' where N=1 or 2. '
        write(*,*)'On unix systems, a path may precede the name.'
        write(*,*)'On vms systems, do not type the .dat'
        go to 10
      endif
      fout=fname(1:i1-1)//fname(i1+i2+1:i1+i2+2)//fname(i1+1:i1+i2-1)
c*****Open input file
      open(unit=nin,file=fname,status='old',err=40)
c*****We search for the line containing "-DEPENDENT" in order to
c*****identify the dependent variable.  It may be preceded by pairs of
c*****separator lines (**********), which are discarded,
c*****and an arbitrary number of "associated scalar quantities" lines.
15    continue
      do l=1,nhead
        iwrite(l)=1
        read(nin,'(a80)',end=60) line(l)
        if(index(line(l),'**********').ne.0  .and.
     #    index(line(l-1),'**********').ne.0)        then
c*****    We do not require that the separator flag start in column 1.
          iwrite(l)=0
          iwrite(l-1)=0
        elseif(index(line(l),'-DEPENDENT').ne.0)    then
c*****    Found line identifying dependent variable.
c*****    Note that column 1 is blank, except for separator lines.
          nh=l
          flab=line(l)(2:21)
          go to 20
        endif
      enddo
      write(*,*)' error, the string -DEPENDENT is not in the header'
      stop
20    lenf=index(fout,' ')-1
      foutnm=fout(1:lenf)//'.'//flab
      open(unit=nout,file=foutnm,status='unknown')
c*****dump lines already read.
      do l=1,nh
        if(iwrite(l).eq.1)write(nout,'(a80)')line(l)
      enddo
c*****write remaining lines 1 by 1, checking for pair of separator lines.
      read(nin,'(a80)')line(1)
      do l=1,1000000
        read(nin,'(a80)',end=50)line(2)
        if(index(line(1),'**********').ne.0  .and.
     #    index(line(2),'**********').ne.0)        then 
c*****    Found end of u-file data.
          go to 30
        else
c*****    Write one line
          write(nout,'(a80)')line(1)
          line(1)=line(2)
        endif
      enddo
      write(*,*)'error, do loop max hit.',flab,' data not finished'
      stop
c*****All done with this U-file.
30    close(nout)
      go to 15
      stop
c*****error exit
40    write(*,*)'Can not open file ',fname
      stop
c*****normal exit
50    write(nout,'(a80)')line(1)
60    stop
      end

