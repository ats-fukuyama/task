!C
!C
!C
!C
!C
!C
!C
!C
!C 
!C 
MODULE T2WRIT
  USE T2CNST, ONLY:&
       i0ikind,i0rkind
  IMPLICIT NONE

  PUBLIC T2WRIT_MAIN
  PUBLIC T2_WRIT_GPT
  PUBLIC T2_WRIT_GP1
  PRIVATE
CONTAINS
 
  SUBROUTINE T2WRIT_MAIN(d1xxx,i0tag,c10tag)
    USE T2COMM,ONLY:&
         i0fnum,i0spcs,i0nmax0,i0nmax3,i0xmax,i0vmax,&
         i1pdn2,i1rdn2,i0emax,d2rzc3,i3enr,&
         d0rmjr,i0lmax,i0stm2,i1mlvl,&
         !C
         d1bp3,d1bt3,d1er3,d1ep3,d1et3,&
         d2n3, d2fr3,d2fb3,d2ft3,&
         d2p3, d2qr3,d2qb3,d2qt3,&
         d1jm1,d1jm2,d1jm3,d1jm4,d1jm5

    REAL(   i0rkind),DIMENSION(1:i0xmax),INTENT(IN)::d1xxx
    CHARACTER(   10),                    INTENT(IN)::c10tag
    INTEGER(i0ikind),                    INTENT(IN)::i0tag
    
    
    INTEGER(i0ikind)::&
         i1,i2,j2,i0nnc, i0csz, i0ctype,&
         i0rdn2,i0pdn2,i0ofst,i0mpt,i0ecnt,&
         i0mlva,i0mlvb,i0mlvc,i0vnumb,i0fnumb
    REAL(   i0rkind)::d0qf

    CHARACTER(40)::c40fname,c40tname,c40lname
    CHARACTER(10)::c10fnumb
    CHARACTER( 2)::c2vnumb

100 FORMAT(A26)
110 FORMAT(A40)
120 FORMAT(A5)
130 FORMAT(A25)
140 FORMAT(A6,I7,1X,A5)
150 FORMAT(3E15.6)
160 FORMAT(A5,I7,I7)
170 FORMAT(4I7)
175 FORMAT(5I7)
176 FORMAT(6I7)
180 FORMAT(A10,I7)
190 FORMAT(I7)
200 FORMAT(A20)
210 FORMAT(E15.6)
       
    WRITE(c10fnumb,'(I10)')i0tag
    
    c10fnumb = ADJUSTL( c10fnumb)
    i0fnumb  = LEN_TRIM(c10fnumb)
    SELECT CASE(i0fnumb)
    CASE(1)
       c10fnumb = '00000'//c10fnumb
       c10fnumb =  ADJUSTL(c10fnumb)
    CASE(2)
       c10fnumb =  '0000'//c10fnumb
       c10fnumb =  ADJUSTL(c10fnumb)
    CASE(3)
       c10fnumb =   '000'//c10fnumb
       c10fnumb =  ADJUSTL(c10fnumb)
    CASE(4)
       c10fnumb =    '00'//c10fnumb
       c10fnumb =  ADJUSTL(c10fnumb)
    CASE(5)
       c10fnumb =     '0'//c10fnumb
       c10fnumb =  ADJUSTL(c10fnumb)
    CASE(6)
       c10fnumb =         c10fnumb
       c10fnumb = ADJUSTL(c10fnumb)
    CASE default
       WRITE(6,*)'ERROR: ILLEAGAL TIMESTEP NUMBER'
       STOP
    END SELECT
    
    !C
    !C WRITE RESULTS
    !C
    
    c40fname = TRIM(c10tag) //'_' // TRIM(c10fnumb) //'.vtk'
    c40tname = 'TASK/T2 RESULTS (VTK format.ver)'

    !C OPEN OUTPUT FILE
    OPEN(i0fnum, file=TRIM(c40fname), status='unknown')

    !C WRITE HEADER
    WRITE(i0fnum, 100)'# vtk DataFile Version 2.0'
    WRITE(i0fnum, 110) c40tname
    WRITE(i0fnum, 120)'ASCII'
    WRITE(i0fnum, 130)'DATASET UNSTRUCTURED_GRID'

    !C WRITE NUMBER OF POINTS AND DATATYPE
    WRITE(i0fnum,140)'POINTS', i0nmax3, 'float'
    
    !C WRITE COODINATES OF EACH POINTS 
    !C 0    , R-R_0      , PHI      ,Z
    !C ................
    !C ................
    WRITE(i0fnum,150)(&
         d2rzc3(i1,1)-d0rmjr,d2rzc3(i1,2),0.D0,i1=1,i0nmax3)
    
    !C WRITE ELEMENT - NODE RELATIONS 
    i0csz= 5*i0emax-i1pdn2(1)+i1pdn2(1)*(2**(i1mlvl(i0lmax)-1)-1)
    
    WRITE(i0fnum,160)'CELLS', i0emax, i0csz
    
    i0ecnt=0
    
    DO i1=1,i0lmax
       i0mlva=i1mlvl(i1-1)
       i0mlvb=i1mlvl(i1)
       i0mlvc=i1mlvl(i1+1)
       i0rdn2=i1rdn2(i1)
       i0pdn2=i1pdn2(i1)
       i0ofst= 0
       IF(i0mlvb.GE.2)THEN
          DO i2=0,i0mlvb-2
             i0ofst = i0ofst+i1pdn2(1)*(2**i2)
          ENDDO
       END IF
       DO i2=1,i0rdn2
       DO j2=1,i0pdn2
          i0ecnt= i0ecnt+1
          IF(    (i0mlva.EQ.0).AND.(i2.EQ.1))THEN
             i0nnc = 3
             WRITE(i0fnum,170)&
                  i0nnc,&
                  i3enr(i0ecnt,2,1)-1,i3enr(i0ecnt,2,2)-1,&
                  i3enr(i0ecnt,2,3)-1
          ELSEIF((i0mlvb.NE.i0mlvc).AND.&
               (i0mlvc.NE.0).AND.&
               (i2.EQ.i0rdn2))THEN
             i0nnc = 5
             i0mpt = i0stm2+i0ofst+j2
             WRITE(i0fnum,176)&
                  i0nnc,&
                  i3enr(i0ecnt,2,1)-1,i3enr(i0ecnt,2,2)-1,i0mpt-1,&
                  i3enr(i0ecnt,2,3)-1,i3enr(i0ecnt,2,4)-1
          ELSE
             i0nnc = 4
             WRITE(i0fnum,175)&
                  i0nnc,&
                  i3enr(i0ecnt,2,1)-1,i3enr(i0ecnt,2,2)-1,&
                  i3enr(i0ecnt,2,3)-1,i3enr(i0ecnt,2,4)-1
          ENDIF
       ENDDO
       ENDDO
    ENDDO
    
    !C WRITE CELL TYPE: i0ctype
    !C  5: VTK TRIANGLE
    !C  9: VTK QUAD
    WRITE(i0fnum,180) 'CELL_TYPES',i0emax
    DO i1=1,i0lmax
       i0rdn2=i1rdn2(i1)
       i0pdn2=i1pdn2(i1)
       DO i2=1,i0rdn2
       DO j2=1,i0pdn2
          IF(    (i0mlva.EQ.0).AND.(i2.EQ.1))THEN
             i0ctype = 5
          ELSEIF((i0mlvb.NE.i0mlvc).AND.&
                 (i0mlvc.NE.0).AND.&
                 (i2.EQ.i0rdn2))THEN
             i0ctype = 7
          ELSE
             i0ctype = 9
          ENDIF
          WRITE(i0fnum,190)i0ctype
       ENDDO
       ENDDO
    ENDDO
    
    !C WRITE SCALAR DATAS
    
    CALL T2WRIT_CONVERT(d1xxx)
    !C WRITE SCALAR
    WRITE(i0fnum,180)'POINT_DATA', i0nmax3

    !C gm1
    c40lname = 'SCALARS SQRT_G double'
    WRITE(i0fnum,110)c40lname
    WRITE(i0fnum,200)'LOOKUP_TABLE default'
    WRITE(i0fnum,210)(d1jm1(i1), i1 = 1, i0nmax3)
    !C gm2
    c40lname = 'SCALARS Grr double'
    WRITE(i0fnum,110)c40lname
    WRITE(i0fnum,200)'LOOKUP_TABLE default'
    WRITE(i0fnum,210)(d1jm2(i1), i1 = 1, i0nmax3)
    !C gm3
    c40lname = 'SCALARS Grp double'
    WRITE(i0fnum,110)c40lname
    WRITE(i0fnum,200)'LOOKUP_TABLE default'
    WRITE(i0fnum,210)(d1jm3(i1), i1 = 1, i0nmax3)
    !C gm4
    c40lname = 'SCALARS Gpp double'
    WRITE(i0fnum,110)c40lname
    WRITE(i0fnum,200)'LOOKUP_TABLE default'
    WRITE(i0fnum,210)(d1jm4(i1), i1 = 1, i0nmax3)
    !C gm5
    c40lname = 'SCALARS Gtt double'
    WRITE(i0fnum,110)c40lname
    WRITE(i0fnum,200)'LOOKUP_TABLE default'
    WRITE(i0fnum,210)(d1jm5(i1), i1 = 1, i0nmax3)
    !C qprof
!    c40lname = 'SCALARS QF double'
!    WRITE(i0fnum,110)c40lname
!    WRITE(i0fnum,200)'LOOKUP_TABLE default'
!    DO i1 = 1,i0nmax3
!       IF(d1bp3(i1).eq.0.D0)THEN
!          d0qf = 1.D0
!       ELSE
!          d0qf = d1jm5(i1)*d1bt3(i1)/d1bp3(i1)
!       ENDIF
!       WRITE(i0fnum,210)d0qf
!    ENDDO
    !C Bp
    c40lname = 'SCALARS Bp_[T] double'
    WRITE(i0fnum,110)c40lname
    WRITE(i0fnum,200)'LOOKUP_TABLE default'
    WRITE(i0fnum,210)(d1bp3(i1), i1 = 1, i0nmax3)
    !C Bt
    c40lname = 'SCALARS Bt_[T] double'
    WRITE(i0fnum,110)c40lname
    WRITE(i0fnum,200)'LOOKUP_TABLE default'
    WRITE(i0fnum,210)(d1bt3(i1), i1 = 1, i0nmax3)
    !C Er
    c40lname = 'SCALARS Er double'
    WRITE(i0fnum,110)c40lname
    WRITE(i0fnum,200)'LOOKUP_TABLE default'
    WRITE(i0fnum,210)(d1er3(i1), i1 = 1, i0nmax3)
    !C Ep
    c40lname = 'SCALARS Ep double'
    WRITE(i0fnum,110)c40lname
    WRITE(i0fnum,200)'LOOKUP_TABLE default'
    WRITE(i0fnum,210)(d1ep3(i1), i1 = 1, i0nmax3)
    !C Et
    c40lname = 'SCALARS Et double'
    WRITE(i0fnum,110)c40lname
    WRITE(i0fnum,200)'LOOKUP_TABLE default'
    WRITE(i0fnum,210)(d1et3(i1), i1 = 1, i0nmax3)

    !C N
    DO i2 = 1, i0spcs
       WRITE(c2vnumb,'(I2)')i2 
       c2vnumb =  ADJUSTL(c2vnumb)
       i0vnumb = LEN_TRIM(c2vnumb)
       SELECT CASE(i0vnumb)
       CASE(1)
          c2vnumb =     '0'//c2vnumb
          c2vnumb =  ADJUSTL(c2vnumb)
       CASE(2)
          c2vnumb =          c2vnumb
          c2vnumb =  ADJUSTL(c2vnumb)
       CASE default
          WRITE(6,*)'ERROR: ILLEAGAL VARIABLE NUMBER'
          STOP
       END SELECT

       c40lname = 'SCALARS N'//'_'//TRIM(c2vnumb)//' '//'double'
       WRITE(i0fnum,110)c40lname
       WRITE(i0fnum,200)'LOOKUP_TABLE default'
       WRITE(i0fnum,210)(d2n3(i1,i2), i1 = 1, i0nmax3)
    ENDDO

    !C Fr
    DO i2 = 1, i0spcs
       WRITE(c2vnumb,'(I2)')i2 
       c2vnumb =  ADJUSTL(c2vnumb)
       i0vnumb = LEN_TRIM(c2vnumb)
       SELECT CASE(i0vnumb)
       CASE(1)
          c2vnumb =     '0'//c2vnumb
          c2vnumb =  ADJUSTL(c2vnumb)
       CASE(2)
          c2vnumb =          c2vnumb
          c2vnumb =  ADJUSTL(c2vnumb)
       CASE default
          WRITE(6,*)'ERROR: ILLEAGAL VARIABLE NUMBER'
          STOP
       END SELECT

       c40lname = 'SCALARS Fr'//'_'//TRIM(c2vnumb)//' '//'double'
       WRITE(i0fnum,110)c40lname
       WRITE(i0fnum,200)'LOOKUP_TABLE default'
       WRITE(i0fnum,210)(d2fr3(i1,i2), i1 = 1, i0nmax3)
    ENDDO

    !C Fb
    DO i2 = 1, i0spcs
       WRITE(c2vnumb,'(I2)')i2 
       c2vnumb =  ADJUSTL(c2vnumb)
       i0vnumb = LEN_TRIM(c2vnumb)
       SELECT CASE(i0vnumb)       
       CASE(1)
          c2vnumb =     '0'//c2vnumb
          c2vnumb =  ADJUSTL(c2vnumb)
       CASE(2)
          c2vnumb =          c2vnumb
          c2vnumb =  ADJUSTL(c2vnumb)
       CASE default
          WRITE(6,*)'ERROR: ILLEAGAL VARIABLE NUMBER'
          STOP
       END SELECT

       c40lname = 'SCALARS Fb'//'_'//TRIM(c2vnumb)//' '//'double'
       WRITE(i0fnum,110)c40lname
       WRITE(i0fnum,200)'LOOKUP_TABLE default'
       WRITE(i0fnum,210)(d2fb3(i1,i2), i1 = 1, i0nmax3)
    ENDDO

    !C Ft
    DO i2 = 1, i0spcs
       WRITE(c2vnumb,'(I2)')i2 
       c2vnumb =  ADJUSTL(c2vnumb)
       i0vnumb = LEN_TRIM(c2vnumb)
       SELECT CASE(i0vnumb)
       CASE(1)
          c2vnumb =     '0'//c2vnumb
          c2vnumb =  ADJUSTL(c2vnumb)
       CASE(2)
          c2vnumb =          c2vnumb
          c2vnumb =  ADJUSTL(c2vnumb)
       CASE default
          WRITE(6,*)'ERROR: ILLEAGAL VARIABLE NUMBER'
          STOP
       END SELECT

       c40lname = 'SCALARS Ft'//'_'//TRIM(c2vnumb)//' '//'double'
       WRITE(i0fnum,110)c40lname
       WRITE(i0fnum,200)'LOOKUP_TABLE default'
       WRITE(i0fnum,210)(d2ft3(i1,i2), i1 = 1, i0nmax3)
    ENDDO

    !C P
    DO i2 = 1, i0spcs
       WRITE(c2vnumb,'(I2)')i2 
       c2vnumb =  ADJUSTL(c2vnumb)
       i0vnumb = LEN_TRIM(c2vnumb)
       SELECT CASE(i0vnumb)
       CASE(1)
          c2vnumb =     '0'//c2vnumb
          c2vnumb =  ADJUSTL(c2vnumb)
       CASE(2)
          c2vnumb =          c2vnumb
          c2vnumb =  ADJUSTL(c2vnumb)
       CASE default
          WRITE(6,*)'ERROR: ILLEAGAL VARIABLE NUMBER'
          STOP
       END SELECT

       c40lname = 'SCALARS P'//'_'//TRIM(c2vnumb)//' '//'double'
       WRITE(i0fnum,110)c40lname
       WRITE(i0fnum,200)'LOOKUP_TABLE default'
       WRITE(i0fnum,210)(d2p3(i1,i2), i1 = 1, i0nmax3)
    ENDDO

    !C Qr
    DO i2 = 1, i0spcs
       WRITE(c2vnumb,'(I2)')i2 
       c2vnumb =  ADJUSTL(c2vnumb)
       i0vnumb = LEN_TRIM(c2vnumb)
       SELECT CASE(i0vnumb)
       CASE(1)
          c2vnumb =     '0'//c2vnumb
          c2vnumb =  ADJUSTL(c2vnumb)
       CASE(2)
          c2vnumb =          c2vnumb
          c2vnumb =  ADJUSTL(c2vnumb)
       CASE default
          WRITE(6,*)'ERROR: ILLEAGAL VARIABLE NUMBER'
          STOP
       END SELECT

       c40lname = 'SCALARS Qr'//'_'//TRIM(c2vnumb)//' '//'double'
       WRITE(i0fnum,110)c40lname
       WRITE(i0fnum,200)'LOOKUP_TABLE default'
       WRITE(i0fnum,210)(d2Qr3(i1,i2), i1 = 1, i0nmax3)
    ENDDO

    !C Qb
    DO i2 = 1, i0spcs
       WRITE(c2vnumb,'(I2)')i2 
       c2vnumb =  ADJUSTL(c2vnumb)
       i0vnumb = LEN_TRIM(c2vnumb)
       SELECT CASE(i0vnumb) 
       CASE(1)
          c2vnumb =     '0'//c2vnumb
          c2vnumb =  ADJUSTL(c2vnumb)
       CASE(2)
          c2vnumb =          c2vnumb
          c2vnumb =  ADJUSTL(c2vnumb)
       CASE default
          WRITE(6,*)'ERROR: ILLEAGAL VARIABLE NUMBER'
          STOP
       END SELECT

       c40lname = 'SCALARS Qb'//'_'//TRIM(c2vnumb)//' '//'double'
       WRITE(i0fnum,110)c40lname
       WRITE(i0fnum,200)'LOOKUP_TABLE default'
       WRITE(i0fnum,210)(d2Qb3(i1,i2), i1 = 1, i0nmax3)
    ENDDO

    !C Qt
    DO i2 = 1, i0spcs
       WRITE(c2vnumb,'(I2)')i2 
       c2vnumb =  ADJUSTL(c2vnumb)
       i0vnumb = LEN_TRIM(c2vnumb)
       SELECT CASE(i0vnumb) 
       CASE(1)
          c2vnumb =     '0'//c2vnumb
          c2vnumb =  ADJUSTL(c2vnumb)
       CASE(2)
          c2vnumb =          c2vnumb
          c2vnumb =  ADJUSTL(c2vnumb)
       CASE default
          WRITE(6,*)'ERROR: ILLEAGAL VARIABLE NUMBER'
          STOP
       END SELECT

       c40lname = 'SCALARS Qt'//'_'//TRIM(c2vnumb)//' '//'double'
       WRITE(i0fnum,110)c40lname
       WRITE(i0fnum,200)'LOOKUP_TABLE default'
       WRITE(i0fnum,210)(d2Qt3(i1,i2), i1 = 1, i0nmax3)
    ENDDO

    CLOSE(i0fnum)

    RETURN

  END SUBROUTINE T2WRIT_MAIN


 
  SUBROUTINE T2WRIT_CONVERT(d1xxx)

    USE T2COMM,ONLY:&
       d0mfcst,d0btcst,d0ercst,d0epcst,d0etcst,&
       d0nncst,d0frcst,d0fbcst,d0ftcst,&
       d0ppcst,d0qrcst,d0qbcst,d0qtcst,&
         i0nmax3,i0vmax,i0xmax,i0spcs,&
         d1bp3,d1bt3,d1er3,d1ep3,d1et3,&
         d2n3, d2fr3,d2fb3,d2ft3,&
         d2p3, d2qr3,d2qb3,d2qt3,&
         !C
         i2crt,i0nmax1,&
         d2jm1,d1jm1,d1jm2,d1jm3,d1jm4,d1jm5

    REAL(   i0rkind),DIMENSION(1:i0xmax),INTENT(IN)::d1xxx
    INTEGER(i0ikind)::i1,i2,i3,i0vid,i0nid1,i0nid3
    REAL(   i0rkind)::d0jm1,d0jm2,d0jm3,d0jm4,d0jm5
    DO i1 = 1,i0nmax1
       i0nid1 = i2crt(i1,1)
       i0nid3 = i2crt(i1,2)

       d0jm1 = d2jm1(1,i0nid1)
       d0jm2 = d2jm1(2,i0nid1)
       d0jm3 = d2jm1(3,i0nid1)
       d0jm4 = d2jm1(4,i0nid1)
       d0jm5 = d2jm1(5,i0nid1)

       d1jm1(i0nid3) = d0jm1
       d1jm2(i0nid3) = d0jm2
       d1jm3(i0nid3) = d0jm3
       d1jm4(i0nid3) = d0jm4
       d1jm5(i0nid3) = d0jm5
       
       i0vid = i0vmax*(i0nid3 - 1)
       
       d1bp3(i0nid3) = d0mfcst*d1xxx(i0vid+1)*SQRT(d0jm2*d0jm5)
       d1bt3(i0nid3) = d0btcst*d1xxx(i0vid+2)*SQRT(d0jm5)
       d1er3(i0nid3) = d0ercst*d1xxx(i0vid+3)*SQRT(d0jm2)
       d1ep3(i0nid3) = d0epcst*d1xxx(i0vid+4)*SQRT(d0jm4)
       d1et3(i0nid3) = d0etcst*d1xxx(i0vid+5)*SQRT(d0jm5)

       DO i2 = 1, i0spcs
          
          i0vid = i0vmax*(i0nid3 - 1) + 5 + 8*(i2 - 1)
          d2n3( i0nid3,i2) = d1xxx(i0vid+1)
          d2fr3(i0nid3,i2) = d1xxx(i0vid+2)
          d2fb3(i0nid3,i2) = d1xxx(i0vid+3)
          d2ft3(i0nid3,i2) = d1xxx(i0vid+4)
          d2p3( i0nid3,i2) = d1xxx(i0vid+5)
          d2qr3(i0nid3,i2) = d1xxx(i0vid+6)
          d2qb3(i0nid3,i2) = d1xxx(i0vid+7)
          d2qt3(i0nid3,i2) = d1xxx(i0vid+8)
       ENDDO
    ENDDO

    !C


    RETURN

  END SUBROUTINE T2WRIT_CONVERT
  
  SUBROUTINE T2_WRIT_GPT(i0dn,i0tag,d1xxx)
    USE T2COMM,ONLY:&
         i0spcs,i0tstp,i0tmax,i0nmax1,i2crt,&
         d2mfc1,i0xmax,i0vmax,d2rzc1,d0rmjr
    INTEGER(i0ikind),                    INTENT(IN)::i0dn
    INTEGER(i0ikind),                    INTENT(IN)::i0tag
    REAL(   i0rkind),DIMENSION(1:i0xmax),INTENT(IN)::d1xxx

    INTEGER(i0ikind)::i1,i2,i0vid

1000 FORMAT(D15.8,100(',',1X,D15.8))
    
    IF(i0tag.EQ.0) OPEN(i0dn)
    
    WRITE(i0dn,*)'!**************************************************'
    WRITE(i0dn,*)'! TIMESTEP=',i0tstp
    WRITE(i0dn,*)'! 0:t 1:R 2:Z 3:x 4:r 5-9:EM 10-:TR                '
    DO i1 = 1,i0nmax1
       i0vid = i0vmax*(i2crt(i1,2) - 1)
       WRITE(i0dn,1000)d2rzc1(i1,1)-d0rmjr,d2rzc1(i1,2),&
            d2mfc1(i1,2),d2mfc1(i1,1),&
            (d1xxx(i0vid+i2),i2=1,i0vmax)
    ENDDO
    WRITE(i0dn,*)
    WRITE(i0dn,*)
    
    IF(i0tag.EQ.i0tmax) CLOSE(i0dn)
    
  END SUBROUTINE T2_WRIT_GPT

  SUBROUTINE T2_WRIT_GP1(i0dn,i0tag,d1xxx)
    USE T2COMM,ONLY:&
         i0spcs,i0tstp,i0tmax,i0nmax4,i1mfc4,&
         d1mfc4,i0xmax,i0vmax,d2rzc1,d0rmjr
    INTEGER(i0ikind),                    INTENT(IN)::i0dn
    INTEGER(i0ikind),                    INTENT(IN)::i0tag
    REAL(   i0rkind),DIMENSION(1:i0xmax),INTENT(IN)::d1xxx

    INTEGER(i0ikind)::i1,i2,i0vid,i0mfc4

1000 FORMAT(D15.8,100(',',1X,D15.8))
    
    IF(i0tag.EQ.0) OPEN(i0dn)
    
    WRITE(i0dn,*)'!**************************************************'
    WRITE(i0dn,*)'! TIMESTEP=',i0tstp
    WRITE(i0dn,*)'! 0:t 1:R 2:Z 3:x 4:r 5-9:EM 10-:TR                '
    print*,'AAA'
    DO i1 = 1,i0nmax4
       i0mfc4 = i1mfc4(i1)
       i0vid  = i0vmax*(i1mfc4(i1) - 1)
       !C
       print*,'i=',i1,i0mfc4,d1mfc4(i1)
       !C
       WRITE(i0dn,1000)d1mfc4(i1),(d1xxx(i0vid+i2),i2=1,3)
    ENDDO
    WRITE(i0dn,*)
    WRITE(i0dn,*)
    
    IF(i0tag.EQ.i0tmax) CLOSE(i0dn)
    
  END SUBROUTINE T2_WRIT_GP1

 

END MODULE T2WRIT
