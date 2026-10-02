    !C
    !C 
    !C ADDITIONAL INTEGRAL ARRAYS FOR SUPG-FEM
    !C
    !C

    !C 
    !C MS: D5IMSS
    !C
    
    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax
    DO i0nidl = 1, i0nmax
    DO i0nidk = 1, i0nmax
       DO i0didi = 1, i0dmax
          d0temp = 0.D0
          DO i0qidj = 1, i0qmax
          DO i0qidi = 1, i0qmax
             d0ifnci = d4ifnc(i0qidi,i0qidj,i0didi,i0nidi) 
             d0ifncj = d4ifnc(i0qidi,i0qidj,0     ,i0nidj) 
             d0ifnck = d4ifnc(i0qidi,i0qidj,0     ,i0nidk) 
             d0ifncl = d4ifnc(i0qidi,i0qidj,0     ,i0nidl) 
             d0wfct  = d2wfct(i0qidi,i0qidj              )
             d0temp  = d0temp + d0ifnci*d0ifncj*d0ifnck*d0ifncl*d0wfct
          ENDDO
          ENDDO
          d5imss(i0didi,i0nidk,i0nidl,i0nidi,i0nidj) = d0temp
       ENDDO
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    
    !C
    !C AV: D6IAVS
    !C
    
    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax
    DO i0nidl = 1, i0nmax
    DO i0nidk = 1, i0nmax
       DO i0didj = 1, i0dmax
       DO i0didi = 1, i0dmax
          d0temp = 0.D0
          DO i0qidj = 1, i0qmax
          DO i0qidi = 1, i0qmax
             d0ifnci = d4ifnc(i0qidi,i0qidj,i0didi,i0nidi)
             d0ifncj = d4ifnc(i0qidi,i0qidj,i0didj,i0nidj)
             d0ifnck = d4ifnc(i0qidi,i0qidj,0     ,i0nidk)
             d0ifncl = d4ifnc(i0qidi,i0qidj,0     ,i0nidl)
             d0wfct  = d2wfct(i0qidi,i0qidj              )
             d0temp  = d0temp + d0ifnci*d0ifncj*d0ifnck*d0ifncl*d0wfct
             d0ifnci = d4ifnc(i0qidi,i0qidj,i0didi,i0nidi)
             d0ifncj = d4ifnc(i0qidi,i0qidj,0     ,i0nidj)
             d0ifnck = d4ifnc(i0qidi,i0qidj,i0didj,i0nidk)
             d0ifncl = d4ifnc(i0qidi,i0qidj,0     ,i0nidl)
             d0wfct  = d2wfct(i0qidi,i0qidj              )
             d0temp  = d0temp + d0ifnci*d0ifncj*d0ifnck*d0ifncl*d0wfct
          ENDDO
          ENDDO
          d6iavs(i0didi,i0didj,i0nidk,i0nidl,i0nidi,i0nidj) = d0temp
       ENDDO
       ENDDO
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    
    !C
    !C AT: D8IATS*
    !C
    !C    if interpolation functions are quadratic, 
    !C    2nd-order derivatives appear
    !C
    
    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax
    DO i0nidm = 1, i0nmax
    DO i0nidl = 1, i0nmax
    DO i0nidk = 1, i0nmax
       DO i0didk = 1, i0dmax
       DO i0didj = 1, i0dmax
       DO i0didi = 1, i0dmax
          d0temp = 0.D0
          DO i0qidj = 1, i0qmax
          DO i0qidi = 1, i0qmax
             d0ifnci = d4ifnc(i0qidi,i0qidj,i0didi,i0nidi)
             d0ifncj = d4ifnc(i0qidi,i0qidj,i0didj,i0nidj)
             d0ifnck = d4ifnc(i0qidi,i0qidj,i0didk,i0nidk)
             d0ifncl = d4ifnc(i0qidi,i0qidj,0     ,i0nidl)
             d0ifncm = d4ifnc(i0qidi,i0qidj,0     ,i0nidm)
             d0wfct  = d2wfct(i0qidi,i0qidj              )
             d0temp  = d0temp&
                  &  + d0ifnci*d0ifncj*d0ifnck*d0ifncl*d0ifncm*d0wfct
             d0ifnci = d4ifnc(i0qidi,i0qidj,i0didi,i0nidi)
             d0ifncj = d4ifnc(i0qidi,i0qidj,0     ,i0nidj)
             d0ifnck = d4ifnc(i0qidi,i0qidj,i0didk,i0nidk)
             d0ifncl = d4ifnc(i0qidi,i0qidj,i0didj,i0nidl)
             d0ifncm = d4ifnc(i0qidi,i0qidj,0     ,i0nidm)
             d0wfct  = d2wfct(i0qidi,i0qidj              )
             d0temp  = d0temp&
                  &  + d0ifnci*d0ifncj*d0ifnck*d0ifncl*d0ifncm*d0wfct
          ENDDO
          ENDDO
          d8iats(i0didi,i0didj,i0didk,&
               & i0nidk,i0nidl,i0nidm,i0nidi,i0nidj) = d0temp
       ENDDO
       ENDDO
       ENDDO
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    ENDDO

    !C
    !C DT: D7IDTS*
    !C
    !C    if interpolation functions are quadratic, 
    !C    2nd-order derivatives appear
    !C
    
    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax
    DO i0nidl = 1, i0nmax
    DO i0nidk = 1, i0nmax
       DO i0didk = 1, i0dmax
       DO i0didj = 1, i0dmax
       DO i0didi = 1, i0dmax
          d0temp = 0.D0
          DO i0qidj = 1, i0qmax
          DO i0qidi = 1, i0qmax
             d0ifnci = d4ifnc(i0qidi,i0qidj,i0didi,i0nidi)
             d0ifncj = d4ifnc(i0qidi,i0qidj,i0didk,i0nidj)
             d0ifnck = d4ifnc(i0qidi,i0qidj,i0didj,i0nidk)
             d0ifncl = d4ifnc(i0qidi,i0qidj,0     ,i0nidl)
             d0wfct  = d2wfct(i0qidi,i0qidj              )
             d0temp  = d0temp + d0ifnci*d0ifncj*d0ifnck*d0ifncl*d0wfct
          ENDDO
          ENDDO
          d7idts(i0didi,i0didj,i0didk,&
               & i0nidk,i0nidl,i0nidi,i0nidj) = d0temp
       ENDDO
       ENDDO
       ENDDO
    ENDDO
    ENDDO
    ENDDO
    ENDDO

    !C
    !C GV: D6IGVS
    !C

    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax
    DO i0nidl = 1, i0nmax
    DO i0nidk = 1, i0nmax
       DO i0didj = 1, i0dmax
       DO i0didi = 1, i0dmax
          d0temp = 0.D0
          DO i0qidj = 1, i0qmax
          DO i0qidi = 1, i0qmax
             d0ifnci = d4ifnc(i0qidi,i0qidj,i0didi,i0nidi)
             d0ifncj = d4ifnc(i0qidi,i0qidj,i0didj,i0nidj)
             d0ifnck = d4ifnc(i0qidi,i0qidj,0     ,i0nidk)
             d0ifncl = d4ifnc(i0qidi,i0qidj,0     ,i0nidl)
             d0wfct  = d2wfct(i0qidi,i0qidj              )
             d0temp  = d0temp + d0ifnci*d0ifncj*d0ifnck*d0ifncl*d0wfct
          ENDDO
          ENDDO
          d6igvs(i0didi,i0didj,i0nidk,i0nidl,i0nidi,i0nidj) = d0temp
       ENDDO
       ENDDO
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    
    !C
    !C GT: D8IGVS
    !C

    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax
    DO i0nidm = 1, i0nmax
    DO i0nidl = 1, i0nmax
    DO i0nidk = 1, i0nmax
       DO i0didk = 1, i0dmax
       DO i0didj = 1, i0dmax
       DO i0didi = 1, i0dmax
          d0temp = 0.D0
          DO i0qidj = 1, i0qmax
          DO i0qidi = 1, i0qmax
             d0ifnci = d4ifnc(i0qidi,i0qidj,i0didi,i0nidi)
             d0ifncj = d4ifnc(i0qidi,i0qidj,i0didk,i0nidj)
             d0ifnck = d4ifnc(i0qidi,i0qidj,i0didj,i0nidk)
             d0ifncl = d4ifnc(i0qidi,i0qidj,0     ,i0nidl)
             d0ifncm = d4ifnc(i0qidi,i0qidj,0     ,i0nidm)
             d0wfct  = d2wfct(i0qidi,i0qidj              )
             d0temp  = d0temp&
                  &  + d0ifnci*d0ifncj*d0ifnck*d0ifncl*d0ifncm*d0wfct
          ENDDO
          ENDDO
          d8igts(i0didi,i0didj,i0didk,&
               & i0nidk,i0nidl,i0nidm,i0nidi,i0nidj) = d0temp
       ENDDO
       ENDDO
       ENDDO
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    
    !C
    !C ES: D5IESS
    !C
    
    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax
    DO i0nidl = 1, i0nmax
    DO i0nidk = 1, i0nmax
       DO i0didi = 1, i0dmax
          d0temp = 0.D0
          DO i0qidj = 1, i0qmax
          DO i0qidi = 1, i0qmax
             d0ifnci = d4ifnc(i0qidi,i0qidj,i0didi,i0nidi)
             d0ifncj = d4ifnc(i0qidi,i0qidj,0     ,i0nidj)
             d0ifnck = d4ifnc(i0qidi,i0qidj,0     ,i0nidk)
             d0ifncl = d4ifnc(i0qidi,i0qidj,0     ,i0nidl)
             d0wfct  = d2wfct(i0qidi,i0qidj              )
             d0temp  = d0temp&
                  &  + d0ifnci*d0ifncj*d0ifnck*d0ifncl*d0wfct
          ENDDO
          ENDDO
          d5iess(i0didi,i0nidk,i0nidl,i0nidi,i0nidj) = d0temp
       ENDDO
    ENDDO
    ENDDO
    ENDDO
    ENDDO
   
    !C
    !C EV: D7IEVS
    !C
    
    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax
    DO i0nidm = 1, i0nmax
    DO i0nidl = 1, i0nmax
    DO i0nidk = 1, i0nmax
       DO i0didj = 1, i0dmax
       DO i0didi = 1, i0dmax
          d0temp = 0.D0
          DO i0qidj = 1, i0qmax
          DO i0qidi = 1, i0qmax
             d0ifnci = d4ifnc(i0qidi,i0qidj,i0didi,i0nidi)
             d0ifncj = d4ifnc(i0qidi,i0qidj,0,     i0nidj)
             d0ifnck = d4ifnc(i0qidi,i0qidj,i0didj,i0nidk)
             d0ifncl = d4ifnc(i0qidi,i0qidj,0,     i0nidl)
             d0ifncm = d4ifnc(i0qidi,i0qidj,0,     i0nidm)
             d0wfct  = d2wfct(i0qidi,i0qidj              )
             d0temp  = d0temp&
                  &  + d0ifnci*d0ifncj*d0ifnck*d0ifncl*d0ifncm*d0wfct
          ENDDO
          ENDDO
          d7ievs(i0didi,i0didj,&
               & i0nidk,i0nidl,i0nidm,i0nidi,i0nidj) = d0temp
       ENDDO
       ENDDO
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    ENDDO
   
    !C
    !C ET: D9IETS
    !C

    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax
    DO i0nidn = 1, i0nmax
    DO i0nidm = 1, i0nmax
    DO i0nidl = 1, i0nmax
    DO i0nidk = 1, i0nmax
       DO i0didk = 1, i0dmax
       DO i0didj = 1, i0dmax
       DO i0didi = 1, i0dmax
          d0temp = 0.D0
          DO i0qidj = 1, i0qmax
          DO i0qidi = 1, i0qmax
             d0ifnci = d4ifnc(i0qidi,i0qidj,i0didi,i0nidi)
             d0ifncj = d4ifnc(i0qidi,i0qidj,0,     i0nidj)
             d0ifnck = d4ifnc(i0qidi,i0qidj,i0didj,i0nidk)
             d0ifncl = d4ifnc(i0qidi,i0qidj,i0didk,i0nidl)
             d0ifncm = d4ifnc(i0qidi,i0qidj,0,     i0nidm)
             d0ifncn = d4ifnc(i0qidi,i0qidj,0,     i0nidn)
             d0wfct  = d2wfct(i0qidi,i0qidj              )
             d0temp  = d0temp&
                  &  + d0ifnci*d0ifncj*d0ifnck*d0ifncl*d0ifncm*d0ifncn*d0wfct
          ENDDO
          ENDDO
          d9iets(i0didi,i0didj,i0didk,&
               & i0nidk,i0nidl,i0nidm,i0nidn,i0nidi,i0nidj) = d0temp
       ENDDO
       ENDDO
       ENDDO
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    
    !C
    !C SS: D4ISSS
    !C
    
    DO i0nidj = 1, i0nmax
    DO i0nidi = 1, i0nmax
    DO i0nidk = 1, i0nmax
       DO i0didi = 1, i0dmax
          d0temp = 0.D0
          DO i0qidj = 1, i0qmax
          DO i0qidi = 1, i0qmax
             d0ifnci = d4ifnc(i0qidi,i0qidj,i0didi,i0nidi)
             d0ifncj = d4ifnc(i0qidi,i0qidj,0     ,i0nidj)
             d0ifnck = d4ifnc(i0qidi,i0qidj,0     ,i0nidk)
             d0wfct  = d2wfct(i0qidi,i0qidj              )
             d0temp  = d0temp + d0ifnci*d0ifncj*d0ifnck*d0wfct
          ENDDO
          ENDDO
          d4isss(i0didi,i0nidk,i0nidi,i0nidj) = d0temp
       ENDDO
    ENDDO
    ENDDO
    ENDDO
