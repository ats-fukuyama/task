  
    !C PARALLEL FRICTION COEFFICIENTS 
    !C          WITH RESPECT TO MOMENTUM AND TOTAL HEAT FLUX
    !C
    !C DEFINITION OF VARIABLRS 
    !C            FOR NEOCLASSICAL FRICTION COEFFICIENTS
    !C
    !C d2nfcf1: f^{ab}_{1}
    !C d2nfcf2: f^{ab}_{2}
    !C d2nfcf3: f^{ab}_{3}
    !C d2nfcf4: f^{ab}_{4}
    
    DO i0sidj = 1, i0smax
    DO i0sidi = 1, i0smax
       
       d0nfcl11_ab = d2nfcl11(i0sidi,i0sidj)
       d0nfcl12_ab = d2nfcl12(i0sidi,i0sidj)
       d0nfcl21_ab = d2nfcl21(i0sidi,i0sidj)
       d0nfcl22_ab = d2nfcl22(i0sidi,i0sidj)
       
       d0ni_b = d1ni(i0sidj)
       d0pi_b = d1pi(i0sidj)
       
       d0wv4  = d1tt(i0sidi)/d1mm(i0sidi)
       
       d0nfcf1_ab =   (d0nfcl11_ab+d0nfcl12_ab)*d0ni_b 
       d0nfcf2_ab = - 0.4D0*d0nfcl12_ab*d0pi_b 
       d0nfcf3_ab = - (d0nfcl21_ab+d0nfcl22_ab)*d0ni_b
       d0nfcf4_ab =   0.4D0*d0nfcl22_ab*d0pi_b
       
       d0nfcf3_ab = d0wv4*(2.5D0*d0nfcf1_ab+d0nfcf3_ab)
       d0nfcf4_ab = d0wv4*(2.5D0*d0nfcf2_ab+d0nfcf4_ab)
       
       d2nfcf1(i0sidi,i0sidj) = d0nfcf1_ab
       d2nfcf2(i0sidi,i0sidj) = d0nfcf2_ab
       d2nfcf3(i0sidi,i0sidj) = d0nfcf3_ab
       d2nfcf4(i0sidi,i0sidj) = d0nfcf4_ab
       
    ENDDO
    ENDDO
    
    !C NEOCLASSICAL VISCOSITY COEFFICIENTS 
    !C          WITH RESPECT TO MOMENTUM AND TOTAL HEAT FLUX
    !C
    !C DEFINITION OF VARIABLRS 
    !C            FOR NEOCLASSICAL VISCOSITY COEFFICIENTS
    !C
    !C d2nfcf1: f^{ab}_{1}
    !C d2nfcf2: f^{ab}_{2}
    !C d2nfcf3: f^{ab}_{3}
    !C d2nfcf4: f^{ab}_{4}
    
    DO i0sidi = 1, i0smax
       
       d0nvcm1_a = d1nvcm1(i0sidi)
       d0nvcm2_a = d1nvcm2(i0sidi)
       d0nvcm3_a = d1nvcm3(i0sidi)
       
       d0ni_a = d1ni(i0sidi)
       d0pi_a = d1pi(i0sidi)
       d0tt_a = d1tt(i0sidi)
       d0mm_a = d1mm(i0sidi)
       
       d0nvcc1_a =      (d0nvcm1_a-d0nvcm2_a)*d0ni_a
       d0nvcc2_a = 0.4D0*d0nvcm2_a*d0pi_a
       d0nvcc3_a =      (d0nvcm2_a-d0nvcm3_a)*d0ni_a
       d0nvcc4_a = 0.4D0*d0nvcm3_a*d0pi_a
       
       d0nvcc3_a = (d0tt_a/d0mm_a)*(2.5D0*d0nvcc1_a + d0nvcc3_a)
       d0nvcc4_a = (d0tt_a/d0mm_a)*(2.5D0*d0nvcc2_a + d0nvcc4_a)
       
       d1nvcc1(i0sidi) = d0nvcc1_a
       d1nvcc2(i0sidi) = d0nvcc2_a
       d1nvcc3(i0sidi) = d0nvcc3_a
       d1nvcc4(i0sidi) = d0nvcc4_a
       
    ENDDO
    
