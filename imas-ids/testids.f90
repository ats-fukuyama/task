! idstest.f90

PROGRAM ids

  USE task_ids_core_profiles

  IMPLICIT NONE
  INTEGER ns,nr

  CALL get_ids_core_profiles

  DO ns=1,nsmax
     WRITE(6,'(I4,2ES12.4)') ns,pa(ns),pz(ns)
  END DO
  DO nr=1,nrmax
     WRITE(6,'(I4,6ES12.4)') nr,rn(nr,1),ru(nr,1),rtpr(nr,1), &
                                rn(nr,2),ru(nr,2),rtpr(nr,2)
  END DO
END PROGRAM ids
     
  
