  subroutine t3d_input_prof(var,nsw)
    implicit none
    character(*), intent(in) :: var
    integer, intent(in) :: nsw

    select case(nsw)
    case(1)
       ! standard format
       print*,'standard files'
     call t3d_read_standard_format(var)
    case(2)
       ! ufile format
    case(0) 
       call t3d_read_standard_format_edge(var)      
    end select
  end subroutine t3d_input_prof

!---------------------------------------------------------------------

  subroutine t3d_read_standard_format(var)
    use TRCOMM, only : RN, RT, QP, Q0, NRMAX, RM, RG
    implicit none
    character(*), intent(in) :: var
    real(8), dimension(:), allocatable :: rho, ne, ni, te, ti, iota, temp_h
    integer :: i, imax
    character :: cha1(1)

    select case(var)
    case('nt')
       open(1,file='profile_nt.txt')
       read(1,*) cha1
       read(1,'(2x,i)') imax
       allocate(rho(imax), ne(imax), ni(imax), te(imax), ti(imax))
       read(1,*) cha1
       do i = 1, imax
          read(1,*) rho(i), ne(i), ni(i), te(i), ti(i)
!          print*, rho(i), ne(i), ni(i), te(i), ti(i)
       end do
       close(1)

       allocate(temp_h(nrmax))
    
       call t3d_convert_mesh(temp_h, rm, nrmax, ne, rho, imax)
       rn(1:nrmax,1) = temp_h(1:nrmax)
    
       call t3d_convert_mesh(temp_h, rm, nrmax, ni, rho, imax)
       rn(1:nrmax,2) = temp_h(1:nrmax)
       
       call t3d_convert_mesh(temp_h, rm, nrmax, te, rho, imax)
       rt(1:nrmax,1) = temp_h(1:nrmax)
       
       call t3d_convert_mesh(temp_h, rm, nrmax, ti, rho, imax)
       rt(1:nrmax,2) = temp_h(1:nrmax)

       deallocate(rho, ne, ni, te, ti, temp_h)
    case('te')
       open(1,file='profile_nt.txt')
       read(1,*) cha1
       read(1,'(2x,i)') imax
       allocate(rho(imax), ne(imax), ni(imax), te(imax), ti(imax))
       read(1,*) cha1
       do i = 1, imax
          read(1,*) rho(i), ne(i), ni(i), te(i), ti(i)
       end do
       close(1)

       allocate(temp_h(nrmax))
    
       call t3d_convert_mesh(temp_h, rm, nrmax, te, rho, imax)
       rt(1:nrmax,1) = temp_h(1:nrmax)

       deallocate(rho, ne, ni, te, ti, temp_h)
    case('ti')
       open(1,file='profile_nt.txt')
       read(1,*) cha1
       read(1,'(2x,i)') imax
       allocate(rho(imax), ne(imax), ni(imax), te(imax), ti(imax))
       read(1,*) cha1
       do i = 1, imax
          read(1,*) rho(i), ne(i), ni(i), te(i), ti(i)
       end do
       close(1)

       allocate(temp_h(nrmax))
    
       call t3d_convert_mesh(temp_h, rm, nrmax, ti, rho, imax)
       rt(1:nrmax,2) = temp_h(1:nrmax)
       do i=1,imax
       print*,i,rt(i,2)
       enddo
       deallocate(rho, ne, ni, te, ti, temp_h)
    case('nn')
       open(1,file='profile_nt.txt')
       read(1,*) cha1
       read(1,'(2x,i)') imax
       allocate(rho(imax), ne(imax), ni(imax), te(imax), ti(imax))
       read(1,*) cha1
       do i = 1, imax
          read(1,*) rho(i), ne(i), ni(i), te(i), ti(i)
       end do
       close(1)

       allocate(temp_h(nrmax))
    
       call t3d_convert_mesh(temp_h, rm, nrmax, ne, rho, imax)
       rn(1:nrmax,1) = temp_h(1:nrmax)
       call t3d_convert_mesh(temp_h, rm, nrmax, ni, rho, imax)
       rn(1:nrmax,2) = temp_h(1:nrmax)

       deallocate(rho, ne, ni, te, ti, temp_h)
    case('iota')
       open(1,file='profile_iota.txt')
       read(1,*) cha1
       read(1,'(2x,i)') imax
       allocate(rho(imax), iota(imax))
       read(1,*) cha1
       do i = 1, imax
          read(1,*) rho(i), iota(i)
       end do
       close(1)

       allocate(temp_h(nrmax))
    
       call t3d_convert_mesh(temp_h, rm, nrmax, iota, rho, imax)
       qp(1:nrmax) = 1.0d0/temp_h(1:nrmax)

       deallocate(rho, iota, temp_h)
    end select
  end subroutine t3d_read_standard_format


  subroutine t3d_read_standard_format_edge(var)
    use TRCOMM, only : RN, RT, QP, Q0, NRMAX, RM, RG
    implicit none
    character(*), intent(in) :: var
    real(8), dimension(:), allocatable :: rho, ne, ni, te, ti, iota, temp_h
    integer :: i, imax
    character :: cha1(1)
    integer :: nredge =5

    print*,'edge read from standard file'
    select case(var)
    case('te')
       open(1,file='profile_nt.txt')
       read(1,*) cha1
       read(1,'(2x,i)') imax
       allocate(rho(imax), ne(imax), ni(imax), te(imax), ti(imax))
       read(1,*) cha1
       do i = 1, imax
          read(1,*) rho(i), ne(i), ni(i), te(i), ti(i)
       end do
       close(1)

       allocate(temp_h(nrmax))
    
       call t3d_convert_mesh(temp_h, rm, nrmax, te, rho, imax)
       do i=nrmax-nredge,nrmax
       rt(i,1) = temp_h(i)
       enddo
       deallocate(rho, ne, ni, te, ti, temp_h)
    case('ti')
       open(1,file='profile_nt.txt')
       read(1,*) cha1
       read(1,'(2x,i)') imax
       allocate(rho(imax), ne(imax), ni(imax), te(imax), ti(imax))
       read(1,*) cha1
       do i = 1, imax
          read(1,*) rho(i), ne(i), ni(i), te(i), ti(i)
       end do
       close(1)

       allocate(temp_h(nrmax))
    
       call t3d_convert_mesh(temp_h, rm, nrmax, ti, rho, imax)
       do i=nrmax-nredge,nrmax
       rt(i,2) = temp_h(i)
       enddo
       deallocate(rho, ne, ni, te, ti, temp_h)
   end select
end subroutine t3d_read_standard_format_edge 
