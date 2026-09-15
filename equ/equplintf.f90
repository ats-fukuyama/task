c=======================================================================
c            data interface program of "TOPICS-EQ"
c                                                         06/09/21
c=======================================================================
      module eqplintf_mod
      contains
c=======================================================================
c            get data from pl
c-----------------------------------------------------------------------
      subroutine eq_getpl
c
      use bpsd_device
      use bpsd_metric1D
      use bpsd_plasmaf
      implicit none
      bpsd_device_type :: device
      bpsd_metric1D_type :: metric1D
      bpsd_plasmaf_type :: plasmaf
c-----------------------------------------------------------------------
      call bpsd_get('device',device)
      call bpsd_get('metric1D',metric1D)
      call bpsd_get('plasmaf',plasmaf)
      return
      end subroutine eq_getpl
c=======================================================================
c            set data to pl
c-----------------------------------------------------------------------
      subroutine eq_putpl
c
      use bpsd_metric1D
      implicit none
      bpsd_metric1D_type :: metric1D
c-----------------------------------------------------------------------
      allocate(metricD%data(nv))
      call bpsd_get('metric1D',metric1D)
      metric1D%nrmax=nv
      do nr=1,nv
         metric1D%data(nr)%rho=sqrt(hiv(nr)/hiv(nv))
         metric1D%data(nr)%psit=hiv(nr)
         metric1D%data(nr)%psip=hiv(nr)
         metric1D%data(nr)%psip=hiv(nr)
         

      call bpsd_put('metric1D',metric1D)
      deallocate(metricD%data)
      return
      end subroutine eq_putpl
c=======================================================================
      end module eqplintf_mod
c=======================================================================
