module FFTW3
      use,intrinsic::iso_c_binding
      include 'fftw3.f03'
end module FFTW3

subroutine fftwf(N,x,iflag)
      use FFTW3  
      implicit none
      integer,intent(in)::N
      complex*16,intent(inout)::x(0:N-1)
      integer,intent(in)::iflag
      !End of interfacial block!
      complex(C_DOUBLE_COMPLEX),dimension(1:N)::inn
      integer(C_SIZE_T)::sz
      complex(C_DOUBLE_COMPLEX),pointer::inp(:),outp(:)
      type(C_PTR)::plan,p1,p2
      integer::i,j,ierr
      integer,parameter::nth=8

      ierr = fftw_init_threads()

      sz=N
      p1=fftw_alloc_complex(sz)
      p2=fftw_alloc_complex(sz)

      call c_f_pointer(p1,inp,[N])
      call c_f_pointer(p2,outp,[N])

      do i=1,N
       inn(i)=x(i-1)
      enddo
      inp=inn

      call fftw_plan_with_nthreads(nth)
      plan=fftw_plan_dft_1d(N,inp,outp,iflag,FFTW_ESTIMATE)
      call fftw_execute_dft(plan,inp,outp)

      do j=1,N
       if (iflag==-1)then
        x(j-1) = outp(j)/N
       else
        x(j-1)=outp(j)
       endif
      enddo


      call fftw_destroy_plan(plan)
      call fftw_free(p1)
      call fftw_free(p2)
      call fftw_cleanup_threads

      return
end subroutine fftwf
