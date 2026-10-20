program benchmark_fmt_reference
  use fundamental_measure_theory, only: fmt_energy_density
  implicit none
  integer,parameter :: points=100,batches=1000
  integer :: model,i,j,k,batch,offset
  integer(kind=8) :: start,finish,rate
  real :: n(points,6),phi(points),gradient(points,6),hessian(points,6,6),v(43),checksum,s
  character(len=20) :: command
  call get_command_argument(1,command)
  do i=1,points
    s=1.+0.01*real(i-1)
    n(i,:)=[.1*s,.2*s,.3*s,.05+.008*real(i-1),.02*s,-.04*s]
  enddo
  do model=1,3
    if (trim(command)=='check') then
      call fmt_energy_density(model,points,6,n,phi,gradient,hessian)
      do k=1,points
        v(1)=phi(k);v(2:7)=gradient(k,:);offset=8
        do i=1,6
          do j=1,6
            v(offset)=hessian(k,i,j);offset=offset+1
          enddo
        enddo
        write(*,'(a,i0,a,i0,a,*(es24.16,:,","))') 'check_',model,'_',k-1,',',v
      enddo
    else
      checksum=0.
      call system_clock(start,rate)
      do batch=1,batches
        call fmt_energy_density(model,points,6,n,phi,gradient,hessian)
        checksum=checksum+sum(phi)+sum(gradient)+sum(hessian)
      enddo
      call system_clock(finish)
      write(*,'(a,i0,a,es24.16,a,es24.16)') 'source_',model,',', &
        real(finish-start)/real(rate)*1e9/(points*batches),',',checksum
    endif
  enddo
end program benchmark_fmt_reference
