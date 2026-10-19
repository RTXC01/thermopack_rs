program benchmark_mbwr_density_reference
  use mbwr, only: eosmbwr,initializeMBWRmodel,makeParam,MBWR_density,MBWR_pressure
  implicit none
  integer,parameter :: iterations=20000
  type(eosmbwr) :: model
  real,allocatable :: param(:)
  real :: t,p,rho,actual,dp,checksum
  integer :: family,j,mode,i,phase,found
  integer(kind=8) :: start,finish,rate
  character(len=30) :: command,workload,preparation,name
  character(len=2) :: identifier
  call get_command_argument(1,command)
  do family=19,32,13
    call initializeMBWRmodel('C1',model,family)
    if (allocated(param)) deallocate(param)
    allocate(param(1+model%bplen+model%belen))
    write(identifier,'(i2)') family
    do j=1,3
      if (j==1) then
        workload='liquid';t=150.0;phase=1
      else if (j==2) then
        workload='vapor';t=150.0;phase=2
      else
        workload='supercritical';t=300.0;phase=2
      endif
      do mode=1,2
        preparation='prepared'
        if (mode==1) preparation='fresh'
        name=identifier//'_'//trim(workload)//'_'//trim(preparation)
        call makeParam(param,t,model,0)
        if (trim(command)=='check') then
          do i=0,99
            p=5e5+real(i)*5e3
            if (mode==1) call makeParam(param,t,model,0)
            rho=MBWR_density(t,p,phase,param,model,found)
            call MBWR_pressure(rho,param,actual,dp)
            write(*,'(a,i0,a,*(es24.16,:,","))') 'check_'//trim(name)//'_',i,',',rho*1000.0,real(found),actual,dp*1e-3
          enddo
        else
          checksum=0.0
          call system_clock(start,rate)
          do i=0,iterations-1
            p=5e5+real(mod(i,100))*5e3
            if (mode==1) call makeParam(param,t,model,0)
            rho=MBWR_density(t,p,phase,param,model,found)
            call MBWR_pressure(rho,param,actual,dp)
            checksum=checksum+rho*1000.0+real(found)+actual*1e-9+dp*1e-8
          enddo
          call system_clock(finish)
          write(*,'(a,es24.16,a,es24.16)') trim(name)//',',real(finish-start)/real(rate)*1e9/iterations,',',checksum
        endif
      enddo
    enddo
  enddo
end program benchmark_mbwr_density_reference
