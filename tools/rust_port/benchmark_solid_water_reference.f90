program benchmark_solid_water_reference
  use h2o_gibbs
  implicit none
  integer,parameter :: iterations=500000
  integer :: i
  integer(kind=8) :: start,finish,rate
  real :: t,p,values(12),checksum
  character(len=20) :: command
  call get_command_argument(1,command)
  call sho2_init(100.0,-10000.0)
  t=250.0
  if (trim(command)=='check') then
    do i=0,99
      p=1e5+real(i)*1e6
      call evaluate(t,p,values)
      write(*,'(a,i0,a,*(es24.16,:,","))') 'check_',i,',',values
    enddo
  else
    checksum=0.0
    call system_clock(start,rate)
    do i=0,iterations-1
      p=1e5+real(mod(i,100))*1e6
      call evaluate(t,p,values)
      checksum=checksum+sum(values)
    enddo
    call system_clock(finish)
    write(*,'(a,es24.16,a,es24.16)') 'source,',real(finish-start)/real(rate)*1e9/iterations,',',checksum
  endif
contains
  subroutine evaluate(t,p,v)
    real,intent(in) :: t,p
    real,intent(out) :: v(12)
    v=[sh2o_gibbs(t,p),sh2o_dgdp(t,p),sh2o_d2gdp2(t,p),sh2o_dgdT(t,p),sh2o_d2gdT2(t,p), &
         sh2o_d2gdtdp(t,p),sh2o_specvol(t,p),sh2o_entropy(t,p),sh2o_enthalpy(t,p), &
         sh2o_energy(t,p),sh2o_helmholtz(t,p),sh2o_heat_capacity(t,p)]
  end subroutine evaluate
end program benchmark_solid_water_reference
