program benchmark_pets_roots_reference
  use eoslibinit, only: init_pets
  use thermopack_var, only: get_active_eos,base_eos_param
  use pets, only: PETS_eos
  use saft_interface, only: pc_saft_nonassoc_volume_solver
  use eostv, only: pressure
  implicit none
  integer, parameter :: iterations=20000
  class(base_eos_param), pointer :: active
  real :: t,p,v,actual,pt,pv,checksum
  integer :: i,j,phase
  integer(kind=8) :: start,finish,rate
  character(len=20) :: command,name
  call get_command_argument(1,command)
  call init_pets()
  active=>get_active_eos()
  select type(eos=>active)
  class is(PETS_eos)
    do j=1,3
      if (j==1) then
        name='liquid'; t=120.0; phase=1
      else if (j==2) then
        name='vapor'; t=120.0; phase=2
      else
        name='supercritical'; t=300.0; phase=2
      endif
      if (trim(command)=='check') then
        do i=0,99
          p=5e5+real(i)*5e3
          call pc_saft_nonassoc_volume_solver(1,eos,t,p,[1.0],phase,v)
          actual=pressure(t,v,[1.0],dpdt=pt,dpdv=pv)
          write(*,'(a,i0,a,*(es24.16,:,","))') 'check_'//trim(name)//'_',i,',',v,actual,pt,pv
        enddo
      else
        checksum=0.0
        call system_clock(start,rate)
        do i=0,iterations-1
          p=5e5+real(mod(i,100))*5e3
          call pc_saft_nonassoc_volume_solver(1,eos,t,p,[1.0],phase,v)
          checksum=checksum+v
        enddo
        call system_clock(finish)
        write(*,'(a,es24.16,a,es24.16)') trim(name)//',',real(finish-start)/real(rate)*1e9/iterations,',',checksum
      endif
    enddo
  class default
    error stop 'PeTS benchmark requires active PeTS model'
  end select
end program benchmark_pets_roots_reference
