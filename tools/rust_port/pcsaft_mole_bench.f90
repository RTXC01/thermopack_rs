subroutine pcsaft_mole_bench(t,v,n,warmup,iterations,fn,checksum,seconds,status) bind(C,name='pcsaft_mole_bench')
  use iso_c_binding, only: c_double,c_int,c_int64_t
  use thermopack_var, only: nc,get_active_eos,base_eos_param
  use pc_saft_nonassoc, only: PCSAFT_eos,sPCSAFT_eos,F_PC_SAFT_TVn,F_sPC_SAFT_TVn
  implicit none
  real(c_double), value :: t,v
  integer(c_int), value :: warmup,iterations
  real(c_double), intent(in) :: n(16)
  real(c_double), intent(out) :: fn(16),checksum,seconds
  integer(c_int), intent(out) :: status
  integer(c_int64_t) :: start,finish,rate
  integer :: k,i
  real :: amounts(nc), result(nc)
  class(base_eos_param), pointer :: eos
  amounts=n(1:nc)
  eos=>get_active_eos()
  fn=0; checksum=0; seconds=0; status=0
  select type(eos)
  class is(PCSAFT_eos)
    do k=1,warmup
      call F_PC_SAFT_TVn(eos,t,v,amounts,F_n=result)
    end do
    call system_clock(start,rate)
    do k=1,iterations
      call F_PC_SAFT_TVn(eos,t,v,amounts,F_n=result)
      do i=1,nc
        checksum=checksum+result(i)*i
      end do
    end do
  class is(sPCSAFT_eos)
    do k=1,warmup
      call F_sPC_SAFT_TVn(eos,t,v,amounts,F_n=result)
    end do
    call system_clock(start,rate)
    do k=1,iterations
      call F_sPC_SAFT_TVn(eos,t,v,amounts,F_n=result)
      do i=1,nc
        checksum=checksum+result(i)*i
      end do
    end do
  class default
    status=1
    return
  end select
  call system_clock(finish)
  seconds=real(finish-start,c_double)/real(rate,c_double)
  fn(1:nc)=result
end subroutine
