subroutine pcsaft_mole_gradient(t,v,n,fn,status) bind(C,name='pcsaft_mole_gradient')
  use iso_c_binding, only: c_double,c_int
  use thermopack_var, only: nc,get_active_eos,base_eos_param
  use pc_saft_nonassoc, only: PCSAFT_eos,sPCSAFT_eos,F_PC_SAFT_TVn,F_sPC_SAFT_TVn
  implicit none
  real(c_double), value :: t,v
  real(c_double), intent(in) :: n(16)
  real(c_double), intent(out) :: fn(16)
  integer(c_int), intent(out) :: status
  class(base_eos_param), pointer :: eos
  real :: amounts(nc), result(nc)
  integer :: i
  amounts=n(1:nc)
  eos=>get_active_eos()
  status=0
  select type(eos)
  class is(PCSAFT_eos)
    call F_PC_SAFT_TVn(eos,t,v,amounts,F_n=result)
  class is(sPCSAFT_eos)
    call F_sPC_SAFT_TVn(eos,t,v,amounts,F_n=result)
  class default
    status=1
    return
  end select
  do i=1,nc
    fn(i)=result(i)
  end do
end subroutine
