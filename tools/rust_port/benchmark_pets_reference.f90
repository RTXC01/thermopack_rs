program benchmark_pets_reference
  use eoslibinit, only: init_pets
  use pets, only: PETS_eos,F_PETS_TVn
  implicit none
  integer, parameter :: iterations=1000000
  type(PETS_eos) :: eos
  real :: v,f(10),checksum(10),n(1)
  integer :: i
  integer(kind=8) :: start,finish,rate
  character(len=20) :: command
  call get_command_argument(1,command)
  call init_pets()
  eos%SIGMA_PETS=3.42e-10
  eos%EPSDIVK_PETS=136.0
  n=[1.0]
  if (trim(command)=='check') then
    do i=0,999
      v=5e-5+real(i)*1e-6
      call evaluate()
      write(*,'(a,i0,a,*(es24.16,:,","))') 'check_',i,',',f
    enddo
    stop
  endif
  checksum=0.0
  call system_clock(start,rate)
  do i=0,iterations-1
    v=5e-5+real(mod(i,1000))*1e-6
    call evaluate()
    checksum=checksum+f
  enddo
  call system_clock(finish)
  write(*,'(a,es24.16,a,*(es24.16,:,","))') 'helmholtz,', &
       real(finish-start)/real(rate)*1e9/iterations,',',checksum
contains
  subroutine evaluate()
    call F_PETS_TVn(eos,300.0,v,n,F=f(1),F_T=f(2),F_V=f(3),F_n=f(4:4), &
         F_TT=f(5),F_TV=f(6),F_Tn=f(7:7),F_VV=f(8),F_Vn=f(9:9),F_nn=f(10:10))
  end subroutine evaluate
end program benchmark_pets_reference
