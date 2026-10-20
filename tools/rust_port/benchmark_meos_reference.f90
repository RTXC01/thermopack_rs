program benchmark_meos_reference
  use eoslibinit, only: init_multiparameter
  use thermopack_var, only: get_active_eos,base_eos_param
  use eos_parameters, only: single_eos
  use multiparameter_base, only: meos
  use multiparameter_lj, only: meos_lj,init_LJTS
  implicit none
  integer,parameter :: iterations=100000
  integer :: profile,i
  integer(kind=8) :: start,finish,rate
  class(base_eos_param),pointer :: active
  class(meos),pointer :: fluid
  type(meos_lj),target :: local
  real :: t,v,values(24),checksum
  character(len=20) :: command
  call get_command_argument(1,command)
  do profile=1,6
    select case(profile)
    case(1)
      call init_multiparameter('C3','NIST_MEOS','DEFAULT')
    case(2)
      call init_multiparameter('R134A','NIST_MEOS','DEFAULT')
    case(3)
      call init_multiparameter('N-H2','NIST_MEOS','DEFAULT')
    case(4)
      call init_multiparameter('AR','LJ_MEOS','DEFAULT')
    case(5)
      call init_multiparameter('AR','LJTS_MEOS','DEFAULT')
    case(6)
      call init_LJTS(local)
    end select
    if (profile==6) then
      fluid=>local
    else
      active=>get_active_eos()
      select type(eos=>active)
      class is(single_eos)
        fluid=>eos%nist(1)%meos
      class default
        error stop 'Expected a pure MEOS'
      end select
    endif
    t=1.4*fluid%tc
    if (trim(command)=='check') then
      do i=0,99
        v=1./(fluid%rc*(.05+.015*real(i)))
        call evaluate(fluid,t,v,values)
        write(*,'(a,i0,a,i0,a,*(es24.16,:,","))') 'check_',profile,'_',i,',',values
      enddo
    else
      checksum=0.
      call system_clock(start,rate)
      do i=0,iterations-1
        v=1./(fluid%rc*(.05+.015*real(mod(i,100))))
        call evaluate(fluid,t,v,values)
        checksum=checksum+sum(values)
      enddo
      call system_clock(finish)
      write(*,'(a,i0,a,es24.16,a,es24.16)') 'source_',profile,',', &
        real(finish-start)/real(rate)*1e9/iterations,',',checksum
    endif
  enddo
contains
  subroutine evaluate(f,t,v,q)
    class(meos) :: f
    real,intent(in) :: t,v
    real,intent(out) :: q(24)
    real :: a0(0:2,0:2),ar(0:2,0:2),delta,tau,p,pr,pt,a,at,av,att,atv,avv,r
    delta=1./v/f%rc;tau=f%tc/t;r=f%Rgas_meos
    call f%alpha0Derivs_taudelta(delta,tau,a0)
    call f%alphaResDerivs_taudelta(delta,tau,ar)
    q(1:6)=[a0(0,0),a0(0,1),a0(0,2),a0(1,0),a0(2,0),a0(1,1)]
    q(7:12)=[ar(0,0),ar(0,1),ar(0,2),ar(1,0),ar(2,0),ar(1,1)]
    call f%mp_pressure(1./v,t,p,pr,pt)
    q(13:15)=[p,pt,-pr/(v*v)]
    call f%alphaDerivs_Tv(t,v,a,at,av,att,atv,avv)
    q(16)=-r*t*t*at;q(17)=-r*(a+t*at);q(18)=q(16)+p*v
    q(19)=r*t*a;q(20)=q(19)+p*v
    q(21)=f%cv(t,v);q(22)=f%cp(t,v)
    q(23)=ar(0,0)+ar(1,0)-log(1.+ar(1,0));q(24)=f%speed_of_sound(t,v)
  end subroutine evaluate
end program benchmark_meos_reference
