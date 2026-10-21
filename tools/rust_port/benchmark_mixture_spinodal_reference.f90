program benchmark_mixture_spinodal_reference
  use eoslibinit,only:init_cubic
  use eos,only:getCriticalParam
  use spinodal,only:map_stability_limit
  implicit none
  integer,parameter :: iterations=10
  integer :: profile,i,j,nc,n,ierr
  integer(kind=8) :: start,finish,rate
  real :: t(1000),v(1000),p(1000),z(3),tmin,tc,pc,omega,sums(4)
  character(len=16) :: command
  call get_command_argument(1,command)
  do profile=1,3
    select case(profile)
    case(1)
      call init_cubic('C1,C2','PR',mixing='vdW',alpha='Classic',vol_shift=.false.)
      nc=2;z(:2)=[.5,.5]
    case(2)
      call init_cubic('C1,C2','SRK',mixing='vdW',alpha='Classic',vol_shift=.false.)
      nc=2;z(:2)=[.5,.5]
    case(3)
      call init_cubic('C1,C2,C3','PR',mixing='vdW',alpha='Classic',vol_shift=.false.)
      nc=3;z=[.2,.3,.5]
    end select
    tmin=0.
    do j=1,nc
      call getCriticalParam(j,tc,pc,omega)
      tmin=tmin+.6*z(j)*tc
    enddo
    if (trim(command)=='check') then
      call solve()
      do i=1,n
        write(*,'(a,i0,a,i0,a,*(es24.16,:,","))') 'check_',profile,'_',i-1,',',t(i),v(i),p(i)
      enddo
    else
      sums=0.
      call system_clock(start,rate)
      do i=1,iterations
        call solve()
        sums=sums+[real(n),sum(t(:n)),sum(v(:n)),sum(p(:n))]
      enddo
      call system_clock(finish)
      write(*,'(a,i0,a,*(es24.16,:,","))') 'source_',profile,',', &
        real(finish-start)/real(rate)*1e9/iterations,sums
    endif
  enddo
contains
  subroutine solve()
    call map_stability_limit(1e5,z(:nc),tmin,t,p,v,n,ierr,dlnv_override=.04)
    if (ierr/=0) error stop 'Source spinodal curve failed'
  end subroutine solve
end program benchmark_mixture_spinodal_reference
