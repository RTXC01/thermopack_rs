program benchmark_metastable_reference
  use eoslibinit,only:init_cubic
  use eos,only:specificvolume
  use eostv,only:pressure,entropy_tv
  use spinodal,only:tv_meta_ps,map_meta_isentrope
  use thermopack_constants,only:LIQPH
  implicit none
  integer :: profile,i,j,k,nc,n,queries,iterations,ierr
  integer(kind=8) :: start,finish,rate
  real :: z(2),ts(25),vs(25),ps(25),ss(25),t(12),v(12),p(12),s(12),sums(4)
  character(len=16) :: command
  call get_command_argument(1,command)
  do profile=1,4
    select case(profile)
    case(1,4)
      call init_cubic('C1,C2','PR',mixing='vdW',alpha='Classic',vol_shift=.false.)
      nc=2;z=[.5,.5]
    case(2)
      call init_cubic('CO2,C3','SRK',mixing='vdW',alpha='Classic',vol_shift=.false.)
      nc=2;z=[.5,.5]
    case(3)
      call init_cubic('C1','PR',mixing='vdW',alpha='Classic',vol_shift=.false.)
      nc=1;z(1)=1.
    end select
    queries=5
    if(profile<=2) queries=25
    do i=1,queries
      select case(profile)
      case(1)
        ts(i)=250.+real(i-1)
      case(2)
        ts(i)=350.+real(i-1)
      case(3)
        ts(i)=180.+.5*real(i-1)
      case(4)
        ts(i)=250.+.5*real(i-1)
      end select
      if(profile<=2) then
        vs(i)=2e-3+2e-5*real(i-1)
        ps(i)=pressure(ts(i),vs(i),z(:nc))
      else
        ps(i)=1e7
        if(profile==3) ps(i)=6e6
        call specificvolume(ts(i),ps(i),z(:nc),LIQPH,vs(i))
      endif
      call entropy_tv(ts(i),vs(i),z(:nc),ss(i))
    enddo
    if(trim(command)=='check') then
      do i=1,queries
        call solve(i)
        do j=1,n
          write(*,'(a,i0,a,i0,a,i0,a,*(es24.16,:,","))') 'check_',profile,'_',i-1,'_',j-1,',',t(j),v(j),p(j),s(j)
        enddo
      enddo
    else
      iterations=20
      if(profile<=2) iterations=10000
      sums=0.
      call system_clock(start,rate)
      do k=1,iterations
        call solve(mod(k-1,queries)+1)
        sums=sums+[sum(t(:n)),sum(v(:n)),sum(p(:n)),sum(s(:n))]
      enddo
      call system_clock(finish)
      write(*,'(a,i0,a,*(es24.16,:,","))') 'source_',profile,',', &
        real(finish-start)/real(rate)*1e9/iterations,sums
    endif
  enddo
contains
  subroutine solve(index)
    integer,intent(in) :: index
    if(profile<=2) then
      n=1;t(1)=1.02*ts(index);v(1)=.98*vs(index)
      call tv_meta_ps(ps(index),ss(index),z(:nc),t(1),v(1),ierr)
      if(ierr/=0) error stop 'Source local PS failed'
      p(1)=pressure(t(1),v(1),z(:nc))
    else
      n=12
      call map_meta_isentrope(ps(index),ss(index),z(:nc),1e5,n,t,v,p,ierr)
      if(ierr/=0) error stop 'Source metastable curve failed'
    endif
    do j=1,n
      call entropy_tv(t(j),v(j),z(:nc),s(j))
    enddo
  end subroutine solve
end program benchmark_metastable_reference
