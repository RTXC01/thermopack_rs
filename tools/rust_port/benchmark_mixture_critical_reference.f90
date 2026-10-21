program benchmark_mixture_critical_reference
  use eoslibinit,only:init_cubic
  use critical,only:calcCriticalTV
  implicit none
  integer,parameter :: iterations=200,points=25
  integer :: profile,i
  integer(kind=8) :: start,finish,rate
  real :: values(3),sums(3)
  character(len=16) :: command
  call get_command_argument(1,command)
  do profile=1,3
    select case(profile)
    case(1)
      call init_cubic('C1,C2','PR',mixing='vdW',alpha='Classic',vol_shift=.false.)
    case(2)
      call init_cubic('C1,C2','SRK',mixing='vdW',alpha='Classic',vol_shift=.false.)
    case(3)
      call init_cubic('C1,C2,C3','PR',mixing='vdW',alpha='Classic',vol_shift=.false.)
    end select
    if (trim(command)=='check') then
      do i=0,points-1
        call solve(profile,i,values)
        write(*,'(a,i0,a,i0,a,*(es24.16,:,","))') 'check_',profile,'_',i,',',values
      enddo
    else
      sums=0.
      call system_clock(start,rate)
      do i=0,iterations-1
        call solve(profile,mod(i,points),values)
        sums=sums+values
      enddo
      call system_clock(finish)
      write(*,'(a,i0,a,*(es24.16,:,","))') 'source_',profile,',', &
        real(finish-start)/real(rate)*1e9/iterations,sums
    endif
  enddo
contains
  subroutine solve(profile,i,result)
    integer,intent(in)::profile,i
    real,intent(out)::result(3)
    real :: t,v,p,x
    integer :: ierr
    t=0.;v=0.;x=real(i)/real(points-1)
    if (profile==3) then
      call calcCriticalTV(t,v,[.1+.6*x,.2,.7-.6*x],ierr,p=p)
    else
      call calcCriticalTV(t,v,[.1+.8*x,.9-.8*x],ierr,p=p)
    endif
    if (ierr/=0) error stop 'Source critical solve failed'
    result=[t,v,p]
  end subroutine solve
end program benchmark_mixture_critical_reference
