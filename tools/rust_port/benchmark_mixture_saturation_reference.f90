program benchmark_mixture_saturation_reference
  use eoslibinit,only:init_cubic
  use eos,only:specificVolume
  use saturation,only:safe_bubP,safe_dewP,safe_bubT,safe_dewT
  implicit none
  integer,parameter :: iterations=1000
  integer :: profile,mode,i
  integer(kind=8) :: start,finish,rate
  real :: values(8),sums(8)
  character(len=16) :: command
  character(len=3),parameter :: eos_names(2)=['PR ','SRK']
  character(len=2),parameter :: methods(4)=['BP','DP','BT','DT']
  call get_command_argument(1,command)
  do profile=1,2
    call init_cubic('C1,C2',trim(eos_names(profile)),mixing='vdW',alpha='Classic',vol_shift=.false.)
    do mode=1,4
      if (trim(command)=='check') then
        do i=0,99
          call solve(mode,i,values)
          write(*,'(a,i0,a,a,a,i0,a,*(es24.16,:,","))') 'check_',profile,'_',methods(mode),'_',i,',',values
        enddo
      else
        sums=0.
        call system_clock(start,rate)
        do i=0,iterations-1
          call solve(mode,mod(i,100),values)
          sums=sums+values
        enddo
        call system_clock(finish)
        write(*,'(a,i0,a,a,a,*(es24.16,:,","))') 'source_',profile,'_',methods(mode),',', &
          real(finish-start)/real(rate)*1e9/iterations,sums
      endif
    enddo
  enddo
contains
  subroutine solve(mode,i,v)
    integer,intent(in)::mode,i
    real,intent(out)::v(8)
    real :: t,p,x(2),y(2),z(2),vl,vv
    integer :: ierr
    z=[.5,.5];t=160.+.4*real(i);p=1e5+3e4*real(i)
    select case(mode)
    case(1)
      x=z;p=safe_bubP(t,x,y,ierr)
    case(2)
      y=z;p=safe_dewP(t,x,y,ierr)
    case(3)
      x=z;t=safe_bubT(p,x,y,ierr)
    case(4)
      y=z;t=safe_dewT(p,x,y,ierr)
    end select
    if (ierr/=0) error stop 'Source saturation solve failed'
    call specificVolume(t,p,x,1,vl);call specificVolume(t,p,y,2,vv)
    v=[t,p,vl,vv,x,y]
  end subroutine solve
end program benchmark_mixture_saturation_reference
