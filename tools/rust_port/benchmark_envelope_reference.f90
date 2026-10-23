program benchmark_envelope_reference
  use eoslibinit,only:init_cubic
  use saturation_curve,only:envelopePlot
  use eos,only:specificvolume
  use thermopack_constants,only:LIQPH,VAPPH
  implicit none
  integer,parameter :: nmax=1000,iterations=10
  integer :: profile,i,j,repeat,nc,n
  integer(kind=8) :: start,finish,rate
  real :: t(nmax),p(nmax),k(nmax,3),beta(nmax),vl(nmax),vv(nmax),x(nmax,3),y(nmax,3),z(3),sums(12)
  character(len=16) :: command
  call get_command_argument(1,command)
  do profile=1,4
    nc=2;z(:2)=[.5,.5]
    select case(profile)
    case(1)
      call init_cubic('C1,C2','PR',mixing='vdW',alpha='Classic',vol_shift=.false.)
    case(2)
      call init_cubic('C1,C2','SRK',mixing='vdW',alpha='Classic',vol_shift=.false.)
    case(3)
      call init_cubic('C1,C2,C3','PR',mixing='vdW',alpha='Classic',vol_shift=.false.)
      nc=3;z=[.2,.3,.5]
    case(4)
      call init_cubic('C1,C3','PR',mixing='UNIFAC',alpha='Classic',vol_shift=.false.)
    end select
    if(trim(command)=='check') then
      call solve()
      do i=1,n
        write(*,'(a,i0,a,i0,a,*(es24.16,:,","))') 'check_',profile,'_',i-1,',', &
          t(i),p(i),beta(i),vl(i),vv(i),x(i,:nc),y(i,:nc)
      enddo
    else
      sums=0.
      call system_clock(start,rate)
      do repeat=1,iterations
        call solve()
        sums(:6)=sums(:6)+[real(n),sum(t(:n)),sum(p(:n)),sum(vl(:n)),sum(vv(:n)),sum(beta(:n))]
        do j=1,nc
          sums(6+j)=sums(6+j)+sum(x(:n,j))
          sums(6+nc+j)=sums(6+nc+j)+sum(y(:n,j))
        enddo
      enddo
      call system_clock(finish)
      write(*,'(a,i0,a,*(es24.16,:,","))') 'source_',profile,',', &
        real(finish-start)/real(rate)*1e9/iterations,sums(:6+2*nc)
    endif
  enddo
contains
  subroutine solve()
    integer :: point
    call envelopePlot(z(:nc),0.,1e5,1,1.,1.5e7,nmax,t,p,k(:,:nc),beta,n,dS_override=.25,Tme=80.)
    if(n<2 .or. abs(p(n)-1e5)>2.) error stop 'Source envelope incomplete'
    do point=1,n
      x(point,:nc)=z(:nc)/(1.-beta(point)+beta(point)*k(point,:nc))
      y(point,:nc)=k(point,:nc)*x(point,:nc)
      x(point,:nc)=x(point,:nc)/sum(x(point,:nc))
      y(point,:nc)=y(point,:nc)/sum(y(point,:nc))
      call specificvolume(t(point),p(point),x(point,:nc),LIQPH,vl(point))
      call specificvolume(t(point),p(point),y(point,:nc),VAPPH,vv(point))
    enddo
  end subroutine solve
end program benchmark_envelope_reference
