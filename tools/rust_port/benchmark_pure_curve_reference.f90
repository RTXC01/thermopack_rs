program benchmark_pure_curve_reference
  use eoslibinit,only:init_cubic,init_quantum_cubic
  use saturation_curve,only:pure_fluid_saturation_wrapper,envelopePlot
  use eos,only:getCriticalParam,specificvolume
  use thermopack_constants,only:LIQPH,VAPPH
  use thermopack_var,only:tpTmin
  implicit none
  integer,parameter :: nmax=100,iterations=20
  integer :: profile,n,i,repeat,point
  integer(kind=8) :: start,finish,rate
  real :: z(1),t(nmax),p(nmax),vl(nmax),vv(nmax),tc,pc,w,value,sums(5),ratios(nmax,1),beta(nmax)
  logical :: from_t,log_grid,envelope
  character(len=16) :: command,mode
  call get_command_argument(1,command)
  call get_command_argument(2,mode)
  envelope=trim(mode)=='envelope'
  z=1.
  do profile=1,4
    select case(profile)
    case(1)
      call init_cubic('C1','PR',mixing='vdW',alpha='Classic',vol_shift=.false.)
    case(2)
      call init_cubic('CO2','SRK',mixing='vdW',alpha='Classic',vol_shift=.false.)
    case(3)
      call init_cubic('H2O','PR',mixing='vdW',alpha='Classic',vol_shift=.true.)
    case(4)
      call init_quantum_cubic('H2','vdW')
    end select
    tpTmin=80.
    if(profile==4) tpTmin=2.
    call getCriticalParam(1,tc,pc,w)
    from_t=(profile==2);log_grid=(profile==2 .or. profile==4)
    value=.05*pc
    if(from_t) value=.65*tc
    if(trim(command)=='check') then
      call solve()
      do i=1,n
        write(*,'(a,i0,a,i0,a,*(es24.16,:,","))') 'check_',profile,'_',i-1,',',t(i),p(i),vl(i),vv(i)
      enddo
    else
      sums=0.
      call system_clock(start,rate)
      do repeat=1,iterations
        call solve()
        sums=sums+[real(n),sum(t(:n)),sum(p(:n)),sum(vl(:n)),sum(vv(:n))]
      enddo
      call system_clock(finish)
      write(*,'(a,i0,a,*(es24.16,:,","))') 'source_',profile,',', &
        real(finish-start)/real(rate)*1e9/iterations,sums
    endif
  enddo
contains
  subroutine solve()
    if(envelope) then
      if(from_t) then
        call envelopePlot(z,value,0.,2,1.,1.5e7,nmax,t,p,ratios,beta,n)
      else
        call envelopePlot(z,0.,value,1,1.,1.5e7,nmax,t,p,ratios,beta,n)
      endif
      if(any(beta(:n)/=0.) .or. any(ratios(:n,:)/=1.)) error stop 'Source pure envelope labels'
      do point=1,n
        call specificvolume(t(point),p(point),z,LIQPH,vl(point))
        call specificvolume(t(point),p(point),z,VAPPH,vv(point))
      enddo
    else
      call pure_fluid_saturation_wrapper(z,value,from_t,2e4,log_grid,t,p,vl,vv,nmax,n)
    endif
    if(n<1 .or. abs(t(n)-tc)>1e-6 .or. abs(p(n)-pc)>.1) error stop 'Incomplete source curve'
    if(.not.envelope .and. n/=nmax) error stop 'Incomplete source grid'
  end subroutine solve
end program benchmark_pure_curve_reference
