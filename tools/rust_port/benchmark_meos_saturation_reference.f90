program benchmark_meos_saturation_reference
  use eoslibinit, only: init_multiparameter
  use saturation, only: safe_bubP,safe_bubT
  use eos, only: specificVolume,enthalpy
  use thermopack_constants, only: LIQPH,VAPPH
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer,parameter :: iterations=1000
  integer :: n,unit,i,profile,mode,index,j,count
  integer,allocatable :: profiles(:),modes(:),selected(:)
  real,allocatable :: inputs(:)
  real :: q(5),sums(5)
  integer(kind=8) :: start,finish,rate
  character(len=1024) :: command,path
  character(len=5) :: label
  call get_command_argument(1,command)
  call get_command_argument(2,path)
  open(newunit=unit,file=trim(path),status='old',action='read')
  read(unit,*) n
  allocate(profiles(n),modes(n),inputs(n),selected(n))
  do i=1,n
    read(unit,*) profiles(i),modes(i),inputs(i)
  enddo
  close(unit)
  do profile=1,5
    select case(profile)
    case(1)
      call init_multiparameter('C3','NIST_MEOS','DEFAULT')
    case(2)
      call init_multiparameter('R134A','NIST_MEOS','DEFAULT')
    case(3)
      call init_multiparameter('KR','LJ_MEOS','DEFAULT')
    case(4)
      call init_multiparameter('AR','LJ_MEOS','DEFAULT')
    case(5)
      call init_multiparameter('KR','LJTS_MEOS','DEFAULT')
    end select
    do mode=1,2
      label='SAT_P'
      if (mode==2) label='SAT_T'
      count=0
      do i=1,n
        if (profiles(i)==profile.and.modes(i)==mode) then
          count=count+1
          selected(count)=i
        endif
      enddo
      if (count==0) error stop 'Missing benchmark profile'
      if (trim(command)=='check') then
        do j=1,count
          call evaluate(mode,inputs(selected(j)),q)
          write(*,'(a,i0,a,a,a,i0,a,*(es24.16,:,","))') 'check_',profile,'_',label,'_',j-1,',',q
        enddo
      else
        sums=0.
        call system_clock(start,rate)
        do i=0,iterations-1
          index=selected(mod(i,count)+1)
          call evaluate(mode,inputs(index),q)
          sums=sums+q
        enddo
        call system_clock(finish)
        write(*,'(a,i0,a,a,a,es24.16,a,*(es24.16,:,","))') 'source_',profile,'_',label,',', &
          real(finish-start)/real(rate)*1e9/iterations,',',sums
      endif
    enddo
  enddo
contains
  subroutine evaluate(mode,input,q)
    integer,intent(in) :: mode
    real,intent(in) :: input
    real,intent(out) :: q(5)
    real :: t,p,x(1),y(1),vl,vv,hl,hv
    integer :: ierr
    x=1.;y=0.
    if (mode==1) then
      t=input;p=safe_bubP(t,x,y,ierr)
    else
      p=input;t=safe_bubT(p,x,y,ierr)
    endif
    if (ierr/=0) error stop 'Source saturation failed'
    call specificVolume(t,p,x,LIQPH,vl)
    call specificVolume(t,p,x,VAPPH,vv)
    call enthalpy(t,p,x,LIQPH,hl)
    call enthalpy(t,p,x,VAPPH,hv)
    q=[t,p,vl,vv,(hv-hl)/(t*(vv-vl))]
    if (.not.all(ieee_is_finite(q)).or.vl>=vv) error stop 'Unqualified source coexistence'
  end subroutine evaluate
end program benchmark_meos_saturation_reference
