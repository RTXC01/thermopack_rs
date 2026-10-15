! Native benchmark: same composition, temperatures, volumes and derivatives
! as rust/benches/cubic.rs. No Python, ctypes or I/O inside timed loops.
program benchmark_reference
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use eoslibinit, only: init_cubic
  use eostv, only: pressure
  use eos, only: specificvolume,twoPhaseEnthalpy,twoPhaseEntropy
  use tp_solver, only: twoPhaseTPflash
  use ph_solver, only: twoPhasePHflash
  use ps_solver, only: twoPhasePSflash
  use thermopack_constants, only: VAPPH
  use thermopack_var, only: get_active_comps,Rgas
  use compdata, only: gendata_pointer
  use ideal, only: Hideal,Sideal_T,CPideal
  implicit none
  integer, parameter :: iterations=1000000
  integer :: i,phase,ierr,mode,j
  integer(kind=8) :: start, finish, rate
  real :: n(3), v, p, pt, pv, checksum(3)
  real :: t,beta,betaL,x(3),y(3),residual,scale
  real :: enthalpies(0:99),entropies(0:99)
  type(gendata_pointer), pointer :: comps(:)
  call init_cubic('C3','PR',mixing='vdW',alpha='Classic',vol_shift=.false.)
  comps=>get_active_comps()
  checksum=0.0
  call system_clock(start,rate)
  do i=0,iterations-1
    t=150.0+real(mod(i,1000))*0.7
    checksum(1)=checksum(1)+Hideal(comps(1)%p_comp,1,t)
    checksum(2)=checksum(2)+Sideal_T(comps(1)%p_comp,1,t)-Rgas*log(1e5)
    checksum(3)=checksum(3)+CPideal(comps(1)%p_comp,1,t)
  enddo
  call system_clock(finish)
  write(*,'(a,4(es24.16,a))') 'ideal,',real(finish-start)/real(rate)*1e9/iterations,',', &
    checksum(1),',',checksum(2),',',checksum(3),''
  n = [0.7, 0.2, 0.1]
  call init_cubic('CO2,C1,N2','PR',mixing='vdW',alpha='Classic',vol_shift=.false.)
  checksum = 0.0
  call system_clock(start,rate)
  do i=0,iterations-1
    v = 0.0001 + real(mod(i,1000))*1e-6
    p = pressure(300.0,v,n,dpdv=pv,dpdt=pt)
    checksum = checksum + [p,pt,pv]
  enddo
  call system_clock(finish)
  write(*,'(a,4(es24.16,a))') 'pressure,',real(finish-start)/real(rate)*1e9/iterations,',', &
    checksum(1),',',checksum(2),',',checksum(3),''
  checksum = 0.0
  call system_clock(start,rate)
  do i=0,iterations-1
    p = 1e6 + real(mod(i,1000))*1000.0
    call specificvolume(300.0,p,n,VAPPH,v)
    checksum(1) = checksum(1)+v
  enddo
  call system_clock(finish)
  write(*,'(a,2(es24.16,a))') 'volume,',real(finish-start)/real(rate)*1e9/iterations,',',checksum(1),''
  n=[0.5,0.3,0.2]
  call init_cubic('C1,C2,C3','PR',mixing='vdW',alpha='Classic',vol_shift=.false.)
  checksum=0.0
  call system_clock(start,rate)
  do i=0,9999
    t=180.0+real(mod(i,100))*0.8
    call twoPhaseTPflash(t,1e6,n,beta,betaL,phase,x,y)
    checksum=checksum+[beta,betaL,real(phase)]
  enddo
  call system_clock(finish)
  write(*,'(a,4(es24.16,a))') 'tp,',real(finish-start)/real(rate)*1e9/10000,',', &
    checksum(1),',',checksum(2),',',checksum(3),''
  do i=0,99
    t=180.0+real(i)*0.8
    call twoPhaseTPflash(t,1e6,n,beta,betaL,phase,x,y)
    enthalpies(i)=twoPhaseEnthalpy(t,1e6,n,x,y,beta,phase,betaL=betaL)
    entropies(i)=twoPhaseEntropy(t,1e6,n,x,y,beta,phase,betaL=betaL)
  enddo
  do mode=1,2
    ! Preflight each state before timing; source success flags alone do not
    ! establish the recovered temperature or caloric specification.
    do i=0,99
      t=0.0
      beta=0.0
      betaL=0.0
      x=0.0
      y=0.0
      phase=0
      if (mode==1) then
        call twoPhasePHflash(t,1e6,n,beta,betaL,x,y,enthalpies(i),phase,ierr)
      else
        call twoPhasePSflash(t,1e6,n,beta,betaL,x,y,entropies(i),phase,ierr)
      endif
      if (mode==1) then
        residual=twoPhaseEnthalpy(t,1e6,n,x,y,beta,phase,betaL=betaL)-enthalpies(i)
        scale=8.31446261815324*t
      else
        residual=twoPhaseEntropy(t,1e6,n,x,y,beta,phase,betaL=betaL)-entropies(i)
        scale=8.31446261815324
      endif
      ! PS wrapper permits -1 (tolerance flag); independently require closure.
      if (.not.ieee_is_finite(t) .or. .not.ieee_is_finite(residual) &
          .or. ierr>0 .or. ierr < -1 .or. (mode==1 .and. ierr/=0) &
          .or. abs(t-(180.0+real(i)*0.8))>1e-5 .or. abs(residual)>1e-8*scale) then
        write(*,*) 'caloric preflight mode, state, status, temperature: ',mode,i,ierr,t
        error stop 'caloric preflight failed'
      endif
    enddo
    checksum=0.0
    call system_clock(start,rate)
    do i=0,999
      j=mod(i,100)
      t=0.0
      beta=0.0
      betaL=0.0
      x=0.0
      y=0.0
      phase=0
      if (mode==1) then
        call twoPhasePHflash(t,1e6,n,beta,betaL,x,y,enthalpies(j),phase,ierr)
      else
        call twoPhasePSflash(t,1e6,n,beta,betaL,x,y,entropies(j),phase,ierr)
      endif
      if (ierr>0 .or. ierr < -1 .or. (mode==1 .and. ierr/=0)) error stop 'caloric benchmark failed'
      checksum=checksum+[t,beta,real(phase)]
    enddo
    call system_clock(finish)
    if (mode==1) then
      write(*,'(a,4(es24.16,a))') 'ph,',real(finish-start)/real(rate)*1e9/1000,',', &
        checksum(1),',',checksum(2),',',checksum(3),''
    else
      write(*,'(a,4(es24.16,a))') 'ps,',real(finish-start)/real(rate)*1e9/1000,',', &
        checksum(1),',',checksum(2),',',checksum(3),''
    endif
  enddo
end program benchmark_reference
