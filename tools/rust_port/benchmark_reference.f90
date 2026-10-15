! Native benchmark: same composition, temperatures, volumes and derivatives
! as rust/benches/cubic.rs. No Python, ctypes or I/O inside timed loops.
program benchmark_reference
  use eoslibinit, only: init_cubic
  use eostv, only: pressure
  use eos, only: specificvolume
  use tp_solver, only: twoPhaseTPflash
  use thermopack_constants, only: VAPPH
  implicit none
  integer, parameter :: iterations=1000000
  integer :: i,phase
  integer(kind=8) :: start, finish, rate
  real :: n(3), v, p, pt, pv, checksum(3)
  real :: t,beta,betaL,x(3),y(3)
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
end program benchmark_reference
