program benchmark_sound_reference
  use eoslibinit, only: init_cubic
  use speed_of_sound, only: speed_of_sound_TV,sound_velocity_2ph
  use tp_solver, only: twoPhaseTPflash
  implicit none
  integer, parameter :: iterations=500000,eq_iterations=10000
  integer :: i,j,phase(0:99)
  integer(kind=8) :: start,finish,rate
  real :: z(3),v,checksum,t,beta(0:99),betal(0:99),x(3,0:99),y(3,0:99)
  character(len=100) :: reference
  call get_command_argument(1,reference)
  if (len_trim(reference)==0) reference='DEFAULT'
  call init_cubic('CO2,C1,N2','PR',mixing='vdW',alpha='Classic',vol_shift=.false., &
                  parameter_reference=trim(reference))
  z=[0.7,0.2,0.1]
  checksum=0.0
  call system_clock(start,rate)
  do i=0,iterations-1
    v=1e-4+real(mod(i,1000))*1e-6
    checksum=checksum+speed_of_sound_TV(300.0,v,z)
  enddo
  call system_clock(finish)
  write(*,'(a,2(es24.16,a))') 'tv,',real(finish-start)/real(rate)*1e9/iterations,',',checksum,''
  call init_cubic('C1,C2,C3','PR',mixing='vdW',alpha='Classic',vol_shift=.false., &
                  parameter_reference=trim(reference))
  z=[0.5,0.3,0.2]
  do i=0,99
    t=180.0+real(i)*0.8
    call twoPhaseTPflash(t,1e6,z,beta(i),betal(i),phase(i),x(:,i),y(:,i))
  enddo
  checksum=0.0
  call system_clock(start,rate)
  do i=0,eq_iterations-1
    j=mod(i,100)
    t=180.0+real(j)*0.8
    checksum=checksum+sound_velocity_2ph(t,1e6,x(:,j),y(:,j),z,beta(j),betal(j),phase(j))
  enddo
  call system_clock(finish)
  write(*,'(a,2(es24.16,a))') 'equilibrium,',real(finish-start)/real(rate)*1e9/eq_iterations,',',checksum,''
end program benchmark_sound_reference
