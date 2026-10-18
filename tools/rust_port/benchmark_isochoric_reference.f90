program benchmark_isochoric_reference
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use eoslibinit, only: init_cubic
  use eos, only: specificVolume,twoPhaseEnthalpy,twoPhaseEntropy
  use tp_solver, only: twoPhaseTPflash
  use uv_solver, only: twoPhaseUVflash
  use sv_solver, only: twoPhaseSVflash
  implicit none
  integer, parameter :: iterations=2000
  integer :: i,j,mode,phase,unit,status
  integer(kind=8) :: start,finish,rate
  real :: z(3),t,p,u,s,v,bv,bl,x(3),y(3),rows(14,0:99),checksum(3),values(11),actual,actualv,scale
  character(len=100) :: reference,command
  character(len=1000) :: path
  call get_command_argument(1,command)
  call get_command_argument(2,path)
  call get_command_argument(3,reference)
  if (len_trim(reference)==0) reference='DEFAULT'
  call init_cubic('C1,C2,C3','PR',mixing='vdW',alpha='Classic',vol_shift=.false., &
                  parameter_reference=trim(reference))
  z=[0.5,0.3,0.2]
  if (trim(command)=='inputs') then
    do i=0,99
      t=180.0+real(i)*0.8
      p=1e6
      call twoPhaseTPflash(t,p,z,bv,bl,phase,x,y)
      call overall_volume(t,p,z,bv,bl,phase,x,y,v)
      u=twoPhaseEnthalpy(t,p,z,x,y,bv,phase,betaL=bl)-p*v
      s=twoPhaseEntropy(t,p,z,x,y,bv,phase,betaL=bl)
      write(*,'(a,*(es24.16,:,","))') 'input,',[t,p,u,s,v,real(phase),bv,bl,x,y]
    enddo
    stop
  endif
  open(newunit=unit,file=trim(path),status='old',action='read')
  do i=0,99
    read(unit,*) rows(:,i)
  enddo
  close(unit)
  do mode=1,2
    do j=0,99
      call solve(mode,j)
      call overall_volume(t,p,z,bv,bl,phase,x,y,actualv)
      if (mode==1) then
        actual=twoPhaseEnthalpy(t,p,z,x,y,bv,phase,betaL=bl)-p*actualv
        scale=8.31446261815324*t
      else
        actual=twoPhaseEntropy(t,p,z,x,y,bv,phase,betaL=bl)
        scale=8.31446261815324
      endif
      if (.not.ieee_is_finite(actual) .or. abs(actual-rows(mode+2,j))>1e-7*scale &
          .or. abs(actualv-rows(5,j))>1e-7*rows(5,j)) error stop 'volume benchmark specification closure'
      if (trim(command)=='check') then
        if (phase==0) then
          values=[t,p,1.0,bv,bl,x,y]
        else
          values=[t,p,0.0,0.0,1.0,z,0.0,0.0,0.0]
        endif
        if (mode==1) then
          write(*,'(a,i0,a,*(es24.16,:,","))') 'check_uv_',j,',',values
        else
          write(*,'(a,i0,a,*(es24.16,:,","))') 'check_sv_',j,',',values
        endif
      endif
    enddo
    if (trim(command)=='check') cycle
    checksum=0.0
    call system_clock(start,rate)
    do i=0,iterations-1
      j=mod(i,100)
      call solve(mode,j)
      checksum(1)=checksum(1)+t
      checksum(2)=checksum(2)+p
      if (phase==0) checksum(3)=checksum(3)+bv
    enddo
    call system_clock(finish)
    if (mode==1) then
      write(*,'(a,4(es24.16,a))') 'uv,',real(finish-start)/real(rate)*1e9/iterations,',', &
        checksum(1),',',checksum(2),',',checksum(3),''
    else
      write(*,'(a,4(es24.16,a))') 'sv,',real(finish-start)/real(rate)*1e9/iterations,',', &
        checksum(1),',',checksum(2),',',checksum(3),''
    endif
  enddo
contains
  subroutine solve(mode,j)
    integer, intent(in) :: mode,j
    if (mode==1) then
      t=0.0;p=0.0;bv=0.0;bl=0.0;x=0.0;y=0.0;phase=0
      call twoPhaseUVflash(t,p,z,bv,bl,x,y,rows(3,j),rows(5,j),phase)
    else
      t=0.9*rows(1,j);p=1.2*rows(2,j)
      phase=nint(rows(6,j));bv=rows(7,j);bl=rows(8,j);x=rows(9:11,j);y=rows(12:14,j)
      call twoPhaseSVflash(t,p,z,bv,bl,x,y,rows(4,j),rows(5,j),phase,status)
      if (status/=0) error stop 'SV benchmark nonzero status'
    endif
  end subroutine solve
  subroutine overall_volume(t,p,z,bv,bl,phase,x,y,v)
    real, intent(in) :: t,p,z(3),bv,bl,x(3),y(3)
    integer, intent(in) :: phase
    real, intent(out) :: v
    real :: vl,vv
    if (phase==0) then
      call specificVolume(t,p,x,1,vl)
      call specificVolume(t,p,y,2,vv)
      v=bl*vl+bv*vv
    else
      call specificVolume(t,p,z,merge(1,phase,phase==4),v)
    endif
  end subroutine overall_volume
end program benchmark_isochoric_reference
