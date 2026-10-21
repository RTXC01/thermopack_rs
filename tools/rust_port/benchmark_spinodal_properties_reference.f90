program benchmark_spinodal_properties_reference
  use eoslibinit,only:init_cubic
  use eos,only:getCriticalParam
  use spinodal,only:locate_spinodal_prop_pure_fluid
  implicit none
  integer,parameter :: iterations=100
  integer :: profile,i,ierr
  integer(kind=8) :: start,finish,rate
  real :: t,v,p,target,tc,pc,omega,sums(3)
  logical :: found
  character(len=16) :: command
  call get_command_argument(1,command)
  do profile=1,2
    if (profile==1) then
      call init_cubic('C1','PR',mixing='vdW',alpha='Classic',vol_shift=.false.)
    else
      call init_cubic('N2','SRK',mixing='vdW',alpha='Classic',vol_shift=.false.)
    endif
    call getCriticalParam(1,tc,pc,omega)
    if (trim(command)=='check') then
      do i=0,24
        call solve(i)
        write(*,'(a,i0,a,i0,a,*(es24.16,:,","))') 'check_',profile,'_',i,',',t,v,p
      enddo
    else
      sums=0.
      call system_clock(start,rate)
      do i=0,iterations-1
        call solve(mod(i,25))
        sums=sums+[t,v,p]
      enddo
      call system_clock(finish)
      write(*,'(a,i0,a,*(es24.16,:,","))') 'source_',profile,',', &
        real(finish-start)/real(rate)*1e9/iterations,sums
    endif
  enddo
contains
  subroutine solve(i)
    integer,intent(in)::i
    target=tc*(.94+.04*real(i)/24.)
    call locate_spinodal_prop_pure_fluid(4,[1.0],target,1,1e6,t,v,p,found,ierr)
    if (ierr/=0 .or. .not.found) error stop 'Source intersection failed'
  end subroutine solve
end program benchmark_spinodal_properties_reference
