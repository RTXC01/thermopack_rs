program benchmark_spinodal_extrema_reference
  use eoslibinit,only:init_cubic
  use spinodal,only:locate_spinodal_prop_min_max_pure_fluid
  implicit none
  integer,parameter :: iterations=100
  integer :: profile,i,prop,ierr
  integer(kind=8) :: start,finish,rate
  integer :: properties(3)=[1,3,6]
  real :: q(6),sums(6),z(1)=[1.]
  character(len=16) :: command
  call get_command_argument(1,command)
  do profile=1,2
    if(profile==1) then
      call init_cubic('C1','PR',mixing='vdW',alpha='Classic',vol_shift=.false.)
    else
      call init_cubic('N2','SRK',mixing='vdW',alpha='Classic',vol_shift=.false.)
    endif
    if(trim(command)=='check') then
      do i=1,3
        call solve(i)
        write(*,'(a,i0,a,i0,a,*(es24.16,:,","))') 'check_',profile,'_',i-1,',',q
      enddo
    else
      sums=0.
      call system_clock(start,rate)
      do i=1,iterations
        call solve(i)
        sums=sums+q
      enddo
      call system_clock(finish)
      write(*,'(a,i0,a,*(es24.16,:,","))') 'source_',profile,',', &
        real(finish-start)/real(rate)*1e9/iterations,sums
    endif
  enddo
contains
  subroutine solve(index)
    integer,intent(in) :: index
    prop=properties(mod(index-1,3)+1)
    call locate_spinodal_prop_min_max_pure_fluid(prop,z,1e6,q(1),q(2),q(3),q(4),q(5),q(6),ierr)
    if(ierr/=0) error stop 'Source extrema search failed'
  end subroutine solve
end program benchmark_spinodal_extrema_reference
