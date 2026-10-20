! Measure the transcendental cost of the reference math runtime.
program bench_transcendental
  implicit none
  integer, parameter :: iterations=1000000
  integer :: i
  real(kind=8) :: t, checksum, start, finish
  checksum = 0.0d0
  call cpu_time(start)
  do i=0,iterations-1
    t = 150.0d0 + dble(mod(i,1000))*0.7d0
    checksum = checksum + exp(-2.0d0*t/300.0d0)
  enddo
  call cpu_time(finish)
  write(*,'(a,1x,f8.2)') 'exp', (finish-start)/dble(iterations)*1.0d9
  call cpu_time(start)
  do i=0,iterations-1
    t = 150.0d0 + dble(mod(i,1000))*0.7d0
    checksum = checksum + log(t)
  enddo
  call cpu_time(finish)
  write(*,'(a,1x,f8.2)') 'log', (finish-start)/dble(iterations)*1.0d9
  call cpu_time(start)
  do i=0,iterations-1
    t = 150.0d0 + dble(mod(i,1000))*0.7d0
    checksum = checksum + sinh(t/300.0d0)
  enddo
  call cpu_time(finish)
  write(*,'(a,1x,f8.2)') 'sinh', (finish-start)/dble(iterations)*1.0d9
  call cpu_time(start)
  do i=0,iterations-1
    t = 150.0d0 + dble(mod(i,1000))*0.7d0
    checksum = checksum + log(1.0d0-exp(-2.0d0*t/300.0d0))
  enddo
  call cpu_time(finish)
  write(*,'(a,1x,f8.2)') 'log1p', (finish-start)/dble(iterations)*1.0d9
  write(*,'(a,1x,es24.16)') 'checksum', checksum
end program bench_transcendental
