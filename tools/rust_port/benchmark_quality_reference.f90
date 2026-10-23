program benchmark_quality_reference
  use eoslibinit,only:init_cubic
  use eos,only:specificVolume,enthalpy,entropy
  use saturation,only:sat_newton
  use saturation_curve,only:newton_extrapolate
  implicit none
  integer,parameter::iterations=1000,capacity=256
  integer::profile,i,n,j,unit,counts(4),groups(capacity),specs(capacity),indexes(capacity),outputs
  integer(kind=8)::start,finish,rate
  real::ts(capacity),ps(capacity),bs(capacity),ss(capacity),ks(2,capacity),v(12),sums(12)
  character(len=1024)::command,input,mode
  call get_command_argument(1,command);call get_command_argument(2,input)
  call get_command_argument(3,mode);outputs=10
  if(trim(mode)=='caloric')outputs=12
  open(newunit=unit,file=trim(input),status='old',action='read');read(unit,*)n
  if(n>capacity)error stop 'quality input capacity'
  counts=0
  do i=1,n
    read(unit,*)groups(i),ts(i),ps(i),bs(i),specs(i),ss(i),ks(:,i)
    counts(groups(i))=counts(groups(i))+1
  enddo
  close(unit)
  do profile=1,4
    select case(profile)
    case(1)
      call init_cubic('C1,C2','PR',mixing='vdW',alpha='Classic',vol_shift=.false.)
    case(2)
      call init_cubic('C1,C2','SRK',mixing='vdW',alpha='Classic',vol_shift=.false.)
    case(3)
      call init_cubic('C1,C3','PR',mixing='vdW',alpha='Classic',vol_shift=.true.)
    case(4)
      call init_cubic('C1,C3','PR',mixing='UNIFAC',alpha='Classic',vol_shift=.false.)
    end select
    j=0
    do i=1,n
      if(groups(i)/=profile)cycle
      j=j+1;indexes(j)=i
    enddo
    if(trim(command)=='check')then
      do i=1,j
        call solve(indexes(i),v)
        write(*,'(a,i0,a,i0,a,*(es24.16,:,","))')'check_',profile,'_',i-1,',',v(:outputs)
      enddo
    else
      sums=0.;call system_clock(start,rate)
      do i=0,iterations-1
        call solve(indexes(mod(i,j)+1),v);sums=sums+v
      enddo
      call system_clock(finish)
      write(*,'(a,i0,a,*(es24.16,:,","))')'source_',profile,',',real(finish-start)/real(rate)*1e9/iterations,sums(:outputs)
    endif
  enddo
contains
  subroutine solve(i,v)
    integer,intent(in)::i
    real,intent(out)::v(12)
    real::t,p,k(2),z(2),x(2),y(2),vl,vv,u(4),tangent(4),hl,hv,sl,sv
    integer::iteration,status
    z=.5;t=ts(i);p=ps(i);k=ks(:,i)
    iteration=sat_newton(z,k,t,p,bs(i),specs(i),ss(i),ierr=status)
    if(status/=0)error stop 'source quality corrector failed'
    x=z/(1.-bs(i)+bs(i)*k);y=k*x;x=x/sum(x);y=y/sum(y)
    call specificVolume(t,p,x,1,vl);call specificVolume(t,p,y,2,vv)
    u=[log(k),log(t),log(p)]
    call newton_extrapolate(z,u,tangent,bs(i),specs(i),ss(i))
    v=0.;v(:10)=[t,p,bs(i),vl,vv,x,y,p*tangent(4)/(t*tangent(3))]
    if(outputs==12)then
      call enthalpy(t,p,x,1,hl);call enthalpy(t,p,y,2,hv)
      call entropy(t,p,x,1,sl);call entropy(t,p,y,2,sv)
      v(11)=(1.-bs(i))*hl+bs(i)*hv;v(12)=(1.-bs(i))*sl+bs(i)*sv
    endif
  end subroutine solve
end program benchmark_quality_reference
