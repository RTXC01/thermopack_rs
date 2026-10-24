program benchmark_meos_flash_reference
  use eoslibinit, only: init_multiparameter
  use tp_solver, only: twoPhaseTPflash
  use ph_solver, only: twoPhasePHflash
  use ps_solver, only: twoPhasePSflash
  use eos, only: specificVolume,enthalpy,entropy
  use eostv, only: enthalpy_tv,entropy_tv
  use thermopack_constants, only: LIQPH,VAPPH,TWOPH,SINGLEPH
  use thermopack_var, only: tpTmin,tpTmax,activate_model,set_tmin,set_tmax,get_active_thermo_model,thermo_model
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer :: n,unit,i,profile,mode,kind,phase,ierr,j,count,idx,model_index,last_status
  integer, allocatable :: profiles(:),modes(:),kinds(:),has_guess(:),selected(:)
  real, allocatable :: values(:,:)
  type(thermo_model),pointer :: active_model
  real :: t,p,spec,guess,tmin,tmax,beta,betal,z(1),x(1),y(1),vl,vv,h,s,hl,hv,sl,sv,h2,s2,q(8),sums(8),dt,eps,vscale,pscale,hscale,sscale
  integer(kind=8) :: start,finish,rate
  character(len=16) :: command
  character(len=1024) :: path
  call get_command_argument(1,command); call get_command_argument(2,path)
  open(newunit=unit,file=trim(path),status='old',action='read'); read(unit,*) n
  allocate(profiles(n),modes(n),kinds(n),has_guess(n),selected(n),values(14,n))
  do i=1,n
    read(unit,*) profiles(i),modes(i),kinds(i),values(1:4,i),has_guess(i),values(5:14,i)
  enddo
  close(unit)
  do profile=0,5
    if (.not.any(profiles==profile)) cycle
    select case(profile)
    case(0); call init_multiparameter('C3','NIST_MEOS','DEFAULT')
    case(1); call init_multiparameter('R134A','NIST_MEOS','DEFAULT')
    case(2); call init_multiparameter('KR','LJ_MEOS','DEFAULT')
    case(3); call init_multiparameter('AR','LJ_MEOS','DEFAULT')
    case(4); call init_multiparameter('KR','LJTS_MEOS','DEFAULT')
    case(5); call init_multiparameter('NE','LJTS_MEOS','DEFAULT')
    end select
    active_model=>get_active_thermo_model();model_index=active_model%model_idx
    do mode=0,5
      count=0
      do i=1,n
        if(profiles(i)==profile.and.modes(i)==mode)then;count=count+1;selected(count)=i;endif
      enddo
      if(count==0)cycle
      if(trim(command)=='check')then
        do j=1,count
          idx=selected(j); call evaluate(idx,q)
          write(*,'(a,i0,a,i0,a,i0,a,*(es24.16,:,","))') 'check_',profile,'_',mode,'_',j-1,',',q
          write(*,'(a,i0,a,i0,a,i0,a,i0)') 'status_',profile,'_',mode,'_',j-1,',',last_status
        enddo
      else
        sums=0.; call system_clock(start,rate)
        do i=0,999
          idx=selected(mod(i,count)+1); call evaluate(idx,q); sums=sums+q
        enddo
        call system_clock(finish)
        write(*,'(a,i0,a,i0,a,es24.16,a,*(es24.16,:,","))') 'source_',profile,'_',mode,',', &
          real(finish-start)/real(rate)*1e9/1000.,',',sums
      endif
    enddo
  enddo
contains
  subroutine evaluate(i,q)
    integer,intent(in)::i
    real,intent(out)::q(8)
    integer::err,profile_case,kind_case
    real::tt,bb,bbl,vliq,vvap,hh,ss
    profile_case=profiles(i); kind_case=kinds(i);last_status=0
    call activate_model(model_index)
    t=values(1,i); p=values(2,i); spec=values(3,i); guess=values(4,i)
    tmin=values(5,i);tmax=values(6,i);z=1.;x=0.;y=0.;beta=0.;betal=0.;phase=TWOPH
    if(kind_case==1.or.kind_case==2)then
      tpTmin=tmin;tpTmax=tmax
    endif
    if(profile_case==5)then
      eps=39.83;vscale=6.02214076e23*(2.8e-10)**3;pscale=8.31446261815324*eps/vscale;hscale=8.31446261815324*eps;sscale=8.31446261815324
      if(kind_case==1.or.kind_case==2)then
        tpTmin=tmin*eps;tpTmax=tmax*eps
      endif
      t=t*eps;p=p*pscale;guess=guess*eps
      if(kind_case==1)spec=spec*hscale
      if(kind_case==2)spec=spec*sscale
    endif
    if(kind_case==1.or.kind_case==2)then
      call set_tmin(tpTmin);call set_tmax(tpTmax)
    endif
    if(kind_case==5)then
      vl=values(11,i);vv=values(12,i);if(profile_case==5)then;vl=vl*vscale;vv=vv*vscale;endif
      call enthalpy_tv(t,vl,z,hl);call entropy_tv(t,vl,z,sl)
      call enthalpy_tv(t,vv,z,hv);call entropy_tv(t,vv,z,s2)
      q=[t,p,real(SINGLEPH),-1.,hl,sl,hv,s2]
      if(profile_case==5)q([1,2,5,6,7,8])=q([1,2,5,6,7,8])/[eps,pscale,hscale,sscale,hscale,sscale]
      return
    endif
    if(kind_case==1.or.kind_case==2) t=guess
    if(kind_case==0)then
      call twoPhaseTPflash(t,p,z,beta,betal,phase,x,y)
    else if(kind_case==1)then
      call twoPhasePHflash(t,p,z,beta,betal,x,y,spec,phase,err)
      if(err>0.or.err < -1)error stop 'source PH error'
      last_status=err
    else
      call twoPhasePSflash(t,p,z,beta,betal,x,y,spec,phase,err)
      if(err>0.or.err < -1)error stop 'source PS error'
      last_status=err
    endif
    if(phase==TWOPH)then
      call specificVolume(t,p,x,LIQPH,vliq);call specificVolume(t,p,y,VAPPH,vvap)
      call enthalpy(t,p,x,LIQPH,hl);call enthalpy(t,p,y,VAPPH,hv)
      call entropy(t,p,x,LIQPH,sl);call entropy(t,p,y,VAPPH,sv)
      hh=(1.-beta)*hl+beta*hv;ss=(1.-beta)*sl+beta*sv
    else if(phase==VAPPH)then
      call specificVolume(t,p,z,VAPPH,vvap);vliq=0.
      call enthalpy(t,p,z,VAPPH,hh);call entropy(t,p,z,VAPPH,ss)
    else
      call specificVolume(t,p,z,LIQPH,vliq);vvap=0.
      call enthalpy(t,p,z,LIQPH,hh);call entropy(t,p,z,LIQPH,ss)
    endif
    q=[t,p,real(phase),beta,vliq,vvap,hh,ss]
    if(profile_case==5)q([1,2,5,6,7,8])=q([1,2,5,6,7,8])/[eps,pscale,vscale,vscale,hscale,sscale]
    if(.not.all(ieee_is_finite(q)))error stop 'nonfinite source flash'
  end subroutine
end program
