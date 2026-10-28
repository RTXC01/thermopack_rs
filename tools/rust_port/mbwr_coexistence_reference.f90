! Exact-source saturation calls and independent SI kernel checks.
! Copyright thermotools and Contributors. Apache-2.0.
program mbwr_coexistence_reference
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_value,ieee_quiet_nan
  use eoslibinit, only: init_multiparameter,redefine_critical_parameters
  use thermopack_var, only: get_active_thermo_model,thermo_model,Rgas,kRgas,tpTmin,tpPmin
  use eos_parameters, only: single_eos
  use eosdata, only: meosMbwr19,meosMbwr32
  use mbwr, only: eosmbwr,allocNIJL,fillIJL19,fillIJL32,computeZCoeff,computeHelmCoeff,makeParam,MBWR_pressure
  use mbwrdata, only: mbwr19Array,mbwr32Array
  use mbwr_additional, only: MBWR_Fres,MBWR_lnphi
  use saturation, only: safe_bubP,safe_bubT
  use eos, only: specificvolume
  use cubic_eos, only: cb_eos
  use cbselect, only: SelectCubicEOS,SelectMixingRules
  implicit none
  type(thermo_model),pointer :: model
  integer :: mode,record,refine,ierr,nroots,j
  character(len=24) :: comp
  real :: roots(64),rootp,rootk,rootphi(1),rootf(3),rootbound,bl,bv
  real :: t,p,vl,vv,params(4),active(3),x(1),y(1),nan,dpdt,fl(3),fv(3),pl,pg,kl,kg,phiL(1),phiV(1)
  read(*,*) mode,record,comp,t,p,vl,vv,params,active,refine
  nan=ieee_value(0.0,ieee_quiet_nan);dpdt=nan
  write(*,'(A)') 'INIT_BEGIN';flush(6)
  if (mode<=2) then
    if (record<=11) then
      call init_multiparameter(trim(comp),'MBWR32','DEFAULT')
    else
      call init_multiparameter(trim(comp),'MBWR19','DEFAULT')
    endif
  else
    call init_multiparameter('C1','MBWR32','DEFAULT')
  endif
  model=>get_active_thermo_model()
  if (mode>=3) then
    if (mode/=5) then
      deallocate(model%comps(1)%p_comp);allocate(model%comps(1)%p_comp)
      call model%comps(1)%p_comp%init_from_name(trim(comp),'DEFAULT',ierr)
      if (ierr/=0) error stop 'explicit component metadata failed'
      model%complist(1)=trim(comp)
      model%comps(1)%p_comp%tc=active(1)
      model%comps(1)%p_comp%pc=active(2)
      model%comps(1)%p_comp%zc=active(2)/(active(3)*Rgas*active(1))
      model%comps(1)%p_comp%acf=params(4)
    endif
    select type(eos=>model%eos(1)%p_eos)
    type is(single_eos)
      call coefficients(record,eos%mbwr_meos(1))
      eos%subeosidx=merge(meosMbwr32,meosMbwr19,record<=11)
      eos%mbwr_meos(1)%tc=params(1);eos%mbwr_meos(1)%pc=params(2)
      eos%mbwr_meos(1)%rc=params(3)*1e-3;eos%mbwr_meos(1)%acf=params(4)
      eos%mbwr_meos(1)%b_SRK=.08664*kRgas*params(1)/params(2)
      eos%mbwr_meos(1)%m_SRK=.480+1.547*params(4)-.176*params(4)**2
      eos%mbwr_meos(1)%a0_SRK=.42747*(kRgas*params(1))**2/params(2)
    class default
      error stop 'exact single_eos expected'
    end select
    if (mode/=5) then
      select type(alt=>model%cubic_eos_alternative(1)%p_eos)
      class is(cb_eos)
        call SelectCubicEOS(1,model%comps,alt,'CLASSIC','DEFAULT')
        call SelectMixingRules(1,model%comps,alt,'VDW','DEFAULT')
      class default
        error stop 'original cubic fallback expected'
      end select
      if (refine==1) then
        write(*,'(A)') 'REFINE_BEGIN';flush(6)
        call redefine_critical_parameters(.true.,Tc_in=active(1:1),vc_in=1/active(3:3))
      endif
    endif
  endif
  write(*,'(A)') 'INIT_RETURNED';flush(6)
  select type(eos=>model%eos(1)%p_eos)
  type is(single_eos)
    write(*,'(A,*(ES26.17E3,1X))') 'META=',model%comps(1)%p_comp%tc,model%comps(1)%p_comp%pc, &
      model%comps(1)%p_comp%zc,model%comps(1)%p_comp%acf,eos%mbwr_meos(1)%tc,eos%mbwr_meos(1)%pc, &
      eos%mbwr_meos(1)%rc*1000.,eos%mbwr_meos(1)%acf,tpTmin,tpPmin
    if (mode==0) stop
    x=1
    if (mode==1.or.mode==3) then
      write(*,'(A)') 'SAT_BEGIN';flush(6)
      p=safe_bubP(t,x,y,ierr)
    else if (mode==2.or.mode==4) then
      write(*,'(A)') 'SAT_BEGIN';flush(6)
      t=safe_bubT(p,x,y,ierr)
    else
      ierr=0
    endif
    write(*,'(A,I0,A,*(ES26.17E3,1X))') 'POINT_STATUS=',ierr,' TP=',t,p
    if (mode/=5.and.ieee_is_finite(t).and.ieee_is_finite(p).and.t>0.and.p>0) then
      call specificvolume(t,p,x,1,vl)
      call specificvolume(t,p,x,2,vv)
    endif
    if (ieee_is_finite(vl).and.ieee_is_finite(vv).and.vl>0.and.vv>0.and.t>0.and.p>0) then
      call kernel(eos%mbwr_meos(1),t,p,vl,pl,kl,phiL,fl,bl)
      call kernel(eos%mbwr_meos(1),t,p,vv,pg,kg,phiV,fv,bv)
      dpdt=(-Rgas*t*t*(fv(2)-fl(2))+p*(vv-vl))/(t*(vv-vl))
      write(*,'(A,*(ES26.17E3,1X))') 'CHECK=',pl,pg,kl,kg,phiL,phiV,fl(3)+log(1/vl),fv(3)+log(1/vv),fl(2),fv(2),bl,bv
    endif
    write(*,'(A,*(ES26.17E3,1X))') 'STATE=',t,p,vl,vv,dpdt
    if (mode==5) then
      read(*,*) nroots
      if (nroots<0.or.nroots>64) error stop 'invalid root count'
      if (nroots>0) read(*,*) roots(1:nroots)
      do j=1,nroots
        call kernel(eos%mbwr_meos(1),t,p,1/roots(j),rootp,rootk,rootphi,rootf,rootbound)
        write(*,'(A,*(ES26.17E3,1X))') 'ROOT=',roots(j),rootp,rootk,rootphi,rootf(3)+log(roots(j)),rootbound
      enddo
    endif
  class default
    error stop 'exact single_eos expected'
  end select
contains
  subroutine kernel(eos,t,p,v,pout,prho,phi,f,bound)
    type(eosmbwr),intent(in)::eos
    real,intent(in)::t,p,v
    real,intent(out)::pout,prho,phi(1),f(3),bound
    real::pr(16),ft(1),rho
    integer::k
    call makeParam(pr(1:1+eos%bplen+eos%belen),t,eos,0)
    call MBWR_pressure(1/(v*1000.),pr(1:1+eos%bplen+eos%belen),pout,prho)
    prho=prho*1e-3
    rho=1/(v*1000.);bound=0
    do k=1,eos%bplen
      bound=bound+abs(pr(k+1))*rho**k
    enddo
    do k=1,eos%belen
      bound=bound+exp(-pr(1)*rho**2)*abs(pr(k+1+eos%bplen))*rho**(2*k+1)
    enddo
    bound=64*epsilon(t)*bound
    call MBWR_Fres(eos,t,v,1.0,Fr=f(1),F_T=f(2),F_n=ft)
    f(3)=ft(1)
    call MBWR_lnphi(eos,t,p,v,1.0,phi)
  end subroutine
  subroutine coefficients(record,eos)
    integer,intent(in)::record
    type(eosmbwr),intent(inout)::eos
    call eos%dealloc()
    if (record>=1 .and. record<=11) then
      eos%EqNo=32; eos%Ipol=19; eos%Iexp=13
      eos%bplen=9; eos%belen=6; eos%Helm_poly_len=22; eos%HelmLength=40
      call allocNIJL(eos%PCoeff_rhoT,32)
      eos%PCoeff_rhoT%N=mbwr32Array(record)%coeff(1:32)
      eos%gamma=mbwr32Array(record)%coeff(33)
      call fillIJL32(eos%PCoeff_rhoT)
    else if (record>=12 .and. record<=65) then
      eos%EqNo=19; eos%Ipol=13; eos%Iexp=6
      eos%bplen=6; eos%belen=2; eos%Helm_poly_len=16; eos%HelmLength=22
      call allocNIJL(eos%PCoeff_rhoT,19)
      eos%PCoeff_rhoT%N=mbwr19Array(record-11)%coeff(1:19)
      eos%gamma=mbwr19Array(record-11)%coeff(20)
      call fillIJL19(eos%PCoeff_rhoT)
    else
      error stop 'invalid MBWR record'
    endif
    call computeZCoeff(eos)
    call computeHelmCoeff(eos)
  end subroutine
end program
