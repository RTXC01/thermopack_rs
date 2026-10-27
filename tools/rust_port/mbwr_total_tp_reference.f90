! Calls unchanged TVP routines and actual classical single_component TP dispatch.
program mbwr_total_tp_reference
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use eoslibinit, only: init_multiparameter
  use thermopack_var, only: get_active_thermo_model, thermo_model, Rgas,kRgas
  use eos_parameters, only: single_eos
  use eosdata, only: meosMbwr19, meosMbwr32
  use mbwr, only: eosmbwr,allocNIJL,fillIJL19,fillIJL32,computeZCoeff,computeHelmCoeff,makeParam,MBWR_density,MBWR_pressure
  use mbwrdata, only: mbwr19Array,mbwr32Array
  use mbwr_additional, only: MBWR_Hres,MBWR_Sres,MBWR_zfac,MBWR_lnphi
  use ideal, only: set_enthalpy_of_formation,set_standard_entropy,Cpideal,Hideal_mix,TP_Sideal_mix
  use eostv, only: Fres,Fideal,pressure,enthalpy_tvp,entropy_tvp,thermo_tvp,chemical_potential_tv,internal_energy_tv,entropy_tv,enthalpy_tv,free_energy_tv
  use speed_of_sound, only: speed_of_sound_TV
  use single_component, only: enthalpy_single,entropy_single
  implicit none
  type(thermo_model), pointer :: model
  integer :: cases,k,record,profile,ierr,mode,phase,found
  real :: t,v,n(1),p,pe,pt,pv,pn(1),tc,pc,rhoc,acf,rho,prho,param(16)
  real :: ht(4),st(4),phi(4),hc(4),sc(4),zc(4),fc(4),hi(4),si(4),nn(1,1)
  real :: fr(10),fi(10),f(10),pp(4),u(4),h(4),s(4),a(4),g(4),mu(4)
  real :: vals(45),cv,cp,rawsound,sound,cpflag,phiflag,soundflag,phival
  character(len=24) :: comp
  call init_multiparameter('C1','MBWR32','DEFAULT')
  model => get_active_thermo_model()
  read(*,*) cases
  do k=1,cases
    read(*,*) record,t,v,n,comp,profile,mode,p,phase,tc,pc,rhoc,acf
    deallocate(model%comps(1)%p_comp)
    allocate(model%comps(1)%p_comp)
    call model%comps(1)%p_comp%init_from_name(trim(comp),'DEFAULT',ierr)
    if (ierr/=0) error stop 'component metadata initialization failed'
    if (profile==1 .or. profile==3) then
      call set_enthalpy_of_formation(1,-12345.678d0)
      call set_standard_entropy(1,210.25d0,'1ATM')
    endif
    if (profile==2 .or. profile==3) then
      model%comps(1)%p_comp%id_cp%cptype=8
      model%comps(1)%p_comp%id_cp%cp=0
      model%comps(1)%p_comp%id_cp%cp(1)=3.5d0
    endif
    select type(eos => model%eos(1)%p_eos)
    type is(single_eos)
      call coefficients(record,eos%mbwr_meos(1))
      eos%subeosidx=merge(meosMbwr32,meosMbwr19,record<=11)
      if (mode==1) then
        eos%mbwr_meos(1)%tc=tc; eos%mbwr_meos(1)%pc=pc
        eos%mbwr_meos(1)%rc=rhoc*1e-3; eos%mbwr_meos(1)%acf=acf
        eos%mbwr_meos(1)%b_SRK=0.08664*kRgas*tc/pc
        eos%mbwr_meos(1)%m_SRK=0.480+1.547*acf-0.176*acf**2
        eos%mbwr_meos(1)%a0_SRK=0.42747*(kRgas*tc)**2/pc
        call makeParam(param(1:1+eos%mbwr_meos(1)%bplen+eos%mbwr_meos(1)%belen),t,eos%mbwr_meos(1),0)
        rho=MBWR_density(t,p,phase,param(1:1+eos%mbwr_meos(1)%bplen+eos%mbwr_meos(1)%belen),eos%mbwr_meos(1),found)
        call MBWR_pressure(rho,param(1:1+eos%mbwr_meos(1)%bplen+eos%mbwr_meos(1)%belen),pe,prho)
        v=n(1)/(rho*1000.)
        write(*,'(A,I0,A,*(ES26.17E3,1X))') 'CASE=',k,' DENSITY=',rho*1000.,real(found),pe,prho*1e-3
        call enthalpy_single(1,model%comps,eos,t,p,n,phase,.false.,hc(1),hc(2),hc(3),hc(4:4))
        call entropy_single(1,model%comps,eos,t,p,n,phase,.false.,sc(1),sc(2),sc(3),sc(4:4))
      else
        ! Component helpers supply independent ideal values. Keep all ideal
        ! derivatives separately; native comparison verifies the dispatcher
        ! defects using mode=1's actual single_component outputs.
        call Hideal_mix(1,model%comps,t,n,hi(1),hi(2),hi(3),hi(4:4))
        call TP_Sideal_mix(1,model%comps,t,p,n,si(1),si(2),si(3),si(4:4))
        call MBWR_Hres(eos%mbwr_meos(1),t,p,v,n(1),hc(1),hc(2),hc(3),hc(4:4))
        call MBWR_Sres(eos%mbwr_meos(1),t,p,v,n(1),sc(1),sc(2),sc(3),sc(4:4))
        hc(1:2)=hc(1:2)+hi(1:2);sc(1:2)=sc(1:2)+si(1:2)
      endif
      call MBWR_zfac(eos%mbwr_meos(1),t,p,v,n(1),zc(1),zc(2),zc(3),zc(4:4))
      call MBWR_lnphi(eos%mbwr_meos(1),t,p,v,n(1),fc(1:1),fc(2:2),fc(3:3),nn)
      fc(4)=nn(1,1)
      write(*,'(A,I0,A,*(ES26.17E3,1X))') 'CASE=',k,' CLASSICAL=',p,v,hc,sc,zc,fc
    class default
      error stop 'expected exact single_eos dynamic type'
    end select
    if (mode==0) then
    call jets(t,v,n,fr,fi)
    write(*,'(A,I0,A,*(ES26.17E3,1X))') 'CASE=',k,' RESIDUAL=',fr
    f=fr+fi
    pp(1)=pressure(t,v,n,dpdt=pp(2),dpdv=pp(3),dpdn=pp(4:4))
    call internal_energy_tv(t,v,n,u(1),u(2),u(3),u(4:4))
    call enthalpy_tv(t,v,n,h(1),h(2),h(3),h(4:4))
    call entropy_tv(t,v,n,s(1),s(2),s(3),s(4:4))
    call free_energy_tv(t,v,n,a(1),a(2),a(3),a(4:4))
    call chemical(t,v,n,mu)
    g=a+v*pp
    g(3)=a(3)+pp(1)+v*pp(3)
    cv=u(2)
    cpflag=0; cp=0
    if (abs(pp(3))>64*epsilon(t)*Rgas*t*n(1)/v**2) then
      cp=cv-t*pp(2)**2/pp(3)
      if (ieee_is_finite(cp)) cpflag=1
    endif
    if (cpflag==0) cp=0
    phiflag=0; phival=0
    if (pp(1)>0) then
      phival=fr(4)-log(pp(1)*v/(n(1)*Rgas*t))
      if (ieee_is_finite(phival)) phiflag=1
    endif
    if (phiflag==0) phival=0
    rawsound=speed_of_sound_TV(t,v,n)
    soundflag=0; sound=0
    if (ieee_is_finite(rawsound) .and. rawsound>0) then
      soundflag=1; sound=rawsound
    endif
    vals=[f,pp,u,h,s,a,g,mu,cv,cpflag,cp,phiflag,phival,soundflag,sound]
    write(*,'(A,I0,A,*(ES26.17E3,1X))') 'CASE=',k,' VALUES=',vals
    write(*,'(A,I0,A,*(ES26.17E3,1X))') 'CASE=',k,' IDEAL=', &
         model%comps(1)%p_comp%href_int,model%comps(1)%p_comp%sref_int, &
         model%comps(1)%p_comp%mw,real(model%comps(1)%p_comp%id_cp%cptype), &
         Cpideal(model%comps(1)%p_comp,1,t),fi
    write(*,'(A,I0,A,ES26.17E3)') 'CASE=',k,' RAW_SOUND=',rawsound
      pe=pressure(t,v,n,dpdt=pt,dpdv=pv,dpdn=pn)
      call enthalpy_tvp(t,v,n,ht(1),ht(2),ht(3),ht(4:4))
      call entropy_tvp(t,v,n,st(1),st(2),st(3),st(4:4))
      call thermo_tvp(t,v,n,phi(1:1),phi(2:2),phi(3:3),nn)
      phi(4)=nn(1,1)
      write(*,'(A,I0,A,*(ES26.17E3,1X))') 'CASE=',k,' TVP=',pe,ht,st,phi
    endif
  enddo
contains
  subroutine jets(t,v,n,fr,fi)
    real,intent(in)::t,v,n(1)
    real,intent(out)::fr(10),fi(10)
    real::rnn(1,1),inn(1,1)
    call Fres(t,v,n,F=fr(1),F_T=fr(2),F_V=fr(3),F_n=fr(4:4),F_TT=fr(5), &
         F_TV=fr(6),F_Tn=fr(7:7),F_VV=fr(8),F_Vn=fr(9:9),F_nn=rnn)
    fr(10)=rnn(1,1)
    call Fideal(t,v,n,F=fi(1),F_T=fi(2),F_V=fi(3),F_n=fi(4:4),F_TT=fi(5), &
         F_TV=fi(6),F_Tn=fi(7:7),F_VV=fi(8),F_Vn=fi(9:9),F_nn=inn)
    fi(10)=inn(1,1)
  end subroutine
  subroutine chemical(t,v,n,mu)
    real,intent(in)::t,v,n(1)
    real,intent(out)::mu(4)
    real::nn(1,1)
    call chemical_potential_tv(t,v,n,mu(1:1),mu(2:2),mu(3:3),nn)
    mu(4)=nn(1,1)
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
