subroutine query_meos_density(t,p,phase,reduced,values) bind(C,name='query_meos_density')
  use, intrinsic :: iso_c_binding, only: c_int,c_double
  use thermopack_var, only: get_active_eos,base_eos_param
  use eos_parameters, only: single_eos
  use multiparameter_base, only: meos
  use multiparameter_lj, only: meos_lj,init_LJTS
  implicit none
  real(c_double),value :: t,p
  integer(c_int),value :: phase,reduced
  real(c_double),intent(out) :: values(5)
  class(base_eos_param),pointer :: active
  class(meos),pointer :: fluid
  type(meos_lj),target :: local
  real :: rho,actual_p,dp
  integer :: found,status
  if (reduced==1) then
    call init_LJTS(local);fluid=>local
  else
    active=>get_active_eos()
    select type(eos=>active)
    class is(single_eos)
      fluid=>eos%nist(1)%meos
    class default
      error stop 'Expected a multiparameter model'
    end select
  endif
  call fluid%densitySolver(t,p,phase,rho,found,status)
  call fluid%mp_pressure(rho,t,actual_p,dp)
  values=[rho,real(found),real(status),actual_p,dp]
end subroutine query_meos_density

! Pure multiparameter kernel and constructor metadata, with source optional
! properties left non-finite instead of dropping their otherwise valid states.
subroutine query_meos(t,v,amount,reduced,values,metadata) bind(C,name='query_meos')
  use, intrinsic :: iso_c_binding, only: c_int,c_double
  use thermopack_var, only: get_active_eos,base_eos_param
  use eos_parameters, only: single_eos
  use multiparameter_base, only: meos
  use multiparameter_lj, only: meos_lj,init_LJTS
  implicit none
  real(c_double),value :: t,v,amount
  integer(c_int),value :: reduced
  real(c_double),intent(out) :: values(46),metadata(18)
  class(base_eos_param),pointer :: active
  class(meos),pointer :: fluid
  type(meos_lj),target :: local
  real :: a0(0:2,0:2),ar(0:2,0:2),rho,tau,delta,p,pr,pt,r
  real :: a,at,av,att,atv,avv
  real :: f,ft,fv,fn(1),ftt,ftv,ftn(1),fvv,fvn(1),fnn(1,1)
  if (reduced==1) then
    call init_LJTS(local)
    fluid=>local
  else
    active=>get_active_eos()
    select type(eos=>active)
    class is(single_eos)
      fluid=>eos%nist(1)%meos
    class default
      error stop 'Expected a multiparameter model'
    end select
  endif
  rho=1/v;tau=fluid%tc/t;delta=rho/fluid%rc;r=fluid%Rgas_meos
  call fluid%alpha0Derivs_taudelta(delta,tau,a0)
  call fluid%alphaResDerivs_taudelta(delta,tau,ar)
  values(1:6)=[a0(0,0),a0(0,1),a0(0,2),a0(1,0),a0(2,0),a0(1,1)]
  values(7:12)=[ar(0,0),ar(0,1),ar(0,2),ar(1,0),ar(2,0),ar(1,1)]
  call fluid%mp_pressure(rho,t,p,pr,pt)
  values(13:15)=[p,pt,-rho*rho*pr]
  call fluid%alphaDerivs_Tv(t,v,a,at,av,att,atv,avv)
  values(16)=-r*t*t*at
  values(17)=-r*(a+t*at)
  values(18)=values(16)+p*v
  values(19)=r*t*a
  values(20)=values(19)+p*v
  values(21)=fluid%cv(t,v)
  values(22)=fluid%cp(t,v)
  values(23)=ar(0,0)+ar(1,0)-log(1.+ar(1,0))
  values(24)=fluid%speed_of_sound(t,v)
  call fluid%calc_F(t,v*amount,[amount],f,ft,fv,fn,ftt,ftv,ftn,fvv,fvn,fnn)
  values(25:34)=[f,ft,fv,fn(1),ftt,ftv,ftn(1),fvv,fvn(1),fnn(1,1)]
  call fluid%calc_Fid(t,v*amount,[amount],f,ft,fv,fn,ftt,ftv,ftn,fvv,fvn,fnn)
  values(35:44)=[f,ft,fv,fn(1),ftt,ftv,ftn(1),fvv,fvn(1),fnn(1,1)]
  values(45:46)=[fluid%satDeltaEstimate(tau,1),fluid%satDeltaEstimate(tau,2)]
  metadata=0.
  metadata(1:11)=[fluid%tc,fluid%pc,fluid%rc,fluid%acf,fluid%molarMass,fluid%t_triple, &
      fluid%p_triple,fluid%rhoLiq_triple,fluid%rhoVap_triple,fluid%Rgas_fit,fluid%Rgas_meos]
  select type(eos=>fluid)
  class is(meos_lj)
    metadata(12:13)=[eos%maxT,eos%maxP]
    if (reduced==0) metadata(14:15)=[eos%sigma,eos%eps_divk]
    metadata(16:18)=[real(eos%upPol),real(eos%upExp),real(eos%upExpExp)]
  end select
end subroutine query_meos

! Direct FMT grid entry point; optional output combinations exercise all source branches.
subroutine query_fmt(model,ngrid,n,flags,values) bind(C,name='query_fmt')
  use, intrinsic :: iso_c_binding, only: c_int,c_double
  use fundamental_measure_theory, only: fmt_energy_density
  implicit none
  integer(c_int), value :: model,ngrid,flags
  real(c_double), intent(in) :: n(ngrid,6)
  real(c_double), intent(out) :: values(43,ngrid)
  real :: phi(ngrid),gradient(ngrid,6),hessian(ngrid,6,6)
  integer :: i,j,k,offset
  phi=0;gradient=0;hessian=0
  select case(flags)
  case(1)
    call fmt_energy_density(model,ngrid,6,n,phi=phi)
  case(2)
    call fmt_energy_density(model,ngrid,6,n,phi_n=gradient)
  case(3)
    call fmt_energy_density(model,ngrid,6,n,phi=phi,phi_n=gradient)
  case(4)
    call fmt_energy_density(model,ngrid,6,n,phi_nn=hessian)
  case(5)
    call fmt_energy_density(model,ngrid,6,n,phi=phi,phi_nn=hessian)
  case(6)
    call fmt_energy_density(model,ngrid,6,n,phi_n=gradient,phi_nn=hessian)
  case(7)
    call fmt_energy_density(model,ngrid,6,n,phi=phi,phi_n=gradient,phi_nn=hessian)
  case default
    error stop 'Unsupported FMT probe output flags'
  end select
  do k=1,ngrid
    values(1,k)=phi(k)
    values(2:7,k)=gradient(k,:)
    offset=8
    do i=1,6
      do j=1,6
        values(offset,k)=hessian(k,i,j)
        offset=offset+1
      enddo
    enddo
  enddo
end subroutine query_fmt

! Access source-model alpha initialization/values for differential qualification.
! Linked to the existing exact-source library; never duplicates its globals.
subroutine query_alpha(index,t,method,params,values) bind(C,name='query_alpha')
  use, intrinsic :: iso_c_binding, only: c_int,c_double
  use thermopack_var, only: get_active_eos,base_eos_param
  use cubic_eos, only: cb_eos,alpha_corr_db,get_alpha_db_idx_from_alpha_idx
  use cbAlpha, only: cbCalcAlphaTerm
  implicit none
  integer(c_int), value :: index
  real(c_double), value :: t
  integer(c_int), intent(out) :: method
  real(c_double), intent(out) :: params(3),values(3)
  class(base_eos_param), pointer :: active
  integer :: nparam,db_index
  active=>get_active_eos()
  select type(eos=>active)
  class is(cb_eos)
    method=eos%single(index)%alphaMethod
    db_index=get_alpha_db_idx_from_alpha_idx(method)
    nparam=alpha_corr_db(db_index)%n_param
    params=0.0
    params(1:nparam)=eos%single(index)%alphaParams(1:nparam)
    call cbCalcAlphaTerm(size(eos%single),eos,t)
    values=[eos%single(index)%alpha,eos%single(index)%dalphadT,eos%single(index)%d2alphadT2]
  class default
    error stop 'reference probe requires a cubic model'
  end select
end subroutine query_alpha

! Isolate tpInitAlphaCorr from eoslibinit's translation/quantum overrides.
subroutine initialize_alpha(name,nname,reference,nreference) bind(C,name='initialize_alpha')
  use, intrinsic :: iso_c_binding, only: c_int,c_char
  use thermopack_var, only: get_active_eos,get_active_comps,base_eos_param
  use cubic_eos, only: cb_eos
  use compdata, only: gendata_pointer
  use cbAlpha, only: tpInitAlphaCorr
  implicit none
  integer(c_int), value :: nname,nreference
  character(kind=c_char), intent(in) :: name(*),reference(*)
  character(len=nname) :: alpha_name
  character(len=nreference) :: parameter_reference
  class(base_eos_param), pointer :: active
  type(gendata_pointer), pointer :: comps(:)
  integer :: i
  do i=1,nname
    alpha_name(i:i)=name(i)
  enddo
  do i=1,nreference
    parameter_reference(i:i)=reference(i)
  enddo
  active=>get_active_eos()
  comps=>get_active_comps()
  select type(eos=>active)
  class is(cb_eos)
    call tpInitAlphaCorr(size(comps),comps,eos,alpha_name,parameter_reference)
  class default
    error stop 'reference probe requires a cubic model'
  end select
end subroutine initialize_alpha

subroutine query_translation(index,t,enabled,kind,params,values,temp_dep) bind(C,name='query_translation')
  use, intrinsic :: iso_c_binding, only: c_int,c_double
  use thermopack_var, only: get_active_eos,get_active_comps,base_eos_param
  use compdata, only: gendata_pointer
  implicit none
  integer(c_int), value :: index
  real(c_double), value :: t
  integer(c_int), intent(out) :: enabled,kind,temp_dep
  real(c_double), intent(out) :: params(6),values(3)
  class(base_eos_param), pointer :: active
  type(gendata_pointer), pointer :: comps(:)
  logical :: dependent
  active=>get_active_eos()
  comps=>get_active_comps()
  enabled=active%volumeShiftId
  associate(ci=>comps(index)%p_comp%cid)
    kind=ci%c_type
    params=[ci%ciA,ci%ciB,ci%ciC,ci%ciDD,ci%ciE,ci%ciF]
    call ci%get_vol_trs_c(t,values(1),values(2),values(3),dependent)
  end associate
  temp_dep=0
  if (dependent) temp_dep=1
end subroutine query_translation

subroutine query_beta(index,t,method,params,values,critical) bind(C,name='query_beta')
  use, intrinsic :: iso_c_binding, only: c_int,c_double
  use thermopack_var, only: get_active_eos,get_active_comps,base_eos_param
  use cubic_eos, only: cb_eos,cbBetaQuantumIdx
  use cbBeta, only: cbCalcBetaTerm
  use compdata, only: gendata_pointer
  implicit none
  integer(c_int), value :: index
  real(c_double), value :: t
  integer(c_int), intent(out) :: method
  real(c_double), intent(out) :: params(2),values(3),critical(3)
  class(base_eos_param), pointer :: active
  type(gendata_pointer), pointer :: comps(:)
  active=>get_active_eos()
  comps=>get_active_comps()
  select type(eos=>active)
  class is(cb_eos)
    method=eos%single(index)%betaMethod
    params=0.0
    if (method==cbBetaQuantumIdx) params=eos%single(index)%betaParams(1:2)
    call cbCalcBetaTerm(size(eos%single),eos,t)
    values=[eos%single(index)%beta,eos%single(index)%dbetadT,eos%single(index)%d2betadT2]
    critical=[comps(index)%p_comp%tc,comps(index)%p_comp%pc,comps(index)%p_comp%mw]
  class default
    error stop 'reference probe requires a cubic model'
  end select
end subroutine query_beta

subroutine query_rr(n,z,k,single_step,beta,betaL,x,y,status) bind(C,name='query_rr')
  use, intrinsic :: iso_c_binding, only: c_int,c_double
  use tp_solver, only: rr_solve
  implicit none
  integer(c_int), value :: n,single_step
  real(c_double), intent(in) :: z(n),k(n)
  real(c_double), intent(inout) :: beta,betaL,x(n),y(n)
  integer(c_int), intent(out) :: status
  logical :: success
  success=rr_solve(n,z,k,beta,x,y,single_step/=0,betaL)
  status=0
  if (success) status=1
end subroutine query_rr

subroutine query_rr_iteration(n,t,p,z,k,single_step,beta,betaL,x,y,kout,fugl,fugv,g,flags) &
    bind(C,name='query_rr_iteration')
  use, intrinsic :: iso_c_binding, only: c_int,c_double
  use tp_solver, only: rr_successive_substitution_iteration
  use thermopack_var, only: nc
  implicit none
  integer(c_int), value :: n,single_step
  real(c_double), value :: t,p
  real(c_double), intent(in) :: z(n),k(n)
  real(c_double), intent(inout) :: beta,betaL,x(n),y(n),kout(n),fugl(n),fugv(n),g
  integer(c_int), intent(out) :: flags(2)
  logical :: converged,success
  if (n/=nc) error stop 'reference RR iteration component mismatch'
  call rr_successive_substitution_iteration(t,p,z,k,single_step/=0,x,y,beta,kout,fugl,fugv, &
       g,converged,success,betaL)
  flags=0
  if (converged) flags(1)=1
  if (success) flags(2)=1
end subroutine query_rr_iteration

subroutine query_symmetric(n,a,b,ipiv,status) bind(C,name='query_symmetric')
  use, intrinsic :: iso_c_binding, only: c_int,c_double
  implicit none
  integer(c_int), value :: n
  real(c_double), intent(inout) :: a(n*n),b(n)
  integer(c_int), intent(out) :: ipiv(n),status
  real(c_double) :: matrix(n,n),work(n)
  integer :: info
  matrix=transpose(reshape(a,[n,n]))
  call dsytf2('U',n,matrix,n,ipiv,info)
  status=0
  if (info==0) then
    call dsytrs2('U',n,1,matrix,n,ipiv,b,n,work,info)
    if (info==0) status=1
  endif
  a=reshape(transpose(matrix),[n*n])
end subroutine query_symmetric

subroutine query_response(n,t,p,z,bv,bl,x,y,phase,values) bind(C,name='query_response')
  use, intrinsic :: iso_c_binding, only: c_int,c_double
  use state_functions, only: dhdt_twoPhase,dvdt_twoPhase,dvdp_twoPhase
  use eos, only: enthalpy,specificVolume
  implicit none
  integer(c_int), value :: n,phase
  real(c_double), value :: t,p,bv,bl
  real(c_double), intent(in) :: z(n),x(n),y(n)
  real(c_double), intent(out) :: values(3)
  real(c_double) :: h,v
  if (phase==0) then
    values(1)=dhdt_twoPhase(t,p,z,bv,bl,x,y)
    values(2)=dvdt_twoPhase(t,p,z,bv,bl,x,y)
    values(3)=dvdp_twoPhase(t,p,z,bv,bl,x,y)
  else
    call enthalpy(t,p,z,phase,h,dhdt=values(1))
    call specificVolume(t,p,z,phase,v,dvdt=values(2),dvdp=values(3))
  endif
end subroutine query_response

subroutine query_sv(n,t,p,z,bv,bl,x,y,s,v,phase,status) bind(C,name='query_sv')
  use, intrinsic :: iso_c_binding, only: c_int,c_double
  use sv_solver, only: twoPhaseSVflash
  implicit none
  integer(c_int), value :: n
  real(c_double), intent(inout) :: t,p,bv,bl,x(n),y(n)
  real(c_double), intent(in) :: z(n),s,v
  integer(c_int), intent(inout) :: phase
  integer(c_int), intent(out) :: status
  call twoPhaseSVflash(t,p,z,bv,bl,x,y,s,v,phase,status)
end subroutine query_sv

subroutine query_pets(t,v,n,sigma,eps,values) bind(C,name='query_pets')
  use, intrinsic :: iso_c_binding, only: c_double
  use thermopack_var, only: get_active_eos,base_eos_param
  use pets, only: PETS_eos,calc_delta_Ac_LJ
  use saft_interface, only: pets_set_pure_params
  use eostv, only: pressure,internal_energy_tv,entropy_tv,free_energy_tv,enthalpy_tv,chemical_potential_tv, &
       virial_coefficients,thermo_tvp,enthalpy_tvp
  use speed_of_sound, only: speed_of_sound_TV
  implicit none
  real(c_double), value :: t,v,n,sigma,eps
  real(c_double), intent(out) :: values(56)
  real(c_double) :: nn(1),h_tvp
  class(base_eos_param), pointer :: active
  nn=[n]
  call pets_set_pure_params(1,[sigma,eps])
  active=>get_active_eos()
  select type(eos=>active)
  class is(PETS_eos)
    call eos%calc_d_pets(t,values(1:1),values(2:2),values(3:3))
    call eos%alpha_PETS(n/v,t,nn,alp=values(4),alp_T=values(5),alp_rho=values(6), &
         alp_TT=values(7),alp_rhoT=values(8),alp_rhorho=values(9))
    call eos%alpha_pets_hs(n/v,t,nn,alp=values(10),alp_T=values(11),alp_rho=values(12), &
         alp_TT=values(13),alp_rhoT=values(14),alp_rhorho=values(15))
    call eos%alpha_disp(n/v,t,nn,alp=values(16),alp_T=values(17),alp_rho=values(18), &
         alp_TT=values(19),alp_rhoT=values(20),alp_rhorho=values(21))
    call eos%F_PETS_TVn(t,v,nn,F=values(22),F_T=values(23),F_V=values(24),F_n=values(25:25), &
         F_TT=values(26),F_TV=values(27),F_Tn=values(28:28),F_VV=values(29),F_Vn=values(30:30), &
         F_nn=values(31:31))
  class default
    error stop 'PeTS probe requires active PeTS model'
  end select
  call calc_delta_Ac_LJ(1,t,v,nn,sigma,eps,a=values(32),a_T=values(33),a_V=values(34), &
       a_n=values(35:35),a_TT=values(36),a_TV=values(37),a_Tn=values(38:38), &
       a_VV=values(39),a_Vn=values(40:40),a_nn=values(41:41))
  values(42)=pressure(t,v,nn,dpdt=values(43),dpdv=values(44),dpdn=values(45:45))
  call internal_energy_tv(t,v,nn,values(46),dudt=values(50))
  call entropy_tv(t,v,nn,values(47))
  call free_energy_tv(t,v,nn,values(48))
  call enthalpy_tv(t,v,nn,values(49))
  call chemical_potential_tv(t,v,nn,values(51:51))
  values(52)=speed_of_sound_TV(t,v,nn)
  call virial_coefficients(t,[1.0_c_double],values(53),values(54))
  values(55)=0.0
  if (values(42)>0.0) call thermo_tvp(t,v,nn,values(55:55))
  call enthalpy_tvp(t,v,nn,h_tvp,dhdt=values(56))
end subroutine query_pets

subroutine query_pets_potential(n,r,sigma,eps,values,forces) bind(C,name='query_pets_potential')
  use, intrinsic :: iso_c_binding, only: c_int,c_double
  use thermopack_var, only: get_active_eos,base_eos_param
  use pets, only: PETS_eos
  use saft_interface, only: pets_set_pure_params
  implicit none
  integer(c_int), value :: n
  real(c_double), value :: sigma,eps
  real(c_double), intent(in) :: r(n)
  real(c_double), intent(out) :: values(n),forces(n)
  class(base_eos_param), pointer :: active
  call pets_set_pure_params(1,[sigma,eps])
  active=>get_active_eos()
  select type(eos=>active)
  class is(PETS_eos)
    call eos%calc_potential_pets(n,r,values,forces)
  class default
    error stop 'PeTS potential probe requires active PeTS model'
  end select
end subroutine query_pets_potential

subroutine set_pets_parameters(sigma,eps) bind(C,name='set_pets_parameters')
  use, intrinsic :: iso_c_binding, only: c_double
  use saft_interface, only: pets_set_pure_params
  implicit none
  real(c_double), value :: sigma,eps
  call pets_set_pure_params(1,[sigma,eps])
end subroutine set_pets_parameters

subroutine query_pets_volume(t,p,sigma,eps,phase,values) bind(C,name='query_pets_volume')
  use, intrinsic :: iso_c_binding, only: c_int,c_double
  use thermopack_var, only: get_active_eos,base_eos_param
  use pets, only: PETS_eos
  use saft_interface, only: pets_set_pure_params,pc_saft_nonassoc_volume_solver
  use eostv, only: pressure
  implicit none
  real(c_double), value :: t,p,sigma,eps
  integer(c_int), value :: phase
  real(c_double), intent(out) :: values(4)
  class(base_eos_param), pointer :: active
  call pets_set_pure_params(1,[sigma,eps])
  active=>get_active_eos()
  select type(eos=>active)
  class is(PETS_eos)
    call pc_saft_nonassoc_volume_solver(1,eos,t,p,[1.0_c_double],phase,values(1))
  class default
    error stop 'PeTS volume probe requires active PeTS model'
  end select
  values(2)=pressure(t,values(1),[1.0_c_double],dpdt=values(3),dpdv=values(4))
end subroutine query_pets_volume

! Load original coefficients without unrelated critical-parameter iteration.
! alphar_derivatives and MBWR_coef use no critical metadata or ideal Cp data.
subroutine initialize_mbwr_record(record,model)
  use mbwr, only: eosmbwr,allocNIJL,fillIJL19,fillIJL32,computeZCoeff,computeHelmCoeff
  use mbwrdata, only: mbwr19Array,mbwr32Array
  implicit none
  integer, intent(in) :: record
  type(eosmbwr), intent(inout) :: model
  call model%dealloc()
  if (record >= 1 .and. record <= 11) then
    model%EqNo=32; model%Ipol=19; model%Iexp=13
    model%bplen=9; model%belen=6; model%Helm_poly_len=22; model%HelmLength=40
    call allocNIJL(model%PCoeff_rhoT,32)
    model%PCoeff_rhoT%N=mbwr32Array(record)%coeff(1:32)
    model%gamma=mbwr32Array(record)%coeff(33)
    call fillIJL32(model%PCoeff_rhoT)
  else if (record >= 12 .and. record <= 65) then
    model%EqNo=19; model%Ipol=13; model%Iexp=6
    model%bplen=6; model%belen=2; model%Helm_poly_len=16; model%HelmLength=22
    call allocNIJL(model%PCoeff_rhoT,19)
    model%PCoeff_rhoT%N=mbwr19Array(record-11)%coeff(1:19)
    model%gamma=mbwr19Array(record-11)%coeff(20)
    call fillIJL19(model%PCoeff_rhoT)
  else
    error stop 'invalid MBWR record'
  endif
  call computeZCoeff(model)
  call computeHelmCoeff(model)
end subroutine initialize_mbwr_record

subroutine query_mbwr_parameters(record,values) bind(C,name='query_mbwr_parameters')
  use, intrinsic :: iso_c_binding, only: c_int,c_double
  use mbwrdata, only: mbwr19Array,mbwr32Array
  implicit none
  integer(c_int),value :: record
  real(c_double),intent(out) :: values(34)
  values=0.0
  if (record<=11) then
    values(1)=32.0; values(2:34)=mbwr32Array(record)%coeff
  else
    values(1)=19.0; values(2:21)=mbwr19Array(record-11)%coeff
  endif
end subroutine query_mbwr_parameters

subroutine query_mbwr(record,t,v,n,values) bind(C,name='query_mbwr')
  use, intrinsic :: iso_c_binding, only: c_int,c_double
  use mbwr, only: eosmbwr,makeParam,MBWR_pressure
  use mbwr_additional, only: alphar_derivatives,MBWR_Fres,MBWR_press
  implicit none
  integer(c_int),value :: record
  real(c_double),value :: t,v,n
  real(c_double),intent(out) :: values(26)
  type(eosmbwr),save :: model
  integer,save :: previous=0
  real :: result(0:2),rho
  real,allocatable :: param(:)
  if (previous/=record) then
    call initialize_mbwr_record(record,model)
    previous=record
  endif
  rho=n/v
  call alphar_derivatives(result,t,rho,0,2,model)
  values(1)=result(0); values(3)=result(1); values(6)=result(2)
  call alphar_derivatives(result,t,rho,1,1,model)
  values(2)=result(0); values(5)=result(1)
  call alphar_derivatives(result,t,rho,2,0,model)
  values(4)=result(0)
  call MBWR_Fres(model,t,v,n,Fr=values(7),F_T=values(8),F_v=values(9),F_n=values(10:10), &
       F_TT=values(11),F_Tv=values(12),F_Tn=values(13:13),F_vv=values(14),F_vn=values(15:15),F_nn=values(16:16))
  call MBWR_press(model,t,v,n,p=values(17),dpdt=values(18),dpdv=values(19))
  values(20)=-values(19)*v/n
  allocate(param(1+model%bplen+model%belen))
  call makeParam(param,t,model,0)
  call MBWR_pressure(rho*1e-3,param,values(21),values(23),values(26))
  values(23)=values(23)*1e-3; values(26)=values(26)*1e-6
  call makeParam(param,t,model,1)
  call MBWR_pressure(rho*1e-3,param,values(22),values(25))
  values(25)=values(25)*1e-3
  call makeParam(param,t,model,2)
  call MBWR_pressure(rho*1e-3,param,values(24))
end subroutine query_mbwr

subroutine query_mbwr_integrated(record,values) bind(C,name='query_mbwr_integrated')
  use, intrinsic :: iso_c_binding, only: c_int,c_double
  use mbwr, only: eosmbwr
  implicit none
  integer(c_int),value :: record
  real(c_double),intent(out) :: values(40)
  type(eosmbwr) :: model
  call initialize_mbwr_record(record,model)
  values=0.0
  values(1:model%HelmLength)=model%redResHelmCoeff_redrhoT%N
end subroutine query_mbwr_integrated

subroutine query_mbwr_properties(record,t,v,n,rule,factor,values) bind(C,name='query_mbwr_properties')
  use, intrinsic :: iso_c_binding, only: c_int,c_double
  use mbwr, only: eosmbwr
  use mbwr_additional, only: MBWR_press,MBWR_Gres,MBWR_Hres,MBWR_Sres,MBWR_zfac,MBWR_lnphi
  implicit none
  integer(c_int),value :: record,rule
  real(c_double),value :: t,v,n,factor
  real(c_double),intent(out) :: values(18)
  type(eosmbwr),save :: model
  integer,save :: previous=0
  real :: p
  if (previous/=record) then
    call initialize_mbwr_record(record,model)
    previous=record
  endif
  call MBWR_press(model,t,v,n,p=p)
  if (rule==1) p=abs(p)*factor+1000.0
  values(1)=p
  call MBWR_Gres(model,t,p,v,n,values(2))
  call MBWR_Hres(model,t,p,v,n,values(3),values(4),values(5),values(6:6))
  call MBWR_Sres(model,t,p,v,n,values(7),values(8),values(9),values(10:10))
  call MBWR_zfac(model,t,p,v,n,values(11),values(12),values(13),values(14:14))
  call MBWR_lnphi(model,t,p,v,n,values(15:15),values(16:16),values(17:17),values(18:18))
end subroutine query_mbwr_properties

subroutine query_mbwr_metadata(record,values) bind(C,name='query_mbwr_metadata')
  use, intrinsic :: iso_c_binding, only: c_int,c_double
  use mbwrdata
  use thermopack_var, only: kRgas
  implicit none
  integer(c_int),value :: record
  real(c_double),intent(out) :: values(8)
  character(len=8) :: id
  integer :: idx
  values=0.0
  if (record<=11) then
    id=mbwr32Array(record)%ComId
    idx=getMBWR32critPropIndex(id)
    if (idx>0) values(2:4)=[tc32_computed(idx),pc32_computed(idx),rc32_computed(idx)]
  else
    id=mbwr19Array(record-11)%ComId
    idx=getMBWR19critPropIndex(id)
    if (idx>0) values(2:4)=[tc19_computed(idx),pc19_computed(idx),rc19_computed(idx)]
  endif
  if (idx>0) then
    values(1)=1.0
    values(5)=values(3)/(values(4)*kRgas*values(2))
    values(4)=values(4)*1000.0
  endif
  idx=getTriplePointPropIndex(id)
  if (idx>0) then
    values(6)=1.0
    values(7:8)=[tTriple_coolprop(idx),pTriple_coolprop(idx)]
  endif
end subroutine query_mbwr_metadata

subroutine query_mbwr_critical_jet(record,t,rho,values) bind(C,name='query_mbwr_critical_jet')
  use, intrinsic :: iso_c_binding, only: c_int,c_double
  use mbwr, only: eosmbwr,makeParam,MBWR_pressure
  implicit none
  integer(c_int),value :: record
  real(c_double),value :: t,rho
  real(c_double),intent(out) :: values(5)
  type(eosmbwr),save :: model
  integer,save :: previous=0
  real,allocatable :: param(:)
  real :: unused
  if (previous/=record) then
    call initialize_mbwr_record(record,model)
    previous=record
  endif
  allocate(param(1+model%bplen+model%belen))
  call makeParam(param,t,model,0)
  call MBWR_pressure(rho*1e-3,param,values(1),values(2),values(3))
  call makeParam(param,t,model,1)
  call MBWR_pressure(rho*1e-3,param,unused,values(4),values(5))
  values(2)=values(2)*1e-3; values(3)=values(3)*1e-6
  values(4)=values(4)*1e-3; values(5)=values(5)*1e-6
end subroutine query_mbwr_critical_jet

! Explicit metadata makes this a probe of the original density routine, without
! triggering readDbParameters' unrelated critical grid search on legacy fluids.
subroutine query_mbwr_density(record,t,p,phase,critical,acf,values) bind(C,name='query_mbwr_density')
  use, intrinsic :: iso_c_binding, only: c_int,c_double
  use mbwr, only: eosmbwr,makeParam,MBWR_density,MBWR_pressure
  use thermopack_var, only: kRgas
  implicit none
  integer(c_int),value :: record,phase
  real(c_double),value :: t,p,acf
  real(c_double),intent(in) :: critical(3)
  real(c_double),intent(out) :: values(4)
  type(eosmbwr) :: model
  integer :: found
  real :: rho
  real,allocatable :: param(:)
  call initialize_mbwr_record(record,model)
  model%tc=critical(1);model%pc=critical(2);model%rc=critical(3)*1e-3;model%acf=acf
  model%b_SRK=0.08664*kRgas*model%tc/model%pc
  model%m_SRK=0.480+1.547*acf-0.176*acf**2
  model%a0_SRK=0.42747*(kRgas*model%tc)**2/model%pc
  allocate(param(1+model%bplen+model%belen))
  call makeParam(param,t,model,0)
  rho=MBWR_density(t,p,phase,param,model,found)
  values(1)=rho*1000.0;values(2)=real(found)
  call MBWR_pressure(rho,param,values(3),values(4))
  values(4)=values(4)*1e-3
end subroutine query_mbwr_density

subroutine query_mbwr_initialization(record,values,status) bind(C,name='query_mbwr_initialization')
  use, intrinsic :: iso_c_binding, only: c_int,c_double
  use mbwr, only: eosmbwr,initializeMBWRmodel
  use mbwrdata, only: mbwr19Array,mbwr32Array
  use compdata_init, only: getCompDBindex
  implicit none
  integer(c_int),value :: record
  real(c_double),intent(out) :: values(8)
  integer(c_int),intent(out) :: status
  type(eosmbwr) :: model
  character(len=8) :: id
  if (record<=11) then
    id=mbwr32Array(record)%ComId
  else
    id=mbwr19Array(record-11)%ComId
  endif
  status=0
  if (getCompDBindex(id)==0) then
    ! The original unguarded compdb(0) access has undefined behavior.
    ! Report it explicitly instead of executing an out-of-bounds initializer.
    status=1;values=0.0;return
  endif
  if (record<=11) then
    call initializeMBWRmodel(mbwr32Array(record)%ComId,model,32)
  else
    call initializeMBWRmodel(mbwr19Array(record-11)%ComId,model,19)
  endif
  values=[model%tc,model%pc,model%rc*1000.0,model%zc,model%acf,model%b_SRK*1e-3,model%m_SRK,model%a0_SRK*1e-6]
end subroutine query_mbwr_initialization

subroutine query_mbwr_critical_grid(record,initial,values) bind(C,name='query_mbwr_critical_grid')
  use, intrinsic :: iso_c_binding, only: c_int,c_double
  use mbwr, only: eosmbwr,MBWR_criticalParameters
  implicit none
  integer(c_int),value :: record
  real(c_double),intent(in) :: initial(3)
  real(c_double),intent(out) :: values(3)
  type(eosmbwr) :: model
  real :: t,p,r
  call initialize_mbwr_record(record,model)
  t=initial(1);p=initial(2);r=initial(3)*1e-3
  call MBWR_criticalParameters(t,p,r,model)
  values=[t,p,r*1000.0]
end subroutine query_mbwr_critical_grid

! Evaluate the original UNIFAC excess Gibbs energy and its analytical derivatives.
! The component list is installed in the source globals for the duration of the
! call and restored afterwards, so the probe never disturbs a live model.
subroutine query_unifac(ncomponents,names,variant,t,composition,values,count) &
     bind(C,name='query_unifac')
  use, intrinsic :: iso_c_binding, only: c_int,c_double,c_char
  use thermopack_var, only: complist,nc
  use unifac, only: init_unifac, unifdb, Ge_UNIFAC_GH_SG
  implicit none
  integer(c_int), value :: ncomponents
  character(kind=c_char), intent(in) :: names(*)
  integer(c_int), value :: variant
  real(c_double), value :: t
  real(c_double), intent(in) :: composition(*)
  real(c_double), intent(out) :: values(128)
  integer(c_int), intent(out) :: count
  ! Locals
  integer, parameter :: nmax=8
  character(len=512) :: namelist
  character(len=20) :: mrulestr
  character(len=20), pointer :: saved(:)
  integer :: saved_nc, i, j, n
  logical :: had_complist
  real :: ge,dgedt,d2gedt2
  real, allocatable :: dgedn(:),d2gedndt(:),z(:)
  real, allocatable :: d2gedndn(:,:)
  n = int(ncomponents)
  if (n < 1 .or. n > nmax) error stop 'query_unifac: unsupported component count'
  ! The source dummies have the component count as an explicit shape. The
  ! local arrays must match it exactly: a larger actual argument is
  ! associated in array element order and would scramble the matrix.
  allocate(dgedn(n),d2gedndt(n),z(n),d2gedndn(n,n))
  do i=1,512
    namelist(i:i)=names(i)
  enddo
  select case (int(variant))
  case (0)
    mrulestr='UNIFAC'
  case (1)
    mrulestr='UMR'
  case (2)
    mrulestr='VTPR'
  case default
    error stop 'query_unifac: unknown variant'
  end select
  do i=1,n
    z(i)=composition(i)
  enddo
  saved_nc=nc
  had_complist=associated(complist)
  if (had_complist) saved=>complist
  nullify(complist)
  allocate(complist(n))
  do i=1,n
    complist(i)=namelist((i-1)*20+1:i*20)
  enddo
  nc=n
  call init_unifac(unifdb,mrulestr)
  call Ge_UNIFAC_GH_SG(z,t,unifdb,ge,dgedt,d2gedt2,dgedn,d2gedndt,d2gedndn)
  deallocate(complist)
  if (had_complist) complist=>saved
  nc=saved_nc
  values(1)=ge
  values(2)=dgedt
  values(3)=d2gedt2
  do i=1,n
    values(3+i)=dgedn(i)
    values(3+n+i)=d2gedndt(i)
  enddo
  do i=1,n
    do j=1,n
      values(3+2*n+(i-1)*n+j)=d2gedndn(i,j)
    enddo
  enddo
  count=3+2*n+n*n
  deallocate(dgedn,d2gedndt,z,d2gedndn)
end subroutine query_unifac

! Evaluate the original melting or sublimation line correlation.
! The component list is installed in the source globals for the call and
! restored afterwards, so the probe never disturbs a live model.
subroutine query_solid_line(kind,component,tlimit,npoints,temperatures,pressures,count) &
     bind(C,name='query_solid_line')
  use, intrinsic :: iso_c_binding, only: c_int,c_double,c_char
  use thermopack_var, only: complist,nc,nce
  use solid_saturation, only: melting_pressure_correlation, sublimation_pressure_correlation
  implicit none
  integer(c_int), value :: kind
  character(kind=c_char), intent(in) :: component(*)
  real(c_double), value :: tlimit
  integer(c_int), value :: npoints
  real(c_double), intent(out) :: temperatures(*)
  real(c_double), intent(out) :: pressures(*)
  integer(c_int), intent(out) :: count
  ! Locals
  integer, parameter :: nmax = 64
  character(len=20) :: component_name
  character(len=20), pointer :: saved(:)
  integer :: saved_nc, saved_nce, i, ierr, n
  logical :: had_complist
  real, dimension(nmax) :: t_line, p_line
  n = int(npoints)
  if (n < 1 .or. n > nmax) error stop 'query_solid_line: unsupported point count'
  do i=1,20
    component_name(i:i)=component(i)
  enddo
  saved_nc = nc
  saved_nce = nce
  had_complist = associated(complist)
  if (had_complist) saved=>complist
  nullify(complist)
  allocate(complist(1))
  complist(1) = component_name
  nc = 1
  nce = 1
  if (int(kind) == 0) then
    call melting_pressure_correlation(tlimit, 1, .false., n, t_line, p_line, ierr)
  else
    call sublimation_pressure_correlation(tlimit, 1, .false., n, t_line, p_line, ierr)
  endif
  deallocate(complist)
  if (had_complist) complist=>saved
  nc = saved_nc
  nce = saved_nce
  do i=1,n
    temperatures(i) = t_line(i)
    pressures(i) = p_line(i)
  enddo
  count = n
end subroutine query_solid_line

! Shift the reference levels of the solid carbon dioxide model.
subroutine query_solid_co2_init(liquid_entropy,liquid_gibbs,temperature,pressure) &
     bind(C,name='query_solid_co2_init')
  use, intrinsic :: iso_c_binding, only: c_double
  use co2_gibbs, only: sco2init
  implicit none
  real(c_double), value :: liquid_entropy
  real(c_double), value :: liquid_gibbs
  real(c_double), value :: temperature
  real(c_double), value :: pressure
  call sco2init(liquid_entropy,liquid_gibbs,temperature,pressure)
end subroutine query_solid_co2_init

! Evaluate the solid carbon dioxide model with its current reference levels.
subroutine query_solid_co2(temperature,pressure,values,count) &
     bind(C,name='query_solid_co2')
  use, intrinsic :: iso_c_binding, only: c_double,c_int
  use co2_gibbs, only: sco2_gibbs, sco2_dgdp, sco2_d2gdp2, sco2_dgdT, sco2_d2gdT2, &
       sco2_d2gdtdp, sco2_specvol, sco2_entropy, sco2_enthalpy, sco2_energy, &
       sco2_helmholtz, sco2_heat_capacity, sco2_speed_of_sound
  implicit none
  real(c_double), value :: temperature
  real(c_double), value :: pressure
  real(c_double), intent(out) :: values(13)
  integer(c_int), intent(out) :: count
  values(1) = sco2_gibbs(temperature,pressure)
  values(2) = sco2_dgdp(temperature,pressure)
  values(3) = sco2_d2gdp2(temperature,pressure)
  values(4) = sco2_dgdT(temperature,pressure)
  values(5) = sco2_d2gdT2(temperature,pressure)
  values(6) = sco2_d2gdtdp(temperature,pressure)
  values(7) = sco2_specvol(temperature,pressure)
  values(8) = sco2_entropy(temperature,pressure)
  values(9) = sco2_enthalpy(temperature,pressure)
  values(10) = sco2_energy(temperature,pressure)
  values(11) = sco2_helmholtz(temperature,pressure)
  values(12) = sco2_heat_capacity(temperature,pressure)
  values(13) = sco2_speed_of_sound(temperature,pressure)
  count = 13
end subroutine query_solid_co2

subroutine query_solid_h2o_init(sl,gl) bind(C,name='query_solid_h2o_init')
  use, intrinsic :: iso_c_binding, only:c_double
  use h2o_gibbs, only:sho2_init
  implicit none
  real(c_double),value :: sl,gl
  call sho2_init(sl,gl)
end subroutine query_solid_h2o_init

subroutine query_solid_h2o(t,p,values) bind(C,name='query_solid_h2o')
  use, intrinsic :: iso_c_binding, only:c_double
  use h2o_gibbs, only:sh2o_gibbs,sh2o_dgdp,sh2o_d2gdp2,sh2o_dgdT,sh2o_d2gdT2,sh2o_d2gdtdp, &
       sh2o_specvol,sh2o_entropy,sh2o_enthalpy,sh2o_energy,sh2o_helmholtz,sh2o_heat_capacity
  implicit none
  real(c_double),value :: t,p
  real(c_double),intent(out) :: values(12)
  values=[sh2o_gibbs(t,p),sh2o_dgdp(t,p),sh2o_d2gdp2(t,p),sh2o_dgdT(t,p),sh2o_d2gdT2(t,p), &
       sh2o_d2gdtdp(t,p),sh2o_specvol(t,p),sh2o_entropy(t,p),sh2o_enthalpy(t,p),sh2o_energy(t,p), &
       sh2o_helmholtz(t,p),sh2o_heat_capacity(t,p)]
end subroutine query_solid_h2o

! Solve the mu-T and lnf-T problems with the original solver.
! The component list and the initial density are installed in the source
! globals for the call and restored afterwards.
subroutine query_mut_solve(ncomponents,kind,components,densities,temperature,mu,result,status) &
     bind(C,name='query_mut_solve')
  use, intrinsic :: iso_c_binding, only: c_int,c_double,c_char
  use thermopack_var, only: complist,nc,nce
  use mut_solver, only: solve_mu_t, solve_lnf_t
  implicit none
  integer(c_int), value :: ncomponents
  integer(c_int), value :: kind
  character(kind=c_char), intent(in) :: components(*)
  real(c_double), intent(in) :: densities(*)
  real(c_double), value :: temperature
  real(c_double), intent(in) :: mu(*)
  real(c_double), intent(out) :: result(*)
  integer(c_int), intent(out) :: status
  ! Locals
  integer, parameter :: nmax=8
  character(len=20) :: component_name
  character(len=20), pointer :: saved(:)
  integer :: saved_nc, saved_nce, i, n, ierr
  logical :: had_complist
  real :: rho(nmax)
  n = int(ncomponents)
  if (n < 1 .or. n > nmax) error stop 'query_mut_solve: unsupported component count'
  do i=1,20
    component_name(i:i)=components(i)
  enddo
  do i=1,n
    rho(i)=densities(i)
  enddo
  saved_nc = nc
  saved_nce = nce
  had_complist = associated(complist)
  if (had_complist) saved=>complist
  nullify(complist)
  allocate(complist(n))
  do i=1,n
    complist(i)=component_name((i-1)*20+1:i*20)
  enddo
  nc = n
  nce = n
  if (int(kind) == 0) then
    call solve_mu_t(mu(1:n),temperature,rho,ierr)
  else
    call solve_lnf_t(mu(1:n),temperature,rho,ierr)
  endif
  deallocate(complist)
  if (had_complist) complist=>saved
  nc = saved_nc
  nce = saved_nce
  do i=1,n
    result(i)=rho(i)
  enddo
  status = ierr
end subroutine query_mut_solve

! Evaluate the original chemical potential or log fugacity coefficients at a
! given density, to build the target of a mu-T or lnf-T solve.
subroutine query_mut_target(ncomponents,temperature,components,densities,kind,result,count) &
     bind(C,name='query_mut_target')
  use, intrinsic :: iso_c_binding, only: c_int,c_double,c_char
  use thermopack_var, only: complist,nc,nce
  use eosTV, only: thermo_tv, chemical_potential_tv
  implicit none
  integer(c_int), value :: ncomponents
  real(c_double), value :: temperature
  character(kind=c_char), intent(in) :: components(*)
  real(c_double), intent(in) :: densities(*)
  integer(c_int), value :: kind
  real(c_double), intent(out) :: result(*)
  integer(c_int), intent(out) :: count
  ! Locals
  integer, parameter :: nmax=8
  character(len=20) :: component_name
  character(len=20), pointer :: saved(:)
  integer :: saved_nc, saved_nce, i, n
  logical :: had_complist
  real :: rho(nmax), lnf(nmax), mu(nmax)
  n = int(ncomponents)
  if (n < 1 .or. n > nmax) error stop 'query_mut_target: unsupported component count'
  do i=1,20
    component_name(i:i)=components(i)
  enddo
  do i=1,n
    rho(i)=densities(i)
  enddo
  saved_nc = nc
  saved_nce = nce
  had_complist = associated(complist)
  if (had_complist) saved=>complist
  nullify(complist)
  allocate(complist(n))
  do i=1,n
    complist(i)=component_name((i-1)*20+1:i*20)
  enddo
  nc = n
  nce = n
  call thermo_tv(temperature,1.0,rho,lnf)
  if (int(kind) == 0) then
    call chemical_potential_tv(temperature,1.0,rho,mu)
    do i=1,n
      result(i)=mu(i)
    enddo
  else
    do i=1,n
      result(i)=lnf(i)
    enddo
  endif
  deallocate(complist)
  if (had_complist) complist=>saved
  nc = saved_nc
  nce = saved_nce
  count = n
end subroutine query_mut_target

! Evaluate the original hyperdual elementary functions.
! The inputs are the value, the three first order derivatives and the three
! mixed derivatives, followed by the overall differentiation order. The result
! holds the eight components of twelve functions.
subroutine query_hyperdual(inputs,values,count) bind(C,name='query_hyperdual')
  use, intrinsic :: iso_c_binding, only: c_int,c_double
  use hyperdual_mod, only: hyperdual, exp, log, sqrt, sin, cos, tan, sinh, cosh, tanh, &
       asin, acos, atan
  implicit none
  real(c_double), intent(in) :: inputs(9)
  real(c_double), intent(out) :: values(96)
  integer(c_int), intent(out) :: count
  ! Locals
  type(hyperdual) :: x
  integer :: i
  x%f0=inputs(1);x%f1=inputs(2);x%f2=inputs(3);x%f3=inputs(4)
  x%f12=inputs(5);x%f13=inputs(6);x%f23=inputs(7);x%f123=inputs(8)
  x%order=int(inputs(9))
  call store(0,exp(x))
  call store(1,log(x))
  call store(2,sqrt(x))
  call store(3,sin(x))
  call store(4,cos(x))
  call store(5,tan(x))
  call store(6,sinh(x))
  call store(7,cosh(x))
  call store(8,tanh(x))
  call store(9,asin(x))
  call store(10,acos(x))
  call store(11,atan(x))
  count=96
contains
  subroutine store(offset,v)
    integer, intent(in) :: offset
    type(hyperdual), intent(in) :: v
    values(offset*8+1)=v%f0
    values(offset*8+2)=v%f1
    values(offset*8+3)=v%f2
    values(offset*8+4)=v%f3
    values(offset*8+5)=v%f12
    values(offset*8+6)=v%f13
    values(offset*8+7)=v%f23
    values(offset*8+8)=v%f123
  end subroutine store
end subroutine query_hyperdual

! Locate the state where the solid and the fluid phase of the original
! coexist. The component list is installed in the source globals for the
! call and restored afterwards.
subroutine query_solid_saturation(component,temperature,pressure,phase, &
     specification,ln_spec,result,status) bind(C,name='query_solid_saturation')
  use, intrinsic :: iso_c_binding, only: c_int,c_double,c_char
  use thermopack_var, only: complist,nc,nce
  use solideos, only: initDryIce
  use solid_saturation, only: solidPointOnEnvelope
  implicit none
  character(kind=c_char), intent(in) :: component(*)
  real(c_double), value :: temperature, pressure
  integer(c_int), value :: phase
  integer(c_int), value :: specification
  real(c_double), value :: ln_spec
  real(c_double), intent(out) :: result(2)
  integer(c_int), intent(out) :: status
  ! Locals
  character(len=20) :: component_name
  character(len=20), pointer :: saved(:)
  integer :: saved_nc, saved_nce, ierr, iterations, i
  logical :: had_complist
  real :: t, p, beta_sol
  real, dimension(1) :: z
  do i=1,20
    component_name(i:i)=component(i)
  enddo
  saved_nc = nc
  saved_nce = nce
  had_complist = associated(complist)
  if (had_complist) saved=>complist
  nullify(complist)
  allocate(complist(1))
  complist(1)=component_name
  nc = 1
  nce = 1
  call initDryIce()
  z = 1.0
  t = temperature
  p = pressure
  beta_sol = 0.0
  iterations = solidPointOnEnvelope(z,t,p,phase,beta_sol,1,specification,ln_spec,ierr)
  deallocate(complist)
  if (had_complist) complist=>saved
  nc = saved_nc
  nce = saved_nce
  result(1)=t
  result(2)=p
  status=ierr
end subroutine query_solid_saturation

! Evaluate the matrix of the second virial coefficients of the original.
subroutine query_second_virial_matrix(temperature,values,count) &
     bind(C,name='query_second_virial_matrix')
  use, intrinsic :: iso_c_binding, only: c_int,c_double
  use thermopack_var, only: nc
  use eosTV, only: secondVirialCoeffMatrix
  implicit none
  real(c_double), value :: temperature
  real(c_double), intent(out) :: values(16)
  integer(c_int), intent(out) :: count
  ! Locals
  real :: Bmat(nc,nc)
  integer :: i,j
  call secondVirialCoeffMatrix(temperature,Bmat)
  count = nc*nc
  do j=1,nc
    do i=1,nc
      values((j-1)*nc+i) = Bmat(i,j)
    end do
  end do
end subroutine query_second_virial_matrix

! Evaluate one polynomial of the spline of the excess energy of the original,
! with its partial derivatives.
subroutine query_ljs_spline(t_star,rho_star,index,values,count) &
     bind(C,name='query_ljs_spline')
  use, intrinsic :: iso_c_binding, only: c_int,c_double
  use lj_splined, only: calc_ai_LJs_ex, ljs_bh_get_pure_params, &
       ljs_bh_get_bh_diameter_div_sigma
  implicit none
  real(c_double), value :: t_star, rho_star
  integer(c_int), value :: index
  real(c_double), intent(out) :: values(6)
  integer(c_int), intent(out) :: count
  ! Locals
  real :: sigma, eps, d_bh, tau
  real :: a, a_r, a_x, a_rr, a_xx, a_rx
  call ljs_bh_get_pure_params(sigma,eps)
  call ljs_bh_get_bh_diameter_div_sigma(t_star,d_bh)
  tau = 1.0 - d_bh
  call calc_ai_LJs_ex(index+1,tau,rho_star,1.0,a,a_r,a_x,a_rr,a_xx,a_rx)
  values(1)=a
  values(2)=a_r
  values(3)=a_x
  values(4)=a_rr
  values(5)=a_xx
  values(6)=a_rx
  count=6
end subroutine query_ljs_spline

! Evaluate the reduced residual Helmholtz energy of the splined model of the
! original.
subroutine query_ljs_helmholtz(temperature,volume,moles,values,count) &
     bind(C,name='query_ljs_helmholtz')
  use, intrinsic :: iso_c_binding, only: c_int,c_double
  use thermopack_var, only: nc, get_active_eos, base_eos_param
  use lj_splined, only: ljs_bh_eos, calcFresLJs_bh
  implicit none
  real(c_double), value :: temperature, volume, moles
  real(c_double), intent(out) :: values(3)
  integer(c_int), intent(out) :: count
  ! Locals
  class(base_eos_param), pointer :: eos
  real :: n(1)
  real :: F, F_T
  real, dimension(1) :: F_n
  n(1) = moles
  F = 0.0
  F_T = 0.0
  F_n = 0.0
  eos => get_active_eos()
  select type(p_eos => eos)
  class is(ljs_bh_eos)
    call calcFresLJs_bh(p_eos,nc,temperature,volume,n,F,F_T=F_T,F_n=F_n)
  end select
  values(1) = F
  values(2) = F_T
  values(3) = F_n(1)
  count = 3
end subroutine query_ljs_helmholtz

! Evaluate the reduced interaction potential of the original.
subroutine query_ljs_potential(count,radii,values,numbers) &
     bind(C,name='query_ljs_potential')
  use, intrinsic :: iso_c_binding, only: c_int,c_double
  use lj_splined, only: ljs_potential_reduced
  implicit none
  integer(c_int), value :: count
  real(c_double), intent(in) :: radii(32)
  real(c_double), intent(out) :: values(64)
  integer(c_int), intent(out) :: numbers
  ! Locals
  real :: pot(32), force(32)
  integer :: i
  pot = 0.0
  force = 0.0
  call ljs_potential_reduced(count,radii,pot,force)
  numbers = 2*count
  do i=1,count
    values(2*i-1) = pot(i)
    values(2*i) = force(i)
  end do
end subroutine query_ljs_potential

! Evaluate one perturbation term of the White Bear family of the original.
subroutine query_ljs_wca_term(index,t_star,rho_star,values,count) &
     bind(C,name='query_ljs_wca_term')
  use, intrinsic :: iso_c_binding, only: c_int,c_double
  use lj_splined, only: calc_LJs_WCA_ai_tr
  implicit none
  integer(c_int), value :: index
  real(c_double), value :: t_star, rho_star
  real(c_double), intent(out) :: values(1)
  integer(c_int), intent(out) :: count
  ! Locals
  real :: ai
  ai = 0.0
  call calc_LJs_WCA_ai_tr(index,t_star,rho_star,ai)
  values(1) = ai
  count = 1
end subroutine query_ljs_wca_term

! Evaluate the cavity terms of the White Bear reference of the original.
subroutine query_white_bear(temperature,eps_divk,eta,r_star,r_star_t,r_star_tt,values,count) &
     bind(C,name='query_white_bear')
  use, intrinsic :: iso_c_binding, only: c_int,c_double
  use hardsphere_wca, only: cavity_distribution_function, lj_boltzmann_factor, &
       minimum_r_cavity_int_LJ
  implicit none
  real(c_double), value :: temperature, eps_divk, eta, r_star, r_star_t, r_star_tt
  real(c_double), intent(out) :: values(12)
  integer(c_int), intent(out) :: count
  ! Locals
  real :: r_min, r_min_t, r_min_tt
  real :: y_hs, y_hs_e, y_hs_ee, y_hs_t, y_hs_tt, y_hs_et
  real :: e_lj, e_lj_t, e_lj_tt
  call minimum_r_cavity_int_LJ(temperature,eps_divk,r_min,r_min_t,r_min_tt)
  call cavity_distribution_function(eta,r_star,r_star_t,r_star_tt,y_hs,y_hs_e,y_hs_ee, &
       y_hs_t,y_hs_tt,y_hs_et)
  call lj_boltzmann_factor(r_star,r_star_t,r_star_tt,temperature,eps_divk,e_lj,e_lj_t,e_lj_tt)
  values(1) = r_min
  values(2) = r_min_t
  values(3) = r_min_tt
  values(4) = y_hs
  values(5) = y_hs_e
  values(6) = y_hs_ee
  values(7) = y_hs_t
  values(8) = y_hs_tt
  values(9) = y_hs_et
  values(10) = e_lj
  values(11) = e_lj_t
  values(12) = e_lj_tt
  count = 12
end subroutine query_white_bear

! Evaluate the nodes and the weights of a quadrature of the original.
subroutine query_quadrature(quad,values,count) bind(C,name='query_quadrature')
  use, intrinsic :: iso_c_binding, only: c_int,c_double
  use quadratures, only: get_quadrature_positions, get_quadrature_weights
  use quadratures, only: max_n_quadrature
  implicit none
  integer(c_int), value :: quad
  real(c_double), intent(out) :: values(2*max_n_quadrature)
  integer(c_int), intent(out) :: count
  ! Locals
  real :: x(max_n_quadrature), w(max_n_quadrature)
  integer :: n, i
  call get_quadrature_positions(quad,x,n)
  call get_quadrature_weights(quad,w,n)
  count = 2*n
  do i=1,n
    values(i) = x(i)
    values(n+i) = w(i)
  end do
end subroutine query_quadrature

! Evaluate the integral and the differential of a test integrand with a
! quadrature of the original.
subroutine query_quadrature_integral(quad,lower,upper,values,count) &
     bind(C,name='query_quadrature_integral')
  use, intrinsic :: iso_c_binding, only: c_int,c_double
  use quadratures, only: integrate_with_differential
  implicit none
  integer(c_int), value :: quad
  real(c_double), value :: lower, upper
  real(c_double), intent(out) :: values(2)
  integer(c_int), intent(out) :: count
  ! Locals
  real :: integral, differential
  real :: param(1)
  param(1) = 0.0
  call integrate_with_differential(quad,test_integrand,1,param,lower,upper,lower,upper, &
       .false.,integral,differential)
  values(1) = integral
  values(2) = differential
  count = 2
contains
  subroutine test_integrand(r,n,param,fun,fun_T,fun_r)
    real, intent(in) :: r
    integer, intent(in) :: n
    real, intent(in) :: param(n)
    real, intent(out) :: fun, fun_T, fun_r
    fun = r*r
    fun_T = 2.0*r
    fun_r = 2.0
  end subroutine test_integrand
end subroutine query_quadrature_integral

! Evaluate the matrix of the third virial coefficients of the original.
subroutine query_binary_third_virial_matrix(temperature,values,count) &
     bind(C,name='query_binary_third_virial_matrix')
  use, intrinsic :: iso_c_binding, only: c_int,c_double
  use thermopack_var, only: nc
  use eosTV, only: binaryThirdVirialCoeffMatrix
  implicit none
  real(c_double), value :: temperature
  real(c_double), intent(out) :: values(4)
  integer(c_int), intent(out) :: count
  ! Locals
  real :: Cmat(nc,nc)
  integer :: i,j
  call binaryThirdVirialCoeffMatrix(temperature,Cmat)
  count = nc*nc
  do j=1,nc
    do i=1,nc
      values((j-1)*nc+i) = Cmat(i,j)
    end do
  end do
end subroutine query_binary_third_virial_matrix
