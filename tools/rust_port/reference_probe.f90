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
