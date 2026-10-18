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
