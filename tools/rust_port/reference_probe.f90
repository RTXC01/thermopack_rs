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
