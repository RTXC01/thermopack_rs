! Isolated exact-source executable. One model per process (Fortran globals).
program saftvrmie_probe
  use eoslibinit, only: init_saftvrmie, init_quantum_saftvrmie
  use eosTV, only: Fres
  use saftvrmie_containers, only: svrm_opt, saftvrmie_param, calcBinaryMieParmeters, svrmie_var
  use saftvrmie_association, only: g_rdf_saftvrmie_ij_TVN
  use saft_interface, only: calc_saft_dispersion
  implicit none
  character(len=128) :: comps, param, variant
  integer :: nc, order, hs, controls, shift, kind
  real :: t,v,cutoff,f,ft,fv,ftt,ftv,fvv
  real :: h, fv_plus, fv_minus, dv_fv
  real, allocatable :: n(:), fn(:), ftn(:), fvn(:), fnn(:,:)
  read(*,'(A)') comps
  read(*,'(A)') param
  read(*,'(A)') variant
  read(*,*) nc,order,hs,controls,cutoff,shift,kind
  read(*,*) t,v
  allocate(n(nc),fn(nc),ftn(nc),fvn(nc),fnn(nc,nc))
  read(*,*) n
  if (order < 0) then
    call init_saftvrmie(trim(comps), trim(param))
  else
    call init_quantum_saftvrmie(trim(comps), order, parameter_reference=trim(param))
  endif
  call svrm_opt%set_hs_reference(hs)
  svrm_opt%enable_hs=btest(controls,0)
  svrm_opt%enable_A1=btest(controls,1)
  svrm_opt%enable_A2=btest(controls,2)
  svrm_opt%enable_A3=btest(controls,3)
  svrm_opt%enable_chain=btest(controls,4)
  select case (trim(variant))
  case ('jaen-kahn')
    svrm_opt%quantum_correction_spec=1
    call calcBinaryMieParmeters(saftvrmie_param,0)
  case ('geometric')
    svrm_opt%use_epsrule_Lafitte=.false.
    call calcBinaryMieParmeters(saftvrmie_param,0)
  case ('arithmetic-cross')
    svrm_opt%exact_crosspot_eff=.false.
  case ('linear')
    svrm_opt%zeta_mixing_rule=2
  case ('bmcsl')
    svrm_opt%Khs_EoS=1
  case ('kolafa')
    svrm_opt%pure_hs_EoS=0
  case ('pair-a3')
    svrm_opt%a3_model=2
  case ('classical-a2')
    svrm_opt%quantum_correct_A2=.false.
  case ('no-extra')
    svrm_opt%enable_hs_extra=.false.
  end select
  if (cutoff>0) then
    call svrm_opt%set_r_cut(cutoff)
    svrm_opt%enable_shift_correction=shift/=0
  endif
  if (kind==0) then
    call Fres(t,v,n,f,ft,fv,fn,ftt,ftv,fvv,ftn,fvn,fnn)
  else if (kind==2) then
    call g_rdf_saftvrmie_ij_TVN(t,v,n,1,nc,svrmie_var,f,fv,fvv,ft,ftt,ftv,fn,fvn,ftn,fnn)
  else
    call calc_saft_dispersion(t,v,n,f,ft,fv,fn,ftt,ftv,fvv,ftn,fvn,fnn)
  endif
  dv_fv=fvv
  if (kind==0 .and. svrm_opt%enable_hs_extra .and. nc>1) then
    h=v*1e-5
    call Fres(t,v+h,n,F_V=fv_plus)
    call Fres(t,v-h,n,F_V=fv_minus)
    dv_fv=(fv_plus-fv_minus)/(2*h)
  endif
  write(*,'(A,*(1X,ES25.16E3))') 'TVN',f,ft,fv,fn,ftt,ftv,ftn,fvv,fvn,transpose(fnn),dv_fv
end program
