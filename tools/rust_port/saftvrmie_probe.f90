! Isolated exact-source executable. One model per process (Fortran globals).
program saftvrmie_probe
  use eoslibinit, only: init_saftvrmie, init_quantum_saftvrmie
  use eosTV, only: Fres
  use saftvrmie_containers, only: svrm_opt
  use saft_interface, only: calc_saft_dispersion
  implicit none
  character(len=128) :: comps, param
  integer :: nc, order, hs, controls, shift, kind
  real :: t,v,cutoff,f,ft,fv,ftt,ftv,fvv
  real :: h, fv_plus, fv_minus, dv_fv
  real, allocatable :: n(:), fn(:), ftn(:), fvn(:), fnn(:,:)
  read(*,'(A)') comps
  read(*,'(A)') param
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
  if (cutoff>0) then
    call svrm_opt%set_r_cut(cutoff)
    svrm_opt%enable_shift_correction=shift/=0
  endif
  if (kind==0) then
    call Fres(t,v,n,f,ft,fv,fn,ftt,ftv,fvv,ftn,fvn,fnn)
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
