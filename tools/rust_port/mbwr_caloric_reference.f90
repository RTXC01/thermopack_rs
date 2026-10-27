! Total TV calls through the unchanged original dispatch. Only coefficient
! allocation is replaced to avoid unrelated fatal critical-grid initialization.
program mbwr_caloric_reference
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use eoslibinit, only: init_multiparameter
  use thermopack_var, only: get_active_thermo_model, get_active_eos, thermo_model, Rgas
  use eos_parameters, only: single_eos
  use eosdata, only: meosMbwr19, meosMbwr32
  use mbwr, only: eosmbwr,allocNIJL,fillIJL19,fillIJL32,computeZCoeff,computeHelmCoeff
  use mbwrdata, only: mbwr19Array,mbwr32Array
  use ideal, only: set_enthalpy_of_formation,set_standard_entropy,Cpideal
  use eostv, only: Fres,Fideal,pressure,internal_energy_tv,entropy_tv,enthalpy_tv,free_energy_tv,chemical_potential_tv
  use speed_of_sound, only: speed_of_sound_TV
  implicit none
  type(thermo_model), pointer :: model
  integer :: cases,k,record,profile,ierr
  real :: t,v,n(1),fr(10),fi(10),f(10),p(4),u(4),h(4),s(4),a(4),g(4),mu(4)
  real :: vals(45),cv,cp,phi,sound,rawsound,cpflag,phiflag,soundflag
  character(len=24) :: comp,mode
  call get_command_argument(1,mode)
  call init_multiparameter('C1','MBWR32','DEFAULT')
  model => get_active_thermo_model()
  read(*,*) cases
  do k=1,cases
    read(*,*) record,t,v,n,comp,profile
    if (trim(mode)=='public') then
      write(*,'(A)') 'PUBLIC_BEGIN'
      flush(6)
      if (record<=11) then
        call init_multiparameter(trim(comp),'MBWR32','DEFAULT')
      else
        call init_multiparameter(trim(comp),'MBWR19','DEFAULT')
      endif
      model => get_active_thermo_model()
      write(*,'(A)') 'PUBLIC_RETURNED'
    else
    ! A real constructor allocates a fresh gendata instance. Reusing a
    ! previously referenced instance would retain its integration constants:
    ! init_from_name copies gendatadb fields without resetting these fields.
    deallocate(model%comps(1)%p_comp)
    allocate(model%comps(1)%p_comp)
    call model%comps(1)%p_comp%init_from_name(trim(comp),'DEFAULT',ierr)
    if (ierr/=0) error stop 'component metadata initialization failed'
    select type(eos => model%eos(1)%p_eos)
    type is(single_eos)
      call coefficients(record,eos%mbwr_meos(1))
      eos%subeosidx=merge(meosMbwr32,meosMbwr19,record<=11)
    class default
      error stop 'expected exact single_eos dynamic type'
    end select
    endif
    if (profile==1 .or. profile==3) then
      call set_enthalpy_of_formation(1,-12345.678d0)
      call set_standard_entropy(1,210.25d0,'1ATM')
    endif
    if (profile==2 .or. profile==3) then
      ! Same fields assigned by the original ideal-Cp setter. Its constants
      ! intentionally remain unchanged until a reference setter is called.
      model%comps(1)%p_comp%id_cp%cptype=8
      model%comps(1)%p_comp%id_cp%cp=0
      model%comps(1)%p_comp%id_cp%cp(1)=3.5d0
    endif
    call jets(t,v,n,fr,fi)
    f=fr+fi
    p(1)=pressure(t,v,n,dpdt=p(2),dpdv=p(3),dpdn=p(4:4))
    call internal_energy_tv(t,v,n,u(1),u(2),u(3),u(4:4))
    call enthalpy_tv(t,v,n,h(1),h(2),h(3),h(4:4))
    call entropy_tv(t,v,n,s(1),s(2),s(3),s(4:4))
    call free_energy_tv(t,v,n,a(1),a(2),a(3),a(4:4))
    call chemical(t,v,n,mu)
    g=a+v*p
    g(3)=a(3)+p(1)+v*p(3)
    cv=u(2)
    cpflag=0; cp=0
    if (abs(p(3))>64*epsilon(t)*Rgas*t*n(1)/v**2) then
      cp=cv-t*p(2)**2/p(3)
      if (ieee_is_finite(cp)) cpflag=1
    endif
    if (cpflag==0) cp=0
    phiflag=0; phi=0
    if (p(1)>0) then
      phi=fr(4)-log(p(1)*v/(n(1)*Rgas*t))
      if (ieee_is_finite(phi)) phiflag=1
    endif
    if (phiflag==0) phi=0
    rawsound=speed_of_sound_TV(t,v,n)
    soundflag=0; sound=0
    if (ieee_is_finite(rawsound) .and. rawsound>0) then
      soundflag=1; sound=rawsound
    endif
    vals=[f,p,u,h,s,a,g,mu,cv,cpflag,cp,phiflag,phi,soundflag,sound]
    write(*,'(A,I0,A,*(ES26.17E3,1X))') 'CASE=',k,' VALUES=',vals
    write(*,'(A,I0,A,*(ES26.17E3,1X))') 'CASE=',k,' IDEAL=', &
         model%comps(1)%p_comp%href_int,model%comps(1)%p_comp%sref_int, &
         model%comps(1)%p_comp%mw,real(model%comps(1)%p_comp%id_cp%cptype), &
         Cpideal(model%comps(1)%p_comp,1,t),fi
    write(*,'(A,I0,A,ES26.17E3)') 'CASE=',k,' RAW_SOUND=',rawsound
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
