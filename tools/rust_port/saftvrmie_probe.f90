! Isolated exact-source executable. One model per process (Fortran globals).
program saftvrmie_probe
  use eoslibinit, only: init_saftvrmie, init_quantum_saftvrmie
  use eosTV, only: Fres
  use eos, only: zfac, enthalpy, entropy, thermo, specificVolume
  use saturation, only: safe_bubP
  use speed_of_sound, only: sound_velocity_2ph
  use thermopack_constants, only: LIQPH, VAPPH
  use thermopack_var, only: get_active_eos, base_eos_param
  use association_var, only: ELLIOT
  use saftvrmie_containers, only: svrm_opt, saftvrmie_param, calcBinaryMieParmeters, svrmie_var
  use saftvrmie_association, only: g_rdf_saftvrmie_ij_TVN
  use saft_interface, only: calc_saft_dispersion, calc_saft_hard_sphere, &
       calc_hard_sphere_diameter, calc_hard_sphere_diameter_ij, &
       sigma_eff_ij, epsilon_eff_ij, alpha, de_Broglie_wavelength, potential
  implicit none
  character(len=128) :: comps, param, variant
  integer :: nc, order, hs, controls, shift, kind, ph, ierr
  real :: t, v, cutoff, f, ft, fv, ftt, ftv, fvv
  real :: h, ftmp, fv_plus, fv_minus, dv_fv, dt_fv
  real :: z, dzdt, dzdp, hr, dhrdt, dhrdp, ht, dhtdt, dhtdp
  real :: sr, dsrdt, dsrdp, st, dstdt, dstdp, sos
  real :: dij, dijt, seff, eeff, lam_db, psat, vl, vv
  real :: r_arr(2), pot_arr(2), force_arr(2)
  real, allocatable :: n(:), x(:), y(:), fn(:), ftn(:), fvn(:), fnn(:,:)
  real, allocatable :: dzdn(:), dhrdn(:), dhtdn(:), dsrdn(:), dstdn(:)
  real, allocatable :: lnphi(:), dlnphidt(:), dlnphidp(:), dlnphidn(:,:)
  real, allocatable :: d_diag(:), dt_diag(:), dtt_diag(:), amat(:,:)
  class(base_eos_param), pointer :: act_eos
  read(*,'(A)') comps
  read(*,'(A)') param
  read(*,'(A)') variant
  read(*,*) nc, order, hs, controls, cutoff, shift, kind
  read(*,*) t, v
  allocate(n(nc), x(nc), y(nc), fn(nc), ftn(nc), fvn(nc), fnn(nc,nc))
  allocate(dzdn(nc), dhrdn(nc), dhtdn(nc), dsrdn(nc), dstdn(nc))
  allocate(lnphi(nc), dlnphidt(nc), dlnphidp(nc), dlnphidn(nc,nc))
  allocate(d_diag(nc), dt_diag(nc), dtt_diag(nc), amat(nc,nc))
  read(*,*) n
  if (order < 0) then
    call init_saftvrmie(trim(comps), trim(param))
  else if (trim(variant) == 'additive-exact-quantum') then
    call init_quantum_saftvrmie(trim(comps), order, additive_hs_ref=.true., parameter_reference=trim(param))
  else
    call init_quantum_saftvrmie(trim(comps), order, parameter_reference=trim(param))
  endif
  if (trim(variant) /= 'additive-exact-quantum') then
    call svrm_opt%set_hs_reference(hs)
  endif
  svrm_opt%enable_hs = btest(controls, 0)
  svrm_opt%enable_A1 = btest(controls, 1)
  svrm_opt%enable_A2 = btest(controls, 2)
  svrm_opt%enable_A3 = btest(controls, 3)
  svrm_opt%enable_chain = btest(controls, 4)
  select case (trim(variant))
  case ('jaen-kahn')
    svrm_opt%quantum_correction_spec = 1
    call calcBinaryMieParmeters(saftvrmie_param, 0)
  case ('geometric')
    svrm_opt%use_epsrule_Lafitte = .false.
    call calcBinaryMieParmeters(saftvrmie_param, 0)
  case ('arithmetic-cross')
    svrm_opt%exact_crosspot_eff = .false.
  case ('linear')
    svrm_opt%zeta_mixing_rule = 2
  case ('linear-then-nonadd')
    svrm_opt%zeta_mixing_rule = 2
    call svrm_opt%set_hs_reference(4)
  case ('bmcsl')
    svrm_opt%Khs_EoS = 1
  case ('kolafa')
    svrm_opt%pure_hs_EoS = 0
  case ('pair-a3')
    svrm_opt%a3_model = 2
  case ('classical-a2')
    svrm_opt%quantum_correct_A2 = .false.
  case ('no-extra')
    svrm_opt%enable_hs_extra = .false.
  case ('elliot')
    act_eos => get_active_eos()
    act_eos%assoc%delta_combrule = ELLIOT
  end select
  if (cutoff > 0) then
    call svrm_opt%set_r_cut(cutoff)
    svrm_opt%enable_shift_correction = (shift /= 0)
  endif
  if (kind == 4) then
    ! Normal TP path: v carries pressure [Pa], shift carries phase (1=LIQ, 2=VAP, 3=MINGIBBS).
    x = n / sum(n)
    ph = shift
    call zfac(t, v, x, ph, z, dzdt, dzdp, dzdn)
    call enthalpy(t, v, x, ph, hr, dhrdt, dhrdp, dhrdn, residual=.true.)
    call enthalpy(t, v, x, ph, ht, dhtdt, dhtdp, dhtdn, residual=.false.)
    call entropy(t, v, x, ph, sr, dsrdt, dsrdp, dsrdn, residual=.true.)
    call entropy(t, v, x, ph, st, dstdt, dstdp, dstdn, residual=.false.)
    call thermo(t, v, x, ph, lnphi, dlnphidt, dlnphidp, dlnphidn)
    sos = sound_velocity_2ph(t, v, x, x, x, 0.0, 1.0, ph)
    write(*,'(A,*(1X,ES25.16E3))') 'TP', z, dzdt, dzdp, dzdn, &
         hr, dhrdt, dhrdp, dhrdn, ht, dhtdt, dhtdp, dhtdn, &
         sr, dsrdt, dsrdp, dsrdn, st, dstdt, dstdp, dstdn, &
         lnphi, dlnphidt, dlnphidp, transpose(dlnphidn), sos
  else if (kind == 5) then
    ! Model-side diameters, effective parameters, alpha, and de Broglie wavelength.
    call calc_hard_sphere_diameter(t, d_diag, dt_diag, dtt_diag)
    call calc_hard_sphere_diameter_ij(1, nc, t, dij, dijt)
    call sigma_eff_ij(1, nc, t, seff)
    call epsilon_eff_ij(1, nc, t, eeff)
    call alpha(t, amat)
    if (order >= 0) then
      call de_Broglie_wavelength(1, t, lam_db)
      write(*,'(A,*(1X,ES25.16E3))') 'HOOKS', d_diag, dt_diag, dtt_diag, &
           dij, dijt, seff, eeff, transpose(amat), lam_db
    else
      write(*,'(A,*(1X,ES25.16E3))', advance='no') 'HOOKS', d_diag, dt_diag, dtt_diag, &
           dij, dijt, seff, eeff, transpose(amat)
      write(*,'(A)') ' UNSUPPORTED_SOURCE'
    endif
  else if (kind == 6) then
    ! Pure-fluid saturation pressure and coexisting liquid/vapor molar volumes.
    x = n / sum(n)
    psat = safe_bubP(t, x, y, ierr)
    if (ierr /= 0) stop 2
    call specificVolume(t, psat, x, LIQPH, vl)
    call specificVolume(t, psat, x, VAPPH, vv)
    write(*,'(A,*(1X,ES25.16E3))') 'SAT', psat, vl, vv
  else if (kind == 7) then
    ! Pair potential u(r)/kB [K] (with default 500*eps cap) and uncapped force -du/dr/kB [K/m].
    ! v carries r / sigma_ij(1, nc). Pass a 2-point array so Fortran truncation handles r > r_cut.
    r_arr(1) = v * saftvrmie_param%sigma_ij(1, nc)
    r_arr(2) = 100.0 * saftvrmie_param%sigma_ij(1, nc)
    call potential(1, nc, 2, r_arr, t, pot_arr, force=force_arr)
    write(*,'(A,*(1X,ES25.16E3))') 'POT', pot_arr(1), force_arr(1)
  else
    if (kind == 0) then
      call Fres(t, v, n, f, ft, fv, fn, ftt, ftv, fvv, ftn, fvn, fnn)
    else if (kind == 2) then
      call g_rdf_saftvrmie_ij_TVN(t, v, n, 1, nc, svrmie_var, f, fv, fvv, ft, ftt, ftv, fn, fvn, ftn, fnn)
    else if (kind == 3) then
      call calc_saft_hard_sphere(t, v, n, f, ft, fv, fn, ftt, ftv, fvv, ftn, fvn, fnn)
    else
      call calc_saft_dispersion(t, v, n, f, ft, fv, fn, ftt, ftv, fvv, ftn, fvn, fnn)
    endif
    dv_fv = fvv
    dt_fv = ftv
    if ((kind == 0 .or. kind == 3) .and. svrm_opt%enable_hs_extra .and. nc > 1) then
      h = v * 1e-5
      if (kind == 0) then
        call Fres(t, v + h, n, F_V=fv_plus)
        call Fres(t, v - h, n, F_V=fv_minus)
      else
        call calc_saft_hard_sphere(t, v + h, n, a=ftmp, a_V=fv_plus)
        call calc_saft_hard_sphere(t, v - h, n, a=ftmp, a_V=fv_minus)
      endif
      dv_fv = (fv_plus - fv_minus) / (2 * h)
      if (kind == 3) then
        h = t * 1e-5
        call calc_saft_hard_sphere(t + h, v, n, a=ftmp, a_V=fv_plus)
        call calc_saft_hard_sphere(t - h, v, n, a=ftmp, a_V=fv_minus)
        dt_fv = (fv_plus - fv_minus) / (2 * h)
      endif
    endif
    write(*,'(A,*(1X,ES25.16E3))') 'TVN', f, ft, fv, fn, ftt, ftv, ftn, fvv, fvn, transpose(fnn), dv_fv, dt_fv
  endif
end program
