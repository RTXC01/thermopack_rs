! Isolated exact-source probe for src/lj_splined.f90 and src/hardsphere_wca.f90.
! One model per process (Fortran module globals).
program ljs_probe
  use eoslibinit, only: init_ljs, init_lj
  use eosTV, only: Fres
  use eos, only: zfac, enthalpy, entropy, thermo, specificVolume
  use speed_of_sound, only: sound_velocity_2ph
  use thermopack_var, only: get_active_eos, base_eos_param
  use saftvrmie_containers, only: svrm_opt, saftvrmie_dhs, saftvrmie_zeta, &
       allocate_saftvrmie_dhs, allocate_saftvrmie_zeta, &
       cleanup_saftvrmie_dhs, cleanup_saftvrmie_zeta
  use lj_splined, only: ljs_bh_eos, ljs_wca_eos, ljx_ux_eos, ljx_ux_eos_constructor, &
       ljs_bh_model_control, ljs_bh_set_pure_params, ljs_wca_model_control, &
       ljs_uv_model_control, ljs_wca_set_pure_params, &
       calc_ljs_dispersion, calc_ljs_hard_sphere, calc_wca_soft_repulsion, &
       calc_uf_wca, calc_uf_wca_tvn, calc_uv_wca_tvn, calc_LJs_a1, calc_LJs_a2
  use hardsphere_wca, only: calc_dhs_WCA, calc_cavity_integral_LJ_Fres_test_TVN, &
       calc_cavity_integral_LJ_Fres_test, calcZetaX_vdW_no_segments
  use saft_interface, only: calc_hard_sphere_diameter, alpha, potential, &
       estimate_critical_parameters, test_fmt_compatibility
  implicit none

  character(len=64) :: model, variant
  integer :: kind, phase_flag
  real :: t, v_or_p, n1, sigma, eps
  real :: n(1), x(1), f, ft, fv, fn(1), ftt, ftv, fvv, ftn(1), fvn(1), fnn(1,1)
  real :: vol, z, dzdt, dzdp, dzdn(1), hr, dhrdt, dhrdp, dhrdn(1)
  real :: ht, dhtdt, dhtdp, dhtdn(1), sr, dsrdt, dsrdp, dsrdn(1)
  real :: st, dstdt, dstdp, dstdn(1), lnphi(1), dlnphidt(1), dlnphidp(1), dlnphidn(1,1), sos
  real :: d_arr(1), dt_arr(1), dtt_arr(1), amat(1,1), tc_est, vc_est
  real :: r_arr(2), pot_arr(2), force_arr(2)
  real :: uf_a, uf_ar, uf_at, uf_arr, uf_att, uf_art
  real :: fe, fee, fte
  real :: n2(2)
  logical :: is_fmt, na_en, is_lj_flag
  class(base_eos_param), pointer :: act_eos
  class(ljs_wca_eos), pointer :: wca_ptr
  type(ljx_ux_eos) :: ctor_ljs_uf, ctor_lj_uf, uf_local
  type(saftvrmie_dhs) :: dhs1, dhs2
  type(saftvrmie_zeta) :: zeta2

  read(*,'(A)') model
  read(*,'(A)') variant
  read(*,*) kind, phase_flag
  read(*,*) t, v_or_p, n1, sigma, eps

  n(1) = n1
  x(1) = 1.0

  ! Standalone non-EOS probes (kinds 6, 7, 8)
  if (trim(model) == 'CAVITY') then
    call calc_cavity_integral_LJ_Fres_test_TVN(1, sigma, eps, t, v_or_p, n, &
         f, ft, fv, fn, ftt, fvv, ftv, ftn, fvn, fnn)
    write(*,'(A,*(1X,ES25.16E3))') 'TVN', f, ft, fv, fn(1), ftt, ftv, ftn(1), fvv, fvn(1), fnn(1,1)
    stop
  end if

  if (trim(model) == 'ZETAX2') then
    ! Binary packing fraction without segments: t=T, v_or_p=V, n1=n(1), sigma=sigma1, eps=eps1
    ! Second component uses fixed (sigma2=3.7e-10, eps2=150.0, n2=0.7) with arithmetic cross diameter.
    call allocate_saftvrmie_dhs(2, dhs2)
    call allocate_saftvrmie_zeta(2, zeta2)
    call allocate_saftvrmie_dhs(1, dhs1)
    call calc_dhs_WCA(1, sigma, eps, t, dhs1)
    dhs2%d(1,1) = dhs1%d(1,1); dhs2%d_T(1,1) = dhs1%d_T(1,1); dhs2%d_TT(1,1) = dhs1%d_TT(1,1)
    call calc_dhs_WCA(1, 3.7e-10, 150.0, t, dhs1)
    dhs2%d(2,2) = dhs1%d(1,1); dhs2%d_T(2,2) = dhs1%d_T(1,1); dhs2%d_TT(2,2) = dhs1%d_TT(1,1)
    dhs2%d(1,2) = 0.5*(dhs2%d(1,1) + dhs2%d(2,2)); dhs2%d(2,1) = dhs2%d(1,2)
    dhs2%d_T(1,2) = 0.5*(dhs2%d_T(1,1) + dhs2%d_T(2,2)); dhs2%d_T(2,1) = dhs2%d_T(1,2)
    dhs2%d_TT(1,2) = 0.5*(dhs2%d_TT(1,1) + dhs2%d_TT(2,2)); dhs2%d_TT(2,1) = dhs2%d_TT(1,2)
    n2 = (/ n1, 0.7 /)
    call calcZetaX_vdW_no_segments(2, t, v_or_p, n2, dhs2, zeta2)
    write(*,'(A,*(1X,ES25.16E3))') 'ZETAX2', zeta2%zx, zeta2%zx_T, zeta2%zx_TT, &
         zeta2%zx_V, zeta2%zx_VV, zeta2%zx_VVV, zeta2%zx_TV, zeta2%zx_VVT, zeta2%zx_VTT, &
         zeta2%zx_n(1), zeta2%zx_n(2), zeta2%zx_Tn(1), zeta2%zx_Tn(2), &
         zeta2%zx_Vn(1), zeta2%zx_Vn(2), zeta2%zx_VTn(1), zeta2%zx_VTn(2), &
         zeta2%zx_nn(1,1), zeta2%zx_nn(1,2), zeta2%zx_nn(2,1), zeta2%zx_nn(2,2), &
         zeta2%zx_Vnn(1,1), zeta2%zx_Vnn(1,2), zeta2%zx_Vnn(2,1), zeta2%zx_Vnn(2,2)
    call cleanup_saftvrmie_dhs(dhs1)
    call cleanup_saftvrmie_dhs(dhs2)
    call cleanup_saftvrmie_zeta(zeta2)
    stop
  end if

  if (trim(model) == 'UF-KERNEL') then
    is_lj_flag = (trim(variant) == 'lj')
    if (kind == 7) then
      call calc_uf_wca(v_or_p, t, is_lj_flag, uf_a, uf_ar, uf_at, uf_arr, uf_att, uf_art)
      write(*,'(A,*(1X,ES25.16E3))') 'UF6', uf_a, uf_ar, uf_at, uf_arr, uf_att, uf_art
    else
      uf_local%lj_potential = is_lj_flag
      uf_local%is_uf_theory = .true.
      call calc_uf_wca_tvn(uf_local, 1, sigma, eps, t, v_or_p, n, &
           f, ft, fv, fn, ftt, fvv, ftv, ftn, fvn, fnn)
      write(*,'(A,*(1X,ES25.16E3))') 'TVN', f, ft, fv, fn(1), ftt, ftv, ftn(1), fvv, fvn(1), fnn(1,1)
    end if
    stop
  end if

  if (trim(model) == 'DEFECTS') then
    if (trim(variant) == 'uf-ctor') then
      ctor_ljs_uf = ljx_ux_eos_constructor('LJS-UF')
      ctor_lj_uf = ljx_ux_eos_constructor('LJ-UF')
      write(*,'(A,4(1X,I2))') 'DEFECT_UF_CTOR', &
           merge(1, 0, ctor_ljs_uf%lj_potential), merge(1, 0, ctor_ljs_uf%is_uf_theory), &
           merge(1, 0, ctor_lj_uf%lj_potential), merge(1, 0, ctor_lj_uf%is_uf_theory)
    else if (trim(variant) == 'bh-set-pure') then
      call init_ljs('BH')
      call ljs_bh_set_pure_params(sigma, eps)
      act_eos => get_active_eos()
      select type (p_bh => act_eos)
      class is (ljs_bh_eos)
        write(*,'(A,*(1X,ES25.16E3))') 'DEFECT_BH_SETPARAM', &
             p_bh%saftvrmie_param%sigma_ij(1,1), p_bh%saftvrmie_var%sigma_eff%d(1,1), &
             p_bh%saftvrmie_param%sigma_ij_cube(1,1), p_bh%saftvrmie_param%eps_divk_ij(1,1)
      end select
    else if (trim(variant) == 'cavity-test') then
      call calc_cavity_integral_LJ_Fres_test(1, sigma, eps, t, v_or_p, f, ft, fe, ftt, fee, fte)
      write(*,'(A,*(1X,ES25.16E3))') 'DEFECT_CAV_TEST', f, ft, fe, ftt, fee, fte
    else if (trim(variant) == 'ljs-uf-init') then
      call init_ljs('UF')
    end if
    stop
  end if

  ! Initialize EOS model via normal public ThermoPack entry point
  select case (trim(model))
  case ('LJS-BH')
    call init_ljs('BH')
  case ('LJS-WCA')
    call init_ljs('WCA')
  case ('LJS-UV')
    call init_ljs('UV')
  case ('LJS-UF')
    ! Upstream init_ljs('UF') aborts in assocIndices_bookkeeping (missing eosLJS_UF)
    ! and ljx_ux_eos_constructor reverses INDEX arguments. Initialize ljx_ux_eos via
    ! init_ljs('UV') (identical Argon/ideal-gas setup) and select UF + LJS flags.
    call init_ljs('UV')
    act_eos => get_active_eos()
    select type (p_ux => act_eos)
    class is (ljx_ux_eos)
      p_ux%is_uf_theory = .true.
      p_ux%lj_potential = .false.
    end select
  case ('LJ-UF')
    call init_lj('UF')
    act_eos => get_active_eos()
    select type (p_ux => act_eos)
    class is (ljx_ux_eos)
      p_ux%is_uf_theory = .true.
      p_ux%lj_potential = .true.
    end select
  case default
    stop 10
  end select

  act_eos => get_active_eos()

  ! Apply custom parameters when requested (sigma > 0 and eps > 0)
  if (sigma > 0.0 .and. eps > 0.0) then
    select type (p_eos => act_eos)
    class is (ljs_bh_eos)
      call p_eos%set_sigma_eps(sigma, eps)
    class is (ljs_wca_eos)
      call ljs_wca_set_pure_params(sigma, eps)
    end select
  end if

  ! Apply model-specific variant controls
  select case (trim(variant))
  case ('default')
    continue
  case ('lafitte')
    call ljs_bh_model_control(.true., .true., .true., .true., .true., .true.)
  case ('no-chi')
    call ljs_bh_model_control(.false., .false., .true., .true., .true., .true.)
  case ('bh-a1-only')
    call ljs_bh_model_control(.false., .true., .false., .true., .false., .false.)
  case ('bh-a2-only')
    call ljs_bh_model_control(.false., .true., .false., .false., .true., .false.)
  case ('bh-a3-only')
    call ljs_bh_model_control(.false., .true., .false., .false., .false., .true.)
  case ('bh-hs-only')
    call ljs_bh_model_control(.false., .true., .true., .false., .false., .false.)
  case ('no-cavity')
    call ljs_wca_model_control(.false., .true., .true., .true., .true., .true.)
  case ('no-hs')
    call ljs_wca_model_control(.true., .false., .true., .true., .true., .true.)
  case ('wca-a1-only')
    call ljs_wca_model_control(.false., .false., .true., .false., .false., .false.)
  case ('wca-a2-only')
    call ljs_wca_model_control(.false., .false., .false., .true., .false., .false.)
  case ('wca-a3-only')
    call ljs_wca_model_control(.false., .false., .false., .false., .true., .false.)
  case ('wca-a4-only')
    call ljs_wca_model_control(.false., .false., .false., .false., .false., .true.)
  case ('temp-u')
    call ljs_uv_model_control(.true.)
  case ('no-virial')
    select type (p_ux => act_eos)
    class is (ljx_ux_eos)
      p_ux%enable_virial_term = .false.
    end select
  case ('kolafa')
    if (trim(model) == 'LJS-BH') then
      call svrm_opt%set_hs_reference(4)
    end if
    svrm_opt%pure_hs_EoS = 0
  case default
    stop 11
  end select

  if (kind == 0) then
    call Fres(t, v_or_p, n, f, ft, fv, fn, ftt, ftv, fvv, ftn, fvn, fnn)
    write(*,'(A,*(1X,ES25.16E3))') 'TVN', f, ft, fv, fn(1), ftt, ftv, ftn(1), fvv, fvn(1), fnn(1,1)
  else if (kind == 1) then
    call calc_ljs_dispersion(act_eos, 1, t, v_or_p, n, f, ft, fv, fn, ftt, fvv, ftv, ftn, fvn, fnn)
    write(*,'(A,*(1X,ES25.16E3))') 'TVN', f, ft, fv, fn(1), ftt, ftv, ftn(1), fvv, fvn(1), fnn(1,1)
  else if (kind == 2) then
    call calc_ljs_hard_sphere(act_eos, 1, t, v_or_p, n, f, ft, fv, fn, ftt, fvv, ftv, ftn, fvn, fnn)
    write(*,'(A,*(1X,ES25.16E3))') 'TVN', f, ft, fv, fn(1), ftt, ftv, ftn(1), fvv, fvn(1), fnn(1,1)
  else if (kind == 3) then
    select type (p_wca => act_eos)
    class is (ljs_wca_eos)
      wca_ptr => p_wca
      call calc_wca_soft_repulsion(wca_ptr, 1, t, v_or_p, n, f, ft, fv, fn, ftt, fvv, ftv, ftn, fvn, fnn)
    class default
      stop 12
    end select
    write(*,'(A,*(1X,ES25.16E3))') 'TVN', f, ft, fv, fn(1), ftt, ftv, ftn(1), fvv, fvn(1), fnn(1,1)
  else if (kind == 4) then
    ! Normal TP state: v_or_p = P [Pa], phase_flag = 1 (LIQPH) or 2 (VAPPH)
    call specificVolume(t, v_or_p, x, phase_flag, vol)
    call zfac(t, v_or_p, x, phase_flag, z, dzdt, dzdp, dzdn)
    call enthalpy(t, v_or_p, x, phase_flag, hr, dhrdt, dhrdp, dhrdn, residual=.true.)
    call enthalpy(t, v_or_p, x, phase_flag, ht, dhtdt, dhtdp, dhtdn, residual=.false.)
    call entropy(t, v_or_p, x, phase_flag, sr, dsrdt, dsrdp, dsrdn, residual=.true.)
    call entropy(t, v_or_p, x, phase_flag, st, dstdt, dstdp, dstdn, residual=.false.)
    call thermo(t, v_or_p, x, phase_flag, lnphi, dlnphidt, dlnphidp, dlnphidn)
    sos = sound_velocity_2ph(t, v_or_p, x, x, x, 0.0, 1.0, phase_flag)
    write(*,'(A,*(1X,ES25.16E3))') 'TP', vol, z, dzdt, dzdp, &
         hr, dhrdt, dhrdp, ht, dhtdt, dhtdp, &
         sr, dsrdt, dsrdp, st, dstdt, dstdp, &
         lnphi(1), dlnphidt(1), dlnphidp(1), sos
  else if (kind == 5) then
    ! Model-side hooks: diameter, alpha, critical estimate, FMT compatibility, potential(r)
    call calc_hard_sphere_diameter(t, d_arr, dt_arr, dtt_arr)
    call alpha(t, amat)
    call estimate_critical_parameters(1, tc_est, vc_est)
    call test_fmt_compatibility(is_fmt, na_en)
    select type (p_eos => act_eos)
    class is (ljs_bh_eos)
      r_arr(1) = v_or_p * p_eos%saftvrmie_param%sigma_ij(1,1)
      r_arr(2) = 10.0 * p_eos%saftvrmie_param%sigma_ij(1,1)
    class is (ljs_wca_eos)
      r_arr(1) = v_or_p * p_eos%sigma
      r_arr(2) = 10.0 * p_eos%sigma
    end select
    call potential(1, 1, 2, r_arr, t, pot_arr, force=force_arr)
    write(*,'(A,*(1X,ES25.16E3))') 'HOOKS', d_arr(1), dt_arr(1), dtt_arr(1), &
         amat(1,1), tc_est, vc_est, &
         merge(1.0, 0.0, is_fmt), merge(1.0, 0.0, na_en), pot_arr(1), force_arr(1)
  else
    stop 13
  end if
end program ljs_probe
