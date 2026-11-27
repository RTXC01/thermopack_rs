! Isolated exact-source Fortran probe for CPA second derivatives/mixing,
! UMR/VTPR and HV/WS cubic mixing, apparent_composition.f90, tuning.f90,
! and complexmodelinit.f90.
! One model per process (ThermoPack uses module globals).
program cpa_cubic_tuning_probe
  use eoslibinit, only: init_cpa, init_cubic, init_thermo
  use complexmodelinit, only: init_VTPR, init_UMR
  use eosTV, only: Fres
  use eos, only: zfac, enthalpy, entropy, thermo, specificVolume
  use speed_of_sound, only: sound_velocity_2ph
  use thermopack_var, only: get_active_eos, base_eos_param, add_eos
  use cubic_eos, only: cb_eos, cpa_eos, cbMixHVCPA, cbMixHVCPA2
  use cbmix, only: calcBmixCubic
  use wong_sandler, only: fidel_alpha
  use association_var, only: ELLIOT, association
  use AssocSchemeUtils, only: assocIndices_bookkeeping, assoc_scheme_2B, assoc_scheme_1ea
  use eosdata, only: eosCPA
  use saft_association, only: solve_for_X_k
  use unifac, only: unifdb, getUNIFACgroupInteraction
  use apparent_compostion, only: apparent_container, apparent_constructor
  implicit none

  character(len=64) :: mode, comps, eos_str, mix_str, ref_str, variant
  integer :: nc, phase_flag, i, j, k, ns, nce, ncsym
  real :: t, v_or_p
  real, allocatable :: n(:), x(:), fn(:), ftn(:), fvn(:), fnn(:,:)
  real, allocatable :: dzdn(:), dhrdn(:), dhtdn(:), dsrdn(:), dstdn(:)
  real, allocatable :: lnphi(:), dlnphidt(:), dlnphidp(:), dlnphidn(:,:)
  real, allocatable :: xk(:)
  real :: f, ft, fv, ftt, ftv, fvv
  real :: vol, z, dzdt, dzdp, hr, dhrdt, dhrdp, ht, dhtdt, dhtdp
  real :: sr, dsrdt, dsrdp, st, dstdt, dstdp, sos
  real :: p_vals(10), bmix, fa
  real :: g1, g2, g3, g4, g5, g6, g7, g8
  integer :: ci_type
  class(base_eos_param), pointer :: act_eos

  ! Apparent composition variables
  type(apparent_container) :: app
  type(association) :: tmp_assoc
  real, allocatable :: app_n(:), app_ne(:), fe_n(:), fe_tn(:), fe_vn(:), fe_nn(:,:)
  real, allocatable :: app_fn(:), app_ftn(:), app_fvn(:), app_fnn(:,:)
  real, allocatable :: lnf_r(:), dt_r(:), dp_r(:), dn_r(:,:)
  real, allocatable :: lnf_a(:), dt_a(:), dp_a(:), dn_a(:,:), mod_lnf(:)
  real :: sumne, v_val

  read(*,'(A)') mode
  read(*,'(A)') variant

  ! ---- 1. Defect probes ----
  if (trim(mode) == 'DEFECT') then
    if (trim(variant) == 'init_umr') then
      k = add_eos()
      call init_UMR(2, 'C1,C3', 2)
    else if (trim(variant) == 'cpa_hv1_init') then
      call init_cpa('H2O,C1', 'SRK', mixing='HV1')
    else if (trim(variant) == 'cpa_hv2_init') then
      call init_cpa('H2O,C1', 'SRK', mixing='HV2')
    else if (trim(variant) == 'cpa_setkijandji') then
      call init_cpa('H2O,CO2', 'SRK')
      call thermopack_setkijandji(1, 2, 0.05)
    else if (trim(variant) == 'cpa_setlijandji') then
      call init_cpa('H2O,CO2', 'SRK')
      call thermopack_setlijandji(1, 2, 0.02)
    else if (trim(variant) == 'cpa_setkij_single') then
      call init_cpa('H2O,CO2', 'SRK')
      call thermopack_setkij(1, 2, 0.075)
      call thermopack_getkij(1, 2, g1)
      call thermopack_getkij(2, 1, g2)
      write(*,'(A,2(1X,ES25.16E3))') 'DEFECT_CPA_SETKIJ', g1, g2
    else if (trim(variant) == 'vtpr_bmix_defect') then
      k = add_eos()
      call init_VTPR(2, 'C1,C3', 2)
      act_eos => get_active_eos()
      allocate(n(2), fn(2), ftn(2), fvn(2), fnn(2,2))
      n = (/ 0.4, 0.6 /)
      select type (p_cb => act_eos)
      class is (cb_eos)
        bmix = calcBmixCubic(2, p_cb, n) * 1.0e-3
        call Fres(310.0, 2.0e-4, n, f, ft, fv, fn, ftt, ftv, fvv, ftn, fvn, fnn)
        write(*,'(A,3(1X,ES25.16E3))') 'DEFECT_VTPR_BMIX', bmix, p_cb%b * 1.0e-3, f
      end select
    else if (trim(variant) == 'apparent_nan') then
      ! H2O + NaCl + CaCl2 with n_CaCl2 = 0 -> xe(Ca2+) = 0 -> NaN in lnfug(NaCl)
      app = apparent_constructor(3, 5)
      app%ncsym = 1
      call app%set_v_stoich_ij(2, 2, 1.0) ! NaCl -> Na+
      call app%set_v_stoich_ij(2, 3, 1.0) ! NaCl -> Cl-
      call app%set_v_stoich_ij(3, 4, 1.0) ! CaCl2 -> Ca2+
      call app%set_v_stoich_ij(3, 3, 2.0) ! CaCl2 -> 2 Cl-
      call app%update_v_sum()
      allocate(app_n(3), app_ne(5), lnf_r(5), lnf_a(3))
      app_n = (/ 0.9, 0.1, 0.0 /)
      call app%apparent_to_real_mole_numbers(app_n, app_ne)
      lnf_r = (/ -0.2, -0.5, -0.6, -0.8, 0.0 /)
      call app%TP_lnfug_apparent(3, app_ne, app_n, 1.0e5, lnf_r, lnf_a)
      write(*,'(A,3(1X,L2))') 'DEFECT_APPARENT_NAN', &
           (lnf_a(1) /= lnf_a(1)), (lnf_a(2) /= lnf_a(2)), (lnf_a(3) /= lnf_a(3))
    else if (trim(variant) == 'scheme_1ea') then
      call assocIndices_bookkeeping(tmp_assoc, 2, eosCPA, (/ assoc_scheme_2B, assoc_scheme_1ea /))
      write(*,'(A,4(1X,I0))') 'DEFECT_1EA', tmp_assoc%numAssocComps, &
           tmp_assoc%numAssocSites, tmp_assoc%comp_vs_sites(1,1), &
           tmp_assoc%comp_vs_sites(1,2)
    else
      stop 20
    end if
    stop
  end if

  ! ---- 2. Standalone fidel_alpha ----
  if (trim(mode) == 'FIDEL') then
    read(*,*) g1, g2
    fa = fidel_alpha(g1, g2)
    write(*,'(A,1X,ES25.16E3)') 'FIDEL', fa
    stop
  end if

  ! ---- 3. Standalone apparent_composition.f90 ----
  if (trim(mode) == 'APPARENT') then
    read(*,*) nc, nce, ncsym
    app = apparent_constructor(nc, nce)
    app%ncsym = ncsym
    read(*,*) k
    do i = 1, k
      read(*,*) phase_flag, j, v_val
      call app%set_v_stoich_ij(phase_flag, j, v_val)
    end do
    call app%update_v_sum()
    allocate(app_n(nc), app_ne(nce))
    allocate(fe_n(nce), fe_tn(nce), fe_vn(nce), fe_nn(nce,nce))
    allocate(app_fn(nc), app_ftn(nc), app_fvn(nc), app_fnn(nc,nc))
    allocate(lnf_r(nce), dt_r(nce), dp_r(nce), dn_r(nce,nce))
    allocate(lnf_a(nc), dt_a(nc), dp_a(nc), dn_a(nc,nc), mod_lnf(nc))
    read(*,*) t, v_or_p
    read(*,*) app_n
    read(*,*) fe_n
    read(*,*) fe_tn
    read(*,*) fe_vn
    do i = 1, nce
      read(*,*) fe_nn(i,:)
    end do
    read(*,*) lnf_r
    read(*,*) dt_r
    read(*,*) dp_r
    do i = 1, nce
      read(*,*) dn_r(i,:)
    end do
    call app%apparent_to_real_mole_numbers(app_n, app_ne)
    call app%real_to_apparent_differentials(fe_n, fe_tn, fe_vn, fe_nn, &
         app_fn, app_ftn, app_fvn, app_fnn)
    call app%TP_lnfug_apparent(nc, app_ne, app_n, v_or_p, lnf_r, lnf_a, &
         dt_r, dp_r, dn_r, dt_a, dp_a, dn_a)
    mod_lnf = lnf_a
    call app%getModFugacity(t, v_or_p, app_n / sum(app_n), mod_lnf, sumne)
    write(*,'(A,*(1X,ES25.16E3))') 'APPARENT', &
         app_ne, app_fn, app_ftn, app_fvn, &
         ((app_fnn(i,j), j=1,nc), i=1,nc), &
         lnf_a, dt_a, dp_a, &
         ((dn_a(i,j), j=1,nc), i=1,nc), &
         mod_lnf, sumne
    stop
  end if

  ! ---- 4. EOS-based modes: CPA_TVN, CPA_TP, CUBIC_TVN, CUBIC_TP, BMIX, TUNING_* ----
  read(*,'(A)') comps
  read(*,'(A)') eos_str
  read(*,'(A)') mix_str
  read(*,'(A)') ref_str
  read(*,*) nc, phase_flag
  read(*,*) t, v_or_p
  allocate(n(nc), x(nc), fn(nc), ftn(nc), fvn(nc), fnn(nc,nc))
  allocate(dzdn(nc), dhrdn(nc), dhtdn(nc), dsrdn(nc), dstdn(nc))
  allocate(lnphi(nc), dlnphidt(nc), dlnphidp(nc), dlnphidn(nc,nc))
  read(*,*) n
  x = n / sum(n)

  if (trim(mode) == 'CPA_TVN' .or. trim(mode) == 'CPA_TP') then
    if (trim(variant) == 'hv1_tuned' .or. trim(variant) == 'hv2_tuned') then
      ! Upstream cbselect.f90 omits mixGE%mGE assignment on HVCPA/HVCPA2; initialize
      ! CPA-SRK Classic and configure mixGE cleanly to validate the HV-CPA kernel.
      call init_cpa(trim(comps), trim(eos_str), mixing='Classic', parameter_reference=trim(ref_str))
      act_eos => get_active_eos()
      select type (p_cpa => act_eos)
      class is (cpa_eos)
        call p_cpa%mixGE%excess_gibbs_allocate_and_init(nc)
        if (trim(variant) == 'hv1_tuned') then
          p_cpa%mruleidx = cbMixHVCPA
          p_cpa%mixGE%mGE = cbMixHVCPA
        else
          p_cpa%mruleidx = cbMixHVCPA2
          p_cpa%mixGE%mGE = cbMixHVCPA2
        end if
      end select
      read(*,*) p_vals(1:8)
      call thermopack_setHVparam(1, 2, p_vals(1), p_vals(2), p_vals(3), p_vals(4), &
           p_vals(5), p_vals(6), p_vals(7), p_vals(8))
    else
      call init_cpa(trim(comps), trim(eos_str), mixing=trim(mix_str), parameter_reference=trim(ref_str))
      act_eos => get_active_eos()
      if (trim(variant) == 'orig_rdf') then
        select type (p_cpa => act_eos)
        class is (cpa_eos)
          p_cpa%useSimplifiedCPA = .false.
        end select
      else if (trim(variant) == 'elliot') then
        act_eos%assoc%delta_combrule = ELLIOT
      else if (trim(variant) == 'ws_tuned') then
        read(*,*) p_vals(1:6)
        call thermopack_setWSparam(1, 2, p_vals(1), p_vals(2), p_vals(3), p_vals(4), &
             p_vals(5), p_vals(6))
      end if
    end if
  else if (trim(mode) == 'BMIX') then
    if (trim(variant) == 'init_vtpr') then
      k = add_eos()
      call init_VTPR(nc, trim(comps), 2)
    else
      call init_thermo(trim(eos_str), trim(mix_str), 'Classic', trim(comps), 2, b_exponent=v_or_p)
    end if
    act_eos => get_active_eos()
    select type (p_cb => act_eos)
    class is (cb_eos)
      bmix = calcBmixCubic(nc, p_cb, n) * 1.0e-3
      write(*,'(A,*(1X,ES25.16E3))') 'BMIX', bmix, &
           ((p_cb%lowcase_bij(i,j) * 1.0e-3, j=1,nc), i=1,nc)
    end select
    stop
  else
    ! Cubic / tuning modes
    if (trim(variant) == 'init_vtpr') then
      k = add_eos()
      call init_VTPR(nc, trim(comps), 2)
    else
      call init_cubic(trim(comps), trim(eos_str), mixing=trim(mix_str), &
           parameter_reference=trim(ref_str))
    end if

    if (trim(variant) == 'hv_ws_tuned') then
      read(*,*) p_vals(1:10) ! alpha_ij, alpha_ji, a_ij, a_ji, b_ij, b_ji, c_ij, c_ji, kij, lij
      call thermopack_setHVparam(1, 2, p_vals(1), p_vals(2), p_vals(3), p_vals(4), &
           p_vals(5), p_vals(6), p_vals(7), p_vals(8))
      call thermopack_setkijandji(1, 2, p_vals(9))
      call thermopack_setlijandji(1, 2, p_vals(10))
    else if (trim(variant) == 'tune_hv') then
      read(*,*) p_vals(1:8)
      call thermopack_setHVparam(1, 2, p_vals(1), p_vals(2), p_vals(3), p_vals(4), &
           p_vals(5), p_vals(6), p_vals(7), p_vals(8))
      call thermopack_getHVparam(1, 2, g1, g2, g3, g4, g5, g6, g7, g8)
      call Fres(t, v_or_p, n, f, ft, fv, fn, ftt, ftv, fvv, ftn, fvn, fnn)
      write(*,'(A,*(1X,ES25.16E3))') 'TUNING', g1, g2, g3, g4, g5, g6, g7, g8, &
           f, ft, fv, fn, ftt, ftv, ftn, fvv, fvn, ((fnn(i,j), j=1,nc), i=1,nc)
      stop
    else if (trim(variant) == 'tune_ws') then
      read(*,*) p_vals(1:6)
      call thermopack_setWSparam(1, 2, p_vals(1), p_vals(2), p_vals(3), p_vals(4), &
           p_vals(5), p_vals(6))
      call thermopack_getWSparam(1, 2, g1, g2, g3, g4, g5, g6)
      call Fres(t, v_or_p, n, f, ft, fv, fn, ftt, ftv, fvv, ftn, fvn, fnn)
      write(*,'(A,*(1X,ES25.16E3))') 'TUNING', g1, g2, g3, g4, g5, g6, &
           f, ft, fv, fn, ftt, ftv, ftn, fvv, fvn, ((fnn(i,j), j=1,nc), i=1,nc)
      stop
    else if (trim(variant) == 'tune_twu_mc_beta') then
      read(*,*) p_vals(1:8) ! twu(1..3) for comp1, mc(1..3) for comp2, qbeta(1..2) for comp2
      call thermopack_setTWUparam(1, p_vals(1), p_vals(2), p_vals(3))
      call thermopack_setMCparam(2, p_vals(4), p_vals(5), p_vals(6))
      call thermopack_set_beta_corr(2, 2, 'QUANTUM', p_vals(7:8))
      call thermopack_getTWUparam(1, g1, g2, g3)
      call thermopack_getMCparam(2, g4, g5, g6)
      call Fres(t, v_or_p, n, f, ft, fv, fn, ftt, ftv, fvv, ftn, fvn, fnn)
      write(*,'(A,*(1X,ES25.16E3))') 'TUNING', g1, g2, g3, g4, g5, g6, &
           f, ft, fv, fn, ftt, ftv, ftn, fvv, fvn, ((fnn(i,j), j=1,nc), i=1,nc)
      stop
    else if (trim(variant) == 'tune_alpha_name') then
      read(*,*) p_vals(1:2)
      call thermopack_set_alpha_corr(1, 1, 'SW', p_vals(1:1))
      call thermopack_set_alpha_corr(1, 2, 'PR78', p_vals(2:2))
      call Fres(t, v_or_p, n, f, ft, fv, fn, ftt, ftv, fvv, ftn, fvn, fnn)
      write(*,'(A,*(1X,ES25.16E3))') 'TUNING', p_vals(1), p_vals(2), &
           f, ft, fv, fn, ftt, ftv, ftn, fvv, fvn, ((fnn(i,j), j=1,nc), i=1,nc)
      stop
    else if (trim(variant) == 'tune_kij_lij') then
      read(*,*) p_vals(1:3) ! kij_sym, kij_asym_override(2,1), lij_sym
      call thermopack_setkijandji(1, 2, p_vals(1))
      call thermopack_setkij(2, 1, p_vals(2))
      call thermopack_setlijandji(1, 2, p_vals(3))
      call thermopack_getkij(1, 2, g1)
      call thermopack_getkij(2, 1, g2)
      call thermopack_getlij(1, 2, g3)
      call thermopack_getlij(2, 1, g4)
      call Fres(t, v_or_p, n, f, ft, fv, fn, ftt, ftv, fvv, ftn, fvn, fnn)
      write(*,'(A,*(1X,ES25.16E3))') 'TUNING', g1, g2, g3, g4, &
           f, ft, fv, fn, ftt, ftv, ftn, fvv, fvn, ((fnn(i,j), j=1,nc), i=1,nc)
      stop
    else if (trim(variant) == 'tune_vshift') then
      read(*,*) ci_type
      read(*,*) p_vals(1:6)
      call thermopack_set_volume_shift_parameters(1, p_vals(1), p_vals(2), p_vals(3), &
           p_vals(4), p_vals(5), p_vals(6), ci_type)
      act_eos => get_active_eos()
      act_eos%volumeShiftId = 1 ! VOSHIFT_PENELOUX
      call thermopack_get_volume_shift_parameters(1, g1, g2, g3, g4, g5, g6, ci_type)
      call Fres(t, v_or_p, n, f, ft, fv, fn, ftt, ftv, fvv, ftn, fvn, fnn)
      write(*,'(A,*(1X,ES25.16E3))') 'TUNING', g1, g2, g3, g4, g5, g6, real(ci_type), &
           f, ft, fv, fn, ftt, ftv, ftn, fvv, fvn, ((fnn(i,j), j=1,nc), i=1,nc)
      stop
    else if (trim(variant) == 'tune_unifac_get') then
      read(*,*) i, j
      call getUNIFACgroupInteraction(unifdb, i, j, g1, g2, g3)
      call Fres(t, v_or_p, n, f, ft, fv, fn, ftt, ftv, fvv, ftn, fvn, fnn)
      write(*,'(A,*(1X,ES25.16E3))') 'TUNING', g1, g2, g3, &
           f, ft, fv, fn, ftt, ftv, ftn, fvv, fvn, ((fnn(i,j), j=1,nc), i=1,nc)
      stop
    end if
  end if

  if (trim(mode) == 'CPA_TVN') then
    call Fres(t, v_or_p, n, f, ft, fv, fn, ftt, ftv, fvv, ftn, fvn, fnn)
    act_eos => get_active_eos()
    ns = act_eos%assoc%numAssocSites
    allocate(xk(ns))
    xk = 0.5
    call act_eos%assoc%state%init(nc, t, v_or_p, n)
    call solve_for_X_k(act_eos, nc, xk, maxit=50, tol=1.0e-13)
    write(*,'(A,*(1X,ES25.16E3))') 'TVN', f, ft, fv, fn, ftt, ftv, ftn, fvv, fvn, &
         ((fnn(i,j), j=1,nc), i=1,nc), xk
  else if (trim(mode) == 'CUBIC_TVN') then
    call Fres(t, v_or_p, n, f, ft, fv, fn, ftt, ftv, fvv, ftn, fvn, fnn)
    write(*,'(A,*(1X,ES25.16E3))') 'TVN', f, ft, fv, fn, ftt, ftv, ftn, fvv, fvn, &
         ((fnn(i,j), j=1,nc), i=1,nc)
  else if (trim(mode) == 'CPA_TP' .or. trim(mode) == 'CUBIC_TP') then
    call specificVolume(t, v_or_p, n, phase_flag, vol)
    call zfac(t, v_or_p, n, phase_flag, z, dzdt, dzdp, dzdn)
    call enthalpy(t, v_or_p, n, phase_flag, hr, dhrdt, dhrdp, dhrdn, residual=.true.)
    call enthalpy(t, v_or_p, n, phase_flag, ht, dhtdt, dhtdp, dhtdn, residual=.false.)
    call entropy(t, v_or_p, n, phase_flag, sr, dsrdt, dsrdp, dsrdn, residual=.true.)
    call entropy(t, v_or_p, n, phase_flag, st, dstdt, dstdp, dstdn, residual=.false.)
    call thermo(t, v_or_p, n, phase_flag, lnphi, dlnphidt, dlnphidp, dlnphidn)
    sos = sound_velocity_2ph(t, v_or_p, x, x, x, 0.0, 1.0, phase_flag)
    write(*,'(A,*(1X,ES25.16E3))') 'TP', vol, z, dzdt, dzdp, dzdn, &
         hr, dhrdt, dhrdp, dhrdn, ht, dhtdt, dhtdp, dhtdn, &
         sr, dsrdt, dsrdp, dsrdn, st, dstdt, dstdp, dstdn, &
         lnphi, dlnphidt, dlnphidp, ((dlnphidn(i,j), j=1,nc), i=1,nc), sos
  else
    stop 21
  end if
end program cpa_cubic_tuning_probe
