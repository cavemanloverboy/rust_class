use crate::precision;
use crate::utils::zero_padded_i8_array;

/// These defaults were obtained from `precisions.h`
impl Default for precision {
    fn default() -> Self {
        precision {
            a_ini_over_a_today_default: 1e-14,
            background_Nloga: 3000,
            background_evolver: 1,
            tol_background_integration: 1e-10,
            background_integration_stepsize: 0.5,
            tol_initial_Omega_r: 1e-4,
            tol_M_ncdm: 1e-7,
            tol_ncdm: 1e-3,
            tol_ncdm_synchronous: 1e-3,
            tol_ncdm_newtonian: 1e-3,
            tol_ncdm_bg: 1e-5,
            tol_ncdm_initial_w: 1e-3,
            tol_tau_eq: 1e-6,
            Omega0_cdm_min_synchronous: 1e-10,
            tol_shooting_deltax: 1e-4,
            tol_shooting_deltaF: 1e-6,
            tol_shooting_deltax_rel: 1e-5,
            tol_fraction_accuracy: 1e-10,
            M_nfsm_threshold: 1e-4,
            sBBN_file: zero_padded_i8_array::<256>("class_public/external/bbn/sBBN_2017.dat"),
            thermo_z_initial: 5e6,
            thermo_z_initial_if_idm: 1e9,
            thermo_z_linear: 1e4,
            thermo_Nz_lin: 20000,
            thermo_Nz_log: 5000,
            thermo_evolver: 1,
            tol_thermo_integration: 1e-6,
            thermo_integration_stepsize: 0.1,
            thermo_rate_smoothing_radius: 50,
            z_end_reco_test: 500.0,
            primordial_black_hole_Nz: 75000,
            noninjection_Nz_log: 10000,
            noninjection_Nk_acc_diss: 500,
            k_min_acc_diss: 0.12,
            k_max_acc_diss: 1e6,
            z_wkb_acc_diss: 1e6,
            recfast_Heswitch: 6,
            recfast_fudge_He: 0.86,
            recfast_Hswitch: 1,
            recfast_fudge_H: 1.14,
            recfast_delta_fudge_H: -0.015,
            recfast_AGauss1: -0.14,
            recfast_AGauss2: 0.079,
            recfast_zGauss1: 7.28,
            recfast_zGauss2: 6.73,
            recfast_wGauss1: 0.18,
            recfast_wGauss2: 0.3,
            recfast_z_He_1: 8000.0,
            recfast_delta_z_He_1: 50.0,
            recfast_z_He_2: 4500.0,
            recfast_delta_z_He_2: 100.0,
            recfast_z_He_3: 3500.0,
            recfast_delta_z_He_3: 50.0,
            recfast_z_early_H_recombination: 2870.0,
            recfast_delta_z_early_H_recombination: 50.0,
            recfast_z_full_H_recombination: 1600.0,
            recfast_delta_z_full_H_recombination: 50.0,
            recfast_delta_z_reio: 2.0,
            recfast_x_He0_trigger: 0.995,
            recfast_x_He0_trigger2: 0.995,
            recfast_x_He0_trigger_delta: 0.05,
            recfast_x_H0_trigger: 0.995,
            recfast_x_H0_trigger2: 0.995,
            recfast_x_H0_trigger_delta: 0.05,
            recfast_z_switch_late: 800.0,
            hyrec_path: zero_padded_i8_array::<256>("class_public/external/HyRec2020/"),
            reionization_z_start_max: 50.0,
            reionization_sampling: 1.5e-2,
            reionization_optical_depth_tol: 1.0e-4,
            reionization_start_factor: 8.0,
            chi_z_Galli: zero_padded_i8_array::<256>("class_public/external/heating/Galli_et_al_2013.dat"),
            z_start_chi_approx: 2.0e3,
            k_min_tau0: 0.1,
            k_max_tau0_over_l_max: 2.4,
            k_step_sub: 0.05,
            k_step_super: 0.002,
            k_step_transition: 0.2,
            k_step_super_reduction: 0.1,
            k_per_decade_for_pk: 10.0,
            idmdr_boost_k_per_decade_for_pk: 1.0,
            k_per_decade_for_bao: 70.0,
            k_bao_center: 3.0,
            k_bao_width: 4.0,
            start_small_k_at_tau_c_over_tau_h: 0.0015,
            start_large_k_at_tau_h_over_tau_k: 0.07,
            tight_coupling_trigger_tau_c_over_tau_h: 0.015,
            tight_coupling_trigger_tau_c_over_tau_k: 0.01,
            tight_coupling_trigger_tau_c_over_tau_dmu_idm_g: 0.01,
            tight_coupling_trigger_tau_c_over_tau_R_idm_b: 0.01,
            start_sources_at_tau_c_over_tau_h: 0.008,
            tight_coupling_approximation: 5, //compromise_CLASS
            idm_dr_tight_coupling_trigger_tau_c_over_tau_k: 0.01,
            idm_dr_tight_coupling_trigger_tau_c_over_tau_h: 0.015,
            l_max_g: 12,
            l_max_pol_g: 10,
            l_max_dr: 17,
            l_max_ur: 17,
            l_max_idr: 17,
            l_max_ncdm: 17,
            l_max_g_ten: 5,
            l_max_pol_g_ten: 5,
            curvature_ini: 1.0,
            entropy_ini: 1.0,
            gw_ini: 1.0,
            perturbations_integration_stepsize: 0.5,
            perturbations_sampling_stepsize: 0.1,
            tol_perturbations_integration: 1e-5,
            c_gamma_k_H_square_max: 1e3,
            tol_tau_approx: 1e-10,
            radiation_streaming_approximation: 2, //rsa_MD_with_reio
            radiation_streaming_trigger_tau_over_tau_k: 45.0,
            radiation_streaming_trigger_tau_c_over_tau: 5.0,
            idr_streaming_approximation: 0, //rsa_idr_none
            idr_streaming_trigger_tau_over_tau_k: 50.0,
            idr_streaming_trigger_tau_c_over_tau: 10.0,
            ur_fluid_approximation: 2, //ufa_CLASS
            ur_fluid_trigger_tau_over_tau_k: 30.0,
            ncdm_fluid_approximation: 2, //ncdmfa_CLASS
            ncdm_fluid_trigger_tau_over_tau_k: 31.0,
            neglect_CMB_sources_below_visibility: 1e-3,
            evolver: 1, //ndf15
            k_per_decade_primordial: 10.0,
            primordial_inflation_ratio_min: 100.0,
            primordial_inflation_ratio_max: 1.0/50.0,
            primordial_inflation_phi_ini_maxit: 10000,
            primordial_inflation_pt_stepsize: 0.01,
            primordial_inflation_bg_stepsize: 0.005,
            primordial_inflation_tol_integration: 1e-3,
            primordial_inflation_attractor_precision_pivot: 0.001,
            primordial_inflation_attractor_precision_initial: 0.1,
            primordial_inflation_attractor_maxit: 10,
            primordial_inflation_tol_curvature: 1e-3,
            primordial_inflation_aH_ini_target: 0.9,
            primordial_inflation_end_dphi: 1e-10,
            primordial_inflation_end_logstep: 10.0, // TODO: Verify
            primordial_inflation_small_epsilon: 0.1,
            primordial_inflation_small_epsilon_tol: 0.01,
            primordial_inflation_extra_efolds: 2.0,
            l_linstep: 40,
            l_logstep: 1.12,
            hyper_x_min: 1e-5,
            hyper_sampling_flat: 8.0,
            hyper_sampling_curved_low_nu: 7.0,
            hyper_sampling_curved_high_nu: 3.0,
            hyper_nu_sampling_step: 1000.0,
            hyper_phi_min_abs: 1e-10,
            hyper_x_tol: 1e-4,
            hyper_flat_approximation_nu: 4000.0,
            q_linstep: 0.45,
            q_logstep_spline: 170.0,
            q_logstep_open: 6.0,
            q_logstep_trapzd: 20.0,
            q_numstep_transition: 250.0,
            transfer_neglect_delta_k_S_t0: 0.15,
            transfer_neglect_delta_k_S_t1: 0.04,
            transfer_neglect_delta_k_S_t2: 0.15,
            transfer_neglect_delta_k_S_e: 0.11,
            transfer_neglect_delta_k_V_t1: 1.0,
            transfer_neglect_delta_k_V_t2: 1.0,
            transfer_neglect_delta_k_V_e: 1.0,
            transfer_neglect_delta_k_V_b: 1.0,
            transfer_neglect_delta_k_T_t2: 0.2,
            transfer_neglect_delta_k_T_e: 0.25,
            transfer_neglect_delta_k_T_b: 0.1,
            transfer_neglect_late_source: 400.0,
            l_switch_limber: 10.0,
            l_switch_limber_for_nc_local_over_z: 100.0,
            l_switch_limber_for_nc_los_over_z: 30.0,
            selection_cut_at_sigma: 5.0,
            selection_sampling: 50.0,
            selection_sampling_bessel: 20.0,
            selection_sampling_bessel_los: 20.0, //in source it is set to selection_sampling_bessel, which was just set to 20
            selection_tophat_edge: 0.1,
            sigma_k_per_decade: 80.0,
            nonlinear_min_k_max: 5.0,
            k_max_for_pk_sigma8_min: 10.0,
            k_max_for_pk_sigma8_max: 100.0,
            halofit_min_k_nonlinear: 1e-4,
            halofit_k_per_decade: 80.0,
            halofit_sigma_precision: 0.05,
            halofit_tol_sigma: 1e-6,
            pk_eq_z_max: 5.0,
            pk_eq_Nzlog: 10,
            pk_eq_tol: 1e-7,
            hmcode_max_k_extra: 1e-6,
            hmcode_tol_sigma: 1e-6,
            n_hmcode_tables: 64,
            rmin_for_sigtab: 1e-5,
            rmax_for_sigtab: 1e-3,
            ainit_for_growtab: 1e-3,
            amax_for_growtab: 1.0,
            nsteps_for_p1h_integral: 256,
            mmin_for_p1h_integral: 1e-3,
            mmax_for_p1h_integral: 1e18,
            accurate_lensing: 0, //false
            num_mu_minus_lmax: 70,
            delta_l_max: 500,
            tol_gauss_legendre: std::f64::EPSILON,
            sd_z_min: 1.02e3,
            sd_z_max: 5.0e6,
            sd_z_size: 400,
            sd_x_min: 1e-2,
            sd_x_max: 5e1,
            sd_x_size: 500,
            tol_sd_detector: 1e-5,
            sd_external_path: zero_padded_i8_array::<256>("class_public/external/distortions"),
            smallest_allowed_variation: std::f64::EPSILON,
            error_message: [0_i8; 2048],
        }
    }
}
