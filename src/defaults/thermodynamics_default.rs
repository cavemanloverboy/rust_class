use crate::{thermodynamics, *};
use crate::utils::*;

impl Default for thermodynamics {
    fn default() -> Self {
        thermodynamics {
            #[doc = "< \\f$ Y_{He} \\f$: primordial helium mass fraction rho_He/(rho_H+rho_He),"]
            #[doc = "close but not exactly equal to the density fraction 4*n_He/(n_H+4*n_He)"]
            YHe: f64::default(),
            #[doc = "< Related to variation of fundamental constants (sensitivity of YHe to alpha)"]
            bbn_alpha_sensitivity: f64::default(),
            #[doc = "< recombination code"]
            recombination: recombination_algorithm::default(),
            #[doc = "< photo-ionization coefficient mode of the recfast algorithm"]
            recfast_photoion_mode: recfast_photoion_modes::default(),
            #[doc = "< reionization scheme"]
            reio_parametrization: reionization_parametrization::default(),
            #[doc = "< is the input parameter the reionization redshift or optical depth?"]
            reio_z_or_tau: reionization_z_or_tau::default(),
            #[doc = "< if above set to tau, input value of reionization optical depth"]
            tau_reio: f64::default(),
            #[doc = "< if above set to z,   input value of reionization redshift"]
            z_reio: f64::default(),
            #[doc = "< do we want to include in computation derivatives of baryon sound speed?"]
            compute_cb2_derivatives: ::std::os::raw::c_short::default(),
            #[doc = "< do we want to compute the simplest analytic approximation to the photon damping (or diffusion) scale?"]
            compute_damping_scale: ::std::os::raw::c_short::default(),
            #[doc = "< Do we have idm with baryons?"]
            has_idm_b: ::std::os::raw::c_short::default(),
            #[doc = "< Do we have idm with photons?"]
            has_idm_g: ::std::os::raw::c_short::default(),
            #[doc = "< Do we have idm with dark radiation?"]
            has_idm_dr: ::std::os::raw::c_short::default(),
            #[doc = "< width of H reionization"]
            reionization_width: f64::default(),
            #[doc = "< shape of H reionization"]
            reionization_exponent: f64::default(),
            #[doc = "< redshift for of helium reionization"]
            helium_fullreio_redshift: f64::default(),
            #[doc = "< width of helium reionization"]
            helium_fullreio_width: f64::default(),
            #[doc = "< with how many bins do we want to describe reionization?"]
            binned_reio_num: ::std::os::raw::c_int::default(),
            #[doc = "< central z value for each bin"]
            binned_reio_z: &mut f64::default(),
            #[doc = "< imposed \\f$ X_e(z)\\f$ value at center of each bin"]
            binned_reio_xe: &mut f64::default(),
            #[doc = "< sharpness of tanh() step interpolating between binned values"]
            binned_reio_step_sharpness: f64::default(),
            #[doc = "< with how many jumps do we want to describe reionization?"]
            many_tanh_num: ::std::os::raw::c_int::default(),
            #[doc = "< central z value for each tanh jump"]
            many_tanh_z: &mut f64::default(),
            #[doc = "< imposed \\f$ X_e(z)\\f$ value at the end of each jump (ie at later times)"]
            many_tanh_xe: &mut f64::default(),
            #[doc = "< sharpness of tanh() steps"]
            many_tanh_width: f64::default(),
            #[doc = "< with how many jumps do we want to describe reionization?"]
            reio_inter_num: ::std::os::raw::c_int::default(),
            #[doc = "< discrete z values"]
            reio_inter_z: &mut f64::default(),
            #[doc = "< discrete \\f$ X_e(z)\\f$ values"]
            reio_inter_xe: &mut f64::default(),
            #[doc = "< true if some exotic mechanism"]
            #[doc = "injects energy and affects the"]
            #[doc = "evolution of ionization and/or"]
            #[doc = "temperature and/or other"]
            #[doc = "thermodynamics variables that are"]
            #[doc = "relevant for the calculation of CMB"]
            #[doc = "anisotropies (and spectral"]
            #[doc = "distorsions if requested)."]
            has_exotic_injection: ::std::os::raw::c_short::default(),
            #[doc = "< structure to store exotic energy injections and their energy deposition"]
            in_: injection::default(),
            #[doc = "< parameter describing CDM annihilation (f <sigma*v> / m_cdm, see e.g. 0905.0003)"]
            annihilation: f64::default(),
            #[doc = "< flag to specify if we want to use the on-the-spot approximation"]
            has_on_the_spot: ::std::os::raw::c_short::default(),
            #[doc = "< parameter describing CDM decay (f/tau, see e.g. 1109.6322)"]
            decay: f64::default(),
            #[doc = "< if this parameter is non-zero,"]
            #[doc = "the function F(z)=(f <sigma*v> /"]
            #[doc = "m_cdm)(z) will be a parabola in"]
            #[doc = "log-log scale between zmin and"]
            #[doc = "zmax, with a curvature given by"]
            #[doc = "annihlation_variation (must be"]
            #[doc = "negative), and with a maximum in"]
            #[doc = "zmax; it will be constant outside"]
            #[doc = "this range"]
            annihilation_variation: f64::default(),
            #[doc = "< if annihilation_variation is non-zero,"]
            #[doc = "this is the value of z at which the"]
            #[doc = "parameter annihilation is defined, i.e."]
            #[doc = "F(annihilation_z)=annihilation"]
            annihilation_z: f64::default(),
            #[doc = "< if annihilation_variation is non-zero,"]
            #[doc = "redshift above which annihilation rate"]
            #[doc = "is maximal"]
            annihilation_zmax: f64::default(),
            #[doc = "< if annihilation_variation is non-zero,"]
            #[doc = "redshift below which annihilation rate"]
            #[doc = "is constant"]
            annihilation_zmin: f64::default(),
            #[doc = "< takes the contribution of DM annihilation in halos into account"]
            annihilation_f_halo: f64::default(),
            #[doc = "< characteristic redshift for DM annihilation in halos"]
            annihilation_z_halo: f64::default(),
            #[doc = "< presence of varying fundamental constants?"]
            has_varconst: ::std::os::raw::c_short::default(),
            #[doc = "< ionization fraction \\f$ x_e \\f$"]
            index_th_xe: ::std::os::raw::c_int::default(),
            #[doc = "< Thomson scattering rate \\f$ d \\kappa / d \\tau\\f$ (units 1/Mpc)"]
            index_th_dkappa: ::std::os::raw::c_int::default(),
            #[doc = "< Baryon drag optical depth"]
            index_th_tau_d: ::std::os::raw::c_int::default(),
            #[doc = "< scattering rate derivative \\f$ d^2 \\kappa / d \\tau^2 \\f$"]
            index_th_ddkappa: ::std::os::raw::c_int::default(),
            #[doc = "< scattering rate second derivative \\f$ d^3 \\kappa / d \\tau^3 \\f$"]
            index_th_dddkappa: ::std::os::raw::c_int::default(),
            #[doc = "< \\f$ exp^{-\\kappa} \\f$"]
            index_th_exp_m_kappa: ::std::os::raw::c_int::default(),
            #[doc = "< visibility function \\f$ g = (d \\kappa / d \\tau) * exp^{-\\kappa} \\f$"]
            index_th_g: ::std::os::raw::c_int::default(),
            #[doc = "< visibility function derivative \\f$ (d g / d \\tau) \\f$"]
            index_th_dg: ::std::os::raw::c_int::default(),
            #[doc = "< visibility function second derivative \\f$ (d^2 g / d \\tau^2) \\f$"]
            index_th_ddg: ::std::os::raw::c_int::default(),
            #[doc = "< idm temperature \\f$ T_idm \\f$"]
            index_th_T_idm: ::std::os::raw::c_int::default(),
            #[doc = "< idm sound speed squared \\f$ c_idm^2 \\f$"]
            index_th_c2_idm: ::std::os::raw::c_int::default(),
            #[doc = "< idr temperature \\f$ T_idr \\f$"]
            index_th_T_idr: ::std::os::raw::c_int::default(),
            #[doc = "< scattering rate of idr with idm_g_dr (i.e. idr opacity to idm_g_dr scattering) (units 1/Mpc)"]
            index_th_dmu_idm_dr: ::std::os::raw::c_int::default(),
            #[doc = "< derivative of the idm_g_dr scattering rate"]
            index_th_ddmu_idm_dr: ::std::os::raw::c_int::default(),
            #[doc = "< second derivative of the idm_g_dr scattering rate"]
            index_th_dddmu_idm_dr: ::std::os::raw::c_int::default(),
            #[doc = "< idr self-interaction rate"]
            index_th_dmu_idr: ::std::os::raw::c_int::default(),
            #[doc = "< optical depth of idm_dr (due to interactions with idr)"]
            index_th_tau_idm_dr: ::std::os::raw::c_int::default(),
            #[doc = "< optical depth of idr (due to self-interactions)"]
            index_th_tau_idr: ::std::os::raw::c_int::default(),
            #[doc = "< visibility function of idm_idr"]
            index_th_g_idm_dr: ::std::os::raw::c_int::default(),
            #[doc = "< idm_g scattering rate \\f$ d \\mu / d \\tau\\f$  (analogous to Thomson scattering) (see 1802.06589 for details)"]
            index_th_dmu_idm_g: ::std::os::raw::c_int::default(),
            #[doc = "< derivative of idm_g scattering, \\f$ d^2 \\mu / d \\tau^2 \\f$"]
            index_th_ddmu_idm_g: ::std::os::raw::c_int::default(),
            #[doc = "< second derivative of idm_g scattering rate, \\f$ d^3 \\mu / d \\tau^3 \\f$"]
            index_th_dddmu_idm_g: ::std::os::raw::c_int::default(),
            #[doc = "< \\f$ exp^{-\\mu} \\f$"]
            index_th_exp_mu_idm_g: ::std::os::raw::c_int::default(),
            #[doc = "< idm_b interaction coefficient"]
            index_th_R_idm_b: ::std::os::raw::c_int::default(),
            #[doc = "< derivative of idm_b interaction coefficient wrt conformal time"]
            index_th_dR_idm_b: ::std::os::raw::c_int::default(),
            #[doc = "< second derivative of ibm_b interaction coefficient wrt conformal time"]
            index_th_ddR_idm_b: ::std::os::raw::c_int::default(),
            #[doc = "< baryon temperature \\f$ T_b \\f$"]
            index_th_Tb: ::std::os::raw::c_int::default(),
            #[doc = "< derivative of baryon temperature"]
            index_th_dTb: ::std::os::raw::c_int::default(),
            #[doc = "< baryon equation of state parameter \\f$ w_b = k_B T_b / \\mu \\f$"]
            index_th_wb: ::std::os::raw::c_int::default(),
            #[doc = "< squared baryon adiabatic sound speed \\f$ c_b^2 \\f$"]
            index_th_cb2: ::std::os::raw::c_int::default(),
            #[doc = "< derivative wrt conformal time of squared baryon sound speed \\f$ d [c_b^2] / d \\tau \\f$ (only computed if some non-minimal tight-coupling schemes is requested)"]
            index_th_dcb2: ::std::os::raw::c_int::default(),
            #[doc = "< second derivative wrt conformal time of squared baryon sound speed  \\f$ d^2 [c_b^2] / d \\tau^2 \\f$ (only computed if some non0-minimal tight-coupling schemes is requested)"]
            index_th_ddcb2: ::std::os::raw::c_int::default(),
            #[doc = "< maximum variation rate of \\f$ exp^{-\\kappa}\\f$, g and \\f$ (d g / d \\tau) \\f$, used for computing integration step in perturbation module"]
            index_th_rate: ::std::os::raw::c_int::default(),
            #[doc = "< simple analytic approximation to the photon comoving damping scale"]
            index_th_r_d: ::std::os::raw::c_int::default(),
            #[doc = "< size of thermodynamics vector"]
            th_size: ::std::os::raw::c_int::default(),
            #[doc = "< number of lines (redshift steps) in the tables"]
            tt_size: ::std::os::raw::c_int::default(),
            #[doc = "< vector z_table[index_z] with values of redshift (vector of size tt_size)"]
            z_table: &mut f64::default(),
            #[doc = "< vector tau_table[index_tau] with values of conformal time (vector of size tt_size)"]
            tau_table: &mut f64::default(),
            #[doc = "< table thermodynamics_table[index_z*pth->tt_size+pba->index_th] with all other quantities (array of size th_size*tt_size)"]
            thermodynamics_table: &mut f64::default(),
            #[doc = "< table d2thermodynamics_dz2_table[index_z*pth->tt_size+pba->index_th] with values of \\f$ d^2 t_i / dz^2 \\f$ (array of size th_size*tt_size)"]
            d2thermodynamics_dz2_table: &mut f64::default(),
            #[doc = "< z at which the visibility reaches its maximum (= recombination redshift)"]
            z_rec: f64::default(),
            #[doc = "< conformal time at which the visibility reaches its maximum (= recombination time)"]
            tau_rec: f64::default(),
            #[doc = "< comoving sound horizon at recombination"]
            rs_rec: f64::default(),
            #[doc = "< physical sound horizon at recombination"]
            ds_rec: f64::default(),
            #[doc = "< conformal angular diameter distance to recombination"]
            ra_rec: f64::default(),
            #[doc = "< physical angular diameter distance to recombination"]
            da_rec: f64::default(),
            #[doc = "< comoving photon damping scale at recombination"]
            rd_rec: f64::default(),
            #[doc = "< redshift at which photon optical depth crosses one"]
            z_star: f64::default(),
            #[doc = "< confirmal time at which photon optical depth crosses one"]
            tau_star: f64::default(),
            #[doc = "< comoving sound horizon at z_star"]
            rs_star: f64::default(),
            #[doc = "< physical sound horizon at z_star"]
            ds_star: f64::default(),
            #[doc = "< conformal angular diameter distance to z_star"]
            ra_star: f64::default(),
            #[doc = "< physical angular diameter distance to z_star"]
            da_star: f64::default(),
            #[doc = "< comoving photon damping scale at z_star"]
            rd_star: f64::default(),
            #[doc = "< baryon drag redshift"]
            z_d: f64::default(),
            #[doc = "< baryon drag time"]
            tau_d: f64::default(),
            #[doc = "< physical sound horizon at baryon drag"]
            ds_d: f64::default(),
            #[doc = "< comoving sound horizon at baryon drag"]
            rs_d: f64::default(),
            #[doc = "< at at which the visibility goes below a fixed fraction of the maximum visibility, used for an approximation in perturbation module"]
            tau_cut: f64::default(),
            #[doc = "< [ratio ra_rec / (tau0-tau_rec)]: gives CMB rescaling in angular space relative to flat model (=1 for curvature K=0)"]
            angular_rescaling: f64::default(),
            #[doc = "< minimum value of tau at which free-streaming approximation can be switched on"]
            tau_free_streaming: f64::default(),
            #[doc = "< trigger for dark radiation free streaming approximation (idm-idr)"]
            tau_idr_free_streaming: f64::default(),
            #[doc = "< decoupling time for idr"]
            tau_idr: f64::default(),
            #[doc = "< decoupling time for idm from idr"]
            tau_idm_dr: f64::default(),
            #[doc = "< initial conformal time at which thermodynamical variables have been be integrated"]
            tau_ini: f64::default(),
            #[doc = "< \\f$ f_{He} \\f$: primordial helium-to-hydrogen nucleon ratio 4*n_He/n_H"]
            fHe: f64::default(),
            #[doc = "< total number density of electrons today (free or not)"]
            n_e: f64::default(),
            #[doc = "< dark matter mass for idm"]
            m_idm: f64::default(),
            #[doc = "< strength of the coupling between interacting dark matter and interacting dark radiation (idm-idr)"]
            a_idm_dr: f64::default(),
            #[doc = "< strength of the self coupling for interacting dark radiation (idr-idr)"]
            b_idr: f64::default(),
            #[doc = "< temperature dependence of the interactions between dark matter and dark radiation"]
            n_index_idm_dr: f64::default(),
            #[doc = "< cross section between interacting dark matter and baryons"]
            cross_idm_b: f64::default(),
            #[doc = "< temperature dependence of the interactions between dark matter and baryons"]
            n_index_idm_b: ::std::os::raw::c_int::default(),
            #[doc = "< numerical n-dependent coefficient for idm_b"]
            n_coeff_idm_b: f64::default(),
            #[doc = "< cross section between interacting dark matter and photons"]
            cross_idm_g: f64::default(),
            #[doc = "< ratio between idm_g cross section and idm mass"]
            u_idm_g: f64::default(),
            #[doc = "< temperature dependence of the interactions between dark matter and photons"]
            n_index_idm_g: ::std::os::raw::c_int::default(),
            #[doc = "< flag for calling thermodynamics_at_z and find position in interpolation table normally"]
            inter_normal: ::std::os::raw::c_short::default(),
            #[doc = "< flag for calling thermodynamics_at_z and find position in interpolation table starting from previous position in previous call"]
            inter_closeby: ::std::os::raw::c_short::default(),
            #[doc = "< flag regulating the amount of information sent to standard output (none if set to zero)"]
            thermodynamics_verbose: ::std::os::raw::c_short::default(),
            #[doc = "< flag regulating the amount of information sent to standard output from hyrec (none if set to zero)"]
            hyrec_verbose: ::std::os::raw::c_short::default(),
            #[doc = "< zone for writing error messages"]
            error_message: zero_padded_i8_array::<2048>(""),
        }
    }
}



impl Default for injection{
    fn default() -> Self {
        injection {
            DM_annihilation_efficiency: f64::default(),
            DM_annihilation_cross_section: f64::default(),
            DM_annihilation_mass: f64::default(),
            DM_annihilation_fraction: f64::default(),
            DM_annihilation_variation: f64::default(),
            DM_annihilation_z: f64::default(),
            DM_annihilation_zmax: f64::default(),
            DM_annihilation_zmin: f64::default(),
            DM_annihilation_f_halo: f64::default(),
            DM_annihilation_z_halo: f64::default(),
            DM_decay_fraction: f64::default(),
            DM_decay_Gamma: f64::default(),
            PBH_evaporation_fraction: f64::default(),
            PBH_evaporation_mass: f64::default(),
            PBH_accretion_fraction: f64::default(),
            PBH_accretion_mass: f64::default(),
            PBH_accretion_recipe: PBH_accretion_approx::default(),
            PBH_accretion_relative_velocities: f64::default(),
            PBH_accretion_ADAF_delta: f64::default(),
            PBH_accretion_eigenvalue: f64::default(),
            f_eff_type: ::std::os::raw::c_int::default(),
            f_eff_file: zero_padded_i8_array::<256>(""),
            chi_type: ::std::os::raw::c_int::default(),
            chi_z_file: zero_padded_i8_array::<256>(""),
            chi_x_file: zero_padded_i8_array::<256>(""),
            Nz_size: ::std::os::raw::c_int::default(),
            z_initial: f64::default(),
            z_start_chi_approx: f64::default(),
            H0: f64::default(),
            T_g0: f64::default(),
            Omega0_b: f64::default(),
            Omega0_cdm: f64::default(),
            rho0_cdm: f64::default(),
            H: f64::default(),
            a: f64::default(),
            t: f64::default(),
            rho_g: f64::default(),
            rho_b: f64::default(),
            rho_cdm: f64::default(),
            fHe: f64::default(),
            N_e0: f64::default(),
            heat_capacity: f64::default(),
            T_b: f64::default(),
            x_e: f64::default(),
            nH: f64::default(),
            z_table: &mut f64::default(),
            z_size: ::std::os::raw::c_int::default(),
            tol_z_table: f64::default(),
            filled_until_index_z: ::std::os::raw::c_int::default(),
            filled_until_z: f64::default(),
            last_index_z_feff: ::std::os::raw::c_int::default(),
            last_index_z_chi: ::std::os::raw::c_int::default(),
            last_index_z_inj: ::std::os::raw::c_int::default(),
            last_index_z: ::std::os::raw::c_int::default(),
            index_z_store: ::std::os::raw::c_int::default(),
            last_index_x_chi: ::std::os::raw::c_int::default(),
            PBH_z_evaporation: f64::default(),
            PBH_QCD_activation: f64::default(),
            PBH_table_z: &mut f64::default(),
            PBH_table_mass: &mut f64::default(),
            PBH_table_mass_dd: &mut f64::default(),
            PBH_table_F: &mut f64::default(),
            PBH_table_F_dd: &mut f64::default(),
            Nz_PBH: ::std::os::raw::c_int::default(),
            injection_table: return_raw_dbl_ptr::<f64>(),
            index_inj_cool: ::std::os::raw::c_int::default(),
            index_inj_diss: ::std::os::raw::c_int::default(),
            index_inj_DM_ann: ::std::os::raw::c_int::default(),
            index_inj_DM_dec: ::std::os::raw::c_int::default(),
            index_inj_PBH_eva: ::std::os::raw::c_int::default(),
            index_inj_PBH_acc: ::std::os::raw::c_int::default(),
            index_inj_tot: ::std::os::raw::c_int::default(),
            inj_size: ::std::os::raw::c_int::default(),
            #[doc = " All contributions + total"]
            f_eff: f64::default(),
            feff_z_size: ::std::os::raw::c_int::default(),
            feff_table: &mut f64::default(),
            chiz_table: &mut f64::default(),
            chiz_size: ::std::os::raw::c_int::default(),
            chix_table: &mut f64::default(),
            chix_size: ::std::os::raw::c_int::default(),
            deposition_table: return_raw_dbl_ptr::<f64>(),
            chi: &mut f64::default(),
            index_dep_heat: ::std::os::raw::c_int::default(),
            index_dep_ionH: ::std::os::raw::c_int::default(),
            index_dep_ionHe: ::std::os::raw::c_int::default(),
            index_dep_lya: ::std::os::raw::c_int::default(),
            index_dep_lowE: ::std::os::raw::c_int::default(),
            dep_size: ::std::os::raw::c_int::default(),
            pvecdeposition: &mut f64::default(),
            has_exotic_injection: ::std::os::raw::c_int::default(),
            has_DM_ann: ::std::os::raw::c_int::default(),
            has_DM_dec: ::std::os::raw::c_int::default(),
            has_PBH_eva: ::std::os::raw::c_int::default(),
            has_PBH_acc: ::std::os::raw::c_int::default(),
            to_store: ::std::os::raw::c_int::default(),
            error_message: zero_padded_i8_array::<2048>(""),
        }
    }
}