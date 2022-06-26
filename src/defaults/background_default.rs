use crate::background;
use crate::utils::*;

impl Default for background {
    fn default() -> Self {
        background {
            #[doc = "< \\f$ H_0 \\f$: Hubble parameter (in fact, [\\f$H_0/c\\f$]) in \\f$ Mpc^{-1} \\f$"]
            H0: f64::default(),
            #[doc = "< reduced Hubble parameter"]
            h: f64::default(),
            #[doc = "< \\f$ \\Omega_{0 \\gamma} \\f$: photons"]
            Omega0_g: f64::default(),
            #[doc = "< \\f$ T_{cmb} \\f$: current CMB temperature in Kelvins"]
            T_cmb: f64::default(),
            #[doc = "< \\f$ \\Omega_{0 b} \\f$: baryons"]
            Omega0_b: f64::default(),
            #[doc = "< \\f$ \\Omega_{0 \\nu r} \\f$: ultra-relativistic neutrinos"]
            Omega0_ur: f64::default(),
            #[doc = "< \\f$ \\Omega_{0 cdm} \\f$: cold dark matter"]
            Omega0_cdm: f64::default(),
            #[doc = "< \\f$ \\Omega_{0 idm} \\f$: interacting dark matter with photons, baryons, and idr"]
            Omega0_idm: f64::default(),
            #[doc = "< \\f$ \\Omega_{0 idr} \\f$: interacting dark radiation"]
            Omega0_idr: f64::default(),
            #[doc = "< \\f$ T_{idr} \\f$: current temperature of interacting dark radiation in Kelvins"]
            T_idr: f64::default(),
            #[doc = "< \\f$ \\Omega_{0 dcdm}+\\Omega_{0 dr} \\f$: decaying cold dark matter (dcdm) decaying to dark radiation (dr)"]
            Omega0_dcdmdr: f64::default(),
            #[doc = "< \\f$ \\Omega_{ini,dcdm} \\f$: rescaled initial value for dcdm density (see 1407.2418 for definitions)"]
            Omega_ini_dcdm: f64::default(),
            #[doc = "< \\f$ \\Gamma_{dcdm} \\f$: decay constant for decaying cold dark matter"]
            Gamma_dcdm: f64::default(),
            tau_dcdm: f64::default(),
            #[doc = "< Number of distinguishable ncdm species"]
            N_ncdm: i32::default(),
            #[doc = "< list of filenames for tabulated p-s-d"]
            ncdm_psd_files: &mut i8::default(),
            #[doc = "< list of flags for each species, set to true if p-s-d is passed through file"]
            got_files: &mut i32::default(),
            #[doc = "< list of parameters for specifying/modifying ncdm p.s.d.'s, to be customized for given model"]
            #[doc = "(could be e.g. mixing angles)"]
            ncdm_psd_parameters: &mut f64::default(),
            #[doc = "< vector of masses of non-cold relic: dimensionless ratios m_ncdm/T_ncdm"]
            M_ncdm: &mut f64::default(),
            #[doc = "< list of ncdm masses in eV (inferred from M_ncdm and other parameters above)"]
            m_ncdm_in_eV: &mut f64::default(),
            #[doc = "< Omega0_ncdm for each species and for the total Omega0_ncdm"]
            Omega0_ncdm: &mut f64::default(),
            #[doc = "< Omega0_ncdm for each species and for the total Omega0_ncdm"]
            Omega0_ncdm_tot: f64::default(),
            #[doc = "< list of 1st parameters in p-s-d of non-cold relics: relative temperature"]
            #[doc = "T_ncdm1/T_gamma; and its default value"]
            T_ncdm: &mut f64::default(),
            #[doc = "< list of 1st parameters in p-s-d of non-cold relics: relative temperature"]
            #[doc = "T_ncdm1/T_gamma; and its default value"]
            T_ncdm_default: f64::default(),
            #[doc = "< list of 2nd parameters in p-s-d of non-cold relics: relative chemical potential"]
            #[doc = "ksi_ncdm1/T_ncdm1; and its default value"]
            ksi_ncdm: &mut f64::default(),
            #[doc = "< list of 2nd parameters in p-s-d of non-cold relics: relative chemical potential"]
            #[doc = "ksi_ncdm1/T_ncdm1; and its default value"]
            ksi_ncdm_default: f64::default(),
            #[doc = "< vector of degeneracy parameters in factor of p-s-d: 1 for one family of neutrinos"]
            #[doc = "(= one neutrino plus its anti-neutrino, total g*=1+1=2, so deg = 0.5 g*); and its"]
            #[doc = "default value"]
            deg_ncdm: &mut f64::default(),
            #[doc = "< vector of degeneracy parameters in factor of p-s-d: 1 for one family of neutrinos"]
            #[doc = "(= one neutrino plus its anti-neutrino, total g*=1+1=2, so deg = 0.5 g*); and its"]
            #[doc = "default value"]
            deg_ncdm_default: f64::default(),
            #[doc = "< Vector of numbers of q bins"]
            ncdm_input_q_size: &mut i32::default(),
            #[doc = "< Vector of maximum value of q"]
            ncdm_qmax: &mut f64::default(),
            #[doc = "< \\f$ \\Omega_{0_k} \\f$: curvature contribution"]
            Omega0_k: f64::default(),
            #[doc = "< \\f$ \\Omega_{0_\\Lambda} \\f$: cosmological constant"]
            Omega0_lambda: f64::default(),
            #[doc = "< \\f$ \\Omega_{0 de} \\f$: fluid"]
            Omega0_fld: f64::default(),
            #[doc = "< \\f$ \\Omega_{0 scf} \\f$: scalar field"]
            Omega0_scf: f64::default(),
            #[doc = "< flag switching on PPF perturbation equations instead of true fluid equations for perturbations. It could have been defined inside"]
            #[doc = "perturbation structure, but we leave it here in such way to have all fld parameters grouped."]
            use_ppf: i16::default(),
            #[doc = "< ppf parameter defined in eq. (16) of 0808.3125 [astro-ph]"]
            c_gamma_over_c_fld: f64::default(),
            #[doc = "< parametrisation scheme for fluid equation of state"]
            fluid_equation_of_state: u32::default(),
            #[doc = "< \\f$ w0_{DE} \\f$: current fluid equation of state parameter"]
            w0_fld: f64::default(),
            #[doc = "< \\f$ wa_{DE} \\f$: fluid equation of state parameter derivative"]
            wa_fld: f64::default(),
            #[doc = "< \\f$ c^2_{s~DE} \\f$: sound speed of the fluid in the frame comoving with the fluid (so, this is"]
            #[doc = "not [delta p/delta rho] in the synchronous or newtonian gauge!)"]
            cs2_fld: f64::default(),
            #[doc = "< \\f$ wa_{DE} \\f$: Early Dark Energy density parameter"]
            Omega_EDE: f64::default(),
            #[doc = "< list of parameters describing the scalar field potential"]
            scf_parameters: &mut f64::default(),
            #[doc = "< whether the scalar field has attractor initial conditions"]
            attractor_ic_scf: i16::default(),
            #[doc = "< index in scf_parameters used for tuning"]
            scf_tuning_index: i32::default(),
            #[doc = "< \\f$ \\phi(t_0) \\f$: scalar field initial value"]
            phi_ini_scf: f64::default(),
            #[doc = "< \\f$ d\\phi(t_0)/d\\tau \\f$: scalar field initial derivative wrt conformal time"]
            phi_prime_ini_scf: f64::default(),
            #[doc = "< size of scf_parameters"]
            scf_parameters_size: i32::default(),
            #[doc = "< finestructure constant for varying fundamental constants"]
            varconst_alpha: f64::default(),
            #[doc = "< electron mass for varying fundamental constants"]
            varconst_me: f64::default(),
            #[doc = "< dependence of the varying fundamental constants as a function of time"]
            varconst_dep: u32::default(),
            #[doc = "< redshift of transition between varied fundamental constants and normal fundamental constants in the 'varconst_instant' case"]
            varconst_transition_redshift: f64::default(),
            #[doc = "< age in Gyears"]
            age: f64::default(),
            #[doc = "< conformal age in Mpc"]
            conformal_age: f64::default(),
            #[doc = "< \\f$ K \\f$: Curvature parameter \\f$ K=-\\Omega0_k*a_{today}^2*H_0^2\\f$;"]
            K: f64::default(),
            #[doc = "< K/|K|: -1, 0 or 1"]
            sgnK: i32::default(),
            #[doc = "< so-called \"effective neutrino number\", computed at earliest time in interpolation table"]
            Neff: f64::default(),
            #[doc = "< \\f$ \\Omega_{0 dcdm} \\f$: decaying cold dark matter"]
            Omega0_dcdm: f64::default(),
            #[doc = "< \\f$ \\Omega_{0 dr} \\f$: decay radiation"]
            Omega0_dr: f64::default(),
            #[doc = "< total non-relativistic matter today"]
            Omega0_m: f64::default(),
            #[doc = "< total ultra-relativistic radiation today"]
            Omega0_r: f64::default(),
            #[doc = "< total dark energy density today, currently defined as 1 - Omega0_m - Omega0_r - Omega0_k"]
            Omega0_de: f64::default(),
            #[doc = "< total non-free-streaming matter, that is, cdm, baryons and wdm"]
            Omega0_nfsm: f64::default(),
            #[doc = "< scale factor at radiation/matter equality"]
            a_eq: f64::default(),
            #[doc = "< Hubble rate at radiation/matter equality [Mpc^-1]"]
            H_eq: f64::default(),
            #[doc = "< redshift at radiation/matter equality"]
            z_eq: f64::default(),
            #[doc = "< conformal time at radiation/matter equality [Mpc]"]
            tau_eq: f64::default(),
            #[doc = "< scale factor (in fact (a/a_0), see"]
            #[doc = "normalisation conventions explained"]
            #[doc = "at beginning of background.c)"]
            index_bg_a: i32::default(),
            #[doc = "< Hubble parameter in \\f$Mpc^{-1}\\f$"]
            index_bg_H: i32::default(),
            #[doc = "< its derivative w.r.t. conformal time"]
            index_bg_H_prime: i32::default(),
            #[doc = "< photon density"]
            index_bg_rho_g: i32::default(),
            #[doc = "< baryon density"]
            index_bg_rho_b: i32::default(),
            #[doc = "< cdm density"]
            index_bg_rho_cdm: i32::default(),
            #[doc = "< idm density"]
            index_bg_rho_idm: i32::default(),
            #[doc = "< cosmological constant density"]
            index_bg_rho_lambda: i32::default(),
            #[doc = "< fluid density"]
            index_bg_rho_fld: i32::default(),
            #[doc = "< fluid equation of state"]
            index_bg_w_fld: i32::default(),
            #[doc = "< density of interacting dark radiation"]
            index_bg_rho_idr: i32::default(),
            #[doc = "< relativistic neutrinos/relics density"]
            index_bg_rho_ur: i32::default(),
            #[doc = "< dcdm density"]
            index_bg_rho_dcdm: i32::default(),
            #[doc = "< dr density"]
            index_bg_rho_dr: i32::default(),
            #[doc = "< scalar field value"]
            index_bg_phi_scf: i32::default(),
            #[doc = "< scalar field derivative wrt conformal time"]
            index_bg_phi_prime_scf: i32::default(),
            #[doc = "< scalar field potential V"]
            index_bg_V_scf: i32::default(),
            #[doc = "< scalar field potential derivative V'"]
            index_bg_dV_scf: i32::default(),
            #[doc = "< scalar field potential second derivative V''"]
            index_bg_ddV_scf: i32::default(),
            #[doc = "< scalar field energy density"]
            index_bg_rho_scf: i32::default(),
            #[doc = "< scalar field pressure"]
            index_bg_p_scf: i32::default(),
            #[doc = "< scalar field pressure"]
            index_bg_p_prime_scf: i32::default(),
            #[doc = "< density of first ncdm species (others contiguous)"]
            index_bg_rho_ncdm1: i32::default(),
            #[doc = "< pressure of first ncdm species (others contiguous)"]
            index_bg_p_ncdm1: i32::default(),
            #[doc = "< another statistical momentum useful in ncdma approximation"]
            index_bg_pseudo_p_ncdm1: i32::default(),
            #[doc = "< Total density"]
            index_bg_rho_tot: i32::default(),
            #[doc = "< Total pressure"]
            index_bg_p_tot: i32::default(),
            #[doc = "< Conf. time derivative of total pressure"]
            index_bg_p_tot_prime: i32::default(),
            #[doc = "< relativistic density fraction (\\f$ \\Omega_{\\gamma} + \\Omega_{\\nu r} \\f$)"]
            index_bg_Omega_r: i32::default(),
            #[doc = "< critical density"]
            index_bg_rho_crit: i32::default(),
            #[doc = "< non-relativistic density fraction (\\f$ \\Omega_b + \\Omega_cdm + \\Omega_{\\nu nr} \\f$)"]
            index_bg_Omega_m: i32::default(),
            #[doc = "< conformal distance (from us) in Mpc"]
            index_bg_conf_distance: i32::default(),
            #[doc = "< angular diameter distance in Mpc"]
            index_bg_ang_distance: i32::default(),
            #[doc = "< luminosity distance in Mpc"]
            index_bg_lum_distance: i32::default(),
            #[doc = "< proper (cosmological) time in Mpc"]
            index_bg_time: i32::default(),
            #[doc = "< comoving sound horizon in Mpc"]
            index_bg_rs: i32::default(),
            #[doc = "< scale independent growth factor D(a) for CDM perturbations"]
            index_bg_D: i32::default(),
            #[doc = "< corresponding velocity growth factor [dlnD]/[dln a]"]
            index_bg_f: i32::default(),
            #[doc = "< value of fine structure constant in varying fundamental constants"]
            index_bg_varc_alpha: i32::default(),
            #[doc = "< value of effective electron mass in varying fundamental constants"]
            index_bg_varc_me: i32::default(),
            #[doc = "< size of background vector in the \"short format\""]
            bg_size_short: i32::default(),
            #[doc = "< size of background vector in the \"normal format\""]
            bg_size_normal: i32::default(),
            #[doc = "< size of background vector in the \"long format\""]
            bg_size: i32::default(),
            #[doc = "< number of lines (i.e. time-steps) in the four following array"]
            bt_size: i32::default(),
            #[doc = "< vector loga_table[index_loga] with values of log(a) (in fact \\f$ log(a/a0) \\f$, logarithm of relative scale factor compared to today)"]
            loga_table: &mut f64::default(),
            #[doc = "< vector tau_table[index_loga] with values of conformal time \\f$ \\tau \\f$ (in fact \\f$ a_0 c tau \\f$, see normalisation conventions explained at beginning of background.c)"]
            tau_table: &mut f64::default(),
            #[doc = "< vector z_table[index_loga] with values of \\f$ z \\f$ (redshift)"]
            z_table: &mut f64::default(),
            #[doc = "< table background_table[index_tau*pba->bg_size+pba->index_bg] with all other quantities (array of size bg_size*bt_size)"]
            background_table: &mut f64::default(),
            #[doc = "< vector d2tau_dz2_table[index_loga] with values of \\f$ d^2 \\tau / dz^2 \\f$ (conformal time)"]
            d2tau_dz2_table: &mut f64::default(),
            #[doc = "< vector d2z_dtau2_table[index_loga] with values of \\f$ d^2 z / d\\tau^2 \\f$ (conformal time)"]
            d2z_dtau2_table: &mut f64::default(),
            #[doc = "< table d2background_dtau2_table[index_loga*pba->bg_size+pba->index_bg] with values of \\f$ d^2 b_i / d\\log(a)^2 \\f$"]
            d2background_dloga2_table: &mut f64::default(),
            #[doc = "< {B} dcdm density"]
            index_bi_rho_dcdm: i32::default(),
            #[doc = "< {B} dr density"]
            index_bi_rho_dr: i32::default(),
            #[doc = "< {B} fluid density"]
            index_bi_rho_fld: i32::default(),
            #[doc = "< {B} scalar field value"]
            index_bi_phi_scf: i32::default(),
            #[doc = "< {B} scalar field derivative wrt conformal time"]
            index_bi_phi_prime_scf: i32::default(),
            #[doc = "< {C} proper (cosmological) time in Mpc"]
            index_bi_time: i32::default(),
            #[doc = "< {C} sound horizon"]
            index_bi_rs: i32::default(),
            #[doc = "< {C} conformal time in Mpc"]
            index_bi_tau: i32::default(),
            #[doc = "< {C} scale independent growth factor D(a) for CDM perturbations."]
            index_bi_D: i32::default(),
            #[doc = "< {C} D satisfies \\f$ [D''(\\tau)=-aHD'(\\tau)+3/2 a^2 \\rho_M D(\\tau) \\f$"]
            index_bi_D_prime: i32::default(),
            #[doc = "< Number of {B} parameters"]
            bi_B_size: i32::default(),
            #[doc = "< Number of {B}+{C} parameters"]
            bi_size: i32::default(),
            #[doc = "< presence of cold dark matter?"]
            has_cdm: i16::default(),
            #[doc = "< presence of interacting dark matter with photons, baryons, and idr"]
            has_idm: i16::default(),
            #[doc = "< presence of decaying cold dark matter?"]
            has_dcdm: i16::default(),
            #[doc = "< presence of relativistic decay radiation?"]
            has_dr: i16::default(),
            #[doc = "< presence of a scalar field?"]
            has_scf: i16::default(),
            #[doc = "< presence of non-cold dark matter?"]
            has_ncdm: i16::default(),
            #[doc = "< presence of cosmological constant?"]
            has_lambda: i16::default(),
            #[doc = "< presence of fluid with constant w and cs2?"]
            has_fld: i16::default(),
            #[doc = "< presence of ultra-relativistic neutrinos/relics?"]
            has_ur: i16::default(),
            #[doc = "< presence of interacting dark radiation?"]
            has_idr: i16::default(),
            #[doc = "< presence of global spatial curvature?"]
            has_curvature: i16::default(),
            #[doc = "< presence of varying fundamental constants?"]
            has_varconst: i16::default(),
            #[doc = "< Vector of integers according to quadrature strategy."]
            ncdm_quadrature_strategy: &mut i32::default(),
            #[doc = "< Pointers to vectors of background sampling in q"]
            q_ncdm_bg: return_raw_dbl_ptr::<f64>(),
            #[doc = "< Pointers to vectors of corresponding quadrature weights w"]
            w_ncdm_bg: return_raw_dbl_ptr::<f64>(),
            #[doc = "< Pointers to vectors of perturbation sampling in q"]
            q_ncdm: return_raw_dbl_ptr::<f64>(),
            #[doc = "< Pointers to vectors of corresponding quadrature weights w"]
            w_ncdm: return_raw_dbl_ptr::<f64>(),
            #[doc = "< Pointers to vectors of logarithmic derivatives of p-s-d"]
            dlnf0_dlnq_ncdm:return_raw_dbl_ptr::<f64>(),
            #[doc = "< Size of the q_ncdm_bg arrays"]
            q_size_ncdm_bg: &mut i32::default(),
            #[doc = "< Size of the q_ncdm arrays"]
            q_size_ncdm: &mut i32::default(),
            #[doc = "< List of normalization factors for calculating energy density etc."]
            factor_ncdm: &mut f64::default(),
            #[doc = "< flag is set to true if shooting failed."]
            shooting_failed: i16::default(),
            #[doc = "< Error message from shooting failed."]
            shooting_error: zero_padded_i8_array::<2048>(""),
            #[doc = "< flag regulating the amount of information sent to standard output (none if set to zero)"]
            background_verbose: i16::default(),
            #[doc = "< zone for writing error messages"]
            error_message: zero_padded_i8_array::<2048>(""),
        }
    }
}