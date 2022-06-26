use crate::{perturbations, tensor_methods, selection_type, possible_gauges};
use crate::utils::*;

impl Default for perturbations {
    fn default() -> Self {
        perturbations {
            #[doc = "< do we need to compute perturbations at all ?"]
            has_perturbations: ::std::os::raw::c_short::default(),
            #[doc = "< do we need any harmonic space spectrum \\f$ C_l \\f$ (and hence Bessel functions, transfer functions, ...)?"]
            has_cls: ::std::os::raw::c_short::default(),
            #[doc = "< do we need scalars?"]
            has_scalars: ::std::os::raw::c_short::default(),
            #[doc = "< do we need vectors?"]
            has_vectors: ::std::os::raw::c_short::default(),
            #[doc = "< do we need tensors?"]
            has_tensors: ::std::os::raw::c_short::default(),
            #[doc = "< do we need adiabatic mode?"]
            has_ad: ::std::os::raw::c_short::default(),
            #[doc = "< do we need isocurvature bi mode?"]
            has_bi: ::std::os::raw::c_short::default(),
            #[doc = "< do we need isocurvature cdi mode?"]
            has_cdi: ::std::os::raw::c_short::default(),
            #[doc = "< do we need isocurvature nid mode?"]
            has_nid: ::std::os::raw::c_short::default(),
            #[doc = "< do we need isocurvature niv mode?"]
            has_niv: ::std::os::raw::c_short::default(),
            #[doc = " Do we want to consider perturbed temperature and ionization fraction?"]
            has_perturbed_recombination: ::std::os::raw::c_short::default(),
            #[doc = "< way to treat neutrinos in tensor perturbations(neglect, approximate as massless, take exact equations)"]
            tensor_method: tensor_methods::default(),
            #[doc = "< will we evolve ur tensor perturbations (either because we have ur species, or we have ncdm species with massless approximation) ?"]
            evolve_tensor_ur: ::std::os::raw::c_short::default(),
            #[doc = "< will we evolve ncdm tensor perturbations (if we have ncdm species and we use the exact method) ?"]
            evolve_tensor_ncdm: ::std::os::raw::c_short::default(),
            #[doc = "< do we need \\f$ C_l \\f$'s for CMB temperature?"]
            has_cl_cmb_temperature: ::std::os::raw::c_short::default(),
            #[doc = "< do we need \\f$ C_l \\f$'s for CMB polarization?"]
            has_cl_cmb_polarization: ::std::os::raw::c_short::default(),
            #[doc = "< do we need \\f$ C_l \\f$'s for CMB lensing potential?"]
            has_cl_cmb_lensing_potential: ::std::os::raw::c_short::default(),
            #[doc = "< do we need \\f$ C_l \\f$'s for galaxy lensing potential?"]
            has_cl_lensing_potential: ::std::os::raw::c_short::default(),
            #[doc = "< do we need \\f$ C_l \\f$'s for density number count?"]
            has_cl_number_count: ::std::os::raw::c_short::default(),
            #[doc = "< do we need matter Fourier spectrum?"]
            has_pk_matter: ::std::os::raw::c_short::default(),
            #[doc = "< do we need to output individual matter density transfer functions?"]
            has_density_transfers: ::std::os::raw::c_short::default(),
            #[doc = "< do we need to output individual matter velocity transfer functions?"]
            has_velocity_transfers: ::std::os::raw::c_short::default(),
            #[doc = "< do we need to output individual transfer functions for scalar metric perturbations?"]
            has_metricpotential_transfers: ::std::os::raw::c_short::default(),
            #[doc = "< should we convert density and velocity transfer functions to Nbody gauge?"]
            has_Nbody_gauge_transfers: ::std::os::raw::c_short::default(),
            #[doc = "< do we want to compute non-linear corrections with an algorithm relying on delta_m (like halofit)?"]
            has_nl_corrections_based_on_delta_m: ::std::os::raw::c_short::default(),
            #[doc = "< in dCl, do we want density terms ?"]
            has_nc_density: ::std::os::raw::c_short::default(),
            #[doc = "< in dCl, do we want redshift space distortion terms ?"]
            has_nc_rsd: ::std::os::raw::c_short::default(),
            #[doc = "< in dCl, do we want lensing terms ?"]
            has_nc_lens: ::std::os::raw::c_short::default(),
            #[doc = "< in dCl, do we want gravity terms ?"]
            has_nc_gr: ::std::os::raw::c_short::default(),
            #[doc = "< maximum l value for CMB scalars \\f$ C_l \\f$'s"]
            l_scalar_max: ::std::os::raw::c_int::default(),
            #[doc = "< maximum l value for CMB vectors \\f$ C_l \\f$'s"]
            l_vector_max: ::std::os::raw::c_int::default(),
            #[doc = "< maximum l value for CMB tensors \\f$ C_l \\f$'s"]
            l_tensor_max: ::std::os::raw::c_int::default(),
            #[doc = "< maximum l value for LSS \\f$ C_l \\f$'s (density and lensing potential in  bins)"]
            l_lss_max: ::std::os::raw::c_int::default(),
            #[doc = "< maximum value of k in 1/Mpc required for the output of P(k,z) and T(k,z)"]
            k_max_for_pk: f64::default(),
            #[doc = "< number of selection functions"]
            #[doc = "(i.e. bins) for matter density \\f$ C_l \\f$'s"]
            selection_num: ::std::os::raw::c_int::default(),
            #[doc = "< type of selection functions"]
            selection: selection_type::default(),
            #[doc = "< centers of selection functions"]
            selection_mean: [0.0; 100usize],
            #[doc = "< widths of selection functions"]
            selection_width: [0.0; 100usize],
            #[doc = "< in temperature calculation, do we want to include the intrinsic temperature + Sachs Wolfe term?"]
            switch_sw: ::std::os::raw::c_int::default(),
            #[doc = "< in temperature calculation, do we want to include the early integrated Sachs Wolfe term?"]
            switch_eisw: ::std::os::raw::c_int::default(),
            #[doc = "< in temperature calculation, do we want to include the late integrated Sachs Wolfe term?"]
            switch_lisw: ::std::os::raw::c_int::default(),
            #[doc = "< in temperature calculation, do we want to include the Doppler term?"]
            switch_dop: ::std::os::raw::c_int::default(),
            #[doc = "< in temperature calculation, do we want to include the polarization-related term?"]
            switch_pol: ::std::os::raw::c_int::default(),
            #[doc = "< at which redshift do we define the cut between eisw and lisw ?"]
            eisw_lisw_split_z: f64::default(),
            #[doc = "< Do we want to store perturbations?"]
            store_perturbations: ::std::os::raw::c_int::default(),
            #[doc = "< Number of perturbation outputs (default=0)"]
            k_output_values_num: ::std::os::raw::c_int::default(),
            #[doc = "< List of k values where perturbation output is requested."]
            k_output_values: [0.0; 30usize],
            #[doc = "< 3 x effective squared sound speed for the ultrarelativistic perturbations"]
            three_ceff2_ur: f64::default(),
            #[doc = "< 3 x effective viscosity parameter for the ultrarelativistic perturbations"]
            three_cvis2_ur: f64::default(),
            #[doc = "< when we compute only the matter spectrum / transfer functions, but not the CMB, we are sometimes interested to sample source functions at very high redshift, way before recombination. This z_max_pk will then fix the initial sampling time of the sources."]
            z_max_pk: f64::default(),
            #[doc = "< Angular contribution to collisional term at l>=2 for idm_fr-idr"]
            alpha_idm_dr: &mut f64::default(),
            #[doc = "< Angular contribution to collisional term at l>=2 for idr-idr"]
            beta_idr: &mut f64::default(),
            #[doc = "< Nature of the interacting dark radiation (free streaming or fluid)"]
            idr_nature: ::std::os::raw::c_int::default(),
            #[doc = "< do we need CMB-related sources (temperature, polarization) ?"]
            has_cmb: ::std::os::raw::c_short::default(),
            #[doc = "< do we need LSS-related sources (lensing potential, ...) ?"]
            has_lss: ::std::os::raw::c_short::default(),
            #[doc = "< do we have idm-dr interactions?"]
            has_idm_dr: ::std::os::raw::c_short::default(),
            #[doc = "< do we need to consider the dark matter sound speed in interaction models?"]
            has_idm_soundspeed: ::std::os::raw::c_short::default(),
            #[doc = "< gauge in which to perform this calculation"]
            gauge: possible_gauges::default(),
            #[doc = "< index value for scalars"]
            index_md_scalars: ::std::os::raw::c_int::default(),
            #[doc = "< index value for tensors"]
            index_md_tensors: ::std::os::raw::c_int::default(),
            #[doc = "< index value for vectors"]
            index_md_vectors: ::std::os::raw::c_int::default(),
            #[doc = "< number of modes included in computation"]
            md_size: ::std::os::raw::c_int::default(),
            #[doc = "< index value for adiabatic"]
            index_ic_ad: ::std::os::raw::c_int::default(),
            #[doc = "< index value for CDM isocurvature"]
            index_ic_cdi: ::std::os::raw::c_int::default(),
            #[doc = "< index value for baryon isocurvature"]
            index_ic_bi: ::std::os::raw::c_int::default(),
            #[doc = "< index value for neutrino density isocurvature"]
            index_ic_nid: ::std::os::raw::c_int::default(),
            #[doc = "< index value for neutrino velocity isocurvature"]
            index_ic_niv: ::std::os::raw::c_int::default(),
            #[doc = "< index value for unique possibility for tensors"]
            index_ic_ten: ::std::os::raw::c_int::default(),
            #[doc = "< for a given mode, ic_size[index_md] = number of initial conditions included in computation"]
            ic_size: &mut ::std::os::raw::c_int::default(),
            #[doc = "< do we need source for CMB temperature?"]
            has_source_t: ::std::os::raw::c_short::default(),
            #[doc = "< do we need source for CMB polarization?"]
            has_source_p: ::std::os::raw::c_short::default(),
            #[doc = "< do we need source for delta of total matter?"]
            has_source_delta_m: ::std::os::raw::c_short::default(),
            #[doc = "< do we ALSO need source for delta of ONLY cdm and baryon?"]
            has_source_delta_cb: ::std::os::raw::c_short::default(),
            #[doc = "< do we need source for delta total?"]
            has_source_delta_tot: ::std::os::raw::c_short::default(),
            #[doc = "< do we need source for delta of gammas?"]
            has_source_delta_g: ::std::os::raw::c_short::default(),
            #[doc = "< do we need source for delta of baryons?"]
            has_source_delta_b: ::std::os::raw::c_short::default(),
            #[doc = "< do we need source for delta of cold dark matter?"]
            has_source_delta_cdm: ::std::os::raw::c_short::default(),
            #[doc = "< do we need source for delta of interacting dark matter"]
            has_source_delta_idm: ::std::os::raw::c_short::default(),
            #[doc = "< do we need source for delta of interacting dark radiation?"]
            has_source_delta_idr: ::std::os::raw::c_short::default(),
            #[doc = "< do we need source for delta of DCDM?"]
            has_source_delta_dcdm: ::std::os::raw::c_short::default(),
            #[doc = "< do we need source for delta of dark energy?"]
            has_source_delta_fld: ::std::os::raw::c_short::default(),
            #[doc = "< do we need source for delta from scalar field?"]
            has_source_delta_scf: ::std::os::raw::c_short::default(),
            #[doc = "< do we need source for delta of decay radiation?"]
            has_source_delta_dr: ::std::os::raw::c_short::default(),
            #[doc = "< do we need source for delta of ultra-relativistic neutrinos/relics?"]
            has_source_delta_ur: ::std::os::raw::c_short::default(),
            #[doc = "< do we need source for delta of all non-cold dark matter species (e.g. massive neutrinos)?"]
            has_source_delta_ncdm: ::std::os::raw::c_short::default(),
            #[doc = "< do we need source for theta of total matter?"]
            has_source_theta_m: ::std::os::raw::c_short::default(),
            #[doc = "< do we ALSO need source for theta of ONLY cdm and baryon?"]
            has_source_theta_cb: ::std::os::raw::c_short::default(),
            #[doc = "< do we need source for theta total?"]
            has_source_theta_tot: ::std::os::raw::c_short::default(),
            #[doc = "< do we need source for theta of gammas?"]
            has_source_theta_g: ::std::os::raw::c_short::default(),
            #[doc = "< do we need source for theta of baryons?"]
            has_source_theta_b: ::std::os::raw::c_short::default(),
            #[doc = "< do we need source for theta of cold dark matter?"]
            has_source_theta_cdm: ::std::os::raw::c_short::default(),
            #[doc = "< do we need source for theta of interacting dark matter"]
            has_source_theta_idm: ::std::os::raw::c_short::default(),
            #[doc = "< do we need source for theta of interacting dark radiation?"]
            has_source_theta_idr: ::std::os::raw::c_short::default(),
            #[doc = "< do we need source for theta of DCDM?"]
            has_source_theta_dcdm: ::std::os::raw::c_short::default(),
            #[doc = "< do we need source for theta of dark energy?"]
            has_source_theta_fld: ::std::os::raw::c_short::default(),
            #[doc = "< do we need source for theta of scalar field?"]
            has_source_theta_scf: ::std::os::raw::c_short::default(),
            #[doc = "< do we need source for theta of ultra-relativistic neutrinos/relics?"]
            has_source_theta_dr: ::std::os::raw::c_short::default(),
            #[doc = "< do we need source for theta of ultra-relativistic neutrinos/relics?"]
            has_source_theta_ur: ::std::os::raw::c_short::default(),
            #[doc = "< do we need source for theta of all non-cold dark matter species (e.g. massive neutrinos)?"]
            has_source_theta_ncdm: ::std::os::raw::c_short::default(),
            #[doc = "< do we need source for metric fluctuation phi?"]
            has_source_phi: ::std::os::raw::c_short::default(),
            #[doc = "< do we need source for metric fluctuation phi'?"]
            has_source_phi_prime: ::std::os::raw::c_short::default(),
            #[doc = "< do we need source for metric fluctuation (phi+psi)?"]
            has_source_phi_plus_psi: ::std::os::raw::c_short::default(),
            #[doc = "< do we need source for metric fluctuation psi?"]
            has_source_psi: ::std::os::raw::c_short::default(),
            #[doc = "< do we need source for metric fluctuation h?"]
            has_source_h: ::std::os::raw::c_short::default(),
            #[doc = "< do we need source for metric fluctuation h'?"]
            has_source_h_prime: ::std::os::raw::c_short::default(),
            #[doc = "< do we need source for metric fluctuation eta?"]
            has_source_eta: ::std::os::raw::c_short::default(),
            #[doc = "< do we need source for metric fluctuation eta'?"]
            has_source_eta_prime: ::std::os::raw::c_short::default(),
            #[doc = "< do we need source for metric fluctuation H_T_Nb'?"]
            has_source_H_T_Nb_prime: ::std::os::raw::c_short::default(),
            #[doc = "< do we need source for metric fluctuation gamma in Nbody gauge?"]
            has_source_k2gamma_Nb: ::std::os::raw::c_short::default(),
            #[doc = "< index value for temperature (j=0 term)"]
            index_tp_t0: ::std::os::raw::c_int::default(),
            #[doc = "< index value for temperature (j=1 term)"]
            index_tp_t1: ::std::os::raw::c_int::default(),
            #[doc = "< index value for temperature (j=2 term)"]
            index_tp_t2: ::std::os::raw::c_int::default(),
            #[doc = "< index value for polarization"]
            index_tp_p: ::std::os::raw::c_int::default(),
            #[doc = "< index value for matter density fluctuation"]
            index_tp_delta_m: ::std::os::raw::c_int::default(),
            #[doc = "< index value for delta cb"]
            index_tp_delta_cb: ::std::os::raw::c_int::default(),
            #[doc = "< index value for total density fluctuation"]
            index_tp_delta_tot: ::std::os::raw::c_int::default(),
            #[doc = "< index value for delta of gammas"]
            index_tp_delta_g: ::std::os::raw::c_int::default(),
            #[doc = "< index value for delta of baryons"]
            index_tp_delta_b: ::std::os::raw::c_int::default(),
            #[doc = "< index value for delta of cold dark matter"]
            index_tp_delta_cdm: ::std::os::raw::c_int::default(),
            #[doc = "< index value for delta of interacting dark matter"]
            index_tp_delta_idm: ::std::os::raw::c_int::default(),
            #[doc = "< index value for delta of DCDM"]
            index_tp_delta_dcdm: ::std::os::raw::c_int::default(),
            #[doc = "< index value for delta of dark energy"]
            index_tp_delta_fld: ::std::os::raw::c_int::default(),
            #[doc = "< index value for delta of scalar field"]
            index_tp_delta_scf: ::std::os::raw::c_int::default(),
            #[doc = "< index value for delta of decay radiation"]
            index_tp_delta_dr: ::std::os::raw::c_int::default(),
            #[doc = "< index value for delta of ultra-relativistic neutrinos/relics"]
            index_tp_delta_ur: ::std::os::raw::c_int::default(),
            #[doc = "< index value for delta of interacting dark radiation"]
            index_tp_delta_idr: ::std::os::raw::c_int::default(),
            #[doc = "< index value for delta of first non-cold dark matter species (e.g. massive neutrinos)"]
            index_tp_delta_ncdm1: ::std::os::raw::c_int::default(),
            #[doc = "< Gas temperature perturbation"]
            index_tp_perturbed_recombination_delta_temp: ::std::os::raw::c_int::default(),
            #[doc = "< Inionization fraction perturbation"]
            index_tp_perturbed_recombination_delta_chi: ::std::os::raw::c_int::default(),
            #[doc = "< index value for matter velocity fluctuation"]
            index_tp_theta_m: ::std::os::raw::c_int::default(),
            #[doc = "< index value for theta cb"]
            index_tp_theta_cb: ::std::os::raw::c_int::default(),
            #[doc = "< index value for total velocity fluctuation"]
            index_tp_theta_tot: ::std::os::raw::c_int::default(),
            #[doc = "< index value for theta of gammas"]
            index_tp_theta_g: ::std::os::raw::c_int::default(),
            #[doc = "< index value for theta of baryons"]
            index_tp_theta_b: ::std::os::raw::c_int::default(),
            #[doc = "< index value for theta of cold dark matter"]
            index_tp_theta_cdm: ::std::os::raw::c_int::default(),
            #[doc = "< index value for theta of DCDM"]
            index_tp_theta_dcdm: ::std::os::raw::c_int::default(),
            #[doc = "< index value for theta of dark energy"]
            index_tp_theta_fld: ::std::os::raw::c_int::default(),
            #[doc = "< index value for theta of scalar field"]
            index_tp_theta_scf: ::std::os::raw::c_int::default(),
            #[doc = "< index value for theta of ultra-relativistic neutrinos/relics"]
            index_tp_theta_ur: ::std::os::raw::c_int::default(),
            #[doc = "< index value for theta of interacting dark radiation"]
            index_tp_theta_idr: ::std::os::raw::c_int::default(),
            #[doc = "< index value for theta of interacting dark matter"]
            index_tp_theta_idm: ::std::os::raw::c_int::default(),
            #[doc = "< index value for F1 of decay radiation"]
            index_tp_theta_dr: ::std::os::raw::c_int::default(),
            #[doc = "< index value for theta of first non-cold dark matter species (e.g. massive neutrinos)"]
            index_tp_theta_ncdm1: ::std::os::raw::c_int::default(),
            #[doc = "< index value for metric fluctuation phi"]
            index_tp_phi: ::std::os::raw::c_int::default(),
            #[doc = "< index value for metric fluctuation phi'"]
            index_tp_phi_prime: ::std::os::raw::c_int::default(),
            #[doc = "< index value for metric fluctuation phi+psi"]
            index_tp_phi_plus_psi: ::std::os::raw::c_int::default(),
            #[doc = "< index value for metric fluctuation psi"]
            index_tp_psi: ::std::os::raw::c_int::default(),
            #[doc = "< index value for metric fluctuation h"]
            index_tp_h: ::std::os::raw::c_int::default(),
            #[doc = "< index value for metric fluctuation h'"]
            index_tp_h_prime: ::std::os::raw::c_int::default(),
            #[doc = "< index value for metric fluctuation eta"]
            index_tp_eta: ::std::os::raw::c_int::default(),
            #[doc = "< index value for metric fluctuation eta'"]
            index_tp_eta_prime: ::std::os::raw::c_int::default(),
            #[doc = "< index value for metric fluctuation H_T_Nb'"]
            index_tp_H_T_Nb_prime: ::std::os::raw::c_int::default(),
            #[doc = "< index value for metric fluctuation gamma times k^2 in Nbody gauge"]
            index_tp_k2gamma_Nb: ::std::os::raw::c_int::default(),
            #[doc = "< number of types tp_size[index_md] included in computation for each mode"]
            tp_size: &mut ::std::os::raw::c_int::default(),
            #[doc = "< k_size_cmb[index_md] number of k values used"]
            #[doc = "for CMB calculations, requiring a fine"]
            #[doc = "sampling in k-space"]
            k_size_cmb: &mut ::std::os::raw::c_int::default(),
            #[doc = "< k_size_cl[index_md] number of k values used"]
            #[doc = "for non-CMB \\f$ C_l \\f$ calculations, requiring a coarse"]
            #[doc = "sampling in k-space."]
            k_size_cl: &mut ::std::os::raw::c_int::default(),
            #[doc = "< number of k values for the P(k,z) and T(k,z) output, not including possible additional values for non-linear corrections"]
            k_size_pk: ::std::os::raw::c_int::default(),
            #[doc = "< k_size[index_md] = total number of k values,"]
            #[doc = "including those needed for all C_l, P(k),"]
            #[doc = "nonlinear corrections"]
            k_size: &mut ::std::os::raw::c_int::default(),
            #[doc = "< k[index_md][index_k] = list of values"]
            k: return_raw_dbl_ptr::<f64>(),
            #[doc = "< minimum value (over all modes)"]
            k_min: f64::default(),
            #[doc = "< maximum value (over all modes)"]
            k_max: f64::default(),
            #[doc = "< array of tau values"]
            tau_sampling: &mut f64::default(),
            #[doc = "< number of values in this array"]
            tau_size: ::std::os::raw::c_int::default(),
            #[doc = "< used in presence of selection functions (for matter density, cosmic shear...)"]
            selection_min_of_tau_min: f64::default(),
            #[doc = "< used in presence of selection functions (for matter density, cosmic shear...)"]
            selection_max_of_tau_max: f64::default(),
            #[doc = "< used in presence of selection functions (for matter density, cosmic shear...)"]
            selection_delta_tau: f64::default(),
            #[doc = "< value of conformal time below which W(tau) is considered to vanish for each bin"]
            selection_tau_min: &mut f64::default(),
            #[doc = "< value of conformal time above which W(tau) is considered to vanish for each bin"]
            selection_tau_max: &mut f64::default(),
            #[doc = "< value of conformal time at the center of each bin"]
            selection_tau: &mut f64::default(),
            #[doc = "< selection function W(tau), normalized to \\f$ \\int W(tau) dtau=1 \\f$, stored in selection_function[bin*ppt->tau_size+index_tau]"]
            selection_function: &mut f64::default(),
            #[doc = "< Pointer towards the source interpolation table"]
            #[doc = "sources[index_md]"]
            #[doc = "[index_ic * ppt->tp_size[index_md] + index_tp]"]
            #[doc = "[index_tau * ppt->k_size + index_k]"]
            sources: return_raw_tpl_ptr::<f64>(),
            #[doc = "< log of the arrau tau_sampling, covering only the"]
            #[doc = "final time range required for the output of"]
            #[doc = "Fourier transfer functions (used for interpolations)"]
            ln_tau: &mut f64::default(),
            #[doc = "< total number of values in this array"]
            ln_tau_size: ::std::os::raw::c_int::default(),
            #[doc = "< first index relevant for output of P(k,z) and T(k,z)"]
            index_ln_tau_pk: ::std::os::raw::c_int::default(),
            #[doc = "< Pointer towards the source interpolation table"]
            #[doc = "late_sources[index_md]"]
            #[doc = "[index_ic * ppt->tp_size[index_md] + index_tp]"]
            #[doc = "[index_tau * ppt->k_size + index_k]"]
            #[doc = "Note that this is not a replication of part of the sources table,"]
            #[doc = "it is just poiting towards the same memory zone, at the place where the late_sources actually start"]
            late_sources: return_raw_tpl_ptr(),
            #[doc = "< Pointer towards the splined source interpolation table with second derivatives with respect to time"]
            #[doc = "ddlate_sources[index_md]"]
            #[doc = "[index_ic * ppt->tp_size[index_md] + index_tp]"]
            #[doc = "[index_tau * ppt->k_size + index_k]"]
            ddlate_sources: return_raw_tpl_ptr(),
            #[doc = "< List of indices corresponding to k-values close to k_output_values for each mode."]
            #[doc = "index_k_output_values[index_md*k_output_values_num+ik]"]
            index_k_output_values: &mut ::std::os::raw::c_int::default(),
            #[doc = "< _DELIMITER_ separated string of titles for scalar perturbation output files."]
            scalar_titles: [i8::default(); 8000],
            #[doc = "< _DELIMITER_ separated string of titles for vector perturbation output files."]
            vector_titles: [i8::default(); 8000],
            #[doc = "< _DELIMITER_ separated string of titles for tensor perturbation output files."]
            tensor_titles: [i8::default(); 8000],
            #[doc = "< number of titles/columns in scalar perturbation output files"]
            number_of_scalar_titles: ::std::os::raw::c_int::default(),
            #[doc = "< number of titles/columns in vector perturbation output files"]
            number_of_vector_titles: ::std::os::raw::c_int::default(),
            #[doc = "< number of titles/columns in tensor perturbation output files"]
            number_of_tensor_titles: ::std::os::raw::c_int::default(),
            #[doc = "< Array of double pointers to perturbation output for scalars"]
            scalar_perturbations_data: [&mut f64::default(); 30usize],
            #[doc = "< Array of double pointers to perturbation output for vectors"]
            vector_perturbations_data: [&mut f64::default(); 30usize],
            #[doc = "< Array of double pointers to perturbation output for tensors"]
            tensor_perturbations_data: [&mut f64::default(); 30usize],
            #[doc = "< Array of sizes of scalar double pointers"]
            size_scalar_perturbation_data: [::std::os::raw::c_int::default(); 30usize],
            #[doc = "< Array of sizes of vector double pointers"]
            size_vector_perturbation_data: [::std::os::raw::c_int::default(); 30usize],
            #[doc = "< Array of sizes of tensor double pointers"]
            size_tensor_perturbation_data: [::std::os::raw::c_int::default(); 30usize],
            #[doc = "< flag regulating the amount of information sent to standard output (none if set to zero)"]
            perturbations_verbose: ::std::os::raw::c_short::default(),
            #[doc = "< zone for writing error messages"]
            error_message: zero_padded_i8_array::<2048>(""),
        }
    }
}