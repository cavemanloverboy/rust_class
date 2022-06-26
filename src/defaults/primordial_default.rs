use crate::{primordial, primordial_spectrum_type, potential_shape, phi_pivot_methods, inflation_module_behavior};
use crate::utils::*;

impl Default for primordial {
    fn default() -> Self {
        primordial {
            #[doc = "< pivot scale in \\f$ Mpc^{-1} \\f$"]
            k_pivot: f64::default(),
            has_k_max_for_primordial_pk: ::std::os::raw::c_int::default(),
            #[doc = "< maximum value of k in 1/Mpc in P(k)"]
            k_max_for_primordial_pk: f64::default(),
            #[doc = "< type of primordial spectrum (simple analytic from, integration of inflationary perturbations, etc.)"]
            primordial_spec_type: primordial_spectrum_type::default(),
            #[doc = "< usual scalar amplitude = curvature power spectrum at pivot scale"]
            A_s: f64::default(),
            #[doc = "< usual scalar tilt = [curvature power spectrum tilt at pivot scale -1]"]
            n_s: f64::default(),
            #[doc = "< usual scalar running"]
            alpha_s: f64::default(),
            #[doc = "< running of running"]
            beta_s: f64::default(),
            #[doc = "< usual tensor to scalar ratio of power spectra, \\f$ r=A_T/A_S=P_h/P_R \\f$"]
            r: f64::default(),
            #[doc = "< usual tensor tilt = [GW power spectrum tilt at pivot scale]"]
            n_t: f64::default(),
            #[doc = "< usual tensor running"]
            alpha_t: f64::default(),
            #[doc = "< baryon isocurvature (BI) entropy-to-curvature ratio \\f$ S_{bi}/R \\f$"]
            f_bi: f64::default(),
            #[doc = "< BI tilt"]
            n_bi: f64::default(),
            #[doc = "< BI running"]
            alpha_bi: f64::default(),
            #[doc = "< CDM isocurvature (CDI) entropy-to-curvature ratio \\f$ S_{cdi}/R \\f$"]
            f_cdi: f64::default(),
            #[doc = "< CDI tilt"]
            n_cdi: f64::default(),
            #[doc = "< CDI running"]
            alpha_cdi: f64::default(),
            #[doc = "< neutrino density isocurvature (NID) entropy-to-curvature ratio \\f$ S_{nid}/R \\f$"]
            f_nid: f64::default(),
            #[doc = "< NID tilt"]
            n_nid: f64::default(),
            #[doc = "< NID running"]
            alpha_nid: f64::default(),
            #[doc = "< neutrino velocity isocurvature (NIV) entropy-to-curvature ratio \\f$ S_{niv}/R \\f$"]
            f_niv: f64::default(),
            #[doc = "< NIV tilt"]
            n_niv: f64::default(),
            #[doc = "< NIV running"]
            alpha_niv: f64::default(),
            #[doc = "< ADxBI cross-correlation at pivot scale, from -1 to 1"]
            c_ad_bi: f64::default(),
            #[doc = "< ADxBI cross-correlation tilt"]
            n_ad_bi: f64::default(),
            #[doc = "< ADxBI cross-correlation running"]
            alpha_ad_bi: f64::default(),
            #[doc = "< ADxCDI cross-correlation at pivot scale, from -1 to 1"]
            c_ad_cdi: f64::default(),
            #[doc = "< ADxCDI cross-correlation tilt"]
            n_ad_cdi: f64::default(),
            #[doc = "< ADxCDI cross-correlation running"]
            alpha_ad_cdi: f64::default(),
            #[doc = "< ADxNID cross-correlation at pivot scale, from -1 to 1"]
            c_ad_nid: f64::default(),
            #[doc = "< ADxNID cross-correlation tilt"]
            n_ad_nid: f64::default(),
            #[doc = "< ADxNID cross-correlation running"]
            alpha_ad_nid: f64::default(),
            #[doc = "< ADxNIV cross-correlation at pivot scale, from -1 to 1"]
            c_ad_niv: f64::default(),
            #[doc = "< ADxNIV cross-correlation tilt"]
            n_ad_niv: f64::default(),
            #[doc = "< ADxNIV cross-correlation running"]
            alpha_ad_niv: f64::default(),
            #[doc = "< BIxCDI cross-correlation at pivot scale, from -1 to 1"]
            c_bi_cdi: f64::default(),
            #[doc = "< BIxCDI cross-correlation tilt"]
            n_bi_cdi: f64::default(),
            #[doc = "< BIxCDI cross-correlation running"]
            alpha_bi_cdi: f64::default(),
            #[doc = "< BIxNIV cross-correlation at pivot scale, from -1 to 1"]
            c_bi_nid: f64::default(),
            #[doc = "< BIxNIV cross-correlation tilt"]
            n_bi_nid: f64::default(),
            #[doc = "< BIxNIV cross-correlation running"]
            alpha_bi_nid: f64::default(),
            #[doc = "< BIxNIV cross-correlation at pivot scale, from -1 to 1"]
            c_bi_niv: f64::default(),
            #[doc = "< BIxNIV cross-correlation tilt"]
            n_bi_niv: f64::default(),
            #[doc = "< BIxNIV cross-correlation running"]
            alpha_bi_niv: f64::default(),
            #[doc = "< CDIxNID cross-correlation at pivot scale, from -1 to 1"]
            c_cdi_nid: f64::default(),
            #[doc = "< CDIxNID cross-correlation tilt"]
            n_cdi_nid: f64::default(),
            #[doc = "< CDIxNID cross-correlation running"]
            alpha_cdi_nid: f64::default(),
            #[doc = "< CDIxNIV cross-correlation at pivot scale, from -1 to 1"]
            c_cdi_niv: f64::default(),
            #[doc = "< CDIxNIV cross-correlation tilt"]
            n_cdi_niv: f64::default(),
            #[doc = "< CDIxNIV cross-correlation running"]
            alpha_cdi_niv: f64::default(),
            #[doc = "< NIDxNIV cross-correlation at pivot scale, from -1 to 1"]
            c_nid_niv: f64::default(),
            #[doc = "< NIDxNIV cross-correlation tilt"]
            n_nid_niv: f64::default(),
            #[doc = "< NIDxNIV cross-correlation running"]
            alpha_nid_niv: f64::default(),
            #[doc = " parameters describing the case primordial_spec_type = inflation_V"]
            potential: potential_shape::default(),
            #[doc = "< one parameter of the function V(phi)"]
            V0: f64::default(),
            #[doc = "< one parameter of the function V(phi)"]
            V1: f64::default(),
            #[doc = "< one parameter of the function V(phi)"]
            V2: f64::default(),
            #[doc = "< one parameter of the function V(phi)"]
            V3: f64::default(),
            #[doc = "< one parameter of the function V(phi)"]
            V4: f64::default(),
            #[doc = "< one parameter of the function H(phi)"]
            H0: f64::default(),
            #[doc = "< one parameter of the function H(phi)"]
            H1: f64::default(),
            #[doc = "< one parameter of the function H(phi)"]
            H2: f64::default(),
            #[doc = "< one parameter of the function H(phi)"]
            H3: f64::default(),
            #[doc = "< one parameter of the function H(phi)"]
            H4: f64::default(),
            #[doc = "< value of inflaton at the end of inflation"]
            phi_end: f64::default(),
            #[doc = "< flag for method used to define and find the pivot scale"]
            phi_pivot_method: phi_pivot_methods::default(),
            #[doc = "< For each of the above methods, critical value to be reached between pivot and end of inflation (N_star, [aH]ratio, etc.)"]
            phi_pivot_target: f64::default(),
            #[doc = "< Specifies if the inflation module computes the primordial spectrum numerically (default) or analytically"]
            behavior: inflation_module_behavior::default(),
            #[doc = "< string with the command for calling 'external_Pk'"]
            command: &mut ::std::os::raw::c_char::default(),
            #[doc = "< one parameter of the primordial computed in 'external_Pk'"]
            custom1: f64::default(),
            #[doc = "< one parameter of the primordial computed in 'external_Pk'"]
            custom2: f64::default(),
            #[doc = "< one parameter of the primordial computed in 'external_Pk'"]
            custom3: f64::default(),
            #[doc = "< one parameter of the primordial computed in 'external_Pk'"]
            custom4: f64::default(),
            #[doc = "< one parameter of the primordial computed in 'external_Pk'"]
            custom5: f64::default(),
            #[doc = "< one parameter of the primordial computed in 'external_Pk'"]
            custom6: f64::default(),
            #[doc = "< one parameter of the primordial computed in 'external_Pk'"]
            custom7: f64::default(),
            #[doc = "< one parameter of the primordial computed in 'external_Pk'"]
            custom8: f64::default(),
            #[doc = "< one parameter of the primordial computed in 'external_Pk'"]
            custom9: f64::default(),
            #[doc = "< one parameter of the primordial computed in 'external_Pk'"]
            custom10: f64::default(),
            #[doc = "< number of modes included in computation"]
            md_size: ::std::os::raw::c_int::default(),
            #[doc = "< for a given mode, ic_size[index_md] = number of initial conditions included in computation"]
            ic_size: &mut ::std::os::raw::c_int::default(),
            #[doc = "< number of ordered pairs of (index_ic1, index_ic2); this number is just N(N+1)/2  where N = ic_size[index_md]"]
            ic_ic_size: &mut ::std::os::raw::c_int::default(),
            #[doc = "< number of ln(k) values"]
            lnk_size: ::std::os::raw::c_int::default(),
            #[doc = "< list of ln(k) values lnk[index_k]"]
            lnk: &mut f64::default(),
            #[doc = "< depends on indices index_md, index_ic1, index_ic2, index_k as:"]
            #[doc = "lnpk[index_md][index_k*ppm->ic_ic_size[index_md]+index_ic1_ic2]"]
            #[doc = "where index_ic1_ic2 labels ordered pairs (index_ic1, index_ic2) (since"]
            #[doc = "the primordial spectrum is symmetric in (index_ic1, index_ic2))."]
            #[doc = "- for diagonal elements (index_ic1 = index_ic2) this arrays contains"]
            #[doc = "ln[P(k)] where P(k) is positive by construction."]
            #[doc = "- for non-diagonal elements this arrays contains the k-dependent"]
            #[doc = "cosine of the correlation angle, namely"]
            #[doc = "P(k )_(index_ic1, index_ic2)/sqrt[P(k)_index_ic1 P(k)_index_ic2]"]
            #[doc = "This choice is convenient since the sign of the non-diagonal cross-correlation"]
            #[doc = "is arbitrary. For fully correlated or anti-correlated initial conditions,"]
            #[doc = "this non -diagonal element is independent on k, and equal to +1 or -1."]
            lnpk: return_raw_dbl_ptr(),
            #[doc = "< second derivative of above array, for spline interpolation. So:"]
            #[doc = "- for index_ic1 = index_ic, we spline ln[P(k)] vs. ln(k), which is"]
            #[doc = "good since this function is usually smooth."]
            #[doc = "- for non-diagonal coefficients, we spline"]
            #[doc = "P(k)_(index_ic1, index_ic2)/sqrt[P(k)_index_ic1 P(k)_index_ic2]"]
            #[doc = "vs. ln(k), which is fine since this quantity is often assumed to be"]
            #[doc = "constant (e.g for fully correlated/anticorrelated initial conditions)"]
            #[doc = "or nearly constant, and with arbitrary sign."]
            ddlnpk: return_raw_dbl_ptr(),
            #[doc = "< is_non_zero[index_md][index_ic1_ic2] set to false if pair"]
            #[doc = "(index_ic1, index_ic2) is uncorrelated"]
            #[doc = "(ensures more precision and saves time with respect to the option"]
            #[doc = "of simply setting P(k)_(index_ic1, index_ic2) to zero)"]
            is_non_zero: return_raw_dbl_ptr(),
            #[doc = "< all amplitudes in matrix form: amplitude[index_md][index_ic1_ic2]"]
            amplitude: return_raw_dbl_ptr(),
            #[doc = "< all tilts in matrix form: tilt[index_md][index_ic1_ic2]"]
            tilt: return_raw_dbl_ptr(),
            #[doc = "< all runnings in matrix form: running[index_md][index_ic1_ic2]"]
            running: return_raw_dbl_ptr(),
            #[doc = "< scale factor"]
            index_in_a: ::std::os::raw::c_int::default(),
            #[doc = "< inflaton vev"]
            index_in_phi: ::std::os::raw::c_int::default(),
            #[doc = "< its time derivative"]
            index_in_dphi: ::std::os::raw::c_int::default(),
            #[doc = "< Mukhanov variable (real part)"]
            index_in_ksi_re: ::std::os::raw::c_int::default(),
            #[doc = "< Mukhanov variable (imaginary part)"]
            index_in_ksi_im: ::std::os::raw::c_int::default(),
            #[doc = "< Mukhanov variable (real part, time derivative)"]
            index_in_dksi_re: ::std::os::raw::c_int::default(),
            #[doc = "< Mukhanov variable (imaginary part, time derivative)"]
            index_in_dksi_im: ::std::os::raw::c_int::default(),
            #[doc = "< tensor perturbation (real part)"]
            index_in_ah_re: ::std::os::raw::c_int::default(),
            #[doc = "< tensor perturbation (imaginary part)"]
            index_in_ah_im: ::std::os::raw::c_int::default(),
            #[doc = "< tensor perturbation (real part, time derivative)"]
            index_in_dah_re: ::std::os::raw::c_int::default(),
            #[doc = "< tensor perturbation (imaginary part, time derivative)"]
            index_in_dah_im: ::std::os::raw::c_int::default(),
            #[doc = "< size of vector of background quantities only"]
            in_bg_size: ::std::os::raw::c_int::default(),
            #[doc = "< full size of vector"]
            in_size: ::std::os::raw::c_int::default(),
            #[doc = "< in inflationary module, value of"]
            #[doc = "phi_pivot (set to 0 for inflation_V,"]
            #[doc = "inflation_H; found by code for"]
            #[doc = "inflation_V_end)"]
            phi_pivot: f64::default(),
            #[doc = "< in inflationary module, value of phi when \\f$ k_{min}=aH \\f$"]
            phi_min: f64::default(),
            #[doc = "< in inflationary module, value of phi when \\f$ k_{max}=aH \\f$"]
            phi_max: f64::default(),
            #[doc = "< in inflationary module, value of phi at the end of inflation"]
            phi_stop: f64::default(),
            #[doc = "< flag regulating the amount of information sent to standard output (none if set to zero)"]
            primordial_verbose: ::std::os::raw::c_short::default(),
            #[doc = "< zone for writing error messages"]
            error_message: zero_padded_i8_array(""),
        }
    }
}