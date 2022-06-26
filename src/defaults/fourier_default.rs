use crate::{fourier, non_linear_method, source_extrapolation, hmcode_baryonic_feedback_model};
use crate::utils::*;

impl Default for fourier {
    fn default() -> fourier {
        fourier {
            #[doc = "< method for computing non-linear corrections (none, Halogit, etc.)"]
            method: non_linear_method::default(),
            #[doc = "< method for analytical extrapolation of sources beyond pre-computed range"]
            extrapolation_method: source_extrapolation::default(),
            feedback: hmcode_baryonic_feedback_model::default(),
            #[doc = " to choose between different baryonic feedback models"]
            #[doc = "in hmcode (dmonly, gas cooling, Agn or supernova feedback)"]
            c_min: f64::default(),
            #[doc = " for HMcode: minimum concentration in Bullock 2001 mass-concentration relation"]
            eta_0: f64::default(),
            #[doc = " for HMcode: halo bloating parameter"]
            z_infinity: f64::default(),
            #[doc = "< flag: in case wa_fld is defined and non-zero, should we use the pk_eq method?"]
            has_pk_eq: ::std::os::raw::c_short::default(),
            #[doc = "< set equal to phr->index_md_scalars"]
            #[doc = "(useful since this module only deals with"]
            #[doc = "scalars)"]
            index_md_scalars: ::std::os::raw::c_int::default(),
            #[doc = "< for a given mode, ic_size[index_md] = number of initial conditions included in computation"]
            ic_size: ::std::os::raw::c_int::default(),
            #[doc = "< for a given mode, ic_ic_size[index_md] = number of pairs of (index_ic1, index_ic2) with index_ic2 >= index_ic1; this number is just N(N+1)/2  where N = ic_size[index_md]"]
            ic_ic_size: ::std::os::raw::c_int::default(),
            #[doc = "< for a given mode, is_non_zero[index_md][index_ic1_ic2] is set to true if the pair of initial conditions (index_ic1, index_ic2) are statistically correlated, or to false if they are uncorrelated"]
            is_non_zero: &mut ::std::os::raw::c_short::default(),
            #[doc = "< do we want spectra for total matter?"]
            has_pk_m: ::std::os::raw::c_short::default(),
            #[doc = "< do we want spectra for cdm+baryons?"]
            has_pk_cb: ::std::os::raw::c_short::default(),
            #[doc = "< index of pk for matter (defined only when has_pk_m is TRUE)"]
            index_pk_m: ::std::os::raw::c_int::default(),
            #[doc = "< index of pk for cold dark matter plus baryons (defined only when has_pk_cb is TRUE"]
            index_pk_cb: ::std::os::raw::c_int::default(),
            #[doc = "< always equal to index_pk_m"]
            #[doc = "(always defined, useful e.g. for weak lensing spectrum)"]
            index_pk_total: ::std::os::raw::c_int::default(),
            #[doc = "< equal to index_pk_cb if it exists, otherwise to index_pk_m"]
            #[doc = "(always defined, useful e.g. for galaxy clustering spectrum)"]
            index_pk_cluster: ::std::os::raw::c_int::default(),
            #[doc = "< k_size = total number of pk"]
            pk_size: ::std::os::raw::c_int::default(),
            #[doc = "< do we need matter Fourier spectrum?"]
            has_pk_matter: ::std::os::raw::c_short::default(),
            #[doc = "< k_size = total number of k values"]
            k_size: ::std::os::raw::c_int::default(),
            #[doc = "< k_size = number of k values for P(k,z) and T(k,z) output)"]
            k_size_pk: ::std::os::raw::c_int::default(),
            #[doc = "< k[index_k] = list of k values"]
            k: &mut f64::default(),
            #[doc = "< ln_k[index_k] = list of log(k) values"]
            ln_k: &mut f64::default(),
            #[doc = "< log(tau) array, only needed if user wants"]
            #[doc = "some output at z>0, instead of only z=0.  This"]
            #[doc = "array only covers late times, used for the"]
            #[doc = "output of P(k) or T(k), and matching the"]
            #[doc = "condition z(tau) < z_max_pk"]
            ln_tau: &mut f64::default(),
            #[doc = "< total number of values in this array"]
            ln_tau_size: ::std::os::raw::c_int::default(),
            #[doc = "< first index relevant for output of P(k,z) and T(k,z)"]
            index_ln_tau_pk: ::std::os::raw::c_int::default(),
            #[doc = "< Matter power spectrum (linear)."]
            #[doc = "Depends on indices index_pk, index_ic1_ic2, index_k, index_tau as:"]
            #[doc = "ln_pk_ic_l[index_pk][(index_tau * pfo->k_size + index_k)* pfo->ic_ic_size + index_ic1_ic2]"]
            #[doc = "where index-pk labels P(k) types (m = total matter, cb = baryons+CDM),"]
            #[doc = "while index_ic1_ic2 labels ordered pairs (index_ic1, index_ic2) (since"]
            #[doc = "the primordial spectrum is symmetric in (index_ic1, index_ic2))."]
            #[doc = "- for diagonal elements (index_ic1 = index_ic2) this arrays contains"]
            #[doc = "ln[P(k)] where P(k) is positive by construction."]
            #[doc = "- for non-diagonal elements this arrays contains the k-dependent"]
            #[doc = "cosine of the correlation angle, namely"]
            #[doc = "P(k)_(index_ic1, index_ic2)/sqrt[P(k)_index_ic1 P(k)_index_ic2]"]
            #[doc = "This choice is convenient since the sign of the non-diagonal cross-correlation"]
            #[doc = "could be negative. For fully correlated or anti-correlated initial conditions,"]
            #[doc = "this non-diagonal element is independent on k, and equal to +1 or -1."]
            ln_pk_ic_l: return_raw_dbl_ptr(),
            #[doc = "< second derivative of above array with respect to log(tau), for spline interpolation. So:"]
            #[doc = "- for index_ic1 = index_ic, we spline ln[P(k)] vs. ln(k), which is"]
            #[doc = "good since this function is usually smooth."]
            #[doc = "- for non-diagonal coefficients, we spline"]
            #[doc = "P(k)_(index_ic1, index_ic2)/sqrt[P(k)_index_ic1 P(k)_index_ic2]"]
            #[doc = "vs. ln(k), which is fine since this quantity is often assumed to be"]
            #[doc = "constant (e.g for fully correlated/anticorrelated initial conditions)"]
            #[doc = "or nearly constant, and with arbitrary sign."]
            ddln_pk_ic_l: return_raw_dbl_ptr(),
            #[doc = "< Total matter power spectrum summed over initial conditions (linear)."]
            #[doc = "Only depends on indices index_pk,index_k, index_tau as:"]
            #[doc = "ln_pk[index_pk][index_tau * pfo->k_size + index_k]"]
            ln_pk_l: return_raw_dbl_ptr(),
            #[doc = "< second derivative of above array with respect to log(tau), for spline interpolation."]
            ddln_pk_l: return_raw_dbl_ptr(),
            #[doc = "< Total matter power spectrum summed over initial conditions (nonlinear)."]
            #[doc = "Only depends on indices index_pk,index_k, index_tau as:"]
            #[doc = "ln_pk[index_pk][index_tau * pfo->k_size + index_k]"]
            ln_pk_nl: return_raw_dbl_ptr(),
            #[doc = "< second derivative of above array with respect to log(tau), for spline interpolation."]
            ddln_pk_nl: return_raw_dbl_ptr(),
            #[doc = "< sigma8[index_pk]"]
            sigma8: &mut f64::default(),
            k_size_extra: ::std::os::raw::c_int::default(),
            #[doc = "< tau_size = number of values"]
            tau_size: ::std::os::raw::c_int::default(),
            #[doc = "< tau[index_tau] = list of time values, covering"]
            #[doc = "all the values of the perturbation module"]
            tau: &mut f64::default(),
            #[doc = "< nl_corr_density[index_pk][index_tau * ppt->k_size + index_k]"]
            nl_corr_density: return_raw_dbl_ptr(),
            #[doc = "< wavenumber at which non-linear corrections become important,"]
            #[doc = "defined differently by different non_linear_method's"]
            k_nl: return_raw_dbl_ptr(),
            #[doc = "< index of smallest value of tau at which nonlinear corrections have been computed"]
            #[doc = "(so, for tau<tau_min_nl, the array nl_corr_density only contains some factors 1"]
            index_tau_min_nl: ::std::os::raw::c_int::default(),
            #[doc = "< index of w in table pk_eq_w_and_Omega"]
            index_pk_eq_w: ::std::os::raw::c_int::default(),
            #[doc = "< index of Omega_m in table pk_eq_w_and_Omega"]
            index_pk_eq_Omega_m: ::std::os::raw::c_int::default(),
            #[doc = "< number of indices in table pk_eq_w_and_Omega"]
            pk_eq_size: ::std::os::raw::c_int::default(),
            #[doc = "< number of times (and raws in table pk_eq_w_and_Omega)"]
            pk_eq_tau_size: ::std::os::raw::c_int::default(),
            #[doc = "< table of time values"]
            pk_eq_tau: &mut f64::default(),
            #[doc = "< table of background quantites"]
            pk_eq_w_and_Omega: &mut f64::default(),
            #[doc = "< table of second derivatives"]
            pk_eq_ddw_and_ddOmega: &mut f64::default(),
            #[doc = "< amount of information written in standard output"]
            fourier_verbose: ::std::os::raw::c_short::default(),
            #[doc = "< zone for writing error messages"]
            error_message: zero_padded_i8_array(""),
        }
    }
}