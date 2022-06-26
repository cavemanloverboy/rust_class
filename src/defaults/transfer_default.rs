use crate::transfer;
use crate::utils::*;

impl Default for transfer {
    fn default() -> transfer {
        transfer {
            #[doc = "< normally set to one, can be used"]
            #[doc = "exceptionally to rescale by hand the CMB"]
            #[doc = "lensing potential"]
            lcmb_rescale: f64::default(),
            #[doc = "< normally set to zero, can be used"]
            #[doc = "exceptionally to tilt by hand the CMB"]
            #[doc = "lensing potential"]
            lcmb_tilt: f64::default(),
            #[doc = "< if lcmb_tilt non-zero, corresponding pivot"]
            #[doc = "scale"]
            lcmb_pivot: f64::default(),
            #[doc = "< light-to-mass bias in the transfer function of density number count"]
            selection_bias: [f64::default(); 100usize],
            #[doc = "< magnification bias in the transfer function of density number count"]
            selection_magnification_bias: [f64::default(); 100usize],
            #[doc = "< Has dN/dz (selection function) input file?"]
            has_nz_file: ::std::os::raw::c_short::default(),
            #[doc = "< Use analytic form for dN/dz (selection function) distribution?"]
            has_nz_analytic: ::std::os::raw::c_short::default(),
            #[doc = "< dN/dz (selection function) input file name"]
            nz_file_name: zero_padded_i8_array(""),
            #[doc = "< number of redshift values in input tabulated selection function"]
            nz_size: ::std::os::raw::c_int::default(),
            #[doc = "< redshift values in input tabulated selection function"]
            nz_z: &mut f64::default(),
            #[doc = "< input tabulated values of selection function"]
            nz_nz: &mut f64::default(),
            #[doc = "< second derivatives in splined selection function"]
            nz_ddnz: &mut f64::default(),
            #[doc = "< Has dN/dz (evolution function) input file?"]
            has_nz_evo_file: ::std::os::raw::c_short::default(),
            #[doc = "< Use analytic form for dN/dz (evolution function) distribution?"]
            has_nz_evo_analytic: ::std::os::raw::c_short::default(),
            #[doc = "< dN/dz (evolution function) input file name"]
            nz_evo_file_name: zero_padded_i8_array(""),
            #[doc = "< number of redshift values in input tabulated evolution function"]
            nz_evo_size: ::std::os::raw::c_int::default(),
            #[doc = "< redshift values in input tabulated evolution function"]
            nz_evo_z: &mut f64::default(),
            #[doc = "< input tabulated values of evolution function"]
            nz_evo_nz: &mut f64::default(),
            #[doc = "< log of tabulated values of evolution function"]
            nz_evo_dlog_nz: &mut f64::default(),
            #[doc = "< second derivatives in splined log of evolution function"]
            nz_evo_dd_dlog_nz: &mut f64::default(),
            #[doc = "< copy of same flag in perturbation structure"]
            has_cls: ::std::os::raw::c_short::default(),
            #[doc = "< number of modes included in computation"]
            md_size: ::std::os::raw::c_int::default(),
            #[doc = "< index for transfer type = temperature (j=0 term)"]
            index_tt_t0: ::std::os::raw::c_int::default(),
            #[doc = "< index for transfer type = temperature (j=1 term)"]
            index_tt_t1: ::std::os::raw::c_int::default(),
            #[doc = "< index for transfer type = temperature (j=2 term)"]
            index_tt_t2: ::std::os::raw::c_int::default(),
            #[doc = "< index for transfer type = E-polarization"]
            index_tt_e: ::std::os::raw::c_int::default(),
            #[doc = "< index for transfer type = B-polarization"]
            index_tt_b: ::std::os::raw::c_int::default(),
            #[doc = "< index for transfer type = CMB lensing"]
            index_tt_lcmb: ::std::os::raw::c_int::default(),
            #[doc = "< index for first bin of transfer type = matter density"]
            index_tt_density: ::std::os::raw::c_int::default(),
            #[doc = "< index for first bin of transfer type = galaxy lensing"]
            index_tt_lensing: ::std::os::raw::c_int::default(),
            #[doc = "< index for first bin of transfer type = redshift space distortion of number count"]
            index_tt_rsd: ::std::os::raw::c_int::default(),
            #[doc = "< index for first bin of transfer type = doppler effect for of number count (j=0 term)"]
            index_tt_d0: ::std::os::raw::c_int::default(),
            #[doc = "< index for first bin of transfer type = doppler effect for of number count (j=1 term)"]
            index_tt_d1: ::std::os::raw::c_int::default(),
            #[doc = "< index for first bin of transfer type = lensing for of number count"]
            index_tt_nc_lens: ::std::os::raw::c_int::default(),
            #[doc = "< index for first bin of transfer type = gravity term G1 for of number count"]
            index_tt_nc_g1: ::std::os::raw::c_int::default(),
            #[doc = "< index for first bin of transfer type = gravity term G2 for of number count"]
            index_tt_nc_g2: ::std::os::raw::c_int::default(),
            #[doc = "< index for first bin of transfer type = gravity term G3 for of number count"]
            index_tt_nc_g3: ::std::os::raw::c_int::default(),
            #[doc = "< index for first bin of transfer type = gravity term G3 for of number count"]
            index_tt_nc_g4: ::std::os::raw::c_int::default(),
            #[doc = "< index for first bin of transfer type = gravity term G3 for of number count"]
            index_tt_nc_g5: ::std::os::raw::c_int::default(),
            #[doc = "< number of requested transfer types tt_size[index_md] for each mode"]
            tt_size: &mut ::std::os::raw::c_int::default(),
            #[doc = "< number of multipole values for which we effectively compute the transfer function,l_size_tt[index_md][index_tt]"]
            l_size_tt: return_raw_dbl_ptr(),
            #[doc = "< number of multipole values for each requested mode, l_size[index_md]"]
            l_size: &mut ::std::os::raw::c_int::default(),
            #[doc = "< greatest of all l_size[index_md]"]
            l_size_max: ::std::os::raw::c_int::default(),
            #[doc = "< list of multipole values l[index_l]"]
            l: &mut ::std::os::raw::c_int::default(),
            #[doc = "< correction between l and k space due to curvature (= comoving angular diameter distance to recombination / comoving radius to recombination)"]
            angular_rescaling: f64::default(),
            #[doc = "< number of wavenumber values"]
            q_size: u64::default(),
            #[doc = "< list of wavenumber values, q[index_q]"]
            q: &mut f64::default(),
            #[doc = "< list of wavenumber values for each requested mode, k[index_md][index_q]. In flat universes k=q. In non-flat universes q and k differ through q2 = k2 + K(1+m), where m=0,1,2 for scalar, vector, tensor. q should be used throughout the transfer module, excepted when interpolating or manipulating the source functions S(k,tau): for a given value of q this should be done in k(q)."]
            k: return_raw_dbl_ptr(),
            #[doc = "< index of the first q value using the flat rescaling approximation"]
            index_q_flat_approximation: ::std::os::raw::c_int::default(),
            #[doc = "< table of transfer functions for each mode, initial condition, type, multipole and wavenumber, with argument transfer[index_md][((index_ic * ptr->tt_size[index_md] + index_tt) * ptr->l_size[index_md] + index_l) * ptr->q_size + index_q]"]
            transfer: return_raw_dbl_ptr(),
            #[doc = "< flag regulating the amount of information sent to standard output (none if set to zero)"]
            transfer_verbose: ::std::os::raw::c_short::default(),
            #[doc = "< zone for writing error messages"]
            error_message:zero_padded_i8_array(""),
        }
    }
}