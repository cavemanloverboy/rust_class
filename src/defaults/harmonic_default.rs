use crate::harmonic;
use crate::utils::*;


impl Default for harmonic {
    fn default() -> harmonic {
        harmonic {
            #[doc = "< sets the number of cross-correlation spectra"]
            #[doc = "that you want to calculate: 0 means only"]
            #[doc = "auto-correlation, 1 means only adjacent bins,"]
            #[doc = "and number of bins minus one means all"]
            #[doc = "correlations"]
            non_diag: ::std::os::raw::c_int::default(),
            #[doc = "< number of modes (scalar, tensor, ...) included in computation"]
            md_size: ::std::os::raw::c_int::default(),
            #[doc = "< index for scalar modes"]
            index_md_scalars: ::std::os::raw::c_int::default(),
            #[doc = "< for a given mode, ic_size[index_md] = number of initial conditions included in computation"]
            ic_size: &mut ::std::os::raw::c_int::default(),
            #[doc = "< for a given mode, ic_ic_size[index_md] = number of pairs of (index_ic1, index_ic2) with index_ic2 >= index_ic1; this number is just N(N+1)/2  where N = ic_size[index_md]"]
            ic_ic_size: &mut ::std::os::raw::c_int::default(),
            #[doc = "< for a given mode, is_non_zero[index_md][index_ic1_ic2] is set to true if the pair of initial conditions (index_ic1, index_ic2) are statistically correlated, or to false if they are uncorrelated"]
            is_non_zero: return_raw_dbl_ptr(),
            #[doc = "< do we want \\f$ C_l^{TT}\\f$? (T = temperature)"]
            has_tt: ::std::os::raw::c_int::default(),
            #[doc = "< do we want \\f$ C_l^{EE}\\f$? (E = E-polarization)"]
            has_ee: ::std::os::raw::c_int::default(),
            #[doc = "< do we want \\f$ C_l^{TE}\\f$?"]
            has_te: ::std::os::raw::c_int::default(),
            #[doc = "< do we want \\f$ C_l^{BB}\\f$? (B = B-polarization)"]
            has_bb: ::std::os::raw::c_int::default(),
            #[doc = "< do we want \\f$ C_l^{\\phi\\phi}\\f$? (\\f$ \\phi \\f$ = CMB lensing potential)"]
            has_pp: ::std::os::raw::c_int::default(),
            #[doc = "< do we want \\f$ C_l^{T\\phi}\\f$?"]
            has_tp: ::std::os::raw::c_int::default(),
            #[doc = "< do we want \\f$ C_l^{E\\phi}\\f$?"]
            has_ep: ::std::os::raw::c_int::default(),
            #[doc = "< do we want \\f$ C_l^{dd}\\f$? (d = density)"]
            has_dd: ::std::os::raw::c_int::default(),
            #[doc = "< do we want \\f$ C_l^{Td}\\f$?"]
            has_td: ::std::os::raw::c_int::default(),
            #[doc = "< do we want \\f$ C_l^{\\phi d}\\f$?"]
            has_pd: ::std::os::raw::c_int::default(),
            #[doc = "< do we want \\f$ C_l^{ll}\\f$? (l = galaxy lensing potential)"]
            has_ll: ::std::os::raw::c_int::default(),
            #[doc = "< do we want \\f$ C_l^{Tl}\\f$?"]
            has_tl: ::std::os::raw::c_int::default(),
            #[doc = "< do we want \\f$ C_l^{dl}\\f$?"]
            has_dl: ::std::os::raw::c_int::default(),
            #[doc = "< index for type \\f$ C_l^{TT} \\f$"]
            index_ct_tt: ::std::os::raw::c_int::default(),
            #[doc = "< index for type \\f$ C_l^{EE} \\f$"]
            index_ct_ee: ::std::os::raw::c_int::default(),
            #[doc = "< index for type \\f$ C_l^{TE} \\f$"]
            index_ct_te: ::std::os::raw::c_int::default(),
            #[doc = "< index for type \\f$ C_l^{BB} \\f$"]
            index_ct_bb: ::std::os::raw::c_int::default(),
            #[doc = "< index for type \\f$ C_l^{\\phi\\phi} \\f$"]
            index_ct_pp: ::std::os::raw::c_int::default(),
            #[doc = "< index for type \\f$ C_l^{T\\phi} \\f$"]
            index_ct_tp: ::std::os::raw::c_int::default(),
            #[doc = "< index for type \\f$ C_l^{E\\phi} \\f$"]
            index_ct_ep: ::std::os::raw::c_int::default(),
            #[doc = "< first index for type \\f$ C_l^{dd} \\f$((d_size*d_size-(d_size-non_diag)*(d_size-non_diag-1)/2) values)"]
            index_ct_dd: ::std::os::raw::c_int::default(),
            #[doc = "< first index for type \\f$ C_l^{Td} \\f$(d_size values)"]
            index_ct_td: ::std::os::raw::c_int::default(),
            #[doc = "< first index for type \\f$ C_l^{pd} \\f$(d_size values)"]
            index_ct_pd: ::std::os::raw::c_int::default(),
            #[doc = "< first index for type \\f$ C_l^{ll} \\f$((d_size*d_size-(d_size-non_diag)*(d_size-non_diag-1)/2) values)"]
            index_ct_ll: ::std::os::raw::c_int::default(),
            #[doc = "< first index for type \\f$ C_l^{Tl} \\f$(d_size values)"]
            index_ct_tl: ::std::os::raw::c_int::default(),
            #[doc = "< first index for type \\f$ C_l^{dl} \\f$(d_size values)"]
            index_ct_dl: ::std::os::raw::c_int::default(),
            #[doc = "< number of bins for which density Cl's are computed"]
            d_size: ::std::os::raw::c_int::default(),
            #[doc = "< number of \\f$ C_l \\f$ types requested"]
            ct_size: ::std::os::raw::c_int::default(),
            #[doc = "< number of multipole values for each requested mode, l_size[index_md]"]
            l_size: &mut ::std::os::raw::c_int::default(),
            #[doc = "< greatest of all l_size[index_md]"]
            l_size_max: ::std::os::raw::c_int::default(),
            #[doc = "< list of multipole values l[index_l]"]
            l: &mut f64::default(),
            #[doc = "< last multipole (given as an input) at which"]
            #[doc = "we want to output \\f$ C_l\\f$'s for a given mode and type;"]
            #[doc = "l[index_md][l_size[index_md]-1] can be larger"]
            #[doc = "than l_max[index_md], in order to ensure a"]
            #[doc = "better interpolation with no boundary effects"]
            l_max_ct: return_raw_dbl_ptr(),
            #[doc = "< last multipole (given as an input) at which"]
            #[doc = "we want to output \\f$ C_l\\f$'s for a given mode (maximized over types);"]
            #[doc = "l[index_md][l_size[index_md]-1] can be larger"]
            #[doc = "than l_max[index_md], in order to ensure a"]
            #[doc = "better interpolation with no boundary effects"]
            l_max: &mut ::std::os::raw::c_int::default(),
            #[doc = "< last multipole (given as an input) at which"]
            #[doc = "we want to output \\f$ C_l\\f$'s (maximized over modes and types);"]
            #[doc = "l[index_md][l_size[index_md]-1] can be larger"]
            #[doc = "than l_max[index_md], in order to ensure a"]
            #[doc = "better interpolation with no boundary effects"]
            l_max_tot: ::std::os::raw::c_int::default(),
            #[doc = "< table of anisotropy spectra for each mode, multipole, pair of initial conditions and types, cl[index_md][(index_l * phr->ic_ic_size[index_md] + index_ic1_ic2) * phr->ct_size + index_ct]"]
            cl: return_raw_dbl_ptr(),
            #[doc = "< second derivatives of previous table with respect to l, in view of spline interpolation"]
            ddcl: return_raw_dbl_ptr(),
            #[doc = "< a pointer to the fourier structure is"]
            #[doc = "stored in the harmonic structure. This odd,"]
            #[doc = "unusual and unelegant feature has been"]
            #[doc = "introduced in v2.8 in order to keep in use"]
            #[doc = "some deprecated functions harmonic_pk_...()"]
            #[doc = "that are now pointing at new function"]
            #[doc = "fourier_pk_...(). In the future, if the"]
            #[doc = "deprecated functions are removed, it will"]
            #[doc = "be possible to remove also this pointer."]
            pfo: &mut crate::fourier::default(),
            #[doc = "< flag regulating the amount of information sent to standard output (none if set to zero)"]
            harmonic_verbose: ::std::os::raw::c_short::default(),
            #[doc = "< zone for writing error messages"]
            error_message: zero_padded_i8_array(""),
        }
    }
}