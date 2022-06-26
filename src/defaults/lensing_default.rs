use crate::lensing;
use crate::utils::*;

impl Default for lensing {
    fn default() -> lensing {
        lensing {
            #[doc = "< do we need to compute lensed \\f$ C_l\\f$'s at all ?"]
            has_lensed_cls: ::std::os::raw::c_short::default(),
            #[doc = "< do we want lensed \\f$ C_l^{TT}\\f$? (T = temperature)"]
            has_tt: ::std::os::raw::c_int::default(),
            #[doc = "< do we want lensed \\f$ C_l^{EE}\\f$? (E = E-polarization)"]
            has_ee: ::std::os::raw::c_int::default(),
            #[doc = "< do we want lensed \\f$ C_l^{TE}\\f$?"]
            has_te: ::std::os::raw::c_int::default(),
            #[doc = "< do we want \\f$ C_l^{BB}\\f$? (B = B-polarization)"]
            has_bb: ::std::os::raw::c_int::default(),
            #[doc = "< do we want \\f$ C_l^{\\phi\\phi}\\f$? (\\f$ \\phi \\f$ = CMB lensing potential)"]
            has_pp: ::std::os::raw::c_int::default(),
            #[doc = "< do we want \\f$ C_l^{T\\phi}\\f$?"]
            has_tp: ::std::os::raw::c_int::default(),
            #[doc = "< do we want \\f$ C_l^{dd}\\f$? (d = matter density)"]
            has_dd: ::std::os::raw::c_int::default(),
            #[doc = "< do we want \\f$ C_l^{Td}\\f$?"]
            has_td: ::std::os::raw::c_int::default(),
            #[doc = "< do we want \\f$ C_l^{ll}\\f$? (l = lensing potential)"]
            has_ll: ::std::os::raw::c_int::default(),
            #[doc = "< do we want \\f$ C_l^{Tl}\\f$?"]
            has_tl: ::std::os::raw::c_int::default(),
            #[doc = "< index for type \\f$ C_l^{TT} \\f$"]
            index_lt_tt: ::std::os::raw::c_int::default(),
            #[doc = "< index for type \\f$ C_l^{EE} \\f$"]
            index_lt_ee: ::std::os::raw::c_int::default(),
            #[doc = "< index for type \\f$ C_l^{TE} \\f$"]
            index_lt_te: ::std::os::raw::c_int::default(),
            #[doc = "< index for type \\f$ C_l^{BB} \\f$"]
            index_lt_bb: ::std::os::raw::c_int::default(),
            #[doc = "< index for type \\f$ C_l^{\\phi\\phi} \\f$"]
            index_lt_pp: ::std::os::raw::c_int::default(),
            #[doc = "< index for type \\f$ C_l^{T\\phi} \\f$"]
            index_lt_tp: ::std::os::raw::c_int::default(),
            #[doc = "< index for type \\f$ C_l^{dd} \\f$"]
            index_lt_dd: ::std::os::raw::c_int::default(),
            #[doc = "< index for type \\f$ C_l^{Td} \\f$"]
            index_lt_td: ::std::os::raw::c_int::default(),
            #[doc = "< index for type \\f$ C_l^{dd} \\f$"]
            index_lt_ll: ::std::os::raw::c_int::default(),
            #[doc = "< index for type \\f$ C_l^{Td} \\f$"]
            index_lt_tl: ::std::os::raw::c_int::default(),
            #[doc = "< number of \\f$ C_l\\f$ types requested"]
            lt_size: ::std::os::raw::c_int::default(),
            #[doc = "< last multipole in all calculations (same as in harmonic module)"]
            l_unlensed_max: ::std::os::raw::c_int::default(),
            #[doc = "< last multipole at which lensed spectra are computed"]
            l_lensed_max: ::std::os::raw::c_int::default(),
            #[doc = "< number of l values"]
            l_size: ::std::os::raw::c_int::default(),
            #[doc = "< last multipole (given as an input) at which"]
            #[doc = "we want to output \\f$ C_l \\f$'s for a given mode and type"]
            l_max_lt: &mut ::std::os::raw::c_int::default(),
            #[doc = "< table of multipole values l[index_l]"]
            l: &mut f64::default(),
            #[doc = "< table of anisotropy spectra for each"]
            #[doc = "multipole and types,"]
            #[doc = "cl[index_l * ple->lt_size + index_lt]"]
            cl_lens: &mut f64::default(),
            #[doc = "< second derivatives for interpolation"]
            ddcl_lens: &mut f64::default(),
            #[doc = "< flag regulating the amount of information sent to standard output (none if set to zero)"]
            lensing_verbose: ::std::os::raw::c_short::default(),
            #[doc = "< zone for writing error messages"]
            error_message: zero_padded_i8_array(""),
        }
    }
}