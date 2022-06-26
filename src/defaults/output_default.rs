use crate::output;
use crate::utils::*;

impl Default for output {
    fn default() -> output {
        output {
            #[doc = "< root for all file names"]
            root: [::std::os::raw::c_char::default(); 224usize],
            #[doc = "< number of redshift at which P(k,z) and T_i(k,z) should be written"]
            z_pk_num: ::std::os::raw::c_int::default(),
            #[doc = "< value(s) of redshift at which P(k,z) and T_i(k,z) should be written"]
            z_pk: [f64::default(); 100usize],
            #[doc = "< flag stating whether we should write a header in output files"]
            write_header: ::std::os::raw::c_short::default(),
            #[doc = "< which format for output files (definitions, order of columns, etc.)"]
            output_format: u32::default(),
            #[doc = "< flag for outputing background evolution in file"]
            write_background: ::std::os::raw::c_short::default(),
            #[doc = "< flag for outputing thermodynamical evolution in file"]
            write_thermodynamics: ::std::os::raw::c_short::default(),
            #[doc = "< flag for outputing perturbations of selected wavenumber(s) in file(s)"]
            write_perturbations: ::std::os::raw::c_short::default(),
            #[doc = "< flag for outputing scalar/tensor primordial spectra in files"]
            write_primordial: ::std::os::raw::c_short::default(),
            #[doc = "< flag for outputing exotic energy injection/deposition in files"]
            write_exotic_injection: ::std::os::raw::c_short::default(),
            #[doc = "< flag for outputing non-injected contributions in files"]
            write_noninjection: ::std::os::raw::c_short::default(),
            #[doc = "< flag for outputing spectral distortions in files"]
            write_distortions: ::std::os::raw::c_short::default(),
            #[doc = "< flag regulating the amount of information sent to standard output (none if set to zero)"]
            output_verbose: ::std::os::raw::c_short::default(),
            #[doc = "< zone for writing error messages"]
            error_message: zero_padded_i8_array(""),
        }
    }
}