use crate::distortions;
use crate::noninjection;
use crate::utils::*;
use crate::reio_approx;

impl Default for distortions {
    fn default() -> Self {
        distortions {
            #[doc = "< Which approximation to use for the branching ratios?"]
            sd_branching_approx: ::std::os::raw::c_int::default(),
            #[doc = "< Number of PCA components for the calculation of residual distortions"]
            sd_PCA_size: ::std::os::raw::c_int::default(),
            #[doc = "< Name of detector list file"]
            sd_detector_file_name: [i8::default(); 612],
            #[doc = "< Name of detector"]
            sd_detector_name: [i8::default(); 100],
            #[doc = "< Minimum frequency of chosen detector"]
            sd_detector_nu_min: f64::default(),
            #[doc = "< Maximum frequency of chosen detector"]
            sd_detector_nu_max: f64::default(),
            #[doc = "< Bin size of chosen detector"]
            sd_detector_nu_delta: f64::default(),
            #[doc = "< Number of frequency bins of chosen detector"]
            sd_detector_bin_number: ::std::os::raw::c_int::default(),
            #[doc = "< Sensitivity of the chosen detector"]
            sd_detector_delta_Ic: f64::default(),
            #[doc = "< Calculation method for Sunyaev Zeldovich contributions from re-ionization"]
            sd_reio_type: reio_approx::default(),
            #[doc = "< Possible additional y contribution (manually) to the SD signal"]
            sd_add_y: f64::default(),
            #[doc = "< Possible additional mu contribution (manually) to the SD signal"]
            sd_add_mu: f64::default(),
            #[doc = "< Redshift of the transition of mu to y era"]
            z_muy: f64::default(),
            #[doc = "< Redshift of the transition from thermal shift to mu era"]
            z_th: f64::default(),
            #[doc = "< Minimum redshift"]
            z_min: f64::default(),
            #[doc = "< Maximum redshift"]
            z_max: f64::default(),
            #[doc = "< Lenght of redshift array"]
            z_size: ::std::os::raw::c_int::default(),
            #[doc = "< Redshift intervals"]
            z_delta: f64::default(),
            #[doc = "< Redshift list z[index_z] = list of values"]
            z: &mut f64::default(),
            #[doc = "< Weights for integration over z"]
            z_weights: &mut f64::default(),
            #[doc = "< Minimum dimentionless frequency"]
            x_min: f64::default(),
            #[doc = "< Maximum dimentionless frequency"]
            x_max: f64::default(),
            #[doc = "< dimentionless frequency intervals"]
            x_delta: f64::default(),
            #[doc = "< Lenght of dimentionless frequency array"]
            x_size: ::std::os::raw::c_int::default(),
            #[doc = "< Dimensionless frequency x[index_x] = list of values"]
            x: &mut f64::default(),
            #[doc = "< Weights for integration over x"]
            x_weights: &mut f64::default(),
            #[doc = "< Conversion factor nu[GHz] = x_to_nu * x"]
            x_to_nu: f64::default(),
            #[doc = "< Conversion from unitless DI to DI[10^26 W m^-2 Hz^-1 sr^-1]"]
            DI_units: f64::default(),
            #[doc = "< Full path of detector noise file"]
            sd_detector_noise_file: [i8::default(); 868],
            #[doc = "< Full path of PCA generator file"]
            sd_PCA_file_generator: [i8::default(); 612],
            #[doc = "< Full path of detector list file"]
            sd_detector_list_file: [i8::default(); 612],
            #[doc = "< Branching ratios br_table[index_type][index_z]"]
            br_table: return_raw_dbl_ptr::<f64>(),
            #[doc = "< Spectral Distortion parameters (g,mu,y,r) sd_parameter_table[index_type]"]
            sd_parameter_table: &mut f64::default(),
            #[doc = "< Spectral Distortion shapes (G,M,Y,R) sd_shape_table[index_type][index_x]"]
            sd_shape_table: return_raw_dbl_ptr::<f64>(),
            #[doc = "< Spectral Distortion Intensities (final deltaI seperated by component) sd_table[index_type][index_x]"]
            sd_table: return_raw_dbl_ptr::<f64>(),
            #[doc = "< temperature shift/g type distortion"]
            index_type_g: ::std::os::raw::c_int::default(),
            #[doc = "< mu type distortion"]
            index_type_mu: ::std::os::raw::c_int::default(),
            #[doc = "< y type distortion"]
            index_type_y: ::std::os::raw::c_int::default(),
            #[doc = "< PCA type distortion (first index)"]
            index_type_PCA: ::std::os::raw::c_int::default(),
            #[doc = "< Number of total components for the type array"]
            type_size: ::std::os::raw::c_int::default(),
            epsilon: f64::default(),
            dQrho_dz_tot: &mut f64::default(),
            Drho_over_rho: f64::default(),
            #[doc = "< DI[index_x] = list of values"]
            DI: &mut f64::default(),
            #[doc = "< Redshift array for reading from file br_exact_z[index_z]"]
            br_exact_z: &mut f64::default(),
            #[doc = "< Number of redshift values for reading from file"]
            br_exact_Nz: ::std::os::raw::c_int::default(),
            #[doc = "< temperature shift/g distortion branching ratio f_g_exact[index_z]"]
            f_g_exact: &mut f64::default(),
            #[doc = "< second derivative of the above ddf_g_exact[index_z]"]
            ddf_g_exact: &mut f64::default(),
            #[doc = "< y distortion branching ratio f_y_exact[index_z]"]
            f_y_exact: &mut f64::default(),
            #[doc = "< second derivative of the above ddf_y_exact[index_z]"]
            ddf_y_exact: &mut f64::default(),
            #[doc = "< mu distortion shape branching ratio f_mu_exact[index_z]"]
            f_mu_exact: &mut f64::default(),
            #[doc = "< second derivative of the above ddf_mu_exact[index_z]"]
            ddf_mu_exact: &mut f64::default(),
            #[doc = "< PCA component E branching ratio for reading from file E_vec[index_e*br_exact_Nz+index_z] with index_e=[1..8]"]
            E_vec: &mut f64::default(),
            #[doc = "< second derivative of the above ddE_vec[index_e*br_exact_Nz+index_z]"]
            ddE_vec: &mut f64::default(),
            #[doc = "< number of PCA component E branching ratios"]
            E_vec_size: ::std::os::raw::c_int::default(),
            #[doc = "< Frquency array for reading from file PCA_nu[index_nu]"]
            PCA_nu: &mut f64::default(),
            #[doc = "< Number of frequency values for reading from file"]
            PCA_Nnu: ::std::os::raw::c_int::default(),
            #[doc = "< temperature shift/g distortion shape PCA_G_T[index_nu]"]
            PCA_G_T: &mut f64::default(),
            #[doc = "< second derivative of the above ddPCA_G_T[index_nu]"]
            ddPCA_G_T: &mut f64::default(),
            #[doc = "< y distortion shape PCA_Y_SZ[index_nu]"]
            PCA_Y_SZ: &mut f64::default(),
            #[doc = "< second derivative of the above ddPCA_Y_SZ[index_nu]"]
            ddPCA_Y_SZ: &mut f64::default(),
            #[doc = "< mu distortion shape PCA_M_mu[index_nu]"]
            PCA_M_mu: &mut f64::default(),
            #[doc = "< second derivative of the above ddPCA_M_mu[index_nu]"]
            ddPCA_M_mu: &mut f64::default(),
            #[doc = "< PCA component S shape for reading from file S_vec[index_s*S_vec_size+index_x] with index_s=[1..8]"]
            S_vec: &mut f64::default(),
            #[doc = "< second derivative of the above ddS_vec[index_s*S_vec_size+index_x]"]
            ddS_vec: &mut f64::default(),
            #[doc = "< number of PCA component S spectral shapes"]
            S_vec_size: ::std::os::raw::c_int::default(),
            #[doc = "< delta_Ic[index_x] for detectors with given sensitivity in each bin"]
            delta_Ic_array: &mut f64::default(),
            #[doc = "< do we need to compute spectral distortions?"]
            has_distortions: ::std::os::raw::c_int::default(),
            #[doc = "< does the user specify their own detector?"]
            has_user_defined_detector: ::std::os::raw::c_int::default(),
            #[doc = "< does the user specify the name of their detector?"]
            has_user_defined_name: ::std::os::raw::c_int::default(),
            #[doc = "< do we have a file for the detector specification?"]
            has_detector_file: ::std::os::raw::c_int::default(),
            #[doc = "< do we include the SZ effect?"]
            has_SZ_effect: ::std::os::raw::c_int::default(),
            #[doc = "< shall we only take exotic injection contributions?"]
            include_only_exotic: ::std::os::raw::c_int::default(),
            #[doc = "< shall we include the g distortion in the total distortion ?"]
            include_g_distortion: ::std::os::raw::c_int::default(),
            #[doc = "< do we have terms that are not injected (like dissipation of acoustic waves)?"]
            has_noninjected: ::std::os::raw::c_int::default(),
            #[doc = "< noninjection file structure"]
            ni: noninjection::default(),
            #[doc = "< flag regulating the amount of information sent to standard output (none if set to zero)"]
            distortions_verbose: ::std::os::raw::c_short::default(),
            #[doc = "< zone for writing error messages"]
            error_message: zero_padded_i8_array::<2048>(""),
        }
    }
}


impl Default for noninjection {
    fn default() -> noninjection {
        noninjection {
            k_min: f64::default(),
            k_max: f64::default(),
            k_size: ::std::os::raw::c_int::default(),
            k: &mut f64::default(),
            k_weights: &mut f64::default(),
            pk_primordial_k: &mut f64::default(),
            integrand_approx: &mut f64::default(),
            z_table_coarse: &mut f64::default(),
            z_size_coarse: ::std::os::raw::c_int::default(),
            logz_max: f64::default(),
            noninjection_table: &mut f64::default(),
            ddnoninjection_table: &mut f64::default(),
            z_size: ::std::os::raw::c_int::default(),
            z_table: &mut f64::default(),
            last_index_z: ::std::os::raw::c_int::default(),
            f_nu_wkb: f64::default(),
            dkD_dz: f64::default(),
            kD: f64::default(),
            H0: f64::default(),
            T_g0: f64::default(),
            Omega0_b: f64::default(),
            Omega0_cdm: f64::default(),
            rho0_cdm: f64::default(),
            fHe: f64::default(),
            N_e0: f64::default(),
            H: f64::default(),
            a: f64::default(),
            rho_g: f64::default(),
            heat_capacity: f64::default(),
            nH: f64::default(),
            T_b: f64::default(),
            T_g: f64::default(),
            x_e: f64::default(),
            photon_dep_table: &mut f64::default(),
            error_message: zero_padded_i8_array::<2048>(""),
        }
    }
}