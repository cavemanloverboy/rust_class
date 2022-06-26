use serde::{Serialize, Deserialize};
use serde_yaml;
use crate::file_content;
use crate::utils::zero_padded_i8_array;

#[derive(Serialize, Deserialize, Debug)]
pub struct FileContent {

    output: Vec<Outputs>,

    lensing: Boolean,

    lcmb_rescale: u16,
    lcmb_tilt: f64,
    lcmb_pivot: f64,

    non_linear: Option<NonLinear>,
    ic: InitialConditions,
    modes: Vec<Modes>,
    l_max_scalars: usize,
    P_k_max_h_RSLASH_Mpc: f64,
    z_pk: Vec<f64>,
    h: f64,
    T_cmb: f64,

    omega_b: f64,
    omega_cdm: f64,
    omega_dcdmdr: f64,
    Omega_k: f64,
    Omega_fld: f64,
    Omega_scf: f64,
    N_ur: f64,

    YHe: HeliumFraction,
    recombination: RecombinationCode,

    z_reio: f64,

    reio_parametrization: ReionizationParameterization,
    reionization_exponent: f64,
    reionization_width: f64,
    helium_fullreio_redshift: f64,
    helium_fullreio_width: f64,

    DM_annihilation_cross_section: f64,
    DM_annihilation_mass: f64,
    DM_decay_fraction: f64,
    DM_decay_Gamma: f64,

    f_eff_type: FEffTypeEnum,
    chi_type: ChiTypeEnum,


    P_k_ini_SPACE_type: PrimordialSpectrum,
    k_pivot: f64,
    A_s: f64,
    n_s: f64,
    alpha_s: f64,


    sd_branching_approx: SdBranchingAppprox,
    sd_PCA_size: i16,
    sd_detector_name: SdDetectorName,

    include_SZ_effect: Boolean,

    overwrite_root: Boolean,
    write_background: Boolean,
    write_thermodynamics: Boolean,
    write_primordial: Boolean,
    write_parameters: Boolean,
    write_warnings: Boolean,

    input_verbose: u16,
    background_verbose: u16,
    thermodynamics_verbose: u16,
    perturbations_verbose: u16,
    transfer_verbose: u16,
    primordial_verbose: u16,
    harmonic_verbose: u16,
    fourier_verbose: u16,
    lensing_verbose: u16,
    output_verbose: u16,

}

#[derive(Serialize, Deserialize, Debug)]
pub enum Boolean {
    yes,
    yeap,
    no
}

#[derive(Serialize, Deserialize, Debug)]
pub enum Outputs {
    tCl,
    pCl,
    lCl,
    nCl,
    sCl,
    mPk,
    dTk,
    vTk,
    Sd,
}

#[derive(Serialize, Deserialize, Debug)]
pub enum NonLinear {
    halofit,
    HMCode,
}

#[derive(Serialize, Deserialize, Debug)]
pub enum InitialConditions {
    ad,
    bi,
    cdi,
    nid,
    nvi
}

#[derive(Serialize, Deserialize, Debug)]
pub enum Modes {
    s,
    v,
    t
}

#[derive(Serialize, Deserialize, Debug)]
pub enum PrimordialSpectrum {
    analytic_Pk,
    inflation_V,
    inflation_H,
    inflation_V_end,
    two_SPACE_scales,
    external_Pk,
}

#[derive(Serialize, Deserialize, Debug)]
pub enum HeliumFraction {
    BBN,
    Other(f64),
}

#[derive(Serialize, Deserialize, Debug)]
pub enum RecombinationCode {
    RECFAST,
    HyRec
}

#[derive(Serialize, Deserialize, Debug)]
pub enum ReionizationParameterization {
    reio_camb
}

#[derive(Serialize, Deserialize, Debug)]
pub enum FEffTypeEnum {
    on_the_spot
}

#[derive(Serialize, Deserialize, Debug)]
pub enum ChiTypeEnum {
    CK_2004,
}

#[derive(Serialize, Deserialize, Debug)]
pub enum SdBranchingAppprox {
    exact,
}


#[derive(Serialize, Deserialize, Debug)]
pub enum SdDetectorName {
    PIXIE,
}



impl Default for FileContent{

    fn default() -> FileContent {
        FileContent {
            output: vec![Outputs::tCl,Outputs::pCl,Outputs::lCl,Outputs::mPk],

            lensing: Boolean::yes,
        
            lcmb_rescale: 1,
            lcmb_tilt: 0.0,
            lcmb_pivot: 0.1,
        
            non_linear: None,
            ic: InitialConditions::ad,
            modes: vec![Modes::s],
            l_max_scalars: 2500,
            P_k_max_h_RSLASH_Mpc: 1.0,
            z_pk: vec![0.0],
            h: 0.67810,
            T_cmb: 2.7255,
        
            omega_b: 0.02238280,
            omega_cdm: 0.1201075,
            omega_dcdmdr: 0.0,
            Omega_k: 0.0,
            Omega_fld: 0.0,
            Omega_scf: 0.0,
            N_ur: 3.044,
        
            YHe: HeliumFraction::BBN,
            recombination: RecombinationCode::HyRec,
        
            z_reio: 7.6711,
        
            reio_parametrization: ReionizationParameterization::reio_camb,
            reionization_exponent: 1.5,
            reionization_width: 0.5,
            helium_fullreio_redshift: 3.5,
            helium_fullreio_width: 0.5,
        
            DM_annihilation_cross_section: 0.0,
            DM_annihilation_mass: 0.0,
            DM_decay_fraction: 0.0,
            DM_decay_Gamma: 0.0,
        
            f_eff_type: FEffTypeEnum::on_the_spot,
            chi_type: ChiTypeEnum::CK_2004,
        
        
            P_k_ini_SPACE_type: PrimordialSpectrum::analytic_Pk,
            k_pivot: 0.05,
            A_s: 2.100549e-09,
            n_s: 0.9660499,
            alpha_s: 0.0,
        
        
            sd_branching_approx: SdBranchingAppprox::exact,
            sd_PCA_size: 2,
            sd_detector_name: SdDetectorName::PIXIE,
        
            include_SZ_effect: Boolean::no,
        
            overwrite_root: Boolean::no,
            write_background: Boolean::no,
            write_thermodynamics: Boolean::no,
            write_primordial: Boolean::no,
            write_parameters: Boolean::yeap,
            write_warnings: Boolean::yes,
        
            input_verbose: 1,
            background_verbose: 1,
            thermodynamics_verbose: 1,
            perturbations_verbose: 1,
            transfer_verbose: 1,
            primordial_verbose: 1,
            harmonic_verbose: 1,
            fourier_verbose: 1,
            lensing_verbose: 1,
            output_verbose: 1,
        }
    }
}

impl FileContent {
    pub fn to_file_content(&self) -> file_content {
        
        // Convert struct to string
        let contents_of_file: String = serde_yaml::to_string(&self)
            .unwrap()
            .replace("\n  -", ",")
            .replace("_SPACE_", " ")
            .replace("Other(", "")
            .replace(")", "")
            .replace("[","")
            .replace("]","")
            .replace("_RSLASH_","/")
            .replace("None", "")
            .replace("---\n","")
            .replace(":"," =")
            .replace("\"","")
            .replace("=,","=")
            .replace("~", "");

        println!("{}",contents_of_file);

        file_content {
            filename: &mut i8::default(),
            size: contents_of_file.as_bytes().len() as i32,
            #[doc = "< list of (size) names"]
            name: {
                let arrays: Vec<&str> = contents_of_file.lines().map(|line| { println!("some_line: {}", line); *line.split(" = ").collect::<Vec<&str>>().get(0).unwrap() }).collect();
                &mut zero_padded_i8_array(arrays.join("").as_str())
            },
            #[doc = "< list of (size) values"]
            value: {
                let arrays: Vec<&str> = contents_of_file.lines().map(|line| { println!("{}", line); *line.split(" = ").collect::<Vec<&str>>().get(1).unwrap() }).collect();
                &mut zero_padded_i8_array(arrays.join(" ").as_str())
            },
            #[doc = "< set to _TRUE_ if this parameter is effectively read"]
            read: &mut 0,
        }

    }
}

//