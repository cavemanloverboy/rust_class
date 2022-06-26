#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]

include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

pub mod defaults;
pub mod utils;
pub mod file_content_rust;
use file_content_rust::FileContent;

pub fn calculate_spectra() {

    // // struct precision pr;        /* for precision parameters */
    // // struct background ba;       /* for cosmological background */
    // // struct thermodynamics th;           /* for thermodynamics */
    // // struct perturbations pt;         /* for source functions */
    // // struct primordial pm;       /* for primordial spectra */
    // // struct fourier fo;        /* for non-linear spectra */
    // // struct transfer tr;        /* for transfer functions */
    // // struct harmonic hr;          /* for output spectra */
    // // struct lensing le;          /* for lensed spectra */
    // struct distortions sd;      /* for spectral distortions */
    // struct output op;           /* for output files */

    // Define precision struct raw ptr
    let mut precision_struct = precision {
        ..
        precision::default()
    };

    // Define background struct
    let mut background_struct = background {
        ..
        background::default()
    };

    // Define thermodynamics struct
    let mut thermodynamics_struct = thermodynamics {
        ..
        thermodynamics::default()
    };

    // define perturbations struct
    let mut perturbations_struct = perturbations {
        ..
        perturbations::default()
    };

    // Define primordial struct
    let mut primordial_struct = primordial {
        ..
        primordial::default()
    };

    // Define fourier struct
    let mut fourier_struct = fourier {
        ..
        fourier::default()
    };

    // Define transfer struct
    let mut transfer_struct = transfer {
        ..
        transfer::default()
    };

    // Define harmonic struct
    let mut harmonic_struct = harmonic {
        ..
        harmonic::default()
    };

    // Define lensing struct
    let mut lensing_struct = lensing {
        ..
        lensing::default()
    };

    // Define distortions struct
    let mut distortions_struct = distortions {
        ..
        distortions::default()
    };

    // Define output struct
    let mut output_struct = output {
        ..
        output::default()
    };

    // Define file content in rust
    let my_file_content = FileContent::default();

    let mut my_file_content = my_file_content.to_file_content();

    // Define error message
    let mut errmsg = i8::default();

    unsafe {
        {
            // input init scope without find_file

            println!("input_read_from_file");
            input_read_from_file(
                &mut my_file_content,
                &mut precision_struct,
                &mut background_struct,
                &mut thermodynamics_struct,
                &mut perturbations_struct,
                &mut transfer_struct,
                &mut primordial_struct,
                &mut harmonic_struct,
                &mut fourier_struct,
                &mut lensing_struct,
                &mut distortions_struct,
                &mut output_struct,
                &mut errmsg
            );
            
            println!("parser_free");
            parser_free(&mut my_file_content);
        }


        println!("background_init");
        background_init(
            &mut precision_struct,
            &mut background_struct,
        );

        thermodynamics_init(
            &mut precision_struct,
            &mut background_struct,
            &mut thermodynamics_struct,
        );

        perturbations_init(
            &mut precision_struct,
            &mut background_struct,
            &mut thermodynamics_struct,
            &mut perturbations_struct,
        );

        primordial_init(
            &mut precision_struct,
            &mut perturbations_struct,
            &mut primordial_struct,
        );

        fourier_init(
            &mut precision_struct,
            &mut background_struct,
            &mut thermodynamics_struct,
            &mut perturbations_struct,
            &mut primordial_struct,
            &mut fourier_struct,
        );

        transfer_init(
            &mut precision_struct,
            &mut background_struct,
            &mut thermodynamics_struct,
            &mut perturbations_struct,
            &mut fourier_struct,
            &mut transfer_struct,
        );

        harmonic_init(
            &mut precision_struct,
            &mut background_struct,
            &mut perturbations_struct,
            &mut primordial_struct,
            &mut fourier_struct,
            &mut transfer_struct,
            &mut harmonic_struct,
        );

        lensing_init(
            &mut precision_struct,
            &mut perturbations_struct,
            &mut harmonic_struct,
            &mut fourier_struct,
            &mut lensing_struct,
        );

        distortions_init(
            &mut precision_struct,
            &mut background_struct,
            &mut thermodynamics_struct,
            &mut perturbations_struct,
            &mut primordial_struct,
            &mut distortions_struct,
        );

        output_init(
            &mut background_struct,
            &mut thermodynamics_struct,
            &mut perturbations_struct,
            &mut primordial_struct,
            &mut transfer_struct,
            &mut harmonic_struct,
            &mut fourier_struct,
            &mut lensing_struct,
            &mut distortions_struct,
            &mut output_struct,
        );
    }
}

// #[cfg(test)]
// mod tests {
//     #[test]
//     fn it_works() {
//         let result = 2 + 2;
//         assert_eq!(result, 4);
//     }

//     #[test]
//     fn test_package() {
//         use crate::calculate_spectra;

//         calculate_spectra()
//     }
// }


