use bytemuck::cast;
use std::io::Write;

/// Generates zero-padded i8 arrays from strings. Required by bindgen-generated bindings.
pub fn zero_padded_i8_array<const L: usize>(string: &str) -> [i8; L]
where
    [u8; L]: bytemuck::Pod,
    [i8; L]: bytemuck::Pod,
{

    let mut buf = [0; L];
    write!(buf.as_mut_slice(), "{}", string);
    cast::<[u8; L], [i8; L]>(buf)
}

pub fn return_raw_dbl_ptr<T: Default>() -> *mut *mut T {
    let mut value = T::default();
    let mut raw_ptr: *mut T = &mut value;
    let raw_dbl_ptr: *mut *mut T = &mut raw_ptr;
    raw_dbl_ptr
}

pub fn return_raw_tpl_ptr<T: Default>() -> *mut *mut *mut T {
    let mut value = T::default();
    let mut raw_ptr: *mut T = &mut value;
    let mut raw_dbl_ptr: *mut *mut T = &mut raw_ptr;
    let raw_tpl_ptr: *mut *mut *mut T = &mut raw_dbl_ptr;
    raw_tpl_ptr
}